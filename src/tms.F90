module tms

  use ace_header 
  use cross_section,  only: cdintegral
  use error,          only: fatal_error, warning 
  use global
  use output,         only: write_message
  use particle_header 
  use random_lcg,     only: prn
  use string,         only: to_str

  implicit none


contains   

!===============================================================================
! TMS_INIT makes all the initializations needed for transport with TMS 
! (find maximum temperature for each nuclide, generate majorants)
!===============================================================================

  subroutine tms_init()
    integer :: n       ! number of TMS materials in problem

    n = tms_set_maximum_temperatures() 
    
    if( n == 0 ) then
       ! if there are no TMS materials in the problem, init nothing
       return
    end if
   
    ! Check for ptables. TMS method is not yet able to handle URR and is 
    ! incompatible with the ptable treatment. To avoid misusage of TMS, 
    ! force ptables = false. 

    ! NOTE: If you know what you are doing, this limitation can be removed with
    ! the consequence that ptables will work correctly in materials in which
    ! TMS is not used, but the materials with TMS will be handled using 
    ! infinitely dilute cross sections. (not tested)

    if ( urr_ptables_on ) then
       message="TMS cannot be used with URR ptable treatment. Please " //& 
            "add <ptables>false</ptables> in settings.xml"
       call fatal_error()
    end if

    call tms_set_thresholds()
    
    ! Print message (this is unimportant and can be removed )
    
    message="Initializing TMS for all nuclides in " // trim(to_str(n)) &
         // " materials..." 
    call write_message(8)

    call tms_calculate_majorants()

  end subroutine tms_init

!===============================================================================
! TMS_SET_MAXIMUM_TEMPERATURES finds the maximum temperature of each nuclide
! for which the TMS treatment is used. Also performs some important checks.
! Returns the number of TMS materials. 
!===============================================================================

  function tms_set_maximum_temperatures() result(n)
    integer :: i           ! loop index for materials
    integer :: u           ! loop index for nuclides 
    integer :: n           ! number of TMS materials 
    integer :: n_different ! number of nuclides for which 
    integer :: nuclide_i   ! index of nuclide
    type(Material),    pointer :: mat => null()

    n = 0;
    
    ! loop over all materials
    
    do i = 1, n_materials

       mat => materials(i) 
             
       if ( mat % tmstemp >= 0.0 ) then

          n = n + 1          

          ! Check that no S(a,b) tables are associated with
          ! the material. TMS is currently (and will be in the near future)
          ! incompatible with S(a,b)

          if ( mat % n_sab > 0) then
             message="S(a,b) tables associated with TMS material (id = "// &
                  trim(to_str(mat % id)) // " )"
             call fatal_error()
          end if          

          ! loop over nuclides in material and set maximum 
          
          n_different = 0

          do u = 1, mat % n_nuclides
                          
             nuclide_i = mat % nuclide(u)
             
             ! If material temperature is larger than max T of nuclide,
             ! update max_kT
             
             if ( nuclides(nuclide_i) % max_kT < mat % tmstemp ) then
                nuclides(nuclide_i) % max_kT = mat % tmstemp
             end if      
             
             if ( mat % tmstemp - nuclides(nuclide_i) % kT & 
                  > 0.001*K_BOLTZMANN) then
                n_different = n_different + 1 
             end if

             
          end do

          ! if there are no nuclides in this material for which the temperature 
          ! difference between the cross sections and tmstemp is significant, 
          ! do not use TMS for this material at all

          if(n_different == 0) then
             
             message="TMS treatment disabled for material " // trim( &
                  to_str(mat % id)) // ": nuclide temperatures " &
                  // "are already equal to tmstemp ( " // trim( &
                  to_str(mat % tmstemp / K_BOLTZMANN)) // " K )"
             call warning()

             mat % tmstemp = -1.000
             n = n - 1
          end if

       end if
    end do       
    
  end function tms_set_maximum_temperatures


!===============================================================================
! TMS_SET_THRESHOLDS finds maximum energies for which TMS will be used 
!===============================================================================

  subroutine tms_set_thresholds()
    integer :: i
    integer :: m           ! loop index for reactions

    real(8) :: low_threshold ! lower threshold energy 
    real(8) :: ures_boundary ! energy of ures boundary
    real(8) :: E_t        ! threshold energy 

    type(Reaction),    pointer :: rxn => null()
    type(Nuclide),     pointer :: nuc => null() ! nuclide pointer

    integer :: q

    do i=1, n_nuclides_total
       
       nuc => nuclides(i)
       
       if( nuc % max_kT > nuc % kT ) then
          
          ! Find lowest energy threshold of a reaction
          low_threshold = 20.00
          
          do m = 1, nuc % n_reaction
             
             rxn => nuc % reactions(m)                                      
             
             if (rxn % threshold > 1 ) then
                E_t = nuc % energy(rxn % threshold )
                
                if(E_t < low_threshold ) &
                     low_threshold = E_t  
                                
             end if
             
          end do
          
          ! Find ures boundary for nuclide (if exists)
          
          ures_boundary=20.0;
                    
          if ( nuc % urr_present) then
       
             ures_boundary=nuc % urr_data % energy(1)
             
          end if
    
          ! According to the NJOY99 manual, the BROADR module of NJOY
          ! Doppler broadens nuclides until one of the following conditions
          ! is met:
          !
          ! - E > 1.00 MeV (this is the default value)
          ! - E > lower energy boundary of the region of unresolved resonances
          ! - E > lowest threshold energy for a threshold reaction
          !
          ! For the results to be in accordance with NJOY and to 
          ! save some calculation work, the same conditions are adopted with TMS
          
          nuc % tms_emax = min(ures_boundary, min(low_threshold, 1.00))
          
       end if
    
    end do
  end subroutine tms_set_thresholds



!===============================================================================
! TMS_CALCULATE_MAJORANTS generates the microscopic majorant cross sections 
! for each nuclide. ("temperature majorant" of total cross section)
!===============================================================================

  subroutine tms_calculate_majorants()

    class(Nuclide), pointer   :: nuc ! pointer to nuclide     
    integer :: u           ! loop index for nuclides 
    integer :: i, ii       ! loop indices over energy grid 
    integer :: ilow, ihigh ! indices used for searching e
    real(8) :: dT       ! delta-T, temperature diff between xs and majorant
    real(8) :: dkT      ! delta-kT, temperature diff between xs and majorant in MeV
    real(8) :: e        ! energy 
    real(8) :: ar       ! awr / dkT 
    real(8) :: df       ! factor used in determining of DBRC-style majorant 
                        ! energy boundaries
    real(8) :: emin, emax ! energy boundaries for majorant
    real(8) :: max_xs   ! maximum value of cross section 
    real(8) :: f        ! factor for cross section interpolation 
    real(8) :: xs       ! temporary storage for cross sections
    
    ! loop over all nuclides and calculate majorant for nuclides with 
    ! max_kT >= 0 (=nuclides to be used with TMS)

    do u = 1, n_nuclides_total

       nuc => nuclides(u)
       
       if ( nuc % max_kT >= 0.0 ) then
          
          ! Calculate dT & ar
          dkT= nuc % max_kT - nuc % kT
          dT = (nuc % max_kT - nuc % kT) / K_BOLTZMANN    
          ar = nuc % awr / dkT;
          
          ! Write message 
          message = "Generating TMS majorant for nuclide " // & 
               trim(nuc % name) // ". Delta-T = " // trim( to_str(dT) ) &
               // " K."   
          call write_message()
          
          ! Allocate memory for majorant cross section 
          
          allocate(nuc % tms_majorant(nuc % n_grid))

          ! initial guesses for the lower and upper boundary

          ilow = 1
          ihigh = 1

          ! loop over energy grid

          do i=1, nuc % n_grid
            
             ! Current energy grid point

             e = nuc % energy(i);
             
             ! If we are above the upper TMS energy limit for this nuclide, 
             ! just copy cross section as-is to majorant and continue 
             ! with the next grid point
             
             if ( e > nuc % tms_emax ) then
                
                nuc % tms_majorant(i) = nuc % total(i) 
                cycle
             end if

             !===========================================!
             ! Determine energy limits corresponding to e!
             !===========================================!

             ! This is done the old "DBRC" way (Becker, Dagan, Rothenstein)

             df = 4.0/sqrt(ar*e);
             
             ! If the neutron energy is small enough compared to the energy of 
             ! the target, the relative collision energy can be arbitrarily small
             ! --> Let's approximate 0.0 eV with the lowest grid point

             if ( e < 16.0 / ar) then
                emin = nuc % energy(1)            
             else
                ! otherwise, calculate the minimum boundary energy normally
                ! min energy corresponds to a parallel collision with mu = 1

                emin = e*(1 - df)*(1 - df)
             end if
              
             ! Calculate maximum energy 
             ! max energy corresponds to a head-on collision with mu = -1
             
             emax = e*(1 + df)*(1 + df);

             ! A few sanity checks
             
             ! This should never happen
             if ( emin < nuc % energy(1) ) then
                emin = nuc % energy(1)
             end if

             ! This happens.
             if ( emax > nuc % energy(nuc % n_grid) ) then
                emax = nuc % energy(nuc  % n_grid)
             end if

             !====================================================!
             ! Find max total xs between the limits emin and emax !
             !====================================================!
             
             ! First find energy grid indices for the lowest and 
             ! highest point that are within the boundaries

             ! As the boundaries are changing monotonously, the boundaries
             ! of the previous point act as a good first guess 

             do while ( nuc % energy(ilow) < emin )
                ilow = ilow + 1                
             end do
             
             if (ihigh < ilow) then
                ihigh = ilow
             end if
             
             do while ( ihigh + 1 < nuc % n_grid ) 
                
                if( nuc % energy(ihigh + 1) >= emax ) exit
                                
                ihigh = ihigh + 1  
                
             end do

             ! look for maximum between ilow and ihigh
             
             max_xs = nuc % total(ilow);
            
             do ii = ilow + 1, ihigh 
                
                if ( nuc % total(ii) > max_xs ) then
                   max_xs = nuc % total(ii)
                end if
                
             end do

             ! interpolate extremes and check 

             ! low:
             if( ilow == 1 ) then
                xs = nuc % total(1)
             else
                f = (emin - nuc % energy(ilow-1))/( nuc % energy(ilow) - nuc % energy(ilow-1))
                xs = nuc % total(ilow-1) + f * (nuc % total(ilow) - nuc % total(ilow-1))
             end if
             
             if( xs > max_xs) max_xs = xs

             ! high:

             if( ihigh >= nuc % n_grid ) then
                xs = nuc % total(nuc % n_grid)
             else
                f = (emax - nuc % energy(ihigh))/( nuc % energy(ihigh + 1) &
                  - nuc % energy(ihigh))
                xs = nuc % total(ihigh) + f * (nuc % total(ihigh+1) - nuc % total(ihigh))
             end if

             if ( xs > max_xs) max_xs = xs
             
             ! Multiply by D-b integral for constant cross section and store 
             
             nuc % tms_majorant(i) = max_xs * cdintegral(e, dkT, nuc % awr)
                          
          end do

       end if

    end do
    
  end subroutine tms_calculate_majorants

      
!===============================================================================
! TMS_SAMPLE_NUCLIDE samples the target nuclide candidate based on majorant 
! cross sections
!===============================================================================
  
      function tms_sample_nuclide(p) result(i_nuclide)
        
        type(Particle), intent(in) :: p ! particle pointer 
        
        type(Material), pointer :: mat ! material pointer 
        integer :: i                   ! loop variable over nuclides in material
        integer :: i_nuclide           ! nuclide index
        real(8) :: xs_nuc              ! auxiliary varaible used in sampling
        real(8) :: E                   ! energy

        ! set energy and material pointer 

        E = p % E        
        mat => materials(p % material)

        ! Sample proportion of macroscopic xs

        xs_nuc = material_xs % tms_majorant * prn()                                  
        
        ! Loop over nuclides until the sampled nuclide is found

        ! NOTE: arranging the nuclide array properly (most probable nuclide
        ! first) would increase the performance when modelling burned 
        ! materials 

        do i = 1, mat % n_nuclides
           
           i_nuclide = mat % nuclide(i) 
                   
           xs_nuc = xs_nuc - mat % atom_density(i) * & 
                micro_xs(i_nuclide) % tms_majorant
           
           if(xs_nuc < 0.0) exit

        end do
        
        ! Check 
        if( xs_nuc > 0.0 ) then
           message = "TMS nuclide sampling failed" 
           call fatal_error()
        end if
        
      end function tms_sample_nuclide

!=============================================================================
! TMS_RESET_XS_SUMS and the two following subroutines are related to storing
! the sampled cross sections during TMS neutron tracking. The stored cross
! sections are used later to calculate macroscopic cross sections to be used 
! in scoring of tallies. TMS_RESET_XS_SUMS is called in the beginning of each
! track. 
!=============================================================================

subroutine tms_reset_xs_sums(mat)
  type(Material),intent(in) :: mat ! material pointer
  
  integer :: n                     ! loop index for nuclides
  integer :: i_nuclide             ! nuclide index 
  type(Nuclide), pointer:: nuc     ! nuclide pointer 

  do n=1,mat % n_nuclides 
     i_nuclide = mat % nuclide(n) 

     micro_xs(i_nuclide) % sum_total = ZERO
     micro_xs(i_nuclide) % sum_elastic = ZERO
     micro_xs(i_nuclide) % sum_absorption = ZERO
     micro_xs(i_nuclide) % sum_fission = ZERO
     micro_xs(i_nuclide) % sum_nu_fission = ZERO
     micro_xs(i_nuclide) % sum_kappa_fission = ZERO

     micro_xs(i_nuclide) % tms_n_samples = 0

  end do

end subroutine tms_reset_xs_sums


!=============================================================================
! TMS_ACCUMULATE_XS_SUMS adds previously sampled microssopic cross sections
! of i_nuclide to the global sums
!=============================================================================

subroutine tms_accumulate_xs_sums(i_nuclide)
  integer, intent(in) :: i_nuclide ! nuclide index
  
  ! add to sums

  micro_xs(i_nuclide) % sum_total = micro_xs(i_nuclide) % sum_total + &
       micro_xs(i_nuclide) % total
  micro_xs(i_nuclide) % sum_elastic = micro_xs(i_nuclide) % sum_elastic + &
       micro_xs(i_nuclide) % elastic
  micro_xs(i_nuclide) % sum_absorption = micro_xs(i_nuclide) % sum_absorption +&
       micro_xs(i_nuclide) % absorption
  micro_xs(i_nuclide) % sum_fission = micro_xs(i_nuclide) % sum_fission + &
       micro_xs(i_nuclide) % fission
  micro_xs(i_nuclide) % sum_nu_fission = micro_xs(i_nuclide) % sum_nu_fission + &
       micro_xs(i_nuclide) % nu_fission
  micro_xs(i_nuclide) % sum_kappa_fission = micro_xs(i_nuclide) % &
       sum_kappa_fission +  micro_xs(i_nuclide) % kappa_fission
  
  ! increase counter
  micro_xs(i_nuclide) % tms_n_samples = micro_xs(i_nuclide) % tms_n_samples + 1 

end subroutine tms_accumulate_xs_sums

end module tms
