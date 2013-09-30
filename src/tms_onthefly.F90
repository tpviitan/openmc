module tms_onthefly 

!===============================================================================
! 
!===============================================================================

  use error,  only: fatal_error
  use global
  use string, only: to_str
  use output, only: write_message
  use particle_header, only: Particle
  use ace_header 

  implicit none

contains 

!===============================================================================
! INIT_TMS makes all the initializations needed for transport with TMS 
! (find maximum temperature for each nuclide, generate majorants)
!===============================================================================

  subroutine init_tms()
    integer :: n

    n = set_maximum_temperatures() 
    
    if( n == 0 ) then
       ! if there are no TMS materials in the problem, init nothing
       return
    end if

! Print message (this is unimportant and can be removed )
    
    message="Initializing TMS for all nuclides in " // trim(to_str(n)) &
         // " materials..." 
    call write_message()

    call calculate_tms_majorants()

  end subroutine init_tms

!===============================================================================
! SET_MAXIMUM_TEMPERATURES finds the maximum temperature of each nuclide
! for which the TMS treatment is used. Returns the number of TMS materials. 
!===============================================================================

  function set_maximum_temperatures() result(n)
    integer :: i        ! loop index for materials
    integer :: u        ! loop index for nuclides 
    integer :: n        ! Number of TMS materials 
    integer :: nuclide_i ! index of nuclide
    type(Material),    pointer :: mat => null()
    type(Nuclide),     pointer :: nuc => null()

    n = 0;
    
    ! loop over materials
    
    do i = 1, n_materials

       mat => materials(i) 
       
       if ( mat % tmstemp >= 0.0 ) then
          n = n + 1

          ! loop over nuclides in material and set maximum 
          
          do u = 1, mat % n_nuclides
             
             nuclide_i = mat % nuclide(u)
             
             ! If material temperature is larger than max T of nuclide,
             ! update max_kT
             
             if ( nuclides(nuclide_i) % max_kT < mat % tmstemp ) then
                nuclides(nuclide_i) % max_kT = mat % tmstemp
             end if
                          
          end do                    
       end if
    end do       
    
  end function set_maximum_temperatures

!===============================================================================
! CALCULATE_TMS_MAJORANTS generates the microscopic majorant cross sections 
! for each nuclide. ("temperature majorant of total cross section")
!===============================================================================

  subroutine calculate_tms_majorants()

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

    
! loop over all nuclides and calculate majorant for nuclides with 
! max_kT >= 0 (nuclides to be used with TMS)

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

             !===========================================!
             ! Determine energy limits corresponding to e!
             !===========================================!

             ! This is done the old "DBRC" way (Becker, Dagan, Rothenstein)
             ! More efficient ways exist but remain to be published

             df = 4.0/sqrt(ar*e);
             
             ! If the neutron energy is small enough compared to the energy of 
             ! the target, the relative collision energ can be arbitrarily small
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
             
             do while ( nuc % energy(ihigh + 1) < emax )
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

             f = (emin - nuc % energy(ilow-1))/( nuc % energy(ilow) - nuc % energy(ilow-1))
             if ( nuc % total(ilow-1) + f * (nuc % total(ilow) - nuc % total(ilow-1)) &
                  > max_xs ) then

                max_xs = nuc % total(ilow-1) + f * (nuc % total(ilow) - nuc % total(ilow-1))

             end if
             
             ! high:
             f = (emax - nuc % energy(ihigh))/( nuc % energy(ihigh + 1) &
                  - nuc % energy(ihigh))
             if ( nuc % total(ihigh) + f * (nuc % total(ihigh+1) - nuc % total(ihigh)) &
                  > max_xs ) then

                max_xs = nuc % total(ihigh) + f * (nuc % total(ihigh+1) - &
                     nuc % total(ihigh))

             end if
             
             ! Multiply by D-b integral for constant cross section and store 
             
             nuc % tms_majorant(i) = max_xs * cdintegral(e, dkT, nuc % awr)
                          
          end do

       end if

    end do
    
    end subroutine calculate_tms_majorants

!===============================================================================
! CDINTEGRAL calculates the Doppler-broadening integral for constant cross 
! section. In the TMS articles [1-x], this is denoted with g(E,T,A)
!===============================================================================


      function cdintegral(E,dkT,awr) result(cdint) 
        real(8), intent(in)      :: E    ! Neutron energy (always L-frame)
        real(8), intent(in)      :: dkT  ! Delta kT (MeV)
        real(8), intent(in)      :: awr  ! Atomic Weight Ratio (AWR)
        
        real(8) :: a, ainv, cdint;
        
        a = sqrt(awr*e/(dkT));
        ainv = 1/a;

! This can be easily optimized, since the value is practically 1.0
! at high energies 

        cdint=(1.0 + 0.5*ainv*ainv)*erf(a) + exp(-a**2)*ainv/SQRTPI;
      
      end function cdintegral
          
!===============================================================================
! TMS_SAMPLE_NUCLIDE samples the target nuclide candidate and stores in nuc 
!===============================================================================
  
      function tms_sample_nuclide(p) result(i_nuclide)

        type(Particle), intent(in) :: p ! particle pointer 

        type(Nuclide), pointer :: nuc ! nuclide pointer 
        type(Material), pointer :: mat ! material pointer 
        integer :: i       ! loop variable over nuclides in material
        integer :: i_grid  ! energy grid index
        integer :: i_nuclide  ! nuclide index
        real(8) :: xs_nuc  ! aux varaible used in TMS nuclide sampling
        real(8) :: E       ! energy
        real(8) :: xs_maj_nuc ! nuclide-wise microscopic majorant

        ! set energy and material pointer 
        E = p % E        
        mat => materials(p % material)

        ! Sample proportion of macroscopic xs

        xs_nuc = material_xs % tms_majorant * prn()                                  
        
        ! init index
        i = 1

        ! arranging the nuclide array properly (most probable nuclide
        ! first) would increase significantly the performance in 
        ! case of burned materials 

        do i = 1, mat % n_nuclides
           
           nuc => nuclides(i_nuclide)
        
           xs_nuc = xs_nuc - mat % atom_density(i) * &
                max( nuc % tms_majorant(i_grid), nuc % tms_majorant(i_grid+1) )
                           
           if(xs_nuc < 0.0) exit

        end do

        if( xs_nuc > 0.0 ) then
           message = "TMS nuclide sampling failed" 
           call fatal_error()
        end if
                
        i_nuclide = mat % nuclide(i) 

      end function tms_sample_nuclide

    
   !======================================================================

   !======================================================================

   subroutine tms_update_majorants(p)   
     
     type(Particle), intent(inout) :: p ! Particle pointer 
     real(8) :: atom_density 
     integer :: i_nuclide
     integer :: i_grid
     integer :: i ! loop variable over nuclides 

     type(Material), pointer :: mat 
     type(Nuclide), pointer :: nuc 
     
     material_xs % tms_majorant = ZERO
     
     mat => materials(p % material)
     
     ! Find energy index on unionized grid          
     if (grid_method == GRID_UNION) call find_energy_index(p % E)
     
     do i = 1, mat % n_nuclides
        
        atom_density = mat % atom_density(i)
        i_nuclide = mat % nuclide(i) 
        nuc => nuclides(i_nuclide)
        
        ! Get i_grid
        
        select case(grid_method)
        case(GRID_UNION)
           i_grid = nuc % grid_index(union_grid_index)          
           
        case(GRID_NUCLIDE)
           
           ! If we're not using the unionized grid, we have to do a binary search on
           ! the nuclide energy grid in order to determine which points to
           ! interpolate between
           
           if (p % E < nuc % energy(1)) then
              i_grid = 1
           elseif (p % E > nuc % energy(nuc % n_grid)) then
              i_grid = nuc % n_grid - 1
           else
              i_grid = binary_search(nuc % energy, nuc % n_grid, p % E)
           end if
           
        end select
        
        micro_xs(i_nuclide) % tms_majorant = &
             max(nuc % tms_majorant(i_grid), nuc % tms_majorant(i_grid + 1))
        
        material_xs % tms_majorant = material_xs % tms_majorant + &
             atom_density * micro_xs(i_nuclide) % tms_majorant
        
     end do
   
   end subroutine tms_update_majorants


!===============================================================================
! TMS_CALCULATE_NUCLIDE_XS 
!===============================================================================

      ! S(a,b) -materiaalitarkistus pitaa saada jonnekin!

      subroutine tms_calculate_nuclide_xs(i_nuclide, kT_mat, E) 
        integer, intent(in) :: i_nuclide ! index of nuclide 
        real(8), intent(in) :: kT_mat  ! material temperature in MeV
        real(8), intent(in) :: E       ! L-frame energy of neutron         
        
        type(Nuclide), pointer :: nuc ! nuclide pointer
        integer :: i_grid ! energy grid index 
        real(8) :: beta_vn ! beta * speed of neutron 
        real(8) :: beta_vt ! beta * speed of target
        real(8) :: beta_vr_sq ! square of beta * relative speed 
        real(8) :: alpha   ! 
        real(8) :: c     !
        real(8) :: Er
        real(8) :: accept_prob ! 
        real(8) :: beta_vt_sq !
        real(8) :: mu 
        real(8) :: r1, r2 
        real(8) :: f ! xs interpolation factor        
        real(8) :: kT      ! temperature difference in MeV
        real(8) :: cdint ! Doppler integral for constant xs
     
        ! set nuclide pointer 

        nuc => nuclides(i_nuclide) 
        
        ! amount of thermal motion to be "added" is the differenct between
        ! the material temperature and the nuclide (xs) temperature 
        
        kT = kT_mat - nuc % kT 
     
        !========================
        ! First some checks 
        !========================
        
        ! This should be checked already in the initialization phase
        ! (is it ?) 
        if ( kT < -0.001*K_BOLTZMANN ) then
           message = "Negative temperature in TMS sampling: T = " &
                // trim(to_str(kT/K_BOLTZMANN))
           call fatal_error()
        else if ( kT < 0.001*K_BOLTZMANN ) then
           
           ! if the temperature difference is neglible skip sampling 
           ! and use Er = E 
           
           Er = E 
        else if ( nuc % urr_present .and. E > nuc % urr_data % energy(1)) then
           ! if the neutron is above lower boundary of ures region,
           ! skip sampling and use Er = E (no Doppler-broadening!)           

           Er = E
        else
           
           !========================
           ! Sample target velocity
           !========================
        
           ! The distribution from which the target velocity is sampled is the same
           ! as in the free gas treatment -> also the sampling procedure is the same 
           ! (decided to copy it here from sample_target_velocity in case DBRC will 
           ! be implemented in near future and to get rid of energy thresholds ) 
  
           ! calculate beta
           beta_vn = sqrt(nuc % awr * E / kT)
           
           alpha = ONE/(ONE + sqrt(pi)*beta_vn/TWO)
        
           do
              ! Sample two random numbers
              r1 = prn()
              r2 = prn()
              
              if (prn() < alpha) then
                 ! With probability alpha, we sample the distribution p(y) =
                 ! y*e^(-y). This can be done with sampling scheme C45 frmo the Monte
                 ! Carlo sampler
                 
                 beta_vt_sq = -log(r1*r2)
                 
              else
                 ! With probability 1-alpha, we sample the distribution p(y) = y^2 *
                 ! e^(-y^2). This can be done with sampling scheme C61 from the Monte
                 ! Carlo sampler
                 
                 c = cos(PI/TWO * prn())
                 beta_vt_sq = -log(r1) - log(r2)*c*c
              end if
              
              ! Determine beta * vt
              beta_vt = sqrt(beta_vt_sq)
              
              ! Sample cosine of angle between neutron and target velocity
              mu = TWO*prn() - ONE
              
              
              beta_vr_sq = beta_vn*beta_vn + beta_vt_sq - 2*beta_vn*beta_vt*mu
              ! Determine rejection probability
              accept_prob = sqrt(beta_vr_sq) /(beta_vn + beta_vt)
              
              ! Perform rejection sampling on vt and mu
              if (prn() < accept_prob) exit
           end do
           
           ! Store "relative energy" = energy corresponding to rel. velocity
           Er = beta_vr_sq * kT / nuc % awr 
           
        end if
     
        !======================================================
        ! Find cross sections corresponding to this energy
        !======================================================
        
        message="E " // to_str(E) // " Er: " // to_str(Er) // " nuc " // nuc % name
        call write_message()
        
        ! Find energy index on unionized grid          
        if (grid_method == GRID_UNION) call find_energy_index(Er)
        
        ! Get i_grid
        
        select case(grid_method)
        case(GRID_UNION)
           i_grid = nuc % grid_index(union_grid_index)          
        
        case(GRID_NUCLIDE)
           
           ! If we're not using the unionized grid, we have to do a binary search on
           ! the nuclide energy grid in order to determine which points to
           ! interpolate between
           
           if (E < nuc % energy(1)) then
              i_grid = 1
           elseif (E > nuc % energy(nuc % n_grid)) then
              i_grid = nuc % n_grid - 1
           else
              i_grid = binary_search(nuc % energy, nuc % n_grid, E)
           end if
           
        end select
        
        ! check for rare case where two energy points are the same
        if (nuc % energy(i_grid) == nuc % energy(i_grid+1)) i_grid = i_grid + 1
        
        ! interpolation factor
        
        f = ( Er - nuc % energy(i_grid)) / & 
             ( nuc % energy( i_grid + 1 ) - nuc % energy(i_grid))
        
        ! get Doppler-integral for constant cross sections only once
        ! since it is quite costly to calculate
        
        cdint=cdintegral(E, kT, nuc % awr )
                        
        !==========================================================================
        ! do the same things as in calculate_nuclide_xs, but multiply cross
        ! sections with cdint.
        !==========================================================================
        
        micro_xs(i_nuclide) % index_grid    = i_grid
        micro_xs(i_nuclide) % interp_factor = f
        
        ! Initialize sab treatment to false
        micro_xs(i_nuclide) % index_sab   = NONE
        micro_xs(i_nuclide) % elastic_sab = ZERO
        micro_xs(i_nuclide) % use_ptable  = .false.
        
        ! Initialize nuclide cross-sections to zero
        micro_xs(i_nuclide) % fission    = ZERO
        micro_xs(i_nuclide) % nu_fission = ZERO
        micro_xs(i_nuclide) % kappa_fission  = ZERO
        
        ! Calculate microscopic nuclide total cross section
        micro_xs(i_nuclide) % total = ( (ONE - f) * nuc % total(i_grid) &
             + f * nuc % total(i_grid+1) ) * cdint 
        
        ! Calculate microscopic nuclide total cross section
        micro_xs(i_nuclide) % elastic = ( (ONE - f) * nuc % elastic(i_grid) &
             + f * nuc % elastic(i_grid+1) ) * cdint
        
        ! Calculate microscopic nuclide absorption cross section
        micro_xs(i_nuclide) % absorption = ( (ONE - f) * nuc % absorption( &
             i_grid) + f * nuc % absorption(i_grid+1) ) * cdint 
        
        if (nuc % fissionable) then
           ! Calculate microscopic nuclide total cross section
           micro_xs(i_nuclide) % fission = ( (ONE - f) * nuc % fission(i_grid) &
                + f * nuc % fission(i_grid+1) ) * cdint 
           
           ! Calculate microscopic nuclide nu-fission cross section
           micro_xs(i_nuclide) % nu_fission = ( (ONE - f) * nuc % nu_fission( &
                i_grid) + f * nuc % nu_fission(i_grid+1) ) * cdint 
           
           ! Calculate microscopic nuclide kappa-fission cross section
           ! The ENDF standard (ENDF-102) states that MT 18 stores
           ! the fission energy as the Q_value (fission(1))
           micro_xs(i_nuclide) % kappa_fission = &
                nuc % reactions(nuc % index_fission(1)) % Q_value * &
                micro_xs(i_nuclide) % fission * cdint 
        end if
        
        ! Do not store energy, because we want to recalculate these 
        ! for each collision 

        micro_xs(i_nuclide) % last_E = ZERO
        
   end subroutine tms_calculate_nuclide_xs
         
 end module tms_onthefly
 
