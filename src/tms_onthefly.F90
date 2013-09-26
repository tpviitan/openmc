module tms_onthefly 

!===============================================================================
! 
!===============================================================================

  use error,  only: fatal_error
  use global
  use string, only: to_str
  use output, only: write_message

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
             
             ! Store maximum cross section as majorant
             
             nuc % tms_majorant(i) = max_xs
 
            
             
!             if (nuc % zaid == 92238 .and. e > 6.0e-6 .and. e < 7.5e-6 ) then
!                message = " " // trim(to_str(e)) // " " // trim(to_str(nuc % tms_majorant(i))) // " " &
!                     // trim(to_str(nuc % total(i))) 
!                call write_message()
!             end if

          end do

       end if

    end do
    
    end subroutine calculate_tms_majorants

!===============================================================================
! CDINTEGRAL calculates the Doppler-broadening integral for constant cross 
! section. In the TMS articles [1-x], this is denoted with g(E,T,A)
!===============================================================================


      function cdintegral(e,dt,awr) result(cdint) 
        real(8), intent(in)      :: e    ! Neutron energy (always L-frame)
        real(8), intent(in)      :: dt   ! Delta T
        real(8), intent(in)      :: awr  ! Atomic Weight Ratio (AWR)
        
        real(8) :: a, ainv, cdint;
        
        a = sqrt(awr*e/(dt*K_BOLTZMANN));
        ainv = 1/a;

! This can be easily optimized, since the value is practically 1.0
! at high energies 

        cdint=(1.0 + 0.5*ainv*ainv)*erf(a) + exp(-a**2)*ainv/SQRTPI;
      
      end function cdintegral

      
      

end module tms_onthefly
