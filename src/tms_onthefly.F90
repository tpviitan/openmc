module tms_onthefly 

!===============================================================================
! 
!===============================================================================

  use error,  only: fatal_error
  use global
  use string, only: to_str
  use output, only: write_message

      ! Global variables here?
#ifdef MPI
      use mpi
#endif

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
       ! if there are no TMS materials in the system, init nothing
       return
    end if

! Print message (this is unimportant and can be removed )
    
    message="Initializing TMS for all nuclides in " // trim(to_str(n)) &
         // " materials." 
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
             
             if ( nuclides(i) % max_kT < mat % tmstemp ) then
                nuclides(i) % max_kT = mat % tmstemp
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
