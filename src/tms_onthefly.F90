module tms_onthefly 

!===============================================================================
! 
!===============================================================================
  use global

      ! Global variables here?
#ifdef MPI
      use mpi
#endif

      implicit none

contains 

!===============================================================================
! CDINTEGRAL calculates the Doppler-broadening integral for constant cross 
! section. In the TMS articles [1-x], this is denoted with g(E,T,A)
!===============================================================================

      function cdintegral(e,dt,awr) result(cdint) 
        real(8), intent(in)      :: e    ! Neutron energy (always L-frame)
        real(8), intent(in)      :: dt   ! Delta T
        real(8), intent(in)     :: awr  ! Atomic Weight Ratio (AWR)
        
        real(8) :: a, ainv, cdint;
        
        a = sqrt(awr*e/(dt*K_BOLTZMANN));
        ainv = 1/a;

! This can be easily optimized, since the value is practically 1.0
! at high energies 

        cdint=(1.0 + 0.5*ainv*ainv)*erf(a) + exp(-a**2)*ainv/SQRTPI;
      
      end function cdintegral
      

end module tms_onthefly
