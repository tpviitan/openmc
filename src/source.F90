module source

  use bank_header,     only: Bank
  use constants
  use error,           only: fatal_error
  use geometry_header, only: BASE_UNIVERSE
  use global
  use math,            only: maxwell_spectrum, watt_spectrum
  use output,          only: write_message
  use particle_header, only: Particle
  use random_lcg,      only: prn, set_particle_seed
  use string,          only: to_str

#ifdef MPI
  use mpi
#endif

  implicit none

contains

!===============================================================================
! INITIALIZE_SOURCE initializes particles in the source bank
!===============================================================================

  subroutine initialize_source()

    integer(8) :: i          ! loop index over bank sites
    integer(8) :: id         ! particle id

    type(Bank), pointer :: src => null() ! source bank site

    message = "Initializing source particles..."
    call write_message(6)

    if (path_source /= '') then
      ! Read the source from a binary file instead of sampling from some
      ! assumed source distribution

      message = 'This feature is currently disabled and will be added back in.'
      call fatal_error()

    else
      ! Generation source sites from specified distribution in user input
      do i = 1, work
        ! Get pointer to source bank site
        src => source_bank(i)

        ! initialize random number seed
        id = work_index(rank) + i
        call set_particle_seed(id)

        ! sample external source distribution
        call sample_external_source(src)
      end do
    end if

  end subroutine initialize_source

!===============================================================================
! SAMPLE_EXTERNAL_SOURCE
!===============================================================================

  subroutine sample_external_source(site)

    type(Bank), pointer :: site ! source site

    integer :: i          ! dummy loop index
    real(8) :: r(3)       ! sampled coordinates
    real(8) :: phi        ! azimuthal angle
    real(8) :: mu         ! cosine of polar angle
    real(8) :: p_min(3)   ! minimum coordinates of source
    real(8) :: p_max(3)   ! maximum coordinates of source
    real(8) :: a          ! Arbitrary parameter 'a'
    real(8) :: b          ! Arbitrary parameter 'b'

    ! Set weight to one by default
    site % wgt = ONE

    ! Sample position
    select case (external_source % type_space)
    case (SRC_SPACE_BOX)
      ! Coordinates sampled uniformly over a box
      p_min = external_source % params_space(1:3)
      p_max = external_source % params_space(4:6)
      r = (/ (prn(), i = 1,3) /)
      site % xyz = p_min + r*(p_max - p_min)

    case (SRC_SPACE_POINT)
      ! Point source
      site % xyz = external_source % params_space

    end select

    ! Sample angle
    select case (external_source % type_angle)
    case (SRC_ANGLE_ISOTROPIC)
      ! Sample isotropic distribution
      phi = TWO*PI*prn()
      mu = TWO*prn() - ONE
      site % uvw(1) = mu
      site % uvw(2) = sqrt(ONE - mu*mu) * cos(phi)
      site % uvw(3) = sqrt(ONE - mu*mu) * sin(phi)

    case (SRC_ANGLE_MONO)
      ! Monodirectional source
      site % uvw = external_source % params_angle

    case default
      message = "No angle distribution specified for external source!"
      call fatal_error()
    end select

    ! Sample energy distribution
    select case (external_source % type_energy)
    case (SRC_ENERGY_MONO)
      ! Monoenergtic source
      site % E = external_source % params_energy(1)

    case (SRC_ENERGY_MAXWELL)
      a = external_source % params_energy(1)
      do
        ! Sample Maxwellian fission spectrum
        site % E = maxwell_spectrum(a)

        ! resample if energy is >= 20 MeV
        if (site % E < 20) exit
      end do

    case (SRC_ENERGY_WATT)
      a = external_source % params_energy(1)
      b = external_source % params_energy(2)
      do
        ! Sample Watt fission spectrum
        site % E = watt_spectrum(a, b)

        ! resample if energy is >= 20 MeV
        if (site % E < 20) exit
      end do

    case default
      message = "No energy distribution specified for external source!"
      call fatal_error()
    end select

  end subroutine sample_external_source

!===============================================================================
! GET_SOURCE_PARTICLE returns the next source particle 
!===============================================================================

  subroutine get_source_particle(p, index_source)

    type(Particle), intent(inout) :: p
    integer(8),     intent(in)    :: index_source

    integer(8) :: particle_seed  ! unique index for particle
    type(Bank), pointer, save :: src => null()
!$omp threadprivate(src)

    ! set defaults
    call p % initialize()

    ! Copy attributes from source to particle
    src => source_bank(index_source)
    call copy_source_attributes(p, src)

    ! set identifier for particle
    p % id = work_index(rank) + index_source

    ! set random number seed
    particle_seed = (overall_gen - 1)*n_particles + p % id
    call set_particle_seed(particle_seed)

    ! set particle trace
    trace = .false.
    if (current_batch == trace_batch .and. current_gen == trace_gen .and. &
         p % id == trace_particle) trace = .true.

  end subroutine get_source_particle

!===============================================================================
! COPY_SOURCE_ATTRIBUTES
!===============================================================================

  subroutine copy_source_attributes(p, src)

    type(Particle), intent(inout) :: p
    type(Bank),     pointer       :: src

    ! copy attributes from source bank site
    p % wgt         = src % wgt
    p % last_wgt    = src % wgt
    p % coord % xyz = src % xyz
    p % coord % uvw = src % uvw
    p % last_xyz    = src % xyz
    p % E           = src % E
    p % last_E      = src % E

  end subroutine copy_source_attributes

end module source
