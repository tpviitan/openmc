module cross_section

  use ace_header,      only: Nuclide, SAlphaBeta, Reaction, UrrData
  use constants
  use error,           only: fatal_error, warning
  use fission,         only: nu_total
  use global
  use material_header, only: Material
  use particle_header, only: Particle
  use random_lcg,      only: prn
  use search,          only: binary_search
  use string,          only: to_str
  use output,          only: write_message

  implicit none
  save

  integer :: union_grid_index
!$omp threadprivate(union_grid_index)

contains

!===============================================================================
! CALCULATE_XS determines the macroscopic cross sections for the material the
! particle is currently traveling through.
!===============================================================================

  subroutine calculate_xs(p)

    type(Particle), intent(inout) :: p

    integer :: i             ! loop index over nuclides
    integer :: i_nuclide     ! index into nuclides array
    integer :: i_sab         ! index into sab_tables array
    integer :: j             ! index in mat % i_sab_nuclides
    real(8) :: atom_density  ! atom density of a nuclide
    logical :: check_sab     ! should we check for S(a,b) table?
    type(Material), pointer, save :: mat => null() ! current material
!$omp threadprivate(mat)

    ! Set all material macroscopic cross sections to zero
    material_xs % total      = ZERO
    material_xs % elastic    = ZERO
    material_xs % absorption = ZERO
    material_xs % fission    = ZERO
    material_xs % nu_fission = ZERO
    material_xs % kappa_fission  = ZERO

    ! Exit subroutine if material is void
    if (p % material == MATERIAL_VOID) return

    mat => materials(p % material)

    ! Find energy index on unionized grid
    if (grid_method == GRID_UNION) call find_energy_index(p % E)

    ! Determine if this material has S(a,b) tables
    check_sab = (mat % n_sab > 0)

    ! Initialize position in i_sab_nuclides
    j = 1

    ! Add contribution from each nuclide in material
    do i = 1, mat % n_nuclides
      ! ========================================================================
      ! CHECK FOR S(A,B) TABLE

      i_sab = 0

      ! Check if this nuclide matches one of the S(a,b) tables specified -- this
      ! relies on i_sab_nuclides being in sorted order
      if (check_sab) then
        if (i == mat % i_sab_nuclides(j)) then
          ! Get index in sab_tables
          i_sab = mat % i_sab_tables(j)

          ! If particle energy is greater than the highest energy for the S(a,b)
          ! table, don't use the S(a,b) table
          if (p % E > sab_tables(i_sab) % threshold_inelastic) i_sab = 0

          ! Increment position in i_sab_nuclides
          j = j + 1

          ! Don't check for S(a,b) tables if there are no more left
          if (j > mat % n_sab) check_sab = .false.
        end if
      end if

      ! ========================================================================
      ! CALCULATE MICROSCOPIC CROSS SECTION

      ! Determine microscopic cross sections for this nuclide
      i_nuclide = mat % nuclide(i)

      ! Copy atom density of nuclide in material
      atom_density = mat % atom_density(i)

      ! Calculate microscopic cross section for this nuclide
      if (p % E /= micro_xs(i_nuclide) % last_E .and. &
           mat % tmstemp < 0.0 ) then
         call calculate_nuclide_xs(i_nuclide, i_sab, p % E)
      else if( mat % tmstemp >= 0.0 ) then

         ! This samples target-at-rest cross sections for i_nuclide
         call tms_sample_nuclide_xs(i_nuclide, mat % tmstemp, p % E)
         
         ! The material macroscopic cross sections, used in scoring
         ! of keff and reaction rates, need to be multiplied by cdint
         ! in TMS mode 
         atom_density = atom_density * micro_xs(i_nuclide) % cdint         

      end if

      ! ========================================================================
      ! ADD TO MACROSCOPIC CROSS SECTION

      ! Add contributions to material macroscopic total cross section
      material_xs % total = material_xs % total + &
           atom_density * micro_xs(i_nuclide) % total

      ! Add contributions to material macroscopic scattering cross section
      material_xs % elastic = material_xs % elastic + &
           atom_density * micro_xs(i_nuclide) % elastic

      ! Add contributions to material macroscopic absorption cross section
      material_xs % absorption = material_xs % absorption + & 
           atom_density * micro_xs(i_nuclide) % absorption

      ! Add contributions to material macroscopic fission cross section
      material_xs % fission = material_xs % fission + &
           atom_density * micro_xs(i_nuclide) % fission

      ! Add contributions to material macroscopic nu-fission cross section
      material_xs % nu_fission = material_xs % nu_fission + &
           atom_density * micro_xs(i_nuclide) % nu_fission 
           
      ! Add contributions to material macroscopic energy release from fission
      material_xs % kappa_fission = material_xs % kappa_fission + &
           atom_density * micro_xs(i_nuclide) % kappa_fission
            
    end do
       
  end subroutine calculate_xs

!===============================================================================
! CALCULATE_NUCLIDE_XS determines microscopic cross sections for a nuclide of a
! given index in the nuclides array at the energy of the given particle
!===============================================================================

  subroutine calculate_nuclide_xs(i_nuclide, i_sab, E)

    integer, intent(in) :: i_nuclide ! index into nuclides array
    integer, intent(in) :: i_sab     ! index into sab_tables array
    real(8), intent(in) :: E         ! energy

    integer :: i_grid ! index on nuclide energy grid
    real(8) :: f      ! interp factor on nuclide energy grid
    type(Nuclide), pointer, save :: nuc => null()
    !$omp threadprivate(nuc)

    ! Set pointer to nuclide
    nuc => nuclides(i_nuclide)

    ! Determine index on nuclide energy grid
    select case (grid_method)
    case (GRID_UNION)
      ! If we're using the unionized grid with pointers, finding the index on
      ! the nuclide energy grid is as simple as looking up the pointer

      i_grid = nuc % grid_index(union_grid_index)

    case (GRID_NUCLIDE)
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

    ! calculate interpolation factor
    f = (E - nuc%energy(i_grid))/(nuc%energy(i_grid+1) - nuc%energy(i_grid))

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
    micro_xs(i_nuclide) % total = (ONE - f) * nuc % total(i_grid) &
         + f * nuc % total(i_grid+1)

    ! Calculate microscopic nuclide total cross section
    micro_xs(i_nuclide) % elastic = (ONE - f) * nuc % elastic(i_grid) &
         + f * nuc % elastic(i_grid+1)

    ! Calculate microscopic nuclide absorption cross section
    micro_xs(i_nuclide) % absorption = (ONE - f) * nuc % absorption( &
         i_grid) + f * nuc % absorption(i_grid+1)

    if (nuc % fissionable) then
      ! Calculate microscopic nuclide total cross section
      micro_xs(i_nuclide) % fission = (ONE - f) * nuc % fission(i_grid) &
           + f * nuc % fission(i_grid+1)

      ! Calculate microscopic nuclide nu-fission cross section
      micro_xs(i_nuclide) % nu_fission = (ONE - f) * nuc % nu_fission( &
           i_grid) + f * nuc % nu_fission(i_grid+1)
           
      ! Calculate microscopic nuclide kappa-fission cross section
      ! The ENDF standard (ENDF-102) states that MT 18 stores
      ! the fission energy as the Q_value (fission(1))
      micro_xs(i_nuclide) % kappa_fission = &
           nuc % reactions(nuc % index_fission(1)) % Q_value * &
           micro_xs(i_nuclide) % fission
    end if

    ! If there is S(a,b) data for this nuclide, we need to do a few
    ! things. Since the total cross section was based on non-S(a,b) data, we
    ! need to correct it by subtracting the non-S(a,b) elastic cross section and
    ! then add back in the calculated S(a,b) elastic+inelastic cross section.

    if (i_sab > 0) call calculate_sab_xs(i_nuclide, i_sab, E)

    ! if the particle is in the unresolved resonance range and there are
    ! probability tables, we need to determine cross sections from the table

    if (urr_ptables_on .and. nuc % urr_present) then
      if (E > nuc % urr_data % energy(1) .and. &
           E < nuc % urr_data % energy(nuc % urr_data % n_energy)) then
        call calculate_urr_xs(i_nuclide, E)
      end if
    end if

    ! Set last evaluated energy -- if we're in S(a,b) region, force
    ! re-calculation of cross-section
    if (i_sab == 0) then
      micro_xs(i_nuclide) % last_E = E
    else
      micro_xs(i_nuclide) % last_E = ZERO
    end if

  end subroutine calculate_nuclide_xs

!===============================================================================
! CALCULATE_SAB_XS determines the elastic and inelastic scattering
! cross-sections in the thermal energy range. These cross sections replace
! whatever data were taken from the normal Nuclide table.
!===============================================================================

  subroutine calculate_sab_xs(i_nuclide, i_sab, E)

    integer, intent(in) :: i_nuclide ! index into nuclides array
    integer, intent(in) :: i_sab     ! index into sab_tables array
    real(8), intent(in) :: E         ! energy

    integer :: i_grid    ! index on S(a,b) energy grid
    real(8) :: f         ! interp factor on S(a,b) energy grid
    real(8) :: inelastic ! S(a,b) inelastic cross section
    real(8) :: elastic   ! S(a,b) elastic cross section
    type(SAlphaBeta), pointer, save :: sab => null()
!$omp threadprivate(sab)

    ! Set flag that S(a,b) treatment should be used for scattering
    micro_xs(i_nuclide) % index_sab = i_sab

    ! Get pointer to S(a,b) table
    sab => sab_tables(i_sab)

    ! Get index and interpolation factor for inelastic grid
    if (E < sab % inelastic_e_in(1)) then
      i_grid = 1
      f = ZERO
    else
      i_grid = binary_search(sab % inelastic_e_in, sab % n_inelastic_e_in, E)
      f = (E - sab%inelastic_e_in(i_grid)) / & 
           (sab%inelastic_e_in(i_grid+1) - sab%inelastic_e_in(i_grid))
    end if

    ! Calculate S(a,b) inelastic scattering cross section
    inelastic = (ONE - f) * sab % inelastic_sigma(i_grid) + &
         f * sab % inelastic_sigma(i_grid + 1)

    ! Check for elastic data
    if (E < sab % threshold_elastic) then
      ! Determine whether elastic scattering is given in the coherent or
      ! incoherent approximation. For coherent, the cross section is
      ! represented as P/E whereas for incoherent, it is simply P

      if (sab % elastic_mode == SAB_ELASTIC_EXACT) then
        if (E < sab % elastic_e_in(1)) then
          ! If energy is below that of the lowest Bragg peak, the elastic
          ! cross section will be zero
          elastic = ZERO
        else
          i_grid = binary_search(sab % elastic_e_in, &
               sab % n_elastic_e_in, E)
          elastic = sab % elastic_P(i_grid) / E
        end if
      else
        ! Determine index on elastic energy grid
        if (E < sab % elastic_e_in(1)) then
          i_grid = 1
        else
          i_grid = binary_search(sab % elastic_e_in, &
               sab % n_elastic_e_in, E)
        end if

        ! Get interpolation factor for elastic grid
        f = (E - sab%elastic_e_in(i_grid))/(sab%elastic_e_in(i_grid+1) - &
             sab%elastic_e_in(i_grid))

        ! Calculate S(a,b) elastic scattering cross section
        elastic = (ONE - f) * sab % elastic_P(i_grid) + &
             f * sab % elastic_P(i_grid + 1)
      end if
    else
      ! No elastic data
      elastic = ZERO
    end if

    ! Correct total and elastic cross sections
    micro_xs(i_nuclide) % total = micro_xs(i_nuclide) % total - &
         micro_xs(i_nuclide) % elastic + inelastic + elastic
    micro_xs(i_nuclide) % elastic = inelastic + elastic

    ! Store S(a,b) elastic cross section for sampling later
    micro_xs(i_nuclide) % elastic_sab = elastic

  end subroutine calculate_sab_xs

!===============================================================================
! CALCULATE_URR_XS determines cross sections in the unresolved resonance range
! from probability tables
!===============================================================================

  subroutine calculate_urr_xs(i_nuclide, E)

    integer, intent(in) :: i_nuclide ! index into nuclides array
    real(8), intent(in) :: E         ! energy

    integer :: i_energy   ! index for energy
    integer :: i_table    ! index for table
    real(8) :: f          ! interpolation factor
    real(8) :: r          ! pseudo-random number
    real(8) :: elastic    ! elastic cross section
    real(8) :: capture    ! (n,gamma) cross section
    real(8) :: fission    ! fission cross section
    real(8) :: inelastic  ! inelastic cross section
    type(UrrData),  pointer, save :: urr => null()
    type(Nuclide),  pointer, save :: nuc => null()
    type(Reaction), pointer, save :: rxn => null()
!$omp threadprivate(urr, nuc, rxn)

    micro_xs(i_nuclide) % use_ptable = .true.

    ! get pointer to probability table
    nuc => nuclides(i_nuclide)
    urr => nuc % urr_data

    ! determine energy table
    i_energy = 1
    do
      if (E < urr % energy(i_energy + 1)) exit
      i_energy = i_energy + 1
    end do

    ! determine interpolation factor on table
    f = (E - urr % energy(i_energy)) / &
         (urr % energy(i_energy + 1) - urr % energy(i_energy))

    ! sample probability table using the cumulative distribution
    r = prn()
    i_table = 1
    do
      if (urr % prob(i_energy, URR_CUM_PROB, i_table) > r) exit
      i_table = i_table + 1
    end do

    ! determine elastic, fission, and capture cross sections from probability
    ! table
    if (urr % interp == LINEAR_LINEAR) then
      elastic = (ONE - f) * urr % prob(i_energy, URR_ELASTIC, i_table) + &
           f * urr % prob(i_energy + 1, URR_ELASTIC, i_table)
      fission = (ONE - f) * urr % prob(i_energy, URR_FISSION, i_table) + &
           f * urr % prob(i_energy + 1, URR_FISSION, i_table)
      capture = (ONE - f) * urr % prob(i_energy, URR_N_GAMMA, i_table) + &
           f * urr % prob(i_energy + 1, URR_N_GAMMA, i_table)
    elseif (urr % interp == LOG_LOG) then
      ! Get logarithmic interpolation factor
      f = log(E / urr % energy(i_energy)) / &
           log(urr % energy(i_energy + 1) / urr % energy(i_energy))

      ! Calculate elastic cross section/factor
      elastic = ZERO
      if (urr % prob(i_energy, URR_ELASTIC, i_table) > ZERO) then
        elastic = exp((ONE - f) * log(urr % prob(i_energy, URR_ELASTIC, &
             i_table)) + f * log(urr % prob(i_energy + 1, URR_ELASTIC, &
             i_table)))
      end if

      ! Calculate fission cross section/factor
      fission = ZERO
      if (urr % prob(i_energy, URR_FISSION, i_table) > ZERO) then
        fission = exp((ONE - f) * log(urr % prob(i_energy, URR_FISSION, &
             i_table)) + f * log(urr % prob(i_energy + 1, URR_FISSION, &
             i_table)))
      end if

      ! Calculate capture cross section/factor
      capture = ZERO
      if (urr % prob(i_energy, URR_N_GAMMA, i_table) > ZERO) then
        capture = exp((ONE - f) * log(urr % prob(i_energy, URR_N_GAMMA, &
             i_table)) + f * log(urr % prob(i_energy + 1, URR_N_GAMMA, &
             i_table)))
      end if
    end if

    ! Determine treatment of inelastic scattering
    inelastic = ZERO
    if (urr % inelastic_flag > 0) then
      ! Get pointer to inelastic scattering reaction
      rxn => nuc % reactions(nuc % urr_inelastic)

      ! Get index on energy grid and interpolation factor
      i_energy = micro_xs(i_nuclide) % index_grid
      f = micro_xs(i_nuclide) % interp_factor

      ! Determine inelastic scattering cross section
      if (i_energy >= rxn % threshold) then
        inelastic = (ONE - f) * rxn % sigma(i_energy - rxn%threshold + 1) + &
             f * rxn % sigma(i_energy - rxn%threshold + 2)
      end if
    end if

    ! Multiply by smooth cross-section if needed
    if (urr % multiply_smooth) then
      elastic = elastic * micro_xs(i_nuclide) % elastic
      capture = capture * (micro_xs(i_nuclide) % absorption - &
           micro_xs(i_nuclide) % fission)
      fission = fission * micro_xs(i_nuclide) % fission
    end if

    ! Set elastic, absorption, fission, and total cross sections. Note that the
    ! total cross section is calculated as sum of partials rather than using the
    ! table-provided value
    micro_xs(i_nuclide) % elastic = elastic
    micro_xs(i_nuclide) % absorption = capture + fission
    micro_xs(i_nuclide) % fission = fission
    micro_xs(i_nuclide) % total = elastic + inelastic + capture + fission

    ! Determine nu-fission cross section
    if (nuc % fissionable) then
      micro_xs(i_nuclide) % nu_fission = nu_total(nuc, E) * &
           micro_xs(i_nuclide) % fission
    end if

  end subroutine calculate_urr_xs

!===============================================================================
! FIND_ENERGY_INDEX determines the index on the union energy grid at a certain
! energy
!===============================================================================

  subroutine find_energy_index(E)

    real(8), intent(in) :: E ! energy of particle

    ! if particle's energy is outside of energy grid range, set to first or last
    ! index. Otherwise, do a binary search through the union energy grid.
    if (E < e_grid(1)) then
      union_grid_index = 1
    elseif (E > e_grid(n_grid)) then
      union_grid_index = n_grid - 1
    else
      union_grid_index = binary_search(e_grid, n_grid, E)
    end if

  end subroutine find_energy_index

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

    if ( urr_ptables_on ) then
       message="TMS cannot be used with URR ptable treatment. Please " //& 
            "add <ptables>false</ptables> in settings.xml"
       call fatal_error()
    end if

    call tms_set_thresholds()
    
! Print message (this is unimportant and can be removed )
    
    message="Initializing TMS for all nuclides in " // trim(to_str(n)) &
         // " materials..." 
    call write_message()

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

          ! Check at this point that no S(a,b) tables are associated with
          ! the material. TMS is currently (and will be in the near future)
          ! not compatible with S(a,b)
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
          ! save some calculation work, the same conditions are adopted in TMS
          
          nuc % tms_emax = min(ures_boundary, min(low_threshold, 1.00))
          
          message="nuc: " // to_str(nuc % zaid) // &
               " emax " // to_str(nuc % tms_emax )
          call write_message()
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
             
             ! If we are above the upper limit for TMS, 
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
! CDINTEGRAL calculates the Doppler-broadening integral for constant cross 
! section. In the TMS articles [1-x], this is denoted with g(E,T,A), also 
! sometimes referred to as the "normalization factor"
!===============================================================================


      function cdintegral(E,dkT,awr) result(cdint) 
        real(8), intent(in)      :: E    ! Neutron energy (always L-frame)
        real(8), intent(in)      :: dkT  ! Delta kT (MeV)
        real(8), intent(in)      :: awr  ! Atomic Weight Ratio (AWR)
        
        real(8) :: a, ainv, cdint;
        
! In case of insignificant temperature change, just return 1 
! (this avoids division-by-zero problems). This should happen only
! if a material with TMS has nuclides with cross sections both at 
! tmstemp and at a temperature T < tmstemp

        if( dkT < 0.001*K_BOLTZMANN) then
           cdint = ONE 
           return
        end if
         
        a = sqrt(awr*E/(dkT));

        ! With higher values of a the routine can be somewhat optimized
        if(a > 250) then
           cdint=ONE
        else if ( a > 2.568 ) then
           ainv = 1/a
           cdint=ONE + 0.5*ainv*ainv
        else
           ainv = 1/a
           cdint=(ONE + 0.5*ainv*ainv)*erf(a) + exp(-a*a)*ainv/SQRTPI
        end if

      end function cdintegral
          
!===============================================================================
! TMS_SAMPLE_NUCLIDE samples the target nuclide candidate based on majorant 
! cross sections
!===============================================================================
  
      function tms_sample_nuclide(p) result(i_nuclide)

        type(Particle), intent(in) :: p ! particle pointer 

        type(Material), pointer :: mat ! material pointer 
        integer :: i       ! loop variable over nuclides in material
        integer :: i_nuclide  ! nuclide index
        real(8) :: xs_nuc  ! aux varaible used in TMS nuclide sampling
        real(8) :: E       ! energy

        ! set energy and material pointer 
        E = p % E        
        mat => materials(p % material)

        ! Sample proportion of macroscopic xs

        xs_nuc = material_xs % tms_majorant * prn()                                  
        
        ! arranging the nuclide array properly (most probable nuclide
        ! first) would increase significantly the performance in 
        ! case of burned materials 

        do i = 1, mat % n_nuclides
           
           i_nuclide = mat % nuclide(i) 
                   
           xs_nuc = xs_nuc - mat % atom_density(i) * & 
                micro_xs(i_nuclide) % tms_majorant
           
           if(xs_nuc < 0.0) exit

        end do

        if( xs_nuc > 0.0 ) then
           message = "TMS nuclide sampling failed" 
           call fatal_error()
        end if
                
      end function tms_sample_nuclide
    
   !======================================================================

   !======================================================================

   subroutine tms_update_majorants(p)   
     
     type(Particle), intent(inout) :: p ! Particle pointer 

     real(8) :: atom_density     
     integer :: i_nuclide               ! index of nuclide
     integer :: i_grid                  ! energy grid index 
     integer :: i                       ! loop variable over nuclides 
     type(Material), pointer :: mat     
     type(Nuclide), pointer  :: nuc 
     
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
! TMS_SAMPLE_NUCLIDE_XS 
!===============================================================================

      ! S(a,b) -materiaalitarkistus pitaa saada jonnekin!

      subroutine tms_sample_nuclide_xs(i_nuclide, kT_mat, E) 
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
     
        ! set nuclide pointer 

        nuc => nuclides(i_nuclide) 
        
        ! amount of thermal motion to be "added" is the differenct between
        ! the material temperature and the nuclide (xs) temperature 
        
        kT = kT_mat - nuc % kT 

        ! init value 

        Er = -1.00 

        !========================
        ! First some checks 
        !========================
        
        ! This should be checked already in the initialization phase
        ! (is it ?) 
        if ( kT < -0.001*K_BOLTZMANN ) then
           message = "Negative temperature in TMS sampling: T = " &
                // trim(to_str(kT/K_BOLTZMANN))
           call fatal_error()
        else if ( kT < 0.001*K_BOLTZMANN) then
           
           ! if the temperature difference is negligible skip sampling 
           ! and use Er = E 
           
           Er = E 
           
        else if ( E > nuc % tms_emax ) then

           ! Skip TMS if E > tms_emax. For more details, see 
           ! tms_set_thresholds()

           Er = E
           
        end if

        if(Er < ZERO ) then
           
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
           
           ! Store "relative energy" = neutron energy corresponding 
           ! to rel. velocity
           Er = beta_vr_sq * kT / nuc % awr 
           
        end if
     
        !======================================================
        ! Find cross sections corresponding to Er
        !======================================================

        ! Check if the previous values can be used and return 
        
        if( Er == micro_xs(i_nuclide) % last_E) then

           ! If the last collision was in a non-tms material
           ! cdint must be updated

           micro_xs(i_nuclide) % cdint = cdintegral(E, kT, nuc % awr )
           return 
        end if
        
        ! This search could be perhaps optimized, since the energy is usually
        ! close to E (binary search might not be the optimal way)

        ! Find energy index on unionized grid          
        if (grid_method == GRID_UNION) call find_energy_index(Er)
        
        ! Get i_grid
        
        select case(grid_method)
        case(GRID_UNION)
           i_grid = nuc % grid_index(union_grid_index)          
        
        case(GRID_NUCLIDE)
                      
           if (Er < nuc % energy(1)) then
              i_grid = 1
           elseif (Er > nuc % energy(nuc % n_grid)) then
              i_grid = nuc % n_grid - 1
           else
              i_grid = binary_search(nuc % energy, nuc % n_grid, Er)
           end if
           
        end select
        
        ! check for rare case where two energy points are the same
        if (nuc % energy(i_grid) == nuc % energy(i_grid+1)) i_grid = i_grid + 1
        
        ! interpolation factor

        f = ( Er - nuc % energy(i_grid)) / & 
             ( nuc % energy( i_grid + 1 ) - nuc % energy(i_grid))
        
        ! get Doppler-integral for constant cross sections
        
        micro_xs(i_nuclide) % cdint = cdintegral(E, kT, nuc % awr )

        !==========================================================================
        ! do the same things as in calculate_nuclide_xs minus URR and S(a,b) 
        ! related stuff (neither works with TMS)
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
        micro_xs(i_nuclide) % total = (ONE - f) * nuc % total(i_grid) &
             + f * nuc % total(i_grid+1) 
        
        ! Calculate microscopic nuclide total cross section
        micro_xs(i_nuclide) % elastic = (ONE - f) * nuc % elastic(i_grid) &
             + f * nuc % elastic(i_grid+1)
        
        ! Calculate microscopic nuclide absorption cross section
        micro_xs(i_nuclide) % absorption = (ONE - f) * nuc % absorption( &
             i_grid) + f * nuc % absorption(i_grid+1) 
        
        if (nuc % fissionable) then
           ! Calculate microscopic nuclide total cross section
           micro_xs(i_nuclide) % fission =  (ONE - f) * nuc % fission(i_grid) &
                + f * nuc % fission(i_grid+1) 
           
           ! Calculate microscopic nuclide nu-fission cross section
           micro_xs(i_nuclide) % nu_fission =  (ONE - f) * nuc % nu_fission( &
                i_grid) + f * nuc % nu_fission(i_grid+1)
           
           ! Calculate microscopic nuclide kappa-fission cross section
           ! The ENDF standard (ENDF-102) states that MT 18 stores
           ! the fission energy as the Q_value (fission(1))
           micro_xs(i_nuclide) % kappa_fission = &
                nuc % reactions(nuc % index_fission(1)) % Q_value * &
                micro_xs(i_nuclide) % fission
        end if
        
        ! This is not used in TMS transport, set to ZERO just in case

        micro_xs(i_nuclide) % last_E = Er
        
   end subroutine tms_sample_nuclide_xs

end module cross_section
