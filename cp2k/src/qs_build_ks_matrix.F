!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2002  CP2K developers group                                 !
!-----------------------------------------------------------------------------!

#include "cp_prep_globals.h"

!!****h* cp2k/qs_build_ks_matrix [1.0] *
!!
!!   NAME
!!     qs_build_ks_matrix
!!
!!   FUNCTION
!!     routines that build the Kohn-Sham matrix (i.e calculate the coulomb
!!     and xc parts
!!
!!   NOTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     05.2002 moved from qs_scf (see there the history) [fawzi]
!!     JGH [30.08.02] multi-grid arrays independent from density and potential
!!
!!   SOURCE
!****************************************************************************
MODULE qs_build_ks_matrix

  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE coefficient_types,               ONLY: coeff_allocate,&
                                             coeff_copy,&
                                             coeff_sumup,&
                                             coeff_transform_space,&
                                             coeff_type,&
                                             coeff_zero
  USE cp_error_handling,               ONLY: cp_assert,&
                                             cp_assertion_failed,&
                                             cp_debug,&
                                             cp_error_dealloc_ref,&
                                             cp_error_get_logger,&
                                             cp_error_init,&
                                             cp_error_message,&
                                             cp_error_type,&
                                             cp_internal_error,&
                                             cp_unimplemented_error
  USE cp_log_handling,                 ONLY: cp_failure_level,&
                                             cp_fatal_level,&
                                             cp_log,&
                                             cp_logger_type,&
                                             cp_note_level,&
                                             cp_to_string,&
                                             cp_warning_level
  USE cube_utils,                      ONLY: cube_info_type,&
                                             destroy_cube_info,&
                                             init_cube_info
  USE dft_types,                       ONLY: dft_control_type
  USE external_potential_types,        ONLY: all_potential_type,&
                                             get_potential,&
                                             gth_potential_type
  USE functionals
  USE gaussian_gridlevels,             ONLY: destroy_gaussian_gridlevel,&
                                             gridlevel_info_type,&
                                             init_gaussian_gridlevel
  USE global_types,                    ONLY: DEBUG,&
                                             global_environment_type
  USE greens_fn,                       ONLY: pw_green_fn_init
  USE hartree,                         ONLY: calculate_hartree
  USE kinds,                           ONLY: wp=>dp
  USE l_utils,                         ONLY: destroy_l_info,&
                                             init_l_info,&
                                             l_info_type
  USE message_passing,                 ONLY: mp_allgather,&
                                             mp_max,&
                                             mp_range,&
                                             mp_sum,&
                                             mp_sync
  USE pw_env_types,                    ONLY: pw_env_get,&
                                             pw_env_type
  USE pw_grid_types,                   ONLY: HALFSPACE,&
                                             pw_grid_type
  USE pw_grids,                        ONLY: pw_find_cutoff,&
                                             pw_grid_construct,&
                                             pw_grid_setup
  USE pw_pool_types,                   ONLY: pw_pool_p_type,&
                                             pw_pool_type
  USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                             COMPLEXDATA3D,&
                                             REALDATA3D,&
                                             REALSPACE,&
                                             RECIPROCALSPACE,&
                                             pw_type
  USE pws,                             ONLY: init_pw_poisson_solver,&
                                             pw_poisson_solver
  USE qs_collocate_density,            ONLY: calculate_rho_core,&
                                             calculate_rho_elec,&
                                             calculate_total_rho
  USE qs_core_energies,                ONLY: calculate_ecore,&
                                             calculate_ecore_overlap,&
                                             calculate_ecore_self
  USE qs_energy_types,                 ONLY: qs_energy_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type,&
                                             set_qs_env
  USE qs_integrate_potential,          ONLY: integrate_v_core_rspace,&
                                             integrate_v_rspace
  USE qs_overlap,                      ONLY: write_blacs_matrix,&
                                             write_sparse_matrix
  USE qs_rho_types,                    ONLY: qs_rho_did_change,&
                                             qs_rho_type
  USE simulation_cell,                 ONLY: cell_type,&
                                             get_cell
  USE sparse_matrix_types,             ONLY: add_matrices,&
                                             allocate_matrix,&
                                             deallocate_matrix,&
                                             first_block_node,&
                                             get_block_node,&
                                             get_matrix_info,&
                                             next_block_node,&
                                             real_block_node_type,&
                                             real_matrix_p_type,&
                                             replicate_matrix_structure,&
                                             set_matrix,&
                                             symmetrise_diagonal_blocks,&
                                             transfer_matrix
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
  USE timings,                         ONLY: timeset,&
                                             timestop
  USE xc,                              ONLY: xc_calculate

  IMPLICIT NONE

  PRIVATE

  LOGICAL, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER :: moduleN='qs_build_ks_matrix'

  PUBLIC :: qs_ks_env_type
  PUBLIC :: qs_ks_init, qs_ks_dealloc_ref, qs_ks_did_change,&
       qs_ks_update_qs_env
!***
!****************************************************************************

!!****s* qs_build_ks_matrix/qs_ks_env_type [1.0] *
!!
!!   NAME
!!     qs_ks_env_type
!!
!!   FUNCTION
!!     calculation environement to calculate the ks matrix,
!!     holds all the needed vars.
!!     assumes that the core hamiltonian and energy are up to date.
!!
!!   NOTES
!!     -
!!
!!   ATTRIBUTES
!!     -
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!   SOURCE
!***************************************************************************
  TYPE qs_ks_env_type

     ! definition of all grid types used (working grids and multigrids)
     TYPE(pw_grid_type), DIMENSION(:), POINTER :: pw_grid

     ! number of the grid used as auxilliary basis for the density
     INTEGER :: auxbas_grid

     LOGICAL :: initialized, &
                s_mstruct_changed, &
                p_changed, &
                forces_up_to_date
     INTEGER :: n_evals, &
                n_electrons
     REAL(kind=wp) :: total_rho_core_rspace, &
                total_rho_rspace, &
                total_rho_gspace

     TYPE(coeff_type), DIMENSION(:), POINTER       :: rho_gspace, &
                                                      v_gspace, &
                                                      v_rspace
     TYPE(coeff_type)                              :: v_rspace_old

     ! the multigrids
     TYPE(coeff_type), DIMENSION(:), POINTER       :: mgrid_gspace, &
                                                      mgrid_rspace

     ! hartree is supposed to contain the hartree potential (for cube ouput)
     TYPE(coeff_type)          :: v_hartree_rspace

     ! the core density
     TYPE(coeff_type)          :: rho_core

     TYPE(real_matrix_p_type)  :: v
  END TYPE qs_ks_env_type
!!***
!****************************************************************************

CONTAINS

!!****f* qs_build_ks_matrix/qs_ks_init [1.0] *
!!
!!   NAME
!!     qs_ks_init
!!
!!   FUNCTION
!!     initializes the given ks_env.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - ks_env: the ks env to be initialized
!!     - qs_env: the qs environment
!!     - global_env: the global environment
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     05.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE qs_ks_init(ks_env, qs_env, global_env, error)
    TYPE(qs_ks_env_type), INTENT(inout) :: ks_env
    TYPE(qs_environment_type), INTENT(in),TARGET :: qs_env
    TYPE(global_environment_type), INTENT(in), TARGET :: global_env
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    INTEGER, PARAMETER :: maxgridpoints = 100
    CHARACTER(len=*), PARAMETER :: routineN='qs_ks_init',&
         routineP=moduleN//':'//routineN
    LOGICAL :: failure
    REAL(wp)          :: rel_cutoff,maxradius
    INTEGER :: ngrid_level
    REAL(wp), DIMENSION(:), POINTER :: cutoff
    INTEGER  :: cmax,cmaxl,igrid_level, istat,maxl
    TYPE(dft_control_type), POINTER :: dft_control
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: h
    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set

    IF (ks_env%initialized) RETURN

    failure=.FALSE.

    NULLIFY(ks_env%rho_gspace, &
         ks_env%mgrid_gspace, ks_env%mgrid_rspace,&
         ks_env%v_gspace, ks_env%v_rspace, &
         ks_env%v%matrix)

    CALL init_functionals()
    CALL get_qs_env(qs_env=qs_env,&
         dft_control=dft_control,atomic_kind_set=atomic_kind_set)

    rel_cutoff = dft_control%qs_control%relative_cutoff
    cutoff => dft_control%qs_control%e_cutoff
    ngrid_level = SIZE(cutoff)

    CALL init_grid(ks_env,qs_env=qs_env,globenv=global_env, &
         ngrid_level=ngrid_level,error=error)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         maxl=maxl, nelectron=ks_env%n_electrons)
    ks_env%n_electrons=ks_env%n_electrons-dft_control%charge

    ks_env%total_rho_core_rspace=-1.0
    ks_env%total_rho_rspace=-1.0
    ks_env%total_rho_gspace=-1.0
    ks_env%s_mstruct_changed=.TRUE.
    ks_env%p_changed=.TRUE.
    ks_env%forces_up_to_date=.FALSE.
    ks_env%n_evals=0
    ks_env%initialized=.TRUE.
  END SUBROUTINE qs_ks_init
!***************************************************************************

!!****f* qs_build_ks_matrix/qs_ks_dealloc_ref [1.0] *
!!
!!   NAME
!!     qs_ks_dealloc_ref
!!
!!   FUNCTION
!!     releases the memory used by the ks_env (but not ks_env)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - ks_env: the ks_env to be released
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     05.2002 created [fawzi]
!!
!!*** **********************************************************************
  SUBROUTINE qs_ks_dealloc_ref(ks_env,qs_env,global_env,error)
    TYPE(qs_ks_env_type), INTENT(inout) :: ks_env
    TYPE(qs_environment_type), INTENT(in),TARGET :: qs_env
    TYPE(global_environment_type), INTENT(in), TARGET :: global_env
    TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error

    LOGICAL :: failure
    CHARACTER(len=*), PARAMETER :: routineN='qs_ks_dealloc_ref',&
         routineP=moduleN//':'//routineN
    INTEGER :: igrid_level,istat
    failure=.FALSE.

    CPPrecondition(ks_env%initialized,cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL cp_unimplemented_error(routineP, "grid deallocation not done in "//&
            CPSourceFileRef,&
            error=error, error_level=cp_warning_level)

       ks_env%initialized=.FALSE.
    END IF
  END SUBROUTINE qs_ks_dealloc_ref

!***************************************************************************

  SUBROUTINE init_grid(ks_env,qs_env, globenv,ngrid_level, error)
    TYPE(qs_ks_env_type), INTENT(inout) :: ks_env
    TYPE(qs_environment_type), INTENT(in), TARGET :: qs_env
    TYPE(global_environment_type), INTENT(in), TARGET :: globenv
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error
    INTEGER, INTENT(in) :: ngrid_level

!   *** Local parameters ***

    CHARACTER(LEN=*), PARAMETER :: routineN ="init_grid",&
         routineP=routineN//':'//moduleN

!   *** Local variables ***

    INTEGER :: handle,i,igrid_level,istat,output_unit ! ,grid_info
    INTEGER :: nspin, auxbas
    LOGICAL :: ionode
    TYPE(cell_type), POINTER          :: cell
    TYPE(pw_grid_type), DIMENSION(:), POINTER :: pw_grid
    TYPE(dft_control_type), POINTER :: dft_control
    TYPE(pw_env_type), POINTER :: pw_env
    TYPE(pw_pool_p_type), DIMENSION(:), POINTER :: pw_pools
    REAL(kind=wp), DIMENSION(:), POINTER :: cutoff

!   ---------------------------------------------------------------------------

    CALL timeset("init_grid","I","",handle)

    ionode = globenv%ionode
    output_unit = globenv%scr
    CALL get_qs_env(qs_env,cell=cell,dft_control=dft_control, pw_env=pw_env)
    CALL pw_env_get(pw_env, pw_pools=pw_pools)
    cutoff => dft_control%qs_control%e_cutoff

    ! the density grid coincideds with one of the multi grids
    ! maybe this could be generalised in the future
    DO i=1,ngrid_level
      IF ( cutoff(i) == dft_control%qs_control%cutoff ) ks_env % auxbas_grid = i
    END DO

    ! to be changed
    nspin = 1
    ALLOCATE (ks_env%rho_gspace(nspin),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"rho_gspace",0)
    ALLOCATE (ks_env%v_rspace(nspin),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"v_rspace",0)
    ALLOCATE (ks_env%v_gspace(nspin),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"v_gspace",0)

!   *** Allocate the multi-grid arrays ***

    ALLOCATE (ks_env%pw_grid(ngrid_level),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"pw_grid",0)
    ALLOCATE (ks_env%mgrid_rspace(ngrid_level),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"mgrid_rspace",0)
    ALLOCATE (ks_env%mgrid_gspace(ngrid_level),STAT=istat)
    IF (istat /= 0) CALL stop_memory(routineP,"mgrid_gspace",0)

!   *** Construct the various coeff_type grids ***

    pw_grid => ks_env%pw_grid
    DO igrid_level=1,ngrid_level
       pw_grid(igrid_level)=pw_pools(igrid_level)%pool%pw_grid

       ! setup the multigrids
       CALL coeff_allocate(ks_env%mgrid_rspace(igrid_level),pw_grid(igrid_level),&
            REALDATA3D)
       ks_env%mgrid_rspace(igrid_level)%pw%in_space = REALSPACE
       CALL coeff_allocate(ks_env%mgrid_gspace(igrid_level),pw_grid(igrid_level),&
            COMPLEXDATA1D)
       ks_env%mgrid_gspace(igrid_level)%pw%in_space = RECIPROCALSPACE

    END DO

!   *** Following quantities are only needed for the auxillary cutoff grid ***

    auxbas = ks_env % auxbas_grid 

    CALL pw_green_fn_init(cell%green,pw_grid(auxbas))

    CALL init_pw_poisson_solver(cell%green)

    CALL coeff_allocate(ks_env%rho_gspace(1),pw_grid(auxbas),COMPLEXDATA1D)
    ks_env%rho_gspace(1)%pw%in_space = RECIPROCALSPACE
    CALL coeff_allocate(ks_env%v_rspace(1),pw_grid(auxbas),REALDATA3D)
    ks_env%v_rspace(1)%pw%in_space = REALSPACE
    CALL coeff_allocate(ks_env%v_gspace(1),pw_grid(auxbas),COMPLEXDATA1D)
    ks_env%v_gspace(1)%pw%in_space = RECIPROCALSPACE

    CALL coeff_allocate(ks_env%v_hartree_rspace,pw_grid(auxbas),REALDATA3D)
    ks_env%v_hartree_rspace%pw%in_space = REALSPACE
    CALL coeff_allocate(ks_env%v_rspace_old,pw_grid(auxbas),REALDATA3D)
    ks_env%v_rspace_old%pw%in_space = REALSPACE
    CALL coeff_allocate(ks_env%rho_core,pw_grid(auxbas),COMPLEXDATA1D)
    ks_env%rho_core%pw%in_space = RECIPROCALSPACE

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE init_grid

! *****************************************************************************

  !!****f* qs_build_ks_matrix/qs_ks_build_kohn_sham_matrix [1.0] *
  !!
  !!   NAME
  !!     qs_ks_build_kohn_sham_matrix
  !!
  !!   FUNCTION
  !!     private routine where the real calculations are made: the
  !!     KS matrix is calculated
  !!
  !!   NOTES
  !!     -
  !!
  !!   INPUTS
  !!     - ks_env: the ks_env that old all the temporary objects that
  !!       the calculation of the KS matrix needs
  !!     - qs_env: the qs_env to update
  !!     - global_env: the global environment
  !!     - rebuild: if true builds again th
  !!     - calculate_forces: if true calculate the quantities needed
  !!       to calculate the forces. Defaults to false.
  !!     - just_energy: if true updates the energies but not the
  !!       ks matrix. Defaults to false
  !!     - ks: the KS matrix that is created
  !!     - p: the actual density matrix, if not given it uses rho%rho_ao
  !!     - rho: the density, rho%rho_r and rho%rho_g must be valid
  !!       if p is not vaild then rho%rho_ao must also be valid
  !!     - energy: the place where energies are stored
  !!     - error: variable to control error logging, stopping,... 
  !!       see module cp_error_handling 
  !!
  !!   AUTHOR
  !!     Matthias Krack
  !!
  !!   MODIFICATION HISTORY
  !!     6.2002 moved from qs_scf to qs_build_ks_matrix, use of ks_env
  !!            new did_change scheme [fawzi]
  !!
  !!*** *********************************************************************
  SUBROUTINE qs_ks_build_kohn_sham_matrix(ks_env,qs_env,globenv,ks,p,&
       rho,energy,rebuild,calculate_forces,just_energy,error)
    TYPE(qs_ks_env_type), INTENT(inout), TARGET :: ks_env
    TYPE(global_environment_type), INTENT(IN) :: globenv
    TYPE(qs_environment_type), INTENT(INOUT)  :: qs_env
    TYPE(real_matrix_p_type), DIMENSION(:), INTENT(inout) :: ks
    TYPE(real_matrix_p_type), DIMENSION(:), INTENT(in) :: p
    TYPE(qs_energy_type), INTENT(inout), TARGET :: energy
    LOGICAL, INTENT(in) :: calculate_forces,rebuild, just_energy
    TYPE(cp_error_type), INTENT(inout), OPTIONAL :: error
    TYPE(qs_rho_type), POINTER :: rho

    !   *** Local variables ***

    REAL(wp)          :: total_rho_elec_rspace
    CHARACTER(len=*), PARAMETER :: routineN='qs_ks_build_kohn_sham_matrix',&
         routineP=moduleN//':'//routineN
    INTEGER :: handle,output_unit
    LOGICAL :: ionode, failure,r_build
    TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: h,s
    TYPE(dft_control_type), POINTER :: dft_control
    TYPE(atomic_kind_type), DIMENSION(:), POINTER :: atomic_kind_set
    TYPE(pw_pool_p_type), DIMENSION(:), POINTER :: pw_pools
    TYPE(pw_env_type), POINTER :: pw_env
    TYPE(qs_rho_type), POINTER :: my_rho
    failure=.FALSE.
    r_build=rebuild
    NULLIFY(my_rho,pw_env,pw_pools,atomic_kind_set,dft_control,h,s)
    !   -----------------------------------------------------------------------

    CALL timeset(routineN,"I","",handle)

    CPPrecondition(ks_env%initialized,cp_failure_level,routineP,error,failure)

    ionode = globenv%ionode
    output_unit = globenv%scr

    CALL get_qs_env(qs_env=qs_env,&
                    dft_control=dft_control,&
                    h=h,&
                    s=s,&
                    pw_env=pw_env,&
                    rho=my_rho)
    my_rho => rho
    CPPrecondition(ASSOCIATED(h),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(s),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(my_rho),cp_failure_level,routineP,error,failure)
    CPPrecondition(ASSOCIATED(pw_env),cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(ks)>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(SIZE(p)>0,cp_failure_level,routineP,error,failure)

    IF (.NOT.failure) THEN
       IF (ks_env%s_mstruct_changed) THEN
          ! update qs_env%k structure?
          
          IF (ASSOCIATED(ks_env%v%matrix)) THEN
             CALL deallocate_matrix(ks_env%v%matrix)
          END IF
          NULLIFY (ks_env%v%matrix)
          CALL replicate_matrix_structure(h(1)%matrix,ks_env%v%matrix,&
               "POTENTIAL MATRIX")
          r_build=.TRUE.

          CALL calculate_rho_core(ks_env%rho_core,pw_env%cube_info(1),&
               pw_env%l_info,&
               ks_env%total_rho_core_rspace,qs_env,error=error)
          ks_env%s_mstruct_changed=.FALSE.
        END IF

! to rm
       CALL calculate_ecore(h(1)%matrix,p(1)%matrix,s(1)%matrix,energy%core,&
                            globenv)

       CALL calculate_rho_elec(p=p(1)%matrix,p_old=my_rho%rho_ao(1)%matrix,&
            rebuild=r_build,&
            rho=my_rho%rho_r(1),&
            rho_gspace=ks_env%rho_gspace(1),&
            rho_gspace_old=my_rho%rho_g(1),&
            mgrid_rspace=ks_env%mgrid_rspace,mgrid_gspace=ks_env%mgrid_gspace,&
            gridlevel_info=pw_env%gridlevel_info,cube_info=pw_env%cube_info,&
            l_info=pw_env%l_info,&
            total_rho=total_rho_elec_rspace,qs_env=qs_env,error=error)

! has this to be done here???
!deb   CALL qs_rho_did_change(my_rho,rho_r_valid=.TRUE.,error=error)
       my_rho%rho_r_valid=.TRUE.
! to rm

       CALL coeff_zero(ks_env%v_rspace(1))

       ks_env%total_rho_rspace=total_rho_elec_rspace+&
            ks_env%total_rho_core_rspace

       CALL calculate_xc_potential(ks_env=ks_env,&
            qs_env=qs_env,&
            globenv=globenv)
       
       ! add the core density to the electronic density
       CALL coeff_sumup(ks_env%rho_core,ks_env%rho_gspace(1))

       ks_env%total_rho_gspace = calculate_total_rho(ks_env%rho_gspace(1))

       IF (ionode.AND.globenv%print%total_densities) THEN
          WRITE (UNIT=output_unit,FMT="(/,(T3,A,T40,2F20.10))")&
               "Total electronic density (r-space): ",&
               total_rho_elec_rspace,total_rho_elec_rspace + &
               REAL(ks_env%n_electrons,wp),&
               "Total core charge density (r-space):",&
               ks_env%total_rho_core_rspace,ks_env%total_rho_core_rspace - &
               REAL(ks_env%n_electrons,wp)
          WRITE (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
               "Total charge density (r-space):     ",ks_env%total_rho_rspace,&
               "Total charge density (g-space):     ",ks_env%total_rho_gspace
       END IF

       CALL calculate_hartree(ks_env%rho_gspace(1),energy%hartree,&
            ks_env%v_gspace(1))

       energy%total = energy%core_overlap + energy%core_self + &
            energy%core + energy%hartree + energy%exc

       IF (ionode.AND.globenv%print%scf_energies) THEN
          WRITE (UNIT=output_unit,FMT="(/,(T3,A,T60,F20.10))")&
            "Core Hamiltonian energy:                       ",energy%core,&
            "Hartree energy:                                ",energy%hartree,&
            "Exchange-correlation energy:                   ",energy%exc
       END IF

       IF (.NOT. just_energy) THEN

          CALL coeff_transform_space(ks_env%v_gspace(1),ks_env%v_hartree_rspace)

          ks_env%v_hartree_rspace%pw%cr3d(:,:,:) =&
               ks_env%v_hartree_rspace%pw%pw_grid%dvol*&
               ks_env%v_hartree_rspace%pw%cr3d(:,:,:)

          IF (calculate_forces) THEN
             CALL integrate_v_core_rspace(ks_env%v_hartree_rspace,&
                  pw_env%cube_info(1),&
                  pw_env%l_info,qs_env)
          END IF

          ks_env%v_rspace(1)%pw%cr3d(:,:,:) =&
               ks_env%v_rspace(1)%pw%pw_grid%dvol*&
               ks_env%v_rspace(1)%pw%cr3d(:,:,:) +&
               ks_env%v_hartree_rspace%pw%cr3d(:,:,:)

          IF (r_build) CALL set_matrix(ks_env%v%matrix,0.0_wp)
          CALL integrate_v_rspace(ks_env%v_rspace(1),ks_env%mgrid_rspace,&
               ks_env%mgrid_gspace, ks_env%auxbas_grid,&
               pw_env%gridlevel_info,ks_env%v_rspace_old,&
               pw_env%cube_info,pw_env%l_info,r_build,p(1),ks_env%v,&
               qs_env,globenv,calculate_forces)
          CALL add_matrices(ks(1)%matrix,1.0_wp,h(1)%matrix,1.0_wp,&
               ks_env%v%matrix)

       END IF

    END IF

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE qs_ks_build_kohn_sham_matrix
! *****************************************************************************

!!****f* qs_build_ks_matrix/qs_ks_update_qs_env [1.0] *
!!
!!   NAME
!!     qs_ks_update_qs_env
!!
!!   FUNCTION
!!     updates the Kohn Sham matrix of the given qs_env (facility method)
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - ks_env: the ks_env that old all the temporary objects that
!!       the calculation of the KS matrix needs
!!     - qs_env: the qs_env to update
!!     - global_env: the global environment
!!     - rebuild: if true builds again th
!!     - calculate_forces: if true calculate the quantities needed
!!       to calculate the forces. Defaults to false.
!!     - just_energy: if true updates the energies but not the
!!       ks matrix. Defaults to false
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE qs_ks_update_qs_env(ks_env,qs_env, global_env,&
     rebuild,calculate_forces, just_energy, error)
  TYPE(qs_ks_env_type), INTENT(inout) :: ks_env
  TYPE(qs_environment_type), INTENT(inout), TARGET :: qs_env
  TYPE(global_environment_type), INTENT(in), TARGET :: global_env
  LOGICAL, INTENT(in), OPTIONAL :: rebuild, calculate_forces, just_energy
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure, r_build, c_forces, energy_only, did_update
  CHARACTER(len=*), PARAMETER :: routineN='qs_ks_update_qs_env',&
        routineP=moduleN//':'//routineN
  TYPE(real_matrix_p_type), DIMENSION(:), POINTER :: ks,p
  TYPE(qs_energy_type), POINTER :: energy
  TYPE(cp_logger_type), POINTER :: logger
  TYPE(qs_rho_type), POINTER :: rho
  failure=.FALSE.
  did_update=.FALSE.
  logger => cp_error_get_logger(error)
  r_build=.FALSE.; c_forces=.FALSE.; energy_only = .FALSE.
  IF (PRESENT(just_energy)) energy_only=just_energy
  IF (PRESENT(calculate_forces)) c_forces=calculate_forces
  IF (PRESENT(rebuild)) r_build=rebuild

  CPPrecondition(ks_env%initialized,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (ks_env%p_changed .OR. ks_env%s_mstruct_changed .OR.&
          (c_forces.AND..NOT.ks_env%forces_up_to_date)) THEN
        did_update=.TRUE.
        CALL get_qs_env(qs_env,k=ks,p=p,rho=rho,energy=energy)
        CALL qs_ks_build_kohn_sham_matrix(ks_env,qs_env=qs_env,&
             globenv=global_env,ks=ks,p=p,&
             rho=rho,energy=energy,rebuild=r_build,&
             calculate_forces=c_forces,just_energy=energy_only,&
             error=error)

        !   *** Keep a copy of the old density matrix ***
        CALL transfer_matrix(p(1)%matrix,rho%rho_ao(1)%matrix)
        
        IF(.NOT.energy_only) THEN
           ks_env%p_changed=.FALSE.
           ks_env%s_mstruct_changed=.FALSE.
           ks_env%forces_up_to_date=.NOT.c_forces
        END IF
     END IF
  END IF
  IF (cp_debug.AND.debug_this_module) THEN
     IF (did_update) THEN
        CALL cp_log(logger,cp_note_level,routineP,&
             "did update")
     ELSE
        CALL cp_log(logger,cp_note_level,routineP,&
             "did NOT update")
     END IF
  END IF
END SUBROUTINE qs_ks_update_qs_env
!***************************************************************************

!!****f* qs_build_ks_matrix/qs_ks_did_change [1.0] *
!!
!!   NAME
!!     qs_ks_did_change
!!
!!   FUNCTION
!!     tells that some of the things relevant to the ks calculation
!!     did change. has to be called when changes happen otherwise
!!     the calculation will give wrong results.
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     - ks_env: the environement that is informed about the changes
!!     - qs_env: the qs_environement
!!     - global_env: the global environment
!!     - s_mstruct_changed: if true it means that the structure of the
!!       overlap matrix has changed
!!       (atoms have moved)
!!     - p_changed: if true it means that the density has changed
!!     - full_reset: if true everything has changed
!!     - error: variable to control error logging, stopping,... 
!!       see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     4.2002 created [fawzi]
!!
!!*** **********************************************************************
SUBROUTINE qs_ks_did_change(ks_env,qs_env,global_env,s_mstruct_changed,&
     p_changed,full_reset,error)
  TYPE(qs_ks_env_type), INTENT(inout) :: ks_env
  TYPE(qs_environment_type), INTENT(inout), TARGET :: qs_env
  TYPE(global_environment_type), INTENT(in), TARGET :: global_env
  LOGICAL, INTENT(in), OPTIONAL :: s_mstruct_changed, p_changed, full_reset
  TYPE(cp_error_type), OPTIONAL, INTENT(inout) :: error
  
  LOGICAL :: failure
  CHARACTER(len=*), PARAMETER :: routineN='qs_ks_did_change',&
        routineP=moduleN//':'//routineN
  failure=.FALSE.
  
  CPPrecondition(ks_env%initialized,cp_failure_level,routineP,error,failure)
  IF (.NOT. failure) THEN
     IF (PRESENT(p_changed)) THEN
        IF (p_changed) ks_env%p_changed=.TRUE.
     END IF
     IF (PRESENT(s_mstruct_changed)) THEN
        IF (s_mstruct_changed) THEN
           ks_env%s_mstruct_changed=.TRUE.
        END IF
     END IF
     IF (PRESENT(full_reset)) THEN
        IF (full_reset) THEN
           ks_env%s_mstruct_changed=.TRUE.
        END IF
     END IF
  END IF
END SUBROUTINE qs_ks_did_change
!***************************************************************************

!!****f* qs_build_ks_matrix/calculate_xc_potential [1.0] *
!!
!!   NAME
!!     calculate_xc_potential
!!
!!   FUNCTION
!!     -
!!
!!   NOTES
!!     -
!!
!!   INPUTS
!!     error: variable to control error logging, stopping,... 
!!            see module cp_error_handling 
!!
!!   AUTHOR
!!     Fawzi Mohamed
!!
!!   MODIFICATION HISTORY
!!     - 05.2002 modified to use the mp_allgather function each pe
!!       computes only part of the grid and this is broadcasted to all 
!!       instead of summed. 
!!       This scales significantly better (e.g. factor 3 on 12 cpus
!!       32 H2O) [Joost VdV]
!!     - moved to qs_build_ks_matrix [fawzi]
!!
!!*** *********************************************************************
  SUBROUTINE calculate_xc_potential(ks_env,qs_env,globenv)
    IMPLICIT NONE
    TYPE(qs_ks_env_type), INTENT(inout),TARGET :: ks_env
    TYPE(qs_environment_type), INTENT(in), TARGET :: qs_env
    TYPE(global_environment_type), INTENT(IN), TARGET :: globenv

!   *** Local variables ***

    INTEGER :: handle

    INTEGER, DIMENSION(3) :: ng

    REAL(wp), DIMENSION(:,:,:), POINTER :: drho,rho,vxc
    TYPE(dft_control_type), POINTER :: dft_control
    TYPE(qs_rho_type), POINTER :: rho_struct
    TYPE(qs_energy_type), POINTER :: energy
    TYPE(cp_error_type) :: error

    REAL(wp) :: ec,ex
    LOGICAL :: failure
    INTEGER :: i
    CHARACTER(len=*), PARAMETER :: routineN='calculate_xc_potential',&
         routineP=moduleN//':'//routineN

!   ---------------------------------------------------------------------------

    CALL timeset( routineN ,"I","",handle)

    failure=.FALSE.

    CALL get_qs_env(qs_env=qs_env,&
                    dft_control=dft_control,&
                    rho=rho_struct,&
                    energy=energy)

    ! test if the real space density is available
    CPPrecondition(ASSOCIATED(rho_struct),cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_struct%ref_count>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(rho_struct%rho_r_valid,cp_failure_level,routineP,error,failure)
    IF ( .NOT. failure ) THEN

      ! for gradient corrected functional we need the density in g space
      IF ( qs_env%dft_control%gradient_functional ) THEN
         IF ( .NOT. rho_struct%rho_g_valid ) THEN
            DO i=1,SIZE(rho_struct%rho_r)
               call coeff_transform_space(rho_struct%rho_r(i),rho_struct%rho_g(i))
            END DO
         END IF
      END IF

      IF ( qs_env%dft_control%nspins == 1) THEN
         call xc_calculate (rho_struct%rho_r(1)%pw,rho_struct%rho_g(1)%pw,&
                            ks_env%v_rspace(1)%pw,&
                            energy%exc,qs_env%dft_control)
         energy%exc = rho_struct%rho_r(1)%pw%pw_grid%dvol*energy%exc
      ELSEIF ( qs_env%dft_control%nspins == 2) THEN
         call stop_program(routineP,"NA")
      ELSE
         call stop_program(routineP,"nspins > 2")
      END IF

    END IF

    CALL timestop(0.0_wp,handle)

  END SUBROUTINE calculate_xc_potential

! *****************************************************************************

END MODULE qs_build_ks_matrix

! *****************************************************************************
