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
!!
!!   SOURCE
!****************************************************************************
module qs_build_ks_matrix
  use cp_log_handling, only: cp_fatal_level, cp_failure_level,&
       cp_warning_level, cp_note_level, cp_to_string, cp_log
  use cp_error_handling, only: cp_debug, cp_error_type, cp_error_init,&
       cp_error_dealloc_ref, cp_error_message, cp_assert,&
       cp_assertion_failed, cp_internal_error, cp_unimplemented_error,&
       cp_error_get_logger
  use kinds, only: wp=>dp
  use global_types, only: global_environment_type,DEBUG
  use timings, only: timeset, timestop
  USE qs_energy_types, ONLY: qs_energy_type
  use qs_environment_types, only: qs_environment_type, get_qs_env,&
       set_qs_env
  use coefficient_types,        only: coeff_allocate,&
       coeff_copy,&
       coeff_transform_space,&
       coeff_type,&
       coeff_zero,&
       coeff_sumup
  use cube_utils,               only: cube_info_type,&
       destroy_cube_info,&
       init_cube_info
  use pw_grid_types,            only: HALFSPACE,pw_grid_type
  use pw_grids,                 only: pw_find_cutoff,&
       pw_grid_construct,&
       pw_grid_setup
  use pw_types,                 only: COMPLEXDATA1D,COMPLEXDATA3D,&
       REALDATA3D,REALSPACE,&
       RECIPROCALSPACE,&
       pw_type
  use pws,                      only: init_pw_poisson_solver,&
       pw_poisson_solver
  use greens_fn,                only: pw_green_fn_init
!JH
  use functionals
  use xc,                       only: xc_calculate
!JH
  use gaussian_gridlevels,      only: gridlevel_info_type,&
       init_gaussian_gridlevel,& 
       destroy_gaussian_gridlevel
  use hartree,                  only: calculate_hartree
  use l_utils,                  only: destroy_l_info,&
       init_l_info,&
       l_info_type
  use qs_collocate_density,     only: calculate_rho_core,&
       calculate_rho_elec,&
       calculate_total_rho
  use qs_core_energies,         only: calculate_ecore,&
       calculate_ecore_overlap,&
       calculate_ecore_self
  use dft_types,                only: dft_control_type
  use external_potential_types, only: all_potential_type,&
       get_potential,&
       gth_potential_type
  use message_passing,          only: mp_max,mp_sum,mp_sync,&
       mp_range,mp_allgather
  use qs_overlap,               only: write_blacs_matrix,&
       write_sparse_matrix
  use simulation_cell,          only: cell_type,&
       get_cell
  use sparse_matrix_types,      only: add_matrices,&
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
  use atomic_kind_types,        only: atomic_kind_type,&
       get_atomic_kind,&
       get_atomic_kind_set
  use qs_integrate_potential,   only: integrate_v_rspace, &
       integrate_v_core_rspace
  use termination, only: stop_memory, stop_program
  use qs_energy_types, only : qs_energy_type
  implicit none
  private

  logical, private, parameter :: debug_this_module=.true.
  character(len=*), private, parameter :: moduleN='qs_build_ks_matrix'

  public :: qs_ks_env_type
  public :: qs_ks_init, qs_ks_dealloc_ref, qs_ks_did_change,&
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
  type qs_ks_env_type
     type(pw_grid_type), dimension(:), pointer :: pw_grid
     logical :: initialized, s_struct_changed, p_changed, forces_up_to_date
     integer :: n_evals, n_electrons
     real(kind=wp) :: total_rho_core_rspace, total_rho_rspace, total_rho_gspace
     type(coeff_type), dimension(:), pointer       :: rho_gspace,rho_rspace,&
          rho_rspace_old,v_gspace,v_rspace
     type(coeff_type)          :: v_rspace_old,work_gspace,work_rspace
     type(gridlevel_info_type) :: gridlevel_info
     type(l_info_type)         :: l_info
     type(real_matrix_p_type)  :: v
     type(cube_info_type), dimension(:), pointer :: cube_info
     type(coeff_type)          :: rho_core
     type(coeff_type), dimension(3) :: nabla_rho_rspace
  end type qs_ks_env_type
!!***
!****************************************************************************

contains

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
  subroutine qs_ks_init(ks_env, qs_env, global_env, error)
    type(qs_ks_env_type), intent(out) :: ks_env
    type(qs_environment_type), intent(in),target :: qs_env
    type(global_environment_type), intent(in), target :: global_env
    type(cp_error_type), optional, intent(inout) :: error

    integer, parameter :: maxgridpoints = 100
    character(len=*), parameter :: routineN='qs_ks_init',&
         routineP=moduleN//':'//routineN
    logical :: failure
    real(wp)          :: rel_cutoff,maxradius
    integer :: ngrid_level
    real(wp), dimension(:), pointer :: cutoff
    integer  :: cmax,cmaxl,igrid_level, istat,maxl
    type(dft_control_type), pointer :: dft_control
    type(real_matrix_p_type), DIMENSION(:), pointer :: h
    type(atomic_kind_type), dimension(:), pointer :: atomic_kind_set
    failure=.false.

    nullify(ks_env%rho_rspace, ks_env%rho_gspace, ks_env%rho_rspace_old,&
         ks_env%v_gspace, ks_env%v_rspace, ks_env%cube_info,&
         ks_env%v%matrix)

    call init_functionals()
    call get_qs_env(qs_env=qs_env,&
         dft_control=dft_control,atomic_kind_set=atomic_kind_set)

    rel_cutoff = dft_control%qs_control%relative_cutoff
    cutoff => dft_control%qs_control%e_cutoff
    ngrid_level = size(cutoff)

    call init_gaussian_gridlevel(ks_env%gridlevel_info,&
         ngrid_level,cutoff,rel_cutoff)

    call init_grid(ks_env,qs_env=qs_env,globenv=global_env, &
         ngrid_level=ngrid_level,error=error)

    allocate (ks_env%cube_info(ngrid_level),STAT=istat)
    if (istat /= 0) call stop_memory(routineP,"cube_info",0)

    cmax = 0
    do igrid_level=1,ngrid_level
       maxradius = maxgridpoints*maxval(ks_env%rho_rspace(igrid_level)%pw%pw_grid%dr)
       call init_cube_info(ks_env%cube_info(igrid_level),&
            ks_env%rho_rspace(igrid_level)%pw%pw_grid%dr(:),&
            maxradius,cmaxl)
       cmax = max(cmax,cmaxl)
    end do

    call get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
         maxl=maxl, nelectron=ks_env%n_electrons)
    ks_env%n_electrons=ks_env%n_electrons-dft_control%charge

    ! maxl+1 to allow for the forces
    call init_l_info(ks_env%l_info,maxl+1,cmax)

    ks_env%total_rho_core_rspace=-1.0
    ks_env%total_rho_rspace=-1.0
    ks_env%total_rho_gspace=-1.0
    ks_env%s_struct_changed=.true.
    ks_env%p_changed=.true.
    ks_env%forces_up_to_date=.false.
    ks_env%n_evals=0
    ks_env%initialized=.true.
  end subroutine qs_ks_init
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
  subroutine qs_ks_dealloc_ref(ks_env,qs_env,global_env,error)
    type(qs_ks_env_type), intent(inout) :: ks_env
    type(qs_environment_type), intent(in),target :: qs_env
    type(global_environment_type), intent(in), target :: global_env
    type(cp_error_type), optional, intent(inout) :: error

    logical :: failure
    character(len=*), parameter :: routineN='qs_ks_dealloc_ref',&
         routineP=moduleN//':'//routineN
    integer :: igrid_level,istat
    failure=.false.

    CPPrecondition(ks_env%initialized,cp_failure_level,routineP,error,failure)
    if (.not. failure) then
       call cp_unimplemented_error(routineP, "grid deallocation not done in "//&
            CPSourceFileRef,&
            error=error, error_level=cp_warning_level)

       do igrid_level=1,size(ks_env%cube_info)
          call destroy_cube_info(ks_env%cube_info(igrid_level))
       end do

       deallocate (ks_env%cube_info,STAT=istat)
       if (istat /= 0) call stop_memory(routineP,"cube_info")

       call destroy_l_info(ks_env%l_info)
       call destroy_gaussian_gridlevel(ks_env%gridlevel_info,global_env)

       ks_env%initialized=.false.
    end if
  end subroutine qs_ks_dealloc_ref
!***************************************************************************


  subroutine init_grid(ks_env,qs_env, globenv,ngrid_level, error)
    type(qs_ks_env_type), intent(inout) :: ks_env
    type(qs_environment_type), intent(in), target :: qs_env
    type(global_environment_type), intent(in), target :: globenv
    type(cp_error_type), intent(inout), optional :: error
    integer, intent(in) :: ngrid_level

!   *** Local parameters ***

    character(LEN=*), parameter :: routineN ="init_grid",&
         routineP=routineN//':'//moduleN

!   *** Local variables ***

    integer :: handle,i,igrid_level,istat,output_unit ! ,grid_info
    logical :: ionode
    type(cell_type), pointer          :: cell
    type(pw_grid_type), dimension(:), pointer :: pw_grid
    type(dft_control_type), pointer :: dft_control
    real(kind=wp), dimension(:), pointer :: cutoff

!   ---------------------------------------------------------------------------

    call timeset("init_grid","I","",handle)

    ionode = globenv%ionode
    output_unit = globenv%scr
    call get_qs_env(qs_env,cell=cell,dft_control=dft_control)
    cutoff => dft_control%qs_control%e_cutoff

!   *** Allocate the multi-grid arrays ***

    allocate (ks_env%pw_grid(ngrid_level),STAT=istat)
    if (istat /= 0) call stop_memory(routineP,"pw_grid",0)
    allocate (ks_env%rho_rspace(ngrid_level),STAT=istat)
    if (istat /= 0) call stop_memory(routineP,"rho_rspace",0)
    allocate (ks_env%rho_gspace(ngrid_level),STAT=istat)
    if (istat /= 0) call stop_memory(routineP,"rho_gspace",0)
    allocate (ks_env%v_rspace(ngrid_level),STAT=istat)
    if (istat /= 0) call stop_memory(routineP,"v_rspace",0)
    allocate (ks_env%v_gspace(ngrid_level),STAT=istat)
    if (istat /= 0) call stop_memory(routineP,"v_gspace",0)

    allocate (ks_env%rho_rspace_old(ngrid_level),STAT=istat)
    if (istat /= 0) call stop_memory(routineP,"rho_rspace_old",0)

!   *** Construct the various coeff_type grids ***

    pw_grid => ks_env%pw_grid
    do igrid_level=1,ngrid_level
       call pw_grid_construct(pw_grid(igrid_level))
       pw_grid(igrid_level)%grid_span = HALFSPACE
       pw_grid(igrid_level)%para%rs_dims ( 1 ) = globenv % num_pe
       pw_grid(igrid_level)%para%rs_dims ( 2 ) = 1

       IF (igrid_level == 1) THEN
         IF (globenv%print%pw_grid_information) THEN
           CALL pw_grid_setup(cell,pw_grid(igrid_level),cutoff(igrid_level),&
                              pe_group=globenv%group,&
                              info=output_unit)
         ELSE
           CALL pw_grid_setup(cell,pw_grid(igrid_level),cutoff(igrid_level),&
                              pe_group=globenv%group)
         END IF
       ELSE
         IF (globenv%print%pw_grid_information) THEN
           CALL pw_grid_setup(cell,pw_grid(igrid_level),cutoff(igrid_level),&
                              pe_group=globenv%group,&
                              info=output_unit,&
                              ref_grid=pw_grid(igrid_level-1))
         ELSE
           CALL pw_grid_setup(cell,pw_grid(igrid_level),cutoff(igrid_level),&
                              pe_group=globenv%group,&
                              ref_grid=pw_grid(igrid_level-1))
         END IF
       END IF

       call coeff_allocate(ks_env%rho_rspace(igrid_level),pw_grid(igrid_level),&
            REALDATA3D)
       ks_env%rho_rspace(igrid_level)%pw%in_space = REALSPACE
       call coeff_allocate(ks_env%rho_gspace(igrid_level),pw_grid(igrid_level),&
            COMPLEXDATA1D)
       ks_env%rho_gspace(igrid_level)%pw%in_space = RECIPROCALSPACE
       call coeff_allocate(ks_env%rho_rspace_old(igrid_level),pw_grid(igrid_level),&
            REALDATA3D)
       ks_env%rho_rspace_old%pw%in_space = REALSPACE
       call coeff_allocate(ks_env%v_rspace(igrid_level),pw_grid(igrid_level),&
            REALDATA3D)
       ks_env%v_rspace(igrid_level)%pw%in_space = REALSPACE
       call coeff_allocate(ks_env%v_gspace(igrid_level),pw_grid(igrid_level),&
            COMPLEXDATA1D)
       ks_env%v_gspace(igrid_level)%pw%in_space = RECIPROCALSPACE

    end do

!   *** Following quantities are only needed for the highest cutoff grid ***

    call pw_green_fn_init(cell%green,pw_grid(1))

    call init_pw_poisson_solver(cell%green)

    call coeff_allocate(ks_env%work_rspace,pw_grid(1),REALDATA3D)
    ks_env%work_rspace%pw%in_space = REALSPACE
    call coeff_allocate(ks_env%work_gspace,pw_grid(1),COMPLEXDATA1D)
    ks_env%work_gspace%pw%in_space = RECIPROCALSPACE
    call coeff_allocate(ks_env%v_rspace_old,pw_grid(1),REALDATA3D)
    ks_env%v_rspace_old%pw%in_space = REALSPACE
    call coeff_allocate(ks_env%rho_core,pw_grid(1),COMPLEXDATA1D)
    ks_env%rho_core%pw%in_space = RECIPROCALSPACE

!JH if (dft_control%gradient_functional) then
!JH    do i=1,3
!JH       call coeff_allocate(ks_env%nabla_rho_rspace(i),pw_grid(1),REALDATA3D)
!JH       ks_env%nabla_rho_rspace(i)%pw%in_space = REALSPACE
!JH    end do
!JH    call coeff_allocate(ks_env%drho_rspace,pw_grid(1),REALDATA3D)
!JH    ks_env%drho_rspace%pw%in_space = REALSPACE
!JH    call coeff_allocate(ks_env%vxcg_rspace,pw_grid(1),REALDATA3D)
!JH    ks_env%vxcg_rspace%pw%in_space = REALSPACE
!JH end if

    call timestop(0.0_wp,handle)

  end subroutine init_grid

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
  !!     - p: the actual density matrix
  !!     - p_old: the old density matrix
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
  subroutine qs_ks_build_kohn_sham_matrix(ks_env,qs_env,globenv,ks,p,&
       p_old,energy,rebuild,&
       calculate_forces,just_energy,error)
    type(qs_ks_env_type), intent(inout), target :: ks_env
    type(global_environment_type), intent(IN) :: globenv
    type(qs_environment_type), intent(INOUT)  :: qs_env
    type(real_matrix_p_type), dimension(:), intent(inout) :: ks
    type(real_matrix_p_type), dimension(:), intent(in) :: p,p_old
    type(qs_energy_type), intent(inout), target :: energy
    logical, intent(in) :: calculate_forces,rebuild, just_energy
    type(cp_error_type), intent(inout), optional :: error

    !   *** Local variables ***

    real(wp)          :: total_rho_elec_rspace
    character(len=*), parameter :: routineN='qs_ks_build_kohn_sham_matrix',&
         routineP=moduleN//':'//routineN
    integer :: handle,output_unit
    logical :: ionode, failure,r_build
    type(real_matrix_p_type), DIMENSION(:), pointer :: h,s
    type(dft_control_type), pointer :: dft_control
    type(atomic_kind_type), dimension(:), pointer :: atomic_kind_set
    failure=.false.
    r_build=rebuild
    !   -----------------------------------------------------------------------

    call timeset(routineN,"I","",handle)

    CPPrecondition(ks_env%initialized,cp_failure_level,routineP,error,failure)

    ionode = globenv%ionode
    output_unit = globenv%scr

    CALL get_qs_env(qs_env=qs_env,&
                    dft_control=dft_control,&
                    h=h,&
                    s=s)

    CPPrecondition(associated(h),cp_failure_level,routineP,error,failure)
    CPPrecondition(associated(s),cp_failure_level,routineP,error,failure)
    CPPrecondition(size(ks)>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(size(p)>0,cp_failure_level,routineP,error,failure)
    CPPrecondition(size(p_old)>0,cp_failure_level,routineP,error,failure)

    if (.not.failure) then
       if (ks_env%s_struct_changed) then
          ! update qs_env%k structure?
          
          if (associated(ks_env%v%matrix)) then
             call deallocate_matrix(ks_env%v%matrix)
          end if
          nullify (ks_env%v%matrix)
          call replicate_matrix_structure(h(1)%matrix,ks_env%v%matrix,"POTENTIAL MATRIX")
          r_build=.true.

          CALL calculate_rho_core(ks_env%rho_core,ks_env%cube_info(1),&
               ks_env%l_info,&
               ks_env%total_rho_core_rspace,qs_env,globenv)
          ks_env%s_struct_changed=.false.
       end if

       CALL calculate_ecore(h(1)%matrix,p(1)%matrix,s(1)%matrix,energy%core,&
                            globenv)

       call calculate_rho_elec(p(1),p_old(1),r_build,ks_env%rho_rspace,&
            ks_env%rho_rspace_old,&
            ks_env%rho_gspace,ks_env%gridlevel_info,ks_env%cube_info,&
            ks_env%l_info,&
            total_rho_elec_rspace,qs_env,globenv)

       call coeff_zero(ks_env%v_rspace(1))

       ks_env%total_rho_rspace=total_rho_elec_rspace+&
            ks_env%total_rho_core_rspace
       if (dft_control%gradient_functional) then
          ! we have to wait for the general xc functions since we changed
          ! to a complex1D data structure for the density
          call stop_program("build_kohn","feature has been removed")
          !MK      CALL calculate_nabla_rho_elec(globenv)
          !MK      CALL coeff_zero(vxcg_rspace)
       end if

       call calculate_xc_potential(ks_env=ks_env,&
            qs_env=qs_env,&
            globenv=globenv)

       ! add the core density to the electronic density
       call coeff_sumup(ks_env%rho_core,ks_env%rho_gspace(1))

       ks_env%total_rho_gspace = calculate_total_rho(ks_env%rho_gspace(1))

       if (ionode.and.globenv%print%total_densities) then
          write (UNIT=output_unit,FMT="(/,(T3,A,T40,2F20.10))")&
               "Total electronic density (r-space): ",&
               total_rho_elec_rspace,total_rho_elec_rspace + &
               real(ks_env%n_electrons,wp),&
               "Total core charge density (r-space):",&
               ks_env%total_rho_core_rspace,ks_env%total_rho_core_rspace - &
               real(ks_env%n_electrons,wp)
          write (UNIT=output_unit,FMT="(T3,A,T60,F20.10)")&
               "Total charge density (r-space):     ",ks_env%total_rho_rspace,&
               "Total charge density (g-space):     ",ks_env%total_rho_gspace
       end if

       call calculate_hartree(ks_env%rho_gspace(1),energy%hartree,&
            ks_env%work_gspace)

       energy%total = energy%core_overlap + energy%core_self + &
            energy%core + energy%hartree + energy%exc

       if (ionode.and.globenv%print%scf_energies) then
          write (UNIT=output_unit,FMT="(/,(T3,A,T60,F20.10))")&
            "Core Hamiltonian energy:                       ",energy%core,&
            "Hartree energy:                                ",energy%hartree,&
            "Exchange-correlation energy:                   ",energy%exc
       end if

       if (.not. just_energy) then

          call coeff_transform_space(ks_env%work_gspace,ks_env%work_rspace)

          ! work_rspace is the hartree potential in rspace

          ks_env%work_rspace%pw%cr3d(:,:,:) =&
               ks_env%work_rspace%pw%pw_grid%dvol*&
               ks_env%work_rspace%pw%cr3d(:,:,:)

          IF (calculate_forces) THEN
             CALL integrate_v_core_rspace(ks_env%work_rspace,&
                  ks_env%cube_info(1),&
                  ks_env%l_info,qs_env)
          END IF

          ks_env%v_rspace(1)%pw%cr3d(:,:,:) =&
               ks_env%v_rspace(1)%pw%pw_grid%dvol*&
               ks_env%v_rspace(1)%pw%cr3d(:,:,:) +&
               ks_env%work_rspace%pw%cr3d(:,:,:)

          if (r_build) call set_matrix(ks_env%v%matrix,0.0_wp)
          call integrate_v_rspace(ks_env%v_rspace,ks_env%v_gspace,&
               ks_env%gridlevel_info,ks_env%v_rspace_old,&
               ks_env%cube_info,ks_env%l_info,r_build,p(1),ks_env%v,&
               qs_env,globenv,calculate_forces)
          call add_matrices(ks(1)%matrix,1.0_wp,h(1)%matrix,1.0_wp,&
               ks_env%v%matrix)

       end if

    end if

    call timestop(0.0_wp,handle)

  end subroutine qs_ks_build_kohn_sham_matrix
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
subroutine qs_ks_update_qs_env(ks_env,qs_env, global_env,&
     rebuild,calculate_forces, just_energy, error)
  type(qs_ks_env_type), intent(inout) :: ks_env
  type(qs_environment_type), intent(inout), target :: qs_env
  type(global_environment_type), intent(in), target :: global_env
  logical, intent(in), optional :: rebuild, calculate_forces, just_energy
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure, r_build, c_forces, energy_only, did_update
  character(len=*), parameter :: routineN='qs_ks_update_qs_env',&
        routineP=moduleN//':'//routineN
  type(real_matrix_p_type), dimension(:), pointer :: ks,p,p_old
  type(qs_energy_type), pointer :: energy
  failure=.false.
  did_update=.false.
  r_build=.false.; c_forces=.false.; energy_only = .false.
  if (present(just_energy)) energy_only=just_energy
  if (present(calculate_forces)) c_forces=calculate_forces
  if (present(rebuild)) r_build=rebuild

  CPPrecondition(ks_env%initialized,cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     if (ks_env%p_changed .or. ks_env%s_struct_changed .or.&
          (c_forces.and..not.ks_env%forces_up_to_date)) then
        did_update=.true.
        call get_qs_env(qs_env,k=ks,p=p,p_old=p_old,energy=energy)
        call qs_ks_build_kohn_sham_matrix(ks_env,qs_env=qs_env,&
             globenv=global_env,ks=ks,p=p,&
             p_old=p_old,energy=energy,rebuild=r_build,&
             calculate_forces=c_forces,just_energy=energy_only,&
             error=error)

        !   *** Keep a copy of the old density matrix ***
        call transfer_matrix(p(1)%matrix,p_old(1)%matrix)
        
        if(.not.energy_only) then
           ks_env%p_changed=.false.
           ks_env%s_struct_changed=.false.
           ks_env%forces_up_to_date=.not.c_forces
        end if
     end if
  end if
  if (cp_debug.and.debug_this_module) then
     if (did_update) then
        call cp_log(cp_error_get_logger(error),cp_note_level,routineP,&
             "did update")
     else
        call cp_log(cp_error_get_logger(error),cp_note_level,routineP,&
             "did NOT update")
     end if
  end if
end subroutine qs_ks_update_qs_env
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
!!     - s_struct_changed: if true it means that the structure of the
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
subroutine qs_ks_did_change(ks_env,qs_env,global_env,s_struct_changed,p_changed,&
     full_reset,error)
  type(qs_ks_env_type), intent(inout) :: ks_env
  type(qs_environment_type), intent(inout), target :: qs_env
  type(global_environment_type), intent(in), target :: global_env
  logical, intent(in), optional :: s_struct_changed, p_changed, full_reset
  type(cp_error_type), optional, intent(inout) :: error
  
  logical :: failure
  character(len=*), parameter :: routineN='qs_ks_did_change',&
        routineP=moduleN//':'//routineN
  failure=.false.
  
  CPPrecondition(ks_env%initialized,cp_failure_level,routineP,error,failure)
  if (.not. failure) then
     if (present(p_changed)) then
        if (p_changed) ks_env%p_changed=.true.
     end if
     if (present(s_struct_changed)) then
        if (s_struct_changed) then
           ks_env%p_changed=.true.
           ks_env%s_struct_changed=.true.
        end if
     end if
     if (present(full_reset)) then
        if (full_reset) then
           ks_env%p_changed=.true.
           ks_env%s_struct_changed=.true.
        end if
     end if
  end if
end subroutine qs_ks_did_change
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
  subroutine calculate_xc_potential(ks_env,qs_env,globenv)
    type(qs_ks_env_type), intent(inout),target :: ks_env
    type(qs_environment_type), intent(in), target :: qs_env
    type(global_environment_type), intent(IN), target :: globenv

!   *** Local variables ***

    integer :: handle

    integer, dimension(3) :: ng

    real(wp), dimension(:,:,:), pointer :: drho,rho,vxc,vxcg
    type(dft_control_type), pointer :: dft_control

    TYPE(qs_energy_type), POINTER :: energy

    REAL(wp) :: ec,ex

!   ---------------------------------------------------------------------------

    call timeset("calculate_xc_potential","I","",handle)

    call get_qs_env(qs_env=qs_env,&
                    dft_control=dft_control,&
                    energy=energy)

!old
if(0==0) then
    ng(:) = ks_env%rho_rspace(1)%pw%pw_grid%npts(3)
    rho => ks_env%rho_rspace(1)%pw%cr3d(:,:,:)
    vxc => ks_env%v_rspace(1)%pw%cr3d(:,:,:)
    ex = 0.0_wp
    ec = 0.0_wp
    call pade(rho(:,:,:),ex,vxc(:,:,:))
    call mp_sum(ex,globenv%group)
    call mp_sum(ec,globenv%group)
    energy%exc = ks_env%rho_rspace(1)%pw%pw_grid%dvol*(ex+ec)
else
!old
!   call xc_calculate ()
endif

    call timestop(0.0_wp,handle)

  end subroutine calculate_xc_potential

! *****************************************************************************

end module qs_build_ks_matrix

! *****************************************************************************
