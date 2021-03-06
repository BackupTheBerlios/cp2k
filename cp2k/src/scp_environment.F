!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief initialize scp environment
!> \author CJM
! *****************************************************************************
MODULE scp_environment
  USE array_types,                     ONLY: array_i1d_obj,&
                                             array_new,&
                                             array_nullify,&
                                             array_release
  USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                             get_atomic_kind,&
                                             get_atomic_kind_set
  USE cp_dbcsr_interface,              ONLY: cp_dbcsr_copy,&
                                             cp_dbcsr_create,&
                                             cp_dbcsr_finalize,&
                                             cp_dbcsr_get_diag,&
                                             cp_dbcsr_get_info,&
                                             cp_dbcsr_init,&
                                             cp_dbcsr_set
  USE cp_dbcsr_operations,             ONLY: cp_dbcsr_add_block_node
  USE cp_dbcsr_types,                  ONLY: cp_dbcsr_type
  USE cp_para_types,                   ONLY: cp_para_env_type
  USE dbcsr_types,                     ONLY: dbcsr_distribution_obj,&
                                             dbcsr_repl_full,&
                                             dbcsr_type_symmetric
  USE dbcsr_util,                      ONLY: convert_offsets_to_sizes
  USE distribution_1d_types,           ONLY: distribution_1d_type
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                             section_vals_type
  USE kinds,                           ONLY: dp
  USE particle_types,                  ONLY: get_particle_set,&
                                             particle_type
  USE qs_environment_types,            ONLY: get_qs_env,&
                                             qs_environment_type
  USE scp_coeff_types,                 ONLY: aux_coeff_set_create,&
                                             aux_coeff_set_initialize,&
                                             aux_coeff_set_release,&
                                             aux_coeff_set_type
  USE scp_energy_types,                ONLY: allocate_scp_energy,&
                                             scp_energy_type
  USE scp_environment_methods,         ONLY: scp_qs_env_setup
  USE scp_environment_types,           ONLY: scp_environment_type,&
                                             set_scp_env
  USE scp_force_types,                 ONLY: allocate_scp_force,&
                                             init_scp_force,&
                                             scp_force_type
  USE semi_empirical_types,            ONLY: get_se_param,&
                                             semi_empirical_type
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: scp_dft_init, scp_nddo_init

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'scp_environment'


CONTAINS

! *****************************************************************************
!> \brief Initializes the SCP environment using the QS environment
! *****************************************************************************
  SUBROUTINE scp_dft_init ( scp_env, qs_env, error )

    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(qs_environment_type), OPTIONAL, &
      POINTER                                :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_dft_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, nkind, stat
    INTEGER, ALLOCATABLE, DIMENSION(:)       :: natom_of_kind
    LOGICAL                                  :: failure
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(aux_coeff_set_type), POINTER        :: aux_coeff_set
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(distribution_1d_type), POINTER      :: local_particles
    TYPE(scp_energy_type), POINTER           :: energy
    TYPE(scp_force_type), POINTER            :: force( : )
    TYPE(section_vals_type), POINTER         :: dft_section, input

    CALL timeset(routineN,handle)

    failure = .FALSE.
    ! Nullifying  pointers
    NULLIFY ( aux_coeff_set )
    NULLIFY ( energy, force )
    NULLIFY ( para_env )
    NULLIFY ( input )
    NULLIFY ( atomic_kind_set )
    NULLIFY ( local_particles )

    IF ( PRESENT ( qs_env ) ) THEN
      CALL get_qs_env (qs_env, para_env = para_env, input = input, &
                       atomic_kind_set = atomic_kind_set, &
                       local_particles = local_particles, error = error )

      dft_section =>  section_vals_get_subs_vals( input, "DFT", error = error )
      ! Set the input section of the SCP
      CALL set_scp_env(scp_env,input=dft_section,&
           error=error)

      ! Allocate and initialize coefs pointer
      CALL aux_coeff_set_create(aux_coeff_set=aux_coeff_set,atomic_kind_set=atomic_kind_set,&
           distribution=local_particles, error=error)
      CALL aux_coeff_set_initialize(aux_coeff_set,para_env,error)
      CALL set_scp_env (scp_env=scp_env,aux_coeff_set = aux_coeff_set, error=error)
      CALL aux_coeff_set_release(aux_coeff_set, error=error)

      ! Allocate the data structure for energies ***
      CALL allocate_scp_energy (energy)

      ! Allocate the data structure for energies ***
      nkind = SIZE(atomic_kind_set)
      ALLOCATE (natom_of_kind(nkind),STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,&
                             natom_of_kind=natom_of_kind)
      CALL allocate_scp_force(force,natom_of_kind)
      DEALLOCATE (natom_of_kind,STAT=stat)
      CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
      CALL init_scp_force(force)

      ! set the scp_env
      CALL set_scp_env (scp_env=scp_env,energy=energy,force=force,error=error)

      !Sets up pw, grids, densities
      CALL scp_qs_env_setup(scp_env,qs_env,error)
    END IF
    CALL timestop(handle)

  END SUBROUTINE scp_dft_init
! *****************************************************************************
!> \brief Initializes the SCP environment using the QS environment
! *****************************************************************************
  SUBROUTINE scp_nddo_init ( scp_env, qs_env, error )

    TYPE(scp_environment_type), POINTER      :: scp_env
    TYPE(qs_environment_type), POINTER       :: qs_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'scp_nddo_init', &
      routineP = moduleN//':'//routineN

    INTEGER :: handle, iatom, ii, ikind, ilist, jj, natom, natorb, &
      neighbor_list_id, nkind, nrow, nsgf_a, nsgf_b, stat
    INTEGER, ALLOCATABLE                     :: first_sgf_a( : ), &
                                                first_sgf_b( : ), &
                                                last_sgf_a( : ), &
                                                last_sgf_b( : )
    INTEGER, DIMENSION(:), POINTER           :: atom_list, cbs, rbs
    LOGICAL                                  :: failure
    REAL(KIND=dp)                            :: app, asp, ass
    REAL(KIND=dp), ALLOCATABLE, DIMENSION(:) :: tmp
    REAL(KIND=dp), DIMENSION(9, 9)           :: a_scp
    REAL(KIND=dp), POINTER                   :: ascp_diag( : ), block( :, : )
    TYPE(array_i1d_obj)                      :: col_blk_sizes, row_blk_sizes
    TYPE(atomic_kind_type), DIMENSION(:), &
      POINTER                                :: atomic_kind_set
    TYPE(atomic_kind_type), POINTER          :: atomic_kind
    TYPE(cp_dbcsr_type), POINTER             :: ks_scp, pscp
    TYPE(cp_para_env_type), POINTER          :: para_env
    TYPE(dbcsr_distribution_obj), POINTER    :: dbcsr_dist
    TYPE(particle_type), DIMENSION(:), &
      POINTER                                :: particle_set
    TYPE(scp_energy_type), POINTER           :: energy
    TYPE(semi_empirical_type), POINTER       :: se_kind

    CALL timeset(routineN,handle)

    failure = .FALSE.
    ! Nullifying  pointers
    NULLIFY ( energy )
    NULLIFY ( particle_set )
    NULLIFY ( para_env )
    NULLIFY ( atomic_kind_set, atomic_kind )
    NULLIFY ( pscp, ks_scp, block, ascp_diag, se_kind, atom_list )

    CALL get_qs_env (qs_env, para_env = para_env, &
                     atomic_kind_set = atomic_kind_set, &
                     particle_set = particle_set, &
                     neighbor_list_id=neighbor_list_id, &
                     dbcsr_dist=dbcsr_dist,&
                     error = error )

    ! Allocate the data structure for energies ***
    CALL allocate_scp_energy (energy)
    ! set the scp_env
    CALL set_scp_env (scp_env=scp_env,energy=energy,error=error)
    nkind = SIZE(atomic_kind_set)
    natom = SIZE(particle_set)

    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,nsgf=nsgf_a)
    CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set,nsgf=nsgf_b)

    ALLOCATE (first_sgf_a(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (last_sgf_a(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (first_sgf_b(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    ALLOCATE (last_sgf_b(natom),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL get_particle_set(particle_set=particle_set,&
                          first_sgf=first_sgf_a,&
                          last_sgf=last_sgf_a,&
                          error=error)
    CALL get_particle_set(particle_set=particle_set,&
                          first_sgf=first_sgf_b,&
                          last_sgf=last_sgf_b,&
                          error=error)

    ! prepare for allocation
    ALLOCATE (rbs(natom),cbs(natom), STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL convert_offsets_to_sizes (first_sgf_a, rbs, last_sgf_a)
    CALL convert_offsets_to_sizes (first_sgf_b, cbs, last_sgf_b)
    CALL array_nullify (row_blk_sizes)
    CALL array_nullify (col_blk_sizes)
    CALL array_new (row_blk_sizes, rbs, gift=.TRUE.)
    CALL array_new (col_blk_sizes, cbs, gift=.TRUE.)

    ALLOCATE(pscp)
    CALL cp_dbcsr_init(pscp,error=error)
    CALL cp_dbcsr_create(matrix=pscp, &
         name="PSCP", &
         dist=dbcsr_dist, matrix_type=dbcsr_type_symmetric,&
         row_blk_size=row_blk_sizes, col_blk_size=col_blk_sizes, &
         nblks=0, nze=0, mutable_work=.TRUE., replication_type=dbcsr_repl_full, &
         error=error)

    CALL array_release (row_blk_sizes)
    CALL array_release (col_blk_sizes)

    DO ikind=1,nkind
       atomic_kind => atomic_kind_set(ikind)
       CALL get_atomic_kind(atomic_kind=atomic_kind,&
                            se_parameter=se_kind, &
                            atom_list=atom_list)
       CALL get_se_param(se_kind, ass=ass,app=app,asp=asp,natorb=natorb)
       IF (natorb>4) THEN
          CALL cp_unimplemented_error(fromWhere=routineP, &
               message="SCP-NDDO not implemened for d-orbitals!", &
               error=error, error_level=cp_failure_level)
       END IF
       a_scp( 1  ,  1) = ass
       a_scp( 1  ,2:4) = asp
       a_scp( 2:4,  1) = asp
       a_scp( 2:4,2:4) = app

       ! Loop over the replicated list of particles, filling the sparse matrix structure to exploit
       DO ilist=1,SIZE ( atom_list )
          iatom = atom_list ( ilist )

          NULLIFY(block)
          CALL cp_dbcsr_add_block_node(matrix=pscp,&
                                    block_row=iatom,&
                                    block_col=iatom,&
                                    BLOCK=BLOCK,error=error)
          DO ii=1,SIZE(block,1)
             DO jj=1,SIZE(block,2)
                block(ii,jj) = a_scp (ii,jj) 
             END DO
          END DO
       END DO ! ilist
    END DO ! ikind 

    CALL cp_dbcsr_finalize(pscp,error=error)

    CALL cp_dbcsr_get_info ( pscp, nfullrows_total=nrow )

    ALLOCATE (ascp_diag(nrow),STAT=stat)
    ALLOCATE (tmp(nrow),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    CALL cp_dbcsr_get_diag( pscp, tmp, error=error )
    ascp_diag=tmp
    CALL set_scp_env ( scp_env=scp_env, ascp_diag=ascp_diag, error = error )
     

    CALL cp_dbcsr_set(pscp,0.0_dp,error=error)
    ALLOCATE(ks_scp)
    CALL cp_dbcsr_init(ks_scp, error=error)
    CALL cp_dbcsr_copy(ks_scp,pscp,'KS_SCP',error=error)

    CALL set_scp_env (scp_env=scp_env,ks_scp=ks_scp,pscp=pscp,error=error)

    DEALLOCATE (first_sgf_a,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (last_sgf_a,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (first_sgf_b,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DEALLOCATE (last_sgf_b,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE scp_nddo_init

END MODULE scp_environment

