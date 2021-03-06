!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief 
!> \par History
!> \author  MI 07.2009
! *****************************************************************************
MODULE gle_system_types
  USE bibliography,                    ONLY: Ceriotti2009,&
                                             Ceriotti2009b,&
                                             cite_reference
  USE extended_system_types,           ONLY: create_map_info_type,&
                                             map_info_type,&
                                             release_map_info_type
  USE f77_blas
  USE input_section_types,             ONLY: section_vals_type,&
                                             section_vals_val_get
  USE kinds,                           ONLY: dp
  USE parallel_rng_types,              ONLY: GAUSSIAN,&
                                             create_rng_stream,&
                                             delete_rng_stream,&
                                             next_rng_seed,&
                                             rng_stream_type
  USE string_utilities,                ONLY: compress
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: gle_dealloc,&
            gle_init, gle_thermo_create,&
            gle_type, gle_thermo_type

!
  TYPE gle_thermo_type
     INTEGER                                 :: degrees_of_freedom
     REAL(KIND=dp)                           :: nkt, kin_energy, thermostat_energy
     REAL(KIND=dp), DIMENSION(:), POINTER    :: s
     TYPE(rng_stream_type), POINTER          :: gaussian_rng_stream
  END TYPE gle_thermo_type

! *****************************************************************************
  TYPE gle_type
     INTEGER                                 :: ndim
     INTEGER                                 :: glob_num_gle, loc_num_gle, region
     INTEGER, DIMENSION(:), POINTER          :: mal
     REAL(dp)                                :: temp, dt, dt_fact
     REAL(dp), POINTER                       :: gle_s(:,:), gle_t(:,:) 
     REAL(dp), POINTER                       :: a_mat(:,:), c_mat(:,:) 
     TYPE(gle_thermo_type), POINTER          :: nvt(:)
     TYPE(map_info_type), POINTER            :: map_info
  END TYPE gle_type

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'gle_system_types'

CONTAINS

! *****************************************************************************
!> \brief 
!> \param 
! *****************************************************************************
  SUBROUTINE gle_init(gle, dt,temp, section, error)
    TYPE(gle_type), POINTER                  :: gle
    REAL(dp), INTENT(IN)                     :: dt, temp
    TYPE(section_vals_type), POINTER         :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'gle_init', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, ir, istat, j, k, n_rep
    LOGICAL                                  :: failure
    REAL(dp), DIMENSION(:), POINTER          :: list
    REAL(KIND=dp)                            :: a_scale

    failure = .FALSE.
    NULLIFY(gle%nvt)
    NULLIFY(gle%gle_s)
    NULLIFY(gle%gle_t)
    NULLIFY(gle%map_info)
    gle%loc_num_gle=0
    gle%glob_num_gle=0
    gle%temp = temp
    gle%dt = dt*0.5_dp

    CALL cite_reference(Ceriotti2009)   
    CALL cite_reference(Ceriotti2009b)   
    CALL section_vals_val_get(section,"NDIM",i_val=gle%ndim,error=error)
    CALL section_vals_val_get(section,"A_SCALE",r_val=a_scale,error=error)

    ALLOCATE(gle%a_mat(gle%ndim,gle%ndim), STAT=istat)
    ALLOCATE(gle%c_mat(gle%ndim,gle%ndim), STAT=istat)
    ALLOCATE(gle%gle_s(gle%ndim,gle%ndim), STAT=istat)
    ALLOCATE(gle%gle_t(gle%ndim,gle%ndim), STAT=istat)

    CALL section_vals_val_get(section,"A_LIST",n_rep_val=n_rep, error=error)

    j = 1
    k = 1
    DO ir = 1,n_rep
      NULLIFY(list)
      CALL section_vals_val_get(section,"A_LIST",&
             i_rep_val=ir,r_vals=list,error=error)

      IF(ASSOCIATED(list)) THEN
         DO i = 1,SIZE(list)
            IF(j>gle%ndim) THEN
               CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                  routineP,"GLE: Too many elements in A_LIST",error,failure)
            END IF
            gle%a_mat(j,k) = list(i)
            k = k + 1
            IF(k>gle%ndim) THEN
              k = 1
              j = j + 1
            END IF
         END DO
      END IF
    END DO ! ir
    IF(j<gle%ndim+1) THEN
       CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
            routineP,"GLE: Too few elements in A_LIST",error,failure)
    END IF
    gle%a_mat = gle%a_mat * a_scale

    CALL section_vals_val_get(section,"C_LIST",n_rep_val=n_rep, error=error)
    IF(n_rep>0) THEN
      j = 1
      k = 1
      DO ir = 1,n_rep
        NULLIFY(list)
        CALL section_vals_val_get(section,"C_LIST",&
             i_rep_val=ir,r_vals=list,error=error)

        IF(ASSOCIATED(list)) THEN
          DO i = 1,SIZE(list)
            IF(j>gle%ndim) THEN
                 CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
                    routineP,"GLE: Too many elements in C_LIST",error,failure)
            END IF
            gle%c_mat(j,k) = list(i)
            k = k + 1
            IF(k>gle%ndim) THEN
                k = 1
                j = j + 1
            END IF
          END DO
        END IF
      END DO ! ir
      IF(j<gle%ndim+1) THEN
         CALL cp_assert(.FALSE.,cp_failure_level,cp_assertion_failed,&
            routineP,"GLE: Too few elements in C_LIST",error,failure)
      END IF
    ELSE
      gle%c_mat = 0.0_dp
      DO i = 1,gle%ndim
          gle%c_mat(i,i) =  gle%temp
      END DO
    END IF
    CALL create_map_info_type(gle%map_info, error)
  END SUBROUTINE gle_init

! *****************************************************************************
!> \brief 
!> \param 
! *****************************************************************************
  SUBROUTINE gle_thermo_create(gle, mal_size, error)
    TYPE(gle_type), POINTER                  :: gle
    INTEGER, INTENT(IN)                      :: mal_size
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'gle_thermo_create', &
      routineP = moduleN//':'//routineN

    CHARACTER(LEN=40)                        :: name
    INTEGER                                  :: i, ithermo, my_index, stat
    LOGICAL                                  :: failure
    REAL(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: seed
    REAL(KIND=dp), DIMENSION(3, 2)           :: initial_seed, my_seed

    failure = .FALSE.
    CPPrecondition(ASSOCIATED(gle),cp_fatal_level,routineP,error,failure)
    CPPrecondition(.NOT.ASSOCIATED(gle%nvt),cp_fatal_level,routineP,error,failure)
    
    ALLOCATE ( gle%nvt(gle%loc_num_gle),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    DO i = 1, gle%loc_num_gle
       NULLIFY(gle%nvt(i)%s)
       ALLOCATE(gle%nvt(i)%s(gle%ndim),stat=stat)
       CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
       gle%nvt(i)%kin_energy = 0.0_dp
       gle%nvt(i)%thermostat_energy = 0.0_dp
    END DO

    ALLOCATE ( gle%mal(mal_size),stat=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
    gle%mal(:) = 0

    ! Initialize the gaussian stream random number
    initial_seed = next_rng_seed(error=error)
    ALLOCATE (seed(3,2,gle%glob_num_gle),STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    seed(:,:,1) = initial_seed
    DO ithermo=2,gle%glob_num_gle
       seed(:,:,ithermo) = next_rng_seed(seed(:,:,ithermo-1),error=error)
    END DO

    ! Update initial seed
    initial_seed = next_rng_seed(seed(:,:,gle%glob_num_gle),error=error)
    DO ithermo = 1, gle%loc_num_gle
       NULLIFY(gle%nvt(ithermo)%gaussian_rng_stream)
       my_index = gle%map_info%index(ithermo)
       my_seed  = seed(:,:,my_index)
       WRITE (UNIT=name,FMT="(A,I8)") "Wiener process for Thermostat #",my_index
       CALL compress(name)
       CALL create_rng_stream(rng_stream=gle%nvt(ithermo)%gaussian_rng_stream,&
            name=name,distribution_type=GAUSSIAN, extended_precision=.TRUE.,&
            seed=my_seed,error=error)
    END DO

    DEALLOCATE (seed,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)
 
  END SUBROUTINE gle_thermo_create

! *****************************************************************************
!> \brief Deallocate type for GLE thermostat
!> \param error variable to control error logging, stopping,...
!>        see module cp_error_handling
! *****************************************************************************
  SUBROUTINE gle_dealloc ( gle, error )
    TYPE(gle_type), POINTER                  :: gle
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'gle_dealloc', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: i, stat
    LOGICAL                                  :: failure

    failure = .FALSE.

    IF (ASSOCIATED(gle)) THEN
       IF ( ASSOCIATED(gle%a_mat)) THEN
          DEALLOCATE (gle%a_mat, STAT=stat)
          CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       END IF
       IF ( ASSOCIATED(gle%c_mat)) THEN
          DEALLOCATE (gle%c_mat, STAT=stat)
          CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       END IF
       IF ( ASSOCIATED(gle%gle_t)) THEN
          DEALLOCATE (gle%gle_t, STAT=stat)
          CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       END IF
       IF ( ASSOCIATED(gle%gle_s)) THEN
          DEALLOCATE (gle%gle_s, STAT=stat)
          CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       END IF
       IF (ASSOCIATED (gle%nvt)) THEN
          DO i = 1,SIZE(gle%nvt)
            DEALLOCATE (gle%nvt(i)%s, STAT=stat)
            IF (ASSOCIATED(gle%nvt(i)%gaussian_rng_stream)) THEN    
              CALL delete_rng_stream(gle%nvt(i)%gaussian_rng_stream,error=error)
            END IF
            CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
          END DO
          DEALLOCATE (gle%nvt, STAT=stat)
          CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       ENDIF
       IF (ASSOCIATED (gle%mal)) THEN
          DEALLOCATE (gle%mal, STAT=stat)
          CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
       END IF

       CALL release_map_info_type(gle%map_info, error)
       DEALLOCATE (gle, STAT=stat)
       CPPrecondition(stat==0,cp_fatal_level,routineP,error,failure)
    ENDIF

  END SUBROUTINE gle_dealloc

END MODULE gle_system_types
