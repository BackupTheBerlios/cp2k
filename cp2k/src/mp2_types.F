!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief Types needed for MP2 calculations
!> \par History
!>       2011.05 created [Mauro Del Ben]
!> \author MDB
! *****************************************************************************
MODULE mp2_types
  USE hfx_types,                       ONLY: pair_list_element_type
  USE input_constants,                 ONLY: mp2_method_direct,&
                                             mp2_method_gpw,&
                                             mp2_method_laplace,&
                                             mp2_method_none
  USE kinds,                           ONLY: default_path_length,&
                                             dp
  USE timings,                         ONLY: timeset,&
                                             timestop
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mp2_types'

  PUBLIC :: mp2_type, mp2_laplace_type, mp2_direct_type, mp2_map_type, &
                                        mp2_map_type_1, &
                                        mp2_map_type_3, &
                                        mp2_method_none,&
                                        mp2_method_laplace,&
                                        mp2_method_direct,&
                                        mp2_method_gpw,&
                                        init_TShPSC_lmax

  PUBLIC :: mp2_env_create, mp2_env_release,mp2_biel_type,&
            pair_list_type_mp2,&
            mp2_potential_type, mp2_gpw_type

  INTEGER, SAVE                                         :: init_TShPSC_lmax = -1

! TYPE definitions

  TYPE mp2_map_type_1
     INTEGER :: iatom,iset,nsgpf
  END TYPE
 
  TYPE mp2_map_type
     INTEGER :: iatom,jatom,iset,jset,nsgpf
  END TYPE

  TYPE mp2_map_type_3
     INTEGER :: iatom,jatom,katom,iset,jset,kset,nsgpf
  END TYPE

  TYPE mp2_biel_type
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: index_table
  END TYPE mp2_biel_type

  TYPE mp2_laplace_type
       INTEGER       :: n_quadrature
       LOGICAL       :: dbcsr
       LOGICAL       :: only_coulomb
       LOGICAL       :: super_p
       LOGICAL       :: sos
       REAL(dp)      :: threshold
  END TYPE

  TYPE mp2_direct_type
    LOGICAL  :: big_send
  END TYPE

  TYPE mp2_potential_type
    INTEGER                                  :: potential_type     !! 1/r, TShPSC
    REAL(dp)                                 :: truncation_radius  !! truncation radius if cutoff potential in use
    CHARACTER(default_path_length)           :: filename
  END TYPE

  TYPE mp2_gpw_type
       REAL(KIND=dp)            :: eps_grid, eps_filter
       INTEGER                  :: print_level
  END TYPE mp2_gpw_type

  TYPE mp2_type
       REAL(KIND=dp)            :: energy
       LOGICAL                  :: compute_reference
       INTEGER                  :: method
!      TYPE(mp2_laplace_type)   :: laplace
       TYPE(mp2_direct_type)    :: direct_canonical
       TYPE(mp2_potential_type) :: potential_parameter
       TYPE(mp2_gpw_type)       :: mp2_gpw 
       REAL(dp) :: mp2_memory
       INTEGER  :: mp2_num_proc
  END TYPE

  TYPE pair_list_type_mp2
    TYPE(pair_list_element_type), DIMENSION(:), ALLOCATABLE :: elements
    INTEGER :: n_element
  END TYPE pair_list_type_mp2

CONTAINS

  SUBROUTINE mp2_env_release(mp2_env,error)
    TYPE(mp2_type), POINTER                  :: mp2_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_env_release', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CALL timeset(routineN,handle)

    CPPostcondition(ASSOCIATED(mp2_env),cp_failure_level,routineP,error,failure)

    DEALLOCATE(mp2_env,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE mp2_env_release

  SUBROUTINE mp2_env_create(mp2_env,error)
    TYPE(mp2_type), POINTER                  :: mp2_env
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(LEN=*), PARAMETER :: routineN = 'mp2_env_create', &
      routineP = moduleN//':'//routineN

    INTEGER                                  :: handle, stat
    LOGICAL                                  :: failure

    failure=.FALSE.

    CALL timeset(routineN,handle)

    CPPostcondition(.NOT.ASSOCIATED(mp2_env),cp_failure_level,routineP,error,failure)

    ALLOCATE(mp2_env,STAT=stat)
    CPPostcondition(stat==0,cp_failure_level,routineP,error,failure)

    CALL timestop(handle)

  END SUBROUTINE mp2_env_create


END MODULE mp2_types


