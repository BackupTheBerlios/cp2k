!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief input section for atomic properties
!> \par History
!>      07.2011 created
!> \author JHU
! *****************************************************************************
MODULE input_cp2k_atprop
  USE cp_output_handling,              ONLY: cp_print_key_section_create
  USE f77_blas
  USE input_constants
  USE input_keyword_types,             ONLY: keyword_create,&
                                             keyword_release,&
                                             keyword_type
  USE input_section_types,             ONLY: section_add_keyword,&
                                             section_add_subsection,&
                                             section_create,&
                                             section_release,&
                                             section_type
#include "cp_common_uses.h"

  IMPLICIT NONE
  PRIVATE

  LOGICAL, PRIVATE, PARAMETER :: debug_this_module=.TRUE.
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_atprop'

  PUBLIC :: create_atprop_section

CONTAINS

! *****************************************************************************
!> \brief Creates the ATOMIC section
!> \param section the section to create
!> \param error variable to control error logging, stopping,... 
!>        see module cp_error_handling 
!> \author JHU
! *****************************************************************************
  SUBROUTINE create_atprop_section(section,error)
    TYPE(section_type), POINTER              :: section
    TYPE(cp_error_type), INTENT(inout)       :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'create_atprop_section', &
      routineP = moduleN//':'//routineN

    LOGICAL                                  :: failure
    TYPE(keyword_type), POINTER              :: keyword
    TYPE(section_type), POINTER              :: print_key

    failure=.FALSE.

    CPPrecondition(.NOT.ASSOCIATED(section),cp_failure_level,routineP,error,failure)
    IF (.NOT. failure) THEN
       CALL section_create(section,name="ATOMIC",&
            description="Controls the calculation and printing of atomic properties ",&
            repeats=.FALSE., required=.FALSE.,&
            error=error)

       NULLIFY(keyword)

       CALL keyword_create(keyword, name="ENERGY",&
            description="Calculate atomic energies ",&
            usage="ENERGY {logical}",&
            repeats=.FALSE.,&
            required=.FALSE.,&
            n_var=1,&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)    

       CALL keyword_create(keyword, name="PRESSURE",&
            description="Calculate atomic pressure tensors ",&
            usage="PRESSURE {logical}",&
            repeats=.FALSE.,&
            required=.FALSE.,&
            n_var=1,&
            default_l_val=.FALSE.,&
            lone_keyword_l_val=.TRUE.,&
            error=error)
       CALL section_add_keyword(section,keyword,error=error)
       CALL keyword_release(keyword,error=error)

       NULLIFY(print_key)
       CALL cp_print_key_section_create(print_key,"PRINT",&
            description="Controls the printing of basic information during the run", &
            print_level=low_print_level,add_last=add_last_numeric,filename="__STD_OUT__",&
            error=error)

! additional print keywords
!      CALL keyword_create(keyword, name="APPEND",&
!           description="append the cube files when they already exist",&
!           default_l_val=.FALSE., lone_keyword_l_val=.TRUE.,error=error)
!      CALL section_add_keyword(print_key,keyword,error=error)
!      CALL keyword_release(keyword,error=error)

       CALL section_add_subsection(section, print_key, error=error)
       CALL section_release(print_key,error=error)

    END IF
  END SUBROUTINE create_atprop_section

END MODULE input_cp2k_atprop
