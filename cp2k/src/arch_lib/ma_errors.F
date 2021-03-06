!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011 Christiane Ribeiro and the CP2K developers group!
!-----------------------------------------------------------------------------!
!> \brief   Machine Hardware support for cp2k
!           Interface to the HWLOC library
!> \author  Christiane Pousa Ribeiro
!> \date    2011-24-10
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - Created 2011-24-10

MODULE ma_errors

  USE ma_kinds

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ma_errors'

  PUBLIC :: ma_error
  PUBLIC :: ma_error_stop

  INTEGER, PARAMETER, PUBLIC  :: ma_error_allocation = 1
  INTEGER, PARAMETER, PUBLIC  :: ma_error_finalize = 2


CONTAINS

! ****************************************************************************
!> \brief To show the library erros
!> \var message         The error message
!> \var unit_num        The output console
! ****************************************************************************
  SUBROUTINE ma_error (message)
    CHARACTER(LEN=*), INTENT(IN), OPTIONAL   :: message

    INTEGER                                  :: unit_num

      unit_num = default_output_unit 

      IF (PRESENT(message)) THEN
       WRITE(unit_num,"(T2,A)") "MACHINE| "//TRIM(message)
       WRITE(unit_num,"(T2,A)") "MACHINE| Abnormal program termination"
      ENDIF
      STOP 1
  END SUBROUTINE ma_error

! ****************************************************************************
!> \brief To show the library erros
!> \var message         The error message
!> \var unit_num        The output console
! ****************************************************************************
  SUBROUTINE ma_error_stop (error_type)
    INTEGER                                  :: error_type

      SELECT CASE (error_type)
       CASE (ma_error_allocation)
        CALL ma_error("Memory allocation error")
       CASE (ma_error_finalize)
        CALL ma_error("Finalization of arch_lib can not be performed")
       END SELECT
      STOP 1
  END SUBROUTINE ma_error_stop



END MODULE ma_errors
