!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \par History
!>      JGH (15-Mar-2001) : Update small grid when cell changes
!>                          with dg_grid_change
! *****************************************************************************
MODULE dgs
  USE cell_types,                      ONLY: cell_create,&
                                             cell_release,&
                                             cell_type,&
                                             get_cell_param,&
                                             init_cell
  USE f77_blas
  USE fft_tools,                       ONLY: BWFFT,&
                                             FFT_RADIX_CLOSEST,&
                                             FFT_RADIX_NEXT_ODD,&
                                             fft3d,&
                                             fft_radix_operations
  USE kinds,                           ONLY: dp,&
                                             dp_size
  USE pw_grid_info,                    ONLY: pw_find_cutoff
  USE pw_grid_types,                   ONLY: HALFSPACE,&
                                             pw_grid_type
  USE pw_grids,                        ONLY: pw_grid_change,&
                                             pw_grid_setup
  USE pw_types,                        ONLY: COMPLEXDATA3D,&
                                             REALDATA3D,&
                                             pw_p_type
  USE realspace_grid_types,            ONLY: realspace_grid_p_type,&
                                             realspace_grid_type
  USE structure_factors,               ONLY: structure_factor_evaluate
  USE termination,                     ONLY: stop_memory,&
                                             stop_program
#include "cp_common_uses.h"

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dgs'

  PUBLIC :: dg_get_patch, dg_get_delta
  PUBLIC :: dg_pme_grid_setup, dg_find_max_cutoff,  &
       dg_sum_patch, dg_sum_patch_force_3d, dg_sum_patch_force_1d,  &
       dg_get_strucfac, dg_grid_change, dg_find_cutoff

  INTERFACE dg_sum_patch
     MODULE PROCEDURE dg_sum_patch_coef, dg_sum_patch_arr
  END INTERFACE

  INTERFACE dg_sum_patch_force_3d
     MODULE PROCEDURE dg_sum_patch_force_coef_3d, dg_sum_patch_force_arr_3d
  END INTERFACE

  INTERFACE dg_sum_patch_force_1d
     MODULE PROCEDURE dg_sum_patch_force_coef_1d, dg_sum_patch_force_arr_1d
  END INTERFACE

  INTERFACE dg_get_patch
     MODULE PROCEDURE dg_get_patch_1, dg_get_patch_2
  END INTERFACE

  INTERFACE dg_add_patch
     MODULE PROCEDURE dg_add_patch_simple, dg_add_patch_folded
  END INTERFACE

  INTERFACE dg_int_patch_3d
     MODULE PROCEDURE dg_int_patch_simple_3d, dg_int_patch_folded_3d
  END INTERFACE

  INTERFACE dg_int_patch_1d
     MODULE PROCEDURE dg_int_patch_simple_1d, dg_int_patch_folded_1d
  END INTERFACE

CONTAINS

! *****************************************************************************
SUBROUTINE dg_pme_grid_setup ( box_b, npts_s, cutoff_radius, grid_s, grid_b, &
           grid_ref, rs_dims, iounit, fft_usage, error )

    TYPE(cell_type), POINTER                 :: box_b
    INTEGER, DIMENSION(:), INTENT(IN)        :: npts_s
    REAL(KIND=dp), INTENT(IN)                :: cutoff_radius
    TYPE(pw_grid_type), POINTER              :: grid_s, grid_b
    TYPE(pw_grid_type), INTENT(IN), OPTIONAL :: grid_ref
    INTEGER, DIMENSION(2), INTENT(in), &
      OPTIONAL                               :: rs_dims
    INTEGER, INTENT(IN), OPTIONAL            :: iounit
    LOGICAL, INTENT(IN), OPTIONAL            :: fft_usage
    TYPE(cp_error_type), INTENT(INOUT)       :: error

    INTEGER, DIMENSION(2, 3)                 :: bo
    REAL(KIND=dp)                            :: cutoff, ecut
    TYPE(cell_type), POINTER                 :: box_s, unit_box

  NULLIFY(box_s, unit_box)
  CALL cell_create(box_s,error=error)
  CALL cell_create(unit_box,error=error)

  CALL dg_find_cutoff ( box_b, npts_s, cutoff_radius, grid_s, grid_b, cutoff ,error=error)

  ecut = 0.5_dp * cutoff * cutoff
  bo = grid_b % bounds
  IF ( PRESENT ( grid_ref ) ) THEN
     CALL pw_grid_setup ( box_b, grid_b, bounds=bo, cutoff=ecut, spherical=.TRUE.,&
          ref_grid = grid_ref, rs_dims=rs_dims, iounit=iounit, fft_usage=fft_usage,&
          error=error)
  ELSE
     CALL pw_grid_setup ( box_b, grid_b, bounds=bo, cutoff=ecut, spherical=.TRUE.,&
          rs_dims=rs_dims, iounit=iounit, fft_usage=fft_usage, error=error)
  ENDIF

  CALL dg_find_basis ( grid_b % npts, box_b, unit_box )

  CALL dg_set_cell ( grid_s % npts, unit_box, box_s )

  bo = grid_s % bounds
  CALL pw_grid_setup ( box_s, grid_s, bounds=bo, cutoff=ecut, iounit=iounit, fft_usage=fft_usage,&
                       error=error)

  CALL cell_release(box_s,error=error)
  CALL cell_release(unit_box,error=error)
END SUBROUTINE dg_pme_grid_setup

! *****************************************************************************
SUBROUTINE dg_find_cutoff ( box_b, npts_s, cutoff_radius, grid_s,  &
                            grid_b, cutoff, error)

    TYPE(cell_type), POINTER                 :: box_b
    INTEGER, DIMENSION(:), INTENT(IN)        :: npts_s
    REAL(KIND=dp), INTENT(IN)                :: cutoff_radius
    TYPE(pw_grid_type), POINTER              :: grid_s, grid_b
    REAL(KIND=dp), INTENT(OUT)               :: cutoff
    TYPE(cp_error_type), INTENT(inout)       :: error

    INTEGER                                  :: nout( 3 )
    REAL(KIND=dp)                            :: cell_lengths( 3 ), dr( 3 )

  CALL fft_radix_operations ( npts_s ( 1 ), nout ( 1 ), &
       operation = FFT_RADIX_NEXT_ODD )
  CALL fft_radix_operations ( npts_s ( 1 ), nout ( 2 ), &
       operation = FFT_RADIX_NEXT_ODD )
  CALL fft_radix_operations ( npts_s ( 1 ), nout ( 3 ), &
       operation = FFT_RADIX_NEXT_ODD )

  CALL get_cell_param ( box_b, cell_lengths , error=error)
  CALL dg_get_spacing ( nout, cutoff_radius, dr )
  CALL dg_find_radix ( dr, cell_lengths, grid_b % npts )

! In-line code to set grid_b % npts = npts_s if necessary
  IF ( nout ( 1 ) > grid_b % npts ( 1 ) ) THEN
    grid_b % npts ( 1 ) = nout ( 1 )
    dr ( 1 ) = cell_lengths ( 1 ) / REAL ( nout ( 1 ),KIND=dp)
  ENDIF
  IF ( nout ( 2 ) > grid_b % npts ( 2 ) ) THEN
    grid_b % npts ( 2 ) = nout ( 2 )
    dr ( 2 ) = cell_lengths ( 2 ) / REAL ( nout ( 2 ),KIND=dp)
  ENDIF
  IF ( nout ( 3 ) > grid_b % npts ( 3 ) ) THEN
    grid_b % npts ( 3 ) = nout ( 3 )
    dr ( 3 ) = cell_lengths ( 3 ) / REAL ( nout ( 3 ),KIND=dp)
  ENDIF

! big grid bounds
  grid_b % bounds ( 1, : ) = - grid_b % npts / 2
  grid_b % bounds ( 2, : ) = + ( grid_b % npts - 1 ) / 2
  grid_b % grid_span = HALFSPACE
! small grid bounds
  grid_s % bounds ( 1, : ) =   -nout ( : ) / 2
  grid_s % bounds ( 2, : ) = ( +nout ( : ) - 1 ) / 2
  grid_s % grid_span = HALFSPACE
  grid_s % npts = nout

  cutoff = pw_find_cutoff ( grid_b % npts, box_b%h_inv, error=error)

END SUBROUTINE  dg_find_cutoff

! *****************************************************************************
SUBROUTINE dg_get_spacing ( npts, cutoff_radius, dr )

    INTEGER, DIMENSION(:), INTENT(IN)        :: npts
    REAL(KIND=dp), INTENT(IN)                :: cutoff_radius
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: dr

  dr ( : ) = cutoff_radius  / ( REAL ( npts ( : ),KIND=dp) / 2.0_dp )

 END SUBROUTINE dg_get_spacing

! *****************************************************************************
SUBROUTINE dg_grid_change ( box_b, grid_b, grid_s, error)
    TYPE(cell_type), POINTER                 :: box_b
    TYPE(pw_grid_type), POINTER              :: grid_b, grid_s
    TYPE(cp_error_type), INTENT(inout)       :: error

    TYPE(cell_type), POINTER                 :: box_s, unit_box

  NULLIFY(box_s, unit_box)
  CALL cell_create(unit_box ,error=error)
  CALL cell_create(box_s,error=error)
  CALL dg_find_basis  ( grid_b % npts, box_b, unit_box )
  CALL dg_set_cell    ( grid_s % npts, unit_box, box_s )
  CALL pw_grid_change ( box_s, grid_s )
  CALL cell_release(unit_box,error=error)
  CALL cell_release(box_s,error=error)
END SUBROUTINE dg_grid_change

! *****************************************************************************
SUBROUTINE dg_find_radix ( dr, cell_lengths, npts )

    REAL(KIND=dp), INTENT(INOUT)             :: dr( 3 )
    REAL(KIND=dp), INTENT(IN)                :: cell_lengths( 3 )
    INTEGER, DIMENSION(:), INTENT(OUT)       :: npts

    INTEGER, DIMENSION(3)                    :: nin

  nin ( : ) = NINT ( cell_lengths ( : ) / dr ( : ) )
  CALL fft_radix_operations ( nin ( 1 ), npts ( 1 ), &
       operation = FFT_RADIX_CLOSEST )
  CALL fft_radix_operations ( nin ( 2 ), npts ( 2 ), &
       operation = FFT_RADIX_CLOSEST )
  CALL fft_radix_operations ( nin ( 3 ), npts ( 3 ), &
       operation = FFT_RADIX_CLOSEST )
  dr ( : ) = cell_lengths ( : ) / REAL ( npts ( : ),KIND=dp)

END SUBROUTINE dg_find_radix

! *****************************************************************************
SUBROUTINE dg_find_basis( npts, box, unit_box )
    INTEGER, DIMENSION(:), INTENT(IN)        :: npts
    TYPE(cell_type), POINTER                 :: box, unit_box

    INTEGER                                  :: i

  DO i = 1, 3
     unit_box % hmat ( :, i ) = box % hmat ( :, i ) / REAL ( npts ( : ),KIND=dp)
  END DO

  CALL init_cell ( unit_box )

END SUBROUTINE dg_find_basis

!! Calculation of the basis on the mesh 'box'

! *****************************************************************************
SUBROUTINE dg_set_cell (npts, unit_box, box )
    INTEGER, DIMENSION(:), INTENT(IN)        :: npts
    TYPE(cell_type), POINTER                 :: unit_box, box

! computing the unit vector along a, b, c and scaling it to length dr:

  box % hmat ( :, 1 ) = unit_box % hmat ( :, 1 ) * npts ( 1 )
  box % hmat ( :, 2 ) = unit_box % hmat ( :, 2 ) * npts ( 2 )
  box % hmat ( :, 3 ) = unit_box % hmat ( :, 3 ) * npts ( 3 )

  CALL init_cell ( box )

END SUBROUTINE dg_set_cell

! *****************************************************************************
SUBROUTINE dg_find_max_cutoff ( npts, grid_b )

    INTEGER, DIMENSION(:), INTENT(IN)        :: npts
    TYPE(pw_grid_type), DIMENSION(:), &
      INTENT(IN)                             :: grid_b

    CHARACTER(LEN=*), PARAMETER :: routineN = 'dg_find_max_cutoff', &
      routineP = moduleN//':'//routineN

  CALL stop_program(routineN,moduleN,__LINE__,"Not yet implemented.")

END SUBROUTINE dg_find_max_cutoff

! *****************************************************************************
SUBROUTINE dg_get_strucfac ( box, r, npts_s, npts_b, centre, lb, ex, ey, ez )

    TYPE(cell_type), INTENT(IN)              :: box
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: r
    INTEGER, DIMENSION(:), INTENT(IN)        :: npts_s, npts_b
    INTEGER, INTENT(OUT)                     :: centre( 3 )
    INTEGER, INTENT(IN)                      :: lb( 3 )
    COMPLEX(KIND=dp), DIMENSION(lb(1):), &
      INTENT(OUT)                            :: ex
    COMPLEX(KIND=dp), DIMENSION(lb(2):), &
      INTENT(OUT)                            :: ey
    COMPLEX(KIND=dp), DIMENSION(lb(3):), &
      INTENT(OUT)                            :: ez

    REAL(KIND=dp)                            :: delta( 3 )

  CALL dg_get_delta ( box, r, npts_s, npts_b, centre, delta)

  CALL structure_factor_evaluate ( delta, npts_s, lb, ex, ey, ez )

END SUBROUTINE dg_get_strucfac

! *****************************************************************************
SUBROUTINE dg_get_delta ( box, r, npts_s, npts_b, centre, delta )

    TYPE(cell_type), INTENT(IN)              :: box
    REAL(KIND=dp), DIMENSION(:), INTENT(IN)  :: r
    INTEGER, DIMENSION(:), INTENT(IN)        :: npts_s, npts_b
    INTEGER, DIMENSION(:), INTENT(OUT)       :: centre
    REAL(KIND=dp), DIMENSION(:), INTENT(OUT) :: delta

    REAL(KIND=dp), DIMENSION(3)              :: grid_i, s

! compute the scaled coordinate of atomi

  s = MATMUL ( box % h_inv, r )
  s = s - NINT ( s )

! find the continuous ``grid'' point (on big grid)
  grid_i ( 1:3 ) = REAL ( npts_b ( 1:3 ),KIND=dp) * s ( 1:3 )

! find the closest grid point (on big grid)
  centre ( : ) = NINT ( grid_i ( : ) )

! find the distance vector
  delta ( : ) = ( grid_i ( : ) - centre ( : ) ) / REAL ( npts_s ( : ),KIND=dp)

  centre ( : ) = centre ( : ) + npts_b(:)/2
  centre ( : ) = MODULO ( centre ( : ), npts_b(:) )
  centre ( : ) = centre ( : ) - npts_b(:)/2

END SUBROUTINE dg_get_delta

! *****************************************************************************
SUBROUTINE dg_sum_patch_coef ( rs, rhos, center )

    TYPE(realspace_grid_type), POINTER       :: rs
    TYPE(pw_p_type), INTENT(IN)              :: rhos
    INTEGER, DIMENSION(3), INTENT(IN)        :: center

    INTEGER                                  :: i, ia, ii
    INTEGER, DIMENSION(3)                    :: nc
    LOGICAL                                  :: folded

  folded = .FALSE.

  DO i =  rhos % pw % pw_grid % bounds( 1, 1 ), rhos % pw % pw_grid % bounds( 2, 1 )
    ia = i - rhos % pw % pw_grid % bounds( 1, 1 ) + 1
    ii = center ( 1 ) + i - rs % lb_local ( 1 )
    IF ( ii < 0 ) THEN
      rs % px ( ia ) = ii + rs % npts_local( 1 ) + 1
      folded = .TRUE.
    ELSEIF ( ii >=  rs % npts_local( 1 ) ) THEN
      rs % px ( ia ) = ii - rs % npts_local( 1 ) + 1
      folded = .TRUE.
    ELSE
      rs % px ( ia ) = ii + 1
    ENDIF
  END DO
  DO i =  rhos % pw % pw_grid % bounds( 1, 2 ), rhos % pw % pw_grid % bounds( 2, 2 )
    ia = i - rhos % pw % pw_grid % bounds( 1, 2 ) + 1
    ii = center ( 2 ) + i - rs % lb_local ( 2 )
    IF ( ii < 0 ) THEN
      rs % py ( ia ) = ii + rs % npts_local( 2 ) + 1
      folded = .TRUE.
    ELSEIF ( ii >=  rs % npts_local( 2 ) ) THEN
      rs % py ( ia ) = ii - rs % npts_local( 2 ) + 1
      folded = .TRUE.
    ELSE
      rs % py ( ia ) = ii + 1
    ENDIF
  END DO
  DO i = rhos % pw % pw_grid % bounds( 1, 3 ), rhos % pw % pw_grid % bounds( 2, 3 )
    ia = i - rhos % pw % pw_grid % bounds( 1, 3 ) + 1
    ii = center ( 3 ) + i - rs % lb_local ( 3 )
    IF ( ii < 0 ) THEN
      rs % pz ( ia ) = ii + rs % npts_local( 3 ) + 1
      folded = .TRUE.
    ELSEIF ( ii >=  rs % npts_local( 3 ) ) THEN
      rs % pz ( ia ) = ii - rs % npts_local( 3 ) + 1
      folded = .TRUE.
    ELSE
      rs % pz ( ia ) = ii + 1
    ENDIF
  END DO

  IF ( folded ) THEN
    CALL dg_add_patch ( rs % r, rhos % pw % cr3d, rhos % pw % pw_grid % npts,&
         rs % px, rs % py, rs % pz )
  ELSE
    nc ( 1 ) = rs % px ( 1 ) - 1
    nc ( 2 ) = rs % py ( 1 ) - 1
    nc ( 3 ) = rs % pz ( 1 ) - 1
    CALL dg_add_patch ( rs % r, rhos % pw % cr3d, rhos % pw % pw_grid % npts , nc )
  END IF

END SUBROUTINE dg_sum_patch_coef

! *****************************************************************************
SUBROUTINE dg_sum_patch_arr ( rs, rhos, center )

    TYPE(realspace_grid_type), POINTER       :: rs
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: rhos
    INTEGER, DIMENSION(3), INTENT(IN)        :: center

    INTEGER                                  :: i, ia, ii
    INTEGER, DIMENSION(3)                    :: lb, nc, ns, ub
    LOGICAL                                  :: folded

  ns ( 1 ) = SIZE ( rhos, 1 )
  ns ( 2 ) = SIZE ( rhos, 2 )
  ns ( 3 ) = SIZE ( rhos, 3 )
  lb = -(ns-1) / 2
  ub = lb + ns - 1
  folded = .FALSE.

  DO i = lb ( 1 ), ub ( 1 )
    ia = i - lb ( 1 ) + 1
    ii = center ( 1 ) + i - rs % lb_local ( 1 )
    IF ( ii < 0 ) THEN
      rs % px ( ia ) = ii + rs % npts_local( 1 ) + 1
      folded = .TRUE.
    ELSEIF ( ii >= rs % npts_local( 1 ) ) THEN
      rs % px ( ia ) = ii - rs % npts_local( 1 ) + 1
      folded = .TRUE.
    ELSE
      rs % px ( ia ) = ii + 1
    ENDIF
  END DO
  DO i = lb ( 2 ), ub ( 2 )
    ia = i - lb ( 2 ) + 1
    ii = center ( 2 ) + i - rs % lb_local ( 2 )
    IF ( ii < 0 ) THEN
      rs % py ( ia ) = ii + rs % npts_local( 2 ) + 1
      folded = .TRUE.
    ELSEIF ( ii >= rs % npts_local( 2 ) ) THEN
      rs % py ( ia ) = ii - rs % npts_local( 2 ) + 1
      folded = .TRUE.
    ELSE
      rs % py ( ia ) = ii + 1
    ENDIF
  END DO
  DO i = lb ( 3 ), ub ( 3 )
    ia = i - lb ( 3 ) + 1
    ii = center ( 3 ) + i - rs % lb_local ( 3 )
    IF ( ii < 0 ) THEN
      rs % pz ( ia ) = ii + rs % npts_local( 3 ) + 1
      folded = .TRUE.
    ELSEIF ( ii >= rs % npts_local( 3 ) ) THEN
      rs % pz ( ia ) = ii - rs % npts_local( 3 ) + 1
      folded = .TRUE.
    ELSE
      rs % pz ( ia ) = ii + 1
    ENDIF
  END DO

  IF ( folded ) THEN
    CALL dg_add_patch ( rs % r, rhos, ns, rs % px, rs % py, rs % pz )
  ELSE
    nc ( 1 ) = rs % px ( 1 ) - 1
    nc ( 2 ) = rs % py ( 1 ) - 1
    nc ( 3 ) = rs % pz ( 1 ) - 1
    CALL dg_add_patch ( rs % r, rhos, ns, nc )
  END IF

END SUBROUTINE dg_sum_patch_arr

! *****************************************************************************
SUBROUTINE dg_sum_patch_force_arr_3d ( drpot, rhos, center, force )

    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: drpot
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: rhos
    INTEGER, DIMENSION(3), INTENT(IN)        :: center
    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: force

    INTEGER                                  :: i, ia, ii
    INTEGER, DIMENSION(3)                    :: lb, nc, ns, ub
    LOGICAL                                  :: folded

  ns ( 1 ) = SIZE ( rhos, 1 )
  ns ( 2 ) = SIZE ( rhos, 2 )
  ns ( 3 ) = SIZE ( rhos, 3 )
  lb = -(ns-1) / 2
  ub = lb + ns - 1
  folded = .FALSE.

  DO i = lb ( 1 ), ub ( 1 )
    ia = i - lb ( 1 ) + 1
    ii = center ( 1 ) + i - drpot(1)%rs_grid%lb_local(1)
    IF ( ii < 0 ) THEN
      drpot(1)%rs_grid%px(ia) = ii + drpot(1)%rs_grid%npts_local(1) + 1
      folded = .TRUE.
    ELSEIF ( ii >= drpot(1)%rs_grid%npts_local(1) ) THEN
      drpot(1)%rs_grid%px(ia) = ii - drpot(1)%rs_grid%npts_local(1) + 1
      folded = .TRUE.
    ELSE
      drpot(1)%rs_grid%px(ia) = ii + 1
    ENDIF
  END DO
  DO i = lb ( 2 ), ub ( 2 )
    ia = i - lb ( 2 ) + 1
    ii = center ( 2 ) + i - drpot(1)%rs_grid%lb_local(2)
    IF ( ii < 0 ) THEN
      drpot(1)%rs_grid%py(ia) = ii + drpot(1)%rs_grid%npts_local(2) + 1
      folded = .TRUE.
    ELSEIF ( ii >=  drpot(1)%rs_grid%npts_local(2) ) THEN
      drpot(1)%rs_grid%py(ia) = ii - drpot(1)%rs_grid%npts_local(2) + 1
      folded = .TRUE.
    ELSE
      drpot(1)%rs_grid%py(ia) = ii + 1
    ENDIF
  END DO
  DO i = lb ( 3 ), ub ( 3 )
    ia = i - lb ( 3 ) + 1
    ii = center ( 3 ) + i - drpot(1)%rs_grid%lb_local(3)
    IF ( ii < 0 ) THEN
      drpot(1)%rs_grid%pz(ia) = ii + drpot(1)%rs_grid%npts_local(3) + 1
      folded = .TRUE.
    ELSEIF ( ii >= drpot(1)%rs_grid%npts_local(3) ) THEN
      drpot(1)%rs_grid%pz(ia) = ii - drpot(1)%rs_grid%npts_local(3) + 1
      folded = .TRUE.
    ELSE
      drpot(1)%rs_grid%pz(ia) = ii + 1
    ENDIF
  END DO

  IF ( folded ) THEN
    CALL dg_int_patch_3d ( drpot(1)%rs_grid%r, drpot(2)%rs_grid%r,&
         drpot(3)%rs_grid%r, rhos, force, ns, &
         drpot(1)%rs_grid%px, drpot(1)%rs_grid%py, drpot(1)%rs_grid%pz )
  ELSE
    nc ( 1 ) = drpot(1)%rs_grid%px(1) - 1
    nc ( 2 ) = drpot(1)%rs_grid%py(1) - 1
    nc ( 3 ) = drpot(1)%rs_grid%pz(1) - 1
    CALL dg_int_patch_3d ( drpot(1)%rs_grid%r, drpot(2)%rs_grid%r,&
         drpot(3)%rs_grid%r, rhos, force, ns, nc )
  END IF

END SUBROUTINE dg_sum_patch_force_arr_3d

! *****************************************************************************
SUBROUTINE dg_sum_patch_force_arr_1d ( drpot, rhos, center, force )

    TYPE(realspace_grid_type), POINTER       :: drpot
    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: rhos
    INTEGER, DIMENSION(3), INTENT(IN)        :: center
    REAL(KIND=dp), INTENT(OUT)               :: force

    INTEGER                                  :: i, ia, ii
    INTEGER, DIMENSION(3)                    :: lb, nc, ns, ub
    LOGICAL                                  :: folded

  ns ( 1 ) = SIZE ( rhos, 1 )
  ns ( 2 ) = SIZE ( rhos, 2 )
  ns ( 3 ) = SIZE ( rhos, 3 )
  lb = -(ns-1) / 2
  ub = lb + ns - 1
  folded = .FALSE.

  DO i = lb ( 1 ), ub ( 1 )
    ia = i - lb ( 1 ) + 1
    ii = center ( 1 ) + i - drpot % lb_local ( 1 )
    IF ( ii < 0 ) THEN
      drpot % px ( ia ) = ii +  drpot % npts_local( 1 ) + 1
      folded = .TRUE.
    ELSEIF ( ii >= drpot % desc % npts( 1 ) ) THEN
      drpot % px ( ia ) = ii - drpot % npts_local( 1 ) + 1
      folded = .TRUE.
    ELSE
      drpot % px ( ia ) = ii + 1
    ENDIF
  END DO
  DO i = lb ( 2 ), ub ( 2 )
    ia = i - lb ( 2 ) + 1
    ii = center ( 2 ) + i - drpot % lb_local ( 2 )
    IF ( ii < 0 ) THEN
      drpot % py ( ia ) = ii + drpot % npts_local( 2 ) + 1
      folded = .TRUE.
    ELSEIF ( ii >= drpot % desc % npts( 2 ) ) THEN
      drpot % py ( ia ) = ii - drpot % npts_local( 2 ) + 1
      folded = .TRUE.
    ELSE
      drpot % py ( ia ) = ii + 1
    ENDIF
  END DO
  DO i = lb ( 3 ), ub ( 3 )
    ia = i - lb ( 3 ) + 1
    ii = center ( 3 ) + i - drpot % lb_local ( 3 )
    IF ( ii < 0 ) THEN
      drpot % pz ( ia ) = ii + drpot % npts_local( 3 ) + 1
      folded = .TRUE.
    ELSEIF ( ii >= drpot % desc % npts( 3 ) ) THEN
      drpot % pz ( ia ) = ii - drpot % npts_local( 3 ) + 1
      folded = .TRUE.
    ELSE
      drpot % pz ( ia ) = ii + 1
    ENDIF
  END DO

  IF ( folded ) THEN
    CALL dg_int_patch_1d ( drpot % r, rhos, force, ns, &
         drpot % px, drpot % py, drpot % pz )
  ELSE
    nc ( 1 ) = drpot % px ( 1 ) - 1
    nc ( 2 ) = drpot % py ( 1 ) - 1
    nc ( 3 ) = drpot % pz ( 1 ) - 1
    CALL dg_int_patch_1d ( drpot % r, rhos, force, ns, nc )
  END IF

END SUBROUTINE dg_sum_patch_force_arr_1d

! *****************************************************************************
SUBROUTINE dg_sum_patch_force_coef_3d ( drpot, rhos, center, force )

    TYPE(realspace_grid_p_type), &
      DIMENSION(:), POINTER                  :: drpot
    TYPE(pw_p_type), INTENT(IN)              :: rhos
    INTEGER, DIMENSION(3), INTENT(IN)        :: center
    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: force

    INTEGER                                  :: i, ia, ii
    INTEGER, DIMENSION(3)                    :: nc
    LOGICAL                                  :: folded

  folded = .FALSE.

  DO i = rhos % pw % pw_grid % bounds( 1, 1 ), rhos % pw % pw_grid % bounds( 2, 1 )
    ia = i - rhos % pw % pw_grid % bounds( 1, 1 ) + 1
    ii = center ( 1 ) + i - drpot ( 1 )%rs_grid % lb_local ( 1 )
    IF ( ii < 0 ) THEN
      drpot ( 1 )%rs_grid % px ( ia ) = ii + drpot ( 1 )%rs_grid % desc % npts( 1 ) + 1
      folded = .TRUE.
    ELSEIF ( ii >= drpot ( 1 )%rs_grid % desc % npts( 1 ) ) THEN
      drpot ( 1 )%rs_grid % px ( ia ) = ii - drpot ( 1 )%rs_grid % desc % npts( 1 ) + 1
      folded = .TRUE.
    ELSE
      drpot ( 1 )%rs_grid % px ( ia ) = ii + 1
    ENDIF
  END DO
  DO i = rhos % pw % pw_grid % bounds( 1, 2 ), rhos % pw % pw_grid % bounds( 2, 2 )
    ia = i - rhos % pw % pw_grid % bounds( 1, 2 ) + 1
    ii = center ( 2 ) + i - drpot ( 1 )%rs_grid % lb_local ( 2 )
    IF ( ii < 0 ) THEN
      drpot ( 1 )%rs_grid % py ( ia ) = ii + drpot ( 1 )%rs_grid % desc % npts( 2 ) + 1
      folded = .TRUE.
    ELSEIF ( ii >= drpot ( 1 )%rs_grid % desc % npts( 2 ) ) THEN
      drpot ( 1 )%rs_grid % py ( ia ) = ii - drpot ( 1 )%rs_grid % desc % npts( 2 ) + 1
      folded = .TRUE.
    ELSE
      drpot ( 1 )%rs_grid % py ( ia ) = ii + 1
    ENDIF
  END DO
  DO i = rhos % pw % pw_grid % bounds( 1, 3 ), rhos % pw % pw_grid % bounds( 2, 3 )
    ia = i - rhos % pw % pw_grid % bounds( 1, 3 ) + 1
    ii = center ( 3 ) + i - drpot ( 1 )%rs_grid % lb_local ( 3 )
    IF ( ii < 0 ) THEN
      drpot ( 1 )%rs_grid % pz ( ia ) = ii + drpot ( 1 )%rs_grid % desc % npts( 3 ) + 1
      folded = .TRUE.
    ELSEIF ( ii >= drpot ( 1 )%rs_grid % desc % npts( 3 ) ) THEN
      drpot ( 1 )%rs_grid % pz ( ia ) = ii - drpot ( 1 )%rs_grid % desc % npts( 3 ) + 1
      folded = .TRUE.
    ELSE
      drpot ( 1 )%rs_grid % pz ( ia ) = ii + 1
    ENDIF
  END DO

  IF ( folded ) THEN
    CALL dg_int_patch_3d ( drpot ( 1 )%rs_grid % r, drpot ( 2 )%rs_grid % r,&
         drpot ( 3 )%rs_grid % r, rhos % pw % cr3d, force, rhos % pw % pw_grid % npts, &
         drpot ( 1 )%rs_grid % px, drpot ( 1 )%rs_grid % py, drpot ( 1 )%rs_grid % pz )
  ELSE
    nc ( 1 ) = drpot ( 1 )%rs_grid % px ( 1 ) - 1
    nc ( 2 ) = drpot ( 1 )%rs_grid % py ( 1 ) - 1
    nc ( 3 ) = drpot ( 1 )%rs_grid % pz ( 1 ) - 1
    CALL dg_int_patch_3d ( drpot ( 1 )%rs_grid % r, drpot ( 2 )%rs_grid % r,&
         drpot ( 3 )%rs_grid % r, rhos % pw % cr3d, force, rhos % pw % pw_grid % npts, nc )
  END IF

END SUBROUTINE dg_sum_patch_force_coef_3d

! *****************************************************************************
SUBROUTINE dg_sum_patch_force_coef_1d ( drpot, rhos, center, force )

    TYPE(realspace_grid_type), POINTER       :: drpot
    TYPE(pw_p_type), INTENT(IN)              :: rhos
    INTEGER, DIMENSION(3), INTENT(IN)        :: center
    REAL(KIND=dp), INTENT(OUT)               :: force

    INTEGER                                  :: i, ia, ii
    INTEGER, DIMENSION(3)                    :: nc
    LOGICAL                                  :: folded

  folded = .FALSE.

  DO i = rhos % pw % pw_grid % bounds( 1, 1 ), rhos % pw % pw_grid % bounds( 2, 1 )
    ia = i - rhos % pw % pw_grid % bounds( 1, 1 ) + 1
    ii = center ( 1 ) + i - drpot % lb_local ( 1 )
    IF ( ii < 0 ) THEN
      drpot % px ( ia ) = ii + drpot % desc % npts( 1 ) + 1
      folded = .TRUE.
    ELSEIF ( ii >= drpot % desc % npts( 1 ) ) THEN
      drpot % px ( ia ) = ii - drpot % desc % npts( 1 ) + 1
      folded = .TRUE.
    ELSE
      drpot % px ( ia ) = ii + 1
    ENDIF
  END DO
  DO i = rhos % pw % pw_grid % bounds( 1, 2 ), rhos % pw % pw_grid % bounds( 2, 2 )
    ia = i - rhos % pw % pw_grid % bounds( 1, 2 ) + 1
    ii = center ( 2 ) + i - drpot % lb_local ( 2 )
    IF ( ii < 0 ) THEN
      drpot % py ( ia ) = ii + drpot % desc % npts( 2 ) + 1
      folded = .TRUE.
    ELSEIF ( ii >= drpot % desc % npts( 2 ) ) THEN
      drpot % py ( ia ) = ii - drpot % desc % npts( 2 ) + 1
      folded = .TRUE.
    ELSE
      drpot % py ( ia ) = ii + 1
    ENDIF
  END DO
  DO i = rhos % pw % pw_grid % bounds( 1, 3 ), rhos % pw % pw_grid % bounds( 2, 3 )
    ia = i - rhos % pw % pw_grid % bounds( 1, 3 ) + 1
    ii = center ( 3 ) + i - drpot % lb_local ( 3 )
    IF ( ii < 0 ) THEN
      drpot % pz ( ia ) = ii + drpot % desc % npts( 3 ) + 1
      folded = .TRUE.
    ELSEIF ( ii >= drpot % desc % npts( 3 ) ) THEN
      drpot % pz ( ia ) = ii - drpot % desc % npts( 3 ) + 1
      folded = .TRUE.
    ELSE
      drpot % pz ( ia ) = ii + 1
    ENDIF
  END DO

  IF ( folded ) THEN
    CALL dg_int_patch_1d ( drpot % r, rhos % pw % cr3d, force,&
         rhos % pw % pw_grid % npts, drpot % px, drpot % py, drpot % pz )
  ELSE
    nc ( 1 ) = drpot % px ( 1 ) - 1
    nc ( 2 ) = drpot % py ( 1 ) - 1
    nc ( 3 ) = drpot % pz ( 1 ) - 1
    CALL dg_int_patch_1d ( drpot % r, rhos % pw % cr3d, force,rhos % pw % pw_grid % npts, nc )
  END IF

END SUBROUTINE dg_sum_patch_force_coef_1d

! *****************************************************************************
SUBROUTINE dg_get_patch_1 ( rho0, rhos1, charge1, ex1, ey1, ez1 )

    TYPE(pw_p_type), INTENT(IN)              :: rho0
    TYPE(pw_p_type), INTENT(INOUT)           :: rhos1
    REAL(KIND=dp), INTENT(IN)                :: charge1
    COMPLEX(KIND=dp), DIMENSION(:), &
      INTENT(IN)                             :: ex1, ey1, ez1

    CHARACTER(len=*), PARAMETER :: routineN = 'dg_get_patch_1', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp)                         :: za, zb
    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:)                           :: zs
    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: cd
    INTEGER                                  :: ierr, n, nd( 3 )

  nd = rhos1 % pw % pw_grid % npts
  n = PRODUCT ( nd )

    ALLOCATE (zs(nd(1)*nd(2)),STAT=ierr)
    IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                    "zs",2*dp_size*nd(1)*nd(2))
    zs=0.0_dp
    ALLOCATE (cd(nd(1),nd(2),nd(3)),STAT=ierr)
    IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                    "cd",2*dp_size*nd(1)*nd(2)*nd(3))
    cd=0.0_dp

  za = CMPLX ( 0.0_dp, 0.0_dp,KIND=dp)
  zb = CMPLX ( charge1, 0.0_dp,KIND=dp)
  CALL rankup ( nd, za, cd, zb, ex1, ey1, ez1, zs )
  IF ( rho0 % pw % in_use == REALDATA3D ) &
      CALL vr_x_vc ( n, rho0 % pw % cr3d, cd )
  IF ( rho0 % pw % in_use == COMPLEXDATA3D ) &
      CALL vc_x_vc ( n, rho0 % pw % cc3d, cd )
  CALL fft3d ( BWFFT, nd, cd )
  CALL copy_cr ( n, cd, rhos1 % pw % cr3d )

  DEALLOCATE (zs,STAT=ierr)
  IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"zs")
  DEALLOCATE (cd,STAT=ierr)
  IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"cd")

END SUBROUTINE dg_get_patch_1

! *****************************************************************************
SUBROUTINE dg_get_patch_2 ( rho0, rhos1, rhos2, charge1, charge2, &
           ex1, ey1, ez1, ex2, ey2, ez2 )

    TYPE(pw_p_type), INTENT(IN)              :: rho0
    TYPE(pw_p_type), INTENT(INOUT)           :: rhos1, rhos2
    REAL(KIND=dp), INTENT(IN)                :: charge1, charge2
    COMPLEX(KIND=dp), DIMENSION(:), &
      INTENT(IN)                             :: ex1, ey1, ez1, ex2, ey2, ez2

    CHARACTER(len=*), PARAMETER :: routineN = 'dg_get_patch_2', &
      routineP = moduleN//':'//routineN

    COMPLEX(KIND=dp)                         :: za, zb
    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:)                           :: zs
    COMPLEX(KIND=dp), ALLOCATABLE, &
      DIMENSION(:, :, :)                     :: cd
    INTEGER                                  :: ierr, n, nd( 3 )

  nd = rhos1 % pw % pw_grid % npts
  n = PRODUCT ( nd )

    ALLOCATE (zs(nd(1)*nd(2)),STAT=ierr)
    IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                    "zs",2*dp_size*nd(1)*nd(2))
    zs=0.0_dp
    ALLOCATE (cd(nd(1),nd(2),nd(3)),STAT=ierr)
    IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,&
                                    "cd",2*dp_size*nd(1)*nd(2)*nd(3))
    cd=0.0_dp

  za = CMPLX ( 0.0_dp, 0.0_dp,KIND=dp)
  zb = CMPLX ( charge2, 0.0_dp,KIND=dp)
  CALL rankup ( nd, za, cd, zb, ex2, ey2, ez2, zs )
  za = CMPLX ( 0.0_dp, 1.0_dp,KIND=dp)
  zb = CMPLX ( charge1, 0.0_dp,KIND=dp)
  CALL rankup ( nd, za, cd, zb, ex1, ey1, ez1, zs )
  IF ( rho0 % pw % in_use == REALDATA3D ) &
      CALL vr_x_vc ( n, rho0 % pw % cr3d, cd )
  IF ( rho0 % pw % in_use == COMPLEXDATA3D ) &
      CALL vc_x_vc ( n, rho0 % pw % cc3d, cd )
  CALL fft3d ( BWFFT, nd, cd )
  CALL copy_cri ( n, cd, rhos1 % pw % cr3d, rhos2 % pw % cr3d )

  DEALLOCATE (zs,stat=ierr)
  IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"zs")
  DEALLOCATE (cd,stat=ierr)
  IF (ierr /= 0) CALL stop_memory(routineN,moduleN,__LINE__,"cd")

END SUBROUTINE dg_get_patch_2

! *****************************************************************************
SUBROUTINE dg_add_patch_simple ( rb, rs, ns, nc )

    REAL(KIND=dp), INTENT(INOUT)             :: rb( :, :, : )
    REAL(KIND=dp), INTENT(IN)                :: rs( :, :, : )
    INTEGER, INTENT(IN)                      :: ns( 3 ), nc( 3 )

    INTEGER                                  :: i, ii, j, jj, k, kk

  DO k = 1, ns ( 3 )
    kk = nc ( 3 ) + k
    DO j = 1, ns ( 2 )
      jj = nc ( 2 ) + j
      DO i = 1, ns ( 1 )
        ii = nc ( 1 ) + i
        rb ( ii, jj, kk ) = rb ( ii, jj, kk ) + rs ( i, j, k )
      END DO
    END DO
  END DO

END SUBROUTINE dg_add_patch_simple

! *****************************************************************************
SUBROUTINE dg_add_patch_folded ( rb, rs, ns, px, py, pz )

    REAL(KIND=dp), INTENT(INOUT)             :: rb( :, :, : )
    REAL(KIND=dp), INTENT(IN)                :: rs( :, :, : )
    INTEGER, INTENT(IN)                      :: ns( : )
    INTEGER, DIMENSION(:), INTENT(IN)        :: px, py, pz

    INTEGER                                  :: i, ii, j, jj, k, kk

  DO k = 1, ns ( 3 )
    kk = pz ( k )
    DO j = 1, ns ( 2 )
      jj = py ( j )
      DO i = 1, ns ( 1 )
        ii =  px ( i )
        rb ( ii, jj, kk ) = rb ( ii, jj, kk ) + rs ( i, j, k )
      END DO
    END DO
  END DO

END SUBROUTINE dg_add_patch_folded

! *****************************************************************************
SUBROUTINE dg_int_patch_simple_3d ( rbx, rby, rbz, rs, f, ns, nc )

    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: rbx, rby, rbz, rs
    REAL(KIND=dp), DIMENSION(3), INTENT(OUT) :: f
    INTEGER, INTENT(IN)                      :: ns( 3 ), nc( 3 )

    INTEGER                                  :: i, ii, j, jj, k, kk
    REAL(KIND=dp)                            :: s

  f = 0.0_dp
  DO k = 1, ns ( 3 )
    kk = nc ( 3 ) + k
    DO j = 1, ns ( 2 )
      jj = nc ( 2 ) + j
      DO i = 1, ns ( 1 )
        ii = nc ( 1 ) + i
        s = rs ( i, j, k )
        f ( 1 ) = f ( 1 ) + s * rbx ( ii, jj, kk )
        f ( 2 ) = f ( 2 ) + s * rby ( ii, jj, kk )
        f ( 3 ) = f ( 3 ) + s * rbz ( ii, jj, kk )
      END DO
    END DO
  END DO

END SUBROUTINE dg_int_patch_simple_3d

! *****************************************************************************
SUBROUTINE dg_int_patch_simple_1d ( rb, rs, f, ns, nc )

    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: rb, rs
    REAL(KIND=dp), INTENT(OUT)               :: f
    INTEGER, INTENT(IN)                      :: ns( 3 ), nc( 3 )

    INTEGER                                  :: i, ii, j, jj, k, kk
    REAL(KIND=dp)                            :: s

  f = 0.0_dp
  DO k = 1, ns ( 3 )
    kk = nc ( 3 ) + k
    DO j = 1, ns ( 2 )
      jj = nc ( 2 ) + j
      DO i = 1, ns ( 1 )
        ii = nc ( 1 ) + i
        s = rs ( i, j, k )
        f = f  + s * rb ( ii, jj, kk )
      END DO
    END DO
  END DO

END SUBROUTINE dg_int_patch_simple_1d

! *****************************************************************************
SUBROUTINE dg_int_patch_folded_3d ( rbx, rby, rbz, rs, f, ns, px, py, pz )

    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: rbx, rby, rbz, rs
    REAL(KIND=dp), DIMENSION(3), &
      INTENT(INOUT)                          :: f
    INTEGER, INTENT(IN)                      :: ns( 3 )
    INTEGER, DIMENSION(:), INTENT(IN)        :: px, py, pz

    INTEGER                                  :: i, ii, j, jj, k, kk
    REAL(KIND=dp)                            :: s

  f = 0.0_dp
  DO k = 1, ns ( 3 )
    kk = pz ( k )
    DO j = 1, ns ( 2 )
      jj = py ( j )
      DO i = 1, ns ( 1 )
        ii = px ( i )
        s = rs ( i, j, k )
        f ( 1 ) = f ( 1 ) + s * rbx ( ii, jj, kk )
        f ( 2 ) = f ( 2 ) + s * rby ( ii, jj, kk )
        f ( 3 ) = f ( 3 ) + s * rbz ( ii, jj, kk )
      END DO
    END DO
  END DO

END SUBROUTINE dg_int_patch_folded_3d

! *****************************************************************************
SUBROUTINE dg_int_patch_folded_1d ( rb, rs, f, ns, px, py, pz )

    REAL(KIND=dp), DIMENSION(:, :, :), &
      INTENT(IN)                             :: rb, rs
    REAL(KIND=dp), INTENT(INOUT)             :: f
    INTEGER, INTENT(IN)                      :: ns( 3 )
    INTEGER, DIMENSION(:), INTENT(IN)        :: px, py, pz

    INTEGER                                  :: i, ii, j, jj, k, kk
    REAL(KIND=dp)                            :: s

  f = 0.0_dp
  DO k = 1, ns ( 3 )
    kk = pz ( k )
    DO j = 1, ns ( 2 )
      jj = py ( j )
      DO i = 1, ns ( 1 )
        ii = px ( i )
        s = rs ( i, j, k )
        f = f + s * rb ( ii, jj, kk )
      END DO
    END DO
  END DO

END SUBROUTINE dg_int_patch_folded_1d

END MODULE dgs
