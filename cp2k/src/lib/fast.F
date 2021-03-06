!-----------------------------------------------------------------------------!
! Some F77 style routines to get full performance
!-----------------------------------------------------------------------------!
! *****************************************************************************
SUBROUTINE rankup ( n, za, cmat, zb, ex, ey, ez, scr )
!
! cmat(i,j,k) <- za * cmat(i,j,k) + ex(i) * ey(j) * ez(k)
!
  USE lib_kinds,                        ONLY: wp

  INTEGER, DIMENSION ( 3 ), INTENT ( IN ) :: n
  COMPLEX ( wp ), INTENT ( IN ) :: za, zb
  COMPLEX ( wp ), INTENT ( IN ) :: ex ( * )
  COMPLEX ( wp ), INTENT ( IN ) :: ey ( * )
  COMPLEX ( wp ), INTENT ( IN ) :: ez ( * )
  COMPLEX ( wp ), INTENT ( INOUT ) :: cmat ( * )
  COMPLEX ( wp ), INTENT ( INOUT ) :: scr ( * )

  INTEGER :: n2, n3
  COMPLEX ( wp ) :: cone
  
  cone = CMPLX ( 1._wp, KIND = wp )
  n2 = n ( 1 ) * n ( 2 )
  n3 = n2 * n ( 3 )
  scr ( 1:n2 ) = CMPLX ( 0._wp, KIND = wp )
#if defined (__SGL) 
  CALL cgeru ( n ( 1 ), n ( 2 ), zb, ex, 1, ey, 1, scr, n ( 1 ) )
  CALL cscal ( n3, za, cmat, 1 )
  CALL cgeru ( n2, n ( 3 ), cone, scr, 1, ez, 1, cmat, n2 )
#else
  CALL zgeru ( n ( 1 ), n ( 2 ), zb, ex, 1, ey, 1, scr, n ( 1 ) )
  CALL zscal ( n3, za, cmat, 1 )
  CALL zgeru ( n2, n ( 3 ), cone, scr, 1, ez, 1, cmat, n2 )
#endif

END SUBROUTINE rankup 

!-----------------------------------------------------------------------------!

! *****************************************************************************
SUBROUTINE vc_x_vc ( n, cvec2, cvec )
!
! cden(i) <- cden(i) * cden2(i)
!
  USE lib_kinds,                        ONLY: wp

  INTEGER, INTENT ( IN ) :: n
  REAL ( wp ), INTENT ( IN ) :: cvec2 (2,*)
  REAL ( wp ), INTENT ( INOUT ) :: cvec (2,*)

  INTEGER :: i
  REAL ( wp ) :: re_c, im_c, re_c2, im_c2

!$OMP PARALLEL DO
  DO i = 1, n
    re_c = cvec ( 1, i )
    im_c = cvec ( 2, i )
    re_c2 = cvec2 ( 1, i )
    im_c2 = cvec2 ( 2, i )
    cvec ( 1, i ) = re_c * re_c2 - im_c * im_c2
    cvec ( 2, i ) = im_c * re_c2 + re_c * im_c2
  END DO

END SUBROUTINE vc_x_vc 

!-----------------------------------------------------------------------------!

! *****************************************************************************
SUBROUTINE vr_x_vc ( n, rvec, cvec )
!
! cden(i) <- cden(i) * rden(i)
!
  USE lib_kinds,                        ONLY: wp

  INTEGER, INTENT ( IN ) :: n
  REAL ( wp ), INTENT ( IN ) :: rvec (*)
  REAL ( wp ), INTENT ( INOUT ) :: cvec (2,*)

  INTEGER :: i

!$OMP PARALLEL DO
  DO i = 1, n
    cvec ( 1, i ) = cvec ( 1, i ) * rvec ( i )
    cvec ( 2, i ) = cvec ( 2, i ) * rvec ( i )
  END DO

END SUBROUTINE vr_x_vc 

!-----------------------------------------------------------------------------!

! *****************************************************************************
SUBROUTINE copy_cri ( n, z, r1, r2 )
!
! r1 = real ( z )
! r2 = imag ( z )
!
  USE lib_kinds,                        ONLY: wp

  INTEGER, INTENT ( IN ) :: n
  REAL ( wp ), INTENT ( OUT ) :: r1(*), r2(*)
  REAL ( wp ), INTENT ( IN ) :: z(2,*)

  INTEGER :: i

!$OMP PARALLEL DO
  DO i = 1, n
    r1 ( i ) = z ( 1, i )
    r2 ( i ) = z ( 2, i )
  END DO

END SUBROUTINE copy_cri

!-----------------------------------------------------------------------------!
! *****************************************************************************
SUBROUTINE copy_cc ( n, in, out )
!
  USE lib_kinds,                        ONLY: wp

  INTEGER, INTENT ( IN ) :: n
  REAL ( wp ), INTENT ( IN ) :: in(n)
  REAL ( wp ), INTENT ( OUT ) :: out(n)
  INTEGER :: i
!$OMP PARALLEL DO
  DO i=1,n
       out(i)=in(i)
  ENDDO

END SUBROUTINE copy_cc

!-----------------------------------------------------------------------------!

! *****************************************************************************
SUBROUTINE copy_cr ( n, z, r1 )
!
! r1 = real ( z )
!
  USE lib_kinds,                        ONLY: wp

  INTEGER, INTENT ( IN ) :: n
  REAL ( wp ), INTENT ( OUT ) :: r1(n)
  REAL ( wp ), INTENT ( IN ) :: z(2,n)
  INTEGER :: i

!$OMP PARALLEL DO
  DO i=1,n
       r1(i)=z(1,i)
  ENDDO

END SUBROUTINE copy_cr

!-----------------------------------------------------------------------------!

! *****************************************************************************
SUBROUTINE copy_rc ( n, r1, z )
!
! z = r1
!
  USE lib_kinds,                        ONLY: wp

  INTEGER, INTENT ( IN ) :: n
  REAL ( wp ), INTENT ( IN ) :: r1(n)
  REAL ( wp ), INTENT ( OUT ) :: z(2,n)
  INTEGER :: i
!$OMP PARALLEL DO
  DO i=1,n
       z(1,i)=r1(i)
       z(2,i)=0._wp
  ENDDO

END SUBROUTINE copy_rc

!-----------------------------------------------------------------------------!

! *****************************************************************************
SUBROUTINE copy_cz ( n, c, z )
!
! z = c
!
  USE lib_kinds,                        ONLY: dp, sp

  INTEGER, INTENT ( IN ) :: n
  COMPLEX ( sp ), INTENT ( IN ) :: c(n)
  COMPLEX ( dp ), INTENT ( OUT ) :: z(n)
  INTEGER :: i
!$OMP PARALLEL DO
  DO i=1,n
    z(i)=c(i)
  ENDDO

END SUBROUTINE copy_cz
! *****************************************************************************
SUBROUTINE copy_zc ( n, z, c )
!
! c = z
!
  USE lib_kinds,                        ONLY: dp, sp

  INTEGER, INTENT ( IN ) :: n
  COMPLEX ( dp ), INTENT ( IN ) :: z(n)
  COMPLEX ( sp ), INTENT ( OUT ) :: c(n)
  INTEGER :: i
!$OMP PARALLEL DO
  DO i=1,n
    c(i)=z(i)
  ENDDO

END SUBROUTINE copy_zc
! *****************************************************************************
SUBROUTINE copy_zz ( n, c, z )
!
! z = c
!
  USE lib_kinds,                        ONLY: wp

  INTEGER, INTENT ( IN ) :: n
  COMPLEX ( wp ), INTENT ( IN ) :: c(n)
  COMPLEX ( wp ), INTENT ( OUT ) :: z(n)
  INTEGER :: i
!$OMP PARALLEL DO
  DO i=1,n
    z(i)=c(i)
  ENDDO

END SUBROUTINE copy_zz

!-----------------------------------------------------------------------------!

