!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief wrappers for the actual blacs calls.
!>      all functionality needed in the code should actually be provide by cp_blacs_env
!>      these functions should be private members of that module
!> \note
!>      http://www.netlib.org/blacs/BLACS/QRef.html 
!> \par History
!>      12.2003 created [Joost]
!> \author Joost VandeVondele
! *****************************************************************************
MODULE cp_blacs_calls
  USE f77_blas
  USE kinds,                           ONLY: dp
#include "cp_common_uses.h"

  IMPLICIT NONE
  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'cp_blacs_calls'
  PRIVATE
  ! setup / info calls
  PUBLIC :: cp_blacs_gridinit, cp_blacs_get, cp_blacs_set, cp_blacs_pnum, &
            cp_blacs_pinfo, cp_blacs_gridexit, cp_blacs_gridinfo
  ! actual message passing
  PUBLIC :: cp_blacs_zgebs2d, cp_blacs_zgebr2d, cp_blacs_dgebs2d, &
            cp_blacs_dgebr2d, cp_blacs_igsum2d

  INTERFACE cp_blacs_zgebs2d
     MODULE PROCEDURE cp_blacs_zgebs2d_array, cp_blacs_zgebs2d_scalar
  END INTERFACE 

  INTERFACE cp_blacs_zgebr2d
     MODULE PROCEDURE cp_blacs_zgebr2d_array, cp_blacs_zgebr2d_scalar
  END INTERFACE 

  INTERFACE cp_blacs_dgebs2d
     MODULE PROCEDURE cp_blacs_dgebs2d_array, cp_blacs_dgebs2d_scalar
  END INTERFACE 

  INTERFACE cp_blacs_dgebr2d
     MODULE PROCEDURE cp_blacs_dgebr2d_array, cp_blacs_dgebr2d_scalar
  END INTERFACE 

  INTERFACE cp_blacs_igsum2d
     MODULE PROCEDURE cp_blacs_igsum2d_array, cp_blacs_igsum2d_scalar
  END INTERFACE 
!***
CONTAINS

! *****************************************************************************
SUBROUTINE cp_blacs_gridinit(context,order,nprow,npcol)
   INTEGER, INTENT(INOUT) :: context
   CHARACTER(len=1), INTENT(IN):: order
   INTEGER, INTENT(IN)    :: nprow, npcol
#if defined(__BLACS)
   CALL blacs_gridinit(context,order,nprow,npcol)
#endif
END SUBROUTINE cp_blacs_gridinit

! *****************************************************************************
SUBROUTINE cp_blacs_gridexit(context)
   INTEGER, INTENT(IN) :: context
#if defined(__BLACS)
   CALL blacs_gridexit(context)
#endif
END SUBROUTINE cp_blacs_gridexit

! *****************************************************************************
SUBROUTINE cp_blacs_gridinfo(context,nprow,npcol,myprow,mypcol)
   INTEGER, INTENT(IN)  :: context
   INTEGER, INTENT(OUT) :: nprow,npcol,myprow,mypcol
#if defined(__BLACS)
   CALL blacs_gridinfo(context,nprow,npcol,myprow,mypcol)
#endif
END SUBROUTINE cp_blacs_gridinfo

! *****************************************************************************
FUNCTION cp_blacs_pnum(context,myprow,mypcol) RESULT(pnum)
   INTEGER, INTENT(IN) :: context,myprow,mypcol
   INTEGER :: pnum
#if defined(__BLACS)
   INTEGER, EXTERNAL :: blacs_pnum
   pnum = blacs_pnum(context,myprow,mypcol)
#else
   pnum = 0
#endif
END FUNCTION cp_blacs_pnum

! WHAT = 0 : Handle indicating default system context;  ! DO NOT USE (i.e. use para_env%group)
! WHAT = 1 : The BLACS message ID range; 
! WHAT = 2 : The BLACS debug level the library was compiled with; 
! WHAT = 10: Handle indicating the system context used to define the BLACS context whose handle is ICONTXT; 
! WHAT = 11: Number of rings multiring topology is presently using; 
! WHAT = 12: Number of branches general tree topology is presently using. 
! *****************************************************************************
SUBROUTINE cp_blacs_get(context,what,val)
   INTEGER, INTENT(IN)  :: context,what
   INTEGER, INTENT(OUT) :: val
#if defined(__BLACS)
   CALL blacs_get(context,what,val)
#endif
END SUBROUTINE cp_blacs_get

! WHAT = 0 : Handle indicating default system context;  ! DO NOT USE (i.e. use para_env%group)
! WHAT = 1 : The BLACS message ID range; 
! WHAT = 2 : The BLACS debug level the library was compiled with; 
! WHAT = 10: Handle indicating the system context used to define the BLACS context whose handle is ICONTXT; 
! WHAT = 11: Number of rings multiring topology is presently using; 
! WHAT = 12: Number of branches general tree topology is presently using. 
! WHAT = 15: If non-zero, makes topology choice for repeatable collectives
! *****************************************************************************
SUBROUTINE cp_blacs_set(context,what,val)
   INTEGER, INTENT(IN)  :: context,what,val
#if defined(__BLACS)
   CALL blacs_set(context,what,val)
#endif
END SUBROUTINE cp_blacs_set

! the total number available to the system, not for a particular context (DO NOT USE)
! *****************************************************************************
SUBROUTINE cp_blacs_pinfo(mype,npe)
   INTEGER, INTENT(OUT) :: mype,npe
#if defined(__BLACS)
   CALL blacs_pinfo(mype,npe)
#endif
END SUBROUTINE cp_blacs_pinfo

! the strange scalar / array interface seems to be needed

! *****************************************************************************
SUBROUTINE cp_blacs_zgebs2d_array(ICONTXT, SCOPE, TOP,             M, N, A, LDA )
  INTEGER, INTENT(IN)     :: ICONTXT
  CHARACTER(len=1), INTENT(IN) :: SCOPE,TOP
  INTEGER, INTENT(IN)     :: M,N,LDA
  COMPLEX(KIND=dp), DIMENSION(LDA,*) :: A
#if defined(__BLACS)
  CALL zgebs2d(ICONTXT, SCOPE, TOP,             M, N, A, LDA )
#endif
END SUBROUTINE
! *****************************************************************************
SUBROUTINE cp_blacs_zgebr2d_array(ICONTXT, SCOPE, TOP,             M, N, A, LDA, RSRC, CSRC )
  INTEGER, INTENT(IN)     :: ICONTXT
  CHARACTER(len=1), INTENT(IN) :: SCOPE,TOP
  INTEGER, INTENT(IN)     :: M,N,LDA
  INTEGER, INTENT(IN)     :: RSRC,CSRC
  COMPLEX(KIND=dp), DIMENSION(LDA,*) :: A
#if defined(__BLACS)
  CALL zgebr2d(ICONTXT, SCOPE, TOP,             M, N, A, LDA, RSRC, CSRC )
#endif
END SUBROUTINE
! *****************************************************************************
SUBROUTINE cp_blacs_zgebs2d_scalar(ICONTXT, SCOPE, TOP,             M, N, A, LDA )
  INTEGER, INTENT(IN)     :: ICONTXT
  CHARACTER(len=1), INTENT(IN) :: SCOPE,TOP
  INTEGER, INTENT(IN)     :: M,N,LDA
  COMPLEX(KIND=dp)            :: A
#if defined(__BLACS)
  CALL zgebs2d(ICONTXT, SCOPE, TOP,             M, N, A, LDA )
#endif
END SUBROUTINE
! *****************************************************************************
SUBROUTINE cp_blacs_zgebr2d_scalar(ICONTXT, SCOPE, TOP,             M, N, A, LDA, RSRC, CSRC )
  INTEGER, INTENT(IN)     :: ICONTXT
  CHARACTER(len=1), INTENT(IN) :: SCOPE,TOP
  INTEGER, INTENT(IN)     :: M,N,LDA
  INTEGER, INTENT(IN)     :: RSRC,CSRC
  COMPLEX(KIND=dp)            :: A
#if defined(__BLACS)
  CALL zgebr2d(ICONTXT, SCOPE, TOP,             M, N, A, LDA, RSRC, CSRC )
#endif
END SUBROUTINE

! *****************************************************************************
SUBROUTINE cp_blacs_dgebs2d_array(ICONTXT, SCOPE, TOP,             M, N, A, LDA )
  INTEGER, INTENT(IN)     :: ICONTXT
  CHARACTER(len=1), INTENT(IN) :: SCOPE,TOP
  INTEGER, INTENT(IN)     :: M,N,LDA
  REAL(KIND=dp), DIMENSION(LDA,*) :: A
#if defined(__BLACS)
  CALL dgebs2d(ICONTXT, SCOPE, TOP,             M, N, A, LDA )
#endif
END SUBROUTINE
! *****************************************************************************
SUBROUTINE cp_blacs_dgebr2d_array(ICONTXT, SCOPE, TOP,             M, N, A, LDA, RSRC, CSRC )
  INTEGER, INTENT(IN)     :: ICONTXT
  CHARACTER(len=1), INTENT(IN) :: SCOPE,TOP
  INTEGER, INTENT(IN)     :: M,N,LDA
  INTEGER, INTENT(IN)     :: RSRC,CSRC
  REAL(KIND=dp), DIMENSION(LDA,*) :: A
#if defined(__BLACS)
  CALL dgebr2d(ICONTXT, SCOPE, TOP,             M, N, A, LDA, RSRC, CSRC )
#endif
END SUBROUTINE
! *****************************************************************************
SUBROUTINE cp_blacs_dgebs2d_scalar(ICONTXT, SCOPE, TOP,             M, N, A, LDA )
  INTEGER, INTENT(IN)     :: ICONTXT
  CHARACTER(len=1), INTENT(IN) :: SCOPE,TOP
  INTEGER, INTENT(IN)     :: M,N,LDA
  REAL(KIND=dp)               :: A
#if defined(__BLACS)
  CALL dgebs2d(ICONTXT, SCOPE, TOP,             M, N, A, LDA )
#endif
END SUBROUTINE
! *****************************************************************************
SUBROUTINE cp_blacs_dgebr2d_scalar(ICONTXT, SCOPE, TOP,             M, N, A, LDA, RSRC, CSRC )
  INTEGER, INTENT(IN)     :: ICONTXT
  CHARACTER(len=1), INTENT(IN) :: SCOPE,TOP
  INTEGER, INTENT(IN)     :: M,N,LDA
  INTEGER, INTENT(IN)     :: RSRC,CSRC
  REAL(KIND=dp)               :: A
#if defined(__BLACS)
  CALL dgebr2d(ICONTXT, SCOPE, TOP,             M, N, A, LDA, RSRC, CSRC )
#endif
END SUBROUTINE

! *****************************************************************************
SUBROUTINE cp_blacs_igsum2d_array(ICONTXT, SCOPE, TOP,             M, N, A, LDA, RDEST, CDEST )
  INTEGER, INTENT(IN)     :: ICONTXT
  CHARACTER(len=1), INTENT(IN) :: SCOPE,TOP
  INTEGER, INTENT(IN)     :: M,N,LDA
  INTEGER, INTENT(IN)     :: RDEST,CDEST
  INTEGER, DIMENSION(LDA,*) :: A
#if defined(__BLACS)
  CALL igsum2d(ICONTXT, SCOPE, TOP,             M, N, A, LDA, RDEST, CDEST )
#endif
END SUBROUTINE
! *****************************************************************************
SUBROUTINE cp_blacs_igsum2d_scalar(ICONTXT, SCOPE, TOP,             M, N, A, LDA, RDEST, CDEST )
  INTEGER, INTENT(IN)     :: ICONTXT
  CHARACTER(len=1), INTENT(IN) :: SCOPE,TOP
  INTEGER, INTENT(IN)     :: M,N,LDA
  INTEGER, INTENT(IN)     :: RDEST,CDEST
  INTEGER                 :: A
#if defined(__BLACS)
  CALL igsum2d(ICONTXT, SCOPE, TOP,             M, N, A, LDA, RDEST, CDEST )
#endif
END SUBROUTINE

END MODULE cp_blacs_calls
