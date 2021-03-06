!-----------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations         !
!   Copyright (C) 2000 - 2011  CP2K developers group                          !
!-----------------------------------------------------------------------------!

! *****************************************************************************
!> \brief   Tools usually found in a standard library.
!> \author  Urban Borstnik
!> \date    2011-11-02
!> \version 1.0
!>
!> <b>Modification history:</b>
!> - Created 2011-11-02
! *****************************************************************************
MODULE dbcsr_toollib

  USE dbcsr_error_handling
  USE dbcsr_kinds,                     ONLY: dp,&
                                             int_4,&
                                             int_8
  USE dbcsr_machine,                   ONLY: m_walltime

  !$ USE OMP_LIB

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'dbcsr_toollib'

  REAL, PARAMETER                      :: default_resize_factor = 1.618034


  PUBLIC :: dbcsr_sleep
  ! Block sizes and arrays
  PUBLIC :: dbcsr_unpack_i8_2i4, make_coordinate_tuple
  PUBLIC :: swap
  ! string routines
  PUBLIC :: uppercase
  ! math routines
  PUBLIC :: gcd, lcm, ceil_log2
  ! utility routines
  PUBLIC :: sort, joaat_hash
  PUBLIC :: partitions

  INTERFACE swap
     MODULE PROCEDURE iswap, bswap
  END INTERFACE

  INTERFACE sort
     MODULE PROCEDURE sort_i4, sort_i8
  END INTERFACE


CONTAINS

  ELEMENTAL FUNCTION make_coordinate_tuple (most, least) RESULT (tuple)
    INTEGER, INTENT(IN)                      :: most, least
    INTEGER(KIND=int_8)                      :: tuple

!tuple = IOR (ISHFT (most, 32), least)

    tuple = most
    tuple = IOR (ISHFT(tuple, 32), INT(least, int_8))
  END FUNCTION make_coordinate_tuple


! *****************************************************************************
!> \brief Swaps two integers
!> \param[in,out] a,b         Integers to swap
! *****************************************************************************
  ELEMENTAL SUBROUTINE iswap(a, b)
    INTEGER, INTENT(INOUT)                   :: a, b

    INTEGER                                  :: tmp

    tmp = a
    a = b
    b = tmp
  END SUBROUTINE iswap


! *****************************************************************************
!> \brief Swaps two logicals
!> \param[in,out] a,b         Logicals to swap
! *****************************************************************************
  ELEMENTAL SUBROUTINE bswap(a, b)
    LOGICAL, INTENT(INOUT)                   :: a, b

    LOGICAL                                  :: tmp

    tmp = a
    a = b
    b = tmp
  END SUBROUTINE bswap


! *****************************************************************************
!> \brief Busy-loop sleep.
!> \param[in] duration        time to sleep/s.
! *****************************************************************************
  SUBROUTINE dbcsr_sleep(duration)
    REAL(KIND=dp), INTENT(IN)                :: duration

    REAL(KIND=dp)                            :: t_now, t_start

    t_start = m_walltime ()
    t_now = t_start
    DO WHILE (t_now .LT. t_start + duration)
       t_now = m_walltime ()
    ENDDO
  END SUBROUTINE dbcsr_sleep


! *****************************************************************************
!> \brief Splits an array of int8 values into two int4 arrays.
!> \param[in] merged          array of merged values
!> \param[out] array_upper    array to fill with the upper bytes of the merged
!>                            values
!> \param[out] array_lower    array to fill with the lower bytes of the merged
!>                            values
! *****************************************************************************
  SUBROUTINE dbcsr_unpack_i8_2i4 (merged, array_upper, array_lower)
    INTEGER(KIND=int_8), DIMENSION(:), &
      INTENT(IN)                             :: merged
    INTEGER(KIND=int_4), DIMENSION(:), &
      INTENT(OUT)                            :: array_upper, array_lower

    INTEGER(KIND=int_8), PARAMETER           :: lmask8 = 4294967295_int_8

    INTEGER                                  :: i

!
!   ---------------------------------------------------------------------------
! Lmask is used to filter in the lower 4 bytes and so its lower 32 bits are
! set to 1: lmask8 = 2^32-1.
! Umask is used to filter in the higher 4 bytes and so its higher 32 bits
! are set to 1: umask8 = 2^32-1 << 32
!lmask8 = 4294967295 ! 2^32-1
!umask8 = 18446744069414584320 ! (2^32-1) * 2^32 = (2^64-1)-(2^32-1)

    DO i = 1, SIZE (merged)
       array_upper(i) = INT (ISHFT (merged(i), -32), KIND=int_4)
       array_lower(i) = INT (IAND (merged(i), lmask8), KIND=int_4)
    ENDDO
  END SUBROUTINE dbcsr_unpack_i8_2i4


  SUBROUTINE uppercase(string)
    CHARACTER(LEN=*), INTENT(INOUT)          :: string

    INTEGER                                  :: i, iascii

    DO i=1,LEN_TRIM(string)
       iascii = ICHAR(string(i:i))
       IF ((iascii >= 97).AND.(iascii <= 122)) THEN
          string(i:i) = CHAR(iascii - 32)
       END IF
    END DO

  END SUBROUTINE uppercase


  ELEMENTAL FUNCTION gcd(a,b)
    INTEGER, INTENT(IN)                      :: a, b
    INTEGER                                  :: gcd

    INTEGER                                  :: aa, ab, l, rem, s

    aa=ABS(a)
    ab=ABS(b)
    IF (aa<ab) THEN
       s=aa
       l=ab
    ELSE
       s=ab
       l=aa
    ENDIF
    IF (s.NE.0) THEN
       DO
          rem=MOD(l,s)
          IF (rem==0) EXIT
          l=s
          s=rem
       ENDDO
       GCD=s
    ELSE
       GCD=l
    ENDIF
  END FUNCTION gcd

  ELEMENTAL FUNCTION lcm(a,b)
    INTEGER, INTENT(IN)                      :: a, b
    INTEGER                                  :: lcm

    INTEGER                                  :: tmp

    tmp = gcd(a,b)
    IF (tmp==0) THEN
       lcm = 0
    ELSE
       ! could still overflow if the true lcm is larger than maxint
       lcm = ABS((a/tmp)*b)
    END IF
  END FUNCTION lcm


! *****************************************************************************
!> \brief Calculates the ceiling of a log base 2 for integers
!> \note Should be careful for highest integer values.
!> \param[in] lv    Value for which to calculate the function
!> \result cl2      Result of function on lv
! *****************************************************************************
  ELEMENTAL FUNCTION ceil_log2 (lv) RESULT (cl2)
    INTEGER, INTENT(IN)                      :: lv
    INTEGER                                  :: cl2

    INTEGER                                  :: bit

    bit = 0
    DO WHILE (ISHFT(1, bit) .LT. lv .AND. bit .LT. BIT_SIZE(lv))
       bit = bit + 1
    ENDDO
    bit = bit-1
    ! Set to zero if lv is too big.
    IF (bit .EQ. BIT_SIZE(lv)-2) bit = 0
    cl2 = bit+1
  END FUNCTION ceil_log2


  SUBROUTINE sort_i4 ( arr, n, index )
    INTEGER, INTENT(IN)                      :: n
    INTEGER(KIND=int_4), INTENT(INOUT)       :: arr(1:n)
    INTEGER, INTENT(OUT)                     :: INDEX(1:n)

    CHARACTER(LEN=*), PARAMETER :: routineN = 'sort_i4', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: m = 7, nstack = 50

    INTEGER                                  :: i, ib, ir, istack(1:nstack), &
                                                itemp, j, jstack, k, l
    INTEGER(kind=int_4)                      :: a, temp

    IF (n==0) RETURN
    DO i = 1, n
       INDEX(i) = i
    END DO
    !
    IF (ALL(arr==arr(1))) RETURN ! Nothing to order..
    jstack = 0
    l = 1
    ir = n
1   IF (ir-l<m) THEN
       DO j = l + 1, ir
          a = arr(j)
          ib = INDEX(j)
          DO i = j - 1, 1, -1
             IF (arr(i)<=a) GO TO 2
             arr(i+1) = arr(i)
             INDEX(i+1) = INDEX(i)
          END DO
          i = 0
2         arr(i+1) = a
          INDEX(i+1) = ib
       END DO
       IF (jstack==0) RETURN
       ir = istack(jstack)
       l = istack(jstack-1)
       jstack = jstack - 2
    ELSE
       k = (l+ir)/2
       temp = arr(k)
       arr(k) = arr(l+1)
       arr(l+1) = temp
       itemp = INDEX(k)
       INDEX(k) = INDEX(l+1)
       INDEX(l+1) = itemp
       IF (arr(l+1)>arr(ir)) THEN
          temp = arr(l+1)
          arr(l+1) = arr(ir)
          arr(ir) = temp
          itemp = INDEX(l+1)
          INDEX(l+1) = INDEX(ir)
          INDEX(ir) = itemp
       END IF
       IF (arr(l)>arr(ir)) THEN
          temp = arr(l)
          arr(l) = arr(ir)
          arr(ir) = temp
          itemp = INDEX(l)
          INDEX(l) = INDEX(ir)
          INDEX(ir) = itemp
       END IF
       IF (arr(l+1)>arr(l)) THEN
          temp = arr(l+1)
          arr(l+1) = arr(l)
          arr(l) = temp
          itemp = INDEX(l+1)
          INDEX(l+1) = INDEX(l)
          INDEX(l) = itemp
       END IF
       i = l + 1
       j = ir
       a = arr(l)
       ib = INDEX(l)
3      CONTINUE
       i = i + 1
       IF (arr(i)<a) GO TO 3
4      CONTINUE
       j = j - 1
       IF (arr(j)>a) GO TO 4
       IF (j<i) GO TO 5
       temp = arr(i)
       arr(i) = arr(j)
       arr(j) = temp
       itemp = INDEX(i)
       INDEX(i) = INDEX(j)
       INDEX(j) = itemp
       GO TO 3
5      arr(l) = arr(j)
       arr(j) = a
       INDEX(l) = INDEX(j)
       INDEX(j) = ib
       jstack = jstack + 2
       IF (jstack>nstack) STOP ' Nstack too small in sortr'
       IF (ir-i+1>=j-l) THEN
          istack(jstack) = ir
          istack(jstack-1) = i
          ir = j - 1
       ELSE
          istack(jstack) = j - 1
          istack(jstack-1) = l
          l = i
       END IF
    END IF
    GO TO 1
  END SUBROUTINE sort_i4


  SUBROUTINE sort_i8 ( arr, n, index )
    INTEGER, INTENT(IN)                      :: n
    INTEGER(KIND=int_8), INTENT(INOUT)       :: arr(1:n)
    INTEGER, INTENT(OUT)                     :: INDEX(1:n)

    CHARACTER(LEN=*), PARAMETER :: routineN = 'sort_i8', &
      routineP = moduleN//':'//routineN
    INTEGER, PARAMETER                       :: m = 7, nstack = 50

    INTEGER                                  :: i, ib, ir, istack(1:nstack), &
                                                itemp, j, jstack, k, l
    INTEGER(kind=int_8)                      :: a, temp

    IF (n==0) RETURN
    DO i = 1, n
       INDEX(i) = i
    END DO
    !
    IF (ALL(arr==arr(1))) RETURN ! Nothing to order..
    jstack = 0
    l = 1
    ir = n
1   IF (ir-l<m) THEN
       DO j = l + 1, ir
          a = arr(j)
          ib = INDEX(j)
          DO i = j - 1, 1, -1
             IF (arr(i)<=a) GO TO 2
             arr(i+1) = arr(i)
             INDEX(i+1) = INDEX(i)
          END DO
          i = 0
2         arr(i+1) = a
          INDEX(i+1) = ib
       END DO
       IF (jstack==0) RETURN
       ir = istack(jstack)
       l = istack(jstack-1)
       jstack = jstack - 2
    ELSE
       k = (l+ir)/2
       temp = arr(k)
       arr(k) = arr(l+1)
       arr(l+1) = temp
       itemp = INDEX(k)
       INDEX(k) = INDEX(l+1)
       INDEX(l+1) = itemp
       IF (arr(l+1)>arr(ir)) THEN
          temp = arr(l+1)
          arr(l+1) = arr(ir)
          arr(ir) = temp
          itemp = INDEX(l+1)
          INDEX(l+1) = INDEX(ir)
          INDEX(ir) = itemp
       END IF
       IF (arr(l)>arr(ir)) THEN
          temp = arr(l)
          arr(l) = arr(ir)
          arr(ir) = temp
          itemp = INDEX(l)
          INDEX(l) = INDEX(ir)
          INDEX(ir) = itemp
       END IF
       IF (arr(l+1)>arr(l)) THEN
          temp = arr(l+1)
          arr(l+1) = arr(l)
          arr(l) = temp
          itemp = INDEX(l+1)
          INDEX(l+1) = INDEX(l)
          INDEX(l) = itemp
       END IF
       i = l + 1
       j = ir
       a = arr(l)
       ib = INDEX(l)
3      CONTINUE
       i = i + 1
       IF (arr(i)<a) GO TO 3
4      CONTINUE
       j = j - 1
       IF (arr(j)>a) GO TO 4
       IF (j<i) GO TO 5
       temp = arr(i)
       arr(i) = arr(j)
       arr(j) = temp
       itemp = INDEX(i)
       INDEX(i) = INDEX(j)
       INDEX(j) = itemp
       GO TO 3
5      arr(l) = arr(j)
       arr(j) = a
       INDEX(l) = INDEX(j)
       INDEX(j) = ib
       jstack = jstack + 2
       IF (jstack>nstack) STOP ' Nstack too small in sortr'
       IF (ir-i+1>=j-l) THEN
          istack(jstack) = ir
          istack(jstack-1) = i
          ir = j - 1
       ELSE
          istack(jstack) = j - 1
          istack(jstack-1) = l
          l = i
       END IF
    END IF
    GO TO 1
  END SUBROUTINE sort_i8


! *****************************************************************************
!> \brief Calculates partition
!> \param[in] total_size      Length to partition
!> \param[in] me              My position (from 0 to num-1)
!> \param[in] num             Number of postitions
!> \param[out] first          First element of my position
!> \param[out] last           Last element of my position
! *****************************************************************************
  pure SUBROUTINE partitions (total_size, me, num, first, last)
    INTEGER, INTENT(IN)                      :: total_size, me, num
    INTEGER, INTENT(OUT)                     :: first, last

    first = INT(&
           (  INT(total_size,kind=int_8) &
            * INT(me,kind=int_8) &
           ) &
         / INT(num,kind=int_8), kind=int_4) + 1
    last = INT(&
           (  INT(total_size,kind=int_8) &
            * INT(me+1,kind=int_8) &
           ) &
         / INT(num,kind=int_8), kind=int_4)
  END SUBROUTINE partitions


! *****************************************************************************
!> \brief generates the hash of a string and the index in the table
!> \param key a string of any length
!> \note
!>       http://en.wikipedia.org/wiki/Hash_table
!>       http://www.burtleburtle.net/bob/hash/doobs.html
!>       However, since fortran doesn't have an unsigned 4 byte int 
!>       we compute it using an integer with the appropriate range 
!>       we return already the index in the table as a final result
!> \par History
!>       09.2006 created [Joost VandeVondele]
!>       2010-11-23 copied from CP2K and converted to take integers
! *****************************************************************************
  FUNCTION joaat_hash(key) RESULT(hash_index)
    INTEGER, DIMENSION(:), INTENT(IN)        :: key
    INTEGER                                  :: hash_index

    INTEGER(KIND=int_8), PARAMETER           :: b32 = 2_int_8**32-1_int_8

    INTEGER                                  :: i, j
    INTEGER(KIND=int_8)                      :: byte, hash

    hash=0_int_8
    DO i=1,SIZE(key)
       DO j = 0, 3
          byte = IAND (ISHFT (key(i), -j*8), 255)
          hash=IAND(hash+byte                           ,b32)
          hash=IAND(     hash+IAND(ISHFT(hash,10),b32)  ,b32)
          hash=IAND(IEOR(hash,IAND(ISHFT(hash,-6),b32)) ,b32)
       ENDDO
    ENDDO
    hash=IAND(     hash+IAND(ISHFT(hash,  3),b32)  ,b32)
    hash=IAND(IEOR(hash,IAND(ISHFT(hash,-11),b32)) ,b32)
    hash=IAND(     hash+IAND(ISHFT(hash, 15),b32)  ,b32)
    hash_index = hash
  END FUNCTION joaat_hash

END MODULE dbcsr_toollib
