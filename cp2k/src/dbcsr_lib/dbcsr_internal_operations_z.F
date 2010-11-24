
! *****************************************************************************
!> \brief Issues actual DGEMM calls.
!>
!> \param[in] params           Stack of DGEMM parameters
!> \param[in] n                Number of parameters
!> \param[in] left_data_a      Left-matrix data
!> \param[in] right_data_a     Right-matrix data
!> \param[in,out] product_data_area  Data for results
!> \param[out] lflop           (optional) Number of FLOPs used by DGEMM
! *****************************************************************************
  SUBROUTINE process_dgemm_stack_z(params,&
       params_a, params_b, params_c,&
       stack_size, stack_size_a, stack_size_b, stack_size_c,&
       left_data_a, right_data_a, product_data_a, use_plasma, lflop, ltime, error)
    INTEGER, INTENT(IN)                       :: stack_size, stack_size_a,&
                                                 stack_size_b, stack_size_c
    TYPE(dgemm_join), DIMENSION(1:stack_size), &
      INTENT(IN)                              :: params
    TYPE(block_parameters), DIMENSION(1:stack_size_a), INTENT(IN) :: params_a
    TYPE(block_parameters), DIMENSION(1:stack_size_b), INTENT(IN) :: params_b
    TYPE(block_parameters), DIMENSION(1:stack_size_c), INTENT(IN) :: params_c
    COMPLEX(kind=real_8), DIMENSION(:), INTENT(IN)         :: left_data_a, &
                                                 right_data_a
    COMPLEX(kind=real_8), DIMENSION(:), INTENT(INOUT)      :: product_data_a
    LOGICAL, INTENT(IN), OPTIONAL             :: use_plasma
    INTEGER(KIND=int_8), INTENT(OUT), &
      OPTIONAL                                :: lflop
    REAL(KIND=real_8), INTENT(OUT), OPTIONAL  :: ltime
    TYPE(dbcsr_error_type), INTENT(inout)     :: error

    CHARACTER(len=*), PARAMETER :: routineN = 'process_dgemm_stack_z', &
      routineP = moduleN//':'//routineN
    REAL, PARAMETER                           :: resize_factor = 1.618034

    INTEGER                                   :: c, plasma_info, lda, ldb, ldc, r, sp,&
                                                 poff
    INTEGER                                   :: sp_c, sp_a, sp_b
    LOGICAL                                   :: my_use_plasma
    LOGICAL, PARAMETER                        :: dbgv = .FALSE.

#ifdef __PLASMA
    INCLUDE 'plasmaf.h'
#endif
    DEBUG_HEADER
!   ---------------------------------------------------------------------------
    DEBUG_BODY

    my_use_plasma = .FALSE.
    IF (PRESENT(use_plasma)) my_use_plasma = use_plasma
    IF (PRESENT (lflop)) lflop = INT(0, int_8)
    IF (PRESENT (ltime)) ltime = 0.0_real_8
    ! Setup encapsulated data area
    DO sp = 1, stack_size
       sp_c = params(sp)%p_c
       sp_a = params(sp)%p_a
       sp_b = params(sp)%p_b
       IF (params_c(sp_c)%nze .EQ. 0 .OR. params(sp)%last_k .LT. 0) THEN
          CYCLE
       ENDIF
       !
       poff = params_c(sp_c)%offset
       !
       ! Avoid nagging compilers
       IF (params(sp)%beta%c_dp .EQ. CMPLX(0.0, 0.0, real_8)) THEN
          IF (.NOT. params_c(sp_c)%tr &
               .OR. params(sp)%last_n .EQ. params_c(sp_c)%logical_cols) THEN
             product_data_a(poff : poff-1+params_c(sp_c)%logical_rows*params(sp)%last_n) = CMPLX(0.0, 0.0, real_8)
             !CALL dbcsr_data_clear (product_data_area,&
             !     lb=poff,&
             !     ub=poff-1+params_c(sp_c)%logical_rows*params(sp)%last_n)
          ELSE
             ! It must be transposed.
             FORALL (r = 1 : params_c(sp_c)%logical_rows,&
                  c = 1 : params(sp)%last_n)
                product_data_a(poff-1+(r-1)*params_c(sp_c)%logical_cols+c) =&
                     CMPLX(0.0, 0.0, real_8)
             END FORALL
             !### optimize out the inner loop (turn into a range)
             !DO r = 1, params_c(sp_c)%logical_rows
             !   DO c = 1, params(sp)%last_n
             !      CALL dbcsr_data_clear(product_data_area,&
             !           lb=poff-1+(r-1)*params_c(sp_c)%logical_cols+c,&
             !           ub=poff-1+(r-1)*params_c(sp_c)%logical_cols+c)
             !   ENDDO
             !ENDDO
             !FORALL (r = 1:params(sp)%lrows_c, c = 1:params(sp)%last_n)
             !   product_data_d((r-1)*params(sp)%lcols_c+c) = 0.0_dp
             !END FORALL
          ENDIF
       ENDIF
       IF (params_a(sp_a)%nze .EQ. 0&
            .OR. params_b(sp_b)%nze .EQ. 0 &
            .OR. params(sp)%last_k .EQ. 0) THEN
          CYCLE
       ENDIF
       ! Now prepare and issue the multiplies.
       lda = select_n_or_t (params_a(sp_a)%tr,&
            params_a(sp_a)%logical_rows, params_a(sp_a)%logical_cols)
       ldb = select_n_or_t (params_b(sp_b)%tr,&
            params_b(sp_b)%logical_rows, params_b(sp_b)%logical_cols)
       ldc = select_n_or_t (params_c(sp_c)%tr,&
            params_c(sp_c)%logical_rows, params_c(sp_c)%logical_cols)
       !write(*,*)"LDx", lda, ldb, ldc
       IF (dbgv) THEN
          WRITE(*,*)"Multiplying left data", sp, params_a(sp_a)%offset,&
               params_a(sp_a)%nze
          WRITE(*,'(10(1X,F9.3))')left_data_a(params_a(sp_a)%offset&
               :params_a(sp_a)%offset+params_a(sp_a)%nze-1)
          WRITE(*,*)"Multiplying right data", sp, &
               params_b(sp_b)%offset, params_b(sp_b)%nze
          WRITE(*,'(10(1X,F9.3))')right_data_a(params_b(sp_b)%offset&
               :params_b(sp_b)%offset+params_b(sp_b)%nze-1)
          WRITE(*,*)"Multiplying product data", sp, poff, params_c(sp_c)%nze,"/",&
               SIZE(product_data_a), poff+params_c(sp_c)%nze-1
          WRITE(*,'(10(1X,F9.3))')product_data_a(poff:&
               poff+params_c(sp_c)%nze-1)
          WRITE(*,*)"again"
          WRITE(*,'(10(1X,F9.3))')product_data_a
       ENDIF
       IF (.NOT. params_c(sp_c)%tr) THEN
          IF (PRESENT (ltime)) ltime = ltime + m_walltime()
          IF(my_use_plasma) THEN
#ifdef __PLASMA
             CALL plasma_ZGEMM(&
                  plasma_mat_type(params_a(sp_a)%tr),&
                  plasma_mat_type(params_b(sp_b)%tr),&
                  params_c(sp_c)%logical_rows, params(sp)%last_n,& !m, n
                  params(sp)%last_k,& ! k
                  params(sp)%alpha%c_dp,&
                  left_data_a(params_a(sp_a)%offset), lda,&
                  right_data_a(params_b(sp_b)%offset), ldb,&
                  params(sp)%beta%c_dp,&
                  product_data_a(poff), ldc,&
                  plasma_info)
             CALL dbcsr_assert( plasma_info, "EQ", 0, dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                  "plasma_gemm failed", __LINE__, error=error)
#else
             CALL dbcsr_assert( .FALSE., dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                  "plasma badly set", __LINE__, error=error)
#endif
          ELSEIF (internal_gemm) THEN
             CALL internal_ZGEMM(&
                  params_a(sp_a)%tr,&
                  params_b(sp_b)%tr,&
                  params_c(sp_c)%logical_rows, params(sp)%last_n,& !m, n
                  params(sp)%last_k,& ! k
                  params(sp)%alpha%c_dp,&
                  left_data_a(params_a(sp_a)%offset:), lda,&
                  right_data_a(params_b(sp_b)%offset:), ldb,&
                  params(sp)%beta%c_dp,&
                  product_data_a(poff:), ldc)
          ELSE
             CALL ZGEMM(&
                  blas_mat_type(params_a(sp_a)%tr),&
                  blas_mat_type(params_b(sp_b)%tr),&
                  params_c(sp_c)%logical_rows, params(sp)%last_n,& !m, n
                  params(sp)%last_k,& ! k
                  params(sp)%alpha%c_dp,&
                  left_data_a(params_a(sp_a)%offset), lda,&
                  right_data_a(params_b(sp_b)%offset), ldb,&
                  params(sp)%beta%c_dp,&
                  product_data_a(poff), ldc)
          ENDIF
          IF (PRESENT (ltime)) ltime = m_walltime() - ltime
          IF (PRESENT (lflop)) &
               lflop = lflop + INT(2, int_8) * &
               INT(params_c(sp_c)%logical_rows, int_8)* &
               INT(params(sp)%last_n, int_8) *&
               INT(params(sp)%last_k, int_8)
       ELSE
          IF (PRESENT (ltime)) ltime = ltime + m_walltime()
          IF(my_use_plasma) THEN
#ifdef __PLASMA
             CALL plasma_ZGEMM(&
                  plasma_mat_type(.NOT.params_b(sp_b)%tr),&
                  plasma_mat_type(.NOT.params_a(sp_a)%tr),&
                  params_c(sp_c)%logical_cols, params_c(sp_c)%logical_rows,& !m, n (switched)
                  params(sp)%last_k,& ! k
                  params(sp)%alpha%c_dp,&
                  right_data_a(params_b(sp_b)%offset), ldb,&
                  left_data_a(params_a(sp_a)%offset), lda,&
                  params(sp)%beta%c_dp,&
                  product_data_a(poff), ldc,&
                  plasma_info)
             CALL dbcsr_assert( plasma_info, "EQ", 0, dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                  "plasma_gemm failed", __LINE__, error=error)
#else
             CALL dbcsr_assert( .FALSE., dbcsr_fatal_level, dbcsr_internal_error, routineN,&
                  "plasma badly set", __LINE__, error=error)
#endif
          ELSEIF (internal_gemm) THEN
             CALL internal_ZGEMM(&
                  .NOT.params_b(sp_b)%tr,&
                  .NOT.params_a(sp_a)%tr,&
                  params_c(sp_c)%logical_cols, params_c(sp_c)%logical_rows,& !m, n (switched)
                  params(sp)%last_k,& ! k
                  params(sp)%alpha%c_dp,&
                  right_data_a(params_b(sp_b)%offset:), ldb,&
                  left_data_a(params_a(sp_a)%offset:), lda,&
                  params(sp)%beta%c_dp,&
                  product_data_a(poff:), ldc)
          ELSE
             CALL ZGEMM(&
                  blas_mat_type(.NOT.params_b(sp_b)%tr),&
                  blas_mat_type(.NOT.params_a(sp_a)%tr),&
                  params_c(sp_c)%logical_cols, params_c(sp_c)%logical_rows,& !m, n (switched)
                  params(sp)%last_k,& ! k
                  params(sp)%alpha%c_dp,&
                  right_data_a(params_b(sp_b)%offset), ldb,&
                  left_data_a(params_a(sp_a)%offset), lda,&
                  params(sp)%beta%c_dp,&
                  product_data_a(poff), ldc)
          ENDIF
          IF (PRESENT (ltime)) ltime = m_walltime() - ltime
          IF (PRESENT (lflop)) &
               lflop = lflop + INT(2, int_8) * &
               INT(params(sp)%last_n, int_8) * &
               INT(params_c(sp_c)%logical_cols, int_8) * &
               INT(params(sp)%last_k, int_8)
       ENDIF
       IF (dbgv) THEN
          WRITE(*,*)"Resulting product data", poff,":",poff+params_c(sp_c)%nze-1,&
               SIZE(product_data_a)
          WRITE(*,'(10(1X,F9.3))')product_data_a(poff:&
               poff+params_c(sp_c)%nze-1)
          WRITE(*,*)"again"
          WRITE(*,'(10(1X,F9.3))')product_data_a
       ENDIF
    ENDDO
  END SUBROUTINE process_dgemm_stack_z


  PURE SUBROUTINE internal_ZGEMM(&
       TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
    LOGICAL, INTENT(IN)                      :: TRANSA, TRANSB
    INTEGER, INTENT(IN)                      :: M, N, K
    INTEGER, INTENT(IN)                      :: LDC
    COMPLEX(kind=real_8), INTENT(INOUT)                   :: C(LDC,*)
    COMPLEX(kind=real_8), INTENT(IN)                      :: BETA
    INTEGER, INTENT(IN)                      :: LDB
    COMPLEX(kind=real_8), INTENT(IN)                      :: B(LDB,*)
    INTEGER, INTENT(IN)                      :: LDA
    COMPLEX(kind=real_8), INTENT(IN)                      :: A(LDA,*), ALPHA
    IF (.NOT. transa .AND. .NOT. transb) THEN
       CALL internal_ZGEMM_nn(M,N,K,ALPHA,A,B,BETA,C)
    ELSEIF (.NOT. transa .AND. transb) THEN
       CALL internal_ZGEMM_nt(M,N,K,ALPHA,A,B,BETA,C)
    ELSEIF (transa .AND. .NOT. transb) THEN
       CALL internal_ZGEMM_tn(M,N,K,ALPHA,A,B,BETA,C)
    ELSEIF (transa .AND. transb) THEN
       CALL internal_ZGEMM_tt(M,N,K,ALPHA,A,B,BETA,C)
    ENDIF
  END SUBROUTINE internal_ZGEMM

  PURE SUBROUTINE internal_ZGEMM_nn(&
       M,N,K,ALPHA,A,B,BETA,C)
    INTEGER, INTENT(IN)                      :: M, N, K
    COMPLEX(kind=real_8), INTENT(INOUT)                   :: C(M,N)
    COMPLEX(kind=real_8), INTENT(IN)                      :: BETA
    COMPLEX(kind=real_8), INTENT(IN)                      :: B(K,N)
    COMPLEX(kind=real_8), INTENT(IN)                      :: A(M,K), ALPHA
    C(:,:) = BETA * C(:,:) &
         + ALPHA * MATMUL (A, B)
  END SUBROUTINE internal_ZGEMM_nn

  PURE SUBROUTINE internal_ZGEMM_nt(&
       M,N,K,ALPHA,A,B,BETA,C)
    INTEGER, INTENT(IN)                      :: M, N, K
    COMPLEX(kind=real_8), INTENT(INOUT)                   :: C(M,N)
    COMPLEX(kind=real_8), INTENT(IN)                      :: BETA
    COMPLEX(kind=real_8), INTENT(IN)                      :: B(N,K)
    COMPLEX(kind=real_8), INTENT(IN)                      :: A(M,K), ALPHA
    C(:,:) = BETA * C(:,:) &
         + ALPHA * MATMUL (A, TRANSPOSE(B))
  END SUBROUTINE internal_ZGEMM_nt

  PURE SUBROUTINE internal_ZGEMM_tn(&
       M,N,K,ALPHA,A,B,BETA,C)
    INTEGER, INTENT(IN)                      :: M, N, K
    COMPLEX(kind=real_8), INTENT(INOUT)                   :: C(M,N)
    COMPLEX(kind=real_8), INTENT(IN)                      :: BETA
    COMPLEX(kind=real_8), INTENT(IN)                      :: B(K,N)
    COMPLEX(kind=real_8), INTENT(IN)                      :: A(K,M), ALPHA
    C(:,:) = BETA * C(:,:) &
         + ALPHA * MATMUL (TRANSPOSE(A), B)
  END SUBROUTINE internal_ZGEMM_tn

  PURE SUBROUTINE internal_ZGEMM_tt(&
       M,N,K,ALPHA,A,B,BETA,C)
    INTEGER, INTENT(IN)                      :: M, N, K
    COMPLEX(kind=real_8), INTENT(INOUT)                   :: C(M,N)
    COMPLEX(kind=real_8), INTENT(IN)                      :: BETA
    COMPLEX(kind=real_8), INTENT(IN)                      :: B(N,K)
    COMPLEX(kind=real_8), INTENT(IN)                      :: A(K,M), ALPHA
    C(:,:) = BETA * C(:,:) &
         + ALPHA * MATMUL (TRANSPOSE(A), TRANSPOSE(B))
  END SUBROUTINE internal_ZGEMM_tt