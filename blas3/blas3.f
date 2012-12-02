      SUBROUTINE CGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
c     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      COMPLEX            ALPHA, BETA
c     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), C( LDC, * )
c     ..
c
c  Purpose
c  =======
c
c  CGEMM  performs one of the matrix-matrix operations
c
c     C := alpha*op( A )*op( B ) + beta*C,
c
c  where  op( X ) is one of
c
c     op( X ) = X   or   op( X ) = X'   or   op( X ) = conjg( X' ),
c
c  alpha and beta are scalars, and A, B and C are matrices, with op( A )
c  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
c
c  Parameters
c  ==========
c
c  TRANSA - CHARACTER*1.
c           On entry, TRANSA specifies the form of op( A ) to be used in
c           the matrix multiplication as follows:
c
c              TRANSA = 'N' or 'n',  op( A ) = A.
c
c              TRANSA = 'T' or 't',  op( A ) = A'.
c
c              TRANSA = 'C' or 'c',  op( A ) = conjg( A' ).
c
c           Unchanged on exit.
c
c  TRANSB - CHARACTER*1.
c           On entry, TRANSB specifies the form of op( B ) to be used in
c           the matrix multiplication as follows:
c
c              TRANSB = 'N' or 'n',  op( B ) = B.
c
c              TRANSB = 'T' or 't',  op( B ) = B'.
c
c              TRANSB = 'C' or 'c',  op( B ) = conjg( B' ).
c
c           Unchanged on exit.
c
c  M      - INTEGER.
c           On entry,  M  specifies  the number  of rows  of the  matrix
c           op( A )  and of the  matrix  C.  M  must  be at least  zero.
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry,  N  specifies the number  of columns of the matrix
c           op( B ) and the number of columns of the matrix C. N must be
c           at least zero.
c           Unchanged on exit.
c
c  K      - INTEGER.
c           On entry,  K  specifies  the number of columns of the matrix
c           op( A ) and the number of rows of the matrix op( B ). K must
c           be at least  zero.
c           Unchanged on exit.
c
c  ALPHA  - COMPLEX         .
c           On entry, ALPHA specifies the scalar alpha.
c           Unchanged on exit.
c
c  A      - COMPLEX          array of DIMENSION ( LDA, ka ), where ka is
c           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
c           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
c           part of the array  A  must contain the matrix  A,  otherwise
c           the leading  k by m  part of the array  A  must contain  the
c           matrix A.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
c           LDA must be at least  max( 1, m ), otherwise  LDA must be at
c           least  max( 1, k ).
c           Unchanged on exit.
c
c  B      - COMPLEX          array of DIMENSION ( LDB, kb ), where kb is
c           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
c           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
c           part of the array  B  must contain the matrix  B,  otherwise
c           the leading  n by k  part of the array  B  must contain  the
c           matrix B.
c           Unchanged on exit.
c
c  LDB    - INTEGER.
c           On entry, LDB specifies the first dimension of B as declared
c           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
c           LDB must be at least  max( 1, k ), otherwise  LDB must be at
c           least  max( 1, n ).
c           Unchanged on exit.
c
c  BETA   - COMPLEX         .
c           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
c           supplied as zero then C need not be set on input.
c           Unchanged on exit.
c
c  C      - COMPLEX          array of DIMENSION ( LDC, n ).
c           Before entry, the leading  m by n  part of the array  C must
c           contain the matrix  C,  except when  beta  is zero, in which
c           case C need not be set on entry.
c           On exit, the array  C  is overwritten by the  m by n  matrix
c           ( alpha*op( A )*op( B ) + beta*C ).
c
c  LDC    - INTEGER.
c           On entry, LDC specifies the first dimension of C as declared
c           in  the  calling  (sub)  program.   LDC  must  be  at  least
c           max( 1, m ).
c           Unchanged on exit.
c
c
c  Level 3 Blas routine.
c
c  -- Written on 8-February-1989.
c     Jack Dongarra, Argonne National Laboratory.
c     Iain Duff, AERE Harwell.
c     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c     Sven Hammarling, Numerical Algorithms Group Ltd.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX
c     .. Local Scalars ..
      LOGICAL            CONJA, CONJB, NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      COMPLEX            TEMP
c     .. Parameters ..
      COMPLEX            ONE
      PARAMETER        ( ONE  = ( 1.0E+0, 0.0E+0 ) )
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
c     ..
c     .. Executable Statements ..
c
c     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
c     conjugated or transposed, set  CONJA and CONJB  as true if  A  and
c     B  respectively are to be  transposed but  not conjugated  and set
c     NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A
c     and the number of rows of  B  respectively.
c
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      CONJA = LSAME( TRANSA, 'C' )
      CONJB = LSAME( TRANSB, 'C' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
c
c     Test the input parameters.
c
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.
     $         ( .NOT.CONJA                ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.
     $         ( .NOT.CONJB                ).AND.
     $         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CGEMM ', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
c
c     And when  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
c
c     Start the operations.
c
      IF( NOTB )THEN
         IF( NOTA )THEN
c
c           Form  C := alpha*A*B + beta*C.
c
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE IF( CONJA )THEN
c
c           Form  C := alpha*conjg( A' )*B + beta*C.
c
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + CONJG( A( L, I ) )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         ELSE
c
c           Form  C := alpha*A'*B + beta*C
c
            DO 150, J = 1, N
               DO 140, I = 1, M
                  TEMP = ZERO
                  DO 130, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  130             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  140          CONTINUE
  150       CONTINUE
         END IF
      ELSE IF( NOTA )THEN
         IF( CONJB )THEN
c
c           Form  C := alpha*A*conjg( B' ) + beta*C.
c
            DO 200, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 160, I = 1, M
                     C( I, J ) = ZERO
  160             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 170, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  170             CONTINUE
               END IF
               DO 190, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*CONJG( B( J, L ) )
                     DO 180, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  180                CONTINUE
                  END IF
  190          CONTINUE
  200       CONTINUE
         ELSE
c
c           Form  C := alpha*A*B'          + beta*C
c
            DO 250, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 210, I = 1, M
                     C( I, J ) = ZERO
  210             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 220, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  220             CONTINUE
               END IF
               DO 240, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 230, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  230                CONTINUE
                  END IF
  240          CONTINUE
  250       CONTINUE
         END IF
      ELSE IF( CONJA )THEN
         IF( CONJB )THEN
c
c           Form  C := alpha*conjg( A' )*conjg( B' ) + beta*C.
c
            DO 280, J = 1, N
               DO 270, I = 1, M
                  TEMP = ZERO
                  DO 260, L = 1, K
                     TEMP = TEMP + CONJG( A( L, I ) )*CONJG( B( J, L ) )
  260             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  270          CONTINUE
  280       CONTINUE
         ELSE
c
c           Form  C := alpha*conjg( A' )*B' + beta*C
c
            DO 310, J = 1, N
               DO 300, I = 1, M
                  TEMP = ZERO
                  DO 290, L = 1, K
                     TEMP = TEMP + CONJG( A( L, I ) )*B( J, L )
  290             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  300          CONTINUE
  310       CONTINUE
         END IF
      ELSE
         IF( CONJB )THEN
c
c           Form  C := alpha*A'*conjg( B' ) + beta*C
c
            DO 340, J = 1, N
               DO 330, I = 1, M
                  TEMP = ZERO
                  DO 320, L = 1, K
                     TEMP = TEMP + A( L, I )*CONJG( B( J, L ) )
  320             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  330          CONTINUE
  340       CONTINUE
         ELSE
c
c           Form  C := alpha*A'*B' + beta*C
c
            DO 370, J = 1, N
               DO 360, I = 1, M
                  TEMP = ZERO
                  DO 350, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  350             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  360          CONTINUE
  370       CONTINUE
         END IF
      END IF
c
      RETURN
c
c     End of CGEMM .
c
      END
      SUBROUTINE CHEMM ( SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
c     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO
      INTEGER            M, N, LDA, LDB, LDC
      COMPLEX            ALPHA, BETA
c     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), C( LDC, * )
c     ..
c
c  Purpose
c  =======
c
c  CHEMM  performs one of the matrix-matrix operations
c
c     C := alpha*A*B + beta*C,
c
c  or
c
c     C := alpha*B*A + beta*C,
c
c  where alpha and beta are scalars, A is an hermitian matrix and  B and
c  C are m by n matrices.
c
c  Parameters
c  ==========
c
c  SIDE   - CHARACTER*1.
c           On entry,  SIDE  specifies whether  the  hermitian matrix  A
c           appears on the  left or right  in the  operation as follows:
c
c              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,
c
c              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
c
c           Unchanged on exit.
c
c  UPLO   - CHARACTER*1.
c           On  entry,   UPLO  specifies  whether  the  upper  or  lower
c           triangular  part  of  the  hermitian  matrix   A  is  to  be
c           referenced as follows:
c
c              UPLO = 'U' or 'u'   Only the upper triangular part of the
c                                  hermitian matrix is to be referenced.
c
c              UPLO = 'L' or 'l'   Only the lower triangular part of the
c                                  hermitian matrix is to be referenced.
c
c           Unchanged on exit.
c
c  M      - INTEGER.
c           On entry,  M  specifies the number of rows of the matrix  C.
c           M  must be at least zero.
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry, N specifies the number of columns of the matrix C.
c           N  must be at least zero.
c           Unchanged on exit.
c
c  ALPHA  - COMPLEX         .
c           On entry, ALPHA specifies the scalar alpha.
c           Unchanged on exit.
c
c  A      - COMPLEX          array of DIMENSION ( LDA, ka ), where ka is
c           m  when  SIDE = 'L' or 'l'  and is n  otherwise.
c           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of
c           the array  A  must contain the  hermitian matrix,  such that
c           when  UPLO = 'U' or 'u', the leading m by m upper triangular
c           part of the array  A  must contain the upper triangular part
c           of the  hermitian matrix and the  strictly  lower triangular
c           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
c           the leading  m by m  lower triangular part  of the  array  A
c           must  contain  the  lower triangular part  of the  hermitian
c           matrix and the  strictly upper triangular part of  A  is not
c           referenced.
c           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of
c           the array  A  must contain the  hermitian matrix,  such that
c           when  UPLO = 'U' or 'u', the leading n by n upper triangular
c           part of the array  A  must contain the upper triangular part
c           of the  hermitian matrix and the  strictly  lower triangular
c           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
c           the leading  n by n  lower triangular part  of the  array  A
c           must  contain  the  lower triangular part  of the  hermitian
c           matrix and the  strictly upper triangular part of  A  is not
c           referenced.
c           Note that the imaginary parts  of the diagonal elements need
c           not be set, they are assumed to be zero.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in the  calling (sub) program. When  SIDE = 'L' or 'l'  then
c           LDA must be at least  max( 1, m ), otherwise  LDA must be at
c           least max( 1, n ).
c           Unchanged on exit.
c
c  B      - COMPLEX          array of DIMENSION ( LDB, n ).
c           Before entry, the leading  m by n part of the array  B  must
c           contain the matrix B.
c           Unchanged on exit.
c
c  LDB    - INTEGER.
c           On entry, LDB specifies the first dimension of B as declared
c           in  the  calling  (sub)  program.   LDB  must  be  at  least
c           max( 1, m ).
c           Unchanged on exit.
c
c  BETA   - COMPLEX         .
c           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
c           supplied as zero then C need not be set on input.
c           Unchanged on exit.
c
c  C      - COMPLEX          array of DIMENSION ( LDC, n ).
c           Before entry, the leading  m by n  part of the array  C must
c           contain the matrix  C,  except when  beta  is zero, in which
c           case C need not be set on entry.
c           On exit, the array  C  is overwritten by the  m by n updated
c           matrix.
c
c  LDC    - INTEGER.
c           On entry, LDC specifies the first dimension of C as declared
c           in  the  calling  (sub)  program.   LDC  must  be  at  least
c           max( 1, m ).
c           Unchanged on exit.
c
c
c  Level 3 Blas routine.
c
c  -- Written on 8-February-1989.
c     Jack Dongarra, Argonne National Laboratory.
c     Iain Duff, AERE Harwell.
c     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c     Sven Hammarling, Numerical Algorithms Group Ltd.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX, REAL
c     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, INFO, J, K, NROWA
      COMPLEX            TEMP1, TEMP2
c     .. Parameters ..
      COMPLEX            ONE
      PARAMETER        ( ONE  = ( 1.0E+0, 0.0E+0 ) )
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
c     ..
c     .. Executable Statements ..
c
c     Set NROWA as the number of rows of A.
c
      IF( LSAME( SIDE, 'L' ) )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      UPPER = LSAME( UPLO, 'U' )
c
c     Test the input parameters.
c
      INFO = 0
      IF(      ( .NOT.LSAME( SIDE, 'L' ) ).AND.
     $         ( .NOT.LSAME( SIDE, 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER              ).AND.
     $         ( .NOT.LSAME( UPLO, 'L' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 7
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 9
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 12
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CHEMM ', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
c
c     And when  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
c
c     Start the operations.
c
      IF( LSAME( SIDE, 'L' ) )THEN
c
c        Form  C := alpha*A*B + beta*C.
c
         IF( UPPER )THEN
            DO 70, J = 1, N
               DO 60, I = 1, M
                  TEMP1 = ALPHA*B( I, J )
                  TEMP2 = ZERO
                  DO 50, K = 1, I - 1
                     C( K, J ) = C( K, J ) + TEMP1*A( K, I )
                     TEMP2     = TEMP2     +
     $                           B( K, J )*CONJG(  A( K, I ) )
   50             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = TEMP1*REAL( A( I, I ) ) +
     $                           ALPHA*TEMP2
                  ELSE
                     C( I, J ) = BETA *C( I, J )         +
     $                           TEMP1*REAL( A( I, I ) ) +
     $                           ALPHA*TEMP2
                  END IF
   60          CONTINUE
   70       CONTINUE
         ELSE
            DO 100, J = 1, N
               DO 90, I = M, 1, -1
                  TEMP1 = ALPHA*B( I, J )
                  TEMP2 = ZERO
                  DO 80, K = I + 1, M
                     C( K, J ) = C( K, J ) + TEMP1*A( K, I )
                     TEMP2     = TEMP2     +
     $                           B( K, J )*CONJG(  A( K, I ) )
   80             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = TEMP1*REAL( A( I, I ) ) +
     $                           ALPHA*TEMP2
                  ELSE
                     C( I, J ) = BETA *C( I, J )         +
     $                           TEMP1*REAL( A( I, I ) ) +
     $                           ALPHA*TEMP2
                  END IF
   90          CONTINUE
  100       CONTINUE
         END IF
      ELSE
c
c        Form  C := alpha*B*A + beta*C.
c
         DO 170, J = 1, N
            TEMP1 = ALPHA*REAL( A( J, J ) )
            IF( BETA.EQ.ZERO )THEN
               DO 110, I = 1, M
                  C( I, J ) = TEMP1*B( I, J )
  110          CONTINUE
            ELSE
               DO 120, I = 1, M
                  C( I, J ) = BETA*C( I, J ) + TEMP1*B( I, J )
  120          CONTINUE
            END IF
            DO 140, K = 1, J - 1
               IF( UPPER )THEN
                  TEMP1 = ALPHA*A( K, J )
               ELSE
                  TEMP1 = ALPHA*CONJG( A( J, K ) )
               END IF
               DO 130, I = 1, M
                  C( I, J ) = C( I, J ) + TEMP1*B( I, K )
  130          CONTINUE
  140       CONTINUE
            DO 160, K = J + 1, N
               IF( UPPER )THEN
                  TEMP1 = ALPHA*CONJG( A( J, K ) )
               ELSE
                  TEMP1 = ALPHA*A( K, J )
               END IF
               DO 150, I = 1, M
                  C( I, J ) = C( I, J ) + TEMP1*B( I, K )
  150          CONTINUE
  160       CONTINUE
  170    CONTINUE
      END IF
c
      RETURN
c
c     End of CHEMM .
c
      END
      SUBROUTINE CHER2K( UPLO, TRANS, N, K, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
c     .. Scalar Arguments ..
      CHARACTER*1        UPLO, TRANS
      INTEGER            N, K, LDA, LDB, LDC
      REAL               BETA
      COMPLEX            ALPHA
c     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), C( LDC, * )
c     ..
c
c  Purpose
c  =======
c
c  CHER2K  performs one of the hermitian rank 2k operations
c
c     C := alpha*A*conjg( B' ) + conjg( alpha )*B*conjg( A' ) + beta*C,
c
c  or
c
c     C := alpha*conjg( A' )*B + conjg( alpha )*conjg( B' )*A + beta*C,
c
c  where  alpha and beta  are scalars with  beta  real,  C is an  n by n
c  hermitian matrix and  A and B  are  n by k matrices in the first case
c  and  k by n  matrices in the second case.
c
c  Parameters
c  ==========
c
c  UPLO   - CHARACTER*1.
c           On  entry,   UPLO  specifies  whether  the  upper  or  lower
c           triangular  part  of the  array  C  is to be  referenced  as
c           follows:
c
c              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
c                                  is to be referenced.
c
c              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
c                                  is to be referenced.
c
c           Unchanged on exit.
c
c  TRANS  - CHARACTER*1.
c           On entry,  TRANS  specifies the operation to be performed as
c           follows:
c
c              TRANS = 'N' or 'n'    C := alpha*A*conjg( B' )          +
c                                         conjg( alpha )*B*conjg( A' ) +
c                                         beta*C.
c
c              TRANS = 'C' or 'c'    C := alpha*conjg( A' )*B          +
c                                         conjg( alpha )*conjg( B' )*A +
c                                         beta*C.
c
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry,  N specifies the order of the matrix C.  N must be
c           at least zero.
c           Unchanged on exit.
c
c  K      - INTEGER.
c           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
c           of  columns  of the  matrices  A and B,  and on  entry  with
c           TRANS = 'C' or 'c',  K  specifies  the number of rows of the
c           matrices  A and B.  K must be at least zero.
c           Unchanged on exit.
c
c  ALPHA  - COMPLEX         .
c           On entry, ALPHA specifies the scalar alpha.
c           Unchanged on exit.
c
c  A      - COMPLEX          array of DIMENSION ( LDA, ka ), where ka is
c           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
c           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
c           part of the array  A  must contain the matrix  A,  otherwise
c           the leading  k by n  part of the array  A  must contain  the
c           matrix A.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
c           then  LDA must be at least  max( 1, n ), otherwise  LDA must
c           be at least  max( 1, k ).
c           Unchanged on exit.
c
c  B      - COMPLEX          array of DIMENSION ( LDB, kb ), where kb is
c           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
c           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
c           part of the array  B  must contain the matrix  B,  otherwise
c           the leading  k by n  part of the array  B  must contain  the
c           matrix B.
c           Unchanged on exit.
c
c  LDB    - INTEGER.
c           On entry, LDB specifies the first dimension of B as declared
c           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
c           then  LDB must be at least  max( 1, n ), otherwise  LDB must
c           be at least  max( 1, k ).
c           Unchanged on exit.
c
c  BETA   - REAL            .
c           On entry, BETA specifies the scalar beta.
c           Unchanged on exit.
c
c  C      - COMPLEX          array of DIMENSION ( LDC, n ).
c           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
c           upper triangular part of the array C must contain the upper
c           triangular part  of the  hermitian matrix  and the strictly
c           lower triangular part of C is not referenced.  On exit, the
c           upper triangular part of the array  C is overwritten by the
c           upper triangular part of the updated matrix.
c           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
c           lower triangular part of the array C must contain the lower
c           triangular part  of the  hermitian matrix  and the strictly
c           upper triangular part of C is not referenced.  On exit, the
c           lower triangular part of the array  C is overwritten by the
c           lower triangular part of the updated matrix.
c           Note that the imaginary parts of the diagonal elements need
c           not be set,  they are assumed to be zero,  and on exit they
c           are set to zero.
c
c  LDC    - INTEGER.
c           On entry, LDC specifies the first dimension of C as declared
c           in  the  calling  (sub)  program.   LDC  must  be  at  least
c           max( 1, n ).
c           Unchanged on exit.
c
c
c  Level 3 Blas routine.
c
c  -- Written on 8-February-1989.
c     Jack Dongarra, Argonne National Laboratory.
c     Iain Duff, AERE Harwell.
c     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c     Sven Hammarling, Numerical Algorithms Group Ltd.
c
c  -- Modified 8-Nov-93 to set C(J,J) to REAL( C(J,J) ) when BETA = 1.
c     Ed Anderson, Cray Research Inc.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX, REAL
c     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, INFO, J, L, NROWA
      COMPLEX            TEMP1, TEMP2
c     .. Parameters ..
      REAL               ONE
      PARAMETER        ( ONE  = 1.0E+0 )
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      IF( LSAME( TRANS, 'N' ) )THEN
         NROWA = N
      ELSE
         NROWA = K
      END IF
      UPPER = LSAME( UPLO, 'U' )
c
      INFO = 0
      IF(      ( .NOT.UPPER               ).AND.
     $         ( .NOT.LSAME( UPLO , 'L' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.LSAME( TRANS, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANS, 'C' ) )      )THEN
         INFO = 2
      ELSE IF( N  .LT.0               )THEN
         INFO = 3
      ELSE IF( K  .LT.0               )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 7
      ELSE IF( LDB.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDC.LT.MAX( 1, N     ) )THEN
         INFO = 12
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CHER2K', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( ( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
c
c     And when  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO )THEN
         IF( UPPER )THEN
            IF( BETA.EQ.REAL( ZERO ) )THEN
               DO 20, J = 1, N
                  DO 10, I = 1, J
                     C( I, J ) = ZERO
   10             CONTINUE
   20          CONTINUE
            ELSE
               DO 40, J = 1, N
                  DO 30, I = 1, J - 1
                     C( I, J ) = BETA*C( I, J )
   30             CONTINUE
                  C( J, J ) = BETA*REAL( C( J, J ) )
   40          CONTINUE
            END IF
         ELSE
            IF( BETA.EQ.REAL( ZERO ) )THEN
               DO 60, J = 1, N
                  DO 50, I = J, N
                     C( I, J ) = ZERO
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 80, J = 1, N
                  C( J, J ) = BETA*REAL( C( J, J ) )
                  DO 70, I = J + 1, N
                     C( I, J ) = BETA*C( I, J )
   70             CONTINUE
   80          CONTINUE
            END IF
         END IF
         RETURN
      END IF
c
c     Start the operations.
c
      IF( LSAME( TRANS, 'N' ) )THEN
c
c        Form  C := alpha*A*conjg( B' ) + conjg( alpha )*B*conjg( A' ) +
c                   C.
c
         IF( UPPER )THEN
            DO 130, J = 1, N
               IF( BETA.EQ.REAL( ZERO ) )THEN
                  DO 90, I = 1, J
                     C( I, J ) = ZERO
   90             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 100, I = 1, J - 1
                     C( I, J ) = BETA*C( I, J )
  100             CONTINUE
                  C( J, J ) = BETA*REAL( C( J, J ) )
               ELSE
                  C( J, J ) = REAL( C( J, J ) )
               END IF
               DO 120, L = 1, K
                  IF( ( A( J, L ).NE.ZERO ).OR.
     $                ( B( J, L ).NE.ZERO )     )THEN
                     TEMP1 = ALPHA*CONJG( B( J, L ) )
                     TEMP2 = CONJG( ALPHA*A( J, L ) )
                     DO 110, I = 1, J - 1
                        C( I, J ) = C( I, J ) + A( I, L )*TEMP1 +
     $                                          B( I, L )*TEMP2
  110                CONTINUE
                     C( J, J ) = REAL( C( J, J ) )         +
     $                           REAL( A( J, L )*TEMP1 +
     $                                 B( J, L )*TEMP2   )
                  END IF
  120          CONTINUE
  130       CONTINUE
         ELSE
            DO 180, J = 1, N
               IF( BETA.EQ.REAL( ZERO ) )THEN
                  DO 140, I = J, N
                     C( I, J ) = ZERO
  140             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 150, I = J + 1, N
                     C( I, J ) = BETA*C( I, J )
  150             CONTINUE
                  C( J, J ) = BETA*REAL( C( J, J ) )
               ELSE
                  C( J, J ) = REAL( C( J, J ) )
               END IF
               DO 170, L = 1, K
                  IF( ( A( J, L ).NE.ZERO ).OR.
     $                ( B( J, L ).NE.ZERO )     )THEN
                     TEMP1 = ALPHA*CONJG( B( J, L ) )
                     TEMP2 = CONJG( ALPHA*A( J, L ) )
                     DO 160, I = J + 1, N
                        C( I, J ) = C( I, J ) + A( I, L )*TEMP1 +
     $                                          B( I, L )*TEMP2
  160                CONTINUE
                     C( J, J ) = REAL( C( J, J ) )         +
     $                           REAL( A( J, L )*TEMP1 +
     $                                 B( J, L )*TEMP2   )
                  END IF
  170          CONTINUE
  180       CONTINUE
         END IF
      ELSE
c
c        Form  C := alpha*conjg( A' )*B + conjg( alpha )*conjg( B' )*A +
c                   C.
c
         IF( UPPER )THEN
            DO 210, J = 1, N
               DO 200, I = 1, J
                  TEMP1 = ZERO
                  TEMP2 = ZERO
                  DO 190, L = 1, K
                     TEMP1 = TEMP1 + CONJG( A( L, I ) )*B( L, J )
                     TEMP2 = TEMP2 + CONJG( B( L, I ) )*A( L, J )
  190             CONTINUE
                  IF( I.EQ.J )THEN
                     IF( BETA.EQ.REAL( ZERO ) )THEN
                        C( J, J ) = REAL(        ALPHA  *TEMP1 +
     $                                    CONJG( ALPHA )*TEMP2   )
                     ELSE
                        C( J, J ) = BETA*REAL( C( J, J ) )         +
     $                              REAL(        ALPHA  *TEMP1 +
     $                                    CONJG( ALPHA )*TEMP2   )
                     END IF
                  ELSE
                     IF( BETA.EQ.REAL( ZERO ) )THEN
                        C( I, J ) = ALPHA*TEMP1 + CONJG( ALPHA )*TEMP2
                     ELSE
                        C( I, J ) = BETA *C( I, J ) +
     $                              ALPHA*TEMP1 + CONJG( ALPHA )*TEMP2
                     END IF
                  END IF
  200          CONTINUE
  210       CONTINUE
         ELSE
            DO 240, J = 1, N
               DO 230, I = J, N
                  TEMP1 = ZERO
                  TEMP2 = ZERO
                  DO 220, L = 1, K
                     TEMP1 = TEMP1 + CONJG( A( L, I ) )*B( L, J )
                     TEMP2 = TEMP2 + CONJG( B( L, I ) )*A( L, J )
  220             CONTINUE
                  IF( I.EQ.J )THEN
                     IF( BETA.EQ.REAL( ZERO ) )THEN
                        C( J, J ) = REAL(        ALPHA  *TEMP1 +
     $                                    CONJG( ALPHA )*TEMP2   )
                     ELSE
                        C( J, J ) = BETA*REAL( C( J, J ) )         +
     $                              REAL(        ALPHA  *TEMP1 +
     $                                    CONJG( ALPHA )*TEMP2   )
                     END IF
                  ELSE
                     IF( BETA.EQ.REAL( ZERO ) )THEN
                        C( I, J ) = ALPHA*TEMP1 + CONJG( ALPHA )*TEMP2
                     ELSE
                        C( I, J ) = BETA *C( I, J ) +
     $                              ALPHA*TEMP1 + CONJG( ALPHA )*TEMP2
                     END IF
                  END IF
  230          CONTINUE
  240       CONTINUE
         END IF
      END IF
c
      RETURN
c
c     End of CHER2K.
c
      END
      SUBROUTINE CHERK ( UPLO, TRANS, N, K, ALPHA, A, LDA,
     $                   BETA, C, LDC )
c     .. Scalar Arguments ..
      CHARACTER*1        UPLO, TRANS
      INTEGER            N, K, LDA, LDC
      REAL               ALPHA, BETA
c     .. Array Arguments ..
      COMPLEX            A( LDA, * ), C( LDC, * )
c     ..
c
c  Purpose
c  =======
c
c  CHERK  performs one of the hermitian rank k operations
c
c     C := alpha*A*conjg( A' ) + beta*C,
c
c  or
c
c     C := alpha*conjg( A' )*A + beta*C,
c
c  where  alpha and beta  are  real scalars,  C is an  n by n  hermitian
c  matrix and  A  is an  n by k  matrix in the  first case and a  k by n
c  matrix in the second case.
c
c  Parameters
c  ==========
c
c  UPLO   - CHARACTER*1.
c           On  entry,   UPLO  specifies  whether  the  upper  or  lower
c           triangular  part  of the  array  C  is to be  referenced  as
c           follows:
c
c              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
c                                  is to be referenced.
c
c              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
c                                  is to be referenced.
c
c           Unchanged on exit.
c
c  TRANS  - CHARACTER*1.
c           On entry,  TRANS  specifies the operation to be performed as
c           follows:
c
c              TRANS = 'N' or 'n'   C := alpha*A*conjg( A' ) + beta*C.
c
c              TRANS = 'C' or 'c'   C := alpha*conjg( A' )*A + beta*C.
c
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry,  N specifies the order of the matrix C.  N must be
c           at least zero.
c           Unchanged on exit.
c
c  K      - INTEGER.
c           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
c           of  columns   of  the   matrix   A,   and  on   entry   with
c           TRANS = 'C' or 'c',  K  specifies  the number of rows of the
c           matrix A.  K must be at least zero.
c           Unchanged on exit.
c
c  ALPHA  - REAL            .
c           On entry, ALPHA specifies the scalar alpha.
c           Unchanged on exit.
c
c  A      - COMPLEX          array of DIMENSION ( LDA, ka ), where ka is
c           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
c           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
c           part of the array  A  must contain the matrix  A,  otherwise
c           the leading  k by n  part of the array  A  must contain  the
c           matrix A.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
c           then  LDA must be at least  max( 1, n ), otherwise  LDA must
c           be at least  max( 1, k ).
c           Unchanged on exit.
c
c  BETA   - REAL            .
c           On entry, BETA specifies the scalar beta.
c           Unchanged on exit.
c
c  C      - COMPLEX          array of DIMENSION ( LDC, n ).
c           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
c           upper triangular part of the array C must contain the upper
c           triangular part  of the  hermitian matrix  and the strictly
c           lower triangular part of C is not referenced.  On exit, the
c           upper triangular part of the array  C is overwritten by the
c           upper triangular part of the updated matrix.
c           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
c           lower triangular part of the array C must contain the lower
c           triangular part  of the  hermitian matrix  and the strictly
c           upper triangular part of C is not referenced.  On exit, the
c           lower triangular part of the array  C is overwritten by the
c           lower triangular part of the updated matrix.
c           Note that the imaginary parts of the diagonal elements need
c           not be set,  they are assumed to be zero,  and on exit they
c           are set to zero.
c
c  LDC    - INTEGER.
c           On entry, LDC specifies the first dimension of C as declared
c           in  the  calling  (sub)  program.   LDC  must  be  at  least
c           max( 1, n ).
c           Unchanged on exit.
c
c
c  Level 3 Blas routine.
c
c  -- Written on 8-February-1989.
c     Jack Dongarra, Argonne National Laboratory.
c     Iain Duff, AERE Harwell.
c     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c     Sven Hammarling, Numerical Algorithms Group Ltd.
c
c  -- Modified 8-Nov-93 to set C(J,J) to REAL( C(J,J) ) when BETA = 1.
c     Ed Anderson, Cray Research Inc.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          CMPLX, CONJG, MAX, REAL
c     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, INFO, J, L, NROWA
      REAL               RTEMP
      COMPLEX            TEMP
c     .. Parameters ..
      REAL               ONE ,         ZERO
      PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      IF( LSAME( TRANS, 'N' ) )THEN
         NROWA = N
      ELSE
         NROWA = K
      END IF
      UPPER = LSAME( UPLO, 'U' )
c
      INFO = 0
      IF(      ( .NOT.UPPER               ).AND.
     $         ( .NOT.LSAME( UPLO , 'L' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.LSAME( TRANS, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANS, 'C' ) )      )THEN
         INFO = 2
      ELSE IF( N  .LT.0               )THEN
         INFO = 3
      ELSE IF( K  .LT.0               )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 7
      ELSE IF( LDC.LT.MAX( 1, N     ) )THEN
         INFO = 10
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CHERK ', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( ( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
c
c     And when  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO )THEN
         IF( UPPER )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 20, J = 1, N
                  DO 10, I = 1, J
                     C( I, J ) = ZERO
   10             CONTINUE
   20          CONTINUE
            ELSE
               DO 40, J = 1, N
                  DO 30, I = 1, J - 1
                     C( I, J ) = BETA*C( I, J )
   30             CONTINUE
                  C( J, J ) = BETA*REAL( C( J, J ) )
   40          CONTINUE
            END IF
         ELSE
            IF( BETA.EQ.ZERO )THEN
               DO 60, J = 1, N
                  DO 50, I = J, N
                     C( I, J ) = ZERO
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 80, J = 1, N
                  C( J, J ) = BETA*REAL( C( J, J ) )
                  DO 70, I = J + 1, N
                     C( I, J ) = BETA*C( I, J )
   70             CONTINUE
   80          CONTINUE
            END IF
         END IF
         RETURN
      END IF
c
c     Start the operations.
c
      IF( LSAME( TRANS, 'N' ) )THEN
c
c        Form  C := alpha*A*conjg( A' ) + beta*C.
c
         IF( UPPER )THEN
            DO 130, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 90, I = 1, J
                     C( I, J ) = ZERO
   90             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 100, I = 1, J - 1
                     C( I, J ) = BETA*C( I, J )
  100             CONTINUE
                  C( J, J ) = BETA*REAL( C( J, J ) )
               ELSE
                  C( J, J ) = REAL( C( J, J ) )
               END IF
               DO 120, L = 1, K
                  IF( A( J, L ).NE.CMPLX( ZERO ) )THEN
                     TEMP = ALPHA*CONJG( A( J, L ) )
                     DO 110, I = 1, J - 1
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  110                CONTINUE
                     C( J, J ) = REAL( C( J, J )      ) +
     $                           REAL( TEMP*A( I, L ) )
                  END IF
  120          CONTINUE
  130       CONTINUE
         ELSE
            DO 180, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 140, I = J, N
                     C( I, J ) = ZERO
  140             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  C( J, J ) = BETA*REAL( C( J, J ) )
                  DO 150, I = J + 1, N
                     C( I, J ) = BETA*C( I, J )
  150             CONTINUE
               ELSE
                  C( J, J ) = REAL( C( J, J ) )
               END IF
               DO 170, L = 1, K
                  IF( A( J, L ).NE.CMPLX( ZERO ) )THEN
                     TEMP      = ALPHA*CONJG( A( J, L ) )
                     C( J, J ) = REAL( C( J, J )      )   +
     $                           REAL( TEMP*A( J, L ) )
                     DO 160, I = J + 1, N
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  160                CONTINUE
                  END IF
  170          CONTINUE
  180       CONTINUE
         END IF
      ELSE
c
c        Form  C := alpha*conjg( A' )*A + beta*C.
c
         IF( UPPER )THEN
            DO 220, J = 1, N
               DO 200, I = 1, J - 1
                  TEMP = ZERO
                  DO 190, L = 1, K
                     TEMP = TEMP + CONJG( A( L, I ) )*A( L, J )
  190             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  200          CONTINUE
               RTEMP = ZERO
               DO 210, L = 1, K
                  RTEMP = RTEMP + CONJG( A( L, J ) )*A( L, J )
  210          CONTINUE
               IF( BETA.EQ.ZERO )THEN
                  C( J, J ) = ALPHA*RTEMP
               ELSE
                  C( J, J ) = ALPHA*RTEMP + BETA*REAL( C( J, J ) )
               END IF
  220       CONTINUE
         ELSE
            DO 260, J = 1, N
               RTEMP = ZERO
               DO 230, L = 1, K
                  RTEMP = RTEMP + CONJG( A( L, J ) )*A( L, J )
  230          CONTINUE
               IF( BETA.EQ.ZERO )THEN
                  C( J, J ) = ALPHA*RTEMP
               ELSE
                  C( J, J ) = ALPHA*RTEMP + BETA*REAL( C( J, J ) )
               END IF
               DO 250, I = J + 1, N
                  TEMP = ZERO
                  DO 240, L = 1, K
                     TEMP = TEMP + CONJG( A( L, I ) )*A( L, J )
  240             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  250          CONTINUE
  260       CONTINUE
         END IF
      END IF
c
      RETURN
c
c     End of CHERK .
c
      END
      SUBROUTINE CSYMM ( SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
c     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO
      INTEGER            M, N, LDA, LDB, LDC
      COMPLEX            ALPHA, BETA
c     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), C( LDC, * )
c     ..
c
c  Purpose
c  =======
c
c  CSYMM  performs one of the matrix-matrix operations
c
c     C := alpha*A*B + beta*C,
c
c  or
c
c     C := alpha*B*A + beta*C,
c
c  where  alpha and beta are scalars, A is a symmetric matrix and  B and
c  C are m by n matrices.
c
c  Parameters
c  ==========
c
c  SIDE   - CHARACTER*1.
c           On entry,  SIDE  specifies whether  the  symmetric matrix  A
c           appears on the  left or right  in the  operation as follows:
c
c              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,
c
c              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
c
c           Unchanged on exit.
c
c  UPLO   - CHARACTER*1.
c           On  entry,   UPLO  specifies  whether  the  upper  or  lower
c           triangular  part  of  the  symmetric  matrix   A  is  to  be
c           referenced as follows:
c
c              UPLO = 'U' or 'u'   Only the upper triangular part of the
c                                  symmetric matrix is to be referenced.
c
c              UPLO = 'L' or 'l'   Only the lower triangular part of the
c                                  symmetric matrix is to be referenced.
c
c           Unchanged on exit.
c
c  M      - INTEGER.
c           On entry,  M  specifies the number of rows of the matrix  C.
c           M  must be at least zero.
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry, N specifies the number of columns of the matrix C.
c           N  must be at least zero.
c           Unchanged on exit.
c
c  ALPHA  - COMPLEX         .
c           On entry, ALPHA specifies the scalar alpha.
c           Unchanged on exit.
c
c  A      - COMPLEX          array of DIMENSION ( LDA, ka ), where ka is
c           m  when  SIDE = 'L' or 'l'  and is n  otherwise.
c           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of
c           the array  A  must contain the  symmetric matrix,  such that
c           when  UPLO = 'U' or 'u', the leading m by m upper triangular
c           part of the array  A  must contain the upper triangular part
c           of the  symmetric matrix and the  strictly  lower triangular
c           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
c           the leading  m by m  lower triangular part  of the  array  A
c           must  contain  the  lower triangular part  of the  symmetric
c           matrix and the  strictly upper triangular part of  A  is not
c           referenced.
c           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of
c           the array  A  must contain the  symmetric matrix,  such that
c           when  UPLO = 'U' or 'u', the leading n by n upper triangular
c           part of the array  A  must contain the upper triangular part
c           of the  symmetric matrix and the  strictly  lower triangular
c           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
c           the leading  n by n  lower triangular part  of the  array  A
c           must  contain  the  lower triangular part  of the  symmetric
c           matrix and the  strictly upper triangular part of  A  is not
c           referenced.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in the  calling (sub) program. When  SIDE = 'L' or 'l'  then
c           LDA must be at least  max( 1, m ), otherwise  LDA must be at
c           least max( 1, n ).
c           Unchanged on exit.
c
c  B      - COMPLEX          array of DIMENSION ( LDB, n ).
c           Before entry, the leading  m by n part of the array  B  must
c           contain the matrix B.
c           Unchanged on exit.
c
c  LDB    - INTEGER.
c           On entry, LDB specifies the first dimension of B as declared
c           in  the  calling  (sub)  program.   LDB  must  be  at  least
c           max( 1, m ).
c           Unchanged on exit.
c
c  BETA   - COMPLEX         .
c           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
c           supplied as zero then C need not be set on input.
c           Unchanged on exit.
c
c  C      - COMPLEX          array of DIMENSION ( LDC, n ).
c           Before entry, the leading  m by n  part of the array  C must
c           contain the matrix  C,  except when  beta  is zero, in which
c           case C need not be set on entry.
c           On exit, the array  C  is overwritten by the  m by n updated
c           matrix.
c
c  LDC    - INTEGER.
c           On entry, LDC specifies the first dimension of C as declared
c           in  the  calling  (sub)  program.   LDC  must  be  at  least
c           max( 1, m ).
c           Unchanged on exit.
c
c
c  Level 3 Blas routine.
c
c  -- Written on 8-February-1989.
c     Jack Dongarra, Argonne National Laboratory.
c     Iain Duff, AERE Harwell.
c     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c     Sven Hammarling, Numerical Algorithms Group Ltd.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, INFO, J, K, NROWA
      COMPLEX            TEMP1, TEMP2
c     .. Parameters ..
      COMPLEX            ONE
      PARAMETER        ( ONE  = ( 1.0E+0, 0.0E+0 ) )
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
c     ..
c     .. Executable Statements ..
c
c     Set NROWA as the number of rows of A.
c
      IF( LSAME( SIDE, 'L' ) )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      UPPER = LSAME( UPLO, 'U' )
c
c     Test the input parameters.
c
      INFO = 0
      IF(      ( .NOT.LSAME( SIDE, 'L' ) ).AND.
     $         ( .NOT.LSAME( SIDE, 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER              ).AND.
     $         ( .NOT.LSAME( UPLO, 'L' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 7
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 9
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 12
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CSYMM ', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
c
c     And when  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
c
c     Start the operations.
c
      IF( LSAME( SIDE, 'L' ) )THEN
c
c        Form  C := alpha*A*B + beta*C.
c
         IF( UPPER )THEN
            DO 70, J = 1, N
               DO 60, I = 1, M
                  TEMP1 = ALPHA*B( I, J )
                  TEMP2 = ZERO
                  DO 50, K = 1, I - 1
                     C( K, J ) = C( K, J ) + TEMP1    *A( K, I )
                     TEMP2     = TEMP2     + B( K, J )*A( K, I )
   50             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = TEMP1*A( I, I ) + ALPHA*TEMP2
                  ELSE
                     C( I, J ) = BETA *C( I, J ) +
     $                           TEMP1*A( I, I ) + ALPHA*TEMP2
                  END IF
   60          CONTINUE
   70       CONTINUE
         ELSE
            DO 100, J = 1, N
               DO 90, I = M, 1, -1
                  TEMP1 = ALPHA*B( I, J )
                  TEMP2 = ZERO
                  DO 80, K = I + 1, M
                     C( K, J ) = C( K, J ) + TEMP1    *A( K, I )
                     TEMP2     = TEMP2     + B( K, J )*A( K, I )
   80             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = TEMP1*A( I, I ) + ALPHA*TEMP2
                  ELSE
                     C( I, J ) = BETA *C( I, J ) +
     $                           TEMP1*A( I, I ) + ALPHA*TEMP2
                  END IF
   90          CONTINUE
  100       CONTINUE
         END IF
      ELSE
c
c        Form  C := alpha*B*A + beta*C.
c
         DO 170, J = 1, N
            TEMP1 = ALPHA*A( J, J )
            IF( BETA.EQ.ZERO )THEN
               DO 110, I = 1, M
                  C( I, J ) = TEMP1*B( I, J )
  110          CONTINUE
            ELSE
               DO 120, I = 1, M
                  C( I, J ) = BETA*C( I, J ) + TEMP1*B( I, J )
  120          CONTINUE
            END IF
            DO 140, K = 1, J - 1
               IF( UPPER )THEN
                  TEMP1 = ALPHA*A( K, J )
               ELSE
                  TEMP1 = ALPHA*A( J, K )
               END IF
               DO 130, I = 1, M
                  C( I, J ) = C( I, J ) + TEMP1*B( I, K )
  130          CONTINUE
  140       CONTINUE
            DO 160, K = J + 1, N
               IF( UPPER )THEN
                  TEMP1 = ALPHA*A( J, K )
               ELSE
                  TEMP1 = ALPHA*A( K, J )
               END IF
               DO 150, I = 1, M
                  C( I, J ) = C( I, J ) + TEMP1*B( I, K )
  150          CONTINUE
  160       CONTINUE
  170    CONTINUE
      END IF
c
      RETURN
c
c     End of CSYMM .
c
      END
      SUBROUTINE CSYR2K( UPLO, TRANS, N, K, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
c     .. Scalar Arguments ..
      CHARACTER*1        UPLO, TRANS
      INTEGER            N, K, LDA, LDB, LDC
      COMPLEX            ALPHA, BETA
c     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), C( LDC, * )
c     ..
c
c  Purpose
c  =======
c
c  CSYR2K  performs one of the symmetric rank 2k operations
c
c     C := alpha*A*B' + alpha*B*A' + beta*C,
c
c  or
c
c     C := alpha*A'*B + alpha*B'*A + beta*C,
c
c  where  alpha and beta  are scalars,  C is an  n by n symmetric matrix
c  and  A and B  are  n by k  matrices  in the  first  case  and  k by n
c  matrices in the second case.
c
c  Parameters
c  ==========
c
c  UPLO   - CHARACTER*1.
c           On  entry,   UPLO  specifies  whether  the  upper  or  lower
c           triangular  part  of the  array  C  is to be  referenced  as
c           follows:
c
c              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
c                                  is to be referenced.
c
c              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
c                                  is to be referenced.
c
c           Unchanged on exit.
c
c  TRANS  - CHARACTER*1.
c           On entry,  TRANS  specifies the operation to be performed as
c           follows:
c
c              TRANS = 'N' or 'n'    C := alpha*A*B' + alpha*B*A' +
c                                         beta*C.
c
c              TRANS = 'T' or 't'    C := alpha*A'*B + alpha*B'*A +
c                                         beta*C.
c
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry,  N specifies the order of the matrix C.  N must be
c           at least zero.
c           Unchanged on exit.
c
c  K      - INTEGER.
c           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
c           of  columns  of the  matrices  A and B,  and on  entry  with
c           TRANS = 'T' or 't',  K  specifies  the number of rows of the
c           matrices  A and B.  K must be at least zero.
c           Unchanged on exit.
c
c  ALPHA  - COMPLEX         .
c           On entry, ALPHA specifies the scalar alpha.
c           Unchanged on exit.
c
c  A      - COMPLEX          array of DIMENSION ( LDA, ka ), where ka is
c           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
c           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
c           part of the array  A  must contain the matrix  A,  otherwise
c           the leading  k by n  part of the array  A  must contain  the
c           matrix A.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
c           then  LDA must be at least  max( 1, n ), otherwise  LDA must
c           be at least  max( 1, k ).
c           Unchanged on exit.
c
c  B      - COMPLEX          array of DIMENSION ( LDB, kb ), where kb is
c           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
c           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
c           part of the array  B  must contain the matrix  B,  otherwise
c           the leading  k by n  part of the array  B  must contain  the
c           matrix B.
c           Unchanged on exit.
c
c  LDB    - INTEGER.
c           On entry, LDB specifies the first dimension of B as declared
c           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
c           then  LDB must be at least  max( 1, n ), otherwise  LDB must
c           be at least  max( 1, k ).
c           Unchanged on exit.
c
c  BETA   - COMPLEX         .
c           On entry, BETA specifies the scalar beta.
c           Unchanged on exit.
c
c  C      - COMPLEX          array of DIMENSION ( LDC, n ).
c           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
c           upper triangular part of the array C must contain the upper
c           triangular part  of the  symmetric matrix  and the strictly
c           lower triangular part of C is not referenced.  On exit, the
c           upper triangular part of the array  C is overwritten by the
c           upper triangular part of the updated matrix.
c           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
c           lower triangular part of the array C must contain the lower
c           triangular part  of the  symmetric matrix  and the strictly
c           upper triangular part of C is not referenced.  On exit, the
c           lower triangular part of the array  C is overwritten by the
c           lower triangular part of the updated matrix.
c
c  LDC    - INTEGER.
c           On entry, LDC specifies the first dimension of C as declared
c           in  the  calling  (sub)  program.   LDC  must  be  at  least
c           max( 1, n ).
c           Unchanged on exit.
c
c
c  Level 3 Blas routine.
c
c  -- Written on 8-February-1989.
c     Jack Dongarra, Argonne National Laboratory.
c     Iain Duff, AERE Harwell.
c     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c     Sven Hammarling, Numerical Algorithms Group Ltd.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, INFO, J, L, NROWA
      COMPLEX            TEMP1, TEMP2
c     .. Parameters ..
      COMPLEX            ONE
      PARAMETER        ( ONE  = ( 1.0E+0, 0.0E+0 ) )
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      IF( LSAME( TRANS, 'N' ) )THEN
         NROWA = N
      ELSE
         NROWA = K
      END IF
      UPPER = LSAME( UPLO, 'U' )
c
      INFO = 0
      IF(      ( .NOT.UPPER               ).AND.
     $         ( .NOT.LSAME( UPLO , 'L' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.LSAME( TRANS, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANS, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( N  .LT.0               )THEN
         INFO = 3
      ELSE IF( K  .LT.0               )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 7
      ELSE IF( LDB.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDC.LT.MAX( 1, N     ) )THEN
         INFO = 12
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CSYR2K', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( ( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
c
c     And when  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO )THEN
         IF( UPPER )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 20, J = 1, N
                  DO 10, I = 1, J
                     C( I, J ) = ZERO
   10             CONTINUE
   20          CONTINUE
            ELSE
               DO 40, J = 1, N
                  DO 30, I = 1, J
                     C( I, J ) = BETA*C( I, J )
   30             CONTINUE
   40          CONTINUE
            END IF
         ELSE
            IF( BETA.EQ.ZERO )THEN
               DO 60, J = 1, N
                  DO 50, I = J, N
                     C( I, J ) = ZERO
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 80, J = 1, N
                  DO 70, I = J, N
                     C( I, J ) = BETA*C( I, J )
   70             CONTINUE
   80          CONTINUE
            END IF
         END IF
         RETURN
      END IF
c
c     Start the operations.
c
      IF( LSAME( TRANS, 'N' ) )THEN
c
c        Form  C := alpha*A*B' + alpha*B*A' + C.
c
         IF( UPPER )THEN
            DO 130, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 90, I = 1, J
                     C( I, J ) = ZERO
   90             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 100, I = 1, J
                     C( I, J ) = BETA*C( I, J )
  100             CONTINUE
               END IF
               DO 120, L = 1, K
                  IF( ( A( J, L ).NE.ZERO ).OR.
     $                ( B( J, L ).NE.ZERO )     )THEN
                     TEMP1 = ALPHA*B( J, L )
                     TEMP2 = ALPHA*A( J, L )
                     DO 110, I = 1, J
                        C( I, J ) = C( I, J ) + A( I, L )*TEMP1 +
     $                                          B( I, L )*TEMP2
  110                CONTINUE
                  END IF
  120          CONTINUE
  130       CONTINUE
         ELSE
            DO 180, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 140, I = J, N
                     C( I, J ) = ZERO
  140             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 150, I = J, N
                     C( I, J ) = BETA*C( I, J )
  150             CONTINUE
               END IF
               DO 170, L = 1, K
                  IF( ( A( J, L ).NE.ZERO ).OR.
     $                ( B( J, L ).NE.ZERO )     )THEN
                     TEMP1 = ALPHA*B( J, L )
                     TEMP2 = ALPHA*A( J, L )
                     DO 160, I = J, N
                        C( I, J ) = C( I, J ) + A( I, L )*TEMP1 +
     $                                          B( I, L )*TEMP2
  160                CONTINUE
                  END IF
  170          CONTINUE
  180       CONTINUE
         END IF
      ELSE
c
c        Form  C := alpha*A'*B + alpha*B'*A + C.
c
         IF( UPPER )THEN
            DO 210, J = 1, N
               DO 200, I = 1, J
                  TEMP1 = ZERO
                  TEMP2 = ZERO
                  DO 190, L = 1, K
                     TEMP1 = TEMP1 + A( L, I )*B( L, J )
                     TEMP2 = TEMP2 + B( L, I )*A( L, J )
  190             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP1 + ALPHA*TEMP2
                  ELSE
                     C( I, J ) = BETA *C( I, J ) +
     $                           ALPHA*TEMP1 + ALPHA*TEMP2
                  END IF
  200          CONTINUE
  210       CONTINUE
         ELSE
            DO 240, J = 1, N
               DO 230, I = J, N
                  TEMP1 = ZERO
                  TEMP2 = ZERO
                  DO 220, L = 1, K
                     TEMP1 = TEMP1 + A( L, I )*B( L, J )
                     TEMP2 = TEMP2 + B( L, I )*A( L, J )
  220             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP1 + ALPHA*TEMP2
                  ELSE
                     C( I, J ) = BETA *C( I, J ) +
     $                           ALPHA*TEMP1 + ALPHA*TEMP2
                  END IF
  230          CONTINUE
  240       CONTINUE
         END IF
      END IF
c
      RETURN
c
c     End of CSYR2K.
c
      END
      SUBROUTINE CSYRK ( UPLO, TRANS, N, K, ALPHA, A, LDA,
     $                   BETA, C, LDC )
c     .. Scalar Arguments ..
      CHARACTER*1        UPLO, TRANS
      INTEGER            N, K, LDA, LDC
      COMPLEX            ALPHA, BETA
c     .. Array Arguments ..
      COMPLEX            A( LDA, * ), C( LDC, * )
c     ..
c
c  Purpose
c  =======
c
c  CSYRK  performs one of the symmetric rank k operations
c
c     C := alpha*A*A' + beta*C,
c
c  or
c
c     C := alpha*A'*A + beta*C,
c
c  where  alpha and beta  are scalars,  C is an  n by n symmetric matrix
c  and  A  is an  n by k  matrix in the first case and a  k by n  matrix
c  in the second case.
c
c  Parameters
c  ==========
c
c  UPLO   - CHARACTER*1.
c           On  entry,   UPLO  specifies  whether  the  upper  or  lower
c           triangular  part  of the  array  C  is to be  referenced  as
c           follows:
c
c              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
c                                  is to be referenced.
c
c              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
c                                  is to be referenced.
c
c           Unchanged on exit.
c
c  TRANS  - CHARACTER*1.
c           On entry,  TRANS  specifies the operation to be performed as
c           follows:
c
c              TRANS = 'N' or 'n'   C := alpha*A*A' + beta*C.
c
c              TRANS = 'T' or 't'   C := alpha*A'*A + beta*C.
c
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry,  N specifies the order of the matrix C.  N must be
c           at least zero.
c           Unchanged on exit.
c
c  K      - INTEGER.
c           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
c           of  columns   of  the   matrix   A,   and  on   entry   with
c           TRANS = 'T' or 't',  K  specifies  the number of rows of the
c           matrix A.  K must be at least zero.
c           Unchanged on exit.
c
c  ALPHA  - COMPLEX         .
c           On entry, ALPHA specifies the scalar alpha.
c           Unchanged on exit.
c
c  A      - COMPLEX          array of DIMENSION ( LDA, ka ), where ka is
c           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
c           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
c           part of the array  A  must contain the matrix  A,  otherwise
c           the leading  k by n  part of the array  A  must contain  the
c           matrix A.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
c           then  LDA must be at least  max( 1, n ), otherwise  LDA must
c           be at least  max( 1, k ).
c           Unchanged on exit.
c
c  BETA   - COMPLEX         .
c           On entry, BETA specifies the scalar beta.
c           Unchanged on exit.
c
c  C      - COMPLEX          array of DIMENSION ( LDC, n ).
c           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
c           upper triangular part of the array C must contain the upper
c           triangular part  of the  symmetric matrix  and the strictly
c           lower triangular part of C is not referenced.  On exit, the
c           upper triangular part of the array  C is overwritten by the
c           upper triangular part of the updated matrix.
c           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
c           lower triangular part of the array C must contain the lower
c           triangular part  of the  symmetric matrix  and the strictly
c           upper triangular part of C is not referenced.  On exit, the
c           lower triangular part of the array  C is overwritten by the
c           lower triangular part of the updated matrix.
c
c  LDC    - INTEGER.
c           On entry, LDC specifies the first dimension of C as declared
c           in  the  calling  (sub)  program.   LDC  must  be  at  least
c           max( 1, n ).
c           Unchanged on exit.
c
c
c  Level 3 Blas routine.
c
c  -- Written on 8-February-1989.
c     Jack Dongarra, Argonne National Laboratory.
c     Iain Duff, AERE Harwell.
c     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c     Sven Hammarling, Numerical Algorithms Group Ltd.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, INFO, J, L, NROWA
      COMPLEX            TEMP
c     .. Parameters ..
      COMPLEX            ONE
      PARAMETER        ( ONE  = ( 1.0E+0, 0.0E+0 ) )
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      IF( LSAME( TRANS, 'N' ) )THEN
         NROWA = N
      ELSE
         NROWA = K
      END IF
      UPPER = LSAME( UPLO, 'U' )
c
      INFO = 0
      IF(      ( .NOT.UPPER               ).AND.
     $         ( .NOT.LSAME( UPLO , 'L' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.LSAME( TRANS, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANS, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( N  .LT.0               )THEN
         INFO = 3
      ELSE IF( K  .LT.0               )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 7
      ELSE IF( LDC.LT.MAX( 1, N     ) )THEN
         INFO = 10
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CSYRK ', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( ( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
c
c     And when  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO )THEN
         IF( UPPER )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 20, J = 1, N
                  DO 10, I = 1, J
                     C( I, J ) = ZERO
   10             CONTINUE
   20          CONTINUE
            ELSE
               DO 40, J = 1, N
                  DO 30, I = 1, J
                     C( I, J ) = BETA*C( I, J )
   30             CONTINUE
   40          CONTINUE
            END IF
         ELSE
            IF( BETA.EQ.ZERO )THEN
               DO 60, J = 1, N
                  DO 50, I = J, N
                     C( I, J ) = ZERO
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 80, J = 1, N
                  DO 70, I = J, N
                     C( I, J ) = BETA*C( I, J )
   70             CONTINUE
   80          CONTINUE
            END IF
         END IF
         RETURN
      END IF
c
c     Start the operations.
c
      IF( LSAME( TRANS, 'N' ) )THEN
c
c        Form  C := alpha*A*A' + beta*C.
c
         IF( UPPER )THEN
            DO 130, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 90, I = 1, J
                     C( I, J ) = ZERO
   90             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 100, I = 1, J
                     C( I, J ) = BETA*C( I, J )
  100             CONTINUE
               END IF
               DO 120, L = 1, K
                  IF( A( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*A( J, L )
                     DO 110, I = 1, J
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  110                CONTINUE
                  END IF
  120          CONTINUE
  130       CONTINUE
         ELSE
            DO 180, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 140, I = J, N
                     C( I, J ) = ZERO
  140             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 150, I = J, N
                     C( I, J ) = BETA*C( I, J )
  150             CONTINUE
               END IF
               DO 170, L = 1, K
                  IF( A( J, L ).NE.ZERO )THEN
                     TEMP      = ALPHA*A( J, L )
                     DO 160, I = J, N
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  160                CONTINUE
                  END IF
  170          CONTINUE
  180       CONTINUE
         END IF
      ELSE
c
c        Form  C := alpha*A'*A + beta*C.
c
         IF( UPPER )THEN
            DO 210, J = 1, N
               DO 200, I = 1, J
                  TEMP = ZERO
                  DO 190, L = 1, K
                     TEMP = TEMP + A( L, I )*A( L, J )
  190             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  200          CONTINUE
  210       CONTINUE
         ELSE
            DO 240, J = 1, N
               DO 230, I = J, N
                  TEMP = ZERO
                  DO 220, L = 1, K
                     TEMP = TEMP + A( L, I )*A( L, J )
  220             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  230          CONTINUE
  240       CONTINUE
         END IF
      END IF
c
      RETURN
c
c     End of CSYRK .
c
      END
      SUBROUTINE CTRMM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   B, LDB )
c     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      COMPLEX            ALPHA
c     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * )
c     ..
c
c  Purpose
c  =======
c
c  CTRMM  performs one of the matrix-matrix operations
c
c     B := alpha*op( A )*B,   or   B := alpha*B*op( A )
c
c  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
c  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
c
c     op( A ) = A   or   op( A ) = A'   or   op( A ) = conjg( A' ).
c
c  Parameters
c  ==========
c
c  SIDE   - CHARACTER*1.
c           On entry,  SIDE specifies whether  op( A ) multiplies B from
c           the left or right as follows:
c
c              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
c
c              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
c
c           Unchanged on exit.
c
c  UPLO   - CHARACTER*1.
c           On entry, UPLO specifies whether the matrix A is an upper or
c           lower triangular matrix as follows:
c
c              UPLO = 'U' or 'u'   A is an upper triangular matrix.
c
c              UPLO = 'L' or 'l'   A is a lower triangular matrix.
c
c           Unchanged on exit.
c
c  TRANSA - CHARACTER*1.
c           On entry, TRANSA specifies the form of op( A ) to be used in
c           the matrix multiplication as follows:
c
c              TRANSA = 'N' or 'n'   op( A ) = A.
c
c              TRANSA = 'T' or 't'   op( A ) = A'.
c
c              TRANSA = 'C' or 'c'   op( A ) = conjg( A' ).
c
c           Unchanged on exit.
c
c  DIAG   - CHARACTER*1.
c           On entry, DIAG specifies whether or not A is unit triangular
c           as follows:
c
c              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
c
c              DIAG = 'N' or 'n'   A is not assumed to be unit
c                                  triangular.
c
c           Unchanged on exit.
c
c  M      - INTEGER.
c           On entry, M specifies the number of rows of B. M must be at
c           least zero.
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry, N specifies the number of columns of B.  N must be
c           at least zero.
c           Unchanged on exit.
c
c  ALPHA  - COMPLEX         .
c           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
c           zero then  A is not referenced and  B need not be set before
c           entry.
c           Unchanged on exit.
c
c  A      - COMPLEX          array of DIMENSION ( LDA, k ), where k is m
c           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
c           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
c           upper triangular part of the array  A must contain the upper
c           triangular matrix  and the strictly lower triangular part of
c           A is not referenced.
c           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
c           lower triangular part of the array  A must contain the lower
c           triangular matrix  and the strictly upper triangular part of
c           A is not referenced.
c           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
c           A  are not referenced either,  but are assumed to be  unity.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
c           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
c           then LDA must be at least max( 1, n ).
c           Unchanged on exit.
c
c  B      - COMPLEX          array of DIMENSION ( LDB, n ).
c           Before entry,  the leading  m by n part of the array  B must
c           contain the matrix  B,  and  on exit  is overwritten  by the
c           transformed matrix.
c
c  LDB    - INTEGER.
c           On entry, LDB specifies the first dimension of B as declared
c           in  the  calling  (sub)  program.   LDB  must  be  at  least
c           max( 1, m ).
c           Unchanged on exit.
c
c
c  Level 3 Blas routine.
c
c  -- Written on 8-February-1989.
c     Jack Dongarra, Argonne National Laboratory.
c     Iain Duff, AERE Harwell.
c     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c     Sven Hammarling, Numerical Algorithms Group Ltd.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX
c     .. Local Scalars ..
      LOGICAL            LSIDE, NOCONJ, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      COMPLEX            TEMP
c     .. Parameters ..
      COMPLEX            ONE
      PARAMETER        ( ONE  = ( 1.0E+0, 0.0E+0 ) )
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOCONJ = LSAME( TRANSA, 'T' )
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
c
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.
     $         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.
     $         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
     $         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CTRMM ', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( N.EQ.0 )
     $   RETURN
c
c     And when  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
c
c     Start the operations.
c
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
c
c           Form  B := alpha*A*B.
c
            IF( UPPER )THEN
               DO 50, J = 1, N
                  DO 40, K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*B( K, J )
                        DO 30, I = 1, K - 1
                           B( I, J ) = B( I, J ) + TEMP*A( I, K )
   30                   CONTINUE
                        IF( NOUNIT )
     $                     TEMP = TEMP*A( K, K )
                        B( K, J ) = TEMP
                     END IF
   40             CONTINUE
   50          CONTINUE
            ELSE
               DO 80, J = 1, N
                  DO 70 K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        TEMP      = ALPHA*B( K, J )
                        B( K, J ) = TEMP
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )*A( K, K )
                        DO 60, I = K + 1, M
                           B( I, J ) = B( I, J ) + TEMP*A( I, K )
   60                   CONTINUE
                     END IF
   70             CONTINUE
   80          CONTINUE
            END IF
         ELSE
c
c           Form  B := alpha*A'*B   or   B := alpha*conjg( A' )*B.
c
            IF( UPPER )THEN
               DO 120, J = 1, N
                  DO 110, I = M, 1, -1
                     TEMP = B( I, J )
                     IF( NOCONJ )THEN
                        IF( NOUNIT )
     $                     TEMP = TEMP*A( I, I )
                        DO 90, K = 1, I - 1
                           TEMP = TEMP + A( K, I )*B( K, J )
   90                   CONTINUE
                     ELSE
                        IF( NOUNIT )
     $                     TEMP = TEMP*CONJG( A( I, I ) )
                        DO 100, K = 1, I - 1
                           TEMP = TEMP + CONJG( A( K, I ) )*B( K, J )
  100                   CONTINUE
                     END IF
                     B( I, J ) = ALPHA*TEMP
  110             CONTINUE
  120          CONTINUE
            ELSE
               DO 160, J = 1, N
                  DO 150, I = 1, M
                     TEMP = B( I, J )
                     IF( NOCONJ )THEN
                        IF( NOUNIT )
     $                     TEMP = TEMP*A( I, I )
                        DO 130, K = I + 1, M
                           TEMP = TEMP + A( K, I )*B( K, J )
  130                   CONTINUE
                     ELSE
                        IF( NOUNIT )
     $                     TEMP = TEMP*CONJG( A( I, I ) )
                        DO 140, K = I + 1, M
                           TEMP = TEMP + CONJG( A( K, I ) )*B( K, J )
  140                   CONTINUE
                     END IF
                     B( I, J ) = ALPHA*TEMP
  150             CONTINUE
  160          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
c
c           Form  B := alpha*B*A.
c
            IF( UPPER )THEN
               DO 200, J = N, 1, -1
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 170, I = 1, M
                     B( I, J ) = TEMP*B( I, J )
  170             CONTINUE
                  DO 190, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*A( K, J )
                        DO 180, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  180                   CONTINUE
                     END IF
  190             CONTINUE
  200          CONTINUE
            ELSE
               DO 240, J = 1, N
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 210, I = 1, M
                     B( I, J ) = TEMP*B( I, J )
  210             CONTINUE
                  DO 230, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*A( K, J )
                        DO 220, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  220                   CONTINUE
                     END IF
  230             CONTINUE
  240          CONTINUE
            END IF
         ELSE
c
c           Form  B := alpha*B*A'   or   B := alpha*B*conjg( A' ).
c
            IF( UPPER )THEN
               DO 280, K = 1, N
                  DO 260, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        IF( NOCONJ )THEN
                           TEMP = ALPHA*A( J, K )
                        ELSE
                           TEMP = ALPHA*CONJG( A( J, K ) )
                        END IF
                        DO 250, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  250                   CONTINUE
                     END IF
  260             CONTINUE
                  TEMP = ALPHA
                  IF( NOUNIT )THEN
                     IF( NOCONJ )THEN
                        TEMP = TEMP*A( K, K )
                     ELSE
                        TEMP = TEMP*CONJG( A( K, K ) )
                     END IF
                  END IF
                  IF( TEMP.NE.ONE )THEN
                     DO 270, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  270                CONTINUE
                  END IF
  280          CONTINUE
            ELSE
               DO 320, K = N, 1, -1
                  DO 300, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        IF( NOCONJ )THEN
                           TEMP = ALPHA*A( J, K )
                        ELSE
                           TEMP = ALPHA*CONJG( A( J, K ) )
                        END IF
                        DO 290, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  290                   CONTINUE
                     END IF
  300             CONTINUE
                  TEMP = ALPHA
                  IF( NOUNIT )THEN
                     IF( NOCONJ )THEN
                        TEMP = TEMP*A( K, K )
                     ELSE
                        TEMP = TEMP*CONJG( A( K, K ) )
                     END IF
                  END IF
                  IF( TEMP.NE.ONE )THEN
                     DO 310, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  310                CONTINUE
                  END IF
  320          CONTINUE
            END IF
         END IF
      END IF
c
      RETURN
c
c     End of CTRMM .
c
      END
      SUBROUTINE CTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   B, LDB )
c     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      COMPLEX            ALPHA
c     .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * )
c     ..
c
c  Purpose
c  =======
c
c  CTRSM  solves one of the matrix equations
c
c     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
c
c  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
c  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
c
c     op( A ) = A   or   op( A ) = A'   or   op( A ) = conjg( A' ).
c
c  The matrix X is overwritten on B.
c
c  Parameters
c  ==========
c
c  SIDE   - CHARACTER*1.
c           On entry, SIDE specifies whether op( A ) appears on the left
c           or right of X as follows:
c
c              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
c
c              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
c
c           Unchanged on exit.
c
c  UPLO   - CHARACTER*1.
c           On entry, UPLO specifies whether the matrix A is an upper or
c           lower triangular matrix as follows:
c
c              UPLO = 'U' or 'u'   A is an upper triangular matrix.
c
c              UPLO = 'L' or 'l'   A is a lower triangular matrix.
c
c           Unchanged on exit.
c
c  TRANSA - CHARACTER*1.
c           On entry, TRANSA specifies the form of op( A ) to be used in
c           the matrix multiplication as follows:
c
c              TRANSA = 'N' or 'n'   op( A ) = A.
c
c              TRANSA = 'T' or 't'   op( A ) = A'.
c
c              TRANSA = 'C' or 'c'   op( A ) = conjg( A' ).
c
c           Unchanged on exit.
c
c  DIAG   - CHARACTER*1.
c           On entry, DIAG specifies whether or not A is unit triangular
c           as follows:
c
c              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
c
c              DIAG = 'N' or 'n'   A is not assumed to be unit
c                                  triangular.
c
c           Unchanged on exit.
c
c  M      - INTEGER.
c           On entry, M specifies the number of rows of B. M must be at
c           least zero.
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry, N specifies the number of columns of B.  N must be
c           at least zero.
c           Unchanged on exit.
c
c  ALPHA  - COMPLEX         .
c           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
c           zero then  A is not referenced and  B need not be set before
c           entry.
c           Unchanged on exit.
c
c  A      - COMPLEX          array of DIMENSION ( LDA, k ), where k is m
c           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
c           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
c           upper triangular part of the array  A must contain the upper
c           triangular matrix  and the strictly lower triangular part of
c           A is not referenced.
c           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
c           lower triangular part of the array  A must contain the lower
c           triangular matrix  and the strictly upper triangular part of
c           A is not referenced.
c           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
c           A  are not referenced either,  but are assumed to be  unity.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
c           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
c           then LDA must be at least max( 1, n ).
c           Unchanged on exit.
c
c  B      - COMPLEX          array of DIMENSION ( LDB, n ).
c           Before entry,  the leading  m by n part of the array  B must
c           contain  the  right-hand  side  matrix  B,  and  on exit  is
c           overwritten by the solution matrix  X.
c
c  LDB    - INTEGER.
c           On entry, LDB specifies the first dimension of B as declared
c           in  the  calling  (sub)  program.   LDB  must  be  at  least
c           max( 1, m ).
c           Unchanged on exit.
c
c
c  Level 3 Blas routine.
c
c  -- Written on 8-February-1989.
c     Jack Dongarra, Argonne National Laboratory.
c     Iain Duff, AERE Harwell.
c     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c     Sven Hammarling, Numerical Algorithms Group Ltd.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          CONJG, MAX
c     .. Local Scalars ..
      LOGICAL            LSIDE, NOCONJ, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      COMPLEX            TEMP
c     .. Parameters ..
      COMPLEX            ONE
      PARAMETER        ( ONE  = ( 1.0E+0, 0.0E+0 ) )
      COMPLEX            ZERO
      PARAMETER        ( ZERO = ( 0.0E+0, 0.0E+0 ) )
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOCONJ = LSAME( TRANSA, 'T' )
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
c
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.
     $         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.
     $         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
     $         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'CTRSM ', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( N.EQ.0 )
     $   RETURN
c
c     And when  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
c
c     Start the operations.
c
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
c
c           Form  B := alpha*inv( A )*B.
c
            IF( UPPER )THEN
               DO 60, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 30, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   30                CONTINUE
                  END IF
                  DO 50, K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 40, I = 1, K - 1
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   40                   CONTINUE
                     END IF
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 100, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 70, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   70                CONTINUE
                  END IF
                  DO 90 K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 80, I = K + 1, M
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   80                   CONTINUE
                     END IF
   90             CONTINUE
  100          CONTINUE
            END IF
         ELSE
c
c           Form  B := alpha*inv( A' )*B
c           or    B := alpha*inv( conjg( A' ) )*B.
c
            IF( UPPER )THEN
               DO 140, J = 1, N
                  DO 130, I = 1, M
                     TEMP = ALPHA*B( I, J )
                     IF( NOCONJ )THEN
                        DO 110, K = 1, I - 1
                           TEMP = TEMP - A( K, I )*B( K, J )
  110                   CONTINUE
                        IF( NOUNIT )
     $                     TEMP = TEMP/A( I, I )
                     ELSE
                        DO 120, K = 1, I - 1
                           TEMP = TEMP - CONJG( A( K, I ) )*B( K, J )
  120                   CONTINUE
                        IF( NOUNIT )
     $                     TEMP = TEMP/CONJG( A( I, I ) )
                     END IF
                     B( I, J ) = TEMP
  130             CONTINUE
  140          CONTINUE
            ELSE
               DO 180, J = 1, N
                  DO 170, I = M, 1, -1
                     TEMP = ALPHA*B( I, J )
                     IF( NOCONJ )THEN
                        DO 150, K = I + 1, M
                           TEMP = TEMP - A( K, I )*B( K, J )
  150                   CONTINUE
                        IF( NOUNIT )
     $                     TEMP = TEMP/A( I, I )
                     ELSE
                        DO 160, K = I + 1, M
                           TEMP = TEMP - CONJG( A( K, I ) )*B( K, J )
  160                   CONTINUE
                        IF( NOUNIT )
     $                     TEMP = TEMP/CONJG( A( I, I ) )
                     END IF
                     B( I, J ) = TEMP
  170             CONTINUE
  180          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
c
c           Form  B := alpha*B*inv( A ).
c
            IF( UPPER )THEN
               DO 230, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 190, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  190                CONTINUE
                  END IF
                  DO 210, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 200, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  200                   CONTINUE
                     END IF
  210             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 220, I = 1, M
                        B( I, J ) = TEMP*B( I, J )
  220                CONTINUE
                  END IF
  230          CONTINUE
            ELSE
               DO 280, J = N, 1, -1
                  IF( ALPHA.NE.ONE )THEN
                     DO 240, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  240                CONTINUE
                  END IF
                  DO 260, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 250, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  250                   CONTINUE
                     END IF
  260             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 270, I = 1, M
                       B( I, J ) = TEMP*B( I, J )
  270                CONTINUE
                  END IF
  280          CONTINUE
            END IF
         ELSE
c
c           Form  B := alpha*B*inv( A' )
c           or    B := alpha*B*inv( conjg( A' ) ).
c
            IF( UPPER )THEN
               DO 330, K = N, 1, -1
                  IF( NOUNIT )THEN
                     IF( NOCONJ )THEN
                        TEMP = ONE/A( K, K )
                     ELSE
                        TEMP = ONE/CONJG( A( K, K ) )
                     END IF
                     DO 290, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  290                CONTINUE
                  END IF
                  DO 310, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        IF( NOCONJ )THEN
                           TEMP = A( J, K )
                        ELSE
                           TEMP = CONJG( A( J, K ) )
                        END IF
                        DO 300, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  300                   CONTINUE
                     END IF
  310             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 320, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  320                CONTINUE
                  END IF
  330          CONTINUE
            ELSE
               DO 380, K = 1, N
                  IF( NOUNIT )THEN
                     IF( NOCONJ )THEN
                        TEMP = ONE/A( K, K )
                     ELSE
                        TEMP = ONE/CONJG( A( K, K ) )
                     END IF
                     DO 340, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  340                CONTINUE
                  END IF
                  DO 360, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        IF( NOCONJ )THEN
                           TEMP = A( J, K )
                        ELSE
                           TEMP = CONJG( A( J, K ) )
                        END IF
                        DO 350, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  350                   CONTINUE
                     END IF
  360             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 370, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  370                CONTINUE
                  END IF
  380          CONTINUE
            END IF
         END IF
      END IF
c
      RETURN
c
c     End of CTRSM .
c
      END
      SUBROUTINE DGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
c     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
c     ..
c
c  Purpose
c  =======
c
c  DGEMM  performs one of the matrix-matrix operations
c
c     C := alpha*op( A )*op( B ) + beta*C,
c
c  where  op( X ) is one of
c
c     op( X ) = X   or   op( X ) = X',
c
c  alpha and beta are scalars, and A, B and C are matrices, with op( A )
c  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
c
c  Parameters
c  ==========
c
c  TRANSA - CHARACTER*1.
c           On entry, TRANSA specifies the form of op( A ) to be used in
c           the matrix multiplication as follows:
c
c              TRANSA = 'N' or 'n',  op( A ) = A.
c
c              TRANSA = 'T' or 't',  op( A ) = A'.
c
c              TRANSA = 'C' or 'c',  op( A ) = A'.
c
c           Unchanged on exit.
c
c  TRANSB - CHARACTER*1.
c           On entry, TRANSB specifies the form of op( B ) to be used in
c           the matrix multiplication as follows:
c
c              TRANSB = 'N' or 'n',  op( B ) = B.
c
c              TRANSB = 'T' or 't',  op( B ) = B'.
c
c              TRANSB = 'C' or 'c',  op( B ) = B'.
c
c           Unchanged on exit.
c
c  M      - INTEGER.
c           On entry,  M  specifies  the number  of rows  of the  matrix
c           op( A )  and of the  matrix  C.  M  must  be at least  zero.
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry,  N  specifies the number  of columns of the matrix
c           op( B ) and the number of columns of the matrix C. N must be
c           at least zero.
c           Unchanged on exit.
c
c  K      - INTEGER.
c           On entry,  K  specifies  the number of columns of the matrix
c           op( A ) and the number of rows of the matrix op( B ). K must
c           be at least  zero.
c           Unchanged on exit.
c
c  ALPHA  - DOUBLE PRECISION.
c           On entry, ALPHA specifies the scalar alpha.
c           Unchanged on exit.
c
c  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
c           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
c           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
c           part of the array  A  must contain the matrix  A,  otherwise
c           the leading  k by m  part of the array  A  must contain  the
c           matrix A.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
c           LDA must be at least  max( 1, m ), otherwise  LDA must be at
c           least  max( 1, k ).
c           Unchanged on exit.
c
c  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
c           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
c           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
c           part of the array  B  must contain the matrix  B,  otherwise
c           the leading  n by k  part of the array  B  must contain  the
c           matrix B.
c           Unchanged on exit.
c
c  LDB    - INTEGER.
c           On entry, LDB specifies the first dimension of B as declared
c           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
c           LDB must be at least  max( 1, k ), otherwise  LDB must be at
c           least  max( 1, n ).
c           Unchanged on exit.
c
c  BETA   - DOUBLE PRECISION.
c           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
c           supplied as zero then C need not be set on input.
c           Unchanged on exit.
c
c  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
c           Before entry, the leading  m by n  part of the array  C must
c           contain the matrix  C,  except when  beta  is zero, in which
c           case C need not be set on entry.
c           On exit, the array  C  is overwritten by the  m by n  matrix
c           ( alpha*op( A )*op( B ) + beta*C ).
c
c  LDC    - INTEGER.
c           On entry, LDC specifies the first dimension of C as declared
c           in  the  calling  (sub)  program.   LDC  must  be  at  least
c           max( 1, m ).
c           Unchanged on exit.
c
c
c  Level 3 Blas routine.
c
c  -- Written on 8-February-1989.
c     Jack Dongarra, Argonne National Laboratory.
c     Iain Duff, AERE Harwell.
c     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c     Sven Hammarling, Numerical Algorithms Group Ltd.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     .. Local Scalars ..
      LOGICAL            NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      DOUBLE PRECISION   TEMP
c     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
c     ..
c     .. Executable Statements ..
c
c     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
c     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
c     and  columns of  A  and the  number of  rows  of  B  respectively.
c
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
c
c     Test the input parameters.
c
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.
     $         ( .NOT.LSAME( TRANSB, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DGEMM ', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
c
c     And if  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
c
c     Start the operations.
c
      IF( NOTB )THEN
         IF( NOTA )THEN
c
c           Form  C := alpha*A*B + beta*C.
c
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE
c
c           Form  C := alpha*A'*B + beta*C
c
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
      ELSE
         IF( NOTA )THEN
c
c           Form  C := alpha*A*B' + beta*C
c
            DO 170, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 130, I = 1, M
                     C( I, J ) = ZERO
  130             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 140, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  140             CONTINUE
               END IF
               DO 160, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 150, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  150                CONTINUE
                  END IF
  160          CONTINUE
  170       CONTINUE
         ELSE
c
c           Form  C := alpha*A'*B' + beta*C
c
            DO 200, J = 1, N
               DO 190, I = 1, M
                  TEMP = ZERO
                  DO 180, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  180             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  190          CONTINUE
  200       CONTINUE
         END IF
      END IF
c
      RETURN
c
c     End of DGEMM .
c
      END
      SUBROUTINE DSYMM ( SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
c     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO
      INTEGER            M, N, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
c     ..
c
c  Purpose
c  =======
c
c  DSYMM  performs one of the matrix-matrix operations
c
c     C := alpha*A*B + beta*C,
c
c  or
c
c     C := alpha*B*A + beta*C,
c
c  where alpha and beta are scalars,  A is a symmetric matrix and  B and
c  C are  m by n matrices.
c
c  Parameters
c  ==========
c
c  SIDE   - CHARACTER*1.
c           On entry,  SIDE  specifies whether  the  symmetric matrix  A
c           appears on the  left or right  in the  operation as follows:
c
c              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,
c
c              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
c
c           Unchanged on exit.
c
c  UPLO   - CHARACTER*1.
c           On  entry,   UPLO  specifies  whether  the  upper  or  lower
c           triangular  part  of  the  symmetric  matrix   A  is  to  be
c           referenced as follows:
c
c              UPLO = 'U' or 'u'   Only the upper triangular part of the
c                                  symmetric matrix is to be referenced.
c
c              UPLO = 'L' or 'l'   Only the lower triangular part of the
c                                  symmetric matrix is to be referenced.
c
c           Unchanged on exit.
c
c  M      - INTEGER.
c           On entry,  M  specifies the number of rows of the matrix  C.
c           M  must be at least zero.
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry, N specifies the number of columns of the matrix C.
c           N  must be at least zero.
c           Unchanged on exit.
c
c  ALPHA  - DOUBLE PRECISION.
c           On entry, ALPHA specifies the scalar alpha.
c           Unchanged on exit.
c
c  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
c           m  when  SIDE = 'L' or 'l'  and is  n otherwise.
c           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of
c           the array  A  must contain the  symmetric matrix,  such that
c           when  UPLO = 'U' or 'u', the leading m by m upper triangular
c           part of the array  A  must contain the upper triangular part
c           of the  symmetric matrix and the  strictly  lower triangular
c           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
c           the leading  m by m  lower triangular part  of the  array  A
c           must  contain  the  lower triangular part  of the  symmetric
c           matrix and the  strictly upper triangular part of  A  is not
c           referenced.
c           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of
c           the array  A  must contain the  symmetric matrix,  such that
c           when  UPLO = 'U' or 'u', the leading n by n upper triangular
c           part of the array  A  must contain the upper triangular part
c           of the  symmetric matrix and the  strictly  lower triangular
c           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
c           the leading  n by n  lower triangular part  of the  array  A
c           must  contain  the  lower triangular part  of the  symmetric
c           matrix and the  strictly upper triangular part of  A  is not
c           referenced.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
c           LDA must be at least  max( 1, m ), otherwise  LDA must be at
c           least  max( 1, n ).
c           Unchanged on exit.
c
c  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
c           Before entry, the leading  m by n part of the array  B  must
c           contain the matrix B.
c           Unchanged on exit.
c
c  LDB    - INTEGER.
c           On entry, LDB specifies the first dimension of B as declared
c           in  the  calling  (sub)  program.   LDB  must  be  at  least
c           max( 1, m ).
c           Unchanged on exit.
c
c  BETA   - DOUBLE PRECISION.
c           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
c           supplied as zero then C need not be set on input.
c           Unchanged on exit.
c
c  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
c           Before entry, the leading  m by n  part of the array  C must
c           contain the matrix  C,  except when  beta  is zero, in which
c           case C need not be set on entry.
c           On exit, the array  C  is overwritten by the  m by n updated
c           matrix.
c
c  LDC    - INTEGER.
c           On entry, LDC specifies the first dimension of C as declared
c           in  the  calling  (sub)  program.   LDC  must  be  at  least
c           max( 1, m ).
c           Unchanged on exit.
c
c
c  Level 3 Blas routine.
c
c  -- Written on 8-February-1989.
c     Jack Dongarra, Argonne National Laboratory.
c     Iain Duff, AERE Harwell.
c     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c     Sven Hammarling, Numerical Algorithms Group Ltd.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, INFO, J, K, NROWA
      DOUBLE PRECISION   TEMP1, TEMP2
c     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
c     ..
c     .. Executable Statements ..
c
c     Set NROWA as the number of rows of A.
c
      IF( LSAME( SIDE, 'L' ) )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      UPPER = LSAME( UPLO, 'U' )
c
c     Test the input parameters.
c
      INFO = 0
      IF(      ( .NOT.LSAME( SIDE, 'L' ) ).AND.
     $         ( .NOT.LSAME( SIDE, 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER              ).AND.
     $         ( .NOT.LSAME( UPLO, 'L' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 7
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 9
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 12
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DSYMM ', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
c
c     And when  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
c
c     Start the operations.
c
      IF( LSAME( SIDE, 'L' ) )THEN
c
c        Form  C := alpha*A*B + beta*C.
c
         IF( UPPER )THEN
            DO 70, J = 1, N
               DO 60, I = 1, M
                  TEMP1 = ALPHA*B( I, J )
                  TEMP2 = ZERO
                  DO 50, K = 1, I - 1
                     C( K, J ) = C( K, J ) + TEMP1    *A( K, I )
                     TEMP2     = TEMP2     + B( K, J )*A( K, I )
   50             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = TEMP1*A( I, I ) + ALPHA*TEMP2
                  ELSE
                     C( I, J ) = BETA *C( I, J ) +
     $                           TEMP1*A( I, I ) + ALPHA*TEMP2
                  END IF
   60          CONTINUE
   70       CONTINUE
         ELSE
            DO 100, J = 1, N
               DO 90, I = M, 1, -1
                  TEMP1 = ALPHA*B( I, J )
                  TEMP2 = ZERO
                  DO 80, K = I + 1, M
                     C( K, J ) = C( K, J ) + TEMP1    *A( K, I )
                     TEMP2     = TEMP2     + B( K, J )*A( K, I )
   80             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = TEMP1*A( I, I ) + ALPHA*TEMP2
                  ELSE
                     C( I, J ) = BETA *C( I, J ) +
     $                           TEMP1*A( I, I ) + ALPHA*TEMP2
                  END IF
   90          CONTINUE
  100       CONTINUE
         END IF
      ELSE
c
c        Form  C := alpha*B*A + beta*C.
c
         DO 170, J = 1, N
            TEMP1 = ALPHA*A( J, J )
            IF( BETA.EQ.ZERO )THEN
               DO 110, I = 1, M
                  C( I, J ) = TEMP1*B( I, J )
  110          CONTINUE
            ELSE
               DO 120, I = 1, M
                  C( I, J ) = BETA*C( I, J ) + TEMP1*B( I, J )
  120          CONTINUE
            END IF
            DO 140, K = 1, J - 1
               IF( UPPER )THEN
                  TEMP1 = ALPHA*A( K, J )
               ELSE
                  TEMP1 = ALPHA*A( J, K )
               END IF
               DO 130, I = 1, M
                  C( I, J ) = C( I, J ) + TEMP1*B( I, K )
  130          CONTINUE
  140       CONTINUE
            DO 160, K = J + 1, N
               IF( UPPER )THEN
                  TEMP1 = ALPHA*A( J, K )
               ELSE
                  TEMP1 = ALPHA*A( K, J )
               END IF
               DO 150, I = 1, M
                  C( I, J ) = C( I, J ) + TEMP1*B( I, K )
  150          CONTINUE
  160       CONTINUE
  170    CONTINUE
      END IF
c
      RETURN
c
c     End of DSYMM .
c
      END
      SUBROUTINE DSYR2K( UPLO, TRANS, N, K, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
c     .. Scalar Arguments ..
      CHARACTER*1        UPLO, TRANS
      INTEGER            N, K, LDA, LDB, LDC
      DOUBLE PRECISION   ALPHA, BETA
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), C( LDC, * )
c     ..
c
c  Purpose
c  =======
c
c  DSYR2K  performs one of the symmetric rank 2k operations
c
c     C := alpha*A*B' + alpha*B*A' + beta*C,
c
c  or
c
c     C := alpha*A'*B + alpha*B'*A + beta*C,
c
c  where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
c  and  A and B  are  n by k  matrices  in the  first  case  and  k by n
c  matrices in the second case.
c
c  Parameters
c  ==========
c
c  UPLO   - CHARACTER*1.
c           On  entry,   UPLO  specifies  whether  the  upper  or  lower
c           triangular  part  of the  array  C  is to be  referenced  as
c           follows:
c
c              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
c                                  is to be referenced.
c
c              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
c                                  is to be referenced.
c
c           Unchanged on exit.
c
c  TRANS  - CHARACTER*1.
c           On entry,  TRANS  specifies the operation to be performed as
c           follows:
c
c              TRANS = 'N' or 'n'   C := alpha*A*B' + alpha*B*A' +
c                                        beta*C.
c
c              TRANS = 'T' or 't'   C := alpha*A'*B + alpha*B'*A +
c                                        beta*C.
c
c              TRANS = 'C' or 'c'   C := alpha*A'*B + alpha*B'*A +
c                                        beta*C.
c
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry,  N specifies the order of the matrix C.  N must be
c           at least zero.
c           Unchanged on exit.
c
c  K      - INTEGER.
c           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
c           of  columns  of the  matrices  A and B,  and on  entry  with
c           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
c           of rows of the matrices  A and B.  K must be at least  zero.
c           Unchanged on exit.
c
c  ALPHA  - DOUBLE PRECISION.
c           On entry, ALPHA specifies the scalar alpha.
c           Unchanged on exit.
c
c  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
c           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
c           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
c           part of the array  A  must contain the matrix  A,  otherwise
c           the leading  k by n  part of the array  A  must contain  the
c           matrix A.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
c           then  LDA must be at least  max( 1, n ), otherwise  LDA must
c           be at least  max( 1, k ).
c           Unchanged on exit.
c
c  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
c           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
c           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
c           part of the array  B  must contain the matrix  B,  otherwise
c           the leading  k by n  part of the array  B  must contain  the
c           matrix B.
c           Unchanged on exit.
c
c  LDB    - INTEGER.
c           On entry, LDB specifies the first dimension of B as declared
c           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
c           then  LDB must be at least  max( 1, n ), otherwise  LDB must
c           be at least  max( 1, k ).
c           Unchanged on exit.
c
c  BETA   - DOUBLE PRECISION.
c           On entry, BETA specifies the scalar beta.
c           Unchanged on exit.
c
c  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
c           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
c           upper triangular part of the array C must contain the upper
c           triangular part  of the  symmetric matrix  and the strictly
c           lower triangular part of C is not referenced.  On exit, the
c           upper triangular part of the array  C is overwritten by the
c           upper triangular part of the updated matrix.
c           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
c           lower triangular part of the array C must contain the lower
c           triangular part  of the  symmetric matrix  and the strictly
c           upper triangular part of C is not referenced.  On exit, the
c           lower triangular part of the array  C is overwritten by the
c           lower triangular part of the updated matrix.
c
c  LDC    - INTEGER.
c           On entry, LDC specifies the first dimension of C as declared
c           in  the  calling  (sub)  program.   LDC  must  be  at  least
c           max( 1, n ).
c           Unchanged on exit.
c
c
c  Level 3 Blas routine.
c
c
c  -- Written on 8-February-1989.
c     Jack Dongarra, Argonne National Laboratory.
c     Iain Duff, AERE Harwell.
c     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c     Sven Hammarling, Numerical Algorithms Group Ltd.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, INFO, J, L, NROWA
      DOUBLE PRECISION   TEMP1, TEMP2
c     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      IF( LSAME( TRANS, 'N' ) )THEN
         NROWA = N
      ELSE
         NROWA = K
      END IF
      UPPER = LSAME( UPLO, 'U' )
c
      INFO = 0
      IF(      ( .NOT.UPPER               ).AND.
     $         ( .NOT.LSAME( UPLO , 'L' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.LSAME( TRANS, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANS, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANS, 'C' ) )      )THEN
         INFO = 2
      ELSE IF( N  .LT.0               )THEN
         INFO = 3
      ELSE IF( K  .LT.0               )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 7
      ELSE IF( LDB.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDC.LT.MAX( 1, N     ) )THEN
         INFO = 12
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DSYR2K', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( ( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
c
c     And when  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO )THEN
         IF( UPPER )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 20, J = 1, N
                  DO 10, I = 1, J
                     C( I, J ) = ZERO
   10             CONTINUE
   20          CONTINUE
            ELSE
               DO 40, J = 1, N
                  DO 30, I = 1, J
                     C( I, J ) = BETA*C( I, J )
   30             CONTINUE
   40          CONTINUE
            END IF
         ELSE
            IF( BETA.EQ.ZERO )THEN
               DO 60, J = 1, N
                  DO 50, I = J, N
                     C( I, J ) = ZERO
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 80, J = 1, N
                  DO 70, I = J, N
                     C( I, J ) = BETA*C( I, J )
   70             CONTINUE
   80          CONTINUE
            END IF
         END IF
         RETURN
      END IF
c
c     Start the operations.
c
      IF( LSAME( TRANS, 'N' ) )THEN
c
c        Form  C := alpha*A*B' + alpha*B*A' + C.
c
         IF( UPPER )THEN
            DO 130, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 90, I = 1, J
                     C( I, J ) = ZERO
   90             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 100, I = 1, J
                     C( I, J ) = BETA*C( I, J )
  100             CONTINUE
               END IF
               DO 120, L = 1, K
                  IF( ( A( J, L ).NE.ZERO ).OR.
     $                ( B( J, L ).NE.ZERO )     )THEN
                     TEMP1 = ALPHA*B( J, L )
                     TEMP2 = ALPHA*A( J, L )
                     DO 110, I = 1, J
                        C( I, J ) = C( I, J ) +
     $                              A( I, L )*TEMP1 + B( I, L )*TEMP2
  110                CONTINUE
                  END IF
  120          CONTINUE
  130       CONTINUE
         ELSE
            DO 180, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 140, I = J, N
                     C( I, J ) = ZERO
  140             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 150, I = J, N
                     C( I, J ) = BETA*C( I, J )
  150             CONTINUE
               END IF
               DO 170, L = 1, K
                  IF( ( A( J, L ).NE.ZERO ).OR.
     $                ( B( J, L ).NE.ZERO )     )THEN
                     TEMP1 = ALPHA*B( J, L )
                     TEMP2 = ALPHA*A( J, L )
                     DO 160, I = J, N
                        C( I, J ) = C( I, J ) +
     $                              A( I, L )*TEMP1 + B( I, L )*TEMP2
  160                CONTINUE
                  END IF
  170          CONTINUE
  180       CONTINUE
         END IF
      ELSE
c
c        Form  C := alpha*A'*B + alpha*B'*A + C.
c
         IF( UPPER )THEN
            DO 210, J = 1, N
               DO 200, I = 1, J
                  TEMP1 = ZERO
                  TEMP2 = ZERO
                  DO 190, L = 1, K
                     TEMP1 = TEMP1 + A( L, I )*B( L, J )
                     TEMP2 = TEMP2 + B( L, I )*A( L, J )
  190             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP1 + ALPHA*TEMP2
                  ELSE
                     C( I, J ) = BETA *C( I, J ) +
     $                           ALPHA*TEMP1 + ALPHA*TEMP2
                  END IF
  200          CONTINUE
  210       CONTINUE
         ELSE
            DO 240, J = 1, N
               DO 230, I = J, N
                  TEMP1 = ZERO
                  TEMP2 = ZERO
                  DO 220, L = 1, K
                     TEMP1 = TEMP1 + A( L, I )*B( L, J )
                     TEMP2 = TEMP2 + B( L, I )*A( L, J )
  220             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP1 + ALPHA*TEMP2
                  ELSE
                     C( I, J ) = BETA *C( I, J ) +
     $                           ALPHA*TEMP1 + ALPHA*TEMP2
                  END IF
  230          CONTINUE
  240       CONTINUE
         END IF
      END IF
c
      RETURN
c
c     End of DSYR2K.
c
      END
      SUBROUTINE DSYRK ( UPLO, TRANS, N, K, ALPHA, A, LDA,
     $                   BETA, C, LDC )
c     .. Scalar Arguments ..
      CHARACTER*1        UPLO, TRANS
      INTEGER            N, K, LDA, LDC
      DOUBLE PRECISION   ALPHA, BETA
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), C( LDC, * )
c     ..
c
c  Purpose
c  =======
c
c  DSYRK  performs one of the symmetric rank k operations
c
c     C := alpha*A*A' + beta*C,
c
c  or
c
c     C := alpha*A'*A + beta*C,
c
c  where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
c  and  A  is an  n by k  matrix in the first case and a  k by n  matrix
c  in the second case.
c
c  Parameters
c  ==========
c
c  UPLO   - CHARACTER*1.
c           On  entry,   UPLO  specifies  whether  the  upper  or  lower
c           triangular  part  of the  array  C  is to be  referenced  as
c           follows:
c
c              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
c                                  is to be referenced.
c
c              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
c                                  is to be referenced.
c
c           Unchanged on exit.
c
c  TRANS  - CHARACTER*1.
c           On entry,  TRANS  specifies the operation to be performed as
c           follows:
c
c              TRANS = 'N' or 'n'   C := alpha*A*A' + beta*C.
c
c              TRANS = 'T' or 't'   C := alpha*A'*A + beta*C.
c
c              TRANS = 'C' or 'c'   C := alpha*A'*A + beta*C.
c
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry,  N specifies the order of the matrix C.  N must be
c           at least zero.
c           Unchanged on exit.
c
c  K      - INTEGER.
c           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
c           of  columns   of  the   matrix   A,   and  on   entry   with
c           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
c           of rows of the matrix  A.  K must be at least zero.
c           Unchanged on exit.
c
c  ALPHA  - DOUBLE PRECISION.
c           On entry, ALPHA specifies the scalar alpha.
c           Unchanged on exit.
c
c  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
c           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
c           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
c           part of the array  A  must contain the matrix  A,  otherwise
c           the leading  k by n  part of the array  A  must contain  the
c           matrix A.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
c           then  LDA must be at least  max( 1, n ), otherwise  LDA must
c           be at least  max( 1, k ).
c           Unchanged on exit.
c
c  BETA   - DOUBLE PRECISION.
c           On entry, BETA specifies the scalar beta.
c           Unchanged on exit.
c
c  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
c           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
c           upper triangular part of the array C must contain the upper
c           triangular part  of the  symmetric matrix  and the strictly
c           lower triangular part of C is not referenced.  On exit, the
c           upper triangular part of the array  C is overwritten by the
c           upper triangular part of the updated matrix.
c           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
c           lower triangular part of the array C must contain the lower
c           triangular part  of the  symmetric matrix  and the strictly
c           upper triangular part of C is not referenced.  On exit, the
c           lower triangular part of the array  C is overwritten by the
c           lower triangular part of the updated matrix.
c
c  LDC    - INTEGER.
c           On entry, LDC specifies the first dimension of C as declared
c           in  the  calling  (sub)  program.   LDC  must  be  at  least
c           max( 1, n ).
c           Unchanged on exit.
c
c
c  Level 3 Blas routine.
c
c  -- Written on 8-February-1989.
c     Jack Dongarra, Argonne National Laboratory.
c     Iain Duff, AERE Harwell.
c     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c     Sven Hammarling, Numerical Algorithms Group Ltd.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, INFO, J, L, NROWA
      DOUBLE PRECISION   TEMP
c     .. Parameters ..
      DOUBLE PRECISION   ONE ,         ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      IF( LSAME( TRANS, 'N' ) )THEN
         NROWA = N
      ELSE
         NROWA = K
      END IF
      UPPER = LSAME( UPLO, 'U' )
c
      INFO = 0
      IF(      ( .NOT.UPPER               ).AND.
     $         ( .NOT.LSAME( UPLO , 'L' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.LSAME( TRANS, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANS, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANS, 'C' ) )      )THEN
         INFO = 2
      ELSE IF( N  .LT.0               )THEN
         INFO = 3
      ELSE IF( K  .LT.0               )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 7
      ELSE IF( LDC.LT.MAX( 1, N     ) )THEN
         INFO = 10
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DSYRK ', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( ( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
c
c     And when  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO )THEN
         IF( UPPER )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 20, J = 1, N
                  DO 10, I = 1, J
                     C( I, J ) = ZERO
   10             CONTINUE
   20          CONTINUE
            ELSE
               DO 40, J = 1, N
                  DO 30, I = 1, J
                     C( I, J ) = BETA*C( I, J )
   30             CONTINUE
   40          CONTINUE
            END IF
         ELSE
            IF( BETA.EQ.ZERO )THEN
               DO 60, J = 1, N
                  DO 50, I = J, N
                     C( I, J ) = ZERO
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 80, J = 1, N
                  DO 70, I = J, N
                     C( I, J ) = BETA*C( I, J )
   70             CONTINUE
   80          CONTINUE
            END IF
         END IF
         RETURN
      END IF
c
c     Start the operations.
c
      IF( LSAME( TRANS, 'N' ) )THEN
c
c        Form  C := alpha*A*A' + beta*C.
c
         IF( UPPER )THEN
            DO 130, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 90, I = 1, J
                     C( I, J ) = ZERO
   90             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 100, I = 1, J
                     C( I, J ) = BETA*C( I, J )
  100             CONTINUE
               END IF
               DO 120, L = 1, K
                  IF( A( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*A( J, L )
                     DO 110, I = 1, J
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  110                CONTINUE
                  END IF
  120          CONTINUE
  130       CONTINUE
         ELSE
            DO 180, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 140, I = J, N
                     C( I, J ) = ZERO
  140             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 150, I = J, N
                     C( I, J ) = BETA*C( I, J )
  150             CONTINUE
               END IF
               DO 170, L = 1, K
                  IF( A( J, L ).NE.ZERO )THEN
                     TEMP      = ALPHA*A( J, L )
                     DO 160, I = J, N
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  160                CONTINUE
                  END IF
  170          CONTINUE
  180       CONTINUE
         END IF
      ELSE
c
c        Form  C := alpha*A'*A + beta*C.
c
         IF( UPPER )THEN
            DO 210, J = 1, N
               DO 200, I = 1, J
                  TEMP = ZERO
                  DO 190, L = 1, K
                     TEMP = TEMP + A( L, I )*A( L, J )
  190             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  200          CONTINUE
  210       CONTINUE
         ELSE
            DO 240, J = 1, N
               DO 230, I = J, N
                  TEMP = ZERO
                  DO 220, L = 1, K
                     TEMP = TEMP + A( L, I )*A( L, J )
  220             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  230          CONTINUE
  240       CONTINUE
         END IF
      END IF
c
      RETURN
c
c     End of DSYRK .
c
      END
      SUBROUTINE DTRMM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   B, LDB )
c     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      DOUBLE PRECISION   ALPHA
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
c     ..
c
c  Purpose
c  =======
c
c  DTRMM  performs one of the matrix-matrix operations
c
c     B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
c
c  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
c  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
c
c     op( A ) = A   or   op( A ) = A'.
c
c  Parameters
c  ==========
c
c  SIDE   - CHARACTER*1.
c           On entry,  SIDE specifies whether  op( A ) multiplies B from
c           the left or right as follows:
c
c              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
c
c              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
c
c           Unchanged on exit.
c
c  UPLO   - CHARACTER*1.
c           On entry, UPLO specifies whether the matrix A is an upper or
c           lower triangular matrix as follows:
c
c              UPLO = 'U' or 'u'   A is an upper triangular matrix.
c
c              UPLO = 'L' or 'l'   A is a lower triangular matrix.
c
c           Unchanged on exit.
c
c  TRANSA - CHARACTER*1.
c           On entry, TRANSA specifies the form of op( A ) to be used in
c           the matrix multiplication as follows:
c
c              TRANSA = 'N' or 'n'   op( A ) = A.
c
c              TRANSA = 'T' or 't'   op( A ) = A'.
c
c              TRANSA = 'C' or 'c'   op( A ) = A'.
c
c           Unchanged on exit.
c
c  DIAG   - CHARACTER*1.
c           On entry, DIAG specifies whether or not A is unit triangular
c           as follows:
c
c              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
c
c              DIAG = 'N' or 'n'   A is not assumed to be unit
c                                  triangular.
c
c           Unchanged on exit.
c
c  M      - INTEGER.
c           On entry, M specifies the number of rows of B. M must be at
c           least zero.
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry, N specifies the number of columns of B.  N must be
c           at least zero.
c           Unchanged on exit.
c
c  ALPHA  - DOUBLE PRECISION.
c           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
c           zero then  A is not referenced and  B need not be set before
c           entry.
c           Unchanged on exit.
c
c  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
c           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
c           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
c           upper triangular part of the array  A must contain the upper
c           triangular matrix  and the strictly lower triangular part of
c           A is not referenced.
c           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
c           lower triangular part of the array  A must contain the lower
c           triangular matrix  and the strictly upper triangular part of
c           A is not referenced.
c           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
c           A  are not referenced either,  but are assumed to be  unity.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
c           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
c           then LDA must be at least max( 1, n ).
c           Unchanged on exit.
c
c  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
c           Before entry,  the leading  m by n part of the array  B must
c           contain the matrix  B,  and  on exit  is overwritten  by the
c           transformed matrix.
c
c  LDB    - INTEGER.
c           On entry, LDB specifies the first dimension of B as declared
c           in  the  calling  (sub)  program.   LDB  must  be  at  least
c           max( 1, m ).
c           Unchanged on exit.
c
c
c  Level 3 Blas routine.
c
c  -- Written on 8-February-1989.
c     Jack Dongarra, Argonne National Laboratory.
c     Iain Duff, AERE Harwell.
c     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c     Sven Hammarling, Numerical Algorithms Group Ltd.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     .. Local Scalars ..
      LOGICAL            LSIDE, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      DOUBLE PRECISION   TEMP
c     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
c
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.
     $         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.
     $         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
     $         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRMM ', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( N.EQ.0 )
     $   RETURN
c
c     And when  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
c
c     Start the operations.
c
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
c
c           Form  B := alpha*A*B.
c
            IF( UPPER )THEN
               DO 50, J = 1, N
                  DO 40, K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*B( K, J )
                        DO 30, I = 1, K - 1
                           B( I, J ) = B( I, J ) + TEMP*A( I, K )
   30                   CONTINUE
                        IF( NOUNIT )
     $                     TEMP = TEMP*A( K, K )
                        B( K, J ) = TEMP
                     END IF
   40             CONTINUE
   50          CONTINUE
            ELSE
               DO 80, J = 1, N
                  DO 70 K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        TEMP      = ALPHA*B( K, J )
                        B( K, J ) = TEMP
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )*A( K, K )
                        DO 60, I = K + 1, M
                           B( I, J ) = B( I, J ) + TEMP*A( I, K )
   60                   CONTINUE
                     END IF
   70             CONTINUE
   80          CONTINUE
            END IF
         ELSE
c
c           Form  B := alpha*A'*B.
c
            IF( UPPER )THEN
               DO 110, J = 1, N
                  DO 100, I = M, 1, -1
                     TEMP = B( I, J )
                     IF( NOUNIT )
     $                  TEMP = TEMP*A( I, I )
                     DO 90, K = 1, I - 1
                        TEMP = TEMP + A( K, I )*B( K, J )
   90                CONTINUE
                     B( I, J ) = ALPHA*TEMP
  100             CONTINUE
  110          CONTINUE
            ELSE
               DO 140, J = 1, N
                  DO 130, I = 1, M
                     TEMP = B( I, J )
                     IF( NOUNIT )
     $                  TEMP = TEMP*A( I, I )
                     DO 120, K = I + 1, M
                        TEMP = TEMP + A( K, I )*B( K, J )
  120                CONTINUE
                     B( I, J ) = ALPHA*TEMP
  130             CONTINUE
  140          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
c
c           Form  B := alpha*B*A.
c
            IF( UPPER )THEN
               DO 180, J = N, 1, -1
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 150, I = 1, M
                     B( I, J ) = TEMP*B( I, J )
  150             CONTINUE
                  DO 170, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*A( K, J )
                        DO 160, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  160                   CONTINUE
                     END IF
  170             CONTINUE
  180          CONTINUE
            ELSE
               DO 220, J = 1, N
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 190, I = 1, M
                     B( I, J ) = TEMP*B( I, J )
  190             CONTINUE
                  DO 210, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*A( K, J )
                        DO 200, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  200                   CONTINUE
                     END IF
  210             CONTINUE
  220          CONTINUE
            END IF
         ELSE
c
c           Form  B := alpha*B*A'.
c
            IF( UPPER )THEN
               DO 260, K = 1, N
                  DO 240, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = ALPHA*A( J, K )
                        DO 230, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  230                   CONTINUE
                     END IF
  240             CONTINUE
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( K, K )
                  IF( TEMP.NE.ONE )THEN
                     DO 250, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  250                CONTINUE
                  END IF
  260          CONTINUE
            ELSE
               DO 300, K = N, 1, -1
                  DO 280, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = ALPHA*A( J, K )
                        DO 270, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  270                   CONTINUE
                     END IF
  280             CONTINUE
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( K, K )
                  IF( TEMP.NE.ONE )THEN
                     DO 290, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  290                CONTINUE
                  END IF
  300          CONTINUE
            END IF
         END IF
      END IF
c
      RETURN
c
c     End of DTRMM .
c
      END
      SUBROUTINE DTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   B, LDB )
c     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      DOUBLE PRECISION   ALPHA
c     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
c     ..
c
c  Purpose
c  =======
c
c  DTRSM  solves one of the matrix equations
c
c     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
c
c  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
c  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
c
c     op( A ) = A   or   op( A ) = A'.
c
c  The matrix X is overwritten on B.
c
c  Parameters
c  ==========
c
c  SIDE   - CHARACTER*1.
c           On entry, SIDE specifies whether op( A ) appears on the left
c           or right of X as follows:
c
c              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
c
c              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
c
c           Unchanged on exit.
c
c  UPLO   - CHARACTER*1.
c           On entry, UPLO specifies whether the matrix A is an upper or
c           lower triangular matrix as follows:
c
c              UPLO = 'U' or 'u'   A is an upper triangular matrix.
c
c              UPLO = 'L' or 'l'   A is a lower triangular matrix.
c
c           Unchanged on exit.
c
c  TRANSA - CHARACTER*1.
c           On entry, TRANSA specifies the form of op( A ) to be used in
c           the matrix multiplication as follows:
c
c              TRANSA = 'N' or 'n'   op( A ) = A.
c
c              TRANSA = 'T' or 't'   op( A ) = A'.
c
c              TRANSA = 'C' or 'c'   op( A ) = A'.
c
c           Unchanged on exit.
c
c  DIAG   - CHARACTER*1.
c           On entry, DIAG specifies whether or not A is unit triangular
c           as follows:
c
c              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
c
c              DIAG = 'N' or 'n'   A is not assumed to be unit
c                                  triangular.
c
c           Unchanged on exit.
c
c  M      - INTEGER.
c           On entry, M specifies the number of rows of B. M must be at
c           least zero.
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry, N specifies the number of columns of B.  N must be
c           at least zero.
c           Unchanged on exit.
c
c  ALPHA  - DOUBLE PRECISION.
c           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
c           zero then  A is not referenced and  B need not be set before
c           entry.
c           Unchanged on exit.
c
c  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
c           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
c           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
c           upper triangular part of the array  A must contain the upper
c           triangular matrix  and the strictly lower triangular part of
c           A is not referenced.
c           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
c           lower triangular part of the array  A must contain the lower
c           triangular matrix  and the strictly upper triangular part of
c           A is not referenced.
c           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
c           A  are not referenced either,  but are assumed to be  unity.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
c           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
c           then LDA must be at least max( 1, n ).
c           Unchanged on exit.
c
c  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
c           Before entry,  the leading  m by n part of the array  B must
c           contain  the  right-hand  side  matrix  B,  and  on exit  is
c           overwritten by the solution matrix  X.
c
c  LDB    - INTEGER.
c           On entry, LDB specifies the first dimension of B as declared
c           in  the  calling  (sub)  program.   LDB  must  be  at  least
c           max( 1, m ).
c           Unchanged on exit.
c
c
c  Level 3 Blas routine.
c
c
c  -- Written on 8-February-1989.
c     Jack Dongarra, Argonne National Laboratory.
c     Iain Duff, AERE Harwell.
c     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c     Sven Hammarling, Numerical Algorithms Group Ltd.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     .. Local Scalars ..
      LOGICAL            LSIDE, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      DOUBLE PRECISION   TEMP
c     .. Parameters ..
      DOUBLE PRECISION   ONE         , ZERO
      PARAMETER        ( ONE = 1.0D+0, ZERO = 0.0D+0 )
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
c
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.
     $         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.
     $         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
     $         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'DTRSM ', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( N.EQ.0 )
     $   RETURN
c
c     And when  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
c
c     Start the operations.
c
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
c
c           Form  B := alpha*inv( A )*B.
c
            IF( UPPER )THEN
               DO 60, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 30, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   30                CONTINUE
                  END IF
                  DO 50, K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 40, I = 1, K - 1
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   40                   CONTINUE
                     END IF
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 100, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 70, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   70                CONTINUE
                  END IF
                  DO 90 K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 80, I = K + 1, M
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   80                   CONTINUE
                     END IF
   90             CONTINUE
  100          CONTINUE
            END IF
         ELSE
c
c           Form  B := alpha*inv( A' )*B.
c
            IF( UPPER )THEN
               DO 130, J = 1, N
                  DO 120, I = 1, M
                     TEMP = ALPHA*B( I, J )
                     DO 110, K = 1, I - 1
                        TEMP = TEMP - A( K, I )*B( K, J )
  110                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  120             CONTINUE
  130          CONTINUE
            ELSE
               DO 160, J = 1, N
                  DO 150, I = M, 1, -1
                     TEMP = ALPHA*B( I, J )
                     DO 140, K = I + 1, M
                        TEMP = TEMP - A( K, I )*B( K, J )
  140                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  150             CONTINUE
  160          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
c
c           Form  B := alpha*B*inv( A ).
c
            IF( UPPER )THEN
               DO 210, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 170, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  170                CONTINUE
                  END IF
                  DO 190, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 180, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  180                   CONTINUE
                     END IF
  190             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 200, I = 1, M
                        B( I, J ) = TEMP*B( I, J )
  200                CONTINUE
                  END IF
  210          CONTINUE
            ELSE
               DO 260, J = N, 1, -1
                  IF( ALPHA.NE.ONE )THEN
                     DO 220, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  220                CONTINUE
                  END IF
                  DO 240, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 230, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  230                   CONTINUE
                     END IF
  240             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 250, I = 1, M
                       B( I, J ) = TEMP*B( I, J )
  250                CONTINUE
                  END IF
  260          CONTINUE
            END IF
         ELSE
c
c           Form  B := alpha*B*inv( A' ).
c
            IF( UPPER )THEN
               DO 310, K = N, 1, -1
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 270, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  270                CONTINUE
                  END IF
                  DO 290, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 280, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  280                   CONTINUE
                     END IF
  290             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 300, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  300                CONTINUE
                  END IF
  310          CONTINUE
            ELSE
               DO 360, K = 1, N
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 320, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  320                CONTINUE
                  END IF
                  DO 340, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 330, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  330                   CONTINUE
                     END IF
  340             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 350, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  350                CONTINUE
                  END IF
  360          CONTINUE
            END IF
         END IF
      END IF
c
      RETURN
c
c     End of DTRSM .
c
      END
      SUBROUTINE SGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
c     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      REAL               ALPHA, BETA
c     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), C( LDC, * )
c     ..
c
c  Purpose
c  =======
c
c  SGEMM  performs one of the matrix-matrix operations
c
c     C := alpha*op( A )*op( B ) + beta*C,
c
c  where  op( X ) is one of
c
c     op( X ) = X   or   op( X ) = X',
c
c  alpha and beta are scalars, and A, B and C are matrices, with op( A )
c  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
c
c  Parameters
c  ==========
c
c  TRANSA - CHARACTER*1.
c           On entry, TRANSA specifies the form of op( A ) to be used in
c           the matrix multiplication as follows:
c
c              TRANSA = 'N' or 'n',  op( A ) = A.
c
c              TRANSA = 'T' or 't',  op( A ) = A'.
c
c              TRANSA = 'C' or 'c',  op( A ) = A'.
c
c           Unchanged on exit.
c
c  TRANSB - CHARACTER*1.
c           On entry, TRANSB specifies the form of op( B ) to be used in
c           the matrix multiplication as follows:
c
c              TRANSB = 'N' or 'n',  op( B ) = B.
c
c              TRANSB = 'T' or 't',  op( B ) = B'.
c
c              TRANSB = 'C' or 'c',  op( B ) = B'.
c
c           Unchanged on exit.
c
c  M      - INTEGER.
c           On entry,  M  specifies  the number  of rows  of the  matrix
c           op( A )  and of the  matrix  C.  M  must  be at least  zero.
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry,  N  specifies the number  of columns of the matrix
c           op( B ) and the number of columns of the matrix C. N must be
c           at least zero.
c           Unchanged on exit.
c
c  K      - INTEGER.
c           On entry,  K  specifies  the number of columns of the matrix
c           op( A ) and the number of rows of the matrix op( B ). K must
c           be at least  zero.
c           Unchanged on exit.
c
c  ALPHA  - REAL            .
c           On entry, ALPHA specifies the scalar alpha.
c           Unchanged on exit.
c
c  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
c           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
c           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
c           part of the array  A  must contain the matrix  A,  otherwise
c           the leading  k by m  part of the array  A  must contain  the
c           matrix A.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
c           LDA must be at least  max( 1, m ), otherwise  LDA must be at
c           least  max( 1, k ).
c           Unchanged on exit.
c
c  B      - REAL             array of DIMENSION ( LDB, kb ), where kb is
c           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
c           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
c           part of the array  B  must contain the matrix  B,  otherwise
c           the leading  n by k  part of the array  B  must contain  the
c           matrix B.
c           Unchanged on exit.
c
c  LDB    - INTEGER.
c           On entry, LDB specifies the first dimension of B as declared
c           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
c           LDB must be at least  max( 1, k ), otherwise  LDB must be at
c           least  max( 1, n ).
c           Unchanged on exit.
c
c  BETA   - REAL            .
c           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
c           supplied as zero then C need not be set on input.
c           Unchanged on exit.
c
c  C      - REAL             array of DIMENSION ( LDC, n ).
c           Before entry, the leading  m by n  part of the array  C must
c           contain the matrix  C,  except when  beta  is zero, in which
c           case C need not be set on entry.
c           On exit, the array  C  is overwritten by the  m by n  matrix
c           ( alpha*op( A )*op( B ) + beta*C ).
c
c  LDC    - INTEGER.
c           On entry, LDC specifies the first dimension of C as declared
c           in  the  calling  (sub)  program.   LDC  must  be  at  least
c           max( 1, m ).
c           Unchanged on exit.
c
c
c  Level 3 Blas routine.
c
c  -- Written on 8-February-1989.
c     Jack Dongarra, Argonne National Laboratory.
c     Iain Duff, AERE Harwell.
c     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c     Sven Hammarling, Numerical Algorithms Group Ltd.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     .. Local Scalars ..
      LOGICAL            NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      REAL               TEMP
c     .. Parameters ..
      REAL               ONE         , ZERO
      PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
c     ..
c     .. Executable Statements ..
c
c     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
c     transposed and set  NROWA, NCOLA and  NROWB  as the number of rows
c     and  columns of  A  and the  number of  rows  of  B  respectively.
c
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
c
c     Test the input parameters.
c
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.
     $         ( .NOT.LSAME( TRANSB, 'C' ) ).AND.
     $         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'SGEMM ', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
c
c     And if  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
c
c     Start the operations.
c
      IF( NOTB )THEN
         IF( NOTA )THEN
c
c           Form  C := alpha*A*B + beta*C.
c
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE
c
c           Form  C := alpha*A'*B + beta*C
c
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         END IF
      ELSE
         IF( NOTA )THEN
c
c           Form  C := alpha*A*B' + beta*C
c
            DO 170, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 130, I = 1, M
                     C( I, J ) = ZERO
  130             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 140, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  140             CONTINUE
               END IF
               DO 160, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 150, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  150                CONTINUE
                  END IF
  160          CONTINUE
  170       CONTINUE
         ELSE
c
c           Form  C := alpha*A'*B' + beta*C
c
            DO 200, J = 1, N
               DO 190, I = 1, M
                  TEMP = ZERO
                  DO 180, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  180             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  190          CONTINUE
  200       CONTINUE
         END IF
      END IF
c
      RETURN
c
c     End of SGEMM .
c
      END
      SUBROUTINE SSYMM ( SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
c     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO
      INTEGER            M, N, LDA, LDB, LDC
      REAL               ALPHA, BETA
c     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), C( LDC, * )
c     ..
c
c  Purpose
c  =======
c
c  SSYMM  performs one of the matrix-matrix operations
c
c     C := alpha*A*B + beta*C,
c
c  or
c
c     C := alpha*B*A + beta*C,
c
c  where alpha and beta are scalars,  A is a symmetric matrix and  B and
c  C are  m by n matrices.
c
c  Parameters
c  ==========
c
c  SIDE   - CHARACTER*1.
c           On entry,  SIDE  specifies whether  the  symmetric matrix  A
c           appears on the  left or right  in the  operation as follows:
c
c              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,
c
c              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
c
c           Unchanged on exit.
c
c  UPLO   - CHARACTER*1.
c           On  entry,   UPLO  specifies  whether  the  upper  or  lower
c           triangular  part  of  the  symmetric  matrix   A  is  to  be
c           referenced as follows:
c
c              UPLO = 'U' or 'u'   Only the upper triangular part of the
c                                  symmetric matrix is to be referenced.
c
c              UPLO = 'L' or 'l'   Only the lower triangular part of the
c                                  symmetric matrix is to be referenced.
c
c           Unchanged on exit.
c
c  M      - INTEGER.
c           On entry,  M  specifies the number of rows of the matrix  C.
c           M  must be at least zero.
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry, N specifies the number of columns of the matrix C.
c           N  must be at least zero.
c           Unchanged on exit.
c
c  ALPHA  - REAL            .
c           On entry, ALPHA specifies the scalar alpha.
c           Unchanged on exit.
c
c  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
c           m  when  SIDE = 'L' or 'l'  and is  n otherwise.
c           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of
c           the array  A  must contain the  symmetric matrix,  such that
c           when  UPLO = 'U' or 'u', the leading m by m upper triangular
c           part of the array  A  must contain the upper triangular part
c           of the  symmetric matrix and the  strictly  lower triangular
c           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
c           the leading  m by m  lower triangular part  of the  array  A
c           must  contain  the  lower triangular part  of the  symmetric
c           matrix and the  strictly upper triangular part of  A  is not
c           referenced.
c           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of
c           the array  A  must contain the  symmetric matrix,  such that
c           when  UPLO = 'U' or 'u', the leading n by n upper triangular
c           part of the array  A  must contain the upper triangular part
c           of the  symmetric matrix and the  strictly  lower triangular
c           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
c           the leading  n by n  lower triangular part  of the  array  A
c           must  contain  the  lower triangular part  of the  symmetric
c           matrix and the  strictly upper triangular part of  A  is not
c           referenced.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
c           LDA must be at least  max( 1, m ), otherwise  LDA must be at
c           least  max( 1, n ).
c           Unchanged on exit.
c
c  B      - REAL             array of DIMENSION ( LDB, n ).
c           Before entry, the leading  m by n part of the array  B  must
c           contain the matrix B.
c           Unchanged on exit.
c
c  LDB    - INTEGER.
c           On entry, LDB specifies the first dimension of B as declared
c           in  the  calling  (sub)  program.   LDB  must  be  at  least
c           max( 1, m ).
c           Unchanged on exit.
c
c  BETA   - REAL            .
c           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
c           supplied as zero then C need not be set on input.
c           Unchanged on exit.
c
c  C      - REAL             array of DIMENSION ( LDC, n ).
c           Before entry, the leading  m by n  part of the array  C must
c           contain the matrix  C,  except when  beta  is zero, in which
c           case C need not be set on entry.
c           On exit, the array  C  is overwritten by the  m by n updated
c           matrix.
c
c  LDC    - INTEGER.
c           On entry, LDC specifies the first dimension of C as declared
c           in  the  calling  (sub)  program.   LDC  must  be  at  least
c           max( 1, m ).
c           Unchanged on exit.
c
c
c  Level 3 Blas routine.
c
c  -- Written on 8-February-1989.
c     Jack Dongarra, Argonne National Laboratory.
c     Iain Duff, AERE Harwell.
c     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c     Sven Hammarling, Numerical Algorithms Group Ltd.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, INFO, J, K, NROWA
      REAL               TEMP1, TEMP2
c     .. Parameters ..
      REAL               ONE         , ZERO
      PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
c     ..
c     .. Executable Statements ..
c
c     Set NROWA as the number of rows of A.
c
      IF( LSAME( SIDE, 'L' ) )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      UPPER = LSAME( UPLO, 'U' )
c
c     Test the input parameters.
c
      INFO = 0
      IF(      ( .NOT.LSAME( SIDE, 'L' ) ).AND.
     $         ( .NOT.LSAME( SIDE, 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER              ).AND.
     $         ( .NOT.LSAME( UPLO, 'L' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 7
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 9
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 12
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'SSYMM ', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
c
c     And when  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
c
c     Start the operations.
c
      IF( LSAME( SIDE, 'L' ) )THEN
c
c        Form  C := alpha*A*B + beta*C.
c
         IF( UPPER )THEN
            DO 70, J = 1, N
               DO 60, I = 1, M
                  TEMP1 = ALPHA*B( I, J )
                  TEMP2 = ZERO
                  DO 50, K = 1, I - 1
                     C( K, J ) = C( K, J ) + TEMP1    *A( K, I )
                     TEMP2     = TEMP2     + B( K, J )*A( K, I )
   50             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = TEMP1*A( I, I ) + ALPHA*TEMP2
                  ELSE
                     C( I, J ) = BETA *C( I, J ) +
     $                           TEMP1*A( I, I ) + ALPHA*TEMP2
                  END IF
   60          CONTINUE
   70       CONTINUE
         ELSE
            DO 100, J = 1, N
               DO 90, I = M, 1, -1
                  TEMP1 = ALPHA*B( I, J )
                  TEMP2 = ZERO
                  DO 80, K = I + 1, M
                     C( K, J ) = C( K, J ) + TEMP1    *A( K, I )
                     TEMP2     = TEMP2     + B( K, J )*A( K, I )
   80             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = TEMP1*A( I, I ) + ALPHA*TEMP2
                  ELSE
                     C( I, J ) = BETA *C( I, J ) +
     $                           TEMP1*A( I, I ) + ALPHA*TEMP2
                  END IF
   90          CONTINUE
  100       CONTINUE
         END IF
      ELSE
c
c        Form  C := alpha*B*A + beta*C.
c
         DO 170, J = 1, N
            TEMP1 = ALPHA*A( J, J )
            IF( BETA.EQ.ZERO )THEN
               DO 110, I = 1, M
                  C( I, J ) = TEMP1*B( I, J )
  110          CONTINUE
            ELSE
               DO 120, I = 1, M
                  C( I, J ) = BETA*C( I, J ) + TEMP1*B( I, J )
  120          CONTINUE
            END IF
            DO 140, K = 1, J - 1
               IF( UPPER )THEN
                  TEMP1 = ALPHA*A( K, J )
               ELSE
                  TEMP1 = ALPHA*A( J, K )
               END IF
               DO 130, I = 1, M
                  C( I, J ) = C( I, J ) + TEMP1*B( I, K )
  130          CONTINUE
  140       CONTINUE
            DO 160, K = J + 1, N
               IF( UPPER )THEN
                  TEMP1 = ALPHA*A( J, K )
               ELSE
                  TEMP1 = ALPHA*A( K, J )
               END IF
               DO 150, I = 1, M
                  C( I, J ) = C( I, J ) + TEMP1*B( I, K )
  150          CONTINUE
  160       CONTINUE
  170    CONTINUE
      END IF
c
      RETURN
c
c     End of SSYMM .
c
      END
      SUBROUTINE SSYR2K( UPLO, TRANS, N, K, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
c     .. Scalar Arguments ..
      CHARACTER*1        UPLO, TRANS
      INTEGER            N, K, LDA, LDB, LDC
      REAL               ALPHA, BETA
c     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), C( LDC, * )
c     ..
c
c  Purpose
c  =======
c
c  SSYR2K  performs one of the symmetric rank 2k operations
c
c     C := alpha*A*B' + alpha*B*A' + beta*C,
c
c  or
c
c     C := alpha*A'*B + alpha*B'*A + beta*C,
c
c  where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
c  and  A and B  are  n by k  matrices  in the  first  case  and  k by n
c  matrices in the second case.
c
c  Parameters
c  ==========
c
c  UPLO   - CHARACTER*1.
c           On  entry,   UPLO  specifies  whether  the  upper  or  lower
c           triangular  part  of the  array  C  is to be  referenced  as
c           follows:
c
c              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
c                                  is to be referenced.
c
c              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
c                                  is to be referenced.
c
c           Unchanged on exit.
c
c  TRANS  - CHARACTER*1.
c           On entry,  TRANS  specifies the operation to be performed as
c           follows:
c
c              TRANS = 'N' or 'n'   C := alpha*A*B' + alpha*B*A' +
c                                        beta*C.
c
c              TRANS = 'T' or 't'   C := alpha*A'*B + alpha*B'*A +
c                                        beta*C.
c
c              TRANS = 'C' or 'c'   C := alpha*A'*B + alpha*B'*A +
c                                        beta*C.
c
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry,  N specifies the order of the matrix C.  N must be
c           at least zero.
c           Unchanged on exit.
c
c  K      - INTEGER.
c           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
c           of  columns  of the  matrices  A and B,  and on  entry  with
c           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
c           of rows of the matrices  A and B.  K must be at least  zero.
c           Unchanged on exit.
c
c  ALPHA  - REAL            .
c           On entry, ALPHA specifies the scalar alpha.
c           Unchanged on exit.
c
c  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
c           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
c           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
c           part of the array  A  must contain the matrix  A,  otherwise
c           the leading  k by n  part of the array  A  must contain  the
c           matrix A.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
c           then  LDA must be at least  max( 1, n ), otherwise  LDA must
c           be at least  max( 1, k ).
c           Unchanged on exit.
c
c  B      - REAL             array of DIMENSION ( LDB, kb ), where kb is
c           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
c           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
c           part of the array  B  must contain the matrix  B,  otherwise
c           the leading  k by n  part of the array  B  must contain  the
c           matrix B.
c           Unchanged on exit.
c
c  LDB    - INTEGER.
c           On entry, LDB specifies the first dimension of B as declared
c           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
c           then  LDB must be at least  max( 1, n ), otherwise  LDB must
c           be at least  max( 1, k ).
c           Unchanged on exit.
c
c  BETA   - REAL            .
c           On entry, BETA specifies the scalar beta.
c           Unchanged on exit.
c
c  C      - REAL             array of DIMENSION ( LDC, n ).
c           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
c           upper triangular part of the array C must contain the upper
c           triangular part  of the  symmetric matrix  and the strictly
c           lower triangular part of C is not referenced.  On exit, the
c           upper triangular part of the array  C is overwritten by the
c           upper triangular part of the updated matrix.
c           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
c           lower triangular part of the array C must contain the lower
c           triangular part  of the  symmetric matrix  and the strictly
c           upper triangular part of C is not referenced.  On exit, the
c           lower triangular part of the array  C is overwritten by the
c           lower triangular part of the updated matrix.
c
c  LDC    - INTEGER.
c           On entry, LDC specifies the first dimension of C as declared
c           in  the  calling  (sub)  program.   LDC  must  be  at  least
c           max( 1, n ).
c           Unchanged on exit.
c
c
c  Level 3 Blas routine.
c
c
c  -- Written on 8-February-1989.
c     Jack Dongarra, Argonne National Laboratory.
c     Iain Duff, AERE Harwell.
c     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c     Sven Hammarling, Numerical Algorithms Group Ltd.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, INFO, J, L, NROWA
      REAL               TEMP1, TEMP2
c     .. Parameters ..
      REAL               ONE         , ZERO
      PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      IF( LSAME( TRANS, 'N' ) )THEN
         NROWA = N
      ELSE
         NROWA = K
      END IF
      UPPER = LSAME( UPLO, 'U' )
c
      INFO = 0
      IF(      ( .NOT.UPPER               ).AND.
     $         ( .NOT.LSAME( UPLO , 'L' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.LSAME( TRANS, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANS, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANS, 'C' ) )      )THEN
         INFO = 2
      ELSE IF( N  .LT.0               )THEN
         INFO = 3
      ELSE IF( K  .LT.0               )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 7
      ELSE IF( LDB.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDC.LT.MAX( 1, N     ) )THEN
         INFO = 12
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'SSYR2K', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( ( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
c
c     And when  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO )THEN
         IF( UPPER )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 20, J = 1, N
                  DO 10, I = 1, J
                     C( I, J ) = ZERO
   10             CONTINUE
   20          CONTINUE
            ELSE
               DO 40, J = 1, N
                  DO 30, I = 1, J
                     C( I, J ) = BETA*C( I, J )
   30             CONTINUE
   40          CONTINUE
            END IF
         ELSE
            IF( BETA.EQ.ZERO )THEN
               DO 60, J = 1, N
                  DO 50, I = J, N
                     C( I, J ) = ZERO
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 80, J = 1, N
                  DO 70, I = J, N
                     C( I, J ) = BETA*C( I, J )
   70             CONTINUE
   80          CONTINUE
            END IF
         END IF
         RETURN
      END IF
c
c     Start the operations.
c
      IF( LSAME( TRANS, 'N' ) )THEN
c
c        Form  C := alpha*A*B' + alpha*B*A' + C.
c
         IF( UPPER )THEN
            DO 130, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 90, I = 1, J
                     C( I, J ) = ZERO
   90             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 100, I = 1, J
                     C( I, J ) = BETA*C( I, J )
  100             CONTINUE
               END IF
               DO 120, L = 1, K
                  IF( ( A( J, L ).NE.ZERO ).OR.
     $                ( B( J, L ).NE.ZERO )     )THEN
                     TEMP1 = ALPHA*B( J, L )
                     TEMP2 = ALPHA*A( J, L )
                     DO 110, I = 1, J
                        C( I, J ) = C( I, J ) +
     $                              A( I, L )*TEMP1 + B( I, L )*TEMP2
  110                CONTINUE
                  END IF
  120          CONTINUE
  130       CONTINUE
         ELSE
            DO 180, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 140, I = J, N
                     C( I, J ) = ZERO
  140             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 150, I = J, N
                     C( I, J ) = BETA*C( I, J )
  150             CONTINUE
               END IF
               DO 170, L = 1, K
                  IF( ( A( J, L ).NE.ZERO ).OR.
     $                ( B( J, L ).NE.ZERO )     )THEN
                     TEMP1 = ALPHA*B( J, L )
                     TEMP2 = ALPHA*A( J, L )
                     DO 160, I = J, N
                        C( I, J ) = C( I, J ) +
     $                              A( I, L )*TEMP1 + B( I, L )*TEMP2
  160                CONTINUE
                  END IF
  170          CONTINUE
  180       CONTINUE
         END IF
      ELSE
c
c        Form  C := alpha*A'*B + alpha*B'*A + C.
c
         IF( UPPER )THEN
            DO 210, J = 1, N
               DO 200, I = 1, J
                  TEMP1 = ZERO
                  TEMP2 = ZERO
                  DO 190, L = 1, K
                     TEMP1 = TEMP1 + A( L, I )*B( L, J )
                     TEMP2 = TEMP2 + B( L, I )*A( L, J )
  190             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP1 + ALPHA*TEMP2
                  ELSE
                     C( I, J ) = BETA *C( I, J ) +
     $                           ALPHA*TEMP1 + ALPHA*TEMP2
                  END IF
  200          CONTINUE
  210       CONTINUE
         ELSE
            DO 240, J = 1, N
               DO 230, I = J, N
                  TEMP1 = ZERO
                  TEMP2 = ZERO
                  DO 220, L = 1, K
                     TEMP1 = TEMP1 + A( L, I )*B( L, J )
                     TEMP2 = TEMP2 + B( L, I )*A( L, J )
  220             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP1 + ALPHA*TEMP2
                  ELSE
                     C( I, J ) = BETA *C( I, J ) +
     $                           ALPHA*TEMP1 + ALPHA*TEMP2
                  END IF
  230          CONTINUE
  240       CONTINUE
         END IF
      END IF
c
      RETURN
c
c     End of SSYR2K.
c
      END
      SUBROUTINE SSYRK ( UPLO, TRANS, N, K, ALPHA, A, LDA,
     $                   BETA, C, LDC )
c     .. Scalar Arguments ..
      CHARACTER*1        UPLO, TRANS
      INTEGER            N, K, LDA, LDC
      REAL               ALPHA, BETA
c     .. Array Arguments ..
      REAL               A( LDA, * ), C( LDC, * )
c     ..
c
c  Purpose
c  =======
c
c  SSYRK  performs one of the symmetric rank k operations
c
c     C := alpha*A*A' + beta*C,
c
c  or
c
c     C := alpha*A'*A + beta*C,
c
c  where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
c  and  A  is an  n by k  matrix in the first case and a  k by n  matrix
c  in the second case.
c
c  Parameters
c  ==========
c
c  UPLO   - CHARACTER*1.
c           On  entry,   UPLO  specifies  whether  the  upper  or  lower
c           triangular  part  of the  array  C  is to be  referenced  as
c           follows:
c
c              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
c                                  is to be referenced.
c
c              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
c                                  is to be referenced.
c
c           Unchanged on exit.
c
c  TRANS  - CHARACTER*1.
c           On entry,  TRANS  specifies the operation to be performed as
c           follows:
c
c              TRANS = 'N' or 'n'   C := alpha*A*A' + beta*C.
c
c              TRANS = 'T' or 't'   C := alpha*A'*A + beta*C.
c
c              TRANS = 'C' or 'c'   C := alpha*A'*A + beta*C.
c
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry,  N specifies the order of the matrix C.  N must be
c           at least zero.
c           Unchanged on exit.
c
c  K      - INTEGER.
c           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
c           of  columns   of  the   matrix   A,   and  on   entry   with
c           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
c           of rows of the matrix  A.  K must be at least zero.
c           Unchanged on exit.
c
c  ALPHA  - REAL            .
c           On entry, ALPHA specifies the scalar alpha.
c           Unchanged on exit.
c
c  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
c           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
c           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
c           part of the array  A  must contain the matrix  A,  otherwise
c           the leading  k by n  part of the array  A  must contain  the
c           matrix A.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
c           then  LDA must be at least  max( 1, n ), otherwise  LDA must
c           be at least  max( 1, k ).
c           Unchanged on exit.
c
c  BETA   - REAL            .
c           On entry, BETA specifies the scalar beta.
c           Unchanged on exit.
c
c  C      - REAL             array of DIMENSION ( LDC, n ).
c           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
c           upper triangular part of the array C must contain the upper
c           triangular part  of the  symmetric matrix  and the strictly
c           lower triangular part of C is not referenced.  On exit, the
c           upper triangular part of the array  C is overwritten by the
c           upper triangular part of the updated matrix.
c           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
c           lower triangular part of the array C must contain the lower
c           triangular part  of the  symmetric matrix  and the strictly
c           upper triangular part of C is not referenced.  On exit, the
c           lower triangular part of the array  C is overwritten by the
c           lower triangular part of the updated matrix.
c
c  LDC    - INTEGER.
c           On entry, LDC specifies the first dimension of C as declared
c           in  the  calling  (sub)  program.   LDC  must  be  at  least
c           max( 1, n ).
c           Unchanged on exit.
c
c
c  Level 3 Blas routine.
c
c  -- Written on 8-February-1989.
c     Jack Dongarra, Argonne National Laboratory.
c     Iain Duff, AERE Harwell.
c     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c     Sven Hammarling, Numerical Algorithms Group Ltd.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, INFO, J, L, NROWA
      REAL               TEMP
c     .. Parameters ..
      REAL               ONE ,         ZERO
      PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      IF( LSAME( TRANS, 'N' ) )THEN
         NROWA = N
      ELSE
         NROWA = K
      END IF
      UPPER = LSAME( UPLO, 'U' )
c
      INFO = 0
      IF(      ( .NOT.UPPER               ).AND.
     $         ( .NOT.LSAME( UPLO , 'L' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.LSAME( TRANS, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANS, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANS, 'C' ) )      )THEN
         INFO = 2
      ELSE IF( N  .LT.0               )THEN
         INFO = 3
      ELSE IF( K  .LT.0               )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 7
      ELSE IF( LDC.LT.MAX( 1, N     ) )THEN
         INFO = 10
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'SSYRK ', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( ( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
c
c     And when  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO )THEN
         IF( UPPER )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 20, J = 1, N
                  DO 10, I = 1, J
                     C( I, J ) = ZERO
   10             CONTINUE
   20          CONTINUE
            ELSE
               DO 40, J = 1, N
                  DO 30, I = 1, J
                     C( I, J ) = BETA*C( I, J )
   30             CONTINUE
   40          CONTINUE
            END IF
         ELSE
            IF( BETA.EQ.ZERO )THEN
               DO 60, J = 1, N
                  DO 50, I = J, N
                     C( I, J ) = ZERO
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 80, J = 1, N
                  DO 70, I = J, N
                     C( I, J ) = BETA*C( I, J )
   70             CONTINUE
   80          CONTINUE
            END IF
         END IF
         RETURN
      END IF
c
c     Start the operations.
c
      IF( LSAME( TRANS, 'N' ) )THEN
c
c        Form  C := alpha*A*A' + beta*C.
c
         IF( UPPER )THEN
            DO 130, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 90, I = 1, J
                     C( I, J ) = ZERO
   90             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 100, I = 1, J
                     C( I, J ) = BETA*C( I, J )
  100             CONTINUE
               END IF
               DO 120, L = 1, K
                  IF( A( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*A( J, L )
                     DO 110, I = 1, J
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  110                CONTINUE
                  END IF
  120          CONTINUE
  130       CONTINUE
         ELSE
            DO 180, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 140, I = J, N
                     C( I, J ) = ZERO
  140             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 150, I = J, N
                     C( I, J ) = BETA*C( I, J )
  150             CONTINUE
               END IF
               DO 170, L = 1, K
                  IF( A( J, L ).NE.ZERO )THEN
                     TEMP      = ALPHA*A( J, L )
                     DO 160, I = J, N
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  160                CONTINUE
                  END IF
  170          CONTINUE
  180       CONTINUE
         END IF
      ELSE
c
c        Form  C := alpha*A'*A + beta*C.
c
         IF( UPPER )THEN
            DO 210, J = 1, N
               DO 200, I = 1, J
                  TEMP = ZERO
                  DO 190, L = 1, K
                     TEMP = TEMP + A( L, I )*A( L, J )
  190             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  200          CONTINUE
  210       CONTINUE
         ELSE
            DO 240, J = 1, N
               DO 230, I = J, N
                  TEMP = ZERO
                  DO 220, L = 1, K
                     TEMP = TEMP + A( L, I )*A( L, J )
  220             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  230          CONTINUE
  240       CONTINUE
         END IF
      END IF
c
      RETURN
c
c     End of SSYRK .
c
      END
      SUBROUTINE STRMM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   B, LDB )
c     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      REAL               ALPHA
c     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * )
c     ..
c
c  Purpose
c  =======
c
c  STRMM  performs one of the matrix-matrix operations
c
c     B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
c
c  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
c  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
c
c     op( A ) = A   or   op( A ) = A'.
c
c  Parameters
c  ==========
c
c  SIDE   - CHARACTER*1.
c           On entry,  SIDE specifies whether  op( A ) multiplies B from
c           the left or right as follows:
c
c              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
c
c              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
c
c           Unchanged on exit.
c
c  UPLO   - CHARACTER*1.
c           On entry, UPLO specifies whether the matrix A is an upper or
c           lower triangular matrix as follows:
c
c              UPLO = 'U' or 'u'   A is an upper triangular matrix.
c
c              UPLO = 'L' or 'l'   A is a lower triangular matrix.
c
c           Unchanged on exit.
c
c  TRANSA - CHARACTER*1.
c           On entry, TRANSA specifies the form of op( A ) to be used in
c           the matrix multiplication as follows:
c
c              TRANSA = 'N' or 'n'   op( A ) = A.
c
c              TRANSA = 'T' or 't'   op( A ) = A'.
c
c              TRANSA = 'C' or 'c'   op( A ) = A'.
c
c           Unchanged on exit.
c
c  DIAG   - CHARACTER*1.
c           On entry, DIAG specifies whether or not A is unit triangular
c           as follows:
c
c              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
c
c              DIAG = 'N' or 'n'   A is not assumed to be unit
c                                  triangular.
c
c           Unchanged on exit.
c
c  M      - INTEGER.
c           On entry, M specifies the number of rows of B. M must be at
c           least zero.
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry, N specifies the number of columns of B.  N must be
c           at least zero.
c           Unchanged on exit.
c
c  ALPHA  - REAL            .
c           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
c           zero then  A is not referenced and  B need not be set before
c           entry.
c           Unchanged on exit.
c
c  A      - REAL             array of DIMENSION ( LDA, k ), where k is m
c           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
c           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
c           upper triangular part of the array  A must contain the upper
c           triangular matrix  and the strictly lower triangular part of
c           A is not referenced.
c           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
c           lower triangular part of the array  A must contain the lower
c           triangular matrix  and the strictly upper triangular part of
c           A is not referenced.
c           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
c           A  are not referenced either,  but are assumed to be  unity.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
c           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
c           then LDA must be at least max( 1, n ).
c           Unchanged on exit.
c
c  B      - REAL             array of DIMENSION ( LDB, n ).
c           Before entry,  the leading  m by n part of the array  B must
c           contain the matrix  B,  and  on exit  is overwritten  by the
c           transformed matrix.
c
c  LDB    - INTEGER.
c           On entry, LDB specifies the first dimension of B as declared
c           in  the  calling  (sub)  program.   LDB  must  be  at  least
c           max( 1, m ).
c           Unchanged on exit.
c
c
c  Level 3 Blas routine.
c
c  -- Written on 8-February-1989.
c     Jack Dongarra, Argonne National Laboratory.
c     Iain Duff, AERE Harwell.
c     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c     Sven Hammarling, Numerical Algorithms Group Ltd.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     .. Local Scalars ..
      LOGICAL            LSIDE, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      REAL               TEMP
c     .. Parameters ..
      REAL               ONE         , ZERO
      PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
c
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.
     $         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.
     $         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
     $         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'STRMM ', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( N.EQ.0 )
     $   RETURN
c
c     And when  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
c
c     Start the operations.
c
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
c
c           Form  B := alpha*A*B.
c
            IF( UPPER )THEN
               DO 50, J = 1, N
                  DO 40, K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*B( K, J )
                        DO 30, I = 1, K - 1
                           B( I, J ) = B( I, J ) + TEMP*A( I, K )
   30                   CONTINUE
                        IF( NOUNIT )
     $                     TEMP = TEMP*A( K, K )
                        B( K, J ) = TEMP
                     END IF
   40             CONTINUE
   50          CONTINUE
            ELSE
               DO 80, J = 1, N
                  DO 70 K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        TEMP      = ALPHA*B( K, J )
                        B( K, J ) = TEMP
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )*A( K, K )
                        DO 60, I = K + 1, M
                           B( I, J ) = B( I, J ) + TEMP*A( I, K )
   60                   CONTINUE
                     END IF
   70             CONTINUE
   80          CONTINUE
            END IF
         ELSE
c
c           Form  B := alpha*A'*B.
c
            IF( UPPER )THEN
               DO 110, J = 1, N
                  DO 100, I = M, 1, -1
                     TEMP = B( I, J )
                     IF( NOUNIT )
     $                  TEMP = TEMP*A( I, I )
                     DO 90, K = 1, I - 1
                        TEMP = TEMP + A( K, I )*B( K, J )
   90                CONTINUE
                     B( I, J ) = ALPHA*TEMP
  100             CONTINUE
  110          CONTINUE
            ELSE
               DO 140, J = 1, N
                  DO 130, I = 1, M
                     TEMP = B( I, J )
                     IF( NOUNIT )
     $                  TEMP = TEMP*A( I, I )
                     DO 120, K = I + 1, M
                        TEMP = TEMP + A( K, I )*B( K, J )
  120                CONTINUE
                     B( I, J ) = ALPHA*TEMP
  130             CONTINUE
  140          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
c
c           Form  B := alpha*B*A.
c
            IF( UPPER )THEN
               DO 180, J = N, 1, -1
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 150, I = 1, M
                     B( I, J ) = TEMP*B( I, J )
  150             CONTINUE
                  DO 170, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*A( K, J )
                        DO 160, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  160                   CONTINUE
                     END IF
  170             CONTINUE
  180          CONTINUE
            ELSE
               DO 220, J = 1, N
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 190, I = 1, M
                     B( I, J ) = TEMP*B( I, J )
  190             CONTINUE
                  DO 210, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*A( K, J )
                        DO 200, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  200                   CONTINUE
                     END IF
  210             CONTINUE
  220          CONTINUE
            END IF
         ELSE
c
c           Form  B := alpha*B*A'.
c
            IF( UPPER )THEN
               DO 260, K = 1, N
                  DO 240, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = ALPHA*A( J, K )
                        DO 230, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  230                   CONTINUE
                     END IF
  240             CONTINUE
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( K, K )
                  IF( TEMP.NE.ONE )THEN
                     DO 250, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  250                CONTINUE
                  END IF
  260          CONTINUE
            ELSE
               DO 300, K = N, 1, -1
                  DO 280, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = ALPHA*A( J, K )
                        DO 270, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  270                   CONTINUE
                     END IF
  280             CONTINUE
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( K, K )
                  IF( TEMP.NE.ONE )THEN
                     DO 290, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  290                CONTINUE
                  END IF
  300          CONTINUE
            END IF
         END IF
      END IF
c
      RETURN
c
c     End of STRMM .
c
      END
      SUBROUTINE STRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   B, LDB )
c     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      REAL               ALPHA
c     .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * )
c     ..
c
c  Purpose
c  =======
c
c  STRSM  solves one of the matrix equations
c
c     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
c
c  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
c  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
c
c     op( A ) = A   or   op( A ) = A'.
c
c  The matrix X is overwritten on B.
c
c  Parameters
c  ==========
c
c  SIDE   - CHARACTER*1.
c           On entry, SIDE specifies whether op( A ) appears on the left
c           or right of X as follows:
c
c              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
c
c              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
c
c           Unchanged on exit.
c
c  UPLO   - CHARACTER*1.
c           On entry, UPLO specifies whether the matrix A is an upper or
c           lower triangular matrix as follows:
c
c              UPLO = 'U' or 'u'   A is an upper triangular matrix.
c
c              UPLO = 'L' or 'l'   A is a lower triangular matrix.
c
c           Unchanged on exit.
c
c  TRANSA - CHARACTER*1.
c           On entry, TRANSA specifies the form of op( A ) to be used in
c           the matrix multiplication as follows:
c
c              TRANSA = 'N' or 'n'   op( A ) = A.
c
c              TRANSA = 'T' or 't'   op( A ) = A'.
c
c              TRANSA = 'C' or 'c'   op( A ) = A'.
c
c           Unchanged on exit.
c
c  DIAG   - CHARACTER*1.
c           On entry, DIAG specifies whether or not A is unit triangular
c           as follows:
c
c              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
c
c              DIAG = 'N' or 'n'   A is not assumed to be unit
c                                  triangular.
c
c           Unchanged on exit.
c
c  M      - INTEGER.
c           On entry, M specifies the number of rows of B. M must be at
c           least zero.
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry, N specifies the number of columns of B.  N must be
c           at least zero.
c           Unchanged on exit.
c
c  ALPHA  - REAL            .
c           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
c           zero then  A is not referenced and  B need not be set before
c           entry.
c           Unchanged on exit.
c
c  A      - REAL             array of DIMENSION ( LDA, k ), where k is m
c           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
c           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
c           upper triangular part of the array  A must contain the upper
c           triangular matrix  and the strictly lower triangular part of
c           A is not referenced.
c           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
c           lower triangular part of the array  A must contain the lower
c           triangular matrix  and the strictly upper triangular part of
c           A is not referenced.
c           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
c           A  are not referenced either,  but are assumed to be  unity.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
c           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
c           then LDA must be at least max( 1, n ).
c           Unchanged on exit.
c
c  B      - REAL             array of DIMENSION ( LDB, n ).
c           Before entry,  the leading  m by n part of the array  B must
c           contain  the  right-hand  side  matrix  B,  and  on exit  is
c           overwritten by the solution matrix  X.
c
c  LDB    - INTEGER.
c           On entry, LDB specifies the first dimension of B as declared
c           in  the  calling  (sub)  program.   LDB  must  be  at  least
c           max( 1, m ).
c           Unchanged on exit.
c
c
c  Level 3 Blas routine.
c
c
c  -- Written on 8-February-1989.
c     Jack Dongarra, Argonne National Laboratory.
c     Iain Duff, AERE Harwell.
c     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c     Sven Hammarling, Numerical Algorithms Group Ltd.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     .. Local Scalars ..
      LOGICAL            LSIDE, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      REAL               TEMP
c     .. Parameters ..
      REAL               ONE         , ZERO
      PARAMETER        ( ONE = 1.0E+0, ZERO = 0.0E+0 )
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
c
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.
     $         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.
     $         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
     $         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'STRSM ', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( N.EQ.0 )
     $   RETURN
c
c     And when  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
c
c     Start the operations.
c
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
c
c           Form  B := alpha*inv( A )*B.
c
            IF( UPPER )THEN
               DO 60, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 30, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   30                CONTINUE
                  END IF
                  DO 50, K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 40, I = 1, K - 1
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   40                   CONTINUE
                     END IF
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 100, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 70, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   70                CONTINUE
                  END IF
                  DO 90 K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 80, I = K + 1, M
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   80                   CONTINUE
                     END IF
   90             CONTINUE
  100          CONTINUE
            END IF
         ELSE
c
c           Form  B := alpha*inv( A' )*B.
c
            IF( UPPER )THEN
               DO 130, J = 1, N
                  DO 120, I = 1, M
                     TEMP = ALPHA*B( I, J )
                     DO 110, K = 1, I - 1
                        TEMP = TEMP - A( K, I )*B( K, J )
  110                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  120             CONTINUE
  130          CONTINUE
            ELSE
               DO 160, J = 1, N
                  DO 150, I = M, 1, -1
                     TEMP = ALPHA*B( I, J )
                     DO 140, K = I + 1, M
                        TEMP = TEMP - A( K, I )*B( K, J )
  140                CONTINUE
                     IF( NOUNIT )
     $                  TEMP = TEMP/A( I, I )
                     B( I, J ) = TEMP
  150             CONTINUE
  160          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
c
c           Form  B := alpha*B*inv( A ).
c
            IF( UPPER )THEN
               DO 210, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 170, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  170                CONTINUE
                  END IF
                  DO 190, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 180, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  180                   CONTINUE
                     END IF
  190             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 200, I = 1, M
                        B( I, J ) = TEMP*B( I, J )
  200                CONTINUE
                  END IF
  210          CONTINUE
            ELSE
               DO 260, J = N, 1, -1
                  IF( ALPHA.NE.ONE )THEN
                     DO 220, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  220                CONTINUE
                  END IF
                  DO 240, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 230, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  230                   CONTINUE
                     END IF
  240             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 250, I = 1, M
                       B( I, J ) = TEMP*B( I, J )
  250                CONTINUE
                  END IF
  260          CONTINUE
            END IF
         ELSE
c
c           Form  B := alpha*B*inv( A' ).
c
            IF( UPPER )THEN
               DO 310, K = N, 1, -1
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 270, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  270                CONTINUE
                  END IF
                  DO 290, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 280, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  280                   CONTINUE
                     END IF
  290             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 300, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  300                CONTINUE
                  END IF
  310          CONTINUE
            ELSE
               DO 360, K = 1, N
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( K, K )
                     DO 320, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  320                CONTINUE
                  END IF
                  DO 340, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        TEMP = A( J, K )
                        DO 330, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  330                   CONTINUE
                     END IF
  340             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 350, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  350                CONTINUE
                  END IF
  360          CONTINUE
            END IF
         END IF
      END IF
c
      RETURN
c
c     End of STRSM .
c
      END
      SUBROUTINE ZGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
c     .. Scalar Arguments ..
      CHARACTER*1        TRANSA, TRANSB
      INTEGER            M, N, K, LDA, LDB, LDC
      COMPLEX*16         ALPHA, BETA
c     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDC, * )
c     ..
c
c  Purpose
c  =======
c
c  ZGEMM  performs one of the matrix-matrix operations
c
c     C := alpha*op( A )*op( B ) + beta*C,
c
c  where  op( X ) is one of
c
c     op( X ) = X   or   op( X ) = X'   or   op( X ) = conjg( X' ),
c
c  alpha and beta are scalars, and A, B and C are matrices, with op( A )
c  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
c
c  Parameters
c  ==========
c
c  TRANSA - CHARACTER*1.
c           On entry, TRANSA specifies the form of op( A ) to be used in
c           the matrix multiplication as follows:
c
c              TRANSA = 'N' or 'n',  op( A ) = A.
c
c              TRANSA = 'T' or 't',  op( A ) = A'.
c
c              TRANSA = 'C' or 'c',  op( A ) = conjg( A' ).
c
c           Unchanged on exit.
c
c  TRANSB - CHARACTER*1.
c           On entry, TRANSB specifies the form of op( B ) to be used in
c           the matrix multiplication as follows:
c
c              TRANSB = 'N' or 'n',  op( B ) = B.
c
c              TRANSB = 'T' or 't',  op( B ) = B'.
c
c              TRANSB = 'C' or 'c',  op( B ) = conjg( B' ).
c
c           Unchanged on exit.
c
c  M      - INTEGER.
c           On entry,  M  specifies  the number  of rows  of the  matrix
c           op( A )  and of the  matrix  C.  M  must  be at least  zero.
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry,  N  specifies the number  of columns of the matrix
c           op( B ) and the number of columns of the matrix C. N must be
c           at least zero.
c           Unchanged on exit.
c
c  K      - INTEGER.
c           On entry,  K  specifies  the number of columns of the matrix
c           op( A ) and the number of rows of the matrix op( B ). K must
c           be at least  zero.
c           Unchanged on exit.
c
c  ALPHA  - COMPLEX*16      .
c           On entry, ALPHA specifies the scalar alpha.
c           Unchanged on exit.
c
c  A      - COMPLEX*16       array of DIMENSION ( LDA, ka ), where ka is
c           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
c           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
c           part of the array  A  must contain the matrix  A,  otherwise
c           the leading  k by m  part of the array  A  must contain  the
c           matrix A.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
c           LDA must be at least  max( 1, m ), otherwise  LDA must be at
c           least  max( 1, k ).
c           Unchanged on exit.
c
c  B      - COMPLEX*16       array of DIMENSION ( LDB, kb ), where kb is
c           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
c           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
c           part of the array  B  must contain the matrix  B,  otherwise
c           the leading  n by k  part of the array  B  must contain  the
c           matrix B.
c           Unchanged on exit.
c
c  LDB    - INTEGER.
c           On entry, LDB specifies the first dimension of B as declared
c           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
c           LDB must be at least  max( 1, k ), otherwise  LDB must be at
c           least  max( 1, n ).
c           Unchanged on exit.
c
c  BETA   - COMPLEX*16      .
c           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
c           supplied as zero then C need not be set on input.
c           Unchanged on exit.
c
c  C      - COMPLEX*16       array of DIMENSION ( LDC, n ).
c           Before entry, the leading  m by n  part of the array  C must
c           contain the matrix  C,  except when  beta  is zero, in which
c           case C need not be set on entry.
c           On exit, the array  C  is overwritten by the  m by n  matrix
c           ( alpha*op( A )*op( B ) + beta*C ).
c
c  LDC    - INTEGER.
c           On entry, LDC specifies the first dimension of C as declared
c           in  the  calling  (sub)  program.   LDC  must  be  at  least
c           max( 1, m ).
c           Unchanged on exit.
c
c
c  Level 3 Blas routine.
c
c  -- Written on 8-February-1989.
c     Jack Dongarra, Argonne National Laboratory.
c     Iain Duff, AERE Harwell.
c     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c     Sven Hammarling, Numerical Algorithms Group Ltd.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
c     .. Local Scalars ..
      LOGICAL            CONJA, CONJB, NOTA, NOTB
      INTEGER            I, INFO, J, L, NCOLA, NROWA, NROWB
      COMPLEX*16         TEMP
c     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER        ( ONE  = ( 1.0D+0, 0.0D+0 ) )
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
c     ..
c     .. Executable Statements ..
c
c     Set  NOTA  and  NOTB  as  true if  A  and  B  respectively are not
c     conjugated or transposed, set  CONJA and CONJB  as true if  A  and
c     B  respectively are to be  transposed but  not conjugated  and set
c     NROWA, NCOLA and  NROWB  as the number of rows and  columns  of  A
c     and the number of rows of  B  respectively.
c
      NOTA  = LSAME( TRANSA, 'N' )
      NOTB  = LSAME( TRANSB, 'N' )
      CONJA = LSAME( TRANSA, 'C' )
      CONJB = LSAME( TRANSB, 'C' )
      IF( NOTA )THEN
         NROWA = M
         NCOLA = K
      ELSE
         NROWA = K
         NCOLA = M
      END IF
      IF( NOTB )THEN
         NROWB = K
      ELSE
         NROWB = N
      END IF
c
c     Test the input parameters.
c
      INFO = 0
      IF(      ( .NOT.NOTA                 ).AND.
     $         ( .NOT.CONJA                ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.NOTB                 ).AND.
     $         ( .NOT.CONJB                ).AND.
     $         ( .NOT.LSAME( TRANSB, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( K  .LT.0               )THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 8
      ELSE IF( LDB.LT.MAX( 1, NROWB ) )THEN
         INFO = 10
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 13
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'ZGEMM ', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
c
c     And when  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
c
c     Start the operations.
c
      IF( NOTB )THEN
         IF( NOTA )THEN
c
c           Form  C := alpha*A*B + beta*C.
c
            DO 90, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 50, I = 1, M
                     C( I, J ) = ZERO
   50             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 60, I = 1, M
                     C( I, J ) = BETA*C( I, J )
   60             CONTINUE
               END IF
               DO 80, L = 1, K
                  IF( B( L, J ).NE.ZERO )THEN
                     TEMP = ALPHA*B( L, J )
                     DO 70, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
   70                CONTINUE
                  END IF
   80          CONTINUE
   90       CONTINUE
         ELSE IF( CONJA )THEN
c
c           Form  C := alpha*conjg( A' )*B + beta*C.
c
            DO 120, J = 1, N
               DO 110, I = 1, M
                  TEMP = ZERO
                  DO 100, L = 1, K
                     TEMP = TEMP + DCONJG( A( L, I ) )*B( L, J )
  100             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  110          CONTINUE
  120       CONTINUE
         ELSE
c
c           Form  C := alpha*A'*B + beta*C
c
            DO 150, J = 1, N
               DO 140, I = 1, M
                  TEMP = ZERO
                  DO 130, L = 1, K
                     TEMP = TEMP + A( L, I )*B( L, J )
  130             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  140          CONTINUE
  150       CONTINUE
         END IF
      ELSE IF( NOTA )THEN
         IF( CONJB )THEN
c
c           Form  C := alpha*A*conjg( B' ) + beta*C.
c
            DO 200, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 160, I = 1, M
                     C( I, J ) = ZERO
  160             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 170, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  170             CONTINUE
               END IF
               DO 190, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*DCONJG( B( J, L ) )
                     DO 180, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  180                CONTINUE
                  END IF
  190          CONTINUE
  200       CONTINUE
         ELSE
c
c           Form  C := alpha*A*B'          + beta*C
c
            DO 250, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 210, I = 1, M
                     C( I, J ) = ZERO
  210             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 220, I = 1, M
                     C( I, J ) = BETA*C( I, J )
  220             CONTINUE
               END IF
               DO 240, L = 1, K
                  IF( B( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*B( J, L )
                     DO 230, I = 1, M
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  230                CONTINUE
                  END IF
  240          CONTINUE
  250       CONTINUE
         END IF
      ELSE IF( CONJA )THEN
         IF( CONJB )THEN
c
c           Form  C := alpha*conjg( A' )*conjg( B' ) + beta*C.
c
            DO 280, J = 1, N
               DO 270, I = 1, M
                  TEMP = ZERO
                  DO 260, L = 1, K
                     TEMP = TEMP +
     $                      DCONJG( A( L, I ) )*DCONJG( B( J, L ) )
  260             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  270          CONTINUE
  280       CONTINUE
         ELSE
c
c           Form  C := alpha*conjg( A' )*B' + beta*C
c
            DO 310, J = 1, N
               DO 300, I = 1, M
                  TEMP = ZERO
                  DO 290, L = 1, K
                     TEMP = TEMP + DCONJG( A( L, I ) )*B( J, L )
  290             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  300          CONTINUE
  310       CONTINUE
         END IF
      ELSE
         IF( CONJB )THEN
c
c           Form  C := alpha*A'*conjg( B' ) + beta*C
c
            DO 340, J = 1, N
               DO 330, I = 1, M
                  TEMP = ZERO
                  DO 320, L = 1, K
                     TEMP = TEMP + A( L, I )*DCONJG( B( J, L ) )
  320             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  330          CONTINUE
  340       CONTINUE
         ELSE
c
c           Form  C := alpha*A'*B' + beta*C
c
            DO 370, J = 1, N
               DO 360, I = 1, M
                  TEMP = ZERO
                  DO 350, L = 1, K
                     TEMP = TEMP + A( L, I )*B( J, L )
  350             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  360          CONTINUE
  370       CONTINUE
         END IF
      END IF
c
      RETURN
c
c     End of ZGEMM .
c
      END
      SUBROUTINE ZHEMM ( SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
c     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO
      INTEGER            M, N, LDA, LDB, LDC
      COMPLEX*16         ALPHA, BETA
c     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDC, * )
c     ..
c
c  Purpose
c  =======
c
c  ZHEMM  performs one of the matrix-matrix operations
c
c     C := alpha*A*B + beta*C,
c
c  or
c
c     C := alpha*B*A + beta*C,
c
c  where alpha and beta are scalars, A is an hermitian matrix and  B and
c  C are m by n matrices.
c
c  Parameters
c  ==========
c
c  SIDE   - CHARACTER*1.
c           On entry,  SIDE  specifies whether  the  hermitian matrix  A
c           appears on the  left or right  in the  operation as follows:
c
c              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,
c
c              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
c
c           Unchanged on exit.
c
c  UPLO   - CHARACTER*1.
c           On  entry,   UPLO  specifies  whether  the  upper  or  lower
c           triangular  part  of  the  hermitian  matrix   A  is  to  be
c           referenced as follows:
c
c              UPLO = 'U' or 'u'   Only the upper triangular part of the
c                                  hermitian matrix is to be referenced.
c
c              UPLO = 'L' or 'l'   Only the lower triangular part of the
c                                  hermitian matrix is to be referenced.
c
c           Unchanged on exit.
c
c  M      - INTEGER.
c           On entry,  M  specifies the number of rows of the matrix  C.
c           M  must be at least zero.
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry, N specifies the number of columns of the matrix C.
c           N  must be at least zero.
c           Unchanged on exit.
c
c  ALPHA  - COMPLEX*16      .
c           On entry, ALPHA specifies the scalar alpha.
c           Unchanged on exit.
c
c  A      - COMPLEX*16       array of DIMENSION ( LDA, ka ), where ka is
c           m  when  SIDE = 'L' or 'l'  and is n  otherwise.
c           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of
c           the array  A  must contain the  hermitian matrix,  such that
c           when  UPLO = 'U' or 'u', the leading m by m upper triangular
c           part of the array  A  must contain the upper triangular part
c           of the  hermitian matrix and the  strictly  lower triangular
c           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
c           the leading  m by m  lower triangular part  of the  array  A
c           must  contain  the  lower triangular part  of the  hermitian
c           matrix and the  strictly upper triangular part of  A  is not
c           referenced.
c           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of
c           the array  A  must contain the  hermitian matrix,  such that
c           when  UPLO = 'U' or 'u', the leading n by n upper triangular
c           part of the array  A  must contain the upper triangular part
c           of the  hermitian matrix and the  strictly  lower triangular
c           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
c           the leading  n by n  lower triangular part  of the  array  A
c           must  contain  the  lower triangular part  of the  hermitian
c           matrix and the  strictly upper triangular part of  A  is not
c           referenced.
c           Note that the imaginary parts  of the diagonal elements need
c           not be set, they are assumed to be zero.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in the  calling (sub) program. When  SIDE = 'L' or 'l'  then
c           LDA must be at least  max( 1, m ), otherwise  LDA must be at
c           least max( 1, n ).
c           Unchanged on exit.
c
c  B      - COMPLEX*16       array of DIMENSION ( LDB, n ).
c           Before entry, the leading  m by n part of the array  B  must
c           contain the matrix B.
c           Unchanged on exit.
c
c  LDB    - INTEGER.
c           On entry, LDB specifies the first dimension of B as declared
c           in  the  calling  (sub)  program.   LDB  must  be  at  least
c           max( 1, m ).
c           Unchanged on exit.
c
c  BETA   - COMPLEX*16      .
c           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
c           supplied as zero then C need not be set on input.
c           Unchanged on exit.
c
c  C      - COMPLEX*16       array of DIMENSION ( LDC, n ).
c           Before entry, the leading  m by n  part of the array  C must
c           contain the matrix  C,  except when  beta  is zero, in which
c           case C need not be set on entry.
c           On exit, the array  C  is overwritten by the  m by n updated
c           matrix.
c
c  LDC    - INTEGER.
c           On entry, LDC specifies the first dimension of C as declared
c           in  the  calling  (sub)  program.   LDC  must  be  at  least
c           max( 1, m ).
c           Unchanged on exit.
c
c
c  Level 3 Blas routine.
c
c  -- Written on 8-February-1989.
c     Jack Dongarra, Argonne National Laboratory.
c     Iain Duff, AERE Harwell.
c     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c     Sven Hammarling, Numerical Algorithms Group Ltd.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX, DBLE
c     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, INFO, J, K, NROWA
      COMPLEX*16         TEMP1, TEMP2
c     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER        ( ONE  = ( 1.0D+0, 0.0D+0 ) )
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
c     ..
c     .. Executable Statements ..
c
c     Set NROWA as the number of rows of A.
c
      IF( LSAME( SIDE, 'L' ) )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      UPPER = LSAME( UPLO, 'U' )
c
c     Test the input parameters.
c
      INFO = 0
      IF(      ( .NOT.LSAME( SIDE, 'L' ) ).AND.
     $         ( .NOT.LSAME( SIDE, 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER              ).AND.
     $         ( .NOT.LSAME( UPLO, 'L' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 7
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 9
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 12
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'ZHEMM ', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
c
c     And when  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
c
c     Start the operations.
c
      IF( LSAME( SIDE, 'L' ) )THEN
c
c        Form  C := alpha*A*B + beta*C.
c
         IF( UPPER )THEN
            DO 70, J = 1, N
               DO 60, I = 1, M
                  TEMP1 = ALPHA*B( I, J )
                  TEMP2 = ZERO
                  DO 50, K = 1, I - 1
                     C( K, J ) = C( K, J ) + TEMP1*A( K, I )
                     TEMP2     = TEMP2     +
     $                           B( K, J )*DCONJG( A( K, I ) )
   50             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = TEMP1*DBLE( A( I, I ) ) +
     $                           ALPHA*TEMP2
                  ELSE
                     C( I, J ) = BETA *C( I, J )         +
     $                           TEMP1*DBLE( A( I, I ) ) +
     $                           ALPHA*TEMP2
                  END IF
   60          CONTINUE
   70       CONTINUE
         ELSE
            DO 100, J = 1, N
               DO 90, I = M, 1, -1
                  TEMP1 = ALPHA*B( I, J )
                  TEMP2 = ZERO
                  DO 80, K = I + 1, M
                     C( K, J ) = C( K, J ) + TEMP1*A( K, I )
                     TEMP2     = TEMP2     +
     $                           B( K, J )*DCONJG( A( K, I ) )
   80             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = TEMP1*DBLE( A( I, I ) ) +
     $                           ALPHA*TEMP2
                  ELSE
                     C( I, J ) = BETA *C( I, J )         +
     $                           TEMP1*DBLE( A( I, I ) ) +
     $                           ALPHA*TEMP2
                  END IF
   90          CONTINUE
  100       CONTINUE
         END IF
      ELSE
c
c        Form  C := alpha*B*A + beta*C.
c
         DO 170, J = 1, N
            TEMP1 = ALPHA*DBLE( A( J, J ) )
            IF( BETA.EQ.ZERO )THEN
               DO 110, I = 1, M
                  C( I, J ) = TEMP1*B( I, J )
  110          CONTINUE
            ELSE
               DO 120, I = 1, M
                  C( I, J ) = BETA*C( I, J ) + TEMP1*B( I, J )
  120          CONTINUE
            END IF
            DO 140, K = 1, J - 1
               IF( UPPER )THEN
                  TEMP1 = ALPHA*A( K, J )
               ELSE
                  TEMP1 = ALPHA*DCONJG( A( J, K ) )
               END IF
               DO 130, I = 1, M
                  C( I, J ) = C( I, J ) + TEMP1*B( I, K )
  130          CONTINUE
  140       CONTINUE
            DO 160, K = J + 1, N
               IF( UPPER )THEN
                  TEMP1 = ALPHA*DCONJG( A( J, K ) )
               ELSE
                  TEMP1 = ALPHA*A( K, J )
               END IF
               DO 150, I = 1, M
                  C( I, J ) = C( I, J ) + TEMP1*B( I, K )
  150          CONTINUE
  160       CONTINUE
  170    CONTINUE
      END IF
c
      RETURN
c
c     End of ZHEMM .
c
      END
      SUBROUTINE ZHER2K( UPLO, TRANS, N, K, ALPHA, A, LDA, B, LDB, BETA,
     $                   C, LDC )
c     .. Scalar Arguments ..
      CHARACTER          TRANS, UPLO
      INTEGER            K, LDA, LDB, LDC, N
      DOUBLE PRECISION   BETA
      COMPLEX*16         ALPHA
c     ..
c     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDC, * )
c     ..
c
c  Purpose
c  =======
c
c  ZHER2K  performs one of the hermitian rank 2k operations
c
c     C := alpha*A*conjg( B' ) + conjg( alpha )*B*conjg( A' ) + beta*C,
c
c  or
c
c     C := alpha*conjg( A' )*B + conjg( alpha )*conjg( B' )*A + beta*C,
c
c  where  alpha and beta  are scalars with  beta  real,  C is an  n by n
c  hermitian matrix and  A and B  are  n by k matrices in the first case
c  and  k by n  matrices in the second case.
c
c  Parameters
c  ==========
c
c  UPLO   - CHARACTER*1.
c           On  entry,   UPLO  specifies  whether  the  upper  or  lower
c           triangular  part  of the  array  C  is to be  referenced  as
c           follows:
c
c              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
c                                  is to be referenced.
c
c              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
c                                  is to be referenced.
c
c           Unchanged on exit.
c
c  TRANS  - CHARACTER*1.
c           On entry,  TRANS  specifies the operation to be performed as
c           follows:
c
c              TRANS = 'N' or 'n'    C := alpha*A*conjg( B' )          +
c                                         conjg( alpha )*B*conjg( A' ) +
c                                         beta*C.
c
c              TRANS = 'C' or 'c'    C := alpha*conjg( A' )*B          +
c                                         conjg( alpha )*conjg( B' )*A +
c                                         beta*C.
c
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry,  N specifies the order of the matrix C.  N must be
c           at least zero.
c           Unchanged on exit.
c
c  K      - INTEGER.
c           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
c           of  columns  of the  matrices  A and B,  and on  entry  with
c           TRANS = 'C' or 'c',  K  specifies  the number of rows of the
c           matrices  A and B.  K must be at least zero.
c           Unchanged on exit.
c
c  ALPHA  - COMPLEX*16         .
c           On entry, ALPHA specifies the scalar alpha.
c           Unchanged on exit.
c
c  A      - COMPLEX*16       array of DIMENSION ( LDA, ka ), where ka is
c           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
c           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
c           part of the array  A  must contain the matrix  A,  otherwise
c           the leading  k by n  part of the array  A  must contain  the
c           matrix A.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
c           then  LDA must be at least  max( 1, n ), otherwise  LDA must
c           be at least  max( 1, k ).
c           Unchanged on exit.
c
c  B      - COMPLEX*16       array of DIMENSION ( LDB, kb ), where kb is
c           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
c           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
c           part of the array  B  must contain the matrix  B,  otherwise
c           the leading  k by n  part of the array  B  must contain  the
c           matrix B.
c           Unchanged on exit.
c
c  LDB    - INTEGER.
c           On entry, LDB specifies the first dimension of B as declared
c           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
c           then  LDB must be at least  max( 1, n ), otherwise  LDB must
c           be at least  max( 1, k ).
c           Unchanged on exit.
c
c  BETA   - DOUBLE PRECISION            .
c           On entry, BETA specifies the scalar beta.
c           Unchanged on exit.
c
c  C      - COMPLEX*16          array of DIMENSION ( LDC, n ).
c           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
c           upper triangular part of the array C must contain the upper
c           triangular part  of the  hermitian matrix  and the strictly
c           lower triangular part of C is not referenced.  On exit, the
c           upper triangular part of the array  C is overwritten by the
c           upper triangular part of the updated matrix.
c           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
c           lower triangular part of the array C must contain the lower
c           triangular part  of the  hermitian matrix  and the strictly
c           upper triangular part of C is not referenced.  On exit, the
c           lower triangular part of the array  C is overwritten by the
c           lower triangular part of the updated matrix.
c           Note that the imaginary parts of the diagonal elements need
c           not be set,  they are assumed to be zero,  and on exit they
c           are set to zero.
c
c  LDC    - INTEGER.
c           On entry, LDC specifies the first dimension of C as declared
c           in  the  calling  (sub)  program.   LDC  must  be  at  least
c           max( 1, n ).
c           Unchanged on exit.
c
c
c  Level 3 Blas routine.
c
c  -- Written on 8-February-1989.
c     Jack Dongarra, Argonne National Laboratory.
c     Iain Duff, AERE Harwell.
c     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c     Sven Hammarling, Numerical Algorithms Group Ltd.
c
c  -- Modified 8-Nov-93 to set C(J,J) to DBLE( C(J,J) ) when BETA = 1.
c     Ed Anderson, Cray Research Inc.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     ..
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCONJG, MAX
c     ..
c     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, INFO, J, L, NROWA
      COMPLEX*16         TEMP1, TEMP2
c     ..
c     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
      COMPLEX*16         ZERO
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ) )
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      IF( LSAME( TRANS, 'N' ) ) THEN
         NROWA = N
      ELSE
         NROWA = K
      END IF
      UPPER = LSAME( UPLO, 'U' )
c
      INFO = 0
      IF( ( .NOT.UPPER ) .AND. ( .NOT.LSAME( UPLO, 'L' ) ) ) THEN
         INFO = 1
      ELSE IF( ( .NOT.LSAME( TRANS, 'N' ) ) .AND.
     $         ( .NOT.LSAME( TRANS, 'C' ) ) ) THEN
         INFO = 2
      ELSE IF( N.LT.0 ) THEN
         INFO = 3
      ELSE IF( K.LT.0 ) THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, NROWA ) ) THEN
         INFO = 7
      ELSE IF( LDB.LT.MAX( 1, NROWA ) ) THEN
         INFO = 9
      ELSE IF( LDC.LT.MAX( 1, N ) ) THEN
         INFO = 12
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHER2K', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( ( N.EQ.0 ) .OR. ( ( ( ALPHA.EQ.ZERO ) .OR. ( K.EQ.0 ) ) .AND.
     $    ( BETA.EQ.ONE ) ) )RETURN
c
c     And when  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO ) THEN
         IF( UPPER ) THEN
            IF( BETA.EQ.DBLE( ZERO ) ) THEN
               DO 20 J = 1, N
                  DO 10 I = 1, J
                     C( I, J ) = ZERO
   10             CONTINUE
   20          CONTINUE
            ELSE
               DO 40 J = 1, N
                  DO 30 I = 1, J - 1
                     C( I, J ) = BETA*C( I, J )
   30             CONTINUE
                  C( J, J ) = BETA*DBLE( C( J, J ) )
   40          CONTINUE
            END IF
         ELSE
            IF( BETA.EQ.DBLE( ZERO ) ) THEN
               DO 60 J = 1, N
                  DO 50 I = J, N
                     C( I, J ) = ZERO
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 80 J = 1, N
                  C( J, J ) = BETA*DBLE( C( J, J ) )
                  DO 70 I = J + 1, N
                     C( I, J ) = BETA*C( I, J )
   70             CONTINUE
   80          CONTINUE
            END IF
         END IF
         RETURN
      END IF
c
c     Start the operations.
c
      IF( LSAME( TRANS, 'N' ) ) THEN
c
c        Form  C := alpha*A*conjg( B' ) + conjg( alpha )*B*conjg( A' ) +
c                   C.
c
         IF( UPPER ) THEN
            DO 130 J = 1, N
               IF( BETA.EQ.DBLE( ZERO ) ) THEN
                  DO 90 I = 1, J
                     C( I, J ) = ZERO
   90             CONTINUE
               ELSE IF( BETA.NE.ONE ) THEN
                  DO 100 I = 1, J - 1
                     C( I, J ) = BETA*C( I, J )
  100             CONTINUE
                  C( J, J ) = BETA*DBLE( C( J, J ) )
               ELSE
                  C( J, J ) = DBLE( C( J, J ) )
               END IF
               DO 120 L = 1, K
                  IF( ( A( J, L ).NE.ZERO ) .OR. ( B( J, L ).NE.ZERO ) )
     $                 THEN
                     TEMP1 = ALPHA*DCONJG( B( J, L ) )
                     TEMP2 = DCONJG( ALPHA*A( J, L ) )
                     DO 110 I = 1, J - 1
                        C( I, J ) = C( I, J ) + A( I, L )*TEMP1 +
     $                              B( I, L )*TEMP2
  110                CONTINUE
                     C( J, J ) = DBLE( C( J, J ) ) +
     $                           DBLE( A( J, L )*TEMP1+B( J, L )*TEMP2 )
                  END IF
  120          CONTINUE
  130       CONTINUE
         ELSE
            DO 180 J = 1, N
               IF( BETA.EQ.DBLE( ZERO ) ) THEN
                  DO 140 I = J, N
                     C( I, J ) = ZERO
  140             CONTINUE
               ELSE IF( BETA.NE.ONE ) THEN
                  DO 150 I = J + 1, N
                     C( I, J ) = BETA*C( I, J )
  150             CONTINUE
                  C( J, J ) = BETA*DBLE( C( J, J ) )
               ELSE
                  C( J, J ) = DBLE( C( J, J ) )
               END IF
               DO 170 L = 1, K
                  IF( ( A( J, L ).NE.ZERO ) .OR. ( B( J, L ).NE.ZERO ) )
     $                 THEN
                     TEMP1 = ALPHA*DCONJG( B( J, L ) )
                     TEMP2 = DCONJG( ALPHA*A( J, L ) )
                     DO 160 I = J + 1, N
                        C( I, J ) = C( I, J ) + A( I, L )*TEMP1 +
     $                              B( I, L )*TEMP2
  160                CONTINUE
                     C( J, J ) = DBLE( C( J, J ) ) +
     $                           DBLE( A( J, L )*TEMP1+B( J, L )*TEMP2 )
                  END IF
  170          CONTINUE
  180       CONTINUE
         END IF
      ELSE
c
c        Form  C := alpha*conjg( A' )*B + conjg( alpha )*conjg( B' )*A +
c                   C.
c
         IF( UPPER ) THEN
            DO 210 J = 1, N
               DO 200 I = 1, J
                  TEMP1 = ZERO
                  TEMP2 = ZERO
                  DO 190 L = 1, K
                     TEMP1 = TEMP1 + DCONJG( A( L, I ) )*B( L, J )
                     TEMP2 = TEMP2 + DCONJG( B( L, I ) )*A( L, J )
  190             CONTINUE
                  IF( I.EQ.J ) THEN
                     IF( BETA.EQ.DBLE( ZERO ) ) THEN
                        C( J, J ) = DBLE( ALPHA*TEMP1+DCONJG( ALPHA )*
     $                              TEMP2 )
                     ELSE
                        C( J, J ) = BETA*DBLE( C( J, J ) ) +
     $                              DBLE( ALPHA*TEMP1+DCONJG( ALPHA )*
     $                              TEMP2 )
                     END IF
                  ELSE
                     IF( BETA.EQ.DBLE( ZERO ) ) THEN
                        C( I, J ) = ALPHA*TEMP1 + DCONJG( ALPHA )*TEMP2
                     ELSE
                        C( I, J ) = BETA*C( I, J ) + ALPHA*TEMP1 +
     $                              DCONJG( ALPHA )*TEMP2
                     END IF
                  END IF
  200          CONTINUE
  210       CONTINUE
         ELSE
            DO 240 J = 1, N
               DO 230 I = J, N
                  TEMP1 = ZERO
                  TEMP2 = ZERO
                  DO 220 L = 1, K
                     TEMP1 = TEMP1 + DCONJG( A( L, I ) )*B( L, J )
                     TEMP2 = TEMP2 + DCONJG( B( L, I ) )*A( L, J )
  220             CONTINUE
                  IF( I.EQ.J ) THEN
                     IF( BETA.EQ.DBLE( ZERO ) ) THEN
                        C( J, J ) = DBLE( ALPHA*TEMP1+DCONJG( ALPHA )*
     $                              TEMP2 )
                     ELSE
                        C( J, J ) = BETA*DBLE( C( J, J ) ) +
     $                              DBLE( ALPHA*TEMP1+DCONJG( ALPHA )*
     $                              TEMP2 )
                     END IF
                  ELSE
                     IF( BETA.EQ.DBLE( ZERO ) ) THEN
                        C( I, J ) = ALPHA*TEMP1 + DCONJG( ALPHA )*TEMP2
                     ELSE
                        C( I, J ) = BETA*C( I, J ) + ALPHA*TEMP1 +
     $                              DCONJG( ALPHA )*TEMP2
                     END IF
                  END IF
  230          CONTINUE
  240       CONTINUE
         END IF
      END IF
c
      RETURN
c
c     End of ZHER2K.
c
      END
      SUBROUTINE ZHERK( UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, C, LDC )
c     .. Scalar Arguments ..
      CHARACTER          TRANS, UPLO
      INTEGER            K, LDA, LDC, N
      DOUBLE PRECISION   ALPHA, BETA
c     ..
c     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), C( LDC, * )
c     ..
c
c  Purpose
c  =======
c
c  ZHERK  performs one of the hermitian rank k operations
c
c     C := alpha*A*conjg( A' ) + beta*C,
c
c  or
c
c     C := alpha*conjg( A' )*A + beta*C,
c
c  where  alpha and beta  are  real scalars,  C is an  n by n  hermitian
c  matrix and  A  is an  n by k  matrix in the  first case and a  k by n
c  matrix in the second case.
c
c  Parameters
c  ==========
c
c  UPLO   - CHARACTER*1.
c           On  entry,   UPLO  specifies  whether  the  upper  or  lower
c           triangular  part  of the  array  C  is to be  referenced  as
c           follows:
c
c              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
c                                  is to be referenced.
c
c              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
c                                  is to be referenced.
c
c           Unchanged on exit.
c
c  TRANS  - CHARACTER*1.
c           On entry,  TRANS  specifies the operation to be performed as
c           follows:
c
c              TRANS = 'N' or 'n'   C := alpha*A*conjg( A' ) + beta*C.
c
c              TRANS = 'C' or 'c'   C := alpha*conjg( A' )*A + beta*C.
c
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry,  N specifies the order of the matrix C.  N must be
c           at least zero.
c           Unchanged on exit.
c
c  K      - INTEGER.
c           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
c           of  columns   of  the   matrix   A,   and  on   entry   with
c           TRANS = 'C' or 'c',  K  specifies  the number of rows of the
c           matrix A.  K must be at least zero.
c           Unchanged on exit.
c
c  ALPHA  - DOUBLE PRECISION            .
c           On entry, ALPHA specifies the scalar alpha.
c           Unchanged on exit.
c
c  A      - COMPLEX*16       array of DIMENSION ( LDA, ka ), where ka is
c           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
c           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
c           part of the array  A  must contain the matrix  A,  otherwise
c           the leading  k by n  part of the array  A  must contain  the
c           matrix A.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
c           then  LDA must be at least  max( 1, n ), otherwise  LDA must
c           be at least  max( 1, k ).
c           Unchanged on exit.
c
c  BETA   - DOUBLE PRECISION.
c           On entry, BETA specifies the scalar beta.
c           Unchanged on exit.
c
c  C      - COMPLEX*16          array of DIMENSION ( LDC, n ).
c           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
c           upper triangular part of the array C must contain the upper
c           triangular part  of the  hermitian matrix  and the strictly
c           lower triangular part of C is not referenced.  On exit, the
c           upper triangular part of the array  C is overwritten by the
c           upper triangular part of the updated matrix.
c           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
c           lower triangular part of the array C must contain the lower
c           triangular part  of the  hermitian matrix  and the strictly
c           upper triangular part of C is not referenced.  On exit, the
c           lower triangular part of the array  C is overwritten by the
c           lower triangular part of the updated matrix.
c           Note that the imaginary parts of the diagonal elements need
c           not be set,  they are assumed to be zero,  and on exit they
c           are set to zero.
c
c  LDC    - INTEGER.
c           On entry, LDC specifies the first dimension of C as declared
c           in  the  calling  (sub)  program.   LDC  must  be  at  least
c           max( 1, n ).
c           Unchanged on exit.
c
c
c  Level 3 Blas routine.
c
c  -- Written on 8-February-1989.
c     Jack Dongarra, Argonne National Laboratory.
c     Iain Duff, AERE Harwell.
c     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c     Sven Hammarling, Numerical Algorithms Group Ltd.
c
c  -- Modified 8-Nov-93 to set C(J,J) to DBLE( C(J,J) ) when BETA = 1.
c     Ed Anderson, Cray Research Inc.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     ..
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     ..
c     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, DCONJG, MAX
c     ..
c     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, INFO, J, L, NROWA
      DOUBLE PRECISION   RTEMP
      COMPLEX*16         TEMP
c     ..
c     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      IF( LSAME( TRANS, 'N' ) ) THEN
         NROWA = N
      ELSE
         NROWA = K
      END IF
      UPPER = LSAME( UPLO, 'U' )
c
      INFO = 0
      IF( ( .NOT.UPPER ) .AND. ( .NOT.LSAME( UPLO, 'L' ) ) ) THEN
         INFO = 1
      ELSE IF( ( .NOT.LSAME( TRANS, 'N' ) ) .AND.
     $         ( .NOT.LSAME( TRANS, 'C' ) ) ) THEN
         INFO = 2
      ELSE IF( N.LT.0 ) THEN
         INFO = 3
      ELSE IF( K.LT.0 ) THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, NROWA ) ) THEN
         INFO = 7
      ELSE IF( LDC.LT.MAX( 1, N ) ) THEN
         INFO = 10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHERK ', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( ( N.EQ.0 ) .OR. ( ( ( ALPHA.EQ.ZERO ) .OR. ( K.EQ.0 ) ) .AND.
     $    ( BETA.EQ.ONE ) ) )RETURN
c
c     And when  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO ) THEN
         IF( UPPER ) THEN
            IF( BETA.EQ.ZERO ) THEN
               DO 20 J = 1, N
                  DO 10 I = 1, J
                     C( I, J ) = ZERO
   10             CONTINUE
   20          CONTINUE
            ELSE
               DO 40 J = 1, N
                  DO 30 I = 1, J - 1
                     C( I, J ) = BETA*C( I, J )
   30             CONTINUE
                  C( J, J ) = BETA*DBLE( C( J, J ) )
   40          CONTINUE
            END IF
         ELSE
            IF( BETA.EQ.ZERO ) THEN
               DO 60 J = 1, N
                  DO 50 I = J, N
                     C( I, J ) = ZERO
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 80 J = 1, N
                  C( J, J ) = BETA*DBLE( C( J, J ) )
                  DO 70 I = J + 1, N
                     C( I, J ) = BETA*C( I, J )
   70             CONTINUE
   80          CONTINUE
            END IF
         END IF
         RETURN
      END IF
c
c     Start the operations.
c
      IF( LSAME( TRANS, 'N' ) ) THEN
c
c        Form  C := alpha*A*conjg( A' ) + beta*C.
c
         IF( UPPER ) THEN
            DO 130 J = 1, N
               IF( BETA.EQ.ZERO ) THEN
                  DO 90 I = 1, J
                     C( I, J ) = ZERO
   90             CONTINUE
               ELSE IF( BETA.NE.ONE ) THEN
                  DO 100 I = 1, J - 1
                     C( I, J ) = BETA*C( I, J )
  100             CONTINUE
                  C( J, J ) = BETA*DBLE( C( J, J ) )
               ELSE
                  C( J, J ) = DBLE( C( J, J ) )
               END IF
               DO 120 L = 1, K
                  IF( A( J, L ).NE.DCMPLX( ZERO ) ) THEN
                     TEMP = ALPHA*DCONJG( A( J, L ) )
                     DO 110 I = 1, J - 1
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  110                CONTINUE
                     C( J, J ) = DBLE( C( J, J ) ) +
     $                           DBLE( TEMP*A( I, L ) )
                  END IF
  120          CONTINUE
  130       CONTINUE
         ELSE
            DO 180 J = 1, N
               IF( BETA.EQ.ZERO ) THEN
                  DO 140 I = J, N
                     C( I, J ) = ZERO
  140             CONTINUE
               ELSE IF( BETA.NE.ONE ) THEN
                  C( J, J ) = BETA*DBLE( C( J, J ) )
                  DO 150 I = J + 1, N
                     C( I, J ) = BETA*C( I, J )
  150             CONTINUE
               ELSE
                  C( J, J ) = DBLE( C( J, J ) )
               END IF
               DO 170 L = 1, K
                  IF( A( J, L ).NE.DCMPLX( ZERO ) ) THEN
                     TEMP = ALPHA*DCONJG( A( J, L ) )
                     C( J, J ) = DBLE( C( J, J ) ) +
     $                           DBLE( TEMP*A( J, L ) )
                     DO 160 I = J + 1, N
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  160                CONTINUE
                  END IF
  170          CONTINUE
  180       CONTINUE
         END IF
      ELSE
c
c        Form  C := alpha*conjg( A' )*A + beta*C.
c
         IF( UPPER ) THEN
            DO 220 J = 1, N
               DO 200 I = 1, J - 1
                  TEMP = ZERO
                  DO 190 L = 1, K
                     TEMP = TEMP + DCONJG( A( L, I ) )*A( L, J )
  190             CONTINUE
                  IF( BETA.EQ.ZERO ) THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  200          CONTINUE
               RTEMP = ZERO
               DO 210 L = 1, K
                  RTEMP = RTEMP + DCONJG( A( L, J ) )*A( L, J )
  210          CONTINUE
               IF( BETA.EQ.ZERO ) THEN
                  C( J, J ) = ALPHA*RTEMP
               ELSE
                  C( J, J ) = ALPHA*RTEMP + BETA*DBLE( C( J, J ) )
               END IF
  220       CONTINUE
         ELSE
            DO 260 J = 1, N
               RTEMP = ZERO
               DO 230 L = 1, K
                  RTEMP = RTEMP + DCONJG( A( L, J ) )*A( L, J )
  230          CONTINUE
               IF( BETA.EQ.ZERO ) THEN
                  C( J, J ) = ALPHA*RTEMP
               ELSE
                  C( J, J ) = ALPHA*RTEMP + BETA*DBLE( C( J, J ) )
               END IF
               DO 250 I = J + 1, N
                  TEMP = ZERO
                  DO 240 L = 1, K
                     TEMP = TEMP + DCONJG( A( L, I ) )*A( L, J )
  240             CONTINUE
                  IF( BETA.EQ.ZERO ) THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  250          CONTINUE
  260       CONTINUE
         END IF
      END IF
c
      RETURN
c
c     End of ZHERK .
c
      END
      SUBROUTINE ZSYMM ( SIDE, UPLO, M, N, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
c     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO
      INTEGER            M, N, LDA, LDB, LDC
      COMPLEX*16         ALPHA, BETA
c     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDC, * )
c     ..
c
c  Purpose
c  =======
c
c  ZSYMM  performs one of the matrix-matrix operations
c
c     C := alpha*A*B + beta*C,
c
c  or
c
c     C := alpha*B*A + beta*C,
c
c  where  alpha and beta are scalars, A is a symmetric matrix and  B and
c  C are m by n matrices.
c
c  Parameters
c  ==========
c
c  SIDE   - CHARACTER*1.
c           On entry,  SIDE  specifies whether  the  symmetric matrix  A
c           appears on the  left or right  in the  operation as follows:
c
c              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,
c
c              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
c
c           Unchanged on exit.
c
c  UPLO   - CHARACTER*1.
c           On  entry,   UPLO  specifies  whether  the  upper  or  lower
c           triangular  part  of  the  symmetric  matrix   A  is  to  be
c           referenced as follows:
c
c              UPLO = 'U' or 'u'   Only the upper triangular part of the
c                                  symmetric matrix is to be referenced.
c
c              UPLO = 'L' or 'l'   Only the lower triangular part of the
c                                  symmetric matrix is to be referenced.
c
c           Unchanged on exit.
c
c  M      - INTEGER.
c           On entry,  M  specifies the number of rows of the matrix  C.
c           M  must be at least zero.
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry, N specifies the number of columns of the matrix C.
c           N  must be at least zero.
c           Unchanged on exit.
c
c  ALPHA  - COMPLEX*16      .
c           On entry, ALPHA specifies the scalar alpha.
c           Unchanged on exit.
c
c  A      - COMPLEX*16       array of DIMENSION ( LDA, ka ), where ka is
c           m  when  SIDE = 'L' or 'l'  and is n  otherwise.
c           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of
c           the array  A  must contain the  symmetric matrix,  such that
c           when  UPLO = 'U' or 'u', the leading m by m upper triangular
c           part of the array  A  must contain the upper triangular part
c           of the  symmetric matrix and the  strictly  lower triangular
c           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
c           the leading  m by m  lower triangular part  of the  array  A
c           must  contain  the  lower triangular part  of the  symmetric
c           matrix and the  strictly upper triangular part of  A  is not
c           referenced.
c           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of
c           the array  A  must contain the  symmetric matrix,  such that
c           when  UPLO = 'U' or 'u', the leading n by n upper triangular
c           part of the array  A  must contain the upper triangular part
c           of the  symmetric matrix and the  strictly  lower triangular
c           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
c           the leading  n by n  lower triangular part  of the  array  A
c           must  contain  the  lower triangular part  of the  symmetric
c           matrix and the  strictly upper triangular part of  A  is not
c           referenced.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in the  calling (sub) program. When  SIDE = 'L' or 'l'  then
c           LDA must be at least  max( 1, m ), otherwise  LDA must be at
c           least max( 1, n ).
c           Unchanged on exit.
c
c  B      - COMPLEX*16       array of DIMENSION ( LDB, n ).
c           Before entry, the leading  m by n part of the array  B  must
c           contain the matrix B.
c           Unchanged on exit.
c
c  LDB    - INTEGER.
c           On entry, LDB specifies the first dimension of B as declared
c           in  the  calling  (sub)  program.   LDB  must  be  at  least
c           max( 1, m ).
c           Unchanged on exit.
c
c  BETA   - COMPLEX*16      .
c           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
c           supplied as zero then C need not be set on input.
c           Unchanged on exit.
c
c  C      - COMPLEX*16       array of DIMENSION ( LDC, n ).
c           Before entry, the leading  m by n  part of the array  C must
c           contain the matrix  C,  except when  beta  is zero, in which
c           case C need not be set on entry.
c           On exit, the array  C  is overwritten by the  m by n updated
c           matrix.
c
c  LDC    - INTEGER.
c           On entry, LDC specifies the first dimension of C as declared
c           in  the  calling  (sub)  program.   LDC  must  be  at  least
c           max( 1, m ).
c           Unchanged on exit.
c
c
c  Level 3 Blas routine.
c
c  -- Written on 8-February-1989.
c     Jack Dongarra, Argonne National Laboratory.
c     Iain Duff, AERE Harwell.
c     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c     Sven Hammarling, Numerical Algorithms Group Ltd.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, INFO, J, K, NROWA
      COMPLEX*16         TEMP1, TEMP2
c     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER        ( ONE  = ( 1.0D+0, 0.0D+0 ) )
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
c     ..
c     .. Executable Statements ..
c
c     Set NROWA as the number of rows of A.
c
      IF( LSAME( SIDE, 'L' ) )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      UPPER = LSAME( UPLO, 'U' )
c
c     Test the input parameters.
c
      INFO = 0
      IF(      ( .NOT.LSAME( SIDE, 'L' ) ).AND.
     $         ( .NOT.LSAME( SIDE, 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER              ).AND.
     $         ( .NOT.LSAME( UPLO, 'L' ) )      )THEN
         INFO = 2
      ELSE IF( M  .LT.0               )THEN
         INFO = 3
      ELSE IF( N  .LT.0               )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 7
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 9
      ELSE IF( LDC.LT.MAX( 1, M     ) )THEN
         INFO = 12
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'ZSYMM ', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( ( M.EQ.0 ).OR.( N.EQ.0 ).OR.
     $    ( ( ALPHA.EQ.ZERO ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
c
c     And when  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO )THEN
         IF( BETA.EQ.ZERO )THEN
            DO 20, J = 1, N
               DO 10, I = 1, M
                  C( I, J ) = ZERO
   10          CONTINUE
   20       CONTINUE
         ELSE
            DO 40, J = 1, N
               DO 30, I = 1, M
                  C( I, J ) = BETA*C( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         RETURN
      END IF
c
c     Start the operations.
c
      IF( LSAME( SIDE, 'L' ) )THEN
c
c        Form  C := alpha*A*B + beta*C.
c
         IF( UPPER )THEN
            DO 70, J = 1, N
               DO 60, I = 1, M
                  TEMP1 = ALPHA*B( I, J )
                  TEMP2 = ZERO
                  DO 50, K = 1, I - 1
                     C( K, J ) = C( K, J ) + TEMP1    *A( K, I )
                     TEMP2     = TEMP2     + B( K, J )*A( K, I )
   50             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = TEMP1*A( I, I ) + ALPHA*TEMP2
                  ELSE
                     C( I, J ) = BETA *C( I, J ) +
     $                           TEMP1*A( I, I ) + ALPHA*TEMP2
                  END IF
   60          CONTINUE
   70       CONTINUE
         ELSE
            DO 100, J = 1, N
               DO 90, I = M, 1, -1
                  TEMP1 = ALPHA*B( I, J )
                  TEMP2 = ZERO
                  DO 80, K = I + 1, M
                     C( K, J ) = C( K, J ) + TEMP1    *A( K, I )
                     TEMP2     = TEMP2     + B( K, J )*A( K, I )
   80             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = TEMP1*A( I, I ) + ALPHA*TEMP2
                  ELSE
                     C( I, J ) = BETA *C( I, J ) +
     $                           TEMP1*A( I, I ) + ALPHA*TEMP2
                  END IF
   90          CONTINUE
  100       CONTINUE
         END IF
      ELSE
c
c        Form  C := alpha*B*A + beta*C.
c
         DO 170, J = 1, N
            TEMP1 = ALPHA*A( J, J )
            IF( BETA.EQ.ZERO )THEN
               DO 110, I = 1, M
                  C( I, J ) = TEMP1*B( I, J )
  110          CONTINUE
            ELSE
               DO 120, I = 1, M
                  C( I, J ) = BETA*C( I, J ) + TEMP1*B( I, J )
  120          CONTINUE
            END IF
            DO 140, K = 1, J - 1
               IF( UPPER )THEN
                  TEMP1 = ALPHA*A( K, J )
               ELSE
                  TEMP1 = ALPHA*A( J, K )
               END IF
               DO 130, I = 1, M
                  C( I, J ) = C( I, J ) + TEMP1*B( I, K )
  130          CONTINUE
  140       CONTINUE
            DO 160, K = J + 1, N
               IF( UPPER )THEN
                  TEMP1 = ALPHA*A( J, K )
               ELSE
                  TEMP1 = ALPHA*A( K, J )
               END IF
               DO 150, I = 1, M
                  C( I, J ) = C( I, J ) + TEMP1*B( I, K )
  150          CONTINUE
  160       CONTINUE
  170    CONTINUE
      END IF
c
      RETURN
c
c     End of ZSYMM .
c
      END
      SUBROUTINE ZSYR2K( UPLO, TRANS, N, K, ALPHA, A, LDA, B, LDB,
     $                   BETA, C, LDC )
c     .. Scalar Arguments ..
      CHARACTER*1        UPLO, TRANS
      INTEGER            N, K, LDA, LDB, LDC
      COMPLEX*16         ALPHA, BETA
c     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDC, * )
c     ..
c
c  Purpose
c  =======
c
c  ZSYR2K  performs one of the symmetric rank 2k operations
c
c     C := alpha*A*B' + alpha*B*A' + beta*C,
c
c  or
c
c     C := alpha*A'*B + alpha*B'*A + beta*C,
c
c  where  alpha and beta  are scalars,  C is an  n by n symmetric matrix
c  and  A and B  are  n by k  matrices  in the  first  case  and  k by n
c  matrices in the second case.
c
c  Parameters
c  ==========
c
c  UPLO   - CHARACTER*1.
c           On  entry,   UPLO  specifies  whether  the  upper  or  lower
c           triangular  part  of the  array  C  is to be  referenced  as
c           follows:
c
c              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
c                                  is to be referenced.
c
c              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
c                                  is to be referenced.
c
c           Unchanged on exit.
c
c  TRANS  - CHARACTER*1.
c           On entry,  TRANS  specifies the operation to be performed as
c           follows:
c
c              TRANS = 'N' or 'n'    C := alpha*A*B' + alpha*B*A' +
c                                         beta*C.
c
c              TRANS = 'T' or 't'    C := alpha*A'*B + alpha*B'*A +
c                                         beta*C.
c
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry,  N specifies the order of the matrix C.  N must be
c           at least zero.
c           Unchanged on exit.
c
c  K      - INTEGER.
c           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
c           of  columns  of the  matrices  A and B,  and on  entry  with
c           TRANS = 'T' or 't',  K  specifies  the number of rows of the
c           matrices  A and B.  K must be at least zero.
c           Unchanged on exit.
c
c  ALPHA  - COMPLEX*16      .
c           On entry, ALPHA specifies the scalar alpha.
c           Unchanged on exit.
c
c  A      - COMPLEX*16       array of DIMENSION ( LDA, ka ), where ka is
c           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
c           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
c           part of the array  A  must contain the matrix  A,  otherwise
c           the leading  k by n  part of the array  A  must contain  the
c           matrix A.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
c           then  LDA must be at least  max( 1, n ), otherwise  LDA must
c           be at least  max( 1, k ).
c           Unchanged on exit.
c
c  B      - COMPLEX*16       array of DIMENSION ( LDB, kb ), where kb is
c           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
c           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
c           part of the array  B  must contain the matrix  B,  otherwise
c           the leading  k by n  part of the array  B  must contain  the
c           matrix B.
c           Unchanged on exit.
c
c  LDB    - INTEGER.
c           On entry, LDB specifies the first dimension of B as declared
c           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
c           then  LDB must be at least  max( 1, n ), otherwise  LDB must
c           be at least  max( 1, k ).
c           Unchanged on exit.
c
c  BETA   - COMPLEX*16      .
c           On entry, BETA specifies the scalar beta.
c           Unchanged on exit.
c
c  C      - COMPLEX*16       array of DIMENSION ( LDC, n ).
c           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
c           upper triangular part of the array C must contain the upper
c           triangular part  of the  symmetric matrix  and the strictly
c           lower triangular part of C is not referenced.  On exit, the
c           upper triangular part of the array  C is overwritten by the
c           upper triangular part of the updated matrix.
c           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
c           lower triangular part of the array C must contain the lower
c           triangular part  of the  symmetric matrix  and the strictly
c           upper triangular part of C is not referenced.  On exit, the
c           lower triangular part of the array  C is overwritten by the
c           lower triangular part of the updated matrix.
c
c  LDC    - INTEGER.
c           On entry, LDC specifies the first dimension of C as declared
c           in  the  calling  (sub)  program.   LDC  must  be  at  least
c           max( 1, n ).
c           Unchanged on exit.
c
c
c  Level 3 Blas routine.
c
c  -- Written on 8-February-1989.
c     Jack Dongarra, Argonne National Laboratory.
c     Iain Duff, AERE Harwell.
c     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c     Sven Hammarling, Numerical Algorithms Group Ltd.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, INFO, J, L, NROWA
      COMPLEX*16         TEMP1, TEMP2
c     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER        ( ONE  = ( 1.0D+0, 0.0D+0 ) )
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      IF( LSAME( TRANS, 'N' ) )THEN
         NROWA = N
      ELSE
         NROWA = K
      END IF
      UPPER = LSAME( UPLO, 'U' )
c
      INFO = 0
      IF(      ( .NOT.UPPER               ).AND.
     $         ( .NOT.LSAME( UPLO , 'L' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.LSAME( TRANS, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANS, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( N  .LT.0               )THEN
         INFO = 3
      ELSE IF( K  .LT.0               )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 7
      ELSE IF( LDB.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDC.LT.MAX( 1, N     ) )THEN
         INFO = 12
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'ZSYR2K', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( ( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
c
c     And when  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO )THEN
         IF( UPPER )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 20, J = 1, N
                  DO 10, I = 1, J
                     C( I, J ) = ZERO
   10             CONTINUE
   20          CONTINUE
            ELSE
               DO 40, J = 1, N
                  DO 30, I = 1, J
                     C( I, J ) = BETA*C( I, J )
   30             CONTINUE
   40          CONTINUE
            END IF
         ELSE
            IF( BETA.EQ.ZERO )THEN
               DO 60, J = 1, N
                  DO 50, I = J, N
                     C( I, J ) = ZERO
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 80, J = 1, N
                  DO 70, I = J, N
                     C( I, J ) = BETA*C( I, J )
   70             CONTINUE
   80          CONTINUE
            END IF
         END IF
         RETURN
      END IF
c
c     Start the operations.
c
      IF( LSAME( TRANS, 'N' ) )THEN
c
c        Form  C := alpha*A*B' + alpha*B*A' + C.
c
         IF( UPPER )THEN
            DO 130, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 90, I = 1, J
                     C( I, J ) = ZERO
   90             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 100, I = 1, J
                     C( I, J ) = BETA*C( I, J )
  100             CONTINUE
               END IF
               DO 120, L = 1, K
                  IF( ( A( J, L ).NE.ZERO ).OR.
     $                ( B( J, L ).NE.ZERO )     )THEN
                     TEMP1 = ALPHA*B( J, L )
                     TEMP2 = ALPHA*A( J, L )
                     DO 110, I = 1, J
                        C( I, J ) = C( I, J ) + A( I, L )*TEMP1 +
     $                                          B( I, L )*TEMP2
  110                CONTINUE
                  END IF
  120          CONTINUE
  130       CONTINUE
         ELSE
            DO 180, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 140, I = J, N
                     C( I, J ) = ZERO
  140             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 150, I = J, N
                     C( I, J ) = BETA*C( I, J )
  150             CONTINUE
               END IF
               DO 170, L = 1, K
                  IF( ( A( J, L ).NE.ZERO ).OR.
     $                ( B( J, L ).NE.ZERO )     )THEN
                     TEMP1 = ALPHA*B( J, L )
                     TEMP2 = ALPHA*A( J, L )
                     DO 160, I = J, N
                        C( I, J ) = C( I, J ) + A( I, L )*TEMP1 +
     $                                          B( I, L )*TEMP2
  160                CONTINUE
                  END IF
  170          CONTINUE
  180       CONTINUE
         END IF
      ELSE
c
c        Form  C := alpha*A'*B + alpha*B'*A + C.
c
         IF( UPPER )THEN
            DO 210, J = 1, N
               DO 200, I = 1, J
                  TEMP1 = ZERO
                  TEMP2 = ZERO
                  DO 190, L = 1, K
                     TEMP1 = TEMP1 + A( L, I )*B( L, J )
                     TEMP2 = TEMP2 + B( L, I )*A( L, J )
  190             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP1 + ALPHA*TEMP2
                  ELSE
                     C( I, J ) = BETA *C( I, J ) +
     $                           ALPHA*TEMP1 + ALPHA*TEMP2
                  END IF
  200          CONTINUE
  210       CONTINUE
         ELSE
            DO 240, J = 1, N
               DO 230, I = J, N
                  TEMP1 = ZERO
                  TEMP2 = ZERO
                  DO 220, L = 1, K
                     TEMP1 = TEMP1 + A( L, I )*B( L, J )
                     TEMP2 = TEMP2 + B( L, I )*A( L, J )
  220             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP1 + ALPHA*TEMP2
                  ELSE
                     C( I, J ) = BETA *C( I, J ) +
     $                           ALPHA*TEMP1 + ALPHA*TEMP2
                  END IF
  230          CONTINUE
  240       CONTINUE
         END IF
      END IF
c
      RETURN
c
c     End of ZSYR2K.
c
      END
      SUBROUTINE ZSYRK ( UPLO, TRANS, N, K, ALPHA, A, LDA,
     $                   BETA, C, LDC )
c     .. Scalar Arguments ..
      CHARACTER*1        UPLO, TRANS
      INTEGER            N, K, LDA, LDC
      COMPLEX*16         ALPHA, BETA
c     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), C( LDC, * )
c     ..
c
c  Purpose
c  =======
c
c  ZSYRK  performs one of the symmetric rank k operations
c
c     C := alpha*A*A' + beta*C,
c
c  or
c
c     C := alpha*A'*A + beta*C,
c
c  where  alpha and beta  are scalars,  C is an  n by n symmetric matrix
c  and  A  is an  n by k  matrix in the first case and a  k by n  matrix
c  in the second case.
c
c  Parameters
c  ==========
c
c  UPLO   - CHARACTER*1.
c           On  entry,   UPLO  specifies  whether  the  upper  or  lower
c           triangular  part  of the  array  C  is to be  referenced  as
c           follows:
c
c              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
c                                  is to be referenced.
c
c              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
c                                  is to be referenced.
c
c           Unchanged on exit.
c
c  TRANS  - CHARACTER*1.
c           On entry,  TRANS  specifies the operation to be performed as
c           follows:
c
c              TRANS = 'N' or 'n'   C := alpha*A*A' + beta*C.
c
c              TRANS = 'T' or 't'   C := alpha*A'*A + beta*C.
c
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry,  N specifies the order of the matrix C.  N must be
c           at least zero.
c           Unchanged on exit.
c
c  K      - INTEGER.
c           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
c           of  columns   of  the   matrix   A,   and  on   entry   with
c           TRANS = 'T' or 't',  K  specifies  the number of rows of the
c           matrix A.  K must be at least zero.
c           Unchanged on exit.
c
c  ALPHA  - COMPLEX*16      .
c           On entry, ALPHA specifies the scalar alpha.
c           Unchanged on exit.
c
c  A      - COMPLEX*16       array of DIMENSION ( LDA, ka ), where ka is
c           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
c           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
c           part of the array  A  must contain the matrix  A,  otherwise
c           the leading  k by n  part of the array  A  must contain  the
c           matrix A.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
c           then  LDA must be at least  max( 1, n ), otherwise  LDA must
c           be at least  max( 1, k ).
c           Unchanged on exit.
c
c  BETA   - COMPLEX*16      .
c           On entry, BETA specifies the scalar beta.
c           Unchanged on exit.
c
c  C      - COMPLEX*16       array of DIMENSION ( LDC, n ).
c           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
c           upper triangular part of the array C must contain the upper
c           triangular part  of the  symmetric matrix  and the strictly
c           lower triangular part of C is not referenced.  On exit, the
c           upper triangular part of the array  C is overwritten by the
c           upper triangular part of the updated matrix.
c           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
c           lower triangular part of the array C must contain the lower
c           triangular part  of the  symmetric matrix  and the strictly
c           upper triangular part of C is not referenced.  On exit, the
c           lower triangular part of the array  C is overwritten by the
c           lower triangular part of the updated matrix.
c
c  LDC    - INTEGER.
c           On entry, LDC specifies the first dimension of C as declared
c           in  the  calling  (sub)  program.   LDC  must  be  at  least
c           max( 1, n ).
c           Unchanged on exit.
c
c
c  Level 3 Blas routine.
c
c  -- Written on 8-February-1989.
c     Jack Dongarra, Argonne National Laboratory.
c     Iain Duff, AERE Harwell.
c     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c     Sven Hammarling, Numerical Algorithms Group Ltd.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          MAX
c     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I, INFO, J, L, NROWA
      COMPLEX*16         TEMP
c     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER        ( ONE  = ( 1.0D+0, 0.0D+0 ) )
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      IF( LSAME( TRANS, 'N' ) )THEN
         NROWA = N
      ELSE
         NROWA = K
      END IF
      UPPER = LSAME( UPLO, 'U' )
c
      INFO = 0
      IF(      ( .NOT.UPPER               ).AND.
     $         ( .NOT.LSAME( UPLO , 'L' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.LSAME( TRANS, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANS, 'T' ) )      )THEN
         INFO = 2
      ELSE IF( N  .LT.0               )THEN
         INFO = 3
      ELSE IF( K  .LT.0               )THEN
         INFO = 4
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 7
      ELSE IF( LDC.LT.MAX( 1, N     ) )THEN
         INFO = 10
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'ZSYRK ', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( ( N.EQ.0 ).OR.
     $    ( ( ( ALPHA.EQ.ZERO ).OR.( K.EQ.0 ) ).AND.( BETA.EQ.ONE ) ) )
     $   RETURN
c
c     And when  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO )THEN
         IF( UPPER )THEN
            IF( BETA.EQ.ZERO )THEN
               DO 20, J = 1, N
                  DO 10, I = 1, J
                     C( I, J ) = ZERO
   10             CONTINUE
   20          CONTINUE
            ELSE
               DO 40, J = 1, N
                  DO 30, I = 1, J
                     C( I, J ) = BETA*C( I, J )
   30             CONTINUE
   40          CONTINUE
            END IF
         ELSE
            IF( BETA.EQ.ZERO )THEN
               DO 60, J = 1, N
                  DO 50, I = J, N
                     C( I, J ) = ZERO
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 80, J = 1, N
                  DO 70, I = J, N
                     C( I, J ) = BETA*C( I, J )
   70             CONTINUE
   80          CONTINUE
            END IF
         END IF
         RETURN
      END IF
c
c     Start the operations.
c
      IF( LSAME( TRANS, 'N' ) )THEN
c
c        Form  C := alpha*A*A' + beta*C.
c
         IF( UPPER )THEN
            DO 130, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 90, I = 1, J
                     C( I, J ) = ZERO
   90             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 100, I = 1, J
                     C( I, J ) = BETA*C( I, J )
  100             CONTINUE
               END IF
               DO 120, L = 1, K
                  IF( A( J, L ).NE.ZERO )THEN
                     TEMP = ALPHA*A( J, L )
                     DO 110, I = 1, J
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  110                CONTINUE
                  END IF
  120          CONTINUE
  130       CONTINUE
         ELSE
            DO 180, J = 1, N
               IF( BETA.EQ.ZERO )THEN
                  DO 140, I = J, N
                     C( I, J ) = ZERO
  140             CONTINUE
               ELSE IF( BETA.NE.ONE )THEN
                  DO 150, I = J, N
                     C( I, J ) = BETA*C( I, J )
  150             CONTINUE
               END IF
               DO 170, L = 1, K
                  IF( A( J, L ).NE.ZERO )THEN
                     TEMP      = ALPHA*A( J, L )
                     DO 160, I = J, N
                        C( I, J ) = C( I, J ) + TEMP*A( I, L )
  160                CONTINUE
                  END IF
  170          CONTINUE
  180       CONTINUE
         END IF
      ELSE
c
c        Form  C := alpha*A'*A + beta*C.
c
         IF( UPPER )THEN
            DO 210, J = 1, N
               DO 200, I = 1, J
                  TEMP = ZERO
                  DO 190, L = 1, K
                     TEMP = TEMP + A( L, I )*A( L, J )
  190             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  200          CONTINUE
  210       CONTINUE
         ELSE
            DO 240, J = 1, N
               DO 230, I = J, N
                  TEMP = ZERO
                  DO 220, L = 1, K
                     TEMP = TEMP + A( L, I )*A( L, J )
  220             CONTINUE
                  IF( BETA.EQ.ZERO )THEN
                     C( I, J ) = ALPHA*TEMP
                  ELSE
                     C( I, J ) = ALPHA*TEMP + BETA*C( I, J )
                  END IF
  230          CONTINUE
  240       CONTINUE
         END IF
      END IF
c
      RETURN
c
c     End of ZSYRK .
c
      END
      SUBROUTINE ZTRMM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   B, LDB )
c     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      COMPLEX*16         ALPHA
c     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * )
c     ..
c
c  Purpose
c  =======
c
c  ZTRMM  performs one of the matrix-matrix operations
c
c     B := alpha*op( A )*B,   or   B := alpha*B*op( A )
c
c  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
c  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
c
c     op( A ) = A   or   op( A ) = A'   or   op( A ) = conjg( A' ).
c
c  Parameters
c  ==========
c
c  SIDE   - CHARACTER*1.
c           On entry,  SIDE specifies whether  op( A ) multiplies B from
c           the left or right as follows:
c
c              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
c
c              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
c
c           Unchanged on exit.
c
c  UPLO   - CHARACTER*1.
c           On entry, UPLO specifies whether the matrix A is an upper or
c           lower triangular matrix as follows:
c
c              UPLO = 'U' or 'u'   A is an upper triangular matrix.
c
c              UPLO = 'L' or 'l'   A is a lower triangular matrix.
c
c           Unchanged on exit.
c
c  TRANSA - CHARACTER*1.
c           On entry, TRANSA specifies the form of op( A ) to be used in
c           the matrix multiplication as follows:
c
c              TRANSA = 'N' or 'n'   op( A ) = A.
c
c              TRANSA = 'T' or 't'   op( A ) = A'.
c
c              TRANSA = 'C' or 'c'   op( A ) = conjg( A' ).
c
c           Unchanged on exit.
c
c  DIAG   - CHARACTER*1.
c           On entry, DIAG specifies whether or not A is unit triangular
c           as follows:
c
c              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
c
c              DIAG = 'N' or 'n'   A is not assumed to be unit
c                                  triangular.
c
c           Unchanged on exit.
c
c  M      - INTEGER.
c           On entry, M specifies the number of rows of B. M must be at
c           least zero.
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry, N specifies the number of columns of B.  N must be
c           at least zero.
c           Unchanged on exit.
c
c  ALPHA  - COMPLEX*16      .
c           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
c           zero then  A is not referenced and  B need not be set before
c           entry.
c           Unchanged on exit.
c
c  A      - COMPLEX*16       array of DIMENSION ( LDA, k ), where k is m
c           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
c           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
c           upper triangular part of the array  A must contain the upper
c           triangular matrix  and the strictly lower triangular part of
c           A is not referenced.
c           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
c           lower triangular part of the array  A must contain the lower
c           triangular matrix  and the strictly upper triangular part of
c           A is not referenced.
c           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
c           A  are not referenced either,  but are assumed to be  unity.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
c           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
c           then LDA must be at least max( 1, n ).
c           Unchanged on exit.
c
c  B      - COMPLEX*16       array of DIMENSION ( LDB, n ).
c           Before entry,  the leading  m by n part of the array  B must
c           contain the matrix  B,  and  on exit  is overwritten  by the
c           transformed matrix.
c
c  LDB    - INTEGER.
c           On entry, LDB specifies the first dimension of B as declared
c           in  the  calling  (sub)  program.   LDB  must  be  at  least
c           max( 1, m ).
c           Unchanged on exit.
c
c
c  Level 3 Blas routine.
c
c  -- Written on 8-February-1989.
c     Jack Dongarra, Argonne National Laboratory.
c     Iain Duff, AERE Harwell.
c     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c     Sven Hammarling, Numerical Algorithms Group Ltd.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
c     .. Local Scalars ..
      LOGICAL            LSIDE, NOCONJ, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      COMPLEX*16         TEMP
c     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER        ( ONE  = ( 1.0D+0, 0.0D+0 ) )
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOCONJ = LSAME( TRANSA, 'T' )
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
c
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.
     $         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.
     $         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
     $         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'ZTRMM ', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( N.EQ.0 )
     $   RETURN
c
c     And when  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
c
c     Start the operations.
c
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
c
c           Form  B := alpha*A*B.
c
            IF( UPPER )THEN
               DO 50, J = 1, N
                  DO 40, K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*B( K, J )
                        DO 30, I = 1, K - 1
                           B( I, J ) = B( I, J ) + TEMP*A( I, K )
   30                   CONTINUE
                        IF( NOUNIT )
     $                     TEMP = TEMP*A( K, K )
                        B( K, J ) = TEMP
                     END IF
   40             CONTINUE
   50          CONTINUE
            ELSE
               DO 80, J = 1, N
                  DO 70 K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        TEMP      = ALPHA*B( K, J )
                        B( K, J ) = TEMP
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )*A( K, K )
                        DO 60, I = K + 1, M
                           B( I, J ) = B( I, J ) + TEMP*A( I, K )
   60                   CONTINUE
                     END IF
   70             CONTINUE
   80          CONTINUE
            END IF
         ELSE
c
c           Form  B := alpha*A'*B   or   B := alpha*conjg( A' )*B.
c
            IF( UPPER )THEN
               DO 120, J = 1, N
                  DO 110, I = M, 1, -1
                     TEMP = B( I, J )
                     IF( NOCONJ )THEN
                        IF( NOUNIT )
     $                     TEMP = TEMP*A( I, I )
                        DO 90, K = 1, I - 1
                           TEMP = TEMP + A( K, I )*B( K, J )
   90                   CONTINUE
                     ELSE
                        IF( NOUNIT )
     $                     TEMP = TEMP*DCONJG( A( I, I ) )
                        DO 100, K = 1, I - 1
                           TEMP = TEMP + DCONJG( A( K, I ) )*B( K, J )
  100                   CONTINUE
                     END IF
                     B( I, J ) = ALPHA*TEMP
  110             CONTINUE
  120          CONTINUE
            ELSE
               DO 160, J = 1, N
                  DO 150, I = 1, M
                     TEMP = B( I, J )
                     IF( NOCONJ )THEN
                        IF( NOUNIT )
     $                     TEMP = TEMP*A( I, I )
                        DO 130, K = I + 1, M
                           TEMP = TEMP + A( K, I )*B( K, J )
  130                   CONTINUE
                     ELSE
                        IF( NOUNIT )
     $                     TEMP = TEMP*DCONJG( A( I, I ) )
                        DO 140, K = I + 1, M
                           TEMP = TEMP + DCONJG( A( K, I ) )*B( K, J )
  140                   CONTINUE
                     END IF
                     B( I, J ) = ALPHA*TEMP
  150             CONTINUE
  160          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
c
c           Form  B := alpha*B*A.
c
            IF( UPPER )THEN
               DO 200, J = N, 1, -1
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 170, I = 1, M
                     B( I, J ) = TEMP*B( I, J )
  170             CONTINUE
                  DO 190, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*A( K, J )
                        DO 180, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  180                   CONTINUE
                     END IF
  190             CONTINUE
  200          CONTINUE
            ELSE
               DO 240, J = 1, N
                  TEMP = ALPHA
                  IF( NOUNIT )
     $               TEMP = TEMP*A( J, J )
                  DO 210, I = 1, M
                     B( I, J ) = TEMP*B( I, J )
  210             CONTINUE
                  DO 230, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        TEMP = ALPHA*A( K, J )
                        DO 220, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  220                   CONTINUE
                     END IF
  230             CONTINUE
  240          CONTINUE
            END IF
         ELSE
c
c           Form  B := alpha*B*A'   or   B := alpha*B*conjg( A' ).
c
            IF( UPPER )THEN
               DO 280, K = 1, N
                  DO 260, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        IF( NOCONJ )THEN
                           TEMP = ALPHA*A( J, K )
                        ELSE
                           TEMP = ALPHA*DCONJG( A( J, K ) )
                        END IF
                        DO 250, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  250                   CONTINUE
                     END IF
  260             CONTINUE
                  TEMP = ALPHA
                  IF( NOUNIT )THEN
                     IF( NOCONJ )THEN
                        TEMP = TEMP*A( K, K )
                     ELSE
                        TEMP = TEMP*DCONJG( A( K, K ) )
                     END IF
                  END IF
                  IF( TEMP.NE.ONE )THEN
                     DO 270, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  270                CONTINUE
                  END IF
  280          CONTINUE
            ELSE
               DO 320, K = N, 1, -1
                  DO 300, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        IF( NOCONJ )THEN
                           TEMP = ALPHA*A( J, K )
                        ELSE
                           TEMP = ALPHA*DCONJG( A( J, K ) )
                        END IF
                        DO 290, I = 1, M
                           B( I, J ) = B( I, J ) + TEMP*B( I, K )
  290                   CONTINUE
                     END IF
  300             CONTINUE
                  TEMP = ALPHA
                  IF( NOUNIT )THEN
                     IF( NOCONJ )THEN
                        TEMP = TEMP*A( K, K )
                     ELSE
                        TEMP = TEMP*DCONJG( A( K, K ) )
                     END IF
                  END IF
                  IF( TEMP.NE.ONE )THEN
                     DO 310, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  310                CONTINUE
                  END IF
  320          CONTINUE
            END IF
         END IF
      END IF
c
      RETURN
c
c     End of ZTRMM .
c
      END
      SUBROUTINE ZTRSM ( SIDE, UPLO, TRANSA, DIAG, M, N, ALPHA, A, LDA,
     $                   B, LDB )
c     .. Scalar Arguments ..
      CHARACTER*1        SIDE, UPLO, TRANSA, DIAG
      INTEGER            M, N, LDA, LDB
      COMPLEX*16         ALPHA
c     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * )
c     ..
c
c  Purpose
c  =======
c
c  ZTRSM  solves one of the matrix equations
c
c     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
c
c  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
c  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
c
c     op( A ) = A   or   op( A ) = A'   or   op( A ) = conjg( A' ).
c
c  The matrix X is overwritten on B.
c
c  Parameters
c  ==========
c
c  SIDE   - CHARACTER*1.
c           On entry, SIDE specifies whether op( A ) appears on the left
c           or right of X as follows:
c
c              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
c
c              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
c
c           Unchanged on exit.
c
c  UPLO   - CHARACTER*1.
c           On entry, UPLO specifies whether the matrix A is an upper or
c           lower triangular matrix as follows:
c
c              UPLO = 'U' or 'u'   A is an upper triangular matrix.
c
c              UPLO = 'L' or 'l'   A is a lower triangular matrix.
c
c           Unchanged on exit.
c
c  TRANSA - CHARACTER*1.
c           On entry, TRANSA specifies the form of op( A ) to be used in
c           the matrix multiplication as follows:
c
c              TRANSA = 'N' or 'n'   op( A ) = A.
c
c              TRANSA = 'T' or 't'   op( A ) = A'.
c
c              TRANSA = 'C' or 'c'   op( A ) = conjg( A' ).
c
c           Unchanged on exit.
c
c  DIAG   - CHARACTER*1.
c           On entry, DIAG specifies whether or not A is unit triangular
c           as follows:
c
c              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
c
c              DIAG = 'N' or 'n'   A is not assumed to be unit
c                                  triangular.
c
c           Unchanged on exit.
c
c  M      - INTEGER.
c           On entry, M specifies the number of rows of B. M must be at
c           least zero.
c           Unchanged on exit.
c
c  N      - INTEGER.
c           On entry, N specifies the number of columns of B.  N must be
c           at least zero.
c           Unchanged on exit.
c
c  ALPHA  - COMPLEX*16      .
c           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
c           zero then  A is not referenced and  B need not be set before
c           entry.
c           Unchanged on exit.
c
c  A      - COMPLEX*16       array of DIMENSION ( LDA, k ), where k is m
c           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
c           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
c           upper triangular part of the array  A must contain the upper
c           triangular matrix  and the strictly lower triangular part of
c           A is not referenced.
c           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
c           lower triangular part of the array  A must contain the lower
c           triangular matrix  and the strictly upper triangular part of
c           A is not referenced.
c           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
c           A  are not referenced either,  but are assumed to be  unity.
c           Unchanged on exit.
c
c  LDA    - INTEGER.
c           On entry, LDA specifies the first dimension of A as declared
c           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
c           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
c           then LDA must be at least max( 1, n ).
c           Unchanged on exit.
c
c  B      - COMPLEX*16       array of DIMENSION ( LDB, n ).
c           Before entry,  the leading  m by n part of the array  B must
c           contain  the  right-hand  side  matrix  B,  and  on exit  is
c           overwritten by the solution matrix  X.
c
c  LDB    - INTEGER.
c           On entry, LDB specifies the first dimension of B as declared
c           in  the  calling  (sub)  program.   LDB  must  be  at  least
c           max( 1, m ).
c           Unchanged on exit.
c
c
c  Level 3 Blas routine.
c
c  -- Written on 8-February-1989.
c     Jack Dongarra, Argonne National Laboratory.
c     Iain Duff, AERE Harwell.
c     Jeremy Du Croz, Numerical Algorithms Group Ltd.
c     Sven Hammarling, Numerical Algorithms Group Ltd.
c
c
c     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
c     .. External Subroutines ..
      EXTERNAL           XERBLA
c     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
c     .. Local Scalars ..
      LOGICAL            LSIDE, NOCONJ, NOUNIT, UPPER
      INTEGER            I, INFO, J, K, NROWA
      COMPLEX*16         TEMP
c     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER        ( ONE  = ( 1.0D+0, 0.0D+0 ) )
      COMPLEX*16         ZERO
      PARAMETER        ( ZERO = ( 0.0D+0, 0.0D+0 ) )
c     ..
c     .. Executable Statements ..
c
c     Test the input parameters.
c
      LSIDE  = LSAME( SIDE  , 'L' )
      IF( LSIDE )THEN
         NROWA = M
      ELSE
         NROWA = N
      END IF
      NOCONJ = LSAME( TRANSA, 'T' )
      NOUNIT = LSAME( DIAG  , 'N' )
      UPPER  = LSAME( UPLO  , 'U' )
c
      INFO   = 0
      IF(      ( .NOT.LSIDE                ).AND.
     $         ( .NOT.LSAME( SIDE  , 'R' ) )      )THEN
         INFO = 1
      ELSE IF( ( .NOT.UPPER                ).AND.
     $         ( .NOT.LSAME( UPLO  , 'L' ) )      )THEN
         INFO = 2
      ELSE IF( ( .NOT.LSAME( TRANSA, 'N' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'T' ) ).AND.
     $         ( .NOT.LSAME( TRANSA, 'C' ) )      )THEN
         INFO = 3
      ELSE IF( ( .NOT.LSAME( DIAG  , 'U' ) ).AND.
     $         ( .NOT.LSAME( DIAG  , 'N' ) )      )THEN
         INFO = 4
      ELSE IF( M  .LT.0               )THEN
         INFO = 5
      ELSE IF( N  .LT.0               )THEN
         INFO = 6
      ELSE IF( LDA.LT.MAX( 1, NROWA ) )THEN
         INFO = 9
      ELSE IF( LDB.LT.MAX( 1, M     ) )THEN
         INFO = 11
      END IF
      IF( INFO.NE.0 )THEN
         CALL XERBLA( 'ZTRSM ', INFO )
         RETURN
      END IF
c
c     Quick return if possible.
c
      IF( N.EQ.0 )
     $   RETURN
c
c     And when  alpha.eq.zero.
c
      IF( ALPHA.EQ.ZERO )THEN
         DO 20, J = 1, N
            DO 10, I = 1, M
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
         RETURN
      END IF
c
c     Start the operations.
c
      IF( LSIDE )THEN
         IF( LSAME( TRANSA, 'N' ) )THEN
c
c           Form  B := alpha*inv( A )*B.
c
            IF( UPPER )THEN
               DO 60, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 30, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   30                CONTINUE
                  END IF
                  DO 50, K = M, 1, -1
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 40, I = 1, K - 1
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   40                   CONTINUE
                     END IF
   50             CONTINUE
   60          CONTINUE
            ELSE
               DO 100, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 70, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
   70                CONTINUE
                  END IF
                  DO 90 K = 1, M
                     IF( B( K, J ).NE.ZERO )THEN
                        IF( NOUNIT )
     $                     B( K, J ) = B( K, J )/A( K, K )
                        DO 80, I = K + 1, M
                           B( I, J ) = B( I, J ) - B( K, J )*A( I, K )
   80                   CONTINUE
                     END IF
   90             CONTINUE
  100          CONTINUE
            END IF
         ELSE
c
c           Form  B := alpha*inv( A' )*B
c           or    B := alpha*inv( conjg( A' ) )*B.
c
            IF( UPPER )THEN
               DO 140, J = 1, N
                  DO 130, I = 1, M
                     TEMP = ALPHA*B( I, J )
                     IF( NOCONJ )THEN
                        DO 110, K = 1, I - 1
                           TEMP = TEMP - A( K, I )*B( K, J )
  110                   CONTINUE
                        IF( NOUNIT )
     $                     TEMP = TEMP/A( I, I )
                     ELSE
                        DO 120, K = 1, I - 1
                           TEMP = TEMP - DCONJG( A( K, I ) )*B( K, J )
  120                   CONTINUE
                        IF( NOUNIT )
     $                     TEMP = TEMP/DCONJG( A( I, I ) )
                     END IF
                     B( I, J ) = TEMP
  130             CONTINUE
  140          CONTINUE
            ELSE
               DO 180, J = 1, N
                  DO 170, I = M, 1, -1
                     TEMP = ALPHA*B( I, J )
                     IF( NOCONJ )THEN
                        DO 150, K = I + 1, M
                           TEMP = TEMP - A( K, I )*B( K, J )
  150                   CONTINUE
                        IF( NOUNIT )
     $                     TEMP = TEMP/A( I, I )
                     ELSE
                        DO 160, K = I + 1, M
                           TEMP = TEMP - DCONJG( A( K, I ) )*B( K, J )
  160                   CONTINUE
                        IF( NOUNIT )
     $                     TEMP = TEMP/DCONJG( A( I, I ) )
                     END IF
                     B( I, J ) = TEMP
  170             CONTINUE
  180          CONTINUE
            END IF
         END IF
      ELSE
         IF( LSAME( TRANSA, 'N' ) )THEN
c
c           Form  B := alpha*B*inv( A ).
c
            IF( UPPER )THEN
               DO 230, J = 1, N
                  IF( ALPHA.NE.ONE )THEN
                     DO 190, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  190                CONTINUE
                  END IF
                  DO 210, K = 1, J - 1
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 200, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  200                   CONTINUE
                     END IF
  210             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 220, I = 1, M
                        B( I, J ) = TEMP*B( I, J )
  220                CONTINUE
                  END IF
  230          CONTINUE
            ELSE
               DO 280, J = N, 1, -1
                  IF( ALPHA.NE.ONE )THEN
                     DO 240, I = 1, M
                        B( I, J ) = ALPHA*B( I, J )
  240                CONTINUE
                  END IF
                  DO 260, K = J + 1, N
                     IF( A( K, J ).NE.ZERO )THEN
                        DO 250, I = 1, M
                           B( I, J ) = B( I, J ) - A( K, J )*B( I, K )
  250                   CONTINUE
                     END IF
  260             CONTINUE
                  IF( NOUNIT )THEN
                     TEMP = ONE/A( J, J )
                     DO 270, I = 1, M
                       B( I, J ) = TEMP*B( I, J )
  270                CONTINUE
                  END IF
  280          CONTINUE
            END IF
         ELSE
c
c           Form  B := alpha*B*inv( A' )
c           or    B := alpha*B*inv( conjg( A' ) ).
c
            IF( UPPER )THEN
               DO 330, K = N, 1, -1
                  IF( NOUNIT )THEN
                     IF( NOCONJ )THEN
                        TEMP = ONE/A( K, K )
                     ELSE
                        TEMP = ONE/DCONJG( A( K, K ) )
                     END IF
                     DO 290, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  290                CONTINUE
                  END IF
                  DO 310, J = 1, K - 1
                     IF( A( J, K ).NE.ZERO )THEN
                        IF( NOCONJ )THEN
                           TEMP = A( J, K )
                        ELSE
                           TEMP = DCONJG( A( J, K ) )
                        END IF
                        DO 300, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  300                   CONTINUE
                     END IF
  310             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 320, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  320                CONTINUE
                  END IF
  330          CONTINUE
            ELSE
               DO 380, K = 1, N
                  IF( NOUNIT )THEN
                     IF( NOCONJ )THEN
                        TEMP = ONE/A( K, K )
                     ELSE
                        TEMP = ONE/DCONJG( A( K, K ) )
                     END IF
                     DO 340, I = 1, M
                        B( I, K ) = TEMP*B( I, K )
  340                CONTINUE
                  END IF
                  DO 360, J = K + 1, N
                     IF( A( J, K ).NE.ZERO )THEN
                        IF( NOCONJ )THEN
                           TEMP = A( J, K )
                        ELSE
                           TEMP = DCONJG( A( J, K ) )
                        END IF
                        DO 350, I = 1, M
                           B( I, J ) = B( I, J ) - TEMP*B( I, K )
  350                   CONTINUE
                     END IF
  360             CONTINUE
                  IF( ALPHA.NE.ONE )THEN
                     DO 370, I = 1, M
                        B( I, K ) = ALPHA*B( I, K )
  370                CONTINUE
                  END IF
  380          CONTINUE
            END IF
         END IF
      END IF
c
      RETURN
c
c     End of ZTRSM .
c
      END
