      subroutine grsvd ( nu, nv, nb, m, n, w, matu, u, matv, v, b, irhs,
     & ierr, rv1 )

c*********************************************************************72
c
cc GRSVD determines the singular value decomposition of a matrix with N <= M.
c
c  Discussion:
c
c     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE SVD,
c     NUM. MATH. 14, 403-420(1970) BY GOLUB AND REINSCH.
c     HANDBOOK FOR AUTO. COMP., VOL II-LINEAR ALGEBRA, 134-151(1971).
c
c     THIS SUBROUTINE DETERMINES THE SINGULAR VALUE DECOMPOSITION
c     A=USV' OF A REAL M BY N RECTANGULAR MATRIX.  HOUSEHOLDER
c     BIDIAGONALIZATION AND A VARIANT OF THE QR ALGORITHM ARE USED.
c     GRSVD ASSUMES THAT A COPY OF THE MATRIX A IS IN THE ARRAY U. IT
c     ALSO ASSUMES M .GE. N.  IF M .LT. N, THEN COMPUTE THE SINGULAR
c     VALUE DECOMPOSITION OF A' . IF A' =UWV'  , THEN A=V'WU  .
c
c     GRSVD CAN ALSO BE USED TO COMPUTE THE MINIMAL LENGTH LEAST SQUARES
c     SOLUTION TO THE OVERDETERMINED LINEAR SYSTEM A*X=B.
c
c     ON INPUT-
c
c        NU MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
c          ARRAY PARAMETERS U AS DECLARED IN THE CALLING PROGRAM
c          DIMENSION STATEMENT.  NOTE THAT NU MUST BE AT LEAST
c          AS LARGE AS M,
c
c        NV MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
c          ARRAY PARAMETER V AS DECLARED IN THE CALLING PROGRAM
c          DIMENSION STATEMENT.  NV MUST BE AT LEAST AS LARGE AS N,
c
c        NB MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
c          ARRAY PARAMETERS B AS DECLARED IN THE CALLING PROGRAM
c          DIMENSION STATEMENT.  NOTE THAT NB MUST BE AT LEAST
c          AS LARGE AS M,
c
c        M IS THE NUMBER OF ROWS OF A (AND U),
c
c        N IS THE NUMBER OF COLUMNS OF A (AND U) AND THE ORDER OF V,
c
c        A CONTAINS THE RECTANGULAR INPUT MATRIX TO BE DECOMPOSED,
c
c        B CONTAINS THE IRHS RIGHT-HAND-SIDES OF THE OVERDETERMINED
c          LINEAR SYSTEM A*X=B.  IF IRHS .GT. 0,  THEN ON OUTPUT,
c          THE FIRST N COMPONENTS OF THESE IRHS COLUMNS OF B
c          WILL CONTAIN U'B.  THUS, TO COMPUTE THE MINIMAL LENGTH LEAST
c                                                +
c          SQUARES SOLUTION, ONE MUST COMPUTE V*W  TIMES THE COLUMNS OF
c                    +                        +
c          B, WHERE W  IS A DIAGONAL MATRIX, W (I)=0 IF W(I) IS
c          NEGLIGIBLE, OTHERWISE IS 1/W(I).  IF IRHS=0, B MAY COINCIDE
c          WITH A OR U AND WILL NOT BE REFERENCED,
c
c        IRHS IS THE NUMBER OF RIGHT-HAND-SIDES OF THE OVERDETERMINED
c          SYSTEM A*X=B.  IRHS SHOULD BE SET TO ZERO IF ONLY THE SINGULA
c          VALUE DECOMPOSITION OF A IS DESIRED,
c
c        MATU SHOULD BE SET TO .TRUE. IF THE U MATRIX IN THE
c          DECOMPOSITION IS DESIRED, AND TO .FALSE. OTHERWISE,
c
c        MATV SHOULD BE SET TO .TRUE. IF THE V MATRIX IN THE
c          DECOMPOSITION IS DESIRED, AND TO .FALSE. OTHERWISE.
c
c     ON OUTPUT-
c
c        W CONTAINS THE N (NON-NEGATIVE) SINGULAR VALUES OF A (THE
c          DIAGONAL ELEMENTS OF S).  THEY ARE UNORDERED.  IF AN
c          ERROR EXIT IS MADE, THE SINGULAR VALUES SHOULD BE CORRECT
c          FOR INDICES IERR+1,IERR+2,...,N,
c
c        U CONTAINS THE MATRIX U (ORTHOGONAL COLUMN VECTORS) OF THE
c          DECOMPOSITION IF MATU HAS BEEN SET TO .TRUE.  OTHERWISE
c          U IS USED AS A TEMPORARY ARRAY.
c          IF AN ERROR EXIT IS MADE, THE COLUMNS OF U CORRESPONDING
c          TO INDICES OF CORRECT SINGULAR VALUES SHOULD BE CORRECT,
c
c        V CONTAINS THE MATRIX V (ORTHOGONAL) OF THE DECOMPOSITION IF
c          MATV HAS BEEN SET TO .TRUE.  OTHERWISE V IS NOT REFERENCED.
c          IF AN ERROR EXIT IS MADE, THE COLUMNS OF V CORRESPONDING TO
c          INDICES OF CORRECT SINGULAR VALUES SHOULD BE CORRECT,
c
c        IERR IS SET TO
c          ZERO       FOR NORMAL RETURN,
c          K          IF THE K-TH SINGULAR VALUE HAS NOT BEEN
c                     DETERMINED AFTER 30 ITERATIONS,
c          -1         IF IRHS .LT. 0 ,
c          -2         IF M .LT. N ,
c          -3         IF NU .LT. M .OR. NB .LT. M,
c          -4         IF NV .LT. N .
c
c        RV1 IS A TEMPORARY STORAGE ARRAY.
c
c        THIS SUBROUTINE HAS BEEN CHECKED BY THE PFORT VERIFIER
c        (RYDER, B.G. 'THE PFORT VERIFIER', SOFTWARE - PRACTICE AND
c        EXPERIENCE, VOL.4, 359-377, 1974) FOR ADHERENCE TO A LARGE,
c        CAREFULLY DEFINED, PORTABLE SUBSET OF AMERICAN NATIONAL STANDAR
c        FORTRAN CALLED PFORT.
c
c        ORIGINAL VERSION OF THIS CODE IS SUBROUTINE SVD IN RELEASE 2 OF
c        EISPACK.
c
c        MODIFIED BY TONY F. CHAN,
c                    COMP. SCI. DEPT, YALE UNIV.,
c                    BOX 2158, YALE STATION,
c                    CT 06520
c        LAST MODIFIED : JANUARY, 1982.
c
c  Modified:
c
c    21 June 2012
c
c  Author:
c
c    Tony Chan
c
c  Reference:
c
c    Tony Chan,
c    An improved algorithm for computing the singular value decomposition,
c    ACM Transactions on Mathematical Software,
c    Volume 8, Number 1, March 1982, pages 72-83.
c
c  Parameters:
c
      implicit none

      INTEGER I, J, K, L, M, N, II, I1, KK, K1, LL, L1, MN, NU, NV, NB,
     & ITS, IERR, IRHS
      REAL W(1), U(NU,1), V(NV,1), B(NB,IRHS), RV1(1)
      REAL C, F, G, H, S, X, Y, Z, EPS, SCALE, SRELPR, DUMMY
      REAL SQRT, AMAX1, ABS, SIGN
      LOGICAL MATU, MATV
c
c  SRELPR IS A MACHINE-DEPENDENT FUNCTION SPECIFYING
c  THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
c
      IERR = 0
      IF (IRHS.GE.0) GO TO 10
      IERR = -1
      RETURN
   10 IF (M.GE.N) GO TO 20
      IERR = -2
      RETURN
   20 IF (NU.GE.M .AND. NB.GE.M) GO TO 30
      IERR = -3
      RETURN
   30 IF (NV.GE.N) GO TO 40
      IERR = -4
      RETURN
   40 CONTINUE
c
c  HOUSEHOLDER REDUCTION TO BIDIAGONAL FORM.
c
      G = 0.0
      SCALE = 0.0
      X = 0.0

      DO 260 I=1,N

        L = I + 1
        RV1(I) = SCALE*G
        G = 0.0
        S = 0.0
        SCALE = 0.0
c
c  COMPUTE LEFT TRANSFORMATIONS THAT ZERO THE SUBDIAGONAL ELEMENTS
c  OF THE I-TH COLUMN.
c
        DO K=I,M
          SCALE = SCALE + ABS(U(K,I))
        end do

        IF (SCALE.EQ.0.0) GO TO 160

        DO K=I,M
          U(K,I) = U(K,I)/SCALE
          S = S + U(K,I)**2
        end do

        F = U(I,I)
        G = -SIGN(SQRT(S),F)
        H = F*G - S
        U(I,I) = F - G
        IF (I.EQ.N) GO TO 100
c
c     APPLY LEFT TRANSFORMATIONS TO REMAINING COLUMNS OF A.
c
        DO J=L,N

          S = 0.0
          DO K=I,M
            S = S + U(K,I)*U(K,J)
          end do

          F = S/H

          DO K=I,M
            U(K,J) = U(K,J) + F*U(K,I)
          end do

        end do
c
c  APPLY LEFT TRANSFORMATIONS TO THE COLUMNS OF B IF IRHS .GT. 0
c
  100   continue

        DO J=1,IRHS

          S = 0.0
          DO K=I,M
            S = S + U(K,I)*B(K,J)
          end do

          F = S/H
          DO K=I,M
            B(K,J) = B(K,J) + F*U(K,I)
          end do

        end do
c
c  COMPUTE RIGHT TRANSFORMATIONS.
c
        DO K=I,M
          U(K,I) = SCALE*U(K,I)
        end do

  160   continue

        W(I) = SCALE*G
        G = 0.0
        S = 0.0
        SCALE = 0.0
        IF (I.GT.M .OR. I.EQ.N) GO TO 250

        DO K=L,N
          SCALE = SCALE + ABS(U(I,K))
        end do

        IF (SCALE.EQ.0.0) GO TO 250

        DO K=L,N
          U(I,K) = U(I,K)/SCALE
          S = S + U(I,K)**2
        end do

        F = U(I,L)
        G = -SIGN(SQRT(S),F)
        H = F*G - S
        U(I,L) = F - G

        DO K=L,N
          RV1(K) = U(I,K)/H
        end do

        IF (I.EQ.M) GO TO 230

        DO J=L,M

          S = 0.0

          DO K=L,N
            S = S + U(J,K)*U(I,K)
          end do

          DO K=L,N
            U(J,K) = U(J,K) + S*RV1(K)
          end do

        end do

  230   DO K=L,N
          U(I,K) = SCALE*U(I,K)
        end do

  250   X = AMAX1(X,ABS(W(I))+ABS(RV1(I)))

  260 CONTINUE
c
c  ACCUMULATION OF RIGHT-HAND TRANSFORMATIONS.
c
      IF (.NOT.MATV) GO TO 350
c
c  FOR I=N STEP -1 UNTIL 1 DO -- .
c
      DO 340 II=1,N

        I = N + 1 - II
        IF (I.EQ.N) GO TO 330

        IF ( G .ne. 0.0 ) then

          DO J=L,N
            V(J,I) = (U(I,J)/U(I,L))/G
          end do

          DO J=L,N

            S = 0.0
            DO K=L,N
              S = S + U(I,K)*V(K,J)
            end do

            DO K=L,N
              V(K,J) = V(K,J) + S*V(K,I)
            end do

          end do

        end if

        DO J=L,N
          V(I,J) = 0.0
          V(J,I) = 0.0
        end do

  330   V(I,I) = 1.0
        G = RV1(I)
        L = I
  340 CONTINUE
c
c  ACCUMULATION OF LEFT-HAND TRANSFORMATIONS.
c
  350 IF (.NOT.MATU) GO TO 470
c
c  FOR I=MIN(M,N) STEP -1 UNTIL 1 DO --.
c
      MN = N
      IF (M.LT.N) MN = M

      DO 460 II=1,MN
        I = MN + 1 - II
        L = I + 1
        G = W(I)
        IF (I.EQ.N) GO TO 370

        DO J=L,N
          U(I,J) = 0.0
        end do

  370   IF (G.EQ.0.0) GO TO 430
        IF (I.EQ.MN) GO TO 410

        DO J=L,N

          S = 0.0
          DO K=L,M
            S = S + U(K,I)*U(K,J)
          end do

          F = (S/U(I,I))/G
          DO K=I,M
            U(K,J) = U(K,J) + F*U(K,I)
          end do

        end do

  410   DO J=I,M
          U(J,I) = U(J,I)/G
        end do

        GO TO 450

  430   DO J=I,M
          U(J,I) = 0.0
        end do

  450   U(I,I) = U(I,I) + 1.0

  460 CONTINUE
c
c  DIAGONALIZATION OF THE BIDIAGONAL FORM.
c
  470 EPS = SRELPR(DUMMY)*X
c
c  FOR K=N STEP -1 UNTIL 1 DO --.
c
      DO 650 KK=1,N
        K1 = N - KK
        K = K1 + 1
        ITS = 0
c
c  TEST FOR SPLITTING.
c  FOR L=K STEP -1 UNTIL 1 DO.
c
  480   DO 490 LL=1,K
          L1 = K - LL
          L = L1 + 1
          IF (ABS(RV1(L)).LE.EPS) GO TO 550
c
c  RV1(1) IS ALWAYS ZERO, SO THERE IS NO EXIT THROUGH THE BOTTOM OF THE LOOP.
c
          IF (ABS(W(L1)).LE.EPS) GO TO 500
  490   CONTINUE
c
c  CANCELLATION OF RV1(L) IF L GREATER THAN 1.
c
  500   C = 0.0
        S = 1.0

        DO 540 I=L,K
          F = S*RV1(I)
          RV1(I) = C*RV1(I)
          IF (ABS(F).LE.EPS) GO TO 550
          G = W(I)
          H = SQRT(F*F+G*G)
          W(I) = H
          C = G/H
          S = -F/H
c
c  APPLY LEFT TRANSFORMATIONS TO B IF IRHS .GT. 0
c
          DO J=1,IRHS
            Y = B(L1,J)
            Z = B(I,J)
            B(L1,J) = Y*C + Z*S
            B(I,J) = -Y*S + Z*C
          end do


          IF (.NOT.MATU) GO TO 540

          DO J=1,M
            Y = U(J,L1)
            Z = U(J,I)
            U(J,L1) = Y*C + Z*S
            U(J,I) = -Y*S + Z*C
          end do

  540   CONTINUE
c
c  TEST FOR CONVERGENCE.
c
  550   Z = W(K)
        IF (L.EQ.K) GO TO 630
c
c  SHIFT FROM BOTTOM 2 BY 2 MINOR.
c
        IF (ITS.EQ.30) GO TO 660
        ITS = ITS + 1
        X = W(L)
        Y = W(K1)
        G = RV1(K1)
        H = RV1(K)
        F = ((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0*H*Y)
        G = SQRT(F*F+1.0)
        F = ((X-Z)*(X+Z)+H*(Y/(F+SIGN(G,F))-H))/X
c
c  NEXT QR TRANSFORMATION.
c
        C = 1.0
        S = 1.0

        DO 620 I1=L,K1
          I = I1 + 1
          G = RV1(I)
          Y = W(I)
          H = S*G
          G = C*G
          Z = SQRT(F*F+H*H)
          RV1(I1) = Z
          C = F/Z
          S = H/Z
          F = X*C + G*S
          G = -X*S + G*C
          H = Y*S
          Y = Y*C
          IF (.NOT.MATV) GO TO 570

          DO J=1,N
            X = V(J,I1)
            Z = V(J,I)
            V(J,I1) = X*C + Z*S
            V(J,I) = -X*S + Z*C
          end do

  570     Z = SQRT(F*F+H*H)
          W(I1) = Z
c
c  ROTATION CAN BE ARBITRARY IF Z IS ZERO.
c
          IF ( Z .ne. 0.0 ) then
            C = F/Z
            S = H/Z
          end if

          F = C*G + S*Y
          X = -S*G + C*Y
c
c  APPLY LEFT TRANSFORMATIONS TO B IF IRHS .GT. 0
c
          DO J=1,IRHS
            Y = B(I1,J)
            Z = B(I,J)
            B(I1,J) = Y*C + Z*S
            B(I,J) = -Y*S + Z*C
          end do

          IF ( MATU ) then

            DO J=1,M
              Y = U(J,I1)
              Z = U(J,I)
              U(J,I1) = Y*C + Z*S
              U(J,I) = -Y*S + Z*C
            end do

          end if

  620   CONTINUE

        RV1(L) = 0.0
        RV1(K) = F
        W(K) = X
        GO TO 480
c
c  CONVERGENCE.
c
  630   IF (Z.GE.0.0) GO TO 650
c
c  W(K) IS MADE NON-NEGATIVE.
c
        W(K) = -Z
        IF (.NOT.MATV) GO TO 650

        DO J=1,N
          V(J,K) = -V(J,K)
        end do

  650 CONTINUE

      GO TO 670
c
c  SET ERROR -- NO CONVERGENCE TO A SINGULAR VALUE AFTER 30 ITERATIONS.
c
  660 IERR = K
  670 RETURN
      END
      subroutine hybsvd ( na, nu, nv, nz, nb, m, n, a, w, matu, u, matv,
     &  v, z, b, irhs, ierr, rv1 )

c*********************************************************************72
c
cc HYBSVD computes the singular value decomposition of a rectangular matrix.
c
c  Discussion:
c
c     THIS ROUTINE IS A MODIFICATION OF THE GOLUB-REINSCH PROCEDURE (1)
c                                                           T
c     FOR COMPUTING THE SINGULAR VALUE DECOMPOSITION A = UWV  OF A
c     REAL M BY N RECTANGULAR MATRIX. U IS M BY MIN(M,N) CONTAINING
c     THE LEFT SINGULAR VECTORS, W IS A MIN(M,N) BY MIN(M,N) DIAGONAL
c     MATRIX CONTAINING THE SINGULAR VALUES, AND V IS N BY MIN(M,N)
c     CONTAINING THE RIGHT SINGULAR VECTORS.
c
c     THE ALGORITHM IMPLEMENTED IN THIS
c     ROUTINE HAS A HYBRID NATURE.  WHEN M IS APPROXIMATELY EQUAL TO N,
c     THE GOLUB-REINSCH ALGORITHM IS USED, BUT WHEN EITHER OF THE RATIOS
c     M/N OR N/M IS GREATER THAN ABOUT 2,
c     A MODIFIED VERSION OF THE GOLUB-REINSCH
c     ALGORITHM IS USED.  THIS MODIFIED ALGORITHM FIRST TRANSFORMS A
c                                                                T
c     INTO UPPER TRIANGULAR FORM BY HOUSEHOLDER TRANSFORMATIONS L
c     AND THEN USES THE GOLUB-REINSCH ALGORITHM TO FIND THE SINGULAR
c     VALUE DECOMPOSITION OF THE RESULTING UPPER TRIANGULAR MATRIX R.
c     WHEN U IS NEEDED EXPLICITLY IN THE CASE M.GE.N (OR V IN THE CASE
c     M.LT.N), AN EXTRA ARRAY Z (OF SIZE AT LEAST
c     MIN(M,N)**2) IS NEEDED, BUT OTHERWISE Z IS NOT REFERENCED
c     AND NO EXTRA STORAGE IS REQUIRED.  THIS HYBRID METHOD
c     SHOULD BE MORE EFFICIENT THAN THE GOLUB-REINSCH ALGORITHM WHEN
c     M/N OR N/M IS LARGE.  FOR DETAILS, SEE (2).
c
c     WHEN M .GE. N,
c     HYBSVD CAN ALSO BE USED TO COMPUTE THE MINIMAL LENGTH LEAST
c     SQUARES SOLUTION TO THE OVERDETERMINED LINEAR SYSTEM A*X=B.
c     IF M .LT. N (I.E. FOR UNDERDETERMINED SYSTEMS), THE RHS B
c     IS NOT PROCESSED.
c
c     NOTICE THAT THE SINGULAR VALUE DECOMPOSITION OF A MATRIX
c     IS UNIQUE ONLY UP TO THE SIGN OF THE CORRESPONDING COLUMNS
c     OF U AND V.
c
c     THIS ROUTINE HAS BEEN CHECKED BY THE PFORT VERIFIER (3) FOR
c     ADHERENCE TO A LARGE, CAREFULLY DEFINED, PORTABLE SUBSET OF
c     AMERICAN NATIONAL STANDARD FORTRAN CALLED PFORT.
c
c     REFERENCES:
c
c     (1) GOLUB,G.H. AND REINSCH,C. (1970) 'SINGULAR VALUE
c         DECOMPOSITION AND LEAST SQUARES SOLUTIONS,'
c         NUMER. MATH. 14,403-420, 1970.
c
c     (2) CHAN,T.F. (1982) 'AN IMPROVED ALGORITHM FOR COMPUTING
c         THE SINGULAR VALUE DECOMPOSITION,' ACM TOMS, VOL.8,
c         NO. 1, MARCH, 1982.
c
c     (3) RYDER,B.G. (1974) 'THE PFORT VERIFIER,' SOFTWARE -
c         PRACTICE AND EXPERIENCE, VOL.4, 359-377, 1974.
c
c     ON INPUT:
c
c        NA MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
c          ARRAY PARAMETER A AS DECLARED IN THE CALLING PROGRAM
c          DIMENSION STATEMENT.  NOTE THAT NA MUST BE AT LEAST
c          AS LARGE AS M.
c
c        NU MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
c          ARRAY U AS DECLARED IN THE CALLING PROGRAM DIMENSION
c          STATEMENT. NU MUST BE AT LEAST AS LARGE AS M.
c
c        NV MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
c          ARRAY PARAMETER V AS DECLARED IN THE CALLING PROGRAM
c          DIMENSION STATEMENT. NV MUST BE AT LEAST AS LARGE AS N.
c
c        NZ MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
c          ARRAY PARAMETER Z AS DECLARED IN THE CALLING PROGRAM
c          DIMENSION STATEMENT.  NOTE THAT NZ MUST BE AT LEAST
c          AS LARGE AS MIN(M,N).
c
c        NB MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL
c          ARRAY PARAMETER B AS DECLARED IN THE CALLING PROGRAM
c          DIMENSION STATEMENT. NB MUST BE AT LEAST AS LARGE AS M.
c
c        M IS THE NUMBER OF ROWS OF A (AND U).
c
c        N IS THE NUMBER OF COLUMNS OF A (AND NUMBER OF ROWS OF V).
c
c        A CONTAINS THE RECTANGULAR INPUT MATRIX TO BE DECOMPOSED.
c
c        B CONTAINS THE IRHS RIGHT-HAND-SIDES OF THE OVERDETERMINED
c         LINEAR SYSTEM A*X=B. IF IRHS .GT. 0 AND M .GE. N,
c         THEN ON OUTPUT, THE FIRST N COMPONENTS OF THESE IRHS COLUMNS
c         WILL CONTAIN U' B. THUS, TO COMPUTE THE MINIMAL LENGTH LEAST
c                                               +
c         SQUARES SOLUTION, ONE MUST COMPUTE V*W  TIMES THE COLUMNS OF
c                   +                        +
c         B, WHERE W  IS A DIAGONAL MATRIX, W (I)=0 IF W(I) IS
c         NEGLIGIBLE, OTHERWISE IS 1/W(I). IF IRHS=0 OR M.LT.N,
c         B IS NOT REFERENCED.
c
c        IRHS IS THE NUMBER OF RIGHT-HAND-SIDES OF THE OVERDETERMINED
c         SYSTEM A*X=B. IRHS SHOULD BE SET TO ZERO IF ONLY THE SINGULAR
c         VALUE DECOMPOSITION OF A IS DESIRED.
c
c        MATU SHOULD BE SET TO .TRUE. IF THE U MATRIX IN THE
c          DECOMPOSITION IS DESIRED, AND TO .FALSE. OTHERWISE.
c
c        MATV SHOULD BE SET TO .TRUE. IF THE V MATRIX IN THE
c          DECOMPOSITION IS DESIRED, AND TO .FALSE. OTHERWISE.
c
c        WHEN HYBSVD IS USED TO COMPUTE THE MINIMAL LENGTH LEAST
c        SQUARES SOLUTION TO AN OVERDETERMINED SYSTEM, MATU SHOULD
c        BE SET TO .FALSE. , AND MATV SHOULD BE SET TO .TRUE.  .
c
c     ON OUTPUT:
c
c        A IS UNALTERED (UNLESS OVERWRITTEN BY U OR V).
c
c        W CONTAINS THE (NON-NEGATIVE) SINGULAR VALUES OF A (THE
c          DIAGONAL ELEMENTS OF W).  THEY ARE SORTED IN DESCENDING
c          ORDER.  IF AN ERROR EXIT IS MADE, THE SINGULAR VALUES
c          SHOULD BE CORRECT AND SORTED FOR INDICES IERR+1,...,MIN(M,N).
c
c        U CONTAINS THE MATRIX U (ORTHOGONAL COLUMN VECTORS) OF THE
c          DECOMPOSITION IF MATU HAS BEEN SET TO .TRUE.  IF MATU IS
c          FALSE, THEN U IS EITHER USED AS A TEMPORARY STORAGE (IF
c          M .GE. N) OR NOT REFERENCED (IF M .LT. N).
c          U MAY COINCIDE WITH A IN THE CALLING SEQUENCE.
c          IF AN ERROR EXIT IS MADE, THE COLUMNS OF U CORRESPONDING
c          TO INDICES OF CORRECT SINGULAR VALUES SHOULD BE CORRECT.
c
c        V CONTAINS THE MATRIX V (ORTHOGONAL) OF THE DECOMPOSITION IF
c          MATV HAS BEEN SET TO .TRUE.  IF MATV IS
c          FALSE, THEN V IS EITHER USED AS A TEMPORARY STORAGE (IF
c          M .LT. N) OR NOT REFERENCED (IF M .GE. N).
c          IF M .GE. N, V MAY ALSO COINCIDE WITH A.  IF AN ERROR
c          EXIT IS MADE, THE COLUMNS OF V CORRESPONDING TO INDICES OF
c          CORRECT SINGULAR VALUES SHOULD BE CORRECT.
c
c        Z CONTAINS THE MATRIX X IN THE SINGULAR VALUE DECOMPOSITION
c          OF R=XSY',  IF THE MODIFIED ALGORITHM IS USED. IF THE
c          GOLUB-REINSCH PROCEDURE IS USED, THEN IT IS NOT REFERENCED.
c          IF MATU HAS BEEN SET TO .FALSE. IN THE CASE M.GE.N (OR
c          MATV SET TO .FALSE. IN THE CASE M.LT.N), THEN Z IS NOT
c          REFERENCED AND NO EXTRA STORAGE IS REQUIRED.
c
c        IERR IS SET TO
c          ZERO       FOR NORMAL RETURN,
c          K          IF THE K-TH SINGULAR VALUE HAS NOT BEEN
c                     DETERMINED AFTER 30 ITERATIONS.
c          -1         IF IRHS .LT. 0 .
c          -2         IF M .LT. 1 .OR. N .LT. 1
c          -3         IF NA .LT. M .OR. NU .LT. M .OR. NB .LT. M.
c          -4         IF NV .LT. N .
c          -5         IF NZ .LT. MIN(M,N).
c
c        RV1 IS A TEMPORARY STORAGE ARRAY OF LENGTH AT LEAST MIN(M,N).
c
c     PROGRAMMED BY : TONY CHAN
c                     BOX 2158, YALE STATION,
c                     COMPUTER SCIENCE DEPT, YALE UNIV.,
c                     NEW HAVEN, CT 06520.
c     LAST MODIFIED : JANUARY, 1982.
c
c     HYBSVD USES THE FOLLOWING FUNCTIONS AND SUBROUTINES.
c       INTERNAL  GRSVD, MGNSVD, SRELPR
c       FORTRAN   MIN0,ABS,SQRT,FLOAT,SIGN,AMAX1
c       BLAS      SSWAP
c
c  Modified:
c
c    21 June 2012
c
c  Author:
c
c    Tony Chan
c
c  Reference:
c
c    Tony Chan,
c    An improved algorithm for computing the singular value decomposition,
c    ACM Transactions on Mathematical Software,
c    Volume 8, Number 1, March 1982, pages 72-83.
c
c  Parameters:
c
      implicit none

      integer irhs
      integer na
      integer nb
      integer nu
      integer nv
      integer nz

      real a(na,*)
      real b(nb,irhs)
      integer i
      integer ierr
      integer j
      integer m
      integer n 
      REAL W(*), U(NU,*), V(NV,*), Z(NZ,*), RV1(*)
      LOGICAL MATU, MATV
c
c     ERROR CHECK.
c
      IERR = 0

      IF (IRHS.lt.0) then
        IERR = -1
        RETURN
      end if

      IF (M.GE.1 .AND. N.GE.1) GO TO 20
      IERR = -2
      RETURN
   20 IF (NA.GE.M .AND. NU.GE.M .AND. NB.GE.M) GO TO 30
      IERR = -3
      RETURN
   30 IF (NV.GE.N) GO TO 40
      IERR = -4
      RETURN
   40 IF (NZ.GE.MIN0(M,N)) GO TO 50
      IERR = -5
      RETURN
   50 CONTINUE
c
c  FIRST COPIES A INTO EITHER U OR V ACCORDING TO WHETHER
c  M .GE. N OR M .LT. N, AND THEN CALLS SUBROUTINE MGNSVD
c  WHICH ASSUMES THAT NUMBER OF ROWS .GE. NUMBER OF COLUMNS.
c
      IF ( n .le. m ) then

        DO I=1,M
          DO J=1,N
            U(I,J) = A(I,J)
          end do
        end do

        CALL MGNSVD(NU, NV, NZ, NB, M, N, W, MATU, U, MATV, V, Z, B,
     &    IRHS, IERR, RV1)

      else
 
        DO I=1,M
          DO J=1,N
            V(J,I) = A(I,J)
          end do
        end do

        CALL MGNSVD(NV, NU, NZ, NB, N, M, W, MATV, V, MATU, U, Z, B, 0,
     &    IERR, RV1)

      end if

      RETURN
      END
      subroutine mgnsvd ( nu, nv, nz, nb, m, n, w, matu, u, matv, v, z,
     &  b, irhs, ierr, rv1 )

c*********************************************************************72
c
cc MGNSVD computes the singular value decomposition for a matrix with N <= M.
c
c  Discussion:
c
c     THE DESCRIPTION OF SUBROUTINE MGNSVD IS ALMOST IDENTICAL
c     TO THAT FOR SUBROUTINE HYBSVD ABOVE, WITH THE EXCEPTION
c     THAT MGNSVD ASSUMES M .GE. N.
c     IT ALSO ASSUMES THAT A COPY OF THE MATRIX A IS IN THE ARRAY U.
c
c  Modified:
c
c    21 June 2012
c
c  Author:
c
c    Tony Chan
c
c  Reference:
c
c    Tony Chan,
c    An improved algorithm for computing the singular value decomposition,
c    ACM Transactions on Mathematical Software,
c    Volume 8, Number 1, March 1982, pages 72-83.
c
c  Parameters:
c
      implicit none

      integer nb
      INTEGER NU, NV, NZ, M, N, IRHS, IERR, IP1, I, J, K, IM1, IBACK
      REAL W(*), U(NU,*), V(NV,*), Z(NZ,*), B(NB,IRHS), RV1(*)
      REAL XOVRPT, C, R, G, SCALE, SIGN, ABS, SQRT, F, S, H
      integer id
      integer ierrp1
      integer krhs
      integer nm1
      real t
      REAL FLOAT
      LOGICAL MATU, MATV
c
c  SET VALUE FOR C. THE VALUE FOR C DEPENDS ON THE RELATIVE
c  EFFICIENCY OF FLOATING POINT MULTIPLICATIONS, FLOATING POINT
c  ADDITIONS AND TWO-DIMENSIONAL ARRAY INDEXINGS ON THE
c  COMPUTER WHERE THIS SUBROUTINE IS TO BE RUN.  C SHOULD
c  USUALLY BE BETWEEN 2 AND 4.  FOR DETAILS ON CHOOSING C, SEE
c  (2).  THE ALGORITHM IS NOT SENSITIVE TO THE VALUE OF C
c  ACTUALLY USED AS LONG AS C IS BETWEEN 2 AND 4.
c
      C = 4.0
c
c  DETERMINE CROSS-OVER POINT
c
      IF (MATU .AND. MATV) XOVRPT = (C+7./3.)/C
      IF (MATU .AND. .NOT.MATV) XOVRPT = (C+7./3.)/C
      IF (.NOT.MATU .AND. MATV) XOVRPT = 5./3.
      IF (.NOT.MATU .AND. .NOT.MATV) XOVRPT = 5./3.
c
c  DETERMINE WHETHER TO USE GOLUB-REINSCH OR THE MODIFIED ALGORITHM.
c
      R = FLOAT(M)/FLOAT(N)
      IF (R.GE.XOVRPT) GO TO 10
c
c  USE GOLUB-REINSCH PROCEDURE
c
      CALL GRSVD(NU, NV, NB, M, N, W, MATU, U, MATV, V, B, IRHS, IERR,
     &  RV1)
      GO TO 330
c
c  USE MODIFIED ALGORITHM
c
   10 CONTINUE
c
c  TRIANGULARIZE U BY HOUSEHOLDER TRANSFORMATIONS, USING
c  W AND RV1 AS TEMPORARY STORAGE.
c
      DO I=1,N

        G = 0.0
        S = 0.0
        SCALE = 0.0
c
c  PERFORM SCALING OF COLUMNS TO AVOID UNNECESSARY OVERFLOW OR UNDERFLOW
c
        DO K=I,M
          SCALE = SCALE + ABS(U(K,I))
        end do

        IF (SCALE.EQ.0.0) GO TO 110

        DO K=I,M
          U(K,I) = U(K,I)/SCALE
          S = S + U(K,I)**2
        end do
c
c  THE VECTOR E OF THE HOUSEHOLDER TRANSFORMATION I + EE'/H
c  WILL BE STORED IN COLUMN I OF U. THE TRANSFORMED ELEMENT
c  U(I,I) WILL BE STORED IN W(I) AND THE SCALAR H IN RV1(I).
c
        F = U(I,I)
        G = -SIGN(SQRT(S),F)
        H = F*G - S
        U(I,I) = F - G
        RV1(I) = H
        W(I) = SCALE*G

        IF (I.EQ.N) GO TO 70
c
c  APPLY TRANSFORMATIONS TO REMAINING COLUMNS OF A
c
        IP1 = I + 1

        DO J=IP1,N

          S = 0.0
          DO K=I,M
            S = S + U(K,I)*U(K,J)
          end do

          F = S/H
          DO K=I,M
            U(K,J) = U(K,J) + F*U(K,I)
          end do

        end do
c
c  APPLY TRANSFORMATIONS TO COLUMNS OF B IF IRHS .GT. 0
c
   70   continue

        DO J=1,IRHS

          S = 0.0
          DO K=I,M
            S = S + U(K,I)*B(K,J)
          end do

          F = S/H
          DO K=I,M
            B(K,J) = B(K,J) + F*U(K,I)
          end do

        end do

110     continue

      end do
c
c  COPY R INTO Z IF MATU = .TRUE.
c
      IF (.NOT.MATU) GO TO 290

      DO I=1,N
        DO J=I,N
          Z(J,I) = 0.0
          Z(I,J) = U(I,J)
        end do
        Z(I,I) = W(I)
      end do
c
c  ACCUMULATE HOUSEHOLDER TRANSFORMATIONS IN U
c
      DO 240 IBACK=1,N
        I = N - IBACK + 1
        IP1 = I + 1
        G = W(I)
        H = RV1(I)
        IF (I.EQ.N) GO TO 150

        DO J=IP1,N
          U(I,J) = 0.0
        end do

  150   IF (H.EQ.0.0) GO TO 210
        IF (I.EQ.N) GO TO 190

        DO 180 J=IP1,N
          S = 0.0
          DO 160 K=IP1,M
            S = S + U(K,I)*U(K,J)
  160     CONTINUE
          F = S/H
          DO 170 K=I,M
            U(K,J) = U(K,J) + F*U(K,I)
  170     CONTINUE
  180   CONTINUE

  190   S = U(I,I)/H
        DO J=I,M
          U(J,I) = U(J,I)*S
        end do
        GO TO 230

  210   DO 220 J=I,M
          U(J,I) = 0.0
  220   CONTINUE
  230   U(I,I) = U(I,I) + 1.0
  240 CONTINUE
c
c  COMPUTE SVD OF R (WHICH IS STORED IN Z)
c
      CALL GRSVD(NZ, NV, NB, N, N, W, MATU, Z, MATV, V, B, IRHS, IERR,
     &  RV1)
c                                      T
c  FORM L*X TO OBTAIN U (WHERE R=XWY ). X IS RETURNED IN Z
c  BY GRSVD. THE MATRIX MULTIPLY IS DONE ONE ROW AT A TIME,
c  USING RV1 AS SCRATCH SPACE.
c
      DO 280 I=1,M
        DO 260 J=1,N
          S = 0.0
          DO K=1,N
            S = S + U(I,K)*Z(K,J)
          end do
          RV1(J) = S
  260   CONTINUE
        DO 270 J=1,N
          U(I,J) = RV1(J)
  270   CONTINUE
  280 CONTINUE
      GO TO 330
c
c  FORM R IN U BY ZEROING THE LOWER TRIANGULAR PART OF R IN U
c
  290 IF (N.EQ.1) GO TO 320
      DO 310 I=2,N
        IM1 = I - 1
        DO 300 J=1,IM1
          U(I,J) = 0.0
  300   CONTINUE
        U(I,I) = W(I)
  310 CONTINUE
  320 U(1,1) = W(1)

      CALL GRSVD(NU, NV, NB, N, N, W, MATU, U, MATV, V, B, IRHS, IERR,
     & RV1)
  330 CONTINUE
      IERRP1 = IERR + 1
      IF (IERR.LT.0 .OR. N.LE.1 .OR. IERRP1.EQ.N) RETURN
c
c  SORT SINGULAR VALUES AND EXCHANGE COLUMNS OF U AND V ACCORDINGLY.
c  SELECTION SORT MINIMIZES SWAPPING OF U AND V.
c
      NM1 = N - 1
      DO 360 I=IERRP1,NM1
c
c  FIND INDEX OF MAXIMUM SINGULAR VALUE
c
        ID = I
        IP1 = I + 1
        DO J=IP1,N
          IF (W(J).GT.W(ID)) ID = J
        end do

c
c  SWAP SINGULAR VALUES AND VECTORS
c
        IF (ID.EQ.I) GO TO 360
        T = W(I)
        W(I) = W(ID)
        W(ID) = T
        IF (MATV) CALL SSWAP(N, V(1,I), 1, V(1,ID), 1)
        IF (MATU) CALL SSWAP(M, U(1,I), 1, U(1,ID), 1)

        DO KRHS=1,IRHS
          T = B(I,KRHS)
          B(I,KRHS) = B(ID,KRHS)
          B(ID,KRHS) = T
        end do

  360 CONTINUE

      RETURN
      END
      subroutine minfit ( nm, m, n, a, w, ip, b, ierr, rv1 )

c*********************************************************************72
c
cc MINFIT determines the singular value decomposition of a real M by N matrix.
c
c  Discussion:
c
c     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE MINFIT,
c     NUM. MATH. 14, 403-420(1970) BY GOLUB AND REINSCH.
c     HANDBOOK FOR AUTO. COMP., VOL II-LINEAR ALGEBRA, 134-151(1971).
c
c     THIS SUBROUTINE DETERMINES, TOWARDS THE SOLUTION OF THE LINEAR
c     SYSTEM AX=B, THE SINGULAR VALUE DECOMPOSITION A=USV' OF A REAL
c     M BY N RECTANGULAR MATRIX, FORMING U'B RATHER THAN U.  HOUSEHOLDER
c     BIDIAGONALIZATION AND A VARIANT OF THE QR ALGORITHM ARE USED.
c
c     ON INPUT-
c
c        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
c          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
c          DIMENSION STATEMENT.  NOTE THAT NM MUST BE AT LEAST
c          AS LARGE AS THE MAXIMUM OF M AND N,
c
c        M IS THE NUMBER OF ROWS OF A AND B,
c
c        N IS THE NUMBER OF COLUMNS OF A AND THE ORDER OF V,
c
c        A CONTAINS THE RECTANGULAR COEFFICIENT MATRIX OF THE SYSTEM,
c
c        IP IS THE NUMBER OF COLUMNS OF B.  IP CAN BE ZERO,
c
c        B CONTAINS THE CONSTANT COLUMN MATRIX OF THE SYSTEM
c          IF IP IS NOT ZERO.  OTHERWISE B IS NOT REFERENCED.
c
c     ON OUTPUT-
c
c        A HAS BEEN OVERWRITTEN BY THE MATRIX V (ORTHOGONAL) OF THE
c          DECOMPOSITION IN ITS FIRST N ROWS AND COLUMNS.  IF AN
c          ERROR EXIT IS MADE, THE COLUMNS OF V CORRESPONDING TO
c          INDICES OF CORRECT SINGULAR VALUES SHOULD BE CORRECT,
c
c        W CONTAINS THE N (NON-NEGATIVE) SINGULAR VALUES OF A (THE
c          DIAGONAL ELEMENTS OF S).  THEY ARE UNORDERED.  IF AN
c          ERROR EXIT IS MADE, THE SINGULAR VALUES SHOULD BE CORRECT
c          FOR INDICES IERR+1,IERR+2,...,N,
c
c        B HAS BEEN OVERWRITTEN BY U'B.  IF AN ERROR EXIT IS MADE,
c          THE ROWS OF U'B CORRESPONDING TO INDICES OF CORRECT
c          SINGULAR VALUES SHOULD BE CORRECT,
c
c        IERR IS SET TO
c          ZERO       FOR NORMAL RETURN,
c          K          IF THE K-TH SINGULAR VALUE HAS NOT BEEN
c                     DETERMINED AFTER 30 ITERATIONS,
c
c        RV1 IS A TEMPORARY STORAGE ARRAY.
c
c     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
c     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
c
      implicit none

      integer n
      integer nm

      INTEGER I, J, K, L, M, II, IP, I1, KK, K1, LL, L1, M1,
     & ITS, IERR
      REAL A(NM,N), W(N), B(NM,IP), RV1(N)
      REAL C, F, G, H, S, X, Y, Z, EPS, SCALE, MACHEP
      REAL SQRT, AMAX1, ABS, SIGN
c
c  MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
c  THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
c
      MACHEP = 2.0**(-26)

      IERR = 0
c
c  HOUSEHOLDER REDUCTION TO BIDIAGONAL FORM.
c
      G = 0.0
      SCALE = 0.0
      X = 0.0

      DO 220 I=1,N

        L = I + 1
        RV1(I) = SCALE*G
        G = 0.0
        S = 0.0
        SCALE = 0.0
        IF (I.GT.M) GO TO 120

        DO K=I,M
          SCALE = SCALE + ABS(A(K,I))
        end do

        IF (SCALE.EQ.0.0) GO TO 120

        DO K=I,M
          A(K,I) = A(K,I)/SCALE
          S = S + A(K,I)**2
        end do

        F = A(I,I)
        G = -SIGN(SQRT(S),F)
        H = F*G - S
        A(I,I) = F - G
        IF (I.EQ.N) GO TO 60

        DO J=L,N

          S = 0.0

          DO K=I,M
            S = S + A(K,I)*A(K,J)
          end do

          F = S/H

          DO K=I,M
            A(K,J) = A(K,J) + F*A(K,I)
          end do

        end do

   60   IF (IP.EQ.0) GO TO 100

        DO J=1,IP

          S = 0.0

          DO K=I,M
            S = S + A(K,I)*B(K,J)
          end do

          F = S/H

          DO K=I,M
            B(K,J) = B(K,J) + F*A(K,I)
          end do

        end do

  100   continue

        DO K=I,M
          A(K,I) = SCALE*A(K,I)
        end do

  120   continue

        W(I) = SCALE*G
        G = 0.0
        S = 0.0
        SCALE = 0.0
        IF (I.GT.M .OR. I.EQ.N) GO TO 210

        DO K=L,N
          SCALE = SCALE + ABS(A(I,K))
        end do

        IF (SCALE.EQ.0.0) GO TO 210

        DO K=L,N
          A(I,K) = A(I,K)/SCALE
          S = S + A(I,K)**2
        end do

        F = A(I,L)
        G = -SIGN(SQRT(S),F)
        H = F*G - S
        A(I,L) = F - G

        DO K=L,N
          RV1(K) = A(I,K)/H
        end do

        IF (I.EQ.M) GO TO 190

        DO J=L,M

          S = 0.0

          DO K=L,N
            S = S + A(J,K)*A(I,K)
          end do

          DO K=L,N
            A(J,K) = A(J,K) + S*RV1(K)
          end do

        end do

  190   continue

        DO K=L,N
          A(I,K) = SCALE*A(I,K)
        end do

  210   X = AMAX1(X,ABS(W(I))+ABS(RV1(I)))

  220 CONTINUE
c
c  ACCUMULATION OF RIGHT-HAND TRANSFORMATIONS.
c  FOR I=N STEP -1 UNTIL 1 DO --.
c
      DO 300 II=1,N

        I = N + 1 - II
        IF (I.EQ.N) GO TO 290
        IF (G.EQ.0.0) GO TO 270

        DO J=L,N
          A(J,I) = (A(I,J)/A(I,L))/G
        end do

        DO J=L,N

          S = 0.0
          DO K=L,N
            S = S + A(I,K)*A(K,J)
          end do

          DO K=L,N
            A(K,J) = A(K,J) + S*A(K,I)
          end do

        end do

  270   DO J=L,N
          A(I,J) = 0.0
          A(J,I) = 0.0
        end do

  290   A(I,I) = 1.0
        G = RV1(I)
        L = I
  300 CONTINUE

      IF (M.GE.N .OR. IP.EQ.0) GO TO 330
      M1 = M + 1

      DO I=M1,N
        DO J=1,IP
          B(I,J) = 0.0
        end do
      end do
c
c  DIAGONALIZATION OF THE BIDIAGONAL FORM.
c
  330 EPS = MACHEP*X
c
c  FOR K=N STEP -1 UNTIL 1 DO --.
c
      DO 460 KK=1,N

        K1 = N - KK
        K = K1 + 1
        ITS = 0
c
c  TEST FOR SPLITTING.
c  FOR L=K STEP -1 UNTIL 1 DO --.
c
  340   DO LL=1,K

          L1 = K - LL
          L = L1 + 1
          IF (ABS(RV1(L)).LE.EPS) GO TO 390
c
c  RV1(1) IS ALWAYS ZERO, SO THERE IS NO EXIT THROUGH THE BOTTOM OF THE LOOP.
c
          IF (ABS(W(L1)).LE.EPS) GO TO 360
        end do
c
c  CANCELLATION OF RV1(L) IF L GREATER THAN 1.
c
  360   C = 0.0
        S = 1.0

        DO I=L,K

          F = S*RV1(I)
          RV1(I) = C*RV1(I)
          IF (ABS(F).LE.EPS) GO TO 390
          G = W(I)
          H = SQRT(F*F+G*G)
          W(I) = H
          C = G/H
          S = -F/H

          DO J=1,IP
            Y = B(L1,J)
            Z = B(I,J)
            B(L1,J) = Y*C + Z*S
            B(I,J) = -Y*S + Z*C
          end do

        end do

390     continue
c
c  TEST FOR CONVERGENCE.
c
        Z = W(K)
        IF (L.EQ.K) GO TO 440
c
c  SHIFT FROM BOTTOM 2 BY 2 MINOR.
c
        IF (ITS.EQ.30) GO TO 470
        ITS = ITS + 1
        X = W(L)
        Y = W(K1)
        G = RV1(K1)
        H = RV1(K)
        F = ((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0*H*Y)
        G = SQRT(F*F+1.0)
        F = ((X-Z)*(X+Z)+H*(Y/(F+SIGN(G,F))-H))/X
c
c  NEXT QR TRANSFORMATION.
c
        C = 1.0
        S = 1.0

        DO I1=L,K1

          I = I1 + 1
          G = RV1(I)
          Y = W(I)
          H = S*G
          G = C*G
          Z = SQRT(F*F+H*H)
          RV1(I1) = Z
          C = F/Z
          S = H/Z
          F = X*C + G*S
          G = -X*S + G*C
          H = Y*S
          Y = Y*C

          DO J=1,N
            X = A(J,I1)
            Z = A(J,I)
            A(J,I1) = X*C + Z*S
            A(J,I) = -X*S + Z*C
          end do

          Z = SQRT(F*F+H*H)
          W(I1) = Z
c
c  ROTATION CAN BE ARBITRARY IF Z IS ZERO.
c
          IF ( Z .ne. 0.0 ) then
            C = F/Z
            S = H/Z
          end if

          F = C*G + S*Y
          X = -S*G + C*Y

          DO J=1,IP
            Y = B(I1,J)
            Z = B(I,J)
            B(I1,J) = Y*C + Z*S
            B(I,J) = -Y*S + Z*C
          end do

        end do

        RV1(L) = 0.0
        RV1(K) = F
        W(K) = X
        GO TO 340
c
c  CONVERGENCE
c.
  440   IF (Z.GE.0.0) GO TO 460
c
c  W(K) IS MADE NON-NEGATIVE.
c
        W(K) = -Z

        DO J=1,N
          A(J,K) = -A(J,K)
        end do

  460 CONTINUE

      GO TO 480
c
c  SET ERROR -- NO CONVERGENCE TO A SINGULAR VALUE AFTER 30 ITERATIONS.
c
  470 IERR = K
  480 RETURN
      END
      function sasum ( n, sx, incx )

c*********************************************************************72
c
cc SASUM takes the sum of the absolute values.
c
c  Discussion:
c
c    This routine uses single precision real arithmetic.
c
c    This routine uses unrolled loops for increment equal to one.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, real X(*), the vector to be examined.
c
c    Input, integer INCX, the increment between successive entries of X.
c    INCX must not be negative.
c
c    Output, real SASUM, the sum of the absolute values of X.
c
      implicit none

      integer i
      integer incx
      integer m
      integer n
      integer nincx
      real sasum
      real stemp
      real sx(*)

      sasum = 0.0e0
      stemp = 0.0e0
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c  code for increment not equal to 1
c
      nincx = n*incx
      do i = 1,nincx,incx
        stemp = stemp + abs(sx(i))
      end do
      sasum = stemp
      return
c
c  code for increment equal to 1
c
c
c  clean-up loop
c
   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        stemp = stemp + abs(sx(i))
      end do
      if( n .lt. 6 ) go to 60

   40 continue

      do i = m+1, n, 6
        stemp = stemp + abs(sx(i)) + abs(sx(i + 1)) + abs(sx(i + 2))
     &  + abs(sx(i + 3)) + abs(sx(i + 4)) + abs(sx(i + 5))
      end do

   60 sasum = stemp

      return
      end
      subroutine saxpy ( n, sa, sx, incx, sy, incy )

c*********************************************************************72
c
cc SAXPY computes constant times a vector plus a vector.
c
c  Discussion:
c
c    This routine uses single precision real arithmetic.
c
c    This routine uses unrolled loop for increments equal to one.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, real SA, the multiplier.
c
c    Input, real X(*), the vector to be scaled and added to Y.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Input/output, real Y(*), the vector to which a multiple of X is to 
c    be added.
c
c    Input, integer INCY, the increment between successive entries of Y.
c
      implicit none

      integer i
      integer incx
      integer incy
      integer ix
      integer iy
      integer m
      integer n
      real sa
      real sx(*)
      real sy(*)

      if(n.le.0)return
      if (sa .eq. 0.0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c  code for unequal increments or equal increments not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        sy(iy) = sy(iy) + sa*sx(ix)
        ix = ix + incx
        iy = iy + incy
      end do
      return
c
c  code for both increments equal to 1
c
c  clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        sy(i) = sy(i) + sa*sx(i)
      end do
      if( n .lt. 4 ) return
   40 continue

      do i = m+1, n, 4
        sy(i) = sy(i) + sa*sx(i)
        sy(i + 1) = sy(i + 1) + sa*sx(i + 1)
        sy(i + 2) = sy(i + 2) + sa*sx(i + 2)
        sy(i + 3) = sy(i + 3) + sa*sx(i + 3)
      end do

      return
      end
      function sdot ( n, sx, incx, sy, incy )

c*********************************************************************72
c
cc SDOT forms the dot product of two vectors.
c
c  Discussion:
c
c    This routine uses single precision real arithmetic.
c
c    This routine uses unrolled loops for increments equal to one.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vectors.
c
c    Input, real X(*), one of the vectors to be multiplied.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Input, real Y(*), one of the vectors to be multiplied.
c
c    Input, integer INCY, the increment between successive elements of Y.
c
c    Output, real SDOT, the dot product of X and Y.
c
      implicit none

      integer i,incx,incy,ix,iy,m,n
      real sdot
      real stemp
      real sx(*)
      real sy(*)

      stemp = 0.0e0
      sdot = 0.0e0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c  code for unequal increments or equal increments not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        stemp = stemp + sx(ix)*sy(iy)
        ix = ix + incx
        iy = iy + incy
      end do
      sdot = stemp
      return
c
c  code for both increments equal to 1
c
c  clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        stemp = stemp + sx(i)*sy(i)
      end do

      if( n .lt. 5 ) go to 60
   40 continue

      do i = m + 1, n, 5
        stemp = stemp + sx(i)*sy(i) + sx(i + 1)*sy(i + 1) +
     &   sx(i + 2)*sy(i + 2) + sx(i + 3)*sy(i + 3) + sx(i + 4)*sy(i + 4)
      end do

   60 sdot = stemp

      return
      end
      function snrm2 ( n, x, incx )

c*********************************************************************72
c
cc SNRM2 returns the euclidean norm of a real vector.
c
c  Discussion:
c
c    This routine uses single precision real arithmetic.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Sven Hammarling
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, real X(*), the vector whose norm is to be computed.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Output, real SNRM2, the Euclidean norm of X.
c
      implicit none

      real absxi
      integer incx
      integer ix
      integer n
      real norm
      real scale
      real snrm2
      real ssq
      real x(*)

      real one         , zero
      parameter           ( one = 1.0e+0, zero = 0.0e+0 )

      intrinsic             abs, sqrt

      if( n.lt.1 .or. incx.lt.1 )then
         norm  = zero
      else if( n.eq.1 )then
         norm  = abs( x( 1 ) )
      else
         scale = zero
         ssq   = one
         do ix = 1, 1 + ( n - 1 )*incx, incx
            if( x( ix ).ne.zero )then
               absxi = abs( x( ix ) )
               if( scale.lt.absxi )then
                  ssq   = one   + ssq*( scale/absxi )**2
                  scale = absxi
               else
                  ssq   = ssq   +     ( absxi/scale )**2
               end if
            end if
        end do
        norm  = scale * sqrt( ssq )
      end if

      snrm2 = norm
      return
      end
      function srelpr ( dummy )

c*********************************************************************72
c
cc SRELPR computes the floating point relative precision.
c
c  Discussion:
c
c    IF TROUBLE WITH AUTOMATIC COMPUTATION OF THESE QUANTITIES,
c    THEY CAN BE SET BY DIRECT ASSIGNMENT STATEMENTS.
c    ASSUME THE COMPUTER HAS
c
c      B = BASE OF ARITHMETIC
c      T = NUMBER OF BASE  B  DIGITS
c
c    THEN
c
c      SRELPR = B^(1-T)
c
c  Modified:
c
c    21 June 2012
c
c  Author:
c
c    Tony Chan
c
c  Reference:
c
c    Tony Chan,
c    An improved algorithm for computing the singular value decomposition,
c    ACM Transactions on Mathematical Software,
c    Volume 8, Number 1, March 1982, pages 72-83.
c
c  Parameters:
c
      implicit none

      real dummy
      real s
      real srelpr

      srelpr = 1.0

   10 continue

      srelpr = srelpr / 2.0

      s = 1.0 + srelpr

      if ( s.gt.1.0 ) then
        go to 10
      end if

      srelpr = 2.0 * srelpr

      return
      end
      subroutine srot ( n, sx, incx, sy, incy, c, s )

c*********************************************************************72
c
cc SROT applies a plane rotation.
c
c  Discussion:
c
c    This routine uses single precision real arithmetic.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vectors.
c
c    Input/output, real X(*), one of the vectors to be rotated.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Input/output, real Y(*), one of the vectors to be rotated.
c
c    Input, integer INCY, the increment between successive elements of Y.
c
c    Input, real C, S, parameters (presumably the cosine and sine of
c    some angle) that define a plane rotation.
c
      implicit none

      real sx(*),sy(*),stemp,c,s
      integer i,incx,incy,ix,iy,n

      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c  code for unequal increments or equal increments not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        stemp = c*sx(ix) + s*sy(iy)
        sy(iy) = c*sy(iy) - s*sx(ix)
        sx(ix) = stemp
        ix = ix + incx
        iy = iy + incy
      end do
      return
c
c  code for both increments equal to 1
c
   20 do i = 1,n
        stemp = c*sx(i) + s*sy(i)
        sy(i) = c*sy(i) - s*sx(i)
        sx(i) = stemp
      end do

      return
      end
      subroutine srotg ( sa, sb, c, s )

c*********************************************************************72
c
cc SROTG constructs a Givens plane rotation.
c
c  Discussion:
c
c    This routine uses single precision real arithmetic.
c
c    Given values A and B, this routine computes
c
c    SIGMA = sign ( A ) if abs ( A ) >  abs ( B )
c          = sign ( B ) if abs ( A ) <= abs ( B );
c
c    R     = SIGMA * ( A * A + B * B );
c
c    C = A / R if R is not 0
c      = 1     if R is 0;
c
c    S = B / R if R is not 0,
c        0     if R is 0.
c
c    The computed numbers then satisfy the equation
c
c    (  C  S ) ( A ) = ( R )
c    ( -S  C ) ( B ) = ( 0 )
c
c    The routine also computes
c
c    Z = S     if abs ( A ) > abs ( B ),
c      = 1 / C if abs ( A ) <= abs ( B ) and C is not 0,
c      = 1     if C is 0.
c
c    The single value Z encodes C and S, and hence the rotation:
c
c    If Z = 1, set C = 0 and S = 1;
c    If abs ( Z ) < 1, set C = sqrt ( 1 - Z * Z ) and S = Z;
c    if abs ( Z ) > 1, set C = 1/ Z and S = sqrt ( 1 - C * C );
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input/output, real SA, SB.  On input, SA and SB are the values
c    A and B.  On output, SA is overwritten with R, and SB is
c    overwritten with Z.
c
c    Output, real C, S, the cosine and sine of the Givens rotation.
c
      implicit none

      real sa,sb,c,s,roe,scale,r,z

      roe = sb
      if( abs(sa) .gt. abs(sb) ) roe = sa
      scale = abs(sa) + abs(sb)

      if( scale .eq. 0.0 ) then
         c = 1.0
         s = 0.0
         r = 0.0
         z = 0.0
      else
        r = scale*sqrt((sa/scale)**2 + (sb/scale)**2)
        r = sign(1.0,roe)*r
        c = sa/r
        s = sb/r
        z = 1.0
        if( abs(sa) .gt. abs(sb) ) z = s
        if( abs(sb) .ge. abs(sa) .and. c .ne. 0.0 ) z = 1.0/c
      end if

      sa = r
      sb = z

      return
      end
      subroutine sscal ( n, sa, sx, incx )

c*********************************************************************72
c
cc SSCAL scales a vector by a constant.
c
c  Discussion:
c
c    This routine uses single precision real arithmetic.
c
c    This routine uses unrolled loops for increment equal to 1.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, real SA, the multiplier.
c
c    Input/output, real X(*), the vector to be scaled.
c
c    Input, integer INCX, the increment between successive entries of X.
c
      implicit none

      real sa,sx(*)
      integer i,incx,m,n,nincx

      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c  code for increment not equal to 1
c
      nincx = n*incx
      do i = 1,nincx,incx
        sx(i) = sa*sx(i)
      end do
      return
c
c  code for increment equal to 1
c
c
c  clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        sx(i) = sa*sx(i)
      end do
      if( n .lt. 5 ) return
   40 continue

      do i = m + 1, n, 5
        sx(i) = sa * sx(i)
        sx(i + 1) = sa * sx(i + 1)
        sx(i + 2) = sa * sx(i + 2)
        sx(i + 3) = sa * sx(i + 3)
        sx(i + 4) = sa * sx(i + 4)
      end do

      return
      end
      subroutine ssvdc ( x, ldx, n, p, s, e, u, ldu, v, ldv, work, job,
     &  info )

c*********************************************************************72
c
cc SSVDC orthogonally reduces a real NxP matrix X to diagonal form.
c
c  Discussion:
c
c     SSVDC IS A SUBROUTINE TO REDUCE A REAL NXP MATRIX X BY
c     ORTHOGONAL TRANSFORMATIONS U AND V TO DIAGONAL FORM.  THE
c     DIAGONAL ELEMENTS S(I) ARE THE SINGULAR VALUES OF X.  THE
c     COLUMNS OF U ARE THE CORRESPONDING LEFT SINGULAR VECTORS,
c     AND THE COLUMNS OF V THE RIGHT SINGULAR VECTORS.
c
c     ON ENTRY
c
c         X         REAL(LDX,P), WHERE LDX.GE.N.
c                   X CONTAINS THE MATRIX WHOSE SINGULAR VALUE
c                   DECOMPOSITION IS TO BE COMPUTED.  X IS
c                   DESTROYED BY SSVDC.
c
c         LDX       INTEGER.
c                   LDX IS THE LEADING DIMENSION OF THE ARRAY X.
c
c         N         INTEGER.
c                   N IS THE NUMBER OF COLUMNS OF THE MATRIX X.
c
c         P         INTEGER.
c                   P IS THE NUMBER OF ROWS OF THE MATRIX X.
c
c         LDU       INTEGER.
c                   LDU IS THE LEADING DIMENSION OF THE ARRAY U.
c                   (SEE BELOW).
c
c         LDV       INTEGER.
c                   LDV IS THE LEADING DIMENSION OF THE ARRAY V.
c                   (SEE BELOW).
c
c         WORK      REAL(N).
c                   WORK IS A SCRATCH ARRAY.
c
c         JOB       INTEGER.
c                   JOB CONTROLS THE COMPUTATION OF THE SINGULAR
c                   VECTORS.  IT HAS THE DECIMAL EXPANSION AB
c                   WITH THE FOLLOWING MEANING
c
c                        A.EQ.0    DO NOT COMPUTE THE LEFT SINGULAR
c                                  VECTORS.
c                        A.EQ.1    RETURN THE N LEFT SINGULAR VECTORS
c                                  IN U.
c                        A.GE.2    RETURN THE FIRST MIN(N,P) SINGULAR
c                                  VECTORS IN U.
c                        B.EQ.0    DO NOT COMPUTE THE RIGHT SINGULAR
c                                  VECTORS.
c                        B.EQ.1    RETURN THE RIGHT SINGULAR VECTORS
c                                  IN V.
c
c     ON RETURN
c
c         S         REAL(MM), WHERE MM=MIN(N+1,P).
c                   THE FIRST MIN(N,P) ENTRIES OF S CONTAIN THE
c                   SINGULAR VALUES OF X ARRANGED IN DESCENDING
c                   ORDER OF MAGNITUDE.
c
c         E         REAL(P).
c                   E ORDINARILY CONTAINS ZEROS.  HOWEVER SEE THE
c                   DISCUSSION OF INFO FOR EXCEPTIONS.
c
c         U         REAL(LDU,K), WHERE LDU.GE.N.  IF JOBA.EQ.1 THEN
c                                   K.EQ.N, IF JOBA.GE.2 THEN
c                                   K.EQ.MIN(N,P).
c                   U CONTAINS THE MATRIX OF RIGHT SINGULAR VECTORS.
c                   U IS NOT REFERENCED IF JOBA.EQ.0.  IF N.LE.P
c                   OR IF JOBA.EQ.2, THEN U MAY BE IDENTIFIED WITH X
c                   IN THE SUBROUTINE CALL.
c
c         V         REAL(LDV,P), WHERE LDV.GE.P.
c                   V CONTAINS THE MATRIX OF RIGHT SINGULAR VECTORS.
c                   V IS NOT REFERENCED IF JOB.EQ.0.  IF P.LE.N,
c                   THEN V MAY BE IDENTIFIED WITH X IN THE
c                   SUBROUTINE CALL.
c
c         INFO      INTEGER.
c                   THE SINGULAR VALUES (AND THEIR CORRESPONDING
c                   SINGULAR VECTORS) S(INFO+1),S(INFO+2),...,S(M)
c                   ARE CORRECT (HERE M=MIN(N,P)).  THUS IF
c                   INFO.EQ.0, ALL THE SINGULAR VALUES AND THEIR
c                   VECTORS ARE CORRECT.  IN ANY EVENT, THE MATRIX
c                   B = TRANS(U)*X*V IS THE BIDIAGONAL MATRIX
c                   WITH THE ELEMENTS OF S ON ITS DIAGONAL AND THE
c                   ELEMENTS OF E ON ITS SUPER-DIAGONAL (TRANS(U)
c                   IS THE TRANSPOSE OF U).  THUS THE SINGULAR
c                   VALUES OF X AND B ARE THE SAME.
c
c     LINPACK. THIS VERSION DATED 03/19/79 .
c     G.W. STEWART, UNIVERSITY OF MARYLAND, ARGONNE NATIONAL LAB.
c
      implicit none

      INTEGER LDX, N, P, LDU, LDV, JOB, INFO
      REAL X(LDX,1), S(1), E(1), U(LDU,1), V(LDV,1), WORK(1)

      INTEGER I, ITER, J, JOBU, K, KASE, KK, L, LL, LLS, LM1, LP1, LS,
     & LU, M, MAXIT, MM, MM1, MP1, NCT, NCTP1, NCU, NRT, NRTP1
      REAL SDOT, T, R
      REAL B, C, CS, EL, EMM1, F, G, SNRM2, SCALE, SHIFT, SL, SM, SN,
     & SMM1, T1, TEST, ZTEST
      LOGICAL WANTU, WANTV
c
c  SET THE MAXIMUM NUMBER OF ITERATIONS.
c
      MAXIT = 30
c
c  DETERMINE WHAT IS TO BE COMPUTED.
c
      WANTU = .FALSE.
      WANTV = .FALSE.
      JOBU = MOD(JOB,100)/10
      NCU = N
      IF (JOBU.GT.1) NCU = MIN0(N,P)
      IF (JOBU.NE.0) WANTU = .TRUE.
      IF (MOD(JOB,10).NE.0) WANTV = .TRUE.
c
c  REDUCE X TO BIDIAGONAL FORM, STORING THE DIAGONAL ELEMENTS
c  IN S AND THE SUPER-DIAGONAL ELEMENTS IN E.
c
      INFO = 0
      NCT = MIN0(N-1,P)
      NRT = MAX0(0,MIN0(P-2,N))
      LU = MAX0(NCT,NRT)
      IF (LU.LT.1) GO TO 170

      DO 160 L=1,LU

        LP1 = L + 1
        IF (L.GT.NCT) GO TO 20
c
c  COMPUTE THE TRANSFORMATION FOR THE L-TH COLUMN AND
c  PLACE THE L-TH DIAGONAL IN S(L).
c
        S(L) = SNRM2(N-L+1,X(L,L),1)
        IF (S(L).EQ.0.0E0) GO TO 10
        IF (X(L,L).NE.0.0E0) S(L) = SIGN(S(L),X(L,L))
        CALL SSCAL(N-L+1, 1.0E0/S(L), X(L,L), 1)
        X(L,L) = 1.0E0 + X(L,L)
   10   CONTINUE
        S(L) = -S(L)
   20   CONTINUE
        IF (P.LT.LP1) GO TO 50
        DO 40 J=LP1,P
          IF (L.GT.NCT) GO TO 30
          IF (S(L).EQ.0.0E0) GO TO 30
c
c  APPLY THE TRANSFORMATION.
c
          T = -SDOT(N-L+1,X(L,L),1,X(L,J),1)/X(L,L)
          CALL SAXPY(N-L+1, T, X(L,L), 1, X(L,J), 1)
   30     CONTINUE
c
c  PLACE THE L-TH ROW OF X INTO  E FOR THE
c  SUBSEQUENT CALCULATION OF THE ROW TRANSFORMATION.
c
          E(J) = X(L,J)
   40   CONTINUE
   50   CONTINUE
        IF (.NOT.WANTU .OR. L.GT.NCT) GO TO 70
c
c  PLACE THE TRANSFORMATION IN U FOR SUBSEQUENT BACK MULTIPLICATION.
c
        DO I=L,N
          U(I,L) = X(I,L)
        end do

   70   CONTINUE
        IF (L.GT.NRT) GO TO 150
c
c  COMPUTE THE L-TH ROW TRANSFORMATION AND PLACE THE
c  L-TH SUPER-DIAGONAL IN E(L).
c
        E(L) = SNRM2(P-L,E(LP1),1)
        IF (E(L).EQ.0.0E0) GO TO 80
        IF (E(LP1).NE.0.0E0) E(L) = SIGN(E(L),E(LP1))
        CALL SSCAL(P-L, 1.0E0/E(L), E(LP1), 1)
        E(LP1) = 1.0E0 + E(LP1)
   80   CONTINUE
        E(L) = -E(L)
        IF (LP1.GT.N .OR. E(L).EQ.0.0E0) GO TO 120
c
c  APPLY THE TRANSFORMATION.
c
        DO I=LP1,N
          WORK(I) = 0.0E0
        end do

        DO J=LP1,P
          CALL SAXPY(N-L, E(J), X(LP1,J), 1, WORK(LP1), 1)
        end do

        DO J=LP1,P
          CALL SAXPY(N-L, -E(J)/E(LP1), WORK(LP1), 1, X(LP1,J), 1)
        end do

  120   CONTINUE
        IF (.NOT.WANTV) GO TO 140
c
c  PLACE THE TRANSFORMATION IN V FOR SUBSEQUENT BACK MULTIPLICATION.
c
        DO I=LP1,P
          V(I,L) = E(I)
        end do

  140   CONTINUE
  150   CONTINUE
  160 CONTINUE
  170 CONTINUE
c
c   SET UP THE FINAL BIDIAGONAL MATRIX OR ORDER M.
c
      M = MIN0(P,N+1)
      NCTP1 = NCT + 1
      NRTP1 = NRT + 1
      IF (NCT.LT.P) S(NCTP1) = X(NCTP1,NCTP1)
      IF (N.LT.M) S(M) = 0.0E0
      IF (NRTP1.LT.M) E(NRTP1) = X(NRTP1,M)
      E(M) = 0.0E0
c
c  IF REQUIRED, GENERATE U.
c
      IF (.NOT.WANTU) GO TO 300

      DO J=NCTP1,NCU
        DO I=1,N
          U(I,J) = 0.0E0
        end do
        U(J,J) = 1.0E0
      end do

      DO 280 LL=1,NCT
        L = NCT - LL + 1
        IF (S(L).EQ.0.0E0) GO TO 250
        LP1 = L + 1
        DO J=LP1,NCU
          T = -SDOT(N-L+1,U(L,L),1,U(L,J),1)/U(L,L)
          CALL SAXPY(N-L+1, T, U(L,L), 1, U(L,J), 1)
        end do
        CALL SSCAL(N-L+1, -1.0E0, U(L,L), 1)
        U(L,L) = 1.0E0 + U(L,L)
        LM1 = L - 1
        IF (LM1.LT.1) GO TO 240
        DO 230 I=1,LM1
          U(I,L) = 0.0E0
  230   CONTINUE
  240   CONTINUE
        GO TO 270
  250   CONTINUE
        DO I=1,N
          U(I,L) = 0.0E0
        end do
        U(L,L) = 1.0E0
  270   CONTINUE
  280 CONTINUE

  300 CONTINUE
c
c  IF IT IS REQUIRED, GENERATE V.
c
      IF (.NOT.WANTV) GO TO 350
      DO 340 LL=1,P
        L = P - LL + 1
        LP1 = L + 1
        IF (L.GT.NRT) GO TO 320
        IF (E(L).EQ.0.0E0) GO TO 320
        DO J=LP1,P
          T = -SDOT(P-L,V(LP1,L),1,V(LP1,J),1)/V(LP1,L)
          CALL SAXPY(P-L, T, V(LP1,L), 1, V(LP1,J), 1)
        end do
  320   CONTINUE
        DO 330 I=1,P
          V(I,L) = 0.0E0
  330   CONTINUE
        V(L,L) = 1.0E0
  340 CONTINUE
  350 CONTINUE
c
c  MAIN ITERATION LOOP FOR THE SINGULAR VALUES.
c
      MM = M
      ITER = 0
  360 CONTINUE
c
c   QUIT IF ALL THE SINGULAR VALUES HAVE BEEN FOUND.
c
      IF (M.EQ.0) GO TO 620
c
c  IF TOO MANY ITERATIONS HAVE BEEN PERFORMED, SET FLAG AND RETURN.
c
      IF (ITER.LT.MAXIT) GO TO 370
      INFO = M
      GO TO 620
  370 CONTINUE
c
c  THIS SECTION OF THE PROGRAM INSPECTS FOR
c  NEGLIGIBLE ELEMENTS IN THE S AND E ARRAYS.  ON
c  COMPLETION THE VARIABLES KASE AND L ARE SET AS FOLLOWS.
c
c  KASE = 1     IF S(M) AND E(L-1) ARE NEGLIGIBLE AND L.LT.M
c  KASE = 2     IF S(L) IS NEGLIGIBLE AND L.LT.M
c  KASE = 3     IF E(L-1) IS NEGLIGIBLE, L.LT.M, AND
c               S(L), ..., S(M) ARE NOT NEGLIGIBLE (QR STEP).
c  KASE = 4     IF E(M-1) IS NEGLIGIBLE (CONVERGENCE).
c
      DO 390 LL=1,M
        L = M - LL
        IF (L.EQ.0) GO TO 400
        TEST = ABS(S(L)) + ABS(S(L+1))
        ZTEST = TEST + ABS(E(L))
        IF (ZTEST.NE.TEST) GO TO 380
        E(L) = 0.0E0
        GO TO 400
  380   CONTINUE
  390 CONTINUE

  400 CONTINUE
      IF (L.NE.M-1) GO TO 410
      KASE = 4
      GO TO 480
  410 CONTINUE
      LP1 = L + 1
      MP1 = M + 1
      DO 430 LLS=LP1,MP1
        LS = M - LLS + LP1
        IF (LS.EQ.L) GO TO 440
        TEST = 0.0E0
        IF (LS.NE.M) TEST = TEST + ABS(E(LS))
        IF (LS.NE.L+1) TEST = TEST + ABS(E(LS-1))
        ZTEST = TEST + ABS(S(LS))
        IF (ZTEST.NE.TEST) GO TO 420
        S(LS) = 0.0E0
        GO TO 440
  420   CONTINUE
  430 CONTINUE
  440 CONTINUE
      IF (LS.NE.L) GO TO 450
      KASE = 3
      GO TO 470
  450 CONTINUE
      IF (LS.NE.M) GO TO 460
      KASE = 1
      GO TO 470
  460 CONTINUE
      KASE = 2
      L = LS
  470 CONTINUE
  480 CONTINUE
      L = L + 1
c
c  PERFORM THE TASK INDICATED BY KASE.
c
      GO TO (490, 520, 540, 570), KASE
c
c  DEFLATE NEGLIGIBLE S(M).
c
  490 CONTINUE
      MM1 = M - 1
      F = E(M-1)
      E(M-1) = 0.0E0
      DO 510 KK=L,MM1
        K = MM1 - KK + L
        T1 = S(K)
        CALL SROTG(T1, F, CS, SN)
        S(K) = T1
        IF (K.EQ.L) GO TO 500
        F = -SN*E(K-1)
        E(K-1) = CS*E(K-1)
  500   CONTINUE
        IF (WANTV) CALL SROT(P, V(1,K), 1, V(1,M), 1, CS, SN)
  510 CONTINUE
      GO TO 610
c
c  SPLIT AT NEGLIGIBLE S(L).
c
  520 CONTINUE
      F = E(L-1)
      E(L-1) = 0.0E0

      DO K=L,M
        T1 = S(K)
        CALL SROTG(T1, F, CS, SN)
        S(K) = T1
        F = -SN*E(K)
        E(K) = CS*E(K)
        IF (WANTU) CALL SROT(N, U(1,K), 1, U(1,L-1), 1, CS, SN)
      end do

      GO TO 610
c
c  PERFORM ONE QR STEP.
c
  540 CONTINUE
c
c  CALCULATE THE SHIFT.
c
      SCALE = AMAX1(ABS(S(M)),ABS(S(M-1)),ABS(E(M-1)),ABS(S(L)),ABS(E(L)
     & ))
      SM = S(M)/SCALE
      SMM1 = S(M-1)/SCALE
      EMM1 = E(M-1)/SCALE
      SL = S(L)/SCALE
      EL = E(L)/SCALE
      B = ((SMM1+SM)*(SMM1-SM)+EMM1**2)/2.0E0
      C = (SM*EMM1)**2
      SHIFT = 0.0E0
      IF (B.EQ.0.0E0 .AND. C.EQ.0.0E0) GO TO 550
      SHIFT = SQRT(B**2+C)
      IF (B.LT.0.0E0) SHIFT = -SHIFT
      SHIFT = C/(B+SHIFT)
  550 CONTINUE
      F = (SL+SM)*(SL-SM) - SHIFT
      G = SL*EL
c
c  CHASE ZEROS.
c
      MM1 = M - 1
      DO 560 K=L,MM1
        CALL SROTG(F, G, CS, SN)
        IF (K.NE.L) E(K-1) = F
        F = CS*S(K) + SN*E(K)
        E(K) = CS*E(K) - SN*S(K)
        G = SN*S(K+1)
        S(K+1) = CS*S(K+1)
        IF (WANTV) CALL SROT(P, V(1,K), 1, V(1,K+1), 1, CS, SN)
        CALL SROTG(F, G, CS, SN)
        S(K) = F
        F = CS*E(K) + SN*S(K+1)
        S(K+1) = -SN*E(K) + CS*S(K+1)
        G = SN*E(K+1)
        E(K+1) = CS*E(K+1)
        IF (WANTU .AND. K.LT.N) CALL SROT(N, U(1,K), 1, U(1,K+1), 1,
     &    CS, SN)
  560 CONTINUE
      E(M-1) = F
      ITER = ITER + 1
      GO TO 610
c
c  CONVERGENCE.
c
  570 CONTINUE
c
c  MAKE THE SINGULAR VALUE  POSITIVE.
c
      IF (S(L).GE.0.0E0) GO TO 580
      S(L) = -S(L)
      IF (WANTV) CALL SSCAL(P, -1.0E0, V(1,L), 1)
  580 CONTINUE
c
c  ORDER THE SINGULAR VALUE.
c
  590 IF (L.EQ.MM) GO TO 600
      IF (S(L).GE.S(L+1)) GO TO 600
      T = S(L)
      S(L) = S(L+1)
      S(L+1) = T
      IF (WANTV .AND. L.LT.P) CALL SSWAP(P, V(1,L), 1, V(1,L+1), 1)
      IF (WANTU .AND. L.LT.N) CALL SSWAP(N, U(1,L), 1, U(1,L+1), 1)
      L = L + 1
      GO TO 590
  600 CONTINUE
      ITER = 0
      M = M - 1
  610 CONTINUE
      GO TO 360
  620 CONTINUE
      RETURN
      END
      subroutine sswap ( n, sx, incx, sy, incy )

c*********************************************************************72
c
cc SSWAP interchanges two vectors.
c
c  Discussion:
c
c    This routine uses single precision real arithmetic.
c
c    This routine uses unrolled loops for increments equal to 1.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vectors.
c
c    Input/output, real X(*), one of the vectors to swap.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Input/output, real Y(*), one of the vectors to swap.
c
c    Input, integer INCY, the increment between successive elements of Y.
c
      implicit none

      real sx(*),sy(*),stemp
      integer i,incx,incy,ix,iy,m,n

      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c  code for unequal increments or equal increments not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        stemp = sx(ix)
        sx(ix) = sy(iy)
        sy(iy) = stemp
        ix = ix + incx
        iy = iy + incy
      end do
      return
c
c  code for both increments equal to 1
c
c
c  clean-up loop
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        stemp = sx(i)
        sx(i) = sy(i)
        sy(i) = stemp
      end do
      if( n .lt. 3 ) return
   40 continue

      do i = m+1, n, 3
        stemp = sx(i)
        sx(i) = sy(i)
        sy(i) = stemp
        stemp = sx(i + 1)
        sx(i + 1) = sy(i + 1)
        sy(i + 1) = stemp
        stemp = sx(i + 2)
        sx(i + 2) = sy(i + 2)
        sy(i + 2) = stemp
      end do

      return
      end
      subroutine svd ( nm, m, n, a, w, matu, u, matv, v, ierr, rv1 )

c*********************************************************************72
c
cc SVD computes the singular value decomposition of a real rectangular matrix.
c
c  Discussion:
c
c     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE SVD,
c     NUM. MATH. 14, 403-420(1970) BY GOLUB AND REINSCH.
c     HANDBOOK FOR AUTO. COMP., VOL II-LINEAR ALGEBRA, 134-151(1971).
c
c     THIS SUBROUTINE DETERMINES THE SINGULAR VALUE DECOMPOSITION
c          T
c     A=USV  OF A REAL M BY N RECTANGULAR MATRIX.  HOUSEHOLDER
c     BIDIAGONALIZATION AND A VARIANT OF THE QR ALGORITHM ARE USED.
c
c     ON INPUT-
c
c        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
c          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
c          DIMENSION STATEMENT.  NOTE THAT NM MUST BE AT LEAST
c          AS LARGE AS THE MAXIMUM OF M AND N,
c
c        M IS THE NUMBER OF ROWS OF A (AND U),
c
c        N IS THE NUMBER OF COLUMNS OF A (AND U) AND THE ORDER OF V,
c
c        A CONTAINS THE RECTANGULAR INPUT MATRIX TO BE DECOMPOSED,
c
c        MATU SHOULD BE SET TO .TRUE. IF THE U MATRIX IN THE
c          DECOMPOSITION IS DESIRED, AND TO .FALSE. OTHERWISE,
c
c        MATV SHOULD BE SET TO .TRUE. IF THE V MATRIX IN THE
c          DECOMPOSITION IS DESIRED, AND TO .FALSE. OTHERWISE.
c
c     ON OUTPUT-
c
c        A IS UNALTERED (UNLESS OVERWRITTEN BY U OR V),
c
c        W CONTAINS THE N (NON-NEGATIVE) SINGULAR VALUES OF A (THE
c          DIAGONAL ELEMENTS OF S).  THEY ARE UNORDERED.  IF AN
c          ERROR EXIT IS MADE, THE SINGULAR VALUES SHOULD BE CORRECT
c          FOR INDICES IERR+1,IERR+2,...,N,
c
c        U CONTAINS THE MATRIX U (ORTHOGONAL COLUMN VECTORS) OF THE
c          DECOMPOSITION IF MATU HAS BEEN SET TO .TRUE.  OTHERWISE
c          U IS USED AS A TEMPORARY ARRAY.  U MAY COINCIDE WITH A.
c          IF AN ERROR EXIT IS MADE, THE COLUMNS OF U CORRESPONDING
c          TO INDICES OF CORRECT SINGULAR VALUES SHOULD BE CORRECT,
c
c        V CONTAINS THE MATRIX V (ORTHOGONAL) OF THE DECOMPOSITION IF
c          MATV HAS BEEN SET TO .TRUE.  OTHERWISE V IS NOT REFERENCED.
c          V MAY ALSO COINCIDE WITH A IF U IS NOT NEEDED.  IF AN ERROR
c          EXIT IS MADE, THE COLUMNS OF V CORRESPONDING TO INDICES OF
c          CORRECT SINGULAR VALUES SHOULD BE CORRECT,
c
c        IERR IS SET TO
c          ZERO       FOR NORMAL RETURN,
c          K          IF THE K-TH SINGULAR VALUE HAS NOT BEEN
c                     DETERMINED AFTER 30 ITERATIONS,
c
c        RV1 IS A TEMPORARY STORAGE ARRAY.
c
c     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
c     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
c
      implicit none

      INTEGER I, J, K, L, M, N, II, I1, KK, K1, LL, L1, MN, NM, ITS,
     & IERR
      REAL A(NM,N), W(N), U(NM,N), V(NM,N), RV1(N)
      REAL C, F, G, H, S, X, Y, Z, EPS, SCALE, MACHEP
      REAL SQRT, AMAX1, ABS, SIGN
      LOGICAL MATU, MATV
c
c  MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING
c  THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.
c
      MACHEP = 2.**(-26)

      IERR = 0

      DO I=1,M
        DO J=1,N
          U(I,J) = A(I,J)
        end do
      end do
c
c  HOUSEHOLDER REDUCTION TO BIDIAGONAL FORM.
c
      G = 0.0
      SCALE = 0.0
      X = 0.0

      DO 200 I=1,N

        L = I + 1
        RV1(I) = SCALE*G
        G = 0.0
        S = 0.0
        SCALE = 0.0
        IF (I.GT.M) GO TO 100

        DO K=I,M
          SCALE = SCALE + ABS(U(K,I))
        end do

        IF (SCALE.EQ.0.0) GO TO 100

        DO K=I,M
          U(K,I) = U(K,I)/SCALE
          S = S + U(K,I)**2
        end do

        F = U(I,I)
        G = -SIGN(SQRT(S),F)
        H = F*G - S
        U(I,I) = F - G

        IF ( I .NE. N ) then

          DO J = L, N

            S = 0.0
            DO K = I, M
              S = S + U(K,I) * U(K,J)
            end do

            F = S / H

            DO K = I, M
              U(K,J) = U(K,J) + F * U(K,I)
            end do

          end do

        end if

        DO K=I,M
          U(K,I) = SCALE*U(K,I)
        end do

  100   continue

        W(I) = SCALE*G
        G = 0.0
        S = 0.0
        SCALE = 0.0

        IF (I.GT.M .OR. I.EQ.N) GO TO 190

        DO K=L,N
          SCALE = SCALE + ABS(U(I,K))
        end do

        IF (SCALE.EQ.0.0) GO TO 190

        DO K=L,N
          U(I,K) = U(I,K)/SCALE
          S = S + U(I,K)**2
        end do

        F = U(I,L)
        G = -SIGN(SQRT(S),F)
        H = F*G - S
        U(I,L) = F - G

        DO K=L,N
          RV1(K) = U(I,K)/H
        end do

        IF (I.EQ.M) GO TO 170

        DO J = L, M

          S = 0.0
          DO K=L,N
            S = S + U(J,K)*U(I,K)
          end do

          DO K=L,N
            U(J,K) = U(J,K) + S*RV1(K)
          end do

        end do

  170   DO K=L,N
          U(I,K) = SCALE*U(I,K)
        end do

  190   X = AMAX1(X,ABS(W(I))+ABS(RV1(I)))

  200 CONTINUE
c
c  ACCUMULATION OF RIGHT-HAND TRANSFORMATIONS.
c
      IF (.NOT.MATV) GO TO 290
c
c  FOR I=N STEP -1 UNTIL 1 DO -- .
c
      DO 280 II=1,N

        I = N + 1 - II

        IF (I.EQ.N) GO TO 270

        IF (G.ne.0.0) then

          DO J=L,N
            V(J,I) = (U(I,J)/U(I,L))/G
          end do

          DO J=L,N

            S = 0.0
            DO K=L,N
              S = S + U(I,K)*V(K,J)
            end do

            DO K=L,N
              V(K,J) = V(K,J) + S*V(K,I)
            end do

          end do

        end if

        DO J=L,N
          V(I,J) = 0.0
          V(J,I) = 0.0
        end do

  270   V(I,I) = 1.0
        G = RV1(I)
        L = I
  280 CONTINUE
c
c  ACCUMULATION OF LEFT-HAND TRANSFORMATIONS.
c
  290 IF (.NOT.MATU) GO TO 410
c
c  FOR I=MIN(M,N) STEP -1 UNTIL 1 DO --.
c
      MN = N
      IF (M.LT.N) MN = M

      DO 400 II=1,MN
        I = MN + 1 - II
        L = I + 1
        G = W(I)
        IF (I.EQ.N) GO TO 310

        DO J=L,N
          U(I,J) = 0.0
        end do

  310   IF (G.EQ.0.0) GO TO 370
        IF (I.EQ.MN) GO TO 350

        DO J=L,N

          S = 0.0
          DO K=L,M
            S = S + U(K,I)*U(K,J)
          end do

          F = (S/U(I,I))/G
          DO K=I,M
            U(K,J) = U(K,J) + F*U(K,I)
          end do

        end do

  350   DO J=I,M
          U(J,I) = U(J,I)/G
        end do

        GO TO 390

  370   DO J=I,M
          U(J,I) = 0.0
        end do

  390   U(I,I) = U(I,I) + 1.0
  400 CONTINUE
c
c  DIAGONALIZATION OF THE BIDIAGONAL FORM.
c
  410 EPS = MACHEP*X
c
c  FOR K=N STEP -1 UNTIL 1 DO --.
c
      DO 550 KK=1,N
        K1 = N - KK
        K = K1 + 1
        ITS = 0
c
c  TEST FOR SPLITTING.
c  FOR L=K STEP -1 UNTIL 1 DO -- .
c
  420   DO 430 LL=1,K
          L1 = K - LL
          L = L1 + 1
          IF (ABS(RV1(L)).LE.EPS) GO TO 470
c
c  RV1(1) IS ALWAYS ZERO, SO THERE IS NO EXIT
c  THROUGH THE BOTTOM OF THE LOOP.
c
          IF (ABS(W(L1)).LE.EPS) GO TO 440
  430   CONTINUE
c
c  CANCELLATION OF RV1(L) IF L GREATER THAN 1.
c
  440   C = 0.0
        S = 1.0

        DO I=L,K

          F = S*RV1(I)
          RV1(I) = C*RV1(I)
          IF (ABS(F).LE.EPS) GO TO 470
          G = W(I)
          H = SQRT(F*F+G*G)
          W(I) = H
          C = G/H
          S = -F/H

          IF ( MATU ) then

            DO J=1,M
              Y = U(J,L1)
              Z = U(J,I)
              U(J,L1) = Y*C + Z*S
              U(J,I) = -Y*S + Z*C
            end do

          end if

       end do

470     continue
c
c  TEST FOR CONVERGENCE.
c
        Z = W(K)
        IF (L.EQ.K) GO TO 530
c
c  SHIFT FROM BOTTOM 2 BY 2 MINOR.
c
        IF (ITS.EQ.30) GO TO 560
        ITS = ITS + 1
        X = W(L)
        Y = W(K1)
        G = RV1(K1)
        H = RV1(K)
        F = ((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0*H*Y)
        G = SQRT(F*F+1.0)
        F = ((X-Z)*(X+Z)+H*(Y/(F+SIGN(G,F))-H))/X
c
c  NEXT QR TRANSFORMATION.
c
        C = 1.0
        S = 1.0

        DO I1=L,K1

          I = I1 + 1
          G = RV1(I)
          Y = W(I)
          H = S*G
          G = C*G
          Z = SQRT(F*F+H*H)
          RV1(I1) = Z
          C = F/Z
          S = H/Z
          F = X*C + G*S
          G = -X*S + G*C
          H = Y*S
          Y = Y*C

          IF ( MATV ) then

            DO J=1,N
              X = V(J,I1)
              Z = V(J,I)
              V(J,I1) = X*C + Z*S
              V(J,I) = -X*S + Z*C
            end do

          end if

          Z = SQRT(F*F+H*H)
          W(I1) = Z
c
c  ROTATION CAN BE ARBITRARY IF Z IS ZERO.
c
          IF ( Z .ne. 0.0 ) then
            C = F/Z
            S = H/Z
          end if

          F = C*G + S*Y
          X = -S*G + C*Y

          IF ( MATU ) then

            DO  J=1,M
              Y = U(J,I1)
              Z = U(J,I)
              U(J,I1) = Y*C + Z*S
              U(J,I) = -Y*S + Z*C
            end do

          end if

        end do

        RV1(L) = 0.0
        RV1(K) = F
        W(K) = X
        GO TO 420
c
c  CONVERGENCE.
c
  530   IF (Z.GE.0.0) GO TO 550
c
c  W(K) IS MADE NON-NEGATIVE.
c
        W(K) = -Z
        IF (.NOT.MATV) GO TO 550

        DO J=1,N
          V(J,K) = -V(J,K)
        end do

  550 CONTINUE

      GO TO 570
c
c  SET ERROR -- NO CONVERGENCE TO A SINGULAR VALUE AFTER 30 ITERATIONS.
c
  560 IERR = K
  570 RETURN
      END
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    None
c
      implicit none

      character * ( 8 ) ampm
      integer d
      character * ( 8 ) date
      integer h
      integer m
      integer mm
      character * ( 9 ) month(12)
      integer n
      integer s
      character * ( 10 ) time
      integer y

      save month

      data month /
     &  'January  ', 'February ', 'March    ', 'April    ', 
     &  'May      ', 'June     ', 'July     ', 'August   ', 
     &  'September', 'October  ', 'November ', 'December ' /

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) y, m, d
      read ( time, '(i2,i2,i2,1x,i3)' ) h, n, s, mm

      if ( h .lt. 12 ) then
        ampm = 'AM'
      else if ( h .eq. 12 ) then
        if ( n .eq. 0 .and. s .eq. 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h .lt. 12 ) then
          ampm = 'PM'
        else if ( h .eq. 12 ) then
          if ( n .eq. 0 .and. s .eq. 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if

      write ( *, 
     &  '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) 
     &  d, month(m), y, h, ':', n, ':', s, '.', mm, ampm

      return
      end
