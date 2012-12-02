      subroutine bndacc ( g, mdg, nb, ip, ir, mt, jt )  

c*********************************************************************72
c
cc BNDACC performs the accumulation phase of a banded least squares solver.
c
c  Discussion:
c
c    SEQUENTIAL ALGORITHM FOR BANDED LEAST SQUARES PROBLEM..  
c    ACCUMULATION PHASE.      FOR SOLUTION PHASE USE BNDSOL.  
c
c    THE CALLING PROGRAM MUST SET IR=1 AND IP=1 BEFORE THE FIRST CALL  
c    TO BNDACC FOR A NEW CASE.     
c   
c    THE SECOND SUBSCRIPT OF G( ) MUST BE DIMENSIONED AT LEAST 
c    NB+1 IN THE CALLING PROGRAM.  
c
c    BNDACC is to be called once for each block of data [ C(I), B(I) ]
c    to be introduced into the problem.  For each block of data, the
c    user must assign values to MT and JT, and copy the MT by (NB+1)
c    array of data [ C(I), B(I) ] into rows IR through IR+MT-1 of
c    the working array G.
c
c  Modified:
c
c    19 October 2008
c
c  Author:
c
c    Charles Lawson, Richard Hanson
c
c  Reference:
c
c    Charles Lawson, Richard Hanson,
c    Solving Least Squares Problems,
c    SIAM, 1995,
c    ISBN: 0898713560,
c    LC: QA275.L38.
c
c  Parameters:
c
c    Input/output, double precision G(MDG,NB+1), the banded matrix being
c    accumulated.
c
c    Input, integer MDG, the leading dimension of G.
c    MDG must be at least equal to the number of rows.
c
c    Input, integer NB, the bandwidth of G, not counting the diagonal.
c
c    Input/output, integer IP, must be set to 1 by the user before
c    the first call.  Thereafter, its value is controlled by the routine.
c
c    Input/output, integer IR, indicates the index of the first row
c    of G into which new data is to be placed.  The user sets IR to
c    1 before the first call, but does not alter it thereafter.
c    Instead, the program keeps this quantity up to date as more
c    data is read in.
c
c    Input, integer MT, is set by the user to indicate the number of
c    new rows of data being introduced by the current call.
c    MT must be at least 0.
c
c    Input, integer JT, set by the user to indicate the column of the
c    submatrix A(I) that is identified with the first column of C(I).
c    This means that JT must be at least 1.
c
      integer i, j, ie, ig, ig1, ig2, ip, ir, jg, jt, k, kh, l, lp1
      integer mdg, mh, mt, mu, nb, nbp1

      double precision g(mdg,nb+1)
      double precision rho
  
      NBP1 = NB + 1
 
      if ( mt .le. 0 ) then
        return
      end if

      if (jt.eq.ip) go to 70

      if (jt.le.ir) go to 30

      do i=1,mt  
        ig1=jt+mt-i     
        ig2=ir+mt-i     
        do j=1,nbp1  
          g(ig1,j)=g(ig2,j) 
        end do
      end do  
 
      ie=jt-ir  
      do i=1,ie  
        ig=ir+i-1   
        do j=1,nbp1  
          g(ig,j)= 0.0D+00
        end do 
      end do  
 
      ir=jt 

   30 mu=min(nb-1,ir-ip-1) 

      do l=1,mu  

        k=min(l,jt-ip) 

        lp1=l+1 
        ig=ip+l 
        do i=lp1,nb  
          jg=i-k
          g(ig,jg)=g(ig,i)  
        end do

        do i=1,k     
          jg=nbp1-i   
          g(ig,jg) = 0.0D+00
        end do

      end do

      ip=jt 

   70 mh=ir+mt-ip   
      kh=min(nbp1,mh)  

      do i=1,kh  
        call h12 (1,i,max(i+1,ir-ip+1),mh,g(ip,i),1,rho,   
     *            g(ip,i+1),1,mdg,nbp1-i)  
      end do 

      ir = ip + kh  
 
      if ( nb + 1 .le. kh ) then
        do i = 1, nb  
          g(ir-1,i) = 0.0D+00
        end do
      end if

      return
      end   
      subroutine bndsol (mode,g,mdg,nb,ip,ir,x,n,rnorm)

c*********************************************************************72  
c
cc BNDSOL performs the solution phase of a banded least squares solver.
c
c  Discussion:
c
c    SEQUENTIAL SOLUTION OF A BANDED LEAST SQUARES PROBLEM..  
c    SOLUTION PHASE.   FOR THE ACCUMULATION PHASE USE BNDACC.     
c
c  Modified:
c
c    19 October 2008
c
c  Author:
c
c    Charles Lawson, Richard Hanson
c
c  Reference:
c
c    Charles Lawson, Richard Hanson,
c    Solving Least Squares Problems,
c    SIAM, 1995,
c    ISBN: 0898713560,
c    LC: QA275.L38.
c
c  Parameters:
c
c     MODE = 1     SOLVE R*X=Y   WHERE R AND Y ARE IN THE G( ) ARRAY    
c                  AND X WILL BE STORED IN THE X( ) ARRAY.  
c            2     SOLVE (R**T)*X=Y   WHERE R IS IN G( ),   
c                  Y IS INITIALLY IN X( ), AND X REPLACES Y IN X( ),    
c            3     SOLVE R*X=Y   WHERE R IS IN G( ).
c                  Y IS INITIALLY IN X( ), AND X REPLACES Y IN X( ).    
c   
c     THE SECOND SUBSCRIPT OF G( ) MUST BE DIMENSIONED AT LEAST 
c     NB+1 IN THE CALLING PROGRAM.  
c
      integer I, I1, I2, IE, II, IP, IR, IX, J, JG, L
      integer MDG, MODE, N, NB, NP1, IRM1
      double precision G(MDG,*), RNORM, RSQ, S, X(N), ZERO
      parameter(ZERO = 0.0d0)

      RNORM=ZERO
      GO TO (10,90,50), MODE
c
c  MODE = 1  
c                                   ALG. STEP 26
   10      DO 20 J=1,N  
   20      X(J)=G(J,NB+1)   
      RSQ=ZERO  
      NP1=N+1   
      IRM1=IR-1 
      IF (NP1.GT.IRM1) GO TO 40     
           DO 30 J=NP1,IRM1 
   30      RSQ=RSQ+G(J,NB+1)**2     
      RNORM=SQRT(RSQ)   
   40 CONTINUE
c
c  MODE = 3  
c                                   ALG. STEP 27
   50      DO 80 II=1,N 
           I=N+1-II     
c                                   ALG. STEP 28
           S=ZERO   
           L=max(0,I-IP)   
c                                   ALG. STEP 29
           IF (I.EQ.N) GO TO 70     
c                                   ALG. STEP 30
           IE=min(N+1-I,NB)
                DO 60 J=2,IE
                JG=J+L  
                IX=I-1+J
   60           S=S+G(I,JG)*X(IX)   
c                                   ALG. STEP 31
   70      continue
           IF (G(I,L+1) .eq. ZERO) go to 130
   80      X(I)=(X(I)-S)/G(I,L+1)   
c                                   ALG. STEP 32
      RETURN
c
c  MODE = 2
c
   90      DO 120 J=1,N 
           S=ZERO   
           IF (J.EQ.1) GO TO 110    
           I1=max(1,J-NB+1)
           I2=J-1   
                DO 100 I=I1,I2  
                L=J-I+1+max(0,I-IP)
  100           S=S+X(I)*G(I,L)     
  110      L=max(0,J-IP)   
           IF (G(J,L+1) .eq. ZERO) go to 130
  120      X(J)=(X(J)-S)/G(J,L+1)   
      RETURN

  130 write (*,'(/a/a,4i6)')' ZERO DIAGONAL TERM IN BNDSOL.',
     *      ' MODE,I,J,L = ',MODE,I,J,L  
      STOP  
      END
      function diff ( x, y )

c*********************************************************************72
c
cc DIFF is used in tests that depend on machine precision.
c
c  Modified:
c
c    19 October 2008
c
c  Author:
c
c    Charles Lawson, Richard Hanson
c
c  Reference:
c
c    Charles Lawson, Richard Hanson,
c    Solving Least Squares Problems,
c    SIAM, 1995,
c    ISBN: 0898713560,
c    LC: QA275.L38.
c
c  Parameters:
c
c    Input, double precision X, Y.
c
c    Output, double precision DIFF, the value of X - Y.
c
      implicit none

      double precision diff
      double precision x
      double precision y

      diff = x - y  

      return
      end   
      subroutine g1 ( a, b, cterm, sterm, sig )

c*********************************************************************72 
c
cc G1 computes an orthogonal rotation matrix.
c
c  Discussion:
c
c    COMPUTE.. MATRIX   (C, S) SO THAT (C, S)(A) = (SQRT(A**2+B**2))   
c                        (-S,C)         (-S,C)(B)   (   0          )    
c    COMPUTE SIG = SQRT(A**2+B**2) 
c       SIG IS COMPUTED LAST TO ALLOW FOR THE POSSIBILITY THAT 
c       SIG MAY BE IN THE SAME LOCATION AS A OR B .
c
c  Modified:
c
c    19 October 2008
c
c  Author:
c
c    Charles Lawson, Richard Hanson
c
c  Reference:
c
c    Charles Lawson, Richard Hanson,
c    Solving Least Squares Problems,
c    SIAM, 1995,
c    ISBN: 0898713560,
c    LC: QA275.L38.
c
c  Parameters:
c
      double precision a, b, cterm, one, sig, sterm, xr, yr
      parameter(one = 1.0d0)

      if (abs(a) .gt. abs(b)) then
         xr=b/a
         yr=sqrt(one+xr**2)
         cterm=sign(one/yr,a)
         sterm=cterm*xr
         sig=abs(a)*yr     
         return
      endif

      if (b .ne. 0.0D+00 ) then
         xr=a/b
         yr=sqrt(one+xr**2)
         sterm=sign(one/yr,b)
         cterm=sterm*xr
         sig=abs(b)*yr     
         return
      endif

      sig= 0.0D+00  
      cterm= 0.0D+00
      sterm=one   
      return
      end   
      subroutine g2 ( cterm, sterm, x, y )

c*********************************************************************72
c
cc G2 applies the rotation computed by G1 to the vector (X,Y).
c
c  Modified:
c
c    19 October 2008
c
c  Author:
c
c    Charles Lawson, Richard Hanson
c
c  Reference:
c
c    Charles Lawson, Richard Hanson,
c    Solving Least Squares Problems,
c    SIAM, 1995,
c    ISBN: 0898713560,
c    LC: QA275.L38.
c
c  Parameters:
c
      implicit none

      double precision cterm, sterm, x, xr, y

      xr=cterm*x+sterm*y    
      y=-sterm*x+cterm*y    
      x=xr

      return
      end 
      function gen ( anoise )

c*********************************************************************72
c
cc GEN generates numbers for construction of test cases.
c
c  Modified:
c
c    22 October 2008
c
c  Author:
c
c    Charles Lawson, Richard Hanson
c
c  Reference:
c
c    Charles Lawson, Richard Hanson,
c    Solving Least Squares Problems,
c    SIAM, 1995,
c    ISBN: 0898713560,
c    LC: QA275.L38.
c
c  Parameters:
c
c    Input, double precision ANOISE, determines the level of
c    "noise" to be added to the data.
c
c    Output, double precision GEN, a random value with noise added.
c
      implicit none

      double precision ai
      double precision aj
      double precision anoise
      double precision gen
      integer i
      integer j
      integer mi
      integer mj

      save

      if ( anoise .lt. 0.0D+00 ) then
        mi = 891
        mj = 457
        i = 5   
        j = 7   
        aj = 0.0D+00
        gen = 0.0D+00
        return
      end if
c   
c  The sequence of values of J is bounded between 1 and 996.
c  If initial J = 1,2,3,4,5,6,7,8, OR 9, the period is 332.
c
      if ( 0.0D+00 .lt. anoise ) then
        j = j * mj
        j = j - 997 * ( j / 997 )   
        aj = j - 498  
      end if
c
c  The sequence of values of I is bounded between 1 and 999.
c  If initial I = 1,2,3,6,7, or 9, the period will be 50;
c  If initial I = 4 or 8           the period will be 25;
c  If initial I = 5                the period will be 10.
c
      i = i * mi
      i = i - 1000 * ( i / 1000 ) 
      ai = i - 500  
      gen = ai + aj * anoise  

      return
      end
      subroutine h12 (mode,lpivot,l1,m,u,iue,up,c,ice,icv,ncv)  

c*********************************************************************72
c
cc H12 constructs and applies a Householder transformation.
c
c  Discussion:
c
c    The transformation has the form Q = I + U*(U**T)/B.
c
c  Modified:
c
c    19 October 2008
c
c  Author:
c
c    Charles Lawson, Richard Hanson
c
c  Reference:
c
c    Charles Lawson, Richard Hanson,
c    Solving Least Squares Problems,
c    SIAM, 1995,
c    ISBN: 0898713560,
c    LC: QA275.L38.
c
c  Parameters:
c
c     MODE   = 1 OR 2   Selects Algorithm H1 to construct and apply a
c            Householder transformation, or Algorithm H2 to apply a
c            previously constructed transformation.
c
c     LPIVOT IS THE INDEX OF THE PIVOT ELEMENT. 
c
c     L1,M   IF L1 .LE. M   THE TRANSFORMATION WILL BE CONSTRUCTED TO   
c            ZERO ELEMENTS INDEXED FROM L1 THROUGH M.   IF L1 GT. M     
c            THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION.
c
c     U(),IUE,UP    On entry with MODE = 1, U() contains the pivot
c            vector.  IUE is the storage increment between elements.  
c            On exit when MODE = 1, U() and UP contain quantities
c            defining the vector U of the Householder transformation.
c            on entry with MODE = 2, U() and UP should contain
c            quantities previously computed with MODE = 1.  These will
c            not be modified during the entry with MODE = 2.   
c
c     C()    ON ENTRY with MODE = 1 or 2, C() CONTAINS A MATRIX WHICH
c            WILL BE REGARDED AS A SET OF VECTORS TO WHICH THE
c            HOUSEHOLDER TRANSFORMATION IS TO BE APPLIED.
c            ON EXIT C() CONTAINS THE SET OF TRANSFORMED VECTORS.
c
c     ICE    STORAGE INCREMENT BETWEEN ELEMENTS OF VECTORS IN C().  
c
c     ICV    STORAGE INCREMENT BETWEEN VECTORS IN C().  
c
c     NCV    NUMBER OF VECTORS IN C() TO BE TRANSFORMED. IF NCV .LE. 0  
c            NO OPERATIONS WILL BE DONE ON C(). 
c
      integer I, I2, I3, I4, ICE, ICV, INCR, IUE, J
      integer L1, LPIVOT, M, MODE, NCV
      double precision B, C(*), CL, CLINV, ONE, SM
      double precision U(IUE,M)
      double precision UP
      parameter(ONE = 1.0d0)

      IF (0.GE.LPIVOT.OR.LPIVOT.GE.L1.OR.L1.GT.M) RETURN    
      CL=abs(U(1,LPIVOT))   
      IF (MODE.EQ.2) GO TO 60   
c
c  CONSTRUCT THE TRANSFORMATION.
c
          DO 10 J=L1,M  
   10     CL=MAX(abs(U(1,J)),CL)  
      IF (CL) 130,130,20
   20 CLINV=ONE/CL  
      SM=(U(1,LPIVOT)*CLINV)**2   
          DO 30 J=L1,M  
   30     SM=SM+(U(1,J)*CLINV)**2 
      CL=CL*SQRT(SM)   
      IF (U(1,LPIVOT)) 50,50,40     
   40 CL=-CL
   50 UP=U(1,LPIVOT)-CL 
      U(1,LPIVOT)=CL    
      GO TO 70
c 
c  APPLY THE TRANSFORMATION  I+U*(U**T)/B  TO C.
c   
   60 IF (CL) 130,130,70
   70 IF (NCV.LE.0) RETURN  
      B= UP*U(1,LPIVOT)
c
c  B  MUST BE NONPOSITIVE HERE.  IF B = 0., RETURN.
c   
      IF (B) 80,130,130 
   80 B=ONE/B   
      I2=1-ICV+ICE*(LPIVOT-1)   
      INCR=ICE*(L1-LPIVOT)  
          DO 120 J=1,NCV
          I2=I2+ICV     
          I3=I2+INCR    
          I4=I3 
          SM=C(I2)*UP
              DO 90 I=L1,M  
              SM=SM+C(I3)*U(1,I)
   90         I3=I3+ICE 
          IF (SM) 100,120,100   
  100     SM=SM*B   
          C(I2)=C(I2)+SM*UP
              DO 110 I=L1,M 
              C(I4)=C(I4)+SM*U(1,I)
  110         I4=I4+ICE 
  120     CONTINUE  
  130 RETURN
      END   
      SUBROUTINE HFTI (A,MDA,M,N,B,MDB,NB,TAU,KRANK,RNORM,H,G,IP)  

c*********************************************************************72 
c
cc HFTI solves least squares problems.
c
c  Discussion:
c
c    The routine uses the HFTI algorithm, that is,
c    "Householder Forward Triangulation with column Interchanges".
c
c  Modified:
c
c    19 October 2008
c
c  Author:
c
c    Charles Lawson, Richard Hanson
c
c  Reference:
c
c    Charles Lawson, Richard Hanson,
c    Solving Least Squares Problems,
c    SIAM, 1995,
c    ISBN: 0898713560,
c    LC: QA275.L38.
c
c  Parameters:
c
      integer I, II, IP1, J, JB, JJ, K, KP1, KRANK
      integer L, LDIAG, LMAX, M, MDA, MDB, N, NB
      integer IP(N)     
      double precision A(MDA,N),B(MDB,NB),H(N),G(N),RNORM(NB)  
      double precision DIFF, FACTOR, HMAX, SM, TAU, TMP, ZERO     
      parameter(FACTOR = 0.001d0, ZERO = 0.0d0)

      K=0   
      LDIAG=min(M,N)   
      IF (LDIAG.LE.0) GO TO 270     
          DO 80 J=1,LDIAG   
          IF (J.EQ.1) GO TO 20  
c   
c  UPDATE SQUARED COLUMN LENGTHS AND FIND LMAX   
c    
          LMAX=J
              DO 10 L=J,N   
              H(L)=H(L)-A(J-1,L)**2 
              IF (H(L).GT.H(LMAX)) LMAX=L   
   10         CONTINUE  
          IF(DIFF(HMAX+FACTOR*H(LMAX),HMAX)) 20,20,50   
c   
c  COMPUTE SQUARED COLUMN LENGTHS AND FIND LMAX  
c
   20     LMAX=J
              DO 40 L=J,N   
              H(L)=0.   
                  DO 30 I=J,M   
   30             H(L)=H(L)+A(I,L)**2   
              IF (H(L).GT.H(LMAX)) LMAX=L   
   40         CONTINUE  
          HMAX=H(LMAX)  
c
c  LMAX HAS BEEN DETERMINED  
c   
c  DO COLUMN INTERCHANGES IF NEEDED. 
c
   50     CONTINUE  
          IP(J)=LMAX    
          IF (IP(J).EQ.J) GO TO 70  
              DO 60 I=1,M   
              TMP=A(I,J)
              A(I,J)=A(I,LMAX)  
   60         A(I,LMAX)=TMP 
          H(LMAX)=H(J)  
c   
c     COMPUTE THE J-TH TRANSFORMATION AND APPLY IT TO A AND B.  
c
   70     CALL H12 (1,J,J+1,M,A(1,J),1,H(J),A(1,J+1),1,MDA,N-J) 
   80     CALL H12 (2,J,J+1,M,A(1,J),1,H(J),B,1,MDB,NB)     
c   
c     DETERMINE THE PSEUDORANK, K, USING THE TOLERANCE, TAU.
c
          DO 90 J=1,LDIAG   
          IF (ABS(A(J,J)).LE.TAU) GO TO 100     
   90     CONTINUE  
      K=LDIAG   
      GO TO 110 
  100 K=J-1 
  110 KP1=K+1   
c   
c     COMPUTE THE NORMS OF THE RESIDUAL VECTORS.
c   
      IF (NB.LE.0) GO TO 140
          DO 130 JB=1,NB
          TMP=ZERO     
          IF (KP1.GT.M) GO TO 130   
              DO 120 I=KP1,M
  120         TMP=TMP+B(I,JB)**2    
  130     RNORM(JB)=SQRT(TMP)   
  140 CONTINUE 
c 
c  SPECIAL FOR PSEUDORANK = 0  
c
      IF (K.GT.0) GO TO 160 
      IF (NB.LE.0) GO TO 270
          DO 150 JB=1,NB
              DO 150 I=1,N  
  150         B(I,JB)=ZERO 
      GO TO 270 
c   
c     IF THE PSEUDORANK IS LESS THAN N COMPUTE HOUSEHOLDER  
c     DECOMPOSITION OF FIRST K ROWS.
c
  160 IF (K.EQ.N) GO TO 180 
          DO 170 II=1,K 
          I=KP1-II  
  170     CALL H12 (1,I,KP1,N,A(I,1),MDA,G(I),A,MDA,1,I-1)  
  180 CONTINUE  
c   
c   
      IF (NB.LE.0) GO TO 270
          DO 260 JB=1,NB
c   
c     SOLVE THE K BY K TRIANGULAR SYSTEM.   
c
              DO 210 L=1,K  
              SM=ZERO  
              I=KP1-L   
              IF (I.EQ.K) GO TO 200 
              IP1=I+1   
                  DO 190 J=IP1,K    
  190             SM=SM+A(I,J)*B(J,JB)
  200         continue
  210         B(I,JB)=(B(I,JB)-SM)/A(I,I)  
c   
c     COMPLETE COMPUTATION OF SOLUTION VECTOR.  
c
          IF (K.EQ.N) GO TO 240     
              DO 220 J=KP1,N
  220         B(J,JB)=ZERO 
              DO 230 I=1,K  
  230         CALL H12 (2,I,KP1,N,A(I,1),MDA,G(I),B(1,JB),1,MDB,1)  
c   
c      RE-ORDER THE SOLUTION VECTOR TO COMPENSATE FOR THE   
c      COLUMN INTERCHANGES. 
c
  240         DO 250 JJ=1,LDIAG     
              J=LDIAG+1-JJ  
              IF (IP(J).EQ.J) GO TO 250 
              L=IP(J)   
              TMP=B(L,JB)   
              B(L,JB)=B(J,JB)   
              B(J,JB)=TMP   
  250         CONTINUE  
  260     CONTINUE  
c
c     THE SOLUTION VECTORS, X, ARE NOW  
c     IN THE FIRST  N  ROWS OF THE ARRAY B(,).  
c   
  270 KRANK=K   
      RETURN
      END  
      SUBROUTINE LDP (G,MDG,M,N,H,X,XNORM,W,INDEX,MODE)

c*********************************************************************72     
c
cc LDP performs least distance programming.
c
c  Modified:
c
c    19 October 2008
c
c  Author:
c
c    Charles Lawson, Richard Hanson
c
c  Reference:
c
c    Charles Lawson, Richard Hanson,
c    Solving Least Squares Problems,
c    SIAM, 1995,
c    ISBN: 0898713560,
c    LC: QA275.L38.
c
c  Parameters:
c
      integer I, IW, IWDUAL, IY, IZ, J, JF, M, MDG, MODE, N, NP1
      integer INDEX(M)  
      double precision G(MDG,N), H(M), X(N), W(*)  
      double precision DIFF, FAC, ONE, RNORM, XNORM, ZERO
      parameter(ONE = 1.0d0, ZERO = 0.0d0)

      IF (N.LE.0) GO TO 120 
          DO 10 J=1,N   
   10     X(J)=ZERO     
      XNORM=ZERO
      IF (M.LE.0) GO TO 110 
c   
c     THE DECLARED DIMENSION OF W() MUST BE AT LEAST (N+1)*(M+2)+2*M.   
c   
c      FIRST (N+1)*M LOCS OF W()   =  MATRIX E FOR PROBLEM NNLS.
c       NEXT     N+1 LOCS OF W()   =  VECTOR F FOR PROBLEM NNLS.
c       NEXT     N+1 LOCS OF W()   =  VECTOR Z FOR PROBLEM NNLS.
c       NEXT       M LOCS OF W()   =  VECTOR Y FOR PROBLEM NNLS.
c       NEXT       M LOCS OF W()   =  VECTOR WDUAL FOR PROBLEM NNLS.    
c     COPY G**T INTO FIRST N ROWS AND M COLUMNS OF E.   
c     COPY H**T INTO ROW N+1 OF E.  
c   
      IW=0  
          DO 30 J=1,M   
              DO 20 I=1,N   
              IW=IW+1   
   20         W(IW)=G(J,I)  
          IW=IW+1   
   30     W(IW)=H(J)    
      JF=IW+1   
c
c  STORE N ZEROS FOLLOWED BY A ONE INTO F.
c
          DO 40 I=1,N   
          IW=IW+1   
   40     W(IW)=ZERO    
      W(IW+1)=ONE   

      NP1=N+1   
      IZ=IW+2   
      IY=IZ+NP1 
      IWDUAL=IY+M   

      CALL NNLS (W,NP1,NP1,M,W(JF),W(IY),RNORM,W(IWDUAL),W(IZ),INDEX,   
     *           MODE)  
c
c                      USE THE FOLLOWING RETURN IF UNSUCCESSFUL IN NNLS.
c
      IF (MODE.NE.1) RETURN 
      IF (RNORM) 130,130,50 
   50 FAC=ONE   
      IW=IY-1   
          DO 60 I=1,M   
          IW=IW+1   
c
c  HERE WE ARE USING THE SOLUTION VECTOR Y.
c
   60     FAC=FAC-H(I)*W(IW)

      IF (DIFF(ONE+FAC,ONE)) 130,130,70 
   70 FAC=ONE/FAC   
          DO 90 J=1,N   
          IW=IY-1   
              DO 80 I=1,M   
              IW=IW+1   
c
c  HERE WE ARE USING THE SOLUTION VECTOR Y.
c
   80         X(J)=X(J)+G(I,J)*W(IW)
   90     X(J)=X(J)*FAC 
          DO 100 J=1,N  
  100     XNORM=XNORM+X(J)**2   
      XNORM=sqrt(XNORM) 
c
c  SUCCESSFUL RETURN.
c
  110 MODE=1
      RETURN
c                             ERROR RETURN.       N .LE. 0. 
  120 MODE=2
      RETURN
c                             RETURNING WITH CONSTRAINTS NOT COMPATIBLE.
  130 MODE=4
      RETURN
      END   
      subroutine MFEOUT (A, MDA, M, N, NAMES, MODE, UNIT, WIDTH)

c*********************************************************************72
c
cc MFEOUT: LABELED MATRIX OUTPUT FOR USE WITH SINGULAR VALUE ANALYSIS.
c
c  Discussion:
c
c    This 1995 version has additional arguments, UNIT and WIDTH,
c    to support user options regarding the output unit and the width of
c    print lines.  Also allows user to choose length of names in NAMES().
c
c  Modified:
c
c    19 October 2008
c
c  Author:
c
c    Charles Lawson, Richard Hanson
c
c  Reference:
c
c    Charles Lawson, Richard Hanson,
c    Solving Least Squares Problems,
c    SIAM, 1995,
c    ISBN: 0898713560,
c    LC: QA275.L38.
c
c  Parameters:
c
c     All are input arguments.  None are modified by this subroutine.
c
c     A(,)     Array containing matrix to be output.
c     MDA      First dimension of the array, A(,).
c     M, N     No. of rows and columns, respectively in the matrix
c              contained in A(,).
c     NAMES()  [character array]  Array of names.
c              If NAMES(1) contains only blanks, the rest of the NAMES()
c              array will be ignored.
c     MODE     =1  Write header for V matrix and use an F format.
c              =2  Write header for  for candidate solutions and use
c                  P format.
c     UNIT  [integer]   Selects output unit.  If UNIT .ge. 0 then UNIT
c              is the output unit number.  If UNIT = -1, output to
c              the '*' unit.
c     WIDTH [integer]   Selects width of output lines.
c              Each output line from this subroutine will have at most
c              max(26,min(124,WIDTH)) characters plus one additional
c              leading character for Fortran "carriage control".  The
c              carriage control character will always be a blank.
c
      integer I, J, J1, J2, KBLOCK, L, LENNAM
      integer M, MAXCOL, MDA, MODE, N, NAMSIZ, NBLOCK, UNIT, WIDTH
      double precision    A(MDA,N)
      character  NAMES(M)*(*)
      character*4  HEAD (2)
      character*26 FMT1(2)
      character*26 FMT2(2)
      logical      BLKNAM, STAR

      data HEAD(1)/' COL'/
      data HEAD(2)/'SOLN'/
      data FMT1 / '(/7x,00x,8(5x,a4,i4,1x)/)',
     *            '(/7x,00x,8(2x,a4,i4,4x)/)'/
      data FMT2 / '(1x,i4,1x,a00,1x,4p8f14.0)',
     *            '(1x,i4,1x,a00,1x,8g14.6  )'/

      if (M .le. 0 .or. N .le. 0) return
      STAR = UNIT .lt. 0
c
c     The LEN function returns the char length of a single element of
c     the NAMES() array.
c
      LENNAM = len(NAMES(1))
      BLKNAM = NAMES(1) .eq. ' '
      NAMSIZ = 1
      if(.not. BLKNAM) then
         do 30 I = 1,M
            do 10 L = LENNAM, NAMSIZ+1, -1
               if(NAMES(I)(L:L) .ne. ' ') then
                  NAMSIZ = L
                  go to 20
               endif
   10       continue
   20       continue
   30    continue
      endif

      write(FMT1(MODE)(6:7),'(i2.2)') NAMSIZ
      write(FMT2(MODE)(12:13),'(i2.2)') NAMSIZ

   70 format(/' V-Matrix of the Singular Value Decomposition of A*D.'/
     *    ' (Elements of V scaled up by a factor of 10**4)')
   80 format(/' Sequence of candidate solutions, X')

      if(STAR) then
         if (MODE .eq. 1) then
            write (*,70)
         else
            write (*,80)
         endif
      else
         if (MODE .eq. 1) then
            write (UNIT,70)
         else
            write (UNIT,80)
         endif
      endif
c
c  With NAMSIZ characters allowed for the "name" and MAXCOL
c  columns of numbers, the total line width, exclusive of a
c  carriage control character, will be 6 + LENNAM + 14 * MAXCOL.
c
      MAXCOL = max(1,min(8,(WIDTH - 6 - NAMSIZ)/14))

      NBLOCK = (N + MAXCOL -1) / MAXCOL
      J2 = 0

      do 50 KBLOCK = 1, NBLOCK
         J1 = J2 + 1
         J2 = min(N, J2 + MAXCOL)
         if(STAR) then
            write (*,FMT1(MODE)) (HEAD(MODE),J,J=J1,J2)
         else
            write (UNIT,FMT1(MODE)) (HEAD(MODE),J,J=J1,J2)
         endif

         do 40 I=1,M
            if(STAR) then
               if(BLKNAM) then
                  write (*,FMT2(MODE)) I,' ',(A(I,J),J=J1,J2)
               else
                  write (*,FMT2(MODE)) I,NAMES(I),(A(I,J),J=J1,J2)
               endif
            else
               if(BLKNAM) then
                  write (UNIT,FMT2(MODE)) I,' ',(A(I,J),J=J1,J2)
               else
                  write (UNIT,FMT2(MODE)) I,NAMES(I),(A(I,J),J=J1,J2)
               endif
            endif
   40    continue
   50 continue
      end
      SUBROUTINE NNLS (A,MDA,M,N,B,X,RNORM,W,ZZ,INDEX,MODE) 

c*********************************************************************72
c
cc NNLS implements the nonnegative least squares algorithm.
c   
c  Discussion:
c
c    GIVEN AN M BY N MATRIX, A, AND AN M-VECTOR, B,  COMPUTE AN
c    N-VECTOR, X, THAT SOLVES THE LEAST SQUARES PROBLEM   
c   
c      A * X = B  SUBJECT TO X .GE. 0   
c
c  Modified:
c
c    19 October 2008
c
c  Author:
c
c    Charles Lawson, Richard Hanson
c
c  Reference:
c
c    Charles Lawson, Richard Hanson,
c    Solving Least Squares Problems,
c    SIAM, 1995,
c    ISBN: 0898713560,
c    LC: QA275.L38.
c
c  Parameters:
c
c     A(),MDA,M,N     MDA IS THE FIRST DIMENSIONING PARAMETER FOR THE   
c                     ARRAY, A().   ON ENTRY A() CONTAINS THE M BY N    
c                     MATRIX, A.           ON EXIT A() CONTAINS 
c                     THE PRODUCT MATRIX, Q*A , WHERE Q IS AN   
c                     M BY M ORTHOGONAL MATRIX GENERATED IMPLICITLY BY  
c                     THIS SUBROUTINE.  
c
c     B()     ON ENTRY B() CONTAINS THE M-VECTOR, B.   ON EXIT B() CON- 
c             TAINS Q*B.
c
c     X()     ON ENTRY X() NEED NOT BE INITIALIZED.  ON EXIT X() WILL   
c             CONTAIN THE SOLUTION VECTOR.  
c
c     RNORM   ON EXIT RNORM CONTAINS THE EUCLIDEAN NORM OF THE  
c             RESIDUAL VECTOR.  
c
c     W()     AN N-ARRAY OF WORKING SPACE.  ON EXIT W() WILL CONTAIN    
c             THE DUAL SOLUTION VECTOR.   W WILL SATISFY W(I) = 0.  
c             FOR ALL I IN SET P  AND W(I) .LE. 0. FOR ALL I IN SET Z   
c
c     ZZ()     AN M-ARRAY OF WORKING SPACE.     
c
c     INDEX()     AN INTEGER WORKING ARRAY OF LENGTH AT LEAST N.
c                 ON EXIT THE CONTENTS OF THIS ARRAY DEFINE THE SETS    
c                 P AND Z AS FOLLOWS..  
c   
c                 INDEX(1)   THRU INDEX(NSETP) = SET P.     
c                 INDEX(IZ1) THRU INDEX(IZ2)   = SET Z.     
c                 IZ1 = NSETP + 1 = NPP1
c                 IZ2 = N   
c
c     MODE    THIS IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING 
c             MEANINGS. 
c             1     THE SOLUTION HAS BEEN COMPUTED SUCCESSFULLY.
c             2     THE DIMENSIONS OF THE PROBLEM ARE BAD.  
c                   EITHER M .LE. 0 OR N .LE. 0.
c             3    ITERATION COUNT EXCEEDED.  MORE THAN 3*N ITERATIONS. 
c   
      integer I, II, IP, ITER, ITMAX, IZ, IZ1, IZ2, IZMAX, J, JJ, JZ, L
      integer M, MDA, MODE,N, NPP1, NSETP, RTNKEY
      integer INDEX(N)  
      double precision A(MDA,N), B(M), W(N), X(N), ZZ(M)   
      double precision ALPHA, ASAVE, CC, DIFF, DUMMY, FACTOR, RNORM
      double precision SM, SS, T, TEMP, TWO, UNORM, UP, WMAX
      double precision ZERO, ZTEST
      parameter(FACTOR = 0.01d0)
      parameter(TWO = 2.0d0, ZERO = 0.0d0)

      MODE=1
      IF (M .le. 0 .or. N .le. 0) then
         MODE=2
         RETURN
      endif
      ITER=0
      ITMAX=3*N 
c   
c  INITIALIZE THE ARRAYS INDEX() AND X(). 
c   
          DO 20 I=1,N   
          X(I)=ZERO     
   20     INDEX(I)=I    

      IZ2=N 
      IZ1=1 
      NSETP=0   
      NPP1=1
c
c  MAIN LOOP BEGINS HERE.
c  
   30 CONTINUE  
c
c  QUIT IF ALL COEFFICIENTS ARE ALREADY IN THE SOLUTION.
c  OR IF M COLS OF A HAVE BEEN TRIANGULARIZED.    
c   
      IF (IZ1 .GT.IZ2.OR.NSETP.GE.M) GO TO 350   
c   
c  COMPUTE COMPONENTS OF THE DUAL (NEGATIVE GRADIENT) VECTOR W().
c   
      DO 50 IZ=IZ1,IZ2  
         J=INDEX(IZ)   
         SM=ZERO   
         DO 40 L=NPP1,M
   40        SM=SM+A(L,J)*B(L)     
         W(J)=SM   
   50 continue
c
c  FIND LARGEST POSITIVE W(J). 
c
   60 continue
      WMAX=ZERO 
      DO 70 IZ=IZ1,IZ2  
         J=INDEX(IZ)   
         IF (W(J) .gt. WMAX) then
            WMAX=W(J)     
            IZMAX=IZ  
         endif
   70 CONTINUE  
c   
c  IF WMAX .LE. 0. GO TO TERMINATION.
c  THIS INDICATES SATISFACTION OF THE KUHN-TUCKER CONDITIONS.
c   
      IF (WMAX .le. ZERO) go to 350
      IZ=IZMAX  
      J=INDEX(IZ)   
c   
c  THE SIGN OF W(J) IS OK FOR J TO BE MOVED TO SET P.    
c  BEGIN THE TRANSFORMATION AND CHECK NEW DIAGONAL ELEMENT TO AVOID  
c  NEAR LINEAR DEPENDENCE.   
c   
      ASAVE=A(NPP1,J)   
      CALL H12 (1,NPP1,NPP1+1,M,A(1,J),1,UP,DUMMY,1,1,0)    
      UNORM=ZERO
      IF (NSETP .ne. 0) then
          DO 90 L=1,NSETP   
   90       UNORM=UNORM+A(L,J)**2     
      endif
      UNORM=sqrt(UNORM) 
      IF (DIFF(UNORM+ABS(A(NPP1,J))*FACTOR,UNORM) .gt. ZERO) then
c   
c  COL J IS SUFFICIENTLY INDEPENDENT.  COPY B INTO ZZ, UPDATE ZZ
c  AND SOLVE FOR ZTEST ( = PROPOSED NEW VALUE FOR X(J) ).    
c   
         DO 120 L=1,M  
  120        ZZ(L)=B(L)    
         CALL H12 (2,NPP1,NPP1+1,M,A(1,J),1,UP,ZZ,1,1,1)   
         ZTEST=ZZ(NPP1)/A(NPP1,J)  
c   
c  SEE IF ZTEST IS POSITIVE  
c   
         IF (ZTEST .gt. ZERO) go to 140
      endif
c   
c  REJECT J AS A CANDIDATE TO BE MOVED FROM SET Z TO SET P.  
c  RESTORE A(NPP1,J), SET W(J)=0., AND LOOP BACK TO TEST DUAL
c  COEFFS AGAIN.     
c   
      A(NPP1,J)=ASAVE   
      W(J)=ZERO 
      GO TO 60  
c   
c  THE INDEX  J=INDEX(IZ)  HAS BEEN SELECTED TO BE MOVED FROM
c  SET Z TO SET P.    UPDATE B,  UPDATE INDICES,  APPLY HOUSEHOLDER  
c  TRANSFORMATIONS TO COLS IN NEW SET Z,  ZERO SUBDIAGONAL ELTS IN   
c  COL J,  SET W(J)=0.   
c   
  140 continue
      DO 150 L=1,M  
  150    B(L)=ZZ(L)    

      INDEX(IZ)=INDEX(IZ1)  
      INDEX(IZ1)=J  
      IZ1=IZ1+1 
      NSETP=NPP1
      NPP1=NPP1+1   

      IF (IZ1 .le. IZ2) then
         DO 160 JZ=IZ1,IZ2 
            JJ=INDEX(JZ)  
            CALL H12 (2,NSETP,NPP1,M,A(1,J),1,UP,A(1,JJ),1,MDA,1)
  160    continue
      endif

      IF (NSETP .ne. M) then
         DO 180 L=NPP1,M   
  180       A(L,J)=ZERO   
      endif

      W(J)=ZERO 
c
c  SOLVE THE TRIANGULAR SYSTEM.   
c  STORE THE SOLUTION TEMPORARILY IN ZZ().
c
      RTNKEY = 1
      GO TO 400 
  200 CONTINUE  
c   
c  SECONDARY LOOP BEGINS HERE. 
c   
c  ITERATION COUNTER.   
c 
  210 continue  
      ITER=ITER+1   
      IF (ITER .gt. ITMAX) then
         MODE=3
         write (*,'(/a)') ' NNLS quitting on iteration count.'
         GO TO 350 
      endif
c   
c                    SEE IF ALL NEW CONSTRAINED COEFFS ARE FEASIBLE.    
c                                  IF NOT COMPUTE ALPHA.    
c   
      ALPHA=TWO 
      DO 240 IP=1,NSETP 
         L=INDEX(IP)   
         IF (ZZ(IP) .le. ZERO) then
            T=-X(L)/(ZZ(IP)-X(L))     
            IF (ALPHA .gt. T) then
               ALPHA=T   
               JJ=IP 
            endif
         endif
  240 CONTINUE  
c   
c          IF ALL NEW CONSTRAINED COEFFS ARE FEASIBLE THEN ALPHA WILL   
c          STILL = 2.    IF SO EXIT FROM SECONDARY LOOP TO MAIN LOOP.   
c   
      IF (ALPHA.EQ.TWO) GO TO 330   
c   
c          OTHERWISE USE ALPHA WHICH WILL BE BETWEEN 0. AND 1. TO   
c          INTERPOLATE BETWEEN THE OLD X AND THE NEW ZZ.    
c   
      DO 250 IP=1,NSETP 
         L=INDEX(IP)   
         X(L)=X(L)+ALPHA*(ZZ(IP)-X(L)) 
  250 continue
c   
c        MODIFY A AND B AND THE INDEX ARRAYS TO MOVE COEFFICIENT I  
c        FROM SET P TO SET Z.   
c   
      I=INDEX(JJ)   
  260 continue
      X(I)=ZERO 
c   
      IF (JJ .ne. NSETP) then
         JJ=JJ+1   
         DO 280 J=JJ,NSETP 
            II=INDEX(J)   
            INDEX(J-1)=II 
            CALL G1 (A(J-1,II),A(J,II),CC,SS,A(J-1,II))   
            A(J,II)=ZERO  
            DO 270 L=1,N  
               IF (L.NE.II) then
c
c  Apply procedure G2 (CC,SS,A(J-1,L),A(J,L))  
c
                  TEMP = A(J-1,L)
                  A(J-1,L) = CC*TEMP + SS*A(J,L)
                  A(J,L)   =-SS*TEMP + CC*A(J,L)
               endif
  270       CONTINUE  
c
c  Apply procedure G2 (CC,SS,B(J-1),B(J))   
c
            TEMP = B(J-1)
            B(J-1) = CC*TEMP + SS*B(J)    
            B(J)   =-SS*TEMP + CC*B(J)    
  280    continue
      endif

      NPP1=NSETP
      NSETP=NSETP-1     
      IZ1=IZ1-1 
      INDEX(IZ1)=I  
c   
c  SEE IF THE REMAINING COEFFS IN SET P ARE FEASIBLE.  THEY SHOULD
c  BE BECAUSE OF THE WAY ALPHA WAS DETERMINED.
c  IF ANY ARE INFEASIBLE IT IS DUE TO ROUND-OFF ERROR.  ANY   
c  THAT ARE NONPOSITIVE WILL BE SET TO ZERO   
c  AND MOVED FROM SET P TO SET Z. 
c   
      DO 300 JJ=1,NSETP 
         I=INDEX(JJ)   
         IF (X(I) .le. ZERO) go to 260
  300 CONTINUE  
c   
c         COPY B( ) INTO ZZ( ).  THEN SOLVE AGAIN AND LOOP BACK.
c   
      DO 310 I=1,M  
  310     ZZ(I)=B(I)    
      RTNKEY = 2
      GO TO 400 
  320 CONTINUE  
      GO TO 210 
c
c  END OF SECONDARY LOOP,
c   
  330 continue
      DO 340 IP=1,NSETP 
          I=INDEX(IP)   
  340     X(I)=ZZ(IP)   
c
c        ALL NEW COEFFS ARE POSITIVE.  LOOP BACK TO BEGINNING.  
c
      GO TO 30  
c   
c  END OF MAIN LOOP
c   
c  COME TO HERE FOR TERMINATION.  
c  COMPUTE THE NORM OF THE FINAL RESIDUAL VECTOR.    
c 
  350 continue  
      SM=ZERO   
      IF (NPP1 .le. M) then
         DO 360 I=NPP1,M   
  360       SM=SM+B(I)**2 
      else
         DO 380 J=1,N  
  380       W(J)=ZERO     
      endif
      RNORM=sqrt(SM)    
      RETURN
c   
c     THE FOLLOWING BLOCK OF CODE IS USED AS AN INTERNAL SUBROUTINE     
c     TO SOLVE THE TRIANGULAR SYSTEM, PUTTING THE SOLUTION IN ZZ().     
c   
  400 continue
      DO 430 L=1,NSETP  
         IP=NSETP+1-L  
         IF (L .ne. 1) then
            DO 410 II=1,IP
               ZZ(II)=ZZ(II)-A(II,JJ)*ZZ(IP+1)   
  410       continue
         endif
         JJ=INDEX(IP)  
         ZZ(IP)=ZZ(IP)/A(IP,JJ)    
  430 continue
      go to (200, 320), RTNKEY
      END     
      SUBROUTINE QRBD (IPASS,Q,E,NN,V,MDV,NRV,C,MDC,NCC)    

c*********************************************************************72
c
cc QRBD uses the QR algorithm for the singular values of a bidiagonal matrix.
c
c  Discussion:
c
c     THE BIDIAGONAL MATRIX 
c   
c                       (Q1,E2,0...    )
c                       (   Q2,E3,0... )
c                D=     (       .      )
c                       (         .   0)
c                       (           .EN)
c                       (          0,QN)
c   
c                 IS PRE AND POST MULTIPLIED BY 
c                 ELEMENTARY ROTATION MATRICES  
c                 RI AND PI SO THAT 
c   
c                 RK...R1*D*P1**(T)...PK**(T) = DIAG(S1,...,SN) 
c   
c                 TO WITHIN WORKING ACCURACY.   
c   
c    1. EI AND QI OCCUPY E(I) AND Q(I) AS INPUT.  
c   
c    2. RM...R1*C REPLACES 'C' IN STORAGE AS OUTPUT.  
c   
c    3. V*P1**(T)...PM**(T) REPLACES 'V' IN STORAGE AS OUTPUT.
c   
c    4. SI OCCUPIES Q(I) AS OUTPUT.   
c   
c    5. THE SI'S ARE NONINCREASING AND NONNEGATIVE.   
c   
c  Modified:
c
c    19 October 2008
c
c  Author:
c
c    Charles Lawson, Richard Hanson
c
c  Reference:
c
c    Gene Golub, Christian Reinsch,
c    Singular Value Decomposition and Least Squares Solutions,
c    Numerische Mathematik,
c    Volume 14, Number 5, April 1970, pages 403-420.
c
c    Charles Lawson, Richard Hanson,
c    Solving Least Squares Problems,
c    SIAM, 1995,
c    ISBN: 0898713560,
c    LC: QA275.L38.
c
      integer MDC, MDV, NCC, NN, NRV
      double precision C(MDC,NCC), E(NN), Q(NN),V(MDV,NN)
      integer I, II, IPASS, J, K, KK, L, LL, LP1, N, N10, NQRS
      double precision CS, DIFF, DNORM, F, G, H, SMALL
      double precision ONE, SN, T, TEMP, TWO, X, Y, Z, ZERO
      
      logical WNTV ,HAVERS,FAIL     
      parameter(ONE = 1.0d0, TWO = 2.0d0, ZERO = 0.0d0)
  
      N=NN  
      IPASS=1   
      IF (N.LE.0) RETURN
      N10=10*N  
      WNTV=NRV.GT.0     
      HAVERS=NCC.GT.0   
      FAIL=.FALSE.  
      NQRS=0
      E(1)=ZERO
 
      DNORM=ZERO
      DO J=1,N  
        DNORM=max(abs(Q(J))+abs(E(J)),DNORM)   
      end do

           DO 200 KK=1,N
           K=N+1-KK     
c   
c  TEST FOR SPLITTING OR RANK DEFICIENCIES.. 
c   FIRST MAKE TEST FOR LAST DIAGONAL TERM, Q(K), BEING SMALL.   
c 
   20       IF(K.EQ.1) GO TO 50     
            IF(DIFF(DNORM+Q(K),DNORM) .ne. ZERO) go to 50
c   
c  SINCE Q(K) IS SMALL WE WILL MAKE A SPECIAL PASS TO    
c  TRANSFORM E(K) TO ZERO.   
c   
           CS=ZERO  
           SN=-ONE  
                DO 40 II=2,K
                I=K+1-II
                F=-SN*E(I+1)
                E(I+1)=CS*E(I+1)    
                CALL G1 (Q(I),F,CS,SN,Q(I))     
c
c   TRANSFORMATION CONSTRUCTED TO ZERO POSITION (I,K).
c   
                IF (.NOT.WNTV) GO TO 40 
                     DO 30 J=1,NRV  
c
c  Apply procedure G2 (CS,SN,V(J,I),V(J,K))  
c
                        TEMP = V(J,I)
                        V(J,I) = CS*TEMP + SN*V(J,K)
                        V(J,K) =-SN*TEMP + CS*V(J,K)
   30                continue
c
c  ACCUMULATE RT. TRANSFORMATIONS IN V. 
c   
   40           CONTINUE
c   
c  THE MATRIX IS NOW BIDIAGONAL, AND OF LOWER ORDER SINCE E(K) .EQ. ZERO..    
c   
   50           DO 60 LL=1,K
                  L=K+1-LL
                  IF(DIFF(DNORM+E(L),DNORM) .eq. ZERO) go to 100
                  IF(DIFF(DNORM+Q(L-1),DNORM) .eq. ZERO) go to 70
   60           CONTINUE
c
c     THIS LOOP CAN'T COMPLETE SINCE E(1) = ZERO.   
c   
           GO TO 100    
c   
c  CANCELLATION OF E(L), L.GT.1. 
c
   70      CS=ZERO  
           SN=-ONE  
                DO 90 I=L,K 
                F=-SN*E(I)  
                E(I)=CS*E(I)
                IF(DIFF(DNORM+F,DNORM) .eq. ZERO) go to 100
                CALL G1 (Q(I),F,CS,SN,Q(I))     
                IF (HAVERS) then
                     DO 80 J=1,NCC  
c
c  Apply procedure G2 ( CS, SN, C(I,J), C(L-1,J)
c
                        TEMP = C(I,J)
                        C(I,J)   = CS*TEMP + SN*C(L-1,J)
                        C(L-1,J) =-SN*TEMP + CS*C(L-1,J)
   80                continue
                endif
   90           CONTINUE
c   
c         TEST FOR CONVERGENCE.
c
  100      Z=Q(K)   
           IF (L.EQ.K) GO TO 170    
c   
c         SHIFT FROM BOTTOM 2 BY 2 MINOR OF B**(T)*B.
c
           X=Q(L)   
           Y=Q(K-1)     
           G=E(K-1)     
           H=E(K)   
           F=((Y-Z)*(Y+Z)+(G-H)*(G+H))/(TWO*H*Y)
           G=sqrt(ONE+F**2) 
           IF (F .ge. ZERO) then
              T=F+G
           else
              T=F-G
           endif
           F=((X-Z)*(X+Z)+H*(Y/T-H))/X  
c   
c  NEXT QR SWEEP.
c
           CS=ONE   
           SN=ONE   
           LP1=L+1  
                DO 160 I=LP1,K  
                G=E(I)  
                Y=Q(I)  
                H=SN*G  
                G=CS*G  
                CALL G1 (F,H,CS,SN,E(I-1))  
                F=X*CS+G*SN 
                G=-X*SN+G*CS
                H=Y*SN  
                Y=Y*CS  
                IF (WNTV) then
c   
c  ACCUMULATE ROTATIONS (FROM THE RIGHT) IN 'V' 
c
                     DO 130 J=1,NRV 
c                    
c  Apply procedure G2 (CS,SN,V(J,I-1),V(J,I))
c
                        TEMP = V(J,I-1)
                        V(J,I-1) = CS*TEMP + SN*V(J,I)
                        V(J,I)   =-SN*TEMP + CS*V(J,I)
  130                continue
                endif
                CALL G1 (F,H,CS,SN,Q(I-1))  
                F=CS*G+SN*Y 
                X=-SN*G+CS*Y
                IF (HAVERS) then
                     DO 150 J=1,NCC 
c
c  Apply procedure G2 (CS,SN,C(I-1,J),C(I,J))
c
                        TEMP = C(I-1,J)
                        C(I-1,J) = CS*TEMP + SN*C(I,J)
                        C(I,J)   =-SN*TEMP + CS*C(I,J)
  150                continue
                endif
c
c              APPLY ROTATIONS FROM THE LEFT TO 
c              RIGHT HAND SIDES IN 'C'..
c   
  160           CONTINUE
           E(L)=ZERO    
           E(K)=F   
           Q(K)=X   
           NQRS=NQRS+1  
           IF (NQRS.LE.N10) GO TO 20
c
c  RETURN TO 'TEST FOR SPLITTING'.  
c 
           SMALL=ABS(E(K))
           I=K     
c
c  IF FAILURE TO CONVERGE SET SMALLEST MAGNITUDE
c  TERM IN OFF-DIAGONAL TO ZERO.  CONTINUE ON.
c
                DO 165 J=L,K 
                TEMP=ABS(E(J))
                IF(TEMP .EQ. ZERO) GO TO 165
                IF(TEMP .LT. SMALL) THEN
                     SMALL=TEMP
                     I=J
                end if  
  165           CONTINUE
           E(I)=ZERO
           NQRS=0
           FAIL=.TRUE.  
           GO TO 20
c
c  CUTOFF FOR CONVERGENCE FAILURE. 'NQRS' WILL BE 2*N USUALLY.
c 
  170      IF (Z.GE.ZERO) GO TO 190 
           Q(K)=-Z  
           IF (WNTV) then
                DO 180 J=1,NRV  
  180           V(J,K)=-V(J,K)  
           endif
  190      CONTINUE     
c
c  CONVERGENCE. Q(K) IS MADE NONNEGATIVE..   
c   
  200      CONTINUE     
      IF (N.EQ.1) RETURN
           DO 210 I=2,N 
           IF (Q(I).GT.Q(I-1)) GO TO 220
  210      CONTINUE     
      IF (FAIL) IPASS=2 
      RETURN
c
c     EVERY SINGULAR VALUE IS IN ORDER..
c
  220      DO 270 I=2,N 
           T=Q(I-1)     
           K=I-1
                DO 230 J=I,N
                IF (T.GE.Q(J)) GO TO 230
                T=Q(J)  
                K=J     
  230           CONTINUE
           IF (K.EQ.I-1) GO TO 270  
           Q(K)=Q(I-1)  
           Q(I-1)=T     
           IF (HAVERS) then
                DO 240 J=1,NCC  
                T=C(I-1,J)  
                C(I-1,J)=C(K,J)     
  240           C(K,J)=T
           endif

  250      IF (WNTV) then
                DO 260 J=1,NRV  
                T=V(J,I-1)  
                V(J,I-1)=V(J,K)     
  260           V(J,K)=T
           endif
  270      CONTINUE     
c
c         END OF ORDERING ALGORITHM.
c   
      IF (FAIL) IPASS=2 
      RETURN
      END   
      subroutine SVA(A,MDA,M,N,MDATA,B,SING,KPVEC,NAMES,ISCALE,D,WORK)

c*********************************************************************72
c
cc SVA carries out a singular value analysis.
c
c  Discussion:
c
c    SINGULAR VALUE ANALYSIS.  COMPUTES THE SINGULAR VALUE
c    DECOMPOSITION OF THE MATRIX OF A LEAST SQUARES PROBLEM, AND
c    PRODUCES A PRINTED REPORT.
c
c    This 1995 version differs from the original 1973 version by the 
c    addition of the arguments KPVEC() and WORK(), and by allowing user to
c    choose the length of names in NAMES().
c    KPVEC() allows the user to exercise options regarding printing.
c    WORK() provides 2*N locations of work space.  Originally SING() was
c    required to have 3*N elements, of which the last 2*N were used for
c    work space.  Now SING() only needs N elements.
c
c  Modified:
c
c    19 October 2008
c
c  Author:
c
c    Charles Lawson, Richard Hanson
c
c  Reference:
c
c    Charles Lawson, Richard Hanson,
c    Solving Least Squares Problems,
c    SIAM, 1995,
c    ISBN: 0898713560,
c    LC: QA275.L38.
c
c  Parameters:
c
c     A(,)     [inout]  On entry, contains the M x N matrix of the least
c              squares problem to be analyzed.  This could be a matrix
c              obtained by preliminary orthogonal transformations
c              applied to the actual problem matrix which may have had
c              more rows (See MDATA below.)
c
c     MDA   [in]  First dimensioning parameter for A(,).  Require
c              MDA .ge. max(M, N).
c
c     M,N   [in]  No. of rows and columns, respectively, in the
c              matrix, A.  Either M > N or M .le. N is permitted.
c              Require M > 0 and N > 0.
c
c     MDATA [in]  No. of rows in actual least squares problem.
c              Generally MDATA .ge. M.  MDATA is used only in computing
c              statistics for the report and is not used as a loop
c              count or array dimension.
c
c     B()   [inout]  On entry, contains the right-side vector, b, of the
c              least squares problem.  This vector is of length, M.
c              On return, contains the vector, g = (U**t)*b, where U
c              comes from the singular value decomposition of A.  The
c              vector , g, is also of length M.
c
c     SING()   [out] On return, contains the singular values of A, in
c              descending order, in locations indexed 1 thru min(M,N).
c              If M < N, locations indexed from M+1 through N will be
c              set to zero.
c
c     KPVEC()  [integer array, in]  Array of integers to select print
c              options.  KPVEC(1) determines whether the rest of
c              the array is to be used or ignored.
c              If KPVEC(1) = 1, the contents of (KPVEC(I), I=2,4)
c              will be used to set internal variables as follows:
c                       PRBLK = KPVEC(2)
c                       UNIT  = KPVEC(3)
c                       WIDTH = KPVEC(4)
c              If KPVEC(1) = 0 default settings will be used.  The user
c              need not dimension KPVEC() greater than 1.  The subr will
c              set PRBLK = 111111, UNIT = -1, and WIDTH = 69.
c
c              The internal variables PRBLK, UNIT, and WIDTH are
c              interpreted as follows:
c
c              PRBLK    The decimal representation of PRBLK must be
c              representable as at most 6 digits, each being 0 or 1.
c              The decimal digits will be interpreted as independant
c              on/off flags for the 6 possible blocks of printed output.
c              Examples:  111111 selects all blocks, 0 suppresses all
c              printing,  101010 selects the 1st, 3rd, and 5th blocks,
c              etc.
c              The six blocks are:
c              1. Header, with size and scaling option parameters.
c              2. V-matrix.  Amount of output depends on M and N.
c              3. Singular values and related quantities.  Amount of
c                 output depends on N.
c              4. Listing of YNORM and RNORM and their logarithms.
c                 Amount of output depends on N.
c              5. Levenberg-Marquart analysis.
c              6. Candidate solutions.  Amount of output depends on
c                 M and N.
c
c              UNIT     Selects the output unit.  If UNIT .ge. 0,
c              UNIT will be used as the output unit number.
c              If UNIT = -1, output will be written to the "*" output
c              unit, i.e., the standard system output unit.
c              The calling program unit is responsible for opening
c              and/or closing the selected output unit if the host
c              system requires these actions.
c
c        WIDTH    Default value is 79.  Determines the width of
c        blocks 2, 3, and 6 of the output report.
c        Block 3 will use 95(+1) cols if WIDTH .ge. 95, and otherwise
c        69(+1) cols.
c        Blocks 2 and 6 are printed by subroutine MFEOUT.  These blocks
c        generally use at most WIDTH(+1) cols, but will use more if
c        the names are so long that more space is needed to print one
c        name and one numeric column.  The (+1)'s above are reminders
c        that in all cases there is one extra initial column for Fortran
c        "carriage control".  The carriage control character will always
c        be a blank.
c        Blocks 1, 4, and 5 have fixed widths of 63(+1), 66(+1) and
c        66(+1), respectively.
c
c     NAMES()  [in]  NAMES(j), for j = 1, ..., N, may contain a
c              name for the jth component of the solution
c              vector.  The declared length of the elements of the
c              NAMES() array is not specifically limited, but a
c              greater length reduces the space available for columns
c              of the matris to be printed.
c              If NAMES(1) contains only blank characters,
c              it will be assumed that no names have been provided,
c              and this subr will not access the NAMES() array beyond
c              the first element.
c
c     ISCALE   [in]  Set by the user to 1, 2, or 3 to select the column
c              scaling option.
c                  1   SUBR WILL USE IDENTITY SCALING AND IGNORE THE D()
c                      ARRAY.
c                  2   SUBR WILL SCALE NONZERO COLS TO HAVE UNIT EUCLID-
c                      EAN LENGTH AND WILL STORE RECIPROCAL LENGTHS OF
c                      ORIGINAL NONZERO COLS IN D().
c                  3   USER SUPPLIES COL SCALE FACTORS IN D(). SUBR
c                      WILL MULT COL J BY D(J) AND REMOVE THE SCALING
c                      FROM THE SOLN AT THE END.
c
c     D()   [ignored or out or in]  Usage of D() depends on ISCALE as
c              described above.  When used, its length must be
c              at least N.
c
c     WORK()   [scratch]   Work space of length at least 2*N.  Used
c              directly in this subr and also in _SVDRS.
c
      integer I, IE, IPASS, ISCALE, J, K, KPVEC(4), M, MDA, MDATA
      integer MINMN, MINMN1, MPASS, N, NSOL
      integer PRBLK, UNIT, WIDTH
      double precision A(MDA,N), A1, A2, A3, A4, ALAMB, ALN10, B(M)
      double precision D(N), DEL, EL, EL2
      double precision ONE, PCOEF, RL, RNORM, RS
      double precision SB, SING(N), SL, TEN, THOU, TWENTY
      double precision WORK(2*N), YL, YNORM, YS, YSQ, ZERO
      character*(*) NAMES(N)
      logical BLK(6), NARROW, STAR
      parameter( ZERO = 0.0d0, ONE = 1.0d0)
      parameter( TEN = 10.0d0, TWENTY = 20.0d0, THOU = 1000.0d0)

  220 format (1X/' INDEX  SING. VAL.    P COEF   ',
     *        '   RECIPROCAL    G COEF         G**2   ',
     *        '   CUMULATIVE  SCALED SQRT'/
     *    31x,'   SING. VAL.',26x,
     *        '  SUM of SQRS  of CUM.S.S.')
  221 format (1X/' INDEX  SING. VAL.    P COEF   ',
     *        '   RECIPROCAL    G COEF     SCALED SQRT'/
     *    31x,'   SING. VAL.',13x,'  of CUM.S.S.')
  222 format (1X/' INDEX  SING. VAL.    G COEF         G**2   ',
     *        '   CUMULATIVE  SCALED SQRT'/
     *    44x,'  SUM of SQRS  of CUM.S.S.')

  230 format (' ',4X,'0',64X,2g13.4)
  231 format (' ',4X,'0',51X, g13.4)
  232 format (' ',4X,'0',38X,2g13.4)

  240 format (' ',i5,g12.4,6g13.4)

  260 format (1X,' M  = ',I6,',   N  =',I4,',   MDATA  =',I8)
  270 format (1X/' Singular Value Analysis of the least squares',
     *        ' problem,  A*X = B,'/
     *        ' scaled as (A*D)*Y = B.')
  280 format (1X/' Scaling option No.',I2,'.  D is a diagonal',
     * ' matrix with the following diagonal elements..'/(5X,10E12.4))
  290 format (1X/' Scaling option No. 1.   D is the identity matrix.'/
     * 1X)
  300 format (1X/' INDEX',12X,'YNORM      RNORM',11X,
     * '      LOG10      LOG10'/
     * 45X,'      YNORM      RNORM'/1X)
  310 format (' ',I5,6X,2E11.3,11X,2F11.3)
  320 format (1X/
     *' Norms of solution and residual vectors for a range of values'/
     *' of the Levenberg-Marquardt parameter, LAMBDA.'//
     *        '      LAMBDA      YNORM      RNORM',
     *         '      LOG10      LOG10      LOG10'/
     *     34X,'     LAMBDA      YNORM      RNORM')
  330 format (1X, 3E11.3, 3F11.3)

      IF (M.LE.0 .OR. N.LE.0) then
        RETURN
      end if

      MINMN = min(M,N)
      MINMN1 = MINMN + 1

      if(KPVEC(1) .eq. 0) then
         PRBLK = 111111
         UNIT = -1
         WIDTH = 79
      else
         PRBLK = KPVEC(2)
         UNIT = KPVEC(3)
         WIDTH = KPVEC(4)
      endif
      STAR = UNIT .lt. 0
c
c  Build logical array BLK() by testing decimal digits of PRBLK.
c
      do I=6, 1, -1
         J = PRBLK/10
         BLK(I) = (PRBLK - 10*J) .gt. 0
         PRBLK = J
      end do
c
c  Optionally print header and M, N, MDATA.
c
      if(BLK(1)) then
         if(STAR) then
            write (*,270)
            write (*,260) M,N,MDATA
         else
            write (UNIT,270)
            write (UNIT,260) M,N,MDATA
         endif
      endif
c
c  Handle scaling as selected by ISCALE.
c
      if( ISCALE .eq. 1) then
         if(BLK(1)) then
            if(STAR) then
               write (*,290)
            else
               write (UNIT,290)
            endif
         endif
      else
c
c  Apply column scaling to A.
c
         DO 52 J = 1,N
            A1 = D(J)
            if( ISCALE .le. 2) then
               SB = ZERO
               DO 30 I = 1,M
   30             SB = SB + A(I,J)**2
               A1 = sqrt(SB)
               IF (A1.EQ.ZERO) A1 = ONE
               A1 = ONE/A1
               D(J) = A1
            endif
            DO 50 I = 1,M
               A(I,J) = A(I,J)*A1
   50       continue
   52    continue
         if(BLK(1)) then
            if(STAR) then
               write (*,280) ISCALE,(D(J),J = 1,N)
            else
               write (UNIT,280) ISCALE,(D(J),J = 1,N)
            endif
         endif
      endif
c
c  Compute the Singular Value Decomposition of the scaled matrix.
c
      call SVDRS (A,MDA,M,N,B,M,1,SING,WORK)
c
c  Determine NSOL.
c
      NSOL = MINMN
      do 60 J = 1,MINMN
         if(SING(J) .eq. ZERO) then
            NSOL = J-1
            go to 65
         endif
   60 continue
   65 continue
c
c  The array B() contains the vector G.
c  Compute cumulative sums of squares of components of
c   G and store them in WORK(I), I = 1,...,MINMN+1
c
      SB = ZERO
      DO 70 I = MINMN1,M
         SB = SB + B(I)**2
   70 CONTINUE
      WORK(MINMN+1) = SB
      DO  75 J = MINMN, 1, -1
         SB = SB + B(J)**2
         WORK(J) = SB
   75 CONTINUE
c
c  PRINT THE V MATRIX.
c
      if(BLK(2)) CALL MFEOUT (A,MDA,N,N,NAMES,1, UNIT, WIDTH)
c
c  REPLACE V BY D*V IN THE ARRAY A()
c
      if (ISCALE .gt.1) then
         do 82 I = 1,N
            do 80 J = 1,N
               A(I,J) = D(I)*A(I,J)
   80       continue
   82    continue
      endif

      if(BLK(3)) then
c
c  Print singular values and other summary results.
c
c  Output will be done using one of two layouts.  The narrow
c  layout uses 69 cols + 1 for carriage control, and makes two passes
c  through the computation.
c  The wide layout uses 95 cols + 1 for carriage control, and makes
c  only one pass through the computation.
c
c  G  NOW IN  B() ARRAY.  V NOW IN A(,) ARRAY.
c
      NARROW = WIDTH .lt. 95
      MPASS = 1
      if(NARROW) MPASS = 2
      do 170 IPASS = 1, MPASS
         if(STAR) then
            if(NARROW) then
               if(IPASS .eq. 1) then
                  write(*,221)
               else
                  write(*,222)
               endif
            else
               write (*,220)
            endif
         else
            if(NARROW) then
               if(IPASS .eq. 1) then
                  write(UNIT,221)
               else
                  write(UNIT,222)
               endif
            else
               write (UNIT,220)
            endif
         endif
c
c The following stmt converts from integer to floating-point.
c
      A3 = WORK(1)
      A4 = sqrt(A3/ max(1,MDATA))
      if(STAR) then
         if(NARROW) then
            if(IPASS .eq. 1) then
               write(*,231) A4
            else
               write(*,232) A3, A4
            endif
         else
            write (*,230) A3,A4
         endif
      else
         if(NARROW) then
            if(IPASS .eq. 1) then
               write(UNIT,231) A4
            else
               write(UNIT,232) A3, A4
            endif
         else
            write (UNIT,230) A3,A4
         endif
      endif

      DO 160 K = 1,MINMN
         if (SING(K).EQ.ZERO) then
            PCOEF = ZERO
            if(STAR) then
               write (*,240) K,SING(K)
            else
               write (UNIT,240) K,SING(K)
            endif
         else
            PCOEF = B(K) / SING(K)
            A1 = ONE / SING(K)
            A2 = B(K)**2
            A3 = WORK(K+1)
            A4 = sqrt(A3/max(1,MDATA-K))
            if(STAR) then
               if(NARROW) then
                  if(IPASS .eq. 1) then
                     write(*,240) K,SING(K),PCOEF,A1,B(K),      A4
                  else
                     write(*,240) K,SING(K),         B(K),A2,A3,A4
                  endif
               else
                  write (*,240) K,SING(K),PCOEF,A1,B(K),A2,A3,A4
               endif
            else
               if(NARROW) then
                  if(IPASS .eq. 1) then
                     write(UNIT,240) K,SING(K),PCOEF,A1,B(K),      A4
                  else
                     write(UNIT,240) K,SING(K),         B(K),A2,A3,A4
                  endif
               else
                  write (UNIT,240) K,SING(K),PCOEF,A1,B(K),A2,A3,A4
               endif

            endif
         endif
  160 continue
  170 continue
      endif

      if( BLK(4) ) then
c
c  Compute and print values of YNORM, RNORM and their logarithms.
c
      if(STAR) then
         write (*,300)
      else
         write (UNIT,300)
      endif
      YSQ = ZERO
      do 180 J = 0, NSOL
         if(J .ne. 0) YSQ = YSQ + (B(J) / SING(J))**2
         YNORM = sqrt(YSQ)
         RNORM = sqrt(WORK(J+1))
         YL = -THOU
         IF (YNORM .GT. ZERO) YL = log10(YNORM)
         RL = -THOU
         IF (RNORM .GT. ZERO) RL = log10(RNORM)
         if(STAR) then
            write (*,310) J,YNORM,RNORM,YL,RL
         else
            write (UNIT,310) J,YNORM,RNORM,YL,RL
         endif
  180 continue
      endif

      if( BLK(5) .and. SING(1) .ne. ZERO ) then
c
c  COMPUTE VALUES OF XNORM AND RNORM FOR A SEQUENCE OF VALUES OF
c  THE LEVENBERG-MARQUARDT PARAMETER.
c
         EL = log10(SING(1)) + ONE
         EL2 = log10(SING(NSOL)) - ONE
         DEL = (EL2-EL) / TWENTY
         ALN10 = log(TEN)
         if(STAR) then
            write (*,320)
         else
            write (UNIT,320)
         endif
         DO 200 IE = 1,21
c
c  COMPUTE        ALAMB = 10.0**EL
c
            ALAMB = EXP(ALN10*EL)
            YS = ZERO
            RS = WORK(NSOL+1)
            DO 190 I = 1,MINMN
               SL = SING(I)**2 + ALAMB**2
               YS = YS + (B(I)*SING(I)/SL)**2
               RS = RS + (B(I)*(ALAMB**2)/SL)**2
  190       CONTINUE
            YNORM = sqrt(YS)
            RNORM = sqrt(RS)
            RL = -THOU
            IF (RNORM.GT.ZERO) RL = log10(RNORM)
            YL = -THOU
            IF (YNORM.GT.ZERO) YL = log10(YNORM)
            if(STAR) then
               write (*,330) ALAMB,YNORM,RNORM,EL,YL,RL
            else
               write (UNIT,330) ALAMB,YNORM,RNORM,EL,YL,RL
            endif
            EL = EL + DEL
  200    CONTINUE
      endif
c
c  Compute and optionally print candidate solutions.
c
      do 215 K = 1,NSOL
         PCOEF = B(K) / SING(K)
         DO 210 I = 1,N
                A(I,K) = A(I,K) * PCOEF
  210           IF (K.GT.1) A(I,K) = A(I,K) + A(I,K-1)
  215 continue
      if (BLK(6) .and. NSOL.GE.1)
     *       CALL MFEOUT (A,MDA,N,NSOL,NAMES,2,UNIT,WIDTH)
      return
      END
      subroutine SVDRS (A, MDA, M1, N1, B, MDB, NB, S, WORK) 

c*********************************************************************72
c
cc SVDRS does a singular value decomposition with a  right side vector.
c
c  Discussion:
c  
c    This 1995 version differs from the original 1974 version by adding
c    the argument WORK().
c    WORK() provides 2*N1 locations of work space.  Originally S() was
c    required to have 3*N1 elements, of which the last 2*N1 were used for
c    work space.  Now S() only needs N1 elements.
c
c    This subroutine computes the singular value decomposition of the
c    given M1 x N1 matrix, A, and optionally applys the transformations
c    from the left to the NB column vectors of the M1 x NB matrix B.
c    Either M1 .ge. N1  or  M1 .lt. N1 is permitted.
c
c    The singular value decomposition of A is of the form
c
c      A  =  U * S * V**t
c
c    where U is M1 x M1 orthogonal, S is M1 x N1 diagonal with the
c    diagonal terms nonnegative and ordered from large to small, and
c    V is N1 x N1 orthogonal.  Note that these matrices also satisfy
c
c      S = (U**t) * A * V
c
c    The matrix V is returned in the leading N1 rows and
c    columns of the array A(,).
c
c    The singular values, i.e. the diagonal terms of the matrix S,
c    are returned in the array S().  If M1 .lt. N1, positions M1+1
c    through N1 of S() will be set to zero.
c
c    The product matrix  G = U**t * B replaces the given matrix B
c    in the array B(,).
c
c    If the user wishes to obtain a minimum length least squares
c    solution of the linear system
c
c      A * X ~=~ B
c
c    the solution X can be constructed, following use of this subroutine,
c    by computing the sum for i = 1, ..., R of the outer products
c
c      (Col i of V) * (1/S(i)) * (Row i of G)
c
c    Here R denotes the pseudorank of A which the user may choose
c    in the range 0 through Min(M1, N1) based on the sizes of the
c    singular values.
c
c
c    This code gives special treatment to rows and columns that are
c    entirely zero.  This causes certain zero sing. vals. to appear as
c    exact zeros rather than as about MACHEPS times the largest sing. val.
c    It similarly cleans up the associated columns of U and V.  
c
c    METHOD..  
c    1. EXCHANGE COLS OF A TO PACK NONZERO COLS TO THE LEFT.  
c       SET N = NO. OF NONZERO COLS.  
c       USE LOCATIONS A(1,N1),A(1,N1-1),...,A(1,N+1) TO RECORD THE    
c       COL PERMUTATIONS. 
c    2. EXCHANGE ROWS OF A TO PACK NONZERO ROWS TO THE TOP.   
c       QUIT PACKING IF FIND N NONZERO ROWS.  MAKE SAME ROW EXCHANGES 
c       IN B.  SET M SO THAT ALL NONZERO ROWS OF THE PERMUTED A   
c       ARE IN FIRST M ROWS.  IF M .LE. N THEN ALL M ROWS ARE 
c       NONZERO.  IF M .GT. N THEN      THE FIRST N ROWS ARE KNOWN    
c       TO BE NONZERO,AND ROWS N+1 THRU M MAY BE ZERO OR NONZERO.     
c    3. APPLY ORIGINAL ALGORITHM TO THE M BY N PROBLEM.   
c    4. MOVE PERMUTATION RECORD FROM A(,) TO S(I),I=N+1,...,N1.   
c    5. BUILD V UP FROM  N BY N  TO  N1 BY N1  BY PLACING ONES ON     
c       THE DIAGONAL AND ZEROS ELSEWHERE.  THIS IS ONLY PARTLY DONE   
c       EXPLICITLY.  IT IS COMPLETED DURING STEP 6.   
c    6. EXCHANGE ROWS OF V TO COMPENSATE FOR COL EXCHANGES OF STEP 2. 
c    7. PLACE ZEROS IN  S(I),I=N+1,N1  TO REPRESENT ZERO SING VALS.   
c 
c  Modified:
c
c    19 October 2008
c
c  Author:
c
c    Charles Lawson, Richard Hanson
c
c  Reference:
c
c    Charles Lawson, Richard Hanson,
c    Solving Least Squares Problems,
c    SIAM, 1995,
c    ISBN: 0898713560,
c    LC: QA275.L38.
c
c  Parameters:
c
c  A(,)     (In/Out)  On input contains the M1 x N1 matrix A.
c           On output contains the N1 x N1 matrix V.
c
c  LDA      (In)  First dimensioning parameter for A(,).
c           Require LDA .ge. Max(M1, N1).
c
c  M1       (In)  No. of rows of matrices A, B, and G.
c           Require M1 > 0.
c
c  N1       (In)  No. of cols of matrix A, No. of rows and cols of
c           matrix V.  Permit M1 .ge. N1  or  M1 .lt. N1.
c           Require N1 > 0.
c
c  B(,)     (In/Out)  If NB .gt. 0 this array must contain an
c           M1 x NB matrix on input and will contain the
c           M1 x NB product matrix, G = (U**t) * B on output.
c
c  LDB      (In)  First dimensioning parameter for B(,).
c           Require LDB .ge. M1.
c
c  NB       (In)  No. of cols in the matrices B and G.
c           Require NB .ge. 0.
c
c  S()      (Out)  Must be dimensioned at least N1.  On return will
c           contain the singular values of A, with the ordering
c                S(1) .ge. S(2) .ge. ... .ge. S(N1) .ge. 0.
c           If M1 .lt. N1 the singular values indexed from M1+1
c           through N1 will be zero.
c           If the given integer arguments are not consistent, this
c           subroutine will return immediately, setting S(1) = -1.0.
c
c  WORK()  (Scratch)  Work space of total size at least 2*N1.
c           Locations 1 thru N1 will hold the off-diagonal terms of
c           the bidiagonal matrix for subroutine QRBD.  Locations N1+1
c           thru 2*N1 will save info from one call to the next of
c           H12.
c
      integer I, IPASS, J, K, L, M, MDA, MDB, M1
      integer N, NB, N1, NP1, NS, NSP1
      double precision A(MDA,N1),B(MDB,NB), S(N1)
      double precision ONE, T, WORK(N1,2), ZERO
      parameter(ONE = 1.0d0, ZERO = 0.0d0)
c
c                             BEGIN.. SPECIAL FOR ZERO ROWS AND COLS.   
c   
c  PACK THE NONZERO COLS TO THE LEFT 
c   
      N=N1  
      IF (N.LE.0.OR.M1.LE.0) RETURN 
      J=N   
   10 CONTINUE  
         DO 20 I=1,M1  
            IF (A(I,J) .ne. ZERO) go to 50
   20    CONTINUE  
c   
c  COL J  IS ZERO. EXCHANGE IT WITH COL N.   
c   
         IF (J .ne. N) then
            DO 30 I=1,M1  
   30       A(I,J)=A(I,N) 
         endif
         A(1,N)=J  
         N=N-1 
   50    CONTINUE  
         J=J-1 
      IF (J.GE.1) GO TO 10  
c
c  IF N=0 THEN A IS ENTIRELY ZERO AND SVD    
c  COMPUTATION CAN BE SKIPPED    
c
      NS=0  
      IF (N.EQ.0) GO TO 240 
c
c  PACK NONZERO ROWS TO THE TOP.
c  QUIT PACKING IF FIND N NONZERO ROWS   
c
      I=1   
      M=M1  
   60 IF (I.GT.N.OR.I.GE.M) GO TO 150   
      IF (A(I,I)) 90,70,90  
   70     DO 80 J=1,N   
          IF (A(I,J)) 90,80,90  
   80     CONTINUE  
      GO TO 100 
   90 I=I+1 
      GO TO 60
c  
c  ROW I IS ZERO.  EXCHANGE ROWS I AND M.
c
  100 IF(NB.LE.0) GO TO 115 
          DO 110 J=1,NB 
          T=B(I,J)  
          B(I,J)=B(M,J) 
  110     B(M,J)=T  
  115     DO 120 J=1,N  
  120     A(I,J)=A(M,J) 
      IF (M.GT.N) GO TO 140 
          DO 130 J=1,N  
  130     A(M,J)=ZERO   
  140 CONTINUE
c
c  EXCHANGE IS FINISHED  
c
      M=M-1 
      GO TO 60  

  150 CONTINUE  
c
c  END.. SPECIAL FOR ZERO ROWS AND COLUMNS   
c  BEGIN.. SVD ALGORITHM 
c
c     METHOD..  
c     (1)     REDUCE THE MATRIX TO UPPER BIDIAGONAL FORM WITH   
c     HOUSEHOLDER TRANSFORMATIONS.  
c          H(N)...H(1)AQ(1)...Q(N-2) = (D**T,0)**T  
c     WHERE D IS UPPER BIDIAGONAL.  
c   
c     (2)     APPLY H(N)...H(1) TO B.  HERE H(N)...H(1)*B REPLACES B    
c     IN STORAGE.   
c   
c     (3)     THE MATRIX PRODUCT W= Q(1)...Q(N-2) OVERWRITES THE FIRST  
c     N ROWS OF A IN STORAGE.   
c   
c     (4)     AN SVD FOR D IS COMPUTED.  HERE K ROTATIONS RI AND PI ARE 
c     COMPUTED SO THAT  
c          RK...R1*D*P1**(T)...PK**(T) = DIAG(S1,...,SM)    
c     TO WORKING ACCURACY.  THE SI ARE NONNEGATIVE AND NONINCREASING.   
c     HERE RK...R1*B OVERWRITES B IN STORAGE WHILE  
c     A*P1**(T)...PK**(T)  OVERWRITES A IN STORAGE. 
c   
c     (5)     IT FOLLOWS THAT,WITH THE PROPER DEFINITIONS,  
c     U**(T)*B OVERWRITES B, WHILE V OVERWRITES THE FIRST N ROW AND     
c     COLUMNS OF A.     
c   
      L=min(M,N)   
c
c  THE FOLLOWING LOOP REDUCES A TO UPPER BIDIAGONAL AND  
c  ALSO APPLIES THE PREMULTIPLYING TRANSFORMATIONS TO B.     
c   
          DO 170 J=1,L  
          IF (J.GE.M) GO TO 160     
          CALL H12 (1,J,J+1,M,A(1,J),1,T,A(1,J+1),1,MDA,N-J)
          CALL H12 (2,J,J+1,M,A(1,J),1,T,B,1,MDB,NB)
  160     IF (J.GE.N-1) GO TO 170   
          CALL H12 (1,J+1,J+2,N,A(J,1),MDA,work(J,2),A(J+1,1),MDA,1,M-J)   
  170     CONTINUE  
c   
c  COPY THE BIDIAGONAL MATRIX INTO S() and WORK() FOR QRBD.   
c  1986 Jan 8. C. L. Lawson. Changed N to L in following 2 statements.
c
      IF (L.EQ.1) GO TO 190 
          DO 180 J=2,L  
          S(J)=A(J,J) 
  180     WORK(J,1)=A(J-1,J)   
  190 S(1)=A(1,1)     

      NS=N  
      IF (M.GE.N) GO TO 200 
      NS=M+1
      S(NS)=ZERO  
      WORK(NS,1)=A(M,M+1)  
  200 CONTINUE  
c   
c  CONSTRUCT THE EXPLICIT N BY N PRODUCT MATRIX, W=Q1*Q2*...*QL*I    
c  IN THE ARRAY A(). 
c   
          DO 230 K=1,N  
          I=N+1-K   
          IF (I .GT. min(M,N-2)) GO TO 210     
          CALL H12 (2,I+1,I+2,N,A(I,1),MDA,WORK(I,2),A(1,I+1),1,MDA,N-I)   
  210         DO 220 J=1,N  
  220         A(I,J)=ZERO   
  230     A(I,I)=ONE    
c   
c  COMPUTE THE SVD OF THE BIDIAGONAL MATRIX 
c   
      CALL QRBD (IPASS,S(1),WORK(1,1),NS,A,MDA,N,B,MDB,NB)   

      if(IPASS .eq. 2) then
         write (*,'(/a)')
     *      ' FULL ACCURACY NOT ATTAINED IN BIDIAGONAL SVD'
      endif

  240 CONTINUE  
      IF (NS.GE.N) GO TO 260
      NSP1=NS+1 
          DO 250 J=NSP1,N   
  250     S(J)=ZERO   
  260 CONTINUE  
      IF (N.EQ.N1) RETURN   
      NP1=N+1
c
c  MOVE RECORD OF PERMUTATIONS AND STORE ZEROS 
c 
          DO 280 J=NP1,N1   
          S(J)=A(1,J) 
              DO 270 I=1,N  
  270         A(I,J)=ZERO   
  280     CONTINUE
c
c  PERMUTE ROWS AND SET ZERO SINGULAR VALUES.
c
          DO 300 K=NP1,N1   
          I=S(K)  
          S(K)=ZERO   
              DO 290 J=1,N1 
              A(K,J)=A(I,J) 
  290         A(I,J)=ZERO   
          A(I,K)=ONE    
  300     CONTINUE  
c
c  END.. SPECIAL FOR ZERO ROWS AND COLUMNS   
c
      RETURN
      END   
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Discussion:
c
c    This FORTRAN77 version is made available for cases where the
c    FORTRAN90 version cannot be used.
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
