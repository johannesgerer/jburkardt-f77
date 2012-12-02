      program main

c*********************************************************************72
c
cc MAIN is the main program for NMS_PRB.
c
c  Discussion:
c
c    NMS_PRB calls the NMS tests.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'NMS_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the NMS library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
      call test06 ( )
      call test07 ( )
      call test08 ( )
      call test09 ( )

      call test10 ( )
      call test11 ( )
      call test12 ( )
      call test13 ( )
      call test14 ( )
      call test15 ( )
      call test16 ( )
      call test17 ( )
      call test18 ( )
      call test19 ( )

      call test20 ( )
      call test21 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'NMS_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests DCFT2D.
c
C EXAMPLE 11.8:  Plot image and transform of 8 by 8 unit source
C in 64 by 64 (otherwise zero) array.
C
      parameter (n=64) 
      double precision  a(n,n),a2(n,n),err,ermax,data
      complex*16  image(n,n),image2(n,n),w(3*n+8)
      logical   forwd

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  DCFT2D computes a 2D Fourier transform.'
      write ( *, '(a)' ) ' '

      WRITE (*,*) 'COMPUTING...'
      LDA = N
c
C Set up data. IMAGE is original. IMAGE2 is scaled by (-1)**(I+J)
C to place Fourier coefficients in the correct place for viewing.
      DO I = 1,N
        DO J = 1,N
          A(I,J) = 0.0D0
          IF (I.GE.(N/2-4) .AND. I.LE.(N/2+4) .AND. J.GE.(N/2-4)
     *         .AND. J.LE.(N/2+4)) A(I,J) = 1.0D0
          IMAGE(I,J) = DCMPLX(A(I,J),0.0D0)
          IMAGE2(I,J) = IMAGE(I,J)*(-1)**(I-1+J-1)
        end do
      end do
C
C Plot image  with axonometric plotting routine.
C
c      CALL AXNPLT (A,LDA,N)
C
C Compute Fourier transform of IMAGE and IMAGE2
C
      FORWD = .TRUE.
      CALL DCFT2D (N,IMAGE,LDA,W,FORWD)
      CALL DCFT2D (N,IMAGE2,LDA,W,FORWD)
C
C Compute magnitude of components of transforms.
C Actual transforms are unscaled and need to be divided by N*N
C to be correct.
C
      DO 2 I = 1,N
        DO 2 J = 1,N
           A(I,J) = ABS(IMAGE(I,J))
           A2(I,J) = ABS(IMAGE2(I,J))
    2 CONTINUE
C
C PLOT MODULUS OF TRANSFORM
C
c    CALL AXNPLT (A,LDA,N)
c      CALL AXNPLT (A2,LDA,N)
C
C Compute inverse transform of image and see if it agrees with original.
C
      FORWD = .FALSE.
      CALL DCFT2D (N,IMAGE,LDA,W,FORWD)
      ERMAX = 0.0D0
      DO 3 I = 1,N
        DO 3 J = 1,N
           DATA = 0.0D0
           IF (I.GE.(N/2-4) .AND. I.LE.(N/2+4) .AND. J.GE.(N/2-4)
     *         .AND. J.LE.(N/2+4)) DATA = 1.0D0
           ERR = ABS( DATA-ABS(IMAGE(I,J))/(N*N) )
           IF (ERR .GT. ERMAX) ERMAX = ERR
    3 CONTINUE
      WRITE (*,*) 'DCFT2D RESULTS (EX 11.8: PLOTS HAVE BEEN SKIPPED)'
      WRITE (*,'(2X,A,1X,D18.10)') 'MAXIMUM ERROR IS ',ERMAX

      WRITE (*,*)
      WRITE (*,*) 'REFERENCE RESULTS FROM IBM PC/AT'
      WRITE (*,*) ' MAXIMUM ERROR IS    0.1885421720E-15'
      RETURN
      END
      SUBROUTINE AXNPLT(A,LDA,N)

c*********************************************************************72
c
cc AXNPLT creates an axonometric plot of a 2d array.
C
C       N: SIZE OF A, N.LE.350
C       LDA: DIMENSION OF FIRST COMPONENT OF A IN CALLING PROGRAM
C
C       ASSUMPTIONS:    0.0 .LE. A(I,J) .LE 1.0
C
      INTEGER LDA,J,K,IDAMAX,N,IRET,I
      double precision A(LDA,*)
      double precision WX(350)
      double precision WY(350),X(350),Y(350),YMX(350),X0,DEL,YMAX
      double precision YM,YX
C
C  FIND THE LARGEST Y SO CAN SCALE
c
      do j=1,n
        k=idamax(n,a(1,j),1)
        ymx(j)=a(k,j)
      end do

      k=idamax(n,ymx,1)
      ymax=ymx(k)
      x0=0.1
      del=.4
c      call agraf0(iret,0)
c      call minmax(iret,0.0, 1.0, 0.0, 1.0, ym, yx)
c      call howplt(iret,0,1,15)
      do 2 i=1,n
        do j=1,n
           x(j)=(j-1.)/(n-1.)*del+x0+(i-1.)/(n-1.)*del
           y(j)=(a(i,j)/ymax)*del+x0+(i-1.)/(n-1.)*del
        end do
        if(i.eq.1)then
c           call nowait
c           call plot1(iret,n,x,y,wx,wy)
        else
c           if(i.eq.n)call wait
c           call howplt(iret,0,1,15)
c           call plot(iret,n,x,y,wx,wy)
        endif
    2 continue
      return
      END
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests DFZERO.
c
c  Modified:
c
c    04 February 2012
c
      double precision b, c, ae, re, test02_f
      external test02_f

      b   = 2.0d0
      c   = 3.0d0
      ae  = 1.0d-06
      re  = 1.0d-06

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02:'
      write ( *, '(a)' ) '  DFZERO finds the root of a function.'
      write ( *, '(a)' ) ' '

      write (*,'(a)') '  Initial interval: '
      write (*,'(8x,2d18.10)') b, c
      write (*,'(a)') '  Tolerances:       '
      write (*,'(8x,2d18.10)') ae, re

      call dfzero ( test02_f, b, c, c, re, ae, iflag)

      write (*,'(a)') ' '
      write (*,'(a)') '  DFZERO results'
      if (iflag .ne. 1) write (*,'(a)') '  Error code =', iflag
      write (*,'(a19,2x,d18.10)') '  Estimate of zero =', b
      write (*,'(a19,2x,d18.10)') '  Function value =  ', test02_f(b)

      write (*,'(a)') ' '
      write (*,'(a)') '  Reference results:'
      write (*,'(a)') '  Estimate of zero =   0.2094551435e+01'    
      write (*,'(a)') '  Function value =    -0.5157851953e-06'

      return
      end
      double precision function test02_f ( x )

c*********************************************************************72
c
cc TEST02_F is a function needed for TEST02_DFZERO.
c
c  Modified:
c
c    04 February 2012
c
      double precision x

      test02_f = x * (x*x-2.0d0) - 5.0d0

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests DGEFS
c
c  Modified:
c
c    04 February 2012
c
      integer   lda
      parameter (lda=10)
      double precision  a(lda,lda), b(lda), work(lda), rcond
      integer   iwork(lda), i, j, n, itask, ind

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03:'
      write ( *, '(a)' ) '  DGEFS factors and solves a linear system.'
      write ( *, '(a)' ) ' '
C
C Set up problem
c
      n      =  3
      itask  =  1
      a(1,1) = 10.0d0
      a(2,1) = -3.0d0
      a(3,1) =  5.0d0
      a(1,2) = -7.0d0
      a(2,2) =  2.0d0
      a(3,2) = -1.0d0
      a(1,3) =  0.0d0
      a(2,3) =  6.0d0
      a(3,3) =  5.0d0
      b(1)   =  7.0d0
      b(2)   =  4.0d0
      b(3)   =  6.0d0
C
C Print problem information
C
      write (*,*) '  Coefficient matrix ='
      do i = 1,n
         write (*,800) (a(i,j), j = 1,n)
      end do
      write (*,*) '  Right-hand side ='
      write (*,800) (b(j), j = 1,n)
C
C Solve linear system
c
      call dgefs (a, lda, n, b, itask, ind, work, iwork, rcond)
C
C Print results
c
      write (*,*)
      write (*,*) '  DGEFS results'
      IF (IND .EQ. -10) THEN
         WRITE (*,*) ' ERROR CODE =', IND
      ELSE IF (IND .LT.0) THEN
         WRITE (*,*) '  ERROR CODE =', IND
         STOP
      ELSE
         WRITE (*,*) '  Number of accurate digits =', IND
      END IF

      WRITE (*,*) '  SOLUTION ='
      WRITE (*,800) (B(J), J = 1,N)

      WRITE (*,*)
      WRITE (*,*) '  REFERENCE RESULTS FROM IBM PC/AT '
      WRITE (*,*) '   NUMBER OF ACCURATE DIGITS =          14'
      WRITE (*,*) '   SOLUTION ='
      WRITE (*,*) '     0.000000000000E+00     -0.100000000000E+01',    
     *'      0.100000000000E+01'

      return
800   FORMAT (4X, D20.12, 4X, D20.12, 4X, D20.12)
      END
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests DUNI
c
c  Modified:
c
c    04 February 2012
c
      DOUBLE PRECISION  U,DUNI,DUSTAR,USEED
      INTEGER ISEED,I

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  DUNI computes uniform random numbers.'
c
C SET INITIAL SEED
c
      ISEED = 305
      USEED = DUSTAR(ISEED)
c
C USTART RETURNS FLOATING ECHO OF ISEED 
c
      WRITE (*,*) '  DUNI RESULTS '
      WRITE (*,'(5X,I5,5X,D15.8)')ISEED, USEED
      DO I = 1,1000
        U = DUNI()
      end do
      WRITE (*,'(12X,D25.15)') U
      WRITE (*,*) '  REFERENCE RESULTS FROM AN IBM P/C'
      WRITE (*,*) '        305      0.30500000E+03' 
      WRITE (*,*) '                 0.241576136599321E+00'

      return
      END
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests DCFFTI and DCFFTF.
c
c  Modified:
c
c    04 February 2012
c
C Using complex discrete Fourier transform, find the approximate Fourier
C coefficients to Runge's function on [-1,1] with N=16 and N=17. 
C 
       DOUBLE PRECISION  WSAVE(150),RUNGE,X,X0,PI,DEL,F,XJ
       COMPLEX*16 COEFF(0:16), SQTM1
C
C Note: Complex Double Precision is not Fortran Standard and many
C       compiler vendors may fail to implement it or implement it
C       in some other manner.  If the COMPLEX*16 declaration does
C       not work on your compiler we recommend that you check
C       with the vendor's instructions and modify the above accordingly.
C
C Note 0 subscript, which makes indexing easier, allowed in Fortran 77.
C
C Arithmetic statement function for Runge's function. 
       RUNGE(X) = 1.0D0/(1.0D0+25.0D0*X*X) 

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  DCFFTI initializes a complex FFT.'
      write ( *, '(a)' ) '  DCFFTF computes it.'

       X0 = -1.0D0 
       PI = ASIN(1.0D0)*2.0D0 
       SQTM1 = CMPLX(0.0D0,-1.0D0)

       DO 10 N = 16,17 
          CALL DCFFTI (N,WSAVE) 
C Function assumed to be periodic on [-1,1], of length 2.
          DEL = 2.0D0/N 
          F = 2.0D0*PI/(N*DEL) 
          DO 1 J = 0,N-1
C First sample point at -1, last at 1-DEL
              XJ = (-1.0D0) + J*DEL 
              COEFF(J) = CMPLX(RUNGE(XJ),0.0D0) 
    1     CONTINUE 
          CALL DCFFTF (N,COEFF,WSAVE) 
C Returned coefficients must be divided by N for correct normaliziation.
C 
C Note repetition after N/2 in original coefficients. Scaling because
C X0 not at origin destroys this to some extent. 
C
          WRITE (*,*) ' DCFFTF RESULTS FOR N = ' ,N 
          WRITE (*,'(A,1X,2D17.8)') ' CZERO = ',COEFF(0)/N*2.0D0
          WRITE (*,*)  
     *' J          OUTPUT FROM DCFFTF,         SCALED COEFFICIENTS'   
          DO 11 J = 1,N-1
              WRITE (*,'(I5,2D17.8,5X,2D17.8)') 
     *              J, COEFF(J), EXP(-SQTM1*J*F*X0) * COEFF(J)/N *2
   11     CONTINUE 
          WRITE (*,*)

   10  CONTINUE 
       WRITE (*,*) 'REFERENCE RESULTS (PARTIAL) FROM IBM PC/AT'
       WRITE (*,*) ' DCFFTF RESULTS FOR N =           17'
       WRITE (*,*) ' CZERO =  (0.54916124E+00, 0.00000000E+00)'
       WRITE (*,*) '    .  '
       WRITE (*,*) '    .  '
       WRITE (*,*) '  10  -0.60345452E-01  -0.37837322E-16 '  
       WRITE (*,*) '  11   0.11248162E+00  -0.47424091E-16 '  
       WRITE (*,*) '  12  -0.23457294E+00   0.73973161E-16 '  
       WRITE (*,*) '  13   0.42197537E+00  -0.25268310E-17 '  
       WRITE (*,*) '  14  -0.82441735E+00  -0.14952867E-16 '
       WRITE (*,*) '  15   0.14919178E+01   0.31242621E-16 '
       WRITE (*,*) '  16  -0.29260639E+01  -0.69144627E-16 '
      return
      END  
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 tests DRNOR.
c
c  Modified:
c
c    06 February 2012
c
C HISTOGRAM FOR DRNOR
C
      DOUBLE PRECISION A,B
      PARAMETER(NBINS=32,A=-3.0D0,B=3.0D0)
      INTEGER ISEED,I,J,H(NBINS),NR,test06_INBIN
      DOUBLE PRECISION  R,DRNOR,DSTART,RSEED,WIDTH

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) '  DRNOR computes normal random numbers.'

      ISEED = 305
      RSEED = DSTART(ISEED)
      WIDTH=(B-A)/(NBINS-2)
CCCCCC  100 WRITE (*,*) 'EX 10.1: ENTER NUMBER OF NORMALS : '
      DO 200 I = 1,NBINS
          H(I)=0
  200 CONTINUE
CCCCC      READ (*,*) NR
      WRITE (*,*) 'RUNNING 10,000 NORMALS INTO 32 BINS...'
      NR = 10000
      IF (NR .LE. 0) return
      DO I = 1,NR 
          R = DRNOR()
          J = test06_INBIN (A,B,NBINS,WIDTH,R)
          H(J) = H(J)+1
      end do

      WRITE (*,*) 'HISTOGRAM FOR DRNOR: NUMBER IN BIN 1,...,32'
      WRITE (*,*) '   (-INFINITY,-3],(-3,-2.8],...,(2.8,3],(3,INFINITY)'
      WRITE (*,*) ' (VALUES ARE SLIGHTLY COMPUTER DEPENDENT)'
      WRITE (*,'(9I8)') (H(I),I=1,NBINS)

      WRITE (*,*)
      WRITE (*,*)
     *   'REFERENCE RESULTS FROM IBM PC/AT '
      WRITE(*,*)'     16      14      21      32      46      96     135
     *      198     292'
      WRITE(*,*)'    334     454     549     611     665     743     801
     *      785     751'
      WRITE(*,*)'    747     634     503     421     351     274     184
     *      124      90'
      WRITE(*,*)'     53      35      19      13       9'
      return
      END 
      INTEGER FUNCTION test06_INBIN (XMIN,XMAX,NBINS,WIDTH,DATA)

c*********************************************************************72
c
cc TEST06_INBIN bins data for TEST06.
C
C THIS FUNCTION TAKES A DOUBLE PRECISION VALUE IN DATA, AND FINDS
C THE CORRECT BIN FOR IT.  VALES BELOW XMIN COME BACK
C IN 1.  VALUES ABOVE XMAX COME BACK IN NBINS.
C
      integer nbins 
      double precision xmin,xmax,width,data

      if (data .lt. xmin) then
         test06_inbin = 1
      else if (data .ge. xmax) then
         test06_inbin = nbins
      else
         test06_inbin = 2 + (data-xmin)/width
      endif

      return
      end
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07 tests DPCHEZ.
c
c  Modified:
c
c    06 February 2012
c
      double precision x(21), f(21), d(21), wk(42), fe(101), 
     *xe(101), fd(101), r, rp, u, error, errord, a, b, q, dpchqa,
     *xans(4), fans(4), erans(4), edans(4), qans
      logical spline
      data xans/-0.3d-01,-0.2d-01,-0.1d-01,0.0d0/
      data fans/0.971920d0,0.986880d0,0.996560d0,1.0d0/
      data erans/-0.6075110024d-02,-0.3219009901d-02,
     *-0.946234414d-03,0.0d0/
      data edans/0.2932883472d0,0.2677039506d0,0.1744906562d0,0.0d0/
      data qans,ians/0.274679262701527d0,0/
C
C Arithmetic statement functions for Runge's function and derivative.
C
      R(U)  = 1.0D0/(1.0D0+25.0D0*U*U)
      RP(U) = -50.0D0*U*R(U)**2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07'
      write ( *, '(a)' ) 
     &  '  DPCHEZ sets up piecewise cubic Hermite splines.'
      write ( *, '(a)' ) ' '
C
C Compute Runge's function at 21 points in [-1,1].
C
      DO 1 I=1,21
          X(I) = -1.0D0 + (I-1.0D0)/10.0D0
          F(I) = R(X(I))
    1 CONTINUE
      N = 21
      NWK = 42 
      SPLINE = .FALSE.
C
C Compute cubic Hermite interpolant because SPLINE is .FALSE.
C
      CALL DPCHEZ (N,X,F,D,SPLINE,WK,NWK,IERR)
      IF (IERR .LT. 0) THEN
          WRITE (*,*) 'AN ERROR CALLING DPCHEZ, IERR= ',IERR
          STOP
      ENDIF

      NE = 101
C
C Evaluate interpolant and derivative at 101 points from -1 to 0.
C
      DO I=1,NE
        XE(I) = -1.0D0 + (I-1.0D0)/(NE-1.0D0)
      end do

      CALL DPCHEV (N,X,F,D,NE,XE,FE,FD,IERR)
      IF (IERR .NE. 0) THEN
          WRITE (*,*) 'AN ERROR CALLING DPCHEV, IERR= ',IERR
          STOP
      ENDIF

      DO 4 I=1,NE
          ERROR = FE(I) - R(XE(I))
          ERRORD = FD(I) - RP(XE(I))
          WRITE (*,3) XE(I),FE(I),ERROR,ERRORD
    3     FORMAT(1X,D17.10,3X,D17.10,3X,D17.10,3X,D17.10)
    4 CONTINUE
C
C Compute integral over the interval [0,1] 
C
      A = 0.0D0
      B = 1.0D0
      Q = DPCHQA (N,X,F,D,A,B,IERR)
      WRITE (*,'(2X,A,D20.12,3X,A,I5)') 'INTEGRAL FROM 0 TO 1 = ',
     *Q,' IERR = ',IERR

      WRITE (*,*) 
      WRITE (*,*) 
     *   ' REFERENCE RESULTS FROM IBM PC/AT OF PRECEDING 5 LINES '
      DO I=1,4
          WRITE (*,5) XANS(I),FANS(I),ERANS(I),EDANS(I)
    5     FORMAT(1X,D17.10,3X,D17.10,3X,D17.10,3X,D17.10)
      end do

      WRITE (*,'(2X,A,D20.12,3X,A,I5)') 'INTEGRAL FROM 0 TO 1 = ',
     *QANS,' IERR = ',IANS           
      return
      END
      subroutine test08 ( )

c*********************************************************************72
c
cc TEST08 tests DFMIN.
c
      DOUBLE PRECISION A
      double precision B, XSTAR, TOL, DFMIN, test08_f
      EXTERNAL test08_f

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08'
      write ( *, '(a)' ) '  DFMIN minimizes a function.'
      write ( *, '(a)' ) ' '

      A     = 0.1D0
      B     = 0.9D0
      TOL   = 1.0D-5
      XSTAR = DFMIN ( A, B, test08_f, TOL )
      WRITE (*,*) 'DFMIN RESULTS'
      WRITE (*,'(2X,A,D18.10)') ' XSTAR =', XSTAR
      WRITE (*,*)
      WRITE (*,*) 'REFERENCE RESULTS FROM IBM PC/AT'
      WRITE (*,*) '  XSTAR =  0.8164961769E+00'    
      return
      END
      DOUBLE PRECISION FUNCTION test08_f(X)

c*********************************************************************72
c
cc TEST08_F is the function to be minimized.
c
      DOUBLE PRECISION X
      test08_f = X*(X*X - 2.0D0) - 5.0D0
      RETURN
      END
      subroutine test09 ( )

c*********************************************************************72
c
cc TEST09 tests DSNQE.
c
      PARAMETER   (N = 2, LW = 19)
      DOUBLE PRECISION   X(N), FVEC(N), W(LW), TOL
      EXTERNAL    test09_f

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST09'
      write ( *, '(a)' ) 
     &  '  DNSQE solves a system of nonlinear equations.'
      write ( *, '(a)' ) ' '

      TOL  = 1.0D-6
      X(1) = 2.0D0
      X(2) = 3.0D0
      WRITE (*,800) X(1), X(2)
      IOPT = 2
      NPRINT = 0
C
C SOLVE NONLINEAR EQUATIONS
C
      CALL DNSQE ( test09_f, test09_f, IOPT, N, X, FVEC, TOL, NPRINT, 
     &  INFO, W, LW)
C
C PRINT RESULTS
C
      WRITE (*,*)
      WRITE (*,*) '  DNSQE RESULTS'
      IF (INFO .NE. 1) WRITE (*,810) INFO
      WRITE (*,820) X(1), X(2)
      WRITE (*,830) FVEC(1), FVEC(2)

      WRITE (*,*)
      WRITE (*,*) 'REFERENCE RESULTS FROM IBM PC/AT'
      WRITE (*,*) ' ESTIMATE OF SOLUTION '
      WRITE (*,*) '     0.199999999997E+01  0.100000000001E+01'
      WRITE (*,*) ' VALUES OF NONLINEAR FUNCTIONS'
      WRITE (*,*) '    -0.415920631269E-10 -0.940900690694E-10'

      return
800   FORMAT (' INITIAL GUESS', /, 4X,2D20.12)
810   FORMAT (' INFO =', I3)
820   FORMAT (' ESTIMATE OF SOLUTION ', /, 4X,2D20.12)
830   FORMAT (' VALUES OF NONLINEAR FUNCTIONS', /, 4X,2D20.12)
      END
      SUBROUTINE test09_f (N, X, FVEC, IFLAG)

c*********************************************************************72
c
cc TEST09 evaluates the nonlinear equations.
c
      DOUBLE PRECISION   X(N), FVEC(N)
C
C COMPUTE NONLINEAR FUNCTIONS
C
      FVEC(1) = X(1)*X(2) - X(2)**3 - 1.0D0
      FVEC(2) = X(1)**2*X(2) + X(2) - 5.0D0

      RETURN
      END
      subroutine test10 ( )

c*********************************************************************72
c
cc TEST10 tests DQRLS.
c
      PARAMETER (MM = 5, NN = 3)
      DOUBLE PRECISION A(MM,NN), B(MM), X(NN), QRAUX(NN), WORK(NN), TOL
      INTEGER   JPVT(NN)
      DATA      B / 1.0D0, 2.3D0, 4.6D0, 3.1D0, 1.2D0/

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST10'
      write ( *, '(a)' ) '  DQRLS solves a least squares problem.'
      write ( *, '(a)' ) ' '
C
C SET UP LEAST-SQUARES PROBLEM
C     QUADRATIC MODEL, EQUALLY-SPACED POINTS
C
      M = 5
      N = 3
      DO I = 1,M
         A(I,1) = 1.0D0
         DO J = 2,N
            A(I,J) = A(I,J-1)*I
         end do
      end do

      TOL = 1.D-15
      ITASK = 1
      WRITE (*,*)   ' COEFFICIENT MATRIX'
      WRITE (*,800) ((A(I,J), J = 1,N), I = 1,M)
      WRITE (*,*)   ' RIGHT-HAND SIDE'
      WRITE (*,810) (B(I), I = 1,M)
C
C
C SOLVE LEAST-SQUARES PROBLEM
C
      CALL DQRLS (A, MM, M, N, TOL, KR, B, X, B, JPVT, QRAUX, WORK,
     *            ITASK, IND)
      IF (IND .NE. 0) THEN
         WRITE (*,*) ' ERROR CODE =', IND
         STOP
      END IF
      WRITE (*,*)    ' RANK OF MATRIX =', KR
C
C PRINT RESULTS
C
      WRITE (*,*)   ' PARAMETERS'
      WRITE (*,800) (X(J), J = 1,N)
      WRITE (*,*)   ' RESIDUALS'
      WRITE (*,810) (B(I), I = 1,M)
      WRITE (*,*)
      WRITE (*,*) 'REFERENCE RESULTS FROM IBM PC/AT'
      WRITE (*,*) ' RANK OF MATRIX =           3'
      WRITE (*,*) ' PARAMETERS'
      WRITE (*,*) '   -0.30200000E+01   0.44914286E+01  -0.72857143E+00'
      WRITE (*,*) ' RESIDUALS'
      WRITE (*,*) 
     *'   0.2571429E+00 -0.7485714E+00  0.7028571E+00 -0.1885714E+00
     *-0.2285714E-01'

      return
800   FORMAT (2X,3D17.8)
810   FORMAT (2X,5D15.7)
      END
      subroutine test11 ( )

c*********************************************************************72
c
cc TEST11 tests DFZERO.
c
      DOUBLE PRECISION     B, C, AE, RE, test11_f
      EXTERNAL test11_f

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST11'
      write ( *, '(a)' ) '  DFZERO solves a nonlinear equation.'
      write ( *, '(a)' ) ' '

      B   = 2.0D0
      C   = 3.0D0
      AE  = 1.0D-06
      RE  = 1.0D-06
      WRITE (*,*) ' INITIAL INTERVAL: '
      WRITE (*,'(8X,2D18.10)') B, C
      WRITE (*,*) ' TOLERANCES:       '
      WRITE (*,'(8X,2D18.10)') AE, RE

      CALL DFZERO ( test11_f, B, C, C, RE, AE, IFLAG)

      WRITE (*,*)
      WRITE (*,*) '  DFZERO RESULTS'
      IF (IFLAG .NE. 1) WRITE (*,*) '  ERROR CODE =', IFLAG
      WRITE (*,'(A19,2X,D18.10)') '  ESTIMATE OF ZERO =', B
      WRITE (*,'(A19,2X,D18.10)') '  FUNCTION VALUE =', test11_f (B)

      WRITE (*,*)
      WRITE (*,*) '  REFERENCE RESULTS FROM IBM PC/AT'
      WRITE (*,*) '   ESTIMATE OF ZERO =   0.2094551435E+01'    
      WRITE (*,*) '     FUNCTION VALUE =  -0.5157851953E-06'

      return
      END
      DOUBLE PRECISION FUNCTION test11_f ( X )

c*********************************************************************72
c
cc TEST11_F evaluates a nonlinear equation.
c
      DOUBLE PRECISION X

      test11_f = X * (X*X-2.0D0) - 5.0D0

      RETURN
      END
      subroutine test12 ( )

c*********************************************************************72
c
cc TEST12 tests UNCMND.
c
C MAIN PROGRAM FOR NONLINEAR LEAST-SQUARES DATA FITTING
C
      PARAMETER       (N = 2, LWORK = N*(N+10), MD = 4)
      DOUBLE PRECISION  X0(N), X(N), F, WORK(LWORK), T(MD), B(MD)
      COMMON /EXPDAT/ T, B, M
      EXTERNAL test12_dcalcf

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST12'
      write ( *, '(a)' ) 
     &  '  UNCMND performs nonlinear least squares data fitting.'
      write ( *, '(a)' ) ' '
C
C DATA FOR DATA FITTING
C
      T(1) =  0.0D0
      T(2) =  1.0D0
      T(3) =  2.0D0
      T(4) =  3.0D0
      B(1) = 20.0D0
      B(2) =  9.0D0
      B(3) =  3.0D0
      B(4) =  1.0D0
C
C SPECIFY INITIAL ESTIMATE OF THE SOLUTION
C
      M     = 4
      X0(1) = 1.0D0
      X0(2) = 1.0D0
C
C MINIMIZE FUNCTION
C
      CALL UNCMND (N, X0, test12_dcalcf, X, F, IERROR, WORK, LWORK)
C
C PRINT RESULTS
C
      WRITE (*,*) 'UNCMND FOR NONLINEAR LEAST SQUARES RESULTS'
      IF (IERROR .NE. 0) WRITE (*,*) ' ERROR CODE =', IERROR
      WRITE (*,'(A,1X,D15.8)') ' F(X*) =', F
      WRITE (*,'(1X,A)') ' X* =' 
      WRITE (*,'(5X,D20.12)') (X(I), I = 1,N)

      WRITE (*,*) 
      WRITE (*,*) 
     * 'REFERENCE RESULTS (PARTIAL-LAST 8 LINES) FROM IBM PC/AT'
      WRITE (*,*) '-0.52716742E+01  -0.20203044E+02  -0.57189483E+01'
      WRITE (*,*) ' 0.20000000E+02  -0.20605341E+02  -0.52716742E+01'
      WRITE (*,*) ' 0.20000001E+02  -0.20605341E+02  -0.52716742E+01'    
      WRITE (*,*) ' 0.20000000E+02  -0.20605340E+02  -0.52716742E+01'
      WRITE (*,*) 'UNCMND WARNING -- INFO = 1: PROBABLY CONVERGED, GRADI
     *ENT SMALL'
      WRITE (*,*) 'UNCMND FOR NONLINEAR LEAST SQUARES RESULTS'
      WRITE (*,*) ' ERROR CODE =           1'
      WRITE (*,*) ' F(X*) =  0.91000000E+02'
      WRITE (*,*) ' X* = '
      WRITE (*,*) '     0.200000001134E+02'
      WRITE (*,*) '    -0.206053408470E+02'    

      return
      END
      SUBROUTINE test12_dcalcf (N, X, F)

c*********************************************************************72
c
cc TEST12_DCALCF evaluates the objective function.
c
      DOUBLE PRECISION  X(N), F, T(4), B(4)
      COMMON /EXPDAT/ T, B, M

      F = 0.0D0
      WRITE (*,'(1X,D15.8,2X,D15.8,2X,D15.8)') X(1),X(2),X(3)
      DO J = 1,M
         F = F + (B(J) - X(1)*EXP(X(2)*T(J)))**2.0D0
      end do

      RETURN
      END
      subroutine test13 ( )

c*********************************************************************72
c
cc TEST13 tests DQK15.
c
C  Compute erf(1), i.e. integral of 2/sqrt(pi) * exp(-x*x) from 0 to 1.0
c
      EXTERNAL test13_func
      DOUBLE PRECISION A,B,RESULT,ABSERR,RESABS,RESASC,PI

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST13'
      write ( *, '(a)' ) '  DQK15 estimates an integral.'
      write ( *, '(a)' ) ' '

      A = 0.0D0
      B = 1.0D0
      PI = 4.0D0*ATAN(1.0D0)
      CALL DQK15(test13_func,A,B,RESULT,ABSERR,RESABS,RESASC)
      WRITE(*,*) ' DQK15 ESTIMATE OF ERF(1)'
      WRITE(*,*) ' 2.0/SQRT(PI)*RESULT'
      WRITE(*,'(3X,D20.12,3X,D20.12)') 2.0D0/SQRT(PI)*RESULT,ABSERR

      WRITE(*,*)
      WRITE(*,*) 'REFERENCE RESULTS COMPUTED ON IBM PC/AT '
      WRITE(*,*) '    0.842700792950E+00     0.828974787422E-14'
      return
      END
      DOUBLE PRECISION FUNCTION test13_func (X)

c*********************************************************************72
c
cc TEST13_FUNC evaluates the integrand.
c
      DOUBLE PRECISION X
      test13_func = EXP(-X*X)
      RETURN
      END
      subroutine test14 ( )

c*********************************************************************72
c
cc TEST14 demonstates DSVDC.
c
C  Compute erf(1), i.e. integral of 2/sqrt(pi) * exp(-x*x) from 0 to 1.0
c
C       CENSUS DATA ANALYSIS WITH THE SVD
C       FROM SECTION 8.1
C     From the book "Numerical Methods and Software"
C          by  D. Kahaner, C. Moler, S. Nash
C               Prentice Hall 1988
C
      INTEGER LDX,N,P,LDU,LDV,JOB,INFO,I,J,JB
      PARAMETER(LDX=8,N=8,P=3,LDU=N,LDV=P,JOB=11)
      DOUBLE PRECISION POP(N),Y(N),X(LDX,P),S(P),E(P),
     *U(LDU,LDU),V(LDV,P),W(N),C(P),YEAR,POP80,
     *TOL,RELERR,SUM,R,RSQ,RI
C
C       C CONTAINS COEFFICIENTS OF POLYNOMIAL C(1)*1+C(2)*T+C(3)*T*T
C         T=YEAR (1900 ETC.)
C
      DATA POP/
     *   75.994575D0,
     *   91.972266D0,
     *  105.710620D0,
     *  122.775046D0,
     *  131.669275D0,
     *  150.697361D0,
     *  179.323175D0,
     *  203.235298D0/

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST14'
      write ( *, '(a)' ) 
     &  '  DSVDC computes the singular value decomposition.'
      write ( *, '(a)' ) ' '

      DO 1 I=1,8
C        Y(I)=(1935.0D0-1900D0*(I-1))/1935.0D0
        Y(I)=1900.0D0+(I-1)*10D0
        X(I,1)=1.0D0
        X(I,2)=Y(I)
        X(I,3)=Y(I)**2
    1 CONTINUE

      CALL DSVDC(X,LDX,N,P,S,E,U,LDU,V,LDV,W,JOB,INFO)
C       NOTE: X IS DESTROYED BY ABOVE CALL!!!
C
C       WRITE SINGULAR VALUES (DESCENDING ORDER)
      WRITE(*,*) 'SINGULAR VALUES ARE: '
      WRITE(*,'(2X,D16.8)') (S(I),I=1,P)
      DO 2 I=1,P
        C(I)=0.0D0
    2 CONTINUE
C
C       RELERR REFLECTS NUMBER OF ACCURATE DIGITS IN DATA
C        E.G. 6 DIGITS ==> RELERR=1.D-6, ETC.
C       MAKING RELERR LARGER INCREASES RESIDUALS
      RELERR=1.D-6
      TOL=RELERR*S(1)
C
C       MULTIPLY U-TRANS * POP, AND SOLVE FOR COEFFICIENTS C(I)
C
      DO 60 J=1,P
         IF(S(J).LE.TOL)GO TO 60
         SUM=0.0D0
         DO I=1,N
            SUM=SUM+U(I,J)*POP(I)
         end do
         SUM=SUM/S(J)
         DO 50 I=1,P
            C(I)=C(I)+SUM*V(I,J)
   50    CONTINUE
   60 CONTINUE
C
      WRITE(*,*) 'COEFFICIENTS (ASSUMING DATA GOOD TO 6 DIGITS) ARE:'
      WRITE(*,'(2X,D16.8)') (C(I),I=1,P)
      WRITE(*,*)
C
C       EVALUATE MODEL (HORNER'S RULE) AND RESIDUALS AT YEAR =1900,...,1980
C
      RSQ=0.0D0
      DO 75 I=1,9
         YEAR=1900.0D0+(I-1)*10.0D0
         POP80=0.0D0
         DO JB=1,P
            J=P+1-JB
            POP80=YEAR*POP80+C(J)
         end do
         IF(I.LT.9) THEN
            RI=POP(I)-POP80
            RSQ=RSQ+RI*RI
            WRITE(*,'(A,I6,A)') 
     *       ' FOR YEAR',IDINT(YEAR),
     *' POP ESTM, MEAS, AND RESIDUAL ARE'
            WRITE(*,'(3D16.8)')POP80,POP(I),RI
         ELSE
            WRITE(*,'(A,I6,A,D16.8)') 
     *        ' FOR YEAR',IDINT(YEAR),' POP ESTMATE IS ',POP80
         ENDIF
   75 CONTINUE

      R=SQRT(RSQ)
      WRITE(*,'(1X,A,D16.8)')
     *'SQUARE ROOT OF RESIDUAL SUM OF SQUARES IS: ',R

      WRITE(*,*)
      WRITE(*,*)'REFERENCE RESULTS (PARTIAL) FROM IBM PC/AT'
      WRITE(*,*)'SINGULAR VALUES ARE: '
      WRITE(*,*)'   0.10594723E+08'
      WRITE(*,*)'   0.64774566E+02'
      WRITE(*,*)'   0.34620247E-03'
      WRITE(*,*)'COEFFICIENTS (ASSUMING DATA GOOD TO 6 DIGITS) ARE:'
      WRITE(*,*)'  -0.16714353E-02'
      WRITE(*,*)'  -0.16169731E+01'
      WRITE(*,*)'   0.87095700E-03'
      WRITE(*,*)'      .  '
      WRITE(*,*)'      .  '
      WRITE(*,*)'      .  '    
      WRITE(*,*)'FOR YEAR  1980 POP ESTMATE IS    0.21289144E+03'
      WRITE(*,*)'SQUARE ROOT OF RESIDUAL SUM OF SQUARES IS:  ',
     *'  0.16096596E+02'
      return
      END
      subroutine test15 ( )

c*********************************************************************72
c
cc TEST15 demonstates DQ1DA.
c 
       DOUBLE PRECISION A, B, EPS, R, E
       double precision f_15
       external f_15

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST15'
      write ( *, '(a)' ) 
     &  '  DQ1DA estimates the integral of a function.'
      write ( *, '(a)' ) ' '

         A = 0.0D0
         B = 1.0D0          
C Set interval endpoints to [0,1]
         EPS = 0.1D-02       
C Set accuracy request for three digits
         CALL DQ1DA (f_15,A,B,EPS,R,E,KF,IFLAG)
         WRITE(*,*)'         INTERVAL: '
         WRITE(*,'(3X,D18.10,3X,D18.10)') A,B
         WRITE(*,*)  
         WRITE(*,*)'         DQ1DA RESULTS: '
         WRITE(*,1) EPS,R,E,KF,IFLAG
    1    FORMAT(3X,D18.10,3X,D18.10,3X,D18.10,3X,I5,2X,I5)
         WRITE(*,*)
         WRITE(*,*)'REFERENCE RESULTS ON IBM PC/AT (CALLED DUNI)'
         WRITE(*,*)'    0.0000000000E+00     0.1000000000E+01',
     *'    0.1000000000E-02     0.4140675161E-01     0.4286574121E-11',
     *'      30      0'
      return
      END
      DOUBLE PRECISION FUNCTION F_15(X)

c*********************************************************************72
c
cc F_15 evaluates the integrand.
c
      double precision x

      f_15 = dsin(2.0d0*x)-dsqrt(x)

      return 
      end
      subroutine test16 ( )

c*********************************************************************72
c
cc TEST16 demonstates DEZFTF.
c 
C Using the D.P. discrete Fourier transform, find the approximate 
C Fourier coefficients to Runge's function on [-1,1] with N=16 and 
C N=17. 
C 
       PARAMETER (MCOEF=17)
       DOUBLE PRECISION  A(MCOEF/2),B(MCOEF/2),R(MCOEF),
     *WSAVE(3*MCOEF+15),TN,ER,PI,DEL,XJ,DFTA(MCOEF/2),
     *DFTB(MCOEF/2),F,C(MCOEF/2),S(MCOEF/2),X,RUNGE,X0,AZERO 
C
C Arithmetic statement function for Runge's function. 
       RUNGE(X) = 1.0D0/(1.0D0+25.0D0*X*X) 

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST16'
      write ( *, '(a)' ) 
     &  '  DEZFTF computes the discrete Fourier transform.'
      write ( *, '(a)' ) ' '

       X0 = -1.0D0 
       PI = ASIN(1.0D0)*2.0D0
C
       DO 10 N = MCOEF-1,MCOEF 
          CALL DEZFTI (N,WSAVE) 
C Function assumed to be periodic on [-1,1], of length 2.
          DEL = 2.0D0/N 
          F = 2.0D0*PI/(N*DEL) 
          DO 1 J = 1,N 
C First sample point at -1, last at 1-DEL
              XJ = (-1.0D0) + (J-1)*DEL 
              R(J) = RUNGE(XJ) 
C Compute sines and cosines to adjust output of DEZFTF to give
C approximate Fourier coefficients.
              IF (J .LE. N/2) THEN 
                  C(J) = COS(J*F*X0) 
                  S(J) = SIN(J*F*X0) 
              END IF 
    1     CONTINUE 
          CALL DEZFTF (N,R,AZERO,A,B,WSAVE) 
C 
C As a convenience this loop can go to N/2. If N is even last B is zero.
c
          DO J = 1,N/2 
              DFTA(J) = A(J)*C(J) - B(J)*S(J) 
              DFTB(J) = A(J)*S(J) + B(J)*C(J)  
          end do

          WRITE (*,'(A,I3,A,D18.10)') ' DEZFTF RESULTS FOR N= ' ,
     *N, ' AZERO = ',AZERO
          WRITE (*,*) '    J          DFTA(J)              DFTB(J) '
          DO 12 J = 1,N/2
              WRITE(*,'(1X,I5,3X,D18.10,3X,D18.10)') J,
     *DFTA(J),DFTB(J)
   12     CONTINUE

          M = 101 
C
          WRITE (*,*) 'FOR BREVITY 101 EVALUATION POINTS OMITTED'
          IF (M.GT.0)GOTO 10
C Evaluate interpolant at 101 points on [-1,1]     
          WRITE (*,*) ' RESULTS FOR N= ',N 
          DO 20 K = 1,M 
              X = -1.0D0 + 2.0D0*(K-1.0D0)/(M-1.0D0) 
              TN = AZERO 
              DO 19 J = 1,N/2 
                  TN = TN + DFTA(J)*COS(J*F*X) + DFTB(J)*SIN(J*F*X) 
   19         CONTINUE 
              ER = TN - RUNGE(X) 
              WRITE (*,'(2X,D15.8,2X,D15.8,2X,D15.8)') X,TN,ER 
   20     CONTINUE 

          WRITE (*,*)

   10  CONTINUE 
       WRITE (*,*)
       WRITE (*,*) 'REFERENCE RESULTS FROM IBM PC/AT'
       WRITE (*,*)
     *  ' DEZFTF RESULTS FOR N=  17    AZERO =    0.274581'    
       WRITE (*,*)'    J         DFTA(J)              DFTB(J) '
       WRITE (*,*)'    1     0.3442428079E+00     0.5030951647E-16'
       WRITE (*,*)'    2     0.1755197395E+00     0.4668315002E-16'
       WRITE (*,*)'    3     0.9699027603E-01     0.3740740688E-16'
       WRITE (*,*)'    4     0.4964416088E-01     0.2403131256E-16'
       WRITE (*,*)'    5     0.2759681609E-01     0.8202373493E-17'
       WRITE (*,*)'    6     0.1323313221E-01     0.4148226049E-17'
       WRITE (*,*)'    7     0.7099464912E-02     0.1053997863E-16'
       WRITE (*,*)'    8     0.1413251200E-02     0.5512226224E-17'

      return
      END
      subroutine test17 ( )

c*********************************************************************72
c
cc TEST17 demonstates DDRIV2.
C
C Notice that NROOT is set to 1
C
      DOUBLE PRECISION  H
      PARAMETER (N=2, H=10.0, NROOT=1, MINT=2, 
     *           LW=N*N+10*N+2*NROOT+204, LIW=23)
      DOUBLE PRECISION  Y(N+1),W(LW),GFUN,EPS,T,TOUT,EWT,MASS
      INTEGER   IW(LIW)
      EXTERNAL  test17_FSUB,test17_GFUN
      DATA      MASS /0.125D0/

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST17'
      write ( *, '(a)' ) 
     &  '  DDRIV2 solves a differential equation.'
      write ( *, '(a)' ) ' '

      EPS = 1.D-5
c
C Set initial point
c
      T = 0.0D0
      TOUT = T
C Set for pure relative error
      EWT = 0.0D0
C Set initial conditions
      Y(1) = H
      Y(2) = 0.0D0
C Set parameter value
      Y(3) = MASS
      MS = 1
      WRITE (*,*) 'DDRIV2 RESULTS'
      WRITE (*,5) T,Y(1),Y(2),MS
    5 FORMAT(2X,D18.10,2X,D18.10,2X,D18.10,1X,I3)
C
   10 CALL DDRIV2 (N,T,Y,test17_FSUB,TOUT,MS,NROOT,EPS,EWT,MINT,W,LW,
     &  IW,LIW,test17_GFUN)

      TOUT = TOUT+0.1D0
      IF (MS .EQ. 5) THEN
        WRITE (*,'(2X,D18.10,2X,D18.10,2X,D18.10,I4)') 
     *        T,Y(1),Y(2),MS
        WRITE (*,'(2X,A,D18.10)') ' <-- Y=0 AT T= ',T
        GOTO 100
      ELSE
        WRITE (*,'(2X,D18.10,2X,D18.10,2X,D18.10,I4)') T,Y(1),Y(2),MS
C Stop if any output code but 1 or 2.
        IF (MS .GT. 2) GOTO 100
      END IF
      GOTO 10
  100 WRITE (*,*)
      WRITE (*,*) 'REFERENCE RESULTS (LAST LINE) FROM IBM PC/AT'
      WRITE (*,*)
     *'   0.2624999995E+01   -0.5648151218E-15   -0.3999999828E+01   5', 
     *'  <-- Y=0 AT T=   0.2624999995E+01'
      return
      END 
      SUBROUTINE test17_FSUB (N,T,Y,YDOT)

c*********************************************************************72
c
cc TEST17_FSUB evaluates right hand sides of equations.
c
        DOUBLE PRECISION T, Y(*), YDOT(*), MASS, G
        DATA G/32.0D0/
c
C Retrieve parameter
c
        MASS = Y(3)

        YDOT(1) = Y(2)
        YDOT(2) = -(G+1.0D0/MASS*Y(2))
        RETURN
      END
      DOUBLE PRECISION FUNCTION test17_GFUN(N,T,Y,IROOT)

c*********************************************************************72
c
cc TEST17_GFUN is used for root finding.
c
C Integration will stop when the function changes sign.
c
      DOUBLE PRECISION Y(*),T
      test17_GFUN = Y(1)
      RETURN
      END
      subroutine test18 ( )

c*********************************************************************72
c
cc TEST18 finds autocorrelation to el Nino data using direct and FFT methods.
C
      INTEGER N
      PARAMETER(N=168)
      DOUBLE PRECISION EL(0:2*N-1),WSAVE(4*(2*N)+15),ACOV(0:N-1),
     *           A(2*N),B(2*N),ACOVR(0:2*N-1),SUM,AZERO
      COMPLEX*16 CEL(0:2*N-1),CORR(0:2*N-1)
      LOGICAL EX

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST18'
      write ( *, '(a)' ) 
     &  '  Find autocorrelation in el Nino data.'
      write ( *, '(a)' ) ' '
C
C  NOTE:  DOUBLE PRECISION COMPLEX IS NOT STANDARD FORTRAN.
C         SOME COMPILER VENDORS MAY FAIL TO IMPLEMENT IT
C         OR MAY IMPLEMENT IN A DIFFERENT MANNER.  IF YOUR
C         COMPILER DOES NOT ACCEPT THE COMPLEX*16 STATEMENT
C         WE RECOMMEND THAT YOU CHECK WITH YOUR VENDOR AND
C         MODIFY THE ABOVE ACCORDINGLY.
C
C needs N locations for el, 2N for acov, cel, corr
C and 4*(2N)+15 for wsave in complex case. N has maximum of 168.
c
      INQUIRE (FILE='elnino.dat',EXIST=EX,ERR=1000)
      IF(.NOT.EX)GOTO 1000
      OPEN (UNIT=8,FILE='elnino.dat',ERR=1000)
C
C Read data, find mean.
c
      SUM = 0.0D0
      DO I = 0,N-1
         READ (8,*) EL(I)
         SUM = SUM + EL(I)
      end do

      DO 101 I = 0,N-1
         EL(I) = EL(I) - SUM/N
C Subtract mean, and add N zeros for either complex or d.p. FFT usage
         CEL(I) = CMPLX(EL(I),0.0D0)
         EL(I+N) = 0.0D0
         CEL(I+N) = 0.0D0
  101 CONTINUE
C
C Direct calculation. Only sum as far as there is data. 
C Simple, but slow.
C      
      DO 110 J = 0,N-1
         ACOV(J) = 0.0D0
         DO 110 M = 0,N-1-J
            ACOV(J) = ACOV(J)+ EL(M) * EL(M+J)
  110 CONTINUE
C
C Write, scaled correlation
      WRITE (*,*) 
     * 'EX 11.6: AUTOCORRELATION (DIRECT) OUTPUT SUPRESSED'
CCCCC      WRITE (*,'(5D14.6)') ( (ACOV(I)/ACOV(0)), I = 0,N-1)
Cc
C
C FFT approach (complex).
C    Compute FFT of data of length 2N.
C    Compute square of magnitude of transform components and place  
C       in complex array as d.p. parts.
C    Compute inverse transform, throwing away second half and
C       imaginary parts (which are zero), and multiply by length of 
C       sequence, 2N.  
c
      CALL DCFFTI (2*N,WSAVE)
      CALL DCFFTF (2*N,CEL,WSAVE)
c
C DCFFTF returns unscaled transforms. Actual transforms are output
C divided by (2N).
c
      DO I = 0,2*N-1
         CORR(I) = ABS(CEL(I) / (2*N)) **2 
      end do
c
C Since we compute transform times its conjugate, must divide by
C (2N) for each, i.e., (2N)**2.
      CALL DCFFTB (2*N,CORR,WSAVE)

      DO I = 0,N-1
         ACOV(I) = DBLE(CORR(I))*(2*N)
      end do
c
C Autocovariance is inverse transform times sequence length, 2N.
C Normally, all the scaling would be done  only once
C by dividing by 2N. We've broken it up for exposition.
c
      WRITE (*,*) 
     *  'EX 11.6: AUTOCORRELATION (COMPLEX FFT) OUTPUT SUPRESSED'
CCCCC      WRITE (*,'(5D14.6)') ( (ACOV(I)/ACOV(0) ), I = 0,N-1)
C
C FFT approach (d.p.).
C    Compute FFT of data of length 2N.
C    DEZFTF produces correctly scaled A's and B's so no extra scaling
C       needed to get transform.     
C    Compute array of square of each frequency component and place
C       in cosine array (A's) to be back transformed. Set B's to 0.
C    There are N A's, and N B's.
C    Note that care must be taken to compute magnitude correctly, 
C       0.5*(A(I)**2+B(I)**2) for I < N, twice that for I=N.
C    Compute back transform throwing away its second half.
C
      CALL DEZFTI (2*N,WSAVE)
      CALL DEZFTF (2*N,EL,AZERO,A,B,WSAVE)
      AZERO = AZERO*AZERO
C
      DO 150 I = 1,N
         IF(I.NE.N) THEN
            A(I) = (A(I)**2 + B(I)**2) / 2.0D0
         ELSE
            A(I) = (A(I)**2 + B(I)**2)
         ENDIF
         B(I) = 0.0D0
  150 CONTINUE
      CALL DEZFTB (2*N,ACOVR,AZERO,A,B,WSAVE)
      WRITE (*,*) 
     * 'EX 11.6: AUTOCORRELATION (DBLE FFT) OUTPUT REDUCED'
      WRITE (*,'(1X,3D18.10)') ( (ACOVR(I)/ACOVR(0) ), I = 0,19)
CCCCC      WRITE (*,'(5D14.6)') ( (ACOVR(I)/ACOVR(0) ), I = 0,N-1)
      WRITE (*,*)
      WRITE (*,*) 
     * 'REFERENCE RESULTS (EX 11.6 PARTIAL-LAST 7 LINES) FROM IBM PC/AT'
      WRITE (*,*)    
     *'  0.1000000000E+01  0.6067398062E+00  0.3538571074E+00'
      WRITE (*,*)
     *'  0.1843027185E+00 -0.1529114172E-01 -0.2198717698E+00'
      WRITE (*,*)
     *' -0.2962181307E+00 -0.2914580887E+00 -0.1550562043E+00'
      WRITE (*,*)
     *'  0.3680803130E-01  0.1735588664E+00  0.2665974779E+00'
      WRITE (*,*)
     *'  0.3049896546E+00  0.2011691615E+00  0.1780362880E-01'
      WRITE (*,*)
     *' -0.2103672514E+00 -0.3773868663E+00 -0.4379599788E+00'
      WRITE (*,*)
     *' -0.4603327597E+00 -0.4255892823E+00'

      return
 1000 WRITE (*,*) 'CANNOT FIND THE DATA FILE: elnino.dat '
      END
      subroutine test19 ( )

c*********************************************************************72
c
cc TEST19 minimizes a function using UNCMND.
c
      PARAMETER (N = 10, LWORK = N*(N+10))
      DOUBLE PRECISION X0(N), X(N), F, WORK(LWORK), S, W
      EXTERNAL test19_calcf

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST19'
      write ( *, '(a)' ) 
     &  '  UNCMND finds the minimum of a function.'
      write ( *, '(a)' ) ' '
C
C SPECIFY INITIAL ESTIMATE OF THE SOLUTION
C
      WRITE (*,*) 'COMPUTING...'
      DO I = 1,N
         X0(I) = I / FLOAT(N+1)
      end do
C
C MINIMIZE FUNCTION
C
      CALL UNCMND (N, X0, test19_calcf, X, F, IERROR, WORK, LWORK)
C
C PRINT RESULTS
C
      WRITE (*,*) 'UNCMND RESULTS'
      IF (IERROR .NE. 0) WRITE (*,*) ' ERROR CODE =', IERROR
      WRITE (*,'(1X,A,1X,D20.12)') ' F(X*) = ', F
      WRITE (*,*) ' X* ='
      WRITE (*,800) (X(I), I = 1,N)

      WRITE (*,*)
      WRITE (*,*) 'REFERENCE RESULTS FROM IBM PC/AT'
      WRITE (*,*) 'UNCMND WARNING -- INFO = 1: PROBABLY CONVERGED, GRADI
     *ENT SMALL'
      WRITE (*,*) 'UNCMND RESULTS'
      WRITE (*,*) ' ERROR CODE =           1'
      WRITE (*,*) ' F(X*) =    0.100000000000E+01'
      S=0.1D+1
      W=0.9999999D0
      WRITE (*,*) ' X* ='
      WRITE (*,799) S,S,S,S,S
      WRITE (*,799) S,S,S,S,W

      return
799   FORMAT (1X,D14.7,2X,D14.7,2X,D14.7,2X,D14.7,2X,D14.7)
800   FORMAT (1X,D14.7,2X,D14.7,2X,D14.7,2X,D14.7,2X,D14.7)
      END
      SUBROUTINE test19_calcf (N, X, F)

c*********************************************************************72
c
cc TEST19_CALCF evaluates the function to be minimized.
c
      DOUBLE PRECISION X(N), F, T1, T2

      T1 = 0.0
      T2 = 0.0
      DO I = 2,N
         T1 = T1 + (X(I)-X(I-1)**2)**2
         T2 = T2 + (1.0-X(I-1))**2
      end do
      F = 1.0 + 100.0*T1 + T2

      RETURN
      END
      subroutine test20 ( )

c*********************************************************************72
c
cc TEST20 tests DQAGI.
c
C          COMPUTE INTEGRAL OF EXP(-X)*COS(X*X)**2 ON [0,INFINITY)
C               RESULT IS 0.70260...
C
      INTEGER LIMIT,LENW
      PARAMETER(LIMIT=100,LENW=LIMIT*4)
      DOUBLE PRECISION BOUND,EPSABS,EPSREL,RESULT,ABSERR,WORK(LENW)
      INTEGER INF,NEVAL,IER,LAST,IWORK(LIMIT)
      EXTERNAL test20_f

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST20'
      write ( *, '(a)' ) 
     &  '  DQAGI estimates an integral.'
      write ( *, '(a)' ) ' '

      BOUND = 0.0D0
      INF = 1
      EPSABS = 0.0D0
      EPSREL = 1.D-5
      CALL DQAGI(test20_f,BOUND,INF,EPSABS,EPSREL,RESULT,ABSERR,NEVAL,
     *    IER,LIMIT,LENW,LAST,IWORK,WORK)
      WRITE(*,*)'DQAGI RESULT, ABSERR, NEVAL, IER: '
      WRITE(*,'(2X,2D20.10,I6,I5)') RESULT, ABSERR, NEVAL, IER

      WRITE(*,*)
      WRITE(*,*)'REFERENCE RESULTS FROM IBM PC/AT'
      WRITE(*,*)'   BE SURE THAT UNDERFLOWS ARE SET TO ZERO...'
      WRITE(*,*) 'DQAGI RESULT, ABSERR, NEVAL, IER: '
      WRITE(*,*)'     0.7026028726E+00    0.6401518223E-05  1005    0'

      return
      END
      DOUBLE PRECISION FUNCTION test20_f(X)

c*********************************************************************72
c
cc TEST20_F evaluates the integrand.
c
      double precision x,exp,cos
      test20_f = exp(-x)*cos(x*x)**2
      return
      end
      subroutine test21 ( )

c*********************************************************************72
c
cc TEST21 runs the REACTOR SHIELDING PROBLEM
C
      DOUBLE PRECISION MU,E,AZM,D,DIST2C,X,Y,Z,EA,ER,ET,SA,SR,ST,THICK
      INTEGER particle_exit,I,NA,NT,NR,NTOT,NPART
      LOGICAL ABSORB
      COMMON THICK

      ntot = 1000
      thick = 2.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST21'
      write ( *, '(a)' ) '  Reactor shielding problem'
      write ( *, '(a,i8)' ) '  Number of particles = ', ntot
      write ( *, '(a,g14.6)' ) '  Slab thickness = ', thick
c
C  INITIALIZE PARTICLE COUNTS AND ENERGY TALLIES
c
      CALL INIT(NA,NT,NR,NPART,EA,ET,ER,SA,SR,ST)
C
C   MAIN LOOP OVER NTOT PARTICLES.
c
    1 CONTINUE

      NPART=NPART+1
c
C FINISHED, COMPUTE AND PRINT AVERAGES, STANDARD DEVIATIONS
c
      IF(NPART.EQ.NTOT)THEN
        CALL OUTPUT(NA,EA,SA,NR,ER,SR,NT,ET,ST,NTOT)
        GOTO 3
      ENDIF
C
C  SOURCE GENERATES A NEW PARTICLE WITH POSITION, DIRECTION, ENERGY
c
      CALL SOURCE(E,MU,AZM,X,Y,Z)

    2 CONTINUE
C
C  COMPUTE-DISTANCE-TO-COLLISION
c
      D=DIST2C(E)
C
C  UPDATE POSITION TO DECIDE IF PARTICLE EXITS OR COLLIDES
c
      CALL UPDATE(MU,AZM,D,X,Y,Z)

      I=particle_exit(X,Y,Z)
c
C  RETURNS -1, 0 , 1 FOR 'OUT ON LEFT', 'IN', 'OUT ON RIGHT'
c
C  EXIT ON LEFT (REFLECTED), TALLY AND GOTO SOURCE
c
      IF(I.LT.0)THEN

           NR=NR+1
           ER=ER+E  
           SR=SR+E*E
           GOTO 1
c
C  EXIT ON RIGHT (TRANSMITTED), TALLY AND GOTO SOURCE
c
      ELSEIF(I.GT.0)THEN

           NT=NT+1
           ET=ET+E
           ST=ST+E*E
           GOTO 1
c
C  COLLISION & ABSORBED, TALLY AND GOTO SOURCE
c
      ELSEIF(I.EQ.0 .AND. ABSORB())THEN

           NA=NA+1
           EA=EA+E
           SA=SA+E*E
           GOTO 1
c
C  COLLISION & SCATTERED, FIND SCATTERING ANGLE AND ENERGY
c
      ELSE
        CALL SCATT(MU,AZM,E)
        GOTO 2
      ENDIF

3     continue

      RETURN
      END
      SUBROUTINE INIT(NA,NT,NR,NPART,EA,ET,ER,SA,SR,ST)

c*********************************************************************72
c
cc INIT initializes the data.
c
      DOUBLE PRECISION EA,ER,ET,SA,SR,ST
      INTEGER NA,NT,NR,NPART
      NA=0
      NT=0
      NR=0
      NPART=0
      EA=0.0D0
      ET=0.0D0
      ER=0.0D0
      SA=0.0D0
      SR=0.0D0
      ST=0.0D0
      RETURN
      END
      SUBROUTINE SOURCE(E,MU,AZM,X,Y,Z)

c*********************************************************************72
c
cc SOURCE models the emitter.
c
C          EMITTER  
C          RETURNS MU UNIFORM ON [0,1] AS THE COSINE OF AN ANGLE, 
C              AND AZM UNIFORM ON [0,2*PI] AS ASIMUTHAL ANGLE 
c
      DOUBLE PRECISION E,MU,AZM,X,Y,Z,PI,ENERGY,DUNI

      PI=4*ATAN(1.0D0)
      MU=DUNI()
      AZM=DUNI()*(2.0D0*PI)
      X=0.0D0
      Y=0.0D0
      Z=0.0D0
c
C  RETURNS ENERGY FROM SOURCE ENERGY DISTRIBUTION
c
      E=ENERGY()
      RETURN
      END
      SUBROUTINE UPDATE(MU,AZM,D,X,Y,Z)

c*********************************************************************72
c
cc UPDATE
c
      DOUBLE PRECISION MU,AZM,D,X,Y,Z,R     
 
      X=X+D*MU
      R=D*SQRT(1.0D0-MU*MU)
      Y=Y+R*COS(AZM)
      Z=Z+R*SIN(AZM)
      RETURN
      END
      SUBROUTINE SCATT(MU,AZM,E)

c*********************************************************************72
c
cc SCATT returns scattering angle and energy.
c
      DOUBLE PRECISION MU,AZM,E,PI,DUNI
C
C  ISOTROPIC SCATTERING, I.E., UNIFORM ON SPHERE 
C
      PI=4.0D0*ATAN(1.0D0)
      MU=-1.0D0+2.0D0*DUNI()
      AZM=DUNI()*2D0*PI
C
C  FIND SCATTERING ENERGY, UNIFORM IN [.3D0*E,E]
C
      E=(DUNI()*0.7D0+0.3D0)*E
      RETURN
      END
      LOGICAL FUNCTION ABSORB()

c*********************************************************************72
c
cc ABSORB returns true if particle is absorbed. 
c
c  OCCURS WITH PROB PA.
c
      DOUBLE PRECISION PA,DUNI
      PARAMETER(PA=0.1D0)

      IF(DUNI().LE.PA)THEN
        ABSORB=.TRUE.
      ELSE
        ABSORB=.FALSE.
      ENDIF

      RETURN
      END
      INTEGER FUNCTION particle_exit (X,Y,Z)

c*********************************************************************72
c
cc PARTICLE_EXIT reports the location of the particle.
c
c  RETURNS -1, 0, +1 AS PARTICLE IS   OUTSIDE ON LEFT, INSIDE, 
C                            OR OUTSIDE ON RIGHT.
c
      DOUBLE PRECISION X,Y,Z,THICK
      COMMON THICK

      IF(X.GT.THICK)THEN
        particle_exit=1
      ELSEIF(X.LT.0.0D0)THEN
        particle_exit=-1
      ELSE
        particle_exit=0
      ENDIF

      RETURN
      END
      DOUBLE PRECISION FUNCTION DIST2C(E)

c*********************************************************************72
c
cc DIST2C returns distance to collision.
c
c  It uses an exponential distribution with parameter  `cross section'
c
      DOUBLE PRECISION LOG,CROSS,E,DUNI

      DIST2C=-LOG(DUNI())/CROSS(E)

      RETURN
      END
      DOUBLE PRECISION FUNCTION ENERGY()

c*********************************************************************72
c
cc ENERGY returns energy, with distribution const/sqrt(energy) over [EMIN,EMAX].
c
C               USE INVERSE FUNCTION APPROACH TO COMPUTE THIS
C
C      ENERGY MIN, MAX IN MEV
c
      DOUBLE PRECISION C,SQRT,EMIN,EMAX,DUNI
      PARAMETER(EMIN=1.0D-3,EMAX=2.5D0)

      C=1.0D0/(2.0D0*(SQRT(EMAX)-SQRT(EMIN)))
      ENERGY=( DUNI()/(2*C)+SQRT(EMIN) )**2
      RETURN
      END
      DOUBLE PRECISION FUNCTION CROSS(E)

c*********************************************************************72
c
cc CROSS returns cross section (fictional) for energy in range [EMIN,EMAX]
c
      DOUBLE PRECISION E,S,ABS,SIN,Y,EXP

      S=ABS(SIN(100.0D0*(EXP(E)-1.0D0))+SIN(18.81D0*(EXP(E)-1.0D0)))
      Y=MAX(0.02D0,S)
      CROSS=10.0D0*EXP(-0.1D0/Y)
      RETURN
      END
      SUBROUTINE OUTPUT(NA,EA,SA,NR,ER,SR,NT,ET,ST,NTOT)

c*********************************************************************72
c
cc OUTPUT prints information.
c
      DOUBLE PRECISION EA,SA,ER,SR,ET,ST
      INTEGER NA,NR,NT,NTOT

      WRITE(*,*) 'TALLIES '

      IF(NA.GT.0)THEN
           EA=EA/NA
           SA=SQRT(SA/NA-EA*EA)
      ENDIF

      WRITE(*,'(1X,A,D10.4,2D20.10)') '% ABSORBED, ENERGY, SD: ',
     *       FLOAT(NA)/NTOT*100,EA,SA

      IF(NR.GT.0)THEN
           ER=ER/NR
           SR=SQRT(SR/NR-ER*ER)
      ENDIF

      WRITE(*,'(1X,A,D10.4,2D20.10)') '% REFLECTED, ENERGY, SD: ',
     *       FLOAT(NR)/NTOT*100,ER,SR

      IF(NT.GT.0)THEN
        ET=ET/NT
        ST=SQRT(ST/NT-ET*ET)
      ENDIF

      WRITE(*,'(1X,A,D10.4,2D20.10)') '% TRANSMITTED, ENERGY, SD: ',
     *       FLOAT(NT)/NTOT*100,ET,ST
      WRITE(*,*)

      RETURN
      END
