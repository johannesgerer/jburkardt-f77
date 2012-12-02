      PROGRAM IQTEST

c*********************************************************************72
c
C     TEST PROGRAM FOR IQPACK: DOUBLE PRECISION VERSION
C
      DOUBLE PRECISION A,ALFA,ALPHA,B,BETA,FRQ,QFSUM,QFSX,T,WF,WTS
      INTEGER I,IER,IK,IND,IWF,KEY,KIND,LI,LO,LU,MLT,MTK,NDX,NIWF,
     1        NT,NTMAX,NWF,NWTS
      DOUBLE PRECISION ZERO,HALF,ONE,TWO,THREE,FOUR,PI
      DIMENSION T(100),MLT(100),WTS(100),NDX(100),WF(10000),IWF(500)
C
      PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0,TWO=2.0D0)
      PARAMETER (THREE=3.0D0,FOUR=4.0D0)
      PARAMETER (PI=3.14159265358979323846264338327950 D0)
CS    PARAMETER (ZERO=0.0E0,HALF=0.5E0,ONE=1.0E0,TWO=2.0E0)
CS    PARAMETER (THREE=3.0E0,FOUR=4.0E0)
CS    PARAMETER (PI=3.14159 26535 89793 23846 26433 83279 50 E0)
C
C     FUNCTIONS AND SUBROUTINES REFERENCED - CEGQF,CEGQFS,CEIQF,
C         CEIQFS,CGQF,CIQF,CIQFS,CLIQF,CLIQFS,EIQFS,F,COS,
C         SIN
      EXTERNAL F
C     LOGICAL UNIT NUMBERS: LO,LU=OUTPUT  LI=INPUT (NOT USED)
      DATA LO/6/,LI/5/,LU/6/

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS655_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS655 library.'

C     SIZE OF KNOTS ARRAY AND WORKSPACE IN DIMENSION STATEMENTS
      NTMAX=100
      NWF=10000
      NIWF=500
C     GENERATE FEJER TYPE RULES - I.E. GAUSS CHEBYSHEV KNOTS AND
C     LEGENDRE WEIGHT FUNCTION
C     SET WEIGHT FUNCTION TYPE AND PARAMETERS
      KIND=1
      A=ONE
      B=FOUR
C     ALPHA, BETA NOT USED IN LEGENDRE WEIGHT FUNCTION BUT SET ANYWAY
      ALPHA=HALF
      BETA=-HALF/TWO
C
C     SET MULTIPLICITY OF KNOTS
      MTK=2
      KEY=1
C     NUMBER OF KNOTS.
      NT=5
C     IF NT IS CHANGED IT MUST BE .LE. NTMAX
C     SET SIZE OF WEIGHTS ARRAY
      NWTS=MTK*NT
      FRQ=PI/TWO/NT
C
C     TEST THE DRIVERS
      DO 110 IND=1,7
C
C        SET UP THE KNOTS FOR NT-POINT FEJER TYPE RULE
C        FOR IND=2,4,6 SCALE THE KNOTS TO GIVEN A,B
C        OTHERWISE - DEFAULT A,B - DON'T SCALE
C
         DO 10 I=1,NT
            T(I)=COS((TWO*I-1)*FRQ)
            IF(MOD(IND,2).EQ.0) THEN
              IF(KIND.EQ.5.OR.KIND.EQ.6.OR.KIND.EQ.8) THEN
                T(I)=T(I)+A
              ELSE
                T(I)=((B-A)*T(I)+(A+B))/TWO
              ENDIF
            ENDIF
10       CONTINUE
C        SET KNOT MULTIPLICITIES
         DO 20 I=1,NT
            MLT(I)=MTK
20       CONTINUE
         WRITE(LO,*)'----------------------------------------'
         GOTO(30,40,50,60,70,80,90),IND
C
30       CONTINUE
         WRITE(LO,*)' BEGINNING TEST OF CIQFS'
         CALL CIQFS(NT,T,MLT,NWTS,WTS,NDX,KEY,KIND,ALPHA,BETA,
     1              LU,NWF,WF,IER)
         WRITE(LO,*)' RETURN FROM CIQFS: ERROR INDICATOR IER= ',IER
         GOTO 110
C
40       CONTINUE
C        CIQF CALLS CIQFS, CGQF CALLS CGQFS
C        SO THIS TESTS BOTH
         WRITE(LO,*)' BEGINNING TEST OF CIQF,CIQFS,CGQF AND CGQFS'
         WRITE(LO,*)' WITH ALL CLASSICAL WEIGHT FUNCTIONS'
         A=-HALF
         B=TWO
         ALFA=-HALF
         DO 45 IK=1,8
            KIND=IK
            BETA=TWO
            IF(KIND.EQ.8)BETA=-FOUR*FOUR
            WRITE(LO,*)' KNOTS AND WEIGHTS OF GAUSS QF'
            WRITE(LO,*)' COMPUTED BY CGQF(S)'
            CALL CGQF(NT,T,WTS,KIND,ALPHA,BETA,A,B,LU,
     1                 NWF,WF,NIWF,IWF,IER)
            WRITE(LO,*)
            WRITE(LO,*)' RETURN FROM CGQF(S): ERROR INDICATOR IER= ',IER
C           NOW COMPUTE THE WEIGHTS FOR THE SAME KNOTS BY CIQF
            WRITE(LO,*)' WEIGHTS OF GAUSS QF COMPUTED FROM THE'
            WRITE(LO,*)' KNOTS BY CIQF(S)'
            CALL CIQF(NT,T,MLT,NWTS,WTS,NDX,KEY,KIND,ALPHA,BETA,A,B,
     1               LU,NWF,WF,IER)
            WRITE(LO,*)
            WRITE(LO,*)' RETURN FROM CIQF(S): ERROR INDICATOR IER= ',IER
 45      CONTINUE
         GOTO 110
C
50       CONTINUE
         KIND=1
         WRITE(LO,*)' BEGINNING TEST OF CEIQFS'
         QFSX=COS(-ONE)-COS(ONE)
C        QFSX IS EXACT VALUE OF INTEGRAL
         CALL CEIQFS(NT,T,MLT,KIND,ALPHA,BETA,F,QFSUM,
     1               NWF,WF,NIWF,IWF,IER)
         WRITE(LO,*)
         WRITE(LO,*)' RETURN FROM CEIQFS: ERROR INDICATOR IER= ',IER
         WRITE(LO,*)' INTERGRAL OF SIN(X) ON [-1,1] BY FEJER TYPE RULE'
         WRITE(LO,*)' WITH ',NT,' POINTS OF MULTIPLICITY ',MTK
         WRITE(LO,*)' QUADRATURE FORMULA:',QFSUM
         WRITE(LO,*)' EXACT VALUE       :',QFSX
         WRITE(LO,*)' ERROR             :',abs ( qfsum - qfsx )
         GOTO 110
C
60       CONTINUE
         WRITE(LO,*)' BEGINNING TEST OF CEIQF'
         KIND=1
         QFSX=COS(A)-COS(B)
C        QFSX IS EXACT VALUE OF INTEGRAL
         CALL CEIQF(NT,T,MLT,KIND,ALPHA,BETA,A,B,F,QFSUM,
     1              NWF,WF,NIWF,IWF,IER)
         WRITE(LO,*)
         WRITE(LO,*)' RETURN FROM CEIQF: ERROR INDICATOR IER= ',IER
         WRITE(LO,*)' INTERGRAL OF SIN(X) FROM ',A,' TO ',B
         WRITE(LO,*)' BY FEJER TYPE RULE WITH ',NT,' POINTS'
         WRITE(LO,*)' OF MULTIPLICITY ',MTK
         WRITE(LO,*)' QUADRATURE FORMULA:',QFSUM
         WRITE(LO,*)' EXACT VALUE       :',QFSX
         WRITE(LO,*)' ERROR             :',abs ( qfsum - qfsx )
         GOTO 110
C
70       CONTINUE
         WRITE(LO,*)' BEGINNING TEST OF CLIQFS'
         CALL CLIQFS(NT,T,WTS,KIND,ALPHA,BETA,LU,NWF,WF,NIWF,IWF,IER)
         WRITE(LO,*)
         WRITE(LO,*)'RETURN FROM CLIQFS. IER= ',IER
         GOTO 110
C
80       CONTINUE
         WRITE(LO,*)' BEGINNING TEST OF CLIQF AND EIQFS'
         CALL CLIQF(NT,T,WTS,KIND,ALPHA,BETA,A,B,LU,NWF,WF,NIWF,IWF,IER)
         WRITE(LO,*)
         WRITE(LO,*)' RETURN FROM CLIQF: ERROR INDICATOR IER= ',IER
         CALL EIQFS(NT,T,WTS,F,QFSUM,IER)
         WRITE(LO,*)' RETURN FROM EIQFS: ERROR INDICATOR IER= ',IER
         QFSX=COS(A)-COS(B)
         WRITE(LO,*)' INTERGRAL OF SIN(X) FROM ',A,' TO ',B
         WRITE(LO,*)' BY FEJER TYPE RULE WITH ',NT,' POINTS'
         WRITE(LO,*)' OF MULTIPLICITY ONE'
         WRITE(LO,*)' QUADRATURE FORMULA:',QFSUM
         WRITE(LO,*)' EXACT VALUE       :',QFSX
         WRITE(LO,*)' ERROR             :',abs ( qfsum - qfsx )
         GOTO 110
90       CONTINUE
         write ( *, '(a)' ) ' '
         write ( *, '(a)' ) 'TEST08'
         WRITE(LO,*)' BEGINNING TEST OF CEGQF'
         NT=12
C
C        EXPONENTIAL WEIGHT FUNCTION ON (A,B) CHOSEN
C
         KIND=7
         ALPHA=ONE
         BETA=ZERO
         A=ONE
         B=FOUR
C
C        COMPUTE AND EVALUATE THE GAUSS QF
         CALL CEGQF(NT,KIND,ALPHA,BETA,A,B,F,QFSUM,
     1              NWF,WF,NIWF,IWF,IER)
         WRITE(LO,*)' RETURN FROM CEGQF: ERROR INDICATOR IER= ',IER
         IF(IER.NE.0) GOTO 110
C
         QFSX=(B-A)*HALF*(COS(A)-COS(B))+SIN(B)+SIN(A)
     1         -TWO*SIN((A+B)/TWO)
         WRITE(LO,*)' INTERGRAL OF X*SIN(X) FROM ',A,' TO ',B
         WRITE(LO,*)' BY GAUSS-EXPONENTIAL RULE WITH ',NT,' POINTS'
         WRITE(LO,*)' QUADRATURE FORMULA:',QFSUM
         WRITE(LO,*)' EXACT VALUE       :',QFSX
         WRITE(LO,*)' ERROR             :',abs ( qfsum - qfsx )
C
100      CONTINUE
         WRITE(LO,*)' BEGINNING TEST OF CEGQFS'
         NT=12
C
C        EXPONENTIAL WEIGHT FUNCTION, DEFAULT VALUES OF A,B
C
         KIND=7
         ALPHA=ONE
         BETA=ZERO
C
C        COMPUTE AND EVALUATE THE GAUSS QF
         CALL CEGQFS(NT,KIND,ALPHA,BETA,F,QFSUM,
     1               NWF,WF,NIWF,IWF,IER)
         WRITE(LO,*)' RETURN FROM CEGQFS: ERROR INDICATOR IER= ',IER
         IF(IER.NE.0) GOTO 110
C
         A=-ONE
         B=ONE
         QFSX=(B-A)*HALF*(COS(A)-COS(B))+SIN(B)+SIN(A)
     1         -TWO*SIN((A+B)/TWO)
         WRITE(LO,*)' INTERGRAL OF X*SIN(X) FROM ',A,' TO ',B
         WRITE(LO,*)' BY GAUSS-EXPONENTIAL RULE WITH ',NT,' POINTS'
         WRITE(LO,*)' QUADRATURE FORMULA:',QFSUM
         WRITE(LO,*)' EXACT VALUE       :',QFSX
         WRITE(LO,*)' ERROR             :',abs ( qfsum - qfsx )
C
110   CONTINUE

      call test09 ( )
c
c  Compute 15 points of an example of each rule.
c
      do kind = 1, 8
        nt = 15
        if ( kind == 8 ) then
          alpha = 1.0D+00
          beta = - alpha - 2 * nt - 2
        else
          alpha = 0.0D+00
          beta = 0.0D+00
        end if
        call test10 ( nt, kind, alpha, beta )
      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS655_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      STOP
      END
      subroutine test09 ( )

c*********************************************************************72
c
cc TEST09 calls CGQFS to compute and print generalized Gauss-Hermite rules.
c
c  Modified:
c
c    09 January 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nt
      parameter ( nt = 15 )
      integer niwf
      parameter ( niwf = 2 * nt )
      integer nwf
      parameter ( nwf = 3 * nt+ 14 )

      double precision alpha
      double precision beta
      integer ier
      integer io
      integer iwf(niwf)
      integer kind
      double precision t(nt)
      double precision wf(nwf)
      double precision wts(nt)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST09'
      write ( *, '(a)' ) '  Call CGQFS for generalized Hermite rules.'

      kind = 6
      alpha = 1.0D+00
      beta = 0.0D+00
      io = -6

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  NT = ', nt
      write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha

      call cgqfs ( nt, t, wts, kind, alpha, beta, io, nwf, wf, niwf,
     &  iwf, ier )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i2)' ) '  IER = ', ier

      return
      end
      subroutine test10 ( nt, kind, alpha, beta )

c*********************************************************************72
c
cc TEST10 calls CDGQF to compute a quadrature formula.
c
c  Modified:
c
c    10 January 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nt_max
      parameter ( nt_max = 15 )
      integer nwf
      parameter ( nwf = 2 * nt_max )

      double precision alpha
      double precision beta
      integer i
      integer ier
      integer kind
      integer nt
      double precision t(nt_max)
      double precision wf(nwf)
      double precision wts(nt)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST10'
      write ( *, '(a)' ) '  Call CDGQF to compute a quadrature formula.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  KIND = ', kind
      write ( *, '(a,g14.6)' ) '  ALPHA = ', alpha
      write ( *, '(a,g14.6)' ) '  BETA  = ', beta

      call cdgqf ( nt, t, wts, kind, alpha, beta, nwf, wf, ier )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) ' Index     Abscissas                 Weights'
      write ( *, '(a)' ) ' '
      do i = 1, nt
        write ( *, '(2x,i4,2x,g24.16,2x,g24.16)' ) i, t(i), wts(i)
      end do

      return
      end
      DOUBLE PRECISION FUNCTION F(X,I)

c*********************************************************************72
c
C     THIS FUNCTION GENERATES VALUES OF THE INTEGRAND
C     AND POSSIBLY ITS DERIVATIVES. SINCE IT APPEARS AS AN ACTUAL
C     PARAMETER IN THE ROUTINE WHICH USES IT, IT SHOULD BE DECLARED
C     IN AN EXTERNAL STATEMENT IN ANY CALLING PROGRAM
      DOUBLE PRECISION X
      INTEGER I,K,L
C     FUNCTIONS AND SUBROUTINES REFERENCED - COS,SIN,MOD
C
C     RETURNS THE I-TH DERIVATIVE OF SIN(T) AT T=X
C     (ZERO-TH DERIVATIVE IS THE FUNCTION VALUE)
      K=MOD(I,4)
      L=MOD(K,2)
      IF(L.EQ.0)F=SIN(X)
      IF(L.EQ.1)F=COS(X)
      IF(K.GT.1)F=-F
      RETURN
      END
