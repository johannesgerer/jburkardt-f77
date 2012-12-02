      program main

c*********************************************************************72
c
cc MAIN is the main program for TOMS793_PRB.
c
c  Modified:
c
c    08 November 2012
c
c  Author:
c
c    Walter Gautschi
c
c  Reference:
c
c    Walter Gautschi,
c    Algorithm 793: GQRAT, Gauss Quadrature for Rational Functions,
c    ACM Transactions on Mathematical Software,
c    Volume 25, Number 2, pages 213-239, June 1999.
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS793_PRB:'
      write ( *, '(a)' ) '  Test the TOMS793 library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test05 ( )
      call test06 ( )
      call test07 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS793_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 carries out the first test.
c
c  Discussion:
c
c    I1(w) = Integral ( -1 <= t <= +1 ) (pi*t/w) / sin(pi*t/w) dt,
c    with 1 < w.
c
c    The integrand has simple real poles at the integer multiples of w.
c
c    We take the m poles closest to, and on either side, of [-1,+1],
c    to determine the space Qm of rational functions.
c
c    In this example, W is taken to be 1.1.
c
c  Modified:
c
c    08 November 2012
c
c  Author:
c
c    Walter Gautschi
c
c  Reference:
c
c    Walter Gautschi,
c    Algorithm 793: GQRAT, Gauss Quadrature for Rational Functions,
c    ACM Transactions on Mathematical Software,
c    Volume 25, Number 2, pages 213-239, June 1999.
c
      implicit none

      INTEGER NMAX,NCAPM
      PARAMETER (NMAX=10,NCAPM=100)
      DOUBLE PRECISION OM
      PARAMETER (OM=1.1D0)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CONST,EPS,EXACT,SGN,SUM
      REAL ERR
      INTEGER IERAB,IERGQ,IERR,IROUT,K,KOUNT,M,MU,N,NCAP
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION A(NCAPM),ALPHA(NMAX+1),B(NCAPM),BE(NMAX+1),
     +                 BETA(NMAX+1),E(NCAPM),P0(NCAPM),P1(NCAPM),
     +                 P2(NCAPM),W(NCAPM),WG(NMAX),X(NCAPM),XII(2*NMAX),
     +                 XIR(2*NMAX),ZG(NMAX)
      INTEGER IS(2*NMAX)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION D1MACH,F1
      EXTERNAL D1MACH,F1
C     ..
C     .. External Subroutines ..
      EXTERNAL DABMOD,DGQRAT,DRECUR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DBLE,REAL
      integer i1mach, nout

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Test ???'

      nout = i1mach(2)

      IROUT = 1
c      irout=0
      EPS = D1MACH(3)*1.D2
      CALL DRECUR(NCAPM,1,0.D0,0.D0,A,B,IERR)
      WRITE(NOUT,FMT=9000)
c
c  Using M = 0 produces ordinary Gauss-Legendre quadrature.
c
      DO 30 N = 1,NMAX
          M = 2*N
c        m=n
c        m=2
c        m=0
          SGN = 1.D0
          DO MU = 1,M
              SGN = -SGN
              XIR(MU) = SGN/ (OM*DBLE((MU+1)/2))
              XII(MU) = 0.D0
              IS(MU) = 1
          end do
          CALL DABMOD(N+1,NCAPM,M,EPS,IROUT,A,B,XIR,XII,IS,ALPHA,BETA,
     &                NCAP,KOUNT,IERAB,BE,X,W,E,P0,P1,P2)
          IF (IERAB.NE.0) WRITE(NOUT,FMT=9010) IERAB
          CALL DGQRAT(N,M,ALPHA,BETA,XIR,XII,IS,ZG,WG,CONST,IERGQ,E)
          IF (IERGQ.NE.0) WRITE(NOUT,FMT=9020) IERGQ
          SUM = 0.D0
          DO K = 1,N
              SUM = SUM + WG(K)*F1(ZG(K),OM)
          end do
          EXACT = 4.467773646387766D0
          ERR = ABS(REAL((SUM-EXACT)/EXACT))
          WRITE(NOUT,FMT=9030) N,M,SUM,ERR,CONST
   30 CONTINUE

      return
 9000 FORMAT (5X,'n',4X,'m',8X,'integral',9X,'rel error',3X,'err const',
     +       /)
 9010 FORMAT (1X,'ierab=',I3)
 9020 FORMAT (1X,'iergq=',I3)
 9030 FORMAT (1X,2I5,D23.14,E12.4,D12.4)
      END
      FUNCTION F1(T,OM)

c*********************************************************************72
c
cc F1 is the integrand for TEST01.
c
c  Modified:
c
c    08 November 2012
c
c  Author:
c
c    Walter Gautschi
c
c  Reference:
c
c    Walter Gautschi,
c    Algorithm 793: GQRAT, Gauss Quadrature for Rational Functions,
c    ACM Transactions on Mathematical Software,
c    Volume 25, Number 2, pages 213-239, June 1999.
c
      implicit none

      DOUBLE PRECISION F1
      double precision OM,T
      DOUBLE PRECISION PI
      INTRINSIC ATAN,SIN

      PI = 4.D0*ATAN(1.D0)

      IF (T.NE.0.D0) then
        F1 = (PI*T/OM)/SIN(PI*T/OM)
      else
        F1 = 1.0D+00
      end if

      RETURN
      END
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 carries out test 2.
c
c  Discussion:
c
c    I2(w) = Integral ( -1 <= t <= +1 ) (pi*t/w) / sin(pi*t/w) dt,
c    with 1 < w.
c
c    The integrand has simple real poles at the integer multiples of w.
c
c    We take the m poles closest to, and on either side, of [-1,+1],
c    to determine the space Qm of rational functions.
c
c    In this example, W is taken to be very close to 1.
c
c  Modified:
c
c    08 November 2012
c
c  Author:
c
c    Walter Gautschi
c
c  Reference:
c
c    Walter Gautschi,
c    Algorithm 793: GQRAT, Gauss Quadrature for Rational Functions,
c    ACM Transactions on Mathematical Software,
c    Volume 25, Number 2, pages 213-239, June 1999.
c
      implicit none

      INTEGER NMAX,NCAPM
      PARAMETER (NMAX=10,NCAPM=100)
      DOUBLE PRECISION OM
      PARAMETER (OM=1.001D0)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CONST,EPS,EXACT,HN,HP,P,SGN,SUM,SUMN,SUMP,X1
      REAL ERR
      INTEGER IERAB,IERGC,IERGQ,IERR,IROUT,K,KOUNT,M,MU,N,NCAP,NU
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION A(NCAPM),A1(NMAX+1),ALPHA(NMAX+1),B(NCAPM),
     +                 B1(NMAX+1),BE(NMAX+1),BETA(NMAX+1),E(NCAPM),
     +                 P0(NCAPM),P1(NCAPM),P2(NCAPM),W(NCAPM),WG(NMAX),
     +                 X(NCAPM),XII(2*NMAX),XIR(2*NMAX),ZG(NMAX)
      INTEGER IS(2*NMAX)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION D1MACH,F2
      EXTERNAL D1MACH,F2
C     ..
C     .. External Subroutines ..
      EXTERNAL DABMOD,DGCHRS,DGQRAT,DRECUR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DBLE,LOG,REAL
      integer i1mach, nout

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Test ???'

      nout = i1mach(2)
      IROUT = 1
c      irout=0
      EPS = D1MACH(3)*1.D2
      CALL DRECUR(NCAPM,1,0.D0,0.D0,A,B,IERR)
      WRITE(NOUT,FMT=9000)

      DO N = 1,NMAX
          M = 2*N
c        m=n
c        m=4
c        m=2
          SGN = 1.D0
          DO NU = 1,M
              SGN = -SGN
              IF (NU.LE.M-2) XIR(NU) = SGN/DBLE((NU+3)/2)
              XII(NU) = 0.D0
              IS(NU) = 1
          end do
          CALL DABMOD(N+1,NCAPM,M-2,EPS,IROUT,A,B,XIR,XII,IS,A1,B1,NCAP,
     &                KOUNT,IERAB,BE,X,W,E,P0,P1,P2)
          IF (IERAB.NE.0) WRITE(NOUT,FMT=9010) IERAB
          XIR(M-1) = 1.D0/OM
          XIR(M) = -1.D0/OM
          X1 = OM
          HP = LOG(ABS((X1+1.D0)/ (X1-1.D0)))
          HN = -HP
          IF (M.GT.2) THEN
              SUMP = 0.D0
              SUMN = 0.D0
              DO NU = 1,M - 2
                  P = 1.D0
                  DO MU = 1,M - 2
                      IF (MU.EQ.NU) GO TO 20
                      P = (1.D0-XIR(MU)/XIR(NU))*P
                  end do
                  SUMP = SUMP + LOG(ABS((X1+1.D0)* (XIR(NU)+1.D0)/ ((X1-
     &                   1.D0)* (XIR(NU)-1.D0))))/ ((1.D0+XIR(NU)*X1)*P)
                  SUMN = SUMN + LOG(ABS((X1-1.D0)* (XIR(NU)+1.D0)/ ((X1+
     &                   1.D0)* (XIR(NU)-1.D0))))/ ((1.D0-XIR(NU)*X1)*P)
              end do
              HP = SUMP
              HN = SUMN
          END IF

          CALL DGCHRS(N+1,2,A1,B1,X1,HP,HN,ALPHA,BETA,IERGC)
          BETA(1) = - (OM**2)*BETA(1)
          CALL DGQRAT(N,M,ALPHA,BETA,XIR,XII,IS,ZG,WG,CONST,IERGQ,E)
          IF (IERGQ.NE.0) WRITE(NOUT,FMT=9020) IERGQ
          SUM = 0.D0
          DO K = 1,N
              SUM = SUM + WG(K)*F2(ZG(K),OM)
          end do
          EXACT = 12.929256850003D0
          ERR = ABS(REAL((SUM-EXACT)/EXACT))
          WRITE(NOUT,FMT=9030) N,M,SUM,ERR,CONST

      end do

      return

 9000 FORMAT (5X,'n',4X,'m',8X,'integral',10X,'rel error',4X,
     +       'err const',/)
 9010 FORMAT (1X,'ierab=',I3)
 9020 FORMAT (1X,'iergq=',I3)
 9030 FORMAT (1X,2I5,D23.14,E13.4,D13.4)
      END
      FUNCTION F2(T,OM)

c*********************************************************************72
c
cc F2 is the integrand for TEST02.
c
c  Modified:
c
c    08 November 2012
c
c  Author:
c
c    Walter Gautschi
c
c  Reference:
c
c    Walter Gautschi,
c    Algorithm 793: GQRAT, Gauss Quadrature for Rational Functions,
c    ACM Transactions on Mathematical Software,
c    Volume 25, Number 2, pages 213-239, June 1999.
c
      implicit none

      DOUBLE PRECISION f2, OM,T
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION PI
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ATAN,SIN
C     ..
      PI = 4.D0*ATAN(1.D0)

      IF (T.NE.0.D0) then
        F2 = (PI*T/OM)/SIN(PI*T/OM)
      else
        F2 = 1.0D+00
      end if

      RETURN
      END
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 carries out test 3.
c
c  Discussion:
c
c    I3(w) = Integral ( 0 <= t <= 1 ) sqrt(t) * gamma(1+t) / (t+w) dt
c    with 0 < w.
c
c    The integrand has poles at -w and the negative natural numbers.
c
c  Modified:
c
c    08 November 2012
c
c  Author:
c
c    Walter Gautschi
c
c  Reference:
c
c    Walter Gautschi,
c    Algorithm 793: GQRAT, Gauss Quadrature for Rational Functions,
c    ACM Transactions on Mathematical Software,
c    Volume 25, Number 2, pages 213-239, June 1999.
c
      implicit none

      INTEGER NMAX,NCAPM
      PARAMETER (NMAX=10,NCAPM=100)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CONST,EPS,OM,SUM
      REAL ERR
      INTEGER IERAB,IERGQ,IERR,IOM,IROUT,K,KOUNT,M,MU,N,NCAP
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION A(NCAPM),ALPHA(NMAX+1),B(NCAPM),BE(NMAX+1),
     +                 BETA(NMAX+1),E(NCAPM),EXACT(9),OOM(9),P0(NCAPM),
     +                 P1(NCAPM),P2(NCAPM),W(NCAPM),WG(NMAX),X(NCAPM),
     +                 XII(2*NMAX),XIR(2*NMAX),ZG(NMAX)
      INTEGER IS(2*NMAX)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION D1MACH,F3
      EXTERNAL D1MACH,F3
C     ..
C     .. External Subroutines ..
      EXTERNAL DABMOD,DGQRAT,DRECUR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,DBLE,REAL,SQRT
      integer i1mach, nout
C     ..
C     .. Data statements ..
      DATA OOM/.1D0,.5D0,1.D0,1.999D0,2.D0,3.D0,3.001D0,10.5D0,50.D0/
      DATA EXACT/7.6724625863706D0,2.5531371574419D0,1.4784672063106D0,
     +     .81771417926000D0,.81735216597811D0,.56726444593607D0,
     +     .56709144671961D0,.17313056604184D0,.037222082318054D0/

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  Test ???'

      nout = i1mach(2)
      WRITE(NOUT,FMT=9000)
      IROUT = 1
c      irout=0
      EPS = D1MACH(3)*1.D2
      CALL DRECUR(NCAPM,6,0.D0,-.5D0,A,B,IERR)
      DO K = 1,NCAPM
          A(K) = .5D0* (1.D0+A(K))
          B(K) = .25D0*B(K)
      end do
      B(1) = SQRT(8.D0)*B(1)
      WRITE(NOUT,FMT=9010)
      DO 50 IOM = 1,9
          OM = OOM(IOM)
          DO 40 N = 1,NMAX
              M = 2*N
c          m=n
c          m=2
c          m=0
              IF (M.GT.0) THEN
                  DO 20 MU = 1,M
                      IF (MU.EQ.1) THEN
                          XIR(MU) = 1.D0/OM

                      ELSE
                          XIR(MU) = 1.D0/DBLE(MU-1)
                      END IF

                      XII(MU) = 0.D0
                      IS(MU) = 1
   20             CONTINUE
              END IF

              CALL DABMOD(N+1,NCAPM,M,EPS,IROUT,A,B,XIR,XII,IS,ALPHA,
     +                    BETA,NCAP,KOUNT,IERAB,BE,X,W,E,P0,P1,P2)
              IF (IERAB.NE.0) WRITE(NOUT,FMT=9020) IERAB
              CALL DGQRAT(N,M,ALPHA,BETA,XIR,XII,IS,ZG,WG,CONST,IERGQ,E)
              IF (IERGQ.NE.0) WRITE(NOUT,FMT=9030) IERGQ
              SUM = 0.D0
              DO K = 1,N
                  SUM = SUM + WG(K)*F3(ZG(K),OM)
              end do
              ERR = ABS(REAL((SUM-EXACT(IOM))/EXACT(IOM)))
              IF (N.EQ.1) THEN
                  WRITE(NOUT,FMT=9040) N,M,SUM,ERR,CONST,REAL(OM)

              ELSE
                  WRITE(NOUT,FMT=9050) N,M,SUM,ERR,CONST
              END IF

   40     CONTINUE
          WRITE(NOUT,FMT=9000)
   50 CONTINUE
      return

 9000 FORMAT (/)
 9010 FORMAT (5X,'n',4X,'m',8X,'integral',10X,'rel error',4X,
     +       'err const',/)
 9020 FORMAT (1X,'ierab=',I3)
 9030 FORMAT (1X,'iergq=',I3)
 9040 FORMAT (1X,2I5,D23.14,E13.4,D13.4,'  om=',F6.3)
 9050 FORMAT (1X,2I5,D23.14,E13.4,D13.4)
      END
      FUNCTION F3(T,OM)

c*********************************************************************72
c
cc F3 is the integrand for TEST03.
c
c  Modified:
c
c    08 November 2012
c
c  Author:
c
c    Walter Gautschi
c
c  Reference:
c
c    Walter Gautschi,
c    Algorithm 793: GQRAT, Gauss Quadrature for Rational Functions,
c    ACM Transactions on Mathematical Software,
c    Volume 25, Number 2, pages 213-239, June 1999.
c
      implicit none

      DOUBLE PRECISION F3, OM,T
C     ..
C     .. Local Scalars ..
      INTEGER IERR
C     ..
C     .. External Functions ..
      DOUBLE PRECISION DGAMMA
      EXTERNAL DGAMMA

      F3 = DGAMMA ( 1.D0 + T, IERR )/ ( T + OM )

      RETURN
      END
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 carries out test 5.
c
c  Discussion:
c
c    I5(k,eta,theta) = Integral ( 0 <= t < +oo ) t^k sqrt (1+0.5*theta*t)
c      / ( exp ( -eta + t ) + 1 ) dt
c    with eta real, and 0 <= theta.
c
c    This is the generalized Fermi-Dirac integral.
c
c    The poles are all simple and occur in complex conjugate pairs
c    on the line Imag(z)=eta at odd integer multiples of pi away
c    from the real axis.
c
c  Modified:
c
c    08 November 2012
c
c  Author:
c
c    Walter Gautschi
c
c  Reference:
c
c    Walter Gautschi,
c    Algorithm 793: GQRAT, Gauss Quadrature for Rational Functions,
c    ACM Transactions on Mathematical Software,
c    Volume 25, Number 2, pages 213-239, June 1999.
c
      implicit none
c
c      parameter(eta=1.d0)
c      parameter(ak=1.5d0)
c      parameter(ak=2.5d0)
C     .. Parameters ..
      INTEGER NMAX,NCAPM
      PARAMETER (NMAX=10,NCAPM=200)
      DOUBLE PRECISION THETA
      PARAMETER (THETA=1.D-4)
      DOUBLE PRECISION ETA
      PARAMETER (ETA=-1.D0)
      DOUBLE PRECISION AK
      PARAMETER (AK=.5D0)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CONST,DEN,EPS,EXACT,PI,SUM
      REAL ERR
      INTEGER IERAB,IERGQ,IERR,IROUT,K,KOUNT,M,MU,N,NCAP
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION A(NCAPM),ALPHA(NMAX+1),B(NCAPM),BE(NMAX+1),
     +                 BETA(NMAX+1),E(NCAPM),P0(NCAPM),P1(NCAPM),
     +                 P2(NCAPM),W(NCAPM),WG(NMAX),X(NCAPM),XII(2*NMAX),
     +                 XIR(2*NMAX),ZG(NMAX)
      INTEGER IS(2*NMAX)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION D1MACH,F5
      EXTERNAL D1MACH,F5
C     ..
C     .. External Subroutines ..
      EXTERNAL DABMOD,DGQRAT,DRECUR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,DBLE,REAL
      integer i1mach, nout

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  Test ???'

      nout = i1mach(2)
C     ..
      IROUT = 1
c      irout=0
      EPS = D1MACH(3)*1.D2
      CALL DRECUR(NCAPM,7,AK,0.D0,A,B,IERR)
      PI = 4.D0*ATAN(1.D0)
      WRITE(NOUT,FMT=9000)
      DO 30 N = 1,NMAX
          M = 2*N
c        m=2*((n+1)/2)
c        m=2
c        m=0
          IF (M.GT.0) THEN
              DO MU = 1,M - 1,2
                  DEN = ETA**2 + DBLE(MU**2)* (PI**2)
                  XIR(MU) = -ETA/DEN
                  XIR(MU+1) = XIR(MU)
                  XII(MU) = DBLE(MU)*PI/DEN
                  XII(MU+1) = -XII(MU)
                  IS(MU) = 1
                  IS(MU+1) = 1
              end do
          END IF

          CALL DABMOD(N+1,NCAPM,M,EPS,IROUT,A,B,XIR,XII,IS,ALPHA,BETA,
     +                NCAP,KOUNT,IERAB,BE,X,W,E,P0,P1,P2)
          IF (IERAB.NE.0) WRITE(NOUT,FMT=9010) IERAB
          CALL DGQRAT(N,M,ALPHA,BETA,XIR,XII,IS,ZG,WG,CONST,IERGQ,E)
          IF (IERGQ.NE.0) WRITE(NOUT,FMT=9020) IERGQ
          SUM = 0.D0
          DO K = 1,N
              SUM = SUM + WG(K)*F5(ZG(K),ETA,THETA)
          end do
          EXACT = .29051241701949D0
c        exact=.46087845417799d0
c        exact=1.18607350107576d0
c        exact=1.39644182034912d0
c        exact=2.66187327910715d0
c        exact=7.62725609565345d0
          ERR = ABS(REAL((SUM-EXACT)/EXACT))
          WRITE(NOUT,FMT=9030) N,M,SUM,ERR,CONST
   30 CONTINUE
      return

 9000 FORMAT (5X,'n',4X,'m',8X,'integral',10X,'rel error',4X,
     +       'err const',/)
 9010 FORMAT (1X,'ierab=',I3)
 9020 FORMAT (1X,'iergq=',I3)
 9030 FORMAT (1X,2I5,D23.14,E13.4,D13.4)
      END
      FUNCTION F5(T,ETA,THETA)

c*********************************************************************72
c
cc F5 is the integrand for TEST01.
c
c  Modified:
c
c    08 November 2012
c
c  Author:
c
c    Walter Gautschi
c
c  Reference:
c
c    Walter Gautschi,
c    Algorithm 793: GQRAT, Gauss Quadrature for Rational Functions,
c    ACM Transactions on Mathematical Software,
c    Volume 25, Number 2, pages 213-239, June 1999.
c
      implicit none

      DOUBLE PRECISION ETA,F5,T,THETA
      INTRINSIC EXP,SQRT

      F5 = SQRT(1.D0+.5D0*THETA*T)/ (EXP(-ETA)+EXP(-T))

      RETURN
      END
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 carries out test 6.
c
c  Discussion:
c
c    I6(eta,theta) = Integral ( 0 <= t < +oo ) t sqrt (1+0.5*theta*t)
c      / ( exp(-eta) - exp(-t) ) * t^(k-1) * exp(-t) dt
c    where eta < 0, 0 <= theta.
c
c    The poles are complex conjugate, along Imag(z)=eta, at a distance
c    of multiples of 2 pi from the real axis.
c
c  Modified:
c
c    08 November 2012
c
c  Author:
c
c    Walter Gautschi
c
c  Reference:
c
c    Walter Gautschi,
c    Algorithm 793: GQRAT, Gauss Quadrature for Rational Functions,
c    ACM Transactions on Mathematical Software,
c    Volume 25, Number 2, pages 213-239, June 1999.
c
      implicit none
c
c      parameter(ak=1.5d0)
c      parameter(ak=2.5d0)
C     .. Parameters ..
      INTEGER NMAX,NCAPM
      PARAMETER (NMAX=10,NCAPM=200)
      DOUBLE PRECISION THETA
      PARAMETER (THETA=1.D-4)
      DOUBLE PRECISION ETA
      PARAMETER (ETA=-1.D0)
      DOUBLE PRECISION AK
      PARAMETER (AK=.5D0)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION CONST,DEN,EPS,EXACT,PI,SUM
      REAL ERR
      INTEGER IERAB,IERGQ,IERR,IROUT,K,KOUNT,M,MU,N,NCAP
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION A(NCAPM),ALPHA(NMAX+1),B(NCAPM),BE(NMAX+1),
     +                 BETA(NMAX+1),E(NCAPM),P0(NCAPM),P1(NCAPM),
     +                 P2(NCAPM),W(NCAPM),WG(NMAX),X(NCAPM),XII(2*NMAX),
     +                 XIR(2*NMAX),ZG(NMAX)
      INTEGER IS(2*NMAX)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION D1MACH,F6
      EXTERNAL D1MACH,F6
C     ..
C     .. External Subroutines ..
      EXTERNAL DABMOD,DGQRAT,DRECUR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,DBLE,REAL
      integer i1mach, nout

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) '  Test ???'

      nout = i1mach(2)
C     ..
      IROUT = 1
c      irout=0
      EPS = D1MACH(3)*1.D2
      CALL DRECUR(NCAPM,7,AK-1.D0,0.D0,A,B,IERR)
      PI = 4.D0*ATAN(1.D0)
      WRITE(NOUT,FMT=9000)
      DO 30 N = 1,NMAX
          M = 2*N - 1
c        m=2*((n+1)/2)-1
c        m=1
c        m=0
          IF (M.GE.3) THEN
              DO MU = 1,M - 2,2
                  DEN = ETA**2 + DBLE((MU+1)**2)* (PI**2)
                  XIR(MU) = -ETA/DEN
                  XIR(MU+1) = XIR(MU)
                  XII(MU) = DBLE(MU+1)*PI/DEN
                  XII(MU+1) = -XII(MU)
                  IS(MU) = 1
                  IS(MU+1) = 1
              end do
          END IF

          IF (M.GT.0) THEN
              XIR(M) = -1.D0/ETA
              XII(M) = 0.D0
              IS(M) = 1
          END IF

          CALL DABMOD(N+1,NCAPM,M,EPS,IROUT,A,B,XIR,XII,IS,ALPHA,BETA,
     &                NCAP,KOUNT,IERAB,BE,X,W,E,P0,P1,P2)
          IF (IERAB.NE.0) WRITE(NOUT,FMT=9010) IERAB
          CALL DGQRAT(N,M,ALPHA,BETA,XIR,XII,IS,ZG,WG,CONST,IERGQ,E)
          IF (IERGQ.NE.0) WRITE(NOUT,FMT=9020) IERGQ
          SUM = 0.D0
          DO K = 1,N
            SUM = SUM + WG(K)*F6(ZG(K),ETA,THETA)
          end do
          EXACT = .379708865998074D0
c        exact=.52608888707965d0
c        exact=1.266569126543118d0
          ERR = ABS(REAL((SUM-EXACT)/EXACT))
          WRITE(NOUT,FMT=9030) N,M,SUM,ERR,CONST
   30 CONTINUE
      return

 9000 FORMAT (5X,'n',4X,'m',8X,'integral',10X,'rel error',4X,
     +       'err const',/)
 9010 FORMAT (1X,'ierab=',I3)
 9020 FORMAT (1X,'iergq=',I3)
 9030 FORMAT (1X,2I5,D23.14,E13.4,D13.4)
      END
      FUNCTION F6(T,ETA,THETA)

c*********************************************************************72
c
cc F6 is the integrand for TEST06.
c
c  Modified:
c
c    08 November 2012
c
c  Author:
c
c    Walter Gautschi
c
c  Reference:
c
c    Walter Gautschi,
c    Algorithm 793: GQRAT, Gauss Quadrature for Rational Functions,
c    ACM Transactions on Mathematical Software,
c    Volume 25, Number 2, pages 213-239, June 1999.
c
      implicit none

      DOUBLE PRECISION ETA,F6,T,THETA
      DOUBLE PRECISION S,S1,TERM,X
      INTEGER L
      INTRINSIC ABS,DBLE,EXP,SQRT
      integer i1mach, nout

      X = T - ETA
      IF (ABS(X).LE.1.D0) THEN
          L = 0
          TERM = X
          S1 = TERM
   10     continue
          S = S1
          L = L + 1
          IF (L.GT.200) GO TO 20
          TERM = X*TERM/DBLE(L+1)
          S1 = S + TERM
          IF (S1.NE.S) GO TO 10
          F6 = T*SQRT(1.D0+.5D0*THETA*T)/ (EXP(-T)*S)
          RETURN

      ELSE
          F6 = T*SQRT(1.D0+.5D0*THETA*T)/ (EXP(-ETA)-EXP(-T))
          RETURN

      END IF

   20 nout = i1mach(2)
      WRITE(NOUT,FMT=9000)
      RETURN

 9000 FORMAT (1X,'exp series does not converge')
      END
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07 carries out test 7.
c
c  Discussion:
c
c    I7(eta,theta) = Integral ( 0 <= t < +oo ) t sqrt (1+0.5*theta*t)
c      / ( exp(-eta) - exp(-t) ) * t^(k-1) * exp(-t) dt
c    where eta < 0, 0 <= theta.
c
c    The poles are complex conjugate, along Imag(z)=eta, at a distance
c    of multiples of 2 pi from the real axis.
c
c    Here, we consider small values of |eta|.
c
c  Modified:
c
c    08 November 2012
c
c  Author:
c
c    Walter Gautschi
c
c  Reference:
c
c    Walter Gautschi,
c    Algorithm 793: GQRAT, Gauss Quadrature for Rational Functions,
c    ACM Transactions on Mathematical Software,
c    Volume 25, Number 2, pages 213-239, June 1999.
c
      implicit none
c
c Be sure that ncapm is greater or equalt to numax
c
c      parameter(ak=1.5d0)
c      parameter(ak=2.5d0)
C     .. Parameters ..
      INTEGER NMAX,NCAPM,NUMAX
      PARAMETER (NMAX=10,NCAPM=100,NUMAX=50)
      DOUBLE PRECISION THETA
      PARAMETER (THETA=1.D-4)
      DOUBLE PRECISION AK
      PARAMETER (AK=.5D0)
C     ..
C     .. Local Scalars ..
      DOUBLE PRECISION C1,C2,CONST,DX,DY,EPS,ETA,EXACT,FMU,FMU2,FNU,
     +                 FNU2,HN,HP,P,PI,SQPI,SUM,SUM0,T,X1,XA
      REAL ERR
      INTEGER IERAB,IERGC,IERGQ,IERK,IERR,IROUT,J,K,KOUNT,M,MU,N,NCAP,
     +        NU,NU0,NU1
C     ..
C     .. Local Arrays ..
      DOUBLE PRECISION A(NCAPM),A1(NMAX+1),ALPHA(NMAX+1),B(NCAPM),
     +                 B1(NMAX+1),BE(NMAX+1),BETA(NMAX+1),E(NCAPM),
     +                 P0(NCAPM),P1(NCAPM),P2(NCAPM),RHOI(1),RHOR(1),
     +                 ROLDI(1),ROLDR(1),W(NCAPM),WG(NMAX),X(NCAPM),
     +                 XII(2*NMAX),XIR(2*NMAX),ZG(NMAX)
      INTEGER IS(2*NMAX)
C     ..
C     .. External Functions ..
      DOUBLE PRECISION D1MACH,F7
      EXTERNAL D1MACH,F7
C     ..
C     .. External Subroutines ..
      EXTERNAL DABMOD,DGCHRS,DGQRAT,DKNUM,DRECUR
C     ..
C     .. Intrinsic Functions ..
      INTRINSIC ABS,ATAN,DBLE,EXP,REAL,SQRT
      integer i1mach, nout

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST07'
      write ( *, '(a)' ) '  Test ???'

      nout = i1mach(2)
C     ..
      PI = 4.D0*ATAN(1.D0)
      SQPI = SQRT(PI)
      C1 = PI
c      c1=-pi
c      c1=pi
      C2 = SQPI
c      c2=.5*sqpi
c      c2=.75*sqpi
      IROUT = 1
c      irout=0
      EPS = D1MACH(3)*1.D3
      CALL DRECUR(NCAPM,7,AK-1.D0,0.D0,A,B,IERR)
      WRITE(NOUT,FMT=9000)
      DO 60 N = 1,NMAX
          M = 2*N - 1
c        m=2*((n+1)/2)-1
c        m=1
          IF (M.GE.3) THEN
              DO MU = 1,M - 2,2
                  XIR(MU) = 0.D0
                  XIR(MU+1) = 0.D0
                  XII(MU) = 1.D0/ (DBLE(MU+1)*PI)
                  XII(MU+1) = -XII(MU)
                  IS(MU) = 1
                  IS(MU+1) = 1
              end do
          END IF

          CALL DABMOD(N+1,NCAPM,M-1,EPS,IROUT,A,B,XIR,XII,IS,A1,B1,NCAP,
     +                KOUNT,IERAB,BE,X,W,E,P0,P1,P2)
          IF (IERAB.NE.0) WRITE(NOUT,FMT=9010) IERAB
          ETA = -.001D0
          XIR(M) = -1.D0/ETA
          XII(M) = 0.D0
          IS(M) = 1
          X1 = ETA
          XA = ABS(ETA)
          T = 1.D0
          SUM0 = 0.D0
          SUM = -1.D0/ (AK-1.D0)
          J = 0
   20     J = J + 1
          IF (J.GT.100) THEN
              WRITE(NOUT,FMT=9020)
              GO TO 60
          END IF

          T = -XA*T/DBLE(J)
          SUM = SUM + T/ (DBLE(J+1)-AK)
          IF (SUM.NE.SUM0) THEN
              SUM0 = SUM
              GO TO 20
          END IF

          HP = -EXP(XA)* (C1* (XA** (AK-1.D0))-C2*SUM)
          HN = 0.D0
          IF (M.GE.3) THEN
              SUM = 0.D0
              DO 40 NU = 1, (M-1)/2
                  FNU = DBLE(NU)
                  FNU2 = FNU**2
                  P = 1.D0
                  DO 30 MU = 1, (M-1)/2
                      IF (MU.EQ.NU) GO TO 30
                      FMU = DBLE(MU)
                      FMU2 = FMU**2
                      P = (FMU2/ (FMU2-FNU2))*P
   30             CONTINUE
                  NU0 = 10
                  DX = 0.D0
                  DY = 2.D0*FNU*PI
                  CALL DKNUM(0,NU0,NUMAX,DX,DY,EPS,A,B,RHOR,RHOI,NU1,
     +                       IERK,ROLDR,ROLDI)
                  IF (IERK.NE.0) WRITE(NOUT,FMT=9030) IERK
                  SUM = SUM + 2.D0*FNU*PI*P*
     +                  (2.D0*FNU*PI*HP- (2.D0*FNU*PI*RHOR(1)+
     +                  X1*RHOI(1)))/ (X1**2+4.D0*FNU2* (PI**2))
   40         CONTINUE
              HP = SUM
          END IF

          CALL DGCHRS(N+1,1,A1,B1,X1,HP,HN,ALPHA,BETA,IERGC)
          IF (IERGC.NE.0) WRITE(NOUT,FMT=9040) IERGC
          BETA(1) = -ETA*BETA(1)
          CALL DGQRAT(N,M,ALPHA,BETA,XIR,XII,IS,ZG,WG,CONST,IERGQ,E)
          IF (IERGQ.NE.0) WRITE(NOUT,FMT=9050) IERGQ
          SUM = 0.D0
          DO 50 K = 1,N
              SUM = SUM + WG(K)*F7(ZG(K),ETA,THETA)
   50     CONTINUE
          EXACT = 2.2171501009112D0
c        exact=1.7800123286291d0
c        exact=3.7403844583705d0
          ERR = ABS(REAL((SUM-EXACT)/EXACT))
          WRITE(NOUT,FMT=9060) N,M,SUM,ERR,CONST
   60 CONTINUE
      return

 9000 FORMAT (5X,'n',4X,'m',8X,'integral',10X,'rel error',4X,
     +       'err const',/)
 9010 FORMAT (1X,'ierab=',I3)
 9020 FORMAT (1X,'power series for gamma does not converge')
 9030 FORMAT (1X,'ierk for complex z =',I3)
 9040 FORMAT (1X,'iergc=',I3)
 9050 FORMAT (1X,'iergq=',I3)
 9060 FORMAT (1X,2I5,D23.14,E13.4,D13.4)
      END
      FUNCTION F7(T,ETA,THETA)

c*********************************************************************72
c
cc F7 is the integrand for TEST07.
c
c  Modified:
c
c    08 November 2012
c
c  Author:
c
c    Walter Gautschi
c
c  Reference:
c
c    Walter Gautschi,
c    Algorithm 793: GQRAT, Gauss Quadrature for Rational Functions,
c    ACM Transactions on Mathematical Software,
c    Volume 25, Number 2, pages 213-239, June 1999.
c
      implicit none

      DOUBLE PRECISION ETA,F7,T,THETA
      DOUBLE PRECISION S,S1,TERM,X
      INTEGER L
      INTRINSIC ABS,DBLE,EXP,SQRT
      integer i1mach, nout

      X = T - ETA

      IF (ABS(X).LE.1.D0) THEN

          L = 0
          TERM = X
          S1 = TERM
   10     continue
          S = S1
          L = L + 1

          IF ( 200 .lt. L ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'F7 - Fatal error!'
            write ( *, '(a)' ) '  Exponential series not converging.'
            stop
          end if

          TERM = X*TERM/DBLE(L+1)
          S1 = S + TERM
          IF (S1.NE.S) then
            GO TO 10
          end if
          F7 = T*SQRT(1.D0+.5D0*THETA*T)/ (EXP(-T)*S)
      ELSE
          F7 = T*SQRT(1.D0+.5D0*THETA*T)/ (EXP(-ETA)-EXP(-T))
      END IF

      return
      END
