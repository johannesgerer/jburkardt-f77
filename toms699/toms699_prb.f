      program main

C      ALGORITHM 699 , COLLECTED ALGORITHMS FROM ACM.
C      THIS WORK PUBLISHED IN TRANSACTIONS ON MATHEMATICAL SOFTWARE,
C      VOL. 17, NO. 4, DECEMBER, 1991, PP. 457-461.
      program TEQUAD
c
c     Test the re-implementation of Patterson's QUAD.
c
c     Parameters
c
c NEPS    is the number of tolerances at which to exercise QUAD.
c NKASES  is the number of integrands to use to exercise QUAD.
c PI      is pi.
c
      integer NEPS, NKASES
      parameter (NEPS=3, NKASES=13)
      real PI
      parameter (PI = 3.1415 92653 58979 32384 62643)
c
c     Arguments of QUAD.
c
      real A, B
      real EPSIL
      real F
      external F
      integer ICHECK, K, NPTS
      real RESULT(8)
c
c     Other variables
c
c ANSWRS  is an array of RESULT(K) found by calling QUAD; indexed by
c         the tolerance index.
c AS, BS  are arrays of lower and upper bounds, indexed by KASE, which
c         is in common /COM/.
c EPSILS  are values to assign to EPSIL, indexed by I.
c I       is a loop inductor and subscript.
c ICHEKS  is an array of ICHECK values returned by QUAD; indexed by
c         the tolerance index.
c NPTSS   is an array of NPTS values returned by QUAD; indexed by
c         the tolerance index.
c
      real ANSWRS(NEPS), AS(NKASES), BS(NKASES), EPSILS(NEPS)
      integer I, ICHEKS(NEPS), NPTSS(NEPS)
c
c     Common Variables
c
c KASE    tells the integrand function what case to work on.
c
      integer KASE
      common /COM/ KASE
c
c     Data
c
      data AS /0.0, -1.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0,
     1         0.0, 0.0/
      data BS /1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 10.0, pi, 1.0,
     1         1.0, 1.0/
      data EPSILS /1.e-3, 1.e-6, 1.e-8/
c
c     Executable Statements
c
      print 10, epsils
10    format (' eps=',5(:'------- ',1pg9.1,' ------ '))
      print 20, (' ',i=1,neps)
20    format (' KASE',5(:' NPTS     RESULT    CHK ',a1))
      do 100 kase = 1, nkases
        a = as(kase)
        b = bs(kase)
        do 50 i = 1, neps
          epsil = epsils(i)
          call quad (a, b, result, k, epsil, npts, icheck, f)
          answrs(i)=result(k)
          nptss(i)=npts
          icheks(i)=icheck
50      continue
        print 70, kase, (nptss(i),answrs(i),icheks(i),i=1,neps)
70      format (i5,5(i4,1pe16.8,i2,3x))
100   continue
      stop
      end
c     =====     Integrand Function     =================================
      real function F (X)
      real X
      real PI
      parameter (PI = 3.1415 92653 58979 32384 62643)
      real T
      integer KASE
      common /COM/ KASE
      go to (1,2,3,4,5,6,7,8,9,10,11,12,13), kase
1     f=sqrt(x)
      return
2     f=0.92*cosh(x)-cos(x)
      return
3     t=x*x
      f=1.0/(t*(t+1.0)+0.9)
      return
4     f=sqrt(x)**3
      return
5     f=1.0/(1.0+x**4)
      return
6     f=1.0/(1.0+0.5*sin(10.0*pi*x))
      return
7     f=x/(exp(x)-1.0)
      return
8     f=sin(100.0*pi*x)/(pi*x)
      return
9     f=50.0/((2500.0*x*x+1.0)*pi)
      return
10    t=x+x
      f=cos(cos(x)+3.0*(sin(x)+cos(3.0*x)+sin(t))+2.0*cos(t))
      return
11    f=log(x)
      return
12    f=4.0*pi*pi*x*sin(20.0*pi*x)*cos(2.0*pi*x)
      return
13    f=1.0/(1.0+(230.0*x-30.0)**2)
      return
      end
