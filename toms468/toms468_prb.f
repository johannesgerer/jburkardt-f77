      program main

c*********************************************************************72
c
cc TOMS468_PRB tests TOMS468.
c
c  Modified:
c
c    12 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS468_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS468 library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS468_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests QUAD.
c
      implicit none

      real a
      real b
      real epsil
      real exact
      external f01
      external f02
      external f03
      external f04
      external f05
      external f06
      external f07
      external f08
      external f09
      external f10
      external f11
      external f12
      external f13
      integer icheck
      integer k
      integer npts
      real pi
      real result(8)

      pi = 3.141592653589793E+00
      epsil = 0.001E+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Test QUAD, for simple quadrature.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Error tolerance EPSIL = ', epsil
      write ( *, '(a)' ) ' '

      write ( *, '(a,a)' )
     &  '      A         B   ICHECK K     NFUNC     ',
     &  'RESULT(K)        EXACT'
      write ( *, '(a)' ) ' '

      a = 0.0E+00
      b = 1.0E+00

      call quad ( a, b, result, k, epsil, npts, icheck, f01 )

      exact = 2.0E+00 / 3.0E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i2,2x,i8,2x,g14.6,2x,g14.6)' )
     &  a, b, icheck, k, npts, result(k), exact

      a = -1.0E+00
      b =  1.0E+00

      call quad ( a, b, result, k, epsil, npts, icheck, f02 )

      exact = 0.4794282267E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i2,2x,i8,2x,g14.6,2x,g14.6)' )
     &  a, b, icheck, k, npts, result(k), exact

      a = -1.0E+00
      b =  1.0E+00

      call quad ( a, b, result, k, epsil, npts, icheck, f03 )

      exact = 1.582232964E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i2,2x,i8,2x,g14.6,2x,g14.6)' )
     &  a, b, icheck, k, npts, result(k), exact

      a = 0.0E+00
      b = 1.0E+00

      call quad ( a, b, result, k, epsil, npts, icheck, f04 )

      exact = 2.0E+00 / 5.0E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i2,2x,i8,2x,g14.6,2x,g14.6)' )
     &  a, b, icheck, k, npts, result(k), exact

      a = 0.0E+00
      b = 1.0E+00

      call quad ( a, b, result, k, epsil, npts, icheck, f05 )

      exact = 0.8669729873E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i2,2x,i8,2x,g14.6,2x,g14.6)' )
     &  a, b, icheck, k, npts, result(k), exact

      a = 0.0E+00
      b = 1.0E+00

      call quad ( a, b, result, k, epsil, npts, icheck, f06 )

      exact = 1.154700669E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i2,2x,i8,2x,g14.6,2x,g14.6)' )
     &  a, b, icheck, k, npts, result(k), exact

      a = 0.0E+00
      b = 1.0E+00

      call quad ( a, b, result, k, epsil, npts, icheck, f07 )

      exact = 0.7775046341E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i2,2x,i8,2x,g14.6,2x,g14.6)' )
     &  a, b, icheck, k, npts, result(k), exact

      a = 0.1E+00
      b = 1.0E+00

      call quad ( a, b, result, k, epsil, npts, icheck, f08 )

      exact = 0.009098645256E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i2,2x,i8,2x,g14.6,2x,g14.6)' )
     &  a, b, icheck, k, npts, result(k), exact

      a =  0.0E+00
      b = 10.0E+00

      call quad ( a, b, result, k, epsil, npts, icheck, f09 )

      exact = 0.4993638029E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i2,2x,i8,2x,g14.6,2x,g14.6)' )
     &  a, b, icheck, k, npts, result(k), exact

      a = 0.0E+00
      b = pi

      call quad ( a, b, result, k, epsil, npts, icheck, f10 )

      exact = 0.8386763234E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i2,2x,i8,2x,g14.6,2x,g14.6)' )
     &  a, b, icheck, k, npts, result(k), exact

      a = 0.0E+00
      b = 1.0E+00

      call quad ( a, b, result, k, epsil, npts, icheck, f11 )

      exact = -1.0E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i2,2x,i8,2x,g14.6,2x,g14.6)' )
     &  a, b, icheck, k, npts, result(k), exact

      a = 0.0E+00
      b = 1.0E+00

      call quad ( a, b, result, k, epsil, npts, icheck, f12 )

      exact = -0.6346651825E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i2,2x,i8,2x,g14.6,2x,g14.6)' )
     &  a, b, icheck, k, npts, result(k), exact

      a = 0.0E+00
      b = 1.0E+00

      call quad ( a, b, result, k, epsil, npts, icheck, f13 )
c
c  The reference lists an exact value of 0.0013492485650E+00 but this is
c  apparently a typo.
c
      exact = 0.013492485650E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i2,2x,i8,2x,g14.6,2x,g14.6)' )
     &  a, b, icheck, k, npts, result(k), exact

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests QSUB.
c
      implicit none

      real a
      real b
      real epsil
      real exact
      external f01
      external f02
      external f03
      external f04
      external f05
      external f06
      external f07
      external f08
      external f09
      external f10
      external f11
      external f12
      external f13
      integer icheck
      integer npts
      real pi
      real qsub
      real relerr
      real result

      pi = 3.141592653589793E+00
      epsil = 0.001E+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Test QSUB, for quadrature with subdivision.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Error tolerance EPSIL = ', epsil
      write ( *, '(a)' ) ' '

      write ( *, '(a,a)' )
     &  '      A         B   ICHECK   NFUNC     ',
     &  '   RESULT        EXACT          RELERR'
      write ( *, '(a)' ) ' '

      a = 0.0E+00
      b = 1.0E+00

      result = qsub ( a, b, epsil, npts, icheck, relerr, f01 )

      exact = 2.0E+00 / 3.0E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' )
     &  a, b, icheck, npts, result, exact, relerr

      a = -1.0E+00
      b =  1.0E+00

      result = qsub ( a, b, epsil, npts, icheck, relerr, f02 )

      exact = 0.4794282267E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' )
     &  a, b, icheck, npts, result, exact, relerr

      a = -1.0E+00
      b =  1.0E+00

      result = qsub ( a, b, epsil, npts, icheck, relerr, f03 )

      exact = 1.582232964E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' )
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.0E+00
      b = 1.0E+00

      result = qsub ( a, b, epsil, npts, icheck, relerr, f04 )

      exact = 2.0E+00 / 5.0E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' )
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.0E+00
      b = 1.0E+00

      result = qsub ( a, b, epsil, npts, icheck, relerr, f05 )

      exact = 0.8669729873E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' )
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.0E+00
      b = 1.0E+00

      result = qsub ( a, b, epsil, npts, icheck, relerr, f06 )

      exact = 1.154700669E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' )
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.0E+00
      b = 1.0E+00

      result = qsub ( a, b, epsil, npts, icheck, relerr, f07 )

      exact = 0.7775046341E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' )
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.1E+00
      b = 1.0E+00

      result = qsub ( a, b, epsil, npts, icheck, relerr, f08 )

      exact = 0.009098645256E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' )
     &  a, b, icheck, npts, result, exact, relerr

      a =  0.0E+00
      b = 10.0E+00

      result = qsub ( a, b, epsil, npts, icheck, relerr, f09 )

      exact = 0.4993638029E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' )
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.0E+00
      b = pi

      result = qsub ( a, b, epsil, npts, icheck, relerr, f10 )

      exact = 0.8386763234E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' )
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.0E+00
      b = 1.0E+00

      result = qsub ( a, b, epsil, npts, icheck, relerr, f11 )

      exact = -1.0E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' )
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.0E+00
      b = 1.0E+00

      result = qsub ( a, b, epsil, npts, icheck, relerr, f12 )

      exact = -0.6346651825E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' )
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.0E+00
      b = 1.0E+00

      result = qsub ( a, b, epsil, npts, icheck, relerr, f13 )
c
c  The reference lists an exact value of 0.0013492485650E+00 but this is
c  apparently a typo.
c
      exact = 0.013492485650E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' )
     &  a, b, icheck, npts, result, exact, relerr

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests QSUBA.
c
      implicit none

      real a
      real b
      real epsil
      real exact
      external f01
      external f02
      external f03
      external f04
      external f05
      external f06
      external f07
      external f08
      external f09
      external f10
      external f11
      external f12
      external f13
      integer icheck
      integer npts
      real pi
      real qsuba
      real relerr
      real result

      pi = 3.141592653589793E+00
      epsil = 0.001E+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  Test QSUBA,'
      write ( *, '(a)' ) '  for adaptive quadrature with subdivision.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Error tolerance EPSIL = ', epsil
      write ( *, '(a)' ) ' '

      write ( *, '(a,a)' )
     &  '      A         B   ICHECK   NFUNC     ',
     &  '   RESULT        EXACT          RELERR'
      write ( *, '(a)' ) ' '

      a = 0.0E+00
      b = 1.0E+00

      result = qsuba ( a, b, epsil, npts, icheck, relerr, f01 )

      exact = 2.0E+00 / 3.0E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' )
     &  a, b, icheck, npts, result, exact, relerr

      a = -1.0E+00
      b =  1.0E+00

      result = qsuba ( a, b, epsil, npts, icheck, relerr, f02 )

      exact = 0.4794282267E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' )
     &  a, b, icheck, npts, result, exact, relerr

      a = -1.0E+00
      b =  1.0E+00

      result = qsuba ( a, b, epsil, npts, icheck, relerr, f03 )

      exact = 1.582232964E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' )
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.0E+00
      b = 1.0E+00

      result = qsuba ( a, b, epsil, npts, icheck, relerr, f04 )

      exact = 2.0E+00 / 5.0E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' )
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.0E+00
      b = 1.0E+00

      result = qsuba ( a, b, epsil, npts, icheck, relerr, f05 )

      exact = 0.8669729873E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' )
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.0E+00
      b = 1.0E+00

      result = qsuba ( a, b, epsil, npts, icheck, relerr, f06 )

      exact = 1.154700669E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' )
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.0E+00
      b = 1.0E+00

      result = qsuba ( a, b, epsil, npts, icheck, relerr, f07 )

      exact = 0.7775046341E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' )
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.1E+00
      b = 1.0E+00

      result = qsuba ( a, b, epsil, npts, icheck, relerr, f08 )

      exact = 0.009098645256E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' )
     &  a, b, icheck, npts, result, exact, relerr

      a =  0.0E+00
      b = 10.0E+00

      result = qsuba ( a, b, epsil, npts, icheck, relerr, f09 )

      exact = 0.4993638029E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' )
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.0E+00
      b = pi

      result = qsuba ( a, b, epsil, npts, icheck, relerr, f10 )

      exact = 0.8386763234E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' )
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.0E+00
      b = 1.0E+00

      result = qsuba ( a, b, epsil, npts, icheck, relerr, f11 )

      exact = -1.0E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' )
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.0E+00
      b = 1.0E+00

      result = qsuba ( a, b, epsil, npts, icheck, relerr, f12 )

      exact = -0.6346651825E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' )
     &  a, b, icheck, npts, result, exact, relerr

      a = 0.0E+00
      b = 1.0E+00

      result = qsuba ( a, b, epsil, npts, icheck, relerr, f13 )
c
c  The reference lists an exact value of 0.0013492485650E+00 but this is
c  apparently a typo.
c
      exact = 0.013492485650E+00

      write ( *,
     &  '(2x,f8.4,2x,f8.4,2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g12.4)' )
     &  a, b, icheck, npts, result, exact, relerr

      return
      end
      function f01 ( x )

c*********************************************************************72
c
cc F01 evaluates test function 1.
c
      implicit none

      real f01
      real x

      f01 = sqrt ( x )

      return
      end
      function f02 ( x )

c*********************************************************************72
c
cc F02 evaluates test function 2.
c
      implicit none

      real f02
      real x

      f02 = 0.92E+00 * cosh ( x ) - cos ( x )

      return
      end
      function f03 ( x )

c*********************************************************************72
c
cc F03 evaluates test function 3.
c
      implicit none

      real f03
      real x

      f03 = 1.0E+00 / ( x**4 + x**2 + 0.9E+00 )

      return
      end
      function f04 ( x )

c*********************************************************************72
c
cc F04 evaluates test function 4.
c
      implicit none

      real f04
      real x

      f04 = sqrt ( x**3 )

      return
      end
      function f05 ( x )

c*********************************************************************72
c
cc F05 evaluates test function 5.
c
      implicit none

      real f05
      real x

      f05 = 1.0E+00 / ( 1.0E+00 + x**4 )

      return
      end
      function f06 ( x )

c*********************************************************************72
c
cc F06 evaluates test function 6.
c
      implicit none

      real f06
      real x

      f06 = 1.0E+00 / ( 1.0E+00 + 0.5E+00 * sin ( 31.4159E+00 * x ) )

      return
      end
      function f07 ( x )

c*********************************************************************72
c
cc F07 evaluates test function 7.
c
c  Modified:
c
c    12 January 2006
c
      implicit none

      real f07
      real x

      if ( x .eq. 0.0E+00 ) then
        f07 = 1.0E+00
      else
        f07 = x / ( exp ( x ) - 1.0E+00 )
      end if

      return
      end
      function f08 ( x )

c*********************************************************************72
c
cc F08 evaluates test function 8.
c
      implicit none

      real f08
      real pi
      real x

      pi = 3.141592653589793E+00

      f08 = sin ( 100.0E+00 * pi * x ) / ( pi * x )

      return
      end
      function f09 ( x )

c*********************************************************************72
c
cc F09 evaluates test function 9.
c
      implicit none

      real f09
      real pi
      real x

      pi = 3.141592653589793E+00

      f09 = 50.0E+00 / ( 2500.0E+00 * x**2 + 1.0E+00 ) / pi

      return
      end
      function f10 ( x )

c*********************************************************************72
c
cc F10 evaluates test function 10.
c
      implicit none

      real arg
      real f10
      real x

      arg =         cos (           x )
     &  + 3.0E+00 * sin (           x )
     &  + 2.0E+00 * cos ( 2.0E+00 * x )
     &  + 3.0E+00 * cos ( 3.0E+00 * x )
     &  + 3.0E+00 * sin ( 2.0E+00 * x )

      f10 = cos ( arg )

      return
      end
      function f11 ( x )

c*********************************************************************72
c
cc F11 evaluates test function 11.
c
      implicit none

      real f11
      real x

      if ( x .le. 0.0E+00 ) then
       f11 = 0.0E+00
      else
        f11 = log ( x )
      end if

      return
      end
      function f12 ( x )

c*********************************************************************72
c
cc F12 evaluates test function 12.
c
      implicit none

      real f12
      real pi
      real x

      pi = 3.141592653589793E+00

      f12 = 4.0E+00 * pi**2 * x * sin ( 20.0E+00 * pi * x )
     &  * cos ( 2.0E+00 * pi * x )

      return
      end
      function f13 ( x )

c*********************************************************************72
c
cc F13 evaluates test function 13.
c
      implicit none

      real f13
      real x

      f13 = 1.0E+00 / ( 1.0E+00 + ( 230.0E+00 * x - 30.0E+00 )**2 )

      return
      end
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
c  Modified:
c
c    16 September 2005
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

      character ( len = 8 ) date
      character ( len = 10 ) time

      call date_and_time ( date, time )

      write ( *, '(a8,2x,a10)' ) date, time

      return
      end
