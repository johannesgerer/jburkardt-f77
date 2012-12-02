      program main

c*********************************************************************72
c
cc TOMS424_PRB tests CCQUAD.
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

      real a
      real b
      external f1d1
      external fbd1
      external fed1
      external fqd1
      external fxd1
      external fx2d1
      external fx3d1

      a = 0.0E+00
      b = 1.0E+00

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS424_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS424 library.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  CCQUAD, the ACM TOMS algorithm 424'
      write ( *, '(a)' ) '  for Clenshaw-Curtis quadrature to estimate'
      write ( *, '(a)' ) '  the integral of F(X) on [A,B].'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g24.16)' ) '  A = ', a
      write ( *, '(a,g24.16)' ) '  B = ', b

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  F(X)=1'

      call test01 ( f1d1, a, b )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  F(X)=X'

      call test01 ( fxd1, a, b )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  F(X)=X*2'

      call test01 ( fx2d1, a, b )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  F(X)=X*3'

      call test01 ( fx3d1, a, b )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  F(X)=EXP(X)'

      call test01 ( fed1, a, b )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  F(X)=SQRT(X)'

      call test01 ( fqd1, a, b )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  F(X)=1/(1+X*X)'

      call test01 ( fbd1, a, b )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS424_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      return
      end
      subroutine test01 ( f, a, b )

c*********************************************************************72
c
cc TEST01 tests CCQUAD with a given function.
c
      implicit none

      integer limit
      parameter ( limit = 1000 )

      real a
      real b
      real ccquad
      real csxfrm(limit)
      real esterr
      external f
      integer i
      real quad
      real tolerr
      integer used

      tolerr = 1.0E-06

      quad = ccquad ( f, a, b, tolerr, limit, esterr, used,
     *  csxfrm )

      write ( *, '(a,g24.16)' ) '  Integral estimate: ', quad
      write ( *, '(a,g24.16)' ) '  Error estimate:    ', esterr
      write ( *, '(a,i6)'     ) '  Evaluations of F = ', used

      if ( used < 20 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Scaled Discrete Cosine Transform:'
        write ( *, '(a)' ) ' '
        do i = 1, used
          write ( *, '(g24.16)' ) csxfrm(i)
        end do
      end if

      return
      end
      function f1d1 ( x )

c*********************************************************************72
c
cc F1D1(X) = 1.
c
      implicit none

      real f1d1
      real x

      f1d1 = 1.0E+00

      return
      end
      function fbd1 ( x )

c*********************************************************************72
c
cc FBD1(X) = 1 / ( 1 + X**2 )
c
      implicit none

      real fbd1
      real x

      fbd1 = 1.0E+00 / ( 1.0E+00 + x**2 )

      return
      end
      function fed1 ( x )

c*********************************************************************72
c
cc FED1(X) = EXP(X).
c
      implicit none

      real fed1
      real x

      fed1 = exp ( x )

      return
      end
      function fqd1 ( x )

c*********************************************************************72
c
cc FQD1(X) = SQRT ( X ).
c
      implicit none

      real fqd1
      real x

      fqd1 = sqrt ( abs ( x ) )

      return
      end
      function fxd1 ( x )

c*********************************************************************72
c
cc FXD1(X) = X.
c
      implicit none

      real fxd1
      real x

      fxd1 = x

      return
      end
      function fx2d1 ( x )

c*********************************************************************72
c
cc FX2D1(X) = X**2.
c
      implicit none

      real fx2d1
      real x

      fx2d1 = x**2

      return
      end
      function fx3d1 ( x )

c*********************************************************************72
c
cc FX3D1(X) = X**3.
c
      implicit none

      real fx3d1
      real x

      fx3d1 = x**3

      return
      end
      function fxsd1 ( x )

c*********************************************************************72
c
cc FXSD1(X) = X / SQRT ( 1 - X**2 )
c
      implicit none

      real fxsd1
      real x

      fxsd1 = x / sqrt ( 1.0E+00 - x**2 )

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
