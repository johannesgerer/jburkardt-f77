      program main

c*********************************************************************72
c
cc TOMS429_PRB tests POLYAN, ACM TOMS algorithm 429.
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
      write ( *, '(a)' ) 'TOMS429_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS429 library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS429_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests POLYAN with a polynomial suggested by Driessen and Hunt.
c
c  Discussion:
c
c    The original version of POLYAN, when given this polynomial, reported:
c
c    Roots are in annulus of inner radius 0.454E+00 and outer radius 0.836E+01,
c    there are no real positive roots,
c    the negative roots (if any) are between -.454E+00 and -0.836D+01,
c    there are no roots with positive real parts.
c
c    But this is incorrect for the given polynomial.
c
c  Reference:
c
c    H B Driessen and E W Hunt,
c    Remark on Algorithm 429,
c    Communications of the ACM,
c    Volume 16, Number 9, page 579, September 1973.
c
      implicit none

      integer n
      parameter ( n = 4 )

      real c(n)
      real cm(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  p(x) = x^4 + 5.6562x^3 + 5.8854x^2'
      write ( *, '(a)' ) '             + 7.3646x   + 6.1354'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Approximate roots:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  -1.001,'
      write ( *, '(a)' ) '  -4.7741, '
      write ( *, '(a)' ) '   0.0089 + 1.1457 i,'
      write ( *, '(a)' ) '   0.0089 - 1.1457 i.'
      write ( *, '(a)' ) ' '

      c(1) = 5.6562E+00
      c(2) = 5.8854E+00
      c(3) = 7.3646E+00
      c(4) = 6.1354E+00

      call polyan ( c, cm, n )

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests POLYAN with a polynomial with a single root.
c
      implicit none

      integer n
      parameter ( n = 4 )

      real c(n)
      real cm(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  p(x) = x^4 - 8 x^3 + 24 x^2 - 32 x + 16'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Exact roots:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  2 (multiplicity 4)'
      write ( *, '(a)' ) ' '

      c(1) =  -8.0E+00
      c(2) =  24.0E+00
      c(3) = -32.0E+00
      c(4) =  16.0E+00

      call polyan ( c, cm, n )

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests POLYAN with a polynomial in the roots of unity.
c
      implicit none

      integer n
      parameter ( n = 5 )

      real c(n)
      real cm(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  p(x) = x^5 - 1'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Exact roots:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  0.3090 + 0.9510 i'
      write ( *, '(a)' ) ' -0.8090 + 0.5877 i'
      write ( *, '(a)' ) ' -0.8090 - 0.5877 i'
      write ( *, '(a)' ) '  0.3090 - 0.9510 i'
      write ( *, '(a)' ) '  1'
      write ( *, '(a)' ) ' '

      c(1) =  0.0E+00
      c(2) =  0.0E+00
      c(3) =  0.0E+00
      c(4) =  0.0E+00
      c(5) = -1.0E+00

      call polyan ( c, cm, n )

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
