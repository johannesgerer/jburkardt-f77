      program main

c*********************************************************************72
c
cc TOMS493_PRB tests RPOLY, ACM TOMS algorithm 493.
c
c  Modified:
c
c    30 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS493_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS493 library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS493_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 uses a polynomial suggested by Driessen and Hunt.
c
c  Modified:
c
c    30 December 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    HB Driessen, EW Hunt,
c    Remark on Algorithm 429,
c    Communications of the ACM,
c    Volume 16, Number 9, page 579, September 1973.
c
      implicit none

      integer degree
      parameter ( degree = 4 )

      double precision c(degree+1)
      logical fail
      integer i
      double precision root_i(degree)
      double precision root_r(degree)

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

      c(1) = 1.0D+00
      c(2) = 5.6562D+00
      c(3) = 5.8854D+00
      c(4) = 7.3646D+00
      c(5) = 6.1354D+00

      call rpoly ( c, degree, root_r, root_i, fail )

      do i = 1, degree

        write ( *, '(2x,g14.6,2x,g14.6)' ) root_r(i), root_i(i)

      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests a polynomial with a single root.
c
c  Modified:
c
c    30 December 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer degree
      parameter ( degree = 4 )

      double precision c(degree+1)
      logical fail
      integer i
      double precision root_i(degree)
      double precision root_r(degree)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  p(x) = x^4 - 8 x^3 + 24 x^2 - 32 x + 16'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Exact roots:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  2 (multiplicity 4)'
      write ( *, '(a)' ) ' '

      c(1) =   1.0D+00
      c(2) =  -8.0D+00
      c(3) =  24.0D+00
      c(4) = -32.0D+00
      c(5) =  16.0D+00

      call rpoly ( c, degree, root_r, root_i, fail )

      do i = 1, degree

        write ( *, '(2x,g14.6,2x,g14.6)' ) root_r(i), root_i(i)

      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests a polynomial in the roots of unity.
c
c  Modified:
c
c    30 December 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer degree
      parameter ( degree = 5 )

      double precision c(degree+1)
      logical fail
      integer i
      double precision root_i(degree)
      double precision root_r(degree)

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

      c(1) =  1.0D+00
      c(2) =  0.0D+00
      c(3) =  0.0D+00
      c(4) =  0.0D+00
      c(5) =  0.0D+00
      c(6) = -1.0D+00

      call rpoly ( c, degree, root_r, root_i, fail )

      do i = 1, degree

        write ( *, '(2x,g14.6,2x,g14.6)' ) root_r(i), root_i(i)

      end do

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
