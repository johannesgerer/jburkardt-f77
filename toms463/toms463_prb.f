      program main

c*********************************************************************72
c
cc TOMS463_PRB tests algorithm 463.
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
      write ( *, '(a)' ) 'TOMS463_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS463 library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS463_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests SCALE1.
c
      implicit none

      integer test_num

      parameter ( test_num = 3 )

      real dist
      integer n
      integer n_actual
      integer n_test(test_num)
      integer test
      real xmax
      real xmax_test(test_num)
      real xmaxp
      real xmin
      real xmin_test(test_num)
      real xminp

      save n_test
      save xmax_test
      save xmin_test

      data n_test / 5, 5, 9 /
      data xmax_test / 11.1, 10.1, -100.0 /
      data xmin_test / -3.1, 5.2, -12000.0 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  SCALE1 chooses a scale for a plot.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '       XMIN       XMAX   N      XMINP',
     &  '      XMAXP       DIST Nactual'
      write ( *, '(a)' ) ' '

      do test = 1, test_num

        n = n_test(test)
        xmin = xmin_test(test)
        xmax = xmax_test(test)

        call scale1 ( xmin, xmax, n, xminp, xmaxp, dist )

        n_actual = nint ( ( xmaxp - xminp ) / dist )

        write ( *,
     &    '(2x,f9.1,2x,f9.1,2x,i2,2x,f9.1,2x,f9.1,2x,f9.1,2x,i6)' )
     &    xmin, xmax, n, xminp, xmaxp, dist, n_actual

      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests SCALE2.
c
      implicit none

      integer test_num

      parameter ( test_num = 3 )

      real dist
      integer n
      integer n_actual
      integer n_test(test_num)
      integer test
      real xmax
      real xmax_test(test_num)
      real xmaxp
      real xmin
      real xmin_test(test_num)
      real xminp

      save n_test
      save xmax_test
      save xmin_test

      data n_test / 5, 5, 9 /
      data xmax_test / 11.1, 10.1, -100.0 /
      data xmin_test / -3.1, 5.2, -12000.0 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  SCALE2 chooses a scale for a plot.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '       XMIN       XMAX   N      XMINP',
     &  '      XMAXP       DIST Nactual'
      write ( *, '(a)' ) ' '

      do test = 1, test_num

        n = n_test(test)
        xmin = xmin_test(test)
        xmax = xmax_test(test)

        call scale2 ( xmin, xmax, n, xminp, xmaxp, dist )

        n_actual = nint ( ( xmaxp - xminp ) / dist )

        write ( *,
     &    '(2x,f9.1,2x,f9.1,2x,i2,2x,f9.1,2x,f9.1,2x,f9.1,2x,i6)' )
     &    xmin, xmax, n, xminp, xmaxp, dist, n_actual

      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests SCALE3.
c
      implicit none

      integer test_num

      parameter ( test_num = 3 )

      real dist
      integer n
      integer n_actual
      integer n_test(test_num)
      integer test
      real xmax
      real xmax_test(test_num)
      real xmaxp
      real xmin
      real xmin_test(test_num)
      real xminp

      save n_test
      save xmax_test
      save xmin_test

      data n_test / 10, 2, 4 /
      data xmax_test / 125.0, 10.0, 1500.0 /
      data xmin_test / 1.8, 0.1, 0.1 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  SCALE3 chooses a logarithmic scale'
      write ( *, '(a)' ) '  for a plot.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '       XMIN       XMAX   N      XMINP',
     &  '      XMAXP       DIST'
      write ( *, '(a)' ) ' '

      do test = 1, test_num

        n = n_test(test)
        xmin = xmin_test(test)
        xmax = xmax_test(test)

        call scale3 ( xmin, xmax, n, xminp, xmaxp, dist )

        write ( *,
     &    '(2x,f9.1,2x,f9.1,2x,i2,2x,f11.3,2x,f11.3,2x,f11.3)' )
     &    xmin, xmax, n, xminp, xmaxp, dist

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
