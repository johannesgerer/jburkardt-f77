      program main

c*********************************************************************72
c
cc TOMS379_PRB tests TOMS379.
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
      write ( *, '(a)' ) 'TOMS379_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS379 library.'

      call test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS379_PRB'
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

      double precision a
      double precision big
      double precision error
      double precision exact
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
      double precision fifth
      integer i
      integer j
      integer no
      double precision pi
      double precision result
      double precision rum
      double precision squank

      pi = 3.141592653589793D+00
      error = 0.1D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Test SQUANK,'
      write ( *, '(a)' ) '  for adaptive Simpson quadrature.'

      do i = 1, 3

        error = error / 100.0D+00

        write ( *, '(a)' ) ' '
        write ( *, '(a,g14.6)' ) '  Error tolerance ERROR = ', error
        write ( *, '(a)' ) ' '

        write ( *, '(a,a)' )
     &  '            NO      SQUANK',
     &  '          EXACT           ERROR       RUM'
        write ( *, '(a)' ) ' '

        j = 1
        a = 0.0D+00
        big = 1.0D+00

        result = squank ( a, big, error, fifth, rum, no, f01 )

        exact = 2.0D+00 / 3.0D+00

        write ( *,
     &    '(2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g10.2,2x,g10.2)' )
     &    j, no, result, exact, abs ( result - exact ), rum

        j = 2
        a = -1.0D+00
        big =  1.0D+00

        result = squank ( a, big, error, fifth, rum, no, f02 )

        exact = 0.4794282267D+00

        write ( *,
     &    '(2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g10.2,2x,g10.2)' )
     &    j, no, result, exact, abs ( result - exact ), rum

        j = 3
        a = -1.0D+00
        big =  1.0D+00

        result = squank ( a, big, error, fifth, rum, no, f03 )

        exact = 1.582232964D+00

        write ( *,
     &    '(2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g10.2,2x,g10.2)' )
     &    j, no, result, exact, abs ( result - exact ), rum

        j = 4
        a = 0.0D+00
        big = 1.0D+00

        result = squank ( a, big, error, fifth, rum, no, f04 )

        exact = 2.0D+00 / 5.0D+00

        write ( *,
     &    '(2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g10.2,2x,g10.2)' )
     &    j, no, result, exact, abs ( result - exact ), rum

        j = 5
        a = 0.0D+00
        big = 1.0D+00

        result = squank ( a, big, error, fifth, rum, no, f05 )

        exact = 0.8669729873D+00

        write ( *,
     &    '(2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g10.2,2x,g10.2)' )
     &    j, no, result, exact, abs ( result - exact ), rum

        j = 6
        a = 0.0D+00
        big = 1.0D+00

        result = squank ( a, big, error, fifth, rum, no, f06 )

        exact = 1.154700669D+00

        write ( *,
     &    '(2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g10.2,2x,g10.2)' )
     &    j, no, result, exact, abs ( result - exact ), rum

        j = 7
        a = 0.0D+00
        big = 1.0D+00

        result = squank ( a, big, error, fifth, rum, no, f07 )

        exact = 0.7775046341D+00

        write ( *,
     &    '(2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g10.2,2x,g10.2)' )
     &    j, no, result, exact, abs ( result - exact ), rum

        j = 8
        a = 0.1D+00
        big = 1.0D+00

        result = squank ( a, big, error, fifth, rum, no, f08 )

        exact = 0.009098645256D+00

        write ( *,
     &    '(2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g10.2,2x,g10.2)' )
     &    j, no, result, exact, abs ( result - exact ), rum

        j = 9
        a =  0.0D+00
        big = 10.0D+00

        result = squank ( a, big, error, fifth, rum, no, f09 )

        exact = 0.4993638029D+00

        write ( *,
     &    '(2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g10.2,2x,g10.2)' )
     &    j, no, result, exact, abs ( result - exact ), rum

        j = 10
        a = 0.0D+00
        big = pi

        result = squank ( a, big, error, fifth, rum, no, f10 )

        exact = 0.8386763234D+00

        write ( *,
     &    '(2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g10.2,2x,g10.2)' )
     &    j, no, result, exact, abs ( result - exact ), rum

        j = 11
        a = 0.0D+00
        big = 1.0D+00

        result = squank ( a, big, error, fifth, rum, no, f11 )

        exact = -1.0D+00

        write ( *,
     &    '(2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g10.2,2x,g10.2)' )
     &    j, no, result, exact, abs ( result - exact ), rum

        j = 12
        a = 0.0D+00
        big = 1.0D+00

        result = squank ( a, big, error, fifth, rum, no, f12 )

        exact = -0.6346651825D+00

        write ( *,
     &    '(2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g10.2,2x,g10.2)' )
     &    j, no, result, exact, abs ( result - exact ), rum

        j = 13
        a = 0.0D+00
        big = 1.0D+00

        result = squank ( a, big, error, fifth, rum, no, f13 )

        exact = 0.013492485650D+00

        write ( *,
     &    '(2x,i2,2x,i8,2x,g14.6,2x,g14.6,2x,g10.2,2x,g10.2)' )
     &    j, no, result, exact, abs ( result - exact ), rum

      end do

      return
      end
      function f01 ( x )

c*********************************************************************72
c
cc F01 evaluates test function 1.
c
      implicit none

      double precision f01
      double precision x

      f01 = sqrt ( x )

      return
      end
      function f02 ( x )

c*********************************************************************72
c
cc F02 evaluates test function 2.
c
      implicit none

      double precision f02
      double precision x

      f02 = 0.92D+00 * cosh ( x ) - cos ( x )

      return
      end
      function f03 ( x )

c*********************************************************************72
c
cc F03 evaluates test function 3.
c
      implicit none

      double precision f03
      double precision x

      f03 = 1.0D+00 / ( x**4 + x**2 + 0.9D+00 )

      return
      end
      function f04 ( x )

c*********************************************************************72
c
cc F04 evaluates test function 4.
c
      implicit none

      double precision f04
      double precision x

      f04 = sqrt ( x**3 )

      return
      end
      function f05 ( x )

c*********************************************************************72
c
cc F05 evaluates test function 5.
c
      implicit none

      double precision f05
      double precision x

      f05 = 1.0D+00 / ( 1.0D+00 + x**4 )

      return
      end
      function f06 ( x )

c*********************************************************************72
c
cc F06 evaluates test function 6.
c
      implicit none

      double precision f06
      double precision pi
      double precision x

      pi = 3.141592653589793D+00

      f06 = 1.0D+00 / ( 1.0D+00 + 0.5D+00 * sin ( 10.0D+00 * pi * x ) )

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

      double precision f07
      double precision x

      if ( x .eq. 0.0D+00 ) then
        f07 = 1.0D+00
      else
        f07 = x / ( exp ( x ) - 1.0D+00 )
      end if

      return
      end
      function f08 ( x )

c*********************************************************************72
c
cc F08 evaluates test function 8.
c
      implicit none

      double precision f08
      double precision pi
      double precision x

      pi = 3.141592653589793D+00

      f08 = sin ( 100.0D+00 * pi * x ) / ( pi * x )

      return
      end
      function f09 ( x )

c*********************************************************************72
c
cc F09 evaluates test function 9.
c
      implicit none

      double precision f09
      double precision pi
      double precision x

      pi = 3.141592653589793D+00

      f09 = 50.0D+00 / ( 2500.0D+00 * x**2 + 1.0D+00 ) / pi

      return
      end
      function f10 ( x )

c*********************************************************************72
c
cc F10 evaluates test function 10.
c
      implicit none

      double precision arg
      double precision f10
      double precision x

      arg =         cos (           x )
     &  + 3.0D+00 * sin (           x )
     &  + 2.0D+00 * cos ( 2.0D+00 * x )
     &  + 3.0D+00 * cos ( 3.0D+00 * x )
     &  + 3.0D+00 * sin ( 2.0D+00 * x )

      f10 = cos ( arg )

      return
      end
      function f11 ( x )

c*********************************************************************72
c
cc F11 evaluates test function 11.
c
c  Modified:
c
c    12 January 2006
c
      implicit none

      double precision f11
      double precision x

      if ( x .le. 0.0D+00 ) then
        f11 = 0.0D+00
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

      double precision f12
      double precision pi
      double precision x

      pi = 3.141592653589793D+00

      f12 = 4.0D+00 * pi**2 * x * sin ( 20.0D+00 * pi * x )
     &  * cos ( 2.0D+00 * pi * x )

      return
      end
      function f13 ( x )

c*********************************************************************72
c
cc F13 evaluates test function 13.
c
      implicit none

      double precision f13
      double precision x

      f13 = 1.0D+00 / ( 1.0D+00 + ( 230.0D+00 * x - 30.0D+00 )**2 )

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
