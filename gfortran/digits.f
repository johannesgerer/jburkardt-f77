      program main

c*********************************************************************72
c
cc MAIN is the main program for DIGITS.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 November 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DIGITS'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  How many digits are appropriate for data?'

      call test01 ( )
      call test02 ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DIGITS:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 looks at data stored as REAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 November 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer digit_max
      parameter ( digit_max = 20 )

      real c
      real d
      integer digit
      real pi(digit_max)
      real s

      save pi

      data pi /
     &  3.0E+00, 
     &  3.1E+00, 
     &  3.14E+00, 
     &  3.141E+00, 
     &  3.1415E+00, 
     &  3.14159E+00, 
     &  3.141592E+00, 
     &  3.1415926E+00, 
     &  3.14159265E+00, 
     &  3.141592653E+00, 
     &  3.1415926535E+00, 
     &  3.14159265358E+00, 
     &  3.141592653589E+00, 
     &  3.1415926535897E+00, 
     &  3.14159265358979E+00, 
     &  3.141592653589793E+00, 
     &  3.1415926535897932E+00, 
     &  3.14159265358979323E+00, 
     &  3.141592653589793238E+00, 
     &  3.1415926535897932384E+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  REAL arithmetic'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Store PI using 1 through 20 digits.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Digits     sin(PI(I))   1+cos(PI(I))     (PI(I)-PI(20))'
      write ( *, '(a)' ) ' '

      do digit = 1, digit_max
        s = sin ( pi(digit) )
        c = 1.0E+00 + cos ( pi(digit) )
        d = pi(digit) - pi(digit_max)
        write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &    digit, s, c, d
      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 looks at data stored as DOUBLE PRECISION.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 November 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer digit_max
      parameter ( digit_max = 20 )

      double precision c
      double precision d
      integer digit
      double precision pi(digit_max)
      double precision s

      save pi

      data pi /
     &  3.0D+00, 
     &  3.1D+00, 
     &  3.14D+00, 
     &  3.141D+00, 
     &  3.1415D+00, 
     &  3.14159D+00, 
     &  3.141592D+00, 
     &  3.1415926D+00, 
     &  3.14159265D+00, 
     &  3.141592653D+00, 
     &  3.1415926535D+00, 
     &  3.14159265358D+00, 
     &  3.141592653589D+00, 
     &  3.1415926535897D+00, 
     &  3.14159265358979D+00, 
     &  3.141592653589793D+00, 
     &  3.1415926535897932D+00, 
     &  3.14159265358979323D+00, 
     &  3.141592653589793238D+00, 
     &  3.1415926535897932384D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  DOUBLE PRECISION arithmetic'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Store PI using 1 through 20 digits.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Digits     sin(PI(I))   1+cos(PI(I))     (PI(I)-PI(20))'
      write ( *, '(a)' ) ' '

      do digit = 1, digit_max
        s = sin ( pi(digit) )
        c = 1.0D+00 + cos ( pi(digit) )
        d = pi(digit) - pi(digit_max)
        write ( *, '(2x,i8,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &    digit, s, c, d
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
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 January 2007
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

      character * ( 8 ) ampm
      integer d
      character * ( 8 ) date
      integer h
      integer m
      integer mm
      character * ( 9 ) month(12)
      integer n
      integer s
      character * ( 10 ) time
      integer y

      save month

      data month /
     &  'January  ', 'February ', 'March    ', 'April    ', 
     &  'May      ', 'June     ', 'July     ', 'August   ', 
     &  'September', 'October  ', 'November ', 'December ' /

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) y, m, d
      read ( time, '(i2,i2,i2,1x,i3)' ) h, n, s, mm

      if ( h .lt. 12 ) then
        ampm = 'AM'
      else if ( h .eq. 12 ) then
        if ( n .eq. 0 .and. s .eq. 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h .lt. 12 ) then
          ampm = 'PM'
        else if ( h .eq. 12 ) then
          if ( n .eq. 0 .and. s .eq. 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if

      write ( *, 
     &  '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) 
     &  d, month(m), y, h, ':', n, ':', s, '.', mm, ampm

      return
      end
