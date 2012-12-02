      program main

c*********************************************************************72
c
cc MAIN is the main program for EXPONENT_FORMAT_OVERFLOW.
c
c  Discussion:
c
c    EXPONENT_FORMAT_OVERFLOW examines the format used to print real values
c    in exponential format, in cases where the exponent has large magnitude.
c
c    It has been observed that, at least for some compilers, the case in
c    which the exponent has three digits is handled in a very bad way
c    that is misleading and liable to result in errors, particularly if
c    one program writes out the data and another program is to read it
c    back in.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EXPONENT_FORMAT_OVERFLOW:'
      write ( *, '(a)' ) '  FORTRAN77 version.'
     
      call test01 ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EXPONENT_FORMAT_OVERFLOW:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 prints some large values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision x1
      double precision x2
      double precision x3
      double precision x4
      double precision x5
      double precision x6
      double precision x7
      double precision x8
      double precision x9
      double precision x10
      double precision x11
      double precision x12

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) 
     &  '  Real numbers can have exponents greater than 99.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  It is unsettling to discover that there is a'
      write ( *, '(a)' ) 
     &  '  common FLAW in certain output formats, which means'
      write ( *, '(a)' ) 
     &  '  that the printing of real numbers with exponents'
      write ( *, '(a)' ) 
     &  '  of magnitude more than 99 is handled poorly.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  The "E" marker, used to indicate scientific notation,'
      write ( *, '(a)' ) 
     &  '  is simply suppressed, as though to make room for'
      write ( *, '(a)' ) '  the leading digit of the exponent.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  But this inconsistency can be deadly.  In particular,'
      write ( *, '(a)' ) 
     &  '  if you write out such data, it CANNOT be read back in'
      write ( *, '(a)' ) 
     &  '  properly.  (THAT example is easy to write, but I will'
      write ( *, '(a)' ) '  do that later.)'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Define some big numbers:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X1 = 1.79769313486231571D+308'
      write ( *, '(a)' ) '  X2 = 0.99D+00 * X1'
      write ( *, '(a)' ) '  X3 = X1 / 100.0D+00'
      write ( *, '(a)' ) '  X4 = 1.0D+101.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Define some small numbers:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X5 = 2.22507385850720138D-308'
      write ( *, '(a)' ) '  X6 = 1.01D+00 * X5'
      write ( *, '(a)' ) '  X7 = X5 * 100.0D+00'
      write ( *, '(a)' ) '  X8 = 1.0D-101.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Define some comparison numbers:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X9  = 1.0D+98'
      write ( *, '(a)' ) '  X10 = 1.0D-98'
      write ( *, '(a)' ) '  X11 = 123456789.0'
      write ( *, '(a)' ) '  X12 = 0.123456789'

      x1 = 1.79769313486231571D+308
      x2 = 0.99D+00 * x1
      x3 = x1 / 10.0D+00
      x4 = 1.0D+101

      x5 = 2.22507385850720138D-308
      x6 = 1.01D+00 * x5
      x7 = 100.0D+00 * x5
      x8 = 1.0D-101

      x9 = 1.0D+98
      x10 = 1.0D-98
      x11 = 123456789.0D+00
      x12 = 0.123456789D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Print with a WRITE(*,*) format:'
      write ( *, '(a)' ) '  This seems to work OK:'
      write ( *, '(a)' ) ' '
      write ( *, * ) '  X1 = ', x1
      write ( *, * ) '  X2 = ', x2
      write ( *, * ) '  X3 = ', x3
      write ( *, * ) '  X4 = ', x4
      write ( *, * ) '  X5 = ', x5
      write ( *, * ) '  X6 = ', x6
      write ( *, * ) '  X7 = ', x7
      write ( *, * ) '  X8 = ', x8
      write ( *, * ) '  X9 = ', x9
      write ( *, * ) '  X10 = ', x10
      write ( *, * ) '  X11 = ', x11
      write ( *, * ) '  X12 = ', x12
      
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Print with a WRITE(*,''(G24.10)'') format:'
      write ( *, '(a)' ) '  Notice the disaster:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g24.10)' ) '  X1 = ', x1
      write ( *, '(a,g24.10)' ) '  X2 = ', x2
      write ( *, '(a,g24.10)' ) '  X3 = ', x3
      write ( *, '(a,g24.10)' ) '  X4 = ', x4
      write ( *, '(a,g24.10)' ) '  X5 = ', x5
      write ( *, '(a,g24.10)' ) '  X6 = ', x6
      write ( *, '(a,g24.10)' ) '  X7 = ', x7
      write ( *, '(a,g24.10)' ) '  X8 = ', x8
      write ( *, '(a,g24.10)' ) '  X9 = ', x9
      write ( *, '(a,g24.10)' ) '  X10 = ', x10
      write ( *, '(a,g24.10)' ) '  X11 = ', x11
      write ( *, '(a,g24.10)' ) '  X12 = ', x12

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Print with a WRITE(*,''(D24.10)'') format:'
      write ( *, '(a)' ) '  Same disaster:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,d24.10)' ) '  X1 = ', x1
      write ( *, '(a,d24.10)' ) '  X2 = ', x2
      write ( *, '(a,d24.10)' ) '  X3 = ', x3
      write ( *, '(a,d24.10)' ) '  X4 = ', x4
      write ( *, '(a,d24.10)' ) '  X5 = ', x5
      write ( *, '(a,d24.10)' ) '  X6 = ', x6
      write ( *, '(a,d24.10)' ) '  X7 = ', x7
      write ( *, '(a,d24.10)' ) '  X8 = ', x8
      write ( *, '(a,d24.10)' ) '  X9 = ', x9
      write ( *, '(a,d24.10)' ) '  X10 = ', x10
      write ( *, '(a,d24.10)' ) '  X11 = ', x11
      write ( *, '(a,d24.10)' ) '  X12 = ', x12

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Print with a WRITE(*,''(E24.10)'') format:'
      write ( *, '(a)' ) '  Same disaster:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,e24.10)' ) '  X1 = ', x1
      write ( *, '(a,e24.10)' ) '  X2 = ', x2
      write ( *, '(a,e24.10)' ) '  X3 = ', x3
      write ( *, '(a,e24.10)' ) '  X4 = ', x4
      write ( *, '(a,e24.10)' ) '  X5 = ', x5
      write ( *, '(a,e24.10)' ) '  X6 = ', x6
      write ( *, '(a,e24.10)' ) '  X7 = ', x7
      write ( *, '(a,e24.10)' ) '  X8 = ', x8
      write ( *, '(a,e24.10)' ) '  X9 = ', x9
      write ( *, '(a,e24.10)' ) '  X10 = ', x10
      write ( *, '(a,e24.10)' ) '  X11 = ', x11
      write ( *, '(a,e24.10)' ) '  X12 = ', x12

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
 
