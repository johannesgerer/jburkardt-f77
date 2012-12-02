      program main

c*********************************************************************72
c
cc MAIN is the main program for DOUBLE_COMPLEX.
c
c  Discussion:
c
c    DOUBLE_COMPLEX runs the double precision complex demonstration examples.
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
      write ( *, '(a)' ) 'DOUBLE_COMPLEX'
      write ( *, '(a)' ) '  FORTRAN77 version.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Demonstrate double precision complex arithmetic.'

      call test01 ( )
      call test02 ( )
      call test03 ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DOUBLE_COMPLEX'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 shows what you get by default.
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

      complex a
      complex b
      complex c

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Use complex values of the default type.'

      a = cmplx ( 1.0E+00, 1.0E+00 )
      b = sqrt ( a )
      c = a - b * b

      write ( *, '(a)' ) ' '
      write ( *, * ) '  A = ', a
      write ( *, * ) '  B = sqrt ( A ) = ', b
      write ( *, * ) '  C = A - B * B = ', c

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 shows what you get with KIND = 8.
c
c  Discussion:
c
c    COMPLEX ( KIND = 8 ) makes sense under the FORTRAN90 standard, but
c    not the FORTRAN77 standard.  Nonetheless, you may be able to use this
c    data type if your FORTRAN77 compiler has been "modernized".
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

      complex ( kind = 8 ) a
      complex ( kind = 8 ) b
      complex ( kind = 8 ) c

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Use complex values of KIND = 8.'
      write ( *, '(a)' ) ' '

      a = cmplx ( 1.0D+00, 1.0D+00, kind = 8 )
      b = sqrt ( a )
      c = a - b * b

      write ( *, '(a)' ) ' '
      write ( *, * ) '  A = ', a
      write ( *, * ) '  B = sqrt ( A ) = ', b
      write ( *, * ) '  C = A - B * B = ', c

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 shows what happens if you ask for "DOUBLE COMPLEX".
c
c  Discussion:
c
c    "Double complex" is a nonstandard data type that was implemented
c    on some machines for FORTRAN77.
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

      double complex a
      double complex b
      double complex c

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  Use "DOUBLE COMPLEX" values.'

      a = cmplx ( 1.0D+00, 1.0D+00 )
      b = sqrt ( a )
      c = a - b * b

      write ( *, '(a)' ) ' '
      write ( *, * ) '  A = ', a
      write ( *, * ) '  B = sqrt ( A ) = ', b
      write ( *, * ) '  C = A - B * B = ', c

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
