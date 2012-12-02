      program main

c*********************************************************************72
c
cc MAIN is the main program for GFORTRAN_QUADMATH.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 April 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GFORTRAN_QUADMATH'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the GFORTRAN quadmath facility.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GFORTRAN_QUADMATH'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 uses REAL ( KIND = 4 ) arithmetic.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 April 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer divs
      real ( kind = 4 ) x
      real ( kind = 4 ) x_old
      real ( kind = 4 ) y
      real ( kind = 4 ) z

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Using REAL ( KIND = 4 ) arithmetic:'
      write ( *, '(a)' ) 
     &  '  Compute smallest 1/2^DIV that can be added to 1.'

      x = 1.0
      z = 1.0
      divs = 0

10    continue
        x_old = x
        x = x / 2.0
        y = 1.0 + x
        if ( y .le. z ) then
          go to 20
        end if
        divs = divs + 1
      go to 10

20    continue

      write ( *, '(a,i4)' ) '  Number of divisions DIV = ', divs
      write ( *, '(a,g14.6)' ) '  1/2^DIV =         ', x_old
      write ( *, '(a,g14.6)' ) '  Machine epsilon = ', epsilon ( x_old )

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 uses REAL ( KIND = 8 ) arithmetic.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 April 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer divs
      real ( kind = 8 ) x
      real ( kind = 8 ) x_old
      real ( kind = 8 ) y
      real ( kind = 8 ) z

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Using REAL ( KIND = 8 ) arithmetic:'
      write ( *, '(a)' ) 
     &  '  Compute smallest 1/2^DIV that can be added to 1.'

      x = 1.0
      z = 1.0
      divs = 0

10    continue
        x_old = x
        x = x / 2.0
        y = 1.0 + x
        if ( y .le. z ) then
          go to 20
        end if
        divs = divs + 1
      go to 10

20    continue

      write ( *, '(a,i4)' ) '  Number of divisions DIV = ', divs
      write ( *, '(a,g14.6)' ) '  1/2^DIV =         ', x_old
      write ( *, '(a,g14.6)' ) '  Machine epsilon = ', epsilon ( x_old )

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 uses REAL ( KIND = 10 ) arithmetic.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 April 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer divs
      real ( kind = 10 ) x
      real ( kind = 10 ) x_old
      real ( kind = 10 ) y
      real ( kind = 10 ) z

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  Using REAL ( KIND = 10 ) arithmetic:'
      write ( *, '(a)' ) 
     &  '  Compute smallest 1/2^DIV that can be added to 1.'

      x = 1.0
      z = 1.0
      divs = 0

10    continue
        x_old = x
        x = x / 2.0
        y = 1.0 + x
        if ( y .le. z ) then
          go to 20
        end if
        divs = divs + 1
      go to 10

20    continue

      write ( *, '(a,i4)' ) '  Number of divisions DIV = ', divs
      write ( *, '(a,g14.6)' ) '  1/2^DIV =         ', x_old
      write ( *, '(a,g14.6)' ) '  Machine epsilon = ', epsilon ( x_old )

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
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
