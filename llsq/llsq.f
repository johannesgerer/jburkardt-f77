      subroutine llsq ( n, x, y, a, b )

c*********************************************************************72
c
cc LLSQ solves a linear least squares problem matching a line to data.
c
c  Discussion:
c
c    A formula for a line of the form Y = A * X + B is sought, which
c    will minimize the root-mean-square error to N data points ( X(I), Y(I) );
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 March 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of data values.
c
c    Input, double precision X(N), Y(N), the coordinates of the data points.
c
c    Output, double precision A, B, the slope and Y-intercept of the 
c    least-squares approximant to the data.
c
      implicit none

      integer n

      double precision a
      double precision b
      double precision bot
      integer i
      double precision top
      double precision x(n)
      double precision xbar
      double precision y(n)
      double precision ybar
c
c  Special case.
c
      if ( n .eq. 1 ) then
        a = 0.0D+00
        b = y(1)
        return
      end if
c
c  Average X and Y.
c
      xbar = 0.0D+00
      ybar = 0.0D+00
      do i = 1, n
        xbar = xbar + x(i)
        ybar = ybar + y(i)
      end do
      xbar = xbar / dble ( n )
      ybar = ybar / dble ( n )
c
c  Compute Beta.
c
      top = 0.0D+00
      bot = 0.0D+00
      do i = 1, n
        top = top + ( x(i) - xbar ) * ( y(i) - ybar )
        bot = bot + ( x(i) - xbar ) * ( x(i) - xbar )
      end do

      a = top / bot

      b = ybar - a * xbar

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
