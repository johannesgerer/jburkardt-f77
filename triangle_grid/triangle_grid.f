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
      subroutine triangle_grid ( n, t, tg )

c*********************************************************************72
c
cc TRIANGLE_GRID computes points on a triangular grid.
c
c  Discussion:
c
c    The grid is defined by specifying the coordinates of an enclosing
c    triangle T, and the number of subintervals each side of the triangle
c    should be divided into.
c
c    Choosing N = 10, for instance, breaks each side into 10 subintervals,
c    and produces a grid of ((10+1)*(10+2))/2 = 66 points.
c
c              X
c             9 X
c            8 9 X
c           7 8 9 X
c          6 7 8 9 X
c         5 6 7 8 9 X
c        4 5 6 7 8 9 X
c       3 4 5 6 7 8 9 X
c      2 3 4 5 6 7 8 9 X
c     1 2 3 4 5 6 7 8 9 X
c    0 1 2 3 4 5 6 7 8 9 X
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 September 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of subintervals.
c
c    Input, double precision T(2,3), the coordinates of the points
c    defining the triangle.
c
c    Output, double precision TG(2,((N+1)*(N+2))/2), the coordinates
c    of the points in the triangle.
c
      implicit none

      integer n

      integer i
      double precision ir
      integer j
      double precision jr
      integer k
      double precision kr
      integer l
      double precision nr
      integer p
      double precision t(2,3)
      double precision tg(2,((n+1)*(n+2))/2)

      p = 0
      nr = dble ( n )

      do i = 0, n
        ir = dble ( i )
        do j = 0, n - i
          jr = dble ( j )
          k = n - i - j
          kr = dble ( k )
          p = p + 1
          do l = 1, 2
            tg(l,p) = ( ir * t(l,1) + jr * t(l,2) + kr * t(l,3) ) / nr
          end do
        end do
      end do

      return
      end
