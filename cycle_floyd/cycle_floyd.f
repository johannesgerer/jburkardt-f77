      subroutine cycle_floyd ( f, x0, lam, mu )

c*********************************************************************72
c
cc CYCLE_FLOYD finds a cycle in an iterated mapping using Floyd's method.
c
c  Discussion:
c
c    Suppose we a repeatedly apply a function f(), starting with the argument
c    x0, then f(x0), f(f(x0)) and so on.  Suppose that the range of f is finite.
c    Then eventually the iteration must reach a cycle.  Once the cycle is reached,
c    succeeding values stay within that cycle.
c
c    Starting at x0, there is a "nearest element" of the cycle, which is
c    reached after MU applications of f.
c
c    Once the cycle is entered, the cycle has a length LAM, which is the number
c    of steps required to first return to a given value.
c
c    This function uses Floyd's method to determine the values of MU and LAM,
c    given F and X0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Knuth,
c    The Art of Computer Programming,
c    Volume 2, Seminumerical Algorithms,
c    Third Edition,
c    Addison Wesley, 1997,
c    ISBN: 0201896842,
c    LC: QA76.6.K64.
c
c  Parameters:
c
c    Input, external integer F(), the name of the function 
c    to be analyzed.
c
c    Input, integer X0, the starting point.
c
c    Output, integer LAM, the length of the cycle.
c
c    Output, integer MU, the index in the sequence starting
c    at X0, of the first appearance of an element of the cycle.
c
      implicit none

      external f
      integer f
      integer hare
      integer lam
      integer mu
      integer tortoise
      integer x0

      tortoise = f ( x0 )
      hare = f ( tortoise )

10    continue

      if ( tortoise .ne. hare ) then
        tortoise = f ( tortoise )
        hare = f ( f ( hare ) )
        go to 10
      end if

      mu = 0
      tortoise = x0

20    continue

      if ( tortoise .ne. hare ) then
        tortoise = f ( tortoise )
        hare = f ( hare )
        mu = mu + 1
        go to 20
      end if

      lam = 1
      hare = f ( tortoise )

30    continue

      if ( tortoise .ne. hare ) then
        hare = f ( hare )
        lam = lam + 1
        go to 30
      end if

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
