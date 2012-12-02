      subroutine cycle_brent ( f, x0, lam, mu )

c*********************************************************************72
c
cc CYCLE_BRENT finds a cycle in an iterated mapping using Brent's method.
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
c    This function uses Brent's method to determine the values of MU and LAM,
c    given F and X0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Richard Brent,
c    An improved Monte Carlo factorization algorithm,
c    BIT,
c    Volume 20, Number 2, 1980, pages 176-184.
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
      integer i
      integer lam
      integer mu
      integer power
      integer tortoise
      integer x0

      power = 1
      lam = 1
      tortoise = x0
      hare = f ( x0 )

10    continue

      if ( tortoise .ne. hare ) then
        if ( power .eq. lam ) then
          tortoise = hare
          power = power * 2
          lam = 0
        end if
        hare = f ( hare )
        lam = lam + 1
        go to 10
      end if
     
      mu = 0
      tortoise = x0
      hare = x0

      do i = 0, lam - 1
        hare = f ( hare )
      end do

20    continue

      if ( tortoise .ne. hare ) then
        tortoise = f ( tortoise )
        hare = f ( hare )
        mu = mu + 1
        go to 20
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
