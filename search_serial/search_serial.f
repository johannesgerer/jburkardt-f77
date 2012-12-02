      program main

c*********************************************************************72
c
cc MAIN is the main program for SEARCH_SERIAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 October 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer a
      integer b
      integer c
      integer f
      integer fj
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer j
      double precision wtime
      double precision wtime1
      double precision wtime2

      a = 1
      b = i4_huge
      c = 45

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SEARCH_SERIAL:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Search the integers from A to B'
      write ( *, '(a)' ) '  for a value J such that F(J) = C.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  A           = ', a
      write ( *, '(a,i12)' ) '  B           = ', b
      write ( *, '(a,i12)' ) '  C           = ', c

      call cpu_time ( wtime1 )

      call search ( a, b, c, j )

      call cpu_time ( wtime2 )
      wtime = wtime2 - wtime

      if ( j .eq. -1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  No solution was found.'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a,i12)' ) '  Found     J = ', j
        write ( *, '(a,i12)' ) '  Verify F(J) = ', f ( j )
      end if

      write ( *, '(a,g14.6)' ) '  Elapsed CPU time is ', wtime
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SEARCH_SERIAL:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine search ( a, b, c, j )

c*********************************************************************72
c
cc SEARCH searches integers in [A,B] for a J so that F(J) = C.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 October 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer A, B, the search range.
c
c    Input, integer C, the desired function value.
c
c    Output, integer J, the computed solution, or -1
c    if no solution was found.
c
      implicit none

      integer a
      integer b
      integer c
      integer f
      integer fi
      integer i
      integer j

      j = -1

      do i = a, b

        fi = f ( i )

        if ( fi .eq. c ) then
          j = i
          return
        end if

      end do

      return
      end
      function f ( i )

c*********************************************************************72
c
cc F is the function we are analyzing.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 October 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the argument.
c
c    Input, integer F, the value.
c
      implicit none

      integer f
      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer j
      integer k
      integer value

      value = i

      do j = 1, 5

        k = value / 127773

        value = 16807 * ( value - k * 127773 ) - k * 2836

        if ( value .le. 0 ) then
          value = value + i4_huge
        end if

      end do

      f = value

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
