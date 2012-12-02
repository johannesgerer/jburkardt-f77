      program main

c*********************************************************************72
c
cc MAIN is the main program for SEARCH_MPI.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 November 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      include 'mpif.h'

      integer a
      integer b
      integer c
      integer f
      integer fj
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer id
      integer ierr
      integer j
      integer p
      double precision wtime

      call MPI_Init ( ierr )

      call MPI_Comm_Rank ( MPI_COMM_WORLD, id, ierr )

      call MPI_Comm_Size ( MPI_COMM_WORLD, p, ierr )

      a = 1
      b = i4_huge
      c = 45

      if ( id .eq. 0 ) then

        call timestamp ( )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SEARCH_MPI:'
        write ( *, '(a)' ) '  FORTRAN77/MPI version'
        write ( *, '(a)' ) '  Search the integers from A to B'
        write ( *, '(a)' ) '  for a value J such that F(J) = C.'
        write ( *, '(a)' ) ' '
        write ( *, '(a,i12)' ) '  A           = ', a
        write ( *, '(a,i12)' ) '  B           = ', b
        write ( *, '(a,i12)' ) '  C           = ', c
      end if

      wtime = MPI_Wtime ( )

      call search ( a, b, c, id, p, j )

      wtime = MPI_Wtime ( )- wtime
c
c  Any process that finds a solution should report it.
c
      if ( j .ne. -1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i1,a,i12)' ) 
     &    '  Process ', id, ' found     J = ', j
        write ( *, '(a,i12)' ) '  Verify F(J) = ', f ( j )
      end if

      if ( id .eq. 0 ) then
        write ( *, '(a,g14.6)' ) '  Elapsed CPU time is ', wtime
      end if
c
c  Terminate.
c
      if ( id .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SEARCH_MPI:'
        write ( *, '(a)' ) '  Normal end of execution.'
        write ( *, '(a)' ) ' '
        call timestamp ( )
      end if

      call MPI_Finalize ( ierr )

      stop
      end
      subroutine search ( a, b, c, id, p, j )

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
c    01 November 2012
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
c    Input, integer ID, the identifier of this process.
c    0 <= ID < P.
c
c    Input, integer P, the number of processes.
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
      integer id
      integer j
      integer p

      j = -1
c
c  We may have to be more careful here, since, when I + P exceeds
c  the maximum integer, it may actually look to be negative.
c
      do i = a + id, b, p

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
c    01 November 2012
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

        if ( value .lt. 0 ) then
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

