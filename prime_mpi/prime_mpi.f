      program main

c*********************************************************************72
c
cc MAIN is the main program for PRIME_MPI.
c
c  Discussion:
c
c    This program calls a version of PRIME_NUMBER that includes
c    MPI calls for parallel processing.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 August 2009
c
c  Author:
c
c    John Burkardt
c
      include 'mpif.h'

      integer i
      integer id
      integer ierr
      integer n
      integer n_factor
      integer n_hi
      integer n_lo
      integer p
      integer primes
      integer primes_part
      double precision wtime

      n_lo = 1
      n_hi = 262144
      n_factor = 2
c
c  Initialize MPI.
c
      call MPI_Init ( ierr )
c
c  Get this process's ID.
c
      call MPI_Comm_rank ( MPI_COMM_WORLD, id, ierr )
c
c  Find out how many processes are available.
c
      call MPI_Comm_size ( MPI_COMM_WORLD, p, ierr )

      if ( id .eq. 0 ) then
        call timestamp ( )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PRIME_MPI'
        write ( *, '(a)' ) '  FORTRAN77/MPI version'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    '  An MPI example program to count the number of primes.'
        write ( *, '(a,i8)' ) '  The number of processes is ', p
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '         N        Pi          Time'
        write ( *, '(a)' ) ' '
      end if

      n = n_lo

10    continue

      if ( n .le. n_hi ) then

        if ( id .eq. 0 ) then
          wtime = MPI_Wtime ( )
        end if

        call MPI_Bcast ( n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, 
     &    ierr )

        call prime_number ( n, id, p, primes_part )

        call MPI_Reduce ( primes_part, primes, 1, MPI_INTEGER, MPI_SUM, 
     &    0, MPI_COMM_WORLD, ierr )

        if ( id .eq. 0 ) then

          wtime = MPI_Wtime ( ) - wtime

          write ( *, '(2x,i8,2x,i8,g14.6)' ) n, primes, wtime

        end if

        n = n * n_factor

        go to 10

      end if
c
c  Finish up.
c
      call MPI_Finalize ( ierr )

      if ( id .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PRIME_MPI - Master process:'
        write ( *, '(a)' ) '  Normal end of execution.'
        write ( *, '(a)' ) ' '
        call timestamp ( )
      end if

      stop
      end
      subroutine prime_number ( n, id, p, total )

c*********************************************************************72
c
cc PRIME_NUMBER returns a part of the number of primes between 1 and N.
c
c  Discussion:
c
c    In order to divide the work up evenly among P processors, processor
c    ID starts at 2+ID and skips by P.
c
c    A naive algorithm is used.
c
c    Mathematica can return the number of primes less than or equal to N
c    by the command PrimePi[N].
c
c                N  PRIME_NUMBER
c
c                1           0
c               10           4
c              100          25
c            1,000         168
c           10,000       1,229
c          100,000       9,592
c        1,000,000      78,498
c       10,000,000     664,579
c      100,000,000   5,761,455
c    1,000,000,000  50,847,534
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the maximum number to check.
c
c    Input, integer ID, the ID of this process,
c    between 0 and P-1.
c
c    Input, integer P, the number of processes.
c
c    Output, integer TOTAL, the number of prime numbers up to N,
c    starting at 2+ID and skipping by P.
c
      implicit none

      integer i
      integer id
      integer j
      integer n
      integer p
      integer prime
      integer total

      total = 0

      do i = 2+id, n, p

        prime = 1

        do j = 2, i - 1
          if ( mod ( i, j ) .eq. 0 ) then
            prime = 0
            go to 10
          end if
        end do

10      continue

        total = total + prime

      end do

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
