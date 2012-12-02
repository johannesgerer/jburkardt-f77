      program main

c*********************************************************************72
c
cc MAIN is the main program for DAY1_EX3.
c
c  Discussion:
c
c    DAY1_EX3 is exercise 3 for first day of the MPI workshop.
c
c    The instructions say:
c
c    Process 1 computes the squares of the first 200 integers.
c    It sends this data to process 3.
c
c    Process 3 should divide the integers between 20 and 119 by 53,
c    getting a real result, and passes this data back to process 1.
c
c    * I presume the first 200 integers are the numbers 0 through 199.
c
c    * The instructions literally mean that process 3 should look
c      at integers whose VALUES are between 20 and 119.  I doubt that
c      is what the instructor meant, but it's more interesting than
c      simply picking the entries with index between 20 and 119,
c      so that's what I'll do.
c
c    * It is also not completely clear whether only the selected data
c      should be sent back, or the entire array.  Again, it is more
c      interesting to send back only part of the data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 October 2005
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    William Gropp, Ewing Lusk, Anthony Skjellum,
c    Using MPI: Portable Parallel Programming with the
c    Message-Passing Interface,
c    Second Edition,
c    MIT Press, 1999,
c    ISBN: 0262571323,
c    LC: QA76.642.G76.
c
c    Marc Snir, Steve Otto, Steven Huss-Lederman, David Walker, 
c    Jack Dongarra,
c    MPI: The Complete Reference,
c    Volume I: The MPI Core,
c    Second Edition,
c    MIT Press, 1998,
c    ISBN: 0-262-69216-3,
c     LC: QA76.642.M65.
c

c
c  Fortran77 include file:
c
      include 'mpif.h'

      integer i_dim
      integer r_dim

      parameter ( i_dim = 200 )
      parameter ( r_dim = 200 )

      integer count
      integer count2
      integer dest
      integer i
      integer i_buffer(i_dim)
      integer ierr
      integer num_procs
      real r_buffer(r_dim)
      integer rank
      integer source
      integer status(MPI_Status_size)
      integer tag
c
c  Initialize MPI.
c
      call MPI_Init ( ierr )
c
c  Determine this process's rank.
c
      call MPI_Comm_rank ( MPI_COMM_WORLD, rank, ierr )
c
c  Have Process 0 say hello.
c
      if ( rank .eq. 0 ) then
        call timestamp ( )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DAY1_EX3:'
        write ( *, '(a)' ) '  FORTRAN77 version'
        write ( *, '(a)' ) '  MPI exercise #3 for day 1.'
        write ( *, '(a,i6)' ) 
     &    '  The number of processes available is ', num_procs

      end if
c
c  Get the number of processes.
c
      call MPI_Comm_size ( MPI_COMM_WORLD, num_procs, ierr )
c
c  If we don't have at least 4 processes, then bail out now.
c
      if ( num_procs .lt. 4 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' ) 'DAY1_EX3 - Process ', rank
        write ( *, '(a)' ) '  Not enough processes for this task!'
        write ( *, '(a)' ) '  Bailing out now!'
        call MPI_Finalize ( ierr )
        stop
      end if
c
c  Process 1 knows that it will generate 200 integers, and may receive no more
c  than 200 reals.
c
      if ( rank .eq. 1 ) then

        count = 200

        do i = 1, count
          i_buffer(i) = ( i - 1 )**2
        end do

        dest = 3
        tag = 1

        call MPI_Send ( i_buffer, count, MPI_INTEGER, dest, tag, 
     &    MPI_COMM_WORLD, ierr )

        write ( *, '(a,i1,a,i3,a,i1)' ) 'P:', rank, ' sent ', 
     &    count, ' integers to process ', dest

        source = 3
        tag = 2

        call MPI_Recv ( r_buffer, r_dim, MPI_REAL, source, tag, 
     &    MPI_COMM_WORLD, status, ierr )

        write ( *, '(a,i1,a,i1)' ) 'P:', rank, 
     &    ' received real values from process 3.'

        call MPI_Get_count ( status, MPI_REAL, count, ierr )

        write ( *, '(a,i1,a,i3,a)' ) 'P:', rank, 
     &    ' Number of real values received is ', count

        write ( *, '(a,i1,a,3f8.4)' ) 'P:', rank, 
     &    ' First 3 values = ', ( r_buffer(i), i = 1, 3 )
c
c  Process 3 receives the integer data from process 1, selects some 
c  of the data, does a real computation on it, and sends that part 
c  back to process 1.
c
      else if ( rank .eq. 3 ) then
 
        source = 1
        tag = 1

        call MPI_Recv ( i_buffer, i_dim, MPI_INTEGER, source, tag, 
     &      MPI_COMM_WORLD, status, ierr )

        write ( *, '(a)' ) ' '
        write ( *, '(a,i1,a)' ) 'P:', rank, 
     &    ' received integer values from process 1.'

        call MPI_Get_count ( status, MPI_INTEGER, count, ierr )

        write ( *, '(a,i1,a,i6)' ) 'P:', rank, 
     &      ' - Number of integers received is ', count
        write ( *, '(a,i1,a,3i8)' ) 'P:', rank, 
     &      ' First 3 values = ', ( i_buffer(i), i = 1, 3 )

        count2 = 0

        do i = 1, count
         
          if ( 20 .le. i_buffer(i) .and. 
     &      i_buffer(i) .le. 119 ) then

            count2 = count2 + 1
            r_buffer(count2) = real ( i_buffer(i) ) / 53.0E+00

            if ( count2 .le. 3 ) then
              write ( *, '(a,i1,a,i6,a,f8.4)' ) 'P:', rank, 
     &            ' Input integer ', i_buffer(i), ' becomes ', 
     &            r_buffer(count2)
            end if

          end if

        end do

        dest = 1
        tag = 2
      
        call MPI_Send ( r_buffer, count2, MPI_REAL, dest, tag, 
     &    MPI_COMM_WORLD, ierr )

        write ( *, '(a,i1,a,i3,a,i1)' ) 'P:', rank, ' sent ', count2, 
     &      ' reals to process ', dest

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a,i1,a)' ) 'P:', rank, ' - MPI has no work for me!'

      end if
c
c  Terminate MPI.
c
      call MPI_Finalize ( ierr )
c
c  Terminate.
c
      if ( rank .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DAY1_EX3:'
        write ( *, '(a)' ) '  Normal end of execution.'
        write ( *, '(a)' ) ' '
        call timestamp ( )
      end if

      stop
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
