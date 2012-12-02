      program main

c*********************************************************************72
c
cc MAIN is the main program for SEARCH.
c
c  Discussion:
c
c    SEARCH demonstrates the use of MPI routines to carry out a search.
c
c    An array of given size is to be searched for occurrences of a
c    specific value.
c
c    The search is done in parallel.  A master process generates the
c    array and the target value, then distributes the information among
c    a set of worker processes, and waits for them to communicate back
c    the (global) index values at which occurrences of the target value
c    were found.
c
c    An interesting feature of this program is the use of allocatable
c    arrays, which allows the master program to set aside just enough
c    memory for the whole array, and for each worker program to set aside
c    just enough memory for its own part of the array.
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
      include 'mpif.h'

      integer n_max
      parameter ( n_max = 600 )

      integer a(n_max)
      integer dest
      real factor
      integer global
      integer i
      integer ierr
      integer master
      integer my_id
      integer n
      integer npart
      integer num_procs
      real r(n_max)
      real r4_uniform_01
      integer seed
      integer source
      integer start
      integer status(MPI_STATUS_SIZE)
      integer tag
      integer tag_target
      integer tag_size
      integer tag_data
      integer tag_found
      integer tag_done
      integer target
      integer workers_done
      integer x

      master = 0
      seed = 123456789
      tag_target = 1
      tag_size = 2
      tag_data = 3
      tag_found = 4
      tag_done = 5
c
c  Initialize MPI.
c
      call MPI_Init ( ierr )
c
c  Get this processes's rank.
c
      call MPI_Comm_rank ( MPI_COMM_WORLD, my_id, ierr )
c
c  Find out how many processes are available.
c
      call MPI_Comm_size ( MPI_COMM_WORLD, num_procs, ierr )

      if ( my_id .eq. master ) then
        call timestamp ( )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SEARCH - Master process:'
        write ( *, '(a)' ) '  FORTRAN77 version'
        write ( *, '(a)' ) '  An MPI program to search an array.'
        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' ) '  The number of processes is ', num_procs
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6,a)' ) 'Process ', my_id, ' is active.'
c
c  Have the master process generate the target and data.  In a more 
c  realistic application, the data might be in a file which the master 
c  process would read.  Here, the master process decides
c
      if ( my_id .eq. master ) then
c
c  Pick the number of data items per process, and set the total.
c
        factor = r4_uniform_01 ( seed )
        npart = 50 + int ( factor * 100.0E+00 )
        n = npart * num_procs

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SEARCH - Master process:'
        write ( *, '(a,i6)' ) 
     &    '  The number of data items per process is ', npart
        write ( *, '(a,i6)' ) 
     &    '  The total number of data items is       ', n
c
c  Now fill A with values, and pick a value for the target.
c
        do i = 1, n
          r(i) = r4_uniform_01 ( seed )
        end do

        factor = real ( n ) / 10
        do i = 1, n
          a(i) = int ( factor * r(i) )
        end do

        target = a(n/2)

        write ( *, '(a,i6)' ) '  The target value is ', target
c
c  The worker processes need to have the target value, the number of data items,
c  and their individual chunk of the data vector.
c
        do i = 1, num_procs-1

          dest = i
          tag = tag_target

          call MPI_Send ( target, 1, MPI_INTEGER, dest, tag, 
     &      MPI_COMM_WORLD, ierr )

          tag = tag_size

          call MPI_Send ( npart, 1, MPI_INTEGER, dest, tag, 
     &      MPI_COMM_WORLD, ierr )

          start = ( i - 1 ) * npart + 1
          tag = tag_data

          call MPI_Send ( a(start), npart, MPI_INTEGER, dest, tag, 
     &      MPI_COMM_WORLD, ierr )

        end do
c
c  Now the master process simply waits for each worker process to report that 
c  it is done.
c
        workers_done = 0

10      continue

          if ( num_procs-1 .le. workers_done ) then
            go to 20
          end if

          call MPI_Recv ( x, 1, MPI_INTEGER, MPI_ANY_SOURCE, 
     &      MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr )

          source = status(MPI_SOURCE)
          tag = status(MPI_TAG)
        
          if ( tag .eq. tag_done ) then

            workers_done = workers_done + 1

          else if ( tag .eq. tag_found ) then

            write ( *, '(a,i2,2i6)' ) 'P', source, x, a(x)

          else

            write ( *, '(a,i6)' ) 
     &        '  Master process received unknown tag = ', tag

          end if

        go to 10

20      continue
c
c  Each worker process expects to receive the target value, the number of data
c  items, and the data vector.
c
      else 

        source = master
        tag = tag_target

        call MPI_Recv ( target, 1, MPI_INTEGER, source, tag, 
     &    MPI_COMM_WORLD, status, ierr )
 
        source = master
        tag = tag_size

        call MPI_Recv ( npart, 1, MPI_INTEGER, source, tag, 
     &    MPI_COMM_WORLD, status, ierr )

        source = master
        tag = tag_data

        call MPI_Recv ( a, npart, MPI_INTEGER, source, tag, 
     &    MPI_COMM_WORLD, status, ierr )
c
c  The worker simply checks each entry to see if it is equal to the target
c  value.
c
        do i = 1, npart

          if ( a(i) .eq. target ) then

            global = ( my_id - 1 ) * npart + i
            dest = master
            tag = tag_found

            call MPI_Send ( global, 1, MPI_INTEGER, dest, tag,
     &       MPI_COMM_WORLD, ierr )

          end if

        end do  
c
c  When the worker is finished with the loop, it sends a dummy data value with
c  the tag "TAG_DONE" indicating that it is done.
c
        dest = master
        tag = tag_done

        call MPI_Send ( target, 1, MPI_INTEGER, dest, tag, 
     &    MPI_COMM_WORLD, ierr )
         
      end if 
c
c  Terminate MPI.
c
      call MPI_Finalize ( ierr )
c
c  Terminate.
c
      if ( my_id .eq. master ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SEARCH - Master process:'
        write ( *, '(a)' ) '  Normal end of execution.'
        write ( *, '(a)' ) ' '
        call timestamp ( )
      end if
 
      stop
      end
      function r4_uniform_01 ( seed )

c*********************************************************************72
c
cc R4_UNIFORM_01 returns a unit real pseudorandom number.
c
c  Discussion:
c
c    The pseudorandom number should be uniformly distributed
c    between 0 and 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 January 2005
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, real R4_UNIFORM_01, a number strictly between 0 and 1.
c
      implicit none

      integer k
      integer seed
      real r4_uniform_01

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + 2147483647
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r4_uniform_01 = real ( dble ( seed ) * 4.656612875D-10 )

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
