      program main

c*********************************************************************72
c
cc MAIN is the main program for COMMUNICATOR_MPI.
c
c  Discussion:
c
c    This program demonstrates how an MPI program can start with the
c    default communicator MPI_COMM_WORLD, and create new communicators
c    referencing a subset of the total number of processes.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 January 2012
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
      include 'mpif.h'

      integer even_comm_id
      integer even_group_id
      integer even_id
      integer even_id_sum
      integer even_p
      integer even_rank(10)
      integer i
      integer id
      integer ierr
      integer j
      integer odd_comm_id
      integer odd_group_id
      integer odd_id
      integer odd_id_sum
      integer odd_p
      integer odd_rank(10)
      integer p
      integer world_group_id
c
c  Initialize MPI.
c
      call MPI_Init ( ierr )
c
c  Get the number of processes.
c
      call MPI_Comm_size ( MPI_COMM_WORLD, p, ierr )
c
c  Get the individual process ID.
c
      call MPI_Comm_rank ( MPI_COMM_WORLD, id, ierr )
c
c  Process 0 prints an introductory message.
c
      if ( id .eq. 0 ) then
        call timestamp ( )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'COMMUNICATOR_MPI - Master process:'
        write ( *, '(a)' ) '  FORTRAN77/MPI version'
        write ( *, '(a)' ) '  An MPI example program.'
        write ( *, '(a)' ) ' '
        write ( *, '(a,i4)' ) '  The number of processes is ', p
        write ( *, '(a)' ) ' '
      end if
c
c  Every process prints a hello.
c
      write ( *, '(a,i4,a)' ) '  Process ', id, ' says "Hello, world!".'
c
c  Get a group identifier for MPI_COMM_WORLD.
c


      call MPI_Comm_group ( MPI_COMM_WORLD, world_group_id, ierr )
c
c  List the even processes, and create their group.
c
      even_p = ( p + 1 ) / 2
      j = 0
      do i = 0, p - 1, 2
        j = j + 1
        even_rank(j) = i
      end do
      call MPI_Group_incl ( world_group_id, even_p, even_rank, 
     &  even_group_id, ierr )

      call MPI_Comm_create ( MPI_COMM_WORLD, even_group_id, 
     &  even_comm_id, ierr )
c
c  List the odd processes, and create their group.
c
      odd_p = p / 2
      j = 0
      do i = 1, p - 1, 2
        j = j + 1
        odd_rank(j) = i
      end do
      call MPI_Group_incl ( world_group_id, odd_p, odd_rank, 
     &  odd_group_id, ierr )

      call MPI_Comm_create ( MPI_COMM_WORLD, odd_group_id, 
     &  odd_comm_id, ierr )
c
c  Try to get ID of each process in both groups.  
c  If a process is not in a communicator, set its ID to -1.
c
      if ( mod ( id, 2 ) .eq. 0 ) then
        call MPI_Comm_rank ( even_comm_id, even_id, ierr )
        odd_id = -1
      else
        call MPI_Comm_rank ( odd_comm_id,  odd_id, ierr )
        even_id = -1
      end if
c
c  Use MPI_Reduce to sum the global ID of each process in the even group.
c  Assuming 4 processes: EVEN_SUM = 0 + 2 = 2
c
      if ( even_id .ne. -1 ) then
        call MPI_Reduce ( id, even_id_sum, 1, MPI_INTEGER, MPI_SUM, 0, 
     &    even_comm_id, ierr )
      end if

      if ( even_id .eq. 0 ) then
        write ( *, '(a,i4)' ) 
     &    '  Number of processes in even communicator = ', even_p
        write ( *, '(a,i4)' ) 
     &    '  Sum of global ID''s in even communicator  = ', even_id_sum
      end if
c
c  Use MPI_Reduce to sum the global ID of each process in the odd group.
c  Assuming 4 processes: ODD_SUM = 1 + 3 = 4
c
      if ( odd_id .ne. -1 ) then
        call MPI_Reduce ( id, odd_id_sum,  1, MPI_INTEGER, MPI_SUM, 0, 
     &    odd_comm_id, ierr )
      end if

      if ( odd_id .eq. 0 ) then
        write ( *, '(a,i4)' ) 
     &    '  Number of processes in odd communicator  = ', odd_p
        write ( *, '(a,i4)' ) 
     &    '  Sum of global ID''s in odd communicator  = ', odd_id_sum
      end if
c
c  Terminate MPI.
c
      call MPI_Finalize ( ierr )
c
c  Terminate
c
      if ( id .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'COMMUNICATOR_MPI:'
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
c  Discussion:
c
c    This FORTRAN77 version is made available for cases where the
c    FORTRAN90 version cannot be used.
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
