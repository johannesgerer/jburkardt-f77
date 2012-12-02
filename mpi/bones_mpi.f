        program main

c*********************************************************************72
c
cc MAIN is the main program for BONES.
c
c  Discussion:
c
c    BONES is a simple demonstration of the use of MPI by an F77 program.
c
c    This program should be run on at least two processes.
c    Only the first two processes will actually get any work to do.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 October 2005
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

      integer count
      real data(0:99)
      integer dest
      integer i
      integer id
      integer ierr
      integer num_procs
      integer status(MPI_Status_size)
      integer tag
      real value(200)
c
c  Initialize MPI.
c
      call MPI_Init ( ierr )
c
c  Determine this process's ID.
c
      call MPI_Comm_rank ( MPI_COMM_WORLD, id, ierr )
c
c  Find out how many processes are active.
c
      call MPI_Comm_size ( MPI_COMM_WORLD, num_procs, ierr )
c
c  Have Process 0 say hello.
c
      if ( id .eq. 0 ) then
        call timestamp ( )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BONES:'
        write ( *, '(a)' ) '  FORTRAN77 version'
        write ( *, '(a)' ) '  A simple MPI test.'
        write ( *, '(a,i6)' ) 
     &    '  The number of processes available is ', num_procs

      end if
c
c  Process 0 expects to receive as much as 200 real values, 
c  from any source.
c
      if ( id .eq. 0 ) then

        tag = 55
        call MPI_Recv ( value, 200, MPI_REAL, MPI_ANY_SOURCE, tag, 
     &    MPI_COMM_WORLD, status, ierr )

        write ( *, '(a,i1,a,i1)' ) 'P:', id, 
     &    ' Got data from processor ', status(MPI_SOURCE)

        call MPI_Get_count ( status, MPI_REAL, count, ierr )

        write ( *, '(a,i1,a,i3,a)' ) 'P:', id, ' Got ', 
     &    count, ' elements.'

        write ( *, '(a,i1,a,g14.6)' ) 'P:', id, 
     &    ' value(5) = ', value(5)
c
c  Process 1 sends 100 real values to process 0.
c
      else if ( id .eq. 1 ) then
 
        write ( *, '(a)' ) ' '
        write ( *, '(a,i1,a)' ) 'P:', id, 
     &    ' - setting up data to send to process 0.'

        do i = 0, 99
          data(i) = real ( i )
        end do

        dest = 0
        tag = 55
        call MPI_Send ( data, 100, MPI_REAL, dest, tag, 
     &    MPI_COMM_WORLD, ierr )
  
      else

        write ( *, '(a)' ) ' '
        write ( *, '(a,i1,a)' ) 'P:', id, 
     &    ' - MPI has no work for me!'

      end if
c
c  Terminate MPI.
c
      call MPI_Finalize ( ierr )
c
c  Terminate.
c
      if ( id .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BONES:'
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
