        program main

c*********************************************************************72
c
cc MAIN is the main program for TYPE.
c
c  Discussion:
c
c    TYPE demonstrates user-defined datatypes in MPI.
c
c    The datatype contains three integers, stored contiguously
c    in memory.
c
c    Process 0 will set up an example of this structure, and send it
c    to Process 1, which will alter it and send it back.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 October 2007
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
      include "mpif.h"

      integer dest
      integer i
      integer ierr
      integer master
      integer my_id
      integer num_procs
      integer point_count
      integer point_type
      integer source
      integer status(MPI_STATUS_SIZE)
      integer tag
      integer x
      integer y
      integer z
c
c  In FORTRAN77, the easiest way to ensure that data is contiguous is to
c  put it all in a vector.  If that is possible, then MPI can handle the
c  data immediately.  But if the items of data are not of the same type,
c  you can't store them all in a single vector.  FORTRAN77 doesn't allow
c  the user to create an appropriate datatype, but it does allow the user
c  to place disparate items in a common block.
c
c  In this example, therefore, we will show how the contents of a common
c  block can be treated as a single data item, that can be sent between
c  processes.  
c
c  In this example, we use MPI_CONTIGUOUS, which assumes that the data
c  items do have the same numeric type.
c
      common /point/ x, y, z

      master = 0
c
c  Initialize MPI.
c
      call MPI_Init ( ierr )
c
c  Get the number of processes.
c
      call MPI_Comm_size ( MPI_COMM_WORLD, num_procs, ierr )
c
c  Get the individual process ID.
c
      call MPI_Comm_rank ( MPI_COMM_WORLD, my_id, ierr )
c
c  Print a message.
c
      if ( my_id .eq. master ) then
        call timestamp ( )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TYPE - Master process:'
        write ( *, '(a)' ) '  FORTRAN77 version'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  An MPI program to set up a datatype.'
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  The number of processes is ', num_procs
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a)' ) '  Process ', my_id, ' is active.'
c
c  Define the new datatype as 3 contiguous integers.
c
      point_count = 3
      call MPI_Type_contiguous ( point_count, MPI_INTEGER, point_type, 
     &  ierr )
c
c  Commit the new datatype.
c  The value POINT_TYPE will now indicate to MPI that we are using
c  an item of the type we defined.
c
      call MPI_Type_commit ( point_type, ierr )
c
c  Use the new datatype.
c
c  The master will send an item of the given type to the workers.
c
      if ( my_id .eq. master ) then

        x = 1
        y = 2
        z = 4
        dest = 1
        tag = 1

        call MPI_Send ( x, 1, point_type, dest, tag, MPI_COMM_WORLD, 
     &    ierr )

        if ( ierr .ne. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a,i8,a,i8)' ) '  Process ', my_id, 
     &      ' got MPI_Send error ', ierr
        end if

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8,a)' ) '  Process ', my_id, 
     &    ' sent an item of type POINT_TYPE'
        write ( *, '(a,i8,2x,i8,2x,i8)' ) '  Value: ', x, y, z

        source = 1
        tag = 2
        call MPI_Recv ( x, 1, point_type, source, tag, MPI_COMM_WORLD, 
     &    status, ierr )

        write ( *, '(a,i8,a)' ) '  Process ', my_id, 
     &    ' received a modified item of type POINT_TYPE' 
        write ( *, '(a,i8,2x,i8,2x,i8)' ) '  Value: ', x, y, z

      else if ( my_id .eq. 1 ) then

        source = 0
        tag = 1

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8,a)' ) '  Process ', my_id, 
     &    ' expecting an item of type POINT_TYPE.'
 
         call MPI_Recv ( x, 1, point_type, source, tag, MPI_COMM_WORLD, 
     &    status, ierr )

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8,a)' ) '  Process ', my_id, 
     &    ' received an item of type POINT_TYPE'
        write ( *, '(a,i8,2x,i8,2x,i8)' ) '  Value: ', x, y, z
 
        i = x
        x = z * 100
        y = y * 10
        z = i
        dest = 0
        tag = 2
        
        call MPI_Send ( x, 1, point_type, dest, tag, MPI_COMM_WORLD, 
     &    ierr )

        write ( *, '(a,i8,a,i8,2x,i8,2x,i8)' ) '  Process ', my_id, 
     &    ' sent a modified item of type POINT_TYPE' 
        write ( *, '(a,i8,2x,i8,2x,i8)' ) '  Value: ', x, y, z

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8,a)' ) '  Process ', my_id, 
     &    ': MPI has nothing for me to do!'

      end if
c
c  Free the memory associated with the datatype.
c
      call MPI_Type_Free ( point_type, ierr )
c
c  Terminate MPI.
c
      call MPI_Finalize ( ierr )
c
c  Terminate.
c
      if ( my_id .eq. master ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TYPE - Master process:'
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
