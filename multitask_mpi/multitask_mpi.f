      program main

c*********************************************************************72
c
cc MAIN is the main program for MPI_MULTITASK.
c
c  Discussion:
c
c    Message tag 1: P0 sends input to P1
c    Message tag 2: P0 sends input to P2
c    Message tag 3: P1 sends output to P0.
c    Message tag 4: P2 sends output to P0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 October 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      include 'mpif.h'

      integer id
      integer ierr
      integer input1
      integer input2
      integer my_id
      integer output1
      integer output2
      integer p
      double precision wtime
c
c  Process 0 is the "monitor".
c  It chooses the inputs, and sends them to the workers.
c  It waits for the outputs.
c  It plots the outputs.
c
      call MPI_Init ( ierr )

      call MPI_Comm_rank ( MPI_COMM_WORLD, id, ierr )

      call MPI_Comm_size ( MPI_COMM_WORLD, p, ierr )
c
c  Make sure we have enough processes.
c
      if ( p .lt. 3 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MPI_MULTITASK - Fatal error!'
        write ( *, '(a)' ) '  Number of processes must be at least 3.'
        call MPI_Finalize ( ierr )
        stop
      end if
c
c  Run program P0 on process 0, and so on.
c
      if ( id .eq. 0 ) then

        call timestamp ( )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MPI_MULTITASK:'
        write ( *, '(a)' ) '  FORTRAN77 / MPI version'

        wtime = MPI_Wtime ( )

        call p0_set_input ( input1, input2 )
        call p0_send_input ( input1, input2 )
        call p0_receive_output ( output1, output2 )
c
c  Get timing.
c
        wtime = MPI_Wtime ( ) - wtime
        write ( *, '(a,g14.6)' ) '  Process 0 time = ', wtime

        call MPI_Finalize ( ierr )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MPI_MULTITASK:'
        write ( *, '(a)' ) '  Normal end of execution.'

        call timestamp ( )

        stop
c
c  Process 1 works on task 1.
c  It receives input from process 0.
c  It computes the output.
c  It sends the output to process 0.
c
      else if ( id .eq. 1 ) then

        wtime = MPI_Wtime ( )
        call p1_receive_input ( input1 )
        call p1_compute_output ( input1, output1 )
        call p1_send_output ( output1 )
        wtime = MPI_Wtime ( ) - wtime
        write ( *, '(a,g14.6)' ) '  Process 1 time = ', wtime
        call MPI_Finalize ( ierr )
        stop
c
c  Process 2 works on task 2.
c  It receives input from process 0.
c  It computes the output.
c  It sends the output to process 0.
c
      else if ( id .eq. 2 ) then

        wtime = MPI_Wtime ( )
        call p2_receive_input ( input2 )
        call p2_compute_output ( input2, output2 )
        call p2_send_output ( output2 )
        wtime = MPI_Wtime ( ) - wtime
        write ( *, '(a,g14.6)' ) '  Process 2 time = ', wtime
        call MPI_Finalize ( ierr )
        stop

      end if

      end
      subroutine p0_set_input ( input1, input2 )

c*********************************************************************72
c
cc P0_SET_INPUT sets input.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer INPUT1, INPUT2, the values of two
c    inputs used by tasks 1 and 2.
c
      implicit none

      integer input1
      integer input2

      input1 = 10000000
      input2 = 100000

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) 'P0_SET_PARAMETERS:'
      write ( *, '(a,i12)' ) '  Set INPUT1 = ', input1
      write ( *, '(a,i12)' ) '      INPUT2 = ', input2

      return
      end
      subroutine p0_send_input ( input1, input2 )

c*********************************************************************72
c
cc P0_SEND_INPUT sends input to processes 1 and 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer INPUT1, INPUT2, the values of two
c    inputs used by tasks 1 and 2.
c
      implicit none

      include 'mpif.h'

      integer id
      integer ierr
      integer input1
      integer input2
      integer tag

      id = 1
      tag = 1
      call MPI_Send ( input1, 1, MPI_INTEGER, id, tag, MPI_COMM_WORLD, 
     &  ierr )

      id = 2
      tag = 2
      call MPI_Send ( input2, 1, MPI_INTEGER, id, tag, MPI_COMM_WORLD, 
     &  ierr )

      return
      end
      subroutine p0_receive_output ( output1, output2 )

c*********************************************************************72
c
cc P0_RECEIVE_OUTPUT receives output from processes 1 and 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer OUTPUT1, OUTPUT2, the values of the
c    outputs of tasks 1 and 2.
c
      implicit none

      include 'mpif.h'

      integer ierr
      integer output
      integer output_received
      integer output1
      integer output2
      integer source
      integer status(MPI_STATUS_SIZE)

      output_received = 0
c
c  Loop until every worker has checked in.
c
10    continue

      if ( output_received .lt. 2 ) then
c
c  Receive the next message that arrives.
c
        call MPI_Recv ( output, 1, MPI_INTEGER, MPI_ANY_SOURCE, 
     &    MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr )
c
c  The actual source of the message is saved in STATUS.
c
        source = status(MPI_SOURCE)
c
c  Save the value in OUTPUT1 or OUTPUT2.
c
        if ( source == 1 ) then
          output1 = output
        else
          output2 = output
        end if

        output_received = output_received + 1

        go to 10

      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Process 1 returned OUTPUT1 = ', output1
      write ( *, '(a,i8)' ) '  Process 2 returned OUTPUT2 = ', output2

      return
      end
      subroutine p1_receive_input ( input1 )

c*********************************************************************72
c
cc P1_RECEIVE_INPUT receives input from process 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer INPUT1, the value of the parameter.
c
      implicit none

      include 'mpif.h'

      integer id
      integer ierr
      integer input1
      integer status(MPI_STATUS_SIZE)
      integer tag

      id = 0
      tag = 1
      call MPI_Recv ( input1, 1, MPI_INTEGER, id, tag, MPI_COMM_WORLD, 
     &  status, ierr )

      return
      end
      subroutine p1_compute_output ( input1, output1 )

c*********************************************************************72
c
cc P1_COMPUTE_OUTPUT carries out computation number 1.
c
c  Discussion:
c
c    No MPI calls occur in this function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer INPUT1, the problem input.
c
c    Output, integer OUTPUT1, the problem output.
c
      implicit none

      integer i
      integer j
      integer k
      integer input1
      integer output1

      output1 = 0

      do i = 2, input1

        j = i
        k = 0

10      continue

        if ( 1 .lt. j ) then

          if ( mod ( j, 2 ) .eq. 0 ) then
            j = j / 2
          else
            j = 3 * j + 1
          end if

          k = k + 1

          go to 10

        end if

        output1 = max ( output1, k )

      end do

      return
      end
      subroutine p1_send_output ( output1 )

c*********************************************************************72
c
cc P1_SEND_OUTPUT sends output to process 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OUTPUT1, the problem output.
c
      implicit none

      include 'mpif.h'

      integer id
      integer ierr
      integer output1
      integer tag

      id = 0
      tag = 3
      call MPI_Send ( output1, 1, MPI_INTEGER, id, tag, MPI_COMM_WORLD, 
     &  ierr )

      return
      end
      subroutine p2_receive_input ( input2 )

c*********************************************************************72
c
cc P2_RECEIVE_INPUT receives input from process 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer INPUT2, the value of the parameter.
c
      implicit none

      include 'mpif.h'

      integer id
      integer ierr
      integer input2
      integer status(MPI_STATUS_SIZE)
      integer tag

      id = 0
      tag = 2
      call MPI_Recv ( input2, 1, MPI_INTEGER, id, tag, MPI_COMM_WORLD, 
     &  status, ierr )

      return
      end
      subroutine p2_compute_output ( input2, output2 )

c*********************************************************************72
c
cc P2_COMPUTE_OUTPUT carries out computation number 2.
c
c  Discussion:
c
c    No MPI calls occur in this function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer INPUT2, the problem input.
c
c    Output, integer OUTPUT2, the problem output.
c
      implicit none

      integer i
      integer j
      integer input2
      integer output2
      logical prime

      output2 = 0

      do i = 2, input2

        prime = .true.
        do j = 2, i - 1
          if ( mod ( i, j ) .eq. 0 ) then
            prime = .false.
            go to 10
          end if
        end do

10      continue

        if ( prime ) then
          output2 = output2 + 1
        end if

      end do

      return
      end
      subroutine p2_send_output ( output2 )

c*********************************************************************72
c
cc P2_SEND_OUTPUT sends output to process 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer OUTPUT2, the problem output.
c
      implicit none

      include 'mpif.h'

      integer id
      integer ierr
      integer output2
      integer tag

      id = 0
      tag = 4
      call MPI_Send ( output2, 1, MPI_INTEGER, id, tag, MPI_COMM_WORLD, 
     &  ierr )

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
