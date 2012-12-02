      program main

c*********************************************************************72
c
cc MAIN is the main program for INTERVALS.
c
c  Discussion:
c
c    INTERVALS uses MPI routines to multiprocess a computational task.
c
c    We have a function F(X), an interval [XMIN,XMAX], 
c    and a value N.
c
c    We define N equally spaced points in the interval,
c
c      X(I) = ( ( N - I     ) * XMIN 
c             + (     I - 1 ) * XMAX ) 
c             / ( N     - 1 )
c
c    We thus have N-1 subintervals.
c
c    We assume we have N processors available.
c
c    Processor 0 is designated the master processor, assigned
c    to estimating the integral of F(X) over the entire
c    interval [ X(1), X(N) ].
c
c    For I = 1 to N-1, processor I is assigned the subinterval
c
c      [ X(I), X(I+1) ]
c
c    and then estimates the integral Q(I) of F(X) over that
c    subinterval.
c
c    COMMUNICATION:
c
c    Processor 0 communicates to processor I the endpoints of 
c    the interval it is assigned, and the number of sample points
c    to use in that interval.
c
c    Processor I communicates to processor 0 the computed value of
c    Q(I).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 September 2009
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
 
      double precision end_time
      double precision f
      double precision h
      integer i
      integer ierr
      integer m
      integer master
      parameter ( master = 0 )
      integer n
      double precision pi
      parameter ( pi = 3.141592653589793238462643D+00 )
      integer process
      integer process_id
      integer process_num
      double precision q_global
      double precision q_local
      integer received
      integer source
      double precision start_time
      integer status(MPI_Status_size)
      integer tag
      integer target
      double precision x
      double precision xb(2)
      double precision x_max
      parameter ( x_max = 1.0D+00 )
      double precision x_min 
      parameter ( x_min = 0.0D+00 )
c
c  Establish the MPI environment.
c
      call MPI_Init ( ierr )
c
c  Get this process's ID.
c
      call MPI_Comm_rank ( MPI_COMM_WORLD, process_id, ierr )
c
c  Find out how many processes are available.
c
      call MPI_Comm_size ( MPI_COMM_WORLD, process_num, ierr )
c
c  Say hello (once), and shut down right away unless we
c  have at least 2 processes available.
c
      if ( process_id .eq. master ) then
        call timestamp ( )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'INTERVALS - Master process:'
        write ( *, '(a)' ) '  FORTRAN77 version'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  An MPI example program.'
        write ( *, '(a)' ) '  A quadrature over an interval is done by'
        write ( *, '(a)' ) '  assigning subintervals to processes.'
        write ( *, '(a,i8)' ) 
     &    '  The number of processes is ', process_num

        start_time = MPI_Wtime ( )

        if ( process_num .le. 1 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'INTERVALS - Master process:'
          write ( *, '(a)' ) '  Need at least 2 processes.'
          call MPI_Finalize ( ierr )
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'INTERVALS - Master process:'
          write ( *, '(a)' ) '  Abnormal end of execution.'
          stop
        end if

      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a)' ) 'Process ', process_id, ': Active.'
c
c  Every process could figure out the endpoints of its interval
c  on its own.  But we want to demonstrate communication.  So we
c  assume that the assignment of processes to intervals is done
c  only by the master process, which then tells each process
c  what job it is to do.
c
      if ( process_id .eq. master ) then

        do process = 1, process_num-1

          xb(1) = ( dble ( process_num - process     ) * x_min   
     &            + dble (               process - 1 ) * x_max ) 
     &            / dble ( process_num           - 1 )

          xb(2) = ( dble ( process_num - process - 1 ) * x_min   
     &            + dble (               process     ) * x_max ) 
     &            / dble ( process_num           - 1 )
 
          target = process
          tag = 1

          call MPI_Send ( xb, 2, MPI_DOUBLE_PRECISION, target, tag, 
     &      MPI_COMM_WORLD, ierr )

        end do

      else

        source = master
        tag = 1

        call MPI_Recv ( xb, 2, MPI_DOUBLE_PRECISION, source, tag, 
     &    MPI_COMM_WORLD, status, ierr )
    
      end if
c
c  Wait here until everyone has gotten their assignment.
c
      call MPI_Barrier ( MPI_COMM_WORLD, ierr )

      if ( process_id .eq. master ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'INTERVALS - Master process:'
        write ( *, '(a)' ) '  Subintervals have been assigned.'
      end if
c
c  Every process needs to be told the number of points to use.
c  Since this is the same value for everybody, we use a broadcast.
c  Again, we are doing it in this roundabout way to emphasize that
c  the choice for M could really be made at runtime, by processor 0,
c  and then sent out to the others.
c
      m = 100
      source = master

      call MPI_Bcast ( m, 1, MPI_INTEGER, source, MPI_COMM_WORLD, 
     &  ierr )
c
c  Now, every process EXCEPT 0 computes its estimate of the 
c  integral over its subinterval, and sends the result back
c  to process 0.
c
      if ( process_id .ne. master ) then

        q_local = 0.0D+00

        do i = 1, m

          x = ( dble ( 2 * m - 2 * i + 1 ) * xb(1)   
     &        + dble (         2 * i - 1 ) * xb(2) ) 
     &        / dble ( 2 * m             )

          q_local = q_local + f ( x )

        end do

        q_local = q_local * ( xb(2) - xb(1) ) / dble ( m )

        target = master
        tag = 2

        call MPI_Send ( q_local, 1, MPI_DOUBLE_PRECISION, target, tag, 
     &    MPI_COMM_WORLD, ierr )
c
c  Process 0 expects to receive N-1 partial results.
c
      else

        received = 0
        q_global = 0.0D+00

10      continue

        if ( received < process_num - 1 ) then

          source = MPI_ANY_SOURCE
          tag = 2

          call MPI_Recv ( q_local, 1, MPI_DOUBLE_PRECISION, source, tag, 
     &      MPI_COMM_WORLD, status, ierr )

          q_global = q_global + q_local
          received = received + 1

          go to 10

        end if

      end if
c
c  The master process prints the answer.
c
      if ( process_id .eq. master ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'INTERVALS - Master process:'
        write ( *, '(a,g14.6)' ) 
     &  '  Estimate for PI is ', q_global
        write ( *, '(a,g14.6)' ) 
     &  '  Error is           ', q_global - pi

        end_time = MPI_Wtime ( )

        write ( *, '(a)' ) ' '
        write ( *, '(a,f14.6)' ) '  Elapsed wall clock seconds = ', 
     &    end_time - start_time

      end if
c
c  Terminate MPI.
c
      call MPI_Finalize ( ierr )
c
c  Terminate.
c
      if ( process_id == master ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'INTERVALS - Master process:'
        write ( *, '(a)' ) '  Normal end of execution.'
        write ( *, '(a)' ) ' '
        call timestamp ( )
      end if

      stop
      end
      function f ( x )

c*********************************************************************72
c
cc F is the function we are integrating.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 February 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the argument of the function.
c
c    Output, double precision F, the value of the function.
c
      implicit none

      double precision f
      double precision x

      f = 4.0D+00 / ( 1.0D+00 + x * x )

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
