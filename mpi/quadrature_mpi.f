      program main

c*********************************************************************72
c
cc MAIN is the main program for QUADRATURE.
c
c  Discussion:
c
c    QUADRATURE estimates an integral using quadrature.
c
c    The integral of F(X) = 4 / ( 1 + X * X ) from 0 to 1 is PI.
c
c    We break up the interval [0,1] into N subintervals, evaluate
c    F(X) at the midpoint of each subinterval, and multiply the
c    sum of these values by N to get an estimate for the integral.
c
c    If we have M processes available because we are using MPI, then
c    we can ask processes 0, 1, 2, ... M-1 to handle the subintervals
c    in the following order:
c
c          0      1       2            M-1  <-- Process numbers begin at 0
c     ------ ------  ------  -----  ------
c          1      2       3    ...       M
c        M+1    M+2     M+3    ...     2*M
c      2*M+1    2*M+2 2*M+3    ...     3*M
c                              
c    and so on up to subinterval N.  The partial sums collected by 
c    each process are then sent to the master process to be added 
c    together to get the estimated integral.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 October 2007
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

      double precision f
      double precision h
      integer i
      integer id
      integer ierr
      integer master
      parameter ( master = 0 )
      integer n
      integer n_part
      integer p
      double precision q
      double precision q_diff
      double precision q_exact
      parameter ( q_exact = 3.141592653589793238462643D+00 )
      double precision q_part
      double precision wtime
      double precision x
c
c  Establish the MPI environment.
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

      if ( id .eq. master ) then
        call timestamp ( )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'QUADRATURE - Master process:'
        write ( *, '(a)' ) '  FORTRAN77 version'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  An MPI program to estimate an integral.'
        write ( *, '(a,i8)' ) '  The number of processes is ', p

        wtime = MPI_Wtime ( )

      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a)' ) 'Process ', id, ' is active.'
c
c  Assume that the master process just got the value of N from the user.
c  Here, we'll use an assignment statement.
c
      if ( id .eq. master ) then
        n = 1000
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'QUADRATURE - Master process:'
        write ( *, '(a,i8)' ) '  Number of intervals is ', n
      end if
c
c  The master process must broadcast the value of N to all processes.
c
      call MPI_Bcast ( n, 1, MPI_INTEGER, master, MPI_COMM_WORLD, ierr )
c
c  Every process, including the master, now adds up its terms.
c
c  There are N evaluation points total.
c  Each evaluation point is in the center of an interval of width 1/N.
c
      h = 1.0D+00 / dble ( n )

      q_part = 0.0D+00
      n_part = 0

      do i = id+1, n, p

        x = dble ( 2 * i - 1 ) 
     &    / dble ( 2 * n     )

        n_part = n_part + 1
        q_part = q_part + f ( x )

      end do

      q_part = q_part * h

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8))' ) 'QUADRATURE - Process ', id
      write ( *, '(a,g24.16)' ) '  My contribution is ', q_part
c
c  All the partial sums are collected and summed by the master process.
c
      call MPI_Reduce ( q_part, q, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 
     &  master, MPI_COMM_WORLD, ierr )
c
c  The master process scales the total and prints the answer.
c
      if ( id .eq. master ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'QUADRATURE - Master process:'
        write ( *, '(a,g24.16)' ) '  Integral estimate  ', q
        write ( *, '(a,g24.16)' ) '  Exact value is     ', q_exact
        q_diff = abs ( q - q_exact )
        write ( *, '(a,g24.16)' ) '  Error is           ', q_diff

        wtime = MPI_Wtime ( ) - wtime

        write ( *, '(a)' ) ' '
        write ( *, '(a,f14.6)' ) '  Elapsed wall clock seconds = ', 
     &    wtime

      end if
c
c  Terminate MPI.
c
      call MPI_Finalize ( ierr )
c
c  Terminate.
c
      if ( id .eq. master ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'QUADRATURE - Master process:'
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
c  Discussion:
c
c    Integral ( 0 <= X <= 1 ) 4/(1+X*X) dX = PI
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
