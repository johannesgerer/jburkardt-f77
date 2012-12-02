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
      double precision wtime_diff
      double precision wtime_end
      double precision wtime_start
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
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'QUADRATURE - Master process:'
        write ( *, '(a)' ) '  FORTRAN77 version'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  An MPI program to estimate an integral.'
        write ( *, '(a,i8)' ) '  The number of processes is ', p

        wtime_start = MPI_Wtime ( )

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

        wtime_end = MPI_Wtime ( )
        wtime_diff = end_time - start_time

        write ( *, '(a)' ) ' '
        write ( *, '(a,f14.6)' ) '  Elapsed wall clock seconds = ', 
     &    wtime_diff

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
