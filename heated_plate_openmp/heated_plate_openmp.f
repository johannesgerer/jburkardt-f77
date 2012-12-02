      program main

c*********************************************************************72
c
cc MAIN is the main program for HEATED_PLATE_OPENMP.
c
c  Discussion:
c
c    This code solves the steady state heat equation on a rectangular region.
c
c    The sequential version of this program needs approximately
c    18/eps iterations to complete. 
c
c
c    The physical region, and the boundary conditions, are suggested
c    by this diagram;
c
c                   W = 0
c             +------------------+
c             |                  |
c    W = 100  |                  | W = 100
c             |                  |
c             +------------------+
c                   W = 100
c
c    The region is covered with a grid of M by N nodes, and an N by N
c    array W is used to record the temperature.  The correspondence between
c    array indices and locations in the region is suggested by giving the
c    indices of the four corners:
c
c                  I = 0
c          [0][0]-------------[0][N-1]
c             |                  |
c      J = 0  |                  |  J = N-1
c             |                  |
c        [M-1][0]-----------[M-1][N-1]
c                  I = M-1
c
c    The steady state solution to the discrete heat equation satisfies the
c    following condition at an interior grid point:
c
c      W[Central] = (1/4) * ( W[North] + W[South] + W[East] + W[West] )
c
c    where "Central" is the index of the grid point, "North" is the index
c    of its immediate neighbor to the "north", and so on.
c   
c    Given an approximate solution of the steady state heat equation, a
c    "better" solution is given by replacing each interior point by the
c    average of its 4 neighbors - in other words, by using the condition
c    as an ASSIGNMENT statement:
c
c      W[Central]  <=  (1/4) * ( W[North] + W[South] + W[East] + W[West] )
c
c    If this process is repeated often enough, the difference between successive 
c    estimates of the solution will go to zero.
c
c    This program carries out such an iteration, using a tolerance specified by
c    the user, and writes the final estimate of the solution to a file that can
c    be used for graphic processing.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 October 2011
c
c  Author:
c
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Michael Quinn,
c    Parallel Programming in C with MPI and OpenMP,
c    McGraw-Hill, 2004,
c    ISBN13: 978-0071232654,
c    LC: QA76.73.C15.Q55.
c
c  Local parameters:
c
c    Local, double precision DIFF, the norm of the change in the solution from 
c    one iteration to the next.
c
c    Local, double precision MEAN, the average of the boundary values, used 
c    to initialize the values of the solution in the interior.
c
c    Local, double precision U(M,N), the solution at the previous iteration.
c
c    Local, double precision W(M,N), the solution computed at the latest 
c    iteration.
c
      implicit none

      include 'omp_lib.h'

      integer m
      parameter ( m = 500 )
      integer n
      parameter ( n = 500 )

      double precision diff
      double precision eps
      parameter ( eps = 0.001D+00 )
      integer i
      integer iterations
      integer iterations_print
      integer j
      double precision mean
      double precision u(m,n)
      double precision w(m,n)
      double precision wtime

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HEATED_PLATE_OPENMP'
      write ( *, '(a)' ) '  FORTRAN77/OpenMP version'
      write ( *, '(a)' ) '  A program to solve for the steady state '
      write ( *, '(a)' ) '  temperature distribution over a rectangular'
      write ( *, '(a)' ) '  plate.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a,i8,a)' ) 
     &  '  Spatial grid of ', m, ' by ', n, ' points.'
      write ( *, '(a,g14.6)' ) 
     &  '  The iteration will repeat until the change is <= ', eps
      write ( *, '(a,i8)' ) 
     &  '  The number of processors available = ', omp_get_num_procs ( )
      write ( *, '(a,i8)' ) 
     &  '  The number of threads = ', omp_get_max_threads ( )
c
c  Set the boundary values, which don't change.
c
c  OpenMP Note:
c  You CANNOT set MEAN to zero inside the parallel region.
c
      mean = 0.0
c
c$omp parallel shared ( w ) private ( i, j )

c$omp do
      do i = 2, m - 1
        w(i,1) = 100.0
        w(i,n) = 100.0
      end do
c$omp end do

c$omp do
      do j = 1, n
        w(m,j) = 100.0
        w(1,j) =   0.0
      end do
c$omp end do
c
c  Average the boundary values, to come up with a reasonable
c  initial value for the interior.
c
c$omp do reduction ( + : mean )
      do i = 2, m - 1
        mean = mean + w(i,1) + w(i,n)
      end do
c$omp end do

c$omp do reduction ( + : mean )
      do j = 1, n
        mean = mean + w(1,j) + w(m,j)
      end do
c$omp end do

c$omp end parallel
c
c  OpenMP note:
c  You cannot normalize MEAN inside the parallel region.  It
c  only gets its correct value once you leave the parallel region.
c  So we interrupt the parallel region, set MEAN, and go back in.
c
      mean = mean / dble ( 2 * m + 2 * n - 4 )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  MEAN = ', mean
c
c  Initialize the interior solution to the mean value.
c
c$omp parallel shared ( mean, w ) private ( i, j )

c$omp do
      do j = 2, n - 1
        do i = 2, m - 1
          w(i,j) = mean
        end do
      end do
c$omp end do

c$omp end parallel
c
c  Iterate until the  new solution W differs from the old solution U
c  by no more than EPS.
c
      iterations = 0
      iterations_print = 1

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) ' Iteration  Change'
      write ( *, '(a)' ) ' '

      wtime = omp_get_wtime ( )

      diff = eps

10    continue

      if ( eps .le. diff ) then
c
c  OpenMP node: You CANNOT set DIFF to 0.0 inside the parallel region.
c
        diff = 0.0

c$omp parallel shared ( u, w ) private ( i, j )

c$omp do
        do j = 1, n
          do i = 1, m
            u(i,j) = w(i,j)
          end do
        end do
c$omp end do

c$omp do
        do j = 2, n - 1
          do i = 2, m - 1
            w(i,j) = 0.25 
     &        * ( u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1) )
          end do
        end do
c$omp end do

c$omp do reduction ( max : diff )
        do j = 1, n
          do i = 1, m
            diff = max ( diff, abs ( u(i,j) - w(i,j) ) )
          end do
        end do
c$omp end do

c$omp end parallel

        iterations = iterations + 1

        if ( iterations .eq. iterations_print ) then
          write ( *, '(2x,i8,2x,g14.6)' ) iterations, diff
          iterations_print = 2 * iterations_print
        end if

        go to 10

      end if

      wtime = omp_get_wtime ( ) - wtime

      write ( *, '(a)' ) ' '
      write ( *, '(2x,i8,2x,g14.6)' ) iterations, diff
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Error tolerance achieved.'
      write ( *, '(a,g14.6)' ) '  Wall clock time = ', wtime
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HEATED_PLATE_OPENMP:'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
