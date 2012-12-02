      program main

c*********************************************************************72
c
cc MAIN is the main program for COMPUTE_PI.
c
c  Discussion:
c
c    COMPUTE_PI estimates the value of PI.
c
c    This program uses Open MP parallelization directives.  
c
c    It should run properly whether parallelization is used or not.
c
c    However, the parallel version computes the sum in a different
c    order than the serial version; some of the quantities added are
c    quite small, and so this will affect the accuracy of the results.
c
c    The single precision code is noticeably less accurate than the
c    double precision code.  Again, the amount of difference depends
c    on whether the computation is done in parallel or not.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 April 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      include 'omp_lib.h'

      integer r8_logn_max
      parameter ( r8_logn_max = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'COMPUTE_PI'
      write ( *, '(a)' ) '  FORTRAN77/OpenMP version'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Estimate PI by summing a series.'

      write ( *, '(a,i8)' ) 
     &  '  The number of processors = ', omp_get_num_procs ( )
      write ( *, '(a,i8)' ) 
     &  '  The number of threads    = ', omp_get_max_threads ( )

      call r8_test ( r8_logn_max )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'COMPUTE_PI'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
      subroutine r8_test ( logn_max )

c*********************************************************************72
c
cc R8_TEST estimates the value of PI using double precision.
c
c  Discussion:
c
c    PI is estimated using N terms.  N is increased from 10^2 to 10^LOGN_MAX.
c    The calculation is repeated using both sequential and Open MP enabled code.
c    Wall clock time is measured by calling SYSTEM_CLOCK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 April 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      include 'omp_lib.h'

      double precision error
      double precision estimate
      integer logn
      integer logn_max
      character ( len = 3 ) mode
      integer n
      double precision r8_pi
      parameter ( r8_pi = 3.141592653589793D+00 ) 
      double precision wtime

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_TEST:'
      write ( *, '(a)' ) '  Estimate the value of PI,'
      write ( *, '(a)' ) '  using double precision arithmetic.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  N = number of terms computed and added;'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  MODE = SEQ for sequential code;'
      write ( *, '(a)' ) '  MODE = OMP for Open MP enabled code;'
      write ( *, '(a)' ) '  (performance depends on if Open MP is used,'
      write ( *, '(a)' ) '  and how many processes are availablec)'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  ESTIMATE = the computed estimate of PI;'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  ERROR = ( the computed estimate - PI );'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  TIME = elapsed wall clock time;'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  You can''t increase N forever, because:'
      write ( *, '(a)' ) '  A) ROUNDOFF starts to be a problem, and'
      write ( *, '(a)' ) '  B) maximum integer size is a problem.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '             N Mode    Estimate        Error           Time'
      write ( *, '(a)' ) ' '

      n = 1

      do logn = 1, logn_max
c
c  Note that when I set N = 10**LOGN directly, rather than using
c  recursion, I got inaccurate values of N when LOGN was "large",
c  that is, for LOGN = 10, despite the fact that N itself was
c  a KIND = 8 integerc  
c 
c  Sequential calculation.
c
        mode = 'SEQ'

        wtime = omp_get_wtime ( )

        call r8_pi_est_seq ( n, estimate )

        wtime = omp_get_wtime ( ) - wtime

        error = abs ( estimate - r8_pi )

        write ( *, '(i14,2x,a3,2x,f14.10,2x,g14.6,2x,g14.6)' ) 
     &    n, mode, estimate, error, wtime
c
c  Open MP enabled calculation.
c
        mode = 'OMP'

        wtime = omp_get_wtime ( )

        call r8_pi_est_omp ( n, estimate )

        wtime = omp_get_wtime ( ) - wtime

        error = abs ( estimate - r8_pi )

        write ( *, '(i14,2x,a3,2x,f14.10,2x,g14.6,2x,g14.6)' ) 
     &    n, mode, estimate, error, wtime

        n = n * 10

      end do

      return
      end
      subroutine r8_pi_est_omp ( n, estimate )

c*********************************************************************72
c
cc R8_PI_EST_OMP estimates the value of PI, using Open MP.
c
c  Discussion:
c
c    The calculation is based on the formula for the indefinite integral:
c
c      Integral 1 / ( 1 + X**2 ) dx = Arctan ( X ) 
c
c    Hence, the definite integral
c
c      Integral ( 0 <= X <= 1 ) 1 / ( 1 + X**2 ) dx 
c      = Arctan ( 1 ) - Arctan ( 0 )
c      = PI / 4.
c
c    A standard way to approximate an integral uses the midpoint rule.
c    If we create N equally spaced intervals of width 1/N, then the
c    midpoint of the I-th interval is 
c
c      X(I) = (2*I-1)/(2*N).  
c
c    The approximation for the integral is then:
c
c      Sum ( 1 <= I <= N ) (1/N) * 1 / ( 1 + X(I)**2 )
c
c    In order to compute PI, we multiply this by 4; we also can pull out
c    the factor of 1/N, so that the formula you see in the program looks like:
c
c      ( 4 / N ) * Sum ( 1 <= I <= N ) 1 / ( 1 + X(I)**2 )
c
c    Until roundoff becomes an issue, greater accuracy can be achieved by 
c    increasing the value of N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 January 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of terms to add up.
c
c    Output, double precision ESTIMATE, the estimated value of pi.
c
      implicit none

      double precision h
      double precision estimate
      integer i
      integer n
      double precision sum2
      double precision x

      h = 1.0D+00 / dble ( 2 * n )

      sum2 = 0.0D+00

c$omp parallel
c$omp& shared ( h, n ) 
c$omp& private ( i, x )

c$omp do reduction ( + : sum2 )
      do i = 1, n
        x = h * dble ( 2 * i - 1 )
        sum2 = sum2 + 1.0D+00 / ( 1.0D+00 + x * x )
      end do
c$omp end do

c$omp end parallel

      estimate = 4.0D+00 * sum2 / dble ( n )

      return
      end
      subroutine r8_pi_est_seq ( n, estimate )

c*********************************************************************72
c
cc R8_PI_EST_SEQ estimates the value of PI, using sequential execution.
c
c  Discussion:
c
c    The calculation is based on the formula for the indefinite integral:
c
c      Integral 1 / ( 1 + X**2 ) dx = Arctan ( X ) 
c
c    Hence, the definite integral
c
c      Integral ( 0 <= X <= 1 ) 1 / ( 1 + X**2 ) dx 
c      = Arctan ( 1 ) - Arctan ( 0 )
c      = PI / 4.
c
c    A standard way to approximate an integral uses the midpoint rule.
c    If we create N equally spaced intervals of width 1/N, then the
c    midpoint of the I-th interval is 
c
c      X(I) = (2*I-1)/(2*N).  
c
c    The approximation for the integral is then:
c
c      Sum ( 1 <= I <= N ) (1/N) * 1 / ( 1 + X(I)**2 )
c
c    In order to compute PI, we multiply this by 4; we also can pull out
c    the factor of 1/N, so that the formula you see in the program looks like:
c
c      ( 4 / N ) * Sum ( 1 <= I <= N ) 1 / ( 1 + X(I)**2 )
c
c    Until roundoff becomes an issue, greater accuracy can be achieved by 
c    increasing the value of N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 January 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer ( kind = 8 ) N, the number of terms to add up.
c
c    Output, double precision ESTIMATE, the estimated value of pi.
c
      implicit none

      double precision h
      double precision estimate
      integer i
      integer n
      double precision sum2
      double precision x

      h = 1.0D+00 / dble ( 2 * n )

      sum2 = 0.0D+00

      do i = 1, n
        x = h * dble ( 2 * i - 1 )
        sum2 = sum2 + 1.0D+00 / ( 1.0D+00 + x**2 )
      end do

      estimate = 4.0D+00 * sum2 / dble ( n )

      return
      end
