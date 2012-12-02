      program main

c*********************************************************************72
c
cc MAIN is the main program for SCHEDULE_OPENMP.
c
c  Discussion:
c
c    This program demonstrates the difference between default,
c    static and dynamic scheduling for a loop parallelized in OpenMP.
c
c    The purpose of scheduling is to deal with loops in which there is
c    known or suspected imbalance in the work load.  In this example,
c    if the work is divided in the default manner between two threads,
c    the second thread has 3 times the work of the first.  
c
c    Both static and dynamic scheduling, if used, even out the work
c    so that both threads have about the same load.  This could be
c    expected to decrease the run time of the loop by about 1/3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      include 'omp_lib.h'

      integer n
      integer n_factor
      integer n_hi
      integer n_lo
      integer primes
      double precision time1
      double precision time2
      double precision time3

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCHEDULE_OPENMP'
      write ( *, '(a)' ) '  FORTRAN77/OpenMP version'
      write ( *, '(a)' ) '  Count the primes from 1 to N.'
      write ( *, '(a)' ) '  This is an unbalanced work load, particular'
      write ( *, '(a)' ) '  for two threads.  Demonstrate default,'
      write ( *, '(a)' ) '  static and dynamic scheduling.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  Number of processors available = ', omp_get_num_procs ( )
      write ( * ,'(a,i8)' ) 
     &  '  Number of threads =              ', omp_get_max_threads ( )

      n_lo = 1
      n_hi = 131072
      n_factor = 2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '                           Default' //
     &  '        Static       Dynamic'
      write ( *, '(a)' ) '         N     Pi(N)          Time' //
     &  '          Time          Time'
      write ( *, '(a)' ) ' '

      n = n_lo

10    continue

      if ( n .le. n_hi ) then

        time1 = omp_get_wtime ( )
        call prime_default ( n, primes )
        time1 = omp_get_wtime ( ) - time1

        time2 = omp_get_wtime ( )
        call prime_static ( n, primes )
        time2 = omp_get_wtime ( ) - time2

        time3 = omp_get_wtime ( )
        call prime_dynamic ( n, primes )
        time3 = omp_get_wtime ( ) - time3

        write ( *, '(2x,i8,2x,i8,g14.6,g14.6,g14.6)' ) 
     &    n, primes, time1, time2, time3

        n = n * n_factor

        go to 10

      end if
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCHEDULE_OPENMP'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
      subroutine prime_default ( n, total )

c*********************************************************************72
c
cc PRIME_DEFAULT counts primes, using default scheduling.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the maximum number to check.
c
c    Output, integer TOTAL, the number of prime numbers up to N.
c
      implicit none

      integer i
      integer j
      integer n
      integer prime
      integer total

      total = 0

c$omp parallel 
c$omp&  shared ( n ) 
c$omp&  private ( i, j, prime )  

c$omp do reduction ( + : total )

      do i = 2, n

        prime = 1

        do j = 2, i - 1
          if ( mod ( i, j ) .eq. 0 ) then
            prime = 0
            go to 10
          end if
        end do

10      continue

        total = total + prime

      end do

c$omp end do

c$omp end parallel

      return
      end
      subroutine prime_static ( n, total )

c*********************************************************************72
c
cc PRIME_STATIC counts primes, using static scheduling.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the maximum number to check.
c
c    Output, integer TOTAL, the number of prime numbers up to N.
c
      implicit none

      integer i
      integer j
      integer n
      integer prime
      integer total

      total = 0

c$omp parallel 
c$omp&  shared ( n ) 
c$omp&  private ( i, j, prime ) 

c$omp do reduction ( + : total ) schedule ( static, 100 )

      do i = 2, n

        prime = 1

        do j = 2, i - 1
          if ( mod ( i, j ) .eq. 0 ) then
            prime = 0
            go to 10
          end if
        end do

10      continue

        total = total + prime

      end do

c$omp end do

c$omp end parallel

      return
      end
      subroutine prime_dynamic ( n, total )

c*********************************************************************72
c
cc PRIME_DYNAMIC counts primes, using dynamic scheduling.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the maximum number to check.
c
c    Output, integer TOTAL, the number of prime numbers up to N.
c
      implicit none

      integer i
      integer j
      integer n
      integer prime
      integer total

      total = 0

c$omp parallel 
c$omp&  shared ( n ) 
c$omp&  private ( i, j, prime ) 

c$omp do reduction ( + : total ) schedule ( dynamic, 10 )

      do i = 2, n

        prime = 1

        do j = 2, i - 1
          if ( mod ( i, j ) .eq. 0 ) then
            prime = 0
            go to 10
          end if
        end do

10      continue

        total = total + prime

      end do

c$omp end do

c$omp end parallel

      return
      end
