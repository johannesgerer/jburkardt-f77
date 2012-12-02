      program main

c*********************************************************************72
c
cc MAIN is the main program for PRIME_NUMBER_OPENMP.
c
c  Discussion:
c
c    This program calls a version of PRIME_NUMBER that includes
c    OpenMP directives for parallel processing.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 May 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      include 'omp_lib.h'

      integer i
      integer j
      integer n
      integer n_factor
      integer n_hi
      integer n_lo
      integer primes
      double precision wtime

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PRIME_NUMBER_OPENMP'
      write ( *, '(a)' ) '  FORTRAN77/OpenMP version'

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  Number of processors available = ', omp_get_num_procs ( )
      write ( * ,'(a,i8)' ) 
     &  '  Number of threads =              ', omp_get_max_threads ( )

      n_lo = 1
      n_hi = 131072
      n_factor = 2

      call prime_number_sweep_openmp ( n_lo, n_hi, n_factor )

      n_lo = 5
      n_hi = 500000
      n_factor = 10

      call prime_number_sweep_openmp ( n_lo, n_hi, n_factor )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PRIME_NUMBER_OPENMP'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
      subroutine prime_number_sweep_openmp ( n_lo, n_hi, n_factor )

c*********************************************************************72
c
cc PRIME_NUMBER_SWEEP_OPENMP does repeated calls to PRIME_NUMBER.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 August 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N_LO, the first value of N.
c
c    Input, integer N_HI, the last value of N.
c
c    Input, integer N_FACTOR, the factor by which to increase N after
c    each iteration.
c
      implicit none

      include 'omp_lib.h'

      integer i
      integer n
      integer n_factor
      integer n_hi
      integer n_lo
      integer primes
      double precision wtime

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         N        Pi          Time'
      write ( *, '(a)' ) ' '

      n = n_lo

10    continue

      if ( n <= n_hi ) then

        wtime = omp_get_wtime ( )

        call prime_number ( n, primes )

        wtime = omp_get_wtime ( ) - wtime

        write ( *, '(2x,i8,2x,i8,g14.6)' ) n, primes, wtime

        n = n * n_factor

        go to 10

      end if
     
      return
      end
      subroutine prime_number ( n, total )

c*********************************************************************72
c
cc PRIME_NUMBER returns the number of primes between 1 and N.
c
c  Discussion:
c
c    A naive algorithm is used.
c
c    Mathematica can return the number of primes less than or equal to N
c    by the command PrimePi[N].
c
c                N  PRIME_NUMBER
c
c                1           0
c               10           4
c              100          25
c            1,000         168
c           10,000       1,229
c          100,000       9,592
c        1,000,000      78,498
c       10,000,000     664,579
c      100,000,000   5,761,455
c    1,000,000,000  50,847,534
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 May 2009
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
          if ( mod ( i, j ) == 0 ) then
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
