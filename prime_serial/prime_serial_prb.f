      program main

c*********************************************************************72
c
cc MAIN is the main program for PRIME_SERIAL_PRB.
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
      implicit none

      integer n_factor
      integer n_hi
      integer n_lo

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PRIME_SERIAL_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the PRIME_SERIAL library.'

      n_lo = 1
      n_hi = 131072
      n_factor = 2

      call prime_number_sweep ( n_lo, n_hi, n_factor )

      n_lo = 5
      n_hi = 500000
      n_factor = 10

      call prime_number_sweep ( n_lo, n_hi, n_factor )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PRIME_SERIAL_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine prime_number_sweep ( n_lo, n_hi, n_factor )

c*********************************************************************72
c
cc PRIME_SERIAL_NUMBER_SWEEP does repeated calls to PRIME_SERIAL_NUMBER.
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

      integer i
      integer  n
      integer  n_factor
      integer n_hi
      integer n_lo
      integer primes
      double precision time1
      double precision time2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' )
     &  '  Call PRIME_SERIAL_NUMBER to count the primes from 1 to N.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         N        Pi          Time'
      write ( *, '(a)' ) ' '

      n = n_lo

10    continue

      if ( n .le. n_hi ) then

        call cpu_time ( time1 )

        call prime_number ( n, primes )

        call cpu_time ( time2 )

        write ( *, '(2x,i8,2x,i8,g14.6)' ) n, primes, time2 - time1

        n = n * n_factor

        go to 10

      end if

      return
      end


