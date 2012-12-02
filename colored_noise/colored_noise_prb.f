      program main

c*********************************************************************72
c
cc MAIN generates colored noise data for a sequence of values of ALPHA.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 June 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision alpha
      integer i
      integer n
      double precision q_d
      integer seed_init

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'COLORED_NOISE_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the COLORED_NOISE library.'

      n = 128
      q_d = 1.0D+00
      alpha = 0.00D+00
      seed_init = 123456789

      do i = 0, 8
        alpha = 0.25D+00 * dble ( i )
        call test01 ( n, q_d, alpha, seed_init )
      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'COLORED_NOISE_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, * ) ' '
      call timestamp ( )

      return
      end
      subroutine test01 ( n, q_d, alpha, seed_init )

c*****************************************************************************80
c
cc TEST01 calls F_ALPHA with particular parameters.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of elements of the sequence 
c    to generate.
c
c    Input, double precision Q_D, the variance of the sequence.
c
c    Input, double precision ALPHA, the exponent of the power law.
c
c    Input, integer SEED_INIT, the initial seed for the 
c    random number generator.
c
      implicit none

      integer n

      double precision alpha
      integer i
      character * ( 255 ) output_filename
      integer output_unit
      double precision q_d
      integer seed
      integer seed_init
      double precision x(n)

      write ( output_filename, '(a,f4.2,a)' ) "alpha_", alpha, '.txt'
c
c  Report parameters.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a,i8,a)' ) '  Generating ', n, ' sample points.'
      write ( *, '(a,g14.6)' ) '  1/F^ALPHA noise has ALPHA = ', alpha
      write ( *, '(a,g14.6)' ) '  Variance is ', q_d
      write ( *, '(a,i12)' ) 
     &  '  Initial random number seed = ', seed_init

      seed = seed_init

      call f_alpha ( n, q_d, alpha, seed, x )
c
c  Print no more than 10 entries of the data.
c
      call r8vec_print_part ( n, x, 10, '  Noise sample:' )
c
c  Write the data to a file.
c
      call get_unit ( output_unit )

      open ( unit = output_unit, file = output_filename, 
     &  status = 'replace', err = 10 )

      do i = 1, n
        write ( output_unit, '(g14.6)' ) x(i)
      end do

      close ( unit = output_unit )

      write ( *, '(a)' ) '  Data written to file "' 
     &  // trim ( output_filename ) // '."'

      return

10    continue

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01 - Fatal errorc'
      write ( *, '(a)' ) '  Could not open the output file.'
      stop

      end
