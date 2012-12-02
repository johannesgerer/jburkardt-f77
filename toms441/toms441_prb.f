      program main

c*********************************************************************72
c
cc TOMS441_PRB tests DIPOLE.
c
c  Modified:
c
c    12 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer sample_num
      integer test_num

      parameter ( sample_num = 1000 )
      parameter ( test_num = 3 )

      real a
      real alpha_test(test_num)
      real alpha
      real b
      real dipole
      integer i
      real mean
      real r
      real r_test(test_num)
      integer seed
      integer test
      real variance
      real x(sample_num)
      real xmax
      real xmin

      save a_test
      save b_test

      data a_test / 0.0, 0.785398163397448, 1.57079632679490 /
      data b_test / 1.0, 0.5, 0.0 /

      seed = 123456789

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS441_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS441 library.'

      do test = 1, test_num

        alpha = alpha_test(test)
        r = r_test(test)

        a = r * cos ( alpha )
        b = r * sin ( alpha )

        write ( *, '(a)' ) ' '
        write ( *, '(a,g14.6)' ) ' A = ', a
        write ( *, '(a,g14.6)' ) ' B = ', b

        xmax = -1.0E+30
        xmin = +1.0E+30
        mean = 0.0E+00

        do i = 1, sample_num
          x(i) = dipole ( a, b, seed )
          xmax = max ( xmax, x(i) )
          xmin = min ( xmin, x(i) )
          mean = mean + x(i)
        end do

        mean = mean / real ( sample_num )

        variance = 0.0E+00
        do i = 1, sample_num
          variance = variance + ( x(i) - mean )**2
        end do
        variance = variance / real ( sample_num - 1 )

        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)'    ) '  Sample size =     ', sample_num
        write ( *, '(a,g14.6)' ) '  Sample mean =     ', mean
        write ( *, '(a,g14.6)' ) '  Sample variance = ', variance
        write ( *, '(a,g14.6)' ) '  Sample maximum =  ', xmax
        write ( *, '(a,g14.6)' ) '  Sample minimum =  ', xmin

      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS441_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      function r11 ( seed )

c*********************************************************************72
c
cc R11 returns a pseudorandom number between -1 and +1.
c
c  Modified:
c
c    06 December 2005
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, real RD11, a new pseudorandom variate, strictly between -1 and 1.
c
      implicit none

      integer k
      real r11
      integer seed

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + 2147483647
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r11 = 2.0 * real ( dble ( seed ) * 4.656612875D-10 ) - 1.0

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Discussion:
c
c    This FORTRAN77 version is made available for cases where the
c    FORTRAN90 version cannot be used.
c
c  Modified:
c
c    16 September 2005
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

      character ( len = 8 ) date
      character ( len = 10 ) time

      call date_and_time ( date, time )

      write ( *, '(a8,2x,a10)' ) date, time

      return
      end
