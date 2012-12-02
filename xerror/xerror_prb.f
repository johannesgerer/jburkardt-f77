      program main

c*********************************************************************72
c
cc MAIN is the main program for XERROR_PRB.
c
c  Discussion:
c
c    XERROR_PRB calls sample problems for the XERROR library.
c
c  Modified:
c
c    04 April 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'XERROR_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the XERROR library.'

      call test01 ( )
      call test02 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'XERROR_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01

c*********************************************************************72
c
cc TEST01 tests XERROR.
c
c  Discussion:
c
c    This simple example just invokes XERROR, which is responsible
c    for printing out the message, recording the error number, and
c    taking appropriate action if the error level is high enough.
c
c  Modified:
c
c    04 April 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer level
      integer nerr
      real r4_uniform_01
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )
      real x
      real y

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  XERROR can be called when an error occurs.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X        Sqrt ( X ) '
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        x = 2.0E+00 * r4_uniform_01 ( seed ) - 1.0E+00

        if ( 0.0E+00 <= x ) then

          y = sqrt ( x )
          write ( *, '(2x,f10.6,2x,f10.6)' ) x, y

        else

          y = 0.0E+00
          nerr = 1
          level = 0

          call xerror ( 'TEST01 - Illegal argument to SQRT', 33,
     &      nerr, level )

        end if

      end do

      return
      end
      subroutine test02

c*********************************************************************72
c
cc TEST02 tests XERROR.
c
c  Discussion:
c
c    In this example, we want to reduce the amount of error output.
c    By setting KONTRL to 0, we suppress printout of everything but
c    fatal errors.
c
c  Modified:
c
c    06 April 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer kontrl
      integer level
      integer nerr
      real r4_uniform_01
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )
      real x
      real y

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  XGETF gets the error control flag.'
      write ( *, '(a)' ) '  XSETF sets the error control flag.'
      write ( *, '(a)' ) '  XERDMP prints an error message summary.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  By setting the error control flag KONTRL'
      write ( *, '(a)' ) '  to 0, we may reduce the number of error'
      write ( *, '(a)' ) '  messages printed out.'

      call xgetf ( kontrl )
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Current value of KONTRL was ', kontrl

      kontrl = 0
      call xsetf ( kontrl )

      write ( *, '(a,i8)' ) '  Resetting KONTRL to         ', kontrl

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X        Sqrt ( X ) '
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        x = 2.0E+00 * r4_uniform_01 ( seed ) - 1.0E+00

        if ( 0.0E+00 <= x ) then

          y = sqrt ( x )
          write ( *, '(2x,f10.6,2x,f10.6)' ) x, y

        else

          y = 0.0E+00
          nerr = 1
          level = 0

          call xerror ( 'TEST01 - Illegal argument to SQRT', 33,
     &      nerr, level )

        end if

      end do
c
c  At the end, we can call XERDMP to see what error messages we avoided.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Call XERDMP to see what we missed.'

      call xerdmp ( )

      return
      end
      function r4_uniform_01 ( seed )

c*********************************************************************72
c
cc R4_UNIFORM_01 returns a unit pseudorandom R4.
c
c  Discussion:
c
c    This routine implements the recursion
c
c      seed = 16807 * seed mod ( 2**31 - 1 )
c      r4_uniform_01 = seed / ( 2**31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      R4_UNIFORM_01
c      SEED      SEED
c
c         12345   207482415  0.096616
c     207482415  1790989824  0.833995
c    1790989824  2035175616  0.947702
c
c  Modified:
c
c    11 August 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Springer Verlag, pages 201-202, 1983.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley Interscience, page 95, 1998.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, pages 362-376, 1986.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, pages 136-143, 1969.
c
c  Parameters:
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, real R4_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      integer k
      integer seed
      real r4_uniform_01

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R4_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + 2147483647
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r4_uniform_01 = real ( dble ( seed ) * 4.656612875D-10 )

      return
      end
