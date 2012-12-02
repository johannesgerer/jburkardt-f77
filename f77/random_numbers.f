      program main

c*********************************************************************72
c
cc MAIN is the main program for RANDOM_NUMBERS.
c
c  Discussion:
c
c    FORTRAN77 does not have a built in random number generator.
c    One way to deal with this is to write your own.  Here is an example.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 April 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RANDOM_NUMBERS'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Generate some random numbers.'

      call test01 ( )
      call test02 ( )
      call test03 ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RANDOM_NUMBERS'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests R8_UNIFORM_01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      double precision r8_uniform_01
      integer seed
      integer seed_init

      seed_init = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  R8_UNIFORM_01 computes pseudorandom values '
      write ( *, '(a)' ) '  in the interval [0,1].'

      seed = seed_init

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
      write ( *, '(a)' ) ' '

      do i = 1, 10
        write ( *, '(2x,i8,2x,g14.6)' ) i, r8_uniform_01 ( seed )
      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests R8_UNIFORM_01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n

      parameter ( n = 1000 )

      integer i
      double precision r8_uniform_01
      integer seed
      integer seed_in
      integer seed_out
      double precision u(n)
      double precision u_avg
      double precision u_var

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  R8_UNIFORM_01 computes a sequence of '
      write ( *, '(a)' ) '  uniformly distributed pseudorandom numbers.'
c
c  Start with a known seed.
c
      seed = 12345

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  Initial SEED = ', seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  First 10 values:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '     I         Input        Output    R8_UNIFORM_01'
      write ( *, '(a)' ) '                SEED          SEED'
      write ( *, '(a)' ) ' '

      do i = 1, 10

        seed_in = seed
        u(i) = r8_uniform_01 ( seed )
        seed_out = seed

        write ( *, '(i8,2x,i12,2x,i12,2x,g14.6)' ) 
     &    i, seed_in, seed_out, u(i)

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,i10,a)' ) '  Now compute ', n, ' elements.'
      write ( *, '(a)' ) ' '

      u_avg = 0.0D+00
      do i = 1, n
        u(i) = r8_uniform_01 ( seed )
        u_avg = u_avg + u(i)
      end do

      u_avg = u_avg / dble ( n )

      u_var = 0.0D+00
      do i = 1, n
        u_var = u_var + ( u(i) - u_avg )**2
      end do

      u_var = u_var / dble ( n - 1 )

      write ( *, '(a)' ) ' '
      write ( *, '(a,f10.6)' ) '  Average value = ', u_avg
      write ( *, '(a,f10.6)' ) '  Expecting       ', 0.5D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,f10.6)' ) '  Variance =      ', u_var
      write ( *, '(a,f10.6)' ) '  Expecting       ', 1.0D+00 / 12.0D+00

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests R8_UNIFORM_01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      double precision r8_uniform_01
      integer seed
      integer seed_in
      integer seed_out
      integer seed_save
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  R8_UNIFORM_01 computes a sequence of '
      write ( *, '(a)' ) '  pseudorandom numbers but all computations'
      write ( *, '(a)' ) '   depend on the seed value.'
      write ( *, '(a)' ) '  In this test, we show how a sequence of '
      write ( *, '(a)' ) '  "random" values can be manipulated by '
      write ( *, '(a)' ) '  accessing the seed.'

      seed = 1066

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  Initial SEED is ', seed
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Call R8_UNIFORM_01 10 times, and watch SEED.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '     I         Input        Output    R8_UNIFORM_01'
      write ( *, '(a)' ) '                SEED          SEED'
      write ( *, '(a)' ) ' '

      do i = 1, 10

        seed_in = seed

        if ( i == 5 ) then
          seed_save = seed
        end if

        x = r8_uniform_01 ( seed )
        seed_out = seed
        write ( *, '(i8,2x,i12,2x,i12,2x,g14.6)' ) 
     &    i, seed_in, seed_out,x

      end do

      seed = seed_save

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) 
     &  '  Reset SEED to its value at step 5, = ', seed
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Now call R8_UNIFORM_01 10 times, and watch'
      write ( *, '(a)' ) 
     &  '  how SEED and R8_UNIFORM_01 restart themselves.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '     I         Input        Output    R8_UNIFORM_01'
      write ( *, '(a)' ) '                SEED          SEED'
      write ( *, '(a)' ) ' '

      do i = 1, 10
        seed_in = seed
        x = r8_uniform_01 ( seed )
        seed_out = seed
        write ( *, '(i8,2x,i12,2x,i12,2x,g14.6)' ) 
     &    i, seed_in, seed_out, x
      end do

      seed = -12345678

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  What happens with a negative SEED?'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '     I         Input        Output    R8_UNIFORM_01'
      write ( *, '(a)' ) '                SEED          SEED'
      write ( *, '(a)' ) ' '

      do i = 1, 10
        seed_in = seed
        x = r8_uniform_01 ( seed )
        seed_out = seed
        write ( *, '(i8,2x,i12,2x,i12,2x,g14.6)' ) 
     &    i, seed_in, seed_out, x
      end do

      return
      end
      function r8_uniform_01 ( seed )

c*********************************************************************72
c
cc R8_UNIFORM_01 returns a unit pseudorandom R8.
c
c  Discussion:
c
c    This routine implements the recursion
c
c      seed = 16807 * seed mod ( 2**31 - 1 )
c      r8_uniform_01 = seed / ( 2**31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      R8_UNIFORM_01
c      SEED      SEED
c
c         12345   207482415  0.096616
c     207482415  1790989824  0.833995
c    1790989824  2035175616  0.947702
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
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
c    Second Edition,
c    Springer, 1987,
c    ISBN: 0387964673,
c    LC: QA76.9.C65.B73.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, December 1986, pages 362-376.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley, 1998,
c    ISBN: 0471134031,
c    LC: T57.62.H37.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, Number 2, 1969, pages 136-143.
c
c  Parameters:
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R8_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer k
      double precision r8_uniform_01
      integer seed

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + i4_huge
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r8_uniform_01 = dble ( seed ) * 4.656612875D-10

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

