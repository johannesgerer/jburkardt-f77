      program main

c*********************************************************************72
c
cc MAIN is the main program for R8LIB_PRB.
c
c  Discussion:
c
c    R8LIB_PRB calls sample problems for the R8LIB library.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 June 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8LIB_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the R8LIB library.'

      call test001 ( )
      call test002 ( )
      call test003 ( )
      call test004 ( )
      call test005 ( )
      call test006 ( )
      call test007 ( )
      call test008 ( )
      call test009 ( )

      call test010 ( )
      call test011 ( )
      call test012 ( )
      call test013 ( )
      call test014 ( )
      call test015 ( )
      call test016 ( )
      call test017 ( )
      call test018 ( )
      call test019 ( )

      call test020 ( )
      call test021 ( )
      call test022 ( )
      call test023 ( )
      call test0235 ( )
      call test024 ( )
      call test025 ( )
      call test026 ( )
      call test027 ( )
      call test028 ( )
      call test029 ( )
      call test0295 ( )

      call test031 ( )
      call test032 ( )
      call test033 ( )
      call test034 ( )
      call test035 ( )
      call test036 ( )
      call test0365 ( )
      call test037 ( )
      call test038 ( )
      call test0383 ( )
      call test0385 ( )
      call test039 ( )
      call test0393 ( )
      call test0395 ( )
      call test0397 ( )

      call test049 ( )

      call test0555 ( )

      call test1143 ( )
      call test1145 ( )
      call test1147 ( )
      call test115 ( )

      call test1251 ( )
      call test1252 ( )
      call test1255 ( )
      call test1256 ( )
      call test1258 ( )

      call test1465 ( )

      call test1504 ( )
      call test1515 ( )
      call test157 ( )
      call test158 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8LIB_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test001 ( )

c*********************************************************************72
c
cc TEST001 tests R8_ABS.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 August 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision r8
      double precision r8_abs
      double precision r8_absolute
      double precision r8_uniform_ab
      double precision r8_hi
      double precision r8_lo
      integer seed
      integer test
      integer test_num

      r8_hi = 5.0D+00
      r8_lo = -3.0D+00
      seed = 123456789
      test_num = 10

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST001'
      write ( *, '(a)' ) '  R8_ABS returns the absolute value of an R8.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X         R8_ABS(X)'
      write ( *, '(a)' ) ' '

      do test = 1, test_num
        r8 = r8_uniform_ab ( r8_lo, r8_hi, seed )
        r8_absolute = r8_abs ( r8 )
        write ( *, '(2x,f10.6,2x,f10.6)' ) r8, r8_absolute
      end do

      return
      end
      subroutine test002 ( )

c*********************************************************************72
c
cc TEST002 tests R8_ATAN.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 8 )

      double precision r8_atan
      integer test
      double precision x
      double precision xtest(test_num)
      double precision y
      double precision ytest(test_num)

      save xtest
      save ytest

      data xtest /
     &   1.0D+00,  1.0D+00,  0.0D+00, -1.0D+00,
     &  -1.0D+00, -1.0D+00,  0.0D+00,  1.0D+00 /
      data ytest /
     &   0.0D+00,  1.0D+00,  1.0D+00,  1.0D+00,
     &   0.0D+00, -1.0D+00, -1.0D+00, -1.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST002'
      write ( *, '(a)' )
     &  '  R8_ATAN computes the arc-tangent given Y and X;'
      write ( *, '(a)' )
     &  '  ATAN2 is the system version of this routine.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '       X             Y          ATAN2(Y,X)    R8_ATAN(Y,X)'
      write ( *, '(a)' ) ' '

      do test = 1, test_num
        x = xtest(test)
        y = ytest(test)
        write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' )
     &    x, y, atan2 ( y, x ), r8_atan ( y, x )
      end do

      return
      end
      subroutine test003 ( )

c*********************************************************************72
c
cc TEST003 tests R8_CAS.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision r8_cas
      double precision r8_pi
      integer test_num
      parameter ( test_num = 12 )
      integer test
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST003'
      write ( *, '(a)' ) '  R8_CAS evaluates the casine of a number.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        X           R8_CAS ( X )'
      write ( *, '(a)' ) ' '
      do test = 0, test_num
        x = r8_pi ( ) * dble ( test ) / dble ( test_num )
        write ( *, '(2x,g14.6,2x,g14.6)' ) x, r8_cas ( x )
      end do

      return
      end
      subroutine test004 ( )

c*********************************************************************72
c
cc TEST004 tests R8_CEILING.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      double precision r8_ceiling
      double precision rval
      double precision rval2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST004'
      write ( *, '(a)' ) '  R8_CEILING rounds a value up.'
      write ( *, '(a)' ) ' '

      do i = -6, 6
        rval = dble ( i ) / 5.0D+00
        rval2 = r8_ceiling ( rval )
        write ( *, '(2x,g14.6,2x,g14.6)' ) rval, rval2
      end do

      return
      end
      subroutine test005 ( )

c*********************************************************************72
c
cc TEST005 tests R8_DIFF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 15 )

      integer ndig
      double precision r8_diff
      integer test
      double precision x
      double precision y_test(test_num)
      double precision y

      save y_test

      data y_test /
     &  0.0625D+00, 0.125D+00, 0.25D+00, 0.50D+00,  0.874D+00,
     &  0.876D+00,  0.90D+00,  0.95D+00, 0.99D+00,  1.0D+00,
     &  1.01D+00,   1.05D+00,  1.10D+00, 3.0D+00,  10.0D+00 /

      ndig = 3
      x = 1.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST005'
      write ( *, '(a)' ) '  R8_DIFF computes a difference X-Y to a'
      write ( *, '(a)' ) '  given number of binary places.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a)' )
     &  '  For this test, we use ', ndig, ' binary places.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       X       Y       X-Y     R8_DIFF(X,Y)'
      write ( *, '(a)' ) ' '
      do test = 1, test_num
        y = y_test(test)
        write ( *, '(4f10.4)' ) x, y, x - y, r8_diff ( x, y, ndig )
      end do

      return
      end
      subroutine test006 ( )

c*********************************************************************72
c
cc TEST006 tests R8_DIGIT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer maxdig
      parameter ( maxdig = 20 )

      integer i
      integer digit(-2:maxdig)
      integer idigit
      double precision r8_pi
      double precision x

      x = r8_pi ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST006'
      write ( *, '(a)' ) '  R8_DIGIT extracts decimal digits.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g24.16)' ) '  Here, we get digits of ', x
      write ( *, '(a)' ) ' '

      do idigit = -2, maxdig
        call r8_digit ( x, idigit, digit(idigit) )
      end do

      write ( *, '(2x,25i3)' ) ( i, i = -2, maxdig )
      write ( *, '(2x,25i3)' ) digit(-2:maxdig)

      return
      end
      subroutine test007 ( )

c*********************************************************************72
c
cc TEST007 tests R8_EPSILON and R8_EPSILON_COMPUTE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 September 2012
c
c  Author:
c
c    John Burkardt
c
      double precision r8_epsilon
      double precision r8_epsilon_compute
      double precision r1
      double precision r2
      double precision s

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST007'
      write ( *, '(a)' ) 
     &  '  R8_EPSILON returns the R8 machine precision.'
      write ( *, '(a)' ) 
     &  '  R8_EPSILON_COMPUTE computes the R8 machine precision.'
      write ( *, '(a)' ) ' '

      r1 = r8_epsilon ( )
      write ( *, '(a,g24.16)' ) '  R1 = R8_EPSILON()         = ', r1

      r2 = r8_epsilon_compute ( )
      write ( *, '(a,g24.16)' ) '  R2 = R8_EPSILON_COMPUTE() = ', r2

      s = ( 1.0D+00 + r2 ) - 1.0D+00
      write ( *, '(a,g24.16)' ) '  ( 1 + R2 ) - 1            = ', s

      s = ( 1.0D+00 + ( r2 / 2.0D+00 ) ) - 1.0D+00
      write ( *, '(a,g14.6)' ) '  ( 1 + (R2/2) ) - 1        = ', s

      return
      end
      subroutine test008 ( )

c*********************************************************************72
c
cc TEST008 tests R8_FRACTIONAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 October 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fractional
      double precision r8
      double precision r8_fractional
      double precision r8_uniform_ab
      double precision r8_hi
      double precision r8_lo
      integer seed
      integer test
      integer test_num

      r8_hi = 5.0D+00
      r8_lo = -3.0D+00
      seed = 123456789
      test_num = 10

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST008'
      write ( *, '(a)' )
     &  '  R8_FRACTIONAL returns the fractional part of an R8.'
      write ( *, '(a)' ) ' '

      do test = 1, test_num
        r8 = r8_uniform_ab ( r8_lo, r8_hi, seed )
        fractional = r8_fractional ( r8 )
        write ( *, '(2x,f10.6,2x,f10.6)' ) r8, fractional
      end do

      return
      end
      subroutine test009 ( )

c*********************************************************************72
c
cc TEST009 tests R8_HUGE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision r8_huge

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST009'
      write ( *, '(a)' ) '  R8_HUGE returns a "huge" R8;'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g24.16)' ) '    R8_HUGE ( ) =      ',
     &  r8_huge ( )

      return
      end
      subroutine test010 ( )

c*********************************************************************72
c
cc TEST010 tests R8_LOG_2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 18 )

      double precision r8_log_2
      integer ( kind = 4 ) test
      double precision x
      double precision x_test(test_num)

      save x_test

      data x_test /
     &  0.0D+00,  1.0D+00,  2.0D+00,   3.0D+00,  9.0D+00,
     &  10.0D+00, 11.0D+00, 99.0D+00, 101.0D+00, -1.0D+00,
     &  -2.0D+00, -3.0D+00, -9.0D+00,   0.5D+00,  0.33D+00,
     &   0.25D+00, 0.20D+00, 0.01D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST010'
      write ( *, '(a)' ) '  R8_LOG_2 computes the logarithm base 2.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X, R8_LOG_2'
      write ( *, '(a)' ) ' '

      do test = 1, test_num
        x = x_test(test)
        write ( *, '( 2g14.6 )' ) x, r8_log_2 ( x )
      end do

      return
      end
      subroutine test011 ( )

c*********************************************************************72
c
cc TEST011 tests R8_LOG_B.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 10 )

      double precision b
      double precision b_test(test_num)
      double precision r8_log_b
      integer test
      double precision x

      save b_test

      data b_test /
     &  2.0D+00, 3.0D+00, 4.0D+00, 5.0D+00, 6.0D+00,
     &  7.0D+00, 8.0D+00, 16.0D+00, 32.0D+00, 256.0D+00 /

      x = 16.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST011'
      write ( *, '(a)' ) '  R8_LOG_B computes the logarithm base B.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X, B, R8_LOG_B'
      write ( *, '(a)' ) ' '

      do test = 1, test_num

        b = b_test(test)

        write ( *, '( 2x,3g14.6, i12 )' ) x, b, r8_log_b ( x, b )

      end do

      return
      end
      subroutine test012 ( )

c*********************************************************************72
c
cc TEST012 tests R8_MANT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 December 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer is
      integer l
      double precision r
      double precision x

      x = -314.159D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST012'
      write ( *, '(a)' ) '  R8_MANT decomposes a value.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Number to be decomposed:'
      write ( *, '(2x,g14.6)' ) x

      call r8_mant ( x, is, r, l )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a,g14.6,a,i8)' )
     &  '  R8_MANT: X = ', is, ' * ', r, ' * 2**', l

      return
      end
      subroutine test013 ( )

!*********************************************************************72
!
!! TEST013 tests R8_MOD.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    14 June 2007
!
!  Author:
!
!    John Burkardt
!
      implicit none

      double precision r8_mod
      double precision r8_uniform_ab
      integer test
      integer test_num
      parameter ( test_num = 10 )
      integer seed
      double precision x
      double precision x_hi
      double precision x_lo
      double precision y
      double precision z1
      double precision z2

      x_hi = 10.0D+00
      x_lo = - 10.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST013'
      write ( *, '(a)' )
     &  '  R8_MOD returns the remainder after division.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X         Y     MOD(X,Y)    R8_MOD(X,Y)'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        x = r8_uniform_ab ( x_lo, x_hi, seed )
        y = r8_uniform_ab ( x_lo, x_hi, seed )

        z1 =    mod ( x, y )
        z2 = r8_mod ( x, y )

        write ( * , '(2x,f8.4,2x,f8.4,2x,f8.4,2x,f8.4)' ) x, y, z1, z2

      end do

      return
      end
      subroutine test014 ( )

c*********************************************************************72
c
cc TEST014 tests R8_MODP.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 August 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision r8_modp
      double precision r8_uniform_ab
      integer test
      integer test_num
      integer seed
      double precision x
      double precision x_hi
      double precision x_lo
      double precision y
      double precision z1
      double precision z2

      test_num = 10
      x_hi = 10.0D+00
      x_lo = - 10.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST014'
      write ( *, '(a)' )
     &  '  R8_MODP returns the remainder after division.'
      write ( *, '(a)' )
     & '  Unlike the FORTRAN MOD, R8_MODP ( X, Y ) is positive.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X       Y      MOD(X,Y)  R8_MODP(X,Y)'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        x = r8_uniform_ab ( x_lo, x_hi, seed )
        y = r8_uniform_ab ( x_lo, x_hi, seed )

        z1 =    mod  ( x, y )
        z2 = r8_modp ( x, y )

        write ( * , '(2x,f8.4,2x,f8.4,2x,f8.4,2x,f8.4)' ) x, y, z1, z2

      end do

      return
      end
      subroutine test015 ( )

c*********************************************************************72
c
cc TEST015 tests R8_NINT
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 August 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision b
      double precision c
      integer r8_nint
      double precision r8_uniform_ab
      integer i
      integer seed
      integer test
      integer test_num
      double precision x

      seed = 123456789
      test_num = 10

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST015'
      write ( *, '(a)' )
     &  '  R8_NINT produces the nearest integer to an R8.'
      write ( *, '(a)' ) ' '

      b = -10.0D+00
      c = +10.0D+00

      do test = 1, test_num
        x = r8_uniform_ab ( b, c, seed )
        write ( *, '(2x,f10.4,2x,i8)' ) x, r8_nint ( x )
      end do

      return;
      end
      subroutine test016 ( )

c*********************************************************************72
c
cc TEST016 tests R8_NORMAL_01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 August 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision r8_normal_01
      integer seed
      integer test
      integer test_num
      double precision x

      seed = 123456789
      test_num = 20

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST016'
      write ( *, '(a)' ) '  R8_NORMAL_01 generates normally distributed'
      write ( *, '(a)' ) '  random values.'
      write ( *, '(a,i12)' )
     &  '  Using initial random number seed = ', seed
      write ( *, '(a)' ) ' '

      do test = 1, test_num

        x = r8_normal_01 ( seed )
        write ( *, '(2x,g14.6)' ) x

      end do

      return
      end
      subroutine test017 ( )

c*********************************************************************72
c
cc TEST017 tests R8_PI.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 August 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision four
      double precision one
      double precision r8_pi
      double precision v1
      double precision v2

      four = real ( 4, kind = 8 )
      one = real ( 1, kind = 8 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST017'
      write ( *, '(a)' ) '  R8_PI returns the value of PI.'
      write ( *, '(a)' ) ' '
      v1 = r8_pi ( )
      write ( *, '(a,g24.16)' ) '  R8_PI =     ', v1
      v2 = four * atan ( one )
      write ( *, '(a,g24.16)' ) '  4*atan(1) = ', v2

      return
      end
      subroutine test018 ( )

c*********************************************************************72
c
cc TEST018 tests R8_POWER.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 August 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision r8_power
      integer i
      integer p
      double precision r
      double precision value

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST018'
      write ( *, '(a)' ) '  R8_POWER computes R**P.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      R          P       R**P'
      write ( *, '(a)' ) ' '

      do i = -5, 5

        r = 2.0D+00
        p = i
        value = r8_power ( r, p )
        write ( *, '(2x,g14.6,i5,g14.6,i5)' ) r, p, value

      end do

      return
      end
      subroutine test019 ( )

c*********************************************************************72
c
cc TEST019 tests R8_POWER_FAST.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 August 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      integer mults
      integer p
      double precision r
      double precision rp

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST019'
      write ( *, '(a)' ) '  R8_POWER_FAST computes R**P, economizing on'
      write ( *, '(a)' ) '  multiplications.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      R          P       R**P       Mults'
      write ( *, '(a)' ) ' '

      do i = -10, 40

        r = 2.0D+00
        p = i
        call r8_power_fast ( r, p, rp, mults )
        write ( *, '(2x,g14.6,i5,g14.6,i5)' ) r, p, rp, mults

      end do

      return
      end
      subroutine test020 ( )

c*********************************************************************72
c
cc TEST020 tests R8_ROUND2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 September 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      integer nplace
      double precision r8_pi
      double precision x
      double precision xround

      x = r8_pi ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST020'
      write ( *, '(a)' ) '  R8_ROUND2 rounds a number to a'
      write ( *, '(a)' ) '  specified number of base 2 digits.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Test effect on PI:'
      write ( *, '(a,g24.16)' ) '  X = ', x
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  NPLACE  XROUND'
      write ( *, '(a)' ) ' '

      do i = 0, 20
        nplace = i
        call r8_round2 ( nplace, x, xround )
        write ( *, '(2x,i8,g24.16)' ) i, xround
      end do

      return
      end
      subroutine test021 ( )

c*********************************************************************72
c
cc TEST021 tests R8_ROUNDB.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 November 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer base
      integer i
      integer nplace
      double precision r8_pi
      double precision x
      double precision xround

      base = 3
      x = r8_pi ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST021'
      write ( *, '(a)' ) '  R8_ROUNDB rounds a number to a '
      write ( *, '(a)' ) '  specified number of base BASE digits.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Here, we will use BASE = ',base
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Test effect on PI:'
      write ( *, '(a,g24.16)' ) '  X = ', x
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  NPLACE  XROUND'
      write ( *, '(a)' ) ' '

      do i = 0, 20
        nplace = i
        call r8_roundb ( base, nplace, x, xround )
        write ( *, '(2x,i8,g24.16)' ) i, xround
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Try with a negative base:'
      x = 121.0D+00
      base = -3
      nplace = 3
      write ( *, '(a)' ) ' '
      write ( *, '(a,g24.16)' ) '  Input quantity is X = ', x
      write ( *, '(a,i8)' ) '  to be rounded in base ', base

      do nplace = 1, 5

        call r8_roundb ( base, nplace, x, xround )

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8,a,g24.16)' ) '  Output value to ', nplace,
     &    ' places is ', xround

      end do

      return
      end
      subroutine test022 ( )

c*********************************************************************72
c
cc TEST022 tests R8_ROUNDX.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 November 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      integer nplace
      double precision r8_pi
      double precision r8_uniform_01
      integer seed
      double precision x
      double precision xround

      seed = 123456789
      x = r8_pi ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST022'
      write ( *, '(a)' ) '  R8_ROUNDX rounds a number to a '
      write ( *, '(a)' ) '  specified number of decimal digits.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Test effect on PI:'
      write ( *, '(a,g24.16)' ) '  X = ', x
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  NPLACE  XROUND'
      write ( *, '(a)' ) ' '

      do i = 0, 10
        nplace = i
        call r8_roundx ( nplace, x, xround )
        write ( *, '(2x,i8,g24.16)' ) i, xround
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Test effect on random values:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  NPLACE  X     XROUND'
      write ( *, '(a)' ) ' '

      do i = 1, 5

        x = r8_uniform_01 ( seed )

        write ( *, '(a)' ) ' '

        do nplace = 0, 10, 2
          call r8_roundx ( nplace, x, xround )
          write ( *, '(2x,i8,2x,g24.16,2x,g24.16)' ) nplace, x, xround
        end do

      end do

      return
      end
      subroutine test023 ( )

c*********************************************************************72
c
cc TEST023 tests R8_SIGN.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 November 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 5 )

      double precision r8_sign
      integer test
      double precision x
      double precision x_test(test_num)

      save x_test

      data x_test /
     & -1.25D+00, -0.25D+00, 0.0D+00, +0.5D+00, +9.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST023'
      write ( *, '(a)' ) '  R8_SIGN returns the sign of a number.'
      write ( *, '(a)' ) ' '

      do test = 1, test_num
        x = x_test(test)
        write ( *, '(2x,f8.4,2x,f8.4)' ) x, r8_sign ( x )
      end do

      return
      end
      subroutine test0235 ( )

c*********************************************************************72
c
cc TEST0235 tests R8_SWAP.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 November 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision x
      double precision y

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0235'
      write ( *, '(a)' ) '  R8_SWAP swaps two reals.'

      x = 1.0D+00
      y = 3.141592653589793D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Before swapping:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '    X = ', x
      write ( *, '(a,g14.6)' ) '    Y = ', y

      call r8_swap ( x, y )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  After swapping:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '    X = ', x
      write ( *, '(a,g14.6)' ) '    Y = ', y

      return
      end
      subroutine test024 ( )

c*********************************************************************72
c
cc TEST024 tests R8_TO_R8_DISCRETE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 November 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision r8_uniform_ab
      integer ndx
      double precision r
      double precision rd
      double precision rhi
      double precision rhi2
      double precision rlo
      double precision rlo2
      integer seed
      integer test
      integer test_num

      ndx = 19
      rlo = 1.0D+00
      rhi = 10.0D+00
      test_num = 15

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST024'
      write ( *, '(a)' )
     &  '  R8_TO_R8_DISCRETE maps numbers to a discrete set'
      write ( *, '(a)' ) '  of equally spaced numbers in an interval.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of discrete values = ', ndx
      write ( *, '(a,2g14.6)' ) '  Real interval: ', rlo, rhi
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  R   RD'
      write ( *, '(a)' ) ' '

      seed = 123456789

      rlo2 = rlo - 2.0D+00
      rhi2 = rhi + 2.0D+00

      do test = 1, test_num
        r = r8_uniform_ab ( rlo2, rhi2, seed )
        call r8_to_r8_discrete ( r, rlo, rhi, ndx, rd )
        write ( *, '(2x,g14.6,g14.6)' ) r, rd
      end do

      return
      end
      subroutine test025 ( )

c*********************************************************************72
c
cc TEST025 tests R8_TO_I4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 November 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer ix
      integer ixmax
      integer ixmin
      double precision x
      double precision xmax
      double precision xmin

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST025'
      write ( *, '(a)' )
     &  '  R8_TO_I4 finds an integer IX in [IXMIN,IXMAX]'
      write ( *, '(a)' ) '  corresponding to X in [XMIN,XMAX].'

      xmin = 2.5D+00
      x = 3.5D+00
      xmax = 5.5D+00

      ixmin = 10
      ixmax = 40

      call r8_to_i4 ( x, xmin, xmax, ixmin, ixmax, ix )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,g14.6,a,g14.6)' )
     &  '  XMIN ',  xmin, '   X = ',  x, '  XMAX = ', xmax
      write ( *, '(a,i14,a,i14,a,i14)' )
     &  ' IXMIN ', ixmin, '  IX = ', ix, ' IXMAX = ', ixmax

      return
      end
      subroutine test026 ( )

c*********************************************************************72
c
cc TEST026 tests R8_UNIFORM_AB.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 June 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision b
      parameter ( b = 10.0D+00 )
      double precision c
      parameter ( c = 20.0D+00 )
      double precision r8_uniform_ab
      integer i
      integer j
      double precision r
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST026'
      write ( *, '(a)' ) 
     &  '  R8_UNIFORM_AB returns random values in a given range:'
      write ( *, '(a)' ) '  [ A, B ]'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  For this problem:'
      write ( *, '(a,g14.6)' ) '  B = ', b
      write ( *, '(a,g14.6)' ) '  C = ', c
      write ( *, '(a)' ) ' '

      seed = 123456789

      do i = 1, 10
        r = r8_uniform_ab ( b, c, seed )
        write ( *, '(2x,g14.6)' ) r
      end do

      return
      end
      subroutine test027 ( )

c*********************************************************************72
c
cc TEST027 tests R8_UNIFORM_01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 June 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      double precision r8_uniform_01
      integer seed
      integer seed_old
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST027'
      write ( *, '(a)' ) 
     &  '  R8_UNIFORM_01 produces a sequence of random values.'

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  Using random seed ', seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  SEED   R8_UNIFORM_01(SEED)'
      write ( *, '(a)' ) ' '
      do i = 1, 10
        seed_old = seed
        x = r8_uniform_01 ( seed )
        write ( *, '(2x,i12,2x,g14.6)' ) seed, x
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Verify that the sequence can be restarted.'
      write ( *, '(a)' ) 
     &  '  Set the seed back to its original value, and see that'
      write ( *, '(a)' ) '  we generate the same sequence.'

      seed = 123456789
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  SEED   R8_UNIFORM_01(SEED)'
      write ( *, '(a)' ) ' '

      do i = 1, 10
        seed_old = seed
        x = r8_uniform_01 ( seed )
        write ( *, '(2x,i12,2x,g14.6)' ) seed, x
      end do

      return
      end
      subroutine test028 ( )

c*********************************************************************72
c
cc TEST028 tests R8_UNIFORM_01
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 June 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      double precision mean
      integer n 
      parameter ( n = 1000 )
      integer seed
      double precision r8_uniform_01
      double precision t
      double precision variance
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST028'
      write ( *, '(a)' ) '  R8_UNIFORM_01 samples a uniform random'
      write ( *, '(a)' ) '  distribution in [0,1].'

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  Starting with seed = ', seed

      do i = 1, n
        x(i) = r8_uniform_01 ( seed )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  First few values:'
      write ( *, '(a)' ) ' '
      do i = 1, 5
        write ( *, '(2x,i8,2x,g14.6)' ) i, x(i)
      end do

      call r8vec_mean ( n, x, mean )

      call r8vec_variance ( n, x, variance )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of values computed was N = ', n
      write ( *, '(a,g14.6)' ) '  Average value was ', mean
      call r8vec_min ( n, x, t )
      write ( *, '(a,g14.6)' ) '  Minimum value was ', t
      call r8vec_max ( n, x, t )
      write ( *, '(a,g14.6)' ) '  Maximum value was ', t
      write ( *, '(a,g14.6)' ) '  Variance was ', variance

      return
      end
      subroutine test029 ( )

c*********************************************************************72
c
cc TEST029 tests R8_WALSH_1D;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 June 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      double precision r8_walsh_1d
      double precision w0
      double precision wm1
      double precision wm2
      double precision wm3
      double precision wp1
      double precision wp2
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST029'
      write ( *, '(a)' ) '  R8_WALSH_1D evaluates 1D Walsh functions:'
      write ( *, '(a)' ) ' '
      write ( *, * ) 'X  W(+2) W(+1) W(0) W(-1) W(-2) W(-3)'
      write ( *, '(a)' ) ' '

      do i = 0, 32

        x = real ( i, kind = 8 ) / 4.0D+00

        wp2 = r8_walsh_1d ( x,  2 )
        wp1 = r8_walsh_1d ( x,  1 )
        w0  = r8_walsh_1d ( x,  0 )
        wm1 = r8_walsh_1d ( x, -1 )
        wm2 = r8_walsh_1d ( x, -2 )
        wm3 = r8_walsh_1d ( x, -3 )

        write ( *, '(2x,f10.6,6f4.1)' ) x, wp2, wp1, w0, wm1, wm2, wm3

      end do

      return
      end
      subroutine test0295 ( )

c*****************************************************************************80
c
cc TEST0295 tests R8_WRAP;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 July 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      parameter ( a = - 2.0D+00 )
      double precision b
      parameter ( b = 12.0D+00 )
      double precision r
      double precision r2
      double precision r8_uniform_ab
      double precision r8_wrap
      double precision rhi
      parameter ( rhi = 6.5D+00 )
      double precision rlo
      parameter ( rlo = 3.0D+00 )
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 20 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0295'
      write ( *, '(a)' ) 
     &  '  R8_WRAP "wraps" an R8 to lie within an interval:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,g14.6)' ) 
     &  '  Wrapping interval is ', rlo, ', ', rhi
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      R      R8_WRAP ( R )'
      write ( *, '(a)' ) ' '
      seed = 123456789

      do test = 1, test_num

        r = r8_uniform_ab ( a, b, seed )
        r2 = r8_wrap ( r, rlo, rhi )
        write ( *, '(2x,g14.6,2x,g14.6)' ) r, r2

      end do

      return
      end
      subroutine test031 ( )

c*********************************************************************72
c
cc TEST031 tests R82POLY2_TYPE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 September 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 12 )

      double precision a
      double precision a_test(test_num)
      double precision b
      double precision b_test(test_num)
      double precision c
      double precision c_test(test_num)
      double precision d
      double precision d_test(test_num)
      double precision e
      double precision e_test(test_num)
      double precision f
      double precision f_test(test_num)
      integer test
      integer type

      save a_test
      save b_test
      save c_test
      save d_test
      save e_test
      save f_test

      data a_test /
     &  9.0D+00, 4.0D+00, 9.0D+00,  1.0D+00, 0.0D+00, 
     &  1.0D+00, 0.0D+00, 0.0D+00,  0.0D+00, 0.0D+00, 
     &  0.0D+00, 0.0D+00 /
      data b_test /
     &  -4.0D+00,   1.0D+00,  16.0D+00,   1.0D+00, 0.0D+00,  
     &   2.0D+00, 1.0D+00,   1.0D+00,  1.0D+00,  0.0D+00, 
     &   0.0D+00, 0.0D+00 /
      data c_test /
     &   0.0D+00,  -4.0D+00,   0.0D+00,   0.0D+00, 1.0D+00,  
     &   0.0D+00, 0.0D+00,   0.0D+00,  0.0D+00,  0.0D+00, 
     &   0.0D+00, 0.0D+00 /
      data d_test /
     &  -36.0D+00,  3.0D+00,  36.0D+00,  -6.0D+00, 3.0D+00, 
     &  -2.0D+00, 0.0D+00,   0.0D+00,  0.0D+00,  2.0D+00, 
     &   0.0D+00, 0.0D+00 /
      data e_test /
     &  -24.0D+00, -4.0D+00, -32.0D+00, -10.0D+00, -1.0D+00, 
     &   16.0D+00, -6.0D+00, -6.0D+00, -2.0D+00, -1.0D+00, 
     &   0.0D+00, 0.0D+00 /
      data f_test /
     &  -36.0D+00,  1.0D+00, -92.0D+00, 115.0D+00, -3.0D+00, 
     &   33.0D+00, +8.0D+00, 10.0D+00,  +1.0D+00,  1.0D+00, 
     &    0.0D+00, 1.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST031'
      write ( *, '(a)' ) '  R82POLY2_TYPE determines the type of a'
      write ( *, '(a)' ) '  second order equation in two variables.'
      write ( *, '(a)' ) ' '

      do test = 1, test_num

        a = a_test(test)
        b = b_test(test)
        c = c_test(test)
        d = d_test(test)
        e = e_test(test)
        f = f_test(test)

        write ( *, '(a)' ) ' '

        call r82poly2_print ( a, b, c, d, e, f )

        call r82poly2_type ( a, b, c, d, e, f, type )

        write ( *, '(a,i8)' ) '  Type = ', type

        call r82poly2_type_print ( type )

      end do

      return
      end
      subroutine test032 ( )

c*********************************************************************72
c
cc TEST032 tests R82VEC_ORDER_TYPE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integern
      parameter ( n = 4 )
      integer test_num
      parameter ( test_num = 10 )

      integer i
      integer j
      integer order
      integer seed
      integer test
      character * ( 40 ) title
      double precision x(2,n)

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST032'
      write ( *, '(a)' ) '  R82VEC_ORDER_TYPE classifies a R8VEC as'
      write ( *, '(a)' ) '  -1: no order'
      write ( *, '(a)' ) '   0: all equal;'
      write ( *, '(a)' ) '   1: ascending;'
      write ( *, '(a)' ) '   2: strictly ascending;'
      write ( *, '(a)' ) '   3: descending;'
      write ( *, '(a)' ) '   4: strictly descending.'
      write ( *, '(a)' ) ' '

      do test = 1, test_num

        call r8mat_uniform_01 ( 2, n, seed, x )

        do j = 1, n
          do i = 1, 2
            x(i,j) = dble ( nint ( 3.0D+00 * x(i,j) ) )
          end do
        end do

        call r82vec_order_type ( n, x, order )

        write ( title, '(a,i8)' ) '  Order type = ', order

        call r82vec_print ( n, x, title )

      end do

      return
      end
      subroutine test033 ( )

c*********************************************************************72
c
cc TEST033 tests R82VEC_PART_QUICK_A.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 June 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 12 )

      double precision a(2,n)
      double precision b
      double precision c
      integer l
      integer r
      integer seed

      b = 0.0D+00
      c = 10.0D+00
      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST033'
      write ( *, '(a)' ) '  R82VEC_PART_QUICK_A reorders an R82VEC'
      write ( *, '(a)' ) '  as part of a quick sort.'
      write ( *, '(a,i12)' ) 
     &  '  Using initial random number seed = ', seed

      call r8mat_uniform_ab ( 2, n, b, c, seed, a )

      call r82vec_print ( n, a, '  Before rearrangement:' )

      call r82vec_part_quick_a ( n, a, l, r )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Rearranged array'
      write ( *, '(a,i8)' ) '  Left index =  ', l
      write ( *, '(a,i8)' ) '  Key index =   ', l+1
      write ( *, '(a,i8)' ) '  Right index = ', r
      write ( *, '(a)' ) ' '

      call r82vec_print ( l,     a(1:2,1:l),   '  Left half:' )
      call r82vec_print ( 1,     a(1:2,l+1),   '  Key:' )
      call r82vec_print ( n-l-1, a(1:2,l+2:n), '  Right half:' )

      return
      end
      subroutine test034 ( )

c*********************************************************************72
c
cc TEST034 tests R82VEC_SORT_HEAP_INDEX_A.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 June 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 12 )

      double precision a(2,n)
      double precision b
      double precision c
      integer i
      integer indx(n)
      integer seed

      b = 0.0D+00
      c = 10.0D+00
      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST034'
      write ( *, '(a)' ) 
     &  '  R82VEC_SORT_HEAP_INDEX_A index sorts an R82VEC'
      write ( *, '(a)' ) '  using heapsort.'
      write ( *, '(a,i12)' ) 
     &  '  Using initial random number seed = ', seed

      call r8mat_uniform_ab ( 2, n, b, c, seed, a )
c
c  Give a few elements the same first component.
c
      a(1,3) = a(1,5)
      a(1,4) = a(1,12)
c
c  Give a few elements the same second component.
c
      a(2,6) = a(2,1)
      a(2,2) = a(2,9)
c
c  Make two entries equal.
c
      a(1:2,7) = a(1:2,11)

      call r82vec_print ( n, a, '  Before rearrangement:' )

      call r82vec_sort_heap_index_a ( n, a, indx )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I Index A(Index)'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,i8,i8,g14.6,g14.6)' ) i, indx(i), a(1:2,indx(i))
      end do

      call r82vec_permute ( n, indx, a )

      call r82vec_print ( n, a, 
     &  '  After rearrangement by R82VEC_PERMUTE:' )

      return
      end
      subroutine test035 ( )

c*********************************************************************72
c
cc TEST035 tests R82VEC_SORT_QUICK_A.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 June 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 12 )

      double precision a(2,n)
      double precision b
      double precision c
      integer seed

      b = 0.0D+00
      c = 10.0D+00
      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST035'
      write ( *, '(a)' ) '  R82VEC_SORT_QUICK_A sorts an R82VEC'
      write ( *, '(a)' ) '  using quick sort.'
      write ( *, '(a,i12)' ) 
     &  '  Using initial random number seed = ', seed

      call r8mat_uniform_ab ( 2, n, b, c, seed, a )
c
c  Give a few elements the same first component.
c
      a(1,3) = a(1,5)
      a(1,4) = a(1,12)
c
c  Give a few elements the same second component.
c
      a(2,6) = a(2,1)
      a(2,2) = a(2,9)
c
c  Make two entries equal.
c
      a(1:2,7) = a(1:2,11)

      call r82vec_print ( n, a, '  Before rearrangement:' )

      call r82vec_sort_quick_a ( n, a )

      call r82vec_print ( n, a, '  Sorted array:' )

      return
      end
      subroutine test036 ( )

c*********************************************************************72
c
cc TEST036 tests R8BLOCK_EXPAND_LINEAR.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 June 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer l
      parameter ( l = 4 )
      integer m
      parameter ( m = 3 )
      integer n
      parameter ( n = 2 )
      integer lfat
      parameter ( lfat = 1 )
      integer mfat
      parameter ( mfat = 2 )
      integer nfat
      parameter ( nfat = 1 )
      integer l2
      parameter ( l2 = ( l - 1 ) * ( lfat + 1 ) + 1 )
      integer m2
      parameter ( m2 = ( m - 1 ) * ( mfat + 1 ) + 1 )
      integer n2
      parameter ( n2 = ( n - 1 ) * ( nfat + 1 ) + 1 )

      double precision x(l,m,n)
      double precision xfat(l2,m2,n2)

      save x

      data x /
     &  1.0D+00,  2.0D+00,  3.0D+00,   4.0D+00,  1.0D+00, 
     &  4.0D+00,  9.0D+00, 16.0D+00,   1.0D+00,  8.0D+00, 
     & 27.0D+00, 64.0D+00,  2.0D+00,   4.0D+00,  6.0D+00, 
     &  8.0D+00,  2.0D+00,  8.0D+00,  18.0D+00, 32.0D+00, 
     &  2.0D+00, 16.0D+00, 54.0D+00, 128.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST036'
      write ( *, '(a)' ) 
     &  '  R8BLOCK_EXPAND_LINEAR linearly interpolates new data'
      write ( *, '(a)' ) '  between old values in a 3D block.'

      call r8block_print ( l, m, n, x, '  Original block:' )
 
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  LFAT = ', lfat
      write ( *, '(a,i8)' ) '  MFAT = ', mfat
      write ( *, '(a,i8)' ) '  NFAT = ', nfat

      call r8block_expand_linear ( l, m, n, x, lfat, mfat, nfat, xfat )

      call r8block_print ( l2, m2, n2, xfat, '  Fattened block:' )

      return
      end
      subroutine test0365 ( )

c*****************************************************************************80
c
cc TEST0365 tests R8BLOCK_PRINT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 June 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer l
      parameter ( l = 4 )
      integer m
      parameter ( m = 3 )
      integer n
      parameter ( n = 2 )

      double precision x(l,m,n)

      save x

      data x /
     &  1.0D+00,  2.0D+00,  3.0D+00,   4.0D+00,  1.0D+00, 
     &  4.0D+00,  9.0D+00, 16.0D+00,   1.0D+00,  8.0D+00, 
     & 27.0D+00, 64.0D+00,  2.0D+00,   4.0D+00,  6.0D+00, 
     &  8.0D+00,  2.0D+00,  8.0D+00,  18.0D+00, 32.0D+00, 
     &  2.0D+00, 16.0D+00, 54.0D+00, 128.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0365'
      write ( *, '(a)' ) '  R8BLOCK_PRINT prints an R8BLOCK.'

      call r8block_print ( l, m, n, x, '  The 3D array:' )

      return
      end
      subroutine test037 ( )

c*********************************************************************72
c
cc TEST037 tests R8COL_FIND.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 June 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )
      integer m
      parameter ( m = 3 )

      double precision dtab(m,n)
      integer i
      integer icol
      integer j
      integer k
      double precision r8vec(m)

      k = 1
      do i = 1, m
        do j = 1, n

          dtab(i,j) = dble ( k )

          if ( j .eq. 3 ) then
            r8vec(i) = dble ( k )
          end if

          k = k + 1

        end do
      end do

      call r8col_find ( m, n, dtab, r8vec, icol )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST037'
      write ( *, '(a)' ) '  For an R8COL;'
      write ( *, '(a)' ) 
     &  '  R8COL_FIND seeks a column matching given data.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  R8COL_FIND returns ICOL = ', icol

      return
      end
      subroutine test038 ( )

c*********************************************************************72
c
cc TEST038 tests R8COL_INSERT and R8COL_SORT_HEAP_A.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 June 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 3 )
      integer n_max
      parameter ( n_max = 10 )

      double precision a(m,n_max)
      integer col
      double precision r8vec1(m)
      double precision r8vec2(m)
      integer n

      save a
      save r8vec1
      save r8vec2

      data a /
     &  2.0D+00,  6.0D+00, 10.0D+00, 
     &  4.0D+00,  8.0D+00, 12.0D+00, 
     &  1.0D+00,  5.0D+00,  9.0D+00, 
     &  3.0D+00,  7.0D+00, 11.0D+00, 
     &  0.0D+00,  0.0D+00,  0.0D+00, 
     &  0.0D+00,  0.0D+00,  0.0D+00, 
     &  0.0D+00,  0.0D+00,  0.0D+00, 
     &  0.0D+00,  0.0D+00,  0.0D+00, 
     &  0.0D+00,  0.0D+00,  0.0D+00, 
     &  0.0D+00,  0.0D+00,  0.0D+00 /
      data r8vec1 / 3.0D+00, 7.0D+00, 11.0D+00 /
      data r8vec2 / 3.0D+00, 4.0D+00, 18.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST038'
      write ( *, '(a)' ) '  For an R8COL;'
      write ( *, '(a)' ) 
     &  '  R8COL_SORT_HEAP_A does an ascending heap sort'
      write ( *, '(a)' ) '  R8COL_INSERT inserts new columns.'

      n = 4

      call r8mat_print ( m, n, a, '  The unsorted matrix:' )

      call r8col_sort_heap_a ( m, n, a )

      call r8mat_print ( m, n, a, '  The sorted matrix:' )

      call r8vec_print ( m, r8vec1, '  New column:' )

      call r8col_insert ( n_max, m, n, a, r8vec1, col )

      if ( col .lt. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) 
     &    '  The data was already in column ', abs ( col )
      else
        call r8mat_print ( m, n, a, '  The updated matrix:' )
      end if

      call r8vec_print ( m, r8vec2, '  New column:' )

      call r8col_insert ( n_max, m, n, a, r8vec2, col )

      if ( col .lt. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) 
     &    '  The data was already in column ', abs ( col )
      else
        call r8mat_print ( m, n, a, '  The updated matrix:' )
      end if

      return
      end
      subroutine test0383 ( )

c*********************************************************************72
c
cc TEST0383 tests R8COL_PART_QUICK_A.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 June 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 2 )
      integer n
      parameter ( n = 8 )

      double precision a(m,n)
      integer l
      integer r

      save a

      data a /
     &   2.0D+00, 4.0D+00, 
     &   8.0D+00, 8.0D+00, 
     &   6.0D+00, 2.0D+00, 
     &   0.0D+00, 2.0D+00, 
     &  10.0D+00, 6.0D+00, 
     &  10.0D+00, 0.0D+00, 
     &   0.0D+00, 6.0D+00, 
     &   5.0D+00, 8.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0383'
      write ( *, '(a)' ) '  For an R8COL;'
      write ( *, '(a)' ) '  R8COL_PART_QUICK_A partitions the matrix.'

      call r8mat_print ( m, n, a, '  The matrix:' )

      l = 2
      r = 4
      call r8col_part_quick_a ( m, n, a, l, r )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  L = ', l
      write ( *, '(a,i4)' ) '  R = ', r

      call r8mat_print ( m, n, a, '  The partitioned matrix:' )

      return
      end
      subroutine test0385 ( )

c*********************************************************************72
c
cc TEST0385 tests R8COL_SORT_HEAP_INDEX_A.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 October 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 3 )
      integer n
      parameter ( n = 15 )

      double precision a(m,n)
      integer indx(n)
      integer j
      integer j2

      save a

      data a /
     &  2.0D+00,  6.0D+00, 10.0D+00, 
     &  4.0D+00,  8.0D+00, 12.0D+00, 
     &  1.0D+00,  5.0D+00,  9.0D+00, 
     &  3.0D+00,  7.0D+00, 11.0D+00, 
     &  2.0D+00,  6.0D+00,  0.0D+00, 
     &  3.0D+00,  4.0D+00, 18.0D+00, 
     &  0.0D+00,  0.0D+00,  0.0D+00, 
     &  0.0D+00,  6.0D+00, 10.0D+00, 
     &  2.0D+00,  6.0D+00, 10.0D+00, 
     &  3.0D+00,  7.0D+00, 11.0D+00, 
     &  2.0D+00,  0.0D+00, 10.0D+00, 
     &  2.0D+00,  6.0D+00, 10.0D+00, 
     &  1.0D+00,  5.0D+00,  9.0D+00, 
     &  1.0D+00,  5.0D+00,  9.1D+00, 
     &  1.0D+00,  5.1D+00,  9.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0385'
      write ( *, '(a)' ) 
     &  '  R8COL_SORT_HEAP_INDEX_A computes an index vector which'
      write ( *, '(a)' ) '  ascending sorts an R8COL.'

      call r8mat_transpose_print ( m, n, a, 
     &  '  The unsorted R8COL (transposed):' )

      call r8col_sort_heap_index_a ( m, n, a, indx )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The implicitly sorted R8COL (transposed)'
      write ( *, '(a)' ) ' '

      do j = 1, n
        j2 = indx(j)
        write ( *, '(2x,i4,a,2x,f10.1,2x,f10.1,2x,f10.1)' ) 
     &    j2, ':', a(1:m,j2)
      end do

      return
      end
      subroutine test039 ( )

c*********************************************************************72
c
cc TEST039 tests R8COL_SORT_QUICK_A.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 October 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 3 )
      integer n
      parameter ( n = 10 )

      double precision a(m,n)
      double precision b
      parameter ( b = 0.0D+00 )
      double precision c
      parameter ( c = 10.0D+00 )
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST039'
      write ( *, '(a)' ) '  For an R8COL;'
      write ( *, '(a)' ) '  R8COL_SORT_QUICK_A does a quicksort.'

      seed = 123456789

      call r8mat_uniform_ab ( m, n, b, c, seed, a )

      call r8mat_print ( m, n, a, '  The unsorted matrix:' )

      call r8col_sort_quick_a ( m, n, a )

      call r8mat_print ( m, n, a, '  The sorted matrix:' )

      return
      end
      subroutine test0393 ( )

c*********************************************************************72
c
cc TEST0393 tests R8COL_SORTED_TOL_UNIQUE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 October 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 3 )
      integer n
      parameter ( n = 22 )

      double precisiona(m,n)
      double precision r8_epsilon
      double precision tol
      integer unique_num

      save a

      data a / 
     &  1.9D+00,  0.0D+00, 10.0D+00, 
     &  2.0D+00,  6.0D+00, 10.0D+00, 
     &  4.0D+00,  8.0D+00, 12.0D+00, 
     &  1.0D+00,  5.0D+00,  9.0D+00, 
     &  3.0D+00,  7.0D+00, 11.0D+00, 
     &  2.0D+00,  6.0D+00,  0.0D+00, 
     &  2.0D+00,  0.0D+00, 10.1D+00, 
     &  2.0D+00,  0.1D+00, 10.0D+00, 
     &  3.0D+00,  4.0D+00, 18.0D+00, 
     &  1.9D+00,  8.0D+00, 10.0D+00, 
     &  0.0D+00,  0.0D+00,  0.0D+00, 
     &  0.0D+00,  6.0D+00, 10.0D+00, 
     &  2.1D+00,  0.0D+00, 10.0D+00, 
     &  2.0D+00,  6.0D+00, 10.0D+00, 
     &  3.0D+00,  7.0D+00, 11.0D+00, 
     &  2.0D+00,  0.0D+00, 10.0D+00, 
     &  2.0D+00,  0.0D+00, 10.0D+00, 
     &  2.0D+00,  6.0D+00, 10.0D+00, 
     &  1.0D+00,  5.0D+00,  9.0D+00, 
     &  2.0D+00,  0.0D+00, 10.1D+00, 
     &  1.0D+00,  5.0D+00,  9.1D+00, 
     &  1.0D+00,  5.1D+00,  9.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0393'
      write ( *, '(a)' ) 
     &  '  R8COL_SORTED_TOL_UNIQUE finds tolerably unique columns'
      write ( *, '(a)' ) '  in a sorted R8COL.'

      call r8mat_transpose_print ( m, n, a, 
     &  '  The unsorted R8COL (transposed):' )

      call r8col_sort_heap_a ( m, n, a )

      call r8mat_transpose_print ( m, n, a, 
     &  '  The sorted R8COL (transposed):' )

      tol = 0.25D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Using tolerance = ', tol

      call r8col_sorted_tol_unique ( m, n, a, tol, unique_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  Number of tolerably unique columns is ', unique_num

      call r8mat_transpose_print ( m, unique_num, a, 
     &  '  The sorted tolerably unique R8COL (transposed):' )

      return
      end
      subroutine test0395 ( )

c*********************************************************************72
c
cc TEST0395 tests R8COL_SORTED_TOL_UNIQUE_COUNT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 October 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 3 )
      integer n
      parameter ( n = 22 )

      double precision a(m,n)
      integer n2
      double precision r8_epsilon
      double precision tol
      integer unique_num

      save a

      data a /
     &  1.9D+00,  0.0D+00, 10.0D+00, 
     &  2.0D+00,  6.0D+00, 10.0D+00, 
     &  4.0D+00,  8.0D+00, 12.0D+00, 
     &  1.0D+00,  5.0D+00,  9.0D+00,
     &  3.0D+00,  7.0D+00, 11.0D+00,
     &  2.0D+00,  6.0D+00,  0.0D+00, 
     &  2.0D+00,  0.0D+00, 10.1D+00, 
     &  2.0D+00,  0.1D+00, 10.0D+00, 
     &  3.0D+00,  4.0D+00, 18.0D+00, 
     &  1.9D+00,  8.0D+00, 10.0D+00, 
     &  0.0D+00,  0.0D+00,  0.0D+00, 
     &  0.0D+00,  6.0D+00, 10.0D+00, 
     &  2.1D+00,  0.0D+00, 10.0D+00, 
     &  2.0D+00,  6.0D+00, 10.0D+00, 
     &  3.0D+00,  7.0D+00, 11.0D+00, 
     &  2.0D+00,  0.0D+00, 10.0D+00, 
     &  2.0D+00,  0.0D+00, 10.0D+00,
     &  2.0D+00,  6.0D+00, 10.0D+00, 
     &  1.0D+00,  5.0D+00,  9.0D+00, 
     &  2.0D+00,  0.0D+00, 10.1D+00, 
     &  1.0D+00,  5.0D+00,  9.1D+00, 
     &  1.0D+00,  5.1D+00,  9.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0395'
      write ( *, '(a)' ) 
     &  '  R8COL_SORTED_TOL_UNIQUE_COUNT counts tolerably '
      write ( *, '(a)' ) '  unique columns in a sorted R8COL.'

      call r8mat_transpose_print ( m, n, a, 
     &  '  The unsorted R8COL (transposed):' )

      call r8col_sort_heap_a ( m, n, a )

      call r8mat_transpose_print ( m, n, a, 
     &  '  The sorted R8COL (transposed):' )

      tol = 0.25D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Using tolerance = ', tol

      call r8col_sorted_tol_unique_count ( m, n, a, tol, unique_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  Number of tolerably unique columns is ', unique_num

      return
      end
      subroutine test0397 ( )

c*********************************************************************72
c
cc TEST0397 tests R8COL_SORTED_TOL_UNDEX.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 October 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 3 )
      integer n
      parameter ( n = 22 )

      double precision a(m,n)
      double precision au(m,n)
      integer i
      integer j
      integer j2
      integer n_unique
      double precision r8_epsilon
      double precision tol
      integer undx(n)
      integer xdnu(n)

      save a

      data a /
     &  1.9D+00,  0.0D+00, 10.0D+00, 
     &  2.0D+00,  6.0D+00, 10.0D+00, 
     &  4.0D+00,  8.0D+00, 12.0D+00, 
     &  1.0D+00,  5.0D+00,  9.0D+00, 
     &  3.0D+00,  7.0D+00, 11.0D+00, 
     &  2.0D+00,  6.0D+00,  0.0D+00, 
     &  2.0D+00,  0.0D+00, 10.1D+00, 
     &  2.0D+00,  0.1D+00, 10.0D+00, 
     &  3.0D+00,  4.0D+00, 18.0D+00, 
     &  1.9D+00,  8.0D+00, 10.0D+00, 
     &  0.0D+00,  0.0D+00,  0.0D+00, 
     &  0.0D+00,  6.0D+00, 10.0D+00, 
     &  2.1D+00,  0.0D+00, 10.0D+00, 
     &  2.0D+00,  6.0D+00, 10.0D+00, 
     &  3.0D+00,  7.0D+00, 11.0D+00, 
     &  2.0D+00,  0.0D+00, 10.0D+00, 
     &  2.0D+00,  0.0D+00, 10.0D+00, 
     &  2.0D+00,  6.0D+00, 10.0D+00, 
     &  1.0D+00,  5.0D+00,  9.0D+00, 
     &  2.0D+00,  0.0D+00, 10.1D+00, 
     &  1.0D+00,  5.0D+00,  9.1D+00, 
     &  1.0D+00,  5.1D+00,  9.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0397'
      write ( *, '(a)' ) 
     &  '  R8COL_SORTED_TOL_UNDEX produces index vectors which '
      write ( *, '(a)' ) 
     &  '  create a sorted list of the tolerably unique columns'
      write ( *, '(a)' ) 
     &  '  of a sorted R8COL,'
      write ( *, '(a)' ) 
     &  '  and a map from the original R8COL to the (implicit)'
      write ( *, '(a)' ) 
     &  '  R8COL of sorted tolerably unique elements.'

      call r8mat_transpose_print ( m, n, a, 
     &  '  The unsorted R8COL (transposed):' )

      call r8col_sort_heap_a ( m, n, a )

      call r8mat_transpose_print ( m, n, a, 
     &  '  The sorted R8COL (transposed):' )

      tol = 0.25D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Tolerance for equality = ', tol

      call r8col_sorted_tol_unique_count ( m, n, a, tol, n_unique )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  Number of tolerably unique columns is ', n_unique

      call r8col_sorted_tol_undex ( m, n, a, n_unique, tol, undx, xdnu )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  XDNU points to the representative for each item.'
      write ( *, '(a)' ) '  UNDX selects the representatives.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I  XDNU  UNDX'
      write ( *, '(a)' ) ' '
      do i = 1, n_unique
        write ( *, '(2x,i4,2x,i4,2x,i4)' ) i, xdnu(i), undx(i)
      end do
      do i = n_unique + 1, n
        write ( *, '(2x,i4,2x,i4)'       ) i, xdnu(i)
      end do

      do j = 1, n_unique
        au(1:m,j) = a(1:m,undx(j))
      end do

      call r8mat_transpose_print ( m, n_unique, au, 
     &  '  The tolerably unique R8COL (transposed):' )

      return
      end
      subroutine test049 ( )

c*********************************************************************72
c
cc TEST049 tests R8MAT_CHOLESKY_FACTOR, R8MAT_CHORESKY_FACTOR and R8MAT_CHOLESKY_SOLVE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 April 2009
c
c  Author:
c
c    John Burkardt
c
c     implicit none
c
      integer n
      parameter ( n = 5 )

      double precision a(n,n)
      double precision b(n)
      double precision d(n,n)
      integer i
      integer ierror
      integer j
      integer k
      double precision l(n,n)
      double precision r(n,n)
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST049'
      write ( *, '(a)' ) '  For a positive definite symmetric matrix,'
      write ( *, '(a)' ) '  R8MAT_CHOLESKY_FACTOR computes the lower'
      write ( *, '(a)' ) '  triangular Cholesky factor.'
      write ( *, '(a)' ) '  R8MAT_CHORESKY_FACTOR computes the upper'
      write ( *, '(a)' ) '  triangular Cholesky factor.'
      write ( *, '(a)' ) '  R8MAT_CHOLESKY_SOLVE solves a linear system'
      write ( *, '(a)' ) '  using the Cholesky factorization.'

      do i = 1, n
        do j = 1, n
          if ( i .eq. j ) then
            a(i,j) = 2.0D+00
          else if ( abs ( i - j ) .eq. 1 ) then
            a(i,j) = -1.0D+00
          else
            a(i,j) = 0.0D+00
          end if
        end do
      end do

      call r8mat_print ( n, n, a, '  Matrix to be factored:' )
c
c  Compute the lower Cholesky factor.
c
      call r8mat_cholesky_factor ( n, a, l, ierror )

      call r8mat_print ( n, n, l, '  Cholesky factor L:' )

      do i = 1, n
        do j = 1, n
          d(i,j) = 0.0D+00
          do k = 1, n
            d(i,j) = d(i,j) + l(i,k) * l(j,k)
          end do
        end do
      end do

      call r8mat_print ( n, n, d, '  Product L * L'':' )
c
c  Compute the upper Cholesky factor.
c
      call r8mat_choresky_factor ( n, a, r, ierror )

      call r8mat_print ( n, n, r, '  Cholesky factor R:' )

      do i = 1, n
        do j = 1, n
          d(i,j) = 0.0D+00
          do k = 1, n
            d(i,j) = d(i,j) + r(i,k) * r(j,k)
          end do
        end do
      end do

      call r8mat_print ( n, n, d, '  Product R * R'':' )
c
c  Solve a linear system.
c
      do i = 1, n - 1
        b(i) = 0.0D+00
      end do
      b(n) = dble ( n + 1 )

      call r8vec_print ( n, b, '  Right hand side:' )

      call r8mat_cholesky_solve ( n, l, b, x )

      call r8vec_print ( n, x, '  Computed solution:' )

      return
      end
      subroutine test0555 ( )

c*********************************************************************72
c
cc TEST0555 tests R8MAT_FSS.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )
      integer nb
      parameter ( nb = 3 )

      double precision a(n,n)
      double precision b(n,nb)
      integer i
      integer info
      integer seed
      double precision x(n)

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0555'
      write ( *, '(a)' ) '  For a matrix in general storage,'
      write ( *, '(a)' ) 
     &  '  R8MAT_FSS factors and solves multiple linear system.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix order N = ', n
c
c  Set the matrix.
c
      call r8mat_uniform_01 ( n, n, seed, a )
c
c  Set the desired solutions.
c
      do i = 1, n
        x(i) = 1.0D+00
      end do
      call r8mat_mv ( n, n, a, x, b(1,1) )

      do i = 1, n
        x(i) = dble ( i )
      end do
      call r8mat_mv ( n, n, a, x, b(1,2) )

      do i = 1, n
        x(i) = dble ( mod ( i - 1, 3 ) + 1 )
      end do
      call r8mat_mv ( n, n, a, x, b(1,3) )
c
c  Factor and solve the system.
c
      call r8mat_fss ( n, a, nb, b, info )
      
      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST0555 - Fatal error!'
        write ( *, '(a)' ) '  R8MAT_FSS reports the matrix is singular.'
        return
      end if

      call r8mat_print ( n, nb, b, '  Solutions:' )

      return
      end
      subroutine test1143 ( )

c*********************************************************************72
c
cc TEST1143 tests R8VEC_BRACKET5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 October 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )
      integer test_num
      parameter ( test_num = 6 )

      integer left
      integer r8vec_bracket5
      integer right
      integer test
      double precision x(n)
      double precision xtest(test_num)
      double precision xval

      save xtest

      data xtest /
     &  -10.0D+00, 1.0D+00, 4.5D+00, 5.0D+00, 10.0D+00, 12.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1143'
      write ( *, '(a)' ) '  R8VEC_BRACKET5 finds a pair of entries in a'
      write ( *, '(a)' ) '  sorted R8VEC which bracket a value.'

      call r8vec_indicator ( n, x )
      x(6) = x(5)

      call r8vec_print ( n, x, '  Sorted array:' )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        LEFT                   RIGHT'
      write ( *, '(a)' ) '      X(LEFT)       XVAL     X(RIGHT)'
      write ( *, '(a)' ) ' ' 

      do test = 1, test_num

        xval = xtest(test)

        left = r8vec_bracket5 ( n, x, xval )

        if ( left .eq. -1 ) then
          write ( *, '(2x,i10)' ) left
          write ( *, '(2x,10x,2x,f10.4,2x,a)' ) xval, '(Not bracketed!)'
        else
          right = left + 1
          write ( *, '(2x,i10,2x,10x,2x,i10)' ) left, right
          write ( *, '(2x,f10.4,2x,f10.4,2x,f10.4)' ) 
     &      x(left), xval, x(right)
        end if

      end do

      return
      end
      subroutine test1145 ( )

c*********************************************************************72
c
cc TEST1145 tests R8VEC_CHEBYSPACE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 June 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 10 )

      integer n
      double precision r(n_max)
      double precision r1
      double precision r2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1145'
      write ( *, '(a)' ) 
     &  '  R8VEC_CHEBYSPACE computes N Chebyshev points in [R1,R2].'

      r1 = -1.0D+00
      r2 = +1.0D+00
      n = 5

      call r8vec_chebyspace ( n, r1, r2, r )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i4,a,g14.6,a,g14.6)' ) 
     &  '  N = ', n, '  R1 = ', r1, '  R2 = ', r2

      call r8vec_print ( n, r, '  Chebyshev points:' )

      r1 =   0.0D+00
      r2 = +10.0D+00
      n = 7

      call r8vec_chebyspace ( n, r1, r2, r )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i4,a,g14.6,a,g14.6)' ) 
     &  '  N = ', n, '  R1 = ', r1, '  R2 = ', r2

      call r8vec_print ( n, r, '  Chebyshev points:' )

      return
      end
      subroutine test1147 ( )

c*********************************************************************72
c
cc TEST1147 tests R8VEC_CONVOLUTION
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 May 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      integer n

      parameter ( m = 4 )
      parameter ( n = 3 )

      double precision x(m)
      double precision y(n)
      double precision z(m+n-1)
      double precision z_correct(m+n-1)

      save x
      save y
      save z_correct

      data x / 1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00 /
      data y / -1.0D+00, 5.0D+00, 3.0D+00 /
      data z_correct  / -1.0D+00, 3.0D+00, 10.0D+00, 17.0D+00, 
     &  29.0D+00, 12.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1147'
      write ( *, '(a)' ) '  R8VEC_CONVOLUTION computes the convolution'
      write ( *, '(a)' ) '  of two vectors.'

      call r8vec_print ( m, x, '  The factor X:' )
      call r8vec_print ( n, y, '  The factor Y:' )

      call r8vec_convolution ( m, x, n, y, z )

      call r8vec_print ( m + n - 1, z, 
     &  '  The convolution z = x star y:' )

      call r8vec_print ( m + n - 1, z_correct, '  Correct answer:' )

      return
      end
      subroutine test115 ( )

c*********************************************************************72
c
cc TEST115 tests R8VEC_CONVOLUTION_CIRC
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 August 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      double precision x(n)
      double precision y(n)
      double precision z(n)
      double precision z_correct(n)

      save x
      save y
      save z_correct

      data x / 1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00 /
      data y / 1.0D+00, 2.0D+00, 4.0D+00, 8.0D+00 /
      data z_correct / 37.0D+00, 44.0D+00, 43.0D+00, 26.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST115'
      write ( *, '(a)' ) 
     &  '  R8VEC_CONVOLUTION_CIRC computes the circular convolution'
      write ( *, '(a)' ) '  of two vectors.'

      call r8vec_print ( n, x, '  The factor X:' )
      call r8vec_print ( n, y, '  The factor Y:' )

      call r8vec_convolution_circ ( n, x, y, z )

      call r8vec_print ( n, z, 
     &  '  The circular convolution z = x CC y:' )

      call r8vec_print ( n, z_correct, '  Correct answer:' )

      return
      end
      subroutine test1251 ( )

c*********************************************************************72
c
cc TEST1251 tests R8VEC_INDEX_SORTED_RANGE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 October 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 20 )

      integer i
      integer i_hi
      integer i_lo
      integer indx(n)
      double precision r(n)
      double precision r_lo
      double precision r_hi
      double precision r8_uniform_01
      integer seed
      double precision t
      integer test

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1251'
      write ( *, '(a)' )
     &  '  R8VEC_INDEX_SORTED_RANGE seeks the range I_LO:I_HI'
      write ( *, '(a)' ) '  of entries of sorted indexed R so that'
      write ( *, '(a)' )
     &  '  R_LO <= R(INDX(I)) <= R_HI for I_LO <= I <= I_HI.'

      seed = 123456789

      do test = 1, 5

        call r8vec_uniform_01 ( n, seed, r )

        call r8vec_print ( n, r, '  Array' )

        call r8vec_sort_heap_index_a ( n, r, indx )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '     I  INDX    R(INDX(I))'
        write ( *, '(a)' ) ' '
        do i = 1, n
          write ( *, '(2x,i4,2x,i4,2x,g14.6)' ) i, indx(i), r(indx(i))
        end do

        r_lo = r8_uniform_01 ( seed )
        r_hi = r8_uniform_01 ( seed )

        if ( r_hi .lt. r_lo ) then
          t = r_lo
          r_lo = r_hi
          r_hi = t
        end if

        call r8vec_index_sorted_range ( n, r, indx, r_lo, r_hi, i_lo,
     &    i_hi )

        write ( *, '(a)' ) ' '
        if ( i_hi .lt. i_lo ) then
          write ( *, '(2x,a4,2x,6x,g14.6)' ) 'R_LO', r_lo
          write ( *, '(2x,a4,2x,6x,g14.6)' ) 'R_HI', r_hi
          write ( *, '(a)' ) '  Empty range in R.'
        else

          write ( *, '(2x,a4,2x,6x,g14.6)' ) 'R_LO', r_lo
          do i = i_lo, i_hi
            write ( *, '(2x,i4,2x,i4,2x,g14.6)' ) i, indx(i), r(indx(i))
          end do
          write ( *, '(2x,a4,2x,6x,g14.6)' ) 'R_HI', r_hi
        end if

      end do

      return
      end
      subroutine test1252 ( )

c*********************************************************************72
c
cc TEST1252 tests R8VEC_INDEXED_HEAP_D;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 August 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 20 )
      integer n
      parameter ( n = 10 )

      double precision a(m)
      integer i
      integer indx(n)

      save a
      save indx

      data a /
     &  101.0D+00, 102.0D+00, 103.0D+00, 104.0D+00, 105.0D+00,
     &  106.0D+00, 107.0D+00, 108.0D+00, 109.0D+00, 110.0D+00,
     &  111.0D+00, 112.0D+00, 113.0D+00, 114.0D+00, 115.0D+00,
     &  116.0D+00, 117.0D+00, 118.0D+00, 119.0D+00, 120.0D+00 /
      data indx /
     &  1, 11, 17, 5, 7, 13, 15, 3, 19, 9 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1252'
      write ( *, '(a)' ) '  R8VEC_INDEXED_HEAP_D creates a descending'
      write ( *, '(a)' ) '  heap from an indexed R8VEC.'
c
c  Print before.
c
      call r8vec_print ( m, a, '  The data vector:' )
      call i4vec_print ( n, indx, '  The index vector:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  A(INDX):'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i4,2x,g14.6)' ) i, a(indx(i))
      end do
c
c  Create the heap.
c
      call r8vec_indexed_heap_d ( n, a, indx )
c
c  Print afterwards.  Only INDX should change.
c
      call r8vec_print ( m, a,
     &  '  The data vector (should NOT change):' )
      call i4vec_print ( n, indx, '  The index vector (may change):' )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  A(INDX) is now a heap:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i4,2x,g14.6)' ) i, a(indx(i))
      end do

      return
      end
      subroutine test1255 ( )

c*********************************************************************72
c
cc TEST1255 tests R8VEC_INDEXED_HEAP_D_EXTRACT and related routines.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 August 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 20 )
      integer n_max
      parameter ( n_max = 20 )

      double precision a(m)
      integer i
      integer indx(n_max)
      integer indx_extract
      integer indx_insert
      integer indx_max
      integer n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1255'
      write ( *, '(a)' ) '  For an indexed R8VEC,'
      write ( *, '(a)' )
     &  '  R8VEC_INDEXED_HEAP_D_INSERT inserts a value into the heap.'
      write ( *, '(a)' )
     &  '  R8VEC_INDEXED_HEAP_D_EXTRACT extracts the maximum value;'
      write ( *, '(a)' )
     &  '  R8VEC_INDEXED_HEAP_D_MAX reports the maximum value.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  These 3 operations are enough to model a priority queue.'
c
c  Set the data array.  To keep things easy, we will use the indicator vector.
c
      call r8vec_indicator ( m, a )
c
c  The index array will initially be a random subset of the numbers 1 to M,
c  in random order.
c
      n = 5
      indx(1)  =  9
      indx(2)  =  2
      indx(3)  =  8
      indx(4)  = 14
      indx(5)  =  5
      indx(6)  =  7
      indx(7)  = 15
      indx(8)  =  1
      indx(9)  = 19
      indx(10) = 20
      indx(11) =  3

      call r8vec_print ( m, a, '  The data vector:' )
      call i4vec_print ( n, indx, '  The index vector:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  A(INDX):'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i4,2x,g14.6)' ) i, a(indx(i))
      end do
c
c  Create the descending heap.
c
      call r8vec_indexed_heap_d ( n, a, indx )

      call i4vec_print ( n, indx, '  The index vector after heaping:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  A(INDX) after heaping:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i4,2x,g14.6)' ) i, a(indx(i))
      end do
c
c  Insert five entries, and monitor the maximum.
c
      do i = 1, 5

        indx_insert = indx(n+1)

        write ( *, '(a)' ) ' '
        write ( *, '(a,g14.6)' ) '  Inserting value ', a(indx_insert)

        call r8vec_indexed_heap_d_insert ( n, a, indx, indx_insert )

        call r8vec_indexed_heap_d_max ( n, a, indx, indx_max )

        write ( *, '(a,g14.6)' ) '  Current maximum is ', a(indx_max)

      end do
      call r8vec_print ( m, a, '  The data vector after insertions:' )
      call i4vec_print ( n, indx,
     &  '  The index vector after insertions:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  A(INDX) after insertions:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i4,2x,g14.6)' ) i, a(indx(i))
      end do
c
c  Extract the first 5 largest elements.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Now extract the maximum several times.'
      write ( *, '(a)' ) ' '

      do i = 1, 5
        call r8vec_indexed_heap_d_extract ( n, a, indx, indx_extract )
        write ( *, '(a,i8,a,g14.6)' ) '  Extracting maximum element A(',
     &    indx_extract,') = ', a(indx_extract)
      end do

      call r8vec_print ( m, a, '  The data vector after extractions:' )
      call i4vec_print ( n, indx,
     & '  The index vector after extractions:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  A(INDX) after extractions:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i4,2x,g14.6)' ) i, a(indx(i))
      end do

      return
      end
      subroutine test1256 ( )

c*********************************************************************72
c
cc TEST1256 tests R8VEC_LEGENDRE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 June 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 7 )

      integer ( kind = 4 ) n
      double precision r(n_max)
      double precision r1
      double precision r2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1256'
      write ( *, '(a)' ) 
     &  '  R8VEC_LEGENDRE computes N Legendre points in [R1,R2].'

      r1 = -1.0D+00
      r2 = +1.0D+00
      n = 5

      call r8vec_legendre ( n, r1, r2, r )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i4,a,g14.6,a,g14.6)' ) 
     &  '  N = ', n, '  R1 = ', r1, '  R2 = ', r2

      call r8vec_print ( n, r, '  Legendre points:' )

      r1 =   0.0D+00
      r2 = +10.0D+00
      n = 7

      call r8vec_legendre ( n, r1, r2, r )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i4,a,g14.6,a,g14.6)' ) 
     &  '  N = ', n, '  R1 = ', r1, '  R2 = ', r2

      call r8vec_print ( n, r, '  Legendre points:' )

      return
      end
      subroutine test1258 ( )

c*********************************************************************72
c
cc TEST1258 tests R8VEC_LINSPACE and R8VEC_MIDSPACE
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 June 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      double precision a
      double precision b
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1258'
      write ( *, '(a)' ) '  For a R8VEC:'
      write ( *, '(a)' ) 
     &  '  R8VEC_LINSPACE: evenly spaced points between A and B;'
      write ( *, '(a)' ) 
     &  '  R8VEC_MIDSPACE: evenly spaced midpoints between A and B'

      a = 10.0D+00
      b = 20.0D+00

      call r8vec_linspace ( n, a, b, x )
      call r8vec_print ( n, x, '  r8vec_linspace ( 5, 10, 20 )' )

      call r8vec_midspace ( n, a, b, x )
      call r8vec_print ( n, x, '  r8vec_midspace ( 5, 10, 20 )' )

      return
      end
      subroutine test1465 ( )

c*********************************************************************72
c
cc TEST1465 tests R8VEC_SORTED_RANGE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 September 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer i
      integer i_hi
      integer i_lo
      double precision r(n)
      double precision r_lo
      double precision r_hi
      double precision r8_uniform_01
      integer seed
      double precision t
      integer test

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1465'
      write ( *, '(a)' )
     &  '  R8VEC_SORTED_RANGE seeks the range of indices'
      write ( *, '(a)' ) '  in a sorted vector R so that'
      write ( *, '(a)' ) '  R_LO <= R(I_LO:I_HI) <= R_HI.'

      seed = 123456789

      do test = 1, 5

        call r8vec_uniform_01 ( n, seed, r )

        call r8vec_sort_heap_a ( n, r )

        call r8vec_print ( n, r, '  Sorted array R:' )

        r_lo = r8_uniform_01 ( seed )
        r_hi = r8_uniform_01 ( seed )

        if ( r_hi .lt. r_lo ) then
          t = r_lo
          r_lo = r_hi
          r_hi = t
        end if

        call r8vec_sorted_range ( n, r, r_lo, r_hi, i_lo, i_hi )

        write ( *, '(a)' ) ' '
        if ( i_hi .lt. i_lo ) then
          write ( *, '(2x,a4,2x,g14.6)' ) 'R_LO', r_lo
          write ( *, '(2x,a4,2x,g14.6)' ) 'R_HI', r_hi
          write ( *, '(2x,a)' ) '  Empty range in R.'
        else

          write ( *, '(2x,a4,2x,g14.6)' ) 'R_LO', r_lo
          do i = i_lo, i_hi
            write ( *, '(2x,i4,2x,g14.6)' ) i, r(i)
          end do
          write ( *, '(2x,a4,2x,g14.6)' ) 'R_HI', r_hi
        end if

      end do

      return
      end
      subroutine test1504 ( )

c*********************************************************************72
c
cc TEST1504 tests R8VEC_TRANSPOSE_PRINT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 November 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 12 )

      integer seed
      double precision x(n)

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1504'
      write ( *, '(a)' )
     &  '  R8VEC_TRANSPOSE_PRINT prints an R8VEC "tranposed",'
      write ( *, '(a)' )
     &  '  that is, placing multiple entries on a line.'

      call r8vec_uniform_01 ( n, seed, x )

      call r8vec_transpose_print ( n, x, '  The vector X:' )

      return
      end
      subroutine test1515 ( )

c*********************************************************************72
c
cc TEST1515 tests R8VEC_UNIFORM_01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 November 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 20 )

      double precision r(n)
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 3 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST1515'
      write ( *, '(a)' ) '  R8VEC_UNIFORM_01 returns a random R8VEC '
      write ( *, '(a)' ) '  with entries in [0,1].'

      seed = 123456789

      do test = 1, test_num

        write ( *, '(a)' ) ' '
        write ( *, '(a,i12)' ) '  Input SEED = ', seed

        call r8vec_uniform_01 ( n, seed, r )

        call r8vec_print_some ( n, r, 1, 10, '  Random vector:' )

      end do

      return
      end
      subroutine test157 ( )

c*********************************************************************72
c
cc TEST157 tests R8VEC2_SUM_MAX_INDEX.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 June 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      double precision a1(n)
      double precision a2(n)
      double precision b
      double precision c
      integer ival
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST157'
      write ( *, '(a)' ) '  For a pair of R8VEC''s:'
      write ( *, '(a)' ) '  R8VEC2_SUM_MAX_INDEX: index of the sum'
      write ( *, '(a)' ) '  vector with maximum value.'

      b = 0.0D+00
      c = 10.0D+00
      seed = 123456789

      call r8vec_uniform_ab ( n, b, c, seed, a1 )

      b = 0.0D+00
      c = 5.0D+00

      call r8vec_uniform_ab ( n, b, c, seed, a2 )

      call r8vec2_print ( n, a1, a2, '  The pair of vectors:' )

      call r8vec2_sum_max_index ( n, a1, a2, ival )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Index of maximum in A+B: ', ival

      return
      end
      subroutine test158 ( )

c*********************************************************************72
c
cc TEST158 tests R8VECS_PRINT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 June 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 5 )
      integer na
      parameter ( na = 15 )

      double precision a(na)
      integer nvec(m+1)

      save a
      save nvec

      data a /
     &  11.0D+00, 12.0D+00, 13.0D+00, 
     &  21.0D+00, 22.0D+00, 
     &  31.0D+00, 32.0D+00, 33.0D+00, 34.0D+00, 35.0D+00, 36.0D+00, 
     &  37.0D+00, 
     &  41.0D+00, 42.0D+00, 
     &  51.0D+00 /

      data nvec / 1, 4, 6, 13, 15, 16 /
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST158'
      write ( *, '(a)' ) '  R8VECS_PRINT prints a packed R8VEC.'

      call r8vecs_print ( m, nvec, na, a, '  Packed R8VEC:' )

      return
      end

