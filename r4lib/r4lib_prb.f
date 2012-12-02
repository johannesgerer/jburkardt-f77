      program main

c*********************************************************************72
c
cc MAIN is the main program for R4LIB_PRB.
c
c  Discussion:
c
c    R4LIB_PRB tests routines from the R4LIB library.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 June 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R4LIB_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the R4LIB library.'

      call test001 ( )
      call test002 ( )
      call test003 ( )
      call test004 ( )
      call test005 ( )
      call test006 ( )
      call test007 ( )
      call test009 ( )

      call test0235 ( )
      call test027 ( )
      call test028 ( )
      call test137 ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R4LIB_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test001 ( )

c*********************************************************************72
c
cc TEST001 tests R4_ABS.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 June 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      real r4
      real r4_abs
      real r4_absolute
      real r4_uniform
      real r4_hi
      parameter ( r4_hi = 5.0E+00 )
      real r4_lo
      parameter ( r4_lo = -3.0E+00 )
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST001'
      write ( *, '(a)' ) '  R4_ABS returns the absolute value of an R4.'
      write ( *, '(a)' ) ' '

      do test = 1, test_num
        r4 = r4_uniform ( r4_lo, r4_hi, seed )
        r4_absolute = r4_abs ( r4 )
        write ( *, '(2x,f10.6,2x,f10.6)' ) r4, r4_absolute
      end do

      return
      end
      subroutine test002 ( )

c*********************************************************************72
c
cc TEST002 tests R4_ATAN.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 September 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 8 )

      real r4_atan
      integer test
      real x
      real xtest(test_num)
      real y
      real ytest(test_num)

      save xtest
      save ytest

      data xtest /
     &   1.0E+00,  1.0E+00,  0.0E+00, -1.0E+00, 
     &  -1.0E+00, -1.0E+00,  0.0E+00,  1.0E+00 /
      data ytest /
     &   0.0E+00,  1.0E+00,  1.0E+00,  1.0E+00, 
     &   0.0E+00, -1.0E+00, -1.0E+00, -1.0E+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST002'
      write ( *, '(a)' ) 
     &  '  R4_ATAN computes the arc-tangent given Y and X;'
      write ( *, '(a)' ) 
     &  '  ATAN2 is the system version of this routine.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '       X             Y          ATAN2(Y,X)    R4_ATAN(Y,X)'
      write ( *, '(a)' ) ' '

      do test = 1, test_num
        x = xtest(test)
        y = ytest(test)
        write ( *, '(2x,4g14.6)' ) 
     &    x, y, atan2 ( y, x ), r4_atan ( y, x )
      end do

      return
      end
      subroutine test003 ( )

c*********************************************************************72
c
cc TEST003 tests R4_CAS.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 August 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      real r4_cas
      real r4_pi
      integer test_num
      parameter ( test_num = 12 )
      integer test
      real x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST003'
      write ( *, '(a)' ) '  R4_CAS evaluates the casine of a number.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        X           R4_CAS ( X )'
      write ( *, '(a)' ) ' '
      do test = 0, test_num
        x = r4_pi ( ) * real ( test ) / real ( test_num )
        write ( *, '(2x,g14.6,2x,g14.6)' ) x, r4_cas ( x )
      end do

      return
      end
      subroutine test004 ( )

c*********************************************************************72
c
cc TEST004 tests R4_CEILING.
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

      integer i
      integer r4_ceiling
      integer ival
      real rval

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST004'
      write ( *, '(a)' ) '  R4_CEILING rounds a value up.'
      write ( *, '(a)' ) ' '

      do i = -6, 6
        rval = real ( i ) / 5.0E+00
        ival = r4_ceiling ( rval )
        write ( *, '(2x,g14.6,2x,i8)' ) rval, ival
      end do

      return
      end
      subroutine test005 ( )

c*********************************************************************72
c
cc TEST005 tests R4_DIFF.
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

      integer test_num
      parameter ( test_num = 15 )

      integer ndig
      real r4_diff
      integer test
      real x
      real y_test(test_num)
      real y

      save y_test

      data y_test /
     &  0.0625E+00, 0.125E+00, 0.25E+00, 0.50E+00,  0.874E+00, 
     &  0.876E+00,  0.90E+00,  0.95E+00, 0.99E+00,  1.0E+00, 
     &  1.01E+00,   1.05E+00,  1.10E+00, 3.0E+00,  10.0E+00 /

      ndig = 3
      x = 1.0E+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST005'
      write ( *, '(a)' ) '  R4_DIFF computes a difference X-Y to a'
      write ( *, '(a)' ) '  given number of binary places.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a)' ) 
     &  '  For this test, we use ', ndig, ' binary places.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       X       Y       X-Y     R4_DIFF(X,Y)'
      write ( *, '(a)' ) ' '
      do test = 1, test_num
        y = y_test(test)
        write ( *, '(4f10.4)' ) x, y, x - y, r4_diff ( x, y, ndig )
      end do

      return
      end
      subroutine test006 ( )

c*********************************************************************72
c
cc TEST006 tests R4_DIGIT.
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

      integer maxdig
      parameter ( maxdig = 20 )

      integer i
      integer digit(-2:maxdig)
      integer idigit
      real r4_pi
      real x

      x = r4_pi ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST006'
      write ( *, '(a)' ) '  R4_DIGIT extracts decimal digits.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g24.16)' ) '  Here, we get digits of ', x
      write ( *, '(a)' ) ' '

      do idigit = -2, maxdig
        call r4_digit ( x, idigit, digit(idigit) )
      end do

      write ( *, '(2x,25i3)' ) ( i, i = -2, maxdig )
      write ( *, '(2x,25i3)' ) digit(-2:maxdig)

      return
      end
      subroutine test007 ( )

c*********************************************************************72
c
cc TEST007 tests R4_EPSILON
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 June 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      real r4_epsilon
      real r
      real s

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST007'
      write ( *, '(a)' ) 
     &  '  R4_EPSILON produces the R4 machine precision.'
      write ( *, '(a)' ) ' '

      r = r4_epsilon ( )
      write ( *, '(a,g14.6)' ) '  R = R4_EPSILON()   = ', r

      s = ( 1.0E+00 + r ) - 1.0E+00
      write ( *, '(a,g14.6)' ) '  ( 1 + R ) - 1      = ', s

      s = ( 1.0E+00 + ( r / 2.0E+00 ) ) - 1.0E+00
      write ( *, '(a,g14.6)' ) '  ( 1 + (R/2) ) - 1  = ', s

      return
      end
      subroutine test009

c*********************************************************************72
c
cc TEST009 tests R4_HUGE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 June 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      real r4_huge

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST009'
      write ( *, '(a)' ) '  R4_HUGE returns a "huge" R4;'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g24.16)' ) '    R4_HUGE ( ) =      ', r4_huge ( )

      return
      end
      subroutine test0235

c*********************************************************************72
c
cc TEST0235 tests R4_SWAP.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 June 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      real x
      real y

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0235'
      write ( *, '(a)' ) '  R4_SWAP swaps two reals.'

      x = 1.0E+00
      y = 3.141592653589793E+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Before swapping:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '    X = ', x
      write ( *, '(a,g14.6)' ) '    Y = ', y

      call r4_swap ( x, y )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  After swapping:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '    X = ', x
      write ( *, '(a,g14.6)' ) '    Y = ', y

      return
      end
      subroutine test027

c*********************************************************************72
c
cc TEST027 tests R4_UNIFORM_01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 June 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      real r4_uniform_01
      integer seed
      integer seed_old
      real x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST027'
      write ( *, '(a)' ) 
     &  '  R4_UNIFORM_01 produces a sequence of random values.'

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  Using random seed ', seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  SEED   R4_UNIFORM_01(SEED)'
      write ( *, '(a)' ) ' '
      do i = 1, 10
        seed_old = seed
        x = r4_uniform_01 ( seed )
        write ( *, '(2x,i12,2x,g14.6)' ) seed, x
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Verify that the sequence can be restarted.'
      write ( *, '(a)' ) 
     &  '  Set the seed back to its original value, and see that'
      write ( *, '(a)' ) '  we generate the same sequence.'

      seed = 123456789
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  SEED   R4_UNIFORM_01(SEED)'
      write ( *, '(a)' ) ' '

      do i = 1, 10
        seed_old = seed
        x = r4_uniform_01 ( seed )
        write ( *, '(2x,i12,2x,g14.6)' ) seed, x
      end do

      return
      end
      subroutine test028

c*********************************************************************72
c
cc TEST028 tests R4_UNIFORM_01
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 June 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      real mean
      integer n
      parameter ( n = 1000 )
      integer seed
      real r4_uniform_01
      real variance
      real x(n)
      real xmax
      real xmin

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST028'
      write ( *, '(a)' ) '  R4_UNIFORM_01 samples a uniform random'
      write ( *, '(a)' ) '  distribution in [0,1].'

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  Starting with seed = ', seed

      do i = 1, n
        x(i) = r4_uniform_01 ( seed )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  First few values:'
      write ( *, '(a)' ) ' '
      do i = 1, 5
        write ( *, '(2x,i8,2x,g14.6)' ) i, x(i)
      end do

      call r4vec_mean ( n, x, mean )

      call r4vec_variance ( n, x, variance )

      call r4vec_max ( n, x, xmax )
      call r4vec_min ( n, x, xmin )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of values computed was N = ', n
      write ( *, '(a,g14.6)' ) '  Average value was ', mean
      write ( *, '(a,g14.6)' ) 
     &  '  Minimum value was ', xmin
      write ( *, '(a,g14.6)' ) 
     &  '  Maximum value was ', xmax
      write ( *, '(a,g14.6)' ) '  Variance was ', variance

      return
      end
      subroutine test137

c*********************************************************************72
c
cc TEST137 tests R4VEC_SORT_BUBBLE_A.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 June 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 20 )

      real a(n)
      real b
      real c
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST137'
      write ( *, '(a)' ) 
     &  '  R4VEC_SORT_BUBBLE_A ascending sorts a R4VEC.'

      b = 0.0E+00
      c = 3.0E+00 * real ( n )
      seed = 123456789

      call r4vec_uniform ( n, b, c, seed, a )

      call r4vec_print_some ( n, a, 10, '  Original array:' )

      call r4vec_sort_bubble_a ( n, a )

      call r4vec_print_some ( n, a, 10, '  Ascending sorted array:' )

      return
      end
