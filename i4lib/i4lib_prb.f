      program main

c*********************************************************************72
c
cc MAIN is the main program for I4LIB_PRB.
c
c  Discussion:
c
c    I4LIB_PRB calls sample problems for the I4LIB library.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 December 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'I4LIB_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the I4LIB library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
      call test06 ( )
      call test07 ( )
      call test08 ( )
      call test09 ( )

      call test10 ( )
      call test11 ( )
      call test12 ( )
      call test13 ( )
      call test14 ( )
      call test15 ( )
      call test16 ( )
      call test17 ( )
      call test18 ( )
      call test19 ( )

      call test20 ( )
      call test21 ( )
      call test22 ( )
      call test23 ( )
      call test24 ( )
      call test245 ( )
      call test25 ( )
      call test26 ( )
      call test27 ( )
      call test28 ( )

      call test40 ( )

      call test50 ( )

      call test602 ( )
      call test605 ( )

      call test73 ( )

      call test80 ( )
      call test81 ( )
      call test82 ( )
      call test84 ( )
      call test85 ( )
      call test87 ( )
      call test88 ( )
      call test89 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'I4LIB_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests I4_BIT_HI1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      integer i4_uniform_ab
      integer i4_bit_hi1
      integer j
      integer seed
      integer test
      integer test_num

      seed = 123456789
      test_num = 10

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' )
     &  '  I4_BIT_HI1 returns the location of the high 1 bit.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       I  I4_BIT_HI1(I)'
      write ( *, '(a)' ) ' '

      do test = 1, test_num
        i = i4_uniform_ab ( 0, 100, seed )
        j = i4_bit_hi1 ( i )
        write ( *, '(2x,i8,2x,i8)' ) i, j
      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests I4_BIT_LO0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      integer i4_uniform_ab
      integer i4_bit_lo0
      integer j
      integer seed
      integer test
      integer test_num

      seed = 123456789
      test_num = 10

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' )
     &  '  I4_BIT_LO0 returns the location of the low 0 bit.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       I  I4_BIT_LO0(I)'
      write ( *, '(a)' ) ' '

      do test = 1, test_num
        i = i4_uniform_ab ( 0, 100, seed )
        j = i4_bit_lo0 ( i )
        write ( *, '(2x,i8,2x,i8)' ) i, j
      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests I4_BIT_LO1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      integer i4_uniform_ab
      integer i4_bit_lo1
      integer j
      integer seed
      integer test
      integer test_num

      test_num = 10
      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' )
     &  '  I4_BIT_LO1 returns the location of the low 1 bit.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       I  I4_BIT_LO1(I)'
      write ( *, '(a)' ) ' '

      do test = 1, test_num
        i = i4_uniform_ab ( 0, 100, seed )
        j = i4_bit_lo1 ( i )
        write ( *, '(2x,i8,2x,i8)' ) i, j
      end do

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests I4_BIT_REVERSE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      integer i_hi
      integer i4_bit_reverse
      integer j
      integer k

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' )
     &  '  I4_BIT_REVERSE bit reverses I with respect to 2^J'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         I         J  I4_BIT_REVERSE(I,J)'
      write ( *, '(a)' ) ' '

      do j = 0, 4
        i_hi = 2**j - 1
        do i = 0, i_hi
          k = i4_bit_reverse ( i, j )
          write ( *, '(2x,i8,2x,i8,2x,i8)' ) i, j, k
        end do
      end do

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests I4_CHARACTERISTIC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      integer i4_characteristic

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' )
     &  '  I4_CHARACTERISTIC computes the characteristic'
      write ( *, '(a)' ) '  of an integer Q, which is  '
      write ( *, '(a)' ) '    Q if Q is prime;'
      write ( *, '(a)' ) '    P, if Q = P**N for some prime P;'
      write ( *, '(a)' )
     &  '    0, if Q is negative, 0, 1, or the product of '
      write ( *, '(a)' ) '      more than 1 distinct prime.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  I  I4_CHARACTERISTIC'
      write ( *, '(a)' ) ' '

      do i = 1, 50
        write ( *, '(2x,i2,13x,i4)' ) i, i4_characteristic ( i )
      end do

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 tests I4_DIV_ROUNDED.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer a
      integer a_hi
      integer a_lo
      integer b
      integer b_hi
      integer b_lo
      double precision c0
      integer c1
      integer c2
      integer c3
      integer c4
      integer i4_div_rounded
      integer i4_uniform_ab
      integer seed
      integer test
      integer test_num

      a_hi = 100
      a_lo = - 100
      b_hi = 10
      b_lo = - 10
      test_num = 20

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' )
     &  '  I4_DIV_ROUNDED performs rounded integer division.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  C0 = real ( a ) / real ( b )'
      write ( *, '(a)' ) '  C1 = I4_DIV_ROUNDED ( A, B )'
      write ( *, '(a)' ) '  C2 = nint ( real ( a ) / real ( b ) )'
      write ( *, '(a)' ) '  C3 = A / B'
      write ( *, '(a)' ) '  C4 = int ( real ( a ) / real ( b ) )'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  C1 and C2 should be equal;'
      write ( *, '(a)' ) '  C3 and C4 should be equal.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '     A     B           C0         C1    C2      C3    C4'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num
        a = i4_uniform_ab ( a_lo, a_hi, seed )
        b = i4_uniform_ab ( b_lo, b_hi, seed )
        if ( b .eq. 0 ) then
          b = 7
        end if
        c0 = dble ( a ) / dble ( b )
        c1 = i4_div_rounded ( a, b )
        c2 = nint ( dble ( a ) / dble ( b ) )
        c3 = a / b
        c4 = int ( dble ( a ) / dble ( b ) )
        write ( *, '(2x,i4,2x,i4,4x,f14.6,2x,i4,2x,i4,4x,i4,2x,i4)' )
     &    a, b, c0, c1, c2, c3, c4
      end do

      return
      end
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07 tests I4_DIVP.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer a
      integer a_hi
      integer a_lo
      integer b
      integer b_hi
      integer b_lo
      integer c
      integer d
      integer i4_divp
      integer i4_uniform_ab
      integer seed
      integer test
      integer test_num

      a_hi = 100
      a_lo = - 100
      b_hi = 10
      b_lo = - 10
      test_num = 20

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07'
      write ( *, '(a)' )
     &  '  I4_DIVP(A,B) returns the smallest multiplier of J'
      write ( *, '(a)' ) '  that is less than I'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     A     B     C     D'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num
        a = i4_uniform_ab ( a_lo, a_hi, seed )
        b = i4_uniform_ab ( b_lo, b_hi, seed )
        if ( b .eq. 0 ) then
          b = 7
        end if
        c = i4_divp ( a, b )
        d = c * b
        write ( *, '(2x,i4,2x,i4,2x,i4,2x,i4)' ) a, b, c, d
      end do

      return
      end
      subroutine test08 ( )

c*********************************************************************72
c
cc TEST08 tests I4_GCD.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 November 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 7 )

      integer i
      integer i4_gcd
      integer i_test(test_num)
      integer j
      integer j_test(test_num)
      integer test

      data i_test / 36, 49, 0, 12, 36, 1, 91 /
      data j_test / 30, -7, 71, 12, 49, 42, 28/

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08'
      write ( *, '(a)' ) '  I4_GCD computes the greatest common factor,'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I     J   I4_GCD'
      write ( *, '(a)' ) ' '

      do test = 1, test_num
        i = i_test(test)
        j = j_test(test)
        write ( *, '(2x,3i8)') i, j, i4_gcd ( i, j )
      end do

      return
      end
      subroutine test09 ( )

c*********************************************************************72
c
cc TEST09 tests I4_HUGE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i4_huge

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST09'
      write ( *, '(a)' ) '  I4_HUGE returns a huge integer.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  I4_HUGE() = ', i4_huge ( )

      return
      end
      subroutine test10 ( )

c***********************************************************************72
c
cc TEST10 tests I4_HUGE_NORMALIZER.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i4
      integer i4_huge
      double precision i4_huge_normalizer
      double precision r8
      double precision value

      i4 = i4_huge ( )
      r8 = i4_huge_normalizer ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST10'
      write ( *, '(a)' ) '  I4_HUGE_NORMALIZER returns 1/(I4_HUGE+1).'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  I4_HUGE() = ', i4
      write ( *, '(a,g14.6)' ) '  I4_HUGE_NORMALIZER() = ', r8

      value = dble ( i4 ) * r8

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' )
     &  '  I4_HUGE * I4_HUGE_NORMALIZER = ', value

      return
      end
      subroutine test11 ( )

c*********************************************************************72
c
cc TEST11 tests I4_IS_PRIME.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      logical i4_is_prime

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST11'
      write ( *, '(a)' )
     &  '  I4_IS_PRIME reports whether an integer is prime.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  I     I4_IS_PRIME(I)'
      write ( *, '(a)' ) ' '

      do i = -2, 25
        write ( *, '(2x,i8,2x,l1)' ) i, i4_is_prime ( i )
      end do

      return
      end
      subroutine test12 ( )

c*********************************************************************72
c
cc TEST12 tests I4_LCM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 7 )

      integer i
      integer i4_lcm
      integer i_test(test_num)
      integer j
      integer j_test(test_num)
      integer test

      save i_test
      save j_test

      data i_test / 36, 49, 0, 12, 36, 1, 91 /
      data j_test / 30, -7, 71, 12, 49, 42, 28 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST12'
      write ( *, '(a)' ) '  I4_LCM computes the least common multiple.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I     J   I4_LCM'
      write ( *, '(a)' ) ' '

      do test = 1, test_num
        i = i_test(test)
        j = j_test(test)
        write ( *, '(2x,4i8)') i, j, i4_lcm ( i, j )
      end do

      return
      end
      subroutine test13 ( )

c*********************************************************************72
c
cc TEST13 tests I4_LOG_10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 13 )

      integer i4_log_10
      integer test
      integer x
      integer x_test(test_num)

      save x_test

      data x_test /
     &  0, 1, 2, 3, 9, 10, 11, 99, 101, -1, -2, -3, -9 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST13'
      write ( *, '(a)' ) '  I4_LOG_10: whole part of log base 10,'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X, I4_LOG_10'
      write ( *, '(a)' ) ' '

      do test = 1, test_num
        x = x_test(test)
        write ( *, '( 2x, i8, i12 )' ) x, i4_log_10 ( x )
      end do

      return
      end
      subroutine test14 ( )

c*********************************************************************72
c
cc TEST14 tests I4_LOG_2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 17 )

      integer i4_log_2
      integer test
      integer x
      integer x_test(test_num)

      save x_test

      data x_test /
     &    0,    1,    2,    3,    9,
     &   10,   11,   99,  101,   -1,
     &   -2,   -3,   -9, 1000, 1023,
     & 1024, 1025 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST14'
      write ( *, '(a)' ) '  I4_LOG_2: whole part of log base 2.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       X     I4_LOG_2'
      write ( *, '(a)' ) ' '

      do test = 1, test_num
        x = x_test(test)
        write ( *, '( 2x, i8, i12 )' ) x, i4_log_2 ( x )
      end do

      return
      end
      subroutine test15 ( )

c*********************************************************************72
c
cc TEST15 tests I4_LOG_I4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i4
      integer i4_log_i4
      integer j4

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST15'
      write ( *, '(a)' ) '  I4_LOG_I4: logarith of I4 base J4,'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        I4        J4 I4_LOG_I4'
      write ( *, '(a)' ) ' '

      do j4 = 2, 5
        do i4 = 0, 10
          write ( *, '(2x, i8, 2x, i8, 2x, i8 )' )
     &      i4, j4, i4_log_i4 ( i4, j4 )
        end do
        write ( *, '(a)' ) ' '
      end do

      return
      end
      subroutine test16 ( )

c*********************************************************************72
c
cc TEST16 tests I4_LOG_R8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 July 2010
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
      integer i4_log_r8
      integer test
      integer x

      save b_test

      data b_test /
     &  2.0D+00, 3.0D+00,  4.0D+00,  5.0D+00,   6.0D+00,
     &  7.0D+00, 8.0D+00, 16.0D+00, 32.0D+00, 256.0D+00 /

      x = 16

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST16'
      write ( *, '(a)' ) '  I4_LOG_R8: whole part of log base B,'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X, B, I4_LOG_R8'
      write ( *, '(a)' ) ' '

      do test = 1, test_num

        b = b_test(test)

        write ( *, '(2x, i8, g14.6, i12 )' ) x, b, i4_log_r8 ( x, b )

      end do

      return
      end
      subroutine test17 ( )

c*********************************************************************72
c
cc TEST17 tests I4_MANT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer is
      integer j
      integer k
      integer l
      double precision x

      x = -314.159D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST17'
      write ( *, '(a)' ) '  I4_MANT decomposes an integer,'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Number to be decomposed is X = ', x

      call i4_mant ( x, is, j, k, l )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a,i24,a,i24,a,i8)' )
     &  '  I4_MANT: X = ', is, ' * (', j, '/', k, ') * 2**', l

      return
      end
      subroutine test18 ( )

c*********************************************************************72
c
cc TEST18 tests I4_MODDIV;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 4 )

      integer ndivid(test_num)
      integer nmult
      integer nrem
      integer number(test_num)
      integer test

      save ndivid
      save number

      data ndivid / 50, -50, 50, -50 /
      data number / 107, 107, -107, -107 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST18'
      write ( *, '(a)' ) '  I4_MODDIV factors a number'
      write ( *, '(a)' ) '    into a multiple and a remainder.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    Number   Divisor  Multiple Remainder'
      write ( *, '(a)' ) ' '

      do test = 1, test_num
        call i4_moddiv ( number(test), ndivid(test), nmult, nrem )
        write ( *, '(2x,4i10)' ) number(test), ndivid(test), nmult, nrem
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Repeat using FORTRAN MOD:'
      write ( *, '(a)' ) ' '

      do test = 1, test_num
        nrem = mod ( number(test), ndivid(test) )
        nmult = number(test) / ndivid(test)
        write ( *, '(2x,4i10)' ) number(test), ndivid(test), nmult, nrem
      end do

      return
      end
      subroutine test19 ( )

c*********************************************************************72
c
cc TEST19 tests I4_MODP.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 4 )

      integer i4_modp
      integer ndivid(test_num)
      integer nmult
      integer nrem
      integer number(test_num)
      integer test

      save ndivid
      save number

      data ndivid / 50, -50, 50, -50 /
      data number / 107, 107, -107, -107 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST19'
      write ( *, '(a)' ) '  I4_MODP factors a number'
      write ( *, '(a)' ) '    into a multiple and a remainder.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    Number   Divisor  Multiple Remainder'
      write ( *, '(a)' ) ' '

      do test = 1, test_num
        nrem = i4_modp ( number(test), ndivid(test) )
        nmult = ( number(test) - nrem ) / ndivid(test)
        write ( *, '(2x,4i10)' ) number(test), ndivid(test), nmult, nrem
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Repeat using FORTRAN MOD:'
      write ( *, '(a)' ) ' '

      do test = 1, test_num
        nrem = mod ( number(test), ndivid(test) )
        nmult = number(test) / ndivid(test)
        write ( *, '(2x,4i10)' ) number(test), ndivid(test), nmult, nrem
      end do

      return
      end
      subroutine test20 ( )

c*********************************************************************72
c
cc TEST20 tests I4_SIGN.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 5 )

      integer i4_sign
      integer test
      integer x
      integer x_test(test_num)

      save x_test

      data x_test / -10, -7, 0, 5, 9 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST20'
      write ( *, '(a)' ) '  I4_SIGN returns the sign of a number.'
      write ( *, '(a)' ) ' '

      do test = 1, test_num
        x = x_test(test)
        write ( *, '(2i8)' ) x, i4_sign ( x )
      end do

      return
      end
      subroutine test21 ( )

c*********************************************************************72
c
cc TEST21 tests I4_SWAP.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 November 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      integer j

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST21'
      write ( *, '(a)' ) '  I4_SWAP swaps two integers.'

      i = 1
      j = 202

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Before swapping: '
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  I = ', i
      write ( *, '(a,i8)' ) '  J = ', j

      call i4_swap ( i, j )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  After swapping: '
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  I = ', i
      write ( *, '(a,i8)' ) '  J = ', j

      return
      end
      subroutine test22 ( )

c*********************************************************************72
c
cc TEST22 tests I4_WALSH_1D;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      integer i4_walsh_1d
      integer w0
      integer wm1
      integer wm2
      integer wm3
      integer wp1
      integer wp2
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST22'
      write ( *, '(a)' ) '  I4_WALSH_1D evaluates 1D Walsh functions:'
      write ( *, '(a)' ) ' '
      write ( *, * ) 'X  W(+2) W(+1) W(0) W(-1) W(-2) W(-3)'
      write ( *, '(a)' ) ' '

      do i = 0, 32

        x = dble ( i ) / 4.0D+00

        wp2 = i4_walsh_1d ( x,  2 )
        wp1 = i4_walsh_1d ( x,  1 )
        w0  = i4_walsh_1d ( x,  0 )
        wm1 = i4_walsh_1d ( x, -1 )
        wm2 = i4_walsh_1d ( x, -2 )
        wm3 = i4_walsh_1d ( x, -3 )

        write ( *, '(2x,f10.6,6i2)' ) x, wp2, wp1, w0, wm1, wm2, wm3

      end do

      return
      end
      subroutine test23 ( )

c*********************************************************************72
c
cc TEST23 tests I4_WRAP.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 July 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      integer i4_wrap
      integer ihi
      integer ilo

      ilo = 4
      ihi = 8

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST23'
      write ( *, '(a)' )
     &  '  I4_WRAP forces an integer to lie within given limits.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  ILO = ', ilo
      write ( *, '(a,i8)' ) '  IHI = ', ihi
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I  I4_WRAP(I)'
      write ( *, '(a)' ) ' '

      do i = -10, 20
        write ( *, '(2x,2i8)' ) i, i4_wrap ( i, ilo, ihi )
      end do

      return
      end
      subroutine test24 ( )

c*********************************************************************72
c
cc TEST24 tests I4_XOR.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 January 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      integer i4_hi
      integer i4_lo
      integer i4_uniform_ab
      integer i4_xor
      integer j
      integer k
      integer l
      integer seed
      integer test
      integer test_num

      i4_hi = 100
      i4_lo = 0
      seed = 123456789
      test_num = 10

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST24'
      write ( *, '(a)' ) '  I4_XOR returns the bitwise exclusive OR of'
      write ( *, '(a)' ) '  two integers.'
      write ( *, '(a)' ) '  Compare with FORTRAN90 intrinsic IEOR.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         I         J    I4_XOR      IEOR'
      write ( *, '(a)' ) ' '

      do test = 1, test_num
        i = i4_uniform_ab ( i4_lo, i4_hi, seed )
        j = i4_uniform_ab ( i4_lo, i4_hi, seed )
        k = i4_xor ( i, j )
        l = ieor ( i, j )
        write ( *, '(2x,i8,2x,i8,2x,i8,2x,i8)' ) i, j, k,l
      end do

      return
      end
      subroutine test245 ( )

c*********************************************************************72
c
cc TEST245 tests I4BLOCK_PRINT.
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
      integer x(l,m,n) 

      save x

      data x /
     &  1,  2,  3,   4,  1, 
     &  4,  9, 16,   1,  8, 
     & 27, 64,  2,   4,  6, 
     &  8,  2,  8,  18, 32, 
     &  2, 16, 54, 128 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST245'
      write ( *, '(a)' ) '  I4BLOCK_PRINT prints an I4BLOCK.'

      call i4block_print ( l, m, n, x, '  The 3D array:' )

      return
      end
      subroutine test25 ( )

c*********************************************************************72
c
cc TEST25 tests I4COL_FIND_ITEM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 December 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      integer n
      integer test_num

      parameter ( m = 5 )
      parameter ( n = 4 )
      parameter ( test_num = 3 )

      integer a(m,n)
      integer col
      integer i
      integer item
      integer item_test(test_num)
      integer j
      integer row
      integer test

      save item_test

      data item_test / 34, 12, 90 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST25'
      write ( *, '(a)' ) '  I4COL_FIND_ITEM finds the first occurrence'
      write ( *, '(a)' ) '  of an item in an integer array of columns.'

      do i = 1, m
        do j = 1, n
          a(i,j) = 10 * i + j
        end do
      end do

      call i4mat_print ( m, n, a, '  The matrix of columns:' )

      do test = 1, test_num

        item = item_test(test)

        call i4col_find_item ( m, n, a, item, row, col )

        write ( *, '(a,i8,a,i8,a,i8)' ) '  Item ', item,
     &    '  occurs in row ', row, ' and column ', col

      end do

      return
      end
      subroutine test26 ( )

c*********************************************************************72
c
cc TEST26 tests I4COL_FIND_PAIR_WRAP.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 December 2005
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 5 )
      integer n
      parameter ( n = 4 )
      integer test_num
      parameter ( test_num = 5 )

      integer a(m,n)
      integer col
      integer i
      integer item1
      integer item1_test(test_num)
      integer item2
      integer item2_test(test_num)
      integer j
      integer row
      integer test

      save item1_test
      save item2_test

      data item1_test / 22, 32, 22, 54, 54 /
      data item2_test / 32, 22, 23, 14, 11 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST26'
      write ( *, '(a)' )
     &  '  I4COL_FIND_PAIR_WRAP finds the first occurrence of'
      write ( *, '(a)' )
     &  '  a pair of item in an integer array of columns.'
      write ( *, '(a)' )
     &  '  Items in the array are ordered by column, and'
      write ( *, '(a)' ) '  wraparound is allowed.'

      do i = 1, m
        do j = 1, n
          a(i,j) = 10 * i + j
        end do
      end do

      call i4mat_print ( m, n, a, '  The matrix of columns:' )

      do test = 1, test_num

        item1 = item1_test(test)
        item2 = item2_test(test)

        call i4col_find_pair_wrap ( m, n, a, item1, item2, row, col )

        write ( *, '(a,i8,a,i8,a,i8,a,i8)' ) '  Item ', item1,
     &    ' followed by item ', item2, ' occurs in row ',
     &    row, ' and column ', col

      end do

      return
      end
      subroutine test27 ( )

c*********************************************************************72
c
cc TEST27 tests I4COL_SORT_A and I4COL_SORT_D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 December 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 5 )
      integer n
      parameter ( n = 4 )

      integer a(m,n)
      integer b
      integer c
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST27'
      write ( *, '(a)' ) '  I4COL_SORT_A ascending sorts an integer'
      write ( *, '(a)' ) '  array as a table of columns.'
      write ( *, '(a)' ) '  I4COL_SORT_D descending sorts an integer'
      write ( *, '(a)' ) '  array as a table of columns.'

      b = 1
      c = 10
      seed = 123456789

      call i4mat_uniform_ab ( m, n, b, c, seed, a )

      call i4mat_print ( m, n, a, '  The original matrix:' )

      call i4col_sort_a ( m, n, a )

      call i4mat_print ( m, n, a, '  Ascending sorted:' )

      call i4col_sort_d ( m, n, a )

      call i4mat_print ( m, n, a, '  Descending sorted:' )

      return
      end
      subroutine test28 ( )

c*********************************************************************72
c
cc TEST28 tests I4COL_SORT2_A;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 December 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 6 )
      integer n
      parameter ( n = 4 )

      integer a(m,n)
      integer b
      integer c
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST28'
      write ( *, '(a)' ) '  For a rectangular integer matrix:'
      write ( *, '(a)' )
     &  '  I4COL_SORT2_D sorts the elements of the columns.'

      b = 0
      c = 20
      seed = 123456789

      call i4mat_uniform_ab ( m, n, b, c, seed, a )

      call i4mat_print ( m, n, a, '  The matrix:' )

      call i4col_sort2_a ( m, n, a )

      call i4mat_print ( m, n, a,
     &  '  The element-sorted column matrix:' )

      return
      end
      subroutine test40 ( )

c*********************************************************************72
c
cc TEST40 tests I4ROW_SWAP;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 December 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 3 )
      integer n
      parameter ( n = 4 )

      integer a(m,n)
      integer i
      integer j
      integer k
      integer row1
      integer row2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST40'
      write ( *, '(a)' ) '  For an integer matrix of rows,'
      write ( *, '(a)' ) '  I4ROW_SWAP swaps two rows;'

      k = 0
      do i = 1, m
        do j = 1, n
          k = k + 1
          a(i,j) = k
        end do
      end do

      call i4mat_print ( m, n, a, '  The matrix:' )

      row1 = 1
      row2 = 3

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a,i8)' ) '  Swap rows ', row1, ' and ', row2
      write ( *, '(a)' ) ' '

      call i4row_swap ( m, n, a, row1, row2 )

      call i4mat_print ( m, n, a, '  The new matrix:' )

      return
      end
      subroutine test50 ( )

c*********************************************************************72
c
cc TEST50 tests I4VEC_CUM and I4VEC_CUM0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 December 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer a(n)
      integer a_cum(n)
      integer a_cum0(0:n)
      integer b
      integer c
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST50'
      write ( *, '(a)' ) '  For an integer vector:'
      write ( *, '(a)' ) '  I4VEC_CUM:   cumulative sum;'
      write ( *, '(a)' ) '  I4VEC_CUM0:  cumulative sum, zero based;'

      seed = 123456789

      b = -n
      c = n

      call i4vec_uniform_ab ( n, b, c, seed, a )

      call i4vec_print ( n, a, '  Input vector:' )

      call i4vec_cum ( n, a, a_cum )

      call i4vec_print ( n, a_cum, '  Cumulative sums:' )

      call i4vec_cum0 ( n, a, a_cum0 )

      call i4vec_print ( n + 1, a_cum0, '  0-based Cumulative sums:' )

      return
      end
      subroutine test602 ( )

c*********************************************************************72
c
cc TEST602 tests I4VEC_INDEXED_HEAP_D;
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

      integer a(m)
      integer i
      integer indx(n)

      save a
      save indx

      data a /
     &  101, 102, 103, 104, 105, 106, 107, 108, 109, 110,
     &  111, 112, 113, 114, 115, 116, 117, 118, 119, 120 /
      data indx /
     &  1, 11, 17, 5, 7, 13, 15, 3, 19, 9 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST602'
      write ( *, '(a)' ) '  I4VEC_INDEXED_HEAP_D creates a descending'
      write ( *, '(a)' ) '  heap from an indexed vector.'
c
c  Print before.
c
      call i4vec_print ( m, a, '  The data vector:' )
      call i4vec_print ( n, indx, '  The index vector:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  A(INDX):'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i4,2x,i4)' ) i, a(indx(i))
      end do
c
c  Heap the data.
c
      call i4vec_indexed_heap_d ( n, a, indx )
c
c  Print afterwards.  Only INDX should change.
c
      call i4vec_print ( m, a,
     &  '  The data vector (should NOT change):' )
      call i4vec_print ( n, indx, '  The index vector (may change):' )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  A(INDX) is now a descending heap:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i4,2x,i4)' ) i, a(indx(i))
      end do

      return
      end
      subroutine test605 ( )

c*********************************************************************72
c
cc TEST605 tests I4VEC_INDEXED_HEAP_D_EXTRACT and related routines.
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

      integer a(m)
      integer i
      integer indx(n_max)
      integer indx_extract
      integer indx_insert
      integer indx_max
      integer n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST605'
      write ( *, '(a)' ) '  For an indexed I4VEC,'
      write ( *, '(a)' )
     &  '  I4VEC_INDEXED_HEAP_D_INSERT inserts a value into the heap.'
      write ( *, '(a)' )
     &  '  I4VEC_INDEXED_HEAP_D_EXTRACT extracts the maximum value;'
      write ( *, '(a)' )
     &  '  I4VEC_INDEXED_HEAP_D_MAX reports the maximum value.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  These 3 operations are enough to model a priority queue.'
c
c  Set the data array.  To keep things easy, we will use the indicator vector.
c
      call i4vec_indicator ( m, a )
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

      call i4vec_print ( m, a, '  The data vector:' )
      call i4vec_print ( n, indx, '  The index vector:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  A(INDX):'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i4,2x,i4)' ) i, a(indx(i))
      end do
c
c  Create a descending heap from the indexed array.
c
      call i4vec_indexed_heap_d ( n, a, indx )

      call i4vec_print ( n, indx, '  The index vector after heaping:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  A(INDX) after heaping:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i4,2x,i4)' ) i, a(indx(i))
      end do
c
c  Insert five entries, and monitor the maximum.
c
      do i = 1, 5

        indx_insert = indx(n+1)

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' )
     &    '  Inserting value          ', a(indx_insert)

        call i4vec_indexed_heap_d_insert ( n, a, indx, indx_insert )

        call i4vec_indexed_heap_d_max ( n, a, indx, indx_max )

        write ( *, '(a,i8)' ) '  Current maximum is ', a(indx_max)

      end do
      call i4vec_print ( m, a, '  The data vector after insertions:' )
      call i4vec_print ( n, indx,
     &  '  The index vector after insertions:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  A(INDX) after insertions:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i4,2x,i4)' ) i, a(indx(i))
      end do
c
c  Extract the first 5 largest elements.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Now extract the maximum several times.'
      write ( *, '(a)' ) ' '

      do i = 1, 5
        call i4vec_indexed_heap_d_extract ( n, a, indx, indx_extract )
        write ( *, '(a,i8,a,i8)' ) '  Extracting maximum element A(',
     &    indx_extract,') = ', a(indx_extract)
      end do

      call i4vec_print ( m, a, '  The data vector after extractions:' )
      call i4vec_print ( n, indx,
     &  '  The index vector after extractions:' )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  A(INDX) after extractions:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i4,2x,i4)' ) i, a(indx(i))
      end do

      return
      end
      subroutine test73 ( )

c*********************************************************************72
c
cc TEST73 tests I4VEC_RUN_COUNT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 20 )

      integer a(n)
      integer i
      integer run_count
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST73'
      write ( *, '(a)' ) '  I4VEC_RUN_COUNT counts runs in an I4VEC'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) ' Run Count        Sequence'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        call i4vec_uniform_ab ( n, 0, 1, seed, a )

        call i4vec_run_count ( n, a, run_count )

        write ( *, '(2x,i8,8x,20i2)' ) run_count, ( a(i), i = 1, n )

      end do

      return
      end
      subroutine test80 ( )

c*********************************************************************72
c
cc TEST80 tests I4VEC_SORT_INSERT_A.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 January 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer a(n)
      integer b
      integer c
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST80'
      write ( *, '(a)' ) '  I4VEC_SORT_INSERT_A sorts an integer array.'

      seed = 123456789

      b = 0
      c = n

      call i4vec_uniform_ab ( n, b, c, seed, a )

      call i4vec_print ( n, a, '  Unsorted array:' )

      call i4vec_sort_insert_a ( n, a )

      call i4vec_print ( n, a, '  Sorted array:' )

      return
      end
      subroutine test81 ( )

c*****************************************************************************80
c
cc TEST81 tests I4VEC_SORT_QUICK_A.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 20 )

      integer a(n)
      integer b
      integer c
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST81'
      write ( *, '(a)' ) '  I4VEC_SORT_QUICK_A sorts an integer vector'
      write ( *, '(a)' ) '  using quick sort.'

      seed = 123456789

      b = 0
      c = 3 * n

      call i4vec_uniform_ab ( n, b, c, seed, a )

      call i4vec_print ( n, a, '  Unsorted array:' )

      call i4vec_sort_quick_a ( n, a )

      call i4vec_print ( n, a, '  Sorted array:' )

      return
      end
      subroutine test82 ( )

c*********************************************************************72
c
cc TEST82 tests I4VEC_SORT_SHELL_A.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 January 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 20 )

      integer a(n)
      integer b
      integer c
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST82'
      write ( *, '(a)' ) '  I4VEC_SORT_SHELL_A sorts an integer vector'
      write ( *, '(a)' ) '  using Shell''s sort.'

      seed = 123456789

      b = 0
      c = 3 * n

      call i4vec_uniform_ab ( n, b, c, seed, a )

      call i4vec_print ( n, a, '  Unsorted array:' )

      call i4vec_sort_shell_a ( n, a )

      call i4vec_print ( n, a, '  Sorted array:' )

      return
      end
      subroutine test84 ( )

c*********************************************************************72
c
cc TEST84 tests I4VEC_SORTED_UNIQUE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 December 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 20 )

      integer a(n)
      integer b
      integer c
      integer seed
      integer unique_num

      b = 0
      c = n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST84'
      write ( *, '(a)' ) '  I4VEC_SORTED_UNIQUE finds unique entries'
      write ( *, '(a)' ) '  in a sorted array.'

      seed = 123456789

      call i4vec_uniform_ab ( n, b, c, seed, a )

      call i4vec_sort_heap_a ( n, a )

      call i4vec_print ( n, a, '  Input vector:' )

      call i4vec_sorted_unique ( n, a, unique_num )

      call i4vec_print ( unique_num, a, '  Unique entries:' )

      return
      end
      subroutine test85 ( )

c*********************************************************************72
c
cc TEST85 tests I4VEC_TRANSPOSE_PRINT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 December 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 12 )

      integer a(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST85'
      write ( *, '(a)' )
     &  '  I4VEC_TRANSPOSE_PRINT prints an integer vector'
      write ( *, '(a)' )
     &  '  with 5 entries to a row, and an optional title.'

      call i4vec_indicator ( n, a )

      call i4vec_print ( n, a, '  Output from I4VEC_PRINT:' )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  Call I4VEC_TRANSPOSE_PRINT with a short title:'

      call i4vec_transpose_print ( n, a, '  My array:  ' )

      return
      end
      subroutine test87 ( )

c*********************************************************************72
c
cc TEST87 tests I4VEC_UNIQUE_INDEX.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c   30 December 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 20 )

      integer a(n)
      integer b
      integer c
      integer i
      integer seed
      integer unique_index(n)

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST87'
      write ( *, '(a)' )
     &  '  I4VEC_UNIQUE_INDEX, for each entry in an I4VEC'
      write ( *, '(a)' ) '  indexes the unique elements.'

      b = 1
      c = 5

      call i4vec_uniform_ab ( n, b, c, seed, a )

      call i4vec_unique_index ( n, a, unique_index )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         I      A(I)    UNIQUE'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,i8,2x,i8,2x,i8)' ) i, a(i), unique_index(i)
      end do

      return
      end
      subroutine test88 ( )

c*********************************************************************72
c
cc TEST88 tests I4VEC_VALUE_INDEX.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 December 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer max_index
      parameter ( max_index = 3 )
      integer n
      parameter ( n = 25 )

      integer a(n)
      integer b
      integer c
      integer n_index
      integer seed
      integer value
      integer value_index(max_index)

      seed = 123456789

      value = 3

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST88'
      write ( *, '(a)' ) '  I4VEC_VALUE_INDEX indexes entries equal to'
      write ( *, '(a)' ) '  a given value.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The desired value is ', value
      write ( *, '(a,i8)' )
     &  '  Maximum number of indices to find is ', max_index

      b = 1
      c = 5

      call i4vec_uniform_ab ( n, b, c, seed, a )

      call i4vec_print ( n, a, '  Input vector A:' )

      call i4vec_value_index ( n, a, value, max_index, n_index,
     &  value_index )

      call i4vec_print ( n_index, value_index,
     &  '  Indices of entries equal to given value: ' )

      return
      end
      subroutine test89 ( )

c*********************************************************************72
c
cc TEST89 tests I4VEC_VARIANCE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 December 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer a(n)
      integer b
      integer c
      integer seed
      double precision variance

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST89'
      write ( *, '(a)' ) '  For an integer vector:'
      write ( *, '(a)' ) '  I4VEC_VARIANCE:      variance.'

      seed = 123456789

      b = -n
      c = n

      call i4vec_uniform_ab ( n, b, c, seed, a )

      call i4vec_print ( n, a, '  Input vector:' )

      call i4vec_variance ( n, a, variance )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Variance:                 ', variance

      return
      end
