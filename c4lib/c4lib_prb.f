      program main

c*********************************************************************72
c
cc MAIN is the main program for C4LIB_PRB.
c
c  Discussion:
c
c    C4LIB_PRB tests routines from the C4LIB library.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'C4LIB_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the C4LIB library.'

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
      call test25 ( )
      call test26 ( )
      call test27 ( )
      call test28 ( )
      call test29 ( )

      call test30 ( )
      call test31 ( )
      call test32 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'C4LIB_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests C4_ABS.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex c1
      real c4_abs
      complex c4_uniform_01
      integer i
      real r2
      real r3
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  C4_ABS computes the absolute value of a C4.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C4_UNIFORM_01        ' // 
     &  '  R2=C4_ABS(C1)             R3=CABS(C1)' 
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do i = 1, 10
        c1 = c4_uniform_01 ( seed )
        r2 = c4_abs ( c1 )
        r3 = cabs ( c1 )
        write ( *, '(2x,2f12.6,2x,f12.6,12x,2x,f12.6)' ) c1, r2, r3
      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests C4_ACOS and C4_COS.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex c1
      complex c2
      complex c3
      complex c4_acos
      complex c4_cos
      complex c4_uniform_01
      integer i
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  C4_ACOS computes the inverse cosine.'
      write ( *, '(a)' ) '  C4_COS computes the cosine.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C4_UNIFORM_01        ' // 
     &  '  C2 = C4_ACOS(C1)           C3 = C4_COS(C2)'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a,a)' ) ' '

      seed = 123456789

      do i = 1, 10
        c1 = c4_uniform_01 ( seed )
        c2 = c4_acos ( c1 )
        c3 = c4_cos ( c2 )
        write ( *, '(2x,2f12.6,2x,2f12.6,2x,2f12.6)' ) c1, c2, c3
      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests C4_ACOSH and C4_COSH.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex c1
      complex c2
      complex c3
      complex c4_acosh
      complex c4_cosh
      complex c4_uniform_01
      integer i
      integer seed


      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) 
     &  '  C4_ACOSH computes the inverse hyperbolic cosine.'
      write ( *, '(a)' ) '  C4_COSH computes the hyperbolic cosine.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C4_UNIFORM_01        ' // 
     &  '  C2 = C4_ACOSH(C1)           C3 = COSH(C2)'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  ' ---------------------     ---------------------'
      write ( *, '(a,a)' ) ' '

      seed = 123456789

      do i = 1, 10
        c1 = c4_uniform_01 ( seed )
        c2 = c4_acosh ( c1 )
        c3 = c4_cosh ( c2 )
        write ( *, '(2x,2f12.6,2x,2f12.6,2x,2f12.6)' ) c1, c2, c3
      end do

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests C4_ADD.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex c1
      complex c2
      complex c3
      complex c4_add
      complex c4_uniform_01
      integer i
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  C4_ADD adds two C4s'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C4_UNIFORM_01        ' // 
     &  '  C2=C4_UNIFORM_01          C3 = C1 + C2'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a,a)' ) ' '

      seed = 123456789

      do i = 1, 10
        c1 = c4_uniform_01 ( seed )
        c2 = c4_uniform_01 ( seed )
        c3 = c4_add ( c1, c2 )
        write ( *, '(2x,2f12.6,2x,2f12.6,2x,2f12.6)' ) c1, c2, c3
      end do

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests C4_ARG.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex c1
      real c4_arg
      complex c4_uniform_01
      real r2
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  C4_ARG computes the argument of a C4.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C4_UNIFORM_01        ' // 
     &  '  R2=C4_ARG(C1)'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c4_uniform_01 ( seed )

        r2 = c4_arg ( c1 )

        write ( *, '(2x,2f12.4,2x,f12.4)' ) c1, r2

      end do

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 tests C4_ASIN and C4_SIN.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex c1
      complex c2
      complex c3
      complex c4_asin
      complex c4_sin
      complex c4_uniform_01
      integer i
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) '  C4_ASIN computes the inverse sine.'
      write ( *, '(a)' ) '  C4_SIN computes the sine.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C4_UNIFORM_01        ' // 
     &  '  C2 = C4_ASIN(C1)          C3 = C4_SIN(C2)'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a,a)' ) ' '

      seed = 123456789

      do i = 1, 10
        c1 = c4_uniform_01 ( seed )
        c2 = c4_asin ( c1 )
        c3 = c4_sin ( c2 )
        write ( *, '(2x,2f12.6,2x,2f12.6,2x,2f12.6)' ) c1, c2, c3
      end do

      return
      end
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07 tests C4_ASINH and C4_SINH.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex c1
      complex c2
      complex c3
      complex c4_asinh
      complex c4_sinh
      complex c4_uniform_01
      integer i
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07'
      write ( *, '(a)' ) 
     &  '  C4_ASINH computes the inverse hyperbolic sine.'
      write ( *, '(a)' ) '  C4_SINH computes the hyperbolic sine.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C4_UNIFORM_01        ' // 
     &  '  C2 = C4_ASINH(C1)         C3 = C4_SINH(C2)'
      write ( *, '(a)' ) '     ---------------------     ' //
     &  '---------------------     ---------------------'
      write ( *, '(a,a)' ) ' '

      seed = 123456789

      do i = 1, 10
        c1 = c4_uniform_01 ( seed )
        c2 = c4_asinh ( c1 )
        c3 = c4_sinh ( c2 )
        write ( *, '(2x,2f12.6,2x,2f12.6,2x,2f12.6)' ) c1, c2, c3
      end do

      return
      end
      subroutine test08 ( )

c*********************************************************************72
c
cc TEST08 tests C4_ATANH and C4_TANH.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex c1
      complex c2
      complex c3
      complex c4_atan
      complex c4_tan
      complex c4_uniform_01
      integer i
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08'
      write ( *, '(a)' ) '  C4_ATAN computes the inverse tangent.'
      write ( *, '(a)' ) '  C4_TAN computes the tangent.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C4_UNIFORM_01        ' // 
     &  '  C2 = C4_ATAN(C1)          C3 = C4_TAN(C2)'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a,a)' ) ' '

      seed = 123456789

      do i = 1, 10
        c1 = c4_uniform_01 ( seed )
        c2 = c4_atan ( c1 )
        c3 = c4_tan ( c2 )
        write ( *, '(2x,2f12.6,2x,2f12.6,2x,2f12.6)' ) c1, c2, c3
      end do

      return
      end
      subroutine test09 ( )

c*********************************************************************72
c
cc TEST09 tests C4_ATANH and C4_TANH.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex c1
      complex c2
      complex c3
      complex c4_atanh
      complex c4_tanh
      complex c4_uniform_01
      integer i
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST09'
      write ( *, '(a)' ) 
     &  '  C4_ATANH computes the inverse hyperbolic tangent.'
      write ( *, '(a)' ) '  C4_TANH computes the hyperbolic tangent.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C4_UNIFORM_01        ' // 
     &  '  C2 = C4_ATANH(C1)         C3 = C4_TANH(C2)'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a,a)' ) ' '

      seed = 123456789

      do i = 1, 10
        c1 = c4_uniform_01 ( seed )
        c2 = c4_atanh ( c1 )
        c3 = c4_tanh ( c2 )
        write ( *, '(2x,2f12.6,2x,2f12.6,2x,2f12.6)' ) c1, c2, c3
      end do

      return
      end
      subroutine test10 ( )

c*********************************************************************72
c
cc TEST10 tests C4_CONJ.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex c1
      complex c2
      complex c3
      complex c4_conj
      complex c4_uniform_01
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST10'
      write ( *, '(a)' ) '  C4_CONJ computes the conjugate of a C4.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C4_UNIFORM_01        ' // 
     &  '  C2=C4_CONJ(C1)            C3=CONJG(C1)'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c4_uniform_01 ( seed )
        c2 = c4_conj ( c1 )
        c3 = conjg ( c1 )

        write ( *, '(2x,2f12.4,2x,2f12.4,2x,2f12.4)' ) c1, c2, c3

      end do

      return
      end
      subroutine test11 ( )

c*********************************************************************72
c
cc TEST11 tests C4_COS.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex c1
      complex c2
      complex c3
      complex c4_cos
      complex c4_uniform_01
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST11'
      write ( *, '(a)' ) '  C4_COS computes the cosine of a C4.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        C1=C4_UNIFORM_01       ' // 
     &  '   C2=C4_COS(C1)             C3=CCOS(C1)'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c4_uniform_01 ( seed )
        c2 = c4_cos ( c1 )
        c3 = ccos ( c1 )

        write ( *, '(2x,2f12.4,2x,2f12.4,2x,2f12.4)' ) c1, c2, c3

      end do

      return
      end
      subroutine test12 ( )

c*********************************************************************72
c
cc TEST12 tests C4_COSH.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex c1
      complex c2
      complex c4_cosh
      complex c4_uniform_01
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST12'
      write ( *, '(a)' ) 
     &  '  C4_COSH computes the hyperbolic cosine of a C4.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C4_UNIFORM_01        ' // 
     &  '  C2=C4_COSH(C1) '
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c4_uniform_01 ( seed )
        c2 = c4_cosh ( c1 )

        write ( *, '(2x,2f12.4,2x,2f12.4,2x,2f12.4)' ) c1, c2

      end do

      return
      end
      subroutine test13 ( )

c*********************************************************************72
c
cc TEST13 tests C4_CUBE_ROOT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex c1
      complex c2
      complex c3
      complex c4_cube_root
      complex c4_uniform_01
      real power
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST13'
      write ( *, '(a)' ) 
     &  '  C4_CUBE_ROOT computes the principal cube root of a C4.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        C1=C4_UNIFORM_01       ' // 
     &  '   C2=C4_CUBE_ROOT(C1)       C3=C2*C2*C2'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c4_uniform_01 ( seed )
        c2 = c4_cube_root ( c1 )
        c3 = c2 * c2 * c2
        write ( *, '(2x,2f12.4,2x,2f12.4,2x,2f12.4)' ) c1, c2, c3

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        C1=C4_UNIFORM_01       ' // 
     &  '   C2=C1**(1.0/3.0)          C3=C2*C2*C2'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      power = 1.0E+00 / 3.0E+00

      do test = 1, test_num

        c1 = c4_uniform_01 ( seed )
        c2 = c1**power
        c3 = c2 * c2 * c2
        write ( *, '(2x,2f12.4,2x,2f12.4,2x,2f12.4)' ) c1, c2, c3

      end do

      return
      end
      subroutine test14 ( )

c*********************************************************************72
c
cc TEST14 tests C4_DIV.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex c1
      complex c2
      complex c3
      complex c4_div
      complex c4_uniform_01
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST14'
      write ( *, '(a)' ) '  C4_DIV computes C3 = C1 / C2.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        C1=C4_UNIFORM_01       ' // 
     &  '   C2=C4_UNIFORM_01          C3=C1/C2'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c4_uniform_01 ( seed )
        c2 = c4_uniform_01 ( seed )
        c3 = c4_div ( c1, c2 )

        write ( *, '(2x,2f12.4,2x,2f12.4,2x,2f12.4)' ) c1, c2, c3

      end do

      return
      end
      subroutine test15 ( )

c*********************************************************************72
c
cc TEST15 tests C4_EXP.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex c1
      complex c2
      complex c3
      complex c4_exp
      complex c4_uniform_01
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST15'
      write ( *, '(a)' ) '  C4_EXP computes exp ( Z ).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        C1=C4_UNIFORM_01       ' // 
     &  '   C2=C4_EXP(C1)             C3=CEXP(C1)'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c4_uniform_01 ( seed )
        c2 = c4_exp ( c1 )
        c3 = cexp ( c1 )

        write ( *, '(2x,2f12.4,2x,2f12.4,2x,2f12.4)' ) c1, c2, c3

      end do

      return
      end
      subroutine test16 ( )

c*********************************************************************72
c
cc TEST16 tests C4_INV.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex c1
      complex c2
      complex c3
      complex c4_inv
      complex c4_uniform_01
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST16'
      write ( *, '(a)' ) '  C4_INV computes C2 = 1 / C1.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        C1=C4_UNIFORM_01       ' // 
     &  '   C2=C4_INV(C1)             C3=C4_INV(C2)'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c4_uniform_01 ( seed )
        c2 = c4_inv ( c1 )
        c3 = c4_inv ( c2 )

        write ( *, '(2x,2f12.4,2x,2f12.4,2x,2f12.4)' ) c1, c2, c3

      end do

      return
      end
      subroutine test17 ( )

c*********************************************************************72
c
cc TEST17 tests C4_LOG.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex c1
      complex c2
      complex c3
      complex c4_exp
      complex c4_log
      complex c4_uniform_01
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST17'
      write ( *, '(a)' ) '  C4_LOG computes log ( Z ).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        C1=C4_UNIFORM_01       ' // 
     &  '   C2=C4_LOG(C1)             C3=C4_EXP(C2)'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c4_uniform_01 ( seed )
        c2 = c4_log ( c1 )
        c3 = c4_exp ( c2 )

        write ( *, '(2x,2f12.4,2x,2f12.4,2x,2f12.4)' ) c1, c2, c3

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        C1=C4_UNIFORM_01       ' // 
     &  '   C2=CLOG(C1)              C3=CEXP(C2)'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c4_uniform_01 ( seed )
        c2 = clog ( c1 )
        c3 = cexp ( c2 )

        write ( *, '(2x,2f12.4,2x,2f12.4,2x,2f12.4)' ) c1, c2, c3

      end do

      return
      end
      subroutine test18 ( )

c*********************************************************************72
c
cc TEST18 tests C4_MAG.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex c1
      real c4_mag
      complex c4_uniform_01
      real r2
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST18'
      write ( *, '(a)' ) '  C4_MAG computes the magnitude of a C4.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C4_UNIFORM_01        ' // 
     &  '  R2=C4_MAG(C1)'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c4_uniform_01 ( seed )
        r2 = c4_mag ( c1 )

        write ( *, '(2x,2f12.4,2x,f12.4)' ) c1, r2

      end do

      return
      end
      subroutine test19 ( )

c*********************************************************************72
c
cc TEST19 tests C4_MUL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex c1
      complex c2
      complex c3
      complex c4_mul
      complex c4_uniform_01
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST19'
      write ( *, '(a)' ) '  C4_MUL computes C3 = C1 * C2.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        C1=C4_UNIFORM_01       ' // 
     &  '   C2=C4_UNIFORM_01          C3=C1*C2'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c4_uniform_01 ( seed )
        c2 = c4_uniform_01 ( seed )
        c3 = c4_mul ( c1, c2 )

        write ( *, '(2x,2f12.4,2x,2f12.4,2x,2f12.4)' ) c1, c2, c3

      end do

      return
      end
      subroutine test20 ( )

c*********************************************************************72
c
cc TEST20 tests C4_NORMAL_01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex c1
      complex c4_normal_01
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 20 )

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST20'
      write ( *, '(a)' ) '  C4_NORMAL_01 generates unit pseudonormal'
      write ( *, '(a)' ) '    complex values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C4_NORMAL_01'
      write ( *, '(a)' ) '     ---------------------'
      write ( *, '(a)' ) ' '

      do test = 1, test_num

        c1 = c4_normal_01 ( seed )
        write ( *, '(2x,2g14.6)' ) c1

      end do

      return
      end
      subroutine test21 ( )

c*********************************************************************72
c
cc TEST21 tests C4_SIN.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex c1
      complex c2
      complex c3
      complex c4_sin
      complex c4_uniform_01
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST21'
      write ( *, '(a)' ) '  C4_SIN computes the sine of a C4.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        C1=C4_UNIFORM_01       ' // 
     &  '   C2=C4_SIN(C1)             C3=CSIN(C1)'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c4_uniform_01 ( seed )
        c2 = c4_sin ( c1 )
        c3 = csin ( c1 )

        write ( *, '(2x,2f12.4,2x,2f12.4,2x,2f12.4)' ) c1, c2, c3

      end do

      return
      end
      subroutine test22 ( )

c*********************************************************************72
c
cc TEST22 tests C4_SINH.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex c1
      complex c2
      complex c4_sinh
      complex c4_uniform_01
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST22'
      write ( *, '(a)' ) 
     &  '  C4_SINH computes the hyperbolic sine of a C4.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C4_UNIFORM_01        ' // 
     &  '  C2=C4_SINH(C1) '
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c4_uniform_01 ( seed )
        c2 = c4_sinh ( c1 )

        write ( *, '(2x,2f12.4,2x,2f12.4)' ) c1, c2

      end do

      return
      end
      subroutine test23 ( )

c*********************************************************************72
c
cc TEST23 tests C4_SQRT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex c1
      complex c2
      complex c3
      complex c4_sqrt
      complex c4_uniform_01
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST23'
      write ( *, '(a)' ) 
     &  '  C4_SQRT computes the principal square root of a C4.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        C1=C4_UNIFORM_01       ' // 
     &  '   C2=C4_SQRT(C1)            C3=C2*C2'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c4_uniform_01 ( seed )
        c2 = c4_sqrt ( c1 )
        c3 = c2 * c2

        write ( *, '(2x,2f12.4,2x,2f12.4,2x,2f12.4)' ) c1, c2, c3

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        C1=C4_UNIFORM_01       ' // 
     &  '   C2=CSQRT(C1)              C3=C2*C2'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c4_uniform_01 ( seed )
        c2 = csqrt ( c1 )
        c3 = c2 * c2

        write ( *, '(2x,2f12.4,2x,2f12.4,2x,2f12.4)' ) c1, c2, c3

      end do

      return
      end
      subroutine test24 ( )

c*********************************************************************72
c
cc TEST24 tests C4_SUB.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex c4_sub
      complex c4_uniform_01
      integer i
      integer seed
      complex c1
      complex c2
      complex c3

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST24'
      write ( *, '(a)' ) '  C4_SUB subtracts two C4s'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C4_UNIFORM_01        ' // 
     &  '  C2=C4_UNIFORM_01          C3 = C1 - C2'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a,a)' ) ' '

      seed = 123456789

      do i = 1, 10
        c1 = c4_uniform_01 ( seed )
        c2 = c4_uniform_01 ( seed )
        c3 = c4_sub ( c1, c2 )
        write ( *, '(2x,2f12.6,2x,2f12.6,2x,2f12.6)' ) c1, c2, c3
      end do

      return
      end
      subroutine test25 ( )

c*********************************************************************72
c
cc TEST25 tests C4_TAN.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex c1
      complex c2
      complex c4_tan
      complex c4_uniform_01
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST25'
      write ( *, '(a)' ) '  C4_TAN computes the tangent of a C4.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C4_UNIFORM_01        ' // 
     &  '  C2=C4_TAN(C1)'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c4_uniform_01 ( seed )
        c2 = c4_tan ( c1 )

        write ( *, '(2x,2f12.4,2x,2f12.4)' ) c1, c2

      end do

      return
      end
      subroutine test26 ( )

c*********************************************************************72
c
cc TEST26 tests C4_TO_CARTESIAN and CARTESIAN_TO_C4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex c1
      complex c3
      complex c4_uniform_01
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )
      real x2
      real y2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST26'
      write ( *, '(a)' ) '  C4_TO_CARTESIAN converts C4 to (X,Y).'
      write ( *, '(a)' ) '  CARTESIAN_TO_C4 converts (X,Y) to C4.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        C1=C4_UNIFORM_01   ' // 
     &  ' (X2,Y2)=C4_TO_CARTESIAN(C1)     C3=CARTESIAN_TO_C4(X2,Y2)'
      write ( *, '(a)' ) '     --------------------- ' // 
     &  '    ---------------------     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c4_uniform_01 ( seed )
        call c4_to_cartesian ( c1, x2, y2 )
        call cartesian_to_c4 ( x2, y2, c3 )

        write ( *, '(2x,2f12.4,2x,2f12.4,2x,2f12.4)' ) c1, x2, y2, c3

      end do

      return
      end
      subroutine test27 ( )

c*********************************************************************72
c
cc TEST27 tests C4_TO_POLAR and POLAR_TO_C4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex c1
      complex c3
      complex c4_uniform_01
      real r2
      integer seed
      real t2
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST27'
      write ( *, '(a)' ) '  C4_TO_POLAR converts C4 to (R,T).'
      write ( *, '(a)' ) '  POLAR_TO_C4 converts (R,T) to C4.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        C1=C4_UNIFORM_01       ' // 
     &  '(R2,T2)=C4_TO_POLAR(C1)     C3=POLAR_TO_C4(R2,T2)'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c4_uniform_01 ( seed )
        call c4_to_polar ( c1, r2, t2 )
        call polar_to_c4 ( r2, t2, c3 )

        write ( *, '(2x,2f12.4,2x,2f12.4,2x,2f12.4)' ) c1, r2, t2, c3

      end do

      return
      end
      subroutine test28 ( )

c*********************************************************************72
c
cc TEST28 tests C4_UNIFORM_01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex c1
      complex c4_uniform_01
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST28'
      write ( *, '(a)' ) 
     &  '  C4_UNIFORM_01 returns a uniformly random "unit" C4.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C4_UNIFORM_01'
      write ( *, '(a)' ) '     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c4_uniform_01 ( seed )

        write ( *, '(2x,2f12.4)' ) c1

      end do

      return
      end
      subroutine test29 ( )

c*********************************************************************72
c
cc TEST29 tests C4MAT_UNIFORM_01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 September 2006
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

      complex a(m,n)
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST29'
      write ( *, '(a)' ) 
     &  '  C4MAT_UNIFORM_01 computes a "random" complex matrix.'

      seed = 123456789

      call c4mat_uniform_01 ( m, n, seed, a )

      call c4mat_print ( m, n, a, '  The matrix:' )

      return
      end
      subroutine test30 ( )

c*********************************************************************72
c
cc TEST30 tests C4VEC_INDICATOR;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 September 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      complex a(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST30'
      write ( *, '(a)' ) 
     &  '  C4VEC_INDICATOR sets A = (1-1i,2-2i,...,N-Ni)'

      call c4vec_indicator ( n, a )

      call c4vec_print ( n, a, '  The "indicator" vector:' )

      return
      end
      subroutine test31 ( )

c*********************************************************************72
c
cc TEST31 tests C4VEC_SPIRAL;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 13 )

      complex c(n)
      complex c1
      complex c2
      integer m

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST31'
      write ( *, '(a)' ) '  C4VEC_SPIRAL returns N points on a spiral'
      write ( *, '(a)' ) '  which includes M complete turns.'

      m = 1
      c1 = cmplx ( 5.0E+00, 0.0E+00 )
      c2 = cmplx ( 3.0E+00, 0.0E+00 )

      call c4vec_spiral ( n, m, c1, c2, c )

      call c4vec_print ( n, c, '  The spiral points:' )

      return
      end
      subroutine test32 ( )

c*********************************************************************72
c
cc TEST32 tests C4VEC_UNITY;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      complex a(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST32'
      write ( *, '(a)' ) '  C4VEC_UNITY returns the N roots of unity'

      call c4vec_unity ( n, a )

      call c4vec_print ( n, a, '  The N roots of unity:' )

      return
      end
