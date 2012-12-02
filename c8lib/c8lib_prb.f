      program main

c*********************************************************************72
c
cc MAIN is the main program for C8LIB_PRB.
c
c  Discussion:
c
c    C8LIB_PRB tests routines from the C8LIB library.
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
      write ( *, '(a)' ) 'C8LIB_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the C8LIB library.'

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
      write ( *, '(a)' ) 'C8LIB_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests C8_ABS.
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

      double complex c1
      double precision c8_abs
      double complex c8_uniform_01
      integer i
      double precision r2
      double precision r3
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  C8_ABS computes the absolute value of a C8.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C8_UNIFORM_01        ' // 
     &  '  R2=C8_ABS(C1)             R3=CDABS(C1)'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do i = 1, 10
        c1 = c8_uniform_01 ( seed )
        r2 = c8_abs ( c1 )
        r3 = cdabs ( c1 )
        write ( *, '(2x,2f12.6,2x,f12.6,12x,2x,f12.6)' ) c1, r2, r3
      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests C8_ACOS and C8_COS.
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

      double complex c1
      double complex c2
      double complex c3
      double complex c8_acos
      double complex c8_cos
      double complex c8_uniform_01
      integer i
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  C8_ACOS computes the inverse cosine.'
      write ( *, '(a)' ) '  C8_COS computes the cosine.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C8_UNIFORM_01        ' // 
     &  '  C2 = C8_ACOS(C1)           C3 = C8_COS(C2)'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a,a)' ) ' '

      seed = 123456789

      do i = 1, 10
        c1 = c8_uniform_01 ( seed )
        c2 = c8_acos ( c1 )
        c3 = c8_cos ( c2 )
        write ( *, '(2x,2f12.6,2x,2f12.6,2x,2f12.6)' ) c1, c2, c3
      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests C8_ACOSH and C8_COSH.
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

      double complex c1
      double complex c2
      double complex c3
      double complex c8_acosh
      double complex c8_cosh
      double complex c8_uniform_01
      integer i
      integer seed


      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) 
     &  '  C8_ACOSH computes the inverse hyperbolic cosine.'
      write ( *, '(a)' ) '  C8_COSH computes the hyperbolic cosine.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C8_UNIFORM_01        ' // 
     &  '  C2 = C8_ACOSH(C1)           C3 = COSH(C2)'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  ' ---------------------     ---------------------'
      write ( *, '(a,a)' ) ' '

      seed = 123456789

      do i = 1, 10
        c1 = c8_uniform_01 ( seed )
        c2 = c8_acosh ( c1 )
        c3 = c8_cosh ( c2 )
        write ( *, '(2x,2f12.6,2x,2f12.6,2x,2f12.6)' ) c1, c2, c3
      end do

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests C8_ADD.
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

      double complex c1
      double complex c2
      double complex c3
      double complex c8_add
      double complex c8_uniform_01
      integer i
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  C8_ADD adds two C8s'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C8_UNIFORM_01        ' // 
     &  '  C2=C8_UNIFORM_01          C3 = C1 + C2'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a,a)' ) ' '

      seed = 123456789

      do i = 1, 10
        c1 = c8_uniform_01 ( seed )
        c2 = c8_uniform_01 ( seed )
        c3 = c8_add ( c1, c2 )
        write ( *, '(2x,2f12.6,2x,2f12.6,2x,2f12.6)' ) c1, c2, c3
      end do

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests C8_ARG.
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

      double complex c1
      double precision c8_arg
      double complex c8_uniform_01
      double precision r2
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  C8_ARG computes the argument of a C8.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C8_UNIFORM_01        ' // 
     &  '  R2=C8_ARG(C1)'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c8_uniform_01 ( seed )

        r2 = c8_arg ( c1 )

        write ( *, '(2x,2f12.4,2x,f12.4)' ) c1, r2

      end do

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 tests C8_ASIN and C8_SIN.
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

      double complex c1
      double complex c2
      double complex c3
      double complex c8_asin
      double complex c8_sin
      double complex c8_uniform_01
      integer i
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) '  C8_ASIN computes the inverse sine.'
      write ( *, '(a)' ) '  C8_SIN computes the sine.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C8_UNIFORM_01        ' // 
     &  '  C2 = C8_ASIN(C1)          C3 = C8_SIN(C2)'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a,a)' ) ' '

      seed = 123456789

      do i = 1, 10
        c1 = c8_uniform_01 ( seed )
        c2 = c8_asin ( c1 )
        c3 = c8_sin ( c2 )
        write ( *, '(2x,2f12.6,2x,2f12.6,2x,2f12.6)' ) c1, c2, c3
      end do

      return
      end
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07 tests C8_ASINH and C8_SINH.
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

      double complex c1
      double complex c2
      double complex c3
      double complex c8_asinh
      double complex c8_sinh
      double complex c8_uniform_01
      integer i
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07'
      write ( *, '(a)' ) 
     &  '  C8_ASINH computes the inverse hyperbolic sine.'
      write ( *, '(a)' ) '  C8_SINH computes the hyperbolic sine.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C8_UNIFORM_01        ' // 
     &  '  C2 = C8_ASINH(C1)         C3 = C8_SINH(C2)'
      write ( *, '(a)' ) '     ---------------------     ' //
     &  '---------------------     ---------------------'
      write ( *, '(a,a)' ) ' '

      seed = 123456789

      do i = 1, 10
        c1 = c8_uniform_01 ( seed )
        c2 = c8_asinh ( c1 )
        c3 = c8_sinh ( c2 )
        write ( *, '(2x,2f12.6,2x,2f12.6,2x,2f12.6)' ) c1, c2, c3
      end do

      return
      end
      subroutine test08 ( )

c*********************************************************************72
c
cc TEST08 tests C8_ATANH and C8_TANH.
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

      double complex c1
      double complex c2
      double complex c3
      double complex c8_atan
      double complex c8_tan
      double complex c8_uniform_01
      integer i
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08'
      write ( *, '(a)' ) '  C8_ATAN computes the inverse tangent.'
      write ( *, '(a)' ) '  C8_TAN computes the tangent.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C8_UNIFORM_01        ' // 
     &  '  C2 = C8_ATAN(C1)          C3 = C8_TAN(C2)'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a,a)' ) ' '

      seed = 123456789

      do i = 1, 10
        c1 = c8_uniform_01 ( seed )
        c2 = c8_atan ( c1 )
        c3 = c8_tan ( c2 )
        write ( *, '(2x,2f12.6,2x,2f12.6,2x,2f12.6)' ) c1, c2, c3
      end do

      return
      end
      subroutine test09 ( )

c*********************************************************************72
c
cc TEST09 tests C8_ATANH and C8_TANH.
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

      double complex c1
      double complex c2
      double complex c3
      double complex c8_atanh
      double complex c8_tanh
      double complex c8_uniform_01
      integer i
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST09'
      write ( *, '(a)' ) 
     &  '  C8_ATANH computes the inverse hyperbolic tangent.'
      write ( *, '(a)' ) '  C8_TANH computes the hyperbolic tangent.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C8_UNIFORM_01        ' // 
     &  '  C2 = C8_ATANH(C1)         C3 = C8_TANH(C2)'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a,a)' ) ' '

      seed = 123456789

      do i = 1, 10
        c1 = c8_uniform_01 ( seed )
        c2 = c8_atanh ( c1 )
        c3 = c8_tanh ( c2 )
        write ( *, '(2x,2f12.6,2x,2f12.6,2x,2f12.6)' ) c1, c2, c3
      end do

      return
      end
      subroutine test10 ( )

c*********************************************************************72
c
cc TEST10 tests C8_CONJ.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    14 December 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double complex c1
      double complex c2
      double complex c3
      double complex c8_conj
      double complex c8_uniform_01
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST10'
      write ( *, '(a)' ) '  C8_CONJ computes the conjugate of a C8.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C8_UNIFORM_01        ' // 
     &  '  C2=C8_CONJ(C1)            C3=DCONJG(C1)'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c8_uniform_01 ( seed )
        c2 = c8_conj ( c1 )
        c3 = dconjg ( c1 )

        write ( *, '(2x,2f12.4,2x,2f12.4,2x,2f12.4)' ) c1, c2, c3

      end do

      return
      end
      subroutine test11 ( )

c*********************************************************************72
c
cc TEST11 tests C8_COS.
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

      double complex c1
      double complex c2
      double complex c3
      double complex c8_cos
      double complex c8_uniform_01
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST11'
      write ( *, '(a)' ) '  C8_COS computes the cosine of a C8.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        C1=C8_UNIFORM_01       ' // 
     &  '   C2=C8_COS(C1)             C3=CDCOS(C1)'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c8_uniform_01 ( seed )
        c2 = c8_cos ( c1 )
        c3 = cdcos ( c1 )

        write ( *, '(2x,2f12.4,2x,2f12.4,2x,2f12.4)' ) c1, c2, c3

      end do

      return
      end
      subroutine test12 ( )

c*********************************************************************72
c
cc TEST12 tests C8_COSH.
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

      double complex c1
      double complex c2
      double complex c8_cosh
      double complex c8_uniform_01
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST12'
      write ( *, '(a)' ) 
     &  '  C8_COSH computes the hyperbolic cosine of a C8.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C8_UNIFORM_01        ' // 
     &  '  C2=C8_COSH(C1) '
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c8_uniform_01 ( seed )
        c2 = c8_cosh ( c1 )

        write ( *, '(2x,2f12.4,2x,2f12.4,2x,2f12.4)' ) c1, c2

      end do

      return
      end
      subroutine test13 ( )

c*********************************************************************72
c
cc TEST13 tests C8_CUBE_ROOT.
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

      double complex c1
      double complex c2
      double complex c3
      double complex c8_cube_root
      double complex c8_uniform_01
      double precision power
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST13'
      write ( *, '(a)' ) 
     &  '  C8_CUBE_ROOT computes the principal cube root of a C8.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        C1=C8_UNIFORM_01       ' // 
     &  '   C2=C8_CUBE_ROOT(C1)       C3=C2*C2*C2'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c8_uniform_01 ( seed )
        c2 = c8_cube_root ( c1 )
        c3 = c2 * c2 * c2
        write ( *, '(2x,2f12.4,2x,2f12.4,2x,2f12.4)' ) c1, c2, c3

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        C1=C8_UNIFORM_01       ' // 
     &  '   C2=C1**(1.0/3.0)          C3=C2*C2*C2'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      power = 1.0D+00 / 3.0D+00

      do test = 1, test_num

        c1 = c8_uniform_01 ( seed )
        c2 = c1**power
        c3 = c2 * c2 * c2
        write ( *, '(2x,2f12.4,2x,2f12.4,2x,2f12.4)' ) c1, c2, c3

      end do

      return
      end
      subroutine test14 ( )

c*********************************************************************72
c
cc TEST14 tests C8_DIV.
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

      double complex c1
      double complex c2
      double complex c3
      double complex c8_div
      double complex c8_uniform_01
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST14'
      write ( *, '(a)' ) '  C8_DIV computes C3 = C1 / C2.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        C1=C8_UNIFORM_01       ' // 
     &  '   C2=C8_UNIFORM_01          C3=C1/C2'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c8_uniform_01 ( seed )
        c2 = c8_uniform_01 ( seed )
        c3 = c8_div ( c1, c2 )

        write ( *, '(2x,2f12.4,2x,2f12.4,2x,2f12.4)' ) c1, c2, c3

      end do

      return
      end
      subroutine test15 ( )

c*********************************************************************72
c
cc TEST15 tests C8_EXP.
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

      double complex c1
      double complex c2
      double complex c3
      double complex c8_exp
      double complex c8_uniform_01
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST15'
      write ( *, '(a)' ) '  C8_EXP computes exp ( Z ).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        C1=C8_UNIFORM_01       ' // 
     &  '   C2=C8_EXP(C1)             C3=CDEXP(C1)'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c8_uniform_01 ( seed )
        c2 = c8_exp ( c1 )
        c3 = cdexp ( c1 )

        write ( *, '(2x,2f12.4,2x,2f12.4,2x,2f12.4)' ) c1, c2, c3

      end do

      return
      end
      subroutine test16 ( )

c*********************************************************************72
c
cc TEST16 tests C8_INV.
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

      double complex c1
      double complex c2
      double complex c3
      double complex c8_inv
      double complex c8_uniform_01
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST16'
      write ( *, '(a)' ) '  C8_INV computes C2 = 1 / C1.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        C1=C8_UNIFORM_01       ' // 
     &  '   C2=C8_INV(C1)             C3=C8_INV(C2)'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c8_uniform_01 ( seed )
        c2 = c8_inv ( c1 )
        c3 = c8_inv ( c2 )

        write ( *, '(2x,2f12.4,2x,2f12.4,2x,2f12.4)' ) c1, c2, c3

      end do

      return
      end
      subroutine test17 ( )

c*********************************************************************72
c
cc TEST17 tests C8_LOG.
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

      double complex c1
      double complex c2
      double complex c3
      double complex c8_exp
      double complex c8_log
      double complex c8_uniform_01
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST17'
      write ( *, '(a)' ) '  C8_LOG computes log ( Z ).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        C1=C8_UNIFORM_01       ' // 
     &  '   C2=C8_LOG(C1)             C3=C8_EXP(C2)'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c8_uniform_01 ( seed )
        c2 = c8_log ( c1 )
        c3 = c8_exp ( c2 )

        write ( *, '(2x,2f12.4,2x,2f12.4,2x,2f12.4)' ) c1, c2, c3

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        C1=C8_UNIFORM_01       ' // 
     &  '   C2=CDLOG(C1)              C3=CDEXP(C2)'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c8_uniform_01 ( seed )
        c2 = cdlog ( c1 )
        c3 = cdexp ( c2 )

        write ( *, '(2x,2f12.4,2x,2f12.4,2x,2f12.4)' ) c1, c2, c3

      end do

      return
      end
      subroutine test18 ( )

c*********************************************************************72
c
cc TEST18 tests C8_MAG.
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

      double complex c1
      double precision c8_mag
      double complex c8_uniform_01
      double precision r2
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST18'
      write ( *, '(a)' ) '  C8_MAG computes the magnitude of a C8.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C8_UNIFORM_01        ' // 
     &  '  R2=C8_MAG(C1)'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c8_uniform_01 ( seed )
        r2 = c8_mag ( c1 )

        write ( *, '(2x,2f12.4,2x,f12.4)' ) c1, r2

      end do

      return
      end
      subroutine test19 ( )

c*********************************************************************72
c
cc TEST19 tests C8_MUL.
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

      double complex c1
      double complex c2
      double complex c3
      double complex c8_mul
      double complex c8_uniform_01
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST19'
      write ( *, '(a)' ) '  C8_MUL computes C3 = C1 * C2.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        C1=C8_UNIFORM_01       ' // 
     &  '   C2=C8_UNIFORM_01          C3=C1*C2'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c8_uniform_01 ( seed )
        c2 = c8_uniform_01 ( seed )
        c3 = c8_mul ( c1, c2 )

        write ( *, '(2x,2f12.4,2x,2f12.4,2x,2f12.4)' ) c1, c2, c3

      end do

      return
      end
      subroutine test20 ( )

c*********************************************************************72
c
cc TEST20 tests C8_NORMAL_01.
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

      double complex c1
      double complex c8_normal_01
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 20 )

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST20'
      write ( *, '(a)' ) '  C8_NORMAL_01 generates unit pseudonormal'
      write ( *, '(a)' ) '    complex values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C8_NORMAL_01'
      write ( *, '(a)' ) '     ---------------------'
      write ( *, '(a)' ) ' '

      do test = 1, test_num

        c1 = c8_normal_01 ( seed )
        write ( *, '(2x,2g14.6)' ) c1

      end do

      return
      end
      subroutine test21 ( )

c*********************************************************************72
c
cc TEST21 tests C8_SIN.
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

      double complex c1
      double complex c2
      double complex c3
      double complex c8_sin
      double complex c8_uniform_01
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST21'
      write ( *, '(a)' ) '  C8_SIN computes the sine of a C8.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        C1=C8_UNIFORM_01       ' // 
     &  '   C2=C8_SIN(C1)             C3=CDSIN(C1)'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c8_uniform_01 ( seed )
        c2 = c8_sin ( c1 )
        c3 = cdsin ( c1 )

        write ( *, '(2x,2f12.4,2x,2f12.4,2x,2f12.4)' ) c1, c2, c3

      end do

      return
      end
      subroutine test22 ( )

c*********************************************************************72
c
cc TEST22 tests C8_SINH.
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

      double complex c1
      double complex c2
      double complex c8_sinh
      double complex c8_uniform_01
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST22'
      write ( *, '(a)' ) 
     &  '  C8_SINH computes the hyperbolic sine of a C8.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C8_UNIFORM_01        ' // 
     &  '  C2=C8_SINH(C1) '
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c8_uniform_01 ( seed )
        c2 = c8_sinh ( c1 )

        write ( *, '(2x,2f12.4,2x,2f12.4)' ) c1, c2

      end do

      return
      end
      subroutine test23 ( )

c*********************************************************************72
c
cc TEST23 tests C8_SQRT.
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

      double complex c1
      double complex c2
      double complex c3
      double complex c8_sqrt
      double complex c8_uniform_01
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST23'
      write ( *, '(a)' ) 
     &  '  C8_SQRT computes the principal square root of a C8.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        C1=C8_UNIFORM_01       ' // 
     &  '   C2=C8_SQRT(C1)            C3=C2*C2'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c8_uniform_01 ( seed )
        c2 = c8_sqrt ( c1 )
        c3 = c2 * c2

        write ( *, '(2x,2f12.4,2x,2f12.4,2x,2f12.4)' ) c1, c2, c3

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        C1=C8_UNIFORM_01       ' // 
     &  '   C2=CDSQRT(C1)             C3=C2*C2'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c8_uniform_01 ( seed )
        c2 = cdsqrt ( c1 )
        c3 = c2 * c2

        write ( *, '(2x,2f12.4,2x,2f12.4,2x,2f12.4)' ) c1, c2, c3

      end do

      return
      end
      subroutine test24 ( )

c*********************************************************************72
c
cc TEST24 tests C8_SUB.
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

      double complex c8_sub
      double complex c8_uniform_01
      integer i
      integer seed
      double complex c1
      double complex c2
      double complex c3

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST24'
      write ( *, '(a)' ) '  C8_SUB subtracts two C8s'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C8_UNIFORM_01        ' // 
     &  '  C2=C8_UNIFORM_01          C3 = C1 - C2'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a,a)' ) ' '

      seed = 123456789

      do i = 1, 10
        c1 = c8_uniform_01 ( seed )
        c2 = c8_uniform_01 ( seed )
        c3 = c8_sub ( c1, c2 )
        write ( *, '(2x,2f12.6,2x,2f12.6,2x,2f12.6)' ) c1, c2, c3
      end do

      return
      end
      subroutine test25 ( )

c*********************************************************************72
c
cc TEST25 tests C8_TAN.
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

      double complex c1
      double complex c2
      double complex c8_tan
      double complex c8_uniform_01
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST25'
      write ( *, '(a)' ) '  C8_TAN computes the tangent of a C8.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C8_UNIFORM_01        ' // 
     &  '  C2=C8_TAN(C1)'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c8_uniform_01 ( seed )
        c2 = c8_tan ( c1 )

        write ( *, '(2x,2f12.4,2x,2f12.4)' ) c1, c2

      end do

      return
      end
      subroutine test26 ( )

c*********************************************************************72
c
cc TEST26 tests C8_TO_CARTESIAN and CARTESIAN_TO_C8.
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

      double complex c1
      double complex c3
      double complex c8_uniform_01
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )
      double precision x2
      double precision y2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST26'
      write ( *, '(a)' ) '  C8_TO_CARTESIAN converts C8 to (X,Y).'
      write ( *, '(a)' ) '  CARTESIAN_TO_C8 converts (X,Y) to C8.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        C1=C8_UNIFORM_01   ' // 
     &  ' (X2,Y2)=C8_TO_CARTESIAN(C1)     C3=CARTESIAN_TO_C8(X2,Y2)'
      write ( *, '(a)' ) '     --------------------- ' // 
     &  '    ---------------------     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c8_uniform_01 ( seed )
        call c8_to_cartesian ( c1, x2, y2 )
        call cartesian_to_c8 ( x2, y2, c3 )

        write ( *, '(2x,2f12.4,2x,2f12.4,2x,2f12.4)' ) c1, x2, y2, c3

      end do

      return
      end
      subroutine test27 ( )

c*********************************************************************72
c
cc TEST27 tests C8_TO_POLAR and POLAR_TO_C8.
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

      double complex c1
      double complex c3
      double complex c8_uniform_01
      double precision r2
      integer seed
      double precision t2
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST27'
      write ( *, '(a)' ) '  C8_TO_POLAR converts C8 to (R,T).'
      write ( *, '(a)' ) '  POLAR_TO_C8 converts (R,T) to C8.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '        C1=C8_UNIFORM_01       ' // 
     &  '(R2,T2)=C8_TO_POLAR(C1)     C3=POLAR_TO_C8(R2,T2)'
      write ( *, '(a)' ) '     ---------------------     ' // 
     &  '---------------------     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c8_uniform_01 ( seed )
        call c8_to_polar ( c1, r2, t2 )
        call polar_to_c8 ( r2, t2, c3 )

        write ( *, '(2x,2f12.4,2x,2f12.4,2x,2f12.4)' ) c1, r2, t2, c3

      end do

      return
      end
      subroutine test28 ( )

c*********************************************************************72
c
cc TEST28 tests C8_UNIFORM_01.
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

      double complex c1
      double complex c8_uniform_01
      integer seed
      integer test
      integer test_num
      parameter ( test_num = 10 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST28'
      write ( *, '(a)' ) 
     &  '  C8_UNIFORM_01 returns a uniformly random "unit" C8.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C1=C8_UNIFORM_01'
      write ( *, '(a)' ) '     ---------------------'
      write ( *, '(a)' ) ' '

      seed = 123456789

      do test = 1, test_num

        c1 = c8_uniform_01 ( seed )

        write ( *, '(2x,2f12.4)' ) c1

      end do

      return
      end
      subroutine test29 ( )

c*********************************************************************72
c
cc TEST29 tests C8MAT_UNIFORM_01.
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

      double complex a(m,n)
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST29'
      write ( *, '(a)' ) 
     &  '  C8MAT_UNIFORM_01 computes a "random" complex matrix.'

      seed = 123456789

      call c8mat_uniform_01 ( m, n, seed, a )

      call c8mat_print ( m, n, a, '  The matrix:' )

      return
      end
      subroutine test30 ( )

c*********************************************************************72
c
cc TEST30 tests C8VEC_INDICATOR;
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

      double complex a(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST30'
      write ( *, '(a)' ) 
     &  '  C8VEC_INDICATOR sets A = (1-1i,2-2i,...,N-Ni)'

      call c8vec_indicator ( n, a )

      call c8vec_print ( n, a, '  The "indicator" vector:' )

      return
      end
      subroutine test31 ( )

c*********************************************************************72
c
cc TEST31 tests C8VEC_SPIRAL;
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

      double complex c(n)
      double complex c1
      double complex c2
      integer m

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST31'
      write ( *, '(a)' ) '  C8VEC_SPIRAL returns N points on a spiral'
      write ( *, '(a)' ) '  which includes M complete turns.'

      m = 1
      c1 = dcmplx ( 5.0D+00, 0.0D+00 )
      c2 = dcmplx ( 3.0D+00, 0.0D+00 )

      call c8vec_spiral ( n, m, c1, c2, c )

      call c8vec_print ( n, c, '  The spiral points:' )

      return
      end
      subroutine test32 ( )

c*********************************************************************72
c
cc TEST32 tests C8VEC_UNITY;
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

      double complex a(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST32'
      write ( *, '(a)' ) '  C8VEC_UNITY returns the N roots of unity'

      call c8vec_unity ( n, a )

      call c8vec_print ( n, a, '  The N roots of unity:' )

      return
      end
