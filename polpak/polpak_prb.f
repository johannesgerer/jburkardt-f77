      program main

c*********************************************************************72
c
cc MAIN is the main program for POLPAK_PRB.
c
c  Discussion:
c
c    POLPAK_PRB calls sample problems for the POLPAK library.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    31 March 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      write ( *, '(a)' ) ' '
      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POLPAK_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the POLPAK library.'

      call test001 ( )
      call test002 ( )
      call test003 ( )
      call test0035 ( )
      call test005 ( )
      call test007 ( )
      call test008 ( )

      call test010 ( )
      call test0102 ( )
      call test0104 ( )
      call test011 ( )
      call test0115 ( )
      call test013 ( )
      call test012 ( )
      call test014 ( )
      call test0141 ( )
      call test0142 ( )
      call test0143 ( )
      call test015 ( )
      call test016 ( )
      call test017 ( )
      call test0175 ( )
      call test018 ( )
      call test0185 ( )
      call test019 ( )

      call test020 ( )
      call test021 ( )
      call test0215 ( )
      call test0216 ( )
      call test0217 ( )
      call test0218 ( )
      call test024 ( )
      call test02405 ( )
      call test0241 ( )
      call test0242 ( )
      call test0243 ( )
      call test01155 ( )
      call test0245 ( )
      call test025 ( )
      call test0255 ( )
      call test026 ( )
      call test0265 ( )
      call test028 ( )
      call test027 ( )

      call test031 ( )
      call test032 ( )
      call test036 ( )
      call test037 ( )
      call test038 ( )

      call test041 ( )
      call test042 ( )
      call test0425 ( )
      call test023 ( )
      call test043 ( )
      call test044 ( )
      call test045 ( )
      call test046 ( )
      call test047 ( )
      call test048 ( )
      call test049 ( )

      call test050 ( )
      call test0505 ( )
      call test051 ( )
      call test052 ( )
      call test054 ( )
      call test055 ( )
      call test0552 ( )
      call test059 ( )
      call test0595 ( )

      call test060 ( )
      call test057 ( )
      call test058 ( )
      call test061 ( )
      call test0615 ( )
      call test062 ( )
      call test0623 ( )
      call test0625 ( )
      call test063 ( )
      call test0635 ( )
      call test064 ( )
      call test065 ( )
      call test066 ( )
      call test0665 ( )
      call test0667 ( )
      call test067 ( )
      call test0675 ( )
      call test004 ( )
      call test006 ( )
      call test022 ( )
      call test0685 ( )
      call test06855 ( )
      call test06856 ( )
      call test069 ( )
      call test0695 ( )
      call test0696 ( )
      call test0697 ( )

      call test070 ( )
      call test071 ( )
      call test072 ( )
      call test073 ( )
      call test074 ( )
      call test076 ( )
      call test077 ( )
      call test0773 ( )
      call test0775 ( )
      call test078 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POLPAK_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test001 ( )

c*********************************************************************72
c
cc TEST001 tests AGM and AGM_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision agm
      double precision b
      double precision fx
      double precision fx2
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST001'
      write ( *, '(a)' ) 
     &  '  AGM computes the arithmetic geometric mean.'
      write ( *, '(a)' ) '  AGM_VALUES returns some exact values.'
      write ( *, '(a)' ) ' ' 
      write ( *, '(a,a)' ) '      A           B          ',
     &  '   AGM                       AGM                   Diff'
      write ( *, '(a,a)' ) '                             ',
     &  '  (Tabulated)                AGM(A,B)'
      write ( *, '(a)' ) ' '
     
      n_data = 0

10    continue

        call agm_values ( n_data, a, b, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = agm ( a, b )

        write ( *, '(2x,f10.6,2x,f10.6,2x,g24.16,2x,g24.16,2x,g10.4)' ) 
     &    a, b, fx, fx2, abs ( fx - fx2 )

      go to 10

20    continue
     
      return
      end
      subroutine test002 ( )

c*********************************************************************72
c
cc TEST002 tests AGUD and GUD.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 October 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision agud
      double precision gamma
      double precision gud
      integer i
      double precision x
      double precision x2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST002'
      write ( *, '(a)' ) '  AGUD computes the inverse Gudermannian;'
      write ( *, '(a)' ) '  GUD computes the Gudermannian.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       X            GUD(X)     AGUD(GUD(X))'
      write ( *, '(a)' ) ' '

      do i = 0, 10
        x = 1.0D+00 + dble ( i ) / 5.0D+00
        gamma = gud ( x )
        x2 = agud ( gamma )
        write ( *, '(2x,3g14.6)' ) x, gamma, x2
      end do

      return
      end
      subroutine test003 ( )

c*********************************************************************72
c
cc TEST003 tests ALIGN_ENUM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 December 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m_max
      parameter ( m_max = 10 )
      integer n_max
      parameter ( n_max = 10 )

      integer align_enum
      integer i
      integer j
      integer table(0:m_max,0:n_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST003'
      write ( *, '(a)' ) '  ALIGN_ENUM counts the number of possible'
      write ( *, '(a)' ) '  alignments of two biological sequences.'

      do i = 0, m_max
        do j = 0, n_max
          table(i,j) = align_enum ( i, j )
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Alignment enumeration table:'
      write ( *, '(a)' ) ' '
      write ( *, '(4x,5i5,6i8)' ) ( i, i = 0, n_max )
      write ( *, '(a)' ) ' '
      do i = 0, m_max
        write ( *, '(2x,i2,5i5,6i8)' ) i, table(i,0:n_max)
      end do

      return
      end
      subroutine test0035 ( )

c*********************************************************************72
c
cc TEST0035 tests ARC_COSINE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 December 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision arc_cosine
      integer ( kind = 4 ) i
      double precision x
      double precision x2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0035'
      write ( *, '(a)' ) '  ARC_COSINE computes the inverse cosine'
      write ( *, '(a)' ) '  of a given value, and chops out of bound '
      write ( *, '(a)' ) '  arguments.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '    X     ARC_COSINE(X)     COS(ARC_COSINE(X))'
      write ( *, '(a)' ) ' '

      do i = -5, 5
        x = 1.0D+00 + dble ( i ) / 5.0D+00
        a = arc_cosine ( x )
        x2 = cos ( a )
        write ( *, '(2x,3g14.6)' ) x, a, x2
      end do

      return
      end
      subroutine test005 ( )

c*********************************************************************72
c
cc TEST005 tests ATAN4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 December 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 8 )

      double precision atan4
      integer test
      double precision x
      double precision x_test(test_num)
      double precision y
      double precision y_test(test_num)

      data x_test /
     &   1.0D+00,  1.0D+00, 0.0D+00, -1.0D+00, 
     &  -1.0D+00, -1.0D+00, 0.0D+00,  1.0D+00 /
      data y_test /
     &  0.0D+00,  1.0D+00,  1.0D+00,  1.0D+00, 
     &  0.0D+00, -1.0D+00, -1.0D+00, -1.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST005'
      write ( *, '(a)' ) 
     &  '  ATAN4 computes the arc-tangent given Y and X;'
      write ( *, '(a)' ) 
     &  '  ATAN2 is the system version of this routine.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '       X             Y          ATAN2(Y,X)   ATAN4(Y,X)'
      write ( *, '(a)' ) ' '

      do test = 1, test_num
        x = x_test(test)
        y = y_test(test)
        write ( *, '(2x,4g14.6)' ) 
     &  x, y, atan2 ( y, x ), atan4 ( y, x )
      end do

      return
      end
      subroutine test007 ( )

c***********************************************************************72
c
cc TEST007 tests BELL and BELL_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 December 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer c
      integer c2(0:10)
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST007'
      write ( *, '(a)' ) '  BELL computes Bell numbers.'
      write ( *, '(a)' ) '  BELL_VALUES returns some exact values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     N  exact C(I)  computed C(I)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bell_values ( n_data, n, c )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        call bell ( n, c2 )

        write ( *, '(2x,i8,2x,i10,2x,i10)' ) n, c, c2(n)

      go to 10

20    continue
 
      return
      end
      subroutine test008 ( )

c*********************************************************************72
c
cc TEST008 tests BENFORD.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 December 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision benford
      integer i

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST008'
      write ( *, '(a)' ) '  BENFORD(I) is the Benford probability of'
      write ( *, '(a)' ) '  the initial digit sequence I.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  I,  BENFORD(I)'
      write ( *, '(a)' ) ' '

      do i = 1, 9
        write ( *, '(2x,i2,2x,g14.6)' )  i, benford(i)
      end do

      return
      end
      subroutine test010 ( )

c*********************************************************************72
c
cc TEST010 tests BERNOULLI_NUMBER and BERNOULLI_NUMBER_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 December 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision c0
      double precision c1(0:30)
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST010'
      write ( *, '(a)' ) 
     &  '  BERNOULLI_NUMBER computes Bernoulli numbers;'
      write ( *, '(a)' ) 
     &  '  BERNOULLI_NUMBER_VALUES returns some exact values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   I      Exact     Bernoulli'
      write ( *, '(a)' ) ' '
  
      n_data = 0

10    continue

        call bernoulli_number_values ( n_data, n, c0 )

        if ( n_data == 0 ) then
          go to 20
        end if

        call bernoulli_number ( n, c1 )

        write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) n, c0, c1(n)

      go to 10

20    continue
 
      return
      end
      subroutine test0102 ( )

c*********************************************************************72
c
cc TEST0102 tests BERNOULLI_NUMBER2 and BERNOULLI_NUMBER_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 December 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision c0
      double precision c1(0:30)
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0102'
      write ( *, '(a)' ) 
     &  '  BERNOULLI_NUMBER2 computes Bernoulli numbers;'
      write ( *, '(a)' ) 
     &  '  BERNOULLI_NUMBER_VALUES returns some exact values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   I      Exact     Bernoulli2'
      write ( *, '(a)' ) ' '
  
      n_data = 0

10    continue

        call bernoulli_number_values ( n_data, n, c0 )

        if ( n_data == 0 ) then
          go to 20
        end if

        call bernoulli_number2 ( n, c1 )
 
        write ( *, '(2x,i4,2g14.6)' ) n, c0, c1(n)

      go to 10

20    continue
 
      return
      end
      subroutine test0104 ( )

c*********************************************************************72
c
cc TEST0104 tests BERNOULLI_NUMBER3 and BERNOULLI_NUMBER_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 December 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision c0
      double precision c1
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0104'
      write ( *, '(a)' ) 
     &  '  BERNOULLI_NUMBER3 computes Bernoulli numbers.'
      write ( *, '(a)' ) 
     &  '  BERNOULLI_NUMBER_VALUES returns some exact values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   I      Exact     BERNOULLI3'
      write ( *, '(a)' ) ' '
  
      n_data = 0

10    continue

        call bernoulli_number_values ( n_data, n, c0 )

        if ( n_data == 0 ) then
          go to 20
        end if

        call bernoulli_number3 ( n, c1 )

        write ( *, '(2x,i4,2g14.6)' ) n, c0, c1

      go to 10

20    continue
 
      return
      end
      subroutine test011 ( )

c*********************************************************************72
c
cc TEST011 tests BERNOULLI_POLY;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision bx
      integer i
      integer n
      parameter ( n = 15 )
      double precision x

      x = 0.2D+00
 
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST011'
      write ( *, '(a)' ) 
     &'  BERNOULLI_POLY evaluates Bernoulli polynomials;'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  X = ', x
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  I          BX'
      write ( *, '(a)' ) ' '
 
      do i = 1, n
        call bernoulli_poly ( i, x, bx )
        write ( *, '(2x,i2,2x,g16.8)' ) i, bx
      end do
 
      return
      end
      subroutine test0115 ( )

c*********************************************************************72
c
cc TEST0115 tests BERNOULLI_POLY2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 July 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision bx
      integer i
      integer n
      parameter ( n = 15 )
      double precision x

      x = 0.2D+00
 
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0115'
      write ( *, '(a)' ) 
     &  '  BERNOULLI_POLY2 evaluates Bernoulli polynomials. '
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  X = ', x
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  I          BX'
      write ( *, '(a)' ) ' '
 
      do i = 1, n
        call bernoulli_poly2 ( i, x, bx )
        write ( *, '(2x,i2,2x,2g16.8)' ) i, bx
      end do
 
      return
      end
      subroutine test013 ( )

c*********************************************************************72
c
cc TEST013 tests BERNSTEIN_POLY and BERNSTEIN_POLY_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision b
      double precision bvec(0:10)
      integer k
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST013:'
      write ( *, '(a)' ) 
     &  '  BERNSTEIN_POLY evaluates the Bernstein polynomials.'
      write ( *, '(a)' ) 
     &  '  BERNSTEIN_POLY_VALUES returns some exact values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   N   K   X   Exact   B(N,K)(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call bernstein_poly_values ( n_data, n, k, x, b )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        call bernstein_poly ( n, x, bvec )

        write ( *, '(2x,i4,i4,f7.4,2g14.6)' ) n, k, x, b, bvec(k)

      go to 10

20    continue

      return
      end
      subroutine test012 ( )

c*********************************************************************72
c
cc TEST012 tests BETA and BETA_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision beta
      double precision fxy
      double precision fxy2
      integer n_data
      double precision x
      double precision y

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST012:'
      write ( *, '(a)' ) '  BETA evaluates the Beta function.'
      write ( *, '(a)' ) '  BETA_VALUES returns some exact values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     X      Y        Exact F       BETA(X,Y)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call beta_values ( n_data, x, y, fxy )

        if ( n_data == 0 ) then
          go to 20
        end if

        fxy2 = beta ( x, y )

        write ( *, '(2x,2f8.4,2g14.6)' ) x, y, fxy, fxy2

      go to 10

20    continue

      return
      end
      subroutine test014 ( )

c*********************************************************************72
c
cc TEST014 tests BPAB.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      double precision a
      double precision b
      double precision bern(0:n)
      integer i
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST014'
      write ( *, '(a)' ) '  BPAB evaluates Bernstein polynomials.'
      write ( *, '(a)' ) ' '

      x = 0.3D+00
      a = 0.0D+00
      b = 1.0D+00
      call bpab ( n, x, a, b, bern )
 
      write ( *, '(a,i4)' ) '  The Bernstein polynomials of degree ', n
      write ( *, '(a,g14.6)' ) '  based on the interval from ', a
      write ( *, '(a,g14.6)' ) '  to ', b
      write ( *, '(a,g14.6)' ) '  evaluated at X = ', x
      write ( *, '(a)' ) ' '
 
      do i = 0, n
        write ( *, '(2x,i4,2x,g14.6)' )  i, bern(i)
      end do
 
      return
      end
      subroutine test0141 ( )

c*********************************************************************72
c
cc TEST0141 tests C4_ACOSH.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 November 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex a
      complex c4_acosh
      complex c4_cosh
      complex c4_uniform_01
      integer i
      integer seed
      complex x
      complex x2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0141'
      write ( *, '(a)' ) '  C4_ACOSH computes the inverse hyperbolic'
      write ( *, '(a)' ) '  cosine of a given value.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) 
     &  '               X                  C4_ACOSH(X)       ',
     &  '       COSH(C4_ACOSH(X))'
      write ( *, '(a,a)' ) 
     &  '     ---------------------     ---------------------',
     &  '     ---------------------'
      write ( *, '(a,a)' ) ' '

      seed = 123456789

      do i = 0, 10
        x = c4_uniform_01 ( seed )
        a = c4_acosh ( x )
        x2 = c4_cosh ( a )
        write ( *, '(2x,2f12.6,2x,2f12.6,2x,2f12.6)' ) x, a, x2
      end do

      return
      end
      subroutine test0142 ( )

c*********************************************************************72
c
cc TEST0142 tests C4_ASINH.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 November 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex a
      complex c4_asinh
      complex c4_sinh
      complex c4_uniform_01
      integer i
      integer seed
      complex x
      complex x2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0142'
      write ( *, '(a)' ) '  C4_ASINH computes the inverse hyperbolic'
      write ( *, '(a)' ) '  sine of a given value.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) 
     &  '               X                  C4_ASINH(X)       ',
     &  '       SINH(C4_ASINH(X))'
      write ( *, '(a,a)' ) 
     &  '     ---------------------     ---------------------',
     &  '     ---------------------'
      write ( *, '(a,a)' ) ' '

      seed = 123456789

      do i = 0, 10
        x = c4_uniform_01 ( seed )
        a = c4_asinh ( x )
        x2 = c4_sinh ( a )
        write ( *, '(2x,2f12.6,2x,2f12.6,2x,2f12.6)' ) x, a, x2
      end do

      return
      end
      subroutine test0143 ( )

c*********************************************************************72
c
cc TEST0143 tests C4_ATANH.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 November 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex a
      complex c4_atanh
      complex c4_tanh
      complex c4_uniform_01
      integer i
      integer seed
      complex x
      complex x2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0143'
      write ( *, '(a)' ) '  C4_ATANH computes the inverse hyperbolic'
      write ( *, '(a)' ) '  tangent of a given value.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) 
     &  '               X                  C4_ATANH(X)       ',
     &  '       TANH(C4_ATANH(X))'
      write ( *, '(a,a)' ) 
     &  '     ---------------------     ---------------------',
     &  '     ---------------------'
      write ( *, '(a,a)' ) ' '

      seed = 123456789

      do i = 0, 10
        x = c4_uniform_01 ( seed )
        a = c4_atanh ( x )
        x2 = c4_tanh ( a )
        write ( *, '(2x,2f12.6,2x,2f12.6,2x,2f12.6)' ) x, a, x2
      end do

      return
      end
      subroutine test015 ( )

c*********************************************************************72
c
cc TEST015 tests CARDAN and CARDAN_POLY_COEF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 10 )

      double precision c(0:n_max)
      double precision cx1
      double precision cx2(0:n_max)
      integer n
      double precision s
      double precision x

      s = 1.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST015'
      write ( *, '(a)' ) '  CARDAN_POLY_COEF returns the coefficients'
      write ( *, '(a)' ) '  of a Cardan polynomial.'
      write ( *, '(a)' ) 
     &  '  CARDAN evaluates a Cardan polynomial directly.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  We use the parameter S = ', s
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Table of polynomial coefficients:'
      write ( *, '(a)' ) ' '

      do n = 0, n_max
        call cardan_poly_coef ( n, s, c )
        write ( *, '(2x,i2,11f7.0)' ) n, c(0:n)
      end do

      s = 0.5D+00
      x = 0.25D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Compare CARDAN_POLY_COEF + R8POLY_VAL_HORNER'
      write ( *, '(a)' ) '  versus CARDAN alone.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Evaluate polynomials at X = ', x
      write ( *, '(a,g14.6)' ) '  We use the parameter S = ', s
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Order, Horner, Direct'
      write ( *, '(a)' ) ' '

      call cardan ( n, x, s, cx2 )

      do n = 0, n_max

        call cardan_poly_coef ( n, s, c )
        call r8poly_val_horner ( n, c, x, cx1 )

        write ( *, '(2x,i2,2g14.6)' ) n, cx1, cx2(n)

      end do

      return
      end
      subroutine test016 ( )

c*********************************************************************72
c
cc TEST016 tests CATALAN and CATALAN_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer c
      integer c2(0:10)
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST016'
      write ( *, '(a)' ) '  CATALAN computes Catalan numbers.'
      write ( *, '(a)' ) '  CATALAN_VALUES returns some exact values.'
      write ( *, '(a)' ) ' ' 
      write ( *, '(a)' ) '  N  exact C(I)  computed C(I)'
      write ( *, '(a)' ) ' '
 
      n_data = 0

10    continue

        call catalan_values ( n_data, n, c )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        call catalan ( n, c2 )
 
        write ( *, '(2x,i4,2i8)' ) n, c, c2(n)

      go to 10
 
20    continue

      return
      end
      subroutine test017 ( )

c*********************************************************************72
c
cc TEST017 tests CATALAN_ROW_NEXT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer c(0:n)
      integer i
      integer ido

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST017'
      write ( *, '(a)' ) 
     &  '  CATALAN_ROW_NEXT computes a row of Catalan''s triangle.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  First, compute row 7:'

      ido = 0
      i = 7
      call catalan_row_next ( ido, i, c )
      write ( *, '(2x,i2,2x,11i6)' ) i, c(0:i)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Now compute rows one at a time:'
      write ( *, '(a)' ) ' '

      ido = 0
 
      do i = 0, n
        call catalan_row_next ( ido, i, c )
        ido = 1
        write ( *, '(2x,i2,2x,11i6)' ) i, c(0:i)
      end do
 
      return
      end
      subroutine test0175 ( )

c*********************************************************************72
c
cc TEST0175 tests CHARLIER.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )
      integer test_num
      parameter ( test_num = 5 )

      double precision a
      double precision a_test(test_num)
      integer i
      integer j
      integer test
      double precision x
      double precision value(0:n)

      save a_test

      data a_test / 0.25D+00, 0.5D+00, 1.0D+00, 2.0D+00, 10.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0175:'
      write ( *, '(a)' ) '  CHARLIER evaluates Charlier polynomials.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N      A         X        P(N,A,X)'
      write ( *, '(a)' ) ' '

      do test = 1, test_num

        a = a_test(test)

        write ( *, '(a)' ) ' '

        do j = 0, 5

          x = dble ( j ) / 2.0D+00

          call charlier ( n, a, x, value )

          write ( *, '(a)' ) ' '

          do i = 0, n

            write ( *, '(2x,i8,2x,f8.4,2x,f8.4,2x,g14.6)' ) 
     &        i, a, x, value(i)

          end do

        end do

      end do

      return
      end
      subroutine test018 ( )

c*********************************************************************72
c
cc TEST018 tests CHEBY_T_POLYNOMIAL and CHEBY_T_POLYNOMIAL_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 12 )

      double precision fx
      double precision fx2(0:n_max)
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST018:'
      write ( *, '(a)' ) 
     &  '  CHEBY_T_POLYNOMIAL evaluates the Chebyshev T polynomial.'
      write ( *, '(a)' ) 
     &  '  CHEBY_T_POLYNOMIAL_VALUES returns some exact values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     N      X        Exact F       T(N)(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call cheby_t_poly_values ( n_data, n, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        call cheby_t_poly ( 1, n, x, fx2 )

        write ( *, '(2x,i8,f8.4,2g14.6)' ) n, x, fx, fx2(n)

      go to 10

20    continue

      return
      end
      subroutine test0185 ( )

c*********************************************************************72
c
cc TEST0185 tests CHEBY_T_POLYNOMIAL_ZERO.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 4 )

      double precision fx(0:n_max)
      integer i
      integer n
      double precision z(n_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0185:'
      write ( *, '(a)' ) 
     &  '  CHEBY_T_POLYNOMIAL_ZERO returns zeroes of the T(N)(X).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N      X        T(N)(X)'
      write ( *, '(a)' ) ' '

      do n = 1, n_max

        call cheby_t_poly_zero ( n, z )

        do i = 1, n

          call cheby_t_poly ( 1, n, z(i), fx )

          write ( *, '(2x,i8,2x,f8.4,2x,g14.6)' ) n, z(i), fx(n)

        end do

        write ( *, '(a)' ) ' '

      end do

      return
      end
      subroutine test019 ( )

c*********************************************************************72
c
cc TEST019 tests CHEBY_T_POLYNOMIAL_COEF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      double precision c(0:n,0:n)
      integer i
      integer  j

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST019'
      write ( *, '(a)' ) '  CHEBY_T_POLYNOMIAL_COEF determines ' // 
     &  'the Chebyshev T polynomial coefficients.'

      call cheby_t_poly_coef ( n, c )
 
      do i = 0, n
        write ( *, '(a)' ) ' '
        write ( *, '(a,i2,a)' ) '  T(', i, ')'
        write ( *, '(a)' ) ' '
        do j = i, 0, -1
          if ( j .eq. 0 ) then
            write ( *, '(2x,g14.6)' ) c(i,j)
          else if ( j .eq. 1 ) then
            write ( *, '(2x,g14.6,a)' ) c(i,j), ' * x'
          else
            write ( *, '(2x,g14.6,a,i2)' ) c(i,j), ' * x**', j
          end if
        end do
      end do
 
      return
      end
      subroutine test020 ( )

c*********************************************************************72
c
cc TEST020 tests CHEBY_U_POLYNOMIAL and CHEBY_U_POLYNOMIAL_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 12 )

      double precision fx
      double precision fx2(0:n_max)
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST020:'
      write ( *, '(a)' ) 
     &  '  CHEBY_U_POLYNOMIAL evaluates the Chebyshev U polynomial.'
      write ( *, '(a)' ) 
     &  '  CHEBY_U_POLYNOMIAL_VALUES returns some exact values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     N      X        Exact F       U(N)(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call cheby_u_poly_values ( n_data, n, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        call cheby_u_poly ( n, x, fx2 )

        write ( *, '(2x,i8,f8.4,2g14.6)' ) n, x, fx, fx2(n)

      go to 10

20    continue

      return
      end
      subroutine test021 ( )

c*********************************************************************72
c
cc TEST021 tests CHEBY_U_POLYNOMIAL_COEF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      double precision c(0:n,0:n)
      integer i
      integer  j

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST021'
      write ( *, '(a)' ) '  CHEBY_U_POLYNOMIAL_COEF determines ' // 
     &  'the Chebyshev U polynomial coefficients.'

      call cheby_u_poly_coef ( n, c )
 
      do i = 0, n
        write ( *, '(a)' ) ' '
        write ( *, '(a,i2,a)' ) '  T(', i, ')'
        write ( *, '(a)' ) ' '
        do j = i, 0, -1
          if ( j .eq. 0 ) then
            write ( *, '(2x,g14.6)' ) c(i,j)
          else if ( j .eq. 1 ) then
            write ( *, '(2x,g14.6,a)' ) c(i,j), ' * x'
          else
            write ( *, '(2x,g14.6,a,i2)' ) c(i,j), ' * x**', j
          end if
        end do
      end do
 
      return
      end
      subroutine test0215 ( )

c*********************************************************************72
c
cc TEST0215 tests CHEBY_U_POLYNOMIAL_ZERO.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 4 )

      double precision fx(0:n_max)
      integer i
      integer n
      double precision z(n_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0215:'
      write ( *, '(a)' ) 
     &  '  CHEBY_U_POLYNOMIAL_ZERO returns zeroes of the U(N)(X).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N      X        U(N)(X)'
      write ( *, '(a)' ) ' '

      do n = 1, n_max

        call cheby_u_poly_zero ( n, z )

        do i = 1, n

          call cheby_u_poly ( n, z(i), fx )

          write ( *, '(2x,i8,2x,f8.4,2x,g14.6)' ) n, z(i), fx(n)

        end do

        write ( *, '(a)' ) ' '

      end do

      return
      end
      subroutine test0216 ( )

c*********************************************************************72
c
cc TEST0216 tests CHEBYSHEV_DISCRETE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      integer i
      integer j
      integer m
      double precision x
      double precision value(0:n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0216:'
      write ( *, '(a)' ) 
     &  '  CHEBYSHEV_DISCRETE evaluates discrete Chebyshev polynomials.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N      M         X        T(N,M,X)'
      write ( *, '(a)' ) ' '

      m = 5

      do j = 0, 5

        x = dble ( j ) / 2.0D+00

        call chebyshev_discrete ( n, m, x, value )

        write ( *, '(a)' ) ' '

        do i = 0, n

          write ( *, '(2x,i8,2x,i8,2x,f8.4,2x,g14.6)' ) 
     &      i, m, x, value(i)

        end do

      end do

      return
      end
      subroutine test0217 ( )

c*********************************************************************72
c
cc TEST0217 tests COLLATZ_COUNT and COLLATZ_COUNT_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer collatz_count
      integer count
      integer count2
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0217:'
      write ( *, '(a)' ) '  COLLATZ_COUNT(N) counts the length of the'
      write ( *, '(a)' ) '  Collatz sequence beginning with N.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N       COUNT(N)     COUNT(N)'
      write ( *, '(a)' ) '              (computed)    (table)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call collatz_count_values ( n_data, n, count )

        if ( n_data .eq. 0 ) then
          go to 20
        end if
 
        count2 = collatz_count ( n )

        write ( *, '(2x,i8,2x,i8,2x,i8)' ) n, count, count2

      go to 10

20    continue

      return
      end
      subroutine test0218 ( )

c*********************************************************************72
c
cc TEST0218 tests COLLATZ_COUNT_MAX.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 April 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i_max
      integer j_max
      integer n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0218:'
      write ( *, '(a)' ) '  COLLATZ_COUNT_MAX(N) returns the length of'
      write ( *, '(a)' ) '  the longest Collatz sequence from 1 to N.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         N     I_MAX     J_MAX'
      write ( *, '(a)' ) ' '

      n = 10

10    continue

      if ( n <= 100000 ) then

        call collatz_count_max ( n, i_max, j_max )

        write ( *, '(2x,i8,2x,i8,2x,i8)' ) n, i_max, j_max

        n = n * 10

        go to 10

      end if

      return
      end
      subroutine test024 ( )

c*********************************************************************72
c
cc TEST024 tests COMB_ROW.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer c(0:n)
      integer i
      integer ido

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST024'
      write ( *, '(a)' ) 
     &  '  COMB_ROW computes a row of Pascal''s triangle.'
      write ( *, '(a)' ) ' '
 
      ido = 0
 
      do i = 0, n
        call comb_row ( ido, i, c )
        ido = 1
        write ( *, '(2x,i2,2x,11i5)' ) i, c(0:i)
      end do
 
      return
      end
      subroutine test02405 ( )

c*********************************************************************72
c
cc TEST02405 tests COMMUL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer factor(4)
      integer i
      integer ncomb
      integer nfactor

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02405'
      write ( *, '(a)' ) '  COMMUL computes a multinomial coefficient.'
      write ( *, '(a)' ) ' '

      n = 8
      nfactor = 2
      factor(1) = 6
      factor(2) = 2

      call commul ( n, nfactor, factor, ncomb ) 

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  N = ', n
      write ( *, '(a,i8)' ) '  Number of factors = ', nfactor
      do i = 1, nfactor
        write ( *, '(2x,i2,2x,i8)' ) i, factor(i)
      end do
      write ( *, '(a,i12)' ) '  Value of coefficient = ', ncomb

      n = 8
      nfactor = 3
      factor(1) = 2
      factor(2) = 2
      factor(3) = 4
      call commul ( n, nfactor, factor, ncomb )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  N = ', n
      write ( *, '(a,i8)' ) '  Number of factors = ', nfactor
      do i = 1, nfactor
        write ( *, '(2x,i2,2x,i8)' ) i, factor(i)
      end do
      write ( *, '(a,i12)' ) '  Value of coefficient = ', ncomb

      n = 13
      nfactor = 4
      factor(1) = 5
      factor(2) = 3
      factor(3) = 3
      factor(4) = 2
      call commul ( n, nfactor, factor, ncomb )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  N = ', n
      write ( *, '(a,i8)' ) '  Number of factors = ', nfactor
      do i = 1, nfactor
        write ( *, '(2x,i2,2x,i8)' ) i, factor(i)
      end do
      write ( *, '(a,i12)' ) '  Value of coefficient = ', ncomb

      return
      end
      subroutine test0241 ( )

c*********************************************************************72
c
cc TEST0241 tests COS_DEG.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision angle
      double precision cos_deg
      integer i

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0241'
      write ( *, '(a)' ) '  COS_DEG computes the cosine of an angle'
      write ( *, '(a)' ) '  given in degrees.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  ANGLE    COS_DEG(ANGLE)'
      write ( *, '(a)' ) ' '
 
      do i = 0, 360, 10
        angle = dble ( i )
        write ( *, '(2x,f8.2,2x,g14.6)' )  angle, cos_deg ( angle )
      end do
 
      return
      end
      subroutine test0242 ( )

c*********************************************************************72
c
cc TEST0242 tests R8_CAS.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision angle
      double precision angle_rad
      double precision r8_cas
      integer i
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0242'
      write ( *, '(a)' ) '  R8_CAS computes the "casine" of an angle.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  ANGLE    R8_CAS(ANGLE)'
      write ( *, '(a)' ) ' '
 
      do i = 0, 360, 10
        angle = dble ( i )
        angle_rad = angle * pi / 180.0D+00
        write ( *, '(2x,f8.2,2x,g14.6)' )  angle, r8_cas ( angle_rad )
      end do
 
      return
      end
      subroutine test0243 ( )

c*********************************************************************72
c
cc TEST0243 tests COS_POWER_INT and COS_POWER_INT_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    31 March 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision b
      double precision cos_power_int
      double precision fx
      double precision fx2
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0243:'
      write ( *, '(a)' ) '  COS_POWER_INT returns values of '
      write ( *, '(a)' ) '  the integral of COS(X)^N from A to B.'
      write ( *, '(a)' ) 
     &  '  COS_POWER_INT_VALUES stores some selected values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '      A         B          N      Exact           Computed'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call cos_power_int_values ( n_data, a, b, n, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = cos_power_int ( a, b, n )

        write ( *, '(2x,f8.4,2x,f8.4,2x,i8,2x,g14.6,2x,g14.6)' ) 
     &    a, b, n, fx, fx2

      go to 10

20    continue

      return
      end
      subroutine test01155 ( )

c*********************************************************************72
c
cc TEST01155 tests DELANNOY.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 8 )
      integer n
      parameter ( n = 8 )

      integer a(0:m,0:n)
      integer i

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01155'
      write ( *, '(a)' ) 
     &  '  DELANNOY computes the Delannoy numbers A(0:M,0:N).'
      write ( *, '(a)' ) 
     &  '  A(M,N) counts the paths from (0,0) to (M,N).'
      write ( *, '(a)' ) ' '

      call delannoy ( m, n, a )

      do i = 0, m
        write ( *, '(2x,i4,2x,5i4,3i8,i10)' )  i, a(i,0:n)
      end do
  
      return
      end
      subroutine test0245 ( )

c*********************************************************************72
c
cc TEST0245 tests R8_FACTORIAL and R8_FACTORIAL_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fn
      double precision fn2
      integer n_data
      integer n
      double precision r8_factorial

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0245:'
      write ( *, '(a)' ) 
     &  '  R8_FACTORIAL evaluates the factorial function.'
      write ( *, '(a)' ) 
     &  '  R8_FACTORIAL_VALUES returns some exact values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     N       Exact F       R8_FACTORIAL(N)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call r8_factorial_values ( n_data, n, fn )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fn2 = r8_factorial ( n )

        write ( *, '(2x,i4,2g14.6)' ) n, fn, fn2

      go to 10

20    continue

      return
      end
      subroutine test025 ( )

c*********************************************************************72
c
cc TEST025 tests ERROR_F and ERF_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision error_f
      double precision fx
      double precision fx2
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST025:'
      write ( *, '(a)' ) '  ERROR_F evaluates the error function.'
      write ( *, '(a)' ) '  ERF_VALUES returns some exact values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     X      Exact F       ERROR_F(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call erf_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = error_f ( x )

        write ( *, '(2x,f8.4,2g14.6)' ) x, fx, fx2

      go to 10

20    continue

      return
      end
      subroutine test0255 ( )

c*********************************************************************72
c
cc TEST0255 tests ERROR_F and ERF_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 August 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision error_f_inverse
      double precision fx
      integer n_data
      double precision x1
      double precision x2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0255:'
      write ( *, '(a)' ) '  ERROR_F_INVERSE inverts the error function.'
      write ( *, '(a)' ) '  ERF_VALUES returns some exact values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    FX            X    ERROR_F_INVERSE(FX)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call erf_values ( n_data, x1, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        x2 = error_f_inverse ( fx )

        write ( *, '(2x,f8.4,2g14.6)' ) fx, x1, x2

      go to 10

20    continue

      return
      end
      subroutine test026 ( )

c*********************************************************************72
c
cc TEST026 tests EULER_NUMBER and EULER_NUMBER_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer c1
      integer c2(0:12)
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST026'
      write ( *, '(a)' ) '  EULER_NUMBER computes Euler numbers.'
      write ( *, '(a)' ) 
     &  '  EULER_NUMBER_VALUES returns some exact values.'
      write ( *, '(a)' ) ' ' 
      write ( *, '(a)' ) '     N       exact   EULER_NUMBER'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call euler_number_values ( n_data, n, c1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        call euler_number ( n, c2 )

        write ( *, '(2x,i4,2i12,g14.6)' ) n, c1, c2(n)

      go to 10

20    continue
 
      return
      end
      subroutine test0265 ( )

c*********************************************************************72
c
cc TEST0265 tests EULER_NUMBER2 and EULER_NUMBER_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer c1
      double precision c2
      double precision euler_number2
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0265'
      write ( *, '(a)' ) 
     &  '  EULER_NUMBER2 computes Euler numbers.'
      write ( *, '(a)' ) 
     &  '  EULER_NUMBER_VALUES returns some exact values.'
      write ( *, '(a)' ) ' ' 
      write ( *, '(a)' ) '     N       exact   EULER_NUMBER2'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call euler_number_values ( n_data, n, c1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        c2 = euler_number2 ( n )

        write ( *, '(2x,i4,i12,g14.6)' ) n, c1, c2

      go to 10

20    continue
 
      return
      end
      subroutine test028 ( )

c*********************************************************************72
c
cc TEST028 tests EULER_POLY.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision euler_poly
      double precision f
      integer i
      integer n
      parameter ( n = 15 )
      double precision x

      x = 0.5D+00
 
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST028'
      write ( *, '(a)' ) '  EULER_POLY evaluates Euler polynomials.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   N      X             F(X)'
      write ( *, '(a)' ) ' '
   
      do i = 0, n
        f = euler_poly ( i, x )
        write ( *, '(2x,i2,2x,2g14.6)' ) i, x, f
      end do
 
      return
      end
      subroutine test027 ( )

c*********************************************************************72
c
cc  TEST027 tests EULERIAN.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 7 )

      integer e(n,n)
      integer i

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST027'
      write ( *, '(a)' ) '  EULERIAN evaluates Eulerian numbers.'
      write ( *, '(a)' ) ' '
 
      call eulerian ( n, e )

      do i = 1, n
        write ( *, '(2x,10i6)' )  e(i,1:n)
      end do
 
      return
      end
      subroutine test031 ( )

c*********************************************************************72
c
cc TEST031 tests FIBONACCI_DIRECT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer f
      integer i
      integer n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST031'
      write ( *, '(a)' ) 
     &  '  FIBONACCI_DIRECT evalutes a Fibonacci number directly.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       I        F(I)'
      write ( *, '(a)' ) ' '

      n = 20
     
      do i = 1, n
        call fibonacci_direct ( i, f )
        write ( *, '(2x,i8,i10)' ) i, f
      end do
     
      return
      end
      subroutine test032 ( )

c*********************************************************************72
c
cc TEST032 tests FIBONACCI_FLOOR.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer f
      integer i
      integer n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST032'
      write ( *, '(a)' ) 
     &  '  FIBONACCI_FLOOR computes the largest Fibonacci number'
      write ( *, '(a)' ) 
     &  '  less than or equal to a given positive integer.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N  Fibonacci  Index'
      write ( *, '(a)' ) ' ' 

      do n = 1, 20
        call fibonacci_floor ( n, f, i )
        write ( *, '(2x,i8,2x,i8,2x,i8)' ) n, f, i
      end do
     
      return
      end
      subroutine test036 ( )

c*********************************************************************72
c
cc TEST036 tests R8_GAMMA_LOG and GAMMA_LOG_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      double precision fx2
      integer n_data
      double precision r8_gamma_log
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST036:'
      write ( *, '(a)' ) '  R8_GAMMA_LOG evaluates the logarithm of the'
      write ( *, '(a)' ) '  Gamma function.'
      write ( *, '(a)' ) 
     &  '  GAMMA_LOG_VALUES returns some exact values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     X       Exact F       GAMMA_LOG(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call gamma_log_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r8_gamma_log ( x )

        write ( *, '(2x,f8.4,2g18.10)' ) x, fx, fx2

      go to 10

20    continue

      return
      end
      subroutine test037 ( )

c*********************************************************************72
c
cc TEST037 tests GEGENBAUER_POLY and GEGENBAUER_POLY_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 10 )

      double precision a
      double precision c(0:n_max)
      double precision fx
      double precision fx2
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST037:'
      write ( *, '(a)' ) '  GEGENBAUER_POLY computes values of '
      write ( *, '(a)' ) '  the Gegenbauer polynomials.'
      write ( *, '(a)' ) '  GEGENBAUER_POLY_VALUES returns values of '
      write ( *, '(a)' ) '  the Gegenbauer polynomials.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '       N        A           X       GPV      GEGENBAUER'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call gegenbauer_poly_values ( n_data, n, a, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        call gegenbauer_poly ( n, a, x, c )
        fx2 = c(n)

        write ( *, '(2x,i8,2x,f10.4,2x,f10.4,2g14.6)' ) 
     &    n, a, x, fx, fx2

      go to 10

20    continue

      return
      end
      subroutine test038 ( )

c*********************************************************************72
c
cc TEST038 tests GUD and GUD_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      double precision fx2
      double precision gud
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST038:'
      write ( *, '(a)' ) '  GUD evaluates the Gudermannian function.'
      write ( *, '(a)' ) '  GUD_VALUES returns some exact values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     X      Exact F       GUD(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call gud_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = gud ( x )

        write ( *, '(2x,f8.4,2g14.6)' ) x, fx, fx2

      go to 10

20    continue

      return
      end
      subroutine test041 ( )

c*********************************************************************72
c
cc TEST041 tests HERMITE_POLY and HERMITE_POLY_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 12 )

      double precision fx
      double precision fx2(0:n_max)
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST041:'
      write ( *, '(a)' ) 
     &  '  HERMITE_POLY evaluates the Hermite polynomial.'
      write ( *, '(a)' ) 
     &  '  HERMITE_POLY_VALUES returns some exact values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         N    X      Exact F       H(N)(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call hermite_poly_values ( n_data, n, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        call hermite_poly ( n, x, fx2 )

        write ( *, '(2x,i8,f8.4,2g14.6)' ) n, x, fx, fx2(n)

      go to 10

20    continue

      return
      end
      subroutine test042 ( )

c*********************************************************************72
c
cc TEST042 tests HERMITE_POLY_COEF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      double precision c(0:n,0:n)
      integer i
      integer j

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST042'
      write ( *, '(a)' ) 
     &  '  HERMITE_POLY_COEF determines' // 
     &  ' the Hermite polynomial coefficients.'

      call hermite_poly_coef ( n, c )
     
      do i = 0, n
        write ( *, '(a)' ) ' '
        write ( *, '(a,i2,a)' ) '  H(', i, ')'
        write ( *, '(a)' ) ' '
        do j = i, 0, -1
          if ( j .eq. 0 ) then
            write ( *, '(2x,g14.6)' ) c(i,j)
          else if ( j .eq. 1 ) then
            write ( *, '(2x,g14.6,a)' ) c(i,j), ' * x'
          else
            write ( *, '(2x,g14.6,a,i2)' ) c(i,j), ' * x**', j
          end if
        end do
      end do
     
      return
      end
      subroutine test0425 ( )

c*********************************************************************72
c
cc TEST0425 tests R8_HYPER_2F1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision b
      double precision c
      double precision fx
      double precision fx2
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0425:'
      write ( *, '(a)' ) 
     &  '  R8_HYPER_2F1 evaluates the hypergeometric 2F1 function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) '      A       B       C       X      ', 
     &  ' 2F1                       2F1                     DIFF'
      write ( *, '(a,a)' ) '                                     ', 
     &  '(tabulated)               (computed)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call hyper_2f1_values ( n_data, a, b, c, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        call r8_hyper_2f1 ( a, b, c, x, fx2 )

        write ( *, 
     &  '(2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,' //
     &  'g24.16,2x,g24.16,2x,g10.4)' ) 
     &  a, b, c, x, fx, fx2, abs ( fx - fx2 )

      go to 10

20    continue

      return
      end
      subroutine test023 ( )

c*********************************************************************72
c
cc TEST023 tests I4_CHOOSE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer cnk
      integer i4_choose
      integer k
      integer n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST023'
      write ( *, '(a)' ) '  I4_CHOOSE evaluates C(N,K).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     N     K    CNK'
      write ( *, '(a)' ) ' '
 
      do n = 0, 4
        do k = 0, n
          cnk = i4_choose ( n, k )
          write ( *, '(2x,i8,2x,i8,2x,i8)' ) n, k, cnk
        end do
      end do
 
      return
      end
      subroutine test043 ( )

c*********************************************************************72
c
cc TEST043 tests I4_FACTORIAL and I4_FACTORIAL_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer fn
      integer fn2
      integer i4_factorial
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST043:'
      write ( *, '(a)' ) 
     &  '  I4_FACTORIAL evaluates the factorial function.'
      write ( *, '(a)' ) 
     &  '  I4_FACTORIAL_VALUES returns some exact values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     X       Exact F       I4_FACTORIAL(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call i4_factorial_values ( n_data, n, fn )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fn2 = i4_factorial ( n )

        write ( *, '(2x,i4,2i12)' ) n, fn, fn2

      go to 10

20    continue

      return
      end
      subroutine test044 ( )

c*********************************************************************72
c
cc TEST044 tests I4_FACTORIAL2 and I4_FACTORIAL2_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer fn
      integer fn2
      integer n
      integer n_data
      integer i4_factorial2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST044:'
      write ( *, '(a)' ) 
     &  '  I4_FACTORIAL2 evaluates the double factorial function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   N   Exact  I4_FACTORIAL2(N)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call i4_factorial2_values ( n_data, n, fn )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fn2 = i4_factorial2 ( n )

        write ( *, '(2x,i4,2i8)' ) n, fn, fn2

      go to 10

20    continue

      return
      end
      subroutine test045 ( )

c*********************************************************************72
c
cc TEST045 tests PARTITION_COUNT_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer c
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST045:'
      write ( *, '(a)' ) '  For the number of partitions of an integer,'
      write ( *, '(a)' ) 
     &  '  PARTITION_COUNT_VALUES returns some exact values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '           N       Exact F'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call partition_count_values ( n_data, n, c )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        write ( *, '(2x,2i10)' ) n, c

      go to 10

20    continue

      return
      end
      subroutine test046 ( )

c*********************************************************************72
c
cc TEST046 tests I4_PARTITION_DISTINCT_COUNT and PARTITION_DISTINCT_COUNT_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer c
      integer c2
      integer n
      integer n_data
      integer n_max
      parameter ( n_max = 20 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST046:'
      write ( *, '(a)' ) '  For the number of partitions of an integer'
      write ( *, '(a)' ) '  into distinct parts,'
      write ( *, '(a)' ) '  I4_PARTITION_DISTINCT_COUNT'
      write ( *, '(a)' ) '  computes any value.'
      write ( *, '(a)' ) '  PARTITION_DISTINCT_COUNT_VALUES '
      write ( *, '(a)' ) '  returns some exact values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '           N       Exact F    Q(N)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call partition_distinct_count_values ( n_data, n, c )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        if ( n_max .lt. n ) then
          go to 10
        end if

        call i4_partition_distinct_count ( n, c2 )

        write ( *, '(2x,3i10)' ) n, c, c2

      go to 10

20    continue

      return
      end
      subroutine test047 ( )

c*********************************************************************72
c
cc TEST047 tests I4_POCHHAMMER.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      integer i4_pochhammer
      integer j
      integer k

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST047:'
      write ( *, '(a)' ) 
     &  '  I4_POCHHAMMER evaluates the Pochhammer function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I   J   I4_Pochhammer(I,J)'
      write ( *, '(a)' ) ' '

      i = 3
      j = 3
      k = i4_pochhammer ( i, j, k )

      write ( *, '(2x,i4,i4,i4)' ) i, j,  k

      i = 3
      j = 4
      k = i4_pochhammer ( i, j, k )
      write ( *, '(2x,i4,i4,i4)' ) i, j,  k

      i = 3
      j = 5
      k = i4_pochhammer ( i, j, k )
      write ( *, '(2x,i4,i4,i4)' ) i, j,  k

      return
      end
      subroutine test048 ( )

c*********************************************************************72
c
cc TEST048 tests I4_IS_TRIANGULAR, I4_TO_TRIANGLE and TRIANGLE_TO_I4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      logical i4_is_triangular
      integer i2
      integer j
      integer k
      logical l

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST048'
      write ( *, '(a)' ) '  I4_TO_TRIANGLE converts a linear index to a'
      write ( *, '(a)' ) '  triangular one.'
      write ( *, '(a)' ) '  TRIANGLE_TO_I4 converts a triangular index'
      write ( *, '(a)' ) '  to a linear one.'
      write ( *, '(a)' ) '  I4_IS_TRIANGULAR returns T or F depending'
      write ( *, '(a)' ) '  on whether I is triangular.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I  =>   J   K  =>   I   T/F'
      write ( *, '(a)' ) ' '

      do i = 0, 20

        call i4_to_triangle ( i, j, k )

        call triangle_to_i4 ( j, k, i2 )

        l = i4_is_triangular ( i )

        write ( *, '(2x,i4,4x,i4,i4,4x,i4,4x,l1)' )  i, j, k, i2, l

      end do
     
      return
      end
      subroutine test049 ( )

c*********************************************************************72
c
cc TEST049 tests JACOBI_POLY and JACOBI_POLY_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 April 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision b
      double precision c(0:6)
      double precision fx
      double precision fx2
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST049:'
      write ( *, '(a)' ) '  JACOBI_POLY computes values of '
      write ( *, '(a)' ) '  the Jacobi polynomial.'
      write ( *, '(a)' ) '  JACOBI_POLY_VALUES returns values of '
      write ( *, '(a)' ) '  the Jacobi polynomial.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '       N       A       B      X       JPV      JACOBI'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call jacobi_poly_values ( n_data, n, a, b, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        call jacobi_poly ( n, a, b, x, c )
        fx2 = c(n)

        write ( *, '(2x,i8,2x,f8.4,2x,f8.4,f10.4,2g14.6)' ) 
     &  n, a, b, x, fx, fx2

      go to 10

20    continue

      return
      end
      subroutine test050 ( )

c*********************************************************************72
c
cc TEST050 tests JACOBI_SYMBOL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 4 )

      integer l
      integer p
      integer p_test(test_num)
      integer q
      integer test

      save p_test

      data p_test / 3, 9, 10, 12 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST050'
      write ( *, '(a)' ) '  JACOBI_SYMBOL computes the Jacobi symbol'
      write ( *, '(a)' ) '  (Q/P), which records if Q is a quadratic '
      write ( *, '(a)' ) '  residue modulo the number P.'

      do test = 1, test_num
        p = p_test(test)
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Jacobi Symbols for P = ', p
        write ( *, '(a)' ) ' '
        do q = 0, p
          call jacobi_symbol ( q, p, l )
          write ( *, '(2x,3i8)' ) p, q, l
        end do
      end do

      return
      end
      subroutine test0505 ( )

c*********************************************************************72
c
cc TEST0505 tests KRAWTCHOUK
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )
      integer test_num
      parameter ( test_num = 2 )

      integer i
      integer j
      integer m
      double precision p
      double precision p_test(test_num)
      integer test
      double precision x
      double precision value(0:n)

      save p_test

      data p_test / 0.25D+00, 0.50D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0505:'
      write ( *, '(a)' ) 
     &  '  KRAWTCHOUK evaluates Krawtchouk polynomials.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '         N      P         X            M    K(N,P,X,M)'
      write ( *, '(a)' ) ' '

      m = 5

      do test = 1, test_num

        p = p_test(test)

        do j = 0, 5

          x = dble ( j ) / 2.0D+00

          call krawtchouk ( n, p, x, m, value )

          write ( *, '(a)' ) ' '

          do i = 0, n

            write ( *, '(2x,i8,2x,f8.4,2x,f8.4,2x,i8,2x,g14.6)' ) 
     &        i, p, x, m, value(i)

          end do

        end do

      end do

      return
      end
      subroutine test051 ( )

c*********************************************************************72
c
cc TEST051 tests LAGUERRE_ASSOCIATED.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 6 )
      integer n
      parameter ( n = 6 )

      double precision c(0:n)
      integer j
      integer m
      integer m_test(test_num)
      integer test
      double precision x
      double precision x_test(test_num)

      save m_test
      save x_test

      data m_test / 0, 0, 1, 2, 3, 1 /
      data x_test /
     &  0.0D+00, 1.0D+00, 0.0D+00, 0.5D+00, 0.5D+00, 0.5D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST051'
      write ( *, '(a)' ) 
     &  '  LAGUERRE_ASSOCIATED evaluates the associated Laguerre'
      write ( *, '(a)' ) '  polynomials.'

      do test = 1, test_num

        m = m_test(test)
        x = x_test(test)

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Table of L(N,M)(X) for'
        write ( *, '(a)' ) ' '
        write ( *, '(a,i4)' ) '  N(max) = ', n
        write ( *, '(a,i4)' ) '  M      = ', m
        write ( *, '(a,g14.6)' ) '  X =      ', x
        write ( *, '(a)' ) ' '
     
        call laguerre_associated ( n, m, x, c )
     
        do j = 0, n
          write ( *, '(2x,i8,g14.6)' ) j, c(j)
        end do
     
      end do

      return
      end
      subroutine test052 ( )

c*********************************************************************72
c
cc TEST052 tests GEN_LAGUERRE_POLY.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 6 )
      integer n
      parameter ( n = 10 )

      double precision alpha
      double precision alpha_test(test_num)
      double precision c(0:n)
      integer j
      integer test
      double precision x
      double precision x_test(test_num)

      save alpha_test
      save x_test

      data alpha_test /
     &  0.0D+00, 0.0D+00, 0.1D+00, 0.1D+00, 0.5D+00, 1.0D+00 /
      data x_test /
     &  0.0D+00, 1.0D+00, 0.0D+00, 0.5D+00, 0.5D+00, 0.5D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST052'
      write ( *, '(a)' ) '  GEN_LAGUERRE_POLY evaluates the generalized'
      write ( *, '(a)' ) '  Laguerre functions.'

      do test = 1, test_num

        x = x_test(test)
        alpha = alpha_test(test)

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Table of L(N,ALPHA)(X) for'
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '    N(max) = ', n
        write ( *, '(a,g14.6)' ) '    ALPHA =  ', alpha
        write ( *, '(a,g14.6)' ) '    X =      ', x
        write ( *, '(a)' ) ' '
      
        call gen_laguerre_poly ( n, alpha, x, c )
     
        do j = 0, n
          write ( *, '(2x,i8,g14.6)' ) j, c(j)
        end do

      end do
     
      return
      end
      subroutine test054 ( )

c*********************************************************************72
c
cc TEST054 tests LAGUERRE_POLY and LAGUERRE_POLYNOMIAL_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 12 )

      double precision fx
      double precision fx2(0:n_max)
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST054:'
      write ( *, '(a)' ) 
     &  '  LAGUERRE_POLY evaluates the Laguerre polynomial.'
      write ( *, '(a)' ) 
     &  '  LAGUERRE_POLYNOMIAL_VALUES returns some exact values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         N    X      Exact F       L(N)(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

      do

        call laguerre_polynomial_values ( n_data, n, x, fx )

        if ( n_data == 0 ) then
          exit
        end if

        call laguerre_poly ( n, x, fx2 )

        write ( *, '(2x,i8,f8.4,2g14.6)' ) n, x, fx, fx2(n)

      end do

      return
      end
      subroutine test055 ( )

c*********************************************************************72
c
cc TEST055 tests LAGUERRE_POLY_COEF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      double precision c(0:n,0:n)
      double precision fact
      integer i
      integer j

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST055'
      write ( *, '(a)' ) '  LAGUERRE_POLY_COEF determines the ' //
     &  'Laguerre polynomial coefficients.'

      call laguerre_poly_coef ( n, c )
     
      do i = 0, n
        write ( *, '(a)' ) ' '
        write ( *, '(a,i2,a)' ) '  L(', i, ')'
        write ( *, '(a)' ) ' '
        do j = i, 0, -1
          if ( j .eq. 0 ) then
            write ( *, '(2x,g14.6)' ) c(i,j)
          else if ( j .eq. 1 ) then
            write ( *, '(2x,g14.6,a)' ) c(i,j), ' * x'
          else
            write ( *, '(2x,g14.6,a,i2)' ) c(i,j), ' * x**', j
          end if
        end do
      end do
     
      fact = 1.0D+00

      do i = 0, n

        if ( 0 .lt. i ) then
          fact = fact * dble ( i )
        end if

        write ( *, '(a)' ) ' '
        write ( *, '(a,i2,a)' ) '  Factorially scaled L(', i, ')'
        write ( *, '(a)' ) ' '

        do j = i, 0, -1
          if ( j == 0 ) then
            write ( *, '(2x,g14.6)' ) fact * c(i,j)
          else if ( j == 1 ) then
            write ( *, '(2x,g14.6,a)' ) fact * c(i,j), ' * x'
          else
            write ( *, '(2x,g14.6,a,i2)' ) fact * c(i,j), ' * x**', j
          end if
        end do
        
      end do

      return
      end
      subroutine test0552 ( )

c*********************************************************************72
c
cc TEST0552 tests LAMBERT_W, LAMBERT_W_CRUDE and LAMBERT_W_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      double precision fx2
      double precision fx3
      double precision lambert_w
      double precision lambert_w_crude
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0552:'
      write ( *, '(a)' ) 
     &  '  LAMBERT_W estimates the Lambert W function.'
      write ( *, '(a)' ) 
     &  '  LAMBERT_W_CRUDE makes a crude estimate of the'
      write ( *, '(a)' ) '    Lambert W function.'
      write ( *, '(a)' ) 
     &  '  LAMBERT_W_VALUES returns some tabulated values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '           X           W(X)        W(X)         W(X)'
      write ( *, '(a)' ) 
     &  '                   Tabulated       Crude     Estimate'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call lambert_w_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = lambert_w_crude ( x )

        fx3 = lambert_w ( x )

        write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &    x, fx, fx2, fx3

      go to 10

20    continue

      return
      end
      subroutine test059 ( )

c*********************************************************************72
c
cc TEST059 tests LEGENDRE_ASSOCIATED and LEGENDRE_ASSOCIATED_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 20 )

      double precision fx2(0:n_max)
      double precision fx
      integer m
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST059:'
      write ( *, '(a)' ) 
     &  '  LEGENDRE_ASSOCIATED evaluates associated Legendre functions.'
      write ( *, '(a)' ) 
     &  '  LEGENDRE_ASSOCIATED_VALUES returns some exact values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '       N       M    X     Exact F       PNM(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call legendre_associated_values ( n_data, n, m, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        call legendre_associated ( n, m, x, fx2 )

        write ( *, '(2x,i8,2x,i8,f8.4,2g14.6)' ) n, m, x, fx, fx2(n)

      go to 10

20    continue

      return
      end
      subroutine test0595 ( )

c*********************************************************************72
c
cc TEST0595 tests LEGENDRE_ASSOCIATED_NORMALIZED and LEGENDRE_ASSOCIATED_NORMALIZED_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 September 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 20 )

      double precision fx2(0:n_max)
      double precision fx
      integer m
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0595:'
      write ( *, '(a)' ) 
     &  '  LEGENDRE_ASSOCIATED_NORMALIZED evaluates normalized ' //
     &  'associated Legendre functions.'
      write ( *, '(a)' ) 
     &  '  LEGENDRE_ASSOCIATED_NORMALIZED_VALUES returns some ' //
     &  'exact values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '       N       M    X     Exact F       PNM(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call legendre_associated_normalized_values ( n_data, n, m, x, 
     &    fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        call legendre_associated_normalized ( n, m, x, fx2 )

        write ( *, '(2x,i8,2x,i8,f8.4,2g14.6)' ) n, m, x, fx, fx2(n)

      go to 10

20    continue

      return
      end
      subroutine test060 ( )

c*********************************************************************72
c
cc TEST060 tests LEGENDRE_FUNCTION_Q and LEGENDRE_FUNCTION_Q_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 12 )

      double precision fx
      double precision fx2(0:n_max)
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST060:'
      write ( *, '(a)' ) 
     &  '  LEGENDRE_FUNCTION_Q evaluates the Legendre Q function.'
      write ( *, '(a)' ) 
     &  '  LEGENDRE_FUNCTION_Q_VALUES returns some exact values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         N    X      Exact F       Q(N)(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call legendre_function_q_values ( n_data, n, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        call legendre_function_q ( n, x, fx2 )

        write ( *, '(2x,i8,f8.4,2g14.6)' ) n, x, fx, fx2(n)

      go to 10

20    continue

      return
      end
      subroutine test057 ( )

c*********************************************************************72
c
cc TEST057 tests LEGENDRE_POLY and LEGENDRE_POLY_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 12 )

      double precision fx
      double precision fp2(0:n_max)
      double precision fx2(0:n_max)
      integer n
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST057:'
      write ( *, '(a)' ) 
     &  '  LEGENDRE_POLY evaluates the Legendre PN function.'
      write ( *, '(a)' ) 
     &  '  LEGENDRE_POLY_VALUES returns some exact values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         N    X      Exact F       P(N)(X)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call legendre_poly_values ( n_data, n, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        call legendre_poly ( n, x, fx2, fp2 )

        write ( *, '(2x,i8,f8.4,2g14.6)' ) n, x, fx, fx2(n)

      go to 10

20    continue

      return
      end
      subroutine test058 ( )

c*********************************************************************72
c
cc TEST058 tests LEGENDRE_POLY_COEF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      double precision c(0:n,0:n)
      integer i
      integer j

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST058'
      write ( *, '(a)' ) 
     &  '  LEGENDRE_POLY_COEF returns Legendre polynomial coefficients.'

      call legendre_poly_coef ( n, c )
     
      do i = 0, n
        write ( *, '(a)' ) ' '
        write ( *, '(a,i2,a)' ) '  P(', i, ')'
        write ( *, '(a)' ) ' '
        do j = i, 0, -1
          if ( j .eq. 0 ) then
            write ( *, '(2x,g14.6)' ) c(i,j)
          else if ( j .eq. 1 ) then
            write ( *, '(2x,g14.6,a)' ) c(i,j), ' * x'
          else
            write ( *, '(2x,g14.6,a,i2)' ) c(i,j), ' * x**', j
          end if
        end do
      end do

      return
      end
      subroutine test061 ( )

c*********************************************************************72
c
cc TEST061 tests LEGENDRE_SYMBOL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 4 )

      integer l
      integer p
      integer p_test(test_num)
      integer q
      integer test

      save p_test

      data p_test / 7, 11, 13, 17 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST061'
      write ( *, '(a)' ) '  LEGENDRE_SYMBOL computes the Legendre'
      write ( *, '(a)' ) '  symbol (Q/P) which records whether Q is '
      write ( *, '(a)' ) '  a quadratic residue modulo the prime P.'

      do test = 1, test_num
        p = p_test(test)
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Legendre Symbols for P = ', p
        write ( *, '(a)' ) ' '
        do q = 0, p
          call legendre_symbol ( q, p, l )
          write ( *, '(2x,3i8)' ) p, q, l
        end do
      end do

      return
      end
      subroutine test0615 ( )

c*********************************************************************72
c
cc TEST0615 tests LERCH and LERCH_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision fx
      double precision fx2
      double precision lerch
      integer n_data
      integer s
      double precision z

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0615'
      write ( *, '(a)' ) '  LERCH computes the Lerch function.'
      write ( *, '(a)' ) '  LERCH_VALUES returns some tabulated values.'
      write ( *, '(a)' ) ' ' 
      write ( *, '(a)' ) 
     &  '       Z       S       A         Lerch           Lerch'
      write ( *, '(a)' ) 
     &  '                             Tabulated        Computed'
      write ( *, '(a)' ) ' '
     
      n_data = 0

10    continue

        call lerch_values ( n_data, z, s, a, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = lerch ( z, s, a )

        write ( *, '(2x,f8.4,2x,i4,2x,f8.4,2x,g14.6,2x,g14.6)' ) 
     &    z, s, a, fx, fx2

      go to 10

20    continue
     
      return
      end
      subroutine test062 ( )

c*********************************************************************72
c
cc TEST062 tests LOCK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer a(0:n)
      integer i

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST062'
      write ( *, '(a)' ) 
     &  '  LOCK counts the combinations on a button lock.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         I      LOCK(I)'
      write ( *, '(a)' ) ' '

      call lock ( n, a )

      do i = 0, n
        write ( *, '(2x,i8,2x,i10)' )  i, a(i)
      end do
     
      return
      end
      subroutine test0623 ( )

c*********************************************************************72
c
cc TEST0623 tests MEIXNER.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )
      integer test_num
      parameter ( test_num = 3 )

      double precision beta
      double precision beta_test(test_num)
      double precision c
      double precision c_test(test_num)
      integer i
      integer j
      integer test
      double precision v(0:n)
      double precision x

      save beta_test
      save c_test

      data beta_test / 0.5D+00, 1.0D+00, 2.0D+00 /
      data c_test / 0.125D+00, 0.25D+00, 0.5D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POLPAK_TEST0623:'
      write ( *, '(a)' ) '  MEIXNER evaluates Meixner polynomials.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '       N      BETA         C         X        M(N,BETA,C,X)'

      do test = 1, test_num

        beta = beta_test(test)
        c = c_test(test)

        do j = 0, 5

          x = dble ( j ) / 2.0D+00

          call meixner ( n, beta, c, x, v )

          write ( *, '(a)' ) ' '

          do i = 0, n

            write ( *, '(2x,i8,2x,f8.4,2x,f8.4,2x,f8.4,2x,g14.6)' ) 
     &        i, beta, c, x, v(i)

          end do

        end do

      end do

      return
      end
      subroutine test0625 ( )

c*********************************************************************72
c
cc TEST0625 tests MERTENS and MERTENS_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer c
      integer c2
      integer mertens
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0625'
      write ( *, '(a)' ) '  MERTENS computes the Mertens function.'
      write ( *, '(a)' ) '  MERTENS_VALUES returns some exact values.'
      write ( *, '(a)' ) ' ' 
      write ( *, '(a)' ) '         N     Exact   MERTENS(N)'
      write ( *, '(a)' ) ' '
     
      n_data = 0

10    continue

        call mertens_values ( n_data, n, c )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        c2 = mertens ( n )

        write ( *, '(2x,i8,2x,i10,2x,i10)' ) n, c, c2

      go to 10

20    continue
     
      return
      end
      subroutine test063 ( )

c*********************************************************************72
c
cc TEST063 tests MOEBIUS and MOEBIUS_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer c
      integer c2
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST063'
      write ( *, '(a)' ) '  MOEBIUS computes the Moebius function.'
      write ( *, '(a)' ) '  MOEBIUS_VALUES returns some exact values.'
      write ( *, '(a)' ) ' ' 
      write ( *, '(a)' ) '         N     Exact   MOEBIUS(N)'
      write ( *, '(a)' ) ' '
     
      n_data = 0

10    continue

        call moebius_values ( n_data, n, c )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        call moebius ( n, c2 )

        write ( *, '(2x,i8,2x,i10,2x,i10)' ) n, c, c2

      go to 10

20    continue
     
      return
      end
      subroutine test0635 ( )

c*********************************************************************72
c
cc TEST0635 tests MOTZKIN.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer a(0:n)
      integer i

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0635'
      write ( *, '(a)' ) 
     &  '  MOTZKIN computes the Motzkin numbers A(0:N).'
      write ( *, '(a)' ) '  A(N) counts the paths from (0,0) to (N,0).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         I         A(I)'
      write ( *, '(a)' ) ' '

      call motzkin ( n, a )

      do i = 0, n
        write ( *, '(2x,i8,2x,i10)' )  i, a(i)
      end do
     
      return
      end
      subroutine test064 ( )

c*********************************************************************72
c
cc TEST064 tests OMEGA and OMEGA_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer c
      integer c2
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST064'
      write ( *, '(a)' ) 
     &  '  OMEGA counts the distinct prime divisors of an integer N.'
      write ( *, '(a)' ) '  OMEGA_VALUES returns some exact values.'
      write ( *, '(a)' ) ' ' 
      write ( *, '(a)' ) '             N      Exact   OMEGA(N)'
      write ( *, '(a)' ) ' '
     
      n_data = 0

10    continue

        call omega_values ( n_data, n, c )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        call omega ( n, c2 )

        write ( *, '(2x,i12,2x,i10,2x,i10)' ) n, c, c2

      go to 10

20    continue
     
      return
      end
      subroutine test065 ( )

c*********************************************************************72
c
cc TEST065 tests PENTAGON_NUM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer p

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST065'
      write ( *, '(a)' ) 
     &  '  PENTAGON_NUM computes the pentagonal numbers.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         I      Pent(I)'
      write ( *, '(a)' ) ' '

      do n = 1, 10
        call pentagon_num ( n, p )
        write ( *, '(2x,i8,2x,i8)' ) n, p
      end do
     
      return
      end
      subroutine test066 ( )

c*********************************************************************72
c
cc TEST066 tests PHI and PHI_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer c
      integer c2
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST066'
      write ( *, '(a)' ) '  PHI computes the PHI function.'
      write ( *, '(a)' ) '  PHI_VALUES returns some exact values.'
      write ( *, '(a)' ) ' ' 
      write ( *, '(a)' ) '     N     Exact     PHI(N)'
      write ( *, '(a)' ) ' '
     
      n_data = 0

10    continue

        call phi_values ( n_data, n, c )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        call phi ( n, c2 )

        write ( *, '(2x,i8,2x,i10,2x,i10)' ) n, c, c2

      go to 10

20    continue
     
      return
      end
      subroutine test0665 ( )

c*********************************************************************72
c
cc TEST0665 tests POLY_BERNOULLI.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer b
      integer k
      integer n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0665'
      write ( *, '(a)' ) 
     &  '  POLY_BERNOULLI computes the poly-Bernoulli numbers'
      write ( *, '(a)' ) '  of negative index, B_n^(-k)'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     N     K    B_N^(-K)'
      write ( *, '(a)' ) ' '

      do k = 0, 6
        write ( *, '(a)' ) ' '
        do n = 0, 6

          call poly_bernoulli ( n, k, b )

          write ( *, '(2x,i4,2x,i4,2x,i12)' ) n, k, b

        end do
      end do

      return
      end
      subroutine test0667 ( )

c*********************************************************************72
c
cc TEST0667 tests POLY_COEF_COUNT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer degree
      integer dim
      integer poly_coef_count

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0667'
      write ( *, '(a)' ) 
     &  '  POLY_COEF_COUNT counts the number of coefficients'
      write ( *, '(a)' ) 
     &  '  in a polynomial of degree DEGREE and dimension DIM'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) ' Dimension    Degree     Count'

      do dim = 1, 10, 3
        write ( *, '(a)' ) ' '
        do degree = 0, 5
          write ( *, '(2x,i8,2x,i8,2x,i8)' ) 
     &     dim, degree, poly_coef_count ( dim, degree )
        end do
      end do
     
      return
      end
      subroutine test067 ( )

c*********************************************************************72
c
cc TEST067 tests PYRAMID_NUM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer pyramid_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST067'
      write ( *, '(a)' ) '  PYRAMID_NUM computes the pyramidal numbers.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         I    PYR(I)'
      write ( *, '(a)' ) ' '

      do n = 1, 10
        write ( *, '(2x,i8,2x,i8)' ) n, pyramid_num ( n )
      end do
     
      return
      end
      subroutine test0675

c*********************************************************************72
c
cc TEST0675 tests R8_ACOSH.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 November 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      integer i
      double precision r8_acosh
      double precision x
      double precision x2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0675'
      write ( *, '(a)' ) '  R8_ACOSH computes the inverse hyperbolic'
      write ( *, '(a)' ) '  cosine of a given value.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '       X        R8_ACOSH(X)     COSH(R8_ACOSH(X))'
      write ( *, '(a)' ) ' '

      do i = 0, 10
        x = 1.0D+00 + dble ( i ) / 5.0D+00
        a = r8_acosh ( x )
        x2 = dcosh ( a )
        write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) x, a, x2
      end do

      return
      end
      subroutine test004

c*********************************************************************72
c
cc TEST004 tests R8_ASINH.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 November 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      integer i
      double precision r8_asinh
      double precision x
      double precision x2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST004'
      write ( *, '(a)' ) '  R8_ASINH computes the inverse hyperbolic'
      write ( *, '(a)' ) '  sine of a given value.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '       X       R8_ASINH(X)     SINH(R8_ASINH(X))'
      write ( *, '(a)' ) ' '

      do i = 0, 10
        x = 1.0D+00 + dble ( i ) / 5.0D+00
        a = r8_asinh ( x )
        x2 = dsinh ( a )
        write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) x, a, x2
      end do

      return
      end
      subroutine test006

c*********************************************************************72
c
cc TEST006 tests R8_ATANH.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 November 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      integer i
      double precision r8_atanh
      double precision x
      double precision x2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST006'
      write ( *, '(a)' ) '  R8_ATANH computes the inverse hyperbolic'
      write ( *, '(a)' ) '  tangent of a given value.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '      X       R8_ATANH(X)     TANH(R8_ATANH(X))'
      write ( *, '(a)' ) ' '

      do i = -2, 9
        x = dble ( i ) / 10.0D+00
        a = r8_atanh ( x )
        x2 = dtanh ( a )
        write ( *, '(2x,g14.6,2x,g14.6,2x,g14.6)' ) x, a, x2
      end do

      return
      end
      subroutine test022 ( )

c*********************************************************************72
c
cc TEST022 tests R8_CHOOSE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision cnk
      integer k
      integer n
      double precision r8_choose

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST022'
      write ( *, '(a)' ) '  R8_CHOOSE evaluates C(N,K).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     N       K      CNK'
      write ( *, '(a)' ) ' '
 
      do n = 0, 4
        do k = 0, n
          cnk = r8_choose ( n, k )
          write ( *, '(2x,i8,2x,i8,2x,g14.6)' ) n, k, cnk
        end do
      end do
 
      return
      end
      subroutine test0685 ( )

c*********************************************************************72
c
cc TEST0685 tests R8_FACTORIAL_LOG and R8_FACTORIAL_LOG_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fn
      double precision fn2
      double precision r8_factorial_log
      integer n_data
      integer n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0685:'
      write ( *, '(a)' ) 
     &  '  R8_FACTORIAL_LOG evaluates the logarithm of the '
      write ( *, '(a)' ) '  factorial function.'
      write ( *, '(a)' ) 
     &  '  R8_FACTORIAL_LOG_VALUES returns some exact values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     N	   Exact F	 R8_FACTORIAL_LOG(N)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call r8_factorial_log_values ( n_data, n, fn )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fn2 = r8_factorial_log ( n )

        write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) n, fn, fn2

      go to 10

20    continue

      return
      end
      subroutine test06855

c*********************************************************************72
c
cc TEST06855 tests R8_GAMMA and GAMMA_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 January 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      double precision fx2
      integer n_data
      double precision r8_gamma
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06855:'
      write ( *, '(a)' ) '  R8_GAMMA evaluates the Gamma function.'
      write ( *, '(a)' ) '  GAMMA_VALUES returns some exact values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) 
     &  '      X         Gamma(X)                   Gamma(X)',
     &  '               DIFF'
      write ( *, '(a)' ) 
     &  '               (Tabulated)                (R8_GAMMA)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call gamma_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r8_gamma ( x )

        write ( *, '(2x,f8.4,2x,g24.16,2x,g24.16,2x,g10.4)' ) 
     &  x, fx, fx2, dabs ( fx - fx2 )

      go to 10

20    continue

      return
      end
      subroutine test06856 ( )

c*********************************************************************72
c
cc TEST06856 tests R8_PSI and PSI_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      double precision fx2
      integer n_data
      double precision r8_psi
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06856:'
      write ( *, '(a)' ) '  R8_PSI evaluates the Psi function.'
      write ( *, '(a)' ) '  PSI_VALUES returns some exact values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X         Psi(X)                     '
     &  // 'Psi(X)                 DIFF'
      write ( * , '(a)' ) 
     &  '               (Tabulated)                (R8_PSI)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call psi_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = r8_psi ( x )

        write ( *, '(2x,f8.4,2x,g24.16,2x,g24.16,2x,g10.4)' ) 
     &    x, fx, fx2, abs ( fx - fx2 )

      go to 10

20    continue

      return
      end
      subroutine test069 ( )

c*********************************************************************72
c
cc TEST069 tests SIGMA and SIGMA_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer c
      integer c2
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST069'
      write ( *, '(a)' ) '  SIGMA computes the SIGMA function.'
      write ( *, '(a)' ) '  SIGMA_VALUES returns some exact values.'
      write ( *, '(a)' ) ' ' 
      write ( *, '(a)' ) '     N     Exact   SIGMA(N)'
      write ( *, '(a)' ) ' '
     
      n_data = 0

10    continue

        call sigma_values ( n_data, n, c )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        call sigma ( n, c2 )

        write ( *, '(2x,i4,2i10)' ) n, c, c2

      go to 10

20    continue
     
      return
      end
      subroutine test0695 ( )

c*********************************************************************72
c
cc TEST0695 tests SIN_POWER_INT and SIN_POWER_INT_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision b
      double precision fx
      double precision fx2
      integer n
      integer n_data
      double precision sin_power_int

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0695:'
      write ( *, '(a)' ) '  SIN_POWER_INT returns values of '
      write ( *, '(a)' ) '  the integral of SIN(X)^N from A to B.'
      write ( *, '(a)' ) 
     &  '  SIN_POWER_INT_VALUES stores some selected values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '      A         B          N      Exact           Computed'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call sin_power_int_values ( n_data, a, b, n, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = sin_power_int ( a, b, n )

        write ( *, '(2x,f8.4,2x,f8.4,2x,i8,2x,g14.6,2x,g14.6)' ) 
     &    a, b, n, fx, fx2

      go to 10

20    continue

      return
      end
      subroutine test0696 ( )

c*********************************************************************72
c
cc TEST0696 tests SLICE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 August 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_max
      parameter ( dim_max = 5 )
      integer slice_max
      parameter ( slice_max = 8 )

      integer dim_num
      integer p(dim_max,slice_max)
      integer piece_num
      integer slice_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0696:'
      write ( *, '(a)' ) 
     &  '  SLICE determines the maximum number of pieces created'
      write ( *, '(a)' ) '  by SLICE_NUM slices in a DIM_NUM space.'

      do dim_num = 1, dim_max
        do slice_num = 1, slice_max
          call slice ( dim_num, slice_num, piece_num )
          p(dim_num,slice_num) = piece_num
        end do
      end do

      call i4mat_print ( dim_max, slice_max, p, '  Slice Array:' )

      return
      end
      subroutine test0697 ( )

c*********************************************************************72
c
cc TEST0697 tests SPHERICAL_HARMONIC and SPHERICAL_HARMONIC_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 20 )

      double precision c(0:n_max)
      integer l
      integer m
      integer n_data
      double precision phi
      double precision s(0:n_max)
      double precision theta
      double precision yi
      double precision yi2
      double precision yr
      double precision yr2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0697:'
      write ( *, '(a)' ) 
     &  '  SPHERICAL_HARMONIC evaluates spherical harmonic'
      write ( *, '(a)' ) '  functions.'
      write ( *, '(a)' ) 
     &  '  SPHERICAL_HARMONIC_VALUES returns some exact values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '       L       M   THETA    PHI     C              S'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call spherical_harmonic_values ( n_data, l, m, theta, phi, 
     &    yr, yi )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        call spherical_harmonic ( l, m, theta, phi, c, s )

        yr2 = c(l)
        yi2 = s(l)

        write ( *, '(2x,i8,2x,i6,2f8.4,2g14.6)' ) 
     &    l, m, theta, phi, yr,  yi
        write ( *, '(2x,8x,2x,6x,16x,  2g14.6)' )
     &                     yr2, yi2

      go to 10

20    continue

      return
      end
      subroutine test070 ( )

c*********************************************************************72
c
cc TEST070 tests STIRLING1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 8 )
      integer n
      parameter ( n = m )

      integer i
      integer s1(m,n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST070'
      write ( *, '(a)' ) '  STIRLING1: Stirling numbers of first kind.'
      write ( *, '(a,i8)' ) '  Get rows 1 through ', m
      write ( *, '(a)' ) ' '
     
      call stirling1 ( m, n, s1 )
     
      do i = 1, m
        write ( *, '(2x,i8,8i8)' ) i, s1(i,1:n)
      end do
     
      return
      end
      subroutine test071 ( )

c*********************************************************************72
c
cc TEST071 tests STIRLING2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 8 )
      integer n
      parameter ( n = m )

      integer i
      integer s2(m,n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST071'
      write ( *, '(a)' ) '  STIRLING2: Stirling numbers of second kind.'
      write ( *, '(a,i4)' ) '  Get rows 1 through ', m
      write ( *, '(a)' ) ' '
     
      call stirling2 ( m, n, s2 )
     
      do i = 1, m
        write ( *, '(2x,i8,8i8)' ) i, s2(i,1:n)
      end do
     
      return
      end
      subroutine test072 ( )

c*********************************************************************72
c
cc TEST072 tests TAU and TAU_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer c
      integer c2
      integer n
      integer n_data

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST072'
      write ( *, '(a)' ) '  TAU computes the Tau function.'
      write ( *, '(a)' ) '  TAU_VALUES returns some exact values.'
      write ( *, '(a)' ) ' ' 
      write ( *, '(a)' ) '         N  exact C(I)  computed C(I)'
      write ( *, '(a)' ) ' '
     
      n_data = 0

10    continue

        call tau_values ( n_data, n, c )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        call tau ( n, c2 )

        write ( *, '(2x,i8,2x,i10,2x,i10)' ) n, c, c2

      go to 10

20    continue
     
      return
      end
      subroutine test073 ( )

c*********************************************************************72
c
cc TEST073 tests TETRAHEDRON_NUM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer tetrahedron_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST073'
      write ( *, '(a)' ) 
     &  '  TETRAHEDRON_NUM computes the tetrahedron numbers.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         I    TETR(I)'
      write ( *, '(a)' ) ' '

      do n = 1, 10
        write ( *, '(2x,i8,2x,i8)' ) n, tetrahedron_num ( n )
      end do
     
      return
      end
      subroutine test074 ( )

c*********************************************************************72
c
cc TEST074 tests TRIANGLE_NUM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer triangle_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST074'
      write ( *, '(a)' ) 
     &  '  TRIANGLE_NUM computes the triangular numbers.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         I    TRI(I)'
      write ( *, '(a)' ) ' '
     
      do n = 1, 10
        write ( *, '(2x,i8,2x,i8)' ) n, triangle_num ( n )
      end do
     
      return
      end
      subroutine test076 ( )

c*********************************************************************72
c
cc TEST076 tests VIBONACCI.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 20 )
      integer n_time
      parameter ( n_time = 3 )

      integer i
      integer j
      integer seed
      integer v(n,n_time)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST076'
      write ( *, '(a)' ) '  VIBONACCI computes a Vibonacci sequence.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  Number of times we compute the series: ', n_time
      write ( *, '(a)' ) ' '

      seed = 123456789

      do j = 1, n_time
        call vibonacci ( n, seed, v(1,j) ) 
      end do

      do i = 1, n
        write ( *, '(2x,i8,2x,3i8)' ) i, v(i,1:n_time)
      end do
     
      return
      end
      subroutine test0773 ( )

c*********************************************************************72
c
cc TEST0773 tests ZERNIKE_POLY and ZERNIKE_POLY_COEF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 5 )

      double precision c(0:n_max)
      double precision cx1
      double precision cx2(0:n_max)
      integer m
      integer n
      double precision rho
      double precision z1
      double precision z2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0773'
      write ( *, '(a)' ) 
     &  '  ZERNIKE_POLY_COEF returns the coefficients of a'
      write ( *, '(a)' ) '  Zernike polynomial.'
      write ( *, '(a)' ) 
     &  '  ZERNIKE_POLY evaluates a Zernike polynomial directly.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Table of polynomial coefficients:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   N   M'
      write ( *, '(a)' ) ' '

      do n = 0, 5

        write ( *, '(a)' ) ' '

        do m = 0, n
          call zernike_poly_coef ( m, n, c )
          write ( *, '(2x,i2,2x,i2,2x,11f7.0)' ) n, m, c(0:n)
        end do

      end do

      rho = 0.987654321D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Z1: Compute polynomial coefficients,'
      write ( *, '(a)' ) '  then evaluate by Horner''s method;'
      write ( *, '(a)' ) '  Z2: Evaluate directly by recursion.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   N   M       Z1              Z2'
      write ( *, '(a)' ) ' '

      do n = 0, 5

        write ( *, '(a)' ) ' '

        do m = 0, n

          call zernike_poly_coef ( m, n, c )
          call r8poly_val_horner ( n, c, rho, z1 )

          call zernike_poly ( m, n, rho, z2 )

          write ( *, '(2x,i2,2x,i2,2x,g16.8,2x,g16.8)' ) n, m, z1, z2

        end do

      end do

      return
      end
      subroutine test077 ( )

c*********************************************************************72
c
cc TEST077 tests ZECKENDORF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m_max
      parameter ( m_max = 20 )

      integer i_list(m_max)
      integer f_list(m_max)
      integer m
      integer n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST077'
      write ( *, '(a)' ) 
     &  '  ZECKENDORF computes the Zeckendorf decomposition of'
      write ( *, '(a)' ) 
     &  '  an integer N into nonconsecutive Fibonacci numbers.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         N Sum M Parts'
      write ( *, '(a)' ) ' '

      do n = 1, 100

        call zeckendorf ( n, m_max, m, i_list, f_list )

        write ( *, '(2x,i8,2x,15i4)' ) n, f_list(1:m)

      end do

      return
      end
      subroutine test0775 ( )

c*********************************************************************72
c
cc TEST0775 tests ZERNIKE_POLY_COEF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      double precision c(0:n)
      integer m

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0775'
      write ( *, '(a)' ) '  ZERNIKE_POLY_COEF determines the Zernike'
      write ( *, '(a)' ) '  polynomial coefficients.'

      do m = 0, n

        call zernike_poly_coef ( m, n, c )
     
        call r8poly_print ( n, c, '  Zernike polynomial' )

      end do

      return
      end
      subroutine test078 ( )

c*********************************************************************72
c
cc TEST078 tests ZETA and ZETA_VALUES.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer n_data
      double precision n_real
      double precision z1
      double precision z2
      double precision zeta

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST078'
      write ( *, '(a)' ) '  ZETA computes the Zeta function.'
      write ( *, '(a)' ) '  ZETA_VALUES returns some exact values.'
      write ( *, '(a)' ) ' ' 
      write ( *, '(a)' ) '       N    exact Zeta    computed Zeta'
      write ( *, '(a)' ) ' '
     
      n_data = 0

10    continue

        call zeta_values ( n_data, n, z1 )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        n_real = dble ( n )

        z2 = zeta ( n_real )

        write ( *, '(2x,i8,2x,g20.12,2x,g20.12)' ) n, z1, z2

      go to 10

20    continue
     
      return
      end
