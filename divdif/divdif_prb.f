      program main

c*********************************************************************72
c
cc MAIN is the main program for DIVDIF_PRB.
c
c  Discussion:
c
c    DIVDIF_PRB calls a set of problems for DIVDIF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 February 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DIVDIF_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the DIVDIF library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
      call test06 ( )
      call test07 ( )
      call test08 ( )
      call test085 ( )
      call test09 ( )
      call test095 ( )

      call test10 ( )
      call test105 ( )
      call test11 ( )
      call test115 ( )
      call test12 ( )
      call test13 ( )
      call test14 ( )
      call test15 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DIVDIF_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests DIF_APPEND, DIF_ANTIDERIV, DIF_DERIV_TABLE, DIF_SHIFT_ZERO, DIF_VAL;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 October 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer maxtab
      parameter ( maxtab = 10 )

      double precision diftab(maxtab)
      double precision diftab2(maxtab)
      double precision diftab3(maxtab)
      integer ntab
      integer ntab2
      integer ntab3
      double precision xtab(maxtab)
      double precision xtab2(maxtab)
      double precision xtab3(maxtab)
      double precision xval
      double precision ytab(maxtab)
      double precision yval

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  For a divided difference polynomial:'
      write ( *, '(a)' )
     &  '  DATA_TO_DIF_DISPLAY sets up a difference table'
      write ( *, '(a)' ) '    and displays intermediate calculations;'
      write ( *, '(a)' ) '  DIF_APPEND appends a new data point;'
      write ( *, '(a)' ) '  DIF_ANTIDERIV computes the antiderivative;'
      write ( *, '(a)' ) '  DIF_DERIV_TABLE computes the derivative;'
      write ( *, '(a)' )
     &  '  DIF_SHIFT_ZERO shifts all the abscissas to 0;'
      write ( *, '(a)' ) '  DIF_VAL evaluates at a point.'
      write ( *, '(a)' ) ' '
c
c  Set XTAB, YTAB to X, X**2.
c
      ntab = 4
      call r8vec_indicator ( ntab, xtab )
      ytab(1:ntab) = xtab(1:ntab)**2

      call data_to_dif_display ( ntab, xtab, ytab, diftab )

      call dif_print ( ntab, xtab, diftab,
     &  '  The divided difference polynomial:' )
c
c  Append (5,25) to the table.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Append the data (5,25) to the table.'
      write ( *, '(a)' ) ' '

      xval = 5.0D+00
      yval = 25.0D+00

      call dif_append ( ntab, xtab, diftab, xval, yval, ntab, xtab,
     &  diftab )

      call dif_print ( ntab, xtab, diftab,
     &  '  The augmented divided difference polynomial:' )
c
c  Evaluate the polynomial at 2.5.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Evaluate the table at a point.'
      write ( *, '(a)' ) ' '

      xval = 2.5D+00

      call dif_val ( ntab, xtab, diftab, xval, yval )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,g14.6)' ) '  P( ', xval, ' ) = ', yval
c
c  Shift the base to zero.
c
      call dif_shift_zero ( ntab, xtab, diftab )

      call dif_print ( ntab, xtab, diftab,
     &  '  The table, rebased at 0:' )
c
c  Compute a table for the derivative.
c
      call dif_deriv_table ( ntab, xtab, diftab, ntab2, xtab2, diftab2 )

      call dif_print ( ntab2, xtab2, diftab2, '  The derivative:' )

      call dif_val ( ntab2, xtab2, diftab2, xval, yval )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,g14.6)' ) '  P''( ', xval, ' ) = ', yval
c
c  Compute the antiderivative.
c
      call dif_antideriv ( ntab, xtab, diftab, ntab3, xtab3, diftab3 )

      call dif_print ( ntab3, xtab3, diftab3, '  The antiderivative:' )

      call dif_val ( ntab3, xtab3, diftab3, xval, yval )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,g14.6)' ) 
     &  '  (Anti)P( ', xval, ' ) = ', yval

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests DATA_TO_DIF and DIF_VAL.
c
c  Discussion:
c
c    This test demonstrates how a divided difference table can be generated,
c    and then used to approximate a function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 February 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer maxtab
      parameter ( maxtab = 16 )

      double precision diftab(maxtab)
      double precision err
      double precision f02
      integer j
      integer ntab
      integer ntest
      parameter ( ntest = 1001 )
      double precision xtab(maxtab)
      double precision xtest
      double precision ytab(maxtab)
      double precision ytest

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  DATA_TO_DIF takes a set of (X,F(X)) data'
      write ( *, '(a)' ) '    and computes a divided difference table.'
      write ( *, '(a)' )
     &  '  DIF_VAL evaluates the corresponding polynomial'
      write ( *, '(a)' ) '    interpolant at an arbitrary point.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  By increasing the number of data points,'
      write ( *, '(a)' )
     &  '  the approximation should improve for a while.'
      write ( *, '(a)' ) '  However, our function is non-differentiable'
      write ( *, '(a)' )
     &  '  at one point, so the approximation begins to'
      write ( *, '(a)' ) '  misbehave rapidly.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Our interval is [-1,1].'
      write ( *, '(a)' ) '  Our function is F(X) = |X| + X/2 + X^2'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  We estimate the interpolation error using'
      write ( *, '(a)' ) '  1001 equally spaced sample points.'
      write ( *, '(a)' ) ' '

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Order  Interpolation Error'
      write ( *, '(a)' ) ' '

      do ntab = 1, maxtab

        call r8vec_even ( ntab, -1.0D+00, +1.0D+00, xtab )

        do j = 1, ntab
          ytab(j) = f02 ( xtab(j) )
        end do

        call data_to_dif ( ntab, xtab, ytab, diftab )

        do j = 1, ntest
          call r8vec_even_select ( ntest, -1.0D+00, +1.0D+00, j, xtest )
          call dif_val ( ntab, xtab, diftab, xtest, ytest )
          err = err + ( ytest - f02 ( xtest ) )**2
        end do

        err = 2.0D+00 * sqrt ( err ) / dble ( ntab )

        write ( *, ' ( 2x, i6, 2x, g14.6 )' ) ntab, err

      end do

      return
      end
      function f02 ( x )

c*********************************************************************72
c
cc F02 evaluates the function F(X) = |X| + X/2 - X^2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 February 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision f02
      double precision x

      f02 = abs ( x ) + x / 2.0D+00 - x * x

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests DIF_BASIS.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 February 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer ntab
      parameter ( ntab = 5 )

      double precision diftab(ntab,ntab)
      integer i
      integer nstep
      parameter ( nstep = 9 )
      double precision xhi
      double precision xlo
      double precision xtab(ntab)
      double precision xval
      double precision yval

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' )
     &  '  DIF_BASIS computes Lagrange basis polynomials'
      write ( *, '(a)' ) '    in difference form.'
      write ( *, '(a)' ) ' '
c
c  Set the base points.
c
      call r8vec_indicator ( ntab, xtab )

      call r8vec_print ( ntab, xtab, '  The base points:' )
c
c  Get the difference tables for the basis polynomials and print them.
c
      call dif_basis ( ntab, xtab, diftab )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  The table of difference vectors defining the basis'
      write ( *, '(a)' )
     &  '  polynomials.  Each column represents a polynomial.'
      write ( *, '(a)' ) ' '
      do i = 1, ntab
        write ( *, '(2x,5g14.6)' ) diftab(i,1:ntab)
      end do
c
c  Evaluate basis polynomial 3 at a set of points.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  Evaluate basis polynomial #3 at a set of points.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X        Y'
      write ( *, '(a)' ) ' '
      xhi = dble ( ntab )
      xlo = 1.0D+00

      do i = 1, nstep

        xval = ( dble ( nstep - i     ) * xlo
     &         + dble (         i - 1 ) * xhi )
     &         / dble ( nstep     - 1 )

        call dif_val ( ntab, xtab, diftab(1,3), xval, yval )

        write ( *, '(2g14.6)' ) xval, yval

      end do

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests DIF_ROOT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 February 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer maxtab
      parameter ( maxtab = 10 )

      double precision abserr
      double precision cubic
      external cubic
      integer iprint
      integer maxstp
      double precision relerr
      double precision xroot
      double precision xtry1
      double precision xtry2
c
c  Seek a root of the function F(X) = (X+3)*(X+1)*(X-1)
c  given starting estimates of 0.5 and 2.0.
c
      abserr = 0.00001D+00
      iprint = 1
      maxstp = 15
      relerr = 0.00001D+00
      xtry1 = 0.5D+00
      xtry2 = 2.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  DIF_ROOT seeks a zero of F(x).'
      write ( *, '(a)' ) '  F(X) = (X+3)*(X+1)*(X-1)'
      write ( *, '(a)' ) ' '

      call dif_root ( abserr, cubic, iprint, maxstp, maxtab,
     &  relerr, xroot, xtry1, xtry2 )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Estimated root = ', xroot
      write ( *, '(a,g14.6)' ) '  F(X) = ', cubic ( xroot )

      return
      end
      function cubic ( x )

c*********************************************************************72
c
cc CUBIC is the function whose root is sought.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 February 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision cubic
      double precision x

      cubic = ( x + 3.0D+00 ) * ( x + 1.0D+00 ) * ( x - 1.0D+00 )

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests DIF_TO_R8POLY and DIF_SHIFT_ZERO.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 February 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer maxtab
      parameter ( maxtab = 10 )

      double precision c(maxtab)
      double precision diftab1(maxtab)
      double precision diftab2(maxtab)
      integer i
      integer ntab
      double precision xtab1(maxtab)
      double precision xtab2(maxtab)
      double precision ytab1(maxtab)
      double precision ytab2(maxtab)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' )
     &  '  DIF_TO_R8POLY converts a difference table to a'
      write ( *, '(a)' ) '    polynomial;'
      write ( *, '(a)' ) '  DIF_SHIFT_ZERO shifts a divided difference '
      write ( *, '(a)' ) '    table to all zero abscissas;'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  These are equivalent operationsc'
      write ( *, '(a)' ) ' '
c
c  Set XTAB, YTAB to X, F(X)
c
      ntab = 4

      call r8vec_indicator ( ntab, xtab1 )

      do i = 1, ntab
        ytab1(i) = xtab1(i)**3 - 2.0D+00 * xtab1(i)**2
     &    + 3.0D+00 * xtab1(i) - 4.0D+00
      end do

      call r8vec_indicator ( ntab, xtab2 )

      do i = 1, ntab
        ytab2(i) = xtab2(i)**3 - 2.0D+00 * xtab2(i)**2
     &    + 3.0D+00 * xtab2(i) - 4.0D+00
      end do
c
c  Compute and display the finite difference table.
c
      call data_to_dif_display ( ntab, xtab1, ytab1, diftab1 )

      call data_to_dif_display ( ntab, xtab2, ytab2, diftab2 )
c
c  Examine corresponding polynomial.
c
      call dif_print ( ntab, xtab1, diftab1,
     &  '  The divided difference polynomial:' )
c
c  Shift to zero.
c
      call dif_shift_zero ( ntab, xtab1, diftab1 )

      call r8poly_print ( ntab, diftab1, '  Using DIF_SHIFT_ZERO' )
c
c  Shift to zero.
c
      call dif_to_r8poly ( ntab, xtab2, diftab2, c )

      call r8poly_print ( ntab, c, '  Using DIF_TO_R8POLY' )

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 tests R8POLY_*.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 February 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      integer i
      double precision poly_cof(n)
      double precision poly_cof2(n+1)
      double precision poly_cof3(n-1)
      double precision xval
      double precision yval
      double precision yval2
      double precision yval3

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' )
     &  '  R8POLY_ANT_COF computes the coefficients of the'
      write ( *, '(a)' ) '    antiderivative of a polynomial;'
      write ( *, '(a)' )
     &  '  R8POLY_ANT_VAL evaluates the antiderivative of'
      write ( *, '(a)' ) '    a polynomial;'
      write ( *, '(a)' )
     &  '  R8POLY_DER_COF computes the coefficients of the'
      write ( *, '(a)' ) '    derivative of a polynomial;'
      write ( *, '(a)' )
     &  '  R8POLY_DER_VAL evaluates the derivative of'
      write ( *, '(a)' ) '    a polynomial;'
      write ( *, '(a)' ) '  R8POLY_PRINT prints a polynomial;'
      write ( *, '(a)' ) '  R8POLY_VAL evaluates a polynomial.'

      do i = 1, n
        poly_cof(i) = dble ( i )
      end do

      call r8poly_print ( n, poly_cof, '  Our initial polynomial:' )

      call r8poly_ant_cof ( n, poly_cof, poly_cof2 )

      call r8poly_print ( n+1, poly_cof2,
     &  '  The antiderivative polynomial:' )

      call r8poly_der_cof ( n, poly_cof, poly_cof3 )

      call r8poly_print ( n-1, poly_cof3,
     &  '  The derivative polynomial:' )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Evaluate the polynomial, antiderivative and'
      write ( *, '(a)' )
     &  '  derivative, using only the original polynomial'
      write ( *, '(a)' ) '  coefficients:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '     X             P(X)         Anti_P(X)     P''(X)'
      write ( *, '(a)' ) ' '

      do i = 0, 2

        xval = dble ( i )

        call r8poly_val_horner ( n, poly_cof, xval, yval )

        call r8poly_ant_val ( n, poly_cof, xval, yval2 )

        call r8poly_der_val ( n, poly_cof, xval, yval3 )

        write ( *, '(4g14.6)' ) xval, yval, yval2, yval3

      end do

      return
      end
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07 tests R8POLY_BASIS.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 February 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer ntab
      parameter ( ntab = 5 )

      double precision polcof(ntab,ntab)
      integer i
      integer j
      integer nstep
      parameter ( nstep = 9 )
      double precision xhi
      double precision xlo
      double precision xtab(ntab)
      double precision xval
      double precision yval

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07'
      write ( *, '(a)' )
     &  '  R8POLY_BASIS computes Lagrange basis polynomials'
      write ( *, '(a)' ) '  in standard form.'
      write ( *, '(a)' ) ' '
c
c  Set the base points.
c
      call r8vec_indicator ( ntab, xtab )
c
c  Get the difference tables for the basis polynomials and print them.
c
      call r8poly_basis ( ntab, xtab, polcof )

      do i = 1, ntab
        write ( *, '(2x,5g14.6)' ) ( polcof(i,j), j = 1, ntab )
      end do
c
c  Print basis polynomial 3 in polynomial form.
c
      call r8poly_print ( ntab, polcof(1,3),
     &  '  Basis polynomial 3 in standard form:' )
c
c  Evaluate basis polynoimial 3 at a set of points.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  Evaluate basis polynomial 3 at a set of points.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      X        Y'
      write ( *, '(a)' ) ' '
      xhi = dble ( ntab )
      xlo = 1.0D+00

      do i = 1, nstep

        xval = ( dble ( nstep - i     ) * xlo
     &         + dble (         i - 1 ) * xhi )
     &         / dble ( nstep     - 1 )

        call r8poly_val_horner ( ntab, polcof(1,3), xval, yval )

        write ( *, '(2x,2g14.6)' ) xval, yval

      end do

      return
      end
      subroutine test08 ( )

c*********************************************************************72
c
cc TEST08 tests R8POLY_SHIFT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 February 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      integer i
      double precision poly_cof(n)
      double precision scale
      double precision shift

      scale = 2.0D+00
      shift = +3.0D+00
      poly_cof(1) = +6.0D+00
      poly_cof(2) = -1.0D+00
      poly_cof(3) =  2.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08'
      write ( *, '(a)' )
     &  '  R8POLY_SHIFT shifts polynomial coefficients.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Polynomial coefficients for argument X'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(i5,g14.6)' ) i, poly_cof(i)
      end do

      call r8poly_shift ( scale, shift, n, poly_cof )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  SCALE = ', scale
      write ( *, '(a,g14.6)' ) '  SHIFT = ', shift
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Polynomial coefficients for argument '
      write ( *, '(a)' ) '    Z = SCALE * X + SHIFT'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(i5,g14.6)' ) i, poly_cof(i)
      end do

      return
      end
      subroutine test085 ( )

c*********************************************************************72
c
cc TEST085 tests LAGRANGE_VAL;
c
c  Discussion:
c
c    This test demonstrates how a divided difference table can be generated,
c    and then used to approximate a function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 February 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer maxtab
      parameter ( maxtab = 16 )

      double precision err
      double precision f02
      integer j
      integer ntab
      integer ntest
      parameter ( ntest = 1001 )
      double precision xtab(maxtab)
      double precision xtest
      double precision ytab(maxtab)
      double precision ytest

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST085'
      write ( *, '(a)' )
     &  '  LAGRANGE_VAL uses naive Lagrange interpolation'
      write ( *, '(a)' )
     &  '  to compute the polynomial interpolant to data.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  By increasing the number of data points,'
      write ( *, '(a)' )
     &  '  the approximation should improve for a while.'
      write ( *, '(a)' ) '  However, our function is non-differentiable'
      write ( *, '(a)' )
     &  '  at one point, so the approximation begins to'
      write ( *, '(a)' ) '  misbehave rapidly.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Our interval is [-1,1].'
      write ( *, '(a)' ) '  Our function is F(X) = |X| + X/2 + X^2'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  We estimate the interpolation error using'
      write ( *, '(a)' ) '  1001 equally spaced sample points.'
      write ( *, '(a)' ) ' '

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Order  Interpolation Error'
      write ( *, '(a)' ) ' '

      do ntab = 1, maxtab

        call r8vec_even ( ntab, -1.0D+00, +1.0D+00, xtab )
c
c  Now set up some data.
c
        do j = 1, ntab
          ytab(j) = f02 ( xtab(j) )
        end do

        err = 0.0D+00

        do j = 1, ntest
          call r8vec_even_select ( ntest, -1.0D+00, +1.0D+00, j, xtest )
          call lagrange_val ( ntab, xtab, ytab, xtest, ytest )
          err = err + ( ytest - f02 ( xtest ) )**2
        end do

        err = 2.0D+00 * sqrt ( err ) / dble ( ntab )

        write ( *, ' ( 2x, i6, 2x, g14.6 )' ) ntab, err

      end do

      return
      end
      subroutine test09 ( )

c*********************************************************************72
c
cc TEST09 tests LAGRANGE_RULE;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 February 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer maxtab
      parameter ( maxtab = 8 )

      integer j
      integer ntab
      double precision wtab(maxtab)
      double precision xtab(maxtab)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST09'
      write ( *, '(a)' )
     &  '  LAGRANGE_RULE computes Lagrange interpolation formulas;'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Lagrange Interpolation Rules on [-1,1]'
      write ( *, '(a)' ) '  using equally spaced abscissas.'

      do ntab = 1, maxtab

        call r8vec_even ( ntab, -1.0D+00, +1.0D+00, xtab )

        call lagrange_rule ( ntab, xtab, wtab )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '      Abscissa       Weight'
        write ( *, '(a)' ) ' '
        do j = 1, ntab
          write ( *, '(2x,i3,g14.6,2x,g14.6)' ) j, xtab(j), wtab(j)
        end do

      end do

      return
      end
      subroutine test095 ( )

c*********************************************************************72
c
cc TEST095 tests LAGRANGE_RULE and LAGRANGE_SUM.
c
c  Discussion:
c
c    This test demonstrates how a divided difference table can be generated,
c    and then used to approximate a function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 February 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer maxtab
      parameter ( maxtab = 16 )

      double precision err
      double precision f02
      integer j
      integer ntab
      integer ntest
      parameter ( ntest = 1001 )
      double precision wtab(maxtab)
      double precision xtab(maxtab)
      double precision xtest
      double precision ytab(maxtab)
      double precision ytest

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST095'
      write ( *, '(a)' )
     &  '  LAGRANGE_RULE sets the weights for a Lagrange rule.'
      write ( *, '(a)' )
     &  '  LAGRANGE_SUM uses the rule to compute the Lagrange'
      write ( *, '(a)' ) '  interpolant to data at a given point.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  For this test, the data abscissas are '
      write ( *, '(a)' ) '  equally spaced.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  By increasing the number of data points,'
      write ( *, '(a)' )
     &  '  the approximation should improve for a while.'
      write ( *, '(a)' ) '  However, our function is non-differentiable'
      write ( *, '(a)' )
     &  '  at one point, so the approximation begins to'
      write ( *, '(a)' ) '  misbehave rapidly.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Our interval is [-1,1].'
      write ( *, '(a)' ) '  Our function is F(X) = |X| + X/2 + X^2'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  We estimate the interpolation error using'
      write ( *, '(a)' ) '  1001 equally spaced sample points.'
      write ( *, '(a)' ) ' '

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Order  Interpolation Error'
      write ( *, '(a)' ) ' '

      do ntab = 1, maxtab

        call r8vec_even ( ntab, -1.0D+00, +1.0D+00, xtab )

        call lagrange_rule ( ntab, xtab, wtab )
c
c  Now set up some data.
c
        do j = 1, ntab
          ytab(j) = f02 ( xtab(j) )
        end do

        err = 0.0D+00
        do j = 1, ntest
          call r8vec_even_select ( ntest, -1.0D+00, +1.0D+00, j, xtest )
          call lagrange_sum ( ntab, xtab, wtab, ytab, xtest, ytest )
          err = err + ( ytest - f02 ( xtest ) )**2
        end do

        err = 2.0D+00 * sqrt ( err ) / dble ( ntab )

        write ( *, ' ( 2x, i6, 2x, g14.6 )' ) ntab, err

      end do

      return
      end
      subroutine test10 ( )

c*********************************************************************72
c
cc TEST10 tests LAGRANGE_RULE;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 February 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer maxtab
      parameter ( maxtab = 8 )

      integer j
      integer ntab
      double precision wtab(maxtab)
      double precision xtab(maxtab)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST10'
      write ( *, '(a)' )
     &  '  LAGRANGE_RULE computes Lagrange interpolation formulas;'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Lagrange Interpolation Rules on [-1,1]'
      write ( *, '(a)' ) '  using Chebyshev T abscissas.'

      do ntab = 1, maxtab

        call cheby_t_zero ( ntab, xtab )

        call lagrange_rule ( ntab, xtab, wtab )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '      Abscissa       Weight'
        write ( *, '(a)' ) ' '
        do j = 1, ntab
          write ( *, '(2x,i3,g14.6,2x,g14.6)' ) j, xtab(j), wtab(j)
        end do

      end do

      return
      end
      subroutine test105 ( )

c*********************************************************************72
c
cc TEST105 tests LAGRANGE_RULE and LAGRANGE_SUM.
c
c  Discussion:
c
c    This test demonstrates how a divided difference table can be generated,
c    and then used to approximate a function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 February 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer maxtab
      parameter ( maxtab = 16 )

      double precision err
      double precision f02
      integer j
      integer ntab
      integer ntest
      parameter ( ntest = 1001 )
      double precision wtab(maxtab)
      double precision xtab(maxtab)
      double precision xtest
      double precision ytab(maxtab)
      double precision ytest

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST105'
      write ( *, '(a)' )
     &  '  LAGRANGE_RULE sets the weights for a Lagrange rule.'
      write ( *, '(a)' )
     &  '  LAGRANGE_SUM uses the rule to compute the Lagrange'
      write ( *, '(a)' ) '  interpolant to data at a given point.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  For this test, the data abscissas are the'
      write ( *, '(a)' ) '  zeroes of the Chebyshev T polynomials.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  By increasing the number of data points,'
      write ( *, '(a)' )
     &  '  the approximation should improve for a while.'
      write ( *, '(a)' ) '  However, our function is non-differentiable'
      write ( *, '(a)' )
     &  '  at one point, so the approximation begins to'
      write ( *, '(a)' ) '  misbehave rapidly.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Our interval is [-1,1].'
      write ( *, '(a)' ) '  Our function is F(X) = |X| + X/2 + X^2'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  We estimate the interpolation error using'
      write ( *, '(a)' ) '  1001 equally spaced sample points.'
      write ( *, '(a)' ) ' '

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Order  Interpolation Error'
      write ( *, '(a)' ) ' '

      do ntab = 1, maxtab

        call cheby_t_zero ( ntab, xtab )

        call lagrange_rule ( ntab, xtab, wtab )
c
c  Now set up some data.
c
        do j = 1, ntab
          ytab(j) = f02 ( xtab(j) )
        end do

        err = 0.0D+00
        do j = 1, ntest
          call r8vec_even_select ( ntest, -1.0D+00, +1.0D+00, j, xtest )
          call lagrange_sum ( ntab, xtab, wtab, ytab, xtest, ytest )
          err = err + ( ytest - f02 ( xtest ) )**2
        end do

        err = 2.0D+00 * sqrt ( err ) / dble ( ntab )

        write ( *, ' ( 2x, i6, 2x, g14.6 )' ) ntab, err

      end do

      return
      end
      subroutine test11 ( )

c*********************************************************************72
c
cc TEST11 tests LAGRANGE_RULE;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 February 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer maxtab
      parameter ( maxtab = 8 )

      integer j
      integer ntab
      double precision wtab(maxtab)
      double precision xtab(maxtab)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST11'
      write ( *, '(a)' )
     &  '  LAGRANGE_RULE computes Lagrange interpolation formulas;'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Lagrange Interpolation Rules on [-1,1]'
      write ( *, '(a)' ) '  using Chebyshev U abscissas.'

      do ntab = 1, maxtab

        call cheby_u_zero ( ntab, xtab )

        call lagrange_rule ( ntab, xtab, wtab )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '      Abscissa       Weight'
        write ( *, '(a)' ) ' '
        do j = 1, ntab
          write ( *, '(2x,i3,g14.6,2x,g14.6)' ) j, xtab(j), wtab(j)
        end do

      end do

      return
      end
      subroutine test115 ( )

c*********************************************************************72
c
cc TEST115 tests LAGRANGE_RULE and LAGRANGE_SUM.
c
c  Discussion:
c
c    This test demonstrates how a divided difference table can be generated,
c    and then used to approximate a function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 February 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer maxtab
      parameter ( maxtab = 16 )

      double precision err
      double precision f02
      integer j
      integer ntab
      integer ntest
      parameter ( ntest = 1001 )
      double precision wtab(maxtab)
      double precision xtab(maxtab)
      double precision xtest
      double precision ytab(maxtab)
      double precision ytest

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST115'
      write ( *, '(a)' )
     &  '  LAGRANGE_RULE sets the weights for a Lagrange rule.'
      write ( *, '(a)' )
     &  '  LAGRANGE_SUM uses the rule to compute the Lagrange'
      write ( *, '(a)' ) '  interpolant to data at a given point.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  For this test, the data abscissas are the'
      write ( *, '(a)' ) '  zeroes of the Chebyshev U polynomials.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  By increasing the number of data points,'
      write ( *, '(a)' )
     &  '  the approximation should improve for a while.'
      write ( *, '(a)' ) '  However, our function is non-differentiable'
      write ( *, '(a)' )
     &  '  at one point, so the approximation begins to'
      write ( *, '(a)' ) '  misbehave rapidly.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Our interval is [-1,1].'
      write ( *, '(a)' ) '  Our function is F(X) = |X| + X/2 + X^2'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  We estimate the interpolation error using'
      write ( *, '(a)' ) '  1001 equally spaced sample points.'
      write ( *, '(a)' ) ' '

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Order  Interpolation Error'
      write ( *, '(a)' ) ' '

      do ntab = 1, maxtab

        call cheby_u_zero ( ntab, xtab )

        call lagrange_rule ( ntab, xtab, wtab )
c
c  Now set up some data.
c
        do j = 1, ntab
          ytab(j) = f02 ( xtab(j) )
        end do

        err = 0.0D+00
        do j = 1, ntest
          call r8vec_even_select ( ntest, -1.0D+00, +1.0D+00, j, xtest )
          call lagrange_sum ( ntab, xtab, wtab, ytab, xtest, ytest )
          err = err + ( ytest - f02 ( xtest ) )**2
        end do

        err = 2.0D+00 * sqrt ( err ) / dble ( ntab )

        write ( *, ' ( 2x, i6, 2x, g14.6 )' ) ntab, err

      end do

      return
      end
      subroutine test12 ( )

c*********************************************************************72
c
cc TEST12 tests NCC_RULE;
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 February 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer norder
      parameter ( norder = 8 )

      integer i
      double precision weight(norder)
      double precision xtab(norder)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST12'
      write ( *, '(a)' )
     &  '  NCC_RULE computes closed Newton Cotes formulas'
      write ( *, '(a)' ) '  for quadrature (approximate integration).'
      write ( *, '(a)') ' '

      call ncc_rule ( norder, xtab, weight )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Newton-Cotes Closed Quadrature Rule:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      Abscissa       Weight'
      write ( *, '(a)' ) ' '

      do i = 1, norder
        write ( *, '(2x,i3,2g14.6)' ) i, xtab(i), weight(i)
      end do

      return
      end
      subroutine test13 ( )

c*********************************************************************72
c
cc TEST13 tests NCO_RULE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 February 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer norder
      parameter ( norder = 8 )

      integer i
      double precision weight(norder)
      double precision xtab(norder)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST13'
      write ( *, '(a)' )
     &  '  NCO_RULE computes open Newton Cotes formulas'
      write ( *, '(a)' ) '  for quadrature (approximate integration).'
      write ( *, '(a)' ) ' '

      call nco_rule ( norder, xtab, weight )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Newton-Cotes Open Quadrature Rule:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      Abscissa       Weight'
      write ( *, '(a)' ) ' '

      do i = 1, norder
        write ( *, '(2x,i3,2g14.6)' ) i, xtab(i), weight(i)
      end do

      return
      end
      subroutine test14 ( )

c*********************************************************************72
c
cc TEST14 tests ROOTS_TO_DIF and DIF_TO_R8POLY.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 February 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer maxroots
      parameter ( maxroots = 4 )

      double precision c(maxroots+1)
      double precision diftab(maxroots+1)
      integer nroots
      integer ntab
      double precision roots(maxroots)
      double precision xtab(maxroots+1)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST14'
      write ( *, '(a)' )
     &  '  ROOTS_TO_DIF computes the divided difference'
      write ( *, '(a)' ) '    polynomial with given roots;'
      write ( *, '(a)' )
     &  '  DIF_TO_R8POLY converts it to a standard form'
      write ( *, '(a)' ) '    polynomial.'
      write ( *, '(a)' ) ' '

      nroots = 1
      roots(1) = 3.0D+00
      call r8vec_print ( nroots, roots, '  The roots:' )

      call roots_to_dif ( nroots, roots, ntab, xtab, diftab )

      call dif_to_r8poly ( ntab, xtab, diftab, c )

      call r8poly_print ( ntab, c, '  The polynomial:' )

      nroots = 2
      roots(1:2) = (/ 3.0D+00, 1.0D+00 /)
      call r8vec_print ( nroots, roots, '  The roots:' )

      call roots_to_dif ( nroots, roots, ntab, xtab, diftab )

      call dif_to_r8poly ( ntab, xtab, diftab, c )

      call r8poly_print ( ntab, c, '  The polynomial:' )

      nroots = 3
      roots(1:3) = (/ 3.0D+00, 1.0D+00, 2.0D+00 /)
      call r8vec_print ( nroots, roots, '  The roots:' )

      call roots_to_dif ( nroots, roots, ntab, xtab, diftab )

      call dif_to_r8poly ( ntab, xtab, diftab, c )

      call r8poly_print ( ntab, c, '  The polynomial:' )

      nroots = 4
      roots(1:4) = (/ 3.0D+00, 1.0D+00, 2.0D+00, 4.0D+00 /)
      call r8vec_print ( nroots, roots, '  The roots:' )

      call roots_to_dif ( nroots, roots, ntab, xtab, diftab )

      call dif_to_r8poly ( ntab, xtab, diftab, c )

      call r8poly_print ( ntab, c, '  The polynomial:' )

      return
      end
      subroutine test15 ( )

c*********************************************************************72
c
cc TEST15 tests ROOTS_TO_R8POLY.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 February 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer maxroot
      parameter ( maxroot = 5 )

      double precision c(maxroot+1)
      integer nc
      integer nroots
      double precision roots(maxroot)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST15'
      write ( *, '(a)' )
     &  '  ROOTS_TO_R8POLY computes polynomial coefficients'
      write ( *, '(a)' ) '    from roots.'
      write ( *, '(a)' ) ' '

      nroots = 1
      roots(1) = 3.0D+00
      call r8vec_print ( nroots, roots, '  The roots:' )

      call roots_to_r8poly ( nroots, roots, nc, c )

      call r8poly_print ( nc, c, '  The polynomial:' )

      nroots = 2
      roots(1) = 3.0D+00
      roots(2) = 1.0D+00
      call r8vec_print ( nroots, roots, '  The roots:' )

      call roots_to_r8poly ( nroots, roots, nc, c )

      call r8poly_print ( nc, c, '  The polynomial:' )

      nroots = 3
      roots(1) = 3.0D+00
      roots(2) = 1.0D+00
      roots(3) = 2.0D+00
      call r8vec_print ( nroots, roots, '  The roots:' )

      call roots_to_r8poly ( nroots, roots, nc, c )

      call r8poly_print ( nc, c, '  The polynomial:' )

      nroots = 4
      roots(1) = 3.0D+00
      roots(2) = 1.0D+00
      roots(3) = 2.0D+00
      roots(4) = 4.0D+00
      call r8vec_print ( nroots, roots, '  The roots:' )

      call roots_to_r8poly ( nroots, roots, nc, c )

      call r8poly_print ( nc, c, '  The polynomial:' )

      return
      end
