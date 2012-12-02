      subroutine cheby_t_zero ( n, z )

c*********************************************************************72
c
cc CHEBY_T_ZERO returns zeroes of the Chebyshev polynomial T(N)(X).
c
c  Discussion:
c
c    The I-th zero of T(N)(X) is cos((2*I-1)*PI/(2*N)), I = 1 to N
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
c  Parameters:
c
c    Input, integer N, the order of the polynomial.
c
c    Output, double precision Z(N), the zeroes of T(N)(X).
c
      implicit none

      integer n

      double precision angle
      integer i
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision z(n)

      do i = 1, n
        angle = dble ( 2 * i - 1 ) * pi
     &        / dble ( 2 * n     )
        z(i) = cos ( angle )
      end do

      return
      end
      subroutine cheby_u_zero ( n, z )

c*********************************************************************72
c
cc CHEBY_U_ZERO returns zeroes of the Chebyshev polynomial U(N)(X).
c
c  Discussion:
c
c    The I-th zero of U(N)(X) is cos((I-1)*PI/(N-1)), I = 1 to N
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
c  Parameters:
c
c    Input, integer N, the order of the polynomial.
c
c    Output, double precision Z(N), the zeroes of U(N)(X).
c
      implicit none

      integer n

      double precision angle
      integer i
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision z(n)

      do i = 1, n
        angle = dble ( i     ) * pi
     &        / dble ( n + 1 )
        z(i) = cos ( angle )
      end do

      return
      end
      subroutine data_to_dif ( ntab, xtab, ytab, diftab )

c*********************************************************************72
c
cc DATA_TO_DIF sets up a divided difference table from raw data.
c
c  Discussion:
c
c    Space can be saved by using a single array for both the DIFTAB and
c    YTAB dummy parameters.  In that case, the divided difference table will
c    overwrite the Y data without interfering with the computation.
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
c  Reference:
c
c    Carl deBoor,
c    A Practical Guide to Splines,
c    Springer, 2001,
c    ISBN: 0387953663,
c    LC: QA1.A647.v27.
c
c  Parameters:
c
c    Input, integer NTAB, the number of pairs of points
c    (XTAB(I),YTAB(I)) which are to be used as data.  The
c    number of entries to be used in DIFTAB, XTAB and YTAB.
c
c    Input, double precision XTAB(NTAB), the X values at which data was taken.
c    These values must be distinct.
c
c    Input, double precision YTAB(NTAB), the corresponding Y values.
c
c    Output, double precision DIFTAB(NTAB), the divided difference coefficients
c    corresponding to the input (XTAB,YTAB).
c
      implicit none

      integer ntab

      double precision diftab(ntab)
      integer i
      integer j
      logical r8vec_distinct
      double precision xtab(ntab)
      double precision ytab(ntab)

      if ( .not. r8vec_distinct ( ntab, xtab ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DATA_TO_DIF - Fatal errorc'
        write ( *, '(a)' ) '  Two entries of XTAB are equalc'
        stop
      end if
c
c  Copy the data values into DIFTAB.
c
      do i = 1, ntab
        diftab(i) = ytab(i)
      end do
c
c  Compute the divided differences.
c
      do i = 2, ntab
        do j = ntab, i, -1

          diftab(j) = ( diftab(j) - diftab(j-1) )
     &              / ( xtab(j)   - xtab(j+1-i) )

        end do
      end do

      return
      end
      subroutine data_to_dif_display ( ntab, xtab, ytab, diftab )

c*********************************************************************72
c
cc DATA_TO_DIF_DISPLAY computes a divided difference table and shows how.
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
c  Parameters:
c
c    Input, integer NTAB, the number of pairs of points
c    (XTAB(I),YTAB(I)) which are to be used as data.  The
c    number of entries to be used in DIFTAB, XTAB and YTAB.
c
c    Input, double precision XTAB(NTAB), the X values at which data was taken.
c
c    Input, double precision YTAB(NTAB), the corresponding Y values.
c
c    Output, double precision DIFTAB(NTAB), the divided difference
c    coefficients corresponding to the input (XTAB,YTAB).
c
      implicit none

      integer ntab

      double precision diftab(ntab)
      integer i
      integer j
      logical r8vec_distinct
      double precision xtab(ntab)
      double precision ytab(ntab)

      if ( .not. r8vec_distinct ( ntab, xtab ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DATA_TO_DIF_DISPLAY - Fatal error!'
        write ( *, '(a)' ) '  Two entries of XTAB are equal.'
        stop
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Divided difference table:'
      write ( *, '(a)' ) ' '
      write ( *, '(6x,5g14.6)' ) ( xtab(j), j = 1, ntab )
      write ( *, '(a)' ) ' '
      write ( *, '(2x,i3,1x,5g14.6)' ) 0, ( ytab(j), j = 1, ntab )
c
c  Copy the data values into DIFTAB.
c
      do i = 1, ntab
        diftab(i) = ytab(i)
      end do
c
c  Compute the divided differences.
c
      do i = 2, ntab
        do j = ntab, i, -1

          diftab(j) = ( diftab(j) - diftab(j-1) )
     &              / ( xtab(j)   - xtab(j+1-i) )

        end do

        write ( *, '(2x,i3,1x,5g14.6)' ) i-1, ( diftab(j), j = i, ntab )

      end do

      return
      end
      subroutine data_to_r8poly ( ntab, xtab, ytab, c )

c*********************************************************************72
c
cc DATA_TO_R8POLY computes the coefficients of a polynomial interpolating data.
c
c  Discussion:
c
c    Space can be saved by using a single array for both the C and
c    YTAB parameters.  In that case, the coefficients will
c    overwrite the Y data without interfering with the computation.
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
c  Reference:
c
c    Carl deBoor,
c    A Practical Guide to Splines,
c    Springer, 2001,
c    ISBN: 0387953663,
c    LC: QA1.A647.v27.
c
c  Parameters:
c
c    Input, integer NTAB, the number of data points.
c
c    Input, double precision XTAB(NTAB), YTAB(NTAB), the data values.
c
c    Output, double precision C(NTAB), the coefficients of the
c    polynomial that passes through the data (XTAB,YTAB).  C(1) is the
c    constant term.
c
      implicit none

      integer ntab

      double precision c(ntab)
      double precision diftab(ntab)
      logical r8vec_distinct
      double precision xtab(ntab)
      double precision ytab(ntab)

      if ( .not. r8vec_distinct ( ntab, xtab ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DATA_TO_POLY - Fatal errorc'
        write ( *, '(a)' ) '  Two entries of XTAB are equalc'
        stop
      end if

      call data_to_dif ( ntab, xtab, ytab, diftab )

      call dif_to_r8poly ( ntab, xtab, diftab, c )

      return
      end
      subroutine dif_antideriv ( ntab, xtab, diftab, ntab2, xtab2,
     &  diftab2 )

c*********************************************************************72
c
cc DIF_ANTIDERIV computes the antiderivative of a divided difference polynomial.
c
c  Discussion:
c
c    The routine uses the divided difference representation
c    of a polynomial to compute the divided difference representation
c    of the antiderivative of the polynomial.
c
c    The antiderivative of a polynomial P(X) is any polynomial Q(X)
c    with the property that d/dX Q(X) = P(X).
c
c    This routine chooses the antiderivative whose constant term is zero.
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
c  Parameters:
c
c    Input, integer NTAB, the size of the difference table.
c
c    Input, double precision XTAB(NTAB), the abscissas of the
c    difference table.
c
c    Input, double precision DIFTAB(NTAB), the difference table.
c
c    Input, integer NTAB2, the size of the difference table for the
c    antiderivative, which will be NTAB+1.
c
c    Input, double precision XTAB2(NTAB2), the abscissas of the
c    difference table for the antiderivative.
c
c    Input, double precision DIFTAB2(NTAB2), the difference table
c    for the antiderivative.
c
      implicit none

      integer ntab

      double precision diftab(ntab)
      double precision diftab1(ntab)
      double precision diftab2(ntab+1)
      integer i
      integer ntab2
      double precision xtab(ntab)
      double precision xtab1(ntab)
      double precision xtab2(ntab+1)
c
c  Copy the input data, and shift the abscissas to zero.
c
      do i = 1, ntab
        xtab1(i) = xtab(i)
        diftab1(i) = diftab(i)
      end do

      call dif_shift_zero ( ntab, xtab1, diftab1 )
c
c  Append a final zero to XTAB.
c
      ntab2 = ntab + 1
      do i = 1, ntab2
        xtab2(i) = 0.0D+00
      end do
c
c  Get the antiderivative of the standard form polynomial.
c
      call r8poly_ant_cof ( ntab, diftab1, diftab2 )

      return
      end
      subroutine dif_append ( ntab, xtab, diftab, xval, yval, ntab2,
     &  xtab2, diftab2 )

c*********************************************************************72
c
cc DIF_APPEND adds a pair of data values to a divided difference table.
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
c  Parameters:
c
c    Input, integer NTAB, the size of the difference table.
c
c    Input, double precision XTAB(NTAB), the X data values.
c
c    Input, double precision DIFTAB(NTAB), the difference table.
c
c    Input, double precision XVAL, the X data value to be inserted as XTAB(1).
c
c    Input, double precision YVAL, the Y data value to be inserted as YTAB(1).
c
c    Output, integer NTAB2, the updated size of the difference table.
c
c    Output, double precision XTAB2(NTAB2), the updated abscissas.
c
c    Output, double precision DIFTAB2(NTAB2), the updated difference table.
c
      implicit none

      integer ntab

      double precision diftab(ntab)
      double precision diftab2(ntab+1)
      integer i
      integer ntab2
      double precision xtab(ntab)
      double precision xtab2(ntab+1)
      double precision xval
      double precision yval
c
c  Move the original data up one index.
c
      ntab2 = ntab + 1

      do i = ntab2, 2, -1
        xtab2(i) = xtab(i-1)
      end do

      do i = ntab2, 2, -1
        diftab2(i) = diftab(i-1)
      end do
c
c  Insert the new data.
c
      xtab2(1) = xval
      diftab2(1) = yval
c
c  Recompute the difference table.
c
      do i = 2, ntab2
        diftab2(i) = ( diftab2(i) - diftab2(i-1) )
     &             / ( xtab2(i)   - xtab2(1) )
      end do

      return
      end
      subroutine dif_basis ( ntab, xtab, diftab )

c*********************************************************************72
c
cc DIF_BASIS computes all Lagrange basis polynomials in divided difference form.
c
c  Discussion:
c
c    The I-th Lagrange basis polynomial for a set of NTAB X values XTAB,
c    L(I,NTAB,XTAB)(X) is a polynomial of order NTAB-1 which is zero at
c    XTAB(J) for J not equal to I, and 1 when J is equal to I.
c
c    The Lagrange basis polynomials have the property that the interpolating
c    polynomial through a set of NTAB data points (XTAB,YTAB) may be
c    represented as
c
c      P(X) = Sum ( 1 <= I <= N ) YTAB(I) * L(I,NTAB,XTAB)(X)
c
c    Higher order interpolation at selected points may be accomplished
c    using repeated X values, and scaled derivative values.
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
c  Reference:
c
c    Carl deBoor,
c    A Practical Guide to Splines,
c    Springer, 2001,
c    ISBN: 0387953663,
c    LC: QA1.A647.v27.
c
c  Parameters:
c
c    Input, integer NTAB, the number of X data points XTAB,
c    and the number of basis polynomials to compute.
c
c    Input, double precision XTAB(NTAB), the X values upon which the
c    Lagrange basis polynomials are to be based.
c
c    Output, double precision DIFTAB(NTAB,NTAB), the set of divided
c    difference tables.  Column I of DIFTAB contains the table for
c    the I-th Lagrange basis polynomial.
c
      implicit none

      integer ntab

      double precision diftab(ntab,ntab)
      integer i
      integer j
      double precision xtab(ntab)
c
c  Initialize DIFTAB to the identity matrix.
c
      do j = 1, ntab
        do i = 1, ntab
          diftab(i,j) = 0.0D+00
        end do
      end do

      do i = 1, ntab
        diftab(i,i) = 1.0D+00
      end do
c
c  Compute each Lagrange basis polynomial.
c
      do i = 1, ntab
        call data_to_dif ( ntab, xtab, diftab(1,i), diftab(1,i) )
      end do

      return
      end
      subroutine dif_basis_i ( ival, ntab, xtab, diftab )

c*********************************************************************72
c
cc DIF_BASIS_I: I-th Lagrange basis polynomial in divided difference form.
c
c  Discussion:
c
c    The I-th Lagrange basis polynomial for a set of NTAB X values XTAB,
c    L(I,NTAB,XTAB)(X) is a polynomial of order NTAB-1 which is zero at
c    XTAB(J) for J not equal to I, and 1 when J is equal to I.
c
c    The Lagrange basis polynomials have the property that the interpolating
c    polynomial through a set of NTAB data points (XTAB,YTAB) may be
c    represented as
c
c      P(X) = Sum ( 1 <= I <= N ) YTAB(I) * L(I,NTAB,XTAB)(X)
c
c    Higher order interpolation at selected points may be accomplished
c    using repeated X values, and scaled derivative values.
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
c  Reference:
c
c    Carl deBoor,
c    A Practical Guide to Splines,
c    Springer, 2001,
c    ISBN: 0387953663,
c    LC: QA1.A647.v27.
c
c  Parameters:
c
c    Input, integer IVAL, the index of the desired Lagrange
c    basis polynomial.  IVAL should be between 1 and NTAB.
c
c    Input, integer NTAB, the number of data points XTAB.
c
c    Input, double precision XTAB(NTAB), the X values upon which the
c    Lagrange basis polynomial is to be based.
c
c    Output, double precision DIFTAB(NTAB), the divided difference table
c    for the IVAL-th Lagrange basis polynomial.
c
      implicit none

      integer ntab

      double precision diftab(ntab)
      integer i
      integer ival
      double precision xtab(ntab)
c
c  Check IVAL.
c
      if ( ival .lt. 1 .or. ntab .lt. ival ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DIF_BASIS_I - Fatal errorc'
        write ( *, '(a,i6)' ) '  IVAL must be between 1 and ', ntab
        write ( *, '(a,i6)' ) '  but your value is ', ival
        stop
      end if
c
c  Initialize DIFTAB to Delta(I,J).
c
      do i = 1, ntab
        diftab(i) = 0.0D+00
      end do
      diftab(ival) = 1.0D+00
c
c  Compute the IVAL-th Lagrange basis polynomial.
c
      call data_to_dif ( ntab, xtab, diftab, diftab )

      return
      end
      subroutine dif_deriv_table ( nd, xd, yd, ndp, xdp, ydp )

c*********************************************************************72
c
cc DIF_DERIV_TABLE computes the divided difference table for a derivative.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Carl deBoor,
c    A Practical Guide to Splines,
c    Springer, 2001,
c    ISBN: 0387953663,
c    LC: QA1.A647.v27.
c
c  Parameters:
c
c    Input, integer ND, the size of the input table.
c
c    Input, double precision XD(ND), the abscissas for the divided
c    difference table.
c
c    Input, double precision YD(ND), the divided difference table.
c
c    Output, integer NDP, the size of the output table,
c    which is ND-1.
c
c    Input, double precision XDP(NDP), the abscissas for the divided
c    difference table for the derivative.
c
c    Output, double precision YDP(NDP), the divided difference
c    table for the derivative.
c
      implicit none

      integer nd

      integer i
      integer ndp
      double precision xd(nd)
      double precision xd_temp(nd)
      double precision xdp(nd-1)
      double precision yd(nd)
      double precision yd_temp(nd)
      double precision ydp(nd-1)
c
c  Using a temporary copy of the difference table, shift the
c  abscissas to zero.
c
      do i = 1, nd
        xd_temp(i) = xd(i)
        yd_temp(i) = yd(i)
      end do

      call dif_shift_zero ( nd, xd_temp, yd_temp )
c
c  Construct the derivative.
c
      ndp = nd - 1

      do i = 1, ndp
        xdp(i) = 0.0D+00
      end do

      do i = 1, ndp
        ydp(i) = dble ( i ) * yd_temp(i+1)
      end do

      return
      end
      subroutine dif_print ( ntab, xtab, diftab, title )

c*********************************************************************72
c
cc DIF_PRINT prints the polynomial represented by a divided difference table.
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
c  Parameters:
c
c    Input, integer NTAB, the dimension of the arrays DIFTAB and XTAB.
c
c    Input, double precision XTAB(NTAB), the X values for the polynomial.
c
c    Input, double precision DIFTAB(NTAB), the divided difference table
c    for the polynomial.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer ntab

      double precision diftab(ntab)
      integer i
      character * ( * ) title
      double precision xtab(ntab)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) title
      write ( *, '(a)' ) ' '
      write ( *, '( ''  p(x) =                           '', g14.6 )' )
     &  diftab(1)

      do i = 2, ntab
        write ( *, '( ''       + ( x - '', g14.6, '') * ( '', g14.6 )' )
     &   xtab(i-1), diftab(i)
      end do

      write ( *, '(80a1)' ) ( '       )', i = 1, ntab-1 )

      return
      end
      subroutine dif_root ( abserr, fxname, iprint, maxstp, maxtab,
     &  relerr, xroot, xtry1, xtry2 )

c*********************************************************************72
c
cc DIF_ROOT seeks a zero of F(X) using divided difference techniques.
c
c  Discussion:
c
c    The method uses the idea of considering the related function
c
c      H(X) = 1 / F(X)
c
c    The iteration begins with two estimates for the root supplied by
c    the user.
c
c    From the most recent approximation to the root, X(K), the next
c    approximation X(K+1) is determined by:
c
c      X(K+1) = X(K) + H(X(K-R),...,X(K-1)) / H(X(K-R),...,X(K-1),X(K))
c
c    where K-R = 1 until the maximal order NTAB is reached.
c
c    Generally, the next iterate X(K+1) is the zero of a rational function
c    which passes through the previous data points.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 November 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    FM Larkin,
c    Root Finding by Divided Differences,
c    Numerische Mathematik,
c    Volume 37, pages 93-104, 1981.
c
c  Parameters:
c
c    Input, double precision ABSERR, a positive absolute error tolerance.
c    If an estimate X for the root is found with ABS ( F(X) ) <= ABSERR,
c    the iteration is stopped.
c
c    Input, external FXNAME, the name of the function routine which evaluates
c    F(X).  The form of FXNAME must be similar to the following function which
c    has F(X) = ( X - 1 ) * ( X + 1 ).
c
c    function parab ( x )
c
c      double precision parab
c      double precision x
c
c      parab = ( x - 1.0D+00 ) * ( x + 1.0D+00 )
c
c      return
c    end
c
c    Input, integer IPRINT, a switch controlling printed output:
c    0, only print error messages.
c    nonzero, also print a table of the iterative process.
c
c    Input, integer MAXSTP, the limit on how many iterations
c    may be tried.
c
c    Input, integer MAXTAB, the limit on how high an order can be
c    used in the divided difference table.  MAXTAB must be at least 2, and
c    probably should not be too large.  Perhaps a value of 5 or 6 is reasonable,
c    20 is too large.
c
c    Input, double precision RELERR, a tolerance on the size of the change
c    in the root estimates.  If a step is taken, and the change in the root
c    estimate is less than RELERR, the iteration will stop.
c
c    Output, double precision XROOT, the point which the program has
c    produced as an approximate root.
c    Either ABS ( F(XROOT) ) <= ABSERR, or the maximum number of steps was
c    reached, or the current estimate of the root could not be significantly
c    improved.
c
c    Input, double precision XTRY1, XTRY2, two initial approximations to
c    the root, supplied by the user, which must be distinct.
c
      implicit none

      integer maxtab

      double precision abserr
      double precision diftab(maxtab)
      double precision froot
      double precision ftemp1
      double precision ftemp2
      double precision, external :: fxname
      integer iprint
      integer istep
      integer maxstp
      integer ntab
      double precision relerr
      double precision xdelt
      double precision xold
      double precision xroot
      double precision xtab(maxtab)
      double precision xtry1
      double precision xtry2
      double precision yval
c
c  Make sure XTRY1 and XTRY2 are not equal.
c
      if ( xtry1 .eq. xtry2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DIF_ROOT - Fatal errorc'
        write ( *, '(a)' ) '  XTRY1 = XTRY2 on input.'
        stop
      end if
c
c  Make sure MAXTAB is at least 2.
c
      if ( maxtab .lt. 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DIF_ROOT - Fatal errorc'
        write ( *, '(a)' ) '  MAXTAB < 2 on inputc'
        stop
      end if

      xtab(1) = xtry1
      xtab(2) = xtry2
      ftemp1 = fxname ( xtry1 )
      ftemp2 = fxname ( xtry2 )

      if ( abs ( ftemp2 ) .lt. abs ( ftemp1 ) ) then
        xtab(1) = xtry2
        xtab(2) = xtry1
        call r8_swap ( ftemp1, ftemp2 )
      end if
c
c  Initialize the number of steps.
c
      istep = 0
c
c  Initialize the number of data points.
c
      ntab = 2

      if ( 0 .lt. iprint ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' )
     &    '   Step  NTAB    XROOT        F(XROOT)      XDELT'
        write ( *, '(a)' ) ' '
      end if
c
c  Initialize the divided difference table data.
c
      diftab(1) = 1.0D+00 / ftemp1
      diftab(2) = 1.0D+00 / ftemp2

      call data_to_dif ( ntab, xtab, diftab, diftab )
c
c  Initialize values used in the iteration.
c
      xroot = xtry1
      froot = ftemp1
      xdelt = xtry1 - xtry2
c
c  Does the starting data already satisfy the function norm
c  error tolerance ABSERR, or the interval norm error tolerance
c  RELERR?
c
10    continue

        if ( 0 .lt. iprint ) then
          write ( *, '(3x,i4,4x,i2, 3g14.6)' ) 
     &      istep, ntab, xroot, froot, xdelt
        end if

        if ( abs ( froot ) <= abserr ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DIF_ROOT - Absolute convergence,'
          write ( *, '(a)' )
     &      '  The function value meets the error tolerance.'
          go to 20
        end if

        if ( abs ( xdelt ) .le. relerr ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DIF_ROOT - Relative convergence.'
          write ( *, '(a)' ) '  The stepsize meets the error tolerance.'
          go to 20
        end if
c
c  Check the number of steps taken.
c
        if ( maxstp .le. istep ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DIF_ROOT - Nonconvergence!'
          write ( *, '(a)' ) '  The maximum number of steps was taken.'
          go to 20
        end if
c
c  Generate the next point, XVAL.
c
        xold = xroot
        istep = istep + 1
c
c  Algorithm breakdown: The divisor DIFTAB(NTAB) is zero.
c
        if ( diftab(ntab) .eq. 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'DIF_ROOT - Fatal error!'
          write ( *, '(a,i6)' )
     &      '  Algorithm using differences of order ', ntab
          write ( *, '(a)' ) '  A zero-divisor was computed.'
          write ( *, '(a)' ) '  The algorithm has broken down.'
          write ( *, '(a)' )
     &      '  Examine the results.  They may be useful.'
          write ( *, '(a)' )
     &      '  Perhaps a lower value of MAXTAB would help.'
          stop
        end if

        xroot = xtab(ntab) + ( diftab(ntab-1) / diftab(ntab) )
        xdelt = xroot - xold
        froot = fxname ( xroot )

        if ( abs ( froot ) .le. abserr ) then
          go to 10
        end if

        yval = 1.0D+00 / froot
c
c  If we are now using MAXTAB points, we have to remove an old
c  one before adding the new one.
c
        if ( maxtab .le. ntab ) then
          ntab = ntab - 1
        end if

        call dif_append ( ntab, xtab, diftab, xroot, yval, ntab, xtab,
     &    diftab )

      go to 10

20    continue

      return
      end
      subroutine dif_shift_x ( nd, xd, yd, xv )

c*********************************************************************72
c
cc DIF_SHIFT_X replaces one abscissa of a divided difference table.
c
c  Discussion:
c
c    The routine shifts the representation of a divided difference polynomial by
c    dropping the last X value in XD, and adding a new X value to the
c    beginning of the XD array, suitably modifying the coefficients stored
c    in YD.
c
c    The representation of the polynomial is changed, but the polynomial itself
c    should be identical.
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
c  Reference:
c
c    Carl deBoor,
c    A Practical Guide to Splines,
c    Springer, 2001,
c    ISBN: 0387953663,
c    LC: QA1.A647.v27.
c
c  Parameters:
c
c    Input, integer ND, the number of divided difference
c    coefficients, and the number of entries in XD.
c
c    Input/output, double precision XD(ND), the X values used in
c    the representation of the divided difference polynomial.
c    After a call to this routine, the last entry of XD has been dropped,
c    the other entries have shifted up one index, and XV has been inserted
c    at the beginning of the array.
c
c    Input/output, double precision YD(ND), the divided difference
c    coefficients corresponding to the XD array.  On output, this
c    array has been adjusted.
c
c    Input, double precision XV, a new X value which is to be used in
c    the representation of the polynomial.  On output, XD(1) equals
c    XV and the representation of the polynomial has been suitably changed.
c    Note that XV does not have to be distinct from any of the original XD
c    values.
c
      implicit none

      integer nd

      integer i
      double precision xd(nd)
      double precision xv
      double precision yd(nd)
c
c  Recompute the divided difference coefficients.
c
      do i = nd - 1, 1, -1
        yd(i) = yd(i) + ( xv - xd(i) ) * yd(i+1)
      end do
c
c  Shift the XD values up one position and insert XV.
c
      do i = nd, 2, -1
        xd(i) = xd(i-1)
      end do

      xd(1) = xv

      return
      end
      subroutine dif_shift_zero ( nd, xd, yd )

c*********************************************************************72
c
cc DIF_SHIFT_ZERO shifts a divided difference table so all abscissas are zero.
c
c  Discussion:
c
c    When the abscissas are changed, the coefficients naturally
c    must also be changed.
c
c    The resulting pair (XD, YD) still represents the
c    same polynomial, but the entries in YD are now the
c    standard polynomial coefficients.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 November 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Carl deBoor,
c    A Practical Guide to Splines,
c    Springer, 2001,
c    ISBN: 0387953663,
c    LC: QA1.A647.v27.
c
c  Parameters:
c
c    Input, integer ND, the length of the table.
c
c    Input/output, double precision XD(ND), the abscissas for the divided
c    difference table.  On output, all entries have been reset to 0, but
c    (XD,YD) can still be regarded as a divided difference table for the input
c    polynomial.
c
c    Input/output, double precision YD(ND).  On input, the divided difference 
c    table for the polynomial.  On output, the divided difference table for
c    the polynomial, which has been rebased at 0.  Hence, YD is also simply
c    the coefficient array for the standard representation of the polynomial.
c
      implicit none

      integer nd

      integer i
      integer j
      double precision xd(nd)
      double precision yd(nd)

      do j = 1, nd

        do i = nd - 1, 1, -1
          yd(i) = yd(i) - xd(i) * yd(i+1)
        end do
c
c  Shift the XD values up one position.
c
        do i = nd, 2, -1
          xd(i) = xd(i-1)
        end do
        xd(1) = 0.0D+00

      end do

      return
      end
      subroutine dif_to_r8poly ( ntab, xtab, diftab, c )

c*********************************************************************72
c
cc DIF_TO_R8POLY converts a divided difference table to a standard polynomial.
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
c  Reference:
c
c    Carl deBoor,
c    A Practical Guide to Splines,
c    Springer, 2001,
c    ISBN: 0387953663,
c    LC: QA1.A647.v27.
c
c  Parameters:
c
c    Input, integer NTAB, the number of coefficients, and abscissas.
c
c    Input, double precision XTAB(NTAB), the X values used in the divided
c    difference representation of the polynomial.
c
c    Input, double precision DIFTAB(NTAB) the divided difference table.
c
c    Output, double precision C(NTAB), the polyomial coefficients.
c    C(1) is the constant term, and C(NTAB) is the coefficient
c    of X**(NTAB-1).
c
      implicit none

      integer ntab

      double precision c(ntab)
      double precision diftab(ntab)
      integer i
      integer j
      double precision xtab(ntab)

      do i = 1, ntab
        c(i) = diftab(i)
      end do
c
c  Recompute the divided difference coefficients.
c
      do j = 1, ntab - 1
        do i = 1, ntab - j
          c(ntab-i) = c(ntab-i) - xtab(ntab-i-j+1) * c(ntab-i+1)
        end do
      end do

      return
      end
      subroutine dif_val ( ntab, xtab, diftab, xv, yv )

c*********************************************************************72
c
cc DIF_VAL evaluates a divided difference polynomial at a point.
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
c  Reference:
c
c    Carl deBoor,
c    A Practical Guide to Splines,
c    Springer, 2001,
c    ISBN: 0387953663,
c    LC: QA1.A647.v27.
c
c  Parameters:
c
c    Input, integer NTAB, the number of divided difference
c    coefficients, and the number of points XTAB.
c
c    Input, double precision XTAB(NTAB), the X values upon which the
c    divided difference polynomial is based.
c
c    Input, double precision DIFTAB(NTAB), the divided difference
c    polynomial coefficients.
c
c    Input, double precision XV, the value where the polynomial
c    is to be evaluated.
c
c    Output, double precision YV, the value of the polynomial at XV.
c
      implicit none

      integer ntab

      double precision diftab(ntab)
      integer i
      double precision xtab(ntab)
      double precision xv
      double precision yv

      yv = diftab(ntab)
      do i = 1, ntab - 1
        yv = diftab(ntab-i) + ( xv - xtab(ntab-i) ) * yv
      end do

      return
      end
      subroutine dif_vals ( nd, xd, yd, nv, xv, yv )

c*********************************************************************72
c
cc DIF_VALS evaluates a divided difference polynomial at a set of points.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Carl deBoor,
c    A Practical Guide to Splines,
c    Springer, 2001,
c    ISBN: 0387953663,
c    LC: QA1.A647.v27.
c
c  Parameters:
c
c    Input, integer ND, the order of the difference table.
c
c    Input, double precision XD(ND), the X values of the difference table.
c
c    Input, double precision YD(ND), the divided differences.
c
c    Input, integer NV, the number of evaluation points.
c
c    Input, double precision XV(NV), the evaluation points.
c
c    Output, double precision YV(NV), the value of the divided difference
c    polynomial at the evaluation points.
c
      implicit none

      integer nd
      integer nv

      integer i
      integer j
      double precision xd(nd)
      double precision xv(nv)
      double precision yd(nd)
      double precision yv(nv)

      do j = 1, nv
        yv(j) = yd(nd)
        do i = 1, nd - 1
          yv(j) = yd(nd-i) + ( xv(j) - xd(nd-i) ) * yv(j)
        end do
      end do

      return
      end
      subroutine lagrange_rule ( ntab, xtab, wtab )

c*********************************************************************72
c
cc LAGRANGE_RULE computes the weights of a Lagrange interpolation rule.
c
c  Discussion:
c
c    Given NTAB abscissas XTAB, an arbitrary function F(X) can be
c    interpolated by a polynomial P(X) of order NTAB (and degree NTAB-1)
c    using weights that depend only on XTAB.
c
c    Standard Lagrange interpolation can be rewritten into this form,
c    which is more economical than evaluating the individual Lagrange
c    basis polynomials.
c
c    If we define
c
c      W(I) = 1 / product ( 1 <= J <= NTAB, J /= I ) ( XTAB(J) - XTAB(I) )
c
c    then
c
c      P(X) = sum ( 1 <= I <= NTAB ) W(I) * F( XTAB(I) ) / ( X - XTAB(I) )
c           / sum ( 1 <= I <= NTAB ) W(I)                / ( X - XTAB(I) )
c
c    except when X = XTAB(J), for some J, when we set:
c
c      P(XTAB(J)) = F(XTAB(J))
c
c  Modified:
c
c    21 February 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Jean-Paul Berrut, Lloyd Trefethen,
c    Barycentric Lagrange Interpolation,
c    SIAM Review,
c    Volume 46, Number 3, September 2004, pages 501-517.
c
c  Parameters:
c
c    Input, integer NTAB, the order of the rule.
c
c    Input, double precision XTAB(NTAB), the abscissas of the rule.
c
c    Output, double precision WTAB(NTAB), the weights of the rule.
c
      implicit none

      integer ntab

      integer i
      integer j
      double precision wtab(ntab)
      double precision xtab(ntab)

      do i = 1, ntab
        wtab(i) = 1.0D+00
      end do

      do i = 1, ntab
        do j = 1, i - 1
          wtab(j) = ( xtab(i) - xtab(j) ) * wtab(j)
        end do
        wtab(i) = 1.0D+00
        do j = 1, i - 1
          wtab(i) = wtab(i) * ( xtab(j) - xtab(i) )
        end do
      end do

      do i = 1, ntab
        wtab(i) = 1.0D+00 / wtab(i)
      end do

      return
      end
      subroutine lagrange_sum ( ntab, xtab, wtab, ftab, xval, fval )

c*********************************************************************72
c
cc LAGRANGE_SUM carries out a Lagrange interpolation rule.
c
c  Discussion:
c
c    It is assumed that LAGRANGE_RULE has already been called to compute
c    the appropriate weights for the given set of abscissas.
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
c  Reference:
c
c    Jean-Paul Berrut, Lloyd Trefethen,
c    Barycentric Lagrange Interpolation,
c    SIAM Review,
c    Volume 46, Number 3, September 2004, pages 501-517.
c
c  Parameters:
c
c    Input, integer NTAB, the order of the rule.
c
c    Input, double precision XTAB(NTAB), the abscissas of the rule.
c
c    Input, double precision WTAB(NTAB), the weights of the rule.
c
c    Input, double precision FTAB(NTAB), the function values at the abscissas.
c
c    Input, double precision XVAL, a point where an interpolated value is
c    needed.
c
c    Output, double precision FVAL, the interpolated function value.
c
      implicit none

      integer ntab

      double precision bot
      double precision ftab(ntab)
      double precision fval
      integer i
      double precision top
      double precision wtab(ntab)
      double precision xtab(ntab)
      double precision xval

      do i = 1, ntab

        if ( xval .eq. xtab(i) ) then
          fval = ftab(i)
          return
        end if

      end do

      top = 0.0D+00
      bot = 0.0D+00

      do i = 1, ntab
        top = top + wtab(i) * ftab(i) / ( xval - xtab(i) )
        bot = bot + wtab(i) / ( xval - xtab(i) )
      end do

      fval = top / bot

      return
      end
      subroutine lagrange_val ( ntab, xtab, ftab, xval, fval )

c*********************************************************************72
c
cc LAGRANGE_VAL applies a naive form of Lagrange interpolation.
c
c  Discussion:
c
c    Given NTAB abscissas XTAB, an arbitrary function F(X) can be
c    interpolated by a polynomial P(X) of order NTAB (and degree NTAB-1)
c    using Lagrange basis polynomials of degree NTAB-1.
c
c    Standard Lagrange interpolation can be rewritten into this form,
c    which is more economical than evaluating the individual Lagrange
c    basis polynomials.
c
c    If we define
c
c      L(I)(X) = product ( 1 <= J <= NTAB, J /= I )
c        ( X - XTAB(J) ) / ( XTAB(I) - XTAB(J) )
c
c    then
c
c      P(X) = sum ( 1 <= I <= NTAB ) F( XTAB(I) ) * L(I)(X)
c
c    Applying this form of the interpolation rule directly involves
c    about NTAB**2 work.  There are more efficient forms of the rule.
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
c  Parameters:
c
c    Input, integer NTAB, the number of data points.
c
c    Input, double precision XTAB(NTAB), the abscissas.
c
c    Input, double precision FTAB(NTAB), the function values at the abscissas.
c
c    Input, double precision XVAL, a point where an interpolated value is
c    needed.
c
c    Output, double precision FVAL, the interpolated function value.
c
      implicit none

      integer ntab

      double precision ftab(ntab)
      double precision fval
      integer i
      integer j
      double precision poly
      double precision xtab(ntab)
      double precision xval

      fval = 0.0D+00

      do i = 1, ntab
        poly = 1.0D+00
        do j = 1, ntab
          if ( j .ne. i ) then
            poly = poly * ( xval - xtab(j) ) / ( xtab(i) - xtab(j) )
          end if
        end do
        fval = fval + ftab(i) * poly
      end do

      return
      end
      subroutine nc_rule ( norder, a, b, xtab, weight )

c*********************************************************************72
c
cc NC_RULE computes the weights of a Newton-Cotes quadrature rule.
c
c  Discussion:
c
c    For the interval [A,B], the Newton-Cotes quadrature rule estimates
c
c      Integral ( A <= X <= B ) F(X) dX
c
c    using NORDER equally spaced abscissas XTAB(I) and a weight vector
c    WEIGHT(I):
c
c      Sum ( 1 <= I <= N ) WEIGHT(I) * F ( XTAB(I) ).
c
c    For the CLOSED rule, the abscissas include the points A and B.
c    For the OPEN rule, the abscissas do not include A and B.
c
c    For the common closed and open rules, the abscissas are equally spaced.
c
c    This routine allows the user to specify any set of abscissas;
c    hence, it can compute the standard open and closed rules, and
c    other variations.
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
c  Parameters:
c
c    Input, integer NORDER, the order of the rule.
c
c    Input, double precision A, B, the left and right endpoints of the interval
c    over which the quadrature rule is to be applied.
c
c    Input, double precision XTAB(NORDER), the abscissas of the rule.
c
c    Output, double precision WEIGHT(NORDER), the weights of the rule.
c
      implicit none

      integer norder

      double precision a
      double precision b
      double precision poly_cof(norder)
      integer i
      double precision weight(norder)
      double precision xtab(norder)
      double precision yvala
      double precision yvalb

      do i = 1, norder
c
c  Compute the Lagrange basis polynomial which is 1 at XTAB(I),
c  and zero at the other nodes.
c
        call r8poly_basis_1 ( i, norder, xtab, poly_cof )
c
c  Evaluate the antiderivative of the polynomial at the left and
c  right endpoints.
c
        call r8poly_ant_val ( norder, poly_cof, a, yvala )

        call r8poly_ant_val ( norder, poly_cof, b, yvalb )

        weight(i) = yvalb - yvala

      end do

      return
      end
      subroutine ncc_rule ( norder, xtab, weight )

c*********************************************************************72
c
cc NCC_RULE computes the coefficients of a Newton-Cotes closed quadrature rule.
c
c  Discussion:
c
c    For the interval [-1,1], the Newton-Cotes quadrature rule estimates
c
c      Integral ( -1 <= X <= 1 ) F(X) dX
c
c    using NORDER equally spaced abscissas XTAB(I) and a weight vector
c    WEIGHT(I):
c
c      Sum ( 1 <= I <= N ) WEIGHT(I) * F ( XTAB(I) ).
c
c    For the CLOSED rule, the abscissas include A and B.
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
c  Parameters:
c
c    Input, integer NORDER, the order of the rule.
c
c    Output, double precision XTAB(NORDER), the abscissas of the rule.
c
c    Output, double precision WEIGHT(NORDER), the weights of the rule.
c
      implicit none

      integer norder

      double precision a
      double precision b
      integer i
      double precision weight(norder)
      double precision xtab(norder)
c
c  Compute a closed quadrature rule.
c
      a = -1.0D+00
      b = 1.0D+00

      do i = 1, norder
        xtab(i) = ( dble ( norder - i     ) * a
     &            + dble (          i - 1 ) * b )
     &            / dble ( norder - 1     )
      end do

      call nc_rule ( norder, a, b, xtab, weight )

      return
      end
      subroutine nco_rule ( norder, xtab, weight )

c*********************************************************************72
c
cc NCO_RULE computes the coefficients of a Newton-Cotes open quadrature rule.
c
c  Discussion:
c
c    For the interval [-1,1], the Newton-Cotes quadrature rule estimates
c
c      Integral ( -1 <= X <= 1 ) F(X) dX
c
c    using NORDER equally spaced abscissas XTAB(I) and a weight vector
c    WEIGHT(I):
c
c      Sum ( 1 <= I <= N ) WEIGHT(I) * F ( XTAB(I) ).
c
c    For the OPEN rule, the abscissas do not include A and B.
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
c  Parameters:
c
c    Input, integer NORDER, the order of the rule.
c
c    Output, double precision XTAB(NORDER), the abscissas of the rule.
c
c    Output, double precision WEIGHT(NORDER), the weights of the  rule.
c
      implicit none

      integer norder

      double precision a
      double precision b
      integer i
      double precision weight(norder)
      double precision xtab(norder)

      a = -1.0D+00
      b = 1.0D+00

      do i = 1, norder
        xtab(i) = ( dble ( norder + 1 - i ) * a
     &            + dble (              i ) * b )
     &            / dble ( norder + 1     )
      end do

      call nc_rule ( norder, a, b, xtab, weight )

      return
      end
      subroutine r8poly_ant_cof ( n, poly_cof, poly_cof2 )

c*********************************************************************72
c
cc R8POLY_ANT_COF integrates a polynomial in standard form.
c
c  Discussion:
c
c    The antiderivative of a polynomial P(X) is any polynomial Q(X)
c    with the property that d/dX Q(X) = P(X).
c
c    This routine chooses the antiderivative whose constant term is zero.
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
c  Parameters:
c
c    Input, integer N, the order of the polynomial.
c
c    Input, double precision POLY_COF(N), the polynomial coefficients.
c    POLY_COF(1) is the constant term, and POLY_COF(N) is the
c    coefficient of X**(N-1).
c
c    Output, double precision POLY_COF2(N+1), the coefficients of
c    the antiderivative polynomial, in standard form.  The constant
c    term is set to zero.
c
      implicit none

      integer n

      double precision poly_cof(n)
      double precision poly_cof2(n+1)
      integer i
c
c  Set the constant term.
c
      poly_cof2(1) = 0.0D+00
c
c  Integrate the polynomial.
c
      do i = 2, n+1
        poly_cof2(i) = poly_cof(i-1) / dble ( i - 1 )
      end do

      return
      end
      subroutine r8poly_ant_val ( n, poly_cof, xval, yval )

c*********************************************************************72
c
cc R8POLY_ANT_VAL evaluates the antiderivative of a polynomial in standard form.
c
c  Discussion:
c
c    The constant term of the antiderivative is taken to be zero.
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
c  Parameters:
c
c    Input, integer N, the order of the polynomial.
c
c    Input, double precision POLY_COF(N), the polynomial coefficients.
c    POLY_COF(1) is the constant term, and POLY_COF(N) is the coefficient of
c    X**(N-1).
c
c    Input, double precision XVAL, the point where the antiderivative is to be
c    evaluated.
c
c    Output, double precision YVAL, the value of the antiderivative of
c    the polynomial at XVAL.
c
      implicit none

      integer n

      integer i
      double precision poly_cof(n)
      double precision xval
      double precision yval

      yval = 0.0D+00

      do i = n, 1, -1
        yval = ( yval + poly_cof(i) / dble ( i ) ) * xval
      end do

      return
      end
      subroutine r8poly_basis ( ntab, xtab, poly_cof )

c*********************************************************************72
c
cc R8POLY_BASIS computes all Lagrange basis polynomials in standard form.
c
c  Discussion:
c
c    The I-th Lagrange basis polynomial for a set of NTAB X values XTAB,
c    L(I,NTAB,XTAB)(X) is a polynomial of order NTAB-1 which is zero at
c    XTAB(J) for J not equal to I, and 1 when J is equal to I.
c
c    The Lagrange basis polynomials have the property that the interpolating
c    polynomial through a set of NTAB data points (XTAB,YTAB) may be
c    represented as
c
c      P(X) = Sum ( 1 <= I <= N ) YTAB(I) * L(I,NTAB,XTAB)(X)
c
c    Higher order interpolation at selected points may be accomplished
c    using repeated X values, and scaled derivative values.
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
c  Parameters:
c
c    Input, integer NTAB, the number of data points XTAB.
c
c    Input, double precision XTAB(NTAB), the X values upon which the
c    Lagrange basis polynomial is to be based.
c
c    Output, double precision POLY_COF(NTAB,NTAB), the polynomial
c    coefficients for the I-th Lagrange basis polynomial are stored
c    in column I.  POLY_COF(1,I) is the constant term, and POLY_COF(1,NTAB)
c    is the coefficient of X**(NTAB-1).
c
      implicit none

      integer ntab

      integer i
      integer j
      double precision poly_cof(ntab,ntab)
      double precision xtab(ntab)
c
c  Initialize POLY_COF to the identity matrix.
c
      do j = 1, ntab
        do i = 1, ntab
          poly_cof(i,j) = 0.0D+00
        end do
      end do

      do i = 1, ntab
        poly_cof(i,i) = 1.0D+00
      end do
c
c  Compute the divided difference table for the IVAL-th Lagrange basis
c  polynomial.
c
      do i = 1, ntab
        call data_to_dif ( ntab, xtab, poly_cof(1,i), poly_cof(1,i) )
      end do
c
c  Convert the divided difference table coefficients to standard polynomial
c  coefficients.
c
      do i = 1, ntab
        call dif_to_r8poly ( ntab, xtab, poly_cof(1,i), poly_cof(1,i) )
      end do

      return
      end
      subroutine r8poly_basis_1 ( ival, ntab, xtab, poly_cof )

c*********************************************************************72
c
cc R8POLY_BASIS_1 computes the I-th Lagrange basis polynomial in standard form.
c
c  Discussion:
c
c    The I-th Lagrange basis polynomial for a set of NTAB X values XTAB,
c    L(I,NTAB,XTAB)(X) is a polynomial of order NTAB-1 which is zero at
c    XTAB(J) for J not equal to I, and 1 when J is equal to I.
c
c    The Lagrange basis polynomials have the property that the interpolating
c    polynomial through a set of NTAB data points (XTAB,YTAB) may be
c    represented as
c
c      P(X) = Sum ( 1 <= I <= N ) YTAB(I) * L(I,NTAB,XTAB)(X)
c
c    Higher order interpolation at selected points may be accomplished
c    using repeated X values, and scaled derivative values.
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
c  Parameters:
c
c    Input, integer IVAL, the index of the desired Lagrange
c    basis polynomial.  IVAL should be between 1 and NTAB.
c
c    Input, integer NTAB, the number of data points XTAB.
c
c    Input, double precision XTAB(NTAB), the X values upon which the
c    Lagrange basis polynomial is to be based.
c
c    Output, double precision POLY_COF(NTAB), the polynomial
c    coefficients for the IVAL-th Lagrange basis polynomial.
c
      implicit none

      integer ntab

      integer i
      integer ival
      double precision poly_cof(ntab)
      double precision xtab(ntab)
c
c  Check IVAL.
c
      if ( ival .lt. 1 .or. ntab .lt. ival ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8POLY_BASE_1 - Fatal errorc'
        write ( *, '(a,i6)' ) '  IVAL must be between 1 and ', ntab
        write ( *, '(a,i6)' ) '  but your value is ', ival
        stop
      end if
c
c  Initialize POLY_COF to the IVAL-th column of the identity matrix.
c
      do i = 1, ntab
        poly_cof(i) = 0.0D+00
      end do
      poly_cof(ival) = 1.0D+00
c
c  Compute the divided difference table for the IVAL-th Lagrange basis
c  polynomial.
c
      call data_to_dif ( ntab, xtab, poly_cof, poly_cof )
c
c  Convert the divided difference table coefficients to standard polynomial
c  coefficients.
c
      call dif_to_r8poly ( ntab, xtab, poly_cof, poly_cof )

      return
      end
      subroutine r8poly_der_cof ( n, poly_cof, poly_cof2 )

c*********************************************************************72
c
cc R8POLY_DER_COF computes the coefficients of the derivative of a polynomial.
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
c  Parameters:
c
c    Input, integer N, the order of the polynomial.
c
c    Input, double precision POLY_COF(N), the coefficients of the
c    polynomial to be differentiated.  POLY_COF(1) is the constant term, and
c    POLY_COF(N) is the coefficient of X**(N-1).
c
c    Output, double precision POLY_COF2(N-1), the coefficients of the
c    derivative of the polynomial.
c
      implicit none

      integer n

      double precision poly_cof(n)
      double precision poly_cof2(n-1)
      integer i

      do i = 1, n-1
        poly_cof2(i) = dble ( i ) * poly_cof(i+1)
      end do

      return
      end
      subroutine r8poly_der_val ( n, poly_cof, xval, yval )

c*********************************************************************72
c
cc R8POLY_DER_VAL evaluates the derivative of a polynomial in standard form.
c
c  Discussion:
c
c    A polynomial in standard form, with coefficients POLY_COF(*),
c    may be written:
c
c      P(X) = POLY_COF(1)
c           + POLY_COF(2) * X
c           ...
c           + POLY_COF(N) * X**(N-1)
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
c  Parameters:
c
c    Input, integer N, the order of the polynomial.
c
c    Input, double precision POLY_COF(N), the polynomial coefficients.
c    POLY_COF(1) is the constant term, and POLY_COF(N) is the coefficient of
c    X**(N-1).
c
c    Input, double precision XVAL, a value where the derivative of the
c    polynomial is to be evaluated.
c
c    Output, double precision YVAL, the value of the derivative of the
c    polynomial at XVAL.
c
      implicit none

      integer n

      integer i
      double precision poly_cof(n)
      double precision xval
      double precision yval

      yval = dble ( n - 1 ) * poly_cof(n)
      do i = n-1, 2, -1
        yval = yval * xval + dble ( i - 1 ) * poly_cof(i)
      end do

      return
      end
      subroutine r8poly_order ( na, a, order )

c*********************************************************************72
c
cc R8POLY_ORDER returns the order of a polynomial.
c
c  Discussion:
c
c    The order of a polynomial is the degree plus 1.
c
c    The order of a constant polynomial is 1.
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
c  Parameters:
c
c    Input, integer NA, the dimension of A.
c
c    Input, double precision A(NA), the coefficients of the polynomials.
c
c    Output, integer ORDER, the order of A.
c
      implicit none

      integer na

      double precision a(na)
      integer order

      order = na

10    continue

      if ( 1 .lt. order ) then

        if ( a(order) .ne. 0.0D+00 ) then
          return
        end if

        order = order - 1

        go to 10

      end if

      return
      end
      subroutine r8poly_print ( n, a, title )

c*********************************************************************72
c
cc R8POLY_PRINT prints out a polynomial.
c
c  Discussion:
c
c    The power sum form is:
c
c      p(x) = a(0) + a(1)*x + ... + a(n-1)*x**(n-1) + a(n)*x**(n)
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
c  Parameters:
c
c    Input, integer N, the order of the polynomial.
c
c    Input, double precision A(N), the polynomial coefficients.
c    A(1) is the constant term and
c    A(N) is the coefficient of X**(N-1).
c
c    Input, character * ( * ) TITLE, an optional title.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      double precision mag
      integer n2
      character plus_minus
      character * ( * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '

      call r8poly_order ( n, a, n2 )

      if ( a(n2) .lt. 0.0D+00 ) then
        plus_minus = '-'
      else
        plus_minus = ' '
      end if

      mag = abs ( a(n2) )

      if ( 3 .le. n2 ) then
        write ( *, '( ''  p(x) = '', a1, g14.6, '' * x ^ '', i3 )' )
     &    plus_minus, mag, n2-1
      else if ( n2 .eq. 2 ) then
        write ( *, '( ''  p(x) = '', a1, g14.6, '' * x'' )' )
     &    plus_minus, mag
      else if ( n2 .eq. 1 ) then
        write ( *, '( ''  p(x) = '', a1, g14.6 )' ) plus_minus, mag
      end if

      do i = n2 - 1, 1, -1

        if ( a(i) .lt. 0.0D+00 ) then
          plus_minus = '-'
        else
          plus_minus = '+'
        end if

        mag = abs ( a(i) )

        if ( mag .ne. 0.0D+00 ) then

          if ( 3 .le. i ) then
            write ( *,
     &        ' ( ''         '', a1, g14.6, '' * x ^ '', i3 )' )
     &        plus_minus, mag, i-1
          else if ( i .eq. 2 ) then
            write ( *, ' ( ''         '', a1, g14.6, '' * x'' )' )
     &        plus_minus, mag
          else if ( i .eq. 1 ) then
            write ( *, ' ( ''         '', a1, g14.6 )' )
     &        plus_minus, mag
          end if
        end if

      end do

      return
      end
      subroutine r8poly_shift ( scale, shift, n, poly_cof )

c*********************************************************************72
c
cc R8POLY_SHIFT adjusts the coefficients of a polynomial for a new argument.
c
c  Discussion:
c
c    Assuming P(X) is a polynomial in the argument X, of the form:
c
c      P(X) = C(1)
c           + C(2) * X
c           + ...
c           + C(N) * X**(N-1)
c
c    and that Z is related to X by the formula:
c
c      Z = SCALE * X + SHIFT
c
c    then this routine computes coefficients C for the polynomial Q(Z):
c
c      Q(Z) = C(1)
c           + C(2) * Z
c           + ...
c           + C(N) * Z**(N-1)
c
c    so that:
c
c      Q(Z(X)) = P(X)
c
c  Example:
c
c    P(X) = 2 * X**2 - X + 6
c
c    Z = 2.0D+00 * X + 3.0D+00
c
c    Q(Z) = 0.5 *         Z**2 -  3.5 * Z + 12
c
c    Q(Z(X)) = 0.5 * ( 4.0D+00 * X**2 + 12.0D+00 * X +  9 )
c            - 3.5 * (               2.0D+00 * X +  3 )
c                                            + 12
c
c            = 2.0D+00         * X**2 -  1.0D+00 * X +  6
c
c            = P(X)
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
c  Parameters:
c
c    Input, double precision SHIFT, SCALE, the shift and scale applied to X,
c    so that Z = SCALE * X + SHIFT.
c
c    Input, integer N, the order of the polynomial
c
c    Input/output, double precision POLY_COF(N).
c    On input, the coefficient array in terms of the X variable.
c    On output, the coefficient array in terms of the Z variable.
c
      implicit none

      integer n

      integer i
      integer j
      double precision poly_cof(n)
      double precision scale
      double precision shift

      do i = 1, n
        do j = i + 1, n
          poly_cof(j) = poly_cof(j) / scale
        end do
      end do

      do i = 1, n
        do j = n-1, i, -1
          poly_cof(j) = poly_cof(j) - shift * poly_cof(j+1)
        end do
      end do

      return
      end
      subroutine r8poly_val_horner ( n, poly_cof, xval, yval )

c*********************************************************************72
c
cc R8POLY_VAL_HORNER evaluates a polynomial in standard form.
c
c  Discussion:
c
c    A polynomial in standard form, with coefficients POLY_COF(*),
c    may be written:
c
c      P(X) = C(1)
c           + C(2) * X
c           ...
c           + C(N) * X**(N-1)
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
c  Parameters:
c
c    Input, integer N, the order of the polynomial.
c
c    Input, double precision POLY_COF(N), the polynomial coefficients.
c    POLY_COF(1) is the constant term, and POLY_COF(N) is the coefficient of
c    X**(N-1).
c
c    Input, double precision XVAL, a value where the polynomial is to
c    be evaluated.
c
c    Output, double precision YVAL, the value of the polynomial at XVAL.
c
      implicit none

      integer n

      integer i
      double precision poly_cof(n)
      double precision xval
      double precision yval

      yval = 0.0D+00

      do i = n, 1, -1
        yval = yval * xval + poly_cof(i)
      end do

      return
      end
      subroutine r8_swap ( x, y )

c*********************************************************************72
c
cc R8_SWAP switches two R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 November 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, double precision X, Y.  On output, the values of X and
c    Y have been interchanged.
c
      implicit none

      double precision x
      double precision y
      double precision z

      z = x
      x = y
      y = z

      return
      end
      function r8vec_distinct ( n, a )

c*********************************************************************72
c
cc R8VEC_DISTINCT is true if the entries in an R8VEC are distinct.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, double precision A(N), the vector to be checked.
c
c    Output, logical R8VEC_DISTINCT is TRUE if the elements of A are distinct.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      integer j
      logical r8vec_distinct

      r8vec_distinct = .false.

      do i = 2, n
        do j = 1, i - 1
          if ( a(i) .eq. a(j) ) then
            return
          end if
        end do
      end do

      r8vec_distinct = .true.

      return
      end
      subroutine r8vec_even ( n, alo, ahi, a )

c*********************************************************************72
c
cc R8VEC_EVEN returns an R8VEC of evenly spaced values.
c
c  Discussion:
c
c    An R8VEC is a vector of R8 values.
c
c    If N is 1, then the midpoint is returned.
c
c    Otherwise, the two endpoints are returned, and N-2 evenly
c    spaced points between them.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 December 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of values.
c
c    Input, double precision ALO, AHI, the low and high values.
c
c    Output, double precision A(N), N evenly spaced values.
c    Normally, A(1) = ALO and A(N) = AHI.
c    However, if N = 1, then A(1) = 0.5*(ALO+AHI).
c
      implicit none

      integer n

      double precision a(n)
      double precision ahi
      double precision alo
      integer i

      if ( n .eq. 1 ) then

        a(1) = 0.5D+00 * ( alo + ahi )

      else

        do i = 1, n
          a(i) = ( dble ( n - i     ) * alo
     &           + dble (     i - 1 ) * ahi )
     &           / dble ( n     - 1 )
        end do

      end if

      return
      end
      subroutine r8vec_even_select ( n, xlo, xhi, ival, xval )

c*********************************************************************72
c
cc R8VEC_EVEN_SELECT returns the I-th of N evenly spaced values in [ XLO, XHI ].
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    XVAL = ( (N-IVAL) * XLO + (IVAL-1) * XHI ) / real ( N - 1 )
c
c    Unless N = 1, X(1) = XLO and X(N) = XHI.
c
c    If N = 1, then X(1) = 0.5*(XLO+XHI).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of values.
c
c    Input, double precision XLO, XHI, the low and high values.
c
c    Input, integer IVAL, the index of the desired point.
c    IVAL is normally between 1 and N, but may be any integer value.
c
c    Output, double precision XVAL, the IVAL-th of N evenly spaced values
c    between XLO and XHI.
c
      implicit none

      integer n

      integer ival
      double precision xhi
      double precision xlo
      double precision xval

      if ( n .eq. 1 ) then

        xval = 0.5D+00 * ( xlo + xhi )

      else

        xval = ( dble ( n - ival     ) * xlo
     &         + dble (     ival - 1 ) * xhi )
     &         / dble ( n        - 1 )

      end if

      return
      end
      subroutine r8vec_indicator ( n, a )

c*********************************************************************72
c
cc R8VEC_INDICATOR sets an R8VEC to the indicator vector.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of elements of A.
c
c    Output, double precision A(N), the array to be initialized.
c
      implicit none

      integer n

      double precision a(n)
      integer i

      do i = 1, n
        a(i) = dble ( i )
      end do

      return
      end
      subroutine r8vec_print ( n, a, title )

c*********************************************************************72
c
cc R8VEC_PRINT prints an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
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
c    Input, integer N, the number of components of the vector.
c
c    Input, double precision A(N), the vector to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      character ( len = * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
      end do

      return
      end
      subroutine roots_to_dif ( nroots, roots, ntab, xtab, diftab )

c*********************************************************************72
c
cc ROOTS_TO_DIF sets a divided difference table for a polynomial from its roots.
c
c  Discussion:
c
c    This turns out to be a simple task, because of two facts:
c
c    * The divided difference polynomial of one smaller degree which
c      passes through the values ( ROOT(I), 0 ) is the zero polynomial,
c      and hence has a zero divided difference table.
c
c    * We want a polynomial of one degree higher, but we don't want it
c      to pass through an addditional point.  Instead, we specify that
c      the polynomial is MONIC.  This means that the divided difference
c      table is almost the same as for the zero polynomial, except that
c      there is one more pair of entries, an arbitrary X value, and
c      a Y value of 1.
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
c  Parameters:
c
c    Input, integer NROOTS, is the number of roots.
c
c    Input, double precision ROOTS(NROOTS), the roots of
c    the polynomial.
c
c    Output, integer NTAB, is equal to NROOTS+1.
c
c    Output, double precision XTAB(NTAB), the abscissas of the divided
c    difference table.
c
c    Output, double precision DIFTAB(NTAB), the divided difference
c    table.
c
      implicit none

      integer nroots

      double precision diftab(nroots+1)
      integer ntab
      double precision roots(nroots)
      double precision xtab(nroots+1)

      ntab = nroots + 1
c
c  Build the appropriate difference table for the polynomial
c  through ( ROOTS(I), 0 ) of order NTAB-1.
c
      diftab(1:ntab-1) = 0.0D+00
c
c  Append the extra data to make a monic polynomial of order NTAB
c  which is zero at the NTAB-1 roots.
c
      xtab(1:ntab-1) = roots(1:ntab-1)
      xtab(ntab) = 0.0D+00

      diftab(ntab) = 1.0D+00

      return
      end
      subroutine roots_to_r8poly ( nroots, roots, nc, c )

c*********************************************************************72
c
cc ROOTS_TO_R8POLY converts polynomial roots to polynomial coefficients.
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
c  Parameters:
c
c    Input, integer NROOTS, the number of roots specified.
c
c    Input, double precision ROOTS(NROOTS), the roots.
c
c    Output, integer NC, the order of the polynomial, which will
c    be NROOTS + 1.
c
c    Output, double precision C(NC), the coefficients of the polynomial.
c
      implicit none

      integer nroots

      double precision c(nroots+1)
      integer i
      integer nc
      double precision roots(nroots)
      double precision xtab(nroots+1)

      nc = nroots + 1
c
c  Initialize C to (0, 0, ..., 0, 1).
c  Essentially, we are setting up a divided difference table.
c
      do i = 1, nroots
        xtab(i) = roots(i)
      end do
      xtab(nc) = 0.0

      do i = 1, nc - 1
        c(i) = 0.0D+00
      end do
      c(nc) = 1.0D+00
c
c  Convert to standard polynomial form by shifting the abscissas
c  of the divided difference table to 0.
c
      call dif_shift_zero ( nc, xtab, c )

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
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
