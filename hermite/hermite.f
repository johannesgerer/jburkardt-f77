      subroutine dif_deriv ( nd, xd, yd, ndp, xdp, ydp )

c*********************************************************************72
c
cc DIF_DERIV computes the derivative of a polynomial in divided difference form.
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
c    Input, integer ND, the length of the XD and YD arrays.
c
c    Input/output, double precision XD(ND), the X values that
c    correspond to the divided difference table.  On output, XD
c    contains only zeroes.
c
c    Input/output, double precision YD(ND), the divided difference table
c    for the polynomial.  On output, YD is also the coefficient array for 
c    the standard representation of the polynomial.
c
      implicit none

      integer nd

      integer i
      double precision xd(nd)
      double precision xv
      double precision yd(nd)

      xv = 0.0D+00

      do i = 1, nd
        call dif_shift_x ( nd, xd, yd, xv )
      end do

      return
      end
      subroutine dif_to_r8poly ( nd, xd, yd, c )

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
c    05 September 2004
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
c    Input, integer ND, the number of coefficients, and abscissas.
c
c    Input, double precision XD(ND), the X values used in the divided
c    difference representation of the polynomial.
c
c    Input, double precision YD(ND) the divided difference table.
c
c    Output, double precision C(ND), the polyomial coefficients.
c    C(1) is the constant term, and C(ND) is the coefficient
c    of X**(ND-1).
c
      implicit none

      integer nd

      double precision c(nd)
      integer i
      integer j
      double precision xd(nd)
      double precision yd(nd)

      do i = 1, nd
        c(i) = yd(i)
      end do
c
c  Recompute the divided difference coefficients.
c
      do j = 1, nd - 1
        do i = 1, nd - j
          c(nd-i) = c(nd-i) - xd(nd-i-j+1) * c(nd-i+1)
        end do
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
      subroutine hermite_basis_0 ( n, x, i, xv, value )

c*********************************************************************72
c
cc HERMITE_BASIS_0 evaluates a zero-order Hermite interpolation basis function.
c
c  Discussion:
c
c    Given ND points XD, with values YD and derivative values YPD, the
c    Hermite interpolant can be written as:
c
c      H(X) = sum ( 1 <= I <= ND ) YD(I)  * H0(I;X)
c           + sum ( 1 <= I <= ND ) YPD(I) * H1(I;X)
c
c    where H0(I;X) is the I-th zero order Hermite interpolation basis function,
c    and H1(I;X) is the I-th first order Hermite interpolation basis function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of abscissas.
c
c    Input, double precision X(N), the abscissas.
c
c    Input, integer I, the index of the first-order basis function.
c
c    Input, double precision XV, the evaluation point.
c
c    Output, double precision VALUE, the value of the function.
c
      implicit none

      integer n

      double precision factor(n)
      integer i
      integer j
      double precision li
      double precision lp
      double precision lpp
      double precision r8vec_product
      double precision value
      double precision x(n)
      double precision xv

      if ( i .lt. 1 .or. n .lt. i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'HERMITE_BASIS_0 - Fatal error!'
        write ( *, '(a)' ) '  I < 1 or N < I.'
        stop
      end if
c
c  L(X) = product ( X - X(1:N) )
c
c
c  L'(X(I)).
c
      do j = 1, n
        factor(j) = x(i) - x(j)
      end do
      factor(i) = 1.0D+00

      lp = r8vec_product ( n, factor )
c
c  LI(X) = L(X) / ( X - X(I) ) / L'(X(I))
c
      do j = 1, n
        factor(j) = xv - x(j)
      end do
      factor(i) = 1.0D+00

      li = r8vec_product ( n, factor ) / lp
c
c  L''(X(I)).
c
      lpp = 0.0D+00
      do j = 1, n
        factor(j) = x(i) - x(j)
      end do
      factor(i) = 1.0D+00

      do j = 1, n
        if ( j .ne. i ) then
          factor(j) = 1.0D+00

          lpp = lpp + 2.0D+00 * r8vec_product ( n, factor )
          factor(j) = x(i) - x(j)
        end if
      end do

      value = ( 1.0D+00 - ( xv - x(i) ) * lpp / lp ) * li * li

      return
      end
      subroutine hermite_basis_1 ( n, x, i, xv, value )

c*********************************************************************72
c
cc HERMITE_BASIS_1 evaluates a first-order Hermite interpolation basis function.
c
c  Discussion:
c
c    Given ND points XD, with values YD and derivative values YPD, the
c    Hermite interpolant can be written as:
c
c      H(X) = sum ( 1 <= I <= ND ) YD(I)  * H0(I;X)
c           + sum ( 1 <= I <= ND ) YPD(I) * H1(I;X)
c
c    where H0(I;X) is the I-th zero order Hermite interpolation basis function,
c    and H1(I;X) is the I-th first order Hermite interpolation basis function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 May 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of abscissas.
c
c    Input, double precision X(N), the abscissas.
c
c    Input, integer I, the index of the first-order basis function.
c
c    Input, double precision XV, the evaluation point.
c
c    Output, double precision VALUE, the value of the function.
c
      implicit none

      integer n

      double precision bot
      double precision factor(n)
      integer i
      integer j
      double precision r8vec_product
      double precision top
      double precision value
      double precision x(n)
      double precision xv

      if ( i .lt. 1 .or. n .lt. i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'HERMITE_BASIS_1 - Fatal error!'
        write ( *, '(a)' ) '  I < 1 or N < I.'
        stop
      end if

      do j = 1, n
        factor(j) = xv - x(j)
      end do
      factor(i) = 1.0D+00
      top = r8vec_product ( n, factor )

      do j = 1, n
        factor(j) = x(i) - x(j)
      end do
      factor(i) = 1.0D+00
      bot = r8vec_product ( n, factor )

      value = ( xv - x(i) ) * top**2 / bot**2

      return
      end
      subroutine hermite_demo ( n, x, y, yp )

c*********************************************************************72
c
cc HERMITE_DEMO computes and prints Hermite interpolant information for data.
c
c  Discussion:
c
c    Given a set of Hermite data, this routine calls HDATA_TO_DIF to determine
c    and print the divided difference table, and then DIF_TO_R8POLY to 
c    determine and print the coefficients of the polynomial in standard form.
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
c  Parameters:
c
c    Input, integer N, the number of data points.
c
c    Input, double precision X(N), the abscissas.
c
c    Input, double precision Y(N), YP(N), the function and derivative
c    values at the abscissas.
c
      implicit none

      integer n

      double precision cd(0:2*n-1)
      integer i
      integer nd
      integer ndp
      integer nv
      double precision x(n)
      double precision xd(2*n)
      double precision xdp(2*n-1)
      double precision xv(n)
      double precision y(n)
      double precision yd(2*n)
      double precision ydp(2*n-1)
      double precision yp(n)
      double precision yv(n)
      double precision yvp(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HERMITE_DEMO'
      write ( *, '(a)' ) 
     &  '  Compute coefficients CD of the Hermite polynomial'
      write ( *, '(a)' ) '  interpolant to given data (x,y,yp).'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Data:'
      write ( *, '(a)' ) '              X           Y           Y'''
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4)' ) 
     &    i, x(i), y(i), yp(i)
      end do

      call hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Difference table for interpolant:'
      write ( *, '(a)' ) '              XD          YD'
      write ( *, '(a)' ) ' '

      nd = 2 * n

      do i = 1, nd
        write ( *, '(2x,i4,2x,f10.4,2x,f10.4)' ) i, xd(i), yd(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Difference table for interpolant derivative:'
      write ( *, '(a)' ) '              XDP         YDP'
      write ( *, '(a)' ) ' '

      ndp = 2 * n - 1

      do i = 1, ndp
        write ( *, '(2x,i4,2x,f10.4,2x,f10.4)' ) i, xdp(i), ydp(i)
      end do

      call dif_to_r8poly ( nd, xd, yd, cd )

      call r8poly_print ( nd - 1, cd, 
     &  '  Hermite interpolating polynomial:' )
c
c  Verify interpolation claim!
c
      nv = n
      do i = 1, n
        xv(i) = x(i)
      end do

      call hermite_interpolant_value ( nd, xd, yd, xdp, ydp, nv, xv, 
     &  yv, yvp )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Data Versus Interpolant:'
      write ( *, '(a)' ) '              X           Y           ' //
     &  'H           YP          HP'
      write ( *, '(a)' ) ' '
      do i = 1, nv
        write ( *, 
     &  '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4)' ) 
     &    i, xv(i), y(i), yv(i), yp(i), yvp(i)
      end do

      return
      end
      subroutine hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp )

c*********************************************************************72
c
cc HERMITE_INTERPOLANT sets up a divided difference table from Hermite data.
c
c  Discussion:
c
c    The polynomial represented by the divided difference table can be
c    evaluated by calling DIF_VALS.
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
c    Input, integer N, the number of items of data 
c    ( X(I), Y(I), YP(I) ).
c
c    Input, double precision X(N), the abscissas.
c    These values must be distinct.
c
c    Input, double precision Y(N), YP(N), the function and 
c    derivative values.
c
c    Output, double precision XD(2*N), YD(2*N), the divided difference table
c    for the interpolant value.
c
c    Output, double precision XDP(2*N-1), YDP(2*N-1), the divided difference 
c    table for the interpolant derivative.
c
      implicit none

      integer n

      integer i
      integer j
      integer nd
      integer ndp
      double precision x(n)
      double precision xd(2*n)
      double precision xdp(2*n-1)
      double precision y(n)
      double precision yd(2*n)
      double precision ydp(2*n-1)
      double precision yp(n)
c
c  Copy the data:
c
      nd = 2 * n
      do i = 1, n
        xd(2*i-1) = x(i)
        xd(2*i  ) = x(i)
      end do
c
c  Carry out the first step of differencing.
c
      yd(1) = y(1)
      do i = 2, n
        yd(2*i-1) = ( y(i) - y(i-1) ) / ( x(i) - x(i-1) )
      end do
      do i = 1, n
        yd(2*i) = yp(i)  
      end do
c
c  Carry out the remaining steps in the usual way.
c
      do i = 3, nd
        do j = nd, i, -1

          yd(j) = ( yd(j) - yd(j-1) ) / ( xd(j) - xd(j+1-i) )

        end do

      end do
c
c  Compute the difference table for the derivative.
c
      call dif_deriv ( nd, xd, yd, ndp, xdp, ydp )

      return
      end
      subroutine hermite_interpolant_rule ( n, a, b, x, w )

c*********************************************************************72
c
cc HERMITE_INTERPOLANT_RULE: quadrature rule for a Hermite interpolant.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of abscissas.
c
c    Input, double precision A, B, the integration limits.
c
c    Input, double precision X(N), the abscissas.
c
c    Output, double precision W(2*N), the quadrature coefficients, given as
c    pairs for function and derivative values at each abscissa.
c
      implicit none

      integer n

      double precision a
      double precision a_value
      double precision b
      double precision b_value
      double precision c(0:2*n-1)
      integer i
      integer j
      integer k
      integer nd
      double precision w(2*n)
      double precision x(n)
      double precision xd(2*n)
      double precision xdp(2*n-1)
      double precision y(n)
      double precision yd(2*n)
      double precision ydp(2*n-1)
      double precision yp(n)

      nd = 2 * n

      k = 0

      do i = 1, n

        k = k + 1
        do j = 1, n
          y(j) = 0.0D+00
        end do
        y(i) = 1.0D+00
        do j = 1, n
          yp(j) = 0.0D+00
        end do
        call hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp )
        call dif_to_r8poly ( nd, xd, yd, c )
        call r8poly_ant_val ( n, c, a, a_value )
        call r8poly_ant_val ( n, c, b, b_value )
        w(k) = b_value - a_value

        k = k + 1
        do j = 1, n
          y(j) = 0.0D+00
        end do
        do j = 1, n
          yp(j) = 0.0D+00
        end do
        yp(i) = 1.0D+00
        call hermite_interpolant ( n, x, y, yp, xd, yd, xdp, ydp )
        call dif_to_r8poly ( nd, xd, yd, c )
        call r8poly_ant_val ( n, c, a, a_value )
        call r8poly_ant_val ( n, c, b, b_value )
        w(k) = b_value - a_value

      end do

      return
      end
      subroutine hermite_interpolant_value ( nd, xd, yd, xdp, ydp, nv, 
     &  xv, yv, yvp )

c*********************************************************************72
c
cc HERMITE_INTERPOLANT_VALUE evaluates the Hermite interpolant polynomial.
c
c  Discussion:
c
c    In fact, this function will evaluate an arbitrary polynomial that is
c    represented by a difference table.
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
c    Input, integer ND, the order of the difference table.
c
c    Input, double precision XD(ND), YD(ND), the difference table for the
c    interpolant value.
c
c    Input, double precision XDP(ND-1), YDP(ND-1), the difference table for
c    the interpolant derivative.
c
c    Input, integer NV, the number of evaluation points.
c
c    Input, double precision XV(NV), the evaluation points.
c
c    Output, double precision YV(NV), YVP(NV), the value of the interpolant and
c    its derivative at the evaluation points.
c
      implicit none

      integer nd
      integer nv

      integer i
      integer j
      integer ndp
      double precision xd(nd)
      double precision xdp(nd-1)
      double precision xv(nv)
      double precision yd(nd)
      double precision ydp(nd-1)
      double precision yv(nv)
      double precision yvp(nv)

      ndp = nd - 1

      do j = 1, nv

        yv(j) = yd(nd)
        do i = nd - 1, 1, -1
          yv(j) = yd(i) + ( xv(j) - xd(i) ) * yv(j)
        end do

        yvp(j) = ydp(ndp)
        do i = ndp - 1, 1, -1
          yvp(j) = ydp(i) + ( xv(j) - xdp(i) ) * yvp(j)
        end do

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
c    X^(N-1).
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
      subroutine r8poly_degree ( na, a, degree )

c*********************************************************************72
c
cc R8POLY_DEGREE returns the degree of a polynomial.
c
c  Discussion:
c
c    The degree of a polynomial is the index of the highest power
c    of X with a nonzero coefficient.
c
c    The degree of a constant polynomial is 0.  The degree of the
c    zero polynomial is debatable, but this routine returns the
c    degree as 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 March 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NA, the dimension of A.
c
c    Input, double precision A(0:NA), the coefficients of the polynomials.
c
c    Output, integer DEGREE, the degree of A.
c
      implicit none

      integer na

      double precision a(0:na)
      integer degree

      degree = na

10    continue

      if ( 0 .lt. degree ) then

        if ( a(degree) .ne. 0.0D+00 ) then
          return
        end if

        degree = degree - 1

        go to 10

      end if

      return
      end
      subroutine r8poly_print ( n, a, title )

c*********************************************************************72
c
cc R8POLY_PRINT prints out a polynomial.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of A.
c
c    Input, double precision A(0:N), the polynomial coefficients.
c    A(0) is the constant term and
c    A(N) is the coefficient of X**N.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n

      double precision a(0:n)
      integer i
      double precision mag
      integer n2
      character plus_minus
      character * ( * )  title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '

      call r8poly_degree ( n, a, n2 )

      if ( n2 .le. 0 ) then
        write ( *, '( ''  p(x) = 0'' )' )
        return
      end if

      if ( a(n2) .lt. 0.0D+00 ) then
        plus_minus = '-'
      else
        plus_minus = ' '
      end if

      mag = abs ( a(n2) )

      if ( 2 .le. n2 ) then
        write ( *, '( ''  p(x) = '', a1, g14.6, '' * x ^ '', i3 )' )
     &    plus_minus, mag, n2
      else if ( n2 .eq. 1 ) then
        write ( *, '( ''  p(x) = '', a1, g14.6, '' * x'' )' )
     &    plus_minus, mag
      else if ( n2 .eq. 0 ) then
        write ( *, '( ''  p(x) = '', a1, g14.6 )' ) plus_minus, mag
      end if

      do i = n2 - 1, 0, -1

        if ( a(i) .lt. 0.0D+00 ) then
          plus_minus = '-'
        else
          plus_minus = '+'
        end if

        mag = abs ( a(i) )

        if ( mag .ne. 0.0D+00 ) then

          if ( 2 .le. i ) then
            write ( *,
     &        ' ( ''         '', a1, g14.6, '' * x ^ '', i3 )' )
     &        plus_minus, mag, i
          else if ( i .eq. 1 ) then
            write ( *,
     &        ' ( ''         '', a1, g14.6, '' * x'' )' )
     &        plus_minus, mag
          else if ( i .eq. 0 ) then
            write ( *, ' ( ''         '', a1, g14.6 )' )
     &        plus_minus, mag
          end if
        end if

      end do

      return
      end
      subroutine r8vec_chebyshev ( n, a_first, a_last, a )

c*********************************************************************72
c
cc R8VEC_CHEBYSHEV creates a vector of Chebyshev spaced values.
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
c    08 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, double precision A_FIRST, A_LAST, the first and last entries.
c
c    Output, double precision A(N), a vector of Chebyshev spaced data.
c
      implicit none

      integer n

      double precision a(n)
      double precision a_first
      double precision a_last
      double precision c
      integer i
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision theta

      if ( n .eq. 1 ) then

        a(1) = ( a_first + a_last ) / 2.0D+00

      else

        do i = 1, n

          theta = dble ( n - i ) * pi / dble ( n - 1 )

          c = cos ( theta )

          if ( mod ( n, 2 ) .eq. 1 ) then
            if ( 2 * i - 1 .eq. n ) then
              c = 0.0D+00
            end if
          end if

          a(i) = ( ( 1.0D+00 - c ) * a_first  
     &           + ( 1.0D+00 + c ) * a_last ) 
     &           /   2.0D+00

        end do

      end if

      return
      end
      subroutine r8vec_linspace ( n, a_first, a_last, a )

c*********************************************************************72
c
cc R8VEC_LINSPACE creates a vector of linearly spaced values.
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
c    14 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, double precision A_FIRST, A_LAST, the first and last entries.
c
c    Output, double precision A(N), a vector of linearly spaced data.
c
      implicit none

      integer n

      double precision a(n)
      double precision a_first
      double precision a_last
      integer i

      if ( n .eq. 1 ) then

        a(1) = ( a_first + a_last ) / 2.0D+00

      else

        do i = 1, n
          a(i) = ( dble ( n - i     ) * a_first 
     &           + dble (     i - 1 ) * a_last )
     &           / dble ( n     - 1 )
        end do

      end if

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
      function r8vec_product ( n, v1 )

c*********************************************************************72
c
cc R8VEC_PRODUCT multiplies the entries of an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8 values.
c
c    In FORTRAN90, the system routine PRODUCT should be called
c    directly.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the vectors.
c
c    Input, double precision V1(N), the vector.
c
c    Output, double precision R8VEC_PRODUCT, the product of the entries.
c
      implicit none

      integer n

      integer i
      double precision r8vec_product
      double precision v1(n)
      double precision value

      value = 1.0D+00
      do i = 1, n
        value = value * v1(i)
      end do

      r8vec_product = value

      return
      end
      subroutine r8vec_uniform_01 ( n, seed, r )

c*********************************************************************72
c
cc R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
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
c    17 July 2006
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
c    Input, integer N, the number of entries in the vector.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R(N), the vector of pseudorandom values.
c
      implicit none

      integer n

      integer i
      integer k
      integer seed
      double precision r(n)

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + 2147483647
        end if

        r(i) = dble ( seed ) * 4.656612875D-10

      end do

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
