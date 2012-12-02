      subroutine cc_compute_points ( n, points )

c*********************************************************************72
c
cc CC_COMPUTE_POINTS: abscissas of a Clenshaw Curtis rule.
c
c  Discussion:
c
c    Our convention is that the abscissas are numbered from left to right.
c
c    The rule is defined on [-1,1].
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order.
c    1 <= N.
c
c    Output, double precision POINTS(N), the abscissas.
c
      implicit none

      integer n

      integer i
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision points(n)

      if ( n .lt. 1 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CC_COMPUTE_POINTS - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal value of N = ', n
        stop

      else if ( n .eq. 1 ) then

        points(1) = 0.0D+00

      else

        do i = 1, n
          points(i) = cos ( dble ( n - i ) * pi 
     &                    / dble ( n - 1 ) )
        end do

        points(1) = -1.0D+00
        if ( mod ( n, 2 ) .eq. 1 ) then
          points((n+1)/2) = 0.0D+00
        end if
        points(n) = +1.0D+00

      end if

      return
      end
      subroutine i4vec_print ( n, a, title )

c*********************************************************************72
c
cc I4VEC_PRINT prints an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of integer values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components of the vector.
c
c    Input, integer A(N), the vector to be printed.
c
c    Input, character*(*) TITLE, a title.
c
      implicit none

      integer n

      integer a(n)
      integer i
      character*(*) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,a,1x,i12)' ) i, ':', a(i)
      end do

      return
      end
      function i4vec_product ( n, a )

c*********************************************************************72
c
cc I4VEC_PRODUCT returns the product of the entries of an I4VEC.
c
c  Discussion:
c
c    An I4VEC is a vector of I4's.
c
c    In FORTRAN90, this facility is offered by the built in
c    PRODUCT function:
c
c      I4VEC_PRODUCT ( N, A ) = PRODUCT ( A(1:N) )
c
c    In MATLAB, this facility is offered by the built in
c    PROD function:
c
c      I4VEC_PRODUCT ( N, A ) = PROD ( A(1:N) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input, integer A(N), the array.
c
c    Output, integer I4VEC_PRODUCT, the product of the entries.
c
      implicit none

      integer n

      integer a(n)
      integer i
      integer i4vec_product

      i4vec_product = 1
      do i = 1, n
        i4vec_product = i4vec_product * a(i)
      end do

      return
      end
      subroutine lagrange_base_1d ( nd, xd, ni, xi, lb )

c*********************************************************************72
c
cc LAGRANGE_BASE_1D evaluates the Lagrange basis polynomials.
c
c  Discussion:
c
c    Given ND distinct abscissas, XD(1:ND),
c    the I-th Lagrange basis polynomial LB(I)(T) is defined as the polynomial of
c    degree ND - 1 which is 1 at XD(I) and 0 at the ND - 1
c    other abscissas.
c
c    A formal representation is:
c
c      LB(I)(T) = Product ( 1 <= J <= ND, I /= J )
c       ( T - T(J) ) / ( T(I) - T(J) )
c
c    This routine accepts a set of NI values at which all the Lagrange
c    basis polynomials should be evaluated.
c
c    Given data values YD at each of the abscissas, the value of the
c    Lagrange interpolating polynomial at each of the interpolation points
c    is then simple to compute by matrix multiplication:
c
c      YI(1:NI) = LB(1:NI,1:ND) * YD(1:ND)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer ND, the number of data points.
c    ND must be at least 1.
c
c    Input, double precision XD(ND), the data points.
c
c    Input, integer NI, the number of interpolation points.
c
c    Input, double precision XI(NI), the interpolation points.
c
c    Output, double precision LB(NI,ND), the values
c    of the Lagrange basis polynomials at the interpolation points.
c
      implicit none

      integer nd
      integer ni

      integer i
      integer j
      integer k
      double precision lb(ni,nd)
      double precision xd(nd)
      double precision xi(ni)
c
c  Evaluate the polynomial.
c
      do j = 1, nd
        do i = 1, ni
          lb(i,j) = 1.0D+00
        end do
      end do

      do k = 1, nd

        do j = 1, nd

          if ( j .ne. k ) then

            do i = 1, ni
              lb(i,k) = lb(i,k) * ( xi(i) - xd(j) ) / ( xd(k) - xd(j) )
            end do

          end if

        end do

      end do

      return
      end
      subroutine lagrange_interp_nd_size ( m, n_1d, nd )

c*********************************************************************72
c
cc LAGRANGE_INTERP_ND_SIZE sizes an M-dimensional Lagrange interpolant.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the spatial dimension.
c
c    Input, integer N_1D(M), the order of the 1D rule to be used
c    in each dimension.
c
c    Output, integer ND, the number of points in the product grid.
c
      implicit none

      integer m

      integer i4vec_product
      integer n_1d(m)
      integer nd
c
c  Determine the number of data points.
c
      nd = i4vec_product ( m, n_1d )

      return
      end
      subroutine lagrange_interp_nd_grid ( m, n_1d, ab, nd, xd )

c*********************************************************************72
c
cc LAGRANGE_INTERP_ND_GRID sets an M-dimensional Lagrange interpolant grid.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the spatial dimension.
c
c    Input, integer N_1D(M), the order of the 1D rule to be used
c    in each dimension.
c
c    Input, double precision AB(M,2), the upper and lower interval endpoints in 
c    each dimension.
c
c    Input, integer ND, the number of points in the product grid.
c
c    Output, double precision XD(M,ND), the points at which data was sampled.
c
      implicit none

      integer m
      integer n_1d_max
      parameter ( n_1d_max = 100 )
      integer nd

      double precision ab(m,2)
      integer i
      integer j
      integer n
      integer n_1d(m)
      double precision x_1d(n_1d_max)
      double precision xd(m,nd)
c
c  Compute the data points.
c
      do j = 1, nd
        do i = 1, m
          xd(i,j) = 0.0D+00
        end do
      end do

      do i = 1, m

        n = n_1d(i)
        call cc_compute_points ( n, x_1d )

        do j = 1, n
          x_1d(j) = 0.5D+00 * ( ( 1.0D+00 - x_1d(j) ) * ab(i,1) 
     &                        + ( 1.0D+00 + x_1d(j) ) * ab(i,2) )
        end do

        call r8vec_direct_product ( i, n, x_1d, m, nd, xd )

      end do

      return
      end
      subroutine lagrange_interp_nd_grid2 ( m, ind, ab, nd, xd )

c*********************************************************************72
c
cc LAGRANGE_INTERP_ND_GRID2 sets an M-dimensional Lagrange interpolant grid.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the spatial dimension.
c
c    Input, integer IND(M), the index or level of the 1D rule 
c    to be used in each dimension.
c
c    Input, double precision AB(M,2), the upper and lower interval endpoints 
c    in each dimension.
c
c    Input, integer ND, the number of points in the product grid.
c
c    Output, double precision XD(M,ND), the points at which data was sampled.
c
c    Output, double precision ZD(ND), the function evaluated at the points XD.
c
      implicit none

      integer m
      integer n_1d_max
      parameter ( n_1d_max = 100 )
      integer nd

      double precision ab(m,2)
      integer i
      integer ind(m)
      integer j
      integer n
      double precision x_1d(n_1d_max)
      double precision xd(m,nd)
c
c  Compute the data points.
c
      do j = 1, nd
        do i = 1, m
          xd(i,j) = 0.0D+00
        end do
      end do

      do i = 1, m

        call order_from_level_135 ( ind(i), n )
  
        call cc_compute_points ( n, x_1d )

        do j = 1, n
          x_1d(j) = 0.5D+00 * ( ( 1.0D+00 - x_1d(j) ) * ab(i,1) 
     &                        + ( 1.0D+00 + x_1d(j) ) * ab(i,2) )
        end do

        call r8vec_direct_product ( i, n, x_1d, m, nd, xd )

      end do

      return
      end
      subroutine lagrange_interp_nd_size2 ( m, ind, nd )

c*********************************************************************72
c
cc LAGRANGE_INTERP_ND_SIZE2 sizes an M-dimensional Lagrange interpolant.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the spatial dimension.
c
c    Input, integer IND(M), the index or level of the 1D rule 
c    to be used in each dimension.
c
c    Output, integer ND, the number of points in the product grid.
c
      implicit none

      integer m

      integer i
      integer ind(m)
      integer n
      integer nd
c
c  Determine the number of data points.
c
      nd = 1
      do i = 1, m
        call order_from_level_135 ( ind(i), n )
        nd = nd * n
      end do

      return
      end
      subroutine lagrange_interp_nd_value ( m, n_1d, ab, nd, zd, ni, 
     &  xi, zi )

c*********************************************************************72
c
cc LAGRANGE_INTERP_ND_VALUE evaluates an ND Lagrange interpolant.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the spatial dimension.
c
c    Input, integer N_1D(M), the order of the 1D rule to be used
c    in each dimension.
c
c    Input, double precision AB(M,2), the upper and lower interval endpoints 
c    in each dimension.
c
c    Input, integer ND, the number of points in the product grid.
c
c    Input, double precision ZD(ND), the function evaluated at the points XD.
c
c    Input, integer NI, the number of points at which the 
c    interpolant is to be evaluated.
c
c    Input, double precision XI(M,NI), the points at which the interpolant is 
c    to be evaluated.
c
c    Output, double precision ZI(NI), the interpolant evaluated at the 
c    points XI.
c
      implicit none

      integer m
      integer n_1d_max
      parameter ( n_1d_max = 100 )
      integer nd
      integer ni

      double precision ab(m,2)
      integer i
      integer j
      integer k
      integer n
      integer n_1d(m)
      double precision r8vec_dot_product
      double precision value(n_1d_max)
      double precision w(nd)
      double precision x_1d(n_1d_max)
      double precision xi(m,ni)
      double precision zd(nd)
      double precision zi(ni)

      do j = 1, ni

        do i = 1, nd
          w(i) = 1.0D+00
        end do

        do i = 1, m

          n = n_1d(i)
 
          call cc_compute_points ( n, x_1d )

          do k = 1, n
            x_1d(k) = 0.5D+00 * ( ( 1.0D+00 - x_1d(k) ) * ab(i,1) 
     &                          + ( 1.0D+00 + x_1d(k) ) * ab(i,2) )
          end do

          call lagrange_base_1d ( n, x_1d, 1, xi(i,j), value )
          call r8vec_direct_product2 ( i, n, value, m, nd, w )

        end do

        zi(j) = r8vec_dot_product ( nd, w, zd )

      end do

      return
      end
      subroutine lagrange_interp_nd_value2 ( m, ind, ab, nd, zd, ni, 
     &  xi, zi )

c*********************************************************************72
c
cc LAGRANGE_INTERP_ND_VALUE2 evaluates an ND Lagrange interpolant.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the spatial dimension.
c
c    Input, integer IND(M), the index or level of the 1D rule 
c    to be used in each dimension.
c
c    Input, double precision AB(M,2), the upper and lower interval endpoints 
c    in each dimension.
c
c    Input, integer ND, the number of points in the product grid.
c
c    Input, double precision ZD(ND), the function evaluated at the points XD.
c
c    Input, integer NI, the number of points at which the 
c    interpolant is to be evaluated.
c
c    Input, double precision XI(M,NI), the points at which the interpolant 
c    is to be evaluated.
c
c    Output, double precision ZI(NI), the interpolant evaluated at the 
c    points XI.
c
      implicit none

      integer m
      integer n_1d_max
      parameter ( n_1d_max = 100 )
      integer nd
      integer ni

      double precision ab(m,2)
      integer i
      integer ind(m)
      integer j
      integer k
      integer n
      double precision r8vec_dot_product
      double precision value(n_1d_max)
      double precision w(nd)
      double precision x_1d(n_1d_max)
      double precision xi(m,ni)
      double precision zd(nd)
      double precision zi(ni)

      do j = 1, ni

        do i = 1, nd
          w(i) = 1.0D+00
        end do

        do i = 1, m

          call order_from_level_135 ( ind(i), n )

          call cc_compute_points ( n, x_1d )

          do k = 1, n
            x_1d(k) = 0.5D+00 * ( ( 1.0D+00 - x_1d(k) ) * ab(i,1) 
     &                          + ( 1.0D+00 + x_1d(k) ) * ab(i,2) )
          end do

          call lagrange_base_1d ( n, x_1d, 1, xi(i,j), value )
          call r8vec_direct_product2 ( i, n, value, m, nd, w )

        end do

        zi(j) = r8vec_dot_product ( nd, w, zd )

      end do

      return
      end
      subroutine order_from_level_135 ( l, n )

c*********************************************************************72
c
cc ORDER_FROM_LEVEL_135 evaluates the 135 level-to-order relationship.
c
c  Discussion:
c
c    Clenshaw Curtis rules, and some others, often use the following
c    scheme:
c
c    L: 0  1  2  3   4   5
c    N: 1  3  5  9  17  33 ... 2^L+1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer L, the level, which should be 0 or greater.
c
c    Output, integer N, the order.
c
      implicit none

      integer l
      integer n

      if ( l .lt. 0 ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'ORDER_FROM_LEVEL_135 - Fatal errorc'
        write ( *, '(a)' ) '  Illegal input value of Lc'
        stop
      else if ( l .eq. 0 ) then
        n = 1
      else
        n = ( 2 ** l ) + 1
      end if

      return
      end
      subroutine r8mat_print ( m, n, a, title )

c*********************************************************************72
c
cc R8MAT_PRINT prints an R8MAT.
c
c  Discussion:
c
c    An R8MAT is an array of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 May 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows in A.
c
c    Input, integer N, the number of columns in A.
c
c    Input, double precision A(M,N), the matrix.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      character ( len = * ) title

      call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

      return
      end
      subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi,
     &  title )

c*********************************************************************72
c
cc R8MAT_PRINT_SOME prints some of an R8MAT.
c
c  Discussion:
c
c    An R8MAT is an array of R8's.
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
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, double precision A(M,N), an M by N matrix to be printed.
c
c    Input, integer ILO, JLO, the first row and column to print.
c
c    Input, integer IHI, JHI, the last row and column to print.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 5 )
      integer m
      integer n

      double precision a(m,n)
      character * ( 14 ) ctemp(incx)
      integer i
      integer i2hi
      integer i2lo
      integer ihi
      integer ilo
      integer inc
      integer j
      integer j2
      integer j2hi
      integer j2lo
      integer jhi
      integer jlo
      character * ( * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )

      if ( m .le. 0 .or. n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  (None)'
        return
      end if

      do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

        j2hi = j2lo + incx - 1
        j2hi = min ( j2hi, n )
        j2hi = min ( j2hi, jhi )

        inc = j2hi + 1 - j2lo

        write ( *, '(a)' ) ' '

        do j = j2lo, j2hi
          j2 = j + 1 - j2lo
          write ( ctemp(j2), '(i7,7x)') j
        end do

        write ( *, '(''  Col   '',5a14)' ) ( ctemp(j), j = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) ' '

        i2lo = max ( ilo, 1 )
        i2hi = min ( ihi, m )

        do i = i2lo, i2hi

          do j2 = 1, inc

            j = j2lo - 1 + j2

            write ( ctemp(j2), '(g14.6)' ) a(i,j)

          end do

          write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

        end do

      end do

      return
      end
      subroutine r8mat_uniform_01 ( m, n, seed, r )

c*********************************************************************72
c
cc R8MAT_UNIFORM_01 returns a unit pseudorandom R8MAT.
c
c  Discussion:
c
c    An R8MAT is an array of R8's.
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
c    Input, integer M, N, the number of rows and columns in the array.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R(M,N), the array of pseudorandom values.
c
      implicit none

      integer m
      integer n

      integer i
      integer j
      integer k
      integer seed
      double precision r(m,n)

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8MAT_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do j = 1, n

        do i = 1, m

          k = seed / 127773

          seed = 16807 * ( seed - k * 127773 ) - k * 2836

          if ( seed .lt. 0 ) then
            seed = seed + 2147483647
          end if

          r(i,j) = dble ( seed ) * 4.656612875D-10

        end do
      end do

      return
      end
      subroutine r8vec_direct_product ( factor_index, factor_order,
     &  factor_value, factor_num, point_num, x )

c*********************************************************************72
c
cc R8VEC_DIRECT_PRODUCT creates a direct product of R8VEC's.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    To explain what is going on here, suppose we had to construct
c    a multidimensional quadrature rule as the product of K rules
c    for 1D quadrature.
c
c    The product rule will be represented as a list of points and weights.
c
c    The J-th item in the product rule will be associated with
c      item J1 of 1D rule 1,
c      item J2 of 1D rule 2,
c      ...,
c      item JK of 1D rule K.
c
c    In particular,
c      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
c    and
c      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
c
c    So we can construct the quadrature rule if we can properly
c    distribute the information in the 1D quadrature rules.
c
c    This routine carries out that task for the abscissas X.
c
c    Another way to do this would be to compute, one by one, the
c    set of all possible indices (J1,J2,...,JK), and then index
c    the appropriate information.  An advantage of the method shown
c    here is that you can process the K-th set of information and
c    then discard it.
c
c  Example:
c
c    Rule 1:
c      Order = 4
c      X(1:4) = ( 1, 2, 3, 4 )
c
c    Rule 2:
c      Order = 3
c      X(1:3) = ( 10, 20, 30 )
c
c    Rule 3:
c      Order = 2
c      X(1:2) = ( 100, 200 )
c
c    Product Rule:
c      Order = 24
c      X(1:24) =
c        ( 1, 10, 100 )
c        ( 2, 10, 100 )
c        ( 3, 10, 100 )
c        ( 4, 10, 100 )
c        ( 1, 20, 100 )
c        ( 2, 20, 100 )
c        ( 3, 20, 100 )
c        ( 4, 20, 100 )
c        ( 1, 30, 100 )
c        ( 2, 30, 100 )
c        ( 3, 30, 100 )
c        ( 4, 30, 100 )
c        ( 1, 10, 200 )
c        ( 2, 10, 200 )
c        ( 3, 10, 200 )
c        ( 4, 10, 200 )
c        ( 1, 20, 200 )
c        ( 2, 20, 200 )
c        ( 3, 20, 200 )
c        ( 4, 20, 200 )
c        ( 1, 30, 200 )
c        ( 2, 30, 200 )
c        ( 3, 30, 200 )
c        ( 4, 30, 200 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer FACTOR_INDEX, the index of the factor being processed.
c    The first factor processed must be factor 1!
c
c    Input, integer FACTOR_ORDER, the order of the factor.
c
c    Input, double precision FACTOR_VALUE(FACTOR_ORDER), the factor values
c    for factor FACTOR_INDEX.
c
c    Input, integer FACTOR_NUM, the number of factors.
c
c    Input, integer POINT_NUM, the number of elements in the direct product.
c
c    Input/output, double precision X(FACTOR_NUM,POINT_NUM), the elements of the
c    direct product, which are built up gradually.  Before the first call,
c    X might be set to 0.  After each factor has been input, X should
c    have the correct value.
c
c  Local Parameters:
c
c    Local, integer START, the first location of a block of values to set.
c
c    Local, integer CONTIG, the number of consecutive values to set.
c
c    Local, integer SKIP, the distance from the current value of START
c    to the next location of a block of values to set.
c
c    Local, integer REP, the number of blocks of values to set.
c
      implicit none

      integer factor_num
      integer factor_order
      integer point_num

      integer contig
      integer factor_index
      double precision factor_value(factor_order)
      integer i
      integer j
      integer k
      integer rep
      integer skip
      integer start
      double precision x(factor_num,point_num)

      save contig
      save rep
      save skip

      data contig / 0 /
      data rep / 0 /
      data skip / 0 /

      if ( factor_index .eq. 1 ) then
        contig = 1
        skip = 1
        rep = point_num
      end if

      rep = rep / factor_order
      skip = skip * factor_order

      do j = 1, factor_order

        start = 1 + ( j - 1 ) * contig

        do k = 1, rep
          do i = start, start+contig-1
            x(factor_index,i) = factor_value(j)
          end do
          start = start + skip
        end do

      end do

      contig = contig * factor_order

      return
      end
      subroutine r8vec_direct_product2 ( factor_index, factor_order,
     &  factor_value, factor_num, point_num, w )

c*********************************************************************72
c
cc R8VEC_DIRECT_PRODUCT2 creates a direct product of R8VEC's.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    To explain what is going on here, suppose we had to construct
c    a multidimensional quadrature rule as the product of K rules
c    for 1D quadrature.
c
c    The product rule will be represented as a list of points and weights.
c
c    The J-th item in the product rule will be associated with
c      item J1 of 1D rule 1,
c      item J2 of 1D rule 2,
c      ...,
c      item JK of 1D rule K.
c
c    In particular,
c      X(J) = ( X(1,J1), X(2,J2), ..., X(K,JK))
c    and
c      W(J) = W(1,J1) * W(2,J2) * ... * W(K,JK)
c
c    So we can construct the quadrature rule if we can properly
c    distribute the information in the 1D quadrature rules.
c
c    This routine carries out the task involving the weights W.
c
c    Another way to do this would be to compute, one by one, the
c    set of all possible indices (J1,J2,...,JK), and then index
c    the appropriate information.  An advantage of the method shown
c    here is that you can process the K-th set of information and
c    then discard it.
c
c  Example:
c
c    Rule 1:
c      Order = 4
c      W(1:4) = ( 2, 3, 5, 7 )
c
c    Rule 2:
c      Order = 3
c      W(1:3) = ( 11, 13, 17 )
c
c    Rule 3:
c      Order = 2
c      W(1:2) = ( 19, 23 )
c
c    Product Rule:
c      Order = 24
c      W(1:24) =
c        ( 2 * 11 * 19 )
c        ( 3 * 11 * 19 )
c        ( 4 * 11 * 19 )
c        ( 7 * 11 * 19 )
c        ( 2 * 13 * 19 )
c        ( 3 * 13 * 19 )
c        ( 5 * 13 * 19 )
c        ( 7 * 13 * 19 )
c        ( 2 * 17 * 19 )
c        ( 3 * 17 * 19 )
c        ( 5 * 17 * 19 )
c        ( 7 * 17 * 19 )
c        ( 2 * 11 * 23 )
c        ( 3 * 11 * 23 )
c        ( 5 * 11 * 23 )
c        ( 7 * 11 * 23 )
c        ( 2 * 13 * 23 )
c        ( 3 * 13 * 23 )
c        ( 5 * 13 * 23 )
c        ( 7 * 13 * 23 )
c        ( 2 * 17 * 23 )
c        ( 3 * 17 * 23 )
c        ( 5 * 17 * 23 )
c        ( 7 * 17 * 23 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer FACTOR_INDEX, the index of the factor being processed.
c    The first factor processed must be factor 1!
c
c    Input, integer FACTOR_ORDER, the order of the factor.
c
c    Input, double precision FACTOR_VALUE(FACTOR_ORDER), the factor values
c    for factor FACTOR_INDEX.
c
c    Input, integer FACTOR_NUM, the number of factors.
c
c    Input, integer POINT_NUM, the number of elements in the direct product.
c
c    Input/output, double precision W(POINT_NUM), the elements of the
c    direct product, which are built up gradually.  Before the first call,
c    W should be set to 1.
c
c  Local Parameters:
c
c    Local, integer START, the first location of a block of values to set.
c
c    Local, integer CONTIG, the number of consecutive values to set.
c
c    Local, integer SKIP, the distance from the current value of START
c    to the next location of a block of values to set.
c
c    Local, integer REP, the number of blocks of values to set.
c
      implicit none

      integer factor_num
      integer factor_order
      integer point_num

      integer contig
      integer factor_index
      double precision factor_value(factor_order)
      integer i
      integer j
      integer k
      integer rep
      integer skip
      integer start
      double precision w(point_num)

      save contig
      save rep
      save skip

      data contig / 0 /
      data rep / 0 /
      data skip / 0 /

      if ( factor_index .eq. 1 ) then
        contig = 1
        skip = 1
        rep = point_num
      end if

      rep = rep / factor_order
      skip = skip * factor_order

      do j = 1, factor_order

        start = 1 + ( j - 1 ) * contig

        do k = 1, rep
          do i = start, start+contig-1
            w(i) = w(i) * factor_value(j)
          end do
          start = start + skip
        end do

      end do

      contig = contig * factor_order

      return
      end
      function r8vec_dot_product ( n, v1, v2 )

c*********************************************************************72
c
cc R8VEC_DOT_PRODUCT finds the dot product of a pair of R8VEC's.
c
c  Discussion:
c
c    An R8VEC is a vector of R8 values.
c
c    In FORTRAN90, the system routine DOT_PRODUCT should be called
c    directly.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the vectors.
c
c    Input, double precision V1(N), V2(N), the vectors.
c
c    Output, double precision R8VEC_DOT_PRODUCT, the dot product.
c
      implicit none

      integer n

      integer i
      double precision r8vec_dot_product
      double precision v1(n)
      double precision v2(n)
      double precision value

      value = 0.0D+00
      do i = 1, n
        value = value + v1(i) * v2(i)
      end do

      r8vec_dot_product = value

      return
      end
      function r8vec_norm_affine ( n, v0, v1 )

c*********************************************************************72
c
cc R8VEC_NORM_AFFINE returns the affine norm of an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    The affine vector L2 norm is defined as:
c
c      R8VEC_NORM_AFFINE(V0,V1)
c        = sqrt ( sum ( 1 <= I <= N ) ( V1(I) - V0(I) )^2 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 October 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the vectors.
c
c    Input, double precision V0(N), the base vector.
c
c    Input, double precision V1(N), the vector whose affine norm is desired.
c
c    Output, double precision R8VEC_NORM_AFFINE, the L2 norm of V1-V0.
c
      implicit none

      integer n

      integer i
      double precision r8vec_norm_affine
      double precision v0(n)
      double precision v1(n)

      r8vec_norm_affine = 0.0D+00
      do i = 1, n
        r8vec_norm_affine = r8vec_norm_affine
     &    + ( v0(i) - v1(i) )**2
      end do
      r8vec_norm_affine = sqrt ( r8vec_norm_affine )

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
