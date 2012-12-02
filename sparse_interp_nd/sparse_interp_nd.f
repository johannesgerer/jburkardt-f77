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
c    08 October 2008
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
      subroutine comp_next ( n, k, a, more, h, t )

c*********************************************************************72
c
cc COMP_NEXT computes the compositions of the integer N into K parts.
c
c  Discussion:
c
c    A composition of the integer N into K parts is an ordered sequence
c    of K nonnegative integers which sum to N.  The compositions (1,2,1)
c    and (1,1,2) are considered to be distinct.
c
c    The routine computes one composition on each call until there are no more.
c    For instance, one composition of 6 into 3 parts is
c    3+2+1, another would be 6+0+0.
c
c    On the first call to this routine, set MORE = FALSE.  The routine
c    will compute the first element in the sequence of compositions, and
c    return it, as well as setting MORE = TRUE.  If more compositions
c    are desired, call again, and again.  Each time, the routine will
c    return with a new composition.
c
c    However, when the LAST composition in the sequence is computed
c    and returned, the routine will reset MORE to FALSE, signaling that
c    the end of the sequence has been reached.
c
c    This routine originally used a SAVE statement to maintain the
c    variables H and T.  I have decided that it is safer
c    to pass these variables as arguments, even though the user should
c    never alter them.  This allows this routine to safely shuffle
c    between several ongoing calculations.
c
c
c    There are 28 compositions of 6 into three parts.  This routine will
c    produce those compositions in the following order:
c
c     I         A
c     -     ---------
c     1     6   0   0
c     2     5   1   0
c     3     4   2   0
c     4     3   3   0
c     5     2   4   0
c     6     1   5   0
c     7     0   6   0
c     8     5   0   1
c     9     4   1   1
c    10     3   2   1
c    11     2   3   1
c    12     1   4   1
c    13     0   5   1
c    14     4   0   2
c    15     3   1   2
c    16     2   2   2
c    17     1   3   2
c    18     0   4   2
c    19     3   0   3
c    20     2   1   3
c    21     1   2   3
c    22     0   3   3
c    23     2   0   4
c    24     1   1   4
c    25     0   2   4
c    26     1   0   5
c    27     0   1   5
c    28     0   0   6
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 July 2008
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    FORTRAN90 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the integer whose compositions are desired.
c
c    Input, integer K, the number of parts in the composition.
c
c    Input/output, integer A(K), the parts of the composition.
c
c    Input/output, logical MORE, set by the user to start the computation,
c    and by the routine to terminate it.
c
c    Input/output, integer H, T, two internal parameters needed
c    for the computation.  The user should allocate space for these in the
c    calling program, include them in the calling sequence, but never alter
c    themc
c
      implicit none

      integer k

      integer a(k)
      integer h
      integer i
      logical more
      integer n
      integer t
c
c  The first computation.
c
      if ( .not. more ) then

        t = n
        h = 0
        a(1) = n
        do i = 2, k
          a(i) = 0
        end do
c
c  The next computation.
c
      else
c
c  If the first entry A(1) is positive, then set H to zero,
c  so that when we increment H, it points to A(1); we will decrement A(1) by 1
c  and increment A(2).
c
        if ( 1 .lt. t ) then
          h = 0
        end if
c
c  Otherwise, A(1) is 0.  Then by H + 1 is the entry we incremented last time.
c  Set H = H + 1, zero A(H), adding all but one of its value to A(1),
c  and incrementing A(H+1) by 1.
c
        h = h + 1
        t = a(h)
        a(h) = 0
        a(1) = t - 1
        a(h+1) = a(h+1) + 1

      end if
c
c  This is the last element of the sequence if all the
c  items are in the last slot.
c
      more = ( a(k) .ne. n )

      return
      end
      function i4_choose ( n, k )

c*********************************************************************72
c
cc I4_CHOOSE computes the binomial coefficient C(N,K).
c
c  Discussion:
c
c    The value is calculated in such a way as to avoid overflow and
c    roundoff.  The calculation is done in integer arithmetic.
c
c    The formula used is:
c
c      C(N,K) = N! / ( K! * (N-K)! )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 June 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    ML Wolfson, HV Wright,
c    Algorithm 160:
c    Combinatorial of M Things Taken N at a Time,
c    Communications of the ACM,
c    Volume 6, Number 4, April 1963, page 161.
c
c  Parameters:
c
c    Input, integer N, K, are the values of N and K.
c
c    Output, integer I4_CHOOSE, the number of combinations of N
c    things taken K at a time.
c
      implicit none

      integer i
      integer i4_choose
      integer k
      integer mn
      integer mx
      integer n
      integer value

      mn = min ( k, n - k )

      if ( mn .lt. 0 ) then

        value = 0

      else if ( mn .eq. 0 ) then

        value = 1

      else

        mx = max ( k, n - k )
        value = mx + 1

        do i = 2, mn
          value = ( value * ( mx + i ) ) / i
        end do

      end if

      i4_choose = value

      return
      end
      function i4_mop ( i )

c*********************************************************************72
c
cc I4_MOP returns the I-th power of -1 as an I4 value.
c
c  Discussion:
c
c    An I4 is an integer value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the power of -1.
c
c    Output, integer I4_MOP, the I-th power of -1.
c
      implicit none

      integer i
      integer i4_mop

      if ( mod ( i, 2 ) .eq. 0 ) then
        i4_mop = 1
      else
        i4_mop = -1
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
      subroutine lagrange_basis_1d ( nd, xd, ni, xi, lb ) 

c*********************************************************************72
c
cc LAGRANGE_BASIS_1D evaluates a 1D Lagrange basis.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 October 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer ND, the number of data points.
c
c    Input, double precision XD(ND), the interpolation nodes.
c
c    Input, integer NI, the number of evaluation points.
c
c    Input, double precision XI(NI), the evaluation points.
c
c    Output, double precision LB(NI,ND), the value, at the I-th point XI, 
c    of the Jth basis function.
c
      implicit none

      integer nd
      integer ni

      integer i
      integer j
      integer k
      double precision lb(ni,nd)
      double precision t
      double precision xd(nd)
      double precision xi(ni)
      
      do i = 1, ni
        do j = 1, nd
          t = 1.0D+00
          do k = 1, j - 1
            t = t * ( xi(i) - xd(k) ) / ( xd(j) - xd(k) )
          end do
          do k = j + 1, nd
            t = t * ( xi(i) - xd(k) ) / ( xd(j) - xd(k) )
          end do
          lb(i,j) = t
        end do
      end do

      return
      end
      subroutine lagrange_interp_nd_grid2 ( m, ind, a, b, nd, xd )

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
c    to be used in each dimension.  Each entry of IND must be
c    no greater than 10.
c
c    Input, double precision A(M), B(M), the lower and upper limits.
c
c    Input, integer ND, the number of points in the product grid.
c
c    Output, double precision XD(M,ND), the points at which data was sampled.
c
      implicit none

      integer m
      integer nd

      integer n_max
      parameter ( n_max = 1025 )

      double precision a(m)
      double precision b(m)
      integer i
      integer ind(m)
      integer j
      integer n
      double precision x_1d(n_max)
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

        if ( n_max .lt. n ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'LAGRANGE_INTERP_ND_GRID2 - Fatal error!'
          write ( *, '(a,i4)' ) '  N_MAX = ', n_max
          write ( *, '(a,i4)' ) '  N = ', n
          stop
        end if

        call cc_compute_points ( n, x_1d )

        do j = 1, n
          x_1d(j) = 0.5D+00 * ( ( 1.0D+00 - x_1d(j) ) * a(i) 
     &                        + ( 1.0D+00 + x_1d(j) ) * b(i) )
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
      subroutine lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, 
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
c    Input, double precision A(M), B(M), the lower and upper limits.
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
      integer nd
      integer ni

      integer n_max
      parameter ( n_max = 1025 )

      double precision a(m)
      double precision b(m)
      integer i
      integer ind(m)
      integer j
      integer k
      integer n
      double precision r8vec_dot_product
      double precision value(n_max)
      double precision w(nd)
      double precision x_1d(n_max)
      double precision xi(m,ni)
      double precision zd(nd)
      double precision zi(ni)

      do j = 1, ni

        do i = 1, nd
          w(i) = 1.0D+00
        end do

        do i = 1, m

          call order_from_level_135 ( ind(i), n )

          if ( n_max .lt. n ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 
     &        'LAGRANGE_INTERP_ND_VALUE2 - Fatal error!'
            write ( *, '(a,i4)' ) '  N_MAX = ', n_max
            write ( *, '(a,i4)' ) '  N = ', n
            stop
          end if

          call cc_compute_points ( n, x_1d )

          do k = 1, n
            x_1d(k) = 0.5D+00 * ( ( 1.0D+00 - x_1d(k) ) * a(i) 
     &                          + ( 1.0D+00 + x_1d(k) ) * b(i) )
          end do

          call lagrange_basis_1d ( n, x_1d, 1, xi(i,j), value )
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
c    18 July 2012
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
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'ORDER_FROM_LEVEL_135 - Fatal error!'
        write ( *, '(a)' ) '  Illegal input value of L!'
        stop
      else if ( l .eq. 0 ) then
        n = 1
      else
        n = ( 2 ** l ) + 1
      end if

      return
      end
      subroutine smolyak_coefficients ( l_max, m, c, w )

c*********************************************************************72
c
cc SMOLYAK_COEFFICIENTS returns the Smolyak coefficients and counts.
c
c  Discussion:
c
c    The Smolyak sparse interpolant can be written as:
c
c      A(L,M)(X) = sum ( L-M+1 <= |L| <= L_max ) 
c        C(|L|) * g(l1)(x1) * g(l2)(x2) * ... * g(lm)(xm).
c
c    where:
c
c    * L=(l1,l2,...,lm) is a vector of M nonnegative integers;
c    * |L| is the sum of the entries of L;
c    * X=(x1,x2,...,xm) is an M-dimensional point in a product space;
c    * g(i)(xj) is the i-th 1-d interpolation function in dimension j;
c
c    Note that:
c
c    * W(|L|) will represent the number of distinct interpolants for which
c      the sublevel, or sum of the L vector entries, is |L|;
c
c    * the coefficients C and counts W will be zero for sublevels 
c      0 through L_MAX - M (and MATLAB indices 1 through L_MAX-M+1).
c
c    * it will be the case that W' * C = 1, essentially because the interpolant
c      to the identity function must be the identity function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 July 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer L_MAX, the (maximum) level.
c    0 <= L_MAX.
c
c    Input, integer M, the spatial dimension.
c    1 <= M.
c
c    Output, integer C(0:L_MAX), the coefficients for objects 
c    at sublevels 0 through L_MAX.
c
c    Output, integer W(0:L_MAX), the number of objects at 
c    sublevels 0 through L_MAX.
c
      implicit none

      integer l_max

      integer c(0:l_max)
      integer i4_choose
      integer i4_mop
      integer l
      integer l_min
      integer m
      integer w(0:l_max)

      l_min = max ( l_max - m + 1, 0 )

      c(0:l_min-1) = 0
      do l = l_min, l_max
        c(l) = i4_mop ( l_max - l ) * i4_choose ( m - 1, l_max - l )
      end do

      w(0:l_min-1) = 0
      do l = l_min, l_max
        w(l) = i4_choose ( l + m - 1, m - 1 )
      end do

      return
      end
