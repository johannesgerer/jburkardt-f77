      subroutine bernstein_matrix ( n, a )

c*********************************************************************72
c
cc BERNSTEIN_MATRIX returns the Bernstein matrix.
c
c  Discussion:
c
c    The Bernstein matrix of order N is an NxN matrix A which can be used to
c    transform a vector of power basis coefficients C representing a polynomial 
c    P(X) to a corresponding Bernstein basis coefficient vector B:
c
c      B = A * C
c
c    The N power basis vectors are ordered as (1,X,X^2,...X^(N-1)) and the N 
c    Bernstein basis vectors as ((1-X)^(N-1), X*(1_X)^(N-2),...,X^(N-1)).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Output, double precision A(N,N), the Bernstein matrix.
c
      implicit none

      integer n

      double precision a(n,n)
      integer i0
      integer j0
      integer n0
      double precision r8_choose
      double precision r8_mop

      n0 = n - 1

      do j0 = 0, n0
        do i0 = 0, j0
          a(i0+1,j0+1) = r8_mop ( j0 - i0 ) 
     &      * r8_choose ( n0 - i0, j0 - i0 ) * r8_choose ( n0, i0 )
        end do
        do i0 = j0 + 1, n0
          a(i0+1,j0+1) = 0.0D+00
        end do
      end do

      return
      end
      subroutine bernstein_matrix_inverse ( n, a )

c*********************************************************************72
c
cc BERNSTEIN_MATRIX_INVERSE returns the inverse Bernstein matrix.
c
c  Discussion:
c
c    The inverse Bernstein matrix of order N is an NxN matrix A which can 
c    be used to transform a vector of Bernstein basis coefficients B
c    representing a polynomial P(X) to a corresponding power basis 
c    coefficient vector C:
c
c      C = A * B
c
c    The N power basis vectors are ordered as (1,X,X^2,...X^(N-1)) and the N 
c    Bernstein basis vectors as ((1-X)^(N-1), X*(1_X)^(N-2),...,X^(N-1)).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Output, double precision A(N,N), the inverse Bernstein matrix.
c
      implicit none

      integer n

      double precision a(n,n)
      integer i0
      integer j0
      integer n0
      double precision r8_choose

      n0 = n - 1

      do j0 = 0, n0
        do i0 = 0, j0
          a(i0+1,j0+1) = r8_choose ( j0, i0 ) / r8_choose ( n0, i0 )
        end do
        do i0 = j0 + 1, n0
          a(i0+1,j0+1) = 0.0D+00
        end do
      end do

      return
      end
      subroutine bernstein_poly ( n, x, bern )

c*********************************************************************72
c
cc BERNSTEIN_POLY evaluates the Bernstein polynomials at a point X.
c
c  Discussion:
c
c    The Bernstein polynomials are assumed to be based on [0,1].
c
c    The formula is:
c
c      B(N,I)(X) = [N!/(I!*(N-I)!)] * (1-X)^(N-I) * X^I
c
c  First values:
c
c    B(0,0)(X) = 1
c
c    B(1,0)(X) =      1-X
c    B(1,1)(X) =                X
c
c    B(2,0)(X) =     (1-X)^2
c    B(2,1)(X) = 2 * (1-X)    * X
c    B(2,2)(X) =                X^2
c
c    B(3,0)(X) =     (1-X)^3
c    B(3,1)(X) = 3 * (1-X)^2 * X
c    B(3,2)(X) = 3 * (1-X)   * X^2
c    B(3,3)(X) =               X^3
c
c    B(4,0)(X) =     (1-X)^4
c    B(4,1)(X) = 4 * (1-X)^3 * X
c    B(4,2)(X) = 6 * (1-X)^2 * X^2
c    B(4,3)(X) = 4 * (1-X)   * X^3
c    B(4,4)(X) =               X^4
c
c  Special values:
c
c    B(N,I)(X) has a unique maximum value at X = I/N.
c
c    B(N,I)(X) has an I-fold zero at 0 and and N-I fold zero at 1.
c
c    B(N,I)(1/2) = C(N,K) / 2^N
c
c    For a fixed X and N, the polynomials add up to 1:
c
c      Sum ( 0 <= I <= N ) B(N,I)(X) = 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the degree of the Bernstein polynomials 
c    to be used.  For any N, there is a set of N+1 Bernstein polynomials,
c    each of degree N, which form a basis for polynomials on [0,1].
c
c    Input, double precision X, the evaluation point.
c
c    Output, double precision BERN(0:N), the values of the N+1 
c    Bernstein polynomials at X.
c
      implicit none

      integer n

      double precision bern(0:n)
      integer i
      integer j
      double precision x

      if ( n .eq. 0 ) then
     
        bern(0) = 1.0D+00
     
      else if ( 0 .lt. n ) then
     
        bern(0) = 1.0D+00 - x
        bern(1) = x
     
        do i = 2, n
          bern(i) = x * bern(i-1)
          do j = i - 1, 1, -1
            bern(j) =             x   * bern(j-1) 
     &              + ( 1.0D+00 - x ) * bern(j)
          end do
          bern(0) = ( 1.0D+00 - x ) * bern(0)
        end do
     
      end if
     
      return
      end
      subroutine bernstein_poly_values ( n_data, n, k, x, b )

c*********************************************************************72
c
cc BERNSTEIN_POLY_VALUES returns some values of the Bernstein polynomials.
c
c  Discussion:
c
c    The Bernstein polynomials are assumed to be based on [0,1].
c
c    The formula for the Bernstein polynomials is
c
c      B(N,I)(X) = [Nc/(Ic*(N-I)c)] * (1-X)**(N-I) * X**I
c
c    In Mathematica, the function can be evaluated by:
c
c      Binomial[n,i] * (1-x)^(n-i) * x^i
c
c  First values:
c
c    B(0,0)(X) = 1
c
c    B(1,0)(X) =      1-X
c    B(1,1)(X) =                X
c
c    B(2,0)(X) =     (1-X)**2
c    B(2,1)(X) = 2 * (1-X)    * X
c    B(2,2)(X) =                X**2
c
c    B(3,0)(X) =     (1-X)**3
c    B(3,1)(X) = 3 * (1-X)**2 * X
c    B(3,2)(X) = 3 * (1-X)    * X**2
c    B(3,3)(X) =                X**3
c
c    B(4,0)(X) =     (1-X)**4
c    B(4,1)(X) = 4 * (1-X)**3 * X
c    B(4,2)(X) = 6 * (1-X)**2 * X**2
c    B(4,3)(X) = 4 * (1-X)    * X**3
c    B(4,4)(X) =                X**4
c
c  Special values:
c
c    B(N,I)(X) has a unique maximum value at X = I/N.
c
c    B(N,I)(X) has an I-fold zero at 0 and and N-I fold zero at 1.
c
c    B(N,I)(1/2) = C(N,K) / 2**N
c
c    For a fixed X and N, the polynomials add up to 1:
c
c      Sum ( 0 <= I <= N ) B(N,I)(X) = 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 March 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Stephen Wolfram,
c    The Mathematica Book,
c    Fourth Edition,
c    Cambridge University Press, 1999,
c    ISBN: 0-521-64314-7,
c    LC: QA76.95.W65.
c
c  Parameters:
c
c    Input/output, integer N_DATA.  The user sets N_DATA to 0 before the
c    first call.  On each call, the routine increments N_DATA by 1, and
c    returns the corresponding data; when there is no more data, the
c    output value of N_DATA will be 0 again.
c
c    Output, integer N, the degree of the polynomial.
c
c    Output, integer K, the index of the polynomial.
c
c    Output, double precision X, the argument of the polynomial.
c
c    Output, double precision B, the value of the polynomial B(N,K)(X).
c
      implicit none

      integer n_max
      parameter ( n_max = 15 )

      double precision b
      double precision b_vec(n_max)
      integer k
      integer k_vec(n_max)
      integer n
      integer n_data
      integer n_vec(n_max)
      double precision x
      double precision x_vec(n_max)

      save b_vec
      save k_vec
      save n_vec
      save x_vec

      data b_vec /
     &  0.1000000000000000D+01,
     &  0.7500000000000000D+00,
     &  0.2500000000000000D+00,
     &  0.5625000000000000D+00,
     &  0.3750000000000000D+00,
     &  0.6250000000000000D-01,
     &  0.4218750000000000D+00,
     &  0.4218750000000000D+00,
     &  0.1406250000000000D+00,
     &  0.1562500000000000D-01,
     &  0.3164062500000000D+00,
     &  0.4218750000000000D+00,
     &  0.2109375000000000D+00,
     &  0.4687500000000000D-01,
     &  0.3906250000000000D-02 /
      data k_vec /
     &  0,
     &  0, 1,
     &  0, 1, 2,
     &  0, 1, 2, 3,
     &  0, 1, 2, 3, 4 /
      data n_vec /
     &  0,
     &  1, 1,
     &  2, 2, 2,
     &  3, 3, 3, 3,
     &  4, 4, 4, 4, 4 /
      data x_vec /
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00,
     &  0.25D+00 /

      if ( n_data .lt. 0 ) then
        n_data = 0
      end if

      n_data = n_data + 1

      if ( n_max .lt. n_data ) then
        n_data = 0
        n = 0
        k = 0
        x = 0.0D+00
        b = 0.0D+00
      else
        n = n_vec(n_data)
        k = k_vec(n_data)
        x = x_vec(n_data)
        b = b_vec(n_data)
      end if

      return
      end
      subroutine bpab ( n, a, b, x, bern )

c*********************************************************************72
c
cc BPAB evaluates at X the Bernstein polynomials based in [A,B].
c
c  Discussion:
c
c    The formula is:
c
c      BERN(N,I)(X) = [N!/(I!*(N-I)!)] * (B-X)^(N-I) * (X-A)^I / (B-A)^N
c
c  First values:
c
c    B(0,0)(X) =   1
c
c    B(1,0)(X) = (      B-X                ) / (B-A)
c    B(1,1)(X) = (                 X-A     ) / (B-A)
c
c    B(2,0)(X) = (     (B-X)^2             ) / (B-A)^2
c    B(2,1)(X) = ( 2 * (B-X)    * (X-A)    ) / (B-A)^2
c    B(2,2)(X) = (                (X-A)^2  ) / (B-A)^2
c
c    B(3,0)(X) = (     (B-X)^3             ) / (B-A)^3
c    B(3,1)(X) = ( 3 * (B-X)^2  * (X-A)    ) / (B-A)^3
c    B(3,2)(X) = ( 3 * (B-X)    * (X-A)^2  ) / (B-A)^3
c    B(3,3)(X) = (                (X-A)^3  ) / (B-A)^3
c
c    B(4,0)(X) = (     (B-X)^4             ) / (B-A)^4
c    B(4,1)(X) = ( 4 * (B-X)^3  * (X-A)    ) / (B-A)^4
c    B(4,2)(X) = ( 6 * (B-X)^2  * (X-A)^2  ) / (B-A)^4
c    B(4,3)(X) = ( 4 * (B-X)    * (X-A)^3  ) / (B-A)^4
c    B(4,4)(X) = (                (X-A)^4  ) / (B-A)^4
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the degree of the Bernstein polynomials 
c    to be used.  For any N, there is a set of N+1 Bernstein polynomials, 
c    each of degree N, which form a basis for polynomials on [A,B].
c
c    Input, double precision A, B, the endpoints of the interval on which the
c    polynomials are to be based.  A and B should not be equal.
c
c    Input, double precision X, the point at which the polynomials 
c    are to be evaluated.
c
c    Output, double precision BERN(0:N), the values of the N+1
c    Bernstein polynomials at X.
c
      implicit none

      integer n

      double precision a
      double precision b
      double precision bern(0:n)
      integer i
      integer j
      double precision x

      if ( b .eq. a ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BPAB - Fatal error!'
        write ( *, '(a,g14.6)' ) '  A = B = ', a
        stop
      end if

      if ( n .eq. 0 ) then
     
        bern(0) = 1.0D+00
     
      else if ( 0 .lt. n ) then
     
        bern(0) = ( b - x ) / ( b - a )
        bern(1) = ( x - a ) / ( b - a )
     
        do i = 2, n
          bern(i) = ( x - a ) * bern(i-1) / ( b - a )
          do j = i - 1, 1, -1
            bern(j) = ( ( b - x     ) * bern(j)     
     &                + (     x - a ) * bern(j-1) ) 
     &                / ( b     - a )
          end do
          bern(0) = ( b - x ) * bern(0) / ( b - a )
        end do
     
      end if
     
      return
      end
      subroutine bpab_approx ( n, a, b, ydata, nval, xval, yval )

c*********************************************************************72
c
cc BPAB_APPROX evaluates the Bernstein polynomial approximant to F(X) on [A,B].
c
c  Formula:
c
c    BPAB(F)(X) = sum ( 0 <= I <= N ) F(X(I)) * B_BASE(I,X)
c
c    where
c
c      X(I) = ( ( N - I ) * A + I * B ) / N
c      B_BASE(I,X) is the value of the I-th Bernstein basis polynomial at X.
c
c  Discussion:
c
c    The Bernstein polynomial BPAB(F) for F(X) over [A,B] is an approximant, 
c    not an interpolant; in other words, its value is not guaranteed to equal
c    that of F at any particular point.  However, for a fixed interval
c    [A,B], if we let N increase, the Bernstein polynomial converges
c    uniformly to F everywhere in [A,B], provided only that F is continuous.
c    Even if F is not continuous, but is bounded, the polynomial converges
c    pointwise to F(X) at all points of continuity.  On the other hand,
c    the convergence is quite slow compared to other interpolation
c    and approximation schemes.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, integer N, the degree of the Bernstein polynomial
c    to be used.  N must be at least 0.
c
c    Input, double precision A, B, the endpoints of the interval on which the
c    approximant is based.  A and B should not be equal.
c
c    Input, double precision YDATA(0:N), the data values at N+1 equally
c    spaced points in [A,B].  If N = 0, then the evaluation point should
c    be 0.5 * ( A + B).  Otherwise, evaluation point I should be
c    ( (N-I)*A + I*B ) / N ).
c
c    Input, integer NVAL, the number of points at which the
c    approximant is to be evaluated.
c
c    Input, double precision XVAL(NVAL), the point at which the Bernstein 
c    polynomial approximant is to be evaluated.  The entries of XVAL do not 
c    have to lie in the interval [A,B].
c
c    Output, double precision YVAL(NVAL), the values of the Bernstein 
c    polynomial approximant for F, based in [A,B], evaluated at XVAL.
c
      implicit none

      integer n
      integer nval

      double precision a
      double precision b
      double precision bvec(0:n)
      integer i
      double precision r8vec_dot_product
      double precision xval(nval)
      double precision ydata(0:n)
      double precision yval(nval)

      do i = 1, nval
c
c  Evaluate the Bernstein basis polynomials at XVAL.
c
        call bpab ( n, a, b, xval(i), bvec )
c
c  Now compute the sum of YDATA(I) * BVEC(I).
c
        yval(i) = r8vec_dot_product ( n + 1, ydata, bvec )

      end do

      return
      end
      function r8_choose ( n, k )

c*********************************************************************72
c
cc R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.
c
c  Discussion:
c
c    The value is calculated in such a way as to avoid overflow and
c    roundoff.  The calculation is done in R8 arithmetic.
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
c    07 June 2008
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
c    Output, double precision R8_CHOOSE, the number of combinations of N
c    things taken K at a time.
c
      implicit none

      integer i
      integer k
      integer mn
      integer mx
      integer n
      double precision r8_choose
      double precision value

      mn = min ( k, n - k )

      if ( mn .lt. 0 ) then

        value = 0.0D+00

      else if ( mn .eq. 0 ) then

        value = 1.0D+00

      else

        mx = max ( k, n - k )
        value = dble ( mx + 1 )

        do i = 2, mn
          value = ( value * dble ( mx + i ) ) / dble ( i )
        end do

      end if

      r8_choose = value

      return
      end
      function r8_mop ( i )

c*********************************************************************72
c
cc R8_MOP returns the I-th power of -1 as an R8 value.
c
c  Discussion:
c
c    An R8 is a double precision real value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the power of -1.
c
c    Output, double precision R8_MOP, the I-th power of -1.
c
      implicit none

      integer i
      double precision r8_mop

      if ( mod ( i, 2 ) .eq. 0 ) then
        r8_mop = + 1.0D+00
      else
        r8_mop = - 1.0D+00
      end if

      return
      end
      function r8_uniform_01 ( seed )

c*********************************************************************72
c
cc R8_UNIFORM_01 returns a unit pseudorandom R8.
c
c  Discussion:
c
c    This routine implements the recursion
c
c      seed = 16807 * seed mod ( 2**31 - 1 )
c      r8_uniform_01 = seed / ( 2**31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      R8_UNIFORM_01
c      SEED      SEED
c
c         12345   207482415  0.096616
c     207482415  1790989824  0.833995
c    1790989824  2035175616  0.947702
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
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley Interscience, page 95, 1998.
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
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R8_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      double precision r8_uniform_01
      integer k
      integer seed

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + 2147483647
      end if

      r8_uniform_01 = dble ( seed ) * 4.656612875D-10

      return
      end
      subroutine r8mat_is_identity ( n, a, error_frobenius )

c*********************************************************************72
c
cc R8MAT_IS_IDENTITY determines if an R8MAT is the identity.
c
c  Discussion:
c
c    An R8MAT is a matrix of double precision values.
c
c    The routine returns the Frobenius norm of A - I.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 November 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision A(N,N), the matrix.
c
c    Output, double precision ERROR_FROBENIUS, the Frobenius norm
c    of the difference matrix A - I, which would be exactly zero
c    if A were the identity matrix.
c
      implicit none

      integer n

      double precision a(n,n)
      double precision error_frobenius
      integer i
      integer j
      double precision value

      error_frobenius = 0.0D+00

      do i = 1, n
        do j = 1, n
          if ( i .eq. j ) then
            error_frobenius = error_frobenius + ( a(i,j) - 1.0D+00 )**2
          else
            error_frobenius = error_frobenius + a(i,j)**2
          end if
        end do 
      end do

      error_frobenius = sqrt ( error_frobenius )

      return
      end
      subroutine r8mat_mm ( n1, n2, n3, a, b, c )

c*********************************************************************72
c
cc R8MAT_MM multiplies two R8MAT's.
c
c  Discussion:
c
c    An R8MAT is an array of R8 values.
c
c    In FORTRAN90, this operation is more efficiently done by the
c    command:
c
c      C(1:N1,1:N3) = MATMUL ( A(1:N1,1;N2), B(1:N2,1:N3) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 December 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, N3, the order of the matrices.
c
c    Input, double precision A(N1,N2), B(N2,N3), the matrices to multiply.
c
c    Output, double precision C(N1,N3), the product matrix C = A * B.
c
      implicit none

      integer n1
      integer n2
      integer n3

      double precision a(n1,n2)
      double precision b(n2,n3)
      double precision c(n1,n3)
      integer i
      integer j
      integer k

      do i = 1, n1
        do j = 1, n3
          c(i,j) = 0.0D+00
          do k = 1, n2
            c(i,j) = c(i,j) + a(i,k) * b(k,j)
          end do
        end do
      end do

      return
      end
      function r8mat_norm_fro ( m, n, a )

c*********************************************************************72
c
cc R8MAT_NORM_FRO returns the Frobenius norm of an R8MAT.
c
c  Discussion:
c
c    An R8MAT is an array of R8's.
c
c    The Frobenius norm is defined as
c
c      R8MAT_NORM_FRO = sqrt (
c        sum ( 1 .le. I .le. M ) sum ( 1 .le. j .le. N ) A(I,J)**2 )
c
c    The matrix Frobenius norm is not derived from a vector norm, but
c    is compatible with the vector L2 norm, so that:
c
c      r8vec_norm_l2 ( A * x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 July 2008
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
c    Input, double precision A(M,N), the matrix whose Frobenius
c    norm is desired.
c
c    Output, double precision R8MAT_NORM_FRO, the Frobenius norm of A.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer i
      integer j
      double precision r8mat_norm_fro
      double precision value

      value = 0.0D+00
      do j = 1, n
        do i = 1, m
          value = value + a(i,j) * a(i,j)
        end do
      end do
      value = sqrt ( value )

      r8mat_norm_fro = value

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
