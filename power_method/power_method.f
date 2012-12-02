      subroutine fibonacci2 ( n, a )

c*********************************************************************72
c
cc FIBONACCI2 returns the FIBONACCI2 matrix.
c
c  Example:
c
c    N = 5
c
c    0 1 0 0 0
c    1 1 0 0 0
c    0 1 1 0 0
c    0 0 1 1 0
c    0 0 0 1 1
c
c  Properties:
c
c    A is generally not symmetric: A' /= A.
c
c    A is tridiagonal.
c
c    Because A is tridiagonal, it has property A (bipartite).
c
c    A is banded, with bandwidth 3.
c
c    A is integral, therefore det ( A ) is integral, and 
c    det ( A ) * inverse ( A ) is integral.
c
c    A is a zero/one matrix.
c
c    If N = 1 then
c      det ( A ) = 0
c    else
c      det ( A ) = -1
c
c    If 1 < N, then A is unimodular.
c
c    When applied to a Fibonacci1 matrix B, the Fibonacci2 matrix
c    A produces the "next" Fibonacci1 matrix C = A*B.
c
c    Let PHI be the golden ratio (1+sqrt(5))/2.
c
c    For 2 <= N, the eigenvalues and eigenvectors are:
c
c    LAMBDA(1)     = PHI,     vector = (1,PHI,PHI^2,...PHI^(N-1));
c    LAMBDA(2:N-1) = 1        vector = (0,0,0,...,0,1);
c    LAMBDA(N)     = 1 - PHI. vector = ((-PHI)^(N-1),(-PHI)^(N-2),...,1)
c
c    Note that there is only one eigenvector corresponding to 1.
c    Hence, for 3 < N, the matrix is defective.  This fact means, 
c    for instance, that the convergence of the eigenvector in the power 
c    method will be very slow.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Output, double precision A(N,N), the matrix.
c
      implicit none

      integer n

      double precision a(n,n)
      integer i
      integer j

      do j = 1, n
        do i = 1, n

          if ( i == 1 ) then

            if ( j == 2 ) then
              a(i,j) = 1.0D+00
            else
              a(i,j) = 0.0D+00
            end if

          else

            if ( j == i - 1 .or. j == i ) then
              a(i,j) = 1.0D+00
            else
              a(i,j) = 0.0D+00
            end if

          end if

        end do
      end do

      return
      end
      subroutine power_method ( n, a, y, it_max, tol, lambda, it_num )

c*********************************************************************72
c
cc POWER_METHOD applies the power method for real eigenvalues.
c
c  Discussion:
c
c    For a given NxN matrix A and an N vector Y, the power method produces
c    a series of estimates for LAMBDA, the largest eigenvalue, and Y,
c    the eigenvector corresponding to LAMBDA.
c
c    The iteration repeats the following steps
c
c      AY     = A * Y
c      LAMBDA = || AY ||
c      Y      = AY / LAMBDA
c
c    If the matrix A has a single real eigenvalue of maximum modulus,
c    then this iteration will generally produce a good estimate for that
c    eigenvalue and its corresponding eigenvector.
c
c    If there are multiple distinct eigenvalues of the same modulus,
c    perhaps two values of opposite sign, or complex eigenvalues, then
c    the situation is more complicated.
c
c    Separate issues:
c
c    * when estimating the value of LAMBDA, we use the Rayleigh quotient,
c    LAMBDA = ( y' * A * y ) / ( y' * y ).  Since we normalize Y, the
c    bottom of the fraction is 1.  Using this estimate allows us to
c    easily capture the sign of LAMDBA.  Using the eucldean norm
c    instead, for instance, would always give a positive value.
c
c    * If the dominant eigenvalue is negative, then the iteration 
c    as given will produce eigenvector iterates that alternate in sign.  
c    
c    * It is worth knowing whether the successive eigenvector estimates
c    are tending to some value.  Since an eigenvector is really a direction,
c    we need to normalize the vectors, and we need to somehow treat both
c    a vector and its negative as holding the same information.  This
c    means that the proper way of measuring the difference between two
c    eigenvector estimates is to normalize them both, and then compute
c    the cosine between them as y1'y2, followed by the sine, which is
c    sqrt ( 1 - ( y1'y2)^2 ).  If this sine is small, the vectors y1 and y2
c    are "close" in the sense of direction.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 July 2008
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
c    Input/output, double precision Y(N), the estimate for the eigenvector.
c
c    Input, integer IT_MAX, the maximum number of iterations to take.
c    1 <= IT_MAX.
c
c    Input, double precision TOL, an error tolerance.
c
c    Output, double precision LAMBDA, the estimate for the eigenvalue.
c
c    Output, integer IT_NUM, the number of iterations taken.
c
      implicit none

      integer n

      double precision a(n,n)
      double precision ay(n)
      double precision cos_y1y2
      logical debug
      parameter ( debug = .false. )
      integer i
      integer it_max
      integer it_num
      integer j
      double precision lambda
      double precision lambda_old
      double precision norm
      double precision r8_epsilon
      double precision sin_y1y2
      double precision tol
      double precision val_dif
      double precision vec_dif
      double precision y(n)
      double precision y_old(n)

      if ( debug ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    '     IT      Lambda          Delta-Lambda    Delta-Y'
        write ( *, '(a)' ) ' '
      end if
c
c  Force Y to be a vector of unit norm.
c
      norm = 0.0
      do i = 1, n
        norm = norm + y(i) * y(i)
      end do
      norm = sqrt ( norm )

      do i = 1, n
        y(i) = y(i) / norm
      end do

      it_num = 0

      do i = 1, n
        y_old(i) = y(i)
      end do

      do i = 1, n
        ay(i) = 0.0D+00
        do j = 1, n
          ay(i) = ay(i) + a(i,j) * y(j)
        end do
      end do

      lambda = 0.0D+00
      do i = 1, n
        lambda = lambda + y(i) * ay(i)
      end do

      norm = 0.0D+00
      do i = 1, n
        norm = norm + ay(i) * ay(i)
      end do
      norm = sqrt ( norm )

      do i = 1, n
        y(i) = ay(i) / norm
      end do

      if ( lambda .lt. 0.0D+00 ) then
        do i = 1, n
          y(i) = - y(i)
        end do
      end if

      val_dif = 0.0D+00

      cos_y1y2 = 0.0D+00
      do i = 1, n
        cos_y1y2 = cos_y1y2 + y(i) * y_old(i)
      end do

      sin_y1y2 = sqrt ( ( 1.0D+00 - cos_y1y2 ) 
     &  * ( 1.0D+00 + cos_y1y2 ) )

      if ( debug ) then
        write ( *, '(2x,i5,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &    it_num, lambda, val_dif, sin_y1y2
      end if

      do it_num = 1, it_max

        lambda_old = lambda
        do i = 1, n
          y_old(i) = y(i)
        end do

        do i = 1, n
          ay(i) = 0.0D+00
          do j = 1, n
            ay(i) = ay(i) + a(i,j) * y(j)
          end do
        end do

        lambda = 0.0D+00
        do i = 1, n
          lambda = lambda + y(i) * ay(i)
        end do

        norm = 0.0D+00
        do i = 1, n
          norm = norm + ay(i) * ay(i)
        end do
        norm = sqrt ( norm )

        do i = 1, n
          y(i) = ay(i) / norm
        end do

        if ( lambda .lt. 0.0D+00 ) then
          do i = 1, n
            y(i) = - y(i)
          end do
        end if

        val_dif = abs ( lambda - lambda_old )

        cos_y1y2 = 0.0D+00
        do i = 1, n
          cos_y1y2 = cos_y1y2 + y(i) * y_old(i)
        end do

        sin_y1y2 = sqrt ( ( 1.0D+00 - cos_y1y2 ) 
     &    * ( 1.0D+00 + cos_y1y2 ) )

        if ( debug ) then
          write ( *, '(2x,i5,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &      it_num, lambda, val_dif, sin_y1y2
        end if

        if ( val_dif .le. tol ) then
          go to 10
        end if

      end do

10    continue

      do i = 1, n
        y(i) = ay(i) / lambda
      end do

      return
      end
      subroutine power_method2 ( n, a, x_init, it_max, tol, lambda, v,
     &  it_num )

c*********************************************************************72
c
cc POWER_METHOD2 applies the power method for possibly complex eigenvalues.
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
c  Reference:
c
c    Eric VanDeVelde,
c    Concurrent Scientific Programming,
c    Springer, 1994,
c    ISBN: 0-387-94195-9,
c    LC: QA76.58.V35.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision A(N,N), the matrix.
c
c    Input, double precision X_INIT(N), the initial estimate for the eigenvector.
c
c    Input, integer IT_MAX, the maximum number of iterations to take.
c    1 <= IT_MAX.
c
c    Input, double precision TOL, an error tolerance.
c
c    Output, double complex LAMBDA, the estimate for the eigenvalue.
c
c    Output, double complex V(N), the estimate for the eigenvector.
c
c    Output, integer IT_NUM, the number of iterations taken.
c
      implicit none

      integer n

      integer n_max
      parameter ( n_max = 100 )

      double precision a(n,n)
      double precision alpha
      double precision beta
      double precision gamma
      integer i
      integer it
      integer it_max
      integer it_num
      integer j
      double complex lambda
      double precision lambda_imag
      double precision lambda_real
      double precision pi_xx
      double precision pi_xy
      double precision pi_xz
      double precision pi_yy
      double precision pi_yz
      double precision pi_zz
      double precision tol
      double complex v(n)
      double precision x(n_max)
      double precision x_init(n_max)
      double precision y(n)
      double precision z(n_max)

      it_num = 0

      if ( n_max < n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'POWER_METHOD2 - Fatal error!'
        write ( *, '(a)' ) '  Maximum N_MAX exceeded by input N.'
        write ( *, '(a,i8)' ) '  Maximum N_MAX is ', n_max
        write ( *, '(a,i8)' ) '  Input N is       ', n
        write ( *, '(a)' ) '  Increase N_MAX, recompile, and try again.'
        stop
      end if
c
c  Compute data necessary to start the iteration.
c
      do i = 1, n
        x(i) = x_init(i)
      end do

      pi_xx = 0.0
      do i = 1, n
        pi_xx = pi_xx + x(i) * x(i)
      end do

      do i = 1, n
        x(i) = x(i) / pi_xx
      end do

      do i = 1, n
        y(i) = 0.0
        do j = 1, n
          y(i) = y(i) + a(i,j) * x(j)
        end do
      end do

      pi_xy = 0.0
      do i = 1, n
        pi_xy = pi_xy + x(i) * y(i)
      end do

      pi_yy = 0.0
      do i = 1, n
        pi_yy = pi_yy + y(i) * y(i)
      end do

      do it = 1, it_max

        if ( pi_yy - pi_xy * pi_xy .lt. tol * tol * pi_yy ) then
          lambda = pi_xy
          do i = 1, n
            v(i) = y(i) / sqrt ( pi_yy )
          end do
          return
        end if

        do i = 1, n
          z(i) = 0.0D+00
          do j = 1, n
            z(i) = z(i) + a(i,j) * y(j)
          end do
        end do

        pi_xz = 0.0D+00
        do i = 1, n
          pi_xz = pi_xz + x(i) * z(i)
        end do

        pi_yz = 0.0D+00
        do i = 1, n
          pi_yz = pi_yz + y(i) * z(i)
        end do

        pi_zz = 0.0D+00
        do i = 1, n
          pi_zz = pi_zz + z(i) * z(i)
        end do

        alpha = - ( pi_yz - pi_xy * pi_xz ) / ( pi_yy - pi_xy * pi_xy )
        beta = ( pi_xy * pi_yz - pi_yy * pi_xz ) 
     &    / ( pi_yy - pi_xy * pi_xy )
        gamma = pi_zz + alpha * alpha * pi_yy + beta * beta 
     &    + 2.0D+00 
     &    * ( alpha * pi_yz + beta * pi_xz + alpha * beta * pi_xy )

        if ( gamma .lt. tol * tol * pi_zz .and. 
     &       alpha * alpha .lt. 4.0D+00 * beta ) then

          lambda_real = - alpha / 2.0D+00
          lambda_imag = sqrt ( 4.0D+00 * beta - alpha * alpha ) 
     &      / 2.0D+00
          lambda = cmplx ( lambda_real, lambda_imag )

          do i = 1, n
            v(i) = ( lambda * y(i) - z(i) ) 
     &        / sqrt ( beta * pi_yy + alpha * pi_yz + pi_zz )
          end do

          return
        end if
     
        do i = 1, n
          x(i) = y(i) / sqrt ( pi_yy )
        end do

        do i = 1, n
          y(i) = z(i) / sqrt ( pi_yy )
        end do

        pi_xy = pi_yz / pi_yy
        pi_yy = pi_zz / pi_yy

        it_num = it

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POWER_METHOD2 - Fatal error!'
      write ( *, '(a)' ) '  Convergence was not reached.'

      stop
      end
      function r8_epsilon ( )

c*********************************************************************72
c
cc R8_EPSILON returns the R8 roundoff unit.
c
c  Discussion:
c
c    The roundoff unit is a number R which is a power of 2 with the
c    property that, to the precision of the computer's arithmetic,
c      1 .lt. 1 + R
c    but
c      1 = ( 1 + R / 2 )
c
c    FORTRAN90 provides the superior library routine
c
c      EPSILON ( X )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 March 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision R8_EPSILON, the R8 roundoff unit.
c
      implicit none

      double precision r8
      double precision r8_epsilon
      double precision r8_test

      r8 = 1.0D+00
      r8_test = 1.0D+00 + ( r8 / 2.0D+00 )

10    continue

      if ( 1.0D+00 .lt. r8_test ) then
        r8 = r8 / 2.0D+00
        r8_test = 1.0D+00 + ( r8 / 2.0D+00 )
        go to 10
      end if

      r8_epsilon = r8

      return
      end
      subroutine r8vec_uniform_01 ( n, seed, r )

c*********************************************************************72
c
cc R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 August 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Second Edition,
c    Springer, 1987,
c    ISBN: 0387964673,
c    LC: QA76.9.C65.B73.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, December 1986, pages 362-376.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley, 1998,
c    ISBN: 0471134031,
c    LC: T57.62.H37.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, Number 2, 1969, pages 136-143.
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
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer k
      integer seed
      double precision r(n)

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8VEC_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + i4_huge
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
      subroutine tris ( m, n, x, y, z, a )

c*********************************************************************72
c
cc TRIS returns the TRIS matrix.
c
c  Discussion:
c
c    The matrix is a tridiagonal matrix defined by three scalars.
c
c    See page 155 of the Todd reference.
c
c  Formula:
c
c    if ( J = I-1 )
c      A(I,J) = X
c    else if ( J = I )
c      A(I,J) = Y
c    else if ( J = I + 1 )
c      A(I,J) = Z
c    else
c      A(I,J) = 0
c
c  Example:
c
c    M = 5, N = 5, X = 1, Y = 2, Z = 3
c
c    2 3 0 0 0
c    1 2 3 0 0
c    0 1 2 3 0
c    0 0 1 2 3
c    0 0 0 1 2
c
c  Properties:
c
c    A is generally not symmetric: A' /= A.
c
c    A is tridiagonal.
c
c    Because A is tridiagonal, it has property A (bipartite).
c
c    A is banded, with bandwidth 3.
c
c    A is Toeplitz: constant along diagonals.
c
c    If Y is not zero, then for A to be singular, it must be the case that
c
c      0.5 * Y / sqrt ( X * Z ) < 1
c
c    and
c
c      cos (K*PI/(N+1)) = - 0.5 * Y / sqrt ( X * Z ) for some 1 <= K <= N.
c
c    If Y is zero, then A is singular when N is odd, or if X or Z is zero.
c
c    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
c
c    A has eigenvalues
c
c      LAMBDA(I) = Y + 2 * sqrt(X*Z) * COS(I*PI/(N+1))
c
c    The eigenvalues will be complex if X * Z < 0.
c
c    If X = Z, the matrix is symmetric.
c
c    As long as X and Z are nonzero, the matrix is irreducible.
c
c    If X = Z = -1, and Y = 2, the matrix is a symmetric, positive
c    definite M matrix, the negative of the second difference matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    John Todd,
c    Basic Numerical Mathematics,
c    Volume 2: Numerical Algebra,
c    Birkhauser, 1980,
c    ISBN: 0817608117,
c    LC: QA297.T58.
c
c  Parameters:
c
c    Input, integer M, N, the order of the matrix.
c
c    Input, double precision X, Y, Z, the scalars that define A.
c
c    Output, double precision A(M,N), the matrix.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer i
      integer j
      double precision x
      double precision y
      double precision z

      do j = 1, n
        do i = 1, m

          if ( j == i - 1 ) then
            a(i,j) = x
          else if ( j == i ) then
            a(i,j) = y
          else if ( j == i + 1 ) then
            a(i,j) = z
          else
            a(i,j) = 0.0D+00
          end if

        end do
      end do

      return
      end
      subroutine tris_eigenvalues ( n, x, y, z, lambda )

c*********************************************************************72
c
cc TRIS_EIGENVALUES returns the eigenvalues of the TRIS matrix.
c
c  Discussion:
c
c    The eigenvalues will be complex if X * Z < 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision X, Y, Z, the scalars that define A.
c
c    Output, double complex LAMBDA(N), the eigenvalues.
c
      implicit none

      integer n

      double precision angle
      double complex arg
      integer i
      double complex lambda(n)
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision x
      double precision y
      double precision z

      do i = 1, n
        angle = dble ( i ) * pi / dble ( n + 1 )
        arg = cmplx ( x * z, 0.0D+00 )
        lambda(i) = y + 2.0D+00 * sqrt ( arg ) * cos ( angle )
      end do

      return
      end
