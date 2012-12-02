      subroutine combin ( alpha, beta, n, a )

c*********************************************************************72
c
cc COMBIN returns the COMBIN matrix.
c
c  Discussion:
c
c    This matrix is known as the combinatorial matrix.
c
c  Formula:
c
c    If ( I = J ) then
c      A(I,J) = ALPHA + BETA
c    else
c      A(I,J) = BETA
c
c  Example:
c
c    N = 5, ALPHA = 2, BETA = 3
c
c    5 3 3 3 3
c    3 5 3 3 3
c    3 3 5 3 3
c    3 3 3 5 3
c    3 3 3 3 5
c
c  Properties:
c
c    A is symmetric: A' = A.
c
c    Because A is symmetric, it is normal.
c
c    Because A is normal, it is diagonalizable.
c
c    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
c
c    A is a circulant matrix: each row is shifted once to get the next row.
c
c    det ( A ) = ALPHA^(N-1) * ( ALPHA + N * BETA ).
c
c    A has constant row sums.
c
c    Because A has constant row sums,
c    it has an eigenvalue with this value,
c    and a (right) eigenvector of ( 1, 1, 1, ..., 1 ).
c
c    A has constant column sums.
c
c    Because A has constant column sums,
c    it has an eigenvalue with this value,
c    and a (left) eigenvector of ( 1, 1, 1, ..., 1 ).
c
c    LAMBDA(1:N-1) = ALPHA,
c    LAMBDA(N) = ALPHA + N * BETA.
c
c    The eigenvector associated with LAMBDA(N) is (1,1,1,...,1)/sqrt(N).
c
c    The other N-1 eigenvectors are simply any (orthonormal) basis
c    for the space perpendicular to (1,1,1,...,1).
c
c    A is nonsingular if ALPHA /= 0 and ALPHA + N * BETA /= 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Gregory, David Karney,
c    A Collection of Matrices for Testing Computational Algorithms,
c    Wiley, 1969,
c    ISBN: 0882756494,
c    LC: QA263.68
c
c    Donald Knuth,
c    The Art of Computer Programming,
c    Volume 1, Fundamental Algorithms, Second Edition,
c    Addison-Wesley, Reading, Massachusetts, 1973, page 36.
c
c  Parameters:
c
c    Input, double precision ALPHA, BETA, scalars that define A.
c
c    Input, integer N, the order of the matrix.
c
c    Output, double precision A(N,N), the matrix.
c
      implicit none

      integer n

      double precision a(n,n)
      double precision alpha
      double precision beta
      integer i
      integer j

      do j = 1, n
        do i = 1, n
          a(i,j) = beta
        end do
      end do

      do i = 1, n
        a(i,i) = a(i,i) + alpha
      end do

      return
      end
      subroutine combin_inverse ( alpha, beta, n, a )

c*********************************************************************72
c
cc COMBIN_INVERSE returns the inverse of the COMBIN matrix.
c
c  Formula:
c
c    if ( I = J )
c      A(I,J) = (ALPHA+(N-1)*BETA) / (ALPHA*(ALPHA+N*BETA))
c    else
c      A(I,J) =             - BETA / (ALPHA*(ALPHA+N*BETA))
c
c  Example:
c
c    N = 5, ALPHA = 2, BETA = 3
c
c           14 -3 -3 -3 -3
c           -3 14 -3 -3 -3
c   1/34 *  -3 -3 14 -3 -3
c           -3 -3 -3 14 -3
c           -3 -3 -3 -3 14
c
c  Properties:
c
c    A is symmetric: A' = A.
c
c    Because A is symmetric, it is normal.
c
c    Because A is normal, it is diagonalizable.
c
c    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
c
c    A is a circulant matrix: each row is shifted once to get the next row.
c
c    A is Toeplitz: constant along diagonals.
c
c    det ( A ) = 1 / (ALPHA^(N-1) * (ALPHA+N*BETA)).
c
c    A is well defined if ALPHA /= 0D+00 and ALPHA+N*BETA /= 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Knuth,
c    The Art of Computer Programming,
c    Volume 1, Fundamental Algorithms, Second Edition,
c    Addison-Wesley, Reading, Massachusetts, 1973, page 36.
c
c  Parameters:
c
c    Input, double precision ALPHA, BETA, scalars that define the matrix.
c
c    Input, integer N, the order of the matrix.
c
c    Output, double precision A(N,N), the matrix.
c
      implicit none

      integer n

      double precision a(n,n)
      double precision alpha
      double precision beta
      double precision bot
      integer i
      integer j

      if ( alpha .eq. 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'COMBIN_INVERSE - Fatal error!'
        write ( *, '(a)' ) '  The entries of the matrix are undefined'
        write ( *, '(a)' ) '  because ALPHA = 0.'
        stop
      else if ( alpha + n * beta .eq. 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'COMBIN_INVERSE - Fatal error!'
        write ( *, '(a)' ) '  The entries of the matrix are undefined'
        write ( *, '(a)' ) '  because ALPHA+N*BETA is zero.'
        stop
      end if

      bot = alpha * ( alpha + dble ( n ) * beta )

      do j = 1, n
        do i = 1, n

          if ( i .eq. j ) then
            a(i,j) = ( alpha + dble ( n - 1 ) * beta ) / bot
          else
            a(i,j) = - beta / bot
          end if

        end do
      end do

      return
      end
      subroutine condition_hager ( n, a, cond )

c*********************************************************************72
c
cc CONDITION_HAGER estimates the L1 condition number of a matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 April 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    William Hager,
c    Condition Estimates,
c    SIAM Journal on Scientific and Statistical Computing,
c    Volume 5, Number 2, June 1984, pages 311-316.
c
c  Parameters:
c
c    Input, integer N, the dimension of the matrix.
c
c    Input, double precision A(N,N), the matrix.
c
c    Output, double precision COND, the estimated L1 condition.
c
      implicit none

      integer n

      double precision a(n,n)
      double precision a_lu(n,n)
      double precision b(n)
      double precision c1
      double precision c2
      double precision cond
      integer i
      integer i1
      integer i2
      integer info
      integer j
      integer job
      integer pivot(n)
      double precision r8_sign
      double precision r8mat_norm_l1

      i1 = -1
      c1 = 0.0D+00
c
c  Factor the matrix.
c
      do j = 1, n
        do i = 1, n
          a_lu(i,j) = a(i,j)
        end do
      end do

      call r8ge_fa ( n, a_lu, pivot, info )

      do i = 1, n
        b(i) = 1.0D+00 / dble ( n )
      end do

10    continue

        job = 0
        call r8ge_sl ( n, a_lu, pivot, b, job )

        c2 = 0.0D+00
        do i = 1, n
          c2 = c2 + abs ( b(i) )
        end do

        do i = 1, n
          b(i) = r8_sign ( b(i) )
        end do

        job = 1
        call r8ge_sl ( n, a_lu, pivot, b, job )

        call r8vec_max_abs_index ( n, b, i2 )

        if ( 1 .le. i1 ) then
          if ( i1 .eq. i2 .or. c2 .le. c1 ) then
            go to 20
          end if
        end if

        i1 = i2
        c1 = c2

        b(1:n) = 0.0D+00
        b(i1) = 1.0D+00

      go to 10

20    continue

      cond = c2 * r8mat_norm_l1 ( n, n, a )

      return
      end
      subroutine condition_linpack ( n, a, pivot, cond, z )

c*********************************************************************72
c
cc CONDITION_LINPACK estimates the L1 condition number of a matrix.
c
c  Discussion:
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c    R8GE storage is used by LINPACK and LAPACK.
c
c    For the system A * X = B, relative perturbations in A and B
c    of size EPSILON may cause relative perturbations in X of size
c    EPSILON * COND.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 March 2004
c
c  Author:
c
c    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix A.
c
c    Input/output, double precision A(N,N).  On input, a matrix to be factored.
c    On output, the LU factorization of the matrix.
c
c    Output, integer PIVOT(N), the pivot indices.
c
c    Output, double precision COND, the estimated L1 condition.
c
c    Output, double precision Z(N), a work vector whose contents are 
c    usually unimportant.  If A is close to a singular matrix, then Z is
c    an approximate null vector in the sense that
c      norm ( A * Z ) = RCOND * norm ( A ) * norm ( Z ).
c
      implicit none

      integer n

      double precision a(n,n)
      double precision anorm
      double precision bnorm
      double precision cond
      double precision ek
      integer i
      integer info
      integer j
      integer k
      integer l
      integer pivot(n)
      double precision r8vec_norm_l1
      double precision s
      double precision sm
      double precision t
      double precision wk
      double precision wkm
      double precision ynorm
      double precision z(n)
c
c  Compute the L1 norm of A.
c
      anorm = 0.0D+00
      do j = 1, n
        bnorm = 0.0D+00
        do i = 1, n
          bnorm = bnorm + abs ( a(i,j) )
        end do
        anorm = max ( anorm, bnorm )
      end do
c
c  Compute the LU factorization.
c
      call r8ge_fa ( n, a, pivot, info )
c
c  COND = norm(A) * (estimate of norm(inverse(A)))
c
c  estimate of norm(inverse(A)) = norm(Z) / norm(Y)
c
c  where
c    A * Z = Y
c  and
c    A' * Y = E
c
c  The components of E are chosen to cause maximum local growth in the
c  elements of W, where U'*W = E.  The vectors are frequently rescaled
c  to avoid overflow.
c
c  Solve U' * W = E.
c
      ek = 1.0D+00

      do i = 1, n
        z(i) = 0.0D+00
      end do

      do k = 1, n

        if ( z(k) .ne. 0.0D+00 ) then
          ek = sign ( ek, -z(k) )
        end if

        if ( abs ( a(k,k) ) .lt. abs ( ek - z(k) ) ) then
          s = abs ( a(k,k) ) / abs ( ek - z(k) )
          do i = 1, n
            z(i) = s * z(i)
          end do
          ek = s * ek
        end if

        wk = ek - z(k)
        wkm = -ek - z(k)
        s = abs ( wk )
        sm = abs ( wkm )

        if ( a(k,k) .ne. 0.0D+00 ) then
          wk = wk / a(k,k)
          wkm = wkm / a(k,k)
        else
          wk = 1.0D+00
          wkm = 1.0D+00
        end if

        if ( k + 1 .le. n ) then

          do j = k + 1, n
            sm = sm + abs ( z(j) + wkm * a(k,j) )
            z(j) = z(j) + wk * a(k,j)
            s = s + abs ( z(j) )
          end do

          if ( s .lt. sm ) then
            t = wkm - wk
            wk = wkm
            do i = k + 1, n
              z(i) = z(i) + t * a(k,i)
            end do
          end if

        end if

        z(k) = wk

      end do

      t = r8vec_norm_l1 ( n, z )
      do i = 1, n
        z(i) = z(i) / t
      end do
c
c  Solve L' * Y = W
c
      do k = n, 1, -1

        t = 0.0D+00
        do i = k + 1, n
          t = t + a(i,k) * z(i)
        end do

        z(k) = z(k) + t

        t = abs ( z(k) )

        if ( 1.0D+00 .lt. t ) then
          do i = 1, n
            z(i) = z(i) / t
          end do
        end if

        l = pivot(k)

        t    = z(l)
        z(l) = z(k)
        z(k) = t

      end do

      t = r8vec_norm_l1 ( n, z )
      do i = 1, n
        z(i) = z(i) / t
      end do

      ynorm = 1.0D+00
c
c  Solve L * V = Y.
c
      do k = 1, n

        l = pivot(k)

        t    = z(l)
        z(l) = z(k)
        z(k) = t

        do i = k + 1, n
          z(i) = z(i) + z(k) * a(i,k)
        end do

        if ( 1.0D+00 .lt. abs ( z(k) ) ) then
          ynorm = ynorm / abs ( z(k) )
          do i = 1, n
            z(i) = z(i) / abs ( z(k) )
          end do
        end if

      end do

      s = r8vec_norm_l1 ( n, z )
      do i = 1, n
        z(i) = z(i) / s
      end do
      ynorm = ynorm / s
c
c  Solve U * Z = V.
c
      do k = n, 1, -1

        if ( abs ( a(k,k) ) .lt. abs ( z(k) ) ) then
          s = abs ( a(k,k) ) / abs ( z(k) )
          do i = 1, n
            z(i) = s * z(i)
          end do
          ynorm = s * ynorm
        end if

        if ( a(k,k) .ne. 0.0D+00 ) then
          z(k) = z(k) / a(k,k)
        else
          z(k) = 1.0D+00
        end if

        do i = 1, k - 1
          z(i) = z(i) - z(k) * a(i,k)
        end do

      end do
c
c  Normalize Z in the L1 norm.
c
      s = 1.0D+00 / r8vec_norm_l1 ( n, z )
      do i = 1, n
        z(i) = s * z(i)
      end do
      ynorm = s * ynorm

      cond = anorm / ynorm

      return
      end
      subroutine condition_sample1 ( n, a, m, cond )

c*********************************************************************72
c
cc CONDITION_SAMPLE1 estimates the L1 condition number of a matrix.
c
c  Discussion:
c
c    A naive sampling method is used.
c
c    Only "forward" sampling is used, that is, we only look at results
c    of the form y=A*x.
c
c    Presumably, solving systems A*y=x would give us a better idea of 
c    the inverse matrix.
c
c    Moreover, a power sequence y1 = A*x, y2 = A*y1, ... and the same for
c    the inverse might work better too.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 April 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the matrix.
c
c    Input, double precision A(N,N), the matrix.
c
c    Input, integer M, the number of samples to use.
c
c    Output, double precision COND, the estimated L1 condition.
c
      implicit none

      integer n

      double precision a(n,n)
      double precision a_norm
      double precision ainv_norm
      double precision ax(n)
      double precision ax_norm
      double precision cond
      integer i
      integer m
      double precision r8vec_norm_l1
      integer seed
      double precision x(n)
      double precision x_norm

      a_norm = 0.0D+00
      ainv_norm = 0.0D+00
      seed = 123456789

      do i = 1, m

        call r8vec_uniform_unit ( n, seed, x )
        x_norm = r8vec_norm_l1 ( n, x )
        ax = matmul ( a, x )
        ax_norm = r8vec_norm_l1 ( n, ax )

        if ( ax_norm .eq. 0.0D+00 ) then
          cond = 0.0D+00
          return
        end if

        a_norm    = max ( a_norm,    ax_norm / x_norm  )
        ainv_norm = max ( ainv_norm, x_norm  / ax_norm )

      end do

      cond = a_norm * ainv_norm

      return
      end
      subroutine conex1 ( alpha, a )

c*********************************************************************72
c
cc CONEX1 returns the CONEX1 matrix.
c
c  Discussion:
c
c    The CONEX1 matrix is a counterexample to the LINPACK condition
c    number estimator RCOND available in the LINPACK routine DGECO.
c
c  Formula:
c
c    1  -1 -2*ALPHA   0
c    0   1    ALPHA    -ALPHA
c    0   1  1+ALPHA  -1-ALPHA
c    0   0  0           ALPHA
c
c  Example:
c
c    ALPHA = 100
c
c    1  -1  -200     0
c    0   1   100  -100
c    0   1   101  -101
c    0   0     0   100
c
c  Properties:
c
c    A is generally not symmetric: A' /= A.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Alan Cline, RK Rew,
c    A set of counterexamples to three condition number estimators,
c    SIAM Journal on Scientific and Statistical Computing,
c    Volume 4, 1983, pages 602-611.
c
c  Parameters:
c
c    Input, double precision ALPHA, the scalar defining A.  
c    A common value is 100.0.
c
c    Output, double precision A(4,4), the matrix.
c
      implicit none

      double precision a(4,4)
      double precision alpha

      a(1,1) = 1.0D+00
      a(2,1) = 0.0D+00
      a(3,1) = 0.0D+00
      a(4,1) = 0.0D+00

      a(1,2) = -1.0D+00
      a(2,2) = 1.0D+00
      a(3,2) = 1.0D+00
      a(4,2) = 0.0D+00

      a(1,3) = -2.0D+00 * alpha
      a(2,3) = alpha
      a(3,3) = 1.0D+00 + alpha
      a(4,3) = 0.0D+00

      a(1,4) = 0.0D+00
      a(2,4) = -alpha
      a(3,4) = -1.0D+00 - alpha
      a(4,4) = alpha

      return
      end
      subroutine conex1_inverse ( alpha, a )

c*********************************************************************72
c
cc CONEX1_INVERSE returns the inverse of the CONEX1 matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ALPHA, the scalar defining A.  
c
c    Output, double precision A(4,4), the matrix.
c
      implicit none

      double precision    a(4,4)
      double precision    alpha

      a(1,1) =  1.0D+00
      a(1,2) =  1.0D+00 - alpha
      a(1,3) =            alpha
      a(1,4) =  2.0D+00

      a(2,1) =  0.0D+00
      a(2,2) =  1.0D+00 + alpha
      a(2,3) =          - alpha
      a(2,4) =  0.0D+00

      a(3,1) =  0.0D+00
      a(3,2) = -1.0D+00
      a(3,3) =  1.0D+00
      a(3,4) =  1.0D+00 / alpha

      a(4,1) = 0.0D+00
      a(4,2) = 0.0D+00
      a(4,3) = 0.0D+00
      a(4,4) = 1.0D+00 / alpha

      return
      end
      subroutine conex2 ( alpha, a )

c*********************************************************************72
c
cc CONEX2 returns the CONEX2 matrix.
c
c  Formula:
c
c    1   1-1/ALPHA^2  -2
c    0   1/ALPHA      -1/ALPHA
c    0   0             1
c
c  Example:
c
c    ALPHA = 100
c
c    1  0.9999  -2
c    0  0.01    -0.01
c    0  0        1
c
c  Properties:
c
c    A is generally not symmetric: A' /= A.
c
c    A is upper triangular.
c
c    det ( A ) = 1 / ALPHA.
c
c    LAMBDA = ( 1, 1/ALPHA, 1 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Alan Cline, RK Rew,
c    A set of counterexamples to three condition number estimators,
c    SIAM Journal on Scientific and Statistical Computing,
c    Volume 4, 1983, pages 602-611.
c
c  Parameters:
c
c    Input, double precision ALPHA, the scalar defining A.  
c    A common value is 100.0.  ALPHA must not be zero.
c
c    Output, double precision A(3,3), the matrix.
c
      implicit none

      double precision a(3,3)
      double precision alpha

      if ( alpha .eq. 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CONEX2 - Fatal error!'
        write ( *, '(a)' ) '  The input value of ALPHA was zero.'
        stop
      end if

      a(1,1) = 1.0D+00
      a(1,2) = ( alpha**2 - 1.0D+00 ) / alpha**2
      a(1,3) = -2.0D+00

      a(2,1) = 0.0D+00
      a(2,2) = 1.0D+00 / alpha
      a(2,3) = -1.0D+00 / alpha

      a(3,1) = 0.0D+00
      a(3,2) = 0.0D+00
      a(3,3) = 1.0D+00

      return
      end
      subroutine conex2_inverse ( alpha, a )

c*********************************************************************72
c
cc CONEX2_INVERSE returns the inverse of the CONEX2 matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ALPHA, the scalar defining A.  
c    A common value is 100.0.  ALPHA must not be zero.
c
c    Output, double precision A(3,3), the matrix.
c
      implicit none

      double precision a(3,3)
      double precision alpha

      if ( alpha .eq. 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CONEX2_INVERSE - Fatal error!'
        write ( *, '(a)' ) '  The input value of ALPHA was zero.'
        stop
      end if

      a(1,1) = 1.0D+00
      a(1,2) = ( 1.0D+00 - alpha**2 ) / alpha
      a(1,3) = ( 1.0D+00 + alpha**2 ) / alpha**2

      a(2,1) = 0.0D+00
      a(2,2) = alpha
      a(2,3) = 1.0D+00

      a(3,1) = 0.0D+00
      a(3,2) = 0.0D+00
      a(3,3) = 1.0D+00

      return
      end
      subroutine conex3 ( n, a )

c*********************************************************************72
c
cc CONEX3 returns the CONEX3 matrix.
c
c  Formula:
c
c    if ( I = J and I < N )
c      A(I,J) =  1.0D+00 for 1<=I<N
c    else if ( I = J = N )
c      A(I,J) = -1.0D+00
c    else if ( J < I )
c      A(I,J) = -1.0D+00
c    else
c      A(I,J) =  0.0D+00
c
c  Example:
c
c    N = 5
c
c     1  0  0  0  0
c    -1  1  0  0  0
c    -1 -1  1  0  0
c    -1 -1 -1  1  0
c    -1 -1 -1 -1 -1
c
c  Properties:
c
c    A is generally not symmetric: A' /= A.
c
c    A is integral, therefore det ( A ) is integral, and 
c    det ( A ) * inverse ( A ) is integral.
c
c    A is lower triangular.
c
c    det ( A ) = -1.
c
c    A is unimodular.
c
c    LAMBDA = ( 1, 1, 1, 1, -1 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Alan Cline, RK Rew,
c    A set of counterexamples to three condition number estimators,
c    SIAM Journal on Scientific and Statistical Computing,
c    Volume 4, 1983, pages 602-611.
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

          if ( j .lt. i ) then
            a(i,j) = -1.0D+00
          else if ( j .eq. i .and. i .ne. n ) then
            a(i,j) = 1.0D+00
          else if ( j .eq. i .and. i .eq. n ) then
            a(i,j) = - 1.0D+00
          else
            a(i,j) = 0.0D+00
          end if

        end do
      end do

      return
      end
      subroutine conex3_inverse ( n, a )

c*********************************************************************72
c
cc CONEX3_INVERSE returns the inverse of the CONEX3 matrix.
c
c  Example:
c
c    N = 5
c
c     1  0  0  0  0
c     1  1  0  0  0
c     2  1  1  0  0
c     4  2  1  1  0
c    -8 -4 -2 -1 -1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Alan Cline, RK Rew,
c    A set of counterexamples to three condition number estimators,
c    SIAM Journal on Scientific and Statistical Computing,
c    Volume 4, 1983, pages 602-611.
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

          if ( i .lt. n ) then
          
            if ( j .lt. i ) then
              a(i,j) = 2.0D+00**( i - j - 1 )
            else if ( i .eq. j ) then
              a(i,j) = 1.0D+00
            else
              a(i,j) = 0.0D+00
            end if
	
          else if ( i .eq. n ) then
          
            if ( j .lt. i ) then
              a(i,j) = - 2.0D+00**( i - j - 1 )
            else
              a(i,j) = -1.0D+00
            end if
	
          end if
          
        end do
      end do

      return
      end
      subroutine conex4 ( a )

c*********************************************************************72
c
cc CONEX4 returns the CONEX4 matrix.
c
c  Discussion:
c
c    7  10   8   7
c    6   8  10   9
c    5   7   9  10
c    5   7   6   5
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 April 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A(4,4), the matrix.
c
      implicit none

      double precision a(4,4)
      double precision a_save(4,4)
      integer i
      integer j

      save a_save

      data a_save /
     &   7.0D+00,  6.0D+00,  5.0D+00,  5.0D+00, 
     &  10.0D+00,  8.0D+00,  7.0D+00,  7.0D+00, 
     &   8.0D+00, 10.0D+00,  9.0D+00,  6.0D+00, 
     &   7.0D+00,  9.0D+00, 10.0D+00,  5.0D+00 /

      do j = 1, 4
        do i = 1, 4
          a(i,j) = a_save(i,j)
        end do
      end do

      return
      end
      subroutine conex4_inverse ( a )

c*********************************************************************72
c
cc CONEX4_INVERSE returns the inverse CONEX4 matrix.
c
c  Discussion:
c
c   -41  -17   10   68
c    25   10   -6  -41
c    10    5   -3  -17
c    -6   -3    2   10
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 April 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A(4,4), the matrix.
c
      implicit none

      double precision a(4,4)
      double precision a_save(4,4)
      integer i
      integer j

      save a_save

      data a_save /
     & -41.0D+00,  25.0D+00,  10.0D+00, -6.0D+00, 
     & -17.0D+00,  10.0D+00,   5.0D+00, -3.0D+00, 
     &  10.0D+00,  -6.0D+00,  -3.0D+00,  2.0D+00, 
     &  68.0D+00, -41.0D+00, -17.0D+00, 10.0D+00 /

      do j = 1, 4
        do i = 1, 4
          a(i,j) = a_save(i,j)
        end do
      end do

      return
      end
      subroutine kahan ( alpha, m, n, a )

c*********************************************************************72
c
cc KAHAN returns the KAHAN matrix.
c
c  Formula:
c
c    if ( I = J )
c      A(I,I) =  sin(ALPHA)^I
c    else if ( I < J )
c      A(I,J) = - sin(ALPHA)^I * cos(ALPHA)
c    else
c      A(I,J) = 0
c
c  Example:
c
c    ALPHA = 0.25, N = 4
c
c    S  -C*S    -C*S      -C*S
c    0     S^2  -C*S^2    -C*S^2
c    0     0       S^3    -C*S^3
c    0     0       0         S^4
c
c    where
c
c      S = sin(ALPHA), C=COS(ALPHA)
c
c  Properties:
c
c    A is upper triangular.
c
c    A = B * C, where B is a diagonal matrix and C is unit upper triangular.
c    For instance, for the case M = 3, N = 4:
c
c    A = | S 0    0    |  * | 1 -C -C  -C |
c        | 0 S^2  0    |    | 0  1 -C  -C |
c        | 0 0    S^3  |    | 0  0  1  -C |
c
c    A is generally not symmetric: A' /= A.
c
c    A has some interesting properties regarding estimation of
c    condition and rank.
c
c    det ( A ) = sin(ALPHA)^(N*(N+1)/2).
c
c    LAMBDA(I) = sin ( ALPHA )^I
c
c    A is nonsingular if and only if sin ( ALPHA ) =/= 0.
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
c  Reference:
c
c    Nicholas Higham,
c    A survey of condition number estimation for triangular matrices,
c    SIAM Review,
c    Volume 9, 1987, pages 575-596.
c
c    W Kahan,
c    Numerical Linear Algebra,
c    Canadian Mathematical Bulletin,
c    Volume 9, 1966, pages 757-801.
c
c  Parameters:
c
c    Input, double precision ALPHA, the scalar that defines A.  A typical
c    value is 1.2.  The "interesting" range of ALPHA is 0 < ALPHA < PI.
c
c    Input, integer M, N, the order of the matrix.
c
c    Output, double precision A(M,N), the matrix.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision alpha
      double precision csi
      integer i
      integer j
      double precision si

      do i = 1, m

        si = sin ( alpha )**i
        csi = - cos ( alpha ) * si

        do j = 1, n

          if ( j .lt. i ) then
            a(i,j) = 0.0D+00
          else if ( j .eq. i ) then
            a(i,j) = si
          else
            a(i,j) = csi
          end if

        end do
      end do

      return
      end
      subroutine kahan_inverse ( alpha, n, a )

c*********************************************************************72
c
cc KAHAN_INVERSE returns the inverse of the KAHAN matrix.
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
c    Input, double precision ALPHA, the scalar that defines A.  A typical 
c    value is 1.2.  The "interesting" range of ALPHA is 0 < ALPHA < PI.
c
c    Input, integer N, the order of the matrix.
c
c    Output, double precision A(N,N), the matrix.
c
      implicit none

      integer n

      double precision a(n,n)
      double precision alpha
      double precision ci
      integer i
      integer j
      double precision si

      ci = cos ( alpha )

      do i = 1, n
        do j = 1, n

          if ( i .eq. j ) then
            a(i,j) = 1.0D+00
          else if ( i .eq. j - 1 ) then
            a(i,j) = ci
          else if ( i .lt. j ) then
            a(i,j) = ci * ( 1.0D+00 + ci )**( j - i - 1 )
          else
            a(i,j) = 0.0D+00
          end if

        end do
      end do
c
c  Scale the columns.
c
      do j = 1, n
        si = sin ( alpha )**j
        do i = 1, n
          a(i,j) = a(i,j) / si
        end do
      end do

      return
      end
      subroutine r8ge_fa ( n, a, pivot, info )

c*********************************************************************72
c
cc R8GE_FA performs a LINPACK style PLU factorization of an R8GE matrix.
c
c  Discussion:
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c    R8GE storage is used by LINPACK and LAPACK.
c
c    R8GE_FA is a simplified version of the LINPACK routine SGEFA.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 March 2009
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input/output, double precision A(N,N), the matrix to be factored.
c    On output, A contains an upper triangular matrix and the multipliers
c    which were used to obtain it.  The factorization can be written
c    A = L * U, where L is a product of permutation and unit lower
c    triangular matrices and U is upper triangular.
c
c    Output, integer PIVOT(N), a vector of pivot indices.
c
c    Output, integer INFO, singularity flag.
c    0, no singularity detected.
c    nonzero, the factorization failed on the INFO-th step.
c 
      implicit none

      integer n

      double precision a(n,n)
      integer i
      integer info
      integer pivot(n)
      integer j
      integer k
      integer l
      double precision t

      info = 0

      do k = 1, n - 1
c
c  Find L, the index of the pivot row.
c
        l = k
        do i = k + 1, n
          if ( abs ( a(l,k) ) .lt. abs ( a(i,k) ) ) then
            l = i
          end if
        end do

        pivot(k) = l
c
c  If the pivot index is zero, the algorithm has failed.
c
        if ( a(l,k) .eq. 0.0D+00 ) then
          info = k
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8GE_FA - Fatal error!'
          write ( *, '(a,i8)' ) '  Zero pivot on step ', info
          return
        end if
c
c  Interchange rows L and K if necessary.
c
        if ( l .ne. k ) then
          t      = a(l,k)
          a(l,k) = a(k,k)
          a(k,k) = t
        end if
c
c  Normalize the values that lie below the pivot entry A(K,K).
c
        do i = k + 1, n
          a(i,k) = - a(i,k) / a(k,k)
        end do
c
c  Row elimination with column indexing.
c
        do j = k + 1, n

          if ( l .ne. k ) then
            t      = a(l,j)
            a(l,j) = a(k,j)
            a(k,j) = t
          end if

          do i = k + 1, n
            a(i,j) = a(i,j) + a(i,k) * a(k,j)
          end do

        end do

      end do

      pivot(n) = n

      if ( a(n,n) .eq. 0.0D+00 ) then
        info = n
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8GE_FA - Fatal error!'
        write ( *, '(a,i8)' ) '  Zero pivot on step ', info
      end if

      return
      end
      subroutine r8ge_inverse ( n, a, pivot )

c*********************************************************************72
c
cc R8GE_INVERSE computes the inverse of a matrix factored by R8GE_FA.
c
c  Discussion:
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c    R8GE storage is used by LINPACK and LAPACK.
c
c    R8GE_INVERSE is a simplified standalone version of the LINPACK routine
c    SGEDI.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix A.
c
c    Input/output, double precision A(N,N).
c    On input, the factor information computed by R8GE_FA.
c    On output, the inverse matrix.
c
c    Input, integer PIVOT(N), the pivot vector from R8GE_FA.
c
      implicit none

      integer n

      double precision a(n,n)
      integer i
      integer pivot(n)
      integer j
      integer k
      double precision t
      double precision work(n)
c
c  Compute Inverse(U).
c
      do k = 1, n

        a(k,k) = 1.0D+00 / a(k,k)
        do i = 1, k - 1
          a(i,k) = -a(i,k) * a(k,k)
        end do

        do j = k + 1, n

          t = a(k,j)
          a(k,j) = 0.0D+00
          do i = 1, k
            a(i,j) = a(i,j) + a(i,k) * t
          end do

        end do

      end do
c
c  Form Inverse(U) * Inverse(L).
c
      do k = n - 1, 1, -1

        do i = k + 1, n
          work(i) = a(i,k)
          a(i,k) = 0.0D+00
        end do

        do j = k + 1, n
          do i = 1, n
            a(i,k) = a(i,k) + a(i,j) * work(j)
          end do
        end do

        if ( pivot(k) .ne. k ) then

          do i = 1, n
            t             = a(i,k)
            a(i,k)        = a(i,pivot(k))
            a(i,pivot(k)) = t
          end do

        end if

      end do

      return
      end
      subroutine r8ge_sl ( n, a_lu, pivot, b, job )

c*********************************************************************72
c
cc R8GE_SL solves a system factored by R8GE_FA.
c
c  Discussion:
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c    R8GE storage is used by LINPACK and LAPACK.
c
c    R8GE_SL is a simplified version of the LINPACK routine SGESL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 March 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input, double precision A_LU(N,N), the LU factors from R8GE_FA.
c
c    Input, integer PIVOT(N), the pivot vector from R8GE_FA.
c
c    Input/output, double precision B(N).
c    On input, the right hand side vector.
c    On output, the solution vector.
c
c    Input, integer JOB, specifies the operation.
c    0, solve A * x = b.
c    nonzero, solve A' * x = b.
c
      implicit none

      integer n

      double precision a_lu(n,n)
      double precision b(n)
      integer i
      integer job
      integer k
      integer l
      integer pivot(n)
      double precision t
c
c  Solve A * x = b.
c
      if ( job .eq. 0 ) then
c
c  Solve PL * Y = B.
c
        do k = 1, n - 1

          l = pivot(k)

          if ( l .ne. k ) then
            t    = b(l)
            b(l) = b(k)
            b(k) = t
          end if

          do i = k + 1, n
            b(i) = b(i) + a_lu(i,k) * b(k)
          end do

        end do
c
c  Solve U * X = Y.
c
        do k = n, 1, -1
          b(k) = b(k) / a_lu(k,k)
          do i = 1, k - 1
            b(i) = b(i) - a_lu(i,k) * b(k)
          end do
        end do
c
c  Solve A' * X = B.
c
      else
c
c  Solve U' * Y = B.
c
        do k = 1, n
          do i = 1, k - 1
            b(k) = b(k) - a_lu(i,k) * b(i)
          end do
          b(k) = b(k) / a_lu(k,k)
        end do
c
c  Solve ( PL )' * X = Y.
c
        do k = n - 1, 1, -1

          do i = k + 1, n
            b(k) = b(k) + a_lu(i,k) * b(i)
          end do

          l = pivot(k)

          if ( l .ne. k ) then
            t    = b(l)
            b(l) = b(k)
            b(k) = t
          end if

        end do

      end if

      return
      end

