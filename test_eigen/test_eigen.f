      subroutine get_unit ( unit )

c*********************************************************************72
c
cc GET_UNIT returns a free FORTRAN unit number.
c
c  Discussion:
c
c    A "free" FORTRAN unit number is an integer between 1 and 99 which
c    is not currently associated with an I/O device.  A free FORTRAN unit
c    number is needed in order to open a file with the OPEN command.
c
c    If UNIT = 0, then no free FORTRAN unit could be found, although
c    all 99 units were checked (except for units 5, 6 and 9, which
c    are commonly reserved for console I/O).
c
c    Otherwise, UNIT is an integer between 1 and 99, representing a
c    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
c    are special, and will never return those values.
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
c    Output, integer UNIT, the free unit number.
c
      implicit none

      integer i
      integer unit

      unit = 0

      do i = 1, 99

        if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

          open ( unit = i, err = 10 )
          close ( unit = i )

          unit = i

          return
        end if

10      continue

      end do

      return
      end
      subroutine r4vec_uniform_01 ( n, seed, r )

c*********************************************************************72
c
cc R4VEC_UNIFORM_01 returns a unit pseudorandom R4VEC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 March 2006
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
c    Output, real R(N), the vector of pseudorandom values.
c
      implicit none

      integer n

      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer k
      integer seed
      real r(n)

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R4VEC_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + i4_huge
        end if

        r(i) = real ( seed ) * 4.656612875E-10

      end do

      return
      end
      function r8_normal_01 ( seed )

c*********************************************************************72
c
cc R8_NORMAL_01 returns a unit pseudonormal R8.
c
c  Discussion:
c
c    Because this routine uses the Box Muller method, it requires pairs
c    of uniform random values to generate a pair of normal random values.
c    This means that on every other call, the code can use the second
c    value that it calculated.
c
c    However, if the user has changed the SEED value between calls,
c    the routine automatically resets itself and discards the saved data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, double precision R8_NORMAL_01, a sample of the standard normal PDF.
c
      implicit none

      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r1
      double precision r2
      double precision r8_normal_01
      double precision r8_uniform_01
      integer seed
      integer seed1
      integer seed2
      integer seed3
      integer used
      double precision v1
      double precision v2

      save seed1
      save seed2
      save seed3
      save used
      save v2

      data seed2 / 0 /
      data used / 0 /
      data v2 / 0.0D+00 /
c
c  If USED is odd, but the input SEED does not match
c  the output SEED on the previous call, then the user has changed
c  the seed.  Wipe out internal memory.
c
      if ( mod ( used, 2 ) == 1 ) then

        if ( seed .ne. seed2 ) then
          used = 0
          seed1 = 0
          seed2 = 0
          seed3 = 0
          v2 = 0.0D+00
        end if

      end if
c
c  If USED is even, generate two uniforms, create two normals,
c  return the first normal and its corresponding seed.
c
      if ( mod ( used, 2 ) .eq. 0 ) then

        seed1 = seed

        r1 = r8_uniform_01 ( seed )

        if ( r1 .eq. 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8_NORMAL_01 - Fatal error!'
          write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
          stop
        end if

        seed2 = seed

        r2 = r8_uniform_01 ( seed )

        seed3 = seed

        v1 = sqrt ( -2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * pi * r2 )
        v2 = sqrt ( -2.0D+00 * log ( r1 ) ) * sin ( 2.0D+00 * pi * r2 )

        r8_normal_01 = v1
        seed = seed2
c
c  If USED is odd (and the input SEED matched the output value from
c  the previous call), return the second normal and its corresponding seed.
c
      else

        r8_normal_01 = v2
        seed = seed3

      end if

      used = used + 1

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
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R8_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer k
      double precision r8_uniform_01
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
        seed = seed + i4_huge
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r8_uniform_01 = dble ( seed ) * 4.656612875D-10

      return
      end
      subroutine r8bin_print ( bin_num, bin, bin_limit, title )

c*********************************************************************72
c
cc R8BIN_PRINT prints the bins of a real vector.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer BIN_NUM, the number of bins.
c
c    Input, integer BIN(0:BIN_NUM+1).
c    BIN(0) counts entries of X less than BIN_LIMIT(0).
c    BIN(BIN_NUM+1) counts entries greater than or equal to BIN_LIMIT(BIN_NUM).
c    For 1 <= I <= BIN_NUM, BIN(I) counts the entries X(J) such that
c      BIN_LIMIT(I-1) <= X(J) < BIN_LIMIT(I).
c    where H is the bin spacing.
c
c    Input, double precision BIN_LIMIT(0:BIN_NUM), the "limits" of the bins.
c    BIN(I) counts the number of entries X(J) such that
c      BIN_LIMIT(I-1) <= X(J) < BIN_LIMIT(I).
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer bin_num

      integer bin(0:bin_num+1)
      double precision bin_limit(0:bin_num)
      integer i
      character ( len = * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '   Index     Lower Limit   Count     Upper Limit'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,i6,2x,14x,2x,i6,2x,g14.6)' ) 
     &  0, bin(0), bin_limit(0)
      do i = 1, bin_num
        write ( *, '(2x,i6,2x,g14.6,2x,i6,2x,g14.6)' ) 
     &    i, bin_limit(i-1), bin(i), bin_limit(i)
      end do
      write ( *, '(2x,i6,2x,g14.6,2x,i6)') 
     &  bin_num + 1, bin_limit(bin_num), bin(bin_num+1)

      return
      end
      subroutine r8mat_house_axh ( n, a, v, ah )

c*********************************************************************72
c
cc R8MAT_HOUSE_AXH computes A*H where H is a compact Householder matrix.
c
c  Discussion:
c
c    The Householder matrix H(V) is defined by
c
c      H(V) = I - 2 * v * v' / ( v' * v )
c
c    This routine is not particularly efficient.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 February 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of A.
c
c    Input, double precision A(N,N), the matrix.
c
c    Input, double precision V(N), a vector defining a Householder matrix.
c
c    Output, double precision AH(N,N), the product A*H.
c
      implicit none

      integer n

      double precision a(n,n)
      double precision ah(n,n)
      double precision av(n)
      integer i
      integer j
      double precision v(n)
      double precision v_normsq

      v_normsq = 0.0D+00
      do i = 1, n
        v_normsq = v_normsq + v(i)**2
      end do
      
      do i = 1, n
        av(i) = 0.0D+00
        do j = 1, n
          av(i) = av(i) + a(i,j) * v(j)
        end do
      end do

      do i = 1, n
        do j = 1, n
          ah(i,j) = a(i,j)
        end do
      end do

      do i = 1, n
        do j = 1, n
          ah(i,j) = ah(i,j) - 2.0D+00 * av(i) * v(j)
        end do
      end do

      do i = 1, n
        do j = 1, n
          ah(i,j) = ah(i,j) / v_normsq
        end do
      end do

      return
      end
      subroutine r8mat_orth_uniform ( n, seed, a )

c*********************************************************************72
c
cc R8MAT_ORTH_UNIFORM returns a random orthogonal R8MAT.
c
c  Discussion:
c
c    An R8MAT is a two dimensional matrix of R8 values.
c
c    Thanks to Eugene Petrov, B I Stepanov Institute of Physics,
c    National Academy of Sciences of Belarus, for convincingly
c    pointing out the severe deficiencies of an earlier version of
c    this routine.
c
c    Essentially, the computation involves saving the Q factor of the
c    QR factorization of a matrix whose entries are normally distributed.
c    However, it is only necessary to generate this matrix a column at
c    a time, since it can be shown that when it comes time to annihilate
c    the subdiagonal elements of column K, these (transformed) elements of
c    column K are still normally distributed random values.  Hence, there
c    is no need to generate them at the beginning of the process and
c    transform them K-1 times.
c
c    For computational efficiency, the individual Householder transformations
c    could be saved, as recommended in the reference, instead of being
c    accumulated into an explicit matrix format.
c
c  Properties:
c
c    The inverse of A is equal to A'.
c
c    A * A'  = A' * A = I.
c
c    Columns and rows of A have unit Euclidean norm.
c
c    Distinct pairs of columns of A are orthogonal.
c
c    Distinct pairs of rows of A are orthogonal.
c
c    The L2 vector norm of A*x = the L2 vector norm of x for any vector x.
c
c    The L2 matrix norm of A*B = the L2 matrix norm of B for any matrix B.
c
c    The determinant of A is +1 or -1.
c
c    All the eigenvalues of A have modulus 1.
c
c    All singular values of A are 1.
c
c    All entries of A are between -1 and 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 November 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Pete Stewart,
c    Efficient Generation of Random Orthogonal Matrices With an Application
c    to Condition Estimators,
c    SIAM Journal on Numerical Analysis,
c    Volume 17, Number 3, June 1980, pages 403-409.
c
c  Parameters:
c
c    Input, integer N, the order of A.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, double precision A(N,N), the orthogonal matrix.
c
      implicit none

      integer n

      double precision a(n,n)
      integer i
      integer j
      double precision r8_normal_01
      integer seed
      double precision v(n)
      double precision x(n)
c
c  Start with A = the identity matrix.
c
      do i = 1, n
        do j = 1, n
          if ( i .eq. j ) then
            a(i,j) = 1.0D+00
          else
            a(i,j) = 0.0D+00
          end if
        end do
      end do
c
c  Now behave as though we were computing the QR factorization of
c  some other random matrix.  Generate the N elements of the first column,
c  compute the Householder matrix H1 that annihilates the subdiagonal elements,
c  and set A := A * H1' = A * H.
c
c  On the second step, generate the lower N-1 elements of the second column,
c  compute the Householder matrix H2 that annihilates them,
c  and set A := A * H2' = A * H2 = H1 * H2.
c
c  On the N-1 step, generate the lower 2 elements of column N-1,
c  compute the Householder matrix HN-1 that annihilates them, and
c  and set A := A * H(N-1)' = A * H(N-1) = H1 * H2 * ... * H(N-1).
c  This is our random orthogonal matrix.
c
      do j = 1, n - 1
c
c  Set the vector that represents the J-th column to be annihilated.
c
        do i = 1, j - 1
          x(i) = 0.0D+00
        end do

        do i = j, n
          x(i) = r8_normal_01 ( seed )
        end do
c
c  Compute the vector V that defines a Householder transformation matrix
c  H(V) that annihilates the subdiagonal elements of X.
c
        call r8vec_house_column ( n, x, j, v )
c
c  Postmultiply the matrix A by H'(V) = H(V).
c
        call r8mat_house_axh ( n, a, v, a )

      end do

      return
      end
      subroutine r8mat_print ( m, n, a, title )

c*********************************************************************72
c
cc R8MAT_PRINT prints an R8MAT.
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

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

        end do

      end do

      return
      end
      subroutine r8symm_test ( n, lambda_mean, lambda_dev, seed, a, 
     &  q, lambda )

c*********************************************************************72
c
cc R8SYMM_TEST determines a symmetric matrix with a certain eigenstructure.
c
c  Discussion:
c
c    An R8SYMM is a real symmetric matrix stored using full storage, and
c    using R8 arithmetic.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 February 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision LAMBDA_MEAN, the mean value of the normal
c    distribution from which the eigenvalues will be chosen.
c
c    Input, double precision LAMBDA_DEV, the standard deviation of the normal
c    distribution from which the eigenvalues will be chosen.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, double precision A(N,N), the test matrix.
c
c    Output, double precision Q(N,N), the eigenvector matrix.
c
c    Output, double precision LAMBDA(N), the eigenvalue vector.
c
      implicit none

      integer n

      double precision a(n,n)
      integer i
      integer j
      integer k
      double precision lambda(n)
      double precision lambda_dev
      double precision lambda_mean
      double precision q(n,n)
      integer seed
c
c  Choose the eigenvalues LAMBDA.
c
      call r8vec_normal ( n, lambda_mean, lambda_dev, seed, lambda )
c
c  Get a random orthogonal matrix Q.
c
      call r8mat_orth_uniform ( n, seed, q )
c
c  Set A = Q * Lambda*I * Q'.
c
      do i = 1, n
        do j = 1, n
          a(i,j) = 0.0D+00
        end do
      end do

      do i = 1, n
        do j = 1, n
          do k = 1, n
            a(i,j) = a(i,j) + q(i,k) * lambda(k) * q(j,k)
          end do
        end do
      end do

      return
      end
      subroutine r8vec_bin ( n, x, bin_num, bin_min, bin_max, bin, 
     &  bin_limit )

c*********************************************************************72
c
cc R8VEC_BIN computes bins based on a given R8VEC.
c
c  Discussion:
c
c    The user specifies minimum and maximum bin values, BIN_MIN and
c    BIN_MAX, and the number of bins, BIN_NUM.  This determines a
c    "bin width":
c
c      H = ( BIN_MAX - BIN_MIN ) / BIN_NUM
c
c    so that bin I will count all entries X(J) such that
c
c      BIN_LIMIT(I-1) <= X(J) < BIN_LIMIT(I).
c
c    The array X does NOT have to be sorted.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 July 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries of X.
c
c    Input, double precision X(N), an (unsorted) array to be binned.
c
c    Input, integer BIN_NUM, the number of bins.  Two extra bins, #0 and
c    #BIN_NUM+1, count extreme values.
c
c    Input, real ( kind = 8 ) BIN_MIN, BIN_MAX, define the range and size
c    of the bins.  BIN_MIN and BIN_MAX must be distinct.
c    Normally, BIN_MIN < BIN_MAX, and the documentation will assume
c    this, but proper results will be computed if BIN_MIN > BIN_MAX.
c
c    Output, integer BIN(0:BIN_NUM+1).
c    BIN(0) counts entries of X less than BIN_MIN.
c    BIN(BIN_NUM+1) counts entries greater than or equal to BIN_MAX.
c    For 1 <= I <= BIN_NUM, BIN(I) counts the entries X(J) such that
c      BIN_LIMIT(I-1) <= X(J) < BIN_LIMIT(I).
c    where H is the bin spacing.
c
c    Output, double precision BIN_LIMIT(0:BIN_NUM), the "limits" of the bins.
c    BIN(I) counts the number of entries X(J) such that
c      BIN_LIMIT(I-1) <= X(J) < BIN_LIMIT(I).
c
      implicit none

      integer n
      integer bin_num

      integer bin(0:bin_num+1)
      double precision bin_limit(0:bin_num)
      double precision bin_max
      double precision bin_min
      integer i
      integer j
      double precision t
      double precision x(n)

      if ( bin_max .eq. bin_min ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8VEC_BIN - Fatal error!'
        write ( *, '(a)' ) '  BIN_MIN = BIN_MAX.'
        stop
      end if

      do i = 0, bin_num + 1
        bin(i) = 0
      end do

      do i = 1, n

        t = ( x(i) - bin_min ) / ( bin_max - bin_min )

        if ( t .lt. 0.0D+00 ) then
          j = 0
        else if ( 1.0D+00 .le. t ) then
          j = bin_num + 1
        else
          j = 1 + int ( dble ( bin_num ) * t )
        end if

        bin(j) = bin(j) + 1

      end do
c
c  Compute the bin limits.
c
      do i = 0, bin_num
        bin_limit(i) = (   dble ( bin_num - i ) * bin_min   
     &                   + dble (           i ) * bin_max ) 
     &                   / dble ( bin_num     )
      end do

      return
      end
      subroutine r8vec_house_column ( n, a, k, v )

c*********************************************************************72
c
cc R8VEC_HOUSE_COLUMN defines a Householder premultiplier that "packs" a column.
c
c  Discussion:
c
c    The routine returns a vector V that defines a Householder
c    premultiplier matrix H(V) that zeros out the subdiagonal entries of
c    column K of the matrix A.
c
c       H(V) = I - 2 * v * v'
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 February 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix A.
c
c    Input, double precision A(N), column K of the matrix A.
c
c    Input, integer K, the column of the matrix to be modified.
c
c    Output, double precision V(N), a vector of unit L2 norm which defines an
c    orthogonal Householder premultiplier matrix H with the property
c    that the K-th column of H*A is zero below the diagonal.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      integer k
      double precision s
      double precision v(n)
      double precision vnorm

      do i = 1, n
        v(i) = 0.0D+00
      end do

      if ( k .lt. 1 .or. n .le. k ) then
        return
      end if

      s = 0.0D+00
      do i = k, n
        s = s + a(i)**2
      end do
      s = sqrt ( s )

      if ( s .eq. 0.0D+00 ) then
        return
      end if

      v(k) = a(k) + sign ( s, a(k) )
      do i = k + 1, n
        v(i) = a(i)
      end do

      vnorm = 0.0D+00
      do i = k, n
        vnorm = vnorm + v(i) * v(i)
      end do
      vnorm = sqrt ( vnorm )

      do i = k, n
        v(i) = v(i) / vnorm
      end do

      return
      end
      subroutine r8vec_normal ( n, a, b, seed, x )

c*********************************************************************72
c
cc R8VEC_NORMAL returns a scaled pseudonormal R8VEC.
c
c  Discussion:
c
c    The standard normal probability distribution function (PDF) has
c    mean 0 and standard deviation 1.
c
c    This routine can generate a vector of values on one call.  It
c    has the feature that it should provide the same results
c    in the same order no matter how we break up the task.
c
c    Before calling this routine, the user may call RANDOM_SEED
c    in order to set the seed of the random number generator.
c
c    The Box-Muller method is used, which is efficient, but
c    generates an even number of values each time.  On any call
c    to this routine, an even number of new values are generated.
c    Depending on the situation, one value may be left over.
c    In that case, it is saved for the next call.
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
c  Parameters:
c
c    Input, integer N, the number of values desired.  If N is negative,
c    then the code will flush its internal memory; in particular,
c    if there is a saved value to be used on the next call, it is
c    instead discarded.  This is useful if the user has reset the
c    random number seed, for instance.
c
c    Input, real ( kind = 8 ) A, B, the mean and standard deviation.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, double precision X(N), a sample of the standard normal PDF.
c
c  Local parameters:
c
c    Local, integer MADE, records the number of values that have
c    been computed.  On input with negative N, this value overwrites
c    the return value of N, so the user can get an accounting of
c    how much work has been done.
c
c    Local, integer SAVED, is 0 or 1 depending on whether there is a
c    single saved value left over from the previous call.
c
c    Local, integer X_LO_INDEX, X_HI_INDEX, records the range of entries of
c    X that we need to compute.  This starts off as 1:N, but is adjusted
c    if we have a saved value that can be immediately stored in X(1),
c    and so on.
c
c    Local, double precision Y, the value saved from the previous call, if
c    SAVED is 1.
c
      implicit none

      integer n

      double precision a
      double precision b
      integer i
      integer m
      integer made
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r(2)
      double precision r8_uniform_01
      integer saved
      integer seed
      double precision x(n)
      integer x_hi_index
      integer x_lo_index
      double precision y

      save made
      save saved
      save y

      data made / 0 /
      data saved / 0 /
      data y / 0.0D+00 /
c
c  I'd like to allow the user to reset the internal data.
c  But this won't work properly if we have a saved value Y.
c  I'm making a crock option that allows the user to signal
c  explicitly that any internal memory should be flushed,
c  by passing in a negative value for N.
c
      if ( n .lt. 0 ) then
        n = made
        made = 0
        saved = 0
        y = 0.0D+00
        return
      else if ( n .eq. 0 ) then
        return
      end if
c
c  Record the range of X we need to fill in.
c
      x_lo_index = 1
      x_hi_index = n
c
c  Use up the old value, if we have it.
c
      if ( saved .eq. 1 ) then
        x(1) = y
        saved = 0
        x_lo_index = 2
      end if
c
c  Maybe we don't need any more values.
c
      if ( x_hi_index - x_lo_index + 1 .eq. 0 ) then
c
c  If we need just one new value, do that here to avoid null arrays.
c
      else if ( x_hi_index - x_lo_index + 1 .eq. 1 ) then

        r(1) = r8_uniform_01 ( seed )

        if ( r(1) .eq. 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8VEC_NORMAL - Fatal error!'
          write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
          stop
        end if

        r(2) = r8_uniform_01 ( seed )

        x(x_hi_index) =
     &           sqrt ( -2.0D+00 * log ( r(1) ) )
     &           * cos ( 2.0D+00 * pi * r(2) )
        y =      sqrt ( -2.0D+00 * log ( r(1) ) )
     &           * sin ( 2.0D+00 * pi * r(2) )

        saved = 1

        made = made + 2
c
c  If we require an even number of values, that's easy.
c
      else if ( mod ( x_hi_index - x_lo_index + 1, 2 ) .eq. 0 ) then

        do i = x_lo_index, x_hi_index, 2

          call r8vec_uniform_01 ( 2, seed, r )

          x(i) =
     &      sqrt ( -2.0D+00 * log ( r(1) ) )
     &      * cos ( 2.0D+00 * pi * r(2) )

          x(i+1) =
     &      sqrt ( -2.0D+00 * log ( r(1) ) )
     &      * sin ( 2.0D+00 * pi * r(2) )

        end do

        made = made + x_hi_index - x_lo_index + 1
c
c  If we require an odd number of values, we generate an even number,
c  and handle the last pair specially, storing one in X(N), and
c  saving the other for later.
c
      else

        do i = x_lo_index, x_hi_index - 1, 2

          call r8vec_uniform_01 ( 2, seed, r )

          x(i) =
     &      sqrt ( -2.0D+00 * log ( r(1) ) )
     &      * cos ( 2.0D+00 * pi * r(2) )

          x(i+1) =
     &      sqrt ( -2.0D+00 * log ( r(1) ) )
     &      * sin ( 2.0D+00 * pi * r(2) )

        end do

        call r8vec_uniform_01 ( 2, seed, r )

        x(n) = sqrt ( -2.0D+00 * log ( r(1) ) )
     &    * cos ( 2.0D+00 * pi * r(1) )

        y = sqrt ( -2.0D+00 * log ( r(2) ) )
     &    * sin ( 2.0D+00 * pi * r(2) )

        saved = 1

        made = made + x_hi_index - x_lo_index + 2

      end if

      do i = 1, n
        x(i) = a + b * x(i)
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
c    An R8VEC is an array of double precision real values.
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
      integer s_len_trim
      character ( len = * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,2x,g16.8)' ) i, a(i)
      end do

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
      subroutine r8vec2_print ( n, a1, a2, title )

c*********************************************************************72
c
cc R8VEC2_PRINT prints an R8VEC2.
c
c  Discussion:
c
c    An R8VEC2 is a dataset consisting of N pairs of R8s, stored
c    as two separate vectors A1 and A2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components of the vector.
c
c    Input, double precision A1(N), A2(N), the vectors to be printed.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer n

      double precision a1(n)
      double precision a2(n)
      integer i
      character ( len = * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,a,1x,g14.6,2x,g14.6)' ) i, ':', a1(i), a2(i)
      end do

      return
      end
      function s_len_trim ( s )

c*********************************************************************72
c
cc S_LEN_TRIM returns the length of a string to the last nonblank.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 March 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) S, a string.
c
c    Output, integer S_LEN_TRIM, the length of the string to the last nonblank.
c
      implicit none

      integer i
      character*(*) s
      integer s_len_trim

      do i = len ( s ), 1, -1

        if ( s(i:i) .ne. ' ' ) then
          s_len_trim = i
          return
        end if

      end do

      s_len_trim = 0

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
