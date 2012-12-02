      program main

c*********************************************************************72
c
cc MAIN is the main program for SVD_DEMO.
c
c  Discussion:
c
c    SVD_DEMO demonstrates the SVD.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modififed:
c
c    19 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Usage:
c
c    svd_demo m n [seed]
c
c  Command Parameters:
c
c    Command parameter, integer M, N, the number of rows and 
c    columns of the matrix.  If M or N is not supplied on the command line,
c    the user is prompted to supply them.
c
c    Command parameter, integer SEED, a seed for the random number
c    generator.  If SEED is not supplied on the command line, a value
c    is generated internally.
c
c  Local Parameters:
c
c    Local, double precision A(M,N), the matrix whose singular value
c    decomposition we are investigating.
c
c    Local, double precision S(M,N), the diagonal factor
c    in the singular value decomposition of A.
c
c    Output, double precision U(M,M), the first orthogonal factor
c    in the singular value decomposition of A.
c
c    Output, double precision V(N,N), the second orthogonal factor
c    in the singular value decomposition of A.
c
      implicit none

      integer arg_num
      integer i
      integer iarg
      integer iargc
      integer ierror
      integer length
      integer lwork
      integer m
      integer n
      integer seed
      character ( len = 80 ) string

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SVD_DEMO'
      write ( *, '(a)' ) '  FORTRAN90 version'
      write ( *, '(a)' ) 
     &  '  Demonstrate the singular value decomposition (SVD)'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  A real MxN matrix A can be factored as:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    A = U * S * V'''
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  where'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    U = MxM orthogonal,'
      write ( *, '(a)' ) '    S = MxN zero except for diagonal,'
      write ( *, '(a)' ) '    V = NxN orthogonal.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  The diagonal of S contains only nonnegative numbers'
      write ( *, '(a)' ) '  and these are arranged in descending order.'
c
c  Get the number of command line arguments.
c
      arg_num = iargc ( )
c
c  If at least one command line argument, it's M.
c
      if ( 1 .le. arg_num ) then

        iarg = 1
        call getarg ( iarg, string )
        call s_to_i4 ( string, m, ierror, length )

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SVD_DEMO:'
        write ( *, '(a)' ) '  Please enter the value of M.'
        write ( *, '(a)' ) '  (Number of rows in matrix A)'
        write ( *, '(a)' ) '  (We prefer M <= 10!).'

        read ( *, * ) m

      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix row order    M = ', m
c
c  If at least two command line argument, it's N.
c
      if ( 2 .le. arg_num ) then

        iarg = 2
        call getarg ( iarg, string )
        call s_to_i4 ( string, n, ierror, length )

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SVD_DEMO:'
        write ( *, '(a)' ) '  Please enter the value of N.'
        write ( *, '(a)' ) '  (Number of columns in matrix A)'
        write ( *, '(a)' ) '  (We prefer N <= 10!).'

        read ( *, * ) n

      end if

      write ( *, '(a,i8)' ) '  Matrix column order N = ', n
c
c  If a third command line argument, it's SEED.
c
      if ( 3 .le. arg_num ) then

        iarg = 3
        call getarg ( iarg, string )
        call s_to_i4 ( string, seed, ierror, length )
        write ( *, '(a,i12)' ) '  Random number SEED    = ', seed
        write ( *, '(a)' ) '  (Chosen by user.)'

      else

        call get_seed ( seed )
        write ( *, '(a,i12)' ) '  Random number SEED    = ', seed
        write ( *, '(a)' ) '  (Chosen by program.)'

      end if
c
c  FORTRAN77 is clumsy about allocating space for arrays.
c  We should be able to do this by calling a subroutine.
c
      call svd_demo_sub ( m, n, seed )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SVD_DEMO:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine svd_demo_sub ( m, n, seed )

c*********************************************************************72
c
cc SVD_SUB carries out the demonstration, after dimensions are known.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modififed:
c
c    18 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and 
c    columns of the matrix.  If M or N is not supplied on the command line,
c    the user is prompted to supply them.
c
c    Input, integer SEED, a seed for the random number
c    generator.  If SEED is not supplied on the command line, a value
c    is generated internally.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision a_pseudo(n,m)
      integer i
      integer ierror
      integer j
      integer length
      integer lwork
      double precision s(m,n)
      double precision s2(m,n)
      integer seed
      double precision u(m,m)
      double precision u2(m,m)
      double precision v(n,n)
      double precision v2(n,n)
c
c  Generate the matrix A.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  We choose a "random" matrix A, with'
      write ( *, '(a)' ) '  integral values between 0 and 10.'

      call r8mat_uniform_01 ( m, n, seed, a )

      do j = 1, n
        do i = 1, m
          a(i,j) = anint ( 10.0D+00 * a(i,j) )
        end do
      end do

      call r8mat_print ( m, n, a, '  The matrix A:' )
c
c  Get the SVD from LAPACK.
c
      lwork = max ( 3 * min ( m, n ) + max ( m, n ), 5 * min ( m, n ) )
      call get_svd_lapack ( m, n, a, u, s, v, lwork )
c
c  Print the SVD.
c
      call r8mat_print ( m, m, u, '  The orthogonal factor U:' )

      call r8mat_print ( m, n, s, '  The diagonal factor S:' )

      call r8mat_print ( n, n, v, '  The orthogonal factor V:' )
c
c  Check that A = U * S * V'.
c
      call svd_product_test ( m, n, a, u, s, v )
c
c  Compute the norm of the difference between A and the successive
c  sums of rank one approximants.
c
      call rank_one_test ( m, n, a, u, s, v )
c
c  Actually print the sums of rank one approximants.
c
      call rank_one_print_test ( m, n, a, u, s, v )
c
c  Compute the pseudoinverse.
c
      call pseudo_inverse ( m, n, u, s, v, a_pseudo )

      call r8mat_print ( n, m, a_pseudo, '  The pseudoinverse of A:' )
c
c  Test A*A+ = I+, A+*A = I+
c
      call pseudo_product_test ( m, n, a, a_pseudo )
c
c  Demonstrate the use of the pseudoinverse for linear systems.
c
      call pseudo_linear_solve_test ( m, n, a, a_pseudo, seed )
c
c  Get the SVD from LINPACK.
c
      call get_svd_linpack ( m, n, a, u2, s2, v2 )

      call compare_linpack_lapack ( m, n, u, s, v, u2, s2, v2 )

      return
      end
      subroutine compare_linpack_lapack ( m, n, u, s, v, u2, s2, v2 )

c*********************************************************************72
c
cc COMPARE_LINPACK_LAPACK compares the SVD's from LINPACK and LAPACK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns
c    in the matrix A.
c
c    Input, double precision U(M,M), S(M,N), V(N,N), the factors
c    that form the singular value decomposition of A, computed by LAPACK.
c
c    Input, double precision U2(M,M), S2(M,N), V2(N,N), the factors
c    that form the singular value decomposition of A, computed by LINPACK.
c
      implicit none

      integer m
      integer n

      double precision r8mat_dif_fro
      double precision s(m,n)
      double precision s2(m,n)
      double precision u(m,m)
      double precision u2(m,m)
      double precision v(n,n)
      double precision v2(n,n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'COMPARE_LINPACK_LAPACK:'
      write ( *, '(a)' ) 
     &  '  While the singular values should be identical,'
      write ( *, '(a)' ) 
     &  '  the orthogonal factors may have some differences.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Frobenius differences:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  U(LAPACK) - U(LINPACK): ', 
     &  r8mat_dif_fro ( m, m, u, u2 )
      write ( *, '(a,g14.6)' ) '  S(LAPACK) - S(LINPACK): ', 
     &  r8mat_dif_fro ( m, n, s, s2 )
      write ( *, '(a,g14.6)' ) '  V(LAPACK) - V(LINPACK): ', 
     &  r8mat_dif_fro ( n, n, v, v2 )

      return
      end
      subroutine get_seed ( seed )

c*********************************************************************72
c
cc GET_SEED returns a seed for the random number generator.
c
c  Discussion:
c
c    The seed depends on the current time, and ought to be (slightly)
c    different every millisecond.  Thus, calling this routine several
c    times in succession will probably return the SAME seed, but
c    calling it a few minutes or days apart will turn a suitably
c    "random" seed.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer SEED, a pseudorandom seed value.
c
      implicit none

      integer day
      integer hour
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer milli
      integer minute
      integer month
      integer second
      integer seed
      double precision temp
      character * ( 10 ) time
      character * ( 8 ) date
      integer year

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) year, month, day
      read ( time, '(i2,i2,i2,1x,i3)' ) hour, minute, second, milli

      temp = 0.0D+00
      temp = temp + dble ( month - 1 ) / 11.0D+00
      temp = temp + dble ( day   - 1 ) / 30.0D+00
      temp = temp + dble ( hour      ) / 23.0D+00
      temp = temp + dble ( minute    ) / 59.0D+00
      temp = temp + dble ( second    ) / 59.0D+00
      temp = temp + dble ( milli     ) / 999.0D+00

      temp = temp / 6.0D+00
c
c  Force 0 < TEMP <= 1.
c
10    continue

      if ( temp .le. 0.0D+00 ) then
        temp = temp + 1.0D+00
        go to 10
      end if

20    continue

      if ( 1.0D+00 .lt. temp ) then
        temp = temp - 1.0D+00
        go to 20
      end if

      seed = int ( dble ( i4_huge ) * temp )
c
c  Never use a seed of 0 or maximum integer.
c
      if ( seed .eq. 0 ) then
        seed = 1
      end if

      if ( seed .eq. i4_huge ) then
        seed = seed - 1
      end if

      return
      end
      subroutine get_svd_lapack ( m, n, a, u, s, vt, lwork )

c*********************************************************************72
c
cc GET_SVD_LAPACK gets the SVD of a matrix using a call to LAPACK.
c
c  Discussion:
c
c    The singular value decomposition of a real MxN matrix A has the form:
c
c      A = U * S * V'
c
c    where
c
c      U is MxM orthogonal,
c      S is MxN, and entirely zero except for the diagonal;
c      V is NxN orthogonal.
c
c    Moreover, the nonzero entries of S are positive, and appear
c    in order, from largest magnitude to smallest.
c
c    This routine calls the LAPACK routine DGESVD to compute the
c    factorization.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns
c    in the matrix A.
c
c    Input, double precision A(M,N), the matrix whose singular value
c    decomposition we are investigating.
c
c    Output, double precision U(M,M), S(M,N), V(N,N), the factors
c    that form the singular value decomposition of A.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision a_copy(m,n)
      integer i
      integer info
      integer j
      integer lda
      integer ldu
      integer ldv
      character jobu
      character jobv
      integer lwork
      double precision sdiag(min(m,n))
      double precision s(m,n)
      double precision u(m,m)
      double precision v(n,n)
      double precision vt(n,n)
      double precision work(lwork)
c
c  Compute the eigenvalues and eigenvectors.
c
      jobu = 'A'
      jobv = 'A'
      lda = m
      ldu = m
      ldv = n
c
c  The input matrix is destroyed by the routine.  Since we need to keep
c  it around, we only pass a copy to the routine.
c
      do j = 1, n
        do i = 1, m
          a_copy(i,j) = a(i,j)
        end do
      end do

      call dgesvd ( jobu, jobv, m, n, a_copy, lda, sdiag, u, ldu, v, 
     &  ldv, work, lwork, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GET_SVD_LAPACK - Failure!'
        write ( *, '(a)' ) '  The SVD could not be calculated.'
        write ( *, '(a)' ) '  LAPACK routine DGESVD returned a nonzero'
        write ( *, '(a,i8)' ) '  value of the error flag, INFO = ', info
        return
      end if
c
c  Make the MxN matrix S from the diagonal values in SDIAG.
c
      do j = 1, n
        do i = 1, m
          s(i,j) = 0.0D+00
        end do
      end do

      do i = 1, min ( m, n )
        s(i,i) = sdiag(i)
      end do
c
c  Transpose V.
c
      call r8mat_transpose ( n, n, v, vt )

      return
      end
      subroutine get_svd_linpack ( m, n, a, u, s, v )

c*********************************************************************72
c
cc GET_SVD_LINPACK gets the SVD of a matrix using a call to LINPACK.
c
c  Discussion:
c
c    The singular value decomposition of a real MxN matrix A has the form:
c
c      A = U * S * V'
c
c    where
c
c      U is MxM orthogonal,
c      S is MxN, and entirely zero except for the diagonal;
c      V is NxN orthogonal.
c
c    Moreover, the nonzero entries of S are positive, and appear
c    in order, from largest magnitude to smallest.
c
c    This routine calls the LINPACK routine DSVDC to compute the
c    factorization.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns
c    in the matrix A.
c
c    Input, double precision A(M,N), the matrix whose singular value
c    decomposition we are investigating.
c
c    Output, double precision U(M,M), S(M,N), V(N,N), the factors
c    that form the singular value decomposition of A.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision a_copy(m,n)
      double precision e(max(m+1,n))
      integer i
      integer info
      integer j
      integer lda
      integer ldu
      integer ldv
      integer job
      integer lwork
      double precision s(m,n)
      double precision sdiag(max(m+1,n))
      double precision u(m,m)
      double precision v(n,n)
      double precision work(m)
c
c  Compute the eigenvalues and eigenvectors.
c
      job = 11
      lda = m
      ldu = m
      ldv = n
c
c  The input matrix is destroyed by the routine.  Since we need to keep
c  it around, we only pass a copy to the routine.
c
      do j = 1, n
        do i = 1, m
          a_copy(i,j) = a(i,j)
        end do
      end do

      call dsvdc ( a_copy, lda, m, n, sdiag, e, u, ldu, v, ldv, work, 
     &  job, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GET_SVD_LINPACK - Failure!'
        write ( *, '(a)' ) '  The SVD could not be calculated.'
        write ( *, '(a)' ) '  LINPACK routine DSVDC returned a nonzero'
        write ( *, '(a,i8)' ) '  value of the error flag, INFO = ', info
        return
      end if
c
c  Make the MxN matrix S from the diagonal values in SDIAG.
c
      do j = 1, n
        do i = 1, m
          s(i,j) = 0.0D+00
        end do
      end do

      do i = 1, min ( m, n )
        s(i,i) = sdiag(i)
      end do
c
c  Note that we do NOT need to transpose the V that comes out of LINPACK!
c
      return
      end
      subroutine pseudo_inverse ( m, n, u, s, v, a_pseudo )

c*********************************************************************72
c
cc PSEUDO_INVERSE computes the pseudoinverse.
c
c  Discussion:
c
c    Given the singular value decomposition of a real MxN matrix A:
c
c      A = U * S * V'
c
c    where
c
c      U is MxM orthogonal,
c      S is MxN, and entirely zero except for the diagonal;
c      V is NxN orthogonal.
c
c    the pseudo inverse is the NxM matrix A+ with the form
c
c      A+ = V * S+ * U'
c
c    where
c
c      S+ is the NxM matrix whose nonzero diagonal elements are
c      the inverses of the corresponding diagonal elements of S.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns
c    in the matrix A.
c
c    Input, double precision A(M,N), the matrix whose singular value
c    decomposition we are investigating.
c
c    Input, double precision U(M,M), S(M,N), V(N,N), the factors
c    that form the singular value decomposition of A.
c
c    Output, double precision A_PSEUDO(N,M), the pseudo_inverse of A.
c
      implicit none

      integer m
      integer n

      double precision a_pseudo(n,m)
      integer i
      integer j
      double precision s(m,n)
      double precision sp(n,m)
      double precision sput(n,m)
      double precision u(m,m)
      double precision ut(m,m)
      double precision v(n,n)

      do j = 1, m
        do i = 1, n
          sp(i,j) = 0.0D+00
        end do
      end do

      do i = 1, min ( m, n )
        if ( s(i,i) .ne. 0.0D+00 ) then
          sp(i,i) = 1.0D+00 / s(i,i)
        end if
      end do

      call r8mat_transpose ( m, m, u, ut )

      call r8mat_mm ( n, m, m, sp, ut, sput )

      call r8mat_mm ( n, n, m, v, sput, a_pseudo )

      return
      end
      subroutine pseudo_linear_solve_test ( m, n, a, a_pseudo, seed )

c*********************************************************************72
c
cc PSEUDO_LINEAR_SOLVE_TEST uses the pseudoinverse for linear systems.
c
c  Discussion:
c
c    Given an MxN matrix A, and its pseudoinverse A+:
c
c      "Solve" A  * x = b by x = A+  * b.
c
c      "Solve" A' * x = b by x = A+' * b.
c
c    When the system is overdetermined, the solution minimizes the
c    L2 norm of the residual.
c
c    When the system is underdetermined, the solution
c    is the solution of minimum L2 norm.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns
c    in the matrix A.
c
c    Input, double precision A(M,N), the matrix whose singular value
c    decomposition we are investigating.
c
c    Input, double precision A_PSEUDO(N,M), the pseudo_inverse of A.
c
c    Input/output, integer SEED, a seed for the random
c    number generator.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision a_pseudo(n,m)
      double precision apt(m,n)
      double precision at(n,m)
      double precision bm(m)
      double precision bn(n)
      integer i
      double precision r8vec_norm_l2
      double precision rm(m)
      double precision rn(n)
      integer seed
      double precision xm1(m)
      double precision xm2(m)
      double precision xn1(n)
      double precision xn2(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PSEUDO_LINEAR_SOLVE_TEST'
c
c  A * x = b, b in range of A.
c
      call r8vec_uniform_01 ( n, seed, xn1 )

      do i = 1, n
        xn1(i) = anint ( 10.0D+00 * xn1(i) )
      end do

      call r8mat_mv ( m, n, a, xn1, bm )
      call r8mat_mv ( n, m, a_pseudo, bm, xn2 )
      call r8mat_mv ( m, n, a, xn2, rm )
      do i = 1, m
        rm(i) = rm(i) - bm(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Given:'
      write ( *, '(a)' ) '    b = A * x1'
      write ( *, '(a)' ) '  so that b is in the range of A, solve'
      write ( *, '(a)' ) '    A * x = b'
      write ( *, '(a)' ) '  using the pseudoinverse:'
      write ( *, '(a)' ) '    x2 = A+ * b.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) 
     &  '  Norm of x1 = ', r8vec_norm_l2 ( n, xn1 )
      write ( *, '(a,g14.6)' ) 
     &  '  Norm of x2 = ', r8vec_norm_l2 ( n, xn2 )
      write ( *, '(a,g14.6)' ) 
     &  '  Norm of residual = ', r8vec_norm_l2 ( m, rm )
c
c  A * x = b, b not in range of A.
c
      if ( n .lt. m ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  For N < M, most systems A*x=b will not be'
        write ( *, '(a)' ) 
     &    '  exactly and uniquely solvable, except in the '
        write ( *, '(a)' ) '  least squares sense.'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Here is an example:'

        call r8vec_uniform_01 ( m, seed, bm )
        call r8mat_mv ( n, m, a_pseudo, bm, xn2 )
        call r8mat_mv ( m, n, a, xn2, rm )
        do i = 1, m
          rm(i) = rm(i) - bm(i)
        end do

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Given b is NOT in the range of A, solve'
        write ( *, '(a)' ) '    A * x = b'
        write ( *, '(a)' ) '  using the pseudoinverse:'
        write ( *, '(a)' ) '    x2 = A+ * b.'
        write ( *, '(a)' ) ' '
        write ( *, '(a,g14.6)' ) 
     &    '  Norm of x2 = ', r8vec_norm_l2 ( n, xn2 )
        write ( *, '(a,g14.6)' ) 
     &    '  Norm of residual = ', r8vec_norm_l2 ( m, rm )
      end if
c
c  A' * x = b, b is in the range of A'.
c
      call r8vec_uniform_01 ( m, seed, xm1 )
      do i = 1, m
        xm1(i) = anint ( 10.0D+00 * xm1(i) )
      end do
      call r8mat_transpose ( m, n, a, at )
      call r8mat_mv ( n, m, at, xm1, bn )
      call r8mat_transpose ( n, m, a_pseudo, apt )
      call r8mat_mv ( m, n, apt, bn, xm2 )
      call r8mat_mv ( n, m, at, xm2, rn )
      do i = 1, n
        rn(i) = rn(i) - bn(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Given:'
      write ( *, '(a)' ) '    b = A'' * x1'
      write ( *, '(a)' ) '  so that b is in the range of A'', solve'
      write ( *, '(a)' ) '    A'' * x = b'
      write ( *, '(a)' ) '  using the pseudoinverse:'
      write ( *, '(a)' ) '    x2 = A+'' * b.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) 
     &  '  Norm of x1 = ', r8vec_norm_l2 ( m, xm1 )
      write ( *, '(a,g14.6)' ) 
     &  '  Norm of x2 = ', r8vec_norm_l2 ( m, xm2 )
      write ( *, '(a,g14.6)' ) 
     &  '  Norm of residual = ', r8vec_norm_l2 ( n, rn )
c
c  A' * x = b, b is not in the range of A'.

      if ( m .lt. n ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    '  For M < N, most systems A''*x=b will not be'
        write ( *, '(a)' ) 
     &    '  exactly and uniquely solvable, except in the'
        write ( *, '(a)' ) '  least squares sense.'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Here is an example:'

        call r8vec_uniform_01 ( n, seed, bn )
        call r8mat_mv ( m, n, apt, bn, xm2 )
        call r8mat_mv ( n, m, at, xm2, rn )
        do i = 1, n
          rn(i) = rn(i) - bn(i)
        end do

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Given b is NOT in the range of A'', solve'
        write ( *, '(a)' ) '    A'' * x = b'
        write ( *, '(a)' ) '  using the pseudoinverse:'
        write ( *, '(a)' ) '    x2 = A+'' * b.'
        write ( *, '(a)' ) ' '
        write ( *, '(a,g14.6)' ) 
     &    '  Norm of x2 = ', r8vec_norm_l2 ( m, xm2 )
        write ( *, '(a,g14.6)' ) 
     &    '  Norm of residual = ', r8vec_norm_l2 ( n, rn )

      end if

      return
      end
      subroutine pseudo_product_test ( m, n, a, a_pseudo )

c*********************************************************************72
c
cc PSEUDO_PRODUCT_TEST examines pseudoinverse products.
c
c  Discussion:
c
c    Given an MxN matrix A, and its pseudoinverse A+, we must have
c
c      A+ * A * A+ = A+
c      A * A+ * A = A
c      ( A * A+ )' = A * A+ (MxM symmetry)
c      ( A+ * A )' = A+ * A (NxN symmetry)
c
c    If M <= N, A * A+ may be "interesting" (equal to or "like" the identity),
c    if N <= M, A+ * A may be "interesting" (equal to or "like" the identity).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns
c    in the matrix A.
c
c    Input, double precision A(M,N), the matrix whose singular value
c    decomposition we are investigating.
c
c    Input, double precision A_PSEUDO(N,M), the pseudo_inverse of A.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision a_pseudo(n,m)
      double precision bmm(m,m)
      double precision bmmt(m,m)
      double precision bmn(m,n)
      double precision bnm(n,m)
      double precision bnn(n,n)
      double precision bnnt(n,n)
      double precision dif1
      double precision dif2
      double precision dif3
      double precision dif4
      integer i
      double precision r8mat_dif_fro

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PSEUDO_PRODUCT_TEST'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The following relations MUST hold:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   A  * A+ * A  = A'
      write ( *, '(a)' ) '   A+ * A  * A+ = A+'
      write ( *, '(a)' ) ' ( A  * A+ ) is MxM symmetric;'
      write ( *, '(a)' ) ' ( A+ * A  ) is NxN symmetric'

      call r8mat_mm ( n, m, n, a_pseudo, a, bnn )
      call r8mat_mm ( m, n, n, a, bnn, bmn )

      call r8mat_mm ( m, n, m, a, a_pseudo, bmm )
      call r8mat_mm ( n, m, m, a_pseudo, bmm, bnm )

      call r8mat_transpose ( m, m, bmm, bmmt )
      call r8mat_transpose ( n, n, bnn, bnnt )

      dif1 = r8mat_dif_fro ( m, n, a, bmn )
      dif2 = r8mat_dif_fro ( n, m, a_pseudo, bnm )
      dif3 = r8mat_dif_fro ( m, m, bmm, bmmt )
      dif4 = r8mat_dif_fro ( n, n, bnn, bnnt )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Here are the Frobenius norms of the errors'
      write ( *, '(a)' ) '  in these relationships:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '   A  * A+ * A  = A            ', dif1
      write ( *, '(a,g14.6)' ) '   A+ * A  * A+ = A+           ', dif2
      write ( *, '(a,g14.6)' ) ' ( A  * A+ ) is MxM symmetric; ', dif3
      write ( *, '(a,g14.6)' ) ' ( A+ * A  ) is NxN symmetric; ', dif4

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  In some cases, the matrix A * A+'
      write ( *, '(a)' ) '  may be interesting (if M <= N, then'
      write ( *, '(a)' ) '  it MIGHT look like the identity.)'
      write ( *, '(a)' ) ' '

      call r8mat_print ( m, m, bmm, '  A * A+:' )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  In some cases, the matrix A+ * A'
      write ( *, '(a)' ) '  may be interesting (if N <= M, then'
      write ( *, '(a)' ) '  it MIGHT look like the identity.)'
      write ( *, '(a)' ) ' '

      call r8mat_print ( n, n, bnn, '  A+ * A:' )

      return
      end
      function r8mat_dif_fro ( m, n, a, b )

c*********************************************************************72
c
cc R8MAT_DIF_FRO returns the Frobenius norm of the difference of two R8MAT's.
c
c  Discussion:
c
c    An R8MAT is a matrix of double precision values.
c
c    The Frobenius norm is defined as
c
c      R8MAT_DIF_FRO = sqrt (
c        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) ( A(I,J) - B(I,J) )^2 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows in A and B.
c
c    Input, integer N, the number of columns in A and B.
c
c    Input, double precision A(M,N), B(M,N), the matrices
c    for which we want the Frobenius norm of the difference.
c
c    Output, double precision R8MAT_DIF_FRO, the Frobenius norm of
c    the difference of A and B.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision b(m,n)
      integer i
      integer j
      double precision r8mat_dif_fro
      double precision value

      value = 0.0D+00
      do j = 1, n
        do i = 1, m
          value = value + ( a(i,j) - b(i,j) )**2
        end do
      end do

      value = sqrt ( value )

      r8mat_dif_fro = value

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
      subroutine r8mat_mv ( m, n, a, x, y )

c*********************************************************************72
c
cc R8MAT_MV multiplies a matrix times a vector.
c
c  Discussion:
c
c    An R8MAT is an array of R8's.
c
c    In FORTRAN90, this operation can be more efficiently carried
c    out by the command
c
c      Y(1:M) = MATMUL ( A(1:M,1:N), X(1:N) )
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
c    Input, integer M, N, the number of rows and columns of the matrix.
c
c    Input, double precision A(M,N), the M by N matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision Y(M), the product A*X.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer i
      integer j
      double precision x(n)
      double precision y(m)

      do i = 1, m
        y(i) = 0.0D+00
        do j = 1, n
          y(i) = y(i) + a(i,j) * x(j)
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
c        sum ( 1 .le. I .le. M ) sum ( 1 .le. j .le. N ) A(I,J)^2 )
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
      subroutine r8mat_transpose ( m, n, a, at )

c*********************************************************************72
c
cc R8MAT_TRANSPOSE makes a transposed copy of a matrix.
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
c    13 June 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns of the matrix A.
c
c    Input, double precision A(N,N), the matrix to be transposed.
c
c    Output, double precision AT(N,M), the matrix to be transposed.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision at(n,m)
      integer i
      integer j

      do j = 1, m
        do i = 1, n
          at(i,j) = a(j,i)
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
      function r8vec_norm_l2 ( n, a )

c*********************************************************************72
c
cc R8VEC_NORM_L2 returns the L2 norm of an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8 values.
c
c    The vector L2 norm is defined as:
c
c      R8VEC_NORM_L2 = sqrt ( sum ( 1 <= I <= N ) A(I)**2 ).
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
c    Input, integer N, the number of entries in A.
c
c    Input, double precision A(N), the vector whose L2 norm is desired.
c
c    Output, double precision R8VEC_NORM_L2, the L2 norm of A.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      double precision r8vec_norm_l2
      double precision value

      value = 0.0D+00
      do i = 1, n
        value = value + a(i) * a(i)
      end do
      value = sqrt ( value )

      r8vec_norm_l2 = value

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
      subroutine rank_one_print_test ( m, n, a, u, s, v )

c*********************************************************************72
c
cc RANK_ONE_PRINT_TEST prints the sums of rank one matrices.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns
c    in the matrix A.
c
c    Input, double precision A(M,N), the matrix whose singular value
c    decomposition we are investigating.
c
c    Input, double precision U(M,M), S(M,N), V(N,N), the factors
c    that form the singular value decomposition of A.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision a_norm
      double precision dif_norm
      integer i
      integer ij
      integer j
      integer r
      double precision r8mat_dif_fro
      double precision r8mat_norm_fro
      double precision s(m,n)
      double precision svt(m,n)
      character ( len = 80 ) title
      double precision u(m,m)
      double precision usvt(m,n)
      double precision v(n,n)

      a_norm = r8mat_norm_fro ( m, n, a )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RANK_ONE_PRINT_TEST:'
      write ( *, '(a)' ) '  Print the sums of R rank one matrices.'

      do r = 0, min ( m, n )

        do j = 1, n
          do i = 1, r
            svt(i,j) = 0.0D+00
            do ij = 1, r
              svt(i,j) = svt(i,j) + s(i,ij) * v(j,ij)
            end do
          end do
        end do

        do j = 1, n
          do i = 1, m
            usvt(i,j) = 0.0D+00
            do ij = 1, r
              usvt(i,j) = usvt(i,j) + u(i,ij) * svt(ij,j)
            end do
          end do
        end do

        write ( title, '(a,i8)' ) '  Rank R = ', r
        call r8mat_print ( m, n, usvt, title )

      end do

      call r8mat_print ( m, n, a, '  Original matrix A:' )

      return
      end
      subroutine rank_one_test ( m, n, a, u, s, v )

c*********************************************************************72
c
cc RANK_ONE_TEST compares A to the sum of rank one matrices.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns
c    in the matrix A.
c
c    Input, double precision A(M,N), the matrix whose singular value
c    decomposition we are investigating.
c
c    Input, double precision U(M,M), S(M,N), V(N,N), the factors
c    that form the singular value decomposition of A.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision a_norm
      double precision dif_norm
      integer i
      integer ij
      integer j
      integer r
      double precision r8mat_dif_fro
      double precision r8mat_norm_fro
      double precision s(m,n)
      double precision svt(n,n)
      double precision u(m,m)
      double precision usvt(m,n)
      double precision v(n,n)

      a_norm = r8mat_norm_fro ( m, n, a )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'RANK_ONE_TEST:'
      write ( *, '(a)' ) 
     &  '  Compare A to the sum of R rank one matrices.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         R    Absolute      Relative'
      write ( *, '(a)' ) '              Error         Error'
      write ( *, '(a)' ) ' '

      do r = 0, min ( m, n )

        do j = 1, n
          do i = 1, r
            svt(i,j) = 0.0D+00
            do ij = 1, r
              svt(i,j) = svt(i,j) + s(i,ij) * v(j,ij)
            end do
          end do
        end do

        do j = 1, n
          do i = 1, m
            usvt(i,j) = 0.0D+00
            do ij = 1, r
              usvt(i,j) = usvt(i,j) + u(i,ij) * svt(ij,j)
            end do
          end do
        end do

        dif_norm = r8mat_dif_fro ( m, n, a, usvt )

        write ( *, '(2x,i8,2x,g14.6,2x, g14.6)' ) 
     &    r, dif_norm, dif_norm / a_norm

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
      subroutine s_to_i4 ( s, ival, ierror, length )

c*********************************************************************72
c
cc S_TO_I4 reads an I4 from a string.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) S, a string to be examined.
c
c    Output, integer IVAL, the integer value read from the string.
c    If the string is blank, then IVAL will be returned 0.
c
c    Output, integer IERROR, an error flag.
c    0, no error.
c    1, an error occurred.
c
c    Output, integer LENGTH, the number of characters of S
c    used to make IVAL.
c
      implicit none

      character c
      integer i
      integer ierror
      integer isgn
      integer istate
      integer ival
      integer length
      character * ( * ) s
      integer s_len_trim

      ierror = 0
      istate = 0
      isgn = 1
      ival = 0

      do i = 1, s_len_trim ( s )

        c = s(i:i)
c
c  Haven't read anything.
c
        if ( istate .eq. 0 ) then

          if ( c .eq. ' ' ) then

          else if ( c .eq. '-' ) then
            istate = 1
            isgn = -1
          else if ( c .eq. '+' ) then
            istate = 1
            isgn = + 1
          else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
            istate = 2
            ival = ichar ( c ) - ichar ( '0' )
          else
            ierror = 1
            return
          end if
c
c  Have read the sign, expecting digits.
c
        else if ( istate .eq. 1 ) then

          if ( c .eq. ' ' ) then

          else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
            istate = 2
            ival = ichar ( c ) - ichar ( '0' )
          else
            ierror = 1
            return
          end if
c
c  Have read at least one digit, expecting more.
c
        else if ( istate .eq. 2 ) then

          if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
            ival = 10 * ival + ichar ( c ) - ichar ( '0' )
          else
            ival = isgn * ival
            length = i - 1
            return
          end if

        end if

      end do
c
c  If we read all the characters in the string, see if we're OK.
c
      if ( istate .eq. 2 ) then
        ival = isgn * ival
        length = s_len_trim ( s )
      else
        ierror = 1
        length = 0
      end if

      return
      end
      subroutine svd_product_test ( m, n, a, u, s, v )

c*********************************************************************72
c
cc SVD_PRODUCT_TEST tests that A = U * S * V'.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns
c    in the matrix A.
c
c    Input, double precision A(M,N), the matrix whose singular value
c    decomposition we are investigating.
c
c    Input, double precision U(M,M), S(M,N), V(N,N), the factors
c    that form the singular value decomposition of A.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision a_norm
      double precision dif_norm
      integer i
      double precision r8mat_dif_fro
      double precision r8mat_norm_fro
      double precision s(m,n)
      double precision svt(m,n)
      double precision u(m,m)
      double precision usvt(m,n)
      double precision v(n,n)
      double precision vt(n,n)

      a_norm = r8mat_norm_fro ( m, n, a )

      call r8mat_transpose ( n, n, v, vt )
      call r8mat_mm ( m, n, n, s, vt, svt )
      call r8mat_mm ( m, m, n, u, svt, usvt )

      call r8mat_print ( m, n, usvt, '  The product U * S * V'':' )

      dif_norm = r8mat_dif_fro ( m, n, a, usvt )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) 
     &  '  Frobenius Norm of A, A_NORM = ', a_norm
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  ABSOLUTE ERROR for A = U*S*V'':'
      write ( *, '(a,g14.6)' ) 
     &  '  Frobenius norm of difference A-U*S*V'' = ', dif_norm
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  RELATIVE ERROR for A = U*S*V'':'
      write ( *, '(a,g14.6)' ) '  Ratio of DIF_NORM / A_NORM = ', 
     &  dif_norm / a_norm

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
