      subroutine c83_cr_fa ( n, a, a_cr )

c*********************************************************************72
c
cc C83_CR_FA decomposes a C83 matrix using cyclic reduction.
c
c  Discussion:
c
c    The D83 storage format is used for a real tridiagonal matrix.
c    The superdiagonal is stored in entries (1,2:N), the diagonal in
c    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
c    original matrix is "collapsed" vertically into the array.
c
c    Once C83_CR_FA has decomposed a matrix A, then C83_CR_SL may be used 
c    to solve linear systems A * x = b.
c
c    C83_CR_FA does not employ pivoting.  Hence, the results can be more
c    sensitive to ill-conditioning than standard Gauss elimination.  In
c    particular, C83_CR_FA will fail if any diagonal element of the matrix
c    is zero.  Other matrices may also cause C83_CR_FA to fail.
c
c    C83_CR_FA can be guaranteed to work properly if the matrix is strictly
c    diagonally dominant, that is, if the absolute value of the diagonal
c    element is strictly greater than the sum of the absolute values of
c    the offdiagonal elements, for each equation.
c
c    The algorithm may be illustrated by the following figures:
c
c    The initial matrix is given by:
c
c          D1 U1
c          L1 D2 U2
c             L2 D3 U3
c                L3 D4 U4
c                   L4 D5 U5
c                      L5 D6
c
c    Rows and columns are permuted in an odd/even way to yield:
c
c          D1       U1
c             D3    L2 U3
c                D5    L4 U5
c          L1 U2    D2
c             L3 U4    D4
c                L5       D6
c
c    A block LU decomposition is performed to yield:
c
c          D1      |U1
c             D3   |L2 U3
c                D5|   L4 U5
c          --------+--------
c                  |D2'F3
c                  |F1 D4'F4
c                  |   F2 D6'
c
c    For large systems, this reduction is repeated on the lower right hand
c    tridiagonal subsystem until a completely upper triangular system
c    is obtained.  The system has now been factored into the product of a
c    lower triangular system and an upper triangular one, and the information
c    defining this factorization may be used by C83_CR_SL to solve linear
c    systems.
c
c  Example:
c
c    Here is how a C83 matrix of order 5 would be stored:
c
c       *  A12 A23 A34 A45
c      A11 A22 A33 A44 A55
c      A21 A32 A43 A54  *
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 May 2009
c
c  Author:
c
c    FORTRAN77 version by John Burkardt
c
c  Reference:
c
c    Roger Hockney,
c    A fast direct solution of Poisson's equation using Fourier Analysis,
c    Journal of the ACM,
c    Volume 12, Number 1, pages 95-113, January 1965.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input, double complex A(3,N), the matrix.
c
c    Output, double complex A_CR(3,0:2*N), factorization information 
c    needed by C83_CR_SL.
c
      implicit none

      integer n

      double complex a(3,n)
      double complex a_cr(3,0:2*n)
      integer iful
      integer ifulp
      integer ihaf
      integer il
      integer ilp
      integer inc
      integer incr
      integer ipnt
      integer ipntp
      integer j

      if ( n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'C83_CR_FA - Fatal error!'
        write ( *, '(a,i8)' ) '  Nonpositive N = ', n
        stop
      end if

      if ( n .eq. 1 ) then
        a_cr(1,0) = 0.0D+00
        a_cr(1,1) = 0.0D+00
        a_cr(1,2) = 0.0D+00
        a_cr(2,0) = 0.0D+00
        a_cr(2,1) = 1.0D+00 / a(2,1)
        a_cr(2,2) = 0.0D+00
        a_cr(3,0) = 0.0D+00
        a_cr(3,1) = 0.0D+00
        a_cr(3,2) = 0.0D+00
        return
      end if
c
c  Zero out the workspace entries.
c
      a_cr(1,0) = 0.0D+00
      do j = 1, n - 1
        a_cr(1,j) = a(1,j+1)
      end do
      do j = n, 2 * n
        a_cr(1,j) = 0.0D+00
      end do

      a_cr(2,0) = 0.0D+00
      do j = 1, n
        a_cr(2,j) = a(2,j)
      end do
      do j = n + 1, 2 * n
        a_cr(2,j) = 0.0D+00
      end do

      a_cr(3,0) = 0.0D+00
      do j = 1, n - 1
        a_cr(3,j) = a(3,j)
      end do
      do j = n, 2 * n
        a_cr(3,j) = 0.0D+00
      end do

      il = n
      ipntp = 0

10    continue

      if ( 1 < il ) then

        ipnt = ipntp
        ipntp = ipntp + il
        if ( mod ( il, 2 ) .eq. 1 ) then
          inc = il + 1
        else
          inc = il
        end if

        incr = inc / 2
        il = il / 2
        ihaf = ipntp + incr + 1
        ifulp = ipnt + inc + 2

        do ilp = incr, 1, -1
          ifulp = ifulp - 2
          iful = ifulp - 1
          ihaf = ihaf - 1
          a_cr(2,iful) = 1.0D+00 / a_cr(2,iful)
          a_cr(3,iful)  = a_cr(3,iful)  * a_cr(2,iful)
          a_cr(1,ifulp) = a_cr(1,ifulp) * a_cr(2,ifulp+1)
          a_cr(2,ihaf)  = a_cr(2,ifulp) - a_cr(1,iful)  * a_cr(3,iful) 
     &                                  - a_cr(1,ifulp) * a_cr(3,ifulp)
          a_cr(3,ihaf) = -a_cr(3,ifulp) * a_cr(3,ifulp+1)
          a_cr(1,ihaf) = -a_cr(1,ifulp) * a_cr(1,ifulp+1)
        end do

        go to 10

      end if

      a_cr(2,ipntp+1) = 1.0D+00 / a_cr(2,ipntp+1)

      return
      end
      subroutine c83_cr_sl ( n, a_cr, b, x )

c*********************************************************************72
c
cc C83_CR_SL solves a linear system factored by C83_CR_FA.
c
c  Discussion:
c
c    The matrix A must be tridiagonal.  C83_CR_FA is called to compute the
c    LU factors of A.  It does so using a form of cyclic reduction.  If
c    the factors computed by C83_CR_FA are passed to C83_CR_SL, then 
c    a linear system involving the matrix A may be solved.
c
c    Note that C83_CR_FA does not perform pivoting, and so the solution 
c    produced by C83_CR_SL may be less accurate than a solution produced 
c    by a standard Gauss algorithm.  However, such problems can be 
c    guaranteed not to occur if the matrix A is strictly diagonally 
c    dominant, that is, if the absolute value of the diagonal coefficient 
c    is greater than the sum of the absolute values of the two off diagonal 
c    coefficients, for each row of the matrix.
c
c  Example:
c
c    Here is how a C83 matrix of order 5 would be stored:
c
c       *  A12 A23 A34 A45
c      A11 A22 A33 A44 A55
c      A21 A32 A43 A54  *
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 May 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Roger Hockney,
c    A fast direct solution of Poisson's equation using Fourier Analysis,
c    Journal of the ACM,
c    Volume 12, Number 1, pages 95-113, January 1965.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input, double complex A_CR(3,0:2*N), factorization information 
c    computed by C83_CR_FA.
c
c    Input, real ( kind = 8 ) B(N), the right hand sides.
c
c    Output, real ( kind = 8 ) X(N), the solutions of the linear systems.
c
      implicit none

      integer n
      integer nb

      double complex a_cr(3,0:2*n)
      double complex b(n)
      integer i
      integer iful
      integer ifulm
      integer ihaf
      integer il
      integer ipnt
      integer ipntp
      integer ndiv
      double complex rhs(0:2*n)
      double complex x(n)

      if ( n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'C83_CR_SLS - Fatal error!'
        write ( *, '(a,i8)' ) '  Nonpositive N = ', n
        stop
      end if

      if ( n .eq. 1 ) then
        x(1) = a_cr(2,1) * b(1)
        return
      end if
c
c  Set up RHS.
c
      rhs(0) = 0.0D+00
      do i = 1, n
        rhs(i) = b(i)
      end do
      do i = n + 1, 2 * n
        rhs(i) = 0.0D+00
      end do

      il = n
      ndiv = 1
      ipntp = 0

10    continue

      if ( 1 .lt. il ) then

        ipnt = ipntp
        ipntp = ipntp + il
        il = il / 2
        ndiv = ndiv * 2
        ihaf = ipntp

        do iful = ipnt + 2, ipntp, 2
          ihaf = ihaf + 1
          rhs(ihaf) = rhs(iful) 
     &      - a_cr(3,iful-1) * rhs(iful-1) 
     &      - a_cr(1,iful)   * rhs(iful+1)
        end do

        go to 10

      end if

      rhs(ihaf) = a_cr(2,ihaf) * rhs(ihaf)

      ipnt = ipntp

20    continue

      if ( 0 .lt. ipnt ) then

        ipntp = ipnt
        ndiv = ndiv / 2
        il = n / ndiv
        ipnt = ipnt - il
        ihaf = ipntp

        do ifulm = ipnt + 1, ipntp, 2
          iful = ifulm + 1
          ihaf = ihaf + 1
          rhs(iful) = rhs(ihaf)
          rhs(ifulm) = a_cr(2,ifulm) 
     &      * (                     rhs(ifulm) 
     &          - a_cr(3,ifulm-1) * rhs(ifulm-1) 
     &          - a_cr(1,ifulm)   * rhs(iful) )
        end do

        go to 20

      end if

      do i = 1, n
        x(i) = rhs(i)
      end do

      return
      end
      subroutine c83_cr_sls ( n, a_cr, nb, b, x )

c*********************************************************************72
c
cc C83_CR_SLS solves several linear systems factored by C83_CR_FA.
c
c  Discussion:
c
c    The matrix A must be tridiagonal.  C83_CR_FA is called to compute the
c    LU factors of A.  It does so using a form of cyclic reduction.  If
c    the factors computed by C83_CR_FA are passed to C83_CR_SLS, then one or many
c    linear systems involving the matrix A may be solved.
c
c    Note that C83_CR_FA does not perform pivoting, and so the solution 
c    produced by C83_CR_SLS may be less accurate than a solution produced 
c    by a standard Gauss algorithm.  However, such problems can be 
c    guaranteed not to occur if the matrix A is strictly diagonally 
c    dominant, that is, if the absolute value of the diagonal coefficient 
c    is greater than the sum of the absolute values of the two off diagonal 
c    coefficients, for each row of the matrix.
c
c  Example:
c
c    Here is how a C83 matrix of order 5 would be stored:
c
c       *  A12 A23 A34 A45
c      A11 A22 A33 A44 A55
c      A21 A32 A43 A54  *
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 May 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Roger Hockney,
c    A fast direct solution of Poisson's equation using Fourier Analysis,
c    Journal of the ACM,
c    Volume 12, Number 1, pages 95-113, January 1965.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input, double complex A_CR(3,0:2*N), factorization information 
c    computed by C83_CR_FA.
c
c    Input, integer NB, the number of right hand sides.
c
c    Input, real ( kind = 8 ) B(N,NB), the right hand sides.
c
c    Output, real ( kind = 8 ) X(N,NB), the solutions of the linear systems.
c
      implicit none

      integer n
      integer nb

      double complex a_cr(3,0:2*n)
      double complex b(n,nb)
      integer i
      integer iful
      integer ifulm
      integer ihaf
      integer il
      integer ipnt
      integer ipntp
      integer j
      integer ndiv
      double complex rhs(0:2*n,nb)
      double complex x(n,nb)

      if ( n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'C83_CR_SLS - Fatal error!'
        write ( *, '(a,i8)' ) '  Nonpositive N = ', n
        stop
      end if

      if ( n .eq. 1 ) then
        do j = 1, nb
          x(1,j) = a_cr(2,1) * b(1,j)
        end do
        return
      end if
c
c  Set up RHS.
c
      do j = 1, nb
        rhs(0,j) = 0.0D+00
      end do
      do j = 1, nb
        do i = 1, n
          rhs(i,j) = b(i,j)
        end do
      end do
      do j = 1, nb
        do i = n + 1, 2 * n
          rhs(i,j) = 0.0D+00
        end do
      end do

      il = n
      ndiv = 1
      ipntp = 0

10    continue

      if ( 1 .lt. il ) then

        ipnt = ipntp
        ipntp = ipntp + il
        il = il / 2
        ndiv = ndiv * 2

        do j = 1, nb 
          ihaf = ipntp
          do iful = ipnt + 2, ipntp, 2
            ihaf = ihaf + 1
            rhs(ihaf,j) = rhs(iful,j) 
     &        - a_cr(3,iful-1) * rhs(iful-1,j) 
     &        - a_cr(1,iful)   * rhs(iful+1,j)
          end do
        end do

        go to 10

      end if

      do j = 1, nb
        rhs(ihaf,j) = a_cr(2,ihaf) * rhs(ihaf,j)
      end do

      ipnt = ipntp

20    continue

      if ( 0 .lt. ipnt ) then

        ipntp = ipnt
        ndiv = ndiv / 2
        il = n / ndiv
        ipnt = ipnt - il
 
        do j = 1, nb
          ihaf = ipntp
          do ifulm = ipnt + 1, ipntp, 2
            iful = ifulm + 1
            ihaf = ihaf + 1
            rhs(iful,j) = rhs(ihaf,j)
            rhs(ifulm,j) = a_cr(2,ifulm) 
     &        * (                     rhs(ifulm,j) 
     &            - a_cr(3,ifulm-1) * rhs(ifulm-1,j) 
     &            - a_cr(1,ifulm)   * rhs(iful,j) )
          end do
        end do

        go to 20

      end if

      do j = 1, nb
        do i = 1, n
          x(i,j) = rhs(i,j)
        end do
      end do

      return
      end
      subroutine c83_indicator ( n, a )

c*********************************************************************72
c
cc C83_INDICATOR sets up a C83 indicator matrix.
c
c  Discussion:
c
c    The C83 storage format is used for a tridiagonal matrix.
c    The superdiagonal is stored in entries (1,2:N), the diagonal in
c    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
c    original matrix is "collapsed" vertically into the array.
c
c  Example:
c
c    Here is how a C83 matrix of order 5 would be stored:
c
c       *  A12 A23 A34 A45
c      A11 A22 A33 A44 A55
c      A21 A32 A43 A54  *
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 May 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be at least 2.
c
c    Output, double complex A(3,N), the indicator matrix.
c
      implicit none

      integer n

      double complex a(3,n)
      integer i
      integer j

      a(1,1) = 0.0D+00
      do j = 2, n
        i = j - 1
        a(1,j) = dcmplx ( i, j )
      end do

      do j = 1, n
        i = j
        a(2,j) = dcmplx ( i, j )
      end do

      do j = 1, n-1
        i = j + 1
        a(3,j) = dcmplx ( i, j )
      end do
      a(3,n) = 0.0D+00

      return
      end
      subroutine c83_mxv ( n, a, x, b )

c*********************************************************************72
c
cc C83_MXV multiplies a C83 matrix times a C8VEC.
c
c  Discussion:
c
c    The C83 storage format is used for a tridiagonal matrix.
c    The superdiagonal is stored in entries (1,2:N), the diagonal in
c    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
c    original matrix is "collapsed" vertically into the array.
c
c  Example:
c
c    Here is how a C83 matrix of order 5 would be stored:
c
c       *  A12 A23 A34 A45
c      A11 A22 A33 A44 A55
c      A21 A32 A43 A54  *
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 May 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the linear system.
c
c    Input, double complex A(3,N), the matrix.
c
c    Input, double complex X(N), the vector to be multiplied by A.
c
c    Output, double complex B(N), the product A * x.
c
      implicit none

      integer n

      double complex a(3,n)
      double complex b(n)
      integer j
      double complex x(n)

      do j = 1, n
        b(j)   =      a(2,j)   * x(j)
      end do
      do j = 1, n - 1
        b(j) = b(j) + a(1,j+1) * x(j+1)
      end do
      do j = 2, n
        b(j) = b(j) + a(3,j-1) * x(j-1)
      end do

      return
      end
      subroutine c83_print ( n, a, title )

c*********************************************************************72
c
cc C83_PRINT prints a C83 matrix.
c
c  Discussion:
c
c    The C83 storage format is used for a complex tridiagonal matrix.
c    The superdiagonal is stored in entries (1,2:N), the diagonal in
c    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
c    original matrix is "collapsed" vertically into the array.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 May 2009
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
c    Input, double complex A(3,N), the C83 matrix.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n

      double complex a(3,n)
      character * ( * )  title

      call c83_print_some ( n, a, 1, 1, n, n, title )

      return
      end
      subroutine c83_print_some ( n, a, ilo, jlo, ihi, jhi, title )

c*********************************************************************72
c
cc C83_PRINT_SOME prints some of a C83 matrix.
c
c  Discussion:
c
c    The C83 storage format is used for a complex tridiagonal matrix.
c    The superdiagonal is stored in entries (1,2:N), the diagonal in
c    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
c    original matrix is "collapsed" vertically into the array.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 May 2009
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
c    Input, double complex A(3,N), the C83 matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column, to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 3 )
      integer n

      double complex a(3,n)
      character * ( 12 ) citemp(incx)
      character * ( 12 ) crtemp(incx)
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
      double precision xi
      double precision xr
      character * ( * )  title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
c
c  Print the columns of the matrix, in strips of 5.
c
      do j2lo = jlo, jhi, incx

        j2hi = j2lo + incx - 1
        j2hi = min ( j2hi, n )
        j2hi = min ( j2hi, jhi )

        inc = j2hi + 1 - j2lo

        write ( *, '(a)' ) ' '

        do j = j2lo, j2hi
          j2 = j + 1 - j2lo
          write ( crtemp(j2), '(i8,6x)' ) j
          write ( citemp(j2), '(i8,6x)' ) j
        end do

        write ( *, '(''  Col:  '',6a12)' ) 
     &    ( crtemp(j2), citemp(j2), j2 = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) '  ---'
c
c  Determine the range of the rows in this strip.
c
        i2lo = max ( ilo, 1 )
        i2lo = max ( i2lo, j2lo - 1 )
        i2hi = min ( ihi, n )
        i2hi = min ( i2hi, j2hi + 1 )

        do i = i2lo, i2hi
c
c  Print out (up to) INCX entries in row I, that lie in the current strip.
c
          do j2 = 1, inc

            j = j2lo - 1 + j2

            if ( 1 .lt. i - j .or. 1 .lt. j - i ) then

              crtemp(j2) = ' '
              citemp(j2) = ' '

            else

              if ( j .eq. i - 1 ) then
                xr = dble ( a(1,i) )
                xi = dimag ( a(1,i) )
              else if ( j .eq. i ) then
                xr = dble ( a(2,i) )
                xi = dimag ( a(2,i) )
              else if ( j .eq. i + 1 ) then
                xr = dble ( a(3,i) )
                xi = dimag ( a(3,i) )
              end if

              if ( xr .eq. 0.0D+00 .and. xi .eq. 0.0D+00 ) then
                crtemp(j2) = '    0.0'
                citemp(j2) = ' '
              else if ( xr .eq. 0.0D+00 .and. xi .ne. 0.0D+00 ) then
                crtemp(j2) = ' '
                write ( citemp(j2), '(g12.5)' ) xi
              else if ( xr .ne. 0.0D+00 .and. xi .eq. 0.0D+00 ) then
                write ( crtemp(j2), '(g12.5)' ) xr
                citemp(j2) = ' '
              else
                write ( crtemp(j2), '(g12.5)' ) xr
                write ( citemp(j2), '(g12.5)' ) xi
              end if

            end if

          end do

          write ( *, '(i5,a,6a12)' ) 
     &      i, ':', ( crtemp(j2), citemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine c8mat_print ( m, n, a, title )

c*********************************************************************72
c
cc C8MAT_PRINT prints a C8MAT.
c
c  Discussion:
c
c    A C8MAT is a matrix of C8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 December 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns 
c    in the matrix.
c
c    Input, double complex A(M,N), the matrix.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer m
      integer n

      double complex a(m,n)
      character * ( * ) title

      call c8mat_print_some ( m, n, a, 1, 1, m, n, title )

      return
      end
      subroutine c8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, 
     &  title )

c*********************************************************************72
c
cc C8MAT_PRINT_SOME prints some of a C8MAT.
c
c  Discussion:
c
c    A C8MAT is a matrix of C8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 March 2005
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns 
c    in the matrix.
c
c    Input, double complex A(M,N), the matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 4 )
      integer m
      integer n

      double complex a(m,n)
      character * ( 20 ) ctemp(incx)
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
      double complex zero

      zero = dcmplx ( 0.0D+00, 0.0D+00 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
c
c  Print the columns of the matrix, in strips of INCX.
c
      do j2lo = jlo, jhi, incx

        j2hi = j2lo + incx - 1
        j2hi = min ( j2hi, n )
        j2hi = min ( j2hi, jhi )

        inc = j2hi + 1 - j2lo

        write ( *, '(a)' ) ' '

        do j = j2lo, j2hi
          j2 = j + 1 - j2lo
          write ( ctemp(j2), '(i10,10x)' ) j
        end do

        write ( *, '(a,4a20)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) '  ---'
c
c  Determine the range of the rows in this strip.
c
        i2lo = max ( ilo, 1 )
        i2hi = min ( ihi, m )

        do i = i2lo, i2hi
c
c  Print out (up to) INCX entries in row I, that lie in the current strip.
c
          do j2 = 1, inc

            j = j2lo - 1 + j2

            if ( a(i,j) .eq. zero ) then
              ctemp(j2) = '       0.0          '
            else if ( dimag ( a(i,j) ) .eq. 0.0D+00 ) then
              write ( ctemp(j2), '(g10.3,10x)' ) dreal ( a(i,j) )
            else
              write ( ctemp(j2), '(2g10.3)' ) a(i,j)
            end if

          end do

          write ( *, '(i5,a,4a20)' ) i, ':', ( ctemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine c8vec_indicator ( n, a )

c*********************************************************************72
c
cc C8VEC_INDICATOR sets a C8VEC to the indicator vector.
c
c  Discussion:
c
c    A C8VEC is a vector of C8's
c
c    X(1:N) = ( 1-1i, 2-2i, 3-3i, 4-4i, ... )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 December 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of elements of A.
c
c    Output, double complex A(N), the array to be initialized.
c
      implicit none

      integer n

      double complex a(n)
      integer i

      do i = 1, n
        a(i) = dcmplx ( i, -i )
      end do

      return
      end
      subroutine c8vec_print ( n, a, title )

c*********************************************************************72
c
cc C8VEC_PRINT prints a C8VEC, with an optional title.
c
c  Discussion:
c
c    A C8VEC is a vector of C8's
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 December 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components of the vector.
c
c    Input, double complex A(N), the vector to be printed.
c
c    Input, character*(*) TITLE, a title.
c
      implicit none

      integer n

      double complex a(n)
      integer i
      character*(*) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,a,1x,2g14.6)' ) i, ':', a(i)
      end do

      return
      end
      subroutine c8vec_print_some ( n, x, i_lo, i_hi, title )

c*********************************************************************72
c
cc C8VEC_PRINT_SOME prints some of a C8VEC.
c
c  Discussion:
c
c    A C8VEC is a vector of C8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 October 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries of the vector.
c
c    Input, double complex X(N), the vector to be printed.
c
c    Input, integer I_LO, I_HI, the first and last entries
c    to print.
c
c    Input, character*(*) TITLE, a title.
c
      implicit none

      integer n

      integer i
      integer i_hi
      integer i_lo
      character*(*) title
      double complex x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '

      do i = max ( 1, i_lo ), min ( n, i_hi )
        write ( *, '(2x,i8,a,1x,2g14.6)' ) i, ':', x(i)
      end do

      return
      end
      subroutine r83_cr_fa ( n, a, a_cr )

c*********************************************************************72
c
cc R83_CR_FA decomposes a real tridiagonal matrix using cyclic reduction.
c
c  Discussion:
c
c    The R83 storage format is used for a real tridiagonal matrix.
c    The superdiagonal is stored in entries (1,2:N), the diagonal in
c    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
c    original matrix is "collapsed" vertically into the array.
c
c    Once R83_CR_FA has decomposed a matrix A, then R83_CR_SL may be used 
c    to solve linear systems A * x = b.
c
c    R83_CR_FA does not employ pivoting.  Hence, the results can be more
c    sensitive to ill-conditioning than standard Gauss elimination.  In
c    particular, R83_CR_FA will fail if any diagonal element of the matrix
c    is zero.  Other matrices may also cause R83_CR_FA to fail.
c
c    R83_CR_FA can be guaranteed to work properly if the matrix is strictly
c    diagonally dominant, that is, if the absolute value of the diagonal
c    element is strictly greater than the sum of the absolute values of
c    the offdiagonal elements, for each equation.
c
c    The algorithm may be illustrated by the following figures:
c
c    The initial matrix is given by:
c
c          D1 U1
c          L1 D2 U2
c             L2 D3 U3
c                L3 D4 U4
c                   L4 D5 U5
c                      L5 D6
c
c    Rows and columns are permuted in an odd/even way to yield:
c
c          D1       U1
c             D3    L2 U3
c                D5    L4 U5
c          L1 U2    D2
c             L3 U4    D4
c                L5       D6
c
c    A block LU decomposition is performed to yield:
c
c          D1      |U1
c             D3   |L2 U3
c                D5|   L4 U5
c          --------+--------
c                  |D2'F3
c                  |F1 D4'F4
c                  |   F2 D6'
c
c    For large systems, this reduction is repeated on the lower right hand
c    tridiagonal subsystem until a completely upper triangular system
c    is obtained.  The system has now been factored into the product of a
c    lower triangular system and an upper triangular one, and the information
c    defining this factorization may be used by R83_CR_SL to solve linear
c    systems.
c
c  Example:
c
c    Here is how an R83 matrix of order 5 would be stored:
c
c       *  A12 A23 A34 A45
c      A11 A22 A33 A44 A55
c      A21 A32 A43 A54  *
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 May 2009
c
c  Author:
c
c    FORTRAN77 version by John Burkardt
c
c  Reference:
c
c    Roger Hockney,
c    A fast direct solution of Poisson's equation using Fourier Analysis,
c    Journal of the ACM,
c    Volume 12, Number 1, pages 95-113, January 1965.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input, double precision A(3,N), the R83 matrix.
c
c    Output, double precision A_CR(3,0:2*N), factorization information 
c    needed by R83_CR_SL.
c
      implicit none

      integer n

      double precision a(3,n)
      double precision a_cr(3,0:2*n)
      integer iful
      integer ifulp
      integer ihaf
      integer il
      integer ilp
      integer inc
      integer incr
      integer ipnt
      integer ipntp

      if ( n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R83_CR_FA - Fatal error!'
        write ( *, '(a,i8)' ) '  Nonpositive N = ', n
        return
      end if

      if ( n .eq. 1 ) then
        a_cr(1,0:2) = 0.0D+00
        a_cr(2,0) = 0.0D+00
        a_cr(2,1) = 1.0D+00 / a(2,1)
        a_cr(2,2) = 0.0D+00
        a_cr(3,0:2) = 0.0D+00
        return
      end if
c
c  Zero out the workspace entries.
c
      a_cr(1,0) = 0.0D+00
      a_cr(1,1:n-1) = a(1,2:n)
      a_cr(1,n:2*n) = 0.0D+00

      a_cr(2,0) = 0.0D+00
      a_cr(2,1:n) = a(2,1:n)
      a_cr(2,n+1:2*n) = 0.0D+00

      a_cr(3,0) = 0.0D+00
      a_cr(3,1:n-1) = a(3,1:n-1)
      a_cr(3,n:2*n) = 0.0D+00

      il = n
      ipntp = 0

10    continue

      if ( 1 .lt. il ) then

        ipnt = ipntp
        ipntp = ipntp + il
        if ( mod ( il, 2 ) .eq. 1 ) then
          inc = il + 1
        else
          inc = il
        end if

        incr = inc / 2
        il = il / 2
        ihaf = ipntp + incr + 1
        ifulp = ipnt + inc + 2

        do ilp = incr, 1, -1
          ifulp = ifulp - 2
          iful = ifulp - 1
          ihaf = ihaf - 1
          a_cr(2,iful) = 1.0D+00 / a_cr(2,iful)
          a_cr(3,iful)  = a_cr(3,iful)  * a_cr(2,iful)
          a_cr(1,ifulp) = a_cr(1,ifulp) * a_cr(2,ifulp+1)
          a_cr(2,ihaf)  = a_cr(2,ifulp) - a_cr(1,iful) * a_cr(3,iful) 
     &                                - a_cr(1,ifulp) * a_cr(3,ifulp)
          a_cr(3,ihaf) = -a_cr(3,ifulp) * a_cr(3,ifulp+1)
          a_cr(1,ihaf) = -a_cr(1,ifulp) * a_cr(1,ifulp+1)
        end do

        go to 10

      end if

      a_cr(2,ipntp+1) = 1.0D+00 / a_cr(2,ipntp+1)

      return
      end
      subroutine r83_cr_sl ( n, a_cr, b, x )

c*********************************************************************72
c
cc R83_CR_SL solves a real linear system factored by R83_CR_FA.
c
c  Discussion:
c
c    The matrix A must be tridiagonal.  R83_CR_FA is called to compute the
c    LU factors of A.  It does so using a form of cyclic reduction.  If
c    the factors computed by R83_CR_FA are passed to R83_CR_SL, then one or many
c    linear systems involving the matrix A may be solved.
c
c    Note that R83_CR_FA does not perform pivoting, and so the solution 
c    produced by R83_CR_SL may be less accurate than a solution produced 
c    by a standard Gauss algorithm.  However, such problems can be 
c    guaranteed not to occur if the matrix A is strictly diagonally 
c    dominant, that is, if the absolute value of the diagonal coefficient 
c    is greater than the sum of the absolute values of the two off diagonal 
c    coefficients, for each row of the matrix.
c
c  Example:
c
c    Here is how an R83 matrix of order 5 would be stored:
c
c       *  A12 A23 A34 A45
c      A11 A22 A33 A44 A55
c      A21 A32 A43 A54  *
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 May 2010
c
c  Author:
c
c    FORTRAN77 version by John Burkardt
c
c  Reference:
c
c    Roger Hockney,
c    A fast direct solution of Poisson's equation using Fourier Analysis,
c    Journal of the ACM,
c    Volume 12, Number 1, pages 95-113, January 1965.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input, double precision A_CR(3,0:2*N), factorization information 
c    computed by R83_CR_FA.
c
c    Input, double precision B(N), the right hand side.
c
c    Output, double precision X(N), the solution of the linear system.
c
      implicit none

      integer n

      double precision a_cr(3,0:2*n)
      double precision b(n)
      integer i
      integer iful
      integer ifulm
      integer ihaf
      integer il
      integer ipnt
      integer ipntp
      integer ndiv
      double precision rhs(0:2*n)
      double precision x(n)

      if ( n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R83_CR_SL - Fatal error!'
        write ( *, '(a,i8)' ) '  Nonpositive N = ', n
        return
      end if

      if ( n .eq. 1 ) then
        x(1) = a_cr(2,1) * b(1)
        return
      end if
c
c  Set up RHS.
c
      rhs(0) = 0.0D+00
      do i = 1, n
        rhs(i) = b(i)
      end do
      do i = n + 1, 2 * n
        rhs(i) = 0.0D+00
      end do

      il = n
      ndiv = 1
      ipntp = 0

10    continue

      if ( 1 .lt. il ) then

        ipnt = ipntp
        ipntp = ipntp + il
        il = il / 2
        ndiv = ndiv * 2
        ihaf = ipntp

        do iful = ipnt + 2, ipntp, 2
          ihaf = ihaf + 1
          rhs(ihaf) = rhs(iful) - a_cr(3,iful-1) * rhs(iful-1) 
     &                          - a_cr(1,iful)   * rhs(iful+1)
        end do

        go to 10

      end if

      rhs(ihaf) = rhs(ihaf) * a_cr(2,ihaf)
      ipnt = ipntp

20    continue

      if ( 0 .lt. ipnt ) then

        ipntp = ipnt
        ndiv = ndiv / 2
        il = n / ndiv
        ipnt = ipnt - il
        ihaf = ipntp

        do ifulm = ipnt + 1, ipntp, 2
          iful = ifulm + 1
          ihaf = ihaf + 1
          rhs(iful) = rhs(ihaf)
          rhs(ifulm) = a_cr(2,ifulm) 
     &      * ( rhs(ifulm) - a_cr(3,ifulm-1) * rhs(ifulm-1) 
     &                     - a_cr(1,ifulm)   * rhs(iful) )
        end do

        go to 20

      end if

      do i = 1, n
        x(i) = rhs(i)
      end do

      return
      end
      subroutine r83_cr_sls ( n, a_cr, nb, b, x )

c*********************************************************************72
c
cc R83_CR_SLS solves several real linear systems factored by R83_CR_FA.
c
c  Discussion:
c
c    The matrix A must be tridiagonal.  R83_CR_FA is called to compute the
c    LU factors of A.  It does so using a form of cyclic reduction.  If
c    the factors computed by R83_CR_FA are passed to R83_CR_SLS, then one or many
c    linear systems involving the matrix A may be solved.
c
c    Note that R83_CR_FA does not perform pivoting, and so the solutions 
c    produced by R83_CR_SLS may be less accurate than a solution produced 
c    by a standard Gauss algorithm.  However, such problems can be 
c    guaranteed not to occur if the matrix A is strictly diagonally 
c    dominant, that is, if the absolute value of the diagonal coefficient 
c    is greater than the sum of the absolute values of the two off diagonal 
c    coefficients, for each row of the matrix.
c
c  Example:
c
c    Here is how an R83 matrix of order 5 would be stored:
c
c       *  A12 A23 A34 A45
c      A11 A22 A33 A44 A55
c      A21 A32 A43 A54  *
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 May 2010
c
c  Author:
c
c    FORTRAN77 version by John Burkardt
c
c  Reference:
c
c    Roger Hockney,
c    A fast direct solution of Poisson's equation using Fourier Analysis,
c    Journal of the ACM,
c    Volume 12, Number 1, pages 95-113, January 1965.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input, double precision A_CR(3,0:2*N), factorization information 
c    computed by R83_CR_FA.
c
c    Input, integer NB, the number of right hand sides.
c
c    Input, double precision B(N,NB), the right hand sides.
c
c    Output, double precision X(N,NB), the solutions of the linear system.
c
      implicit none

      integer n
      integer nb

      double precision a_cr(3,0:2*n)
      double precision b(n,nb)
      integer i
      integer iful
      integer ifulm
      integer ihaf
      integer il
      integer ipnt
      integer ipntp
      integer j
      integer ndiv
      double precision rhs(0:2*n,nb)
      double precision x(n,nb)

      if ( n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R83_CR_SLS - Fatal error!'
        write ( *, '(a,i8)' ) '  Nonpositive N = ', n
        return
      end if

      if ( n .eq. 1 ) then
        do j = 1, nb
          x(1,j) = a_cr(2,1) * b(1,j)
        end do
        return
      end if
c
c  Set up RHS.
c
      do j = 1, nb
        rhs(0,j) = 0.0D+00
        do i = 1, n
          rhs(i,j) = b(i,j)
        end do
        do i = n + 1, 2 * n
          rhs(i,j) = 0.0D+00
        end do
      end do

      il = n
      ndiv = 1
      ipntp = 0

10    continue

      if ( 1 .lt. il ) then

        ipnt = ipntp
        ipntp = ipntp + il
        il = il / 2
        ndiv = ndiv * 2

        do j = 1, nb
          ihaf = ipntp
          do iful = ipnt + 2, ipntp, 2
            ihaf = ihaf + 1
            rhs(ihaf,j) = rhs(iful,j) 
     &        - a_cr(3,iful-1) * rhs(iful-1,j) 
     &        - a_cr(1,iful)   * rhs(iful+1,j)
          end do
        end do

        go to 10

      end if

      do j = 1, nb
        rhs(ihaf,j) = a_cr(2,ihaf) * rhs(ihaf,j)
      end do

      ipnt = ipntp

20    continue

      if ( 0 .lt. ipnt ) then

        ipntp = ipnt
        ndiv = ndiv / 2
        il = n / ndiv
        ipnt = ipnt - il

        do j = 1, nb
          ihaf = ipntp
          do ifulm = ipnt + 1, ipntp, 2
            iful = ifulm + 1
            ihaf = ihaf + 1
            rhs(iful,j) = rhs(ihaf,j)
            rhs(ifulm,j) = a_cr(2,ifulm) 
     &        * (                     rhs(ifulm,j) 
     &            - a_cr(3,ifulm-1) * rhs(ifulm-1,j) 
     &            - a_cr(1,ifulm)   * rhs(iful,j) )
          end do

        end do

        go to 20

      end if

      do j = 1, nb
        do i = 1, n
          x(i,j) = rhs(i,j)
        end do
      end do

      return
      end
      subroutine r83_gs_sl ( n, a, b, x, tol, it_max, job, it, diff )

c*********************************************************************72
c
cc R83_GS_SL solves an R83 system using Gauss-Seidel iteration.
c
c  Discussion:
c
c    The R83 storage format is used for a tridiagonal matrix.
c    The superdiagonal is stored in entries (1,2:N), the diagonal in
c    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
c    original matrix is "collapsed" vertically into the array.
c
c    This routine simply applies a given number of steps of the
c    iteration to an input approximate solution.  On first call, you can
c    simply pass in the zero vector as an approximate solution.  If
c    the returned value is not acceptable, you may call again, using
c    it as the starting point for additional iterations.
c
c  Example:
c
c    Here is how an R83 matrix of order 5 would be stored:
c
c       *  A12 A23 A34 A45
c      A11 A22 A33 A44 A55
c      A21 A32 A43 A54  *
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be at least 2.
c
c    Input, double precision A(3,N), the R83 matrix.
c
c    Input, double precision B(N), the right hand side of the linear system.
c
c    Input/output, double precision X(N), an approximate solution to 
c    the system.
c
c    Input, double precision TOL, a tolerance.  If the maximum change in
c    the solution is less than TOL, the iteration is terminated early.
c
c    Input, integer IT_MAX, the maximum number of iterations.
c
c    Input, integer JOB, specifies the system to solve.
c    0, solve A * x = b.
c    nonzero, solve A' * x = b.
c
c    Output, integer IT, the number of iterations taken.
c
c    Output, double precision DIFF, the maximum change in the solution
c    on the last iteration.
c
      implicit none

      integer n

      double precision a(3,n)
      double precision b(n)
      double precision diff
      integer i
      integer it
      integer it_max
      integer it_num
      integer job
      double precision tol
      double precision x(n)
      double precision x_norm
      double precision x_old(n)
c
c  No diagonal matrix entry can be zero.
c
      do i = 1, n
        if ( a(2,i) .eq. 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R83_GS_SL - Fatal error!'
          write ( *, '(a,i8)' ) '  Zero diagonal entry, index = ', i
          return
        end if
      end do

      if ( job .eq. 0 ) then

        do it_num = 1, it_max

          it = it_num

          x_old(1:n) = x(1:n)

          x(1) =   ( b(1)                   - a(3,1) * x(2)   ) / a(2,1)
          do i = 2, n - 1
            x(i) = ( b(i) - a(1,i) * x(i-1) - a(3,i) * x(i+1) ) / a(2,i)
          end do
          x(n) =   ( b(n) - a(1,n) * x(n-1)                   ) / a(2,n)

          x_norm = abs ( x(1) )
          do i = 2, n
            x_norm = max ( x_norm, abs ( x(i) ) )
          end do

          diff = abs ( x(1) - x_old(1) )
          do i = 2, n
            diff = max ( diff, abs ( x(i) - x_old(i) ) )
          end do

          if ( diff .le. tol * ( x_norm + 1.0D+00 ) ) then
            exit
          end if

        end do

      else

        do it_num = 1, it_max

          it = it_num

          x_old(1:n) = x(1:n)

          x(1) =   
     &      ( b(1)                     - a(1,2) * x(2)     ) / a(2,1)
          do i = 2, n - 1
            x(i) = 
     &      ( b(i) - a(3,i-1) * x(i-1) - a(1,i+1) * x(i+1) ) / a(2,i)
          end do
          x(n) =   
     &      ( b(n) - a(3,n-1) * x(n-1)                     ) / a(2,n)

          x_norm = abs ( x(1) )
          do i = 2, n
            x_norm = max ( x_norm, abs ( x(i) ) )
          end do

          diff = abs ( x(1) - x_old(1) )
          do i = 2, n
            diff = max ( diff, abs ( x(i) - x_old(i) ) )
          end do

          if ( diff .le. tol * ( x_norm + 1.0D+00 ) ) then
            exit
          end if
       
        end do

      end if

      return
      end
      subroutine r83_mxv ( n, a, x, b )

c*********************************************************************72
c
cc R83_MXV multiplies an R83 matrix times an R8VEC.
c
c  Discussion:
c
c    The R83 storage format is used for a tridiagonal matrix.
c    The superdiagonal is stored in entries (1,2:N), the diagonal in
c    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
c    original matrix is "collapsed" vertically into the array.
c
c  Example:
c
c    Here is how an R83 matrix of order 5 would be stored:
c
c       *  A12 A23 A34 A45
c      A11 A22 A33 A44 A55
c      A21 A32 A43 A54  *
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 May 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the linear system.
c
c    Input, double precision A(3,N), the R83 matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(N), the product A * x.
c
      implicit none

      integer n

      double precision a(3,n)
      double precision b(n)
      integer i
      double precision x(n)

      if ( n == 1 ) then

        b(1) = a(2,1) * x(1)

      else if ( 1 .lt. n ) then

        b(1) = a(2,1) * x(1) + a(1,2) * x(2)

        do i = 2, n - 1
          b(i) = a(3,i-1) * x(i-1) + a(2,i) * x(i) + a(1,i+1) * x(i+1)
        end do

        b(n) = a(3,n-1) * x(n-1) + a(2,n) * x(n)

      end if

      return
      end
      subroutine r83_print ( n, a, title )

c*********************************************************************72
c
cc R83_PRINT prints an R83 matrix.
c
c  Discussion:
c
c    The R83 storage format is used for a tridiagonal matrix.
c    The superdiagonal is stored in entries (1,2:N), the diagonal in
c    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
c    original matrix is "collapsed" vertically into the array.
c
c  Example:
c
c    Here is how an R83 matrix of order 5 would be stored:
c
c       *  A12 A23 A34 A45
c      A11 A22 A33 A44 A55
c      A21 A32 A43 A54  *
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 May 2009
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
c    Input, double precision A(3,N), the R83 matrix.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n

      double precision a(3,n)
      character * ( * ) title

      call r83_print_some ( n, a, 1, 1, n, n, title )

      return
      end
      subroutine r83_print_some ( n, a, ilo, jlo, ihi, jhi, title )

c*********************************************************************72
c
cc R83_PRINT_SOME prints some of an R83 matrix.
c
c  Discussion:
c
c    The R83 storage format is used for a tridiagonal matrix.
c    The superdiagonal is stored in entries (1,2:N), the diagonal in
c    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
c    original matrix is "collapsed" vertically into the array.
c
c  Example:
c
c    Here is how an R83 matrix of order 5 would be stored:
c
c       *  A12 A23 A34 A45
c      A11 A22 A33 A44 A55
c      A21 A32 A43 A54  *
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 May 2009
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
c    Input, double precision A(3,N), the R83 matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column, to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer, parameter :: incx = 5
      integer n

      double precision a(3,n)
      character * ( 14 ) ctemp(incx)
      logical r8_is_int
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
c
c  Print the columns of the matrix, in strips of 5.
c
      do j2lo = jlo, jhi, incx

        j2hi = j2lo + incx - 1
        j2hi = min ( j2hi, n )
        j2hi = min ( j2hi, jhi )

        inc = j2hi + 1 - j2lo

        write ( *, '(a)' ) ' '

        do j = j2lo, j2hi
          j2 = j + 1 - j2lo
          write ( ctemp(j2), '(i7,7x)' ) j
        end do

        write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) '  ---'
c
c  Determine the range of the rows in this strip.
c
        i2lo = max ( ilo, 1 )
        i2lo = max ( i2lo, j2lo - 1 )
        i2hi = min ( ihi, n )
        i2hi = min ( i2hi, j2hi + 1 )

        do i = i2lo, i2hi
c
c  Print out (up to) 5 entries in row I, that lie in the current strip.
c
          do j2 = 1, inc

            j = j2lo - 1 + j2

            if ( 1 .lt. i - j .or. 1 .lt. j - i ) then
              ctemp(j2) = '              '
            else if ( j .eq. i + 1 ) then
              write ( ctemp(j2), '(g14.6)' ) a(1,j)
            else if ( j .eq. i ) then
              write ( ctemp(j2), '(g14.6)' ) a(2,j)
            else if ( j .eq. i-1 ) then
              write ( ctemp(j2), '(g14.6)' ) a(3,j)
            end if

          end do

          write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j2), j2 = 1, inc )

        end do

      end do

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
      subroutine r8vec_print_some ( n, a, max_print, title )

c*********************************************************************72
c
cc R8VEC_PRINT_SOME prints "some" of an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    The user specifies MAX_PRINT, the maximum number of lines to print.
c
c    If N, the size of the vector, is no more than MAX_PRINT, then
c    the entire vector is printed, one entry per line.
c
c    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
c    followed by a line of periods suggesting an omission,
c    and the last entry.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries of the vector.
c
c    Input, double precision A(N), the vector to be printed.
c
c    Input, integer MAX_PRINT, the maximum number of lines to print.
c
c    Input, character*(*) TITLE, a title.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      integer max_print
      character*(*) title

      if ( max_print .le. 0 ) then
        return
      end if

      if ( n .le. 0 ) then
        return
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '

      if ( n .le. max_print ) then

        do i = 1, n
          write ( *, '(i6,a,1x,g14.6)' ) i, ':', a(i)
        end do

      else if ( 3 .le. max_print ) then

        do i = 1, max_print - 2
          write ( *, '(i6,a,1x,g14.6)' ) i, ':', a(i)
        end do

        write ( *, '(a)' ) '......  ..............'
        i = n

        write ( *, '(i6,a,1x,g14.6)' ) i, ':', a(i)

      else

        do i = 1, max_print - 1
          write ( *, '(i6,a,1x,g14.6)' ) i, ':', a(i)
        end do

        i = max_print

        write ( *, '(i6,a,1x,g14.6,a)' ) 
     &    i, ':', a(i), '...more entries...'

      end if

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