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
        write ( *, '(2x,i8,2x,2g14.6)' ) i, a(i)
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
        write ( *, '(2x,i8,2x,2g14.6)' ) i, x(i)
      end do

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
c    26 October 2008
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

        if ( i .ne. 5 .and. i .ne. 6 .and. i .ne. 9 ) then

          open ( unit = i, err = 10, status = 'scratch' )
          close ( unit = i )

          unit = i

          return
        end if

10      continue

      end do

      return
      end
      subroutine hilbert_inverse ( n, a )

c*****************************************************************************80
c
cc HILBERT_INVERSE returns the inverse of the Hilbert matrix.
c
c  Formula:
c
c    A(I,J) =  (-1)**(I+J) * (N+I-1)c * (N+J-1)c /
c           [ (I+J-1) * ((I-1)c*(J-1)c)**2 * (N-I)c * (N-J)c ]
c
c  Example:
c
c    N = 5
c
c       25    -300     1050    -1400     630
c     -300    4800   -18900    26880  -12600
c     1050  -18900    79380  -117600   56700
c    -1400   26880  -117600   179200  -88200
c      630  -12600    56700   -88200   44100
c
c  Properties:
c
c    A is symmetric.
c
c    Because A is symmetric, it is normal, so diagonalizable.
c
c    A is almost impossible to compute accurately by general routines
c    that compute the inverse.
c
c    A is integral.
c
c    The sum of the entries of A is N**2.
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
c    Input, integer N, the order of A.
c
c    Output, double precision A(N,N), the inverse Hilbert matrix.
c
      implicit none

      integer n

      double precision a(n,n)
      integer i
      integer j
c
c  Set the (1,1) entry.
c
      a(1,1) = dble ( n * n )
c
c  Define Row 1, Column J by recursion on Row 1 Column J-1
c
      i = 1
      do j = 2, n
        a(i,j) = - a(i,j-1) 
     &    * dble ( ( n + j - 1 ) * ( i + j - 2 ) * ( n + 1 - j )) 
     &    / dble ( ( i + j - 1 ) * ( j - 1 ) * ( j - 1 ) )
      end do
c
c  Define Row I by recursion on row I-1
c
      do i = 2, n
        do j = 1, n

          a(i,j) = - a(i-1,j) 
     &      * dble ( (n+i-1) * (i+j-2) * (n+1-i) ) 
     &      / dble ( (i+j-1) * ( i - 1 ) * ( i - 1 ) )

        end do
      end do

      return
      end
      function i4_log_10 ( i )

c*********************************************************************72
c
cc I4_LOG_10 returns the integer part of the logarithm base 10 of ABS(X).
c
c  Discussion:
c
c    I4_LOG_10 ( I ) + 1 is the number of decimal digits in I.
c
c  Example:
c
c        I  I4_LOG_10
c    -----  --------
c        0    0
c        1    0
c        2    0
c        9    0
c       10    1
c       11    1
c       99    1
c      100    2
c      101    2
c      999    2
c     1000    3
c     1001    3
c     9999    3
c    10000    4
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the number whose logarithm base 10 is desired.
c
c    Output, integer I4_LOG_10, the integer part of the logarithm base 10 of
c    the absolute value of X.
c
      implicit none

      integer i
      integer i_abs
      integer i4_log_10
      integer ten_pow

      if ( i .eq. 0 ) then

        i4_log_10 = 0

      else

        i4_log_10 = 0
        ten_pow = 10

        i_abs = abs ( i )

10      continue

        if ( ten_pow .le. i_abs ) then
          i4_log_10 = i4_log_10 + 1
          ten_pow = ten_pow * 10
          go to 10
        end if

      end if

      return
      end
      subroutine i4_swap ( i, j )

c*********************************************************************72
c
cc I4_SWAP switches two I4's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 January 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer I, J.  On output, the values of I and
c    J have been interchanged.
c
      implicit none

      integer i
      integer j
      integer k

      k = i
      i = j
      j = k

      return
      end
      function i4_uniform ( a, b, seed )

c*********************************************************************72
c
cc I4_UNIFORM returns a scaled pseudorandom I4.
c
c  Discussion:
c
c    An I4 is an integer value.
c
c    The pseudorandom number should be uniformly distributed
c    between A and B.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 November 2006
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
c    Peter Lewis, Allen Goodman, James Miller
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, pages 136-143, 1969.
c
c  Parameters:
c
c    Input, integer A, B, the limits of the interval.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, integer I4_UNIFORM, a number between A and B.
c
      implicit none

      integer a
      integer b
      integer i4_uniform
      integer k
      real r
      integer seed
      integer value

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_UNIFORM - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + 2147483647
      end if

      r = real ( seed ) * 4.656612875E-10
c
c  Scale R to lie between A-0.5 and B+0.5.
c
      r = ( 1.0E+00 - r ) * ( real ( min ( a, b ) ) - 0.5E+00 )
     &  +             r   * ( real ( max ( a, b ) ) + 0.5E+00 )
c
c  Use rounding to convert R to an integer between A and B.
c
      value = nint ( r )

      value = max ( value, min ( a, b ) )
      value = min ( value, max ( a, b ) )

      i4_uniform = value

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
      integer s_len_trim
      character*(*) title
      integer title_length

      title_length = s_len_trim ( title )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) title(1:title_length)

      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,2x,i12)' ) i, a(i)
      end do

      return
      end
      subroutine i4vec_search_binary_a ( n, a, b, indx )

c*********************************************************************72
c
cc I4VEC_SEARCH_BINARY_A searches an ascending sorted I4VEC for a value.
c
c  Discussion:
c
c    An I4VEC is a vector of I4 values.
c
c    Binary search is used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    31 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Donald Kreher, Douglas Simpson,
c    Algorithm 1.9,
c    Combinatorial Algorithms,
c    CRC Press, 1998, page 26.
c
c  Parameters:
c
c    Input, integer N, the number of elements in the vector.
c
c    Input, integer A(N), the array to be searched.  A must
c    be sorted in ascending order.
c
c    Input, integer B, the value to be searched for.
c
c    Output, integer INDX, the result of the search.
c    -1, B does not occur in A.
c    I, A(I) = B.
c
      implicit none

      integer n

      integer a(n)
      integer b
      integer high
      integer indx
      integer low
      integer mid

      indx = - 1

      low = 1
      high = n

10    continue

      if ( low .le. high ) then

        mid = ( low + high ) / 2

        if ( a(mid) .eq. b ) then
          indx = mid
          go to 20
        else if ( a(mid) .lt. b ) then
          low = mid + 1
        else if ( b .lt. a(mid) ) then
          high = mid - 1
        end if

        go to 10

      end if

20    continue

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
      function r8_uniform ( a, b, seed )

c*********************************************************************72
c
cc R8_UNIFORM returns a scaled pseudorandom R8.
c
c  Discussion:
c
c    The pseudorandom number should be uniformly distributed
c    between A and B.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 January 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A, B, the limits of the interval.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R8_UNIFORM, a number strictly between A and B.
c
      implicit none

      double precision a
      double precision b
      integer k
      double precision r8_uniform
      integer seed

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_UNIFORM - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + 2147483647
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r8_uniform = a + ( b - a ) * dble ( seed ) * 4.656612875D-10

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
      rhs(1:n) = b(1:n)
      rhs(n+1:2*n) = 0.0D+00

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

      x(1:n) = rhs(1:n)

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
      subroutine r83_indicator ( n, a )

c*********************************************************************72
c
cc R83_INDICATOR sets up an R83 indicator matrix.
c
c  Discussion:
c
c    The "indicator matrix" simply has a value like I*10+J at every
c    entry of a dense matrix, or at every entry of a compressed storage
c    matrix for which storage is allocated. 
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
c    Here are the values as stored in an indicator matrix:
c
c      00 12 23 34 45
c      11 22 33 44 55
c      21 32 43 54 00
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
c    Output, double precision A(3,N), the R83 indicator matrix.
c
      implicit none

      integer n

      double precision a(3,n)
      integer fac
      integer i
      integer i4_log_10
      integer j

      fac = 10 ** ( i4_log_10 ( n ) + 1 )

      a(1,1) = 0.0D+00
      do j = 2, n
        i = j - 1
        a(1,j) = dble ( fac * i + j )
      end do

      do j = 1, n
        i = j
        a(2,j) = dble ( fac * i + j )
      end do

      do j = 1, n-1
        i = j + 1
        a(3,j) = dble ( fac * i + j )
      end do
      a(3,n) = 0.0D+00

      return
      end
      subroutine r83_jac_sl ( n, a, b, x, tol, it_max, job, it, diff )

c*********************************************************************72
c
cc R83_JAC_SL solves an R83 system using Jacobi iteration.
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
c    Input/output, double precision X(N), an approximate solution 
c    to the system.
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
      double precision x_new(n)
      double precision x_norm
c
c  No diagonal matrix entry can be zero.
c
      do i = 1, n
        if ( a(2,i) .eq. 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R83_JAC_SL - Fatal error!'
          write ( *, '(a,i8)' ) '  Zero diagonal entry, index = ', i
          return
        end if
      end do

      if ( job .eq. 0 ) then

        do it_num = 1, it_max

          it = it_num

          x_new(1) =   b(1)                   - a(3,1) * x(2)
          do i = 2, n - 1
            x_new(i) = b(i) - a(1,i) * x(i-1) - a(3,i) * x(i+1)
          end do
          x_new(n) =   b(n) - a(1,n) * x(n-1)
c
c  Divide by diagonal terms.
c
          x_new(1:n) = x_new(1:n) / a(2,1:n)
c
c  Measure norms of solution and change:
c
          x_norm = abs ( x(1) )
          do i = 2, n
            x_norm = max ( x_norm, abs ( x(i) ) )
          end do

          diff = abs ( x_new(1) - x(1) )
          do i = 2, n
            diff = max ( diff, abs ( x_new(i) - x(i) ) )
          end do
c
c  Update.
c
          x(1:n) = x_new(1:n)
c
c  Test for early termination.
c
          if ( diff .le. tol * ( x_norm + 1.0D+00 ) ) then
            exit
          end if

        end do

      else

        do it_num = 1, it_max

          it = it_num

          x_new(1) =   b(1)                     - a(1,2) * x(2)
          do i = 2, n - 1
            x_new(i) = b(i) - a(3,i-1) * x(i-1) - a(1,i+1) * x(i+1)
          end do
          x_new(n) =   b(n) - a(3,n-1) * x(n-1)
c
c  Divide by diagonal terms.
c
          x_new(1:n) = x_new(1:n) / a(2,1:n)
c
c  Measure norms of solution and change:
c
          x_norm = abs ( x(1) )
          do i = 2, n
            x_norm = max ( x_norm, abs ( x(i) ) )
          end do

          diff = abs ( x_new(1) - x(1) )
          do i = 2, n
            diff = max ( diff, abs ( x_new(i) - x(i) ) )
          end do
c
c  Update:
c
          x(1:n) = x_new(1:n)
c
c  Test for early termination.
c
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
cc R83_MXV multiplies an R83 matrix times a vector.
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
      double precision x(n)

      b(1:n)   =            a(2,1:n)   * x(1:n)
      b(1:n-1) = b(1:n-1) + a(1,2:n)   * x(2:n)
      b(2:n)   = b(2:n)   + a(3,1:n-1) * x(1:n-1)

      return
      end
      subroutine r83_np_det ( n, a_lu, det )

c*********************************************************************72
c
cc R83_NP_DET returns the determinant of an R83 system factored by R83_NP_FA.
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
c    N must be at least 2.
c
c    Input, double precision A_LU(3,N), the LU factors computed by R83_NP_FA.
c
c    Output, double precision DET, the determinant of the matrix.
c
      implicit none

      integer n

      double precision a_lu(3,n)
      double precision det
      integer j

      det = 1.0D+00
      do j = 1, n
        det = det * a_lu(2,j)
      end do

      return
      end
      subroutine r83_np_fa ( n, a, info )

c*********************************************************************72
c
cc R83_NP_FA factors an R83 matrix without pivoting.
c
c  Discussion:
c
c    The R83 storage format is used for a tridiagonal matrix.
c    The superdiagonal is stored in entries (1,2:N), the diagonal in
c    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
c    original matrix is "collapsed" vertically into the array.
c
c    Because this routine does not use pivoting, it can fail even when
c    the matrix is not singular, and it is liable to make larger
c    errors.
c
c    R83_NP_FA and R83_NP_SL may be preferable to the corresponding
c    LINPACK routine SGTSL for tridiagonal systems, which factors and solves
c    in one step, and does not save the factorization.
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
c    Input/output, double precision A(3,N).
c    On input, the tridiagonal matrix.  On output, factorization information.
c
c    Output, integer INFO, singularity flag.
c    0, no singularity detected.
c    nonzero, the factorization failed on the INFO-th step.
c
      implicit none

      integer n

      double precision a(3,n)
      integer i
      integer info

      info = 0

      do i = 1, n-1

        if ( a(2,i) .eq. 0.0D+00 ) then
          info = i
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R83_NP_FA - Fatal error!'
          write ( *, '(a,i8)' ) '  Zero pivot on step ', info
          return
        end if
c
c  Store the multiplier in L.
c
        a(3,i) = a(3,i) / a(2,i)
c
c  Modify the diagonal entry in the next column.
c
        a(2,i+1) = a(2,i+1) - a(3,i) * a(1,i+1)

      end do

      if ( a(2,n) .eq. 0.0D+00 ) then
        info = n
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R83_NP_FA - Fatal error!'
        write ( *, '(a,i8)' ) '  Zero pivot on step ', info
        return
      end if

      return
      end
      subroutine r83_np_fs ( n, a, b, x )

c*********************************************************************72
c
cc R83_NP_FS factors and solves an R83 system.
c
c  Discussion:
c
c    The R83 storage format is used for a tridiagonal matrix.
c    The superdiagonal is stored in entries (1,2:N), the diagonal in
c    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
c    original matrix is "collapsed" vertically into the array.
c
c    This algorithm requires that each diagonal entry be nonzero.
c    It does not use pivoting, and so can fail on systems that
c    are actually nonsingular.
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
c    Input, integer N, the order of the linear system.
c
c    Input/output, double precision A(3,N).
c    On input, the tridiagonal matrix.
c    On output, the data in these vectors has been overwritten
c    by factorization information.
c
c    Input, double precision B(N), the right hand side of the linear system.
c
c    Output, double precision X(N), the solution of the linear system.
c
      implicit none

      integer n

      double precision a(3,n)
      double precision b(n)
      integer i
      double precision x(n)
      double precision xmult
c
c  The diagonal entries can't be zero.
c
      do i = 1, n
        if ( a(2,i) .eq. 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R83_NP_FS - Fatal error!'
          write ( *, '(a,i8,a)' ) '  A(2,', i, ') = 0.'
          return
        end if
      end do

      do i = 1, n
        x(i) = b(i)
      end do


      do i = 2, n
        xmult = a(3,i-1) / a(2,i-1)
        a(2,i) = a(2,i) - xmult * a(1,i)
        x(i)   = x(i)   - xmult * x(i-1)
      end do

      x(n) = x(n) / a(2,n)
      do i = n-1, 1, -1
        x(i) = ( x(i) - a(1,i+1) * x(i+1) ) / a(2,i)
      end do

      return
      end
      subroutine r83_np_fss ( n, a, nb, b, x )

c*********************************************************************72
c
cc R83_NP_FSS factors and solves multiple R83 systems.
c
c  Discussion:
c
c    The R83 storage format is used for a tridiagonal matrix.
c    The superdiagonal is stored in entries (1,2:N), the diagonal in
c    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
c    original matrix is "collapsed" vertically into the array.
c
c    This algorithm requires that each diagonal entry be nonzero.
c    It does not use pivoting, and so can fail on systems that
c    are actually nonsingular.
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
c    29 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the linear system.
c
c    Input/output, double precision A(3,N).
c    On input, the tridiagonal matrix.
c    On output, the data in these vectors has been overwritten
c    by factorization information.
c
!    Input, integer NB, the number of right hand sides.
!
c    Input, double precision B(N,NB), the right hand sides of the linear system.
c
c    Output, double precision X(N,NB), the solutions of the linear system.
c
      implicit none

      integer n
      integer nb

      double precision a(3,n)
      double precision b(n,nb)
      integer i
      integer j
      double precision x(n,nb)
      double precision xmult
c
c  The diagonal entries can't be zero.
c
      do i = 1, n
        if ( a(2,i) .eq. 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R83_NP_FSS - Fatal error!'
          write ( *, '(a,i8,a)' ) '  A(2,', i, ') = 0.'
          return
        end if
      end do

      do j = 1, nb
        do i = 1, n
          x(i,j) = b(i,j)
        end do
      end do

      do i = 2, n
        xmult = a(3,i-1) / a(2,i-1)
        a(2,i) = a(2,i) - xmult * a(1,i)
        do j = 1, nb
          x(i,j)   = x(i,j) - xmult * x(i-1,j)
        end do
      end do

      do j = 1, nb
        x(n,j) = x(n,j) / a(2,n)
      end do

      do i = n-1, 1, -1
        do j = 1, nb
          x(i,j) = ( x(i,j) - a(1,i+1) * x(i+1,j) ) / a(2,i)
        end do
      end do

      return
      end
      subroutine r83_np_ml ( n, a_lu, x, b, job )

c*********************************************************************72
c
cc R83_NP_ML computes A * x or x * A, where A has been factored by R83_NP_FA.
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
c    N must be at least 2.
c
c    Input, double precision A_LU(3,N), the LU factors from R83_FA.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(N), the product A*x or A'*x.
c
c    Input, integer JOB, specifies the product to find.
c    0, compute A * x.
c    nonzero, compute A' * x.
c
      implicit none

      integer n

      double precision a_lu(3,n)
      double precision b(n)
      integer i
      integer job
      double precision x(n)

      do i = 1, n
        b(i) = x(i)
      end do

      if ( job .eq. 0 ) then
c
c  Compute X := U * X
c
        do i = 1, n

          b(i) = a_lu(2,i) * b(i)

          if ( i .lt. n ) then
            b(i) = b(i) + a_lu(1,i+1) * b(i+1)
          end if

        end do
c
c  Compute X: = L * X.
c
        do i = n, 2, -1
          b(i) = b(i) + a_lu(3,i-1) * b(i-1)
        end do

      else
c
c  Compute X: = L' * X.
c
        do i = 1, n-1
          b(i) = b(i) + a_lu(3,i) * b(i+1)
        end do
c
c  Compute X: = U' * X.
c
        do i = n, 2, -1
          b(i) = a_lu(2,i) * b(i)
          b(i) = b(i) + a_lu(1,i) * b(i-1)
        end do
        b(1) = a_lu(2,1) * b(1)

      end if

      return
      end
      subroutine r83_np_sl ( n, a_lu, b, job )

c*********************************************************************72
c
cc R83_NP_SL solves an R83 system factored by R83_NP_FA.
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
c    N must be at least 2.
c
c    Input, double precision A_LU(3,N), the LU factors from R83_NP_FA.
c
c    Input/output, double precision B(N).
c    On input, B contains the right hand side of the linear system.
c    On output, B contains the solution of the linear system.
c
c    Input, integer JOB, specifies the system to solve.
c    0, solve A * x = b.
c    nonzero, solve A' * x = b.
c
      implicit none

      integer n

      double precision a_lu(3,n)
      double precision b(n)
      integer i
      integer job

      if ( job .eq. 0 ) then
c
c  Solve L * Y = B.
c
        do i = 2, n
          b(i) = b(i) - a_lu(3,i-1) * b(i-1)
        end do
c
c  Solve U * X = Y.
c
        do i = n, 1, -1
          b(i) = b(i) / a_lu(2,i)
          if ( 1 .lt. i ) then
            b(i-1) = b(i-1) - a_lu(1,i) * b(i)
          end if
        end do

      else
c
c  Solve U' * Y = B
c
        do i = 1, n
          b(i) = b(i) / a_lu(2,i)
          if ( i .lt. n ) then
            b(i+1) = b(i+1) - a_lu(1,i+1) * b(i)
          end if
        end do
c
c  Solve L' * X = Y.
c
        do i = n-1, 1, -1
          b(i) = b(i) - a_lu(3,i) * b(i+1)
        end do

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

            if ( 1 .lt. i-j .or. 1 .lt. j-i ) then
              ctemp(j2) = '              '
            else if ( j .eq. i+1 ) then
              if ( r8_is_int ( a(1,j) ) ) then
                write ( ctemp(j2), '(f8.0,6x)' ) a(1,j)
              else
                write ( ctemp(j2), '(g14.6)' ) a(1,j)
              end if
            else if ( j .eq. i ) then
              if ( r8_is_int ( a(2,j) ) ) then
                write ( ctemp(j2), '(f8.0,6x)' ) a(2,j)
              else
                write ( ctemp(j2), '(g14.6)' ) a(2,j)
              end if
            else if ( j .eq. i-1 ) then
              if ( r8_is_int ( a(3,j) ) ) then
                write ( ctemp(j2), '(f8.0,6x)' ) a(3,j)
              else
                write ( ctemp(j2), '(g14.6)' ) a(3,j)
              end if
            end if

          end do

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine r83_random ( n, seed, a )

c*********************************************************************72
c
cc R83_RANDOM randomizes an R83 matrix.
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
c    Input, integer N, the order of the linear system.
c
c    Input/output, integer SEED, a seed for the random 
c    number generator.
c
c    Output, double precision A(3,N), the R83 matrix.
c
      implicit none

      integer n

      double precision a(3,n)
      integer, parameter :: i4_1 = 1
      integer seed

      a(1,1) = 0.0D+00
      call r8vec_uniform_01 ( n - i4_1, seed, a(1,2:n) )

      call r8vec_uniform_01 ( n,        seed, a(2,1:n) )

      call r8vec_uniform_01 ( n - i4_1, seed, a(3,1:n-1) )
      a(3,n) = 0.0D+00

      return
      end
      subroutine r83_to_r8ge ( n, a, b )

c*********************************************************************72
c
cc R83_TO_R8GE copies an R83 matrix to an R8GE matrix.
c
c  Discussion:
c
c    The R83 storage format is used for a tridiagonal matrix.
c    The superdiagonal is stored in entries (1,2:N), the diagonal in
c    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
c    original matrix is "collapsed" vertically into the array.
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
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
c    Output, double precision B(N,N), the R8GE matrix.
c
      implicit none

      integer n

      double precision a(3,n)
      double precision b(n,n)
      integer i
      integer j

      do i = 1, n
        do j = 1, n

          if ( j .eq. i-1 ) then
            b(i,j) = a(3,j)
          else if ( j .eq. i ) then
            b(i,j) = a(2,j)
          else if ( j .eq. i+1 ) then
            b(i,j) = a(1,j)
          else
            b(i,j) = 0.0D+00
          end if

        end do
      end do

      return
      end
      subroutine r83_vxm ( n, a, x, b )

c*********************************************************************72
c
cc R83_VXM multiplies a vector by an R83 matrix.
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
c    Input, integer N, the order of the linear system.
c
c    Input, double precision A(3,N), the R83 matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A'.
c
c    Output, double precision B(N), the product A' * x.
c
      implicit none

      integer n

      double precision a(3,n)
      double precision b(n)
      double precision x(n)

      b(1:n)   =            a(2,1:n)   * x(1:n)
      b(1:n-1) = b(1:n-1) + a(3,1:n-1) * x(2:n) 
      b(2:n)   = b(2:n)   + a(1,2:n)   * x(1:n-1)

      return
      end
      subroutine r83np_fs ( n, a, b, x )

c*********************************************************************72
c
cc R83NP_FS factors and solves an R83NP system.
c
c  Discussion:
c
c    The R83NP storage format is used for a tridiagonal matrix.
c    The subdiagonal   is in entries (1,2:N), 
c    the diagonal      is in entries (2,1:N), 
c    the superdiagonal is in entries (3,1:N-1). 
c
c    This algorithm requires that each diagonal entry be nonzero.
c    It does not use pivoting, and so can fail on systems that
c    are actually nonsingular.
c
c    The "R83NP" format used for this routine is different from the R83 format.
c    Here, we insist that the nonzero entries
c    for a given row now appear in the corresponding column of the
c    packed array.
c
c  Example:
c
c    Here is how an R83NP matrix of order 5 would be stored:
c
c       *  A21 A32 A43 A54
c      A11 A22 A33 A44 A55
c      A12 A23 A34 A45  *
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
c    Input, integer N, the order of the linear system.
c
c    Input/output, double precision A(3,N).
c    On input, the tridiagonal matrix.
c    On output, the data in these vectors has been overwritten
c    by factorization information.
c
c    Input, double precision B(N), the right hand side of the linear system.
c
c    Output, double precision X(N), the solution of the linear system.
c
      implicit none

      integer n

      double precision a(3,n)
      double precision b(n)
      integer i
      double precision x(n)
      double precision xmult
c
c  The diagonal entries can't be zero.
c
      do i = 1, n
        if ( a(2,i) .eq. 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R83NP_FS - Fatal error!'
          write ( *, '(a,i8,a)' ) '  A(2,', i, ') = 0.'
          return
        end if
      end do

      do i = 1, n
        x(i) = b(i)
      end do

      do i = 2, n
        a(2,i) = a(2,i) - a(3,i-1) * a(1,i) / a(2,i-1)
        x(i)   = x(i)   - x(i-1)   * a(1,i) / a(2,i-1)
      end do

      x(n) = x(n) / a(2,n)
      do i = n-1, 1, -1
        x(i) = ( x(i) - a(3,i) * x(i+1) ) / a(2,i)
      end do

      return
      end
      subroutine r83p_det ( n, a_lu, work4, det )

c*********************************************************************72
c
cc R83P_DET computes the determinant of a matrix factored by R83P_FA.
c
c  Discussion:
c
c    The R83P storage format stores a periodic tridiagonal matrix as 
c    a 3 by N array, in which each row corresponds to a diagonal, and 
c    column locations are preserved.  The matrix value 
c    A(1,N) is stored as the array entry A(3,N), and the matrix value
c    A(N,1) is stored as the array entry A(1,1).
c
c  Example:
c
c    Here is how an R83P matrix of order 5 would be stored:
c
c      A51 A12 A23 A34 A45
c      A11 A22 A33 A44 A55
c      A21 A32 A43 A54 A15
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be at least 3.
c
c    Input, double precision A_LU(3,N), LU factors from R83P_FA.
c
c    Input, double precision WORK4, factorization information from R83P_FA.
c
c    Output, double precision DET, the determinant of the matrix.
c
      implicit none

      integer n

      double precision a_lu(3,n)
      double precision det
      integer j
      double precision work4

      det = work4

      do j = 1, n - 1
        det = det * a_lu(2,j)
      end do

      return
      end
      subroutine r83p_fa ( n, a, info, work2, work3, work4 )

c*********************************************************************72
c
cc R83P_FA factors an R83P matrix.
c
c  Discussion:
c
c    The R83P storage format stores a periodic tridiagonal matrix as 
c    a 3 by N array, in which each row corresponds to a diagonal, and 
c    column locations are preserved.  The matrix value 
c    A(1,N) is stored as the array entry A(3,N), and the matrix value
c    A(N,1) is stored as the array entry A(1,1).
c
c    Once the matrix has been factored by R83P_FA, R83P_SL may be called
c    to solve linear systems involving the matrix.
c
c    The logical matrix has a form which is suggested by this diagram:
c
c      D1 U1          L1
c      L2 D2 U2
c         L3 D3 U3
c            L4 D4 U4
c               L5 D5 U5
c      U6          L6 D6
c
c    The algorithm treats the matrix as a border banded matrix:
c
c      ( A1  A2 )
c      ( A3  A4 )
c
c    where:
c
c      D1 U1          | L1
c      L2 D2 U2       |  0
c         L3 D3  U3    |  0
c            L4 D4 U4 |  0
c               L5 D5 | U5
c      ---------------+---
c      U6  0  0  0 L6 | D6
c
c  Example:
c
c    Here is how an R83P matrix of order 5 would be stored:
c
c      A51 A12 A23 A34 A45
c      A11 A22 A33 A44 A55
c      A21 A32 A43 A54 A15
c
c  Method:
c
c    The algorithm rewrites the system as:
c
c         X1 + inverse(A1) A2 X2 = inverse(A1) B1
c
c      A3 X1 +             A4 X2 = B2
c
c    The first equation can be "solved" for X1 in terms of X2:
c
c         X1 = - inverse(A1) A2 X2 + inverse(A1) B1
c
c    allowing us to rewrite the second equation for X2 explicitly:
c
c      ( A4 - A3 inverse(A1) A2 ) X2 = B2 - A3 inverse(A1) B1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be at least 3.
c
c    Input/output, double precision A(3,N).
c    On input, the periodic tridiagonal matrix.  
c    On output, the arrays have been modified to hold information
c    defining the border-banded factorization of submatrices A1
c    and A3.
c
c    Output, integer INFO, singularity flag.
c    0, no singularity detected.
c    nonzero, the factorization failed on the INFO-th step.
c
c    Output, double precision WORK2(N-1), WORK3(N-1), WORK4, 
c    factorization information.
c
      implicit none

      integer n

      double precision a(3,n)
      integer i
      integer info
      integer job
      double precision work2(n-1)
      double precision work3(n-1)
      double precision work4
c
c  Compute inverse(A1):
c
      call r83_np_fa ( n - 1, a, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R83P_FA - Fatal error!'
        write ( *, '(a,i8)' ) '  R83_NP_FA returned INFO = ', info
        write ( *, '(a)' ) '  Factoring failed for column INFO.'
        write ( *, '(a)' ) '  The tridiagonal matrix A1 is singular.'
        write ( *, '(a)' ) '  This algorithm cannot continue.'
        return
      end if
c
c  WORK2 := inverse(A1) * A2.
c
      work2(1) = a(3,n)
      do i = 2, n - 2
        work2(i) = 0.0D+00
      end do
      work2(n-1) = a(1,n)

      job = 0
      call r83_np_sl ( n - 1, a, work2, job )
c
c  WORK3 := inverse ( A1' ) * A3'.
c
      work3(1) = a(1,1)
      do i = 2, n - 2
        work3(i) = 0.0D+00
      end do
      work3(n-1) = a(3,n-1)

      job = 1
      call r83_np_sl ( n - 1, a, work3, job )
c
c  A4 := ( A4 - A3 * inverse(A1) * A2 )
c
      work4 = a(2,n) - a(1,1) * work2(1) - a(3,n-1) * work2(n-1)

      if ( work4 .eq. 0.0D+00 ) then
        info = n
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R83P_FA - Fatal error!'
        write ( *, '(a)' ) '  The factored A4 submatrix is zero.'
        write ( *, '(a)' ) '  This algorithm cannot continue.'
        return
      end if

      return
      end
      subroutine r83p_indicator ( n, a )

c*********************************************************************72
c
cc R83P_INDICATOR sets up an R83P indicator matrix.
c
c  Discussion:
c
c    The "indicator matrix" simply has a value like I*10+J at every
c    entry of a dense matrix, or at every entry of a compressed storage
c    matrix for which storage is allocated. 
c
c    The R83P storage format stores a periodic tridiagonal matrix as 
c    a 3 by N array, in which each row corresponds to a diagonal, and 
c    column locations are preserved.  The matrix value 
c    A(1,N) is stored as the array entry A(3,N), and the matrix value
c    A(N,1) is stored as the array entry A(1,1).
c
c  Example:
c
c    Here is how an R83P matrix of order 5 would be stored:
c
c      A51 A12 A23 A34 A45
c      A11 A22 A33 A44 A55
c      A21 A32 A43 A54 A15
c
c    Here are the values as stored in an indicator matrix:
c
c      51 12 23 34 45
c      11 22 33 44 55
c      21 32 43 54 15
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 May 2009
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
c    Output, double precision A(3,N), the R83P indicator matrix.
c
      implicit none

      integer n

      double precision a(3,n)
      integer fac
      integer i
      integer i4_log_10
      integer j

      fac = 10 ** ( i4_log_10 ( n ) + 1 )

      i = n
      j = 1
      a(1,j) = dble ( fac * i + j )
      do j = 2, n
        i = j - 1
        a(1,j) = dble ( fac * i + j )
      end do

      do j = 1, n
        i = j
        a(2,j) = dble ( fac * i + j )
      end do

      do j = 1, n - 1
        i = j + 1
        a(3,j) = dble ( fac * i + j )
      end do
      i = 1
      j = n
      a(3,j) = dble ( fac * i + j )

      return
      end
      subroutine r83p_ml ( n, a_lu, x, b, job )

c*********************************************************************72
c
cc R83P_ML computes A * x or x * A, where A has been factored by R83P_FA.
c
c  Discussion:
c
c    The R83P storage format stores a periodic tridiagonal matrix as 
c    a 3 by N array, in which each row corresponds to a diagonal, and 
c    column locations are preserved.  The matrix value 
c    A(1,N) is stored as the array entry A(3,N), and the matrix value
c    A(N,1) is stored as the array entry A(1,1).
c
c  Example:
c
c    Here is how an R83P matrix of order 5 would be stored:
c
c      A51 A12 A23 A34 A45
c      A11 A22 A33 A44 A55
c      A21 A32 A43 A54 A15
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be at least 3.
c
c    Input, double precision A_LU(3,N), the LU factors ffrom R83P_FA.
c
c    Input, double precision X(N), the vector to be multiplied by the matrix.
c
c    Output, double precision B(N), the result of the multiplication.
c
c    Input, integer JOB, indicates what product should be computed.
c    0, compute A * x.
c    nonzero, compute A' * x.
c
      implicit none

      integer n

      double precision a_lu(3,n)
      double precision b(n)
      integer job
      double precision x(n)
c
c  Multiply A(1:N-1,1:N-1) and X(1:N-1).
c
      call r83_np_ml ( n - 1, a_lu, x, b, job )
c
c  Add terms from the border.
c
      if ( job .eq. 0 ) then
        b(1) = b(1) + a_lu(3,n) * x(n)
        b(n-1) = b(n-1) + a_lu(1,n) * x(n)
        b(n) = a_lu(1,1) * x(1) + a_lu(3,n-1) * x(n-1) 
     &    + a_lu(2,n) * x(n)
      else
        b(1) = b(1) + a_lu(1,1) * x(n)
        b(n-1) = b(n-1) + a_lu(3,n-1) * x(n)
        b(n) = a_lu(3,n) * x(1) + a_lu(1,n) * x(n-1) 
     &    + a_lu(2,n) * x(n)
      end if

      return
      end
      subroutine r83p_mxv ( n, a, x, b )

c*********************************************************************72
c
cc R83P_MXV multiplies an R83P matrix times a vector.
c
c  Discussion:
c
c    The R83P storage format stores a periodic tridiagonal matrix as 
c    a 3 by N array, in which each row corresponds to a diagonal, and 
c    column locations are preserved.  The matrix value 
c    A(1,N) is stored as the array entry A(3,N), and the matrix value
c    A(N,1) is stored as the array entry A(1,1).
c
c  Example:
c
c    Here is how an R83P matrix of order 5 would be stored:
c
c      A51 A12 A23 A34 A45
c      A11 A22 A33 A44 A55
c      A21 A32 A43 A54 A15
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be at least 3.
c
c    Input, double precision A(3,N), the R83P matrix.
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

      b(1) =   a(3,n)   * x(n)   + a(2,1) * x(1) + a(1,2)   * x(2)

      do i = 2, n-1
        b(i) = a(3,i-1) * x(i-1) + a(2,i) * x(i) + a(1,i+1) * x(i+1)
      end do

      b(n) =   a(3,n-1) * x(n-1) + a(2,n) * x(n) + a(1,1)   * x(1)

      return
      end
      subroutine r83p_print ( n, a, title )

c*********************************************************************72
c
cc R83P_PRINT prints an R83P matrix.
c
c  Discussion:
c
c    The R83P storage format stores a periodic tridiagonal matrix as 
c    a 3 by N array, in which each row corresponds to a diagonal, and 
c    column locations are preserved.  The matrix value 
c    A(1,N) is stored as the array entry A(3,N), and the matrix value
c    A(N,1) is stored as the array entry A(1,1).
c
c  Example:
c
c    Here is how an R83P matrix of order 5 would be stored:
c
c      A51 A12 A23 A34 A45
c      A11 A22 A33 A44 A55
c      A21 A32 A43 A54 A15
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 May 2009
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
c    Input, double precision A(3,N), the R83P matrix.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n

      double precision a(3,n)
      character * ( * ) title

      call r83p_print_some ( n, a, 1, 1, n, n, title )

      return
      end
      subroutine r83p_print_some ( n, a, ilo, jlo, ihi, jhi, title )

c*********************************************************************72
c
cc R83P_PRINT_SOME prints some of an R83P matrix.
c
c  Discussion:
c
c    The R83P storage format stores a periodic tridiagonal matrix as 
c    a 3 by N array, in which each row corresponds to a diagonal, and 
c    column locations are preserved.  The matrix value 
c    A(1,N) is stored as the array entry A(3,N), and the matrix value
c    A(N,1) is stored as the array entry A(1,1).
c
c  Example:
c
c    Here is how an R83P matrix of order 5 would be stored:
c
c      A51 A12 A23 A34 A45
c      A11 A22 A33 A44 A55
c      A21 A32 A43 A54 A15
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 May 2009
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
c    Input, double precision A(3,N), the R83P matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column, to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 5 )
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
c  Determine the column range of the rows in this strip.
c
        i2lo = max ( ilo, 1 )

        if ( 1 .lt. i2lo .or. j2hi .lt. n ) then
          i2lo = max ( i2lo, j2lo - 1 )
        end if

        i2hi = min ( ihi, n )

        if ( i2hi .lt. n .or. 1 .lt. j2lo ) then
          i2hi = min ( i2hi, j2hi + 1 )
        end if

        do i = i2lo, i2hi
c
c  Print out (up to) 5 entries in row I, that lie in the current strip.
c
          do j2 = 1, inc

            j = j2lo - 1 + j2

            if ( i .eq. n .and. j .eq. 1 ) then
              if ( r8_is_int ( a(1,j) ) ) then
                write ( ctemp(j2), '(f8.0,6x)' ) a(1,j)
              else
                write ( ctemp(j2), '(g14.6)' ) a(1,j)
              end if
            else if ( i .eq. 1 .and. j .eq. n ) then
              if ( r8_is_int ( a(3,j) ) ) then
                write ( ctemp(j2), '(f8.0,6x)' ) a(3,j)
              else
                write ( ctemp(j2), '(g14.6)' ) a(3,j)
              end if
            else if ( 1 .lt. i-j .or. 1 .lt. j-i ) then
              ctemp(j2) = '              '
            else if ( j .eq. i+1 ) then
              if ( r8_is_int ( a(1,j) ) ) then
                write ( ctemp(j2), '(f8.0,6x)' ) a(1,j)
              else
                write ( ctemp(j2), '(g14.6)' ) a(1,j)
              end if
            else if ( j .eq. i ) then
              if ( r8_is_int ( a(2,j) ) ) then
                write ( ctemp(j2), '(f8.0,6x)' ) a(2,j)
              else
                write ( ctemp(j2), '(g14.6)' ) a(2,j)
              end if
            else if ( j .eq. i-1 ) then
              if ( r8_is_int ( a(3,j) ) ) then
                write ( ctemp(j2), '(f8.0,6x)' ) a(3,j)
              else
                write ( ctemp(j2), '(g14.6)' ) a(3,j)
              end if
            end if

          end do

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine r83p_random ( n, seed, a )

c*********************************************************************72
c
cc R83P_RANDOM randomizes an R83P matrix.
c
c  Discussion:
c
c    The R83P storage format stores a periodic tridiagonal matrix as 
c    a 3 by N array, in which each row corresponds to a diagonal, and 
c    column locations are preserved.  The matrix value 
c    A(1,N) is stored as the array entry A(3,N), and the matrix value
c    A(N,1) is stored as the array entry A(1,1).
c
c  Example:
c
c    Here is how an R83P matrix of order 5 would be stored:
c
c      A51 A12 A23 A34 A45
c      A11 A22 A33 A44 A55
c      A21 A32 A43 A54 A15
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be at least 3.
c
c    Input/output, integer SEED, a seed for the random
c    number generator.
c
c    Output, double precision A(3,N), the R83P matrix.
c
      implicit none

      integer n

      double precision a(3,n)
      double precision r8_uniform_01
      integer i
      integer j
      integer seed

      do j = 1, n
        do i = 1, 3
          a(i,j) = r8_uniform_01 ( seed )
        end do
      end do

      return
      end
      subroutine r83p_sl ( n, a_lu, b, x, job, work2, work3, work4 )

c*********************************************************************72
c
cc R83P_SL solves an R83P system.
c
c  Discussion:
c
c    The R83P storage format stores a periodic tridiagonal matrix as 
c    a 3 by N array, in which each row corresponds to a diagonal, and 
c    column locations are preserved.  The matrix value 
c    A(1,N) is stored as the array entry A(3,N), and the matrix value
c    A(N,1) is stored as the array entry A(1,1).
c
c    The linear system must have been factored by R83P_FA.
c
c  Example:
c
c    Here is how an R83P matrix of order 5 would be stored:
c
c      A51 A12 A23 A34 A45
c      A11 A22 A33 A44 A55
c      A21 A32 A43 A54 A15
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be at least 3.
c
c    Input, double precision A_LU(3,N), the LU factors from R83P_FA.
c
c    Input, double precision B(N), the right hand side of the linear system.
c
c    Output, double precision X(N), the solution to the linear system.
c
c    Input, integer JOB, specifies the system to solve.
c    0, solve A * x = b.
c    nonzero, solve A' * x = b.
c
c    Input, double precision WORK2(N-1), WORK3(N-1), WORK4, 
c    factor data from R83P_FA.
c
      implicit none

      integer n

      double precision a_lu(3,n)
      double precision b(n)
      integer i
      integer job
      double precision work2(n-1)
      double precision work3(n-1)
      double precision work4
      double precision x(n)

      do i = 1, n
        x(i) = b(i)
      end do

      if ( job .eq. 0 ) then
c
c  Solve A1 * X1 = B1.
c
        call r83_np_sl ( n - 1, a_lu, x, job )
c
c  X2 = B2 - A3 * X1
c
        x(n) = x(n) - a_lu(1,1) * x(1) - a_lu(3,n-1) * x(n-1)
c
c  Solve A4 * X2 = X2
c
        x(n) = x(n) / work4
c
c  X1 := X1 - inverse ( A1 ) * A2 * X2.
c
        do i = 1, n - 1
          x(i) = x(i) - work2(i) * x(n)
        end do

      else
c
c  Solve A1' * X1 = B1.
c
        call r83_np_sl ( n - 1, a_lu, x, job )
c
c  X2 := X2 - A2' * B1
c
        x(n) = x(n) - a_lu(3,n) * x(1) - a_lu(1,n) * x(n-1)
c
c  Solve A4 * X2 = X2.
c
        x(n) = x(n) / work4
c
c  X1 := X1 - transpose ( inverse ( A1 ) * A3 ) * X2.
c
        do i = 1, n - 1
          x(i) = x(i) - work3(i) * x(n)
        end do

      end if

      return
      end
      subroutine r83p_to_r8ge ( n, a, b )

c*********************************************************************72
c
cc R83P_TO_R8GE copies an R83P matrix to an R8GE matrix.
c
c  Discussion:
c
c    The R83P storage format stores a periodic tridiagonal matrix as 
c    a 3 by N array, in which each row corresponds to a diagonal, and 
c    column locations are preserved.  The matrix value 
c    A(1,N) is stored as the array entry A(3,N), and the matrix value
c    A(N,1) is stored as the array entry A(1,1).
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c  Example:
c
c    Here is how an R83P matrix of order 5 would be stored:
c
c      A51 A12 A23 A34 A45
c      A11 A22 A33 A44 A55
c      A21 A32 A43 A54 A15
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be at least 3.
c
c    Input, double precision A(3,N), the R83P matrix.
c
c    Output, double precision B(N,N), the R8GE matrix.
c
      implicit none

      integer n

      double precision a(3,n)
      double precision b(n,n)
      integer i
      integer j

      do i = 1, n
        do j = 1, n
          if ( i .eq. j ) then
            b(i,j) = a(2,j)
          else if ( j .eq. i-1 ) then
            b(i,j) = a(3,j)
          else if ( j .eq. i+1 ) then
            b(i,j) = a(1,j)
          else if ( i .eq. 1 .and. j .eq. n ) then
            b(i,j) = a(3,j)
          else if ( i .eq. n .and. j .eq. 1 ) then
            b(i,j) = a(1,j)
          else
            b(i,j) = 0.0D+00
          end if
        end do
      end do

      return
      end
      subroutine r83p_vxm ( n, a, x, b )

c*********************************************************************72
c
cc R83P_VXM multiplies a vector by an R83P matrix.
c
c  Discussion:
c
c    The R83P storage format stores a periodic tridiagonal matrix as 
c    a 3 by N array, in which each row corresponds to a diagonal, and 
c    column locations are preserved.  The matrix value 
c    A(1,N) is stored as the array entry A(3,N), and the matrix value
c    A(N,1) is stored as the array entry A(1,1).
c
c  Example:
c
c    Here is how an R83P matrix of order 5 would be stored:
c
c      A51 A12 A23 A34 A45
c      A11 A22 A33 A44 A55
c      A21 A32 A43 A54 A15
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be at least 3.
c
c    Input, double precision A(3,N), the R83P matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(N), the product X * A.
c
      implicit none

      integer n

      double precision a(3,n)
      double precision b(n)
      integer i
      double precision x(n)

      b(1) =   a(1,1)   * x(n)   + a(2,1) * x(1) + a(3,1)   * x(2)

      do i = 2, n-1
        b(i) = a(1,i)   * x(i-1) + a(2,i) * x(i) + a(3,i)   * x(i+1)
      end do

      b(n) =   a(1,n)   * x(n-1) + a(2,n) * x(n) + a(3,n)   * x(1)

      return
      end
      subroutine r85_indicator ( n, a )

c*********************************************************************72
c
cc R85_INDICATOR sets up an R85 indicator matrix.
c
c  Discussion:
c
c    The "indicator matrix" simply has a value like I*10+J at every
c    entry of a dense matrix, or at every entry of a compressed storage
c    matrix for which storage is allocated. 
c
c    The R85 storage format represents a pentadiagonal matrix as a 5 
c    by N array, in which each row corresponds to a diagonal, and 
c    column locations are preserved.  Thus, the original matrix is 
c    "collapsed" vertically into the array.
c
c  Example:
c
c    Here is how an R85 matrix of order 6 would be stored:
c
c       *   *  A13 A24 A35 A46
c       *  A12 A23 A34 A45 A56
c      A11 A22 A33 A44 A55 A66
c      A21 A32 A43 A54 A65  *
c      A31 A42 A53 A64  *   *
c
c    Here are the values as stored in an indicator matrix:
c
c      00 00 13 24 35 46
c      00 12 23 34 45 56
c      11 22 33 44 55 66
c      21 32 43 54 65 00
c      31 42 53 64 00 00
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
c    Input, integer N, the order of the matrix.
c    N must be at least 3.
c
c    Output, double precision A(5,N), the indicator matrix.
c
      implicit none

      integer n

      double precision a(5,n)
      integer fac
      integer i
      integer i4_log_10
      integer j

      fac = 10 ** ( i4_log_10 ( n ) + 1 )

      a(1,1) = 0.0D+00
      a(1,2) = 0.0D+00
      do j = 3, n
        i = j - 2
        a(1,j) = dble ( fac * i + j )
      end do

      a(2,1) = 0.0D+00
      do j = 2, n
        i = j - 1
        a(2,j) = dble ( fac * i + j )
      end do

      do j = 1, n
        i = j
        a(3,j) = dble ( fac * i + j )
      end do

      do j = 1, n-1
        i = j + 1
        a(4,j) = dble ( fac * i + j )
      end do
      a(4,n) = 0.0D+00

      do j = 1, n-2
        i = j + 2
        a(5,j) = dble ( fac * i + j )
      end do
      a(5,n-1) = 0.0D+00
      a(5,n) = 0.0D+00

      return
      end
      subroutine r85_np_fs ( n, a, b, x )

c*********************************************************************72
c
cc R85_NP_FS factors and solves an R85 linear system.
c
c  Discussion:
c
c    The R85 storage format represents a pentadiagonal matrix as a 5 
c    by N array, in which each row corresponds to a diagonal, and 
c    column locations are preserved.  Thus, the original matrix is 
c    "collapsed" vertically into the array.
c
c    The factorization algorithm requires that each diagonal entry be nonzero.
c
c    No pivoting is performed, and therefore the algorithm may fail
c    in simple cases where the matrix is not singular.
c
c  Example:
c
c    Here is how an R85 matrix of order 6 would be stored:
c
c       *   *  A13 A24 A35 A46
c       *  A12 A23 A34 A45 A56
c      A11 A22 A33 A44 A55 A66
c      A21 A32 A43 A54 A65  *
c      A31 A42 A53 A64  *   *
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
c    Original FORTRAN77 version by Ward Cheney, David Kincaid.
c    This version by John Burkardt.
c
c  Reference:
c
c    Ward Cheney, David Kincaid,
c    Numerical Mathematics and Computing,
c    Brooks-Cole Publishing, 2004,
c    ISBN: 0534201121.
c
c  Parameters:
c
c    Input, integer N, the order of the linear system.
c
c    Input/output, double precision A(5,N),
c    On input, the pentadiagonal matrix.
c    On output, the data has been overwritten by factorization information.
c
c    Input/output, double precision B(N).
c    On input, B contains the right hand side of the linear system.
c    On output, B has been overwritten by factorization information.
c
c    Output, double precision X(N), the solution of the linear system.
c
      implicit none

      integer n

      double precision a(5,n)
      double precision b(n)
      integer i
      double precision x(n)
      double precision xmult

      do i = 1, n
        if ( a(3,i) .eq. 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R85_NP_FS - Fatal error!'
          write ( *, '(a,i8,a)' ) '  A(3,', i, ') = 0.'
          return
        end if
      end do

      do i = 2, n-1

        xmult = a(2,i) / a(3,i-1)
        a(3,i) = a(3,i) - xmult * a(4,i-1)
        a(4,i) = a(4,i) - xmult * a(5,i-1)

        b(i) = b(i) - xmult * b(i-1)

        xmult = a(1,i+1) / a(3,i-1)
        a(2,i+1) = a(2,i+1) - xmult * a(4,i-1)
        a(3,i+1) = a(3,i+1) - xmult * a(5,i-1)

        b(i+1) = b(i+1) - xmult * b(i-1)

      end do

      xmult = a(2,n) / a(3,n-1)
      a(3,n) = a(3,n) - xmult * a(4,n-1)

      x(n) = ( b(n) - xmult * b(n-1) ) / a(3,n)
      x(n-1) = ( b(n-1) - a(4,n-1) * x(n) ) / a(3,n-1)

      do i = n-2, 1, -1
        x(i) = ( b(i) - a(4,i) * x(i+1) - a(5,i) * x(i+2) ) / a(3,i)
      end do

      return
      end
      subroutine r85_mxv ( n, a, x, b )

c*********************************************************************72
c
cc R85_MXV multiplies an R85 matrix times a vector.
c
c  Discussion:
c
c    The R85 storage format represents a pentadiagonal matrix as a 5 
c    by N array, in which each row corresponds to a diagonal, and 
c    column locations are preserved.  Thus, the original matrix is 
c    "collapsed" vertically into the array.
c
c  Example:
c
c    Here is how an R85 matrix of order 6 would be stored:
c
c       *   *  A13 A24 A35 A46
c       *  A12 A23 A34 A45 A56
c      A11 A22 A33 A44 A55 A66
c      A21 A32 A43 A54 A65  *
c      A31 A42 A53 A64  *   *
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
c    Input, integer N, the order of the linear system.
c
c    Input, double precision A(5,N), the matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(N), the product A * x.
c
      implicit none

      integer n

      double precision a(5,n)
      double precision b(n)
      integer i
      double precision x(n)

      do i = 1, n
        b(i)   =            a(3,i)   * x(i)
      end do

      do i = 3, n
        b(i)   = b(i)   + a(1,i)   * x(i-2)
      end do

      do i = 2, n
        b(i)   = b(i)   + a(2,i)   * x(i-1)
      end do

      do i = 1, n - 1
        b(i) = b(i) + a(4,i) * x(i+1)
      end do

      do i = 1, n - 2
        b(i) = b(i) + a(5,i) * x(i+2)
      end do

      return
      end
      subroutine r85_print ( n, a, title )

c*********************************************************************72
c
cc R85_PRINT prints an R85 matrix.
c
c  Discussion:
c
c    The R85 storage format represents a pentadiagonal matrix as a 5 
c    by N array, in which each row corresponds to a diagonal, and 
c    column locations are preserved.  Thus, the original matrix is 
c    "collapsed" vertically into the array.
c
c  Example:
c
c    Here is how an R85 matrix of order 6 would be stored:
c
c       *   *  A13 A24 A35 A46
c       *  A12 A23 A34 A45 A56
c      A11 A22 A33 A44 A55 A66
c      A21 A32 A43 A54 A65  *
c      A31 A42 A53 A64  *   *
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
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input, double precision A(5,N), the matrix.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n

      double precision a(5,n)
      character * ( * ) title

      call r85_print_some ( n, a, 1, 1, n, n, title )

      return
      end
      subroutine r85_print_some ( n, a, ilo, jlo, ihi, jhi, title )

c*********************************************************************72
c
cc R85_PRINT_SOME prints some of an R85 matrix.
c
c  Discussion:
c
c    The R85 storage format represents a pentadiagonal matrix as a 5 
c    by N array, in which each row corresponds to a diagonal, and 
c    column locations are preserved.  Thus, the original matrix is 
c    "collapsed" vertically into the array.
c
c  Example:
c
c    Here is how an R85 matrix of order 6 would be stored:
c
c       *   *  A13 A24 A35 A46
c       *  A12 A23 A34 A45 A56
c      A11 A22 A33 A44 A55 A66
c      A21 A32 A43 A54 A65  *
c      A31 A42 A53 A64  *   *
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
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input, double precision A(5,N), the matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column, to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer, parameter :: incx = 5
      integer n

      double precision a(5,n)
      character * 14 ctemp(incx)
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

        write ( *, '(''  Col:  '',5a14)' ) ( ctemp(j2), j2 = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) '  ---'
c
c  Determine the range of the rows in this strip.
c
        i2lo = max ( ilo, 1 )
        i2lo = max ( i2lo, j2lo - 2 )

        i2hi = min ( ihi, n )
        i2hi = min ( i2hi, j2hi + 2 )

        do i = i2lo, i2hi
c
c  Print out (up to) 5 entries in row I, that lie in the current strip.
c
          do j2 = 1, inc

            j = j2lo - 1 + j2

            if ( 2 .lt. i-j .or. 2 .lt. j-i ) then
              ctemp(j2) = '              '
            else if ( j .eq. i+2 ) then
              if ( r8_is_int ( a(1,j) ) ) then
                write ( ctemp(j2), '(f8.0,6x)' ) a(1,j)
              else
                write ( ctemp(j2), '(g14.6)' ) a(1,j)
              end if
            else if ( j .eq. i+1 ) then
              if ( r8_is_int ( a(2,j) ) ) then
                write ( ctemp(j2), '(f8.0,6x)' ) a(2,j)
              else
                write ( ctemp(j2), '(g14.6)' ) a(2,j)
              end if
            else if ( j .eq. i ) then
              if ( r8_is_int ( a(3,j) ) ) then
                write ( ctemp(j2), '(f8.0,6x)' ) a(3,j)
              else
                write ( ctemp(j2), '(g14.6)' ) a(3,j)
              end if
            else if ( j .eq. i-1 ) then
              if ( r8_is_int ( a(4,j) ) ) then
                write ( ctemp(j2), '(f8.0,6x)' ) a(4,j)
              else
                write ( ctemp(j2), '(g14.6)' ) a(4,j)
              end if
            else if ( j .eq. i-2 ) then
              if ( r8_is_int ( a(5,j) ) ) then
                write ( ctemp(j2), '(f8.0,6x)' ) a(5,j)
              else
                write ( ctemp(j2), '(g14.6)' ) a(5,j)
              end if
            end if

          end do

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine r85_random ( n, seed, a )

c*********************************************************************72
c
cc R85_RANDOM randomizes an R85 matrix.
c
c  Discussion:
c
c    The R85 storage format represents a pentadiagonal matrix as a 5 
c    by N array, in which each row corresponds to a diagonal, and 
c    column locations are preserved.  Thus, the original matrix is 
c    "collapsed" vertically into the array.
c
c  Example:
c
c    Here is how an R85 matrix of order 6 would be stored:
c
c       *   *  A13 A24 A35 A46
c       *  A12 A23 A34 A45 A56
c      A11 A22 A33 A44 A55 A66
c      A21 A32 A43 A54 A65  *
c      A31 A42 A53 A64  *   *
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
c    Input, integer N, the order of the linear system.
c
c    Input/output, integer SEED, a seed for the random 
c    number generator.
c
c    Output, double precision A(5,N), the matrix.
c
      implicit none

      integer n

      double precision a(5,n)
      double precision r8_uniform_01
      integer j
      integer seed

      a(1,1) = 0.0D+00
      a(1,2) = 0.0D+00
      do j = 3, n
        a(1,j) = r8_uniform_01 ( seed )
      end do

      a(2,1) = 0.0D+00
      do j = 2, n
        a(2,j) = r8_uniform_01 ( seed )
      end do

      do j = 1, n
        a(3,j) = r8_uniform_01 ( seed )
      end do

      do j = 1, n-1
        a(4,j) = r8_uniform_01 ( seed )
      end do

      do j = 1, n-2
        a(5,j) = r8_uniform_01 ( seed )
      end do

      return
      end
      subroutine r85_to_r8ge ( n, a, b )

c*********************************************************************72
c
cc R85_TO_R8GE copies an R85 matrix into an R8GE matrix.
c
c  Discussion:
c
c    The R85 storage format represents a pentadiagonal matrix as a 5 
c    by N array, in which each row corresponds to a diagonal, and 
c    column locations are preserved.  Thus, the original matrix is 
c    "collapsed" vertically into the array.
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c  Example:
c
c    Here is how an R85 matrix of order 6 would be stored:
c
c       *   *  A13 A24 A35 A46
c       *  A12 A23 A34 A45 A56
c      A11 A22 A33 A44 A55 A66
c      A21 A32 A43 A54 A65  *
c      A31 A42 A53 A64  *   *
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
c    Input, integer N, the order of the matrix.
c    N must be at least 3.
c
c    Input, double precision A(5,N), the R85 matrix.
c
c    Output, double precision A(N,N), the R8GE matrix.
c
      implicit none

      integer n

      double precision a(5,n)
      double precision b(n,n)
      integer i
      integer j

      do i = 1, n
        do j = 1, n

          if ( j .eq. i-2 ) then
            b(i,j) = a(1,i)
          else if ( j .eq. i-1 ) then
            b(i,j) = a(2,i)
          else if ( i .eq. j ) then
            b(i,j) = a(3,i)
          else if ( j .eq. i+1 ) then
            b(i,j) = a(4,i)
          else if ( j .eq. i+2 ) then
            b(i,j) = a(5,i)
          else
            b(i,j) = 0.0D+00
          end if

        end do
      end do

      return
      end
      subroutine r85_vxm ( n, a, x, b )

c*********************************************************************72
c
cc R85_VXM multiplies a vector by an R85 matrix.
c
c  Discussion:
c
c    The R85 storage format represents a pentadiagonal matrix as a 5 
c    by N array, in which each row corresponds to a diagonal, and 
c    column locations are preserved.  Thus, the original matrix is 
c    "collapsed" vertically into the array.
c
c  Example:
c
c    Here is how an R85 matrix of order 6 would be stored:
c
c       *   *  A13 A24 A35 A46
c       *  A12 A23 A34 A45 A56
c      A11 A22 A33 A44 A55 A66
c      A21 A32 A43 A54 A65  *
c      A31 A42 A53 A64  *   *
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
c    Input, integer N, the order of the linear system.
c
c    Input, double precision A(5,N), the matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A'.
c
c    Output, double precision B(N), the product A' * x.
c
      implicit none

      integer n

      double precision a(5,n)
      double precision b(n)
      integer i
      double precision x(n)

      do i = 1, n
        b(i)   =            a(3,i)   * x(i)
      end do

      do i = 2, n
        b(i)   = b(i)   + a(4,i-1) * x(i-1)
      end do

      do i = 3, n
        b(i)   = b(i)   + a(5,i-2) * x(i-2)
      end do

      do i = 1, n - 1
        b(i) = b(i) + a(2,i+1)   * x(i+1)
      end do

      do i = 1, n - 2
        b(i) = b(i) + a(1,i+2)   * x(i+2)
      end do

      return
      end
      subroutine r8bb_add ( n1, n2, ml, mu, a, i, j, value )

c*********************************************************************72
c
cc R8BB_ADD adds a value to an entry in an R8BB matrix.
c
c  Discussion:
c
c    The R8BB storage format is for a border banded matrix.  Such a
c    matrix has the logical form:
c
c      A1 | A2
c      ---+---
c      A3 | A4
c
c    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
c    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
c    respectively.
c
c    A should be defined as a vector.  The user must then store
c    the entries of the four blocks of the matrix into the vector A.
c    Each block is stored by columns.
c
c    A1, the banded portion of the matrix, is stored in
c    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
c    general band format.  The reason for the factor of 2 in front of
c    ML is to allocate space that may be required if pivoting occurs.
c
c    The following formulas should be used to determine how to store
c    the entry corresponding to row I and column J in the original matrix:
c
c    Entries of A1:
c
c      1 .le. I .le. N1, 1 .le. J .le. N1, (J-I) .le. MU and (I-J) .le. ML.
c
c      Store the I, J entry into location
c      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
c
c    Entries of A2:
c
c      1 .le. I .le. N1, N1+1 .le. J .le. N1+N2.
c
c      Store the I, J entry into location
c      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
c
c    Entries of A3:
c
c      N1+1 .le. I .le. N1+N2, 1 .le. J .le. N1.
c
c      Store the I, J entry into location
c      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c
c    Entries of A4:
c
c      N1+1 .le. I .le. N1+N2, N1+1 .le. J .le. N1+N2
c
c      Store the I, J entry into location
c      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c      (same formula used for A3).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, the order of the banded and dense 
c    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than N1-1.
c
c    Input/output, double precision A((2*ML+MU+1)*N1+2*N1*N2+N2*N2),
c    the R8BB matrix.
c
c    Input, integer I, J, the row and column of the entry to be
c    incremented.  Some combinations of I and J are illegal.
c
c    Input, double precision VALUE, the value to be added to the 
c    (I,J)-th entry.
c
      implicit none

      integer ml
      integer mu
      integer n1
      integer n2

      double precision a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
      integer i
      integer ij
      integer j
      double precision value

      if ( value .eq. 0.0D+00 ) then
        return
      end if

      if ( i .le. 0 .or. n1 + n2 .lt. i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8BB_ADD - Fatal error!'
        write ( *, '(a,i8)' ) 
     &  '  Illegal input value of row index I = ', i
        stop
      end if

      if ( j .le. 0 .or. n1 + n2 .lt. j ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8BB_ADD - Fatal error!'
        write ( *, '(a,i8)' ) 
     &  '  Illegal input value of column index J = ', j
        stop
      end if
c
c  The A1 block of the matrix.
c
c  Check for out of band problems.
c
c  Normally, we would check the condition MU .lt. (J-I), but the storage
c  format requires extra entries be set aside in case of pivoting, which
c  means that the condition becomes MU+ML .lt. (J-I).
c
      if ( i .le. n1 .and. j .le. n1 ) then
        if ( (mu+ml) .lt. (j-i) .or. ml .lt. (i-j) ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8BB_ADD - Warning!'
          write ( *, '(a,i8,a,i8,a)' ) 
     &      '  Unable to add to entry (', i, ',', j, ').'
        else
          ij = (i-j+ml+mu+1)+(j-1)*(2*ml+mu+1)
        end if
c
c  The A2 block of the matrix.
c
      else if ( i .le. n1 .and. n1 .lt. j ) then
        ij = (2*ml+mu+1)*n1+(j-n1-1)*n1 + i
c
c  The A3 and A4 blocks of the matrix.
c
      else if ( n1 .lt. i ) then
        ij = (2*ml+mu+1)*n1+n2*n1+(j-1)*n2 + (i-n1)
      end if

      a(ij) = a(ij) + value

      return
      end
      subroutine r8bb_fa ( n1, n2, ml, mu, a, pivot, info )

c*********************************************************************72
c
cc R8BB_FA factors an R8BB matrix.
c
c  Discussion:
c
c    The R8BB storage format is for a border banded matrix.  Such a
c    matrix has the logical form:
c
c      A1 | A2
c      ---+---
c      A3 | A4
c
c    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
c    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
c    respectively.
c
c    A should be defined as a vector.  The user must then store
c    the entries of the four blocks of the matrix into the vector A.
c    Each block is stored by columns.
c
c    A1, the banded portion of the matrix, is stored in
c    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
c    general band format.  The reason for the factor of 2 in front of
c    ML is to allocate space that may be required if pivoting occurs.
c
c    The following formulas should be used to determine how to store
c    the entry corresponding to row I and column J in the original matrix:
c
c    Entries of A1:
c
c      1 .le. I .le. N1, 1 .le. J .le. N1, (J-I) .le. MU and (I-J) .le. ML.
c
c      Store the I, J entry into location
c      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
c
c    Entries of A2:
c
c      1 .le. I .le. N1, N1+1 .le. J .le. N1+N2.
c
c      Store the I, J entry into location
c      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
c
c    Entries of A3:
c
c      N1+1 .le. I .le. N1+N2, 1 .le. J .le. N1.
c
c      Store the I, J entry into location
c      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c
c    Entries of A4:
c
c      N1+1 .le. I .le. N1+N2, N1+1 .le. J .le. N1+N2
c
c      Store the I, J entry into location
c      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c      (same formula used for A3).
c
c    Once the matrix has been factored by R8BB_FA, R8BB_SL may be called
c    to solve linear systems involving the matrix.
c
c    R8BB_FA uses LINPACK routines to carry out the factorization.
c
c
c    The linear system must be border banded, of the form:
c
c      ( A1 A2 ) (X1) = (B1)
c      ( A3 A4 ) (X2)   (B2)
c
c    where A1 is a (usually big) banded square matrix, A2 and A3 are
c    column and row strips which may be nonzero, and A4 is a dense
c    square matrix.
c
c    The algorithm rewrites the system as:
c
c         X1 + inv(A1) A2 X2 = inv(A1) B1
c
c      A3 X1 +         A4 X2 = B2
c
c    and then rewrites the second equation as
c
c      ( A4 - A3 inv(A1) A2 ) X2 = B2 - A3 inv(A1) B1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, the order of the banded and dense
c    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative and no greater than N1-1.
c
c    Input/output, double precision A( (2*ML+MU+1)*N1 + 2*N1*N2 + N2*N2 ).
c    On input, the border-banded matrix to be factored.
c    On output, information describing a partial factorization
c    of the original coefficient matrix.  This information is required
c    by R8BB_SL in order to solve linear systems associated with that
c    matrix.
c
c    Output, integer PIVOT(N1+N2), contains pivoting information.
c
c    Output, integer INFO, singularity flag.
c    0, no singularity detected.
c    nonzero, the factorization failed on the INFO-th step.
c
      implicit none

      integer ml
      integer mu
      integer n1
      integer n2

      double precision a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
      integer i
      integer ij
      integer ik
      integer info
      integer pivot(n1+n2)
      integer j
      integer jk
      integer job
      integer k
      integer nband

      nband = (2*ml+mu+1) * n1
c
c  Factor the A1 band matrix, overwriting A1 by its factors.
c
      if ( 0 .lt. n1 ) then

        call r8gb_fa ( n1, ml, mu, a, pivot, info )

        if ( info .ne. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8BB_FA - Fatal error!'
          write ( *, '(a,i8)' ) '  R8GB_FA returned INFO = ', info
          write ( *, '(a)' ) '  Factoring failed for column INFO.'
          write ( *, '(a)' ) '  The band matrix A1 is singular.'
          write ( *, '(a)' ) '  This algorithm cannot continue!'
          stop
        end if

      end if

      if ( 0 .lt. n1 .and. 0 .lt. n2 ) then
c
c  Solve A1 * x = -A2 for x, and overwrite A2 by the results.
c
        do i = nband+1, nband+n1*n2
          a(i) = -a(i)
        end do

        job = 0
        do j = 1, n2
          call r8gb_sl ( n1, ml, mu, a, pivot, a(nband+(j-1)*n1+1), 
     &      job )
        end do
c
c  A4 := A4 + A3 * A2.
c
        do i = 1, n2
          do j = 1, n1
            ij = nband + n1*n2 + (j-1)*n2 + i
            do k = 1, n2
              ik = nband + 2*n1*n2 + (k-1)*n2 + i
              jk = nband + (k-1) * n1 + j
              a(ik) = a(ik) + a(ij) * a(jk)
            end do
          end do
        end do

      end if
c
c  Factor A4.
c
      if ( 0 .lt. n2 ) then

        call r8ge_fa ( n2, a(nband+2*n1*n2+1), pivot(n1+1), info )

        if ( info .ne. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8BB_FA - Fatal error!'
          write ( *, '(a,i8)' ) '  R8GE_FA returned INFO = ',info
          write ( *, '(a)' ) 
     &      '  This indicates singularity in column INFO.'
          write ( *, '(a,i8)' ) 
     &      '  of the A4 submatrix, which is column ', n1+info
          write ( *, '(a)' ) '  of the full matrix.'
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 
     &      '  It is possible that the full matrix is '
          write ( *, '(a)' ) 
     &      '  nonsingular, but the algorithm R8BB_FA may'
          write ( *, '(a)' ) '  not be used for this matrix.'
          stop
        end if
      end if

      return
      end
      subroutine r8bb_get ( n1, n2, ml, mu, a, i, j, value )

c*********************************************************************72
c
cc R8BB_GET returns an entry of an R8BB matrix.
c
c  Discussion:
c
c    The R8BB storage format is for a border banded matrix.  Such a
c    matrix has the logical form:
c
c      A1 | A2
c      ---+---
c      A3 | A4
c
c    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
c    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
c    respectively.
c
c    A should be defined as a vector.  The user must then store
c    the entries of the four blocks of the matrix into the vector A.
c    Each block is stored by columns.
c
c    A1, the banded portion of the matrix, is stored in
c    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
c    general band format.  The reason for the factor of 2 in front of
c    ML is to allocate space that may be required if pivoting occurs.
c
c    The following formulas should be used to determine how to store
c    the entry corresponding to row I and column J in the original matrix:
c
c    Entries of A1:
c
c      1 .le. I .le. N1, 1 .le. J .le. N1, (J-I) .le. MU and (I-J) .le. ML.
c
c      Store the I, J entry into location
c      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
c
c    Entries of A2:
c
c      1 .le. I .le. N1, N1+1 .le. J .le. N1+N2.
c
c      Store the I, J entry into location
c      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
c
c    Entries of A3:
c
c      N1+1 .le. I .le. N1+N2, 1 .le. J .le. N1.
c
c      Store the I, J entry into location
c      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c
c    Entries of A4:
c
c      N1+1 .le. I .le. N1+N2, N1+1 .le. J .le. N1+N2
c
c      Store the I, J entry into location
c      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c      (same formula used for A3).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, the order of the banded and dense
c    blocks. N1 and N2 must be nonnegative, and at least one must be positive.
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than N1-1.
c
c    Input, double precision A((2*ML+MU+1)*N1+2*N1*N2+N2*N2), the R8BB matrix.
c
c    Input, integer I, J, the row and column of the entry to 
c    be retrieved.
c
c    Output, double precision VALUE, the value of the (I,J) entry.
c
      implicit none

      integer ml
      integer mu
      integer n1
      integer n2

      double precision a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
      integer i
      integer ij
      integer j
      double precision value

      if ( i .le. 0 .or. n1+n2 .lt. i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8BB_GET - Fatal error!'
        write ( *, '(a,i8)' ) 
     &    '  Illegal input value of row index I = ', i
        stop
      end if

      if ( j .le. 0 .or. n1+n2 .lt. j ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8BB_GET - Fatal error!'
        write ( *, '(a,i8)' ) 
     &    '  Illegal input value of column index J = ', j
        stop
      end if
c
c  The A1 block of the matrix.
c
c  Check for out of band problems.
c
c  Normally, we would check the condition MU .lt. (J-I), but the storage
c  format requires extra entries be set aside in case of pivoting, which
c  means that the condition becomes MU+ML .lt. (J-I).
c
      if ( i .le. n1 .and. j .le. n1 ) then
        if ( mu+ml .lt. (j-i) .or. ml .lt. (i-j) ) then
          value = 0.0D+00
          return
        else
          ij = (i-j+ml+mu+1)+(j-1)*(2*ml+mu+1)
        end if
c
c  The A2 block of the matrix.
c
      else if ( i .le. n1 .and. n1 .lt. j ) then
        ij = (2*ml+mu+1)*n1+(j-n1-1)*n1+i
c
c  The A3 and A4 blocks of the matrix.
c
      else if ( n1 .lt. i ) then
        ij = (2*ml+mu+1)*n1+n2*n1+(j-1)*n2+(i-n1)
      end if

      value = a(ij)

      return
      end
      subroutine r8bb_indicator ( n1, n2, ml, mu, a )

c*********************************************************************72
c
cc R8BB_INDICATOR sets up an R8BB indicator matrix.
c
c  Discussion:
c
c    The "indicator matrix" simply has a value like I*10+J at every
c    entry of a dense matrix, or at every entry of a compressed storage
c    matrix for which storage is allocated. 
c
c    The R8BB storage format is for a border banded matrix.  Such a
c    matrix has the logical form:
c
c      A1 | A2
c      ---+---
c      A3 | A4
c
c    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
c    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
c    respectively.
c
c  Example:
c
c    With N1 = 4, N2 = 1, ML = 1, MU = 2, the matrix entries would be:
c
c       00
c       00  00
c       00  00  00 --- ---
c      A11 A12 A13  00 ---  A16 A17
c      A21 A22 A23 A24  00  A26 A27
c      --- A32 A33 A34 A35  A36 A37
c      --- --- A43 A44 A45  A46 A47
c      --- --- --- A54 A55  A56 A57
c                       00
c
c      A61 A62 A63 A64 A65  A66 A67
c      A71 A72 A73 A74 A75  A76 A77
c
c    The matrix is actually stored as a vector, and we will simply suggest
c    the structure and values of the indicator matrix as:
c
c      00 00 00 00 00
c      00 00 13 24 35     16 17     61 62 63 64 65     66 67
c      00 12 23 34 45  +  26 27  +  71 72 73 74 75  +  76 77
c      11 22 33 44 55     36 37     
c      21 32 43 54 00     46 47     
c                         56 57     
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, the order of the banded and dense 
c    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative and no greater than N1-1.
c
c    Output, double precision A((2*ML+MU+1)*N1+2*N1*N2+N2*N2), the R8BB matrix.
c
      implicit none

      integer ml
      integer mu
      integer n1
      integer n2

      double precision a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
      integer base
      integer fac
      integer i
      integer i4_log_10
      integer j
      integer row

      fac = 10 ** ( i4_log_10 ( n1 + n2 ) + 1 )
c
c  Set the banded matrix A1.
c
      do j = 1, n1
        do row = 1, 2 * ml + mu + 1
          i = row + j - ml - mu - 1
          if ( ml .lt. row .and. 1 .le. i .and. i .le. n1 ) then
            a(row+(j-1)*(2*ml+mu+1)) = dble ( fac * i + j )
          else
            a(row+(j-1)*(2*ml+mu+1)) = 0.0D+00
          end if
        end do
      end do
c
c  Set the N1 by N2 rectangular strip A2.
c
      base = ( 2 * ml + mu + 1 ) * n1

      do i = 1, n1
        do j = n1 + 1, n1 + n2
          a(base + i + (j-n1-1)*n1 ) = dble ( fac * i + j )
        end do
      end do
c
c  Set the N2 by N1 rectangular strip A3.
c
      base = ( 2 * ml + mu + 1 ) * n1 + n1 * n2

      do i = n1 + 1, n1 + n2
        do j = 1, n1    
          a(base + i-n1 + (j-1)*n2 ) = dble ( fac * i + j )
        end do
      end do
c
c  Set the N2 by N2 square A4.
c
      base = ( 2 * ml + mu + 1 ) * n1 + n1 * n2 + n2 * n1

      do i = n1 + 1, n1 + n2
        do j = n1 + 1, n1 + n2
          a(base + i-n1 + (j-n1-1)*n2 ) = dble ( fac * i + j )
        end do
      end do

      return
      end
      subroutine r8bb_mxv ( n1, n2, ml, mu, a, x, b )

c*********************************************************************72
c
cc R8BB_MXV multiplies an R8BB matrix by an R8VEC.
c
c  Discussion:
c
c    The R8BB storage format is for a border banded matrix.  Such a
c    matrix has the logical form:
c
c      A1 | A2
c      ---+---
c      A3 | A4
c
c    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
c    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
c    respectively.
c
c    A should be defined as a vector.  The user must then store
c    the entries of the four blocks of the matrix into the vector A.
c    Each block is stored by columns.
c
c    A1, the banded portion of the matrix, is stored in
c    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
c    general band format.  The reason for the factor of 2 in front of
c    ML is to allocate space that may be required if pivoting occurs.
c
c    The following formulas should be used to determine how to store
c    the entry corresponding to row I and column J in the original matrix:
c
c    Entries of A1:
c
c      1 .le. I .le. N1, 1 .le. J .le. N1, (J-I) .le. MU and (I-J) .le. ML.
c
c      Store the I, J entry into location
c      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
c
c    Entries of A2:
c
c      1 .le. I .le. N1, N1+1 .le. J .le. N1+N2.
c
c      Store the I, J entry into location
c      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
c
c    Entries of A3:
c
c      N1+1 .le. I .le. N1+N2, 1 .le. J .le. N1.
c
c      Store the I, J entry into location
c      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c
c    Entries of A4:
c
c      N1+1 .le. I .le. N1+N2, N1+1 .le. J .le. N1+N2
c
c      Store the I, J entry into location
c      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c      (same formula used for A3).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, the order of the banded and dense
c    blocks  N1 and N2 must be nonnegative, and at least one must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative and no greater than N1-1.
c
c    Input, double precision A((2*ML+MU+1)*N1+2*N1*N2+N2*N2), the R8BB matrix.
c
c    Input, double precision X(N1+N2), the vector to be multiplied by A.
c
c    Output, double precision B(N1+N2), the result of multiplying A by X.
c
      implicit none

      integer ml
      integer mu
      integer n1
      integer n2

      double precision a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
      double precision b(n1+n2)
      integer ihi
      integer ij
      integer ilo
      integer j
      double precision x(n1+n2)
c
c  Initialize B.
c
      b(1:n1+n2) = 0.0D+00
c
c  Multiply by A1.
c
      do j = 1, n1

        ilo = max ( 1, j - mu - ml )
        ihi = min ( n1, j + ml )
        ij = (j-1) * (2*ml+mu+1) - j + ml + mu + 1

        b(ilo:ihi) = b(ilo:ihi) + a(ij+ilo:ij+ihi) * x(j)

      end do
c
c  Multiply by A2.
c
      do j = n1+1, n1+n2
        ij = (2*ml+mu+1)*n1+(j-n1-1)*n1

        b(1:n1) = b(1:n1) + a(ij+1:ij+n1) * x(j)

      end do
c
c  Multiply by A3 and A4.
c
      do j = 1, n1+n2
        ij = (2*ml+mu+1)*n1+n1*n2+(j-1)*n2-n1

        b(n1+1:n1+n2) = b(n1+1:n1+n2) + a(ij+n1+1:ij+n1+n2) * x(j)

      end do

      return
      end
      subroutine r8bb_print ( n1, n2, ml, mu, a, title )

c*********************************************************************72
c
cc R8BB_PRINT prints an R8BB matrix.
c
c  Discussion:
c
c    The R8BB storage format is for a border banded matrix.  Such a
c    matrix has the logical form:
c
c      A1 | A2
c      ---+---
c      A3 | A4
c
c    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
c    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
c    respectively.
c
c    A should be defined as a vector.  The user must then store
c    the entries of the four blocks of the matrix into the vector A.
c    Each block is stored by columns.
c
c    A1, the banded portion of the matrix, is stored in
c    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
c    general band format.  The reason for the factor of 2 in front of
c    ML is to allocate space that may be required if pivoting occurs.
c
c    The following formulas should be used to determine how to store
c    the entry corresponding to row I and column J in the original matrix:
c
c    Entries of A1:
c
c      1 .le. I .le. N1, 1 .le. J .le. N1, (J-I) .le. MU and (I-J) .le. ML.
c
c      Store the I, J entry into location
c      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
c
c    Entries of A2:
c
c      1 .le. I .le. N1, N1+1 .le. J .le. N1+N2.
c
c      Store the I, J entry into location
c      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
c
c    Entries of A3:
c
c      N1+1 .le. I .le. N1+N2, 1 .le. J .le. N1.
c
c      Store the I, J entry into location
c      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c
c    Entries of A4:
c
c      N1+1 .le. I .le. N1+N2, N1+1 .le. J .le. N1+N2
c
c      Store the I, J entry into location
c      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c      (same formula used for A3).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, the order of the banded and dense blocks.
c    N1 and N2 must be nonnegative, and at least one must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than N1-1.
c
c    Input, double precision A((2*ML+MU+1)*N1+2*N1*N2+N2*N2), the R8BB matrix.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer ml
      integer mu
      integer n1
      integer n2

      double precision a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
      character * ( * )  title

      call r8bb_print_some ( n1, n2, ml, mu, a, 1, 1, n1+n2, n1+n2, 
     &  title )

      return
      end
      subroutine r8bb_print_some ( n1, n2, ml, mu, a, ilo, jlo, ihi, 
     &  jhi, title )

c*********************************************************************72
c
cc R8BB_PRINT_SOME prints some of an R8BB matrix.
c
c  Discussion:
c
c    The R8BB storage format is for a border banded matrix.  Such a
c    matrix has the logical form:
c
c      A1 | A2
c      ---+---
c      A3 | A4
c
c    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
c    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
c    respectively.
c
c    A should be defined as a vector.  The user must then store
c    the entries of the four blocks of the matrix into the vector A.
c    Each block is stored by columns.
c
c    A1, the banded portion of the matrix, is stored in
c    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
c    general band format.  The reason for the factor of 2 in front of
c    ML is to allocate space that may be required if pivoting occurs.
c
c    The following formulas should be used to determine how to store
c    the entry corresponding to row I and column J in the original matrix:
c
c    Entries of A1:
c
c      1 .le. I .le. N1, 1 .le. J .le. N1, (J-I) .le. MU and (I-J) .le. ML.
c
c      Store the I, J entry into location
c      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
c
c    Entries of A2:
c
c      1 .le. I .le. N1, N1+1 .le. J .le. N1+N2.
c
c      Store the I, J entry into location
c      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
c
c    Entries of A3:
c
c      N1+1 .le. I .le. N1+N2, 1 .le. J .le. N1.
c
c      Store the I, J entry into location
c      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c
c    Entries of A4:
c
c      N1+1 .le. I .le. N1+N2, N1+1 .le. J .le. N1+N2
c
c      Store the I, J entry into location
c      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c      (same formula used for A3).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, the order of the banded and dense 
c    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than N1-1.
c
c    Input, double precision A((2*ML+MU+1)*N1+2*N1*N2+N2*N2), the R8BB matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 5 )
      integer ml
      integer mu
      integer n1
      integer n2

      double precision a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
      double precision aij
      character * ( 14 ) ctemp(incx)
      logical r8_is_int
      integer i
      integer i2hi
      integer i2lo
      integer ihi
      integer ij
      integer ilo
      integer inc
      integer j
      integer j2
      integer j2hi
      integer j2lo
      integer jhi
      integer jlo
      character * ( * )  title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
c
c  Print the columns of the matrix, in strips of 5.
c
      do j2lo = jlo, jhi, incx

        j2hi = j2lo + incx - 1
        j2hi = min ( j2hi, n1+n2 )
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
        i2hi = min ( ihi, n1+n2 )

        do i = i2lo, i2hi
c
c  Print out (up to) 5 entries in row I, that lie in the current strip.
c
          do j2 = 1, inc

            j = j2lo - 1 + j2

            aij = 0.0D+00

            if ( i .le. n1 .and. j .le. n1 ) then
              if ( (j-i) .le. mu+ml .and. (i-j) .le. ml ) then
                ij = (i-j+ml+mu+1)+(j-1)*(2*ml+mu+1)
                aij = a(ij)
              end if
            else if ( i .le. n1 .and. n1 .lt. j ) then
              ij = (2*ml+mu+1)*n1+(j-n1-1)*n1+i
              aij = a(ij)
            else if ( n1 .lt. i ) then
              ij = (2*ml+mu+1)*n1+n2*n1+(j-1)*n2+(i-n1)
              aij = a(ij)
            end if

            if ( r8_is_int ( aij ) ) then
              write ( ctemp(j2), '(f8.0,6x)' ) aij
            else
              write ( ctemp(j2), '(g14.6)' ) aij
            end if

          end do

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine r8bb_random ( n1, n2, ml, mu, seed, a )

c*********************************************************************72
c
cc R8BB_RANDOM randomizes an R8BB matrix.
c
c  Discussion:
c
c    The R8BB storage format is for a border banded matrix.  Such a
c    matrix has the logical form:
c
c      A1 | A2
c      ---+---
c      A3 | A4
c
c    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
c    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
c    respectively.
c
c    A should be defined as a vector.  The user must then store
c    the entries of the four blocks of the matrix into the vector A.
c    Each block is stored by columns.
c
c    A1, the banded portion of the matrix, is stored in
c    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
c    general band format.  The reason for the factor of 2 in front of
c    ML is to allocate space that may be required if pivoting occurs.
c
c    The following formulas should be used to determine how to store
c    the entry corresponding to row I and column J in the original matrix:
c
c    Entries of A1:
c
c      1 .le. I .le. N1, 1 .le. J .le. N1, (J-I) .le. MU and (I-J) .le. ML.
c
c      Store the I, J entry into location
c      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
c
c    Entries of A2:
c
c      1 .le. I .le. N1, N1+1 .le. J .le. N1+N2.
c
c      Store the I, J entry into location
c      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
c
c    Entries of A3:
c
c      N1+1 .le. I .le. N1+N2, 1 .le. J .le. N1.
c
c      Store the I, J entry into location
c      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c
c    Entries of A4:
c
c      N1+1 .le. I .le. N1+N2, N1+1 .le. J .le. N1+N2
c
c      Store the I, J entry into location
c      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c      (same formula used for A3).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, the order of the banded and dense blocks.
c    N1 and N2 must be nonnegative, and at least one must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative and no greater than N1-1.
c
c    Input/output, integer SEED, a seed for the random 
c    number generator.
c
c    Output, double precision A((2*ML+MU+1)*N1+2*N1*N2+N2*N2), the R8BB matrix.
c
      implicit none

      integer ml
      integer mu
      integer n1
      integer n2

      double precision a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
      double precision r8_uniform_01
      integer i
      integer ilo
      integer j
      double precision r
      integer row
      integer seed
c
c  Randomize the banded matrix A1.
c  We still believe that the "junk" entries should be set to 0.
c
      do j = 1, n1
        do row = 1, 2*ml+mu+1
          i = row + j - ml - mu - 1
          if ( ml .lt. row .and. 1 .le. i .and. i .le. n1 ) then
            r = r8_uniform_01 ( seed )
          else
            r = 0.0D+00
          end if
          a(row+(j-1)*(2*ml+mu+1)) = r
        end do
      end do
c
c  Randomize the rectangular strips A2+A3+A4.
c
      ilo = (2*ml+mu+1) * n1 + 1

      call r8vec_uniform_01 ( n1*n2+n2*n1+n2*n2, seed, a(ilo:) )

      return
      end
      subroutine r8bb_set ( n1, n2, ml, mu, a, i, j, value )

c*********************************************************************72
c
cc R8BB_SET sets an entry of an R8BB matrix.
c
c  Discussion:
c
c    The R8BB storage format is for a border banded matrix.  Such a
c    matrix has the logical form:
c
c      A1 | A2
c      ---+---
c      A3 | A4
c
c    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
c    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
c    respectively.
c
c    A should be defined as a vector.  The user must then store
c    the entries of the four blocks of the matrix into the vector A.
c    Each block is stored by columns.
c
c    A1, the banded portion of the matrix, is stored in
c    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
c    general band format.  The reason for the factor of 2 in front of
c    ML is to allocate space that may be required if pivoting occurs.
c
c    The following formulas should be used to determine how to store
c    the entry corresponding to row I and column J in the original matrix:
c
c    Entries of A1:
c
c      1 .le. I .le. N1, 1 .le. J .le. N1, (J-I) .le. MU and (I-J) .le. ML.
c
c      Store the I, J entry into location
c      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
c
c    Entries of A2:
c
c      1 .le. I .le. N1, N1+1 .le. J .le. N1+N2.
c
c      Store the I, J entry into location
c      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
c
c    Entries of A3:
c
c      N1+1 .le. I .le. N1+N2, 1 .le. J .le. N1.
c
c      Store the I, J entry into location
c      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c
c    Entries of A4:
c
c      N1+1 .le. I .le. N1+N2, N1+1 .le. J .le. N1+N2
c
c      Store the I, J entry into location
c      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c      (same formula used for A3).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, the order of the banded and dense blocks.
c    N1 and N2 must be nonnegative, and at least one must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than N1-1.
c
c    Input/output, double precision A((2*ML+MU+1)*N1+2*N1*N2+N2*N2),
c    the R8BB matrix.
c
c    Input, integer I, J, the row and column of the entry to be set.
c
c    Input, double precision VALUE, the value to be assigned to the
c    (I,J) entry.
c
      implicit none

      integer ml
      integer mu
      integer n1
      integer n2

      double precision a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
      integer i
      integer ij
      integer j
      double precision value

      if ( i .le. 0 .or. n1+n2 .lt. i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8BB_SET - Fatal error!'
        write ( *, '(a,i8)' ) 
     &    '  Illegal input value of row index I = ', i
        stop
      end if

      if ( j .le. 0 .or. n1+n2 .lt. j ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8BB_SET - Fatal error!'
        write ( *, '(a,i8)' ) 
     &    '  Illegal input value of column index J = ', j
        stop
      end if
c
c  The A1 block of the matrix.
c
c  Check for out of band problems.
c
c  Normally, we would check the condition MU .lt. (J-I), but the storage
c  format requires extra entries be set aside in case of pivoting, which
c  means that the condition becomes MU+ML .lt. (J-I).
c
      if ( i .le. n1 .and. j .le. n1 ) then
        if ( mu+ml .lt. (j-i) .or. ml .lt. (i-j) ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8BB_SET - Warning!'
          write ( *, '(a,i8,a,i8,a)' ) 
     &      '  Unable to set entry (', i, ',', j, ').'
          return
        else
          ij = (i-j+ml+mu+1)+(j-1)*(2*ml+mu+1)
        end if
c
c  The A2 block of the matrix.
c
      else if ( i .le. n1 .and. n1 .lt. j ) then
        ij = (2*ml+mu+1)*n1+(j-n1-1)*n1+i
c
c  The A3 and A4 blocks of the matrix.
c
      else if ( n1 .lt. i ) then
        ij = (2*ml+mu+1)*n1+n2*n1+(j-1)*n2+(i-n1)
      end if

      a(ij) = value

      return
      end
      subroutine r8bb_sl ( n1, n2, ml, mu, a_lu, pivot, b )

c*********************************************************************72
c
cc R8BB_SL solves an R8BB system factored by R8BB_FA.
c
c  Discussion:
c
c    The R8BB storage format is for a border banded matrix.  Such a
c    matrix has the logical form:
c
c      A1 | A2
c      ---+---
c      A3 | A4
c
c    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
c    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
c    respectively.
c
c    A should be defined as a vector.  The user must then store
c    the entries of the four blocks of the matrix into the vector A.
c    Each block is stored by columns.
c
c    A1, the banded portion of the matrix, is stored in
c    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
c    general band format.  The reason for the factor of 2 in front of
c    ML is to allocate space that may be required if pivoting occurs.
c
c    The following formulas should be used to determine how to store
c    the entry corresponding to row I and column J in the original matrix:
c
c    Entries of A1:
c
c      1 .le. I .le. N1, 1 .le. J .le. N1, (J-I) .le. MU and (I-J) .le. ML.
c
c      Store the I, J entry into location
c      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
c
c    Entries of A2:
c
c      1 .le. I .le. N1, N1+1 .le. J .le. N1+N2.
c
c      Store the I, J entry into location
c      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
c
c    Entries of A3:
c
c      N1+1 .le. I .le. N1+N2, 1 .le. J .le. N1.
c
c      Store the I, J entry into location
c      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c
c    Entries of A4:
c
c      N1+1 .le. I .le. N1+N2, N1+1 .le. J .le. N1+N2
c
c      Store the I, J entry into location
c      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c      (same formula used for A3).
c
c    The linear system A * x = b is decomposable into the block system:
c
c      ( A1 A2 ) * (X1) = (B1)
c      ( A3 A4 )   (X2)   (B2)
c
c    All the arguments except B are input quantities only, which are
c    not changed by the routine.  They should have exactly the same values
c    they had on exit from R8BB_FA.
c
c    If more than one right hand side is to be solved, with the same matrix,
c    R8BB_SL should be called repeatedly.  However, R8BB_FA only needs to be
c    called once to create the factorization.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, the order of the banded and dense 
c    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative and no greater than N1-1.
c
c    Input, double precision A_LU( (2*ML+MU+1)*N1 + 2*N1*N2 + N2*N2), 
c    the LU factors from R8BB_FA.
c
c    Input, integer PIVOT(N1+N2), the pivoting information 
c    from R8BB_FA.
c
c    Input/output, double precision B(N1+N2).
c    On input, B contains the right hand side of the linear system.
c    On output, B contains the solution.
c
      implicit none

      integer ml
      integer mu
      integer n1
      integer n2

      double precision a_lu((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
      double precision b(n1+n2)
      integer i
      integer ij
      integer pivot(n1+n2)
      integer j
      integer job
      integer nband

      nband = (2*ml+mu+1)*n1
c
c  Set B1 := inverse(A1) * B1.
c
      if ( 0 .lt. n1 ) then
        job = 0
        call r8gb_sl ( n1, ml, mu, a_lu, pivot, b, job )
      end if
c
c  Modify the right hand side of the second linear subsystem.
c  Set B2 := B2 - A3*B1.
c
      do i = 1, n2
        do j = 1, n1
          ij = nband + n1*n2 + (j-1)*n2 + i
          b(n1+i) = b(n1+i) - a_lu(ij) * b(j)
        end do
      end do
c
c  Set B2 := inverse(A4) * B2.
c
      if ( 0 .lt. n2 ) then
        job = 0
        call r8ge_sl ( n2, a_lu(nband+2*n1*n2+1), pivot(n1+1), 
     &    b(n1+1), job )
      end if
c
c  Modify the first subsolution.
c  Set B1 := B1 + A2*B2.
c
      do i = 1, n1
        do j = 1, n2
          ij = nband + (j-1)*n1 + i
          b(i) = b(i) + a_lu(ij) * b(n1+j)
        end do
      end do

      return
      end
      subroutine r8bb_to_r8ge ( n1, n2, ml, mu, a, b )

c*********************************************************************72
c
cc R8BB_TO_R8GE copies an R8BB matrix to an R8GE matrix.
c
c  Discussion:
c
c    The R8BB storage format is for a border banded matrix.  Such a
c    matrix has the logical form:
c
c      A1 | A2
c      ---+---
c      A3 | A4
c
c    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
c    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
c    respectively.
c
c    A should be defined as a vector.  The user must then store
c    the entries of the four blocks of the matrix into the vector A.
c    Each block is stored by columns.
c
c    A1, the banded portion of the matrix, is stored in
c    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
c    general band format.  The reason for the factor of 2 in front of
c    ML is to allocate space that may be required if pivoting occurs.
c
c    The following formulas should be used to determine how to store
c    the entry corresponding to row I and column J in the original matrix:
c
c    Entries of A1:
c
c      1 .le. I .le. N1, 1 .le. J .le. N1, (J-I) .le. MU and (I-J) .le. ML.
c
c      Store the I, J entry into location
c      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
c
c    Entries of A2:
c
c      1 .le. I .le. N1, N1+1 .le. J .le. N1+N2.
c
c      Store the I, J entry into location
c      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
c
c    Entries of A3:
c
c      N1+1 .le. I .le. N1+N2, 1 .le. J .le. N1.
c
c      Store the I, J entry into location
c      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c
c    Entries of A4:
c
c      N1+1 .le. I .le. N1+N2, N1+1 .le. J .le. N1+N2
c
c      Store the I, J entry into location
c      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c      (same formula used for A3).
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, the order of the banded and dense blocks.
c    N1 and N2 must be nonnegative, and at least one must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than N1-1.
c
c    Input, double precision A((2*ML+MU+1)*N1+2*N1*N2+N2*N2), the R8BB matrix.
c
c    Output, double precision B(N1+N2,N1+N2), the R8GE matrix.
c
      implicit none

      integer ml
      integer mu
      integer n1
      integer n2

      double precision a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
      double precision b(n1+n2,n1+n2)
      integer i
      integer ij
      integer j

      do i = 1, n1
        do j = 1, n1

          if ( mu+ml .lt. (j-i) .or. ml .lt. (i-j) ) then
            b(i,j) = 0.0D+00
          else
            ij = (i-j+ml+mu+1)+(j-1)*(2*ml+mu+1)
            b(i,j) = a(ij)
          end if

        end do
      end do

      do i = 1, n1
        do j = n1+1, n2
          ij = (2*ml+mu+1)*n1+(j-n1-1)*n1+i
          b(i,j) = a(ij)
        end do
      end do

      do i = n1+1, n2
        do j = 1, n1+n2
          ij = (2*ml+mu+1)*n1+n2*n1+(j-1)*n2+(i-n1)
          b(i,j) = a(ij)
        end do
      end do

      return
      end
      subroutine r8bb_vxm ( n1, n2, ml, mu, a, x, b )

c*********************************************************************72
c
cc R8BB_VXM multiplies an R8VEC by an R8BB matrix.
c
c  Discussion:
c
c    The R8BB storage format is for a border banded matrix.  Such a
c    matrix has the logical form:
c
c      A1 | A2
c      ---+---
c      A3 | A4
c
c    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
c    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
c    respectively.
c
c    A should be defined as a vector.  The user must then store
c    the entries of the four blocks of the matrix into the vector A.
c    Each block is stored by columns.
c
c    A1, the banded portion of the matrix, is stored in
c    the first (2*ML+MU+1)*N1 entries of A, using standard LINPACK
c    general band format.  The reason for the factor of 2 in front of
c    ML is to allocate space that may be required if pivoting occurs.
c
c    The following formulas should be used to determine how to store
c    the entry corresponding to row I and column J in the original matrix:
c
c    Entries of A1:
c
c      1 .le. I .le. N1, 1 .le. J .le. N1, (J-I) .le. MU and (I-J) .le. ML.
c
c      Store the I, J entry into location
c      (I-J+ML+MU+1)+(J-1)*(2*ML+MU+1).
c
c    Entries of A2:
c
c      1 .le. I .le. N1, N1+1 .le. J .le. N1+N2.
c
c      Store the I, J entry into location
c      (2*ML+MU+1)*N1+(J-N1-1)*N1+I.
c
c    Entries of A3:
c
c      N1+1 .le. I .le. N1+N2, 1 .le. J .le. N1.
c
c      Store the I, J entry into location
c      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c
c    Entries of A4:
c
c      N1+1 .le. I .le. N1+N2, N1+1 .le. J .le. N1+N2
c
c      Store the I, J entry into location
c      (2*ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c      (same formula used for A3).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, the order of the banded and dense blocks
c    N1 and N2 must be nonnegative, and at least one must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative and no greater than N1-1.
c
c    Input, double precision A((2*ML+MU+1)*N1 + 2*N1*N2 + N2*N2),
c    the R8BB matrix.
c
c    Input, double precision X(N1+N2), the vector to multiply A.
c
c    Output, double precision B(N1+N2), the product X times A.
c
      implicit none

      integer ml
      integer mu
      integer n1
      integer n2

      double precision a((2*ml+mu+1)*n1+2*n1*n2+n2*n2)
      double precision b(n1+n2)
      integer i
      integer ihi
      integer ij
      integer ilo
      integer j
      double precision x(n1+n2)
c
c  Initialize B.
c
      b(1:n1+n2) = 0.0D+00
c
c  Multiply by A1.
c
      do j = 1, n1
        ilo = max ( 1, j - mu - ml )
        ihi = min ( n1, j + ml )
        ij = (j-1) * (2*ml+mu+1) - j + ml + mu + 1
        do i = ilo, ihi
          b(j) = b(j) + x(i) * a(ij+i)
        end do
      end do
c
c  Multiply by A2.
c
      do j = n1+1, n1+n2
        ij = (2*ml+mu+1)*n1+(j-n1-1)*n1
        do i = 1, n1
          b(j) = b(j) + x(i) * a(ij+i)
        end do
      end do
c
c  Multiply by A3 and A4.
c
      do j = 1, n1+n2
        ij = (2*ml+mu+1)*n1+n1*n2+(j-1)*n2-n1
        do i = n1+1, n1+n2
          b(j) = b(j) + x(i) * a(ij+i)
        end do
      end do

      return
      end
      subroutine r8blt_det ( n, ml, a, det )

c*********************************************************************72
c
cc R8BLT_DET computes the determinant of an R8BLT matrix.
c
c  Discussion:
c
c    The R8BLT storage format is used for a banded lower triangular matrix.
c    The matrix is assumed to be zero below the ML-th subdiagonal.
c    The matrix is stored in an ML+1 by N array, in which the diagonal
c    appears in the first row, followed by successive subdiagonals.
c    Columns are preserved.
c
c  Example:
c
c    N = 5, ML = 2
c
c    A11   0   0   0   0
c    A21 A22   0   0   0
c    A31 A32 A33   0   0
c      0 A42 A43 A44   0
c      0   0 A53 A54 A55
c                --- ---
c                    ---
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, integer ML, the lower bandwidth.
c
c    Input, double precision A(ML+1,N), the R8BLT matrix.
c
c    Output, double precision DET, the determinant of A.
c
      implicit none

      integer ml
      integer n

      double precision a(ml+1,n)
      double precision det
      double precision r8vec_product

      det = r8vec_product ( n, a(1,1:n) )

      return
      end
      subroutine r8blt_indicator ( n, ml, a )

c*********************************************************************72
c
cc R8BLT_INDICATOR sets up an R8BLT indicator matrix.
c
c  Discussion:
c
c    The "indicator matrix" simply has a value like I*10+J at every
c    entry of a dense matrix, or at every entry of a compressed storage
c    matrix for which storage is allocated. 
c
c    The R8BLT storage format is used for a banded lower triangular matrix.
c    The matrix is assumed to be zero below the ML-th subdiagonal.
c    The matrix is stored in an ML+1 by N array, in which the diagonal
c    appears in the first row, followed by successive subdiagonals.
c    Columns are preserved.
c
c  Example:
c
c    N = 5, ML = 2
c
c    A11   0   0   0   0
c    A21 A22   0   0   0
c    A31 A32 A33   0   0
c      0 A42 A43 A44   0
c      0   0 A53 A54 A55
c                --- ---
c                    ---
c
c    The indicator matrix is stored as:
c
c      11  22  33  44  55
c      21  32  43  54   0
c      31  42  53   0   0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of columns of the matrix.
c
c    Input, integer ML, the lower bandwidth.
c
c    Output, double precision A(ML+1,N), the R8BLT matrix.
c
      implicit none

      integer ml
      integer n

      double precision a(ml+1,n)
      integer fac
      integer i
      integer i4_log_10
      integer j

      fac = 10 ** ( i4_log_10 ( n ) + 1 )

      do i = 1, n
        do j = max ( 1, i - ml ), i
          a(i-j+1,j) = dble ( fac * i + j )
        end do
      end do

      do i = n+1, n+ml
        do j = i-ml, n
          a(i-j+1,j) = 0.0D+00
        end do
      end do

      return
      end
      subroutine r8blt_mxv ( n, ml, a, x, b )

c*********************************************************************72
c
cc R8BLT_MXV multiplies an R8BLT matrix by an R8VEC.
c
c  Discussion:
c
c    The R8BLT storage format is used for a banded lower triangular matrix.
c    The matrix is assumed to be zero below the ML-th subdiagonal.
c    The matrix is stored in an ML+1 by N array, in which the diagonal
c    appears in the first row, followed by successive subdiagonals.
c    Columns are preserved.
c
c  Example:
c
c    N = 5, ML = 2
c
c    A11   0   0   0   0
c    A21 A22   0   0   0
c    A31 A32 A33   0   0
c      0 A42 A43 A44   0
c      0   0 A53 A54 A55
c                --- ---
c                    ---
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, integer ML, the lower bandwidth.
c
c    Input, double precision A(ML+1,N), the R8BLT matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(N), the product A * x.
c
      implicit none
c
      integer ml
      integer n

      double precision a(ml+1,n)
      double precision b(n)
      integer i
      integer j
      integer jhi
      integer jlo
      double precision x(n)

      do i = 1, n
        b(i) = 0.0D+00
        jlo = max ( 1, i - ml )
        jhi = i
        do j = jlo, jhi
          b(i) = b(i) + a(i-j+1,j) * x(j)
        end do
      end do

      return
      end
      subroutine r8blt_print ( n, ml, a, title )

c*********************************************************************72
c
cc R8BLT_PRINT prints an R8BLT matrix.
c
c  Discussion:
c
c    The R8BLT storage format is used for a banded lower triangular matrix.
c    The matrix is assumed to be zero below the ML-th subdiagonal.
c    The matrix is stored in an ML+1 by N array, in which the diagonal
c    appears in the first row, followed by successive subdiagonals.
c    Columns are preserved.
c
c  Example:
c
c    N = 5, ML = 2
c
c    A11   0   0   0   0
c    A21 A22   0   0   0
c    A31 A32 A33   0   0
c      0 A42 A43 A44   0
c      0   0 A53 A54 A55
c                --- ---
c                    ---
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, integer ML, the lower bandwidth.
c
c    Input, double precision A(ML+1,N), the R8BLT matrix.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer ml
      integer n

      double precision a(ml+1,n)
      character * ( * )  title

      call r8blt_print_some ( n, ml, a, 1, 1, n, n, title )

      return
      end
      subroutine r8blt_print_some ( n, ml, a, ilo, jlo, ihi, jhi, 
     &  title )

c*********************************************************************72
c
cc R8BLT_PRINT_SOME prints some of an R8BLT matrix.
c
c  Discussion:
c
c    The R8BLT storage format is used for a banded lower triangular matrix.
c    The matrix is assumed to be zero below the ML-th subdiagonal.
c    The matrix is stored in an ML+1 by N array, in which the diagonal
c    appears in the first row, followed by successive subdiagonals.
c    Columns are preserved.
c
c  Example:
c
c    N = 5, ML = 2
c
c    A11   0   0   0   0
c    A21 A22   0   0   0
c    A31 A32 A33   0   0
c      0 A42 A43 A44   0
c      0   0 A53 A54 A55
c                --- ---
c                    ---
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, integer ML, the lower bandwidth.
c
c    Input, double precision A(ML+1,N), the R8BLT matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 5 )
      integer ml
      integer n

      double precision a(ml+1,n)
      double precision aij
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
      character ( len = * )  title

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
        i2lo = max ( i2lo, j2lo )
        i2hi = min ( ihi, n )
        i2hi = min ( i2hi, j2hi + ml )

        do i = i2lo, i2hi
c
c  Print out (up to) 5 entries in row I, that lie in the current strip.
c
          do j2 = 1, inc

            j = j2lo - 1 + j2

            if ( j .le. i .and. i .le. j + ml ) then
              aij = a(i-j+1,j)
            else
              aij = 0.0D+00
            end if

            if ( ml .lt. i-j .or. 0 .lt. j-i ) then
              ctemp(j2) = '              '
            else if ( r8_is_int ( aij ) ) then
              write ( ctemp(j2), '(f8.0,6x)' ) aij
            else
              write ( ctemp(j2), '(g14.6)' ) aij
            end if

          end do

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine r8blt_random ( n, ml, seed, a )

c*********************************************************************72
c
cc R8BLT_RANDOM randomizes an R8BLT matrix.
c
c  Discussion:
c
c    The R8BLT storage format is used for a banded lower triangular matrix.
c    The matrix is assumed to be zero below the ML-th subdiagonal.
c    The matrix is stored in an ML+1 by N array, in which the diagonal
c    appears in the first row, followed by successive subdiagonals.
c    Columns are preserved.
c
c  Example:
c
c    N = 5, ML = 2
c
c    A11   0   0   0   0
c    A21 A22   0   0   0
c    A31 A32 A33   0   0
c      0 A42 A43 A44   0
c      0   0 A53 A54 A55
c                --- ---
c                    ---
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of columns of the matrix.
c
c    Input, integer ML, the lower bandwidth.
c
c    Input/output, integer SEED, a seed for the random 
c    number generator.
c
c    Output, double precision A(ML+1,N), the R8BLT matrix.
c
      implicit none

      integer ml
      integer n

      double precision a(ml+1,n)
      double precision r8_uniform_01
      integer i
      integer j
      integer seed

      do i = 1, n
        do j = max ( 1, i - ml ), i
          a(i-j+1,j) = r8_uniform_01 ( seed )
        end do
      end do
c
c  The junk entries can be thought of as corresponding to
c  elements of a phantom portion of the matrix.
c
      do i = n+1, n+ml
        do j = i-ml, n
          a(i-j+1,j) = 0.0D+00
        end do
      end do

      return
      end
      subroutine r8blt_sl ( n, ml, a, b, job )

c*********************************************************************72
c
cc R8BLT_SL solves an R8BLT system.
c
c  Discussion:
c
c    The R8BLT storage format is used for a banded lower triangular matrix.
c    The matrix is assumed to be zero below the ML-th subdiagonal.
c    The matrix is stored in an ML+1 by N array, in which the diagonal
c    appears in the first row, followed by successive subdiagonals.
c    Columns are preserved.
c
c    No factorization of the lower triangular matrix is required.
c
c  Example:
c
c    N = 5, ML = 2
c
c    A11   0   0   0   0
c    A21 A22   0   0   0
c    A31 A32 A33   0   0
c      0 A42 A43 A44   0
c      0   0 A53 A54 A55
c                --- ---
c                    ---
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, integer ML, the lower bandwidth.
c
c    Input, double precision A(ML+1,N), the R8BLT matrix.
c
c    Input/output, double precision B(N).
c    On input, the right hand side.
c    On output, the solution vector.
c
c    Input, integer JOB, is 0 to solve the untransposed system,
c    nonzero to solve the transposed system.
c
      implicit none

      integer ml
      integer n

      double precision a(ml+1,n)
      double precision b(n)
      integer i
      integer ihi
      integer ilo
      integer j
      integer job

      if ( job .eq. 0 ) then

        do j = 1, n
          b(j) = b(j) / a(1,j)
          ihi = min ( j + ml, n )
          do i = j+1, ihi
            b(i) = b(i) - a(i-j+1,j) * b(j)
          end do
        end do

      else

        do j = n, 1, -1
          b(j) = b(j) / a(1,j)
          ilo = max ( j - ml, 1 )
          do i = ilo, j-1
            b(i) = b(i) - a(j-i+1,i) * b(j)
          end do
        end do

      end if

      return
      end
      subroutine r8blt_to_r8ge ( n, ml, a, b )

c*********************************************************************72
c
cc R8BLT_TO_R8GE copies an R8BLT matrix to an R8GE matrix.
c
c  Discussion:
c
c    The R8BLT storage format is used for a banded lower triangular matrix.
c    The matrix is assumed to be zero below the ML-th subdiagonal.
c    The matrix is stored in an ML+1 by N array, in which the diagonal
c    appears in the first row, followed by successive subdiagonals.
c    Columns are preserved.
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c  Example:
c
c    N = 5, ML = 2
c
c    A11   0   0   0   0
c    A21 A22   0   0   0
c    A31 A32 A33   0   0
c      0 A42 A43 A44   0
c      0   0 A53 A54 A55
c                --- ---
c                    ---
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrices.
c    N must be positive.
c
c    Input, integer ML, the lower bandwidth of A.
c    ML must be nonnegative, and no greater than N-1.
c
c    Input, double precision A(ML+1,N), the R8BLT matrix.
c
c    Output, double precision B(N,N), the R8GE matrix.
c
      implicit none

      integer ml
      integer n

      double precision a(ml+1,n)
      double precision b(n,n)
      integer i
      integer j

      do i = 1, n
        do j = 1, n
          if ( j .le. i .and. i .le. j + ml ) then
            b(i,j) = a(i-j+1,j)
          else
            b(i,j) = 0.0D+00
          end if
        end do
      end do

      return
      end
      subroutine r8blt_vxm ( n, ml, a, x, b )

c*********************************************************************72
c
cc R8BLT_VXM multiplies an R8VEC by an R8BLT matrix.
c
c  Discussion:
c
c    The R8BLT storage format is used for a banded lower triangular matrix.
c    The matrix is assumed to be zero below the ML-th subdiagonal.
c    The matrix is stored in an ML+1 by N array, in which the diagonal
c    appears in the first row, followed by successive subdiagonals.
c    Columns are preserved.
c
c  Example:
c
c    N = 5, ML = 2
c
c    A11   0   0   0   0
c    A21 A22   0   0   0
c    A31 A32 A33   0   0
c      0 A42 A43 A44   0
c      0   0 A53 A54 A55
c                --- ---
c                    ---
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, integer ML, the lower bandwidth.
c
c    Input, double precision A(ML+1,N), the R8BLT matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(N), the product X*A.
c
      implicit none

      integer ml
      integer n

      double precision a(ml+1,n)
      double precision b(n)
      integer i
      integer j
      integer jhi
      integer jlo
      double precision x(n)

      b(1:n) = 0.0D+00

      do i = 1, n
        jlo = max ( 1, i - ml )
        jhi = i
        do j = jlo, jhi
          b(j) = b(j) + x(i) * a(i-j+1,j)
        end do
      end do

      return
      end
      subroutine r8bto_indicator ( m, l, a )

c*********************************************************************72
c
cc R8BTO_INDICATOR sets up an R8BTO indicator matrix.
c
c  Discussion:
c
c    The "indicator matrix" simply has a value like I*10+J at every
c    entry of a dense matrix, or at every entry of a compressed storage
c    matrix for which storage is allocated. 
c
c    The R8BTO storage format is for a block Toeplitz matrix. The matrix
c    can be regarded as an L by L array of blocks, each of size M by M.
c    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
c    that is, along its diagonal, the blocks repeat.
c
c    Storage for the matrix consists of the L blocks of the first row,
c    followed by the L-1 blocks of the first column (skipping the first row).
c    These items are stored in the natural way in an (M,M,2*L-1) array.
c
c  Example:
c
c    M = 2, L = 3
c
c    1 2 | 3 4 | 5 6
c    5 5 | 6 6 | 7 7
c    ----+-----+-----
c    7 8 | 1 2 | 3 4
c    8 8 | 5 5 | 6 6
c    ----+-----+-----
c    9 0 | 7 8 | 1 2
c    9 9 | 8 8 | 5 5
c
c    X = (/ 1, 2, 3, 4, 5, 6 /)
c
c    B = (/ 91, 134, 73, 125, 97, 129 /)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the order of the blocks of the matrix A.
c
c    Input, integer L, the number of blocks in a row or column of A.
c
c    Output, double precision A(M,M,2*L-1), the R8BTO matrix.
c
      implicit none

      integer l
      integer m

      double precision a(m,m,2*l-1)
      integer fac
      integer i
      integer i4_log_10
      integer i2
      integer j
      integer j2
      integer k

      fac = 10 ** ( i4_log_10 ( m * l ) + 1 )
c
c  Blocks 1 to L form the first row.
c
      j = 0

      do k = 1, l

        do j2 = 1, m
          j = j + 1
          do i = 1, m
            a(i,j2,k) = dble ( fac * i + j )
          end do
        end do
      end do
c
c  Blocks L+1 through 2*L-1 form the remainder of the first column.
c
      i = m

      do k = l+1, 2*l-1

        do i2 = 1, m
          i = i + 1
          do j = 1, m
            a(i2,j,k) = dble ( fac * i + j )
          end do
        end do

      end do

      return
      end
      subroutine r8bto_mxv ( m, l, a, x, b )

c*********************************************************************72
c
cc R8BTO_MXV multiplies an R8BTO matrix by an R8VEC.
c
c  Discussion:
c
c    The R8BTO storage format is for a block Toeplitz matrix. The matrix
c    can be regarded as an L by L array of blocks, each of size M by M.
c    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
c    that is, along its diagonal, the blocks repeat.
c
c    Storage for the matrix consists of the L blocks of the first row,
c    followed by the L-1 blocks of the first column (skipping the first row).
c    These items are stored in the natural way in an (M,M,2*L-1) array.
c
c  Example:
c
c    M = 2, L = 3
c
c    1 2 | 3 4 | 5 6
c    5 5 | 6 6 | 7 7
c    ----+-----+-----
c    7 8 | 1 2 | 3 4
c    8 8 | 5 5 | 6 6
c    ----+-----+-----
c    9 0 | 7 8 | 1 2
c    9 9 | 8 8 | 5 5
c
c    X = (/ 1, 2, 3, 4, 5, 6 /)
c
c    B = (/ 91, 134, 73, 125, 79, 138 /)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the order of the blocks of the matrix A.
c
c    Input, integer L, the number of blocks in a row or column of A.
c
c    Input, double precision A(M,M,2*L-1), the R8BTO matrix.
c
c    Input, double precision X(M*L), the vector to be multiplied.
c
c    Output, double precision B(M*L), the product A * X.
c
      implicit none

      integer l
      integer m

      double precision a(m,m,2*l-1)
      double precision b(m,l)
      integer i
      integer ii
      integer j
      integer k
      double precision x(m,l)
c
c  Construct the right hand side by blocks.
c
      do i = 1, l

        do ii = 1, m
          b(ii,i) = 0.0D+00
        end do

        do j = 1, i-1
          do ii = 1, m
            do k = 1, m
              b(ii,i) = b(ii,i) + a(ii,k,l+i-j) * x(k,j)
            end do
          end do
        end do

        do j = i, l
          do ii = 1, m
            do k = 1, m
              b(ii,i) = b(ii,i) + a(ii,k,j+1-i) * x(k,j)
            end do
          end do
        end do

      end do

      return
      end
      subroutine r8bto_print ( m, l, a, title )

c*********************************************************************72
c
cc R8BTO_PRINT prints an R8BTO matrix.
c
c  Discussion:
c
c    The R8BTO storage format is for a block Toeplitz matrix. The matrix
c    can be regarded as an L by L array of blocks, each of size M by M.
c    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
c    that is, along its diagonal, the blocks repeat.
c
c    Storage for the matrix consists of the L blocks of the first row,
c    followed by the L-1 blocks of the first column (skipping the first row).
c    These items are stored in the natural way in an (M,M,2*L-1) array.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the order of the blocks of the matrix A.
c
c    Input, integer L, the number of blocks in a row or column of A.
c
c    Input, double precision A(M,M,2*L-1), the R8BTO matrix.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer l
      integer m

      double precision a(m,m,2*l-1)
      character * ( * )  title

      call r8bto_print_some ( m, l, a, 1, 1, m*l, m*l, title )

      return
      end
      subroutine r8bto_print_some ( m, l, a, ilo, jlo, ihi, jhi, title )

c*********************************************************************72
c
cc R8BTO_PRINT_SOME prints some of an R8BTO matrix.
c
c  Discussion:
c
c    The R8BTO storage format is for a block Toeplitz matrix. The matrix
c    can be regarded as an L by L array of blocks, each of size M by M.
c    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
c    that is, along its diagonal, the blocks repeat.
c
c    Storage for the matrix consists of the L blocks of the first row,
c    followed by the L-1 blocks of the first column (skipping the first row).
c    These items are stored in the natural way in an (M,M,2*L-1) array.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the order of the blocks of the matrix A.
c
c    Input, integer L, the number of blocks in a row or column of A.
c
c    Input, double precision A(M,M,2*L-1), the R8BTO matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 5 )
      integer l
      integer m

      double precision a(m,m,2*l-1)
      double precision aij
      character * ( 14 ) ctemp(incx)
      logical r8_is_int
      integer i
      integer i1
      integer i2
      integer i3hi
      integer i3lo
      integer ihi
      integer ilo
      integer inc
      integer j
      integer j1
      integer j2
      integer j3
      integer j3hi
      integer j3lo
      integer jhi
      integer jlo
      integer n
      character * ( * )  title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )

      n = m * l
c
c  Print the columns of the matrix, in strips of 5.
c
      do j3lo = jlo, jhi, incx

        j3hi = j3lo + incx - 1
        j3hi = min ( j3hi, n )
        j3hi = min ( j3hi, jhi )

        inc = j3hi + 1 - j3lo

        write ( *, '(a)' ) ' '

        do j = j3lo, j3hi
          j3 = j + 1 - j3lo
          write ( ctemp(j3), '(i7,7x)' ) j
        end do

        write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j3), j3 = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) '  ---'
c
c  Determine the range of the rows in this strip.
c
        i3lo = max ( ilo, 1 )
        i3hi = min ( ihi, n )

        do i = i3lo, i3hi
c
c  Print out (up to) 5 entries in row I, that lie in the current strip.
c
          do j3 = 1, inc

            j = j3lo - 1 + j3
c
c  i = M * ( i1 - 1 ) + i2
c  j = M * ( j1 - 1 ) + j2
c
            i1 = ( i - 1 ) / m + 1
            i2 = i - m * ( i1 - 1 )
            j1 = ( j - 1 ) / m + 1
            j2 = j - m * ( j1 - 1 )

            if ( i1 .le. j1 ) then
              aij = a(i2,j2,j1+1-i1)
            else
              aij = a(i2,j2,l+i1-j1)
            end if

            if ( r8_is_int ( aij ) ) then
              write ( ctemp(j3), '(f8.0,6x)' ) aij
            else
              write ( ctemp(j3), '(g14.6)' ) aij
            end if

          end do

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j3), j3 = 1, inc )

        end do

      end do

      return
      end
      subroutine r8bto_random ( m, l, seed, a )

c*********************************************************************72
c
cc R8BTO_RANDOM randomizes an R8BTO matrix.
c
c  Discussion:
c
c    The R8BTO storage format is for a block Toeplitz matrix. The matrix
c    can be regarded as an L by L array of blocks, each of size M by M.
c    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
c    that is, along its diagonal, the blocks repeat.
c
c    Storage for the matrix consists of the L blocks of the first row,
c    followed by the L-1 blocks of the first column (skipping the first row).
c    These items are stored in the natural way in an (M,M,2*L-1) array.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the order of the blocks of the matrix A.
c
c    Input, integer L, the number of blocks in a row or column of A.
c
c    Input/output, integer SEED, a seed for the random number 
c    generator.
c
c    Output, double precision A1(M,M,2*L-1), the R8BTO matrix.
c
      implicit none

      integer m
      integer l

      double precision a(m,m,2*l-1)
      double precision r8_uniform_01
      integer i
      integer j
      integer k
      integer seed

      do i = 1, m
        do j = 1, m
          do k = 1, 2 * l - 1 
            a(i,j,k) = r8_uniform_01 ( seed )
          end do
        end do
      end do

      return
      end
      subroutine r8bto_sl ( m, l, a, b, x )

c*********************************************************************72
c
cc R8BTO_SL solves an R8BTO system.
c
c  Discussion:
c
c    The R8BTO storage format is for a block Toeplitz matrix. The matrix
c    can be regarded as an L by L array of blocks, each of size M by M.
c    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
c    that is, along its diagonal, the blocks repeat.
c
c    Storage for the matrix consists of the L blocks of the first row,
c    followed by the L-1 blocks of the first column (skipping the first row).
c    These items are stored in the natural way in an (M,M,2*L-1) array.
c
c  Example:
c
c    M = 2, L = 3
c
c    1 2 | 3 4 | 5 6
c    5 5 | 6 6 | 7 7
c    ----+-----+-----
c    7 8 | 1 2 | 3 4
c    8 8 | 5 5 | 6 6
c    ----+-----+-----
c    9 0 | 7 8 | 1 2
c    9 9 | 8 8 | 5 5
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 June 2009
c
c  Author:
c
c    FORTRAN77 version by John Burkardt
c
c  Parameters:
c
c    Input, integer M, the order of the blocks of the matrix A.
c
c    Input, integer L, the number of blocks in a row or column
c    of A.
c
c    Input, double precision A(M*M,2*L-1), the R8BTO matrix.
c
c    Input, double precision B(M*L), the right hand side vector.
c
c    Output, double precision X(M*L), the solution vector.  X and B
c    may share storage.
c
      implicit none

      integer l
      integer m

      double precision a(m*m,2*l-1)
      double precision b(m,l)
      double precision c1(m,m,l-1)
      double precision c2(m,m,l-1)
      integer i
      integer i1
      integer i2
      integer i3
      integer info
      integer j
      integer job
      integer k
      integer n
      integer pivot(m)
      double precision r1(m,m)
      double precision r2(m,m)
      double precision r3(m,m)
      double precision r5(m,m)
      double precision r6(m,m)
      double precision x(m,l)
c
c  Solve the system with the principal minor of order M.
c
      i3 = 1
      do j = 1, m
        do i = 1, m
          c1(i,j,1) = a(i3,1)
          r1(i,j) = a(i3,1)
          i3 = i3 + 1
        end do
      end do

      do j = 1, m
        do i = 1, m
          r3(i,j) = r1(i,j)
        end do
      end do

      do i = 1, m
        x(i,1) = b(i,1)
      end do

      call r8ge_fa ( m, r3, pivot, info )

      job = 0
      call r8ge_sl ( m, r3, pivot, x(1,1), job )

      if ( l .eq. 1 ) then
        return
      end if
c
c  Recurrent process for solving the system
c  with the block Toeplitz matrix for N = 2 through L.
c
      do n = 2, l
c
c  Compute multiples of the first and last block columns of
c  the inverse of the principal minor of order M*N.
c
        i3 = 1
        do j = 1, m
          do i = 1, m
            r5(i,j) = a(i3,l+n-1)
            r6(i,j) = a(i3,n)
            i3 = i3 + 1
          end do
        end do

        if ( 2 .lt. n ) then

          do j = 1, m
            do i = 1, m
              c1(i,j,n-1) = r2(i,j)
            end do
          end do

          do i1 = 1, n-2
            i2 = n - i1
            do j = 1, m
              i3 = 1
c
c  My apologies, but we have an I1+L, followed in the next line by I1+1.
c
              do i = 1, m
                call daxpy ( m, c1(i,j,i2), a(i3,i1+l), 1, r5(1,j), 1 )
                call daxpy ( m, c2(i,j,i1), a(i3,i1+1), 1, r6(1,j), 1 )
                i3 = i3 + m
              end do
            end do
          end do

        end if

        do j = 1, m
          do i = 1, m
            r2(i,j) = - r5(i,j)
          end do
          job = 0
          call r8ge_sl ( m, r3, pivot, r2(1,j), job )
        end do

        do j = 1, m
          do i = 1, m
            r3(i,j) = r6(i,j)
          end do
        end do

        do j = 1, m
          do i = 1, m
            r6(i,j) = -c1(i,j,1)
          end do
       end do

        do j = 1, m
          do i = 1, m
            do k = 1, m
              c1(k,j,1) = c1(k,j,1) + r2(i,j) * r3(k,i)
            end do
          end do
        end do

        call r8ge_fa ( m, r6, pivot, info )

        do j = 1, m
          call r8ge_sl ( m, r6, pivot, r3(1,j), job )
          do i = 1, m
            do k = 1, m
              r1(k,j) = r1(k,j) + r3(i,j) * r5(k,i)
            end do
          end do
        end do

        if ( 2 .lt. n ) then

          do j = 1, m
            do i = 1, m
              r6(i,j) = c2(i,j,1)
            end do
          end do

          do i1 = 2, n-1

            if ( i1 .ne. n-1 ) then
              do j = 1, m
                do i = 1, m
                  r5(i,j) = c2(i,j,i1)
                end do
              end do
            end if

            do j = 1, m
              do i = 1, m
                c2(i,j,i1) = r6(i,j)
              end do
              do i = 1, m
                call daxpy ( m, r3(i,j), c1(1,i,i1), 1, c2(1,j,i1), 1 )
              end do
            end do

            do j = 1, m
              do i = 1, m
                call daxpy ( m, r2(i,j), r6(1,i), 1, c1(1,j,i1), 1 )
              end do
            end do

            do j = 1, m
              do i = 1, m
                r6(i,j) = r5(i,j)
              end do
            end do

          end do

        end if

        do j = 1, m
          do i = 1, m
            c2(i,j,1) = r3(i,j)
          end do
        end do
c
c  Compute the solution of the system with the principal minor of order M*N.
c
        do j = 1, m
          do i = 1, m
            r3(i,j) = r1(i,j)
          end do
        end do

        do i = 1, m
          x(i,n) = b(i,n)
        end do

        do i1 = 1, n-1
          i2 = n - i1
          i3 = 1
          do i = 1, m
            call daxpy ( m, -x(i,i2), a(i3,i1+l), 1, x(1,n), 1 )
            i3 = i3 + m
          end do
        end do

        call r8ge_fa ( m, r3, pivot, info )

        call r8ge_sl ( m, r3, pivot, x(1,n), job )

        do i1 = 1, n-1
          do i = 1, m
            call daxpy ( m, x(i,n), c2(1,i,i1), 1, x(1,i1), 1 )
          end do
        end do

      end do

      return
      end
      subroutine r8bto_to_r8ge ( m, l, a, n, b )

c*********************************************************************72
c
cc R8BTO_TO_R8GE copies an R8BTO matrix to an R8GE matrix.
c
c  Discussion:
c
c    The R8BTO storage format is for a block Toeplitz matrix. The matrix
c    can be regarded as an L by L array of blocks, each of size M by M.
c    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
c    that is, along its diagonal, the blocks repeat.
c
c    Storage for the matrix consists of the L blocks of the first row,
c    followed by the L-1 blocks of the first column (skipping the first row).
c    These items are stored in the natural way in an (M,M,2*L-1) array.
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the order of the blocks of the R8BTO matrix.
c
c    Input, integer L, the number of blocks in a row or column
c    of the R8BTO matrix.
c
c    Input, double precision A(M,M,2*L-1), the R8BTO matrix.
c
c    Output, integer N, the order of the matrix, N = M*L.
c
c    Output, double precision B(N,N), the R8GE matrix.
c
      implicit none

      integer l
      integer m

      double precision a(m,m,2*l-1)
      double precision b(m*l,m*l)
      integer i
      integer i1
      integer i2
      integer j
      integer j1
      integer j2
      integer n

      n = m * l

      do i = 1, n

        i1 = ( i - 1 ) / m + 1
        i2 = i - m * ( i1 - 1 )

        do j = 1, n

          j1 = ( j - 1 ) / m + 1
          j2 = j - m * ( j1 - 1 )

          if ( i1 .le. j1 ) then
            b(i,j) = a(i2,j2,j1+1-i1)
          else
            b(i,j) = a(i2,j2,l+i1-j1)
          end if

        end do

      end do

      return
      end
      subroutine r8bto_vxm ( m, l, a, x, b )

c*********************************************************************72
c
cc R8BTO_VXM multiplies an R8VEC by an R8BTO matrix.
c
c  Discussion:
c
c    The R8BTO storage format is for a block Toeplitz matrix. The matrix
c    can be regarded as an L by L array of blocks, each of size M by M.
c    The full matrix has order N = M * L.  The L by L matrix is Toeplitz,
c    that is, along its diagonal, the blocks repeat.
c
c    Storage for the matrix consists of the L blocks of the first row,
c    followed by the L-1 blocks of the first column (skipping the first row).
c    These items are stored in the natural way in an (M,M,2*L-1) array.
c
c  Example:
c
c    M = 2, L = 3
c
c    1 2 | 3 4 | 5 6
c    5 5 | 6 6 | 7 7
c    ----+-----+-----
c    7 8 | 1 2 | 3 4
c    8 8 | 5 5 | 6 6
c    ----+-----+-----
c    9 0 | 7 8 | 1 2
c    9 9 | 8 8 | 5 5
c
c    X = (/ 1, 2, 3, 4, 5, 6 /)
c
c    B = (/ 163, 122, 121, 130, 87, 96 /)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the order of the blocks of the matrix A.
c
c    Input, integer L, the number of blocks in a row or 
c    column of A.
c
c    Input, double precision A(M,M,2*L-1), the R8BTO matrix.
c
c    Input, double precision X(M*L), the vector to be multiplied.
c
c    Output, double precision B(M*L), the product X * A.
c
      implicit none

      integer l
      integer m

      double precision a(m,m,2*l-1)
      double precision b(m,l)
      integer i
      integer ii
      integer j
      integer k
      double precision x(m,l)
c
c  Construct the right hand side by blocks.
c
      do i = 1, l

        do ii = 1, m
          b(ii,i) = 0.0D+00
        end do

        do j = 1, i
          do ii = 1, m
            do k = 1, m
              b(ii,i) = b(ii,i) + a(k,ii,i+1-j) * x(k,j)
            end do
          end do
        end do

        do j = i+1, l
          do ii = 1, m
            do k = 1, m
              b(ii,i) = b(ii,i) + a(k,ii,l+j-i) * x(k,j)
            end do
          end do
        end do

      end do

      return
      end
      subroutine r8but_det ( n, mu, a, det )

c*********************************************************************72
c
cc R8BUT_DET computes the determinant of an R8BUT matrix.
c
c  Discussion:
c
c    The R8BUT storage format is used for a banded upper triangular matrix.
c    The matrix is assumed to be zero above the MU-th superdiagonal.
c    The matrix is stored in an MU+1 by N array.
c    Columns are preserved.
c
c    The diagonal is stored in row MU+1 of the array.
c    The first superdiagonal in row MU, columns 2 through N.
c    The second superdiagonal in row MU-1, columns 3 through N.
c    The MU-th superdiagonal in row 1, columns MU+1 through N.
c
c  Example:
c
c    N = 5, MU = 2
c
c    A11 A12 A13   0   0
c      0 A22 A23 A24   0
c      0   0 A33 A34 A35
c      0   0   0 A44 A45
c      0   0   0   0 A55
c                --- ---
c                    ---
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, integer MU, the upper bandwidth.
c
c    Input, double precision A(MU+1,N), the R8BUT matrix.
c
c    Output, double precision DET, the determinant of A.
c
      implicit none

      integer mu
      integer n

      double precision a(mu+1,n)
      double precision det
      integer j

      det = 1.0D+00
      do j = 1, n
        det = det * a(mu+1,j)
      end do

      return
      end
      subroutine r8but_indicator ( n, mu, a )

c*********************************************************************72
c
cc R8BUT_INDICATOR sets up an R8BUT indicator matrix.
c
c  Discussion:
c
c    The "indicator matrix" simply has a value like I*10+J at every
c    entry of a dense matrix, or at every entry of a compressed storage
c    matrix for which storage is allocated. 
c
c    The R8BUT storage format is used for a banded upper triangular matrix.
c    The matrix is assumed to be zero above the MU-th superdiagonal.
c    The matrix is stored in an MU+1 by N array.
c    Columns are preserved.
c
c    The diagonal is stored in row MU+1 of the array.
c    The first superdiagonal in row MU, columns 2 through N.
c    The second superdiagonal in row MU-1, columns 3 through N.
c    The MU-th superdiagonal in row 1, columns MU+1 through N.
c
c  Example:
c
c    N = 5, MU = 2
c
c    A11 A12 A13   0   0
c      0 A22 A23 A24   0
c      0   0 A33 A34 A35
c      0   0   0 A44 A45
c      0   0   0   0 A55
c                --- ---
c                    ---
c
c    The indicator matrix is stored as:
c
c       0   0  13  24  35
c       0  12  23  34  45
c      11  22  33  44  55
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of columns of the matrix.
c
c    Input, integer MU, the upper bandwidth.
c
c    Output, double precision A(MU+1,N), the R8BUT matrix.
c
      implicit none

      integer mu
      integer n

      double precision a(mu+1,n)
      integer fac
      integer i
      integer i4_log_10
      integer j

      fac = 10 ** ( i4_log_10 ( n ) + 1 )

      do i = 1, n
        do j = i, min ( n, i + mu )
          a(i-j+mu+1,j) = dble ( fac * i + j )
        end do
      end do

      do i = 1, mu
        do j = 1, mu+1-i
          a(i,j) = 0.0D+00
        end do
      end do

      return
      end
      subroutine r8but_mxv ( n, mu, a, x, b )

c*********************************************************************72
c
cc R8BUT_MXV multiplies an R8BUT matrix by an R8VEC.
c
c  Discussion:
c
c    The R8BUT storage format is used for a banded upper triangular matrix.
c    The matrix is assumed to be zero above the MU-th superdiagonal.
c    The matrix is stored in an MU+1 by N array.
c    Columns are preserved.
c
c    The diagonal is stored in row MU+1 of the array.
c    The first superdiagonal in row MU, columns 2 through N.
c    The second superdiagonal in row MU-1, columns 3 through N.
c    The MU-th superdiagonal in row 1, columns MU+1 through N.
c
c  Example:
c
c    N = 5, MU = 2
c
c    A11 A12 A13   0   0
c      0 A22 A23 A24   0
c      0   0 A33 A34 A35
c      0   0   0 A44 A45
c      0   0   0   0 A55
c                --- ---
c                    ---
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, integer MU, the upper bandwidth.
c
c    Input, double precision A(MU+1,N), the R8BUT matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(N), the product A * x.
c
      implicit none

      integer mu
      integer n

      double precision a(mu+1,n)
      double precision b(n)
      integer i
      integer j
      double precision x(n)

      do i = 1, n
        b(i) = 0.0D+00
        do j = i, min ( n, i + mu )
          b(i) = b(i) + a(i-j+mu+1,j) * x(j)
        end do
      end do

      return
      end
      subroutine r8but_print ( n, mu, a, title )

c*********************************************************************72
c
cc R8BUT_PRINT prints an R8BUT matrix.
c
c  Discussion:
c
c    The R8BUT storage format is used for a banded upper triangular matrix.
c    The matrix is assumed to be zero above the MU-th superdiagonal.
c    The matrix is stored in an MU+1 by N array.
c    Columns are preserved.
c
c    The diagonal is stored in row MU+1 of the array.
c    The first superdiagonal in row MU, columns 2 through N.
c    The second superdiagonal in row MU-1, columns 3 through N.
c    The MU-th superdiagonal in row 1, columns MU+1 through N.
c
c  Example:
c
c    N = 5, MU = 2
c
c    A11 A12 A13   0   0
c      0 A22 A23 A24   0
c      0   0 A33 A34 A35
c      0   0   0 A44 A45
c      0   0   0   0 A55
c                --- ---
c                    ---
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, integer MU, the upper bandwidth.
c
c    Input, double precision A(MU+1,N), the R8BUT matrix.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer mu
      integer n

      double precision a(mu+1,n)
      character * ( * )  title

      call r8but_print_some ( n, mu, a, 1, 1, n, n, title )

      return
      end
      subroutine r8but_print_some ( n, mu, a, ilo, jlo, ihi, jhi, 
     &  title )

c*********************************************************************72
c
cc R8BUT_PRINT_SOME prints some of an R8BUT matrix.
c
c  Discussion:
c
c    The R8BUT storage format is used for a banded upper triangular matrix.
c    The matrix is assumed to be zero above the MU-th superdiagonal.
c    The matrix is stored in an MU+1 by N array.
c    Columns are preserved.
c
c    The diagonal is stored in row MU+1 of the array.
c    The first superdiagonal in row MU, columns 2 through N.
c    The second superdiagonal in row MU-1, columns 3 through N.
c    The MU-th superdiagonal in row 1, columns MU+1 through N.
c
c  Example:
c
c    N = 5, MU = 2
c
c    A11 A12 A13   0   0
c      0 A22 A23 A24   0
c      0   0 A33 A34 A35
c      0   0   0 A44 A45
c      0   0   0   0 A55
c                --- ---
c                    ---
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, integer MU, the upper bandwidth.
c
c    Input, double precision A(MU+1,N), the R8BUT matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer incx 
      parameter ( incx = 5 )
      integer mu
      integer n

      double precision a(mu+1,n)
      double precision aij
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
          write ( ctemp(j2), '(i7,7x)' ) j
        end do

        write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) '  ---'
c
c  Determine the range of the rows in this strip.
c
        i2lo = max ( ilo, 1 )
        i2lo = max ( i2lo, j2lo )
        i2hi = min ( ihi, n )
        i2hi = min ( i2hi, j2hi + mu )

        do i = i2lo, i2hi
c
c  Print out (up to) 5 entries in row I, that lie in the current strip.
c
          do j2 = 1, inc

            j = j2lo - 1 + j2

            if ( i .le. j .and. j .le. i + mu ) then
              aij = a(i-j+mu+1,j)
              if ( r8_is_int ( aij ) ) then
                write ( ctemp(j2), '(f8.0,6x)' ) aij
              else
                write ( ctemp(j2), '(g14.6)' ) aij
              end if
            else
              ctemp(j2) = '              '
            end if

          end do

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine r8but_random ( n, mu, seed, a )

c*********************************************************************72
c
cc R8BUT_RANDOM randomizes an R8BUT matrix.
c
c  Discussion:
c
c    The R8BUT storage format is used for a banded upper triangular matrix.
c    The matrix is assumed to be zero above the MU-th superdiagonal.
c    The matrix is stored in an MU+1 by N array.
c    Columns are preserved.
c
c    The diagonal is stored in row MU+1 of the array.
c    The first superdiagonal in row MU, columns 2 through N.
c    The second superdiagonal in row MU-1, columns 3 through N.
c    The MU-th superdiagonal in row 1, columns MU+1 through N.
c
c  Example:
c
c    N = 5, MU = 2
c
c    A11 A12 A13   0   0
c      0 A22 A23 A24   0
c      0   0 A33 A34 A35
c      0   0   0 A44 A45
c      0   0   0   0 A55
c                --- ---
c                    ---
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 March 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of columns of the matrix.
c
c    Input, integer MU, the upper bandwidth.
c
c    Input/output, integer SEED, a seed for the random number 
c    generator.
c
c    Output, double precision A(MU+1,N), the R8BUT matrix.
c
      implicit none

      integer mu
      integer n

      double precision a(mu+1,n)
      double precision r8_uniform_01
      integer i
      integer j
      integer seed

      do i = 1, mu + 1

        do j = 1, mu + 1 - i
          a(i,j) = 0.0D+00
        end do

        do j = max ( 1, mu + 2 - i ), n
          a(i,j) = r8_uniform_01 ( seed )
        end do

      end do

      return
      end
      subroutine r8but_sl ( n, mu, a, b, job )

c*********************************************************************72
c
cc R8BUT_SL solves an R8BUT system.
c
c  Discussion:
c
c    The R8BUT storage format is used for a banded upper triangular matrix.
c    The matrix is assumed to be zero above the MU-th superdiagonal.
c    The matrix is stored in an MU+1 by N array.
c    Columns are preserved.
c
c    The diagonal is stored in row MU+1 of the array.
c    The first superdiagonal in row MU, columns 2 through N.
c    The second superdiagonal in row MU-1, columns 3 through N.
c    The MU-th superdiagonal in row 1, columns MU+1 through N.
c
c  Example:
c
c    N = 5, MU = 2
c
c    A11 A12 A13   0   0
c      0 A22 A23 A24   0
c      0   0 A33 A34 A35
c      0   0   0 A44 A45
c      0   0   0   0 A55
c                --- ---
c                    ---
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 March 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, integer MU, the upper bandwidth.
c
c    Input, double precision A(MU+1,N), the R8BUT matrix.
c
c    Input/output, double precision B(N).
c    On input, the right hand side.
c    On output, the solution vector.
c
c    Input, integer JOB, is 0 to solve the untransposed system,
c    nonzero to solve the transposed system.
c
      implicit none

      integer mu
      integer n

      double precision a(mu+1,n)
      double precision b(n)
      integer i
      integer ihi
      integer j
      integer jlo
      integer job

      if ( job .eq. 0 ) then

        do j = n, 1, -1
          b(j) = b(j) / a(j-j+mu+1,j)
          jlo = max ( 1, j - mu )
          do i = jlo, j - 1
            b(i) = b(i) - a(i-j+mu+1,j) * b(j)
          end do
        end do

      else

        do j = 1, n
          b(j) = b(j) / a(j-j+mu+1,j)
          ihi = min ( n, j + mu )
          do i = j + 1, ihi
            b(i) = b(i) - a(j-i+mu+1,i) * b(j)
          end do
        end do

      end if

      return
      end
      subroutine r8but_to_r8ge ( n, mu, a, b )

c*********************************************************************72
c
cc R8BUT_TO_R8GE copies an R8BUT matrix to an R8GE matrix.
c
c  Discussion:
c
c    The R8BUT storage format is used for a banded upper triangular matrix.
c    The matrix is assumed to be zero above the MU-th superdiagonal.
c    The matrix is stored in an MU+1 by N array.
c    Columns are preserved.
c
c    The diagonal is stored in row MU+1 of the array.
c    The first superdiagonal in row MU, columns 2 through N.
c    The second superdiagonal in row MU-1, columns 3 through N.
c    The MU-th superdiagonal in row 1, columns MU+1 through N.
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c  Example:
c
c    N = 5, MU = 2
c
c    A11 A12 A13   0   0
c      0 A22 A23 A24   0
c      0   0 A33 A34 A35
c      0   0   0 A44 A45
c      0   0   0   0 A55
c                --- ---
c                    ---
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 March 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrices.
c    N must be positive.
c
c    Input, integer MU, the upper bandwidth of A.
c    MU must be nonnegative, and no greater than N-1.
c
c    Input, double precision A(MU+1,N), the R8BUT matrix.
c
c    Output, double precision B(N,N), the R8GE matrix.
c
      implicit none

      integer mu
      integer n

      double precision a(mu+1,n)
      double precision b(n,n)
      integer i
      integer j

      do i = 1, n
        do j = 1, n
          if ( i .le. j .and. j .le. i+mu ) then
            b(i,j) = a(i-j+mu+1,j)
          else
            b(i,j) = 0.0D+00
          end if
        end do
      end do

      return
      end
      subroutine r8but_vxm ( n, mu, a, x, b )

c*********************************************************************72
c
cc R8BUT_VXM multiplies an R8VECr by an R8BUT matrix.
c
c  Discussion:
c
c    The R8BUT storage format is used for a banded upper triangular matrix.
c    The matrix is assumed to be zero above the MU-th superdiagonal.
c    The matrix is stored in an MU+1 by N array.
c    Columns are preserved.
c
c    The diagonal is stored in row MU+1 of the array.
c    The first superdiagonal in row MU, columns 2 through N.
c    The second superdiagonal in row MU-1, columns 3 through N.
c    The MU-th superdiagonal in row 1, columns MU+1 through N.
c
c  Example:
c
c    N = 5, MU = 2
c
c    A11 A12 A13   0   0
c      0 A22 A23 A24   0
c      0   0 A33 A34 A35
c      0   0   0 A44 A45
c      0   0   0   0 A55
c                --- ---
c                    ---
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 March 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, integer MU, the upper bandwidth.
c
c    Input, double precision A(MU+1,N), the R8BUT matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(N), the product X*A.
c
      implicit none

      integer mu
      integer n

      double precision a(mu+1,n)
      double precision b(n)
      integer i
      integer ilo
      integer j
      double precision x(n)

      do i = 1, n
        ilo = max ( 1, i - mu )
        b(i) = 0.0D+00
        do j = ilo, i
          b(i) = b(i) + x(j) * a(j-i+mu+1,i)
        end do
      end do

      return
      end
      subroutine r8cb_det ( n, ml, mu, a_lu, det )

c*********************************************************************72
c
cc R8CB_DET computes the determinant of an R8CB matrix factored by R8CB_NP_FA.
c
c  Discussion:
c
c    The R8CB storage format is used for a compact banded matrix.
c    It is assumed that the matrix has lower and upper bandwidths ML and MU,
c    respectively.  The matrix is stored in a way similar to that used
c    by LINPACK and LAPACK for a general banded matrix, except that in
c    this mode, no extra rows are set aside for possible fillin during pivoting.
c    Thus, this storage mode is suitable if you do not intend to factor
c    the matrix, or if you can guarantee that the matrix can be factored
c    without pivoting.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 June 2009
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
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than N-1.
c
c    Input, double precision A_LU(ML+MU+1,N), the LU factors from R8CB_NP_FA.
c
c    Output, double precision DET, the determinant of the matrix.
c
      implicit none

      integer ml
      integer mu
      integer n

      double precision a_lu(ml+mu+1,n)
      double precision det
      double precision r8vec_product

      det = r8vec_product ( n, a_lu(mu+1,1) )

      return
      end
      subroutine r8cb_indicator ( m, n, ml, mu, a )

c*********************************************************************72
c
cc R8CB_INDICATOR sets up an R8CB indicator matrix.
c
c  Discussion:
c
c    The "indicator matrix" simply has a value like I*10+J at every
c    entry of a dense matrix, or at every entry of a compressed storage
c    matrix for which storage is allocated. 
c
c    The R8CB storage format is used for a compact banded matrix.
c    It is assumed that the matrix has lower and upper bandwidths ML and MU,
c    respectively.  The matrix is stored in a way similar to that used
c    by LINPACK and LAPACK for a general banded matrix, except that in
c    this mode, no extra rows are set aside for possible fillin during pivoting.
c    Thus, this storage mode is suitable if you do not intend to factor
c    the matrix, or if you can guarantee that the matrix can be factored
c    without pivoting.
c
c    The original M by N matrix is "collapsed" downward, so that diagonals
c    become rows of the storage array, while columns are preserved.  The
c    collapsed array is logically ML+MU+1 by N.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than min(M,N)-1.
c
c    Output, double precision A(ML+MU+1,N), the R8CB matrix.
c
      implicit none

      integer ml
      integer mu
      integer m
      integer n

      double precision a(ml+mu+1,n)
      integer diag
      integer fac
      integer i
      integer i4_log_10
      integer j
      integer k
      double precision value

      fac = 10 ** ( i4_log_10 ( n ) + 1 )
      k = 0

      do j = 1, n
        do diag = 1, ml + mu + 1

          i = diag + j - mu - 1

          if ( 1 .le. i .and. i .le. m .and. 
     &         i - ml .le. j .and. j .le. i + mu ) then
            value = dble ( fac * i + j )
          else
            k = k + 1
            value = - dble ( k )
          end if

          a(diag,j) = value

        end do
      end do

      return
      end
      subroutine r8cb_ml ( n, ml, mu, a_lu, x, b, job )

c*********************************************************************72
c
cc R8CB_ML computes A * x or A' * X, using R8CB_NP_FA factors.
c
c  Discussion:
c
c    The R8CB storage format is used for a compact banded matrix.
c    It is assumed that the matrix has lower and upper bandwidths ML and MU,
c    respectively.  The matrix is stored in a way similar to that used
c    by LINPACK and LAPACK for a general banded matrix, except that in
c    this mode, no extra rows are set aside for possible fillin during pivoting.
c    Thus, this storage mode is suitable if you do not intend to factor
c    the matrix, or if you can guarantee that the matrix can be factored
c    without pivoting.
c
c    It is assumed that R8CB_NP_FA has overwritten the original matrix
c    information by LU factors.  R8CB_ML is able to reconstruct the
c    original matrix from the LU factor data.
c
c    R8CB_ML allows the user to check that the solution of a linear
c    system is correct, without having to save an unfactored copy
c    of the matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 June 2009
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
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than N-1.
c
c    Input, double precision A_LU(ML+MU+1,N), the LU factors from R8CB_NP_FA.
c
c    Input, double precision X(N), the vector to be multiplied.
c
c    Output, double precision B(N), the result of the multiplication.
c
c    Input, integer JOB, specifies the operation to be done:
c    JOB = 0, compute A * x.
c    JOB nonzero, compute A' * x.
c
      implicit none

      integer ml
      integer mu
      integer n

      double precision a_lu(ml+mu+1,n)
      double precision b(n)
      integer i
      integer ihi
      integer ilo
      integer j
      integer jhi
      integer job
      double precision x(n)

      do i = 1, n
        b(i) = x(i)
      end do

      if ( job .eq. 0 ) then
c
c  Y = U * X.
c
        do j = 1, n
          ilo = max ( 1, j - mu )
          do i = ilo, j - 1
            b(i) = b(i) + a_lu(i-j+mu+1,j) * b(j)
          end do
          b(j) = a_lu(j-j+mu+1,j) * b(j)
        end do
c
c  B = PL * Y = PL * U * X = A * x.
c
        do j = n-1, 1, -1
          ihi = min ( n, j + ml )
          do i = j + 1, ihi
            b(i) = b(i) - a_lu(mu+1+i-j,j) * b(j)
          end do
        end do

      else
c
c  Y = ( PL )' * X.
c
        do j = 1, n-1

          ihi = min ( n, j + ml )
          do i = j+1, ihi
            b(j) = b(j) - b(i) * a_lu(i-j+mu+1,j)
          end do

        end do
c
c  B = U' * Y = ( PL * U )' * X = A' * X.
c
        do i = n, 1, -1
          jhi = min ( n, i + mu )
          do j = i+1, jhi
            b(j) = b(j) + b(i) * a_lu(i-j+mu+1,j)
          end do
          b(i) = b(i) * a_lu(i-i+mu+1,i)
        end do

      end if

      return
      end
      subroutine r8cb_mxv ( n, ml, mu, a, x, b )

c*********************************************************************72
c
cc R8CB_MXV multiplies an R8CB matrix by an R8VEC.
c
c  Discussion:
c
c    The R8CB storage format is used for a compact banded matrix.
c    It is assumed that the matrix has lower and upper bandwidths ML and MU,
c    respectively.  The matrix is stored in a way similar to that used
c    by LINPACK and LAPACK for a general banded matrix, except that in
c    this mode, no extra rows are set aside for possible fillin during pivoting.
c    Thus, this storage mode is suitable if you do not intend to factor
c    the matrix, or if you can guarantee that the matrix can be factored
c    without pivoting.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 June 2009
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
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than N-1.
c
c    Input, double precision A(ML+MU+1,N), the R8CB matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(N), the product A * x.
c
      implicit none

      integer ml
      integer mu
      integer n

      double precision a(ml+mu+1,n)
      double precision b(n)
      integer i
      integer j
      integer jhi
      integer jlo
      double precision x(n)

      do i = 1, n
        b(i) = 0.0D+00
        jlo = max ( 1, i - ml )
        jhi = min ( n, i + mu )
        do j = jlo, jhi
          b(i) = b(i) + a(i-j+mu+1,j) * x(j)
        end do
      end do

      return
      end
      subroutine r8cb_np_fa ( n, ml, mu, a, info )

c*********************************************************************72
c
cc R8CB_NP_FA factors an R8CB matrix by Gaussian elimination.
c
c  Discussion:
c
c    The R8CB storage format is appropriate for a compact banded matrix.
c    It is assumed that the matrix has lower and upper bandwidths ML and MU,
c    respectively.  The matrix is stored in a way similar to that used
c    by LINPACK and LAPACK for a general banded matrix, except that in
c    this mode, no extra rows are set aside for possible fillin during pivoting.
c    Thus, this storage mode is suitable if you do not intend to factor
c    the matrix, or if you can guarantee that the matrix can be factored
c    without pivoting.
c
c    R8CB_NP_FA is a version of the LINPACK routine R8GBFA, modifed to use
c    no pivoting, and to be applied to the R8CB compressed band matrix storage
c    format.  It will fail if the matrix is singular, or if any zero
c    pivot is encountered.
c
c    If R8CB_NP_FA successfully factors the matrix, R8CB_NP_SL may be called
c    to solve linear systems involving the matrix.
c
c    The matrix is stored in a compact version of LINPACK general
c    band storage, which does not include the fill-in entires.
c    The following program segment will store the entries of a banded
c    matrix in the compact format used by this routine:
c
c      m = mu+1
c      do j = 1, n
c        i1 = max ( 1, j-mu )
c        i2 = min ( n, j+ml )
c        do i = i1, i2
c          k = i-j+m
c          a(k,j) = afull(i,j)
c        end do
c      end do
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 June 2009
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
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than N-1.
c
c    Input/output, double precision A(ML+MU+1,N), the compact band matrix.
c    On input, the coefficient matrix of the linear system.
c    On output, the LU factors of the matrix.
c
c    Output, integer INFO, singularity flag.
c    0, no singularity detected.
c    nonzero, the factorization failed on the INFO-th step.
c
      implicit none

      integer ml
      integer mu
      integer n

      double precision a(ml+mu+1,n)
      integer i
      integer info
      integer j
      integer ju
      integer k
      integer lm
      integer m
      integer mm
c
c  The value of M is MU + 1 rather than ML + MU + 1.
c
      m = mu + 1
      info = 0
      ju = 0

      do k = 1, n-1
c
c  If our pivot entry A(MU+1,K) is zero, then we must give up.
c
        if ( a(m,k) .eq. 0.0D+00 ) then
          info = k
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8CB_NP_FA - Fatal error!'
          write ( *, '(a,i8)' ) '  Zero pivot on step ', info
          stop
        end if
c
c  LM counts the number of nonzero elements that lie below the current
c  diagonal entry, A(K,K).
c
c  Multiply the LM entries below the diagonal by -1/A(K,K), turning
c  them into the appropriate "multiplier" terms in the L matrix.
c
        lm = min ( ml, n-k )
        do i = m + 1, m + lm
          a(i,k) = - a(i,k) / a(m,k)
        end do
c
c  MM points to the row in which the next entry of the K-th row is, A(K,J).
c  We then add L(I,K)*A(K,J) to A(I,J) for rows I = K+1 to K+LM.
c
        ju = max ( ju, mu + k )
        ju = min ( ju, n )
        mm = m

        do j = k+1, ju
          mm = mm - 1
          do i = 1, lm
            a(mm+i,j) = a(mm+i,j) + a(mm,j) * a(m+i,k)
          end do
        end do

      end do

      if ( a(m,n) .eq. 0.0D+00 ) then
        info = n
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8CB_NP_FA - Fatal error!'
        write ( *, '(a,i8)' ) '  Zero pivot on step ', info
        stop
      end if

      return
      end
      subroutine r8cb_np_sl ( n, ml, mu, a_lu, b, job )

c*********************************************************************72
c
cc R8CB_NP_SL solves an R8CB system factored by R8CB_NP_FA.
c
c  Discussion:
c
c    The R8CB storage format is used for a compact banded matrix.
c    It is assumed that the matrix has lower and upper bandwidths ML and MU,
c    respectively.  The matrix is stored in a way similar to that used
c    by LINPACK and LAPACK for a general banded matrix, except that in
c    this mode, no extra rows are set aside for possible fillin during pivoting.
c    Thus, this storage mode is suitable if you do not intend to factor
c    the matrix, or if you can guarantee that the matrix can be factored
c    without pivoting.
c
c    R8CB_NP_SL can also solve the related system A' * x = b.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 June 2009
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
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than N-1.
c
c    Input, double precision A_LU(ML+MU+1,N), the LU factors from R8CB_NP_FA.
c
c    Input/output, double precision B(N).
c    On input, B contains the right hand side of the linear system, B.
c    On output, B contains the solution of the linear system, X.
c
c    Input, integer JOB.
c    If JOB is zero, the routine will solve A * x = b.
c    If JOB is nonzero, the routine will solve A' * x = b.
c
      implicit none

      integer ml
      integer mu
      integer n

      double precision a_lu(ml+mu+1,n)
      double precision b(n)
      integer i
      integer job
      integer k
      integer la
      integer lb
      integer lm
      integer m
c
c  The value of M is ML + 1, rather than MU + ML + 1.
c
      m = mu + 1
c
c  Solve A * x = b.
c
      if ( job .eq. 0 ) then
c
c  Solve PL * Y = B.
c
        if ( 0 < ml ) then
          do k = 1, n-1
            lm = min ( ml, n-k )
            do i = 1, lm
              b(k+i) = b(k+i) + b(k) * a_lu(m+i,k)
            end do
          end do
        end if
c
c  Solve U * X = Y.
c
        do k = n, 1, -1

          b(k) = b(k) / a_lu(m,k)
          lm = min ( k, m ) - 1
          la = m - lm
          lb = k - lm
          do i = 0, lm - 1
            b(lb+i) = b(lb+i) - b(k) * a_lu(la+i,k)
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
          lm = min ( k, m ) - 1
          la = m - lm
          lb = k - lm
          do i = 0, lm - 1
            b(k) = b(k) - a_lu(la+i,k) * b(lb+i)
          end do
          b(k) = b(k) / a_lu(m,k)
        end do
c
c  Solve ( PL )' * X = Y.
c
        if ( 0 < ml ) then

          do k = n-1, 1, -1
            lm = min ( ml, n-k )
            do i = 1, lm
              b(k) = b(k) + a_lu(m+i,k) * b(k+i)
            end do
          end do

        end if

      end if

      return
      end
      subroutine r8cb_print ( m, n, ml, mu, a, title )

c*********************************************************************72
c
cc R8CB_PRINT prints an R8CB matrix.
c
c  Discussion:
c
c    The R8CB storage format is used for a compact banded matrix.
c    It is assumed that the matrix has lower and upper bandwidths ML and MU,
c    respectively.  The matrix is stored in a way similar to that used
c    by LINPACK and LAPACK for a general banded matrix, except that in
c    this mode, no extra rows are set aside for possible fillin during pivoting.
c    Thus, this storage mode is suitable if you do not intend to factor
c    the matrix, or if you can guarantee that the matrix can be factored
c    without pivoting.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns of 
c    the matrix.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than min(M,N)-1.
c
c    Input, double precision A(ML+MU+1,N), the R8CB matrix.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer ml
      integer mu
      integer n

      double precision a(ml+mu+1,n)
      integer m
      character * ( * )  title

      call r8cb_print_some ( m, n, ml, mu, a, 1, 1, m, n, title )

      return
      end
      subroutine r8cb_print_some ( m, n, ml, mu, a, ilo, jlo, ihi, 
     &  jhi, title )

c*********************************************************************72
c
cc R8CB_PRINT_SOME prints some of an R8CB matrix.
c
c  Discussion:
c
c    The R8CB storage format is used for a compact banded matrix.
c    It is assumed that the matrix has lower and upper bandwidths ML and MU,
c    respectively.  The matrix is stored in a way similar to that used
c    by LINPACK and LAPACK for a general banded matrix, except that in
c    this mode, no extra rows are set aside for possible fillin during pivoting.
c    Thus, this storage mode is suitable if you do not intend to factor
c    the matrix, or if you can guarantee that the matrix can be factored
c    without pivoting.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns of 
c    the matrix.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than min(M,N)-1.
c
c    Input, double precision A(ML+MU+1,N), the R8CB matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 5 )
      integer ml
      integer mu
      integer n

      double precision a(ml+mu+1,n)
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
      integer m
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
          write ( ctemp(j2), '(i7,7x)' ) j
        end do

        write ( *, '(a,5a14)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) '  ---'
c
c  Determine the range of the rows in this strip.
c
        i2lo = max ( ilo, 1 )
        i2lo = max ( i2lo, j2lo - mu )
        i2hi = min ( ihi, m )
        i2hi = min ( i2hi, j2hi + ml )

        do i = i2lo, i2hi
c
c  Print out (up to) 5 entries in row I, that lie in the current strip.
c
          do j2 = 1, inc

            j = j2lo - 1 + j2

            if ( ml .lt. i - j .or. mu .lt. j - i ) then
              ctemp(j2) = '              '
            else if ( r8_is_int ( a(i-j+mu+1,j) ) ) then
              write ( ctemp(j2), '(f8.0,6x)' ) a(i-j+mu+1,j)
            else
              write ( ctemp(j2), '(g14.6)' ) a(i-j+mu+1,j)
            end if

          end do

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine r8cb_random ( n, ml, mu, seed, a )

c*********************************************************************72
c
cc R8CB_RANDOM randomizes an R8CB matrix.
c
c  Discussion:
c
c    The R8CB storage format is used for a compact banded matrix.
c    It is assumed that the matrix has lower and upper bandwidths ML and MU,
c    respectively.  The matrix is stored in a way similar to that used
c    by LINPACK and LAPACK for a general banded matrix, except that in
c    this mode, no extra rows are set aside for possible fillin during pivoting.
c    Thus, this storage mode is suitable if you do not intend to factor
c    the matrix, or if you can guarantee that the matrix can be factored
c    without pivoting.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 June 2009
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
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than N-1.
c
c    Input/output, integer SEED, a seed for the random 
c    number generator.
c
c    Output, double precision A(ML+MU+1,N), the R8CB matrix.
c
      implicit none

      integer ml
      integer mu
      integer n

      double precision a(ml+mu+1,n)
      double precision r8_uniform_01
      integer i
      integer ihi
      integer ilo
      integer j
      integer seed
c
c  Set the entries that correspond to matrix elements.
c
      do j = 1, n

        ilo = max ( 1, j - mu )
        ihi = min ( n, j + ml )

        do i = j - mu, 0
          a(i-j+mu+1,j) = 0.0D+00
        end do

        do i = ilo, ihi
          a(i-j+mu+1,j) = r8_uniform_01 ( seed )
        end do

        do i = n + 1, j + ml
          a(i-j+mu+1,j) = 0.0D+00 
        end do

      end do

      return
      end
      subroutine r8cb_to_r8vec ( m, n, ml, mu, a, x )

c*********************************************************************72
c
cc R8CB_TO_R8VEC copies an R8CB matrix to an R8VEC.
c
c  Discussion:
c
c    In C++ and FORTRAN, this routine is not really needed.  In MATLAB,
c    a data item carries its dimensionality implicitly, and so cannot be
c    regarded sometimes as a vector and sometimes as an array.
c
c    The R8CB storage format is used for a compact banded matrix.
c    It is assumed that the matrix has lower and upper bandwidths ML and MU,
c    respectively.  The matrix is stored in a way similar to that used
c    by LINPACK and LAPACK for a general banded matrix, except that in
c    this mode, no extra rows are set aside for possible fillin during pivoting.
c    Thus, this storage mode is suitable if you do not intend to factor
c    the matrix, or if you can guarantee that the matrix can be factored
c    without pivoting.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns in 
c    the array.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c
c    Input, double precision A(ML+MU+1,N), the array to be copied.
c
c    Output, double precision X((ML+MU+1)*N), the vector.
c
      implicit none

      integer m
      integer ml
      integer mu
      integer n

      double precision a(ml+mu+1,n)
      integer i
      integer ihi
      integer ilo
      integer j
      double precision x((ml+mu+1)*n)

      do j = 1, n

        ihi = min ( mu, mu + 1 - j )
        do i = 1, ihi
          x(i+(j-1)*(ml+mu+1)) = 0.0D+00
        end do

        ilo = max ( ihi + 1, 1 )
        ihi = min ( ml+mu+1, mu+1+m-j )
        do i = ilo, ihi
          x(i+(j-1)*(ml+mu+1)) = a(i,j)
        end do

        ilo = ihi + 1
        ihi = ml+mu+1
        do i = ilo, ihi
          x(i+(j-1)*(ml+mu+1)) = 0.0D+00
        end do

      end do

      return
      end
      subroutine r8cb_to_r8ge ( n, ml, mu, a, b )

c*********************************************************************72
c
cc R8CB_TO_R8GE copies an R8CB matrix to an R8GE matrix.
c
c  Discussion:
c
c    The R8CB storage format is used for a compact banded matrix.
c    It is assumed that the matrix has lower and upper bandwidths ML and MU,
c    respectively.  The matrix is stored in a way similar to that used
c    by LINPACK and LAPACK for a general banded matrix, except that in
c    this mode, no extra rows are set aside for possible fillin during pivoting.
c    Thus, this storage mode is suitable if you do not intend to factor
c    the matrix, or if you can guarantee that the matrix can be factored
c    without pivoting.
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrices.
c    N must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths of A.
c    ML and MU must be nonnegative, and no greater than N-1.
c
c    Input, double precision A(ML+MU+1,N), the R8CB matrix.
c
c    Output, double precision B(N,N), the R8GE matrix.
c
      implicit none

      integer ml
      integer mu
      integer n

      double precision a(ml+mu+1,n)
      double precision b(n,n)
      integer i
      integer j

      do i = 1, n
        do j = 1, n
          if ( j-mu .le. i .and. i .le. j+ml ) then
            b(i,j) = a(mu+1+i-j,j)
          else
            b(i,j) = 0.0D+00
          end if
        end do
      end do

      return
      end
      subroutine r8cb_vxm ( n, ml, mu, a, x, b )

c*********************************************************************72
c
cc R8CB_VXM multiplies an R8VECr by an R8CB matrix.
c
c  Discussion:
c
c    The R8CB storage format is used for a compact banded matrix.
c    It is assumed that the matrix has lower and upper bandwidths ML and MU,
c    respectively.  The matrix is stored in a way similar to that used
c    by LINPACK and LAPACK for a general banded matrix, except that in
c    this mode, no extra rows are set aside for possible fillin during pivoting.
c    Thus, this storage mode is suitable if you do not intend to factor
c    the matrix, or if you can guarantee that the matrix can be factored
c    without pivoting.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 June 2009
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
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than N-1.
c
c    Input, double precision A(ML+MU+1,N), the R8CB matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(N), the product X*A.
c
      implicit none

      integer ml
      integer mu
      integer n

      double precision a(ml+mu+1,n)
      double precision b(n)
      integer i
      integer j
      integer jhi
      integer jlo
      double precision x(n)

      do i = 1, n
        b(i) = 0.0D+00
      end do

      do i = 1, n
        jlo = max ( 1, i - ml )
        jhi = min ( n, i + mu )
        do j = jlo, jhi
          b(j) = b(j) + x(i) * a(i-j+mu+1,j)
        end do
      end do

      return
      end
      subroutine r8cbb_add ( n1, n2, ml, mu, a, i, j, value )

c*********************************************************************72
c
cc R8CBB_ADD adds a value to an entry of an R8CBB matrix.
c
c  Discussion:
c
c    The R8CBB storage format is for a compressed border banded matrix.  
c    Such a matrix has the logical form:
c
c      A1 | A2
c      ---+---
c      A3 | A4
c
c    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
c    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
c    respectively.  
c
c    The R8CBB format is the same as the R8BB format, except that the banded
c    matrix A1 is stored in compressed band form rather than standard
c    banded form.  In other words, we do not include the extra room
c    set aside for fill in during pivoting.
c
c    A should be defined as a vector.  The user must then store
c    the entries of the four blocks of the matrix into the vector A.
c    Each block is stored by columns.
c
c    A1, the banded portion of the matrix, is stored in
c    the first (ML+MU+1)*N1 entries of A, using the obvious variant
c    of the LINPACK general band format.
c
c    The following formulas should be used to determine how to store
c    the entry corresponding to row I and column J in the original matrix:
c
c    Entries of A1:
c
c      1 .le. I .le. N1, 1 .le. J .le. N1, (J-I) .le. MU and (I-J) .le. ML.
c
c      Store the I, J entry into location
c      (I-J+MU+1)+(J-1)*(ML+MU+1).
c
c    Entries of A2:
c
c      1 .le. I .le. N1, N1+1 .le. J .le. N1+N2.
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+(J-N1-1)*N1+I.
c
c    Entries of A3:
c
c      N1+1 .le. I .le. N1+N2, 1 .le. J .le. N1.
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c
c    Entries of A4:
c
c      N1+1 .le. I .le. N1+N2, N1+1 .le. J .le. N1+N2
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c      (same formula used for A3).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, the order of the banded and dense 
c    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than N1-1.
c
c    Input/output, double precision A((ML+MU+1)*N1 + 2*N1*N2 + N2*N2),
c    the R8CBB matrix.
c
c    Input, integer I, J, the indices of the entry to be 
c    incremented.
c
c    Input, double precision VALUE, the value to be added to the (I,J) entry.
c
      implicit none

      integer ml
      integer mu
      integer n1
      integer n2

      double precision a((ml+mu+1)*n1+2*n1*n2+n2*n2)
      integer i
      integer ij
      integer j
      double precision value

      if ( value .eq. 0.0D+00 ) then
        return
      end if
c
c  Check for I or J out of bounds.
c
      if ( i .le. 0 .or. n1+n2 .lt. i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8CBB_ADD - Fatal error!'
        write ( *, '(a,i8)' ) 
     &    '  Illegal input value of row index I = ', i
        stop
      end if

      if ( j .le. 0 .or. n1+n2 .lt. j ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8CBB_ADD - Fatal error!'
        write ( *, '(a,i8)' ) 
     &    '  Illegal input value of column index J = ',j
        stop
      end if
c
c  The A1 block of the matrix.
c
c  Check for out of band problems.
c
      if ( i .le. n1 .and. j .le. n1 ) then
        if ( mu .lt. (j-i) .or. ml .lt. (i-j) ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8CBB_ADD - Warning!'
          write ( *, '(a,i8,a,i8,a)' ) 
     &      '  Unable to add to entry (', i, ',', j, ').'
          return
        else
          ij = (i-j+mu+1)+(j-1)*(ml+mu+1)
        end if
c
c  The A2 block of the matrix:
c
      else if ( i .le. n1 .and. n1 .lt. j ) then
        ij = (ml+mu+1)*n1+(j-n1-1)*n1 + i
c
c  The A3 and A4 blocks of the matrix.
c
      else if ( n1 .lt. i ) then
        ij = (ml+mu+1)*n1+n2*n1+(j-1)*n2 + (i-n1)
      end if

      a(ij) = a(ij) + value

      return
      end
      subroutine r8cbb_error ( n1, n2, ml, mu, ierror )

c*********************************************************************72
c
cc R8CBB_ERROR checks the dimensions of an R8CBB matrix.
c
c  Discussion:
c
c    The R8CBB storage format is for a compressed border banded matrix.  
c    Such a matrix has the logical form:
c
c      A1 | A2
c      ---+---
c      A3 | A4
c
c    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
c    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
c    respectively.  
c
c    The R8CBB format is the same as the R8BB format, except that the banded
c    matrix A1 is stored in compressed band form rather than standard
c    banded form.  In other words, we do not include the extra room
c    set aside for fill in during pivoting.
c
c    A should be defined as a vector.  The user must then store
c    the entries of the four blocks of the matrix into the vector A.
c    Each block is stored by columns.
c
c    A1, the banded portion of the matrix, is stored in
c    the first (ML+MU+1)*N1 entries of A, using the obvious variant
c    of the LINPACK general band format.
c
c    The following formulas should be used to determine how to store
c    the entry corresponding to row I and column J in the original matrix:
c
c    Entries of A1:
c
c      1 .le. I .le. N1, 1 .le. J .le. N1, (J-I) .le. MU and (I-J) .le. ML.
c
c      Store the I, J entry into location
c      (I-J+MU+1)+(J-1)*(ML+MU+1).
c
c    Entries of A2:
c
c      1 .le. I .le. N1, N1+1 .le. J .le. N1+N2.
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+(J-N1-1)*N1+I.
c
c    Entries of A3:
c
c      N1+1 .le. I .le. N1+N2, 1 .le. J .le. N1.
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c
c    Entries of A4:
c
c      N1+1 .le. I .le. N1+N2, N1+1 .le. J .le. N1+N2
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c      (same formula used for A3).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, the order of the banded and dense 
c    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative and no greater than N1 - 1.
c
c    Output, integer IERROR, reports whether errors were detected.
c    IERROR is set to 0 before the checks are made, and then:
c    IERROR = IERROR + 1 if ML is illegal;
c    IERROR = IERROR + 2 if MU is illegal;
c    IERROR = IERROR + 4 if N1 is illegal;
c    IERROR = IERROR + 8 if N2 is illegal;
c    IERROR = IERROR + 16 if neither N1 nor N2 is positive.
c
      implicit none

      integer ierror
      integer ml
      integer mu
      integer n1
      integer n2

      ierror = 0

      if ( ml .lt. 0 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) 'R8CBB_ERROR - Illegal ML = ', ml
        write ( *, '(a)' ) 
     &    '  but ML must be greater than or equal to 0.'
      else if ( max ( n1 - 1, 0 ) .lt. ml ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) 'R8CBB_ERROR - Illegal ML = ', ml
        write ( *, '(a,i8)' ) 
     &  '  but ML must be .le. Max ( N1 - 1, 0 ) = ', 
     &    max ( n1 - 1, 0 )
      end if

      if ( mu .lt. 0  ) then
        ierror = ierror + 2
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) 'R8CBB_ERROR - Illegal MU = ', mu
        write ( *, '(a)' ) 
     &    '  but MU must be greater than or equal to 0.'
      else if ( max ( n1 - 1, 0 ) .lt. ml ) then
        ierror = ierror + 2
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) 'R8CBB_ERROR - Illegal MU = ', mu
        write ( *, '(a,i8)' ) 
     &  '  but MU must be .le. Max ( N1 - 1, 0 ) = ', 
     &    max ( n1 - 1, 0 )
      end if

      if ( n1 .lt. 0 ) then
        ierror = ierror + 4
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) 'R8CBB_ERROR - Illegal N1 = ', n1
      end if

      if ( n2 .lt. 0 ) then
        ierror = ierror + 8
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) 'R8CBB_ERROR - Illegal N2 = ', n2
      end if

      if ( n1 + n2 .le. 0 ) then
        ierror = ierror + 16
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) 'R8CBB_ERROR - Illegal N1+N2 = ', n1+n2
      end if

      return
      end
      subroutine r8cbb_fa ( n1, n2, ml, mu, a, info )

c*********************************************************************72
c
cc R8CBB_FA factors an R8CBB matrix.
c
c  Discussion:
c
c    The R8CBB storage format is for a compressed border banded matrix.  
c    Such a matrix has the logical form:
c
c      A1 | A2
c      ---+---
c      A3 | A4
c
c    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
c    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
c    respectively.  
c
c    The R8CBB format is the same as the R8BB format, except that the banded
c    matrix A1 is stored in compressed band form rather than standard
c    banded form.  In other words, we do not include the extra room
c    set aside for fill in during pivoting.
c
c    A should be defined as a vector.  The user must then store
c    the entries of the four blocks of the matrix into the vector A.
c    Each block is stored by columns.
c
c    A1, the banded portion of the matrix, is stored in
c    the first (ML+MU+1)*N1 entries of A, using the obvious variant
c    of the LINPACK general band format.
c
c    The following formulas should be used to determine how to store
c    the entry corresponding to row I and column J in the original matrix:
c
c    Entries of A1:
c
c      1 .le. I .le. N1, 1 .le. J .le. N1, (J-I) .le. MU and (I-J) .le. ML.
c
c      Store the I, J entry into location
c      (I-J+MU+1)+(J-1)*(ML+MU+1).
c
c    Entries of A2:
c
c      1 .le. I .le. N1, N1+1 .le. J .le. N1+N2.
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+(J-N1-1)*N1+I.
c
c    Entries of A3:
c
c      N1+1 .le. I .le. N1+N2, 1 .le. J .le. N1.
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c
c    Entries of A4:
c
c      N1+1 .le. I .le. N1+N2, N1+1 .le. J .le. N1+N2
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c      (same formula used for A3).
c
c
c    Once the matrix has been factored by SCCB_FA, SCCB_SL may be called
c    to solve linear systems involving the matrix.
c
c    SCCB_FA uses special non-pivoting versions of LINPACK routines to
c    carry out the factorization.  The special version of the banded
c    LINPACK solver also results in a space saving, since no entries
c    need be set aside for fill in due to pivoting.
c
c    The linear system must be border banded, of the form:
c
c      ( A1 A2 ) (X1) = (B1)
c      ( A3 A4 ) (X2)   (B2)
c
c    where A1 is a (usually big) banded square matrix, A2 and A3 are
c    column and row strips which may be nonzero, and A4 is a dense
c    square matrix.
c
c    The algorithm rewrites the system as:
c
c         X1 + inverse(A1) A2 X2 = inverse(A1) B1
c
c      A3 X1 +             A4 X2 = B2
c
c    and then rewrites the second equation as
c
c      ( A4 - A3 inverse(A1) A2 ) X2 = B2 - A3 inverse(A1) B1
c
c    The algorithm will certainly fail if the matrix A1 is singular,
c    or requires pivoting.  The algorithm will also fail if the A4 matrix,
c    as modified during the process, is singular, or requires pivoting.
c    All these possibilities are in addition to the failure that will
c    if the total matrix A is singular.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, the order of the banded and dense 
c    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than N1-1.
c
c    Input/output, double precision A( (ML+MU+1)*N1 + 2*N1*N2 + N2*N2).
c    On input, A contains the compact border-banded coefficient matrix.
c    On output, A contains information describing a partial factorization
c    of the original coefficient matrix.  
c
c    Output, integer INFO, singularity flag.
c    0, no singularity detected.
c    nonzero, the factorization failed on the INFO-th step.
c
      implicit none

      integer ml
      integer mu
      integer n1
      integer n2

      double precision a((ml+mu+1)*n1+2*n1*n2+n2*n2)
      integer i
      integer ij
      integer ik
      integer info
      integer j
      integer jk
      integer job
      integer k
      integer nband

      nband = (ml+mu+1)*n1
c
c  Factor the A1 band matrix, overwriting A1 by its factors.
c
      if ( 0 .lt. n1 ) then

        call r8cb_np_fa ( n1, ml, mu, a, info )

        if ( info .ne. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8CBB_FA - Fatal error!'
          write ( *, '(a,i8)' ) '  R8CB_NP_FA returned INFO = ', info
          write ( *, '(a)' ) '  Factoring failed for column INFO.'
          write ( *, '(a)' ) '  The band matrix A1 is singular.'
          write ( *, '(a)' ) '  This algorithm cannot continue!'
          stop
        end if

      end if

      if ( 0 .lt. n1 .and. 0 .lt. n2 ) then
c
c  Set A2 := -inverse(A1) * A2.
c
        do i = nband + 1, nband + n1 * n2
          a(i) = - a(i)
        end do

        job = 0

        do j = 1, n2
          call r8cb_np_sl ( n1, ml, mu, a, a(nband+(j-1)*n1+1), job )
        end do
c
c  Set A4 := A4 + A3*A2
c
        do i = 1, n2
          do j = 1, n1
            ij = nband + n1*n2 + (j-1)*n2 + i
            do k = 1, n2
              ik = nband + 2*n1*n2 + (k-1)*n2 + i
              jk = nband + (k-1)*n1 + j
              a(ik) = a(ik) + a(ij) * a(jk)
            end do
          end do
        end do

      end if
c
c  Factor A4.
c
      if ( 0 .lt. n2 ) then

        call r8ge_np_fa ( n2, a(nband+2*n1*n2+1), info )

        if ( info .ne. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8CBB_FA - Fatal error!'
          write ( *, '(a,i8)' ) '  R8GE_NP_FA returned INFO = ',info
          write ( *, '(a)' ) 
     &      '  This indicates singularity in column INFO'
          info = n1 + info
          write ( *, '(a,i8)' ) 
     &      '  of the A4 submatrix, which is column ',info
          write ( *, '(a)' ) '  of the full matrix.'
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  It is possible that the full matrix is '
          write ( *, '(a)' ) 
     &      '  nonsingular, but the algorithm R8CBB_FA may'
          write ( *, '(a)' ) '  not be used for this matrix.'
          stop
        end if

      end if

      return
      end
      subroutine r8cbb_get ( n1, n2, ml, mu, a, i, j, value )

c*********************************************************************72
c
cc R8CBB_GET returns the value of an entry of an R8CBB matrix.
c
c  Discussion:
c
c    The R8CBB storage format is for a compressed border banded matrix.  
c    Such a matrix has the logical form:
c
c      A1 | A2
c      ---+---
c      A3 | A4
c
c    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
c    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
c    respectively.  
c
c    The R8CBB format is the same as the R8BB format, except that the banded
c    matrix A1 is stored in compressed band form rather than standard
c    banded form.  In other words, we do not include the extra room
c    set aside for fill in during pivoting.
c
c    A should be defined as a vector.  The user must then store
c    the entries of the four blocks of the matrix into the vector A.
c    Each block is stored by columns.
c
c    A1, the banded portion of the matrix, is stored in
c    the first (ML+MU+1)*N1 entries of A, using the obvious variant
c    of the LINPACK general band format.
c
c    The following formulas should be used to determine how to store
c    the entry corresponding to row I and column J in the original matrix:
c
c    Entries of A1:
c
c      1 .le. I .le. N1, 1 .le. J .le. N1, (J-I) .le. MU and (I-J) .le. ML.
c
c      Store the I, J entry into location
c      (I-J+MU+1)+(J-1)*(ML+MU+1).
c
c    Entries of A2:
c
c      1 .le. I .le. N1, N1+1 .le. J .le. N1+N2.
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+(J-N1-1)*N1+I.
c
c    Entries of A3:
c
c      N1+1 .le. I .le. N1+N2, 1 .le. J .le. N1.
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c
c    Entries of A4:
c
c      N1+1 .le. I .le. N1+N2, N1+1 .le. J .le. N1+N2
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c      (same formula used for A3).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, the order of the banded and dense 
c    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than N1-1.
c
c    Input, double precision A((ML+MU+1)*N1 + 2*N1*N2 + N2*N2),
c    the R8CBB matrix.
c
c    Input, integer I, J, the row and column of the entry to 
c    retrieve.
c
c    Output, double precision VALUE, the value of the (I,J) entry.
c
      implicit none

      integer ml
      integer mu
      integer n1
      integer n2

      double precision a((ml+mu+1)*n1+2*n1*n2+n2*n2)
      integer i
      integer ij
      integer j
      double precision value
c
c  Check for I or J out of bounds.
c
      if ( i .le. 0 .or. n1+n2 .lt. i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8CBB_GET - Fatal error!'
        write ( *, '(a,i8)' ) 
     &    '  Illegal input value of row index I = ', i
        stop
      end if

      if ( j .le. 0 .or. n1+n2 .lt. j ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8CBB_GET - Fatal error!'
        write ( *, '(a,i8)' ) 
     &    '  Illegal input value of column index J = ', j
        stop
      end if
c
c  The A1 block of the matrix.
c
c  Check for out of band problems.
c
      if ( i .le. n1 .and. j .le. n1 ) then
        if ( mu .lt. (j-i) .or. ml .lt. (i-j) ) then
          value = 0.0D+00
          return
        else
          ij = (i-j+mu+1)+(j-1)*(ml+mu+1)
        end if
c
c  The A2 block of the matrix:
c
      else if ( i .le. n1 .and. n1 .lt. j ) then
        ij = (ml+mu+1)*n1+(j-n1-1)*n1+i
c
c  The A3 and A4 blocks of the matrix.
c
      else if ( n1 .lt. i ) then
        ij = (ml+mu+1)*n1+n2*n1+(j-1)*n2+(i-n1)
      end if

      value = a(ij)

      return
      end
      subroutine r8cbb_indicator ( n1, n2, ml, mu, a )

c*********************************************************************72
c
cc R8CBB_INDICATOR sets up an R8CBB indicator matrix.
c
c  Discussion:
c
c    The "indicator matrix" simply has a value like I*10+J at every
c    entry of a dense matrix, or at every entry of a compressed storage
c    matrix for which storage is allocated. 
c
c    The R8CBB storage format is for a compressed border banded matrix.  
c    Such a matrix has the logical form:
c
c      A1 | A2
c      ---+---
c      A3 | A4
c
c    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
c    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
c    respectively.  
c
c    The R8CBB format is the same as the R8BB format, except that the banded
c    matrix A1 is stored in compressed band form rather than standard
c    banded form.  In other words, we do not include the extra room
c    set aside for fill in during pivoting.
c
c    A should be defined as a vector.  The user must then store
c    the entries of the four blocks of the matrix into the vector A.
c    Each block is stored by columns.
c
c    A1, the banded portion of the matrix, is stored in
c    the first (ML+MU+1)*N1 entries of A, using the obvious variant
c    of the LINPACK general band format.
c
c    The following formulas should be used to determine how to store
c    the entry corresponding to row I and column J in the original matrix:
c
c    Entries of A1:
c
c      1 .le. I .le. N1, 1 .le. J .le. N1, (J-I) .le. MU and (I-J) .le. ML.
c
c      Store the I, J entry into location
c      (I-J+MU+1)+(J-1)*(ML+MU+1).
c
c    Entries of A2:
c
c      1 .le. I .le. N1, N1+1 .le. J .le. N1+N2.
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+(J-N1-1)*N1+I.
c
c    Entries of A3:
c
c      N1+1 .le. I .le. N1+N2, 1 .le. J .le. N1.
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c
c    Entries of A4:
c
c      N1+1 .le. I .le. N1+N2, N1+1 .le. J .le. N1+N2
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c      (same formula used for A3).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, the order of the banded and dense 
c    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative and no greater than N1-1.
c
c    Output, double precision A((ML+MU+1)*N1+2*N1*N2+N2*N2), the R8CBB matrix.
c
      implicit none

      integer ml
      integer mu
      integer n1
      integer n2

      double precision a((ml+mu+1)*n1+2*n1*n2+n2*n2)
      integer base
      integer fac
      integer i
      integer i4_log_10
      integer j
      integer row

      fac = 10 ** ( i4_log_10 ( n1 + n2 ) + 1 )
c
c  Set the banded matrix A1.
c
      do j = 1, n1
        do row = 1, ml + mu + 1
          i = row + j - mu - 1
          if ( 1 .le. i .and. i .le. n1 ) then
            a(row+(j-1)*(ml+mu+1)) = dble ( fac * i + j )
          else
            a(row+(j-1)*(ml+mu+1)) = 0.0D+00
          end if
        end do
      end do
c
c  Set the N1 by N2 rectangular strip A2.
c
      base = ( ml + mu + 1 ) * n1

      do i = 1, n1
        do j = n1 + 1, n1 + n2
          a(base + i + (j-n1-1)*n1 ) = dble ( fac * i + j )
        end do
      end do
c
c  Set the N2 by N1 rectangular strip A3.
c
      base = ( ml + mu + 1 ) * n1 + n1 * n2

      do i = n1 + 1, n1 + n2
        do j = 1, n1    
          a(base + i-n1 + (j-1)*n2 ) = dble ( fac * i + j )
        end do
      end do
c
c  Set the N2 by N2 square A4.
c
      base = ( ml + mu + 1 ) * n1 + n1 * n2 + n2 * n1

      do i = n1 + 1, n1 + n2
        do j = n1 + 1, n1 + n2
          a(base + i-n1 + (j-n1-1)*n2 ) = dble ( fac * i + j )
        end do
      end do

      return
      end
      subroutine r8cbb_mxv ( n1, n2, ml, mu, a, x, b )

c*********************************************************************72
c
cc R8CBB_MXV multiplies an R8CBB matrix by an R8VEC.
c
c  Discussion:
c
c    The R8CBB storage format is for a compressed border banded matrix.  
c    Such a matrix has the logical form:
c
c      A1 | A2
c      ---+---
c      A3 | A4
c
c    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
c    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
c    respectively.  
c
c    The R8CBB format is the same as the R8BB format, except that the banded
c    matrix A1 is stored in compressed band form rather than standard
c    banded form.  In other words, we do not include the extra room
c    set aside for fill in during pivoting.
c
c    A should be defined as a vector.  The user must then store
c    the entries of the four blocks of the matrix into the vector A.
c    Each block is stored by columns.
c
c    A1, the banded portion of the matrix, is stored in
c    the first (ML+MU+1)*N1 entries of A, using the obvious variant
c    of the LINPACK general band format.
c
c    The following formulas should be used to determine how to store
c    the entry corresponding to row I and column J in the original matrix:
c
c    Entries of A1:
c
c      1 .le. I .le. N1, 1 .le. J .le. N1, (J-I) .le. MU and (I-J) .le. ML.
c
c      Store the I, J entry into location
c      (I-J+MU+1)+(J-1)*(ML+MU+1).
c
c    Entries of A2:
c
c      1 .le. I .le. N1, N1+1 .le. J .le. N1+N2.
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+(J-N1-1)*N1+I.
c
c    Entries of A3:
c
c      N1+1 .le. I .le. N1+N2, 1 .le. J .le. N1.
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c
c    Entries of A4:
c
c      N1+1 .le. I .le. N1+N2, N1+1 .le. J .le. N1+N2
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c      (same formula used for A3).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than N1-1.
c
c    Input, integer N1, N2, the order of the banded and dense 
c    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
c
c    Input, double precision A((ML+MU+1)*N1 + 2*N1*N2 + N2*N2),
c    the R8CBB matrix.
c
c    Input, double precision X(N1+N2), the vector to be multiplied by A.
c
c    Output, double precision B(N1+N2), the result of multiplying A by X.
c
      implicit none

      integer ml
      integer mu
      integer n1
      integer n2

      double precision a((ml+mu+1)*n1+2*n1*n2+n2*n2)
      double precision b(n1+n2)
      integer i
      integer ihi
      integer ij
      integer ilo
      integer j
      double precision x(n1+n2)
c
c  Set B to zero.
c
      do i = 1, n1 + n2
        b(i) = 0.0D+00
      end do
c
c  Multiply by A1.
c
      do j = 1, n1
        ilo = max ( 1, j-mu )
        ihi = min ( n1, j+ml )
        ij = (j-1)*(ml+mu+1)-j+mu+1
        do i = ilo, ihi
          b(i) = b(i) + a(ij+i) * x(j)
        end do
      end do
c
c  Multiply by A2.
c
      do j = n1+1, n1+n2
        ij = (ml+mu+1)*n1+(j-n1-1)*n1
        do i = 1, n1
          b(i) = b(i) + a(ij+i) * x(j)
        end do
      end do
c
c  Multiply by A3 and A4.
c
      do j = 1, n1+n2
        ij = (ml+mu+1)*n1+n1*n2+(j-1)*n2-n1
        do i = n1 + 1, n1 + n2
          b(i) = b(i) + a(ij+i) * x(j)
        end do
      end do

      return
      end
      subroutine r8cbb_print ( n1, n2, ml, mu, a, title )

c*********************************************************************72
c
cc R8CBB_PRINT prints an R8CBB matrix.
c
c  Discussion:
c
c    The R8CBB storage format is for a compressed border banded matrix.  
c    Such a matrix has the logical form:
c
c      A1 | A2
c      ---+---
c      A3 | A4
c
c    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
c    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
c    respectively.  
c
c    The R8CBB format is the same as the R8BB format, except that the banded
c    matrix A1 is stored in compressed band form rather than standard
c    banded form.  In other words, we do not include the extra room
c    set aside for fill in during pivoting.
c
c    A should be defined as a vector.  The user must then store
c    the entries of the four blocks of the matrix into the vector A.
c    Each block is stored by columns.
c
c    A1, the banded portion of the matrix, is stored in
c    the first (ML+MU+1)*N1 entries of A, using the obvious variant
c    of the LINPACK general band format.
c
c    The following formulas should be used to determine how to store
c    the entry corresponding to row I and column J in the original matrix:
c
c    Entries of A1:
c
c      1 .le. I .le. N1, 1 .le. J .le. N1, (J-I) .le. MU and (I-J) .le. ML.
c
c      Store the I, J entry into location
c      (I-J+MU+1)+(J-1)*(ML+MU+1).
c
c    Entries of A2:
c
c      1 .le. I .le. N1, N1+1 .le. J .le. N1+N2.
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+(J-N1-1)*N1+I.
c
c    Entries of A3:
c
c      N1+1 .le. I .le. N1+N2, 1 .le. J .le. N1.
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c
c    Entries of A4:
c
c      N1+1 .le. I .le. N1+N2, N1+1 .le. J .le. N1+N2
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c      (same formula used for A3).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, the order of the banded and dense 
c    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than N1-1.
c
c    Input, double precision A((ML+MU+1)*N1+2*N1*N2+N2*N2), the R8CBB matrix.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer ml
      integer mu
      integer n1
      integer n2

      double precision a((ml+mu+1)*n1+2*n1*n2+n2*n2)
      character * ( * )  title

      call r8cbb_print_some ( n1, n2, ml, mu, a, 1, 1, n1+n2, n1+n2, 
     &  title )

      return
      end
      subroutine r8cbb_print_some ( n1, n2, ml, mu, a, ilo, jlo, ihi, 
     &  jhi, title )

c*********************************************************************72
c
cc R8CBB_PRINT_SOME prints some of an R8CBB matrix.
c
c  Discussion:
c
c    The R8CBB storage format is for a compressed border banded matrix.  
c    Such a matrix has the logical form:
c
c      A1 | A2
c      ---+---
c      A3 | A4
c
c    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
c    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
c    respectively.  
c
c    The R8CBB format is the same as the R8BB format, except that the banded
c    matrix A1 is stored in compressed band form rather than standard
c    banded form.  In other words, we do not include the extra room
c    set aside for fill in during pivoting.
c
c    A should be defined as a vector.  The user must then store
c    the entries of the four blocks of the matrix into the vector A.
c    Each block is stored by columns.
c
c    A1, the banded portion of the matrix, is stored in
c    the first (ML+MU+1)*N1 entries of A, using the obvious variant
c    of the LINPACK general band format.
c
c    The following formulas should be used to determine how to store
c    the entry corresponding to row I and column J in the original matrix:
c
c    Entries of A1:
c
c      1 .le. I .le. N1, 1 .le. J .le. N1, (J-I) .le. MU and (I-J) .le. ML.
c
c      Store the I, J entry into location
c      (I-J+MU+1)+(J-1)*(ML+MU+1).
c
c    Entries of A2:
c
c      1 .le. I .le. N1, N1+1 .le. J .le. N1+N2.
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+(J-N1-1)*N1+I.
c
c    Entries of A3:
c
c      N1+1 .le. I .le. N1+N2, 1 .le. J .le. N1.
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c
c    Entries of A4:
c
c      N1+1 .le. I .le. N1+N2, N1+1 .le. J .le. N1+N2
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c      (same formula used for A3).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, the order of the banded and dense 
c    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than N1-1.
c
c    Input, double precision A((ML+MU+1)*N1+2*N1*N2+N2*N2), the R8CBB matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 5 )
      integer ml
      integer mu
      integer n1
      integer n2

      double precision a((ml+mu+1)*n1+2*n1*n2+n2*n2)
      double precision aij
      character ( len = 14 ) ctemp(incx)
      logical r8_is_int
      integer i
      integer i2hi
      integer i2lo
      integer ihi
      integer ij
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
        j2hi = min ( j2hi, n1+n2 )
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
        i2hi = min ( ihi, n1+n2 )

        do i = i2lo, i2hi
c
c  Print out (up to) 5 entries in row I, that lie in the current strip.
c
          do j2 = 1, inc

            j = j2lo - 1 + j2

            aij = 0.0D+00

            if ( i .le. n1 .and. j .le. n1 ) then
              if ( j - i .le. mu .and. i - j .le. ml ) then
                ij = (i-j+mu+1)+(j-1)*(ml+mu+1)
                aij = a(ij)
              end if
            else if ( i .le. n1 .and. n1 .lt. j ) then
              ij = (ml+mu+1)*n1+(j-n1-1)*n1+i
              aij = a(ij)
            else if ( n1 .lt. i ) then
              ij = (ml+mu+1)*n1+n2*n1+(j-1)*n2+(i-n1)
              aij = a(ij)
            end if

            if ( r8_is_int ( aij ) ) then
              write ( ctemp(j2), '(f8.0,6x)' ) aij
            else
              write ( ctemp(j2), '(g14.6)' ) aij
            end if

          end do

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine r8cbb_random ( n1, n2, ml, mu, seed, a )

c*********************************************************************72
c
cc R8CBB_RANDOM randomizes an R8CBB matrix.
c
c  Discussion:
c
c    The R8CBB storage format is for a compressed border banded matrix.  
c    Such a matrix has the logical form:
c
c      A1 | A2
c      ---+---
c      A3 | A4
c
c    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
c    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
c    respectively.  
c
c    The R8CBB format is the same as the R8BB format, except that the banded
c    matrix A1 is stored in compressed band form rather than standard
c    banded form.  In other words, we do not include the extra room
c    set aside for fill in during pivoting.
c
c    A should be defined as a vector.  The user must then store
c    the entries of the four blocks of the matrix into the vector A.
c    Each block is stored by columns.
c
c    A1, the banded portion of the matrix, is stored in
c    the first (ML+MU+1)*N1 entries of A, using the obvious variant
c    of the LINPACK general band format.
c
c    The following formulas should be used to determine how to store
c    the entry corresponding to row I and column J in the original matrix:
c
c    Entries of A1:
c
c      1 .le. I .le. N1, 1 .le. J .le. N1, (J-I) .le. MU and (I-J) .le. ML.
c
c      Store the I, J entry into location
c      (I-J+MU+1)+(J-1)*(ML+MU+1).
c
c    Entries of A2:
c
c      1 .le. I .le. N1, N1+1 .le. J .le. N1+N2.
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+(J-N1-1)*N1+I.
c
c    Entries of A3:
c
c      N1+1 .le. I .le. N1+N2, 1 .le. J .le. N1.
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c
c    Entries of A4:
c
c      N1+1 .le. I .le. N1+N2, N1+1 .le. J .le. N1+N2
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c      (same formula used for A3).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, the order of the banded and dense 
c    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative and no greater than N1-1.
c
c    Input/output, integer SEED, a seed for the random 
c    number generator.
c
c    Output, double precision A((ML+MU+1)*N1 + 2*N1*N2 + N2*N2), 
c    the R8CBB matrix.
c
      implicit none

      integer ml
      integer mu
      integer n1
      integer n2

      double precision a((ml+mu+1)*n1+2*n1*n2+n2*n2)
      double precision r8_uniform_01
      integer i
      integer ilo
      integer j
      double precision r
      integer row
      integer seed
c
c  Randomize the banded matrix A1.
c  We still believe that the "junk" entries should be set to 0.
c
      do j = 1, n1
        do row = 1, ml+mu+1
          i = row + j - mu - 1
          if ( 1 .le. i .and. i .le. n1 ) then
            r = r8_uniform_01 ( seed )
          else
            r = 0.0D+00
          end if
          a(row+(j-1)*(ml+mu+1)) = r
        end do
      end do
c
c  Randomize the rectangular strips A2+A3+A4.
c
      ilo = (ml+mu+1) * n1 + 1

      call r8vec_uniform_01 ( n1*n2+n2*n1+n2*n2, seed, a(ilo) )

      return
      end
      subroutine r8cbb_set ( n1, n2, ml, mu, a, i, j, value )

c*********************************************************************72
c
cc R8CBB_SET sets the value of an entry in an R8CBB matrix.
c
c  Discussion:
c
c    The R8CBB storage format is for a compressed border banded matrix.  
c    Such a matrix has the logical form:
c
c      A1 | A2
c      ---+---
c      A3 | A4
c
c    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
c    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
c    respectively.  
c
c    The R8CBB format is the same as the R8BB format, except that the banded
c    matrix A1 is stored in compressed band form rather than standard
c    banded form.  In other words, we do not include the extra room
c    set aside for fill in during pivoting.
c
c    A should be defined as a vector.  The user must then store
c    the entries of the four blocks of the matrix into the vector A.
c    Each block is stored by columns.
c
c    A1, the banded portion of the matrix, is stored in
c    the first (ML+MU+1)*N1 entries of A, using the obvious variant
c    of the LINPACK general band format.
c
c    The following formulas should be used to determine how to store
c    the entry corresponding to row I and column J in the original matrix:
c
c    Entries of A1:
c
c      1 .le. I .le. N1, 1 .le. J .le. N1, (J-I) .le. MU and (I-J) .le. ML.
c
c      Store the I, J entry into location
c      (I-J+MU+1)+(J-1)*(ML+MU+1).
c
c    Entries of A2:
c
c      1 .le. I .le. N1, N1+1 .le. J .le. N1+N2.
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+(J-N1-1)*N1+I.
c
c    Entries of A3:
c
c      N1+1 .le. I .le. N1+N2, 1 .le. J .le. N1.
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c
c    Entries of A4:
c
c      N1+1 .le. I .le. N1+N2, N1+1 .le. J .le. N1+N2
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c      (same formula used for A3).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, the order of the banded and dense 
c    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than N1-1.
c
c    Input/output, double precision A((ML+MU+1)*N1 + 2*N1*N2 + N2*N2),
c    the R8CBB matrix.
c
c    Input, integer I, J, the row and column of the entry to set.
c
c    Input, double precision VALUE, the value to be assigned to the
c    (I,J) entry.
c
      implicit none

      integer ml
      integer mu
      integer n1
      integer n2

      double precision a((ml+mu+1)*n1+2*n1*n2+n2*n2)
      integer i
      integer ij
      integer j
      double precision value
c
c  Check for I or J out of bounds.
c
      if ( i .le. 0 .or. n1+n2 .lt. i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8CBB_SET - Fatal error!'
        write ( *, '(a,i8)' ) 
     &    '  Illegal input value of row index I = ', i
        stop
      end if

      if ( j .le. 0 .or. n1+n2 .lt. j ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8CBB_SET - Fatal error!'
        write ( *, '(a,i8)' ) 
     &    '  Illegal input value of column index J = ', j
        stop
      end if
c
c  The A1 block of the matrix.
c
c  Check for out of band problems.
c
      if ( i .le. n1 .and. j .le. n1 ) then
        if ( mu .lt. (j-i) .or. ml .lt. (i-j) ) then
          if ( value .ne. 0.0D+00 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'R8CBB_SET - Warning!'
            write ( *, '(a,i8,a,i8,a)' ) 
     &        '  Unable to set entry (', i, ',', j, ').'
          end if
          return
        else
          ij = (i-j+mu+1)+(j-1)*(ml+mu+1)
        end if
c
c  The A2 block of the matrix:
c
      else if ( i .le. n1 .and. n1 .lt. j ) then
        ij = (ml+mu+1)*n1 + (j-n1-1)*n1 + i
c
c  The A3 and A4 blocks of the matrix.
c
      else if ( n1 .lt. i ) then
        ij = (ml+mu+1)*n1 + n2*n1 + (j-1)*n2 + (i-n1)
      end if

      a(ij) = value

      return
      end
      subroutine r8cbb_sl ( n1, n2, ml, mu, a_lu, b )

c*********************************************************************72
c
cc R8CBB_SL solves an R8CBB system factored by R8CBB_FA.
c
c  Discussion:
c
c    The R8CBB storage format is for a compressed border banded matrix.  
c    Such a matrix has the logical form:
c
c      A1 | A2
c      ---+---
c      A3 | A4
c
c    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
c    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
c    respectively.  
c
c    The R8CBB format is the same as the R8BB format, except that the banded
c    matrix A1 is stored in compressed band form rather than standard
c    banded form.  In other words, we do not include the extra room
c    set aside for fill in during pivoting.
c
c    A should be defined as a vector.  The user must then store
c    the entries of the four blocks of the matrix into the vector A.
c    Each block is stored by columns.
c
c    A1, the banded portion of the matrix, is stored in
c    the first (ML+MU+1)*N1 entries of A, using the obvious variant
c    of the LINPACK general band format.
c
c    The following formulas should be used to determine how to store
c    the entry corresponding to row I and column J in the original matrix:
c
c    Entries of A1:
c
c      1 .le. I .le. N1, 1 .le. J .le. N1, (J-I) .le. MU and (I-J) .le. ML.
c
c      Store the I, J entry into location
c      (I-J+MU+1)+(J-1)*(ML+MU+1).
c
c    Entries of A2:
c
c      1 .le. I .le. N1, N1+1 .le. J .le. N1+N2.
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+(J-N1-1)*N1+I.
c
c    Entries of A3:
c
c      N1+1 .le. I .le. N1+N2, 1 .le. J .le. N1.
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c
c    Entries of A4:
c
c      N1+1 .le. I .le. N1+N2, N1+1 .le. J .le. N1+N2
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c      (same formula used for A3).
c
c
c    The linear system A * x = b is decomposable into the block system:
c
c      ( A1 A2 ) * (X1) = (B1)
c      ( A3 A4 )   (X2)   (B2)
c
c    where A1 is a (usually big) banded square matrix, A2 and A3 are
c    column and row strips which may be nonzero, and A4 is a dense
c    square matrix.
c
c    All the arguments except B are input quantities only, which are
c    not changed by the routine.  They should have exactly the same values
c    they had on exit from R8CBB_FA.
c
c    If more than one right hand side is to be solved, with the same
c    matrix, R8CBB_SL should be called repeatedly.  However, R8CBB_FA only
c    needs to be called once to create the factorization.
c
c    See the documentation of R8CBB_FA for details on the matrix storage.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, the order of the banded and dense 
c    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than N1-1.
c
c    Input, double precision A_LU( (ML+MU+1)*N1 + 2*N1*N2 + N2*N2).
c    the LU factors from R8CBB_FA.
c
c    Input/output, double precision B(N1+N2).
c    On input, B contains the right hand side of the linear system.
c    On output, B contains the solution.
c
      implicit none

      integer ml
      integer mu
      integer n1
      integer n2

      double precision a_lu((ml+mu+1)*n1+2*n1*n2+n2*n2)
      double precision b(n1+n2)
      integer i
      integer ij
      integer j
      integer job
      integer nband

      nband = (ml+mu+1)*n1
c
c  Set B1 := inverse(A1) * B1.
c
      if ( 0 .lt. n1 ) then
        job = 0
        call r8cb_np_sl ( n1, ml, mu, a_lu, b, job )
      end if
c
c  Modify the right hand side of the second linear subsystem.
c  Replace B2 by B2-A3*B1.
c
      do j = 1, n1
        do i = 1, n2
          ij = nband + n1*n2 + (j-1)*n2 + i
          b(n1+i) = b(n1+i) - a_lu(ij) * b(j)
        end do
      end do
c
c  Solve A4*B2 = B2.
c
      if ( 0 .lt. n2 ) then
        job = 0
        call r8ge_np_sl ( n2, a_lu(nband+2*n1*n2+1), b(n1+1), job )
      end if
c
c  Modify the first subsolution.
c  Set B1 = B1+A2*B2.
c
      do i = 1, n1
        do j = 1, n2
          ij = nband + (j-1)*n1 + i
          b(i) = b(i) + a_lu(ij) * b(n1+j)
        end do
      end do

      return
      end
      subroutine r8cbb_to_r8ge ( n1, n2, ml, mu, a, b )

c*********************************************************************72
c
cc R8CBB_TO_R8GE copies an R8CBB matrix to an R8GE matrix.
c
c  Discussion:
c
c    The R8CBB storage format is for a compressed border banded matrix.  
c    Such a matrix has the logical form:
c
c      A1 | A2
c      ---+---
c      A3 | A4
c
c    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
c    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
c    respectively.  
c
c    The R8CBB format is the same as the R8BB format, except that the banded
c    matrix A1 is stored in compressed band form rather than standard
c    banded form.  In other words, we do not include the extra room
c    set aside for fill in during pivoting.
c
c    A should be defined as a vector.  The user must then store
c    the entries of the four blocks of the matrix into the vector A.
c    Each block is stored by columns.
c
c    A1, the banded portion of the matrix, is stored in
c    the first (ML+MU+1)*N1 entries of A, using the obvious variant
c    of the LINPACK general band format.
c
c    The following formulas should be used to determine how to store
c    the entry corresponding to row I and column J in the original matrix:
c
c    Entries of A1:
c
c      1 .le. I .le. N1, 1 .le. J .le. N1, (J-I) .le. MU and (I-J) .le. ML.
c
c      Store the I, J entry into location
c      (I-J+MU+1)+(J-1)*(ML+MU+1).
c
c    Entries of A2:
c
c      1 .le. I .le. N1, N1+1 .le. J .le. N1+N2.
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+(J-N1-1)*N1+I.
c
c    Entries of A3:
c
c      N1+1 .le. I .le. N1+N2, 1 .le. J .le. N1.
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c
c    Entries of A4:
c
c      N1+1 .le. I .le. N1+N2, N1+1 .le. J .le. N1+N2
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c      (same formula used for A3).
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, the order of the banded and dense 
c    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than N1-1.
c
c    Input, double precision A((ML+MU+1)*N1+2*N1*N2+N2*N2), the R8CBB matrix.
c
c    Output, double precision B(N1+N2,N1+N2), the R8GE matrix.
c
      implicit none

      integer ml
      integer mu
      integer n1
      integer n2

      double precision a((ml+mu+1)*n1+2*n1*n2+n2*n2)
      double precision b(n1+n2,n1+n2)
      integer i
      integer ij
      integer j

      do i = 1, n1
        do j = 1, n1

          if ( mu+ml .lt. (j-i) .or. ml .lt. (i-j) ) then
            b(i,j) = 0.0D+00
          else
            ij = (i-j+mu+1)+(j-1)*(ml+mu+1)
            b(i,j) = a(ij)
          end if

        end do
      end do

      do i = 1, n1
        do j = n1+1, n2
          ij = (ml+mu+1)*n1+(j-n1-1)*n1+i
          b(i,j) = a(ij)
        end do
      end do

      do i = n1+1, n2
        do j = 1, n1+n2
          ij = (ml+mu+1)*n1+n2*n1+(j-1)*n2+(i-n1)
          b(i,j) = a(ij)
        end do
      end do

      return
      end
      subroutine r8cbb_vxm ( n1, n2, ml, mu, a, x, b )

c*********************************************************************72
c
cc R8CBB_VXM multiplies an R8VEC by an R8CBB matrix.
c
c  Discussion:
c
c    The R8CBB storage format is for a compressed border banded matrix.  
c    Such a matrix has the logical form:
c
c      A1 | A2
c      ---+---
c      A3 | A4
c
c    with A1 a (usually large) N1 by N1 banded matrix, while A2, A3 and A4
c    are dense rectangular matrices of orders N1 by N2, N2 by N1, and N2 by N2,
c    respectively.  
c
c    The R8CBB format is the same as the R8BB format, except that the banded
c    matrix A1 is stored in compressed band form rather than standard
c    banded form.  In other words, we do not include the extra room
c    set aside for fill in during pivoting.
c
c    A should be defined as a vector.  The user must then store
c    the entries of the four blocks of the matrix into the vector A.
c    Each block is stored by columns.
c
c    A1, the banded portion of the matrix, is stored in
c    the first (ML+MU+1)*N1 entries of A, using the obvious variant
c    of the LINPACK general band format.
c
c    The following formulas should be used to determine how to store
c    the entry corresponding to row I and column J in the original matrix:
c
c    Entries of A1:
c
c      1 .le. I .le. N1, 1 .le. J .le. N1, (J-I) .le. MU and (I-J) .le. ML.
c
c      Store the I, J entry into location
c      (I-J+MU+1)+(J-1)*(ML+MU+1).
c
c    Entries of A2:
c
c      1 .le. I .le. N1, N1+1 .le. J .le. N1+N2.
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+(J-N1-1)*N1+I.
c
c    Entries of A3:
c
c      N1+1 .le. I .le. N1+N2, 1 .le. J .le. N1.
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c
c    Entries of A4:
c
c      N1+1 .le. I .le. N1+N2, N1+1 .le. J .le. N1+N2
c
c      Store the I, J entry into location
c      (ML+MU+1)*N1+N1*N2+(J-1)*N2+(I-N1).
c      (same formula used for A3).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 July 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than N1-1.
c
c    Input, integer N1, N2, the order of the banded and dense 
c    blocks.  N1 and N2 must be nonnegative, and at least one must be positive.
c
c    Input, double precision A((ML+MU+1)*N1 + 2*N1*N2 + N2*N2),
c    the R8CBB matrix.
c
c    Input, double precision X(N1+N2), the vector to multiply the matrix.
c
c    Output, double precision B(N1+N2), the product X * A.
c
      implicit none

      integer ml
      integer mu
      integer n1
      integer n2

      double precision a((ml+mu+1)*n1+2*n1*n2+n2*n2)
      double precision b(n1+n2)
      integer i
      integer ihi
      integer ij
      integer ilo
      integer j
      double precision x(n1+n2)
c
c  Set B to zero.
c
      do i = 1, n1 + n2
        b(i) = 0.0D+00
      end do
c
c  Multiply by A1.
c
      do j = 1, n1
        ilo = max ( 1, j-mu )
        ihi = min ( n1, j+ml )
        ij = (j-1)*(ml+mu+1)-j+mu+1
        do i = ilo, ihi
          b(j) = b(j) + x(i) * a(ij+i)
        end do
      end do
c
c  Multiply by A2.
c
      do j = n1+1, n1+n2
        ij = (ml+mu+1)*n1+(j-n1-1)*n1
        do i = 1, n1
          b(j) = b(j) + x(i) * a(ij+i)
        end do
      end do
c
c  Multiply by A3 and A4.
c
      do j = 1, n1+n2
        ij = (ml+mu+1)*n1+n1*n2+(j-1)*n2-n1
        do i = n1+1, n1+n2
          b(j) = b(j) + x(i) * a(ij+i)
        end do
      end do

      return
      end
      subroutine r8cc_get ( m, n, nz_num, col, row, a, i, j, aij )

c*********************************************************************72
c
cc R8CC_GET gets a value of an R8CC matrix.
c
c  Discussion:
c
c    It is legal to request entries of the matrix for which no storage
c    was set aside.  In that case, a zero value will be returned.
c
c    The R8CC format is the double precision sparse compressed column
c    format.  Associated with this format, we have an M by N matrix
c    with NZ_NUM nonzero entries.  We construct the column pointer
c    vector COL of length N+1, such that entries of column J will be
c    stored in positions COL(J) through COL(J+1)-1.  This indexing
c    refers to both the ROW and A vectors, which store the row indices
c    and the values of the nonzero entries.  The entries of the
c    ROW vector corresponding to each column are assumed to be
c    ascending sorted.
c
c    The R8CC format is equivalent to the MATLAB "sparse" format,
c    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    31 August 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Iain Duff, Roger Grimes, John Lewis,
c    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
c    October 1992
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c
c    Input, integer N, the number of columns of the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero entries.
c
c    Input, integer COL(N+1), indicate where each column's data
c    begins.
c
c    Input, integer ROW(NZ_NUM), the row indices.
c
c    Input, double precision A(NZ_NUM), the nonzero entries.
c
c    Input, integer I, J, the indices of the value to retrieve.
c
c    Output, double precision AIJ, the value of A(I,J).
c
      implicit none

      integer n
      integer nz_num

      double precision a(nz_num)
      double precision aij
      integer col(n+1)
      integer i
      integer j
      integer k
      integer m
      integer row(nz_num)
c
c  Seek sparse index K corresponding to full index (I,J).
c
      call r8cc_ijk ( m, n, nz_num, col, row, i, j, k )
c
c  If no K was found, then be merciful, and simply return 0.
c
      if ( k .eq. -1 ) then
        aij = 0.0D+00
      else
        aij = a(k)
      end if

      return
      end
      subroutine r8cc_ijk ( m, n, nz_num, col, row, i, j, k )

c*********************************************************************72
c
cc R8CC_IJK seeks K, the sparse index of (I,J), the full index of an R8CC matrix.
c
c  Discussion:
c
c    The R8CC format is the double precision sparse compressed column
c    format.  Associated with this format, we have an M by N matrix
c    with NZ_NUM nonzero entries.  We construct the column pointer
c    vector COL of length N+1, such that entries of column J will be
c    stored in positions COL(J) through COL(J+1)-1.  This indexing
c    refers to both the ROW and A vectors, which store the row indices
c    and the values of the nonzero entries.  The entries of the
c    ROW vector corresponding to each column are assumed to be
c    ascending sorted.
c
c    The R8CC format is equivalent to the MATLAB "sparse" format,
c    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    31 August 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Iain Duff, Roger Grimes, John Lewis,
c    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
c    October 1992
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c
c    Input, integer N, the number of columns of the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero entries.
c
c    Input, integer COL(N+1), indicate where each column's 
c    data begins.
c
c    Input, integer ROW(NZ_NUM), the row indices.
c
c    Input, integer I, J, the indices of the value to retrieve.
c
c    Output, integer K, the index of the sparse matrix in which 
c    entry (I,J) is stored, or -1 if no such entry exists.
c
      implicit none

      integer n
      integer nz_num

      integer col(n+1)
      integer i
      integer j
      integer k
      integer k1
      integer k2
      integer m
      integer row(nz_num)
c
c  Determine the part of ROW containing row indices of entries
c  in column J.
c
      k1 = col(j)
      k2 = col(j+1)-1
c
c  Seek the location K for which ROW(K) = I.
c  
      call i4vec_search_binary_a ( k2+1-k1, row(k1:k2), i, k )

      if ( k .ne. -1 ) then
        k = k + k1 - 1
      end if

      return
      end
      subroutine r8cc_inc ( m, n, nz_num, col, row, a, i, j, aij )

c*********************************************************************72
c
cc R8CC_INC increments a value of an R8CC matrix.
c
c  Discussion:
c
c    The R8CC format is the double precision sparse compressed column
c    format.  Associated with this format, we have an M by N matrix
c    with NZ_NUM nonzero entries.  We construct the column pointer
c    vector COL of length N+1, such that entries of column J will be
c    stored in positions COL(J) through COL(J+1)-1.  This indexing
c    refers to both the ROW and A vectors, which store the row indices
c    and the values of the nonzero entries.  The entries of the
c    ROW vector corresponding to each column are assumed to be
c    ascending sorted.
c
c    The R8CC format is equivalent to the MATLAB "sparse" format,
c    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    31 August 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Iain Duff, Roger Grimes, John Lewis,
c    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
c    October 1992
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c
c    Input, integer N, the number of columns of the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero entries.
c
c    Input, integer COL(N+1), indicate where each column's 
c    data begins.
c
c    Input, integer ROW(NZ_NUM), the row indices.
c
c    Input/output, double precision A(NZ_NUM), the nonzero entries.
c    On output, entry (I,J) has been incremented.
c
c    Input, integer I, J, the indices of the value to retrieve.
c
c    Input, double precision AIJ, the value to be added to A(I,J).
c
      implicit none

      integer n
      integer nz_num

      double precision a(nz_num)
      double precision aij
      integer col(n+1)
      integer i
      integer j
      integer k
      integer m
      integer row(nz_num)
c
c  Seek sparse index K corresponding to full index (I,J).
c
      call r8cc_ijk ( m, n, nz_num, col, row, i, j, k )
c
c  If no K was found, we fail.
c
      if ( k .eq. -1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8CC_INC - Fatal error!'
        write ( *, '(a)' ) '  R8CC_IJK could not find the entry.'
        write ( *, '(a,i8)' ) '  Row I = ', i
        write ( *, '(a,i8)' ) '  Col J = ', j
        stop
      end if

      a(k) = a(k) + aij

      return
      end
      subroutine r8cc_indicator ( m, n, nz_num, col, row, a )

c*********************************************************************72
c
cc R8CC_INDICATOR sets up an R8CC indicator matrix.
c
c  Discussion:
c
c    The "indicator matrix" simply has a value like I*10+J at every
c    entry of a dense matrix, or at every entry of a compressed storage
c    matrix for which storage is allocated. 
c
c    The R8CC format is the double precision sparse compressed column
c    format.  Associated with this format, we have an M by N matrix
c    with NZ_NUM nonzero entries.  We construct the column pointer
c    vector COL of length N+1, such that entries of column J will be
c    stored in positions COL(J) through COL(J+1)-1.  This indexing
c    refers to both the ROW and A vectors, which store the row indices
c    and the values of the nonzero entries.  The entries of the
c    ROW vector corresponding to each column are assumed to be
c    ascending sorted.
c
c    The R8CC format is equivalent to the MATLAB "sparse" format,
c    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    31 August 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Iain Duff, Roger Grimes, John Lewis,
c    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
c    October 1992
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c
c    Input, integer N, the number of columns of the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero elements in A.
c
c    Input, integer COL(N+1), points to the first element of 
c    each column.
c
c    Input, integer ROW(NZ_NUM), contains the row indices 
c    of the elements.
c
c    Output, double precision A(NZ_NUM), the R8CC matrix.
c
      implicit none

      integer n
      integer nz_num

      double precision a(nz_num)
      integer col(n+1)
      integer fac
      integer i
      integer i4_log_10
      integer j
      integer k
      integer m
      integer row(nz_num)

      fac = 10 ** ( i4_log_10 ( n ) + 1 )

      do j = 1, n
        do k = col(j), col(j+1) - 1
          i = row(k)
          a(k) = dble ( fac * i + j )
        end do
      end do

      return
      end
      subroutine r8cc_kij ( m, n, nz_num, col, row, k, i, j )

c*********************************************************************72
c
cc R8CC_KIJ seeks (I,J), the full index of K, the sparse index of an R8CC matrix.
c
c  Discussion:
c
c    The R8CC format is the double precision sparse compressed column
c    format.  Associated with this format, we have an M by N matrix
c    with NZ_NUM nonzero entries.  We construct the column pointer
c    vector COL of length N+1, such that entries of column J will be
c    stored in positions COL(J) through COL(J+1)-1.  This indexing
c    refers to both the ROW and A vectors, which store the row indices
c    and the values of the nonzero entries.  The entries of the
c    ROW vector corresponding to each column are assumed to be
c    ascending sorted.
c
c    The R8CC format is equivalent to the MATLAB "sparse" format,
c    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    31 August 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Iain Duff, Roger Grimes, John Lewis,
c    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
c    October 1992
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c
c    Input, integer N, the number of columns of the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero entries.
c
c    Input, integer COL(N+1), indicate where each column's data 
c    begins.
c
c    Input, integer ROW(NZ_NUM), the row indices.
c
c    Input, integer K, the sparse index of an entry of the matrix.
c    1 .le. K .le. NZ_NUM.
c
c    Output, integer I, J, the full indices corresponding to the 
c    sparse index K.
c
      implicit none

      integer n
      integer nz_num

      integer col(n+1)
      integer i
      integer j
      integer jj
      integer k
      integer k1
      integer k2
      integer m
      integer row(nz_num)

      i = -1
      j = -1

      if ( k .lt. 1 .or. nz_num .lt. k ) then
        return
      end if
c
c  The row index is easy.
c
      i = row(k)
c
c  Determine the column by bracketing in COL.
c
      do jj = 1, n
        k1 = col(jj)
        k2 = col(jj+1)-1
        if ( k1 .le. k .and. k .le. k2 ) then
          j = jj
          exit
        end if
      end do

      if ( j .eq. -1 ) then
        return
      end if

      return
      end
      subroutine r8cc_mxv ( m, n, nz_num, col, row, a, x, b )

c*********************************************************************72
c
cc R8CC_MXV multiplies an R8CC matrix by an R8VEC.
c
c  Discussion:
c
c    The R8CC format is the double precision sparse compressed column
c    format.  Associated with this format, we have an M by N matrix
c    with NZ_NUM nonzero entries.  We construct the column pointer
c    vector COL of length N+1, such that entries of column J will be
c    stored in positions COL(J) through COL(J+1)-1.  This indexing
c    refers to both the ROW and A vectors, which store the row indices
c    and the values of the nonzero entries.  The entries of the
c    ROW vector corresponding to each column are assumed to be
c    ascending sorted.
c
c    The R8CC format is equivalent to the MATLAB "sparse" format,
c    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    31 August 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Iain Duff, Roger Grimes, John Lewis,
c    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
c    October 1992
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c
c    Input, integer N, the number of columns of the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero elements in A.
c
c    Input, integer COL(N+1), points to the first element of 
c    each column.
c
c    Input, integer ROW(NZ_NUM), contains the row indices of 
c    the elements.
c
c    Input, double precision A(NZ_NUM), the R8CC matrix.
c
c    Input, double precision X(N), the vector to be multiplied.
c
c    Output, double precision B(M), the product A*X.
c
      implicit none

      integer m
      integer n
      integer nz_num

      double precision a(nz_num)
      double precision b(m)
      integer col(n+1)
      integer i
      integer j
      integer k
      integer row(nz_num)
      double precision x(n)

      b(1:m) = 0.0D+00

      do j = 1, n
        do k = col(j), col(j+1) - 1
          i = row(k)
          b(i) = b(i) + a(k) * x(j)
        end do
      end do

      return
      end
      subroutine r8cc_print ( m, n, nz_num, col, row, a, title )

c*********************************************************************72
c
cc R8CC_PRINT prints an R8CC matrix.
c
c  Discussion:
c
c    The R8CC format is the double precision sparse compressed column
c    format.  Associated with this format, we have an M by N matrix
c    with NZ_NUM nonzero entries.  We construct the column pointer
c    vector COL of length N+1, such that entries of column J will be
c    stored in positions COL(J) through COL(J+1)-1.  This indexing
c    refers to both the ROW and A vectors, which store the row indices
c    and the values of the nonzero entries.  The entries of the
c    ROW vector corresponding to each column are assumed to be
c    ascending sorted.
c
c    The R8CC format is equivalent to the MATLAB "sparse" format,
c    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    31 August 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Iain Duff, Roger Grimes, John Lewis,
c    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
c    October 1992
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c
c    Input, integer N, the number of columns of the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero elements in A.
c
c    Input, integer COL(N+1), points to the first element of 
c    each column.
c
c    Input, integer ROW(NZ_NUM), contains the row indices of 
c    the elements.
c
c    Input, double precision A(NZ_NUM), the R8CC matrix.
c
c    Input, character ( len = * ) TITLE, a title.
c 
      implicit none

      integer n
      integer nz_num

      double precision a(nz_num)
      integer col(n+1)
      integer m
      integer row(nz_num)
      character ( len = * )  title

      call r8cc_print_some ( m, n, nz_num, col, row, a, 1, 1, n, n, 
     &  title )

      return
      end
      subroutine r8cc_print_some ( m, n, nz_num, col, row, a, ilo, 
     &  jlo, ihi, jhi, title )

c*********************************************************************72
c
cc R8CC_PRINT_SOME prints some of an R8CC matrix.
c
c  Discussion:
c
c    The R8CC format is the double precision sparse compressed column
c    format.  Associated with this format, we have an M by N matrix
c    with NZ_NUM nonzero entries.  We construct the column pointer
c    vector COL of length N+1, such that entries of column J will be
c    stored in positions COL(J) through COL(J+1)-1.  This indexing
c    refers to both the ROW and A vectors, which store the row indices
c    and the values of the nonzero entries.  The entries of the
c    ROW vector corresponding to each column are assumed to be
c    ascending sorted.
c
c    The R8CC format is equivalent to the MATLAB "sparse" format,
c    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    31 August 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Iain Duff, Roger Grimes, John Lewis,
c    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
c    October 1992
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c
c    Input, integer N, the number of columns of the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero elements in A.
c
c    Input, integer COL(N+1), points to the first element of 
c    each column.
c
c    Input, integer ROW(NZ_NUM), contains the row indices of 
c    the elements.
c
c    Input, double precision A(NZ_NUM), the R8CC matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer, parameter :: incx = 5
      integer n
      integer nz_num

      double precision a(nz_num)
      double precision aij
      integer col(n+1)
      character ( len = 14 ) ctemp(incx)
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
      integer k
      integer m
      integer row(nz_num)
      character ( len = * ) title

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

        write ( *, '(''  Col:  '',5a14)' ) ( ctemp(j2), j2 = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) '  ---'
c
c  Determine the range of the rows in this strip.
c
        i2lo = max ( ilo, 1 )
        i2hi = min ( ihi, m )

        do i = i2lo, i2hi
c
c  Print out (up to) 5 entries in row I, that lie in the current strip.
c
c  1) Assume everything is zero.
c
          aij = 0.0D+00
          do j2 = 1, inc
            write ( ctemp(j2), '(f8.0,6x)' ) aij
          end do
c
c  2) Now consider each column J in J2LO to J2HI, 
c     and look at every nonzero, and check if it occurs in row I.
c
          do j = j2lo, j2hi
            do k = col(j), col(j+1)-1
              if ( row(k) .eq. i ) then
                j2 = j - j2lo + 1
                if ( r8_is_int ( a(k) ) ) then
                  write ( ctemp(j2), '(f8.0,6x)' ) a(k)
                else
                  write ( ctemp(j2), '(g14.6)' ) a(k)
                end if
              end if
            end do
          end do

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine r8cc_random ( m, n, nz_num, col, row, a, seed )

c*********************************************************************72
c
cc R8CC_RANDOM randomizes an R8CC matrix.
c
c  Discussion:
c
c    The R8CC format is the double precision sparse compressed column
c    format.  Associated with this format, we have an M by N matrix
c    with NZ_NUM nonzero entries.  We construct the column pointer
c    vector COL of length N+1, such that entries of column J will be
c    stored in positions COL(J) through COL(J+1)-1.  This indexing
c    refers to both the ROW and A vectors, which store the row indices
c    and the values of the nonzero entries.  The entries of the
c    ROW vector corresponding to each column are assumed to be
c    ascending sorted.
c
c    The R8CC format is equivalent to the MATLAB "sparse" format,
c    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    31 August 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Iain Duff, Roger Grimes, John Lewis,
c    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
c    October 1992
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c
c    Input, integer N, the number of columns of the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero elements in A.
c
c    Input, integer COL(N+1), points to the first element of each 
c    column.
c
c    Input, integer ROW(NZ_NUM), contains the row indices of the 
c    elements.
c
c    Output, double precision A(NZ_NUM), the R8CC matrix.
c
c    Input/output, integer SEED, a seed for the random number 
c    generator.
c
      implicit none

      integer n
      integer nz_num

      double precision a(nz_num)
      integer col(n+1)
      double precision r8_uniform_01
      integer j
      integer k
      integer m
      integer row(nz_num)
      integer seed

      do j = 1, n
        do k = col(j), col(j+1) - 1
          a(k) = r8_uniform_01 ( seed )
        end do
      end do

      return
      end
      subroutine r8cc_read ( col_file, row_file, a_file, m, n, 
     &  nz_num, col, row, a )

c*********************************************************************72
c
cc R8CC_READ reads an R8CC matrix from three files.
c
c  Discussion:
c
c    This routine needs the values of M, N, and NZ_NUM, which can be 
c    determined by a call to R8CC_READ_SIZE.
c
c    The R8CC format is the double precision sparse compressed column
c    format.  Associated with this format, we have an M by N matrix
c    with NZ_NUM nonzero entries.  We construct the column pointer
c    vector COL of length N+1, such that entries of column J will be
c    stored in positions COL(J) through COL(J+1)-1.  This indexing
c    refers to both the ROW and A vectors, which store the row indices
c    and the values of the nonzero entries.  The entries of the
c    ROW vector corresponding to each column are assumed to be
c    ascending sorted.
c
c    The R8CC format is equivalent to the MATLAB "sparse" format,
c    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 August 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Iain Duff, Roger Grimes, John Lewis,
c    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
c    October 1992
c
c  Parameters:
c
c    Input, character ( len = * ) COL_FILE, ROW_FILE, A_FILE, the names of the 
c    files containing the column pointers, row indices, and matrix entries.
c
c    Input, integer M, N, the number of rows and columns in the 
c    matrix.
c
c    Input, integer NZ_NUM, the number of nonzero elements in 
c    the matrix.
c
c    Output, integer COL(N+1), the column pointers.
c
c    Output, integer ROW(NZ_NUM), the row indices.
c
c    Output, double precision A(NZ_NUM), the nonzero elements 
c    of the matrix.
c
      implicit none

      integer n
      integer nz_num

      double precision a(nz_num)
      character ( len = * )  a_file
      integer col(n+1)
      character ( len = * )  col_file
      integer input_unit
      integer ios
      integer k
      integer m
      integer row(nz_num)
      character ( len = * )  row_file

      call get_unit ( input_unit )
c
c  Read the column information.
c
      open ( unit = input_unit, file = col_file, status = 'old', 
     &  iostat = ios )

      if ( ios .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8CC_READ - Fatal error!'
        write ( *, '(a)' ) '  Could not open the input file "' 
     &    // trim ( col_file ) // '".'
        stop
      end if

      do k = 1, n + 1

        read ( input_unit, *, iostat = ios ) col(k)

        if ( ios .ne. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8CC_READ - Fatal error!'
          write ( *, '(a,i8,a)' ) 
     &      '  I/O error while reading record ', k, 
     &      ' of "' // trim ( col_file ) // '".'
          stop
        end if

      end do

      close ( unit = input_unit )
c
c  Read the row information.
c
      open ( unit = input_unit, file = row_file, status = 'old', 
     &  iostat = ios )

      if ( ios .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8CC_READ - Fatal error!'
        write ( *, '(a)' ) '  Could not open the input file "' 
     &    // trim ( row_file ) // '".'
        stop
      end if

      do k = 1, nz_num

        read ( input_unit, *, iostat = ios ) row(k)

        if ( ios .ne. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8CC_READ - Fatal error!'
          write ( *, '(a,i8,a)' ) 
     &      '  I/O error while reading record ', k, 
     &      ' of "' // trim ( row_file ) // '".'
          stop
        end if

      end do

      close ( unit = input_unit )
c
c  Read the value information.
c
      open ( unit = input_unit, file = a_file, status = 'old', 
     &  iostat = ios )

      if ( ios .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8CC_READ - Fatal error!'
        write ( *, '(a)' ) '  Could not open the input file "' 
     &    // trim ( a_file ) // '".'
        stop
      end if

      do k = 1, nz_num

        read ( input_unit, *, iostat = ios ) a(k)

        if ( ios .ne. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8CC_READ - Fatal error!'
          write ( *, '(a,i8,a)' ) 
     &      '  I/O error while reading record ', k, 
     &      ' of "' // trim ( a_file ) // '".'
          stop
        end if

      end do

      close ( unit = input_unit )

      return
      end
      subroutine r8cc_read_size ( col_file, row_file, m, n, nz_num, 
     &  base )

c*********************************************************************72
c
cc R8CC_READ_SIZE reads the sizes of an R8CC sparse matrix from a file.
c
c  Discussion:
c
c    The value of M is "guessed" to be the largest value that occurs in
c    the ROW file.  However, if a row index of 0 is encountered, then
c    the value of M is incremented by 1.
c
c    The value of N is the number of records in the COL file minus 1.
c
c    The value of NZ_NUM is simply the number of records in the ROW file.
c
c    The value of BASE is 0 or 1, depending on whether the program
c    "guesses" that the row and column indices are 0-based or 1-based.
c    Although the first entry of the COL array might be used as evidence,
c    this program makes its determination based on whether it encounters
c    a 0 index in the ROW file.
c
c    The R8CC format is the double precision sparse compressed column
c    format.  Associated with this format, we have an M by N matrix
c    with NZ_NUM nonzero entries.  We construct the column pointer
c    vector COL of length N+1, such that entries of column J will be
c    stored in positions COL(J) through COL(J+1)-1.  This indexing
c    refers to both the ROW and A vectors, which store the row indices
c    and the values of the nonzero entries.  The entries of the
c    ROW vector corresponding to each column are assumed to be
c    ascending sorted.
c
c    The R8CC format is equivalent to the MATLAB "sparse" format,
c    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Iain Duff, Roger Grimes, John Lewis,
c    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
c    October 1992
c
c  Parameters:
c
c    Input, character ( len = * ) COL_FILE, ROW_FILE, the names of the 
c    column and row files that describe the structure of the matrix.
c
c    Output, integer M, N, the inferred number of rows and columns 
c    in the sparse matrix.
c
c    Output, integer NZ_NUM, the number of nonzero entries in the
c    sparse matrix.
c
c    Output, integer BASE, is 0 if the row indexing is believed
c    to be 0-based, and 1 if the row-index is believed to be
c    1-based.  In uncertain cases, BASE = 1 is the default.
c
      implicit none

      integer base
      integer col
      character ( len = * ) col_file
      integer input_unit
      integer ios
      integer m
      integer n
      integer nz_num
      integer row
      character ( len = * ) row_file
c
c  Default values.
c
      m = -1
      n = -1
      nz_num = -1
      base = -1
c
c  Check the COL file first.
c
      call get_unit ( input_unit )

      open ( unit = input_unit, file = col_file, status = 'old', 
     &  iostat = ios )

      if ( ios .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8CC_READ_SIZE - Fatal error!'
        write ( *, '(a)' ) '  Could not open  the input file "' 
     &    // trim ( col_file ) // '".'
        stop
      end if

      n = -1

      do

        read ( input_unit, *, iostat = ios ) col

        if ( ios .ne. 0 ) then
          exit
        end if

        n = n + 1

      end do

      close ( unit = input_unit )
c
c  Check the ROW file.
c
      call get_unit ( input_unit )

      open ( unit = input_unit, file = row_file, status = 'old', 
     &  iostat = ios )

      if ( ios .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8CC_READ_SIZE - Fatal error!'
        write ( *, '(a)' ) '  Could not open the input file "' 
     &    // trim ( row_file ) // '".'
        stop
      end if

      base = 1
      m = 0
      nz_num = 0
      
      do

        read ( input_unit, *, iostat = ios ) row

        if ( ios .ne. 0 ) then
          exit
        end if

        nz_num = nz_num + 1
        m = max ( m, row )
        if ( row .eq. 0 ) then
          base = 0
        end if

      end do

      close ( unit = input_unit )

      return
      end
      subroutine r8cc_set ( m, n, nz_num, col, row, a, i, j, aij )

c*********************************************************************72
c
cc R8CC_SET sets a value of an R8CC matrix.
c
c  Discussion:
c
c    The R8CC format is the double precision sparse compressed column
c    format.  Associated with this format, we have an M by N matrix
c    with NZ_NUM nonzero entries.  We construct the column pointer
c    vector COL of length N+1, such that entries of column J will be
c    stored in positions COL(J) through COL(J+1)-1.  This indexing
c    refers to both the ROW and A vectors, which store the row indices
c    and the values of the nonzero entries.  The entries of the
c    ROW vector corresponding to each column are assumed to be
c    ascending sorted.
c
c    The R8CC format is equivalent to the MATLAB "sparse" format,
c    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    31 August 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Iain Duff, Roger Grimes, John Lewis,
c    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
c    October 1992
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c
c    Input, integer N, the number of columns of the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero entries.
c
c    Input, integer COL(N+1), indicate where each column's data 
c    begins.
c
c    Input, integer ROW(NZ_NUM), the row indices.
c
c    Input/output, double precision A(NZ_NUM), the nonzero entries.
c    On output, the entry of A corresponding to (I,J) has been reset.
c
c    Input, integer I, J, the indices of the value to retrieve.
c
c    Input, double precision AIJ, the new value of A(I,J).
c
      implicit none

      integer n
      integer nz_num

      double precision a(nz_num)
      double precision aij
      integer col(n+1)
      integer i
      integer j
      integer k
      integer m
      integer row(nz_num)
c
c  Seek sparse index K corresponding to full index (I,J).
c
      call r8cc_ijk ( m, n, nz_num, col, row, i, j, k )
c
c  If no K was found, we fail.
c
      if ( k .eq. -1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8CC_SET - Fatal error!'
        write ( *, '(a)' ) '  R8CC_IJK could not find the entry.'
        write ( *, '(a,i8)' ) '  Row I = ', i
        write ( *, '(a,i8)' ) '  Col J = ', j
        stop
      end if

      a(k) = aij

      return
      end
      subroutine r8cc_to_r8ge ( m, n, nz_num, col, row, a, b )

c*********************************************************************72
c
cc R8CC_TO_R8GE converts an R8CC matrix to an R8GE matrix.
c
c  Discussion:
c
c    The R8CC format is the double precision sparse compressed column
c    format.  Associated with this format, we have an M by N matrix
c    with NZ_NUM nonzero entries.  We construct the column pointer
c    vector COL of length N+1, such that entries of column J will be
c    stored in positions COL(J) through COL(J+1)-1.  This indexing
c    refers to both the ROW and A vectors, which store the row indices
c    and the values of the nonzero entries.  The entries of the
c    ROW vector corresponding to each column are assumed to be
c    ascending sorted.
c
c    The R8CC format is equivalent to the MATLAB "sparse" format,
c    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    31 August 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Iain Duff, Roger Grimes, John Lewis,
c    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
c    October 1992
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c
c    Input, integer N, the number of columns of the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero elements in A.
c
c    Input, integer COL(N+1), points to the first element of 
c    each column.
c
c    Input, integer ROW(NZ_NUM), contains the row indices 
c    of the elements.
c
c    Input, double precision A(NZ_NUM), the R8CC matrix.
c
c    Output, double precision B(M,N), the R8GE matrix.
c
      implicit none

      integer m
      integer n
      integer nz_num

      double precision a(nz_num)
      double precision b(m,n)
      integer col(n+1)
      integer j
      integer k
      integer row(nz_num)

      b(1:m,1:n) = 0.0D+00

      do j = 1, n
        do k = col(j), col(j+1)-1
          b(row(k),j) = a(k)
        end do
      end do

      return
      end
      subroutine r8cc_vxm ( m, n, nz_num, col, row, a, x, b )

c*********************************************************************72
c
cc R8CC_VXM multiplies an R8VEC times an R8CC matrix.
c
c  Discussion:
c
c    The R8CC format is the double precision sparse compressed column
c    format.  Associated with this format, we have an M by N matrix
c    with NZ_NUM nonzero entries.  We construct the column pointer
c    vector COL of length N+1, such that entries of column J will be
c    stored in positions COL(J) through COL(J+1)-1.  This indexing
c    refers to both the ROW and A vectors, which store the row indices
c    and the values of the nonzero entries.  The entries of the
c    ROW vector corresponding to each column are assumed to be
c    ascending sorted.
c
c    The R8CC format is equivalent to the MATLAB "sparse" format,
c    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    31 August 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Iain Duff, Roger Grimes, John Lewis,
c    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
c    October 1992
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c
c    Input, integer N, the number of columns of the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero elements in A.
c
c    Input, integer COL(N+1), points to the first element 
c    of each column.
c
c    Input, integer ROW(NZ_NUM), contains the row indices 
c    of the elements.
c
c    Input, double precision A(NZ_NUM), the R8CC matrix.
c
c    Input, double precision X(M), the vector to be multiplied.
c
c    Output, double precision B(N), the product A'*X.
c
      implicit none

      integer m
      integer n
      integer nz_num

      double precision a(nz_num)
      double precision b(n)
      integer col(n+1)
      integer i
      integer j
      integer k
      integer row(nz_num)
      double precision x(m)

      b(1:n) = 0.0D+00

      do j = 1, n
        do k = col(j), col(j+1) - 1
          i = row(k)
          b(j) = b(j) + a(k) * x(i)
        end do
      end do

      return
      end
      subroutine r8cc_write ( col_file, row_file, a_file, m, n, 
     &  nz_num, col, row, a )

c*********************************************************************72
c
cc R8CC_WRITE writes an R8CC matrix to three files.
c
c  Discussion:
c
c    The R8CC format is the double precision sparse compressed column
c    format.  Associated with this format, we have an M by N matrix
c    with NZ_NUM nonzero entries.  We construct the column pointer
c    vector COL of length N+1, such that entries of column J will be
c    stored in positions COL(J) through COL(J+1)-1.  This indexing
c    refers to both the ROW and A vectors, which store the row indices
c    and the values of the nonzero entries.  The entries of the
c    ROW vector corresponding to each column are assumed to be
c    ascending sorted.
c
c    The R8CC format is equivalent to the MATLAB "sparse" format,
c    and the Harwell Boeing "real unsymmetric assembled" (RUA) format.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 September 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Iain Duff, Roger Grimes, John Lewis,
c    User's Guide for the Harwell-Boeing Sparse Matrix Collection,
c    October 1992
c
c  Parameters:
c
c    Input, character ( len = * ) COL_FILE, ROW_FILE, A_FILE, the names of the 
c    files containing the column pointers, row entries, and matrix entries.
c
c    Input, integer M, N, the number of rows and columns 
c    in the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero elements 
c    in the matrix.
c
c    Input, integer COL(N+1), the column pointers.
c
c    Input, integer ROW(NZ_NUM), the row indices.
c
c    Input, double precision A(NZ_NUM), the nonzero elements 
c    of the matrix.
c
      implicit none

      integer n
      integer nz_num

      double precision a(nz_num)
      character ( len = * ) a_file
      integer col(n+1)
      character ( len = * ) col_file
      integer ios
      integer k
      integer m
      integer output_unit
      integer row(nz_num)
      character ( len = * ) row_file

      call get_unit ( output_unit )
c
c  Write the column information.
c
      open ( unit = output_unit, file = col_file, status = 'replace', 
     &  iostat = ios )

      if ( ios .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8CC_WRITE - Fatal error!'
        write ( *, '(a)' ) '  Could not open the output file "' 
     &    // trim ( col_file ) // '".'
        stop
      end if

      do k = 1, n + 1

        write ( output_unit, '(i8)' ) col(k)

      end do

      close ( unit = output_unit )
c
c  Write the row information.
c
      open ( unit = output_unit, file = row_file, status = 'replace', 
     &  iostat = ios )

      if ( ios .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8CC_WRITE - Fatal error!'
        write ( *, '(a)' ) '  Could not open the output file "' 
     &    // trim ( row_file ) // '".'
        stop
      end if

      do k = 1, nz_num

        write ( output_unit, '(i8)' ) row(k)

      end do

      close ( unit = output_unit )
c
c  Write the value information.
c
      open ( unit = output_unit, file = a_file, status = 'replace', 
     &  iostat = ios )

      if ( ios .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8CC_WRITE - Fatal error!'
        write ( *, '(a)' ) '  Could not open the output file "' 
     &    // trim ( a_file ) // '".'
        stop
      end if

      do k = 1, nz_num

        write ( output_unit, '(g14.6)' ) a(k)

      end do

      close ( unit = output_unit )

      return
      end
      subroutine r8ci_eval ( n, a, lambda )

c*********************************************************************72
c
cc R8CI_EVAL returns the eigenvalues of an R8CI matrix.
c
c  Discussion:
c
c    The R8CI storage format is used for a real N by N circulant matrix.
c    An N by N circulant matrix A has the property that the entries on
c    row I appear again on row I+1, shifted one position to the right,
c    with the final entry of row I appearing as the first of row I+1.
c    The R8CI format simply records the first row of the matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2005
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Philip Davis,
c    Circulant Matrices,
c    Second Edition,
c    Chelsea, 1994,
c    ISBN: 0828403384,
c    LC: QA188.D37.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision A(N), the R8CI matrix.
c
c    Output, double complex LAMBDA(N), the eigenvalues.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      double complex lambda(n)
      double complex w(n)

      call c8vec_unity ( n, w )

      lambda(1:n) = dcmplx ( a(n), 0.0D+00 )
      do i = n-1, 1, -1
        lambda(1:n) = lambda(1:n) * w(1:n) + dcmplx ( a(i), 0.0D+00 )
      end do

      call c8vec_sort_a2 ( n, lambda )

      return
      end
      subroutine r8ci_indicator ( n, a )

c*********************************************************************72
c
cc R8CI_INDICATOR sets up an R8CI indicator matrix.
c
c  Discussion:
c
c    The "indicator matrix" simply has a value like I*10+J at every
c    entry of a dense matrix, or at every entry of a compressed storage
c    matrix for which storage is allocated. 
c
c    The R8CI storage format is used for a real N by N circulant matrix.
c    An N by N circulant matrix A has the property that the entries on
c    row I appear again on row I+1, shifted one position to the right,
c    with the final entry of row I appearing as the first of row I+1.
c    The R8CI format simply records the first row of the matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2004
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
c    Output, double precision A(N), the R8CI matrix.
c
      implicit none

      integer n

      double precision a(n)
      integer fac
      integer i
      integer i4_log_10
      integer j

      fac = 10 ** ( i4_log_10 ( n ) + 1 )

      i = 1

      do j = 1, n
        a(j) = dble ( fac * i + j )
      end do

      return
      end
      subroutine r8ci_mxv ( n, a, x, b )

c*********************************************************************72
c
cc R8CI_MXV multiplies an R8CI matrix by an R8VEC.
c
c  Discussion:
c
c    The R8CI storage format is used for a real N by N circulant matrix.
c    An N by N circulant matrix A has the property that the entries on
c    row I appear again on row I+1, shifted one position to the right,
c    with the final entry of row I appearing as the first of row I+1.
c    The R8CI format simply records the first row of the matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 November 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision A(N), the R8CI matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(N), the product A * x.
c
      implicit none

      integer n

      double precision a(n)
      double precision b(n)
      integer i
      double precision x(n)

      do i = 1, n
        b(i) = sum ( a(n+2-i:n) * x(1:i-1) ) 
     &       + sum ( a(1:n+1-i) * x(i:n) )
      end do

      return
      end
      subroutine r8ci_print ( n, a, title )

c*********************************************************************72
c
cc R8CI_PRINT prints an R8CI matrix.
c
c  Discussion:
c
c    The R8CI storage format is used for a real N by N circulant matrix.
c    An N by N circulant matrix A has the property that the entries on
c    row I appear again on row I+1, shifted one position to the right,
c    with the final entry of row I appearing as the first of row I+1.
c    The R8CI format simply records the first row of the matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 May 2000
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
c    Input, double precision A(N), the R8CI matrix.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer n

      double precision a(n)
      character ( len = * ) title

      call r8ci_print_some ( n, a, 1, 1, n, n, title )

      return
      end
      subroutine r8ci_print_some ( n, a, ilo, jlo, ihi, jhi, title )

c*********************************************************************72
c
cc R8CI_PRINT_SOME prints some of an R8CI matrix.
c
c  Discussion:
c
c    The R8CI storage format is used for a real N by N circulant matrix.
c    An N by N circulant matrix A has the property that the entries on
c    row I appear again on row I+1, shifted one position to the right,
c    with the final entry of row I appearing as the first of row I+1.
c    The R8CI format simply records the first row of the matrix.
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
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input, double precision A(N), the R8CI matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 5 )
      integer n

      double precision a(n)
      double precision aij
      character ( len = 14 ) ctemp(incx)
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
      character ( len = * ) title

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
        i2hi = min ( ihi, n )

        do i = i2lo, i2hi
c
c  Print out (up to) 5 entries in row I, that lie in the current strip.
c
          do j2 = 1, inc

            j = j2lo - 1 + j2

            if ( i .le. j ) then
              aij = a(j+1-i)
            else
              aij = a(n+j+1-i)
            end if

            if ( r8_is_int ( aij ) ) then
              write ( ctemp(j2), '(f8.0,6x)' ) aij
            else
              write ( ctemp(j2), '(g14.6)' ) aij
            end if

          end do

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine r8ci_random ( n, seed, a )

c*********************************************************************72
c
cc R8CI_RANDOM randomizes an R8CI matrix.
c
c  Discussion:
c
c    The R8CI storage format is used for a real N by N circulant matrix.
c    An N by N circulant matrix A has the property that the entries on
c    row I appear again on row I+1, shifted one position to the right,
c    with the final entry of row I appearing as the first of row I+1.
c    The R8CI format simply records the first row of the matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 May 2003
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
c    Input/output, integer SEED, a seed for the random number 
c    generator.
c
c    Output, double precision A(N), the R8CI matrix.
c
      implicit none

      integer n

      double precision a(n)
      double precision r8_uniform_01
      integer i
      integer seed

      do i = 1, n
        a(i) = r8_uniform_01 ( seed )
      end do

      return
      end
      subroutine r8ci_sl ( n, a, b, x, job )

c*********************************************************************72
c
cc R8CI_SL solves an R8CI system.
c
c  Discussion:
c
c    The R8CI storage format is used for a real N by N circulant matrix.
c    An N by N circulant matrix A has the property that the entries on
c    row I appear again on row I+1, shifted one position to the right,
c    with the final entry of row I appearing as the first of row I+1.
c    The R8CI format simply records the first row of the matrix.
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
c    Input, integer N, the order of the matrix.
c
c    Input, double precision A(N), the R8CI matrix.
c
c    Input, double precision B(N), the right hand side.
c
c    Output, double precision X(N), the solution of the linear system.
c
c    Input, integer JOB, specifies the system to solve.
c    0, solve A * x = b.
c    nonzero, solve A' * x = b.
c
      implicit none

      integer n

      double precision a(n)
      double precision b(n)
      integer i
      integer job
      integer nsub
      double precision r1
      double precision r2
      double precision r3
      double precision r5
      double precision r6
      double precision work(2*n-2)
      double precision x(n)

      if ( job .eq. 0 ) then
c
c  Solve the system with the principal minor of order 1.
c
        r1 = a(1)
        x(1) = b(1) / r1

        r2 = 0.0D+00
c
c  Recurrent process for solving the system.
c
        do nsub = 2, n
c
c  Compute multiples of the first and last columns of
c  the inverse of the principal minor of order N.
c
          r5 = a(n+2-nsub)
          r6 = a(nsub)

          if ( 2 .lt. nsub ) then

            work(nsub-1) = r2

            do i = 1, nsub-2
              r5 = r5 + a(n+1-i) * work(nsub-i)
              r6 = r6 + a(i+1) * work(n-1+i)
            end do

          end if

          r2 = - r5 / r1
          r3 = - r6 / r1
          r1 = r1 + r5 * r3

          if ( 2 .lt. nsub ) then

            r6 = work(n)
            work(n-1+nsub-1) = 0.0D+00
            do i = 2, nsub-1
              r5 = work(n-1+i)
              work(n-1+i) = work(i) * r3 + r6
              work(i) = work(i) + r6 * r2
              r6 = r5
            end do

          end if

          work(n) = r3
c
c  Compute the solution of the system with the principal minor of order NSUB.
c
          r5 = 0.0D+00
          do i = 1, nsub-1
            r5 = r5 + a(n+1-i) * x(nsub-i)
          end do

          r6 = ( b(nsub) - r5 ) / r1
          x(1:nsub-1) = x(1:nsub-1) + work(n:n+nsub-2) * r6
          x(nsub) = r6

        end do

      else
c
c  Solve the system with the principal minor of order 1.
c
        r1 = a(1)
        x(1) = b(1) / r1

        r2 = 0.0D+00
c
c  Recurrent process for solving the system.
c
        do nsub = 2, n
c
c  Compute multiples of the first and last columns of
c  the inverse of the principal minor of order N.
c
          r5 = a(nsub)
          r6 = a(n+2-nsub)

          if ( 2 .lt. nsub ) then

            work(nsub-1) = r2

            do i = 1, nsub-2
              r5 = r5 + a(i+1) * work(nsub-i)
              r6 = r6 + a(n+1-i) * work(n-1+i)
            end do

          end if

          r2 = - r5 / r1
          r3 = - r6 / r1
          r1 = r1 + r5 * r3

          if ( 2 .lt. nsub ) then

            r6 = work(n)
            work(n-1+nsub-1) = 0.0D+00
            do i = 2, nsub-1
              r5 = work(n-1+i)
              work(n-1+i) = work(i) * r3 + r6
              work(i) = work(i) + r6 * r2
              r6 = r5
            end do

          end if

          work(n) = r3
c
c  Compute the solution of the system with the principal minor of order NSUB.
c
          r5 = 0.0D+00
          do i = 1, nsub-1
            r5 = r5 + a(i+1) * x(nsub-i)
          end do

          r6 = ( b(nsub) - r5 ) / r1
          do i = 1, nsub-1
            x(i) = x(i) + work(n-1+i) * r6
          end do

          x(nsub) = r6

        end do

      end if

      return
      end
      subroutine r8ci_to_r8ge ( n, a, b )

c*********************************************************************72
c
cc R8CI_TO_R8GE copies an R8CI matrix into an R8GE matrix.
c
c  Discussion:
c
c    The R8CI storage format is used for a real N by N circulant matrix.
c    An N by N circulant matrix A has the property that the entries on
c    row I appear again on row I+1, shifted one position to the right,
c    with the final entry of row I appearing as the first of row I+1.
c    The R8CI format simply records the first row of the matrix.
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 November 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision A(N), the R8CI matrix.
c
c    Output, double precision B(N,N), the R8GE matrix.
c
      implicit none

      integer n

      double precision a(n)
      double precision b(n,n)
      integer i

      do i = 1, n
        b(i,1:i-1) = a(n+2-i:n+2*1-2*i)
        b(i,i:n) = a(1:n+1-i)
      end do

      return
      end
      subroutine r8ci_vxm ( n, a, x, b )

c*********************************************************************72
c
cc R8CI_VXM multiplies an R8VEC by an R8CI matrix.
c
c  Discussion:
c
c    The R8CI storage format is used for a real N by N circulant matrix.
c    An N by N circulant matrix A has the property that the entries on
c    row I appear again on row I+1, shifted one position to the right,
c    with the final entry of row I appearing as the first of row I+1.
c    The R8CI format simply records the first row of the matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 March 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision A(N), the R8CI matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(N), the product A' * X.
c
      implicit none

      integer n

      double precision a(n)
      double precision b(n)
      integer i
      double precision x(n)

      do i = 1, n
        b(i) = sum ( a(i:1:-1) * x(1:i) ) 
     &       + sum ( a(n:i+1:-1) * x(i+1:n) )
      end do

      return
      end
      subroutine r8col_swap ( m, n, a, j1, j2 )

c*********************************************************************72
c
cc R8COL_SWAP swaps columns I and J of an R8COL.
c
c  Discussion:
c
c    An R8COL is an M by N array of R8's, regarded as an array of N columns,
c    each of length M.
c
c  Example:
c
c    Input:
c
c      M = 3, N = 4, J1 = 2, J2 = 4
c
c      A = (
c        1.  2.  3.  4.
c        5.  6.  7.  8.
c        9. 10. 11. 12. )
c
c    Output:
c
c      A = (
c        1.  4.  3.  2.
c        5.  8.  7.  6.
c        9. 12. 11. 10. )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 March 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input/output, double precision A(M,N), the M by N array.
c
c    Input, integer J1, J2, the columns to be swapped.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer i
      integer j1
      integer j2
      double precision temp

      if ( j1 .lt. 1 .or. n .lt. j1 .or. j2 .lt. 1 .or. n .lt. j2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8COL_SWAP - Fatal error!'
        write ( *, '(a)' ) '  J1 or J2 is out of bounds.'
        write ( *, '(a,i8)' ) '  J1 =    ', j1
        write ( *, '(a,i8)' ) '  J2 =    ', j2
        write ( *, '(a,i8)' ) '  NCOL = ', n
        stop
      end if

      if ( j1 .eq. j2 ) then
        return
      end if

      do i = 1, m
        temp = a(i,j1)
        a(i,j1) = a(i,j2)
        a(i,j2) = temp
      end do

      return
      end
      subroutine r8gb_det ( n, ml, mu, a_lu, pivot, det )

c*********************************************************************72
c
cc R8GB_DET computes the determinant of a matrix factored by R8GB_FA or R8GB_TRF.
c
c  Discussion:
c
c    The R8GB storage format is for an M by N banded matrix, with lower 
c    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
c    extra superdiagonals, which may be required to store nonzero entries 
c    generated during Gaussian elimination.
c
c    The original M by N matrix is "collapsed" downward, so that diagonals
c    become rows of the storage array, while columns are preserved.  The
c    collapsed array is logically 2*ML+MU+1 by N.  
c
c    R8GB storage is used by LINPACK and LAPACK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 February 2001
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
c    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
c    Sven Hammarling, Alan McKenney, Danny Sorensen,
c    LAPACK User's Guide,
c    Second Edition,
c    SIAM, 1995.
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
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than N-1.
c
c    Input, double precision A_LU(2*ML+MU+1,N), the LU factors from 
c    R8GB_FA or R8GB_TRF.
c
c    Input, integer PIVOT(N), the pivot vector, as computed 
c    by R8GB_FA or R8GB_TRF.
c
c    Output, double precision DET, the determinant of the matrix.
c
      implicit none

      integer ml
      integer mu
      integer n

      double precision a_lu(2*ml+mu+1,n)
      double precision det
      integer i
      integer pivot(n)
      double precision r8vec_product

      det = r8vec_product ( n, a_lu(ml+mu+1,1:n) )

      do i = 1, n
        if ( pivot(i) .ne. i ) then
          det = -det
        end if
      end do

      return
      end
      subroutine r8gb_fa ( n, ml, mu, a, pivot, info )

c*********************************************************************72
c
cc R8GB_FA performs a LINPACK-style PLU factorization of an R8GB matrix.
c
c  Discussion:
c
c    The R8GB storage format is for an M by N banded matrix, with lower 
c    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
c    extra superdiagonals, which may be required to store nonzero entries 
c    generated during Gaussian elimination.
c
c    The original M by N matrix is "collapsed" downward, so that diagonals
c    become rows of the storage array, while columns are preserved.  The
c    collapsed array is logically 2*ML+MU+1 by N.  
c
c    This routine is based on the LINPACK routine SGBFA.
c
c    R8GB storage is used by LINPACK and LAPACK.
c
c    The following program segment will set up the input.
c
c      m = ml + mu + 1
c      do j = 1, n
c        i1 = max ( 1, j-mu )
c        i2 = min ( n, j+ml )
c        do i = i1, i2
c          k = i - j + m
c          a(k,j) = afull(i,j)
c        end do
c      end do
c
c    This uses rows ML+1 through 2*ML+MU+1 of the array A.
c    In addition, the first ML rows in the array are used for
c    elements generated during the triangularization.
c
c    The ML+MU by ML+MU upper left triangle and the
c    ML by ML lower right triangle are not referenced.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 March 1999
c
c  Author:
c
c    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
c    This version by John Burkardt
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
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than N-1.
c
c    Input/output, double precision A(2*ML+MU+1,N), on input, 
c    the matrix in band storage, on output, information about 
c    the LU factorization.
c
c    Output, integer PIVOT(N), the pivot vector.
c
c    Output, integer INFO, singularity flag.
c    0, no singularity detected.
c    nonzero, the factorization failed on the INFO-th step.
c
      implicit none

      integer ml
      integer mu
      integer n

      double precision a(2*ml+mu+1,n)
      integer i0
      integer info
      integer pivot(n)
      integer j
      integer j0
      integer j1
      integer ju
      integer jz
      integer k
      integer l
      integer lm
      integer m
      integer mm
      double precision t

      m = ml + mu + 1
      info = 0
c
c  Zero out the initial fill-in columns.
c
      j0 = mu + 2
      j1 = min ( n, m ) - 1

      do jz = j0, j1
        i0 = m + 1 - jz
        a(i0:ml,jz) = 0.0D+00
      end do

      jz = j1
      ju = 0

      do k = 1, n-1
c
c  Zero out the next fill-in column.
c
        jz = jz + 1
        if ( jz .le. n ) then
          a(1:ml,jz) = 0.0D+00
        end if
c
c  Find L = pivot index.
c
        lm = min ( ml, n-k )

        l = m
        do j = m+1, m+lm
          if ( abs ( a(l,k) ) .lt. abs ( a(j,k) ) ) then
            l = j
          end if
        end do

        pivot(k) = l + k - m
c
c  Zero pivot implies this column already triangularized.
c
        if ( a(l,k) .eq. 0.0D+00 ) then
          info = k
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8GB_FA - Fatal error!'
          write ( *, '(a,i8)' ) '  Zero pivot on step ', info
          stop
        end if
c
c  Interchange if necessary.
c
        t      = a(l,k)
        a(l,k) = a(m,k)
        a(m,k) = t
c
c  Compute multipliers.
c
        a(m+1:m+lm,k) = - a(m+1:m+lm,k) / a(m,k)
c
c  Row elimination with column indexing.
c
        ju = max ( ju, mu+pivot(k) )
        ju = min ( ju, n )
        mm = m

        do j = k+1, ju

          l = l - 1
          mm = mm - 1

          if ( l .ne. mm ) then
            t       = a(l,j)
            a(l,j)  = a(mm,j)
            a(mm,j) = t
          end if

          a(mm+1:mm+lm,j) = a(mm+1:mm+lm,j) + a(mm,j) * a(m+1:m+lm,k)

        end do

      end do

      pivot(n) = n

      if ( a(m,n) .eq. 0.0D+00 ) then
        info = n
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8GB_FA - Fatal error!'
        write ( *, '(a,i8)' ) '  Zero pivot on step ', info
        stop
      end if

      return
      end
      subroutine r8gb_indicator ( m, n, ml, mu, a )

c*********************************************************************72
c
cc R8GB_INDICATOR sets up an R8GB indicator matrix.
c
c  Discussion:
c
c    The "indicator matrix" simply has a value like I*10+J at every
c    entry of a dense matrix, or at every entry of a compressed storage
c    matrix for which storage is allocated. 
c
c    Note that the R8GB storage format includes extra room for
c    fillin entries that occur during Gauss elimination.  These entries
c    are not normally seen or used by the user.  This routine will
c    set those values to zero.
c
c    The R8GB storage format is for an M by N banded matrix, with lower 
c    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
c    extra superdiagonals, which may be required to store nonzero entries 
c    generated during Gaussian elimination.
c
c    The original M by N matrix is "collapsed" downward, so that diagonals
c    become rows of the storage array, while columns are preserved.  The
c    collapsed array is logically 2*ML+MU+1 by N.  
c
c    R8GB storage is used by LINPACK and LAPACK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 March 2005
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than min(M,N)-1.
c
c    Output, double precision A(2*ML+MU+1,N), the R8GB matrix.
c
      implicit none

      integer ml
      integer mu
      integer m
      integer n

      double precision a(2*ml+mu+1,n)
      integer diag
      integer fac
      integer i
      integer i4_log_10
      integer j
      integer k
      double precision value

      fac = 10 ** ( i4_log_10 ( n ) + 1 )
      k = 0

      do j = 1, n
        do diag = 1, 2 * ml + mu + 1

          i = diag + j - ml - mu - 1

          if ( 1 .le. i .and. i .le. m .and. i - ml .le. j 
     &      .and. j .le. i + mu ) then
            value = dble ( fac * i + j )
          else if ( 1 .le. i .and. i .le. m .and. 
     &      i - ml .le. j .and. j .le. i + mu + ml ) then
            value = 0.0D+00
          else
            k = k + 1
            value = - dble ( k )
          end if

          a(diag,j) = value

        end do
      end do

      return
      end
      subroutine r8gb_ml ( n, ml, mu, a_lu, pivot, x, b, job )

c*********************************************************************72
c
cc R8GB_ML computes A * x or A' * X, using R8GB_FA factors.
c
c  Discussion:
c
c    The R8GB storage format is for an M by N banded matrix, with lower 
c    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
c    extra superdiagonals, which may be required to store nonzero entries 
c    generated during Gaussian elimination.
c
c    The original M by N matrix is "collapsed" downward, so that diagonals
c    become rows of the storage array, while columns are preserved.  The
c    collapsed array is logically 2*ML+MU+1 by N.  
c
c    R8GB storage is used by LINPACK and LAPACK.
c
c    It is assumed that R8GB_FA has overwritten the original matrix
c    information by LU factors.  R8GB_ML is able to reconstruct the
c    original matrix from the LU factor data.
c
c    R8GB_ML allows the user to check that the solution of a linear
c    system is correct, without having to save an unfactored copy
c    of the matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    31 October 1998
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
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than N-1.
c
c    Input, double precision A_LU(2*ML+MU+1,N), the LU factors from R8GB_FA.
c
c    Input, integer PIVOT(N), the pivot vector computed by R8GB_FA.
c
c    Input, double precision X(N), the vector to be multiplied.
c
c    Output, double precision B(N), the result of the multiplication.
c
c    Input, integer JOB, specifies the operation to be done:
c    JOB = 0, compute A * x.
c    JOB nonzero, compute A' * X.
c
      implicit none

      integer ml
      integer mu
      integer n

      double precision a_lu(2*ml+mu+1,n)
      double precision b(n)
      integer i
      integer ihi
      integer ilo
      integer pivot(n)
      integer j
      integer jhi
      integer job
      integer k
      double precision x(n)

      b(1:n) = x(1:n)

      if ( job .eq. 0 ) then
c
c  Y = U * X.
c
        do j = 1, n
          ilo = max ( 1, j - ml - mu )
          do i = ilo, j - 1
            b(i) = b(i) + a_lu(i-j+ml+mu+1,j) * b(j)
          end do
          b(j) = a_lu(j-j+ml+mu+1,j) * b(j)
        end do
c
c  B = PL * Y = PL * U * X = A * x.
c
        do j = n-1, 1, -1

          ihi = min ( n, j + ml )
          do i = j+1, ihi
            b(i) = b(i) - a_lu(i-j+ml+mu+1,j) * b(j)
          end do

          k = pivot(j)

          if ( k .ne. j ) then
            call r8_swap ( b(k), b(j) )
          end if

        end do

      else
c
c  Y = ( PL )' * X.
c
        do j = 1, n-1

          k = pivot(j)

          if ( k .ne. j ) then
            call r8_swap ( b(k), b(j) )
          end if

          jhi = min ( n, j + ml )
          do i = j+1, jhi
            b(j) = b(j) - b(i) * a_lu(i-j+ml+mu+1,j)
          end do

        end do
c
c  B = U' * Y = ( PL * U )' * X = A' * X.
c
        do i = n, 1, -1

          jhi = min ( n, i + ml + mu )
          do j = i+1, jhi
            b(j) = b(j) + b(i) * a_lu(i-j+ml+mu+1,j)
          end do
          b(i) = b(i) * a_lu(i-i+ml+mu+1,i)
        end do

      end if

      return
      end
      subroutine r8gb_mu ( n, ml, mu, a_lu, pivot, x, b, job )

c*********************************************************************72
c
cc R8GB_MU computes A * x or A' * X, using R8GB_TRF factors.
c
c  Warning:
c
c    This routine needs to be updated to allow for rectangular matrices.
c
c  Discussion:
c
c    The R8GB storage format is for an M by N banded matrix, with lower 
c    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
c    extra superdiagonals, which may be required to store nonzero entries 
c    generated during Gaussian elimination.
c
c    The original M by N matrix is "collapsed" downward, so that diagonals
c    become rows of the storage array, while columns are preserved.  The
c    collapsed array is logically 2*ML+MU+1 by N.  
c
c    R8GB storage is used by LINPACK and LAPACK.
c
c    It is assumed that R8GB_TRF has overwritten the original matrix
c    information by LU factors.  R8GB_MU is able to reconstruct the
c    original matrix from the LU factor data.
c
c    R8GB_MU allows the user to check that the solution of a linear
c    system is correct, without having to save an unfactored copy
c    of the matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
c    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
c    Sven Hammarling, Alan McKenney, Danny Sorensen,
c    LAPACK User's Guide,
c    Second Edition,
c    SIAM, 1995.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than N-1.
c
c    Input, double precision A_LU(2*ML+MU+1,N), the LU factors from R8GB_TRF.
c
c    Input, integer PIVOT(N), the pivot vector computed by R8GB_TRF.
c
c    Input, double precision X(N), the vector to be multiplied.
c
c    Output, double precision B(N), the result of the multiplication.
c
c    Input, integer JOB, specifies the operation to be done:
c    JOB = 0, compute A * x.
c    JOB nonzero, compute A' * X.
c
      implicit none

      integer ml
      integer mu
      integer n

      double precision a_lu(2*ml+mu+1,n)
      double precision b(n)
      integer i
      integer ihi
      integer ilo
      integer pivot(n)
      integer j
      integer jhi
      integer job
      integer k
      double precision x(n)

      b(1:n) = x(1:n)

      if ( job .eq. 0 ) then
c
c  Y = U * X.
c
        do j = 1, n
          ilo = max ( 1, j - ml - mu )
          do i = ilo, j - 1
            b(i) = b(i) + a_lu(i-j+ml+mu+1,j) * b(j)
          end do
          b(j) = a_lu(j-j+ml+mu+1,j) * b(j)
        end do
c
c  B = PL * Y = PL * U * X = A * x.
c
        do j = n-1, 1, -1

          ihi = min ( n, j + ml )
          do i = j+1, ihi
            b(i) = b(i) + a_lu(i-j+ml+mu+1,j) * b(j)
          end do

          k = pivot(j)

          if ( k .ne. j ) then
            call r8_swap ( b(k), b(j) )
          end if

        end do

      else
c
c  Y = ( PL )' * X.
c
        do j = 1, n-1

          k = pivot(j)

          if ( k .ne. j ) then
            call r8_swap ( b(k), b(j) )
          end if

          jhi = min ( n, j + ml )
          do i = j+1, jhi
            b(j) = b(j) + b(i) * a_lu(i-j+ml+mu+1,j)
          end do

        end do
c
c  B = U' * Y = ( PL * U )' * X = A' * X.
c
        do i = n, 1, -1

          jhi = min ( n, i + ml + mu )
          do j = i+1, jhi
            b(j) = b(j) + b(i) * a_lu(i-j+ml+mu+1,j)
          end do
          b(i) = b(i) * a_lu(i-i+ml+mu+1,i)
        end do

      end if

      return
      end
      subroutine r8gb_mxv ( m, n, ml, mu, a, x, b )

c*********************************************************************72
c
cc R8GB_MXV multiplies an R8GB matrix by an R8VEC.
c
c  Discussion:
c
c    The R8GB storage format is for an M by N banded matrix, with lower 
c    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
c    extra superdiagonals, which may be required to store nonzero entries 
c    generated during Gaussian elimination.
c
c    The original M by N matrix is "collapsed" downward, so that diagonals
c    become rows of the storage array, while columns are preserved.  The
c    collapsed array is logically 2*ML+MU+1 by N.  
c
c    R8GB storage is used by LINPACK and LAPACK.
c
c    LINPACK and LAPACK storage of general band matrices requires
c    an extra ML upper diagonals for possible fill in entries during
c    Gauss elimination.  This routine does not access any entries
c    in the fill in diagonals, because it assumes that the matrix
c    has NOT had Gauss elimination applied to it.  If the matrix
c    has been Gauss eliminated, then the routine R8GB_MU must be
c    used instead.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 January 1999
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
c    Input, integer M, the number of rows of the matrix.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than min(M,N)-1.
c
c    Input, double precision A(2*ML+MU+1,N), the R8GB matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(M), the product A * x.
c
      implicit none

      integer ml
      integer mu
      integer m
      integer n

      double precision a(2*ml+mu+1,n)
      double precision b(m)
      integer i
      integer j
      integer jhi
      integer jlo
      double precision x(n)

      do i = 1, m
        b(i) = 0.0D+00
        jlo = max ( 1, i - ml )
        jhi = min ( n, i + mu )
        do j = jlo, jhi
          b(i) = b(i) + a(i-j+ml+mu+1,j) * x(j)
        end do
      end do

      return
      end
      subroutine r8gb_nz_num ( m, n, ml, mu, a, nz_num )

c*********************************************************************72
c
cc R8GB_NZ_NUM counts the nonzeroes in an R8GB matrix.
c
c  Discussion:
c
c    The R8GB storage format is for an M by N banded matrix, with lower 
c    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
c    extra superdiagonals, which may be required to store nonzero entries 
c    generated during Gaussian elimination.
c
c    The original M by N matrix is "collapsed" downward, so that diagonals
c    become rows of the storage array, while columns are preserved.  The
c    collapsed array is logically 2*ML+MU+1 by N.  
c
c    R8GB storage is used by LINPACK and LAPACK.
c
c    LINPACK and LAPACK band storage requires that an extra ML
c    superdiagonals be supplied to allow for fillin during Gauss
c    elimination.  Even though a band matrix is described as
c    having an upper bandwidth of MU, it effectively has an
c    upper bandwidth of MU+ML.  This routine will examine
c    values it finds in these extra bands, so that both unfactored
c    and factored matrices can be handled.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 October 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than min(M,N)-1.
c
c    Input, double precision A(2*ML+MU+1,N), the R8GB matrix.
c
c    Output, integer NZ_NUM, the number of nonzero entries in A.
c
      implicit none

      integer ml
      integer mu
      integer m
      integer n

      double precision a(2*ml+mu+1,n)
      integer i
      integer j
      integer jhi
      integer jlo
      integer nz_num

      nz_num = 0

      do i = 1, m
        jlo = max ( 1, i - ml )
        jhi = min ( n, i + mu + ml )
        do j = jlo, jhi
          if ( a(i-j+ml+mu+1,j) .ne. 0.0D+00 ) then
            nz_num = nz_num + 1
          end if
        end do
      end do

      return
      end
      subroutine r8gb_print ( m, n, ml, mu, a, title )

c*********************************************************************72
c
cc R8GB_PRINT prints an R8GB matrix.
c
c  Discussion:
c
c    The R8GB storage format is for an M by N banded matrix, with lower 
c    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
c    extra superdiagonals, which may be required to store nonzero entries 
c    generated during Gaussian elimination.
c
c    The original M by N matrix is "collapsed" downward, so that diagonals
c    become rows of the storage array, while columns are preserved.  The
c    collapsed array is logically 2*ML+MU+1 by N.  
c
c    R8GB storage is used by LINPACK and LAPACK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than min(M,N)-1..
c
c    Input, double precision A(2*ML+MU+1,N), the R8GB matrix.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer ml
      integer mu
      integer n

      double precision a(2*ml+mu+1,n)
      integer m
      character ( len = * )  title

      call r8gb_print_some ( m, n, ml, mu, a, 1, 1, m, n, title )

      return
      end
      subroutine r8gb_print_some ( m, n, ml, mu, a, ilo, jlo, ihi, 
     &  jhi, title )

c*********************************************************************72
c
cc R8GB_PRINT_SOME prints some of an R8GB matrix.
c
c  Discussion:
c
c    The R8GB storage format is for an M by N banded matrix, with lower 
c    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
c    extra superdiagonals, which may be required to store nonzero entries 
c    generated during Gaussian elimination.
c
c    The original M by N matrix is "collapsed" downward, so that diagonals
c    become rows of the storage array, while columns are preserved.  The
c    collapsed array is logically 2*ML+MU+1 by N.  
c
c    R8GB storage is used by LINPACK and LAPACK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 April 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than min(M,N)-1..
c
c    Input, double precision A(2*ML+MU+1,N), the R8GB matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 5 )
      integer ml
      integer mu
      integer n

      double precision a(2*ml+mu+1,n)
      character ( len = 14 ) ctemp(incx)
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
      integer m
      character ( len = * )  title

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
        i2lo = max ( i2lo, j2lo - mu - ml )
        i2hi = min ( ihi, m )
        i2hi = min ( i2hi, j2hi + ml )

        do i = i2lo, i2hi
c
c  Print out (up to) 5 entries in row I, that lie in the current strip.
c
          do j2 = 1, inc

            j = j2lo - 1 + j2

            if ( i .lt. j - ml - mu  .or. j + ml .lt. i ) then
              ctemp(j2) = '              '
            else
              write ( ctemp(j2), '(g14.6)' ) a(i-j+ml+mu+1,j)
            end if

          end do

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine r8gb_random ( m, n, ml, mu, seed, a )

c*********************************************************************72
c
cc R8GB_RANDOM randomizes an R8GB matrix.
c
c  Discussion:
c
c    The R8GB storage format is for an M by N banded matrix, with lower 
c    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
c    extra superdiagonals, which may be required to store nonzero entries 
c    generated during Gaussian elimination.
c
c    The original M by N matrix is "collapsed" downward, so that diagonals
c    become rows of the storage array, while columns are preserved.  The
c    collapsed array is logically 2*ML+MU+1 by N.  
c
c    R8GB storage is used by LINPACK and LAPACK.
c
c    LINPACK and LAPACK band storage requires that an extra ML
c    superdiagonals be supplied to allow for fillin during Gauss
c    elimination.  Even though a band matrix is described as
c    having an upper bandwidth of MU, it effectively has an
c    upper bandwidth of MU+ML.  This routine assumes it is setting
c    up an unfactored matrix, so it only uses the first MU upper bands,
c    and does not place nonzero values in the fillin bands.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 October 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than min(M,N)-1.
c
c    Input/output, integer SEED, a seed for the random number 
c    generator.
c
c    Output, double precision A(2*ML+MU+1,N), the R8GB matrix.
c
      implicit none

      integer ml
      integer mu
      integer n

      double precision a(2*ml+mu+1,n)
      double precision r8_uniform_01
      integer i
      integer j
      integer m
      integer row
      integer seed

      do j = 1, n
        do row = 1, 2*ml+mu+1
          i = row + j - ml - mu - 1
          if ( ml .lt. row .and. 1 .le. i .and. i .le. m ) then
            a(row,j) = r8_uniform_01 ( seed )
          else
            a(row,j) = 0.0D+00
          end if
        end do
      end do

      return
      end
      subroutine r8gb_sl ( n, ml, mu, a_lu, pivot, b, job )

c*********************************************************************72
c
cc R8GB_SL solves a system factored by R8GB_FA.
c
c  Discussion:
c
c    The R8GB storage format is for an M by N banded matrix, with lower 
c    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
c    extra superdiagonals, which may be required to store nonzero entries 
c    generated during Gaussian elimination.
c
c    The original M by N matrix is "collapsed" downward, so that diagonals
c    become rows of the storage array, while columns are preserved.  The
c    collapsed array is logically 2*ML+MU+1 by N.  
c
c    R8GB storage is used by LINPACK and LAPACK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 March 1999
c
c  Author:
c
c    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
c    This version by John Burkardt.
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
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than N-1.
c
c    Input, double precision A_LU(2*ML+MU+1,N), the LU factors from R8GB_FA.
c
c    Input, integer PIVOT(N), the pivot vector from R8GB_FA.
c
c    Input/output, double precision B(N).
c    On input, the right hand side vector.
c    On output, the solution.
c
c    Input, integer JOB.
c    0, solve A * x = b.
c    nonzero, solve A' * x = b.
c
      implicit none

      integer ml
      integer mu
      integer n

      double precision a_lu(2*ml+mu+1,n)
      double precision b(n)
      integer pivot(n)
      integer job
      integer k
      integer l
      integer la
      integer lb
      integer lm
      integer m
      double precision t

      m = mu + ml + 1
c
c  Solve A * x = b.
c
      if ( job .eq. 0 ) then
c
c  Solve L * Y = B.
c
        if ( 1 .le. ml ) then

          do k = 1, n-1

            lm = min ( ml, n-k )
            l = pivot(k)

            if ( l .ne. k ) then
              t    = b(l)
              b(l) = b(k)
              b(k) = t
            end if

            b(k+1:k+lm) = b(k+1:k+lm) + b(k) * a_lu(m+1:m+lm,k)

          end do
        end if
c
c  Solve U * X = Y.
c
        do k = n, 1, -1

          b(k) = b(k) / a_lu(m,k)
          lm = min ( k, m ) - 1
          la = m - lm
          lb = k - lm

          b(lb:lb+lm-1) = b(lb:lb+lm-1) - b(k) * a_lu(la:la+lm-1,k)

        end do
c
c  Solve A' * X = B.
c
      else
c
c  Solve U' * Y = B.
c
        do k = 1, n
          lm = min ( k, m ) - 1
          la = m - lm
          lb = k - lm
          b(k) = ( b(k) - sum ( a_lu(la:la+lm-1,k) * b(lb:lb+lm-1) ) ) 
     &      / a_lu(m,k)
        end do
c
c  Solve L' * X = Y.
c
        if ( 1 .le. ml ) then

          do k = n-1, 1, -1

            lm = min ( ml, n-k )
            b(k) = b(k) + sum ( a_lu(m+1:m+lm,k) * b(k+1:k+lm) )
            l = pivot(k)

            if ( l .ne. k ) then
              t    = b(l)
              b(l) = b(k)
              b(k) = t
            end if

          end do

        end if

      end if

      return
      end
      subroutine r8gb_to_r8s3 ( m, n, ml, mu, a, nz_num, isym, row, 
     &  col, b )

c*********************************************************************72
c
cc R8GB_TO_R8S3 copies an R8GB matrix to an R8S3 matrix.
c
c  Discussion:
c
c    The R8GB storage format is for an M by N banded matrix, with lower 
c    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
c    extra superdiagonals, which may be required to store nonzero entries 
c    generated during Gaussian elimination.
c
c    The original M by N matrix is "collapsed" downward, so that diagonals
c    become rows of the storage array, while columns are preserved.  The
c    collapsed array is logically 2*ML+MU+1 by N.  
c
c    R8GB storage is used by LINPACK and LAPACK.
c
c    LINPACK and LAPACK band storage requires that an extra ML
c    superdiagonals be supplied to allow for fillin during Gauss
c    elimination.  Even though a band matrix is described as
c    having an upper bandwidth of MU, it effectively has an
c    upper bandwidth of MU+ML.  This routine will copy nonzero
c    values it finds in these extra bands, so that both unfactored
c    and factored matrices can be handled.
c
c    The R8S3 storage format corresponds to the DLAP/SLAP Triad format.
c
c    The R8S3 storage format stores the row, column and value of each nonzero
c    entry of a sparse matrix.  The entries may be given in any order.  No
c    check is made for the erroneous case in which a given matrix entry is
c    specified more than once.
c
c    There is a symmetry option for square matrices.  If the symmetric storage
c    option is used, the format specifies that only nonzeroes on the diagonal
c    and lower triangle are stored.  However, this routine makes no attempt
c    to enforce this.  The only thing it does is to "reflect" any nonzero
c    offdiagonal value.  Moreover, no check is made for the erroneous case
c    in which both A(I,J) and A(J,I) are specified, but with different values.
c
c    This routine reorders the entries of A so that the first N entries
c    are exactly the diagonal entries of the matrix, in order.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 August 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrices.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrices.
c    N must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths of A1.
c    ML and MU must be nonnegative, and no greater than min(M,N)-1.
c
c    Input, double precision A(2*ML+MU+1,N), the R8GB matrix.
c
c    Input, integer NZ_NUM, the number of nonzero entries in A.
c    This number can be obtained by calling R8GB_NZ_NUM.
c
c    Output, integer ISYM, is 0 if the matrix is not symmetric,
c    and 1 if the matrix is symmetric.  If the matrix is symmetric, then
c    only the nonzeroes on the diagonal and in the lower triangle are stored.
c    For this routine, ISYM is always output 0.
c
c    Output, integer ROW(NZ_NUM), the row indices.
c
c    Output, integer COL(NZ_NUM), the column indices.
c
c    Output, double precision B(NZ_NUM), the R8S3 matrix.
c
      implicit none

      integer m
      integer ml
      integer mu
      integer n
      integer nz_num

      double precision a(2*ml+mu+1,n)
      double precision b(nz_num)
      integer col(nz_num)
      integer i
      integer isym
      integer j
      integer nz
      integer row(nz_num)

      isym = 0
      nz = 0

      do i = 1, m
        do j = 1, n
          if ( i - ml .le. j .and. j .le. i + mu + ml ) then
            if ( a(ml+mu+1+i-j,j) .ne. 0.0D+00 ) then

              if ( nz_num .le. nz ) then
                write ( *, '(a)' ) ' '
                write ( *, '(a)' ) 'R8GB_TO_R8S3 - Fatal error!'
                write ( *, '(a,i8)' ) '  NZ_NUM = ', nz_num
                write ( *, '(a)' ) 
     &            '  But the matrix has more nonzeros than that.'
                stop
              end if

              nz = nz + 1
              row(nz) = i
              col(nz) = j
              b(nz) = a(ml+mu+1+i-j,j)

            end if
          end if
        end do
      end do

      if ( nz .lt. nz_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8GB_TO_R8S3 - Warningc'
        write ( *, '(a,i8)' ) '  NZ_NUM = ', nz_num
        write ( *, '(a,i8)' ) '  But the number of nonzeros is ', nz
      end if

      return
      end
      subroutine r8gb_to_r8sp ( m, n, ml, mu, a, nz_num, row, col, b )

c*********************************************************************72
c
cc R8GB_TO_R8SP copies an R8GB matrix to an R8SP matrix.
c
c  Discussion:
c
c    The R8GB storage format is for an M by N banded matrix, with lower 
c    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
c    extra superdiagonals, which may be required to store nonzero entries 
c    generated during Gaussian elimination.
c
c    The original M by N matrix is "collapsed" downward, so that diagonals
c    become rows of the storage array, while columns are preserved.  The
c    collapsed array is logically 2*ML+MU+1 by N.  
c
c    R8GB storage is used by LINPACK and LAPACK.
c
c    LINPACK and LAPACK band storage requires that an extra ML
c    superdiagonals be supplied to allow for fillin during Gauss
c    elimination.  Even though a band matrix is described as
c    having an upper bandwidth of MU, it effectively has an
c    upper bandwidth of MU+ML.  This routine will copy nonzero
c    values it finds in these extra bands, so that both unfactored
c    and factored matrices can be handled.
c
c    The R8SP storage format stores the row, column and value of each nonzero
c    entry of a sparse matrix.
c
c    It is possible that a pair of indices (I,J) may occur more than
c    once.  Presumably, in this case, the intent is that the actual value
c    of A(I,J) is the sum of all such entries.  This is not a good thing
c    to do, but I seem to have come across this in MATLAB.
c
c    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP 
c    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    14 September 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrices.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrices.
c    N must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths of A1.
c    ML and MU must be nonnegative, and no greater than min(M,N)-1.
c
c    Input, double precision A(2*ML+MU+1,N), the R8GB matrix.
c
c    Input, integer NZ_NUM, the number of nonzero entries in A.
c    This number can be obtained by calling R8GB_NZ_NUM.
c
c    Output, integer ROW(NZ_NUM), the row indices.
c
c    Output, integer COL(NZ_NUM), the column indices.
c
c    Output, double precision B(NZ_NUM), the R8SP matrix.
c
      implicit none

      integer m
      integer ml
      integer mu
      integer n
      integer nz_num

      double precision a(2*ml+mu+1,n)
      double precision b(nz_num)
      integer col(nz_num)
      integer i
      integer j
      integer jhi
      integer jlo
      integer nz
      integer row(nz_num)

      nz = 0

      do i = 1, m

        jlo = max ( 1, i - ml )
        jhi = min ( n, i + mu + ml )

        do j = jlo, jhi

          if ( a(ml+mu+1+i-j,j) .eq. 0.0D+00 ) then
            cycle
          end if

          if ( nz_num .le. nz ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'R8GB_TO_R8SP - Fatal error!'
            write ( *, '(a,i8)' ) '  NZ_NUM = ', nz_num
            write ( *, '(a)' ) 
     &        '  But the matrix has more nonzeros than that.'
            stop
          end if

          nz = nz + 1
          row(nz) = i
          col(nz) = j
          b(nz) = a(ml+mu+1+i-j,j)

        end do
      end do

      if ( nz .lt. nz_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8GB_TO_R8SP - Warning!'
        write ( *, '(a,i8)' ) '  NZ_NUM = ', nz_num
        write ( *, '(a,i8)' ) '  But the number of nonzeros is ', nz
      end if

      return
      end
      subroutine r8gb_to_r8vec ( m, n, ml, mu, a, x )

c*********************************************************************72
c
cc R8GB_TO_R8VEC copies an R8GB matrix to an R8VEC.
c
c  Discussion:
c
c    In C++ and FORTRAN, this routine is not really needed.  In MATLAB,
c    a data item carries its dimensionality implicitly, and so cannot be
c    regarded sometimes as a vector and sometimes as an array.
c
c    The R8GB storage format is for an M by N banded matrix, with lower 
c    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
c    extra superdiagonals, which may be required to store nonzero entries 
c    generated during Gaussian elimination.
c
c    The original M by N matrix is "collapsed" downward, so that diagonals
c    become rows of the storage array, while columns are preserved.  The
c    collapsed array is logically 2*ML+MU+1 by N.  
c
c    R8GB storage is used by LINPACK and LAPACK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 March 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns in 
c    the array.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c
c    Input, double precision A(2*ML+MU+1,N), the array to be copied.
c
c    Output, double precision X((2*ML+MU+1)*N), the vector.
c
      implicit none

      integer m
      integer ml
      integer mu
      integer n

      double precision a(2*ml+mu+1,n)
      integer i
      integer ihi
      integer ilo
      integer j
      double precision x((2*ml+mu+1)*n)

      do j = 1, n

        ihi = min ( ml + mu, ml + mu + 1 - j )
        do i = 1, ihi
          x(i+(j-1)*(2*ml+mu+1)) = 0.0D+00
        end do

        ilo = max ( ihi + 1, 1 )
        ihi = min ( 2*ml+mu+1, ml+mu+m+1-j )
        do i = ilo, ihi
          x(i+(j-1)*(2*ml+mu+1)) = a(i,j)
        end do

        ilo = ihi + 1
        ihi = 2*ml+mu+1
        do i = ilo, ihi
          x(i+(j-1)*(2*ml+mu+1)) = 0.0D+00
        end do

      end do

      return
      end
      subroutine r8gb_to_r8ge ( m, n, ml, mu, a, b )

c*********************************************************************72
c
cc R8GB_TO_R8GE copies an R8GB matrix to an R8GE matrix.
c
c  Discussion:
c
c    The R8GB storage format is for an M by N banded matrix, with lower 
c    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
c    extra superdiagonals, which may be required to store nonzero entries 
c    generated during Gaussian elimination.
c
c    The original M by N matrix is "collapsed" downward, so that diagonals
c    become rows of the storage array, while columns are preserved.  The
c    collapsed array is logically 2*ML+MU+1 by N.  
c
c    R8GB storage is used by LINPACK and LAPACK.
c
c    LINPACK and LAPACK band storage requires that an extra ML
c    superdiagonals be supplied to allow for fillin during Gauss
c    elimination.  Even though a band matrix is described as
c    having an upper bandwidth of MU, it effectively has an
c    upper bandwidth of MU+ML.  This routine will copy nonzero
c    values it finds in these extra bands, so that both unfactored
c    and factored matrices can be handled.
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 March 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrices.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrices.
c    N must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths of A1.
c    ML and MU must be nonnegative, and no greater than min(M,N)-1.
c
c    Input, double precision A(2*ML+MU+1,N), the R8GB matrix.
c
c    Output, double precision B(M,N), the R8GE matrix.
c
      implicit none

      integer m
      integer ml
      integer mu
      integer n

      double precision a(2*ml+mu+1,n)
      double precision b(m,n)
      integer i
      integer j

      do i = 1, m
        do j = 1, n
          if ( i - ml .le. j .and. j .le. i + mu + ml ) then
            b(i,j) = a(ml+mu+1+i-j,j)
          else
            b(i,j) = 0.0D+00
          end if
        end do
      end do

      return
      end
      subroutine r8gb_trf ( m, n, ml, mu, a, pivot, info )

c*********************************************************************72
c
cc R8GB_TRF performs a LAPACK-style PLU factorization of an R8GB matrix.
c
c  Discussion:
c
c    The R8GB storage format is for an M by N banded matrix, with lower 
c    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
c    extra superdiagonals, which may be required to store nonzero entries 
c    generated during Gaussian elimination.
c
c    The original M by N matrix is "collapsed" downward, so that diagonals
c    become rows of the storage array, while columns are preserved.  The
c    collapsed array is logically 2*ML+MU+1 by N.  
c
c    R8GB storage is used by LINPACK and LAPACK.
c
c    This is a simplified, standalone version of the LAPACK
c    routine R8GBTRF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 January 1999
c
c  Author:
c
c    Original FORTRAN77 version by the LAPACK group.
c    This version by John Burkardt.
c
c  Reference:
c
c    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
c    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
c    Sven Hammarling, Alan McKenney, Danny Sorensen,
c    LAPACK User's Guide,
c    Second Edition,
c    SIAM, 1995.
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix A.  
c    0 .le. M.
c
c    Input, integer N, the number of columns of the matrix A.  
c    0 .le. N.
c
c    Input, integer ML, the number of subdiagonals within the
c    band of A.  0 .le. ML.
c
c    Input, integer MU, the number of superdiagonals within 
c    the band of A.  0 .le. MU.
c
c    Input/output, double precision A(2*ML+MU+1,N).  On input, the matrix A
c    in band storage, and on output, information about the PLU factorization.
c
c    Output, integer PIVOT(min(M,N)), the pivot indices;
c    for 1 .le. i .le. min(M,N), row i of the matrix was interchanged with
c    row IPIV(i).
c
c    Output, integer INFO, error flag.
c    = 0: successful exit;
c    .lt. 0: an input argument was illegal;
c    > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
c         has been completed, but the factor U is exactly
c         singular, and division by zero will occur if it is used
c         to solve a system of equations.
c
      implicit none

      integer ml
      integer mu
      integer m
      integer n

      double precision a(2*ml+mu+1,n)
      integer i
      integer info
      integer pivot(*)
      integer j
      integer jp
      integer ju
      integer k
      integer km
      integer kv
      double precision piv
c
      info = 0
c
c  KV is the number of superdiagonals in the factor U, allowing for fill-in.
c
      kv = mu + ml
c
c  Set fill-in elements in columns MU+2 to KV to zero.
c
      do j = mu + 2, min ( kv, n )
        do i = kv - j + 2, ml
          a(i,j) = 0.0D+00
        end do
      end do
c
c  JU is the index of the last column affected by the current stage
c  of the factorization.
c
      ju = 1

      do j = 1, min ( m, n )
c
c  Set the fill-in elements in column J+KV to zero.
c
        if ( j + kv .le. n ) then
          a(1:ml,j+kv) = 0.0D+00
        end if
c
c  Find the pivot and test for singularity.
c  KM is the number of subdiagonal elements in the current column.
c
        km = min ( ml, m-j )

        piv = abs ( a(kv+1,j) )
        jp = kv + 1

        do i = kv + 2, kv + km + 1
          if ( piv .lt. abs ( a(i,j) ) ) then
            piv = abs ( a(i,j) )
            jp = i
          end if
        end do

        jp = jp - kv

        pivot(j) = jp + j - 1

        if ( a(kv+jp,j) .ne. 0.0D+00 ) then

          ju = max ( ju, min ( j+mu+jp-1, n ) )
c
c  Apply interchange to columns J to JU.
c
          if ( jp .ne. 1 ) then

            do i = 0, ju - j
              call r8_swap ( a(kv+jp-i,j+i), a(kv+1-i,j+i) )
            end do

          end if
c
c  Compute the multipliers.
c
          if ( 0 .lt. km ) then

            a(kv+2:kv+km+1,j) = a(kv+2:kv+km+1,j) / a(kv+1,j)
c
c  Update the trailing submatrix within the band.
c
            if ( j .lt. ju ) then

              do k = 1, ju-j

                if ( a(kv+1-k,j+k) .ne. 0.0D+00 ) then

                  do i = 1, km
                    a(kv+i+1-k,j+k) = 
     &                a(kv+i+1-k,j+k) - a(kv+i+1,j) * a(kv+1-k,j+k)
                  end do

                end if

              end do

            end if

          end if

        else
c
c  If pivot is zero, set INFO to the index of the pivot
c  unless a zero pivot has already been found.
c
          if ( info .eq. 0 ) then
            info = j
          end if

        end if

      end do

      return
      end
      subroutine r8gb_trs ( n, ml, mu, nrhs, trans, a, pivot, b, info )

c*********************************************************************72
c
cc R8GB_TRS solves an R8GB linear system factored by R8GB_TRF.
c
c  Discussion:
c
c    The R8GB storage format is for an M by N banded matrix, with lower 
c    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
c    extra superdiagonals, which may be required to store nonzero entries 
c    generated during Gaussian elimination.
c
c    The original M by N matrix is "collapsed" downward, so that diagonals
c    become rows of the storage array, while columns are preserved.  The
c    collapsed array is logically 2*ML+MU+1 by N.  
c
c    R8GB storage is used by LINPACK and LAPACK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 January 1999
c
c  Author:
c
c    Original FORTRAN77 version by the LAPACK group.
c    This version by John Burkardt.
c
c  Reference:
c
c    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
c    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
c    Sven Hammarling, Alan McKenney, Danny Sorensen,
c    LAPACK User's Guide,
c    Second Edition,
c    SIAM, 1995.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix A.
c    N must be positive.
c
c    Input, integer ML, the number of subdiagonals within the 
c    band of A.  ML must be at least 0, and no greater than N - 1.
c
c    Input, integer MU, the number of superdiagonals within the 
c    band of A.  MU must be at least 0, and no greater than N - 1.
c
c    Input, integer NRHS, the number of right hand sides and the 
c    number of columns of the matrix B.  NRHS must be positive.
c
c    Input, character TRANS, specifies the form of the system.
c    'N':  A * x = b  (No transpose)
c    'T':  A'* X = B  (Transpose)
c    'C':  A'* X = B  (Conjugate transpose = Transpose)
c
c    Input, double precision A(2*ML+MU+1,N), the LU factorization of the 
c    band matrix A, computed by R8GB_TRF.  
c
c    Input, integer PIVOT(N), the pivot indices; for 1 .le. I .le. N, 
c    row I of the matrix was interchanged with row PIVOT(I).
c
c    Input/output, double precision B(N,NRHS),
c    On entry, the right hand side vectors B.
c    On exit, the solution vectors, X.
c
c    Output, integer INFO, error flag.
c    = 0:  successful exit
c    .lt. 0: if INFO = -K, the K-th argument had an illegal value
c
      implicit none

      integer ml
      integer mu
      integer n
      integer nrhs

      double precision a(2*ml+mu+1,n)
      double precision b(n,nrhs)
      integer i
      integer info
      integer pivot(*)
      integer j
      integer k
      integer kd
      integer l
      integer lm
      double precision temp
      character trans
c
c  Test the input parameters.
c
      info = 0

      if ( trans .ne. 'N' .and. trans .ne. 'n' .and. 
     &     trans .ne. 'T' .and. trans .ne. 't' .and. 
     &     trans .ne. 'C' .and. trans .ne. 'c' ) then
        info = -1
      else if ( n .le. 0 ) then
        info = -2
      else if ( ml .lt. 0 ) then
        info = -3
      else if ( mu .lt. 0 ) then
        info = -4
      else if ( nrhs .le. 0 ) then
        info = -5
      end if

      if ( info .ne. 0 ) then
        return
      end if

      kd = mu + ml + 1
c
c  Solve A * x = b.
c
c  Solve L * x = b, overwriting b with x.
c
c  L is represented as a product of permutations and unit lower
c  triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
c  where each transformation L(i) is a rank-one modification of
c  the identity matrix.
c
      if ( trans .eq. 'N' .or. trans .eq. 'n' ) then

        if ( 0 .lt. ml ) then

          do j = 1, n - 1

            lm = min ( ml, n-j )
            l = pivot(j)

            do i = 1, nrhs
              call r8_swap ( b(l,i), b(j,i) )
            end do

            do k = 1, nrhs
              if ( b(j,k) .ne. 0.0D+00 ) then
                b(j+1:j+lm,k) = b(j+1:j+lm,k) - a(kd+1:kd+lm,j) * b(j,k)
              end if
            end do

          end do

        end if
c
c  Solve U * x = b, overwriting b with x.
c
        do i = 1, nrhs

          do j = n, 1, -1
            if ( b(j,i) .ne. 0.0D+00 ) then
              l = ml + mu + 1 - j
              b(j,i) = b(j,i) / a(ml+mu+1,j)
              do k = j - 1, max ( 1, j - ml - mu ), -1
                b(k,i) = b(k,i) - a(l+k,j) * b(j,i)
              end do
            end if
          end do

        end do

      else
c
c  Solve A' * x = b.
c
c  Solve U' * x = b, overwriting b with x.
c
        do i = 1, nrhs

          do j = 1, n
            temp = b(j,i)
            l = ml + mu + 1 - j
            do k = max ( 1, j - ml - mu ), j - 1
              temp = temp - a(l+k,j) * b(k,i)
            end do
            b(j,i) = temp / a(ml+mu+1,j)
          end do

        end do
c
c  Solve L' * x = b, overwriting b with x.
c
        if ( 0 .lt. ml ) then

          do j = n - 1, 1, -1

            lm = min ( ml, n-j )

            do k = 1, nrhs
              b(j,k) = b(j,k) - sum ( b(j+1:j+lm,k) * a(kd+1:kd+lm,j) )
            end do

            l = pivot(j)

            do i = 1, nrhs
              call r8_swap ( b(l,i), b(j,i) )
            end do

          end do

        end if

      end if

      return
      end
      subroutine r8gb_vxm ( m, n, ml, mu, a, x, b )

c*********************************************************************72
c
cc R8GB_VXM multiplies an R8VEC by an R8GB matrix.
c
c  Discussion:
c
c    The R8GB storage format is for an M by N banded matrix, with lower 
c    bandwidth ML and upper bandwidth MU.  Storage includes room for ML 
c    extra superdiagonals, which may be required to store nonzero entries 
c    generated during Gaussian elimination.
c
c    The original M by N matrix is "collapsed" downward, so that diagonals
c    become rows of the storage array, while columns are preserved.  The
c    collapsed array is logically 2*ML+MU+1 by N.  
c
c    R8GB storage is used by LINPACK and LAPACK.
c
c    LINPACK and LAPACK storage of general band matrices requires
c    an extra ML upper diagonals for possible fill in entries during
c    Gauss elimination.  This routine does not access any entries
c    in the fill in diagonals, because it assumes that the matrix
c    has NOT had Gauss elimination applied to it.  If the matrix
c    has been Gauss eliminated, then the routine R8GB_MU must be
c    used instead.
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
c    Input, integer M, the number of rows of the matrix.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c    ML and MU must be nonnegative, and no greater than min(M,N)-1.
c
c    Input, double precision A(2*ML+MU+1,N), the R8GB matrix.
c
c    Input, double precision X(M), the vector to be multiplied by A.
c
c    Output, double precision B(N), the product X*A.
c
      implicit none

      integer ml
      integer mu
      integer m
      integer n

      double precision a(2*ml+mu+1,n)
      double precision b(n)
      integer i
      integer ihi
      integer ilo
      integer j
      double precision x(m)

      do j = 1, n
        b(j) = 0.0D+00
        ilo = max ( 1, j - mu )
        ihi = min ( m, j + ml )
        do i = ilo, ihi
          b(j) = b(j) + x(i) * a(i-j+ml+mu+1,j)
        end do
      end do

      return
      end
      subroutine r8gd_error ( n, ndiag, ierror )

c*********************************************************************72
c
cc R8GD_ERROR checks the dimensions of an R8GD matrix.
c
c  Discussion:
c
c    The R8GD storage format is suitable for matrices whose only nonzero entries
c    occur along a few diagonals, but for which these diagonals are not all
c    close enough to the main diagonal for band storage to be efficient.
c
c    In that case, we assign the main diagonal the offset value 0.
c    Each successive superdiagonal gets an offset value 1 higher, until
c    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
c    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
c
c    Now, assuming that only a few of these diagonals contain nonzeros,
c    then for the I-th diagonal to be saved, we stored its offset in
c    OFFSET(I), and its entries in column I of the matrix.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 October 1998
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
c    Input, integer NDIAG, the number of diagonals of the matrix
c    that are stored in the array.
c    NDIAG must be at least 1, and no more than 2 * N - 1.
c
c    Output, integer IERROR, reports whether any errors were 
c    detected.
c    IERROR is set to 0 before the checks are made, and then:
c    IERROR = IERROR + 1 if N is illegal;
c    IERROR = IERROR + 2 if NDIAG is illegal.
c
      implicit none

      integer ierror
      integer n
      integer ndiag

      ierror = 0

      if ( n .lt. 1 ) then
        ierror = ierror + 1
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) 'R8GD_ERROR - Illegal N = ', n
      end if

      if ( ndiag .lt. 1 .or. 2 * n - 1 .lt. ndiag ) then
        ierror = ierror + 2
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) 'R8GD_ERROR - Illegal NDIAG = ', ndiag
      end if

      return
      end
      subroutine r8gd_indicator ( n, ndiag, offset, a )

c*********************************************************************72
c
cc R8GD_INDICATOR sets up an R8GD indicator matrix.
c
c  Discussion:
c
c    The "indicator matrix" simply has a value like I*10+J at every
c    entry of a dense matrix, or at every entry of a compressed storage
c    matrix for which storage is allocated. 
c
c    The R8GD storage format is suitable for matrices whose only nonzero entries
c    occur along a few diagonals, but for which these diagonals are not all
c    close enough to the main diagonal for band storage to be efficient.
c
c    In that case, we assign the main diagonal the offset value 0.
c    Each successive superdiagonal gets an offset value 1 higher, until
c    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
c    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
c
c    Now, assuming that only a few of these diagonals contain nonzeros,
c    then for the I-th diagonal to be saved, we stored its offset in
c    OFFSET(I), and its entries in column I of the matrix.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 January 2004
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
c    Input, integer NDIAG, the number of diagonals of the matrix
c    that are stored in the array.
c    NDIAG must be at least 1, and no more than 2 * N - 1.
c
c    Input, integer OFFSET(NDIAG), the offsets for the diagonal
c    storage.
c
c    Output, double precision A(N,NDIAG), the R8GD matrix.
c
      implicit none

      integer n
      integer ndiag

      double precision a(n,ndiag)
      integer diag
      integer fac
      integer i
      integer i4_log_10
      integer j
      integer offset(ndiag)

      fac = 10 ** ( i4_log_10 ( n ) + 1 )

      do i = 1, n
        do diag = 1, ndiag
          j = i + offset(diag)
          if ( 1 .le. j .and. j .le. n ) then
            a(i,diag) = dble ( fac * i + j )
          else
            a(i,diag) = 0.0D+00
          end if
        end do
      end do

      return
      end
      subroutine r8gd_mxv ( n, ndiag, offset, a, x, b )

c*********************************************************************72
c
cc R8GD_MXV multiplies an R8GD matrix by an R8VEC.
c
c  Discussion:
c
c    The R8GD storage format is suitable for matrices whose only nonzero entries
c    occur along a few diagonals, but for which these diagonals are not all
c    close enough to the main diagonal for band storage to be efficient.
c
c    In that case, we assign the main diagonal the offset value 0.
c    Each successive superdiagonal gets an offset value 1 higher, until
c    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
c    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
c
c    Now, assuming that only a few of these diagonals contain nonzeros,
c    then for the I-th diagonal to be saved, we stored its offset in
c    OFFSET(I), and its entries in column I of the matrix.  
c
c  Example:
c
c    The "offset" value is printed near the first entry of each diagonal
c    in the original matrix, and above the columns in the new matrix.
c
c    Original matrix               New Matrix
c
c      0    1   2   3   4   5        -3  -2   0   1   3   5
c       \    \   \   \   \   \
c        11  12   0  14   0  16      --  --  11  12  14  16
c   -1 =  0  22  23   0  25   0      --  --  22  23  25  --
c   -2 = 31   0  33  34   0  36      --  31  33  34  36  --
c   -3 = 41  42   0  44  45   0      41  42  44  45  --  --
c   -4 =  0  52  53   0  55  56      52  53  55  56  --  --
c   -5 =  0   0  63  64  65  66      63  64  66  --  --  --
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 January 2004
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
c    Input, integer NDIAG, the number of diagonals of the matrix
c    that are stored in the array.
c    NDIAG must be at least 1, and no more than 2 * N - 1.
c
c    Input, integer OFFSET(NDIAG), the offsets for the diagonal 
c    storage.
c
c    Input, double precision A(N,NDIAG), the R8GD matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(N), the product A * x.
c
      implicit none

      integer n
      integer ndiag

      double precision a(n,ndiag)
      double precision b(n)
      integer diag
      integer i
      integer j
      integer offset(ndiag)
      double precision x(n)

      b(1:n) = 0.0D+00

      do i = 1, n
        do diag = 1, ndiag
          j = i + offset(diag)
          if ( 1 .le. j .and. j .le. n ) then
            b(i) = b(i) + a(i,diag) * x(j)
          end if
        end do
      end do

      return
      end
      subroutine r8gd_print ( n, ndiag, offset, a, title )

c*********************************************************************72
c
cc R8GD_PRINT prints an R8GD matrix.
c
c  Discussion:
c
c    The R8GD storage format is suitable for matrices whose only nonzero entries
c    occur along a few diagonals, but for which these diagonals are not all
c    close enough to the main diagonal for band storage to be efficient.
c
c    In that case, we assign the main diagonal the offset value 0.
c    Each successive superdiagonal gets an offset value 1 higher, until
c    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
c    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
c
c    Now, assuming that only a few of these diagonals contain nonzeros,
c    then for the I-th diagonal to be saved, we stored its offset in
c    OFFSET(I), and its entries in column I of the matrix.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 May 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input, integer NDIAG, the number of diagonals of the matrix
c    that are stored in the array.
c    NDIAG must be at least 1, and no more than 2 * N - 1.
c
c    Input, integer OFFSET(NDIAG), the offsets for the 
c    diagonal storage.
c
c    Input, double precision A(N,NDIAG), the R8GD matrix.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer n
      integer ndiag

      double precision a(n,ndiag)
      integer offset(ndiag)
      character ( len = * )  title

      call r8gd_print_some ( n, ndiag, offset, a, 1, 1, n, n, title )

      return
      end
      subroutine r8gd_print_some ( n, ndiag, offset, a, ilo, jlo, ihi, 
     &  jhi, title )

c*********************************************************************72
c
cc R8GD_PRINT_SOME prints some of an R8GD matrix.
c
c  Discussion:
c
c    The R8GD storage format is suitable for matrices whose only nonzero entries
c    occur along a few diagonals, but for which these diagonals are not all
c    close enough to the main diagonal for band storage to be efficient.
c
c    In that case, we assign the main diagonal the offset value 0.
c    Each successive superdiagonal gets an offset value 1 higher, until
c    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
c    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
c
c    Now, assuming that only a few of these diagonals contain nonzeros,
c    then for the I-th diagonal to be saved, we stored its offset in
c    OFFSET(I), and its entries in column I of the matrix.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 January 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input, integer NDIAG, the number of diagonals of the matrix
c    that are stored in the array.
c    NDIAG must be at least 1, and no more than 2 * N - 1.
c
c    Input, integer OFFSET(NDIAG), the offsets for the diagonal 
c    storage.
c
c    Input, double precision A(N,NDIAG), the R8GD matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer, parameter :: incx = 5
      integer n
      integer ndiag

      double precision a(n,ndiag)
      double precision aij
      character ( len = 14 ) ctemp(incx)
      logical r8_is_int
      integer diag
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
      integer off
      integer offset(ndiag)
      character ( len = * ) title

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
        i2hi = min ( ihi, n )

        do i = i2lo, i2hi
c
c  Print out (up to) 5 entries in row I, that lie in the current strip.
c
          do j2 = 1, inc

            j = j2lo - 1 + j2

            aij = 0.0D+00
            off = j - i
            do diag = 1, ndiag
              if ( off .eq. offset(diag) ) then
                aij = a(i,diag)
              end if
            end do

            if ( r8_is_int ( aij ) ) then
              write ( ctemp(j2), '(f8.0,6x)' ) aij
            else
              write ( ctemp(j2), '(g14.6)' ) aij
            end if

          end do

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine r8gd_random ( n, ndiag, offset, seed, a )

c*********************************************************************72
c
cc R8GD_RANDOM randomizes an R8GD matrix.
c
c  Discussion:
c
c    The R8GD storage format is suitable for matrices whose only nonzero entries
c    occur along a few diagonals, but for which these diagonals are not all
c    close enough to the main diagonal for band storage to be efficient.
c
c    In that case, we assign the main diagonal the offset value 0.
c    Each successive superdiagonal gets an offset value 1 higher, until
c    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
c    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
c
c    Now, assuming that only a few of these diagonals contain nonzeros,
c    then for the I-th diagonal to be saved, we stored its offset in
c    OFFSET(I), and its entries in column I of the matrix.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 October 2004
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
c    Input, integer NDIAG, the number of diagonals of the matrix
c    that are stored in the array.
c    NDIAG must be at least 1, and no more than 2 * N - 1.
c
c    Input, integer OFFSET(NDIAG), the offsets for the diagonal 
c    storage.
c
c    Input/output, integer SEED, a seed for the random number 
c    generator.
c
c    Output, double precision A(N,NDIAG), the R8GD matrix.
c
      implicit none

      integer n
      integer ndiag

      double precision a(n,ndiag)
      double precision r8_uniform_01
      integer diag
      integer i
      integer j
      integer offset(ndiag)
      integer seed

      do i = 1, n
        do diag = 1, ndiag
          j = i + offset(diag)
          if ( 1 .le. j .and. j .le. n ) then
            a(i,diag) = r8_uniform_01 ( seed )
          else
            a(i,diag) = 0.0D+00
          end if
        end do
      end do

      return
      end
      subroutine r8gd_to_r8ge ( n, ndiag, offset, a, b )

c*********************************************************************72
c
cc R8GD_TO_R8GE copies an R8GD matrix to an R8GE matrix.
c
c  Discussion:
c
c    The R8GD storage format is suitable for matrices whose only nonzero entries
c    occur along a few diagonals, but for which these diagonals are not all
c    close enough to the main diagonal for band storage to be efficient.
c
c    In that case, we assign the main diagonal the offset value 0.
c    Each successive superdiagonal gets an offset value 1 higher, until
c    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
c    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
c
c    Now, assuming that only a few of these diagonals contain nonzeros,
c    then for the I-th diagonal to be saved, we stored its offset in
c    OFFSET(I), and its entries in column I of the matrix.  
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 January 2004
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
c    Input, integer NDIAG, the number of diagonals of the matrix
c    that are stored in the array.
c    NDIAG must be at least 1, and no more than 2 * N - 1.
c
c    Input, integer OFFSET(NDIAG), the offsets for the diagonal 
c    storage.
c
c    Input, double precision A(N,NDIAG), the R8GD matrix.
c
c    Output, double precision B(N,N), the R8GE matrix.
c
      implicit none

      integer n
      integer ndiag

      double precision a(n,ndiag)
      double precision b(n,n)
      integer diag
      integer i
      integer j
      integer offset(ndiag)

      b(1:n,1:n) = 0.0D+00

      do i = 1, n
        do diag = 1, ndiag
          j = i + offset(diag)
          if ( 1 .le. j .and. j .le. n ) then
            b(i,j) = a(i,diag)
          end if
        end do
      end do

      return
      end
      subroutine r8gd_vxm ( n, ndiag, offset, a, x, b )

c*********************************************************************72
c
cc R8GD_VXM multiplies an R8VEC by an R8GD matrix.
c
c  Discussion:
c
c    The R8GD storage format is suitable for matrices whose only nonzero entries
c    occur along a few diagonals, but for which these diagonals are not all
c    close enough to the main diagonal for band storage to be efficient.
c
c    In that case, we assign the main diagonal the offset value 0.
c    Each successive superdiagonal gets an offset value 1 higher, until
c    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
c    Similarly, the subdiagonals are assigned offsets of -1 through -(N-1).
c
c    Now, assuming that only a few of these diagonals contain nonzeros,
c    then for the I-th diagonal to be saved, we stored its offset in
c    OFFSET(I), and its entries in column I of the matrix.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 January 2004
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
c    Input, integer NDIAG, the number of diagonals of the matrix
c    that are stored in the array.
c    NDIAG must be at least 1, and no more than 2 * N - 1.
c
c    Input, integer OFFSET(NDIAG), the offsets for the diagonal 
c    storage.
c
c    Input, double precision A(N,NDIAG), the R8GD matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(N), the product X*A.
c
      implicit none

      integer n
      integer ndiag

      double precision a(n,ndiag)
      double precision b(n)
      integer diag
      integer i
      integer j
      integer offset(ndiag)
      double precision x(n)

      b(1:n) = 0.0D+00

      do i = 1, n
        do diag = 1, ndiag
          j = i + offset(diag)
          if ( 1 .le. j .and. j .le. n ) then
            b(j) = b(j) + x(i) * a(i,diag)
          end if
        end do
      end do

      return
      end
      subroutine r8ge_det ( n, a_lu, pivot, det )

c*********************************************************************72
c
cc R8GE_DET computes the determinant of a matrix factored by R8GE_FA or R8GE_TRF.
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
c    Input, double precision A_LU(N,N), the LU factors from R8GE_FA or R8GE_TRF.
c
c    Input, integer PIVOT(N), as computed by R8GE_FA or R8GE_TRF.
c
c    Output, double precision DET, the determinant of the matrix.
c
      implicit none

      integer n

      double precision a_lu(n,n)
      double precision det
      integer i
      integer pivot(n)

      det = 1.0D+00

      do i = 1, n
        det = det * a_lu(i,i)
        if ( pivot(i) .ne. i ) then
          det = - det
        end if
      end do

      return
      end
      subroutine r8ge_dilu ( m, n, a, d )

c*********************************************************************72
c
cc R8GE_DILU produces the diagonal incomplete LU factor of an R8GE matrix.
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
c    Input, integer M, the number of rows in A.
c
c    Input, integer N, the number of columns in A.
c
c    Input, double precision A(M,N), the R8GE matrix.
c
c    Output, double precision D(M), the DILU factor.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision d(m)
      integer i
      integer j

      do i = 1, m
        if ( i .le. n ) then
          d(i) = a(i,i)
        else
          d(i) = 0.0D+00
        end if
      end do

      do i = 1, min ( m, n )
        d(i) = 1.0D+00 / d(i)
        do j = i+1, min ( m, n )
          d(j) = d(j) - a(j,i) * d(i) * a(i,j)
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
      subroutine r8ge_fs ( n, a, b, info )

c*********************************************************************72
c
cc R8GE_FS factors and solves an R8GE system.
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
c    This routine does not save the LU factors of the matrix, and hence cannot
c    be used to efficiently solve multiple linear systems, or even to
c    factor A at one time, and solve a single linear system at a later time.
c
c    This routine uses partial pivoting, but no pivot vector is required.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 November 2008
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
c    Input/output, double precision A(N,N).
c    On input, A is the coefficient matrix of the linear system.
c    On output, A is in unit upper triangular form, and
c    represents the U factor of an LU factorization of the
c    original coefficient matrix.
c
c    Input/output, double precision B(N).
c    On input, B is the right hand side of the linear system.
c    On output, B is the solution of the linear system.
c
c    Output, integer INFO, singularity flag.
c    0, no singularity detected.
c    nonzero, the factorization failed on the INFO-th step.
c
      implicit none

      integer n

      double precision a(n,n)
      double precision b(n)
      integer i
      integer info
      integer ipiv
      integer j
      integer jcol
      double precision piv
      double precision temp

      info = 0

      do jcol = 1, n
c
c  Find the maximum element in column I.
c
        piv = abs ( a(jcol,jcol) )
        ipiv = jcol
        do i = jcol + 1, n
          if ( piv .lt. abs ( a(i,jcol) ) ) then
            piv = abs ( a(i,jcol) )
            ipiv = i
          end if
        end do

        if ( piv .eq. 0.0D+00 ) then
          info = jcol
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8GE_FS - Fatal error!'
          write ( *, '(a,i8)' ) '  Zero pivot on step ', info
          return
        end if
c
c  Switch rows JCOL and IPIV, and B.
c
        if ( jcol .ne. ipiv ) then

          do j = 1, n
            temp = a(jcol,j)
            a(jcol,j) = a(ipiv,j)
            a(ipiv,j) = temp
          end do

          temp = b(jcol)
          b(jcol) = b(ipiv)
          b(ipiv) = temp

        end if
c
c  Scale the pivot row.
c
        do j = jcol + 1, n
          a(jcol,j) = a(jcol,j) / a(jcol,jcol)
        end do
        b(jcol) = b(jcol) / a(jcol,jcol)
        a(jcol,jcol) = 1.0D+00
c
c  Use the pivot row to eliminate lower entries in that column.
c
        do i = jcol + 1, n
          if ( a(i,jcol) .ne. 0.0D+00 ) then
            temp = - a(i,jcol)
            a(i,jcol) = 0.0D+00
            do j = jcol + 1, n
              a(i,j) = a(i,j) + temp * a(jcol,j)
            end do
            b(i) = b(i) + temp * b(jcol)
          end if
        end do

      end do
c
c  Back solve.
c
      do j = n, 2, -1
        do i = 1, j - 1
          b(i) = b(i) - a(i,j) * b(j)
        end do
      end do

      return
      end
      subroutine r8ge_fss ( n, a, nb, b, info )

c*********************************************************************72
c
cc R8GE_FSS factors and solves multiple R8GE systems.
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
c    This routine does not save the LU factors of the matrix, and hence cannot
c    be used to efficiently solve multiple linear systems, or even to
c    factor A at one time, and solve a single linear system at a later time.
c
c    This routine uses partial pivoting, but no pivot vector is required.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 June 2009
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
c    Input/output, double precision A(N,N).
c    On input, A is the coefficient matrix of the linear system.
c    On output, A is in unit upper triangular form, and
c    represents the U factor of an LU factorization of the
c    original coefficient matrix.
c
c    Input, integer NB, the number of right hand sides.
c
c    Input/output, double precision B(N,NB).
c    On input, B is the right hand side of the linear system.
c    On output, B is the solution of the linear system.
c
c    Output, integer INFO, singularity flag.
c    0, no singularity detected.
c    nonzero, the factorization failed on the INFO-th step.
c
      implicit none

      integer n
      integer nb

      double precision a(n,n)
      double precision b(n,nb)
      integer i
      integer info
      integer ipiv
      integer j
      integer jcol
      integer k
      double precision piv
      double precision temp

      info = 0

      do jcol = 1, n
c
c  Find the maximum element in column I.
c
        piv = abs ( a(jcol,jcol) )
        ipiv = jcol
        do i = jcol + 1, n
          if ( piv .lt. abs ( a(i,jcol) ) ) then
            piv = abs ( a(i,jcol) )
            ipiv = i
          end if
        end do

        if ( piv .eq. 0.0D+00 ) then
          info = jcol
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8GE_FSS - Fatal error!'
          write ( *, '(a,i8)' ) '  Zero pivot on step ', info
          return
        end if
c
c  Switch rows JCOL and IPIV, and B.
c
        if ( jcol .ne. ipiv ) then

          do j = 1, n
            temp = a(jcol,j)
            a(jcol,j) = a(ipiv,j)
            a(ipiv,j) = temp
          end do

          do j = 1, nb
            temp = b(jcol,j)
            b(jcol,j) = b(ipiv,j)
            b(ipiv,j) = temp
          end do

        end if
c
c  Scale the pivot row.
c
        do j = jcol + 1, n
          a(jcol,j) = a(jcol,j) / a(jcol,jcol)
        end do
        do j = 1, nb
          b(jcol,j) = b(jcol,j) / a(jcol,jcol)
        end do
        a(jcol,jcol) = 1.0D+00
c
c  Use the pivot row to eliminate lower entries in that column.
c
        do i = jcol + 1, n
          if ( a(i,jcol) .ne. 0.0D+00 ) then
            temp = - a(i,jcol)
            a(i,jcol) = 0.0D+00
            do j = jcol + 1, n
              a(i,j) = a(i,j) + temp * a(jcol,j)
            end do
            do j = 1, nb
              b(i,j) = b(i,j) + temp * b(jcol,j)
            end do
          end if
        end do

      end do
c
c  Back solve.
c
      do k = n, 2, -1
        do i = 1, k - 1
          do j = 1, nb
            b(i,j) = b(i,j) - a(i,k) * b(k,j)
          end do
        end do
      end do

      return
      end
      subroutine r8ge_identity ( n, a )

c*********************************************************************72
c
cc R8GE_IDENTITY copies the identity matrix to an R8GE matrix.
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
c    Input, integer N, the order of A.
c
c    Output, double precision A(N,N), the N by N identity matrix.
c
      implicit none

      integer n

      double precision a(n,n)
      integer i
      integer j

      do i = 1, n
        do j = 1, n
          if ( i .eq. j ) then
            a(i,j) = 1.0D+00
          else
            a(i,j) = 0.0D+00
          end if
        end do
      end do

      return
      end
      subroutine r8ge_ilu ( m, n, a, l, u )

c*********************************************************************72
c
cc R8GE_ILU produces the incomplete LU factors of an R8GE matrix.
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
c    The incomplete LU factors of the M by N matrix A are:
c
c      L, an M by M unit lower triangular matrix,
c      U, an M by N upper triangular matrix
c
c    with the property that L and U are computed in the same way as
c    the usual LU factors, except that, whenever an off diagonal element
c    of the original matrix is zero, then the corresponding value of
c    U is forced to be zero.
c
c    This condition means that it is no longer the case that A = L*U.
c
c    On the other hand, L and U will have a simple sparsity structure
c    related to that of A.  The incomplete LU factorization is generally
c    used as a preconditioner in iterative schemes applied to sparse
c    matrices.  It is presented here merely for illustration.
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
c    Input, integer M, the number of rows in A.
c
c    Input, integer N, the number of columns in A.
c
c    Input, double precision A(M,N), the R8GE matrix.
c
c    Output, double precision L(M,M), the M by M unit lower triangular factor.
c
c    Output, double precision U(M,N), the M by N upper triangular factor.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer i
      integer j
      integer k
      double precision l(m,m)
      double precision u(m,n)
c
c  Initialize:
c
c    L := M by M Identity
c    U := A
c
      call r8ge_identity ( m, l )

      call r8mat_copy ( m, n, a, u )

      do j = 1, min ( m - 1, n )
c
c  Zero out the entries in column J, from row J+1 to M.
c
        do i = j + 1, m

          if ( u(i,j) .ne. 0.0D+00 ) then

            l(i,j) = u(i,j) / u(j,j)
            u(i,j) = 0.0D+00

            do k = j + 1, n
              if ( u(i,k) .ne. 0.0D+00 ) then
                u(i,k) = u(i,k) - l(i,j) * u(j,k)
              end if
            end do

          end if

        end do

      end do
 
      return
      end
      subroutine r8ge_indicator ( m, n, a )

c*********************************************************************72
c
cc R8GE_INDICATOR sets up an R8GE indicator matrix.
c
c  Discussion:
c
c    The "indicator matrix" simply has a value like I*10+J at every
c    entry of a dense matrix, or at every entry of a compressed storage
c    matrix for which storage is allocated. 
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
c    Input, integer M, the number of rows of the matrix.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Output, double precision A(M,N), the R8GE matrix.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer fac
      integer i
      integer i4_log_10
      integer j

      fac = 10 ** ( i4_log_10 ( n ) + 1 )

      do i = 1, m
        do j = 1, n
          a(i,j) = dble ( fac * i + j )
        end do
      end do

      return
      end
      subroutine r8ge_mxm ( n1, n2, n3, a, b, c )

c*********************************************************************72
c
cc R8GE_MXM multiplies two R8GE matrices.
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
c    Input, integer N1, N2, N3, the order of the matrices.
c    N must be positive.
c
c    Input, double precision A(N1,N2), B(N2,N3), the R8GE factor matrices.
c
c    Output, double precision C(N1,N3), the R8GE product matrix.
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
      subroutine r8ge_mxv ( m, n, a, x, b )

c*********************************************************************72
c
cc R8GE_MXV multiplies an R8GE matrix times a vector.
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
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 November 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input, double precision A(M,N), the R8GE matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(M), the product A * x.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision b(m)
      integer i
      integer j
      double precision x(n)

      do i = 1, m
        b(i) = 0.0D+00
        do j = 1, n
          b(i) = b(i) + a(i,j) * x(j)
        end do
      end do

      return
      end
      subroutine r8ge_print ( m, n, a, title )

c*********************************************************************72
c
cc R8GE_PRINT prints an R8GE matrix.
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
c    Input, integer M, the number of rows of the matrix.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input, double precision A(M,N), the R8GE matrix.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      character * ( * ) title

      call r8ge_print_some ( m, n, a, 1, 1, m, n, title )

      return
      end
      subroutine r8ge_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

c*********************************************************************72
c
cc R8GE_PRINT_SOME prints some of an R8GE matrix.
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
c    Input, integer M, the number of rows of the matrix.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input, double precision A(M,N), the R8GE matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character * ( * ) TITLE, a title.
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
      write ( *, '(a)' ) title
!
!  Print the columns of the matrix, in strips of 5.
!
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

        write ( *, '(''  Col:  '',5a14)' ) ( ctemp(j2), j2 = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) '  ---'
!
!  Determine the range of the rows in this strip.
!
        i2lo = max ( ilo, 1 )
        i2hi = min ( ihi, m )

        do i = i2lo, i2hi
!
!  Print out (up to) 5 entries in row I, that lie in the current strip.
!
          do j2 = 1, inc

            j = j2lo - 1 + j2

            write ( ctemp(j2), '(g14.6)' ) a(i,j)

          end do

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine r8ge_random ( m, n, seed, a )

c*********************************************************************72
c
cc R8GE_RANDOM randomizes an R8GE matrix.
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
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 November 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input/output, integer SEED, a seed for the random number
c    generator.
c
c    Output, double precision A(M,N), the R8GE matrix.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer i
      integer j
      double precision r8_uniform_01
      integer seed

      do j = 1, n
        do i = 1, m
          a(i,j) = r8_uniform_01 ( seed )
        end do
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
      subroutine r8ge_to_r8vec ( m, n, a, x )

c*********************************************************************72
c
cc R8GE_TO_R8VEC copies an R8GE matrix to an R8VEC.
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
c    In C++ and FORTRAN, this routine is not really needed.  In MATLAB,
c    a data item carries its dimensionality implicitly, and so cannot be
c    regarded sometimes as a vector and sometimes as an array.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns in 
c    the array.
c
c    Input, double precision A(M,N), the array to be copied.
c
c    Output, double precision X(M*N), the vector.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer i
      integer j
      integer k
      double precision x(m*n)

      k = 0
      do j = 1, n
        do i = 1, m
          k = k + 1
          x(k) = a(i,j)
        end do
      end do

      return
      end
      subroutine r8lt_det ( n, a, det )

c*********************************************************************72
c
cc R8LT_DET computes the determinant of an R8LT matrix.
c
c  Discussion:
c
c    The R8LT storage format is used for an M by N lower triangular matrix,
c    and sets aside storage even for the entries that must be zero.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 May 2003
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
c    Input, double precision A(N,N), the R8LT matrix.
c
c    Output, double precision DET, the determinant of the matrix.
c
      implicit none

      integer n

      double precision a(n,n)
      double precision det
      double precision diag(n)

      call r8mat_diag_get_vector ( n, a, diag )

      det = product ( diag(1:n) )

      return
      end
      subroutine r8lt_indicator ( m, n, a )

c*********************************************************************72
c
cc R8LT_INDICATOR sets up an R8LT indicator matrix.
c
c  Discussion:
c
c    The "indicator matrix" simply has a value like I*10+J at every
c    entry of a dense matrix, or at every entry of a compressed storage
c    matrix for which storage is allocated. 
c
c    The R8LT storage format is used for an M by N lower triangular matrix,
c    and sets aside storage even for the entries that must be zero.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    1 February 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns of 
c    the matrix.  M and N must be positive.
c
c    Output, double precision A(M,N), the R8LT matrix.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer fac
      integer i
      integer i4_log_10
      integer j

      fac = 10 ** ( i4_log_10 ( n ) + 1 )

      do i = 1, m
        do j = 1, min ( i, n )
          a(i,j) = dble ( fac * i + j )
        end do
        do j = i+1, n
          a(i,j) = 0.0D+00
        end do
      end do

      return
      end
      subroutine r8lt_inverse ( n, a )

c*********************************************************************72
c
cc R8LT_INVERSE computes the inverse of an R8LT matrix.
c
c  Discussion:
c
c    The R8LT storage format is used for an M by N lower triangular matrix,
c    and sets aside storage even for the entries that must be zero.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 May 2003
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms,
c    Second edition,
c    Academic Press, 1978,
c    ISBN 0-12-519260-6
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input/output, double precision A(N,N).
c
c    On input, the lower triangular matrix to be inverted.
c    On output, the inverse of the lower triangular matrix.
c
      implicit none

      integer n

      double precision a(n,n)
      integer i
      integer j
c
c  Check.
c
      do i = 1, n
        if ( a(i,i) .eq. 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8LT_INVERSE - Fatal error!'
          write ( *, '(a)' ) '  Zero diagonal element.'
          stop
        end if
      end do

      do j = 1, n

        do i = 1, n

          if ( i .lt. j ) then

            a(i,j) = 0.0D+00

          else if ( i .eq. j ) then

            a(i,j) = 1.0D+00 / a(i,j)

          else if ( j .lt. i ) then

            a(i,j) = - sum ( a(i,j:i-1) * a(j:i-1,j) ) / a(i,i)

          end if

        end do
      end do

      return
      end
      subroutine r8lt_mxm ( n, a, b, c )

c*********************************************************************72
c
cc R8LT_MXM multiplies two R8LT matrices.
c
c  Discussion:
c
c    The R8LT storage format is used for an M by N lower triangular matrix,
c    and sets aside storage even for the entries that must be zero.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 May 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrices.
c    N must be positive.
c
c    Input, double precision A(N,N), B(N,N), the R8LT factor matrices.
c
c    Output, double precision C(N,N), the R8LT product matrix.
c
      implicit none

      integer n

      double precision a(n,n)
      double precision b(n,n)
      double precision c(n,n)
      integer i
      integer j

      c(1:n,1:n) = 0.0D+00

      do i = 1, n
        do j = 1, i
          c(i,j) = sum ( a(i,j:i) * b(j:i,j) )
        end do
      end do

      return
      end
      subroutine r8lt_mxv ( m, n, a, x, b )

c*********************************************************************72
c
cc R8LT_MXV multiplies an R8LT matrix by an R8VEC.
c
c  Discussion:
c
c    The R8LT storage format is used for an M by N lower triangular matrix,
c    and sets aside storage even for the entries that must be zero.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 May 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input, double precision A(M,N), the R8LT matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(M), the product A * x.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision b(m)
      integer i
      integer jmax
      double precision x(n)

      do i = 1, m
        jmax = min ( i, n )
        b(i) = sum ( a(i,1:jmax) * x(1:jmax) )
      end do

      return
      end
      subroutine r8lt_print ( m, n, a, title )

c*********************************************************************72
c
cc R8LT_PRINT prints an R8LT matrix.
c
c  Discussion:
c
c    The R8LT storage format is used for an M by N lower triangular matrix,
c    and sets aside storage even for the entries that must be zero.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 May 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input, double precision A(M,N), the R8LT matrix.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      character ( len = * )  title

      call r8lt_print_some ( m, n, a, 1, 1, m, n, title )

      return
      end
      subroutine r8lt_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

c*********************************************************************72
c
cc R8LT_PRINT_SOME prints some of an R8LT matrix.
c
c  Discussion:
c
c    The R8LT storage format is used for an M by N lower triangular matrix,
c    and sets aside storage even for the entries that must be zero.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 May 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input, double precision A(M,N), the R8LT matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer, parameter :: incx = 5
      integer m
      integer n

      double precision a(m,n)
      character ( len = 14 ) ctemp(incx)
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
      character ( len = * )  title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )

      if ( ilo .lt. jlo ) then
        return
      end if
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
        i2lo = max ( i2lo, j2lo )
        i2hi = min ( ihi, m )

        do i = i2lo, i2hi
c
c  Print out (up to) 5 entries in row I, that lie in the current strip.
c
          do j2 = 1, inc

            j = j2lo - 1 + j2

            if ( i .lt. j ) then
              ctemp(j2) = '              '
            else if ( r8_is_int ( a(i,j) ) ) then
              write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
            else
              write ( ctemp(j2), '(g14.6)' ) a(i,j)
            end if

          end do

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine r8lt_random ( m, n, seed, a )

c*********************************************************************72
c
cc R8LT_RANDOM randomizes an R8LT matrix.
c
c  Discussion:
c
c    The R8LT storage format is used for an M by N lower triangular matrix,
c    and sets aside storage even for the entries that must be zero.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 October 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns of
c    the matrix.  M and N must be positive.
c
c    Input/output, integer SEED, a seed for the random number 
c    generator.
c
c    Output, double precision A(M,N), the R8LT matrix.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision r8_uniform_01
      integer i
      integer j
      integer seed

      do i = 1, m
        do j = 1, min ( i, n )
          a(i,j) = r8_uniform_01 ( seed )
        end do
        do j = i+1, n
          a(i,j) = 0.0D+00
        end do
      end do

      return
      end
      subroutine r8lt_sl ( n, a, b, job )

c*********************************************************************72
c
cc R8LT_SL solves an R8LT system.
c
c  Discussion:
c
c    The R8LT storage format is used for an M by N lower triangular matrix,
c    and sets aside storage even for the entries that must be zero.
c
c    No factorization of the lower triangular matrix is required.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 September 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision A(N,N), the R8LT matrix.
c
c    Input/output, double precision B(N).
c
c    On input, the right hand side.
c    On output, the solution vector.
c
c    Input, integer JOB, is 0 to solve the untransposed system,
c    nonzero to solve the transposed system.
c
      implicit none

      integer n

      double precision a(n,n)
      double precision b(n)
      integer j
      integer job

      if ( job .eq. 0 ) then

        do j = 1, n
          b(j) = b(j) / a(j,j)
          b(j+1:n) = b(j+1:n) - a(j+1:n,j) * b(j)
        end do

      else

        do j = n, 1, -1
          b(j) = b(j) / a(j,j)
          b(1:j-1) = b(1:j-1) - a(j,1:j-1) * b(j)
        end do

      end if

      return
      end
      subroutine r8lt_vxm ( m, n, a, x, b )

c*********************************************************************72
c
cc R8LT_VXM multiplies an R8VEC by an R8LT matrix.
c
c  Discussion:
c
c    The R8LT storage format is used for an M by N lower triangular matrix,
c    and sets aside storage even for the entries that must be zero.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 May 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input, double precision A(M,N), the R8LT matrix.
c
c    Input, double precision X(M), the vector to be multiplied by A.
c
c    Output, double precision B(N), the product A * x.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision b(n)
      integer i
      double precision x(m)

      do i = 1, n
        b(i) = sum ( x(i:m) * a(i:m,i) )
      end do

      return
      end
      subroutine r8mat_copy ( m, n, a, b )

c*********************************************************************72
c
cc R8MAT_COPY copies an R8MAT.
c
c  Discussion:
c
c    An R8MAT is an array of R8 values.
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
c    Input, integer M, N, the order of the matrix.
c
c    Input, double precision A(M,N), the matrix to be copied.
c
c    Output, double precision B(M,N), a copy of the matrix.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision b(m,n)
      integer i
      integer j

      do j = 1, n
        do i = 1, m
          b(i,j) = a(i,j)
        end do
      end do

      return
      end
      subroutine r8mat_diag_add_scalar ( n, a, s )

c*********************************************************************72
c
cc R8MAT_DIAG_ADD_SCALAR adds a scalar to the diagonal of a matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of rows and columns of the matrix.
c
c    Input/output, double precision A(N,N), the N by N matrix to be modified.
c
c    Input, double precision S, the value to be added to the diagonal 
c    of the matrix.
c
      implicit none

      integer n

      double precision a(n,n)
      integer i
      double precision s

      do i = 1, n
        a(i,i) = a(i,i) + s
      end do

      return
      end
      subroutine r8mat_diag_get_vector ( n, a, v )

c*********************************************************************72
c
cc R8MAT_DIAG_GET_VECTOR gets the value of the diagonal of a matrix.
c
c  Licensing:
c
c   This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of rows and columns of the matrix.
c
c    Input, double precision A(N,N), the N by N matrix.
c
c    Output, double precision V(N), the diagonal entries of the matrix.
c
      implicit none

      integer n

      double precision a(n,n)
      integer i
      double precision v(n)

      do i = 1, n
        v(i) = a(i,i)
      end do

      return
      end
      subroutine r8mat_diag_set_scalar ( n, a, s )

c*********************************************************************72
c
cc R8MAT_DIAG_SET_SCALAR sets the diagonal of a matrix to a scalar value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of rows and columns of the matrix.
c
c    Input/output, double precision A(N,N), the N by N matrix to be modified.
c
c    Input, double precision S, the value to be assigned to the 
c    diagonal of the matrix.
c
      implicit none

      integer n

      double precision a(n,n)
      integer i
      double precision s

      do i = 1, n
        a(i,i) = s
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
      subroutine r8ncf_indicator ( m, n, nz_num, rowcol, a )

c*********************************************************************72
c
cc R8NCF_INDICATOR sets up an R8NCF indicator matrix.
c
c  Discussion:
c
c    The "indicator matrix" simply has a value like I*10+J at every
c    entry of a dense matrix, or at every entry of a compressed storage
c    matrix for which storage is allocated. 
c
c    The R8NCF storage format stores NZ_NUM, the number of nonzeros, 
c    a real array containing the nonzero values, a 2 by NZ_NUM integer 
c    array storing the row and column of each nonzero entry.
c
c    The R8NCF format is used by NSPCG.  NSPCG requires that the information
c    for the diagonal entries of the matrix must come first.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 August 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns in
c    the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero entries.
c
c    Input, integer ROWCOL(2,NZ_NUM), the coordinates of 
c    the nonzero entries.
c
c    Output, double precision A(NZ_NUM), the indicator matrix.
c
      implicit none

      integer nz_num

      double precision a(nz_num)
      integer fac
      integer i
      integer i4_log_10
      integer isym
      integer j
      integer k
      integer m
      integer n
      integer rowcol(2,nz_num)

      fac = 10 ** ( i4_log_10 ( n ) + 1 )

      do k = 1, nz_num
        i = rowcol(1,k)
        j = rowcol(2,k)
        a(k) = dble ( fac * i + j )
      end do

      return
      end
      subroutine r8ncf_print ( m, n, nz_num, rowcol, a, title )

c*********************************************************************72
c
cc R8NCF_PRINT prints an R8NCF matrix.
c
c  Discussion:
c
c    The R8NCF storage format stores NZ_NUM, the number of nonzeros, 
c    a real array containing the nonzero values, a 2 by NZ_NUM integer 
c    array storing the row and column of each nonzero entry.
c
c    The R8NCF format is used by NSPCG.  NSPCG requires that the information
c    for the diagonal entries of the matrix must come first.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 August 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns of 
c    the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero elements in 
c    the matrix.
c
c    Input, integer ROWCOL(2,NZ_NUM), the row and column indices
c    of the nonzero elements.
c
c    Input, double precision A(NZ_NUM), the nonzero elements of the matrix.
c
c    Input, character ( len = * ) TITLE, a title.
c 
      implicit none

      integer nz_num

      double precision a(nz_num)
      integer m
      integer n
      integer rowcol(2,nz_num)
      character ( len = * ) title

      call r8ncf_print_some ( m, n, nz_num, rowcol, a, 1, 1, m, n, 
     &  title )

      return
      end
      subroutine r8ncf_print_some ( m, n, nz_num, rowcol, a, ilo, jlo, 
     &  ihi, jhi, title )

c*********************************************************************72
c
cc R8NCF_PRINT_SOME prints some of an R8NCF matrix.
c
c  Discussion:
c
c    The R8NCF storage format stores NZ_NUM, the number of nonzeros, 
c    a real array containing the nonzero values, a 2 by NZ_NUM integer 
c    array storing the row and column of each nonzero entry.
c
c    The R8NCF format is used by NSPCG.  NSPCG requires that the information
c    for the diagonal entries of the matrix must come first.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 August 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns of 
c    the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero elements 
c    in the matrix.
c
c    Input, integer ROWCOL(2,NZ_NUM), the row and column indices
c    of the nonzero elements.
c
c    Input, double precision A(NZ_NUM), the nonzero elements of the matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 5 )
      integer nz_num

      double precision a(nz_num)
      double precision aij
      character ( len = 14 ) ctemp(incx)
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
      integer k
      integer m
      integer n
      logical nonzero
      integer rowcol(2,nz_num)
      character ( len = * ) title

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

        write ( *, '(''  Col:  '',5a14)' ) ( ctemp(j2), j2 = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) '  ---'
c
c  Determine the range of the rows in this strip.
c
        i2lo = max ( ilo, 1 )
        i2hi = min ( ihi, m )

        do i = i2lo, i2hi
c
c  Print out (up to) 5 entries in row I, that lie in the current strip.
c
          nonzero = .false.

          aij = 0.0D+00
          do j2 = 1, inc
            write ( ctemp(j2), '(f8.0,6x)' ) aij
          end do

          do k = 1, nz_num

            if ( 
     &        i .eq. rowcol(1,k) .and. 
     &        j2lo .le. rowcol(2,k) .and. 
     &        rowcol(2,k) .le. j2hi ) then 

              j2 = rowcol(2,k) - j2lo + 1
              aij = a(k)

              if ( aij .eq. 0.0D+00 ) then
                cycle
              end if

              nonzero = .true.

              if ( r8_is_int ( aij ) ) then
                write ( ctemp(j2), '(f8.0,6x)' ) aij
              else
                write ( ctemp(j2), '(g14.6)' ) aij
              end if
            end if

          end do

          if ( nonzero ) then
            write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )
          end if

        end do

      end do

      return
      end
      subroutine r8pbl_det ( n, mu, a_lu, det )

c*********************************************************************72
c
cc R8PBL_DET computes the determinant of a matrix factored by R8PBL_FA.
c
c  Discussion:
c
c    The R8PBL storage format is for a symmetric positive definite band matrix.
c
c    To save storage, only the diagonal and lower triangle of A is stored,
c    in a compact diagonal format that preserves columns.
c
c    The diagonal is stored in row 1 of the array.
c    The first subdiagonal in row 2, columns 1 through MU.
c    The second subdiagonal in row 3, columns 1 through MU-1.
c    The MU-th subdiagonal in row MU+1, columns 1 through 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 October 2003
c
c  Author:
c
c    Original FORTRAN77 version Dongarra, Bunch, Moler, Stewart.
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
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input, integer MU, the upper (and lower) bandwidth.
c    MU must be nonnegative, and no greater than N-1.
c
c    Input, double precision A_LU(MU+1,N), the LU factors from R8PBL_FA.
c
c    Output, double precision DET, the determinant of the matrix.
c
      implicit none

      integer mu
      integer n

      double precision a_lu(mu+1,n)
      double precision det

      det = product ( a_lu(1,1:n)**2 )

      return
      end
      subroutine r8pbl_indicator ( n, mu, a )

c*********************************************************************72
c
cc R8PBL_INDICATOR sets up an R8PBL indicator matrix.
c
c  Discussion:
c
c    The "indicator matrix" simply has a value like I*10+J at every
c    entry of a dense matrix, or at every entry of a compressed storage
c    matrix for which storage is allocated. 
c
c    The R8PBL storage format is for a symmetric positive definite band matrix.
c
c    To save storage, only the diagonal and lower triangle of A is stored,
c    in a compact diagonal format that preserves columns.
c
c    The diagonal is stored in row 1 of the array.
c    The first subdiagonal in row 2, columns 1 through MU.
c    The second subdiagonal in row 3, columns 1 through MU-1.
c    The MU-th subdiagonal in row MU+1, columns 1 through 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 January 2004
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
c    Input, integer MU, the number of subdiagonals in the matrix.
c    MU must be at least 0 and no more than N-1.
c
c    Output, double precision A(MU+1,N), the R8PBL matrix.
c
      implicit none

      integer mu
      integer n

      double precision a(mu+1,n)
      integer fac
      integer i
      integer i4_log_10
      integer j

      fac = 10 ** ( i4_log_10 ( n ) + 1 )
c
c  Zero out the "junk" entries.
c
      do j = n+1-mu, n
        do i = n+1, j+mu
          a(i-j+1,j) = 0.0D+00
        end do
      end do
c
c  Set the meaningful values.
c
      do i = 1, n
        do j = max ( 1, i - mu ), i
          a(i-j+1,j) = dble ( fac * i + j )
        end do
      end do

      return
      end
      subroutine r8pbl_print ( n, mu, a, title )

c*********************************************************************72
c
cc R8PBL_PRINT prints an R8PBL matrix.
c
c  Discussion:
c
c    The R8PBL storage format is for a symmetric positive definite band matrix.
c
c    To save storage, only the diagonal and lower triangle of A is stored,
c    in a compact diagonal format that preserves columns.
c
c    The diagonal is stored in row 1 of the array.
c    The first subdiagonal in row 2, columns 1 through MU.
c    The second subdiagonal in row 3, columns 1 through MU-1.
c    The MU-th subdiagonal in row MU+1, columns 1 through 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 May 2003
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
c    Input, integer MU, the upper (and lower) bandwidth.
c    MU must be nonnegative, and no greater than N-1.
c
c    Input, double precision A(MU+1,N), the R8PBL matrix.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer mu
      integer n

      double precision a(mu+1,n)
      character ( len = * ) title

      call r8pbl_print_some ( n, mu, a, 1, 1, n, n, title )

      return
      end
      subroutine r8pbl_print_some ( n, mu, a, ilo, jlo, ihi, jhi, 
     &  title )

c*********************************************************************72
c
cc R8PBL_PRINT_SOME prints some of an R8PBL matrix.
c
c  Discussion:
c
c    The R8PBL storage format is for a symmetric positive definite band matrix.
c
c    To save storage, only the diagonal and lower triangle of A is stored,
c    in a compact diagonal format that preserves columns.
c
c    The diagonal is stored in row 1 of the array.
c    The first subdiagonal in row 2, columns 1 through MU.
c    The second subdiagonal in row 3, columns 1 through MU-1.
c    The MU-th subdiagonal in row MU+1, columns 1 through 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 May 2003
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
c    Input, integer MU, the upper (and lower) bandwidth.
c    MU must be nonnegative, and no greater than N-1.
c
c    Input, double precision A(MU+1,N), the R8PBL matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 5 )
      integer mu
      integer n

      double precision a(mu+1,n)
      double precision aij
      character ( len = 14 ) ctemp(incx)
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
      character ( len = * ) title

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
        i2lo = max ( i2lo, j2lo - mu )
        i2hi = min ( ihi, n )
        i2hi = min ( i2hi, j2hi + mu )

        do i = i2lo, i2hi
c
c  Print out (up to) 5 entries in row I, that lie in the current strip.
c
          do j2 = 1, inc

            j = j2lo - 1 + j2

            if ( i .le. j .and. j .le. i + mu ) then
              aij = a(j-i+1,i)
            else if ( j .le. i .and. i .le. j + mu ) then
              aij = a(i-j+1,j)
            else
              aij = 0.0D+00
            end if

            if ( mu .lt. i-j .or. mu .lt. j-i ) then
              ctemp(j2) = '              '
            else if ( r8_is_int ( aij ) ) then
              write ( ctemp(j2), '(f8.0,6x)' ) aij
            else
              write ( ctemp(j2), '(g14.6)' ) aij
            end if

          end do

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine r8pbl_random ( n, mu, seed, a )

c*********************************************************************72
c
cc R8PBL_RANDOM randomizes an R8PBL matrix.
c
c  Discussion:
c
c    The R8PBL storage format is for a symmetric positive definite band matrix.
c
c    To save storage, only the diagonal and lower triangle of A is stored,
c    in a compact diagonal format that preserves columns.
c
c    The diagonal is stored in row 1 of the array.
c    The first subdiagonal in row 2, columns 1 through MU.
c    The second subdiagonal in row 3, columns 1 through MU-1.
c    The MU-th subdiagonal in row MU+1, columns 1 through 1.
c
c    The matrix returned will be positive definite, but of limited
c    randomness.  The off diagonal elements are random values between
c    0 and 1, and the diagonal element of each row is selected to
c    ensure strict diagonal dominance.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 October 2004
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
c    Input, integer MU, the number of subdiagonals in the matrix.
c    MU must be at least 0 and no more than N-1.
c
c    Input/output, integer SEED, a seed for the random number 
c    generator.
c
c    Output, double precision A(MU+1,N), the R8PBL matrix.
c
      implicit none

      integer mu
      integer n

      double precision a(mu+1,n)
      double precision r8_uniform_01
      integer i
      integer j
      integer jhi
      integer jlo
      double precision r
      integer seed
      double precision sum2
c
c  Zero out the "junk" entries.
c
      do j = n+1-mu, n
        a(2+n-j:mu+1,j) = 0.0D+00
      end do
c
c  Set the off diagonal values.
c
      do i = 1, n
        do j = max ( 1, i - mu ), i - 1
          a(i-j+1,j) = r8_uniform_01 ( seed )
        end do
      end do
c
c  Set the diagonal values.
c
      do i = 1, n

        sum2 = 0.0D+00

        jlo = max ( 1, i - mu )
        do j = jlo, i-1
          sum2 = sum2 + abs ( a(i-j+1,j) )
        end do

        jhi = min ( i + mu, n )
        do j = i+1, jhi
          sum2 = sum2 + abs ( a(j-i+1,i) )
        end do

        r = r8_uniform_01 ( seed )

        a(1,i) = ( 1.0D+00 + r ) * ( sum2 + 0.01D+00 )

      end do

      return
      end
      subroutine r8pbl_to_r8ge ( n, mu, a, b )

c*********************************************************************72
c
cc R8PBL_TO_R8GE copies an R8PBL matrix to an R8GE matrix.
c
c  Discussion:
c
c    The R8PBL storage format is for a symmetric positive definite band matrix.
c
c    To save storage, only the diagonal and lower triangle of A is stored,
c    in a compact diagonal format that preserves columns.
c
c    The diagonal is stored in row 1 of the array.
c    The first subdiagonal in row 2, columns 1 through MU.
c    The second subdiagonal in row 3, columns 1 through MU-1.
c    The MU-th subdiagonal in row MU+1, columns 1 through 1.
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    14 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrices.
c    N must be positive.
c
c    Input, integer MU, the upper bandwidth of A.
c    MU must be nonnegative, and no greater than N-1.
c
c    Input, double precision A(MU+1,N), the R8PBL matrix.
c
c    Output, double precision B(N,N), the R8GE matrix.
c
      implicit none

      integer mu
      integer n

      double precision a(mu+1,n)
      double precision b(n,n)
      integer i
      integer j

      do i = 1, n
        do j = 1, n
          if ( i .le. j .and. j .le. i + mu ) then
            b(i,j) = a(j-i+1,i)
          else if ( i - mu .le. j .and. j .lt. i ) then
            b(i,j) = a(i-j+1,j)
          else
            b(i,j) = 0.0D+00
          end if
        end do
      end do

      return
      end
      subroutine r8pbu_cg ( n, mu, a, b, x )

c*********************************************************************72
c
cc R8PBU_CG uses the conjugate gradient method on an R8PBU system.
c
c  Discussion:
c
c    The R8PBU storage format is for a symmetric positive definite band matrix.
c
c    To save storage, only the diagonal and upper triangle of A is stored,
c    in a compact diagonal format that preserves columns.
c
c    The diagonal is stored in row MU+1 of the array.
c    The first superdiagonal in row MU, columns 2 through N.
c    The second superdiagonal in row MU-1, columns 3 through N.
c    The MU-th superdiagonal in row 1, columns MU+1 through N.
c
c    The matrix A must be a positive definite symmetric band matrix.
c
c    The method is designed to reach the solution after N computational
c    steps.  However, roundoff may introduce unacceptably large errors for
c    some problems.  In such a case, calling the routine again, using
c    the computed solution as the new starting estimate, should improve
c    the results.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 October 1998
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    FS Beckman,
c    The Solution of Linear Equations by the Conjugate Gradient Method,
c    Mathematical Methods for Digital Computers, pages 62-72.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input, integer MU, the number of superdiagonals.
c    MU must be at least 0, and no more than N-1.
c
c    Input, double precision A(MU+1,N), the R8PBU matrix.
c
c    Input, double precision B(N), the right hand side vector.
c
c    Input/output, double precision X(N).
c    On input, an estimate for the solution, which may be 0.
c    On output, the approximate solution vector.
c
      implicit none

      integer mu
      integer n

      double precision a(mu+1,n)
      double precision alpha
      double precision ap(n)
      double precision b(n)
      double precision beta
      integer it
      double precision p(n)
      double precision pap
      double precision pr
      double precision r(n)
      double precision rap
      double precision x(n)
c
c  Initialize
c    AP = A * x,
c    R  = b - A * x,
c    P  = b - A * x.
c
      call r8pbu_mxv ( n, mu, a, x, ap )

      r(1:n) = b(1:n) - ap(1:n)
      p(1:n) = b(1:n) - ap(1:n)
c
c  Do the N steps of the conjugate gradient method.
c
      do it = 1, n
c
c  Compute the matrix*vector product AP=A*P.
c
        call r8pbu_mxv ( n, mu, a, p, ap )
c
c  Compute the dot products
c    PAP = P*AP,
c    PR  = P*R
c  Set
c    ALPHA = PR / PAP.
c
        pap = sum ( p(1:n) * ap(1:n) )
        pr = sum ( p(1:n) * r(1:n) )

        if ( pap .eq. 0.0D+00 ) then
          return
        end if

        alpha = pr / pap
c
c  Set
c    X = X + ALPHA * P
c    R = R - ALPHA * AP.
c
        x(1:n) = x(1:n) + alpha * p(1:n)
        r(1:n) = r(1:n) - alpha * ap(1:n)
c
c  Compute the vector dot product
c    RAP = R*AP
c  Set
c    BETA = - RAP / PAP.
c
        rap = sum ( r(1:n) * ap(1:n) )

        beta = - rap / pap
c
c  Update the perturbation vector
c    P = R + BETA * P.
c
        p(1:n) = r(1:n) + beta * p(1:n)

      end do

      return
      end
      subroutine r8pbu_det ( n, mu, a_lu, det )

c*********************************************************************72
c
cc R8PBU_DET computes the determinant of a matrix factored by R8PBU_FA.
c
c  Discussion:
c
c    The R8PBU storage format is for a symmetric positive definite band matrix.
c
c    To save storage, only the diagonal and upper triangle of A is stored,
c    in a compact diagonal format that preserves columns.
c
c    The diagonal is stored in row MU+1 of the array.
c    The first superdiagonal in row MU, columns 2 through N.
c    The second superdiagonal in row MU-1, columns 3 through N.
c    The MU-th superdiagonal in row 1, columns MU+1 through N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 October 1998
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
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input, integer MU, the number of superdiagonals of the matrix.
c    MU must be at least 0 and no more than N-1.
c
c    Input, double precision A_LU(MU+1,N), the LU factors from R8PBU_FA.
c
c    Output, double precision DET, the determinant of the matrix.
c
      implicit none

      integer mu
      integer n

      double precision a_lu(mu+1,n)
      double precision det

      det = product ( a_lu(mu+1,1:n)**2 )

      return
      end
      subroutine r8pbu_fa ( n, mu, a, info )

c*********************************************************************72
c
cc R8PBU_FA factors an R8PBU matrix.
c
c  Discussion:
c
c    The R8PBU storage format is for a symmetric positive definite band matrix.
c
c    To save storage, only the diagonal and upper triangle of A is stored,
c    in a compact diagonal format that preserves columns.
c
c    The diagonal is stored in row MU+1 of the array.
c    The first superdiagonal in row MU, columns 2 through N.
c    The second superdiagonal in row MU-1, columns 3 through N.
c    The MU-th superdiagonal in row 1, columns MU+1 through N.
c
c    The matrix A must be a positive definite symmetric band matrix.
c
c    Once factored, linear systems A*x=b involving the matrix can be solved
c    by calling R8PBU_SL.  No pivoting is performed.  Pivoting is not necessary
c    for positive definite symmetric matrices.  If the matrix is not positive
c    definite, the algorithm may behave correctly, but it is also possible
c    that an illegal divide by zero will occur.
c
c  Modified:
c
c    31 October 1998
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Author:
c
c    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
c    This FORTRAN77 version by John Burkardt
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
c    Input, integer MU, the number of superdiagonals of the matrix.
c    MU must be at least 0, and no more than N-1.
c
c    Input/output, double precision A(MU+1,N), the N by N matrix, stored 
c    in LINPACK positive definite symmetric band matrix storage.
c    On output, A contains information describing a factored form
c    of the matrix, that can be used to solve linear systems
c    A*x=b, using R8PBU_SL.
c
c    Output, integer INFO, singularity flag.
c    0, the matrix is nonsingular.
c    nonzero, the matrix is singular.
c
      implicit none

      integer mu
      integer n

      double precision a(mu+1,n)
      integer ik
      integer info
      integer j
      integer jk
      integer k
      integer mm
      double precision s

      info = 0

      do j = 1, n

        ik = mu + 1
        jk = max ( j - mu, 1 )
        mm = max ( mu + 2 - j, 1 )

        s = 0.0D+00

        do k = mm, mu

          a(k,j) = ( a(k,j) - sum ( a(ik:ik+k-mm-1,jk) * a(mm:k-1,j) ) ) 
     &      / a(mu+1,jk)

          s = s + a(k,j) * a(k,j)

          ik = ik - 1
          jk = jk + 1

        end do

        s = a(mu+1,j) - s

        if ( s .le. 0.0D+00 ) then
          info = j
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8PBU_FA - Fatal error!'
          write ( *, '(a,i8)' ) '  Nonpositive pivot on step ', info
          stop
        end if

        a(mu+1,j) = sqrt ( s )

      end do

      return
      end
      subroutine r8pbu_indicator ( n, mu, a )

c*********************************************************************72
c
cc R8PBU_INDICATOR sets up an R8PBU indicator matrix.
c
c  Discussion:
c
c    The "indicator matrix" simply has a value like I*10+J at every
c    entry of a dense matrix, or at every entry of a compressed storage
c    matrix for which storage is allocated. 
c
c    The R8PBU storage format is for a symmetric positive definite band matrix.
c
c    To save storage, only the diagonal and upper triangle of A is stored,
c    in a compact diagonal format that preserves columns.
c
c    The diagonal is stored in row MU+1 of the array.
c    The first superdiagonal in row MU, columns 2 through N.
c    The second superdiagonal in row MU-1, columns 3 through N.
c    The MU-th superdiagonal in row 1, columns MU+1 through N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    14 January 2004
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
c    Input, integer MU, the number of superdiagonals in the matrix.
c    MU must be at least 0 and no more than N-1.
c
c    Output, double precision A(MU+1,N), the R8PBU matrix.
c
      implicit none

      integer mu
      integer n

      double precision a(mu+1,n)
      integer fac
      integer i
      integer i4_log_10
      integer j

      fac = 10 ** ( i4_log_10 ( n ) + 1 )
c
c  Zero out the "junk" entries.
c
      do j = 1, mu
        do i = 1, mu+1-j
          a(i,j) = 0.0D+00
        end do
      end do
c
c  Set the meaningful values.
c
      do i = 1, n
        do j = i, min ( i+mu, n )
          a(mu+1+i-j,j) = dble ( fac * i + j )
        end do
      end do

      return
      end
      subroutine r8pbu_ml ( n, mu, a_lu, x, b )

c*********************************************************************72
c
cc R8PBU_ML multiplies an R8VEC times a matrix that was factored by R8PBU_FA.
c
c  Discussion:
c
c    The R8PBU storage format is for a symmetric positive definite band matrix.
c
c    To save storage, only the diagonal and upper triangle of A is stored,
c    in a compact diagonal format that preserves columns.
c
c    The diagonal is stored in row MU+1 of the array.
c    The first superdiagonal in row MU, columns 2 through N.
c    The second superdiagonal in row MU-1, columns 3 through N.
c    The MU-th superdiagonal in row 1, columns MU+1 through N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    31 October 1998
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
c    Input, integer MU, the number of superdiagonals of the matrix.
c    MU must be at least 0 and no more than N-1.
c
c    Input, double precision A_LU(MU+1,N), the LU factors from R8PBU_FA.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(N), the product A * x.
c
      implicit none

      integer mu
      integer n

      double precision a_lu(mu+1,n)
      double precision b(n)
      integer i
      integer ilo
      integer j
      integer jhi
      integer k
      double precision x(n)

      b(1:n) = x(1:n)
c
c  Multiply U * X = Y.
c
      do k = 1, n

        ilo = max ( 1, k - mu )
        do i = ilo, k - 1
          b(i) = b(i) + a_lu(mu+1+i-k,k) * b(k)
        end do

        b(k) = a_lu(mu+1,k) * b(k)

      end do
c
c  Multiply L * Y = B.
c
      do k = n, 1, -1

        jhi = min ( k + mu, n )
        do j = k + 1, jhi
          b(j) = b(j) + a_lu(mu+1+k-j,j) * b(k)
        end do

        b(k) = a_lu(mu+1,k) * b(k)

      end do

      return
      end
      subroutine r8pbu_mxv ( n, mu, a, x, b )

c*********************************************************************72
c
cc R8PBU_MXV multiplies an R8PBU matrix by an R8VEC.
c
c  Discussion:
c
c    The R8PBU storage format is for a symmetric positive definite band matrix.
c
c    To save storage, only the diagonal and upper triangle of A is stored,
c    in a compact diagonal format that preserves columns.
c
c    The diagonal is stored in row MU+1 of the array.
c    The first superdiagonal in row MU, columns 2 through N.
c    The second superdiagonal in row MU-1, columns 3 through N.
c    The MU-th superdiagonal in row 1, columns MU+1 through N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 October 1998
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
c    Input, integer MU, the number of superdiagonals in the matrix.
c    MU must be at least 0 and no more than N-1.
c
c    Input, double precision A(MU+1,N), the R8PBU matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(N), the result vector A * x.
c
      implicit none

      integer mu
      integer n

      double precision a(mu+1,n)
      double precision b(n)
      integer i
      integer ieqn
      integer j
      double precision x(n)
c
c  Multiply X by the diagonal of the matrix.
c
      b(1:n) = a(mu+1,1:n) * x(1:n)
c
c  Multiply X by the superdiagonals of the matrix.
c
      do i = mu, 1, -1
        do j = mu+2-i, n
          ieqn = i + j - mu - 1
          b(ieqn) = b(ieqn) + a(i,j) * x(j)
          b(j) = b(j) + a(i,j) * x(ieqn)
        end do
      end do

      return
      end
      subroutine r8pbu_print ( n, mu, a, title )

c*********************************************************************72
c
cc R8PBU_PRINT prints an R8PBU matrix.
c
c  Discussion:
c
c    The R8PBU storage format is for a symmetric positive definite band matrix.
c
c    To save storage, only the diagonal and upper triangle of A is stored,
c    in a compact diagonal format that preserves columns.
c
c    The diagonal is stored in row MU+1 of the array.
c    The first superdiagonal in row MU, columns 2 through N.
c    The second superdiagonal in row MU-1, columns 3 through N.
c    The MU-th superdiagonal in row 1, columns MU+1 through N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 May 2000
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
c    Input, integer MU, the upper (and lower) bandwidth.
c    MU must be nonnegative, and no greater than N-1.
c
c    Input, double precision A(MU+1,N), the R8PBU matrix.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer mu
      integer n

      double precision a(mu+1,n)
      character ( len = * ) title

      call r8pbu_print_some ( n, mu, a, 1, 1, n, n, title )

      return
      end
      subroutine r8pbu_print_some ( n, mu, a, ilo, jlo, ihi, jhi, 
     &  title )

c*********************************************************************72
c
cc R8PBU_PRINT_SOME prints some of an R8PBU matrix.
c
c  Discussion:
c
c    The R8PBU storage format is for a symmetric positive definite band matrix.
c
c    To save storage, only the diagonal and upper triangle of A is stored,
c    in a compact diagonal format that preserves columns.
c
c    The diagonal is stored in row MU+1 of the array.
c    The first superdiagonal in row MU, columns 2 through N.
c    The second superdiagonal in row MU-1, columns 3 through N.
c    The MU-th superdiagonal in row 1, columns MU+1 through N.
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
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input, integer MU, the upper (and lower) bandwidth.
c    MU must be nonnegative, and no greater than N-1.
c
c    Input, double precision A(MU+1,N), the R8PBU matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 5 )
      integer mu
      integer n

      double precision a(mu+1,n)
      double precision aij
      character ( len = 14 ) ctemp(incx)
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
      character ( len = * ) title

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
        i2lo = max ( i2lo, j2lo - mu )
        i2hi = min ( ihi, n )
        i2hi = min ( i2hi, j2hi + mu )

        do i = i2lo, i2hi
c
c  Print out (up to) 5 entries in row I, that lie in the current strip.
c
          do j2 = 1, inc

            j = j2lo - 1 + j2

            if ( i .le. j .and. j .le. i + mu ) then
              aij = a(mu+1+i-j,j)
            else if ( i - mu .le. j .and. j .le. i ) then
              aij = a(mu+1+j-i,i)
            else
              aij = 0.0D+00
            end if

            if ( mu .lt. i-j .or. mu .lt. j-i ) then
              ctemp(j2) = '              '
            else if ( r8_is_int ( aij ) ) then
              write ( ctemp(j2), '(f8.0,6x)' ) aij
            else
              write ( ctemp(j2), '(g14.6)' ) aij
            end if

          end do

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine r8pbu_random ( n, mu, seed, a )

c*********************************************************************72
c
cc R8PBU_RANDOM randomizes an R8PBU matrix.
c
c  Discussion:
c
c    The R8PBU storage format is for a symmetric positive definite band matrix.
c
c    To save storage, only the diagonal and upper triangle of A is stored,
c    in a compact diagonal format that preserves columns.
c
c    The diagonal is stored in row MU+1 of the array.
c    The first superdiagonal in row MU, columns 2 through N.
c    The second superdiagonal in row MU-1, columns 3 through N.
c    The MU-th superdiagonal in row 1, columns MU+1 through N.
c
c    The matrix returned will be positive definite, but of limited
c    randomness.  The off diagonal elements are random values between
c    0 and 1, and the diagonal element of each row is selected to
c    ensure strict diagonal dominance.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 May 2003
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
c    Input, integer MU, the number of superdiagonals in the matrix.
c    MU must be at least 0 and no more than N-1.
c
c    Input/output, integer SEED, a seed for the random number 
c    generator.
c
c    Output, double precision A(MU+1,N), the R8PBU matrix.
c
      implicit none

      integer mu
      integer n

      double precision a(mu+1,n)
      double precision r8_uniform_01
      integer i
      integer j
      integer jhi
      integer jlo
      double precision r
      integer seed
      double precision sum2
c
c  Zero out the "junk" entries.
c
      do j = 1, mu
        a(1:mu+1-j,j) = 0.0D+00
      end do
c
c  Set the off diagonal values.
c
      do i = 1, n
        do j = i+1, min ( i+mu, n )
          a(mu+1+i-j,j) = r8_uniform_01 ( seed )
        end do
      end do
c
c  Set the diagonal values.
c
      do i = 1, n

        sum2 = 0.0D+00

        jlo = max ( 1, i - mu )
        do j = jlo, i-1
          sum2 = sum2 + abs ( a(mu+1+j-i,i) )
        end do

        jhi = min ( i + mu, n )
        do j = i+1, jhi
          sum2 = sum2 + abs ( a(mu+1+i-j,j) )
        end do

        r = r8_uniform_01 ( seed )

        a(mu+1,i) = ( 1.0D+00 + r ) * ( sum2 + 0.01D+00 )

      end do

      return
      end
      subroutine r8pbu_sl ( n, mu, a_lu, b )

c*********************************************************************72
c
cc R8PBU_SL solves an R8PBU system factored by R8PBU_FA.
c
c  Discussion:
c
c    The R8PBU storage format is for a symmetric positive definite band matrix.
c
c    To save storage, only the diagonal and upper triangle of A is stored,
c    in a compact diagonal format that preserves columns.
c
c    The diagonal is stored in row MU+1 of the array.
c    The first superdiagonal in row MU, columns 2 through N.
c    The second superdiagonal in row MU-1, columns 3 through N.
c    The MU-th superdiagonal in row 1, columns MU+1 through N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    31 October 1998
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
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input, integer MU, the number of superdiagonals of the matrix.
c    MU must be at least 0 and no more than N-1.
c
c    Input, double precision A_LU(MU+1,N), the LU factors from R8PBU_FA.
c
c    Input/output, double precision B(N).
c    On input, B contains the right hand side of the linear system
c    to be solved.
c    On output, B contains X, the solution vector.
c
      implicit none

      integer mu
      integer n

      double precision a_lu(mu+1,n)
      double precision b(n)
      integer i
      integer ilo
      integer k
c
c  Solve L * Y = B.
c
      do k = 1, n
        ilo = max ( 1, k - mu )
        b(k) = ( b(k) - sum ( b(ilo:k-1) * a_lu(mu+1+ilo-k:mu,k) ) ) 
     &    / a_lu(mu+1,k)
      end do
c
c  Solve U * X = Y.
c
      do k = n, 1, -1

        b(k) = b(k) / a_lu(mu+1,k)

        ilo = max ( 1, k - mu )
        do i = ilo, k - 1
          b(i) = b(i) - b(k) * a_lu(mu+1+i-k,k)
        end do

      end do

      return
      end
      subroutine r8pbu_sor ( n, mu, a, b, eps, itchk, itknt, itmax, 
     &  omega, x )

c*********************************************************************72
c
cc R8PBU_SOR uses SOR iteration to solve an R8PBU linear system.
c
c  Discussion:
c
c    The R8PBU storage format is for a symmetric positive definite band matrix.
c
c    To save storage, only the diagonal and upper triangle of A is stored,
c    in a compact diagonal format that preserves columns.
c
c    The diagonal is stored in row MU+1 of the array.
c    The first superdiagonal in row MU, columns 2 through N.
c    The second superdiagonal in row MU-1, columns 3 through N.
c    The MU-th superdiagonal in row 1, columns MU+1 through N.
c
c    The matrix A must be a positive definite symmetric band matrix.
c
c    A relaxation factor OMEGA may be used.
c
c    The iteration will proceed until a convergence test is met,
c    or the iteration limit is reached.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 October 1998
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
c    Input, integer MU, the number of superdiagonals in the 
c    matrix.  MU must be at least 0, and no more than N-1.
c
c    Input, double precision A(MU+1,N), the R8PBU matrix.
c
c    Input, double precision B(N), the right hand side of the system.
c
c    Input, double precision EPS, convergence tolerance for the system. 
c    The vector b - A * x is computed every ITCHK iterations, and if the 
c    maximum entry of this vector is of norm less than EPS, the program
c    will return.
c
c    Input, integer ITCHK, the interval between convergence checks.
c    ITCHK steps will be taken before any check is made on whether the iteration
c    has converged.  ITCHK should be at least 1 and no greater
c    than ITMAX.
c
c    Output, integer ITKNT, the number of iterations taken.
c
c    Input, integer ITMAX, the maximum number of iterations 
c    allowed.  The program will return to the user if this many iterations 
c    are taken without convergence.
c
c    Input, double precision OMEGA, the relaxation factor.  OMEGA must be
c    strictly between 0 and 2.  Use OMEGA = 1 for no relaxation, classical
c    Jacobi iteration.
c
c    Input/output, double precision X(N).
c    On input, a starting vector for the iteration.
c    On output, the current approximation to the solution.
c
      implicit none

      integer mu
      integer n

      double precision a(mu+1,n)
      double precision b(n)
      double precision eps
      double precision err
      integer it
      integer itchk
      integer itknt
      integer itmax
      double precision omega
      double precision x(n)
      double precision xtemp(n)

      if ( itchk .le. 0 .or. itmax .lt. itchk ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8PBU_SOR - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal ITCHK= ', itchk
        stop
      end if

      if ( itmax .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8PBU_SOR - Fatal error!'
        write ( *, '(a,i8)' ) '  Nonpositive ITMAX =', itmax
        stop
      end if

      if ( omega .le. 0.0D+00 .or. 2.0D+00 .le. omega ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8PBU_SOR - Fatal error!'
        write ( *, '(a,g14.6)' ) '  Illegal value of OMEGA = ', omega
        stop
      end if

      itknt = 0
c
c  Take ITCHK steps of the iteration before doing a convergence check.
c
      do while ( itknt .le. itmax )

        do it = 1, itchk
c
c  Compute XTEMP(I) = B(I) + A(I,I) * X(I) - SUM ( J=1 to N ) A(I,J) * X(J).
c
          call r8pbu_mxv ( n, mu, a, x, xtemp )

          xtemp(1:n) = x(1:n) + ( b(1:n) - xtemp(1:n) ) / a(mu+1,1:n)
c
c  Compute the next iterate as a weighted combination of the
c  old iterate and the just computed standard Jacobi iterate.
c
          if ( omega .ne. 1.0D+00 ) then
            xtemp(1:n) = ( 1.0D+00 - omega ) * x(1:n) 
     &        + omega * xtemp(1:n)
          end if

          itknt = itknt + 1
c
c  Copy the new result into the old result vector.
c
          x(1:n) = xtemp(1:n)

        end do
c
c  Compute the maximum residual, the greatest entry in the vector
c  RESID(I) = B(I) - A(I,J) * X(J).
c
        call r8pbu_mxv ( n, mu, a, x, xtemp )

        err = maxval ( abs ( b(1:n) - xtemp(1:n) ) )
c
c  Test to see if we can quit because of convergence,
c
        if ( err .le. eps ) then
          return
        end if

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8PBU_SOR - Warning!'
      write ( *, '(a)' ) '  The iteration did not converge.'

      return
      end
      subroutine r8pbu_to_r8ge ( n, mu, a, b )

c*********************************************************************72
c
cc R8PBU_TO_R8GE copies an R8PBU matrix to an R8GE matrix.
c
c  Discussion:
c
c    The R8PBU storage format is for a symmetric positive definite band matrix.
c
c    To save storage, only the diagonal and upper triangle of A is stored,
c    in a compact diagonal format that preserves columns.
c
c    The diagonal is stored in row MU+1 of the array.
c    The first superdiagonal in row MU, columns 2 through N.
c    The second superdiagonal in row MU-1, columns 3 through N.
c    The MU-th superdiagonal in row 1, columns MU+1 through N.
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    14 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrices.
c    N must be positive.
c
c    Input, integer MU, the upper bandwidth of A1.
c    MU must be nonnegative, and no greater than N-1.
c
c    Input, double precision A(MU+1,N), the R8PBU matrix.
c
c    Output, double precision B(N,N), the R8GE matrix.
c
      implicit none

      integer mu
      integer n

      double precision a(mu+1,n)
      double precision b(n,n)
      integer i
      integer j

      do i = 1, n
        do j = 1, n
          if ( i .le. j .and. j .le. i+mu ) then
            b(i,j) = a(mu+1+i-j,j)
          else if ( i-mu .le. j .and. j .lt. i ) then
            b(i,j) = a(mu+1+j-i,i)
          else
            b(i,j) = 0.0D+00
          end if
        end do
      end do

      return
      end
      subroutine r8po_det ( n, a_lu, det )

c*********************************************************************72
c
cc R8PO_DET computes the determinant of a matrix factored by R8PO_FA.
c
c  Discussion:
c
c    The R8PO storage format is used for a symmetric positive definite 
c    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
c    upper triangular matrix, so it will be in R8GE storage format.)
c
c    Only the diagonal and upper triangle of the square array are used.
c    This same storage scheme is used when the matrix is factored by
c    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
c    is set to zero.
c
c    R8PO storage is used by LINPACK and LAPACK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 October 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision A_LU(N,N), the LU factors from R8PO_FA.
c
c    Output, double precision DET, the determinant of A.
c
      implicit none

      integer n

      double precision a_lu(n,n)
      double precision det
      integer i

      det = 1.0D+00

      do i = 1, n
        det = det * a_lu(i,i)
      end do

      return
      end
      subroutine r8po_fa ( n, a, info )

c*********************************************************************72
c
cc R8PO_FA factors an R8PO matrix.
c
c  Discussion:
c
c    The R8PO storage format is used for a symmetric positive definite 
c    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
c    upper triangular matrix, so it will be in R8GE storage format.)
c
c    Only the diagonal and upper triangle of the square array are used.
c    This same storage scheme is used when the matrix is factored by
c    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
c    is set to zero.
c
c    R8PO storage is used by LINPACK and LAPACK.
c
c    The positive definite symmetric matrix A has a Cholesky factorization
c    of the form:
c
c      A = R' * R
c
c    where R is an upper triangular matrix with positive elements on
c    its diagonal.  This routine overwrites the matrix A with its
c    factor R.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2003
c
c  Author:
c
c    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
c    FORTRAN90 version by John Burkardt.
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
c
c    Input/output, double precision A(N,N).
c    On input, the matrix in R8PO storage.
c    On output, the Cholesky factor R in R8GE storage.
c
c    Output, integer INFO, error flag.
c    0, normal return.
c    K, error condition.  The principal minor of order K is not
c    positive definite, and the factorization was not completed.
c
      implicit none

      integer n

      double precision a(n,n)
      integer i
      integer info
      integer j
      integer k
      double precision s

      do j = 1, n

        do k = 1, j - 1
          a(k,j) = ( a(k,j) - sum ( a(1:k-1,k) * a(1:k-1,j) ) ) / a(k,k)
        end do

        s = a(j,j) - sum ( a(1:j-1,j)**2 )

        if ( s .le. 0.0D+00 ) then
          info = j
          return
        end if

        a(j,j) = sqrt ( s )

      end do

      info = 0
c
c  Since the Cholesky factor is stored in R8GE format, be sure to
c  zero out the lower triangle.
c
      do i = 1, n
        do j = 1, i-1
          a(i,j) = 0.0D+00
        end do
      end do

      return
      end
      subroutine r8po_indicator ( n, a )

c*********************************************************************72
c
cc R8PO_INDICATOR sets up an R8PO indicator matrix.
c
c  Discussion:
c
c    The "indicator matrix" simply has a value like I*10+J at every
c    entry of a dense matrix, or at every entry of a compressed storage
c    matrix for which storage is allocated. 
c
c    The R8PO storage format is used for a symmetric positive definite 
c    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
c    upper triangular matrix, so it will be in R8GE storage format.)
c
c    Only the diagonal and upper triangle of the square array are used.
c    This same storage scheme is used when the matrix is factored by
c    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
c    is set to zero.
c
c    R8PO storage is used by LINPACK and LAPACK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 January 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of rows and columns of 
c    the matrix.  N must be positive.
c
c    Output, double precision A(N,N), the R8PO matrix.
c
      implicit none

      integer n

      double precision a(n,n)
      integer fac
      integer i
      integer i4_log_10
      integer j

      fac = 10 ** ( i4_log_10 ( n ) + 1 )

      do i = 1, n
        do j = 1, i-1
          a(i,j) = 0.0D+00
        end do
        do j = i, n
          a(i,j) = dble ( fac * i + j )
        end do
      end do

      return
      end
      subroutine r8po_inverse ( n, a )

c*********************************************************************72
c
cc R8PO_INVERSE computes the inverse of a matrix factored by R8PO_FA.
c
c  Discussion:
c
c    The R8PO storage format is used for a symmetric positive definite 
c    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
c    upper triangular matrix, so it will be in R8GE storage format.)
c
c    Only the diagonal and upper triangle of the square array are used.
c    This same storage scheme is used when the matrix is factored by
c    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
c    is set to zero.
c
c    R8PO storage is used by LINPACK and LAPACK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 March 1999
c
c  Author:
c
c    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
c    FORTRAN90 version by John Burkardt.
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
c
c    Input/output, double precision A(N,N).
c    On input, the Cholesky factor, in R8GE storage, returned by R8PO_FA.
c    On output, the inverse matrix, in R8PO storage.
c
      implicit none

      integer n

      double precision a(n,n)
      integer j
      integer k
      double precision t
c
c  Compute Inverse ( R ).
c
      do k = 1, n

        a(k,k) = 1.0D+00 / a(k,k)
        a(1:k-1,k) = -a(1:k-1,k) * a(k,k)

        do j = k + 1, n
          t = a(k,j)
          a(k,j) = 0.0D+00
          a(1:k,j) = a(1:k,j) + t * a(1:k,k)
        end do

      end do
c
c  Compute Inverse ( R ) * ( Inverse ( R ) )'.
c
      do j = 1, n

        do k = 1, j - 1
          t = a(k,j)
          a(1:k,k) = a(1:k,k) + t * a(1:k,j)
        end do

        a(1:j,j) = a(1:j,j) * a(j,j)

      end do

      return
      end
      subroutine r8po_ml ( n, a_lu, x, b )

c*********************************************************************72
c
cc R8PO_ML computes A * x = b after A has been factored by R8PO_FA.
c
c  Discussion:
c
c    The R8PO storage format is used for a symmetric positive definite 
c    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
c    upper triangular matrix, so it will be in R8GE storage format.)
c
c    Only the diagonal and upper triangle of the square array are used.
c    This same storage scheme is used when the matrix is factored by
c    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
c    is set to zero.
c
c    R8PO storage is used by LINPACK and LAPACK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision A_LU(N,N), the Cholesky factor from R8PO_FA.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(N), the product A * x.
c
      implicit none

      integer n

      double precision a_lu(n,n)
      double precision b(n)
      integer i
      double precision x(n)
c
c  Compute R * x = y.
c
      do i = 1, n
        b(i) = a_lu(i,i) * x(i) + sum ( a_lu(i,i+1:n) * x(i+1:n) )
      end do
c
c  Compute R' * y = b.
c
      do i = n, 1, -1
        b(i) = a_lu(i,i) * b(i) + sum ( b(1:i-1) * a_lu(1:i-1,i) )
      end do

      return
      end
      subroutine r8po_mxm ( n, a, b, c )

c*********************************************************************72
c
cc R8PO_MXM multiplies two R8PO matrices.
c
c  Discussion:
c
c    The R8PO storage format is used for a symmetric positive definite 
c    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
c    upper triangular matrix, so it will be in R8GE storage format.)
c
c    Only the diagonal and upper triangle of the square array are used.
c    This same storage scheme is used when the matrix is factored by
c    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
c    is set to zero.
c
c    R8PO storage is used by LINPACK and LAPACK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 May 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrices.
c    N must be positive.
c
c    Input, double precision A(N,N), B(N,N), the R8PO factor matrices.
c
c    Output, double precision C(N,N), the R8PO product matrix.
c
      implicit none

      integer n

      double precision a(n,n)
      double precision aik
      double precision b(n,n)
      double precision bkj
      double precision c(n,n)
      integer i
      integer j
      integer k

      c(1:n,1:n) = 0.0D+00

      do i = 1, n

        do j = i, n
          do k = 1, n

            if ( i .le. k ) then
              aik = a(i,k)
            else
              aik = a(k,i)
            end if

            if ( k .le. j ) then
              bkj = b(k,j)
            else
              bkj = b(j,k)
            end if

            c(i,j) = c(i,j) + aik * bkj

          end do
        end do

      end do

      return
      end
      subroutine r8po_mxv ( n, a, x, b )

c*********************************************************************72
c
cc R8PO_MXV multiplies an R8PO matrix by an R8VEC.
c
c  Discussion:
c
c    The R8PO storage format is used for a symmetric positive definite 
c    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
c    upper triangular matrix, so it will be in R8GE storage format.)
c
c    Only the diagonal and upper triangle of the square array are used.
c    This same storage scheme is used when the matrix is factored by
c    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
c    is set to zero.
c
c    R8PO storage is used by LINPACK and LAPACK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 October 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision A(N,N), the R8PO matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(N), the product A * x.
c
      implicit none

      integer n

      double precision a(n,n)
      double precision b(n)
      integer i
      integer j
      double precision x(n)

      do i = 1, n
        b(i) = 0.0D+00
        do j = 1, i-1
          b(i) = b(i) + a(j,i) * x(j)
        end do
        do j = i, n
          b(i) = b(i) + a(i,j) * x(j)
        end do
      end do

      return
      end
      subroutine r8po_print ( n, a, title )

c*********************************************************************72
c
cc R8PO_PRINT prints an R8PO matrix.
c
c  Discussion:
c
c    The R8PO storage format is used for a symmetric positive definite 
c    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
c    upper triangular matrix, so it will be in R8GE storage format.)
c
c    Only the diagonal and upper triangle of the square array are used.
c    This same storage scheme is used when the matrix is factored by
c    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
c    is set to zero.
c
c    R8PO storage is used by LINPACK and LAPACK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 October 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision A(N,N), the R8PO matrix.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer n

      double precision a(n,n)
      character ( len = * ) title

      call r8po_print_some ( n, a, 1, 1, n, n, title )

      return
      end
      subroutine r8po_print_some ( n, a, ilo, jlo, ihi, jhi, title )

c*********************************************************************72
c
cc R8PO_PRINT_SOME prints some of an R8PO matrix.
c
c  Discussion:
c
c    The R8PO storage format is used for a symmetric positive definite 
c    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
c    upper triangular matrix, so it will be in R8GE storage format.)
c
c    Only the diagonal and upper triangle of the square array are used.
c    This same storage scheme is used when the matrix is factored by
c    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
c    is set to zero.
c
c    R8PO storage is used by LINPACK and LAPACK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 October 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision A(N,N), the R8PO matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer, parameter :: incx = 5
      integer n

      double precision a(n,n)
      double precision aij
      character ( len = 14 ) ctemp(incx)
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
      character ( len = * ) title

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

        write ( *, '(''  Col:  '',5a14)' ) ( ctemp(j2), j2 = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) '  ---'
c
c  Determine the range of the rows in this strip.
c
        i2lo = max ( ilo, 1 )
        i2hi = min ( ihi, n )

        do i = i2lo, i2hi
c
c  Print out (up to) 5 entries in row I, that lie in the current strip.
c
          do j2 = 1, inc

            j = j2lo - 1 + j2

            if ( i .le. j ) then
              aij = a(i,j)
            else
              aij = a(j,i)
            end if

            if ( r8_is_int ( aij ) ) then
              write ( ctemp(j2), '(f8.0,6x)' ) aij
            else
              write ( ctemp(j2), '(g14.6)' ) aij
            end if

          end do

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine r8po_random ( n, seed, a )

c*********************************************************************72
c
cc R8PO_RANDOM randomizes an R8PO matrix.
c
c  Discussion:
c
c    The R8PO storage format is used for a symmetric positive definite 
c    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
c    upper triangular matrix, so it will be in R8GE storage format.)
c
c    Only the diagonal and upper triangle of the square array are used.
c    This same storage scheme is used when the matrix is factored by
c    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
c    is set to zero.
c
c    R8PO storage is used by LINPACK and LAPACK.
c
c    The matrix computed here is not simply a set of random numbers in
c    the nonzero slots of the R8PO array.  It is also a positive definite
c    matrix.  It is computed by setting a "random" upper triangular
c    Cholesky factor R, and then computing A = R'*R.
c    The randomness is limited by the fact that all the entries of
c    R will be between 0 and 1.  A truly random R is only required
c    to have positive entries on the diagonal.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 October 2004
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
c    Input/output, integer SEED, a seed for the random 
c    number generator.
c
c    Output, double precision A(N,N), the R8PO matrix.
c
      implicit none

      integer n

      double precision a(n,n)
      double precision r8_uniform_01
      integer i
      integer j
      integer seed
c
c  Set the whole matrix to zero.
c
      a(1:n,1:n) = 0.0D+00

      do i = n, 1, -1
c
c  Set row I of R.
c
        do j = i, n
          a(i,j) = r8_uniform_01 ( seed )
        end do
c
c  Consider element J of row I, last to first.
c
        do j = n, i, -1
c
c  Add multiples of row I to lower elements of column J.
c
          a(i+1:j,j) = a(i+1:j,j) + a(i,i+1:j) * a(i,j)
c
c  Reset element J.
c
          a(i,j) = a(i,i) * a(i,j)

        end do
      end do

      return
      end
      subroutine r8po_sl ( n, a_lu, b )

c*********************************************************************72
c
cc R8PO_SL solves an R8PO system factored by R8PO_FA.
c
c  Discussion:
c
c    The R8PO storage format is used for a symmetric positive definite 
c    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
c    upper triangular matrix, so it will be in R8GE storage format.)
c
c    Only the diagonal and upper triangle of the square array are used.
c    This same storage scheme is used when the matrix is factored by
c    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
c    is set to zero.
c
c    R8PO storage is used by LINPACK and LAPACK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 March 1999
c
c  Author:
c
c    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
c    FORTRAN90 version by John Burkardt.
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
c
c    Input, double precision A_LU(N,N), the Cholesky factor from R8PO_FA.
c
c    Input/output, double precision B(N).
c    On input, the right hand side.
c    On output, the solution vector.
c
      implicit none

      integer n

      double precision a_lu(n,n)
      double precision b(n)
      integer k
c
c  Solve R' * y = b.
c
      do k = 1, n
        b(k) = ( b(k) - sum ( b(1:k-1) * a_lu(1:k-1,k) ) ) / a_lu(k,k)
      end do
c
c  Solve R * x = y.
c
      do k = n, 1, -1
        b(k) = b(k) / a_lu(k,k)
        b(1:k-1) = b(1:k-1) - a_lu(1:k-1,k) * b(k)
      end do

      return
      end
      subroutine r8po_to_r8ge ( n, a, b )

c*********************************************************************72
c
cc R8PO_TO_R8GE copies an R8PO matrix to an R8GE matrix.
c
c  Discussion:
c
c    The R8PO storage format is used for a symmetric positive definite 
c    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
c    upper triangular matrix, so it will be in R8GE storage format.)
c
c    Only the diagonal and upper triangle of the square array are used.
c    This same storage scheme is used when the matrix is factored by
c    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
c    is set to zero.
c
c    R8PO storage is used by LINPACK and LAPACK.
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 October 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision A(N,N), the R8PO matrix.
c
c    Output, double precision B(N,N), the R8GE matrix.
c
      implicit none

      integer n

      double precision a(n,n)
      double precision b(n,n)
      integer i
      integer j

      do i = 1, n
        do j = 1, n
          if ( i .le. j ) then
            b(i,j) = a(i,j)
          else
            b(i,j) = a(j,i)
          end if
        end do
      end do

      return
      end
      subroutine r8pp_det ( n, a_lu, det )

c*********************************************************************72
c
cc R8PP_DET computes the determinant of an R8PP matrix factored by R8PP_FA.
c
c  Discussion:
c
c    The R8PP storage format is used for a symmetric positive
c    definite matrix.  Only the upper triangle of the matrix is stored,
c    by successive partial columns, in an array of length (N*(N+1))/2,
c    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
c
c    R8PP storage is used by LINPACK and LAPACK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 October 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision A_LU((N*(N+1))/2), the LU factors from R8PO_FA.
c
c    Output, double precision DET, the determinant of A.
c
      implicit none

      integer n

      double precision a_lu((n*(n+1))/2)
      double precision det
      integer i
      integer k

      det = 1.0D+00

      k = 0
      do i = 1, n
        k = k + i
        det = det * a_lu(k)
      end do

      return
      end
      subroutine r8pp_fa ( n, a, info )

c*********************************************************************72
c
cc R8PP_FA factors an R8PP matrix.
c
c  Discussion:
c
c    The R8PP storage format is appropriate for a symmetric positive
c    definite matrix.  Only the upper triangle of the matrix is stored,
c    by successive partial columns, in an array of length (N*(N+1))/2,
c    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
c
c    R8PP storage is used by LINPACK and LAPACK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 November 2008
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
c    SIAM, (Society for Industrial and Applied Mathematics),
c    3600 University City Science Center,
c    Philadelphia, PA, 19104-2688.
c    ISBN 0-89871-172-X
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input/output, double precision A((N*(N+1))/2).  On input, an R8PP matrix A.
c    On output, an upper triangular matrix R, stored in packed form, 
c    so that A = R'*R.
c
c    Output, integer INFO, error flag.
c    0, for normal return.
c    K, if the leading minor of order K is not positive definite.
c
      implicit none

      integer n

      double precision a((n*(n+1))/2)
      integer i
      integer info
      integer j
      integer jj
      integer k
      integer kj
      integer kk
      double precision s
      double precision t

      info = 0
      jj = 0

      do j = 1, n

        s = 0.0D+00
        kj = jj
        kk = 0

        do k = 1, j-1

          kj = kj + 1
          t = a(kj)
          do i = 1, k-1
            t = t - a(kk+i) * a(jj+i)
          end do
          kk = kk + k
          t = t / a(kk)
          a(kj) = t
          s = s + t * t

        end do

        jj = jj + j
        s = a(jj) - s

        if ( s .le. 0.0D+00 ) then
          info = j
          return
        end if

        a(jj) = sqrt ( s )

      end do

      return
      end
      subroutine r8pp_indicator ( n, a )

c*********************************************************************72
c
cc R8PP_INDICATOR sets up an R8PP indicator matrix.
c
c  Discussion:
c
c    The "indicator matrix" simply has a value like I*10+J at every
c    entry of a dense matrix, or at every entry of a compressed storage
c    matrix for which storage is allocated. 
c
c    The R8PP storage format is appropriate for a symmetric positive
c    definite matrix.  Only the upper triangle of the matrix is stored,
c    by successive partial columns, in an array of length (N*(N+1))/2,
c    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
c
c    R8PP storage is used by LINPACK and LAPACK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 February 2004
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
c    Output, double precision A((N*(N+1))/2), the R8PP matrix.
c
      implicit none

      integer n

      double precision a((n*(n+1))/2)
      integer fac
      integer i
      integer i4_log_10
      integer j
      integer k

      fac = 10 ** ( i4_log_10 ( n ) + 1 )

      k = 0
      do j = 1, n
        do i = 1, j
          k = k + 1
          a(k) = dble ( fac * i + j )
        end do
      end do

      return
      end
      subroutine r8pp_mxv ( n, a, x, b )

c*********************************************************************72
c
cc R8PP_MXV multiplies an R8PP matrix by an R8VEC.
c
c  Discussion:
c
c    The R8PP storage format is appropriate for a symmetric positive
c    definite matrix.  Only the upper triangle of the matrix is stored,
c    by successive partial columns, in an array of length (N*(N+1))/2,
c    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
c
c    R8PP storage is used by LINPACK and LAPACK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 October 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision A((N*(N+1))/2), the R8PP matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(N), the product A * x.
c
      implicit none

      integer n

      double precision a((n*(n+1))/2)
      double precision b(n)
      integer i
      integer j
      integer k
      double precision x(n)

      do i = 1, n
        b(i) = 0.0D+00
        do j = 1, i-1
          k = j + ( i * ( i - 1 ) ) / 2
          b(i) = b(i) + a(k) * x(j)
        end do
        do j = i, n
          k = i + ( j * ( j - 1 ) ) / 2
          b(i) = b(i) + a(k) * x(j)
        end do
      end do

      return
      end
      subroutine r8pp_print ( n, a, title )

c*********************************************************************72
c
cc R8PP_PRINT prints an R8PP matrix.
c
c  Discussion:
c
c    The R8PP storage format is appropriate for a symmetric positive
c    definite matrix.  Only the upper triangle of the matrix is stored,
c    by successive partial columns, in an array of length (N*(N+1))/2,
c    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
c
c    R8PP storage is used by LINPACK and LAPACK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 May 2000
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
c    Input, double precision A((N*(N+1))/2), the R8PP matrix.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer n

      double precision a((n*(n+1))/2)
      character ( len = * ) title

      call r8pp_print_some ( n, a, 1, 1, n, n, title )

      return
      end
      subroutine r8pp_print_some ( n, a, ilo, jlo, ihi, jhi, title )

c*********************************************************************72
c
cc R8PP_PRINT_SOME prints some of an R8PP matrix.
c
c  Discussion:
c
c    The R8PP storage format is appropriate for a symmetric positive
c    definite matrix.  Only the upper triangle of the matrix is stored,
c    by successive partial columns, in an array of length (N*(N+1))/2,
c    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
c
c    R8PP storage is used by LINPACK and LAPACK.
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
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input, double precision A((N*(N+1))/2), the R8PP matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer, parameter :: incx = 5
      integer n

      double precision a((n*(n+1))/2)
      double precision aij
      character ( len = 14 ) ctemp(incx)
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
      character ( len = * ) title

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
        i2hi = min ( ihi, n )

        do i = i2lo, i2hi
c
c  Print out (up to) 5 entries in row I, that lie in the current strip.
c
          do j2 = 1, inc

            j = j2lo - 1 + j2

            if ( i .le. j ) then
              aij = a(i+(j*(j-1))/2)
            else
              aij = a(j+(i*(i-1))/2)
            end if

            if ( r8_is_int ( aij ) ) then
              write ( ctemp(j2), '(f8.0,6x)' ) aij
            else
              write ( ctemp(j2), '(g14.6)' ) aij
            end if

          end do

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine r8pp_random ( n, seed, a )

c*********************************************************************72
c
cc R8PP_RANDOM randomizes an R8PP matrix.
c
c  Discussion:
c
c    The R8PP storage format is appropriate for a symmetric positive
c    definite matrix.  Only the upper triangle of the matrix is stored,
c    by successive partial columns, in an array of length (N*(N+1))/2,
c    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
c
c    R8PP storage is used by LINPACK and LAPACK.
c
c    The matrix is computed by setting a "random" upper triangular
c    Cholesky factor R, and then computing A = R'*R.
c    The randomness is limited by the fact that all the entries of
c    R will be between 0 and 1.  A truly random R is only required
c    to have positive entries on the diagonal.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 October 1998
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
c    Input/output, integer SEED, a seed for the random 
c    number generator.
c
c    Output, double precision A((N*(N+1))/2), the R8PP matrix.
c
      implicit none

      integer n

      double precision a((n*(n+1))/2)
      double precision r8_uniform_01
      integer i
      integer ii
      integer ij
      integer ik
      integer j
      integer k
      integer kj
      integer seed

      a(1:(n*(n+1))/2) = 0.0D+00

      do i = n, 1, -1
c
c  Set row I of R.
c
        do j = i, n
          ij = i + ( j * ( j - 1 ) ) / 2
          a(ij) = r8_uniform_01 ( seed )
        end do
c
c  Consider element J of row I, last to first.
c
        do j = n, i, -1
c
c  Add multiples of row I to lower elements of column J.
c
          ij = i + ( j * ( j - 1 ) ) / 2

          do k = i+1, j
            kj = k + ( j * ( j - 1 ) ) / 2
            ik = i + ( k * ( k - 1 ) ) / 2
            a(kj) = a(kj) + a(ik) * a(ij)
          end do
c
c  Reset element J.
c
          ii = i + ( i * ( i - 1 ) ) / 2
          a(ij) = a(ii) * a(ij)

        end do
      end do

      return
      end
      subroutine r8pp_sl ( n, a_lu, b )

c*********************************************************************72
c
cc R8PP_SL solves an R8PP system factored by R8PP_FA.
c
c  Discussion:
c
c    The R8PP storage format is appropriate for a symmetric positive
c    definite matrix.  Only the upper triangle of the matrix is stored,
c    by successive partial columns, in an array of length (N*(N+1))/2,
c    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
c
c    R8PP storage is used by LINPACK and LAPACK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 November 2008
c
c  Author:
c
c    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
c    FORTRAN90 version by John Burkardt.
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, (Society for Industrial and Applied Mathematics),
c    3600 University City Science Center,
c    Philadelphia, PA, 19104-2688.
c    ISBN 0-89871-172-X
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision A_LU((N*(N+1))/2), the LU factors from R8PP_FA.
c
c    Input/output, double precision B(N).  On input, the right hand side.
c    On output, the solution.
c
      implicit none

      integer n

      double precision a_lu((n*(n+1))/2)
      double precision b(n)
      integer i
      integer k
      integer kk
      double precision t

      kk = 0

      do k = 1, n
        t = 0.0D+00
        do i = 1, k-1
          t = t + a_lu(kk+i) * b(i)
        end do
        kk = kk + k
        b(k) = ( b(k) - t ) / a_lu(kk)
      end do

      do k = n, 1, -1
        b(k) = b(k) / a_lu(kk)
        kk = kk - k
        t = -b(k)
        do i = 1, k-1
          b(i) = b(i) + t * a_lu(kk+i)
        end do
      end do

      return
      end
      subroutine r8pp_to_r8ge ( n, a, b )

c*********************************************************************72
c
cc R8PP_TO_R8GE copies an R8PP matrix to an R8GE matrix.
c
c  Discussion:
c
c    The R8PP storage format is appropriate for a symmetric positive
c    definite matrix.  Only the upper triangle of the matrix is stored,
c    by successive partial columns, in an array of length (N*(N+1))/2,
c    which contains (A11,A12,A22,A13,A23,A33,A14,...,ANN)  
c
c    R8PP storage is used by LINPACK and LAPACK.
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 September 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision A((N*(N+1))/2), the R8PP matrix.
c
c    Output, double precision B(N,N), the R8GE matrix.
c
      implicit none

      integer n

      double precision a((n*(n+1))/2)
      double precision b(n,n)
      integer i
      integer j

      do i = 1, n
        do j = 1, n
          if ( i .le. j ) then
            b(i,j) = a(i+(j*(j-1))/2)
          else
            b(i,j) = a(j+(i*(i-1))/2)
          end if
        end do
      end do

      return
      end
      subroutine r8row_swap ( m, n, a, i1, i2 )

c*********************************************************************72
c
cc R8ROW_SWAP swaps two rows of an R8ROW.
c
c  Discussion:
c
c    An R8ROW is an M by N array of R8 values, regarded
c    as an array of M rows of length N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 December 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input/output, double precision A(M,N), the M by N array.
c
c    Input, integer I1, I2, the two rows to swap.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer i1
      integer i2
      double precision row(n)

      if ( i1 .lt. 1 .or. m .lt. i1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8ROW_SWAP - Fatal error!'
        write ( *, '(a)' ) '  I1 is out of range.'
        write ( *, '(a,i8)' ) '  I1 = ', i1
        stop
      end if

      if ( i2 .lt. 1 .or. m .lt. i2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8ROW_SWAP - Fatal error!'
        write ( *, '(a)' ) '  I2 is out of range.'
        write ( *, '(a,i8)' ) '  I2 = ', i2
        stop
      end if

      if ( i1 .eq. i2 ) then
        return
      end if

      row(1:n) = a(i1,1:n)
      a(i1,1:n) = a(i2,1:n)
      a(i2,1:n) = row(1:n)

      return
      end
      subroutine r8s3_diagonal ( n, nz_num, isym, row, col, a )

c*********************************************************************72
c
cc R8S3_DIAGONAL reorders a square R8S3 matrix so the diagonal entries are first.
c
c  Discussion:
c
c    The R8S3 storage format corresponds to the DLAP/SLAP Triad format.
c
c    The R8S3 storage format stores the row, column and value of each nonzero
c    entry of a sparse matrix.  The entries may be given in any order.  No
c    check is made for the erroneous case in which a given matrix entry is
c    specified more than once.
c
c    There is a symmetry option for square matrices.  If the symmetric storage
c    option is used, the format specifies that only nonzeroes on the diagonal
c    and lower triangle are stored.  However, this routine makes no attempt
c    to enforce this.  The only thing it does is to "reflect" any nonzero
c    offdiagonal value.  Moreover, no check is made for the erroneous case
c    in which both A(I,J) and A(J,I) are specified, but with different values.
c
c    This routine reorders the entries of A so that the first N entries
c    are exactly the diagonal entries of the matrix, in order.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 November 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero elements in 
c    the matrix.
c
c    Input, integer ISYM, is 0 if the matrix is not symmetric, 
c    and 1 if the matrix is symmetric.  If the matrix is symmetric, then
c    only the nonzeroes on the diagonal and in the lower triangle are stored.
c
c    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and 
c    column indices of the nonzero elements.
c
c    Input/output, double precision A(NZ_NUM), the nonzero elements 
c    of the matrix.
c
      implicit none

      integer n
      integer nz_num

      double precision a(nz_num)
      integer col(nz_num)
      integer found
      integer i
      integer isym
      integer k
      integer row(nz_num)

      found = 0

      do k = 1, nz_num

        do while ( row(k) .eq. col(k) )

          if ( row(k) .eq. k ) then
            found = found + 1
            exit
          end if

          i = row(k)

          call i4_swap ( row(i), row(k) )
          call i4_swap ( col(i), col(k) )
          call r8_swap ( a(i), a(k) )
     
          found = found + 1

          if ( n .le. found ) then
            exit
          end if
         
        end do

        if ( n .le. found ) then
          exit
        end if

      end do

      if ( found .lt. n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8S3_DIAGONAL - Warning!'
        write ( *, '(a,i8)' ) 
     &    '  Number of diagonal entries expected was ', n
        write ( *, '(a,i8)' ) '  Number found was ', found
      end if

      return
      end
      subroutine r8s3_indicator ( n, nz_num, isym, row, col, a )

c*********************************************************************72
c
cc R8S3_INDICATOR sets up an R8S3 indicator matrix.
c
c  Discussion:
c
c    The "indicator matrix" simply has a value like I*10+J at every
c    entry of a dense matrix, or at every entry of a compressed storage
c    matrix for which storage is allocated. 
c
c    The R8S3 storage format corresponds to the DLAP/SLAP Triad format.
c
c    The R8S3 storage format stores the row, column and value of each nonzero
c    entry of a sparse matrix.  The entries may be given in any order.  No
c    check is made for the erroneous case in which a given matrix entry is
c    specified more than once.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    14 January 2004
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
c    Input, integer NZ_NUM, the number of nonzero entries.
c
c    Input, integer ISYM, is 0 if the matrix is not symmetric, 
c    and 1 if the matrix is symmetric.  If the matrix is symmetric, then
c    only the nonzeroes on the diagonal and in the lower triangle are stored.
c
c    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and column
c    indices of the nonzero elements.
c
c    Output, double precision A(NZ_NUM), the indicator matrix.
c
      implicit none

      integer n
      integer nz_num

      double precision a(nz_num)
      integer col(nz_num)
      integer fac
      integer i
      integer i4_log_10
      integer isym
      integer j
      integer k
      integer row(nz_num)

      fac = 10 ** ( i4_log_10 ( n ) + 1 )

      do k = 1, nz_num
        i = row(k)
        j = col(k)
        a(k) = dble ( fac * i + j )
      end do

      return
      end
      subroutine r8s3_jac_sl ( n, nz_num, isym, row, col, a, b, x, 
     &  tol, it_max, job, it, diff )

c*********************************************************************72
c
cc R8S3_JAC_SL solves an R8S3 system using Jacobi iteration.
c
c  Discussion:
c
c    The R8S3 storage format corresponds to the DLAP/SLAP Triad format.
c
c    The R8S3 storage format stores the row, column and value of each nonzero
c    entry of a sparse matrix.  The entries may be given in any order.  No
c    check is made for the erroneous case in which a given matrix entry is
c    specified more than once.
c
c    There is a symmetry option for square matrices.  If the symmetric storage
c    option is used, the format specifies that only nonzeroes on the diagonal
c    and lower triangle are stored.  However, this routine makes no attempt
c    to enforce this.  The only thing it does is to "reflect" any nonzero
c    offdiagonal value.  Moreover, no check is made for the erroneous case
c    in which both A(I,J) and A(J,I) are specified, but with different values.
c
c    This routine REQUIRES that the matrix be square, that the matrix
c    have nonzero diagonal entries, and that the first N entries of
c    the array A be exactly the diagonal entries of the matrix, in order.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 November 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero elements in 
c    the matrix.
c
c    Input, integer ISYM, is 0 if the matrix is not symmetric, 
c    and 1 if the matrix is symmetric.  If the matrix is symmetric, then
c    only the nonzeroes on the diagonal and in the lower triangle are stored.
c
c    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and column 
c    indices of the nonzero elements.
c
c    Input, double precision A(NZ_NUM), the nonzero elements of the matrix.
c
c    Input, double precision B(N), the right hand side of the linear system.
c
c    Input/output, double precision X(N), an approximate solution 
c    to the system.
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
      integer nz_num

      double precision a(nz_num)
      double precision b(n)
      integer col(nz_num)
      double precision diff
      integer i
      integer isym
      integer it
      integer it_max
      integer it_num
      integer j
      integer job
      integer k
      integer row(nz_num)
      double precision tol
      double precision x(n)
      double precision x_new(n)
      double precision x_norm

      if ( job .eq. 0 ) then

        do it_num = 1, it_max

          it = it_num
c
c  Initialize to right hand side.
c
          x_new(1:n) = b(1:n)
c
c  Subtract off-diagonal terms.
c
          do k = n+1, nz_num
            i = row(k)
            j = col(k)
            x_new(i) = x_new(i) - a(k) * x(j)
            if ( isym .eq. 1 ) then
              x_new(j) = x_new(j) - a(k) * x(i)
            end if
          end do
c
c  Divide by diagonal terms.
c
          x_new(1:n) = x_new(1:n) / a(1:n)
c
c  Measure change:
c
          x_norm = maxval ( abs ( x(1:n) ) )
          diff = maxval ( abs ( x_new(1:n) - x(1:n) ) )
c
c  Update.
c
          x(1:n) = x_new(1:n)
c
c  Test for early termination.
c
          if ( diff .le. tol * ( x_norm + 1.0D+00 ) ) then
            exit
          end if

        end do

      else

        do it_num = 1, it_max

          it = it_num
c
c  Initialize to right hand side.
c
          x_new(1:n) = b(1:n)
c
c  Subtract off-diagonal terms.
c
          do k = n+1, nz_num
            i = col(k)
            j = row(k)
            x_new(i) = x_new(i) - a(k) * x(j)
            if ( isym .eq. 1 ) then
              x_new(j) = x_new(j) - a(k) * x(i)
            end if
          end do
c
c  Divide by diagonal terms.
c
          x_new(1:n) = x_new(1:n) / a(1:n)
c
c  Measure change:
c
          x_norm = maxval ( abs ( x(1:n) ) )
          diff = maxval ( abs ( x_new(1:n) - x(1:n) ) )
c
c  Update.
c
          x(1:n) = x_new(1:n)
c
c  Test for early termination.
c
          if ( diff .le. tol * ( x_norm + 1.0D+00 ) ) then
            exit
          end if

        end do

      end if

      return
      end
      subroutine r8s3_mxv ( m, n, nz_num, isym, row, col, a, x, b )

c*********************************************************************72
c
cc R8S3_MXV multiplies an R8S3 matrix by an R8VEC.
c
c  Discussion:
c
c    The R8S3 storage format corresponds to the DLAP/SLAP Triad format.
c
c    The R8S3 storage format stores the row, column and value of each nonzero
c    entry of a sparse matrix.  The entries may be given in any order.  No
c    check is made for the erroneous case in which a given matrix entry is
c    specified more than once.
c
c    There is a symmetry option for square matrices.  If the symmetric storage
c    option is used, the format specifies that only nonzeroes on the diagonal
c    and lower triangle are stored.  However, this routine makes no attempt
c    to enforce this.  The only thing it does is to "reflect" any nonzero
c    offdiagonal value.  Moreover, no check is made for the erroneous case
c    in which both A(I,J) and A(J,I) are specified, but with different values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 November 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c
c    Input, integer N, the number of columns of the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero elements in 
c    the matrix.
c
c    Input, integer ISYM, is 0 if the matrix is not symmetric, and 1
c    if the matrix is symmetric.  If the matrix is symmetric, then
c    only the nonzeroes on the diagonal and in the lower triangle are stored.
c
c    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and column 
c    indices of the nonzero elements.
c
c    Input, double precision A(NZ_NUM), the nonzero elements of the matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(M), the product A * x.
c
      implicit none

      integer m
      integer n
      integer nz_num

      double precision a(nz_num)
      double precision b(m)
      integer col(nz_num)
      integer i
      integer isym
      integer j
      integer k
      integer row(nz_num)
      double precision x(n)

      b(1:m) = 0.0D+00

      do k = 1, nz_num
        i = row(k)
        j = col(k)
        b(i) = b(i) + a(k) * x(j)
      end do
c
c  Handle the symmetric option.
c
      if ( isym .eq. 1 .and. m .eq. n ) then
        do k = 1, nz_num
          i = col(k)
          j = row(k)
          if ( i .ne. j ) then
            b(i) = b(i) + a(k) * x(j)
          end if
        end do
      end if

      return
      end
      subroutine r8s3_print ( m, n, nz_num, isym, row, col, a, title )

c*********************************************************************72
c
cc R8S3_PRINT prints an R8S3 matrix.
c
c  Discussion:
c
c    The R8S3 storage format corresponds to the DLAP/SLAP Triad format.
c
c    The R8S3 storage format stores the row, column and value of each nonzero
c    entry of a sparse matrix.  The entries may be given in any order.  No
c    check is made for the erroneous case in which a given matrix entry is
c    specified more than once.
c
c    There is a symmetry option for square matrices.  If the symmetric storage
c    option is used, the format specifies that only nonzeroes on the diagonal
c    and lower triangle are stored.  However, this routine makes no attempt
c    to enforce this.  The only thing it does is to "reflect" any nonzero
c    offdiagonal value.  Moreover, no check is made for the erroneous case
c    in which both A(I,J) and A(J,I) are specified, but with different values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns of 
c    the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero elements in 
c    the matrix.
c
c    Input, integer ISYM, is 0 if the matrix is not symmetric, 
c    and 1 if the matrix is symmetric.  The symmetric case only makes sense
c    if the matrix is also square, that is, M = N.  In this case, only
c    the nonzeroes on the diagonal and in the lower triangle are stored.
c
c    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and column
c    indices of the nonzero elements.
c
c    Input, double precision A(NZ_NUM), the nonzero elements of the matrix.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer nz_num

      double precision a(nz_num)
      integer col(nz_num)
      integer isym
      integer m
      integer n
      integer row(nz_num)
      character ( len = * ) title

      call r8s3_print_some ( m, n, nz_num, isym, row, col, a, 1, 1, 
     &  m, n, title )

      return
      end
      subroutine r8s3_print_some ( m, n, nz_num, isym, row, col, a, 
     &  ilo, jlo, ihi, jhi, title )

c*********************************************************************72
c
cc R8S3_PRINT_SOME prints some of an R8S3 matrix.
c
c  Discussion:
c
c    The R8S3 storage format corresponds to the DLAP/SLAP Triad format.
c
c    The R8S3 storage format stores the row, column and value of each nonzero
c    entry of a sparse matrix.  The entries may be given in any order.  No
c    check is made for the erroneous case in which a given matrix entry is
c    specified more than once.
c
c    There is a symmetry option for square matrices.  If the symmetric storage
c    option is used, the format specifies that only nonzeroes on the diagonal
c    and lower triangle are stored.  However, this routine makes no attempt
c    to enforce this.  The only thing it does is to "reflect" any nonzero
c    offdiagonal value.  Moreover, no check is made for the erroneous case
c    in which both A(I,J) and A(J,I) are specified, but with different values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns of
c    the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero elements in 
c    the matrix.
c
c    Input, integer ISYM, is 0 if the matrix is not symmetric, and 1
c    if the matrix is symmetric.  The symmetric case only makes sense
c    if the matrix is also square, that is, M = N.  In this case, only
c    the nonzeroes on the diagonal and in the lower triangle are stored.
c
c    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and column
c    indices of the nonzero elements.
c
c    Input, double precision A(NZ_NUM), the nonzero elements of the matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer, parameter :: incx = 5
      integer nz_num

      double precision a(nz_num)
      double precision aij
      integer col(nz_num)
      character ( len = 14 ) ctemp(incx)
      logical r8_is_int
      integer i
      integer i2hi
      integer i2lo
      integer ihi
      integer ilo
      integer inc
      integer isym
      integer j
      integer j2
      integer j2hi
      integer j2lo
      integer jhi
      integer jlo
      integer k
      integer m
      integer n
      logical nonzero
      integer row(nz_num)
      character ( len = * ) title

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

        write ( *, '(''  Col:  '',5a14)' ) ( ctemp(j2), j2 = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) '  ---'
c
c  Determine the range of the rows in this strip.
c
        i2lo = max ( ilo, 1 )
        i2hi = min ( ihi, m )

        do i = i2lo, i2hi
c
c  Print out (up to) 5 entries in row I, that lie in the current strip.
c
          nonzero = .false.

          aij = 0.0D+00
          do j2 = 1, inc
            write ( ctemp(j2), '(f8.0,6x)' ) aij
          end do

          do k = 1, nz_num

            if ( i .eq. row(k) .and. j2lo .le. col(k) .and. 
     &        col(k) .le. j2hi ) then

              j2 = col(k) - j2lo + 1
              aij = a(k)

              if ( aij .eq. 0.0D+00 ) then
                cycle
              end if

              nonzero = .true.

              if ( r8_is_int ( aij ) ) then
                write ( ctemp(j2), '(f8.0,6x)' ) aij
              else
                write ( ctemp(j2), '(g14.6)' ) aij
              end if

            else if ( isym .eq. 1 .and. m .eq. n .and. 
     &        i .eq. col(k) .and. j2lo .le. row(k) .and. 
     &        row(k) .le. j2hi ) then

              j2 = row(k) - j2lo + 1
              aij = a(k)

              if ( aij .eq. 0.0D+00 ) then
                cycle
              end if

              nonzero = .true.

              if ( r8_is_int ( aij ) ) then
                write ( ctemp(j2), '(f8.0,6x)' ) aij
              else
                write ( ctemp(j2), '(g14.6)' ) aij
              end if

            end if

          end do

          if ( nonzero ) then
            write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )
          end if

        end do

      end do

      return
      end
      subroutine r8s3_read ( input_file, n, nz_num, row, col, a )

c*********************************************************************72
c
cc R8S3_READ reads a square R8S3 matrix from a file.
c
c  Discussion:
c
c    This routine needs the value of NZ_NUM, which can be determined
c    by a call to R8S3_READ_SIZE.
c
c    The R8S3 storage format corresponds to the DLAP/SLAP Triad format.
c
c    The R8S3 storage format stores the row, column and value of each nonzero
c    entry of a sparse matrix.  The entries may be given in any order.  No
c    check is made for the erroneous case in which a given matrix entry is
c    specified more than once.
c
c    There is a symmetry option for square matrices.  If the symmetric storage
c    option is used, the format specifies that only nonzeroes on the diagonal
c    and lower triangle are stored.  However, this routine makes no attempt
c    to enforce this.  The only thing it does is to "reflect" any nonzero
c    offdiagonal value.  Moreover, no check is made for the erroneous case
c    in which both A(I,J) and A(J,I) are specified, but with different values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    31 July 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character ( len = * ) INPUT_FILE, the name of the file to be read.
c
c    Unused, integer N, the order of the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero elements in 
c    the matrix.
c
c    Output, integer ROW(NZ_NUM), COL(NZ_NUM), the row and 
c    column indices of the nonzero elements.
c
c    Output, double precision A(NZ_NUM), the nonzero elements 
c    of the matrix.
c
      implicit none

      integer nz_num

      double precision a(nz_num)
      integer col(nz_num)
      character ( len = * ) input_file
      integer input_unit
      integer ios
      integer isym
      integer k
      integer n
      integer row(nz_num)

      call get_unit ( input_unit )

      open ( unit = input_unit, file = input_file, status = 'old', 
     &  iostat = ios )

      if ( ios .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8S3_READ - Fatal error!'
        write ( *, '(a)' ) '  Could not open the input file "' 
     &    // trim ( input_file ) // '".'
        stop
      end if

      do k = 1, nz_num

        read ( input_unit, *, iostat = ios ) row(k), col(k), a(k)

        if ( ios .ne. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8S3_READ - Fatal error!'
          write ( *, '(a,i8)' ) '  I/O error while reading record ', k
          stop
        end if

      end do

      close ( unit = input_unit )

      return
      end
      subroutine r8s3_read_size ( input_file, n, nz_num )

c*********************************************************************72
c
cc R8S3_READ_SIZE reads the size of a square R8S3 matrix from a file.
c
c  Discussion:
c
c    The value of NZ_NUM is simply the number of records in the input file.
c
c    The value of N is determined as the maximum entry in the row and column
c    vectors.
c
c    The R8S3 storage format corresponds to the DLAP/SLAP Triad format.
c
c    The R8S3 storage format stores the row, column and value of each nonzero
c    entry of a sparse matrix.  The entries may be given in any order.  No
c    check is made for the erroneous case in which a given matrix entry is
c    specified more than once.
c
c    There is a symmetry option for square matrices.  If the symmetric storage
c    option is used, the format specifies that only nonzeroes on the diagonal
c    and lower triangle are stored.  However, this routine makes no attempt
c    to enforce this.  The only thing it does is to "reflect" any nonzero
c    offdiagonal value.  Moreover, no check is made for the erroneous case
c    in which both A(I,J) and A(J,I) are specified, but with different values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    31 July 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character ( len = * ) INPUT_FILE, the name of the file to 
c    be read.
c
c    Output, integer N, the order of the matrix.
c
c    Output, integer NZ_NUM, the number of nonzero elements 
c    in the matrix.
c
      implicit none

      double precision a_k
      integer col_k
      character ( len = * ) input_file
      integer input_unit
      integer ios
      integer k
      integer n
      integer nz_num
      integer row_k

      call get_unit ( input_unit )

      open ( unit = input_unit, file = input_file, status = 'old', 
     &  iostat = ios )

      if ( ios .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8S3_READ_SIZE - Fatal error!'
        write ( *, '(a)' ) '  Could not open the input file "' 
     &    // trim ( input_file ) // '".'
        stop
      end if

      nz_num = 0
      n = 0

      do

        read ( input_unit, *, iostat = ios ) row_k, col_k, a_k

        if ( ios .ne. 0 ) then
          exit
        end if

        nz_num = nz_num + 1
        n = max ( n, row_k )
        n = max ( n, col_k )

      end do

      close ( unit = input_unit )

      return
      end
      subroutine r8s3_to_r8ge ( m, n, nz_num, isym, row, col, a, b )

c*********************************************************************72
c
cc R8S3_TO_R8GE copies an R8S3 matrix to an R8GE matrix.
c
c  Discussion:
c
c    The R8S3 storage format corresponds to the DLAP/SLAP Triad format.
c
c    The R8S3 storage format stores the row, column and value of each nonzero
c    entry of a sparse matrix.  The entries may be given in any order.  No
c    check is made for the erroneous case in which a given matrix entry is
c    specified more than once.
c
c    There is a symmetry option for square matrices.  If the symmetric storage
c    option is used, the format specifies that only nonzeroes on the diagonal
c    and lower triangle are stored.  However, this routine makes no attempt
c    to enforce this.  The only thing it does is to "reflect" any nonzero
c    offdiagonal value.  Moreover, no check is made for the erroneous case
c    in which both A(I,J) and A(J,I) are specified, but with different values.
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 November 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns of 
c    the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero elements in 
c    the matrix.
c
c    Input, integer ISYM, is 0 if the matrix is not symmetric, 
c    and 1 if the matrix is symmetric.  The symmetric case only makes sense
c    if the matrix is also square, that is, M = N.  In this case, only
c    the nonzeroes on the diagonal and in the lower triangle are stored.
c
c    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and column
c    indices of the nonzero elements.
c
c    Input, double precision A(NZ_NUM), the nonzero elements of the matrix.
c
c    Output, double precision B(M,N), the R8GE matrix.
c
      implicit none

      integer m
      integer n
      integer nz_num

      double precision a(nz_num)
      double precision b(m,n)
      integer col(nz_num)
      integer i
      integer isym
      integer j
      integer k
      integer row(nz_num)

      b(1:m,1:n) = 0.0D+00

      do k = 1, nz_num
        i = row(k)
        j = col(k)
        b(i,j) = b(i,j) + a(k)
        if ( isym .eq. 1 .and. m .eq. n .and. i .ne. j ) then
          b(j,i) = b(j,i) + a(k)
        end if
      end do

      return
      end
      subroutine r8s3_vxm ( m, n, nz_num, isym, row, col, a, x, b )

c*********************************************************************72
c
cc R8S3_VXM multiplies an R8VEC times an R8S3 matrix.
c
c  Discussion:
c
c    The R8S3 storage format corresponds to the DLAP/SLAP Triad format.
c
c    The R8S3 storage format stores the row, column and value of each nonzero
c    entry of a sparse matrix.  The entries may be given in any order.  No
c    check is made for the erroneous case in which a given matrix entry is
c    specified more than once.
c
c    There is a symmetry option for square matrices.  If the symmetric storage
c    option is used, the format specifies that only nonzeroes on the diagonal
c    and lower triangle are stored.  However, this routine makes no attempt
c    to enforce this.  The only thing it does is to "reflect" any nonzero
c    offdiagonal value.  Moreover, no check is made for the erroneous case
c    in which both A(I,J) and A(J,I) are specified, but with different values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 November 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c
c    Input, integer N, the number of columns of the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero elements in
c    the matrix.
c
c    Input, integer ISYM, is 0 if the matrix is not symmetric, and 1
c    if the matrix is symmetric.  If the matrix is symmetric, then
c    only the nonzeroes on the diagonal and in the lower triangle are stored.
c
c    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and column
c    indices of the nonzero elements.
c
c    Input, double precision A(NZ_NUM), the nonzero elements of the matrix.
c
c    Input, double precision X(M), the vector to be multiplied by A'.
c
c    Output, double precision B(N), the product A' * x.
c
      implicit none

      integer m
      integer n
      integer nz_num

      double precision a(nz_num)
      double precision b(n)
      integer col(nz_num)
      integer i
      integer isym
      integer j
      integer k
      integer row(nz_num)
      double precision x(m)

      b(1:n) = 0.0D+00

      do k = 1, nz_num
        i = col(k)
        j = row(k)
        b(i) = b(i) + a(k) * x(j)
      end do
c
c  Handle the symmetric option.
c
      if ( isym .eq. 1 .and. m .eq. n ) then
        do k = 1, nz_num
          i = row(k)
          j = col(k)
          if ( i .ne. j ) then
            b(i) = b(i) + a(k) * x(j)
          end if
        end do
      end if

      return
      end
      subroutine r8s3_write ( n, nz_num, isym, row, col, a, 
     &  output_file )

c*********************************************************************72
c
cc R8S3_WRITE writes a square R8S3 matrix to a file.
c
c  Discussion:
c
c    The R8S3 storage format corresponds to the DLAP/SLAP Triad format.
c
c    The R8S3 storage format stores the row, column and value of each nonzero
c    entry of a sparse matrix.  The entries may be given in any order.  No
c    check is made for the erroneous case in which a given matrix entry is
c    specified more than once.
c
c    There is a symmetry option for square matrices.  If the symmetric storage
c    option is used, the format specifies that only nonzeroes on the diagonal
c    and lower triangle are stored.  However, this routine makes no attempt
c    to enforce this.  The only thing it does is to "reflect" any nonzero
c    offdiagonal value.  Moreover, no check is made for the erroneous case
c    in which both A(I,J) and A(J,I) are specified, but with different values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 September 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero elements in 
c    the matrix.
c
c    Input, integer ISYM, is 0 if the matrix is not symmetric, and 1
c    if the matrix is symmetric.  If the matrix is symmetric, then
c    only the nonzeroes on the diagonal and in the lower triangle are stored.
c
c    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and column 
c    indices of the nonzero elements.
c
c    Input, double precision A(NZ_NUM), the nonzero elements 
c    of the matrix.
c
c    Input, character ( len = * ) OUTPUT_FILE, the name of the file to which
c    the information is to be written.
c
      implicit none

      integer nz_num

      double precision a(nz_num)
      integer col(nz_num)
      integer ios
      integer isym
      integer k
      integer n
      character ( len = * ) output_file
      integer output_unit
      integer row(nz_num)

      call get_unit ( output_unit )

      open ( unit = output_unit, file = output_file, status = 'replace', 
     &  iostat = ios )

      if ( ios .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8S3_WRITE - Fatal error!'
        write ( *, '(a)' ) '  Could not open the output file "' 
     &    // trim ( output_file ) // '".'
        stop
      end if

      do k = 1, nz_num
        write ( output_unit, '(2x,i8,2x,i8,2x,g16.8)' ) 
     &    row(k), col(k), a(k)
      end do

      close ( unit = output_unit )

      return
      end
      subroutine r8sd_cg ( n, ndiag, offset, a, b, x )

c*********************************************************************72
c
cc R8SD_CG uses the conjugate gradient method on an R8SD linear system.
c
c  Discussion:
c
c    The R8SD storage format is for symmetric matrices whose only nonzero entries
c    occur along a few diagonals, but for which these diagonals are not all
c    close enough to the main diagonal for band storage to be efficient.
c
c    In that case, we assign the main diagonal the offset value 0, and 
c    each successive superdiagonal gets an offset value 1 higher, until
c    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
c
c    Assuming there are NDIAG nonzero diagonals (ignoring subdiagonalsc),
c    we then create an array B that has N rows and NDIAG columns, and simply
c    "collapse" the matrix A to the left:
c
c    For the conjugate gradient method to be applicable, the matrix A must 
c    be a positive definite symmetric matrix.
c
c    The method is designed to reach the solution to the linear system
c      A * x = b
c    after N computational steps.  However, roundoff may introduce
c    unacceptably large errors for some problems.  In such a case,
c    calling the routine a second time, using the current solution estimate
c    as the new starting guess, should result in improved results.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 October 1998
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    FS Beckman,
c    The Solution of Linear Equations by the Conjugate Gradient Method,
c    Mathematical Methods for Digital Computers, pages 62-72.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input, integer NDIAG, the number of diagonals that are stored.
c    NDIAG must be at least 1 and no more than N.
c
c    Input, integer OFFSET(NDIAG), the offsets for the diagonal
c    storage.
c
c    Input, double precision A(N,NDIAG), the R8SD matrix.
c
c    Input, double precision B(N), the right hand side vector.
c
c    Input/output, double precision X(N).
c    On input, an estimate for the solution, which may be 0.
c    On output, the approximate solution vector.  Note that repeated
c    calls to this routine, using the value of X output on the previous
c    call, MAY improve the solution.
c
      implicit none

      integer n
      integer ndiag

      double precision a(n,ndiag)
      double precision alpha
      double precision ap(n)
      double precision b(n)
      double precision beta
      integer it
      integer offset(ndiag)
      double precision p(n)
      double precision pap
      double precision pr
      double precision r(n)
      double precision rap
      double precision x(n)
c
c  Initialize
c    AP = A * x,
c    R  = b - A * x,
c    P  = b - A * x.
c
      call r8sd_mxv ( n, ndiag, offset, a, x, ap )

      r(1:n) = b(1:n) - ap(1:n)
      p(1:n) = b(1:n) - ap(1:n)
c
c  Do the N steps of the conjugate gradient method.
c
      do it = 1, n
c
c  Compute the matrix*vector product AP = A*P.
c
        call r8sd_mxv ( n, ndiag, offset, a, p, ap )
c
c  Compute the dot products
c    PAP = P*AP,
c    PR  = P*R
c  Set
c    ALPHA = PR / PAP.
c
        pap = sum ( p(1:n) * ap(1:n) )
        pr = sum ( p(1:n) * r(1:n) )

        if ( pap .eq. 0.0D+00 ) then
          return
        end if

        alpha = pr / pap
c
c  Set
c    X = X + ALPHA * P
c    R = R - ALPHA * AP.
c
        x(1:n) = x(1:n) + alpha * p(1:n)
        r(1:n) = r(1:n) - alpha * ap(1:n)
c
c  Compute the vector dot product
c    RAP = R*AP
c  Set
c    BETA = - RAP / PAP.
c
        rap = sum ( r(1:n) * ap(1:n) )

        beta = - rap / pap
c
c  Update the perturbation vector
c    P = R + BETA * P.
c
        p(1:n) = r(1:n) + beta * p(1:n)

      end do

      return
      end
      subroutine r8sd_indicator ( n, ndiag, offset, a )

c*********************************************************************72
c
cc R8SD_INDICATOR sets up an R8SD indicator matrix.
c
c  Discussion:
c
c    The "indicator matrix" simply has a value like I*10+J at every
c    entry of a dense matrix, or at every entry of a compressed storage
c    matrix for which storage is allocated. 
c
c    The R8SD storage format is for symmetric matrices whose only nonzero 
c    entries occur along a few diagonals, but for which these diagonals are not 
c    all close enough to the main diagonal for band storage to be efficient.
c
c    In that case, we assign the main diagonal the offset value 0, and 
c    each successive superdiagonal gets an offset value 1 higher, until
c    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
c
c    Assuming there are NDIAG nonzero diagonals (ignoring subdiagonalsc),
c    we then create an array B that has N rows and NDIAG columns, and simply
c    "collapse" the matrix A to the left:
c
c  Example:
c
c    The "offset" value is printed above each column.
c
c    Original matrix               New Matrix
c
c       0   1   2   3   4   5       0   1   3   5
c
c      11  12   0  14   0  16      11  12  14  16
c      21  22  23   0  25   0      22  23  25  --
c       0  32  33  34   0  36      33  34  36  --
c      41   0  43  44  45   0      44  45  --  --
c       0  52   0  54  55  56      55  56  --  --
c      61   0  63   0  65  66      66  --  --  --
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    14 January 2004
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
c    Input, integer NDIAG, the number of diagonals that are stored.
c    NDIAG must be at least 1 and no more than N.
c
c    Input, integer OFFSET(NDIAG), the offsets for the diagonal
c    storage.
c
c    Output, double precision A(N,NDIAG), the R8SD matrix.
c
      implicit none

      integer n
      integer ndiag

      double precision a(n,ndiag)
      integer diag
      integer fac
      integer i
      integer i4_log_10
      integer j
      integer offset(ndiag)

      fac = 10 ** ( i4_log_10 ( n ) + 1 )

      do i = 1, n
        do diag = 1, ndiag
          j = i + offset(diag)
          if ( 1 .le. j .and. j .le. n ) then
            a(i,diag) = dble ( fac * i + j )
          else
            a(i,diag) = 0.0D+00
          end if
        end do
      end do

      return
      end
      subroutine r8sd_mxv ( n, ndiag, offset, a, x, b )

c*********************************************************************72
c
cc R8SD_MXV multiplies an R8SD matrix by an R8VEC.
c
c  Discussion:
c
c    The R8SD storage format is for symmetric matrices whose only nonzero 
c    entries occur along a few diagonals, but for which these diagonals are not 
c    all close enough to the main diagonal for band storage to be efficient.
c
c    In that case, we assign the main diagonal the offset value 0, and 
c    each successive superdiagonal gets an offset value 1 higher, until
c    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
c
c    Assuming there are NDIAG nonzero diagonals (ignoring subdiagonalsc),
c    we then create an array B that has N rows and NDIAG columns, and simply
c    "collapse" the matrix A to the left:
c
c  Example:
c
c    The "offset" value is printed above each column.
c
c    Original matrix               New Matrix
c
c       0   1   2   3   4   5       0   1   3   5
c
c      11  12   0  14   0  16      11  12  14  16
c      21  22  23   0  25   0      22  23  25  --
c       0  32  33  34   0  36      33  34  36  --
c      41   0  43  44  45   0      44  45  --  --
c       0  52   0  54  55  56      55  56  --  --
c      61   0  63   0  65  66      66  --  --  --
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 October 1998
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
c    Input, integer NDIAG, the number of diagonals that are stored.
c    NDIAG must be at least 1 and no more than N.
c
c    Input, integer OFFSET(NDIAG), the offsets for the diagonal
c    storage.
c
c    Input, double precision A(N,NDIAG), the R8SD matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(N), the product A * x.
c
      implicit none

      integer n
      integer ndiag

      double precision a(n,ndiag)
      double precision b(n)
      integer i
      integer j
      integer jdiag
      integer offset(ndiag)
      double precision x(n)

      b(1:n) = 0.0D+00

      do i = 1, n
        do jdiag = 1, ndiag
          j = i + offset(jdiag)
          if ( 1 .le. j .and. j .le. n ) then
            b(i) = b(i) + a(i,jdiag) * x(j)
            if ( offset(jdiag) .ne. 0 ) then
              b(j) = b(j) + a(i,jdiag) * x(i)
            end if
          end if
        end do
      end do

      return
      end
      subroutine r8sd_print ( n, ndiag, offset, a, title )

c*********************************************************************72
c
cc R8SD_PRINT prints an R8SD matrix.
c
c  Discussion:
c
c    The R8SD storage format is for symmetric matrices whose only nonzero 
c    entries occur along a few diagonals, but for which these diagonals are not 
c    all close enough to the main diagonal for band storage to be efficient.
c
c    In that case, we assign the main diagonal the offset value 0, and 
c    each successive superdiagonal gets an offset value 1 higher, until
c    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
c
c    Assuming there are NDIAG nonzero diagonals (ignoring subdiagonalsc),
c    we then create an array B that has N rows and NDIAG columns, and simply
c    "collapse" the matrix A to the left:
c
c  Example:
c
c    The "offset" value is printed above each column.
c
c    Original matrix               New Matrix
c
c       0   1   2   3   4   5       0   1   3   5
c
c      11  12   0  14   0  16      11  12  14  16
c      21  22  23   0  25   0      22  23  25  --
c       0  32  33  34   0  36      33  34  36  --
c      41   0  43  44  45   0      44  45  --  --
c       0  52   0  54  55  56      55  56  --  --
c      61   0  63   0  65  66      66  --  --  --
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 May 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input, integer NDIAG, the number of diagonals of the matrix
c    that are stored in the array.
c    NDIAG must be at least 1, and no more than N.
c
c    Input, integer OFFSET(NDIAG), the offsets for the 
c    diagonal storage.
c
c    Input, double precision A(N,NDIAG), the R8SD matrix.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer n
      integer ndiag

      double precision a(n,ndiag)
      integer offset(ndiag)
      character ( len = * ) title

      call r8sd_print_some ( n, ndiag, offset, a, 1, 1, n, n, title )

      return
      end
      subroutine r8sd_print_some ( n, ndiag, offset, a, ilo, jlo, 
     &  ihi, jhi, title )

c*********************************************************************72
c
cc R8SD_PRINT_SOME prints some of an R8SD matrix.
c
c  Discussion:
c
c    The R8SD storage format is for symmetric matrices whose only nonzero 
c    entries occur along a few diagonals, but for which these diagonals are not 
c    all close enough to the main diagonal for band storage to be efficient.
c
c    In that case, we assign the main diagonal the offset value 0, and 
c    each successive superdiagonal gets an offset value 1 higher, until
c    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
c
c    Assuming there are NDIAG nonzero diagonals (ignoring subdiagonalsc),
c    we then create an array B that has N rows and NDIAG columns, and simply
c    "collapse" the matrix A to the left:
c
c  Example:
c
c    The "offset" value is printed above each column.
c
c    Original matrix               New Matrix
c
c       0   1   2   3   4   5       0   1   3   5
c
c      11  12   0  14   0  16      11  12  14  16
c      21  22  23   0  25   0      22  23  25  --
c       0  32  33  34   0  36      33  34  36  --
c      41   0  43  44  45   0      44  45  --  --
c       0  52   0  54  55  56      55  56  --  --
c      61   0  63   0  65  66      66  --  --  --
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
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input, integer NDIAG, the number of diagonals of the matrix
c    that are stored in the array.
c    NDIAG must be at least 1, and no more than N.
c
c    Input, integer OFFSET(NDIAG), the offsets for the diagonal
c    storage.
c
c    Input, double precision A(N,NDIAG), the R8SD matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer, parameter :: incx = 5
      integer n
      integer ndiag

      double precision a(n,ndiag)
      double precision aij
      character ( len = 14 ) ctemp(incx)
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
      integer jdiag
      integer jhi
      integer jlo
      integer off
      integer offset(ndiag)
      character ( len = * ) title

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
        i2hi = min ( ihi, n )

        do i = i2lo, i2hi
c
c  Print out (up to) 5 entries in row I, that lie in the current strip.
c
          do j2 = 1, inc

            j = j2lo - 1 + j2

            aij = 0.0D+00
            off = j - i
            do jdiag = 1, ndiag
              if ( off .eq. offset(jdiag) ) then
                aij = a(i,jdiag)
              else if ( off .eq. - offset(jdiag) ) then
                aij = a(j,jdiag)
              end if
            end do

            if ( r8_is_int ( aij ) ) then
              write ( ctemp(j2), '(f8.0,6x)' ) aij
            else
              write ( ctemp(j2), '(g14.6)' ) aij
            end if

          end do

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine r8sd_random ( n, ndiag, offset, seed, a )

c*********************************************************************72
c
cc R8SD_RANDOM randomizes an R8SD matrix.
c
c  Discussion:
c
c    The R8SD storage format is for symmetric matrices whose only nonzero 
c    entries occur along a few diagonals, but for which these diagonals are not 
c    all close enough to the main diagonal for band storage to be efficient.
c
c    In that case, we assign the main diagonal the offset value 0, and 
c    each successive superdiagonal gets an offset value 1 higher, until
c    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
c
c    Assuming there are NDIAG nonzero diagonals (ignoring subdiagonalsc),
c    we then create an array B that has N rows and NDIAG columns, and simply
c    "collapse" the matrix A to the left:
c
c  Example:
c
c    The "offset" value is printed above each column.
c
c    Original matrix               New Matrix
c
c       0   1   2   3   4   5       0   1   3   5
c
c      11  12   0  14   0  16      11  12  14  16
c      21  22  23   0  25   0      22  23  25  --
c       0  32  33  34   0  36      33  34  36  --
c      41   0  43  44  45   0      44  45  --  --
c       0  52   0  54  55  56      55  56  --  --
c      61   0  63   0  65  66      66  --  --  --
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 October 1998
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
c    Input, integer NDIAG, the number of diagonals that are stored.
c    NDIAG must be at least 1 and no more than N.
c
c    Input, integer OFFSET(NDIAG), the offsets for the diagonal
c    storage.
c
c    Input/output, integer SEED, a seed for the random number
c    generator.
c
c    Output, double precision A(N,NDIAG), the R8SD matrix.
c
      implicit none

      integer n
      integer ndiag

      double precision a(n,ndiag)
      double precision r8_uniform_01
      integer i
      integer j
      integer jj
      integer offset(ndiag)
      integer seed

      do i = 1, n
        do j = 1, ndiag
          jj = i + offset(j)
          if ( 1 .le. jj .and. jj .le. n ) then
            a(i,j) = r8_uniform_01 ( seed )
          else
            a(i,j) = 0.0D+00
          end if
        end do
      end do

      return
      end
      subroutine r8sd_to_r8ge ( n, ndiag, offset, a, b )

c*********************************************************************72
c
cc R8SD_TO_R8GE copies an R8SD matrix to an R8GE matrix.
c
c  Discussion:
c
c    The R8SD storage format is for symmetric matrices whose only nonzero 
c    entries occur along a few diagonals, but for which these diagonals are not 
c    all close enough to the main diagonal for band storage to be efficient.
c
c    In that case, we assign the main diagonal the offset value 0, and 
c    each successive superdiagonal gets an offset value 1 higher, until
c    the highest superdiagonal (the A(1,N) entry) is assigned the offset N-1.
c
c    Assuming there are NDIAG nonzero diagonals (ignoring subdiagonalsc),
c    we then create an array B that has N rows and NDIAG columns, and simply
c    "collapse" the matrix A to the left:
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c  Example:
c
c    The "offset" value is printed above each column.
c
c    Original matrix               New Matrix
c
c       0   1   2   3   4   5       0   1   3   5
c
c      11  12   0  14   0  16      11  12  14  16
c      21  22  23   0  25   0      22  23  25  --
c       0  32  33  34   0  36      33  34  36  --
c      41   0  43  44  45   0      44  45  --  --
c       0  52   0  54  55  56      55  56  --  --
c      61   0  63   0  65  66      66  --  --  --
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 October 1998
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
c    Input, integer NDIAG, the number of diagonals that are stored.
c    NDIAG must be at least 1 and no more than N.
c
c    Input, integer OFFSET(NDIAG), the offsets for the diagonal
c    storage.
c
c    Input, double precision A(N,NDIAG), the R8SD matrix.
c
c    Output, double precision B(N,N), the R8GE matrix.
c
      implicit none

      integer n
      integer ndiag

      double precision a(n,ndiag)
      double precision b(n,n)
      integer i
      integer j
      integer jj
      integer offset(ndiag)

      b(1:n,1:n) = 0.0D+00

      do i = 1, n
        do j = 1, ndiag
          jj = i + offset(j)
          if ( 1 .le. jj .and. jj .le. n ) then
            b(i,jj) = a(i,j)
            if ( i .ne. jj ) then
              b(jj,i) = a(i,j)
            end if
          end if
        end do
      end do

      return
      end
      subroutine r8sm_ml ( n, a_lu, u, v, pivot, x, b, job )

c*********************************************************************72
c
cc R8SM_ML multiplies a factored square R8SM matrix by an R8VEC.
c
c  Discussion:
c
c    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
c    which is defined by an M by N matrix A, an M vector U, and
c    an N vector V, by B = A - U * V'
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 January 2004
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
c    Input, double precision U(N), V(N), the Sherman Morrison vectors.
c
c    Input, integer PIVOT(N), the pivot vector computed by R8GE_FA.
c
c    Input, double precision X(N), the vector to be multiplied.
c
c    Output, double precision B(N), the result of the multiplication.
c
c    Input, integer JOB, specifies the operation to be done:
c    JOB = 0, compute (A-u*v') * x.
c    JOB nonzero, compute (A-u*v')' * x.
c
      implicit none

      integer n

      double precision a_lu(n,n)
      double precision b(n)
      integer pivot(n)
      integer job
      double precision u(n)
      double precision v(n)
      double precision x(n)

      call r8ge_ml ( n, a_lu, pivot, x, b, job )

      if ( job .eq. 0 ) then

        b(1:n) = b(1:n) - u(1:n) * sum ( v(1:n) * x(1:n) )

      else

        b(1:n) = b(1:n) - v(1:n) * sum ( u(1:n) * x(1:n) )

      end if

      return
      end
      subroutine r8sm_mxv ( m, n, a, u, v, x, b )

c*********************************************************************72
c
cc R8SM_MXV multiplies an R8SM matrix by an R8VEC.
c
c  Discussion:
c
c    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
c    which is defined by an M by N matrix A, an M vector U, and
c    an N vector V, by B = A - U * V'
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 January 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns 
c    of the matrix.
c
c    Input, double precision A(M,N), the R8SM matrix.
c
c    Input, double precision U(M), V(N), the R8SM vectors U and V.
c
c    Input, double precision X(N), the vector to be multiplied by (A-u*v').
c
c    Output, double precision B(M), the product (A-u*v') * x.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision b(m)
      double precision u(m)
      double precision v(n)
      double precision x(n)

      b(1:m) = matmul ( a(1:m,1:n), x(1:n) ) 
     &  - u(1:m) * sum ( v(1:n) * x(1:n) )

      return
      end
      subroutine r8sm_print ( m, n, a, u, v, title )

c*********************************************************************72
c
cc R8SM_PRINT prints an R8SM matrix.
c
c  Discussion:
c
c    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
c    which is defined by an M by N matrix A, an M vector U, and
c    an N vector V, by B = A - U * V'
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 January 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns of 
c    the matrix.
c
c    Input, double precision A(M,N), the R8SM matrix.
c
c    Input, double precision U(M), V(N), the R8SM vectors.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      character ( len = * ) title
      double precision u(m)
      double precision v(n)

      call r8sm_print_some ( m, n, a, u, v, 1, 1, m, n, title )

      return
      end
      subroutine r8sm_print_some ( m, n, a, u, v, ilo, jlo, ihi, jhi, 
     &  title )

c*********************************************************************72
c
cc R8SM_PRINT_SOME prints some of an R8SM matrix.
c
c  Discussion:
c
c    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
c    which is defined by an M by N matrix A, an M vector U, and
c    an N vector V, by B = A - U * V'
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 January 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns 
c    of the matrix.
c
c    Input, double precision A(M,N), the R8SM matrix.
c
c    Input, double precision U(M), V(N), the R8SM vectors.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer, parameter :: incx = 5
      integer m
      integer n

      double precision a(m,n)
      double precision aij
      character ( len = 14 ) ctemp(incx)
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
      double precision u(n)
      double precision v(n)
      character ( len = * ) title

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
        i2hi = min ( ihi, m )

        do i = i2lo, i2hi
c
c  Print out (up to) 5 entries in row I, that lie in the current strip.
c
          do j2 = 1, inc

            j = j2lo - 1 + j2

            aij = a(i,j) - u(i) * v(j)

            if ( r8_is_int ( aij ) ) then
              write ( ctemp(j2), '(f8.0,6x)' ) aij
            else
              write ( ctemp(j2), '(g14.6)' ) aij
            end if

          end do

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine r8sm_random ( m, n, seed, a, u, v )

c*********************************************************************72
c
cc R8SM_RANDOM randomizes an R8SM matrix.
c
c  Discussion:
c
c    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
c    which is defined by an M by N matrix A, an M vector U, and
c    an N vector V, by B = A - U * V'
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 October 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns 
c    of the matrix.
c
c    Input/output, integer SEED, a seed for the random 
c    number generator.
c
c    Output, double precision A(M,N), the R8SM matrix.
c
c    Output, double precision U(M), V(N), the R8SM vectors.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision r8_uniform_01
      integer i
      integer j
      integer seed
      double precision u(m)
      double precision v(n)

      do j = 1, n
        do i = 1, m
          a(i,j) = r8_uniform_01 ( seed )
        end do
      end do

      do i = 1, m
        u(i) = r8_uniform_01 ( seed )
      end do

      do j = 1, n
        v(j) = r8_uniform_01 ( seed )
      end do

      return
      end
      subroutine r8sm_sl ( n, a_lu, u, v, b, ierror, pivot, job )

c*********************************************************************72
c
cc R8SM_SL solves a square R8SM system that has been factored.
c
c  Discussion:
c
c    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
c    which is defined by an M by N matrix A, an M vector U, and
c    an N vector V, by B = A - U * V'
c
c    It is assumed that A has been decomposed into its LU factors
c    by R8GE_FA.  The Sherman Morrison formula allows
c    us to solve linear systems involving (A-u*v') by solving linear
c    systems involving A and adjusting the results.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 November 1998
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    David Kahaner, Cleve Moler, Steven Nash
c    Numerical Methods and Software,
c    Prentice Hall, 1989
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input, double precision A_LU(N,N), the LU factors from R8GE_FA.
c
c    Input, double precision U(N), V(N), the R8SM vectors U and V.
c
c    Input/output, double precision B(N).
c    On input, the right hand side vector.
c    On output, the solution vector.
c
c    Output, integer IERROR, an error flag.
c    0, no error occurred.  The solution was successfully computed.
c    1, an error occurred.  1 - v' * Inverse(A) * u = 0.
c    The solution was not computed.
c
c    Input, integer PIVOT(N), the pivot vector produced by R8GE_FA.
c
c    Input, integer JOB, specifies the system to solve.
c    0, solve (A-u*v') * X = B.
c    nonzero, solve (A-u*v') * X = B.
c
      implicit none

      integer n

      double precision a_lu(n,n)
      double precision alpha
      double precision b(n)
      double precision beta
      integer ierror
      integer pivot(n)
      integer job
      integer job_local
      double precision u(n)
      double precision v(n)
      double precision w(n)

      ierror = 0

      if ( job .eq. 0 ) then
c
c  Solve A' * w = v.
c
        w(1:n) = v(1:n)

        job_local = 1
        call r8ge_sl ( n, a_lu, pivot, w, job_local )
c
c  Set beta = w' * b.
c
        beta = sum ( w(1:n) * b(1:n) )
c
c  Solve A * b = b.
c
        job_local = 0
        call r8ge_sl ( n, a_lu, pivot, b, job_local )
c
c  Solve A * w = u.
c
        w(1:n) = u(1:n)

        job_local = 0
        call r8ge_sl ( n, a_lu, pivot, w, job_local )
c
c  Set alpha = 1 / ( 1 - v' * w ).
c
        alpha = 1.0D+00 - sum ( v(1:n) * w(1:n) )

      else
c
c  Solve A * w = u.
c
        w(1:n) = u(1:n)

        job_local = 0
        call r8ge_sl ( n, a_lu, pivot, w, job_local )
c
c  Set beta = w' * b.
c
        beta = sum ( w(1:n) * b(1:n) )
c
c  Solve A' * b = b.
c
        job_local = 1
        call r8ge_sl ( n, a_lu, pivot, b, job_local )
c
c  Solve A' * w = v.
c
        w(1:n) = v(1:n)

        job_local = 1
        call r8ge_sl ( n, a_lu, pivot, w, job_local )
c
c  Set alpha = 1 / ( 1 - u' * w ).
c
        alpha = 1.0D+00 - sum ( u(1:n) * w(1:n) )

      end if

      if ( alpha .eq. 0.0D+00 ) then
        ierror = 1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8SM_SL - Fatal error!'
        write ( *, '(a)' ) '  The divisor ALPHA is zero.'
        stop
      end if

      alpha = 1.0D+00 / alpha
c
c  Set b = b + alpha * beta * w.
c
      b(1:n) = b(1:n) + alpha * beta * w(1:n)

      return
      end
      subroutine r8sm_to_r8ge ( m, n, a, u, v, b )

c*********************************************************************72
c
cc R8SM_TO_R8GE copies an R8SM matrix to an R8GE matrix.
c
c  Discussion:
c
c    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
c    which is defined by an M by N matrix A, an M vector U, and
c    an N vector V, by B = A - U * V'
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 January 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns 
c    of the matrix.
c
c    Input, double precision A(M,N), the R8SM matrix.
c
c    Input, double precision U(M), V(N), the R8SM vectors.
c
c    Output, double precision B(M,N), the R8GE matrix.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision b(m,n)
      integer i
      double precision u(m)
      double precision v(n)

      do i = 1, m
        b(i,1:n) = a(i,1:n) - u(i) * v(1:n)
      end do

      return
      end
      subroutine r8sm_vxm ( m, n, a, u, v, x, b )

c*********************************************************************72
c
cc R8SM_VXM multiplies an R8VEC by an R8SM matrix.
c
c  Discussion:
c
c    The R8SM storage format is used for an M by N Sherman Morrison matrix B,
c    which is defined by an M by N matrix A, an M vector U, and
c    an N vector V, by B = A - U * V'
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 January 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns 
c    of the matrix.
c
c    Input, double precision A(M,N), the R8SM matrix.
c
c    Input, double precision U(M), V(N), the R8SM vectors.
c
c    Input, double precision X(M), the vector to be multiplied.
c
c    Output, double precision B(N), the product (A-u*v')' * X.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision b(n)
      double precision u(m)
      double precision v(n)
      double precision x(m)

      b(1:n) = matmul ( transpose ( a(1:m,1:n) ), x(1:m) ) 
     &  - v(1:n) * sum ( u(1:m) * x(1:m) )

      return
      end
      subroutine r8sp_check ( m, n, nz_num, row, col, check )

c*********************************************************************72
c
cc R8SP_CHECK checks that an R8SP matrix data structure is properly sorted.
c
c  Discussion:
c
c    This routine assumes that the data structure has been sorted,
c    so that the entries of ROW are ascending sorted, and that the
c    entries of COL are ascending sorted, within the group of entries
c    that have a common value of ROW.
c
c    The R8SP storage format stores the row, column and value of each nonzero
c    entry of a sparse matrix.
c
c    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP 
c    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns of 
c    the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero elements in
c    the matrix.
c
c    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and 
c    column indices of the nonzero elements.
c
c    Output, logical CHECK, is TRUE if the matrix is properly defined.
c
      implicit none

      integer nz_num

      logical check
      integer col(nz_num)
      integer i
      integer j
      integer k
      integer m
      integer n
      integer row(nz_num)

      check = .true.
c
c  Check 1 .le. ROW(*) .le. M.
c
      do k = 1, nz_num

        if ( row(k) .lt. 1 .or. m .lt. row(k) ) then
          check = .false.
          return
        end if

      end do
c
c  Check 1 .le. COL(*) .le. N.
c
      do k = 1, nz_num

        if ( col(k) .lt. 1 .or. n .lt. col(k) ) then
          check = .false.
          return
        end if

      end do
c
c  Check that ROW(K) .le. ROW(K+1).
c
      do k = 1, nz_num - 1

        if ( row(k+1) .lt. row(k) ) then
          check = .false.
          return
        end if

      end do
c
c  Check that, if ROW(K) .eq. ROW(K+1), that COL(K) .lt. COL(K+1).
c
      do k = 1, nz_num - 1

        if ( row(k) .eq. row(k+1) ) then
          if ( col(k+1) .le. col(k) ) then
            check = .false.
            return
          end if
        end if

      end do

      return
      end
      subroutine r8sp_ij_to_k ( nz_num, row, col, i, j, k )

c*********************************************************************72
c
cc R8SP_IJ_TO_K seeks the compressed index of the (I,J) entry of A.
c
c  Discussion:
c
c    If A(I,J) is nonzero, then its value is stored in location K.
c
c    This routine searches the R8SP storage structure for the index K
c    corresponding to (I,J), returning -1 if no such entry was found.
c
c    This routine assumes that the data structure has been sorted,
c    so that the entries of ROW are ascending sorted, and that the
c    entries of COL are ascending sorted, within the group of entries
c    that have a common value of ROW.
c
c    The R8SP storage format stores the row, column and value of each nonzero
c    entry of a sparse matrix.
c
c    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP 
c    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 July 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NZ_NUM, the number of nonzero elements in
c    the matrix.
c
c    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and 
c    column indices of the nonzero elements.
c
c    Input, integer I, J, the row and column indices of the
c    matrix entry.
c
c    Output, integer K, the R8SP index of the (I,J) entry.
c
      implicit none

      integer nz_num

      integer col(nz_num)
      integer hi
      integer i
      integer j
      integer k
      integer lo
      integer md
      integer row(nz_num)

      lo = 1
      hi = nz_num

      do

        if ( hi .lt. lo ) then
          k = -1
          exit
        end if

        md = ( lo + hi ) / 2

        if ( row(md) .lt. i .or. ( row(md) .eq. i .and. 
     &    col(md) .lt. j ) ) then
          lo = md + 1
        else if ( i .lt. row(md) .or. ( row(md) .eq. i .and. 
     &    j .lt. col(md) ) ) then
          hi = md - 1
        else
          k = md
          exit
        end if

      end do

      return
      end
      subroutine r8sp_indicator ( m, n, nz_num, row, col, a )

c*********************************************************************72
c
cc R8SP_INDICATOR sets up an R8SP indicator matrix.
c
c  Discussion:
c
c    The "indicator matrix" simply has a value like I*10+J at every
c    entry of a dense matrix, or at every entry of a compressed storage
c    matrix for which storage is allocated. 
c
c    The "indicator matrix" simply has a value like I*10+J at every
c    entry of a dense matrix, or at every entry of a compressed storage
c    matrix for which storage is allocated. 
c
c    The R8SP storage format stores the row, column and value of each nonzero
c    entry of a sparse matrix.
c
c    It is possible that a pair of indices (I,J) may occur more than
c    once.  Presumably, in this case, the intent is that the actual value
c    of A(I,J) is the sum of all such entries.  This is not a good thing
c    to do, but I seem to have come across this in MATLAB.
c
c    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP 
c    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 January 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns of 
c    the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero elements in
c    the matrix.
c
c    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and 
c    column indices of the nonzero elements.
c
c    Output, double precision A(NZ_NUM), the nonzero elements of the matrix.
c
      implicit none

      integer nz_num

      double precision a(nz_num)
      integer col(nz_num)
      integer fac
      integer i
      integer i4_log_10
      integer j
      integer k
      integer m
      integer n
      integer row(nz_num)

      fac = 10 ** ( i4_log_10 ( n ) + 1 )

      do k = 1, nz_num

        i = row(k)
        j = col(k)
        a(k) = dble ( fac * i + j )

      end do

      return
      end
      subroutine r8sp_mxv ( m, n, nz_num, row, col, a, x, b )

c*********************************************************************72
c
cc R8SP_MXV multiplies an R8SP matrix by an R8VEC.
c
c  Discussion:
c
c    The R8SP storage format stores the row, column and value of each nonzero
c    entry of a sparse matrix.
c
c    It is possible that a pair of indices (I,J) may occur more than
c    once.  Presumably, in this case, the intent is that the actual value
c    of A(I,J) is the sum of all such entries.  This is not a good thing
c    to do, but I seem to have come across this in MATLAB.
c
c    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP 
c    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 January 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns of 
c    the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero elements in
c    the matrix.
c
c    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and 
c    column indices of the nonzero elements.
c
c    Input, double precision A(NZ_NUM), the nonzero elements of the matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(M), the product vector A*X.
c
      implicit none

      integer m
      integer n
      integer nz_num

      double precision a(nz_num)
      double precision b(m)
      integer col(nz_num)
      integer i
      integer j
      integer k
      integer row(nz_num)
      double precision x(n)

      b(1:m) = 0.0D+00

      do k = 1, nz_num

        i = row(k)
        j = col(k)
        b(i) = b(i) + a(k) * x(j)

      end do

      return
      end
      subroutine r8sp_print ( m, n, nz_num, row, col, a, title )

c*********************************************************************72
c
cc R8SP_PRINT prints an R8SP matrix.
c
c  Discussion:
c
c    This version of R8SP_PRINT has been specifically modified to allow,
c    and correctly handle, the case in which a single matrix location
c    A(I,J) is referenced more than once by the sparse matrix structure.
c    In such cases, the routine prints out the sum of all the values.
c
c    The R8SP storage format stores the row, column and value of each nonzero
c    entry of a sparse matrix.
c
c    It is possible that a pair of indices (I,J) may occur more than
c    once.  Presumably, in this case, the intent is that the actual value
c    of A(I,J) is the sum of all such entries.  This is not a good thing
c    to do, but I seem to have come across this in MATLAB.
c
c    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP 
c    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 January 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns of 
c    the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero elements in 
c    the matrix.
c
c    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and column 
c    indices of the nonzero elements.
c
c    Input, double precision A(NZ_NUM), the nonzero elements of the matrix.
c
c    Input, character ( len = * ) TITLE, a title.
c 
      implicit none

      integer nz_num

      double precision a(nz_num)
      integer col(nz_num)
      integer m
      integer n
      integer row(nz_num)
      character ( len = * ) title

      call r8sp_print_some ( m, n, nz_num, row, col, a, 1, 1, 
     &  m, n, title )

      return
      end
      subroutine r8sp_print_some ( m, n, nz_num, row, col, a, ilo, 
     &  jlo, ihi, jhi, title )

c*********************************************************************72
c
cc R8SP_PRINT_SOME prints some of an R8SP matrix.
c
c  Discussion:
c
c    This version of R8SP_PRINT_SOME has been specifically modified to allow,
c    and correctly handle, the case in which a single matrix location
c    A(I,J) is referenced more than once by the sparse matrix structure.
c    In such cases, the routine prints out the sum of all the values.
c
c    The R8SP storage format stores the row, column and value of each nonzero
c    entry of a sparse matrix.
c
c    It is possible that a pair of indices (I,J) may occur more than
c    once.  Presumably, in this case, the intent is that the actual value
c    of A(I,J) is the sum of all such entries.  This is not a good thing
c    to do, but I seem to have come across this in MATLAB.
c
c    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP 
c    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 September 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns 
c    of the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero elements
c    in the matrix.
c
c    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and column
c    indices of the nonzero elements.
c
c    Input, double precision A(NZ_NUM), the nonzero elements of the matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer, parameter :: incx = 5
      integer nz_num

      double precision a(nz_num)
      double precision aij(incx)
      integer col(nz_num)
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
      integer k
      integer m
      integer n
      integer row(nz_num)
      character ( len = * ) title

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
        write ( *, '(''  Col:  '',5(i7,7x))' ) ( j, j = j2lo, j2hi )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) '  ---'
c
c  Determine the range of the rows in this strip.
c
        i2lo = max ( ilo, 1 )
        i2hi = min ( ihi, m )

        do i = i2lo, i2hi
c
c  Print out (up to) 5 entries in row I, that lie in the current strip.
c
          aij(1:inc) = 0.0D+00
c
c  Is matrix entry K actually the value of A(I,J), with J2LO .le. J .le. J2HI?
c  Because MATLAB seems to allow for multiple (I,J,A) entries, we have
c  to sum up what we find.
c 
          do k = 1, nz_num

            if ( i .eq. row(k) .and. 
     &           j2lo .le. col(k) .and. 
     &           col(k) .le. j2hi ) then 

              j2 = col(k) - j2lo + 1
              aij(j2) = aij(j2) + a(k)

            end if

          end do

          if ( any ( aij(1:inc) .ne. 0.0D+00 ) ) then
            write ( *, '(i5,1x,5g14.6)' ) i, aij(1:inc)
          end if

        end do

      end do

      return
      end
      subroutine r8sp_random ( m, n, nz_num, row, col, seed, a )

c*********************************************************************72
c
cc R8SP_RANDOM sets a random R8SP matrix.
c
c  Discussion:
c
c    The R8SP storage format stores the row, column and value of each nonzero
c    entry of a sparse matrix.
c
c    It is possible that a pair of indices (I,J) may occur more than
c    once.  Presumably, in this case, the intent is that the actual value
c    of A(I,J) is the sum of all such entries.  This is not a good thing
c    to do, but I seem to have come across this in MATLAB.
c
c    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP 
c    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 January 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns of 
c    the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero elements 
c    in the matrix.
c
c    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and column
c    indices of the nonzero elements.
c
c    Input/output, integer SEED, a seed for the random number 
c    generator.
c
c    Output, double precision A(NZ_NUM), the nonzero elements of the matrix.
c
      implicit none

      integer nz_num

      double precision a(nz_num)
      integer col(nz_num)
      integer m
      integer n
      integer row(nz_num)
      integer seed

      call r8vec_uniform_01 ( nz_num, seed, a )

      return
      end
      subroutine r8sp_read ( input_file, m, n, nz_num, row, col, a )

c*********************************************************************72
c
cc R8SP_READ reads an R8SP matrix from a file.
c
c  Discussion:
c
c    This routine needs the value of NZ_NUM, which can be determined
c    by a call to R8SP_READ_SIZE.
c
c    The R8SP storage format stores the row, column and value of each nonzero
c    entry of a sparse matrix.
c
c    It is possible that a pair of indices (I,J) may occur more than
c    once.  Presumably, in this case, the intent is that the actual value
c    of A(I,J) is the sum of all such entries.  This is not a good thing
c    to do, but I seem to have come across this in MATLAB.
c
c    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP 
c    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 August 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character ( len = * ) INPUT_FILE, the name of the file to be read.
c
c    Unused, integer M, N, the number of rows and columns of the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero elements 
c    in the matrix.
c
c    Output, integer ROW(NZ_NUM), COL(NZ_NUM), the row and 
c    column indices of the nonzero elements.
c
c    Output, double precision A(NZ_NUM), the nonzero elements 
c    of the matrix.
c
      implicit none

      integer nz_num

      double precision a(nz_num)
      integer col(nz_num)
      character ( len = * ) input_file
      integer input_unit
      integer ios
      integer k
      integer m
      integer n
      integer row(nz_num)

      call get_unit ( input_unit )

      open ( unit = input_unit, file = input_file, status = 'old', 
     &  iostat = ios )

      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8SP_READ - Fatal error!'
        write ( *, '(a)' ) '  Could not open the input file "' 
     &    // trim ( input_file ) // '".'
        stop
      end if

      do k = 1, nz_num

        read ( input_unit, *, iostat = ios ) row(k), col(k), a(k)

        if ( ios /= 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8SP_READ - Fatal error!'
          write ( *, '(a,i8)' ) '  I/O error while reading record ', k
          stop
        end if

      end do

      close ( unit = input_unit )

      return
      end
      subroutine r8sp_read_size ( input_file, m, n, nz_num )

c*********************************************************************72
c
cc R8SP_READ_SIZE reads the size of an R8SP matrix from a file.
c
c  Discussion:
c
c    The value of NZ_NUM is simply the number of records in the input file.
c
c    The values of M and N are determined as the maximum entry in the row 
c    and column vectors.
c
c    The R8SP storage format stores the row, column and value of each nonzero
c    entry of a sparse matrix.
c
c    It is possible that a pair of indices (I,J) may occur more than
c    once.  Presumably, in this case, the intent is that the actual value
c    of A(I,J) is the sum of all such entries.  This is not a good thing
c    to do, but I seem to have come across this in MATLAB.
c
c    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP 
c    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 August 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character ( len = * ) INPUT_FILE, the name of the file to 
c    be read.
c
c    Output, integer M, N, the number of rows and columns
c    of the matrix.
c
c    Output, integer NZ_NUM, the number of nonzero elements 
c    in the matrix.
c
      implicit none

      double precision a_k
      integer col_k
      character ( len = * ) input_file
      integer input_unit
      integer ios
      integer k
      integer m
      integer n
      integer nz_num
      integer row_k

      call get_unit ( input_unit )

      open ( unit = input_unit, file = input_file, status = 'old', 
     &  iostat = ios )

      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8SP_READ_SIZE - Fatal error!'
        write ( *, '(a)' ) '  Could not open the input file "' 
     &    // trim ( input_file ) // '".'
        stop
      end if

      m = 0
      n = 0
      nz_num = 0

      do

        read ( input_unit, *, iostat = ios ) row_k, col_k, a_k

        if ( ios /= 0 ) then
          exit
        end if

        nz_num = nz_num + 1
        m = max ( m, row_k )
        n = max ( n, col_k )

      end do

      close ( unit = input_unit )

      return
      end
      subroutine r8sp_to_r8ge ( m, n, nz_num, row, col, a, b )

c*********************************************************************72
c
cc R8SP_TO_R8GE converts an R8SP matrix to an R8GE matrix.
c
c  Discussion:
c
c    The R8SP storage format stores the row, column and value of each nonzero
c    entry of a sparse matrix.
c
c    It is possible that a pair of indices (I,J) may occur more than
c    once.  Presumably, in this case, the intent is that the actual value
c    of A(I,J) is the sum of all such entries.  This is not a good thing
c    to do, but I seem to have come across this in MATLAB.
c
c    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP 
c    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 January 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns of 
c    the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero elements
c    in the matrix.
c
c    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and column
c    indices of the nonzero elements.
c
c    Input, double precision A(NZ_NUM), the nonzero elements of the matrix.
c
c    Output, double precision B(M,N), the R8GE matrix.
c
      implicit none

      integer m
      integer n
      integer nz_num

      double precision a(nz_num)
      double precision b(m,n)
      integer col(nz_num)
      integer k
      integer row(nz_num)

      b(1:m,1:n) = 0.0D+00

      do k = 1, nz_num
        b(row(k),col(k)) = a(k)
      end do

      return
      end
      subroutine r8sp_to_r8ncf ( m, n, nz_num, row, col, a, rowcol )

c*********************************************************************72
c
cc R8SP_TO_R8NCF converts an R8SP matrix to an R8NCF matrix.
c
c  Discussion:
c
c    The R8SP and R8NCF formats are essentially identical, except that
c    R8SP keeps separate ROW and COLUMN vectors, while R8NCF uses a single
c    ROWCOL array.  Therefore, the input values NZ_NUM and A used in
c    the R8SP representation can be regarded as part of the output
c    values used for the R8NCF representation.
c
c    The R8SP storage format stores the row, column and value of each nonzero
c    entry of a sparse matrix.
c
c    It is possible that a pair of indices (I,J) may occur more than
c    once.  Presumably, in this case, the intent is that the actual value
c    of A(I,J) is the sum of all such entries.  This is not a good thing
c    to do, but I seem to have come across this in MATLAB.
c
c    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP 
c    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
c
c    The R8NCF storage format stores NZ_NUM, the number of nonzeros, 
c    a real array containing the nonzero values, a 2 by NZ_NUM integer 
c    array storing the row and column of each nonzero entry.
c
c    The R8NCF format is used by NSPCG.  NSPCG requires that the information
c    for the diagonal entries of the matrix must come first.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 August 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Unused, integer M, N, the number of rows and columns of the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero elements 
c    in the matrix.
c
c    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and column
c    indices of the nonzero elements.
c
c    Input, double precision A(NZ_NUM), the nonzero elements of the matrix.
c
c    Output, integer ROWCOL(2,NZ_NUM), the R8NCF row and column
c    index vector.
c
      implicit none

      integer nz_num

      double precision a(nz_num)
      integer col(nz_num)
      integer m
      integer n
      integer row(nz_num)
      integer rowcol(2,nz_num)

      rowcol(1,1:nz_num) = row(1:nz_num)
      rowcol(2,1:nz_num) = col(1:nz_num)

      return
      end
      subroutine r8sp_vxm ( m, n, nz_num, row, col, a, x, b )

c*********************************************************************72
c
cc R8SP_VXM multiplies an R8VEC times an R8SP matrix.
c
c  Discussion:
c
c    The R8SP storage format stores the row, column and value of each nonzero
c    entry of a sparse matrix.
c
c    It is possible that a pair of indices (I,J) may occur more than
c    once.  Presumably, in this case, the intent is that the actual value
c    of A(I,J) is the sum of all such entries.  This is not a good thing
c    to do, but I seem to have come across this in MATLAB.
c
c    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP 
c    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
c
c  Modified:
c
c    21 January 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns of 
c    the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero elements in 
c    the matrix.
c
c    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and column 
c    indices of the nonzero elements.
c
c    Input, double precision A(NZ_NUM), the nonzero elements of the matrix.
c
c    Input, double precision X(M), the vector to be multiplied by A.
c
c    Output, double precision B(N), the product vector A'*X.
c
      implicit none

      integer m
      integer n
      integer nz_num

      double precision a(nz_num)
      double precision b(n)
      integer col(nz_num)
      integer i
      integer j
      integer k
      integer row(nz_num)
      double precision x(m)

      b(1:n) = 0.0D+00

      do k = 1, nz_num

        i = row(k)
        j = col(k)
        b(j) = b(j) + a(k) * x(i)

      end do

      return
      end
      subroutine r8sp_write ( m, n, nz_num, row, col, a, output_file )

c*********************************************************************72
c
cc R8SP_WRITE writes a square R8SP matrix to a file.
c
c  Discussion:
c
c    The R8SP storage format stores the row, column and value of each nonzero
c    entry of a sparse matrix.
c
c    It is possible that a pair of indices (I,J) may occur more than
c    once.  Presumably, in this case, the intent is that the actual value
c    of A(I,J) is the sum of all such entries.  This is not a good thing
c    to do, but I seem to have come across this in MATLAB.
c
c    The R8SP format is used by CSPARSE ("sparse triplet"), DLAP/SLAP 
c    ("nonsymmetric SLAP triad"), by MATLAB, and by SPARSEKIT ("COO" format).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 September 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns of 
c    the matrix.
c
c    Input, integer NZ_NUM, the number of nonzero elements in 
c    the matrix.
c
c    Input, integer ROW(NZ_NUM), COL(NZ_NUM), the row and column
c    indices of the nonzero elements.
c
c    Input, double precision A(NZ_NUM), the nonzero elements 
c    of the matrix.
c
c    Input, character ( len = * ) OUTPUT_FILE, the name of the file to which
c    the information is to be written.
c
      implicit none

      integer nz_num

      double precision a(nz_num)
      integer col(nz_num)
      integer ios
      integer k
      integer m
      integer n
      character ( len = * ) output_file
      integer output_unit
      integer row(nz_num)

      call get_unit ( output_unit )

      open ( unit = output_unit, file = output_file, status = 'replace', 
     &  iostat = ios )

      if ( ios /= 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8SP_WRITE - Fatal error!'
        write ( *, '(a)' ) '  Could not open the output file "' 
     &    // trim ( output_file ) // '".'
        stop
      end if

      do k = 1, nz_num
        write ( output_unit, '(2x,i8,2x,i8,2x,g16.8)' ) 
     &    row(k), col(k), a(k)
      end do

      close ( unit = output_unit )

      return
      end
      subroutine r8sr_indicator ( n, nz, row, col, diag, off )

c*********************************************************************72
c
cc R8SR_INDICATOR sets up an R8SR indicator matrix.
c
c  Discussion:
c
c    The "indicator matrix" simply has a value like I*10+J at every
c    entry of a dense matrix, or at every entry of a compressed storage
c    matrix for which storage is allocated. 
c
c    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
c    The off-diagonal entries of row I are stored in entries ROW(I)
c    through ROW(I+1)-1 of OFF.  COL(J) records the column index
c    of the entry in A(J).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 January 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, integer NZ, the number of offdiagonal nonzero elements
c    in the matrix.
c
c    Input, integer ROW(N+1).  The nonzero offdiagonal elements 
c    of row I of A are contained in A(ROW(I)) through A(ROW(I+1)-1).
c
c    Input, integer COL(NZ), contains the column index of the 
c    element in the corresponding position in A.
c
c    Output, double precision DIAG(N), the diagonal elements of A.
c
c    Output, double precision OFF(NZ), the off-diagonal elements of A.
c
      implicit none

      integer n
      integer nz

      integer col(nz)
      double precision diag(n)
      integer fac
      integer i
      integer i4_log_10
      integer j
      integer k
      double precision off(nz)
      integer row(n+1)

      fac = 10 ** ( i4_log_10 ( n ) + 1 )

      do i = 1, n

        j = i
        diag(i) = dble ( fac * i + j )

        do k = row(i), row(i+1)-1
          j = col(k)
          off(k) = dble ( fac * i + j )
        end do

      end do

      return
      end
      subroutine r8sr_mxv ( n, nz, row, col, diag, off, x, b )

c*********************************************************************72
c
cc R8SR_MXV multiplies an R8SR matrix by an R8VEC.
c
c  Discussion:
c
c    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
c    The off-diagonal entries of row I are stored in entries ROW(I)
c    through ROW(I+1)-1 of OFF.  COL(J) records the column index
c    of the entry in A(J).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 October 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, integer NZ, the number of offdiagonal nonzero 
c    elements in the matrix.
c
c    Input, integer ROW(N+1).  The nonzero offdiagonal elements
c    of row I of A are contained in A(ROW(I)) through A(ROW(I+1)-1).
c
c    Input, integer COL(NZ), contains the column index of the 
c    element in the corresponding position in A.
c
c    Input, double precision DIAG(N), the diagonal elements of the matrix.
c
c    Input, double precision OFF(NZ), the off-diagonal elements of the matrix.
c
c    Input, double precision X(N), the vector to be multiplied by the matrix.
c
c    Output, double precision B(N), the product A * X.
c
      implicit none

      integer n
      integer nz

      double precision b(n)
      integer col(nz)
      double precision diag(n)
      integer i
      integer j
      integer k
      double precision off(nz)
      integer row(n+1)
      double precision x(n)

      b(1:n) = diag(1:n) * x(1:n)

      do i = 1, n
        do k = row(i), row(i+1)-1
          j = col(k)
          b(i) = b(i) + off(k) * x(j)
        end do
      end do

      return
      end
      subroutine r8sr_print ( n, nz, row, col, diag, off, title )

c*********************************************************************72
c
cc R8SR_PRINT prints an R8SR matrix.
c
c  Discussion:
c
c    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
c    The off-diagonal entries of row I are stored in entries ROW(I)
c    through ROW(I+1)-1 of OFF.  COL(J) records the column index
c    of the entry in A(J).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 November 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, integer NZ, the number of offdiagonal nonzero elements
c    in A.
c
c    Input, integer ROW(N+1).  The nonzero offdiagonal elements
c    of row I of A are contained in A(ROW(I)) through A(ROW(I+1)-1).
c
c    Input, integer COL(NZ), contains the column index of
c    the element in the corresponding position in A.
c
c    Input, double precision DIAG(N), the diagonal elements of A.
c
c    Input, double precision OFF(NZ), the off-diagonal elements of A.
c
c    Input, character ( len = * ) TITLE, a title.
c 
      implicit none

      integer n
      integer nz

      integer col(nz)
      double precision diag(n)
      double precision off(nz)
      integer row(n+1)
      character ( len = * ) title

      call r8sr_print_some ( n, nz, row, col, diag, off, 1, 1, 
     &  n, n, title )

      return
      end
      subroutine r8sr_print_some ( n, nz, row, col, diag, off, ilo, 
     &  jlo, ihi, jhi, title )

c*********************************************************************72
c
cc R8SR_PRINT_SOME prints some of an R8SR matrix.
c
c  Discussion:
c
c    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
c    The off-diagonal entries of row I are stored in entries ROW(I)
c    through ROW(I+1)-1 of OFF.  COL(J) records the column index
c    of the entry in A(J).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 November 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, integer NZ, the number of offdiagonal nonzero elements
c    in A.
c
c    Input, integer ROW(N+1).  The nonzero offdiagonal elements 
c    of row I of A are contained in A(ROW(I)) through A(ROW(I+1)-1).
c
c    Input, integer COL(NZ), contains the column index of the 
c    element in the corresponding position in A.
c
c    Input, double precision DIAG(N), the diagonal elements of A.
c
c    Input, double precision OFF(NZ), the off-diagonal elements of A.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer, parameter :: incx = 5
      integer n
      integer nz

      double precision aij
      integer col(nz)
      character ( len = 14 ) ctemp(incx)
      logical r8_is_int
      double precision diag(n)
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
      integer k
      double precision off(nz)
      integer row(n+1)
      character ( len = * ) title

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

        write ( *, '(''  Col:  '',5a14)' ) ( ctemp(j2), j2 = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) '  ---'
c
c  Determine the range of the rows in this strip.
c
        i2lo = max ( ilo, 1 )
        i2hi = min ( ihi, n )

        do i = i2lo, i2hi
c
c  Print out (up to) 5 entries in row I, that lie in the current strip.
c
c  1) Assume everything is zero.
c
          aij = 0.0D+00
          do j2 = 1, inc
            write ( ctemp(j2), '(f8.0,6x)' ) aij
          end do
c
c  2) Insert the diagonal entry, if appropriate.
c
          if ( j2lo .le. i .and. i .le. j2hi ) then
            j2 = i - j2lo + 1
            aij = diag(i)
            if ( r8_is_int ( aij ) ) then
              write ( ctemp(j2), '(f8.0,6x)' ) aij
            else
              write ( ctemp(j2), '(g14.6)' ) aij
            end if
          end if
c
c  3) Now examine all the offdiagonal entries.
c
          do k = row(i), row(i+1)-1
            if ( j2lo .le. col(k) .and. col(k) .le. j2hi ) then 
              j2 = col(k) - j2lo + 1
              aij = off(k)
              if ( r8_is_int ( aij ) ) then
                write ( ctemp(j2), '(f8.0,6x)' ) aij
              else
                write ( ctemp(j2), '(g14.6)' ) aij
              end if
            end if
          end do

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine r8sr_random ( n, nz, row, col, diag, off, seed )

c*********************************************************************72
c
cc R8SR_RANDOM randomizes an R8SR matrix.
c
c  Discussion:
c
c    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
c    The off-diagonal entries of row I are stored in entries ROW(I)
c    through ROW(I+1)-1 of OFF.  COL(J) records the column index
c    of the entry in A(J).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 October 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, integer NZ, the number of offdiagonal nonzero elements 
c    in A.
c
c    Input, integer ROW(N+1).  The nonzero offdiagonal elements 
c    of row I of A are contained in A(ROW(I)) through A(ROW(I+1)-1).
c
c    Input, integer COL(NZ), contains the column index of the 
c    element in the corresponding position in A.
c
c    Output, double precision DIAG(N), the diagonal elements of A.
c
c    Output, double precision OFF(NZ), the off-diagonal elements of A.
c
c    Input/output, integer SEED, a seed for the random number 
c    generator.
c
      implicit none

      integer n
      integer nz

      integer col(nz)
      double precision r8_uniform_01
      double precision diag(n)
      integer i
      integer j
      double precision off(nz)
      integer row(n+1)
      integer seed

      do i = 1, n
        diag(i) = r8_uniform_01 ( seed )
        do j = row(i), row(i+1)-1
          off(j) = r8_uniform_01 ( seed )
        end do
      end do

      return
      end
      subroutine r8sr_to_r8ge ( n, nz, row, col, diag, off, b )

c*********************************************************************72
c
cc R8SR_TO_R8GE converts an R8SR matrix to an R8GE matrix.
c
c  Discussion:
c
c    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
c    The off-diagonal entries of row I are stored in entries ROW(I)
c    through ROW(I+1)-1 of OFF.  COL(J) records the column index
c    of the entry in A(J).
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 April 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, integer NZ, the number of offdiagonal nonzero 
c    elements in A.
c
c    Input, integer ROW(N+1).  The nonzero offdiagonal elements 
c    of row I of A are contained in A(ROW(I)) through A(ROW(I+1)-1).
c
c    Input, integer COL(NZ), contains the column index of the 
c    element in the corresponding position in A.
c
c    Input, double precision DIAG(N), the diagonal elements of A.
c
c    Input, double precision OFF(NZ), the off-diagonal elements of A.
c
c    Output, double precision B(N,N), the R8GE matrix.
c
      implicit none

      integer n
      integer nz

      double precision b(n,n)
      integer col(nz)
      double precision diag(n)
      integer i
      integer j
      double precision off(nz)
      integer row(n+1)

      if ( n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8SR_TO_R8GE - Fatal error!'
        write ( *, '(a,i8)' ) 
     &    '  N is less than or equal to zero, N = ', n
        stop
      end if

      b(1:n,1:n) = 0.0D+00

      do i = 1, n
        b(i,i) = diag(i)
      end do

      do i = 1, n
        do j = row(i), row(i+1)-1
          b(i,col(j)) = off(j)
        end do
      end do

      return
      end
      subroutine r8sr_vxm ( n, nz, row, col, diag, off, x, b )

c*********************************************************************72
c
cc R8SR_VXM multiplies an R8VEC times an R8SR matrix.
c
c  Discussion:
c
c    The R8SR storage format stores the diagonal of a sparse matrix in DIAG.
c    The off-diagonal entries of row I are stored in entries ROW(I)
c    through ROW(I+1)-1 of OFF.  COL(J) records the column index
c    of the entry in A(J).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 October 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, integer NZ, the number of offdiagonal nonzero 
c    elements in A.
c
c    Input, integer ROW(N+1).  The nonzero offdiagonal elements 
c    of row I of A are contained in A(ROW(I)) through A(ROW(I+1)-1).
c
c    Input, integer COL(NZ), contains the column index of the 
c    element in the corresponding position in A.
c
c    Input, double precision DIAG(N), the diagonal elements of A.
c
c    Input, double precision OFF(NZ), the off-diagonal elements of A.
c
c    Input, double precision X(N), the vector to be multiplies by A.
c
c    Output, double precision B(N), the product A' * X.
c
      implicit none

      integer n
      integer nz

      double precision b(n)
      integer col(nz)
      double precision diag(n)
      integer i
      integer j
      integer k
      double precision off(nz)
      integer row(n+1)
      double precision x(n)

      b(1:n) = diag(1:n) * x(1:n)

      do i = 1, n
        do k = row(i), row(i+1)-1
          j = col(k)
          b(j) = b(j) + off(k) * x(i)
        end do
      end do

      return
      end
      subroutine r8ss_error ( n, na, diag, ierror )

c*********************************************************************72
c
cc R8SS_ERROR checks dimensions for an R8SS matrix.
c
c  Discussion:
c
c    The R8SS storage format is used for real symmetric skyline matrices.
c    This storage is appropriate when the nonzero entries of the
c    matrix are generally close to the diagonal, but the number
c    of nonzeroes above each diagonal varies in an irregular fashion.
c
c    In this case, the strategy is essentially to assign column J
c    its own bandwidth, and store the strips of nonzeros one after
c    another.  Note that what's important is the location of the
c    furthest nonzero from the diagonal.  A slot will be set up for
c    every entry between that and the diagonal, whether or not
c    those entries are zero.
c
c    A skyline matrix can be Gauss-eliminated without disrupting
c    the storage scheme, as long as no pivoting is required.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 October 1998
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
c    Input, integer NA, the dimension of the array A.
c    NA must be at least N.
c
c    Input, integer DIAG(N), the indices in A of the N diagonal 
c    elements.
c
c    Output, integer IERROR, error indicator.
c    0, no error.
c    1, N is less than 1.
c    2, NA is less than N.
c    3, DIAG(1) is not 1.
c    4, the elements of DIAG are not strictly increasing.
c    5, DIAG(N) is greater than NA.
c
      implicit none

      integer n

      integer diag(n)
      integer i
      integer ierror
      integer na

      ierror = 0

      if ( n .lt. 1 ) then
        ierror = 1
        return
      end if

      if ( na .lt. n ) then
        ierror = 2
        return
      end if

      if ( diag(1) .ne. 1 ) then
        ierror = 3
        return
      end if

      do i = 1, n-1
        if ( diag(i+1) .le. diag(i) ) then
          ierror = 4
          return
        end if
      end do

      if ( na .lt. diag(n) ) then
        ierror = 5
        return
      end if

      return
      end
      subroutine r8ss_indicator ( n, na, diag, a )

c*********************************************************************72
c
cc R8SS_INDICATOR sets up an R8SS indicator matrix.
c
c  Discussion:
c
c    The "indicator matrix" simply has a value like I*10+J at every
c    entry of a dense matrix, or at every entry of a compressed storage
c    matrix for which storage is allocated. 
c
c    The R8SS storage format is used for real symmetric skyline matrices.
c    This storage is appropriate when the nonzero entries of the
c    matrix are generally close to the diagonal, but the number
c    of nonzeroes above each diagonal varies in an irregular fashion.
c
c    In this case, the strategy is essentially to assign column J
c    its own bandwidth, and store the strips of nonzeros one after
c    another.   Note that what's important is the location of the
c    furthest nonzero from the diagonal.  A slot will be set up for
c    every entry between that and the diagonal, whether or not
c    those entries are zero.
c
c    A skyline matrix can be Gauss-eliminated without disrupting
c    the storage scheme, as long as no pivoting is required.
c
c    The user must set aside ( N * ( N + 1 ) ) / 2 entries for the array,
c    although the actual storage needed will generally be about half of
c    that.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 January 2004
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
c    Output, integer NA, the dimension of the array A, which for
c    this special case will be the maximum, ( N * ( N + 1 ) ) / 2
c
c    Output, integer DIAG(N), the indices in A of the N diagonal
c    elements.
c
c    Output, double precision A((N*(N+1))/2), the R8SS matrix.
c
      implicit none

      integer n

      double precision a((n*(n+1))/2)
      integer fac
      integer diag(n)
      integer i
      integer i4_log_10
      integer j
      integer na

      fac = 10 ** ( i4_log_10 ( n ) + 1 )

      na = 0

      do j = 1, n

        do i = 1, j
          na = na + 1
          a(na) = dble ( fac * i + j )
        end do

        diag(j) = na

      end do

      return
      end
      subroutine r8ss_mxv ( n, na, diag, a, x, b )

c*********************************************************************72
c
cc R8SS_MXV multiplies an R8SS matrix by an R8VEC.
c
c  Discussion:
c
c    The R8SS storage format is used for real symmetric skyline matrices.
c    This storage is appropriate when the nonzero entries of the
c    matrix are generally close to the diagonal, but the number
c    of nonzeroes above each diagonal varies in an irregular fashion.
c
c    In this case, the strategy is essentially to assign column J
c    its own bandwidth, and store the strips of nonzeros one after
c    another.  Note that what's important is the location of the
c    furthest nonzero from the diagonal.  A slot will be set up for
c    every entry between that and the diagonal, whether or not
c    those entries are zero.
c
c    A skyline matrix can be Gauss-eliminated without disrupting
c    the storage scheme, as long as no pivoting is required.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 October 1998
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
c    Input, integer NA, the dimension of the array A.
c    NA must be at least N.
c
c    Input, integer DIAG(N), the indices in A of the N diagonal
c    elements.
c
c    Input, double precision A(NA), the R8SS matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(N), the product A*x.
c
      implicit none

      integer n
      integer na

      double precision a(na)
      double precision b(n)
      integer diag(n)
      integer diagold
      integer i
      integer ilo
      integer j
      integer k
      double precision x(n)

      b(1:n) = 0.0D+00

      diagold = 0
      k = 0

      do j = 1, n

        ilo = j + 1 + diagold - diag(j)

        do i = ilo, j-1
          k = k + 1
          b(i) = b(i) + a(k) * x(j)
          b(j) = b(j) + a(k) * x(i)
        end do

        k = k + 1
        b(j) = b(j) + a(k) * x(j)

        diagold = diag(j)

      end do

      return
      end
      subroutine r8ss_print ( n, na, diag, a, title )

c*********************************************************************72
c
cc R8SS_PRINT prints an R8SS matrix.
c
c  Discussion:
c
c    The R8SS storage format is used for real symmetric skyline matrices.
c    This storage is appropriate when the nonzero entries of the
c    matrix are generally close to the diagonal, but the number
c    of nonzeroes above each diagonal varies in an irregular fashion.
c
c    In this case, the strategy is essentially to assign column J
c    its own bandwidth, and store the strips of nonzeros one after
c    another.  Note that what's important is the location of the
c    furthest nonzero from the diagonal.  A slot will be set up for
c    every entry between that and the diagonal, whether or not
c    those entries are zero.
c
c    A skyline matrix can be Gauss-eliminated without disrupting
c    the storage scheme, as long as no pivoting is required.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 May 2000
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
c    Input, integer NA, the dimension of the array A.
c
c    Input, integer DIAG(N), the indices in A of the N 
c    diagonal elements.
c
c    Input, double precision A(NA), the R8SS matrix.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer na
      integer n

      double precision a(na)
      integer diag(n)
      character ( len = * ) title

      call r8ss_print_some ( n, na, diag, a, 1, 1, n, n, title )

      return
      end
      subroutine r8ss_print_some ( n, na, diag, a, ilo, jlo, ihi, 
     &  jhi, title )

c*********************************************************************72
c
cc R8SS_PRINT_SOME prints some of an R8SS matrix.
c
c  Discussion:
c
c    The R8SS storage format is used for real symmetric skyline matrices.
c    This storage is appropriate when the nonzero entries of the
c    matrix are generally close to the diagonal, but the number
c    of nonzeroes above each diagonal varies in an irregular fashion.
c
c    In this case, the strategy is essentially to assign column J
c    its own bandwidth, and store the strips of nonzeros one after
c    another.  Note that what's important is the location of the
c    furthest nonzero from the diagonal.  A slot will be set up for
c    every entry between that and the diagonal, whether or not
c    those entries are zero.
c
c    A skyline matrix can be Gauss-eliminated without disrupting
c    the storage scheme, as long as no pivoting is required.
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
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input, integer NA, the dimension of the array A.
c
c    Input, integer DIAG(N), the indices in A of the N diagonal
c    elements.
c
c    Input, double precision A(NA), the R8SS matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer, parameter :: incx = 5
      integer na
      integer n

      double precision a(na)
      double precision aij
      character ( len = 14 ) ctemp(incx)
      logical r8_is_int
      integer diag(n)
      integer i
      integer i2hi
      integer i2lo
      integer ihi
      integer ij
      integer ijm1
      integer ilo
      integer inc
      integer j
      integer j2
      integer j2hi
      integer j2lo
      integer jhi
      integer jlo
      character ( len = * ) title

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
        i2hi = min ( ihi, n )

        do i = i2lo, i2hi
c
c  Print out (up to) 5 entries in row I, that lie in the current strip.
c
          do j2 = 1, inc

            j = j2lo - 1 + j2

            aij = 0.0D+00

            if ( j .lt. i ) then
              if ( i .eq. 1 ) then
                ijm1 = 0
              else
                ijm1 = diag(i-1)
              end if
              ij = diag(i)
              if ( ijm1 .lt. ij+j-i ) then
                aij = a(ij+j-i)
              end if
            else if ( j .eq. i ) then
              ij = diag(j)
              aij = a(ij)
            else if ( i .lt. j ) then
              if ( j .eq. 1 ) then
                ijm1 = 0
              else
                ijm1 = diag(j-1)
              end if
              ij = diag(j)
              if ( ijm1 .lt. ij+i-j ) then
                aij = a(ij+i-j)
              end if
            end if

            if ( r8_is_int ( aij ) ) then
              write ( ctemp(j2), '(f8.0,6x)' ) aij
            else
              write ( ctemp(j2), '(g14.6)' ) aij
            end if

          end do

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine r8ss_random ( n, na, diag, a, seed )

c*********************************************************************72
c
cc R8SS_RANDOM randomizes an R8SS matrix.
c
c  Discussion:
c
c    The R8SS storage format is used for real symmetric skyline matrices.
c    This storage is appropriate when the nonzero entries of the
c    matrix are generally close to the diagonal, but the number
c    of nonzeroes above each diagonal varies in an irregular fashion.
c
c    In this case, the strategy is essentially to assign column J
c    its own bandwidth, and store the strips of nonzeros one after
c    another.  Note that what's important is the location of the
c    furthest nonzero from the diagonal.  A slot will be set up for
c    every entry between that and the diagonal, whether or not
c    those entries are zero.
c
c    A skyline matrix can be Gauss-eliminated without disrupting
c    the storage scheme, as long as no pivoting is required.
c
c    The user must set aside ( N * ( N + 1 ) ) / 2 entries for the array,
c    although the actual storage needed will generally be about half of
c    that.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 October 2004
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
c    Output, integer NA, the dimension of the array A.
c    NA will be at least N and no greater than ( N * ( N + 1 ) ) / 2.
c
c    Output, integer DIAG(N), the indices in A of the N diagonal
c    elements.
c
c    Output, double precision A((N*(N+1))/2), the R8SS matrix.
c
c    Input/output, integer SEED, a seed for the random number
c    generator.
c
      implicit none

      integer n
      integer na

      double precision a((n*(n+1))/2)
      double precision r8_uniform_01
      integer diag(n)
      integer diagold
      integer i
      integer i4_uniform
      integer ilo
      integer j
      integer k
      integer seed
c
c  Set the values of DIAG.
c
      diag(1) = 1
      na = 1
      do i = 2, n
        k = i4_uniform ( 1, i, seed )
        diag(i) = diag(i-1) + k
        na = na + k
      end do
c
c  Now set the values of A.
c
      diagold = 0
      k = 0

      do j = 1, n

        ilo = j + 1 + diagold - diag(j)

        do i = ilo, j
          k = k + 1
          a(k) = r8_uniform_01 ( seed )
        end do

        diagold = diag(j)

      end do

      return
      end
      subroutine r8ss_to_r8ge ( n, na, diag, a, b  )

c*********************************************************************72
c
cc R8SS_TO_R8GE copies an R8SS matrix to an R8GE matrix.
c
c  Discussion:
c
c    The R8SS storage format is used for real symmetric skyline matrices.
c    This storage is appropriate when the nonzero entries of the
c    matrix are generally close to the diagonal, but the number
c    of nonzeroes above each diagonal varies in an irregular fashion.
c
c    In this case, the strategy is essentially to assign column J
c    its own bandwidth, and store the strips of nonzeros one after
c    another.  Note that what's important is the location of the
c    furthest nonzero from the diagonal.  A slot will be set up for
c    every entry between that and the diagonal, whether or not
c    those entries are zero.
c
c    A skyline matrix can be Gauss-eliminated without disrupting
c    the storage scheme, as long as no pivoting is required.
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c  Example:
c
c    11   0  13  0 15
c     0  22  23  0  0
c    31  32  33 34  0
c     0   0  43 44  0
c    51   0   0  0 55
c
c    A = ( 11 | 22 | 13, 23, 33 | 34, 44 | 15, 0, 0, 0, 55 )
c    NA = 12
c    DIAG = ( 1, 2, 5, 7, 12 )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 October 1998
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
c    Input, integer NA, the dimension of the array A.
c    NA must be at least N.
c
c    Input, integer DIAG(N), the indices in A of the N diagonal
c    elements.
c
c    Input, double precision A(NA), the R8SS matrix.
c
c    Output, double precision B(N,N), the R8GE matrix.
c
      implicit none

      integer n
      integer na

      double precision a(na)
      double precision b(n,n)
      integer diag(n)
      integer diagold
      integer i
      integer ilo
      integer j
      integer k

      diagold = 0
      k = 0

      do j = 1, n

        ilo = j + 1 + diagold - diag(j)

        b(1:ilo-1,j) = 0.0D+00
        b(j,1:ilo-1) = 0.0D+00

        do i = ilo, j-1
          k = k + 1
          b(i,j) = a(k)
          b(j,i) = a(k)
        end do

        k = k + 1
        b(j,j) = a(k)

        diagold = diag(j)

      end do

      return
      end
      subroutine r8sto_indicator ( n, a )

c*********************************************************************72
c
cc R8STO_INDICATOR sets up an R8STO indicator matrix.
c
c  Discussion:
c
c    The "indicator matrix" simply has a value like I*10+J at every
c    entry of a dense matrix, or at every entry of a compressed storage
c    matrix for which storage is allocated. 
c
c    The R8STO storage format is used for a symmetric Toeplitz matrix.
c    It stores the N elements of the first row.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 January 2004
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
c    Output, double precision A(N), the R8STO matrix.
c
      implicit none

      integer n

      double precision a(n)
      integer fac
      integer i
      integer i4_log_10
      integer j
      integer k

      fac = 10 ** ( i4_log_10 ( n ) + 1 )

      i = 1
      k = 0
      do j = 1, n
        k = k + 1
        a(k) = dble ( fac * i + j )
      end do
      
      return
      end
      subroutine r8sto_inverse ( n, a, b )

c*********************************************************************72
c
cc R8STO_INVERSE computes the inverse of an R8STO matrix.
c
c  Discussion:
c
c    The R8STO storage format is used for a symmetric Toeplitz matrix.
c    It stores the N elements of the first row.
c
c    For this routine, the matrix is also required to be positive definite.
c
c    The original implementation of the algorithm assumed that the
c    diagonal element was 1.  The algorithm has been modified so that
c    this is no longer necessary.
c
c    The inverse matrix is NOT guaranteed to be a Toeplitz matrix.  
c    It is guaranteed to be symmetric and persymmetric.
c    The inverse matrix is returned in general storage, that is,
c    as an "SGE" matrix.
c
c  Example:
c
c    To compute the inverse of
c
c     1.0 0.5 0.2
c     0.5 1.0 0.5
c     0.2 0.5 1.0
c
c    we input:
c
c      N = 3
c      A = (/ 1.0, 0.5, 0.2 /)
c
c    with output:
c
c      B(1:3,1:3) = (/ 75, -40,   5,
c                     -40,  96, -40,
c                       5, -40,  75 /) / 56
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 May 2002
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Gene Golub, Charles Van Loan,
c    Section 4.7.3, "Computing the Inverse",
c    Matrix Computations,
c    Third Edition,
c    Johns Hopkins, 1996.
c
c  Parameters:
c
c    Input, integer N, the order of the system.
c
c    Input, double precision A(N), the R8STO matrix.
c
c    Output, double precision B(N,N), the inverse of the matrix, in R8GE format.
c
      implicit none

      integer n

      double precision a(n)
      double precision a2(n-1)
      double precision b(n,n)
      integer i
      integer j
      double precision v(n)

      a2(1:n-1) = a(2:n) / a(1)

      call r8sto_yw_sl ( n - 1, a2, v )
c
c  Compute the N-th entry of V.
c
      v(n) = 1.0D+00 / ( 1.0D+00 + sum ( a2(1:n-1) * v(1:n-1) ) )
c
c  Reverse and scale entries 1 through N-1.
c
      v(1:n-1) = v(n-1:1:-1)

      v(1:n-1) = v(n) * v(1:n-1)
c
c  Set the boundaries of B.
c
      b(1,1:n) = v(n:1:-1)
      b(n,1:n) = v(1:n)
      b(2:n-1,1) = v(n-1:2:-1)
      b(2:n-1,n) = v(2:n-1)
c
c  Fill the interior.
c
      do i = 2, 1+((n-1)/2)
        do j = i, n-i+1
          b(i,j) = b(i-1,j-1) 
     &      + ( v(n+1-j) * v(n+1-i) - v(i-1) * v(j-1) ) / v(n)
          b(j,i) = b(i,j)
          b(n+1-i,n+1-j) = b(i,j)
          b(n+1-j,n+1-i) = b(i,j)
        end do
      end do
c
c  Scale B.
c
      b(1:n,1:n) = b(1:n,1:n) / a(1)

      return
      end
      subroutine r8sto_mxv ( n, a, x, b )

c*********************************************************************72
c
cc R8STO_MXV multiplies an R8STO matrix by an R8VEC.
c
c  Discussion:
c
c    The R8STO storage format is used for a symmetric Toeplitz matrix.
c    It stores the N elements of the first row.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision A(N), the R8STO matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(N), the product A * x.
c
      implicit none

      integer n

      double precision a(n)
      double precision b(n)
      integer i
      double precision x(n)

      do i = 1, n
        b(i) = sum ( a(i:2:-1) * x(1:i-1) ) 
     &    + sum ( a(1:n+1-i) * x(i:n) )
      end do

      return
      end
      subroutine r8sto_print ( n, a, title )

c*********************************************************************72
c
cc R8STO_PRINT prints an R8STO matrix.
c
c  Discussion:
c
c    The R8STO storage format is used for a symmetric Toeplitz matrix.
c    It stores the N elements of the first row.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 May 2000
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
c    Input, double precision A(N), the R8STO matrix.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer n

      double precision a(n)
      character ( len = * ) title

      call r8sto_print_some ( n, a, 1, 1, n, n, title )

      return
      end
      subroutine r8sto_print_some ( n, a, ilo, jlo, ihi, jhi, title )

c*********************************************************************72
c
cc R8STO_PRINT_SOME prints some of an R8STO matrix.
c
c  Discussion:
c
c    The R8STO storage format is used for a symmetric Toeplitz matrix.
c    It stores the N elements of the first row.
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
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input, double precision A(N), the R8STO matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer, parameter :: incx = 5
      integer n

      double precision a(n)
      double precision aij
      character ( len = 14 ) ctemp(incx)
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
      character ( len = * ) title

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
        i2hi = min ( ihi, n )

        do i = i2lo, i2hi
c
c  Print out (up to) 5 entries in row I, that lie in the current strip.
c
          do j2 = 1, inc

            j = j2lo - 1 + j2

            if ( i .le. j ) then
              aij = a(1+j-i)
            else
              aij = a(1+i-j)
            end if

            if ( r8_is_int ( aij ) ) then
              write ( ctemp(j2), '(f8.0,6x)' ) aij
            else
              write ( ctemp(j2), '(g14.6)' ) aij
            end if

          end do

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine r8sto_random ( n, seed, a )

c*********************************************************************72
c
cc R8STO_RANDOM randomizes an R8STO matrix.
c
c  Discussion:
c
c    The R8STO storage format is used for a symmetric Toeplitz matrix.
c    It stores the N elements of the first row.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 May 2003
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
c    Input/output, integer SEED, a seed for the random number
c    generator.
c
c    Output, double precision A(N), the R8STO matrix.
c
      implicit none

      integer n

      double precision a(n)
      integer seed

      call r8vec_uniform_01 ( n, seed, a )

      return
      end
      subroutine r8sto_sl ( n, a, b, x )

c*********************************************************************72
c
cc R8STO_SL solves an R8STO system.
c
c  Discussion:
c
c    The R8STO storage format is used for a symmetric Toeplitz matrix.
c    It stores the N elements of the first row.
c
c    For this routine, the matrix is also required to be positive definite.
c
c    This implementation of the algorithm assumes that the diagonal element
c    is 1.
c
c    The real symmetric Toeplitz matrix can be described by N numbers, which,
c    for convenience, we will label A(0:N-1).
c
c    Note that there is a typographical error in the presentation
c    of this algorithm in the reference, and another in the presentation
c    of a sample problem.  Both involve sign errors.  A minor error
c    makes the algorithm incorrect for the case N = 1.
c
c  Example:
c
c    To solve
c
c     1.0 0.5 0.2    x1    4.0
c     0.5 1.0 0.5 *  x2 = -1.0
c     0.2 0.5 1.0    x3    3.0
c
c    we input:
c
c      N = 3
c      A(0:N-1) = (/ 1.0, 0.5, 0.2 /)
c      B(1:3) = (/ 4.0, -1.0, 3.0 /)
c
c    with output:
c
c      X(1:3) = (/ 355, -376, 285 /) / 56
c             = (/ 6.339, -6.714, 5.089 /)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 September 2003
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Gene Golub, Charles Van Loan,
c    Section 4.7.3, "The General Right Hand Side Problem",
c    Matrix Computations,
c    Third Edition,
c    Johns Hopkins, 1996.
c
c  Parameters:
c
c    Input, integer N, the order of the system.
c
c    Input, double precision A(N), the R8STO matrix.
c
c    Input, double precision B(N), the right hand side of the linear system.
c
c    Output, double precision X(N), the solution of the linear system.
c
      implicit none

      integer n

      double precision a(n)
      double precision b(n)
      double precision beta
      integer k
      double precision x(n)
      double precision y(n)

      k = 0

      beta = 1.0D+00
      x(k+1) = b(k+1) / beta

      if ( k .lt. n-1 ) then
        y(k+1) = -a(k+2) / beta
      end if

      do k = 1, n-1

        beta = ( 1.0D+00 - y(k) * y(k) ) * beta

        x(k+1) = ( b(k+1) - sum ( a(2:k+1) * x(k:1:-1) ) ) / beta

        x(1:k) = x(1:k) + x(k+1) * y(k:1:-1)

        if ( k .lt. n - 1 ) then
          y(k+1) = ( -a(k+2) - sum ( a(2:k+1) * y(k:1:-1) ) ) / beta
          y(1:k) = y(1:k) + y(k+1) * y(k:1:-1)
        end if

      end do

      return
      end
      subroutine r8sto_to_r8ge ( n, a, b )

c*********************************************************************72
c
cc R8STO_TO_R8GE copies an R8STO matrix to an R8GE matrix.
c
c  Discussion:
c
c    The R8STO storage format is used for a symmetric Toeplitz matrix.
c    It stores the N elements of the first row.
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision A(N), the R8STO matrix.
c
c    Output, double precision B(N,N), the R8GE matrix.
c
      implicit none

      integer n

      double precision a(n)
      double precision b(n,n)
      integer i

      do i = 1, n
        b(i,1:i-1) = a(i:2:-1)
        b(i,i:n) = a(1:n-i+1)
      end do

      return
      end
      subroutine r8sto_yw_sl ( n, b, x )

c*********************************************************************72
c
cc R8STO_YW_SL solves the Yule-Walker equations for an R8STO matrix.
c
c  Discussion:
c
c    The R8STO storage format is used for a symmetric Toeplitz matrix.
c    It stores the N elements of the first row.
c
c    The matrix is also required to be positive definite.
c
c    This implementation of the algorithm assumes that the diagonal element
c    is 1.
c
c    The real symmetric Toeplitz matrix can be described by N numbers, which,
c    for convenience, we will label B(0:N-1).  We assume there is one more
c    number, B(N).  If we let A be the symmetric Toeplitz matrix whose first
c    row is B(0:N-1), then the Yule-Walker equations are:
c
c      A * X = -B(1:N)
c
c  Example:
c
c    To solve
c
c     1.0 0.5 0.2    x1   0.5
c     0.5 1.0 0.5 *  x2 = 0.2
c     0.2 0.5 1.0    x3   0.1
c
c    we input:
c
c      N = 3
c      B(1:3) = (/ 0.5, 0.2, 0.1 /)
c
c    with output:
c
c      X(1:3) = (/ -75, 12, -5 /) / 140
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 March 2001
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Gene Golub, Charles Van Loan,
c    Section 4.7.2, "Solving the Yule-Walker Equations",
c    Matrix Computations,
c    Third Edition,
c    Johns Hopkins, 1996.
c
c  Parameters:
c
c    Input, integer N, the order of the system.
c
c    Input, double precision B(N), defines the linear system.  The first
c    entry of A is a 1, followed by B(1) through B(N-1).  The right hand
c    side of the system is -B(1:N).
c
c    Output, double precision X(N), the solution of the linear system.
c
      implicit none

      integer n

      double precision alpha
      double precision b(n)
      double precision beta
      integer i
      double precision x(n)

      x(1) = -b(1)
      beta = 1.0D+00
      alpha = -b(1)

      do i = 1, n-1
        beta = ( 1.0D+00 - alpha * alpha ) * beta
        alpha = - ( b(i+1) + sum ( b(i:1:-1) * x(1:i) ) ) / beta
        x(1:i) = x(1:i) + alpha * x(i:1:-1)
        x(i+1) = alpha
      end do

      return
      end
      subroutine r8to_indicator ( n, a )

c*********************************************************************72
c
cc R8TO_INDICATOR sets up an R8TO indicator matrix.
c
c  Discussion:
c
c    The "indicator matrix" simply has a value like I*10+J at every
c    entry of a dense matrix, or at every entry of a compressed storage
c    matrix for which storage is allocated. 
c
c    The R8TO storage format is used for a real Toeplitz matrix, which 
c    is constant along diagonals.  Thus, in an N by N Toeplitz matrix, 
c    there are at most 2*N-1 distinct entries.  The format stores the 
c    N elements of the first row, followed by the N-1 elements of the 
c    first column (skipping the entry in the first row).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 January 2004
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
c    Output, double precision A(2*N-1), the R8TO matrix.
c
      implicit none

      integer n

      double precision a(2*n-1)
      integer fac
      integer i
      integer i4_log_10
      integer j
      integer k

      fac = 10 ** ( i4_log_10 ( n ) + 1 )

      i = 1
      k = 0
      do j = 1, n
        k = k + 1
        a(k) = dble ( fac * i + j )
      end do

      j = 1
      do i = 2, n
        k = k + 1
        a(k) = dble ( fac * i + j )
      end do
      
      return
      end
      subroutine r8to_mxv ( n, a, x, b )

c*********************************************************************72
c
cc R8TO_MXV multiplies an R8TO matrix by an R8VEC.
c
c  Discussion:
c
c    The R8TO storage format is used for a Toeplitz matrix, which is constant
c    along diagonals.  Thus, in an N by N Toeplitz matrix, there are at most 
c    2*N-1 distinct entries.  The format stores the N elements of the first
c    row, followed by the N-1 elements of the first column (skipping the
c    entry in the first row).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 March 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision A(2*N-1), the R8TO matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(N), the product A * x.
c
      implicit none

      integer n

      double precision a(2*n-1)
      double precision b(n)
      integer i
      double precision x(n)

      b(1) = sum ( a(1:n) * x(1:n) )

      do i = 2, n
        b(i) = sum ( a(n+i-1:n+1:-1) * x(1:i-1) ) 
     &       + sum ( a(1:n+1-i) * x(i:n) )
      end do

      return
      end
      subroutine r8to_print ( n, a, title )

c*********************************************************************72
c
cc R8TO_PRINT prints an R8TO matrix.
c
c  Discussion:
c
c    The R8TO storage format is used for a Toeplitz matrix, which is constant
c    along diagonals.  Thus, in an N by N Toeplitz matrix, there are at most 
c    2*N-1 distinct entries.  The format stores the N elements of the first
c    row, followed by the N-1 elements of the first column (skipping the
c    entry in the first row).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 May 2000
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
c    Input, double precision A(2*N-1), the R8TO matrix.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer n

      double precision a(2*n-1)
      character ( len = * ) title

      call r8to_print_some ( n, a, 1, 1, n, n, title )

      return
      end
      subroutine r8to_print_some ( n, a, ilo, jlo, ihi, jhi, title )

c*********************************************************************72
c
cc R8TO_PRINT_SOME prints some of an R8TO matrix.
c
c  Discussion:
c
c    The R8TO storage format is used for a Toeplitz matrix, which is constant
c    along diagonals.  Thus, in an N by N Toeplitz matrix, there are at most 
c    2*N-1 distinct entries.  The format stores the N elements of the first
c    row, followed by the N-1 elements of the first column (skipping the
c    entry in the first row).
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
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input, double precision A(2*N-1), the R8TO matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer, parameter :: incx = 5
      integer n

      double precision a(2*n-1)
      double precision aij
      character ( len = 14 ) ctemp(incx)
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
      character ( len = * )  title

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
        i2hi = min ( ihi, n )

        do i = i2lo, i2hi
c
c  Print out (up to) 5 entries in row I, that lie in the current strip.
c
          do j2 = 1, inc

            j = j2lo - 1 + j2

            if ( i .le. j ) then
              aij = a(j+1-i)
            else
              aij = a(n+i-j)
            end if

            if ( r8_is_int ( aij ) ) then
              write ( ctemp(j2), '(f8.0,6x)' ) aij
            else
              write ( ctemp(j2), '(g14.6)' ) aij
            end if

          end do

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine r8to_random ( n, seed, a )

c*********************************************************************72
c
cc R8TO_RANDOM randomizes an R8TO matrix.
c
c  Discussion:
c
c    The R8TO storage format is used for a Toeplitz matrix, which is constant
c    along diagonals.  Thus, in an N by N Toeplitz matrix, there are at most 
c    2*N-1 distinct entries.  The format stores the N elements of the first
c    row, followed by the N-1 elements of the first column (skipping the
c    entry in the first row).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 May 2003
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
c    Input/output, integer SEED, a seed for the random number
c    generator.
c
c    Output, double precision A(2*N-1), the R8TO matrix.
c
      implicit none

      integer n

      double precision a(2*n-1)
      integer n2
      integer seed

      n2 = 2 * n - 1

      call r8vec_uniform_01 ( n2, seed, a )

      return
      end
      subroutine r8to_sl ( n, a, b, x, job )

c*********************************************************************72
c
cc R8TO_SL solves an R8TO system.
c
c  Discussion:
c
c    The R8TO storage format is used for a Toeplitz matrix, which is constant
c    along diagonals.  Thus, in an N by N Toeplitz matrix, there are at most 
c    2*N-1 distinct entries.  The format stores the N elements of the first
c    row, followed by the N-1 elements of the first column (skipping the
c    entry in the first row).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 March 2001
c
c  Author:
c
c    FORTRAN90 version by John Burkardt.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision A(2*N-1), the R8TO matrix.
c
c    Input, double precision B(N) the right hand side vector.
c
c    Output, double precision X(N), the solution vector.  X and B may share the
c    same storage.
c
c    Input, integer JOB,
c    0 to solve A*X=B,
c    nonzero to solve A'*X=B.
c
      implicit none

      integer n

      double precision a(2*n-1)
      double precision b(n)
      double precision c1(n-1)
      double precision c2(n-1)
      integer i
      integer job
      integer nsub
      double precision r1
      double precision r2
      double precision r3
      double precision r5
      double precision r6
      double precision x(n)

      if ( n .lt. 1 ) then
        return
      end if
c
c  Solve the system with the principal minor of order 1.
c
      r1 = a(1)
      x(1) = b(1) / r1

      if ( n .eq. 1 ) then
        return
      end if
c
c  Recurrent process for solving the system with the Toeplitz matrix.
c
      do nsub = 2, n
c
c  Compute multiples of the first and last columns of the inverse of
c  the principal minor of order NSUB.
c
        if ( job .eq. 0 ) then
          r5 = a(n+nsub-1)
          r6 = a(nsub)
        else
          r5 = a(nsub)
          r6 = a(n+nsub-1)
        end if

        if ( 2 .lt. nsub ) then

          c1(nsub-1) = r2

          do i = 1, nsub-2
            if ( job .eq. 0 ) then
              r5 = r5 + a(n+i) * c1(nsub-i)
              r6 = r6 + a(i+1) * c2(i)
            else
              r5 = r5 + a(i+1) * c1(nsub-i)
              r6 = r6 + a(n+i) * c2(i)
            end if
          end do

        end if

        r2 = -r5 / r1
        r3 = -r6 / r1
        r1 = r1 + r5 * r3

        if ( 2 .lt. nsub ) then

          r6 = c2(1)
          c2(nsub-1) = 0.0D+00

          do i = 2, nsub-1
            r5 = c2(i)
            c2(i) = c1(i) * r3 + r6
            c1(i) = c1(i) + r6 * r2
            r6 = r5
          end do

        end if

        c2(1) = r3
c
c  Compute the solution of the system with the principal minor of order NSUB.
c
        if ( job .eq. 0 ) then
          r5 = sum ( a(n+1:n+nsub-1) * x(nsub-1:1:-1) )
        else
          r5 = sum ( a(2:nsub) * x(nsub-1:1:-1) )
        end if

        r6 = ( b(nsub) - r5 ) / r1

        x(1:nsub-1) = x(1:nsub-1) + c2(1:nsub-1) * r6
        x(nsub) = r6

      end do

      return
      end
      subroutine r8to_to_r8ge ( n, a, b )

c*********************************************************************72
c
cc R8TO_TO_R8GE copies an R8TO matrix to an R8GE matrix.
c
c  Discussion:
c
c    The R8TO storage format is used for a Toeplitz matrix, which is constant
c    along diagonals.  Thus, in an N by N Toeplitz matrix, there are at most 
c    2*N-1 distinct entries.  The format stores the N elements of the first
c    row, followed by the N-1 elements of the first column (skipping the
c    entry in the first row).
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 March 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision A(2*N-1), the R8TO matrix.
c
c    Output, double precision B(N,N), the R8GE matrix.
c
      implicit none

      integer n

      double precision a(2*n-1)
      double precision b(n,n)
      integer i

      do i = 1, n
        b(i,1:i-1) = a(n+i-1:n+1:-1)
        b(i,i:n) = a(1:n-i+1)
      end do

      return
      end
      subroutine r8to_vxm ( n, a, x, b )

c*********************************************************************72
c
cc R8TO_VXM multiplies an R8VEC by an R8TO matrix.
c
c  Discussion:
c
c    The R8TO storage format is used for a Toeplitz matrix, which is constant
c    along diagonals.  Thus, in an N by N Toeplitz matrix, there are at most 
c    2*N-1 distinct entries.  The format stores the N elements of the first
c    row, followed by the N-1 elements of the first column (skipping the
c    entry in the first row).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 November 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision A(2*N-1), the R8TO matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(N), the product A' * X.
c
      implicit none

      integer n

      double precision a(2*n-1)
      double precision b(n)
      integer i
      double precision x(n)

      do i = 1, n

        b(i) = sum ( a(i:1:-1) * x(1:i) ) + 
     &         sum ( a(n+1:2*n-i) * x(i+1:n) )

      end do

      return
      end
      subroutine r8ut_det ( n, a, det )

c*********************************************************************72
c
cc R8UT_DET computes the determinant of an R8UT matrix.
c
c  Discussion:
c
c    The R8UT storage format is used for an M by N upper triangular 
c    matrix.  The format stores all M*N entries, even those which are zero.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 August 1999
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
c    Input, double precision A(N,N), the R8UT matrix.
c
c    Output, double precision DET, the determinant of the matrix.
c
      implicit none

      integer n

      double precision a(n,n)
      double precision det
      double precision diag(n)

      call r8mat_diag_get_vector ( n, a, diag )

      det = product ( diag(1:n) )

      return
      end
      subroutine r8ut_indicator ( m, n, a )

c*********************************************************************72
c
cc R8UT_INDICATOR sets up an R8UT indicator matrix.
c
c  Discussion:
c
c    The "indicator matrix" simply has a value like I*10+J at every
c    entry of a dense matrix, or at every entry of a compressed storage
c    matrix for which storage is allocated. 
c
c    The R8UT storage format is used for an M by N upper triangular 
c    matrix.  The format stores all M*N entries, even those which are zero.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 January 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns of 
c    the matrix.  M and N must be positive.
c
c    Output, double precision A(M,N), the R8UT matrix.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer fac
      integer i
      integer i4_log_10
      integer j

      fac = 10 ** ( i4_log_10 ( n ) + 1 )

      do i = 1, m
        do j = 1, min ( i-1, n )
          a(i,j) = 0.0D+00
        end do
        do j = i, n
          a(i,j) = dble ( fac * i + j )
        end do
      end do

      return
      end
      subroutine r8ut_inverse ( n, a )

c*********************************************************************72
c
cc R8UT_INVERSE computes the inverse of an R8UT matrix.
c
c  Discussion:
c
c    The R8UT storage format is used for an M by N upper triangular 
c    matrix.  The format stores all M*N entries, even those which are zero.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms,
c    Academic Press, 1978, second edition,
c    ISBN 0-12-519260-6
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input/output, double precision A(N,N).
c    On input, the upper triangular matrix to be inverted.
c    On output, the inverse of the upper triangular matrix.
c
      implicit none

      integer n

      double precision a(n,n)
      integer i
      integer j
c
c  Check.
c
      do i = 1, n
        if ( a(i,i) .eq. 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8UT_INVERSE - Fatal error!'
          write ( *, '(a)' ) '  Zero diagonal element.'
          stop
        end if
      end do

      do j = n, 1, -1

        do i = n, 1, -1

          if ( j .lt. i ) then

            a(i,j) = 0.0D+00

          else if ( i .eq. j ) then

            a(i,j) = 1.0D+00 / a(i,j)

          else if ( i .lt. j ) then

            a(i,j) = - sum ( a(i,i+1:j) * a(i+1:j,j) ) / a(i,i)

          end if

        end do
      end do

      return
      end
      subroutine r8ut_mxm ( n, a, b, c )

c*********************************************************************72
c
cc R8UT_MXM multiplies two R8UT matrices.
c
c  Discussion:
c
c    The R8UT storage format is used for an M by N upper triangular 
c    matrix.  The format stores all M*N entries, even those which are zero.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 September 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrices.
c    N must be positive.
c
c    Input, double precision A(N,N), B(N,N), the R8UT factor matrices.
c
c    Output, double precision C(N,N), the R8UT product matrix.
c
      implicit none

      integer n

      double precision a(n,n)
      double precision b(n,n)
      double precision c(n,n)
      integer i
      integer j

      do i = 1, n
        do j = 1, i-1
          c(i,j) = 0.0D+00
        end do
        do j = i, n
          c(i,j) = sum ( a(i,i:j) * b(i:j,j) )
        end do
      end do

      return
      end
      subroutine r8ut_mxv ( m, n, a, x, b )

c*********************************************************************72
c
cc R8UT_MXV multiplies an R8UT matrix by an R8VEC.
c
c  Discussion:
c
c    The R8UT storage format is used for an M by N upper triangular 
c    matrix.  The format stores all M*N entries, even those which are zero.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input, double precision A(M,N), the R8UT matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(M), the product A * x.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision b(m)
      integer i
      double precision x(n)

      do i = 1, m
        b(i) = sum ( a(i,i:n) * x(i:n) )
      end do

      return
      end
      subroutine r8ut_print ( m, n, a, title )

c*********************************************************************72
c
cc R8UT_PRINT prints an R8UT matrix.
c
c  Discussion:
c
c    The R8UT storage format is used for an M by N upper triangular 
c    matrix.  The format stores all M*N entries, even those which are zero.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 May 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input, double precision A(M,N), the R8UT matrix.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      character ( len = * )  title

      call r8ut_print_some ( m, n, a, 1, 1, m, n, title )

      return
      end
      subroutine r8ut_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

c*********************************************************************72
c
cc R8UT_PRINT_SOME prints some of an R8UT matrix.
c
c  Discussion:
c
c    The R8UT storage format is used for an M by N upper triangular 
c    matrix.  The format stores all M*N entries, even those which are zero.
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
c    Input, integer M, the number of rows of the matrix.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input, double precision A(M,N), the R8UT matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer, parameter :: incx = 5
      integer m
      integer n

      double precision a(m,n)
      character ( len = 14 ) ctemp(incx)
      logical                r8_is_int
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
      character ( len = * )  title

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
        i2hi = min ( ihi, m )
        i2hi = min ( i2hi, j2hi )

        do i = i2lo, i2hi
c
c  Print out (up to) 5 entries in row I, that lie in the current strip.
c
          do j2 = 1, inc

            j = j2lo - 1 + j2

            if ( j .lt. i ) then
              ctemp(j2) = '              '
            else if ( r8_is_int ( a(i,j) ) ) then
              write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
            else
              write ( ctemp(j2), '(g14.6)' ) a(i,j)
            end if

          end do

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine r8ut_random ( m, n, seed, a )

c*********************************************************************72
c
cc R8UT_RANDOM randomizes an R8UT matrix.
c
c  Discussion:
c
c    The R8UT storage format is used for an M by N upper triangular 
c    matrix.  The format stores all M*N entries, even those which are zero.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 October 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns of 
c    the matrix.  M and N must be positive.
c
c    Input/output, integer SEED, a seed for the random number 
c    generator.
c
c    Output, double precision A(M,N), the R8UT matrix.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision r8_uniform_01
      integer i
      integer j
      integer seed

      do i = 1, m
        do j = 1, min ( i-1, n )
          a(i,j) = 0.0D+00
        end do
        do j = i, n
          a(i,j) = r8_uniform_01 ( seed )
        end do
      end do

      return
      end
      subroutine r8ut_sl ( n, a, b, job )

c*********************************************************************72
c
cc R8UT_SL solves an R8UT system.
c
c  Discussion:
c
c    The R8UT storage format is used for an M by N upper triangular 
c    matrix.  The format stores all M*N entries, even those which are zero.
c
c    No factorization of the upper triangular matrix is required.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 September 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision A(N,N), the R8UT matrix.
c
c    Input/output, double precision B(N).
c    On input, the right hand side.
c    On output, the solution vector.
c
c    Input, integer JOB, is 0 to solve the untransposed system,
c    nonzero to solve the transposed system.
c
      implicit none

      integer n

      double precision a(n,n)
      double precision b(n)
      integer j
      integer job

      if ( job .eq. 0 ) then

        do j = n, 1, -1
          b(j) = b(j) / a(j,j)
          b(1:j-1) = b(1:j-1) - a(1:j-1,j) * b(j)
        end do

      else

        do j = 1, n
          b(j) = b(j) / a(j,j)
          b(j+1:n) = b(j+1:n) - a(j,j+1:n) * b(j)
        end do

      end if

      return
      end
      subroutine r8ut_vxm ( m, n, a, x, b )

c*********************************************************************72
c
cc R8UT_VXM multiplies an R8VEC by an R8UT matrix.
c
c  Discussion:
c
c    The R8UT storage format is used for an M by N upper triangular 
c    matrix.  The format stores all M*N entries, even those which are zero.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c    M must be positive.
c
c    Input, integer N, the number of columns of the matrix.
c    N must be positive.
c
c    Input, double precision A(M,N), the R8UT matrix.
c
c    Input, double precision X(M), the vector to be multiplied by A.
c
c    Output, double precision B(N), the product A' * x.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision b(n)
      integer i
      integer jhi
      double precision x(m)

      do i = 1, n
        jhi = min ( i, m )
        b(i) = sum ( x(1:jhi) * a(1:jhi,i) )
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
c    An R8VEC is an array of double precision real values.
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
      subroutine r8vec_print_some ( n, a, max_print, title )

c*********************************************************************72
c
cc R8VEC_PRINT_SOME prints "some" of an R8VEC.
c
c  Discussion:
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
      integer s_len_trim
      character*(*) title

      if ( max_print .le. 0 ) then
        return
      end if

      if ( n .le. 0 ) then
        return
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) title
      write ( *, '(a)' ) ' '

      if ( n .le. max_print ) then

        do i = 1, n
          write ( *, '(i6,2x,g14.6)' ) i, a(i)
        end do

      else if ( 3 .le. max_print ) then

        do i = 1, max_print-2
          write ( *, '(i6,2x,g14.6)' ) i, a(i)
        end do

        write ( *, '(a)' ) '......  ..............'
        i = n

        write ( *, '(i6,2x,g14.6)' ) i, a(i)

      else

        do i = 1, max_print-1
          write ( *, '(i6,2x,g14.6)' ) i, a(i)
        end do

        i = max_print

        write ( *, '(i6,2x,g14.6,a)' ) i, a(i), '...more entries...'

      end if

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
      subroutine r8vec_read ( input_file, n, r )

c*********************************************************************72
c
cc R8VEC_READ reads an R8VEC from a file.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    31 July 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character ( len = * ) INPUT_FILE, the name of the file to be read.
c
c    Input, integer N, the size of the vector.
c
c    Output, double precision R(N), the vector.
c
      implicit none

      integer n

      character ( len = * ) input_file
      integer input_unit
      integer ios
      integer k
      double precision r(n)

      call get_unit ( input_unit )

      open ( unit = input_unit, file = input_file, status = 'old', 
     &  iostat = ios )

      if ( ios .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8VEC_READ - Fatal errorc'
        write ( *, '(a)' ) '  Could not open the input file "' 
     &    // trim ( input_file ) // '".'
        stop
      end if

      do k = 1, n

        read ( input_unit, *, iostat = ios ) r(k)

        if ( ios .ne. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8VEC_READ - Fatal errorc'
          write ( *, '(a,i8)' ) 
     &      '  I/O error while reading record ', k
          stop
        end if

      end do

      close ( unit = input_unit )

      return
      end
      subroutine r8vec_read_size ( input_file, n )

c*********************************************************************72
c
cc R8VEC_READ_SIZE reads the size of an R8VEC from a file.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    31 July 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character ( len = * ) INPUT_FILE, the name of the file to 
c    be read.
c
c    Output, integer N, the size of the vector.
c
      implicit none

      character ( len = * )  input_file
      integer input_unit
      integer ios
      integer n
      double precision r

      call get_unit ( input_unit )

      open ( unit = input_unit, file = input_file, status = 'old', 
     &  iostat = ios )

      if ( ios .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8VEC_READ_SIZE - Fatal errorc'
        write ( *, '(a)' ) '  Could not open the input file "' 
     &    // trim ( input_file ) // '".'
        stop
      end if

      n = 0

      do

        read ( input_unit, *, iostat = ios ) r

        if ( ios .ne. 0 ) then
          exit
        end if

        n = n + 1

      end do

      close ( unit = input_unit )

      return
      end
      subroutine r8vec_to_r8cb ( m, n, ml, mu, x, a )

c*********************************************************************72
c
cc R8VEC_TO_R8CB copies an R8VEC into an R8CB matrix.
c
c  Discussion:
c
c    In C++ and FORTRAN, this routine is not really needed.  In MATLAB,
c    a data item carries its dimensionality implicitly, and so cannot be
c    regarded sometimes as a vector and sometimes as an array.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 March 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns 
c    in the array.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c
c    Input, double precision X((ML+MU+1)*N), the vector to be copied
c    into the array.
c
c    Output, double precision A(ML+MU+1,N), the array.
c
      implicit none

      integer m
      integer ml
      integer mu
      integer n

      double precision a(ml+mu+1,n)
      integer i
      integer j
      double precision x((ml+mu+1)*n)

      do j = 1, n
        do i = 1, ml + mu + 1

          if ( 1 .le. i + j - mu - 1 .and. i + j - mu - 1 .le. m ) then
            a(i,j) = x(i+(ml+mu+1)*(j-1))
          else
            a(i,j) = 0.0D+00
          end if

        end do
      end do

      return
      end
      subroutine r8vec_to_r8gb ( m, n, ml, mu, x, a )

c*********************************************************************72
c
cc R8VEC_TO_R8GB copies an R8VEC into an R8GB matrix.
c
c  Discussion:
c
c    In C++ and FORTRAN, this routine is not really needed.  In MATLAB,
c    a data item carries its dimensionality implicitly, and so cannot be
c    regarded sometimes as a vector and sometimes as an array.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 March 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns 
c    in the array.
c
c    Input, integer ML, MU, the lower and upper bandwidths.
c
c    Input, double precision X((2*ML+MU+1)*N), the vector to be copied
c    into the array.
c
c    Output, double precision A(2*ML+MU+1,N), the array.
c
      implicit none

      integer m
      integer ml
      integer mu
      integer n

      double precision a(2*ml+mu+1,n)
      integer i
      integer j
      double precision x((2*ml+mu+1)*n)

      do j = 1, n
        do i = 1, 2 * ml + mu + 1

          if ( 1 .le. i + j - ml - mu - 1 .and. 
     &      i + j - ml - mu - 1 .le. m ) then
            a(i,j) = x(i+(2*ml+mu+1)*(j-1))
          else
            a(i,j) = 0.0D+00
          end if

        end do
      end do

      return
      end
      subroutine r8vec_to_r8ge ( m, n, x, a )

c*********************************************************************72
c
cc R8VEC_TO_R8GE copies an R8VEC into an R8GE matrix.
c
c  Discussion:
c
c    In C++ and FORTRAN, this routine is not really needed.  In MATLAB,
c    a data item carries its dimensionality implicitly, and so cannot be
c    regarded sometimes as a vector and sometimes as an array.
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns 
c    in the array.
c
c    Input, double precision X(M*N), the vector to be copied into the array.
c
c    Output, double precision A(M,N), the array.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer i
      integer j
      integer k
      double precision x(m*n)

      k = 0
      do j = 1, n
        do i = 1, m
          k = k + 1
          a(i,j) = x(k)
        end do
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
      subroutine r8vec_write ( n, r, output_file )

c*********************************************************************72
c
cc R8VEC_WRITE writes an R8VEC to a file.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 September 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision R(N), the vector to be written.
c
c    Input, character ( len = * ) OUTPUT_FILE, the name of the file to which
c    the information is to be written.
c
      implicit none

      integer n

      integer i
      integer ios
      character ( len = * ) output_file
      integer output_unit
      double precision r(n)

      call get_unit ( output_unit )

      open ( unit = output_unit, file = output_file, status = 'replace', 
     &  iostat = ios )

      if ( ios .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8VEC_WRITE - Fatal errorc'
        write ( *, '(a)' ) '  Could not open the output file "' 
     &    // trim ( output_file ) // '".'
        stop
      end if

      do i = 1, n
        write ( output_unit, '(2x,g16.8)' ) r(i)
      end do

      close ( unit = output_unit )

      return
      end
      subroutine r8vec2_print_some ( n, x1, x2, max_print, title )

c*********************************************************************72
c
cc R8VEC2_PRINT_SOME prints "some" of a pair of R8VEC's.
c
c  Discussion:
c
c    The user specifies MAX_PRINT, the maximum number of lines to print.
c
c    If N, the size of the vectors, is no more than MAX_PRINT, then
c    the entire vectors are printed, one entry of each per line.
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
c    17 December 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries of the vectors.
c
c    Input, double precision X1(N), X2(N), the vector to be printed.
c
c    Input, integer MAX_PRINT, the maximum number of lines 
c    to print.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer n

      integer i
      integer max_print
      character ( len = * )  title
      double precision x1(n)
      double precision x2(n)

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
          write ( *, '(i8,2x,g14.6,2x,g14.6)' ) i, x1(i), x2(i)
        end do

      else if ( 3 .le. max_print ) then

        do i = 1, max_print-2
          write ( *, '(i8,2x,g14.6,2x,g14.6)' ) i, x1(i), x2(i)
        end do
        write ( *, '(a)' ) '......  ..............  ..............'
        i = n
        write ( *, '(i8,2x,g14.6,2x,g14.6)' ) i, x1(i), x2(i)

      else

        do i = 1, max_print - 1
          write ( *, '(i8,2x,g14.6,2x,g14.6)' ) i, x1(i), x2(i)
        end do
        i = max_print
        write ( *, '(i8,2x,g14.6,2x,g14.6,2x,a)' ) i, x1(i), x2(i), 
     &    '...more entries...'

      end if

      return
      end
      subroutine r8vm_det ( n, a, det )

c*********************************************************************72
c
cc R8VM_DET computes the determinant of an R8VM matrix.
c
c  Discussion:
c
c    The R8VM storage format is used for an M by N Vandermonde matrix.
c    An M by N Vandermonde matrix is defined by the values in its second
c    row, which will be written here as X(1:N).  The matrix has a first 
c    row of 1's, a second row equal to X(1:N), a third row whose entries
c    are the squares of the X values, up to the M-th row whose entries
c    are the (M-1)th powers of the X values.  The matrix can be stored
c    compactly by listing just the values X(1:N).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 November 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of rows and columns of 
c    the matrix.
c
c    Input, double precision A(N), the R8VM matrix.
c
c    Output, double precision DET, the determinant of the matrix.
c
      implicit none

      integer n

      double precision a(n)
      double precision det
      integer i
      integer j

      det = 1.0D+00
      do j = 1, n
        do i = j+1, n
          det = det * ( a(i) - a(j) )
        end do
      end do

      return
      end
      subroutine r8vm_mxv ( m, n, a, x, b )

c*********************************************************************72
c
cc R8VM_MXV multiplies an R8VM matrix by an R8VEC.
c
c  Discussion:
c
c    The R8VM storage format is used for an M by N Vandermonde matrix.
c    An M by N Vandermonde matrix is defined by the values in its second
c    row, which will be written here as X(1:N).  The matrix has a first 
c    row of 1's, a second row equal to X(1:N), a third row whose entries
c    are the squares of the X values, up to the M-th row whose entries
c    are the (M-1)th powers of the X values.  The matrix can be stored
c    compactly by listing just the values X(1:N).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 September 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns of 
c    the matrix.
c
c    Input, double precision A(N), the R8VM matrix.
c
c    Input, double precision X(N), the vector to be multiplied by A.
c
c    Output, double precision B(M), the product A * x.
c
      implicit none

      integer m
      integer n

      double precision a(n)
      double precision b(m)
      integer i
      integer j
      double precision x(n)

      do i = 1, m
        b(i) = 0.0D+00
        do j = 1, n
          if ( i .eq. 1 ) then
            b(i) = b(i) + x(j)
          else
            b(i) = b(i) + a(j)**(i-1) * x(j)
          end if
        end do
      end do

      return
      end
      subroutine r8vm_print ( m, n, a, title )

c*********************************************************************72
c
cc R8VM_PRINT prints an R8VM matrix.
c
c  Discussion:
c
c    The R8VM storage format is used for an M by N Vandermonde matrix.
c    An M by N Vandermonde matrix is defined by the values in its second
c    row, which will be written here as X(1:N).  The matrix has a first 
c    row of 1's, a second row equal to X(1:N), a third row whose entries
c    are the squares of the X values, up to the M-th row whose entries
c    are the (M-1)th powers of the X values.  The matrix can be stored
c    compactly by listing just the values X(1:N).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 September 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns of 
c    the matrix.
c
c    Input, double precision A(N), the R8VM matrix.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer n

      double precision a(n)
      integer m
      character ( len = * )  title

      call r8vm_print_some ( m, n, a, 1, 1, n, n, title )

      return
      end
      subroutine r8vm_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

c*********************************************************************72
c
cc R8VM_PRINT_SOME prints some of an R8VM matrix.
c
c  Discussion:
c
c    The R8VM storage format is used for an M by N Vandermonde matrix.
c    An M by N Vandermonde matrix is defined by the values in its second
c    row, which will be written here as X(1:N).  The matrix has a first 
c    row of 1's, a second row equal to X(1:N), a third row whose entries
c    are the squares of the X values, up to the M-th row whose entries
c    are the (M-1)th powers of the X values.  The matrix can be stored
c    compactly by listing just the values X(1:N).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 September 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns of 
c    the matrix.
c
c    Input, double precision A(N), the R8VM matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer, parameter :: incx = 5
      integer n

      double precision a(n)
      double precision aij
      character ( len = 14 ) ctemp(incx)
      logical                r8_is_int
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
      integer m
      character ( len = * )  title

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
        i2hi = min ( ihi, m )

        do i = i2lo, i2hi
c
c  Print out (up to) 5 entries in row I, that lie in the current strip.
c
          do j2 = 1, inc

            j = j2lo - 1 + j2

            if ( i .eq. 1 ) then
              aij = 1.0D+00
            else
              aij = a(j)**(i-1)
            end if

            if ( r8_is_int ( aij ) ) then
              write ( ctemp(j2), '(f8.0,6x)' ) aij
            else
              write ( ctemp(j2), '(g14.6)' ) aij
            end if

          end do

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine r8vm_random ( m, n, seed, a )

c*********************************************************************72
c
cc R8VM_RANDOM randomizes an R8VM matrix.
c
c  Discussion:
c
c    The R8VM storage format is used for an M by N Vandermonde matrix.
c    An M by N Vandermonde matrix is defined by the values in its second
c    row, which will be written here as X(1:N).  The matrix has a first 
c    row of 1's, a second row equal to X(1:N), a third row whose entries
c    are the squares of the X values, up to the M-th row whose entries
c    are the (M-1)th powers of the X values.  The matrix can be stored
c    compactly by listing just the values X(1:N).
c
c    The parameter M is not actually needed by this routine.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 September 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns of 
c    the matrix.
c
c    Input/output, integer SEED, a seed for the random number 
c    generator.
c
c    Output, double precision A(N), the R8VM matrix.
c
      implicit none

      integer n

      double precision a(n)
      integer m
      integer seed

      call r8vec_uniform_01 ( n, seed, a )

      return
      end
      subroutine r8vm_sl ( n, a, b, x, job, info )

c*********************************************************************72
c
cc R8VM_SL solves an R8VM linear system.
c
c  Discussion:
c
c    The R8VM storage format is used for an M by N Vandermonde matrix.
c    An M by N Vandermonde matrix is defined by the values in its second
c    row, which will be written here as X(1:N).  The matrix has a first 
c    row of 1's, a second row equal to X(1:N), a third row whose entries
c    are the squares of the X values, up to the M-th row whose entries
c    are the (M-1)th powers of the X values.  The matrix can be stored
c    compactly by listing just the values X(1:N).
c
c    Vandermonde systems are very close to singularity.  The singularity
c    gets worse as N increases, and as any pair of values defining
c    the matrix get close.  Even a system as small as N = 10 will
c    involve the 9th power of the defining values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 September 2003
c
c  Author:
c
c    John Burkardt.
c
c  Reference:
c
c    Gene Golub, Charles Van Loan,
c    Matrix Computations,
c    Third Edition,
c    Johns Hopkins, 1996.
c
c  Parameters:
c
c    Input, integer N, the number of rows and columns of 
c    the matrix.
c
c    Input, double precision A(N), the R8VM matrix.
c
c    Input, double precision B(N), the right hand side.
c
c    Output, double precision X(N), the solution of the linear system.
c
c    Input, integer JOB, specifies the system to solve.
c    0, solve A * x = b.
c    nonzero, solve A' * x = b.
c
c    Output, integer INFO.
c    0, no error.
c    nonzero, at least two of the values in A are equal.
c
      implicit none

      integer n

      double precision a(n)
      double precision b(n)
      integer i
      integer info
      integer j
      integer job
      double precision x(n)
c
c  Check for explicit singularity.
c
      info = 0

      do j = 1, n - 1
        do i = j+1, n
          if ( a(i) .eq. a(j) ) then
            info = 1
            return
          end if
        end do
      end do

      x(1:n) = b(1:n)

      if ( job .eq. 0 ) then

        do j = 1, n-1
          do i = n, j+1, -1
            x(i) = x(i) - a(j) * x(i-1)
          end do
        end do

        do j = n-1, 1, -1

          do i = j+1, n
            x(i) = x(i) / ( a(i) - a(i-j) )
          end do

          do i = j, n-1
            x(i) = x(i) - x(i+1)
          end do

        end do

      else

        do j = 1, n-1
          do i = n, j+1, -1
            x(i) = ( x(i) - x(i-1) ) / ( a(i) - a(i-j) )
          end do
        end do

        do j = n-1, 1, -1
          do i = j, n-1
            x(i) = x(i) - x(i+1) * a(j)
          end do
        end do

      end if

      return
      end
      subroutine r8vm_to_r8ge ( m, n, a, b )

c*********************************************************************72
c
cc R8VM_TO_R8GE copies an R8VM matrix to an R8GE matrix.
c
c  Discussion:
c
c    The R8VM storage format is used for an M by N Vandermonde matrix.
c    An M by N Vandermonde matrix is defined by the values in its second
c    row, which will be written here as X(1:N).  The matrix has a first 
c    row of 1's, a second row equal to X(1:N), a third row whose entries
c    are the squares of the X values, up to the M-th row whose entries
c    are the (M-1)th powers of the X values.  The matrix can be stored
c    compactly by listing just the values X(1:N).
c
c    The R8GE storage format is used for a general M by N matrix.  A storage 
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns of 
c    the matrix.
c
c    Input, double precision A(N), the R8VM matrix.
c
c    Output, double precision B(M,N), the R8GE matrix.
c
      implicit none

      integer m
      integer n

      double precision a(n)
      double precision b(m,n)
      integer i
      integer j

      do i = 1, m
        do j = 1, n
          if ( i .eq. 1 ) then
            b(i,j) = 1.0D+00
          else
            b(i,j) = b(i-1,j) * a(j)
          end if
        end do
      end do

      return
      end
      subroutine r8vm_vxm ( m, n, a, x, b )

c*********************************************************************72
c
cc R8VM_VXM multiplies an R8VEC by an R8VM matrix.
c
c  Discussion:
c
c    The R8VM storage format is used for an M by N Vandermonde matrix.
c    An M by N Vandermonde matrix is defined by the values in its second
c    row, which will be written here as X(1:N).  The matrix has a first 
c    row of 1's, a second row equal to X(1:N), a third row whose entries
c    are the squares of the X values, up to the M-th row whose entries
c    are the (M-1)th powers of the X values.  The matrix can be stored
c    compactly by listing just the values X(1:N).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 September 2003
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns of 
c    the matrix.
c
c    Input, double precision A(N), the R8VM matrix.
c
c    Input, double precision X(M), the vector to be multiplied by A.
c
c    Output, double precision B(N), the product A' * x.
c
      implicit none

      integer m
      integer n

      double precision a(n)
      double precision b(n)
      integer i
      integer j
      double precision x(m)

      do j = 1, n
        b(j) = 0.0D+00
        do i = 1, m
          if ( i .eq. 1 ) then
            b(j) = b(j) + x(i)
          else
            b(j) = b(j) + a(j)**(i-1) * x(i)
          end if
        end do
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
      subroutine sort_heap_external ( n, indx, i, j, isgn )

c*********************************************************************72
c
cc SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
c
c  Discussion:
c
c    The actual list of data is not passed to the routine.  Hence this
c    routine may be used to sort integers, reals, numbers, names,
c    dates, shoe sizes, and so on.  After each call, the routine asks
c    the user to compare or interchange two items, until a special
c    return value signals that the sorting is completed.
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
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
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
c    Input, integer N, the number of items to be sorted.
c
c    Input/output, integer INDX, the main communication signal.
c
c    The user must set INDX to 0 before the first call.
c    Thereafter, the user should not change the value of INDX until
c    the sorting is done.
c
c    On return, if INDX is
c
c      greater than 0,
c      * interchange items I and J;
c      * call again.
c
c      less than 0,
c      * compare items I and J;
c      * set ISGN = -1 if I .lt. J, ISGN = +1 if J .lt. I;
c      * call again.
c
c      equal to 0, the sorting is done.
c
c    Output, integer I, J, the indices of two items.
c    On return with INDX positive, elements I and J should be interchanged.
c    On return with INDX negative, elements I and J should be compared, and
c    the result reported in ISGN on the next call.
c
c    Input, integer ISGN, results of comparison of elements I and J.
c    (Used only when the previous call returned INDX less than 0).
c    ISGN .le. 0 means I is less than or equal to J;
c    0 .le. ISGN means I is greater than or equal to J.
c
      implicit none

      integer i
      integer i_save
      integer indx
      integer isgn
      integer j
      integer j_save
      integer k
      integer k1
      integer n
      integer n1

      save i_save
      save j_save
      save k
      save k1
      save n1

      data i_save / 0 /
      data j_save / 0 /
      data k / 0 /
      data k1 / 0 /
      data n1 / 0 /
c
c  INDX = 0: This is the first call.
c
      if ( indx .eq. 0 ) then

        i_save = 0
        j_save = 0
        k = n / 2
        k1 = k
        n1 = n
c
c  INDX .lt. 0: The user is returning the results of a comparison.
c
      else if ( indx .lt. 0 ) then

        if ( indx .eq. -2 ) then

          if ( isgn .lt. 0 ) then
            i_save = i_save + 1
          end if

          j_save = k1
          k1 = i_save
          indx = -1
          i = i_save
          j = j_save
          return

        end if

        if ( 0 .lt. isgn ) then
          indx = 2
          i = i_save
          j = j_save
          return
        end if

        if ( k .le. 1 ) then

          if ( n1 .eq. 1 ) then
            i_save = 0
            j_save = 0
            indx = 0
          else
            i_save = n1
            n1 = n1 - 1
            j_save = 1
            indx = 1
          end if

          i = i_save
          j = j_save
          return

        end if

        k = k - 1
        k1 = k
c
c  0 .lt. INDX, the user was asked to make an interchange.
c
      else if ( indx .eq. 1 ) then

        k1 = k

      end if

10    continue

        i_save = 2 * k1

        if ( i_save .eq. n1 ) then
          j_save = k1
          k1 = i_save
          indx = -1
          i = i_save
          j = j_save
          return
        else if ( i_save .le. n1 ) then
          j_save = i_save + 1
          indx = -2
          i = i_save
          j = j_save
          return
        end if

        if ( k .le. 1 ) then
          go to 20
        end if

        k = k - 1
        k1 = k

      go to 10

20    continue

      if ( n1 .eq. 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
        i = i_save
        j = j_save
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
        i = i_save
        j = j_save
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
