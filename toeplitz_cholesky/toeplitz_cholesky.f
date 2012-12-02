      subroutine r8mat_mmt ( n1, n2, n3, a, b, c )

c*********************************************************************72
c
cc R8MAT_MMT multiplies computes C = A * B' for two R8MAT's.
c
c  Discussion:
c
c    An R8MAT is an array of R8 values.
c
c    In FORTRAN90, this operation is more efficiently done by the
c    command:
c
c      C(1:N1,1:N3) = matmul ( A(1:N1,1;N2), transpose ( B(1:N3,1:N2) ) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, N3, the order of the matrices.
c
c    Input, double precision A(N1,N2), B(N3,N2), the matrices to multiply.
c
c    Output, double precision C(N1,N3), the product matrix C = A * B'.
c
      implicit none

      integer n1
      integer n2
      integer n3

      double precision a(n1,n2)
      double precision b(n3,n2)
      double precision c(n1,n3)
      integer i
      integer j
      integer k

      do i = 1, n1
        do j = 1, n3
          c(i,j) = 0.0D+00
          do k = 1, n2
            c(i,j) = c(i,j) + a(i,k) * b(j,k)
          end do
        end do
      end do

      return
      end
      subroutine r8mat_mtm ( n1, n2, n3, a, b, c )

c*********************************************************************72
c
cc R8MAT_MTM multiplies computes C = A' * B for two R8MAT's.
c
c  Discussion:
c
c    An R8MAT is an array of R8 values.
c
c    In FORTRAN90, this operation is more efficiently done by the
c    command:
c
c      C(1:N1,1:N3) = matmul ( transpose ( A(1:N2,1;N1) ), B(1:N2,1:N3) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, N3, the order of the matrices.
c
c    Input, double precision A(N2,N1), B(N2,N3), the matrices to multiply.
c
c    Output, double precision C(N1,N3), the product matrix C = A' * B.
c
      implicit none

      integer n1
      integer n2
      integer n3

      double precision a(n2,n1)
      double precision b(n2,n3)
      double precision c(n1,n3)
      integer i
      integer j
      integer k

      do i = 1, n1
        do j = 1, n3
          c(i,j) = 0.0D+00
          do k = 1, n2
            c(i,j) = c(i,j) + a(k,i) * b(k,j)
          end do
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
      subroutine toep_cholesky_lower ( n, g, l )

c*********************************************************************72
c
cc TOEP_CHOLESKY_LOWER: lower Cholesky factor of a Toeplitz matrix.
c
c  Discussion:
c
c    The Toeplitz matrix A is supplied in a compressed form G.
c
c    The Toeplitz matrix must be positive semi-definite.
c
c    After factorization, A = L * L'.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Michael Stewart,
c    Cholesky factorization of semi-definite Toeplitz matrices.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision G(2,N), the compressed Toeplitz matrix.
c    G(1,1:N) contains the first row.
c    G(2,2:N) contains the first column.
c
c    Output, double precision L(N,N), the lower Cholesky factor.
c
      implicit none

      integer n

      double precision div
      double precision g(2,n)
      double precision g1j
      double precision g2j
      integer i
      integer j
      double precision l(n,n)
      double precision rho

      do j = 1, n
        do i = 1, n
          l(i,j) = 0.0D+00
        end do
      end do

      do i = 1, n
        l(i,1) = g(1,i)
      end do

      do j = n, 2, -1
        g(1,j) = g(1,j-1)
      end do

      g(1,1) = 0.0D+00

      do i = 2, n

        rho = - g(2,i) / g(1,i)
        div = sqrt ( ( 1.0D+00 - rho ) * ( 1.0D+00 + rho ) )

        do j = i, n
          g1j = g(1,j)
          g2j = g(2,j)
          g(1,j) = (       g1j + rho * g2j ) / div
          g(2,j) = ( rho * g1j +       g2j ) / div
        end do

        do j = i, n
          l(j,i) = g(1,j)
        end do

        do j = n, i + 1, -1
          g(1,j) = g(1,j-1)
        end do

        g(1,i) = 0.0D+00

      end do

      return
      end
      subroutine toep_cholesky_upper ( n, g, r )

c*********************************************************************72
c
cc TOEP_CHOLESKY_UPPER: upper Cholesky factor of a Toeplitz matrix.
c
c  Discussion:
c
c    The Toeplitz matrix A is supplied in a compressed form G.
c
c    The Toeplitz matrix must be positive semi-definite.
c
c    After factorization, A = R' * R.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Michael Stewart,
c    Cholesky factorization of semi-definite Toeplitz matrices.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision G(2,N), the compressed Toeplitz matrix.
c    G(1,1:N) contains the first row.
c    G(2,2:N) contains the first column.
c
c    Output, double precision R(N,N), the upper Cholesky factor.
c
      implicit none

      integer n

      double precision div
      double precision g(2,n)
      double precision g1j
      double precision g2j
      integer i
      integer j
      double precision r(n,n)
      double precision rho

      do j = 1, n
        do i = 1, n
          r(i,j) = 0.0D+00
        end do
      end do

      do j = 1, n
        r(1,j) = g(1,j)
      end do

      do j = n, 2, -1
        g(1,j) = g(1,j-1)
      end do

      g(1,1) = 0.0D+00

      do i = 2, n

        rho = - g(2,i) / g(1,i)
        div = sqrt ( ( 1.0D+00 - rho ) * ( 1.0D+00 + rho ) )

        do j = i, n
          g1j = g(1,j)
          g2j = g(2,j)
          g(1,j) = (       g1j + rho * g2j ) / div
          g(2,j) = ( rho * g1j +       g2j ) / div
        end do

        do j = i, n
          r(i,j) = g(1,j)
        end do
        do j = n, i + 1, -1
          g(1,j) = g(1,j-1)
        end do
        g(1,i) = 0.0D+00
      end do

      return
      end
      subroutine toeplitz_cholesky_lower ( n, a, l )

c*********************************************************************72
c
cc TOEPLITZ_CHOLESKY_LOWER: lower Cholesky factor of a Toeplitz matrix.
c
c  Discussion:
c
c    The Toeplitz matrix must be positive semi-definite.
c
c    After factorization, A = L * L'.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Michael Stewart,
c    Cholesky factorization of semi-definite Toeplitz matrices.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision A(N,N), the Toeplitz matrix.
c
c    Output, double precision L(N,N), the lower Cholesky factor.
c
      implicit none

      integer n

      double precision a(n,n)
      double precision div
      double precision g(2,n)
      double precision g1j
      double precision g2j
      integer i
      integer j
      double precision l(n,n)
      double precision rho

      do j = 1, n
        do i = 1, n
          l(i,j) = 0.0D+00
        end do
      end do

      do j = 1, n
        g(1,j) = a(1,j)
      end do
      g(2,1) = 0.0D+00
      do j = 2, n
        g(2,j) = a(j,1)
      end do 

      do i = 1, n
        l(i,1) = g(1,i)
      end do
      do j = n, 2, -1
        g(1,j) = g(1,j-1)
      end do
      g(1,1) = 0.0D+00

      do i = 2, n

        rho = - g(2,i) / g(1,i)
        div = sqrt ( ( 1.0D+00 - rho ) * ( 1.0D+00 + rho ) )

        do j = i, n
          g1j = g(1,j)
          g2j = g(2,j)
          g(1,j) = (       g1j + rho * g2j ) / div
          g(2,j) = ( rho * g1j +       g2j ) / div
        end do

        do j = i, n
          l(j,i) = g(1,j)
        end do
        do j = n, i + 1, -1
          g(1,j) = g(1,j-1)
        end do
        g(1,i) = 0.0D+00

      end do

      return
      end
      subroutine toeplitz_cholesky_upper ( n, a, r )

c*********************************************************************72
c
cc TOEPLITZ_CHOLESKY_UPPER: upper Cholesky factor of a Toeplitz matrix.
c
c  Discussion:
c
c    The Toeplitz matrix must be positive semi-definite.
c
c    After factorization, A = R' * R.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Michael Stewart,
c    Cholesky factorization of semi-definite Toeplitz matrices.
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision A(N,N), the Toeplitz matrix.
c
c    Output, double precision R(N,N), the upper Cholesky factor.
c
      implicit none

      integer n

      double precision a(n,n)
      double precision div
      double precision g(2,n)
      double precision g1j
      double precision g2j
      integer i
      integer j
      double precision r(n,n)
      double precision rho

      do j = 1, n
        do i = 1, n
          r(i,j) = 0.0D+00
        end do
      end do

      do j = 1, n
        g(1,j) = a(1,j)
      end do

      g(2,1) = 0.0D+00
      do j = 2, n
        g(2,j) = a(j,1) 
      end do

      do j = 1, n
        r(1,j) = g(1,j)
      end do

      do j = n, 2, -1
        g(1,j) = g(1,j-1)
      end do
      g(1,1) = 0.0D+00

      do i = 2, n

        rho = - g(2,i) / g(1,i)
        div = sqrt ( ( 1.0D+00 - rho ) * ( 1.0D+00 + rho ) )
        do j = i, n
          g1j = g(1,j)
          g2j = g(2,j)
          g(1,j) = (       g1j + rho * g2j ) / div
          g(2,j) = ( rho * g1j +       g2j ) / div
        end do

        do j = i, n
          r(i,j) = g(1,j)
        end do

        do j = n, i + 1, -1
          g(1,j) = g(1,j-1)
        end do

        g(1,i) = 0.0D+00

      end do

      return
      end

