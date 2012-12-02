      program band_qr_prb

c*********************************************************************72
c
cc BAND_QR_PRB tests routines from BAND_QR.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BAND_QR_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the BAND_QR library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BAND_QR_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests DGEQRF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 March 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer lda
      integer ldb
      integer lwork
      integer m
      integer n
      integer nrhs
      integer tau_size

      parameter ( m = 5 )
      parameter ( n = 5 )
      parameter ( nrhs = 1 )

      parameter ( lda = m )
      parameter ( ldb = max ( m, n ) )
      parameter ( lwork = n )
      parameter ( tau_size = n )

      double precision a(lda,n)
      double precision b(ldb,nrhs)
      integer i
      integer info
      integer j
      double precision tau(tau_size)
      double precision work(lwork)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Call the standard LAPACK routine'
      write ( *, '(a)' ) '  DGEQRF to get the QR factorization of'
      write ( *, '(a)' ) '  a matrix A stored in GE format, but which'
      write ( *, '(a)' ) '  is actually banded.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Then solve the linear system A*X=B using an'
      write ( *, '(a)' ) '  explicit version of DGEQRS_TWO.'
      write ( *, '(a)' ) '  DGEQRS_TWO is a version of '
      write ( *, '(a)' ) '  DGEQRS with no calls to subroutines.'

      do j = 1, n
        do i = 1, m
          if ( j == i - 1 ) then
            a(i,j) = - 1.0D+00
          else if ( j == i ) then
            a(i,j) = 2.0D+00
          else if ( j == i + 1 ) then
            a(i,j) = -1.0D+00
          else
            a(i,j) = 0.0D+00
          end if
        end do
      end do

      call r8ge_print ( m, n, a, '  Input matrix:' )
c
c  QR factor the matrix.
c
      call dgeqrf ( m, n, a, lda, tau, work, lwork, info )

      if ( info == 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  DGEQRF called successfully.'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  DGEQRF returned error flag.'
        write ( *, '(a,i8)' ) '  INFO = ', info
        return
      end if

      call r8ge_print ( m, n, a, '  Factored matrix:' )

      call r8vec_print ( tau_size, tau, '  Tau:' );
c
c  Set up and solve a linear system using DQEQRS_TWO.
c
      do i = 1, n - 1
        b(i,1) = 0.0D+00
      end do
      b(n,1) = dble ( n + 1 )

      call dgeqrs_two (
     &  m, n, nrhs, a, lda, tau, b, ldb, work, lwork, info )

      if ( info == 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  DGEQRS_TWO called successfully.'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  DGEQRS_TWO returned error flag.'
        write ( *, '(a,i8)' ) '  INFO = ', info
        return
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solution from DGEQRS_TWO:'
      write ( *, '(a)' ) ' '
      do i = 1, m
        write ( *, '(2x,g14.6)' ) b(i,1)
      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests DGEBQR2.
c
c  Discussion:
c
c    The same operation is being carried out as in Test #1, except
c    that DGEBQR2 knows that the matrix A is banded.  The matrix A
c    is still stored as a full storage "GE" matrix, though.
C
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer lda
      integer ldb
      integer lwork
      integer m
      integer ml
      integer mu
      integer n
      integer nrhs
      integer tau_size

      parameter ( m = 5 )
      parameter ( ml = 1 )
      parameter ( mu = 1 )
      parameter ( n = 5 )
      parameter ( nrhs = 1 )

      parameter ( lda = m )
      parameter ( ldb = max ( m, n ) )
      parameter ( tau_size = n )
      parameter ( lwork = min ( ml + mu, n ) )

      double precision a(lda,n)
      double precision b(ldb,nrhs)
      integer i
      integer info
      integer j
      double precision tau(tau_size)
      double precision work(lwork)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Call the LAPACK-like routine'
      write ( *, '(a)' ) '  DGEBQR2 to get the QR factorization of'
      write ( *, '(a)' ) '  a matrix A stored in GE format, but which'
      write ( *, '(a)' ) '  is actually banded.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Then solve the linear system A*X=B using an'
      write ( *, '(a)' ) '  explicit version of DGEQRS_TWO.'
      write ( *, '(a)' ) '  DGEQRS_TWO is a version of '
      write ( *, '(a)' ) '  DGEQRS with no calls to subroutines.'

      do j = 1, n
        do i = 1, m
          if ( j == i - 1 ) then
            a(i,j) = - 1.0D+00
          else if ( j == i ) then
            a(i,j) = 2.0D+00
          else if ( j == i + 1 ) then
            a(i,j) = -1.0D+00
          else
            a(i,j) = 0.0D+00
          end if
        end do
      end do

      call r8ge_print ( m, n, a, '  Input matrix:' )
c
c  QR factor the matrix.
c
      call dgebqr2 ( m, n, ml, mu, a, lda, tau, work, info )

      if ( info == 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  DGEBQR2 called successfully.'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  DGEBQR2 returned error flag.'
        write ( *, '(a,i8)' ) '  INFO = ', info
        return
      end if

      call r8ge_print ( m, n, a, '  Factored matrix:' )

      call r8vec_print ( tau_size, tau, '  Tau:' );
c
c  Set up and solve a linear system using DQEQRS_TWO.
c
      do i = 1, n - 1
        b(i,1) = 0.0D+00
      end do
      b(n,1) = dble ( n + 1 )

      call dgeqrs_two (
     &  m, n, nrhs, a, lda, tau, b, ldb, work, lwork, info )

      if ( info == 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  DGEQRS_TWO called successfully.'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  DGEQRS_TWO returned error flag.'
        write ( *, '(a,i8)' ) '  INFO = ', info
        return
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solution from DGEQRS_TWO:'
      write ( *, '(a)' ) ' '
      do i = 1, m
        write ( *, '(2x,g14.6)' ) b(i,1)
      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests DGBBQR2.
c
c  Discussion:
c
c    The same calculation is carried out as in Test #2.  However,
c    DGBBQR2 knows that the matrix A is banded, and that it is stored
c    in the "GB" format.  The QR factorization data overwrites the
c    values of A, and fits in the GB format.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 April 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer lda
      integer ldb
      integer lwork
      integer m
      integer ml
      integer mu
      integer n
      integer nrhs
      integer tau_size

      parameter ( m = 5 )
      parameter ( ml = 1 )
      parameter ( mu = 1 )
      parameter ( n = 5 )
      parameter ( nrhs = 1 )

      parameter ( lda = 2 * ml + mu + 1 )
      parameter ( ldb = max ( m, n ) )
      parameter ( tau_size = n )
      parameter ( lwork = min ( ml + mu, n ) )

      double precision a(lda,n)
      double precision b(ldb,nrhs)
      integer i
      integer info
      integer j
      double precision tau(tau_size)
      double precision work(lwork)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  Call the LAPACK-like routine'
      write ( *, '(a)' ) '  DGBBQR2 to get the QR factorization of'
      write ( *, '(a)' ) '  a matrix A stored in GB format.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Call DGBBQRS to solve the linear system'
      write ( *, '(a)' ) '  A*X=B.'
c
c  Zero out the matrix.
c
      do i = 1, lda
        do j = 1, n
          a(i,j) = 0.0D+00
        end do
      end do
c
c  Superdiagonal,
c  Diagonal,
c  Subdiagonal.
c
      do j = 2, n
        a(ml + mu + 1 - 1,j) = -1.0D+00
      end do

      do j = 1, n
        a(ml + mu + 1,j) = 2.0D+00
      end do

      do j = 1, n - 1
        a(ml + mu + 1 + 1,j) = -1.0D+00
      end do

      call r8gb_print ( m, n, ml, mu, a, '  Input matrix:' )
c
c  QR factor the matrix.
c
      call dgbbqr2 ( m, n, ml, mu, a, lda, tau, work, info )

      if ( info == 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  DGBBQR2 called successfully.'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  DGBBQR2 returned error flag.'
        write ( *, '(a,i8)' ) '  INFO = ', info
        return
      end if

      call r8gb_print ( m, n, ml, mu, a, '  Factored matrix:' )

      call r8vec_print ( tau_size, tau, '  Tau:' );
c
c  Set up and solve a linear system using DQBBQRS.
c
      do i = 1, n - 1
        b(i,1) = 0.0D+00
      end do
      b(n,1) = dble ( n + 1 )

      call dgbbqrs ( m, n, ml, mu, nrhs, a, lda, tau, b, ldb, work,
     &  lwork, info )

      if ( info == 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  DGBBQRS called successfully.'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  DGBBQRS returned error flag.'
        write ( *, '(a,i8)' ) '  INFO = ', info
        return
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solution from DGBBQRS:'
      write ( *, '(a)' ) ' '
      do i = 1, m
        write ( *, '(2x,g14.6)' ) b(i,1)
      end do

      return
      end

