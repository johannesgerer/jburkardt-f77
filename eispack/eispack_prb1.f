      program main

c*********************************************************************72
c
cc  MAIN is the main program for EISPACK_PRB1.
c
c  Discussion:
c
c    EISPACK_PRB1 calls the EISPACK sample programs.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 November 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EISPACK_PRB1'
      write ( *, '(a)' ) '  FORTRAN77 version.'
      write ( *, '(a)' ) '  Test the EISPACK library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
      call test06 ( )
      call test065 ( )
      call test07 ( )
      call test08 ( )
      call test09 ( )
      call test10 ( )

      call test11 ( )
      call test12 ( )
      call test13 ( )
      call test14 ( )
      call test15 ( )
      call test16 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EISPACK_PRB1'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc  TEST01 tests CG.
c
c  Discussion:
c
c    CG is for the eigenvalues of a complex general matrix.
c
c    eigenvalues and eigenvectors of a complex general matrix
c    note that the eigenvalues of such a matrix are in general complex.
c    however, we will use the same example we used before, namely
c    a hermitian matrix, so the eigenvalues will in fact be real.
c
c    (3     1     0     0+2i)
c    (1     3     0-2i  0   )
c    (0     0+2i  1     1   )
c    (0-2i  0     1     1   )
c
c    The eigenvalues are 2+2*sqrt(2), 2-2*sqrt(2), 4 and 0
c
c    The eigenvector matrix is
c
c    (  1+sqrt(2),  1,                -1,          1)
c    (  1+sqrt(2),  1,                 1,         -1)
c    (     i,       -(1+sqrt(2))*i,    i,          i)
c    (    -i,        (1+sqrt(2))*i,    i,          i)
c
c    Note that the actual eigenvector matrix from EISPACK could
c    be scaled by a real value, or by i, and the columns may
c    appear in any order.
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
      implicit none

      integer n
      parameter ( n = 4 )

      double precision ai(n,n)
      double precision ar(n,n)
      double precision fv1(n)
      double precision fv2(n)
      double precision fv3(n)
      integer i
      integer ierr
      integer j
      integer matz
      double precision wi(n)
      double precision wr(n)
      double precision xi(n,n)
      double precision xr(n,n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  CG computes the eigenvalues and'
      write ( *, '(a)' ) '  eigenvectors of a complex general matrix.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix order = ', n
c
c  Set the values of the matrix.
c
      ar(1,1) = 3.0D+00
      ar(1,2) = 1.0D+00
      ar(1,3) = 0.0D+00
      ar(1,4) = 0.0D+00

      ar(2,1) = 1.0D+00
      ar(2,2) = 3.0D+00
      ar(2,3) = 0.0D+00
      ar(2,4) = 0.0D+00

      ar(3,1) = 0.0D+00
      ar(3,2) = 0.0D+00
      ar(3,3) = 1.0D+00
      ar(3,4) = 1.0D+00

      ar(4,1) = 0.0D+00
      ar(4,2) = 0.0D+00
      ar(4,3) = 1.0D+00
      ar(4,4) = 1.0D+00

      ai(1,1) = 0.0D+00
      ai(1,2) = 0.0D+00
      ai(1,3) = 0.0D+00
      ai(1,4) = 2.0D+00

      ai(2,1) = 0.0D+00
      ai(2,2) = 0.0D+00
      ai(2,3) = -2.0D+00
      ai(2,4) = 0.0D+00

      ai(3,1) = 0.0D+00
      ai(3,2) = 2.0D+00
      ai(3,3) = 0.0D+00
      ai(3,4) = 0.0D+00

      ai(4,1) = -2.0D+00
      ai(4,2) = -0.0D+00
      ai(4,3) = -0.0D+00
      ai(4,4) = 0.0D+00
c
c  matz = 0 for eigenvalues only,
c  matz = 1 for eigenvalues and eigenvectors.
c
      matz = 1
      call cg ( n, n, ar, ai, wr, wi, matz, xr, xi, fv1, fv2, fv3,
     &  ierr )

      if ( ierr .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Warning!'
        write ( *, '(a,i8)' ) '  The error return flag IERR = ',ierr
      end if

      call r8vec2_print ( n, wr, wi,
     &  '  Real and imaginary parts of eigenvalues:' )

      if ( matz .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  The eigenvectors are:'
        do i = 1, n
          write ( *, '(a)' ) ' '
          write ( *, '(a,i8)' ) '  Eigenvector ', i
          write ( *, '(a)' ) ' '
          do j = 1, n
            write ( *, '(2g14.6)' ) xr(i,j), xi(i,j)
          end do
        end do
      end if

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc  TEST02 tests CH.
c
c  Discussion:
c
c    CH is for the eigenvalues of a complex hermitian matrix.
c
c    Eigenvalues and eigenvectors of a complex hermitian matrix
c
c    Note that the eigenvalues (though not the eigenvectors) of
c    a hermitian matrix are real.
c
c    (3     1     0     0+2i)
c    (1     3     0-2i  0   )
c    (0     0+2i  1     1   )
c    (0-2i  0     1     1   )
c
c    The eigenvalues are 2+2*sqrt(2), 2-2*sqrt(2), 4 and 0
c
c    The eigenvector matrix is
c
c    (  1+sqrt(2),  1,                -1,          1)
c    (  1+sqrt(2),  1,                 1,         -1)
c    (     i,       -(1+sqrt(2))*i,    i,          i)
c    (    -i,        (1+sqrt(2))*i,    i,          i)
c
c    Note that the actual eigenvector matrix from EISPACK could
c    be scaled by a real value, or by i, and the columns may
c    appear in any order.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      double precision ar(n,n)
      double precision ai(n,n)
      double precision fm1(2,n)
      double precision fv1(n)
      double precision fv2(n)
      integer i
      integer ierr
      integer j
      integer matz
      double precision w(n)
      double precision xr(n,n)
      double precision xi(n,n)
c
c  Set the values of the matrix.
c
      ar(1,1) = 3.0D+00
      ar(1,2) = 1.0D+00
      ar(1,3) = 0.0D+00
      ar(1,4) = 0.0D+00

      ar(2,1) = 1.0D+00
      ar(2,2) = 3.0D+00
      ar(2,3) = 0.0D+00
      ar(2,4) = 0.0D+00

      ar(3,1) = 0.0D+00
      ar(3,2) = 0.0D+00
      ar(3,3) = 1.0D+00
      ar(3,4) = 1.0D+00

      ar(4,1) = 0.0D+00
      ar(4,2) = 0.0D+00
      ar(4,3) = 1.0D+00
      ar(4,4) = 1.0D+00

      ai(1,1) = 0.0D+00
      ai(1,2) = 0.0D+00
      ai(1,3) = 0.0D+00
      ai(1,4) = 2.0D+00

      ai(2,1) = 0.0D+00
      ai(2,2) = 0.0D+00
      ai(2,3) = -2.0D+00
      ai(2,4) = 0.0D+00

      ai(3,1) = 0.0D+00
      ai(3,2) = 2.0D+00
      ai(3,3) = 0.0D+00
      ai(3,4) = 0.0D+00

      ai(4,1) = -2.0D+00
      ai(4,2) = -0.0D+00
      ai(4,3) = -0.0D+00
      ai(4,4) = 0.0D+00

      matz = 1

      call ch ( n, n, ar, ai, w, matz, xr, xi, fv1, fv2, fm1, ierr )

      if ( ierr .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Warning!'
        write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  CH computes the eigenvalues and'
      write ( *, '(a)' ) '  eigenvectors of a complex hermitian matrix.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix order = ', n
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Error flag = ', ierr
      write ( *, '(a)' ) ' '

      call r8vec_print ( n, w, '  The eigenvalues Lambda:' )

      write ( *, '(a)' ) ' '

      if ( matz .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Eigenvectors are:'
        do i = 1, n
          write ( *, '(a)' ) ' '
          write ( *, '(a,i8)' ) '  Eigenvector ', i
          write ( *, '(a)' ) ' '
          do j = 1, n
            write ( *, '(2g14.6)' ) xr(i,j), xi(i,j)
          end do
        end do
      end if

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests MINFIT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 5 )
      integer nb
      parameter ( nb = 1 )
      integer n
      parameter ( n = 2 )

      double precision a(m,n)
      double precision acopy(m,n)
      double precision b(m,nb)
      integer i
      integer ierr
      integer j
      double precision r(m)
      double precision rv1(n)
      double precision w(n)
      double precision x(n)

      a(1,1) =   1.00D+00
      a(2,1) =   2.05D+00
      a(3,1) =   3.06D+00
      a(4,1) = - 1.02D+00
      a(5,1) =   4.08D+00

      a(1,2) =   1.00D+00
      a(2,2) = - 1.00D+00
      a(3,2) =   1.00D+00
      a(4,2) =   2.00D+00
      a(5,2) = - 1.00D+00

      do j = 1, n
        do i = 1, m
          acopy(i,j) = a(i,j)
        end do
      end do

      b(1,1) = 1.98D+00
      b(2,1) = 0.95D+00
      b(3,1) = 3.98D+00
      b(4,1) = 0.92D+00
      b(5,1) = 2.90D+00

      do i = 1, m
        r(i) = - b(i,1)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' )
     &  '  MINFIT solves an overdetermined linear system'
      write ( *, '(a)' ) '  using least squares methods.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix rows = ', m
      write ( *, '(a,i8)' ) '  Matrix columns = ', n

      call r8mat_print ( m, n, a, '  The matrix A:' )

      call r8mat_print ( m, nb, b, '  The right hand side B:' )

      call minfit ( m, m, n, a, w, nb, b, ierr, rv1 )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  MINFIT error code IERR = ', ierr

      call r8vec_print ( n, w, '  The singular values:' )
c
c  B now contains U' * B.
c  We need to divide by the singular values, and multiply by V.
c
      do i = 1, n
        b(i,1) = b(i,1) / w(i)
      end do

      do i = 1, n
        x(i) = 0.0D+00
        do j = 1, n
          x(i) = x(i) + a(i,j) * b(j,1)
        end do
      end do

      call r8vec_print ( n, x, '  The least squares solution X:' )

      do i = 1, m
        do j = 1, n
          r(i) = r(i) + acopy(i,j) * x(j)
        end do
      end do

      call r8vec_print ( m, r, '  The residual A * X - B:' )

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests RG.
c
c  Discussion:
c
c    RG is for the eigenvalues of a general real matrix.
c
c    The matrix A is nonsymmetric.  The eigenvalues may therefore be
c    complex numbers.
c
c    ( 33  16  72)
c    (-24 -10 -57)
c    ( -8  -4 -17)
c
c    The eigenvalues of A are (1,2,3)
c
c    The eigenvectors of A are
c
c    (-15 -16 -4)
c    ( 12  13  3)
c    (  4   4  1)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      double precision a(n,n)
      double precision acopy(n,n)
      double precision fv1(n)
      integer i
      integer ierr
      integer iv1(n)
      integer j
      integer k
      integer matz
      double precision sum3
      double precision sum1
      double precision sum2
      double precision wi(n)
      double precision wr(n)
      double precision x(n,n)
c
c  Set the values of the matrix.
c
      a(1,1) = 33.0D+00
      a(1,2) = 16.0D+00
      a(1,3) = 72.0D+00

      a(2,1) = -24.0D+00
      a(2,2) = -10.0D+00
      a(2,3) = -57.0D+00

      a(3,1) = -8.0D+00
      a(3,2) = -4.0D+00
      a(3,3) = -17.0D+00
c
c  RG overwrites A with garbage, so save a copy nowc
c
      do j = 1, n
        do i = 1, n
          acopy(i,j) = a(i,j)
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' )
     &  '  RG computes the eigenvalues and eigenvectors of'
      write ( *, '(a)' ) '  a real general matrix.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix order = ', n

      call r8mat_print ( n, n, a, '  The matrix A:' )

      matz = 1

      call rg ( n, n, a, wr, wi, matz, x, iv1, fv1, ierr )

      if ( ierr .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST04 - Warning!'
        write ( *, '(a,i8)' )
     &  '  The error return flag was IERR = ', ierr
      end if

      call r8vec2_print ( n, wr, wi,
     &  '  Real and imaginary parts of eigenvalues:' )

      if ( matz .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  The eigenvectors may be complex:'
        do j = 1, n
          write ( *, '(a)' ) ' '
          write ( *, '(a,i8)' ) '  Eigenvector ', j
          write ( *, '(a)' ) ' '
          do i = 1, n
            if ( wi(j) .eq. 0.0D+00 ) then
              write ( *, '(g14.6)' ) x(i,j)
            else if ( 0.0D+00 .lt. wi(j) ) then
              write ( *, '(2g14.6)' ) x(i,j), x(i,j+1)
            else if ( wi(j) .lt. 0.0D+00 ) then
              write ( *, '(2g14.6)' ) x(i,j-1), -x(i,j)
            end if
          end do
        end do
c
c  Check.
c  First, restore the original values of A.
c
        do j = 1, n
          do i = 1, n
            a(i,j) = acopy(i,j)
          end do
        end do

        do k = 1, n

          write ( *, '(a)' ) ' '
          write ( *, '(a,i8)' )
     &  '  Residuals (A*x-Lambda*x) for eigenvalue ', k
          write ( *, '(a)' ) ' '

          if ( wi(k) .eq. 0.0D+00 ) then

            do i = 1, n
              sum3 = 0.0D+00
              do j = 1, n
                sum3 = sum3 + a(i,j) * x(j,k)
              end do
              sum3 = sum3 - wr(k) * x(i,k)
              write ( *, '(g14.6)' ) sum3
            end do

          else if ( 0.0D+00 .lt. wi(k) ) then

            do i = 1, n
              sum1 = 0.0D+00
              sum2 = 0.0D+00
              do j = 1, n
                sum1 = sum1 + a(i,j) * x(j,k)
                sum2 = sum2 + a(i,j) * x(j,k+1)
              end do
              sum1 = sum1 - wr(k) * x(i,k) + wi(k) * x(i,k+1)
              sum2 = sum2 - wi(k) * x(i,k) - wr(k) * x(i,k+1)
              write ( *, '(2g14.6)' ) sum1, sum2
            end do

          else if ( wi(k) .lt. 0.0D+00 ) then

            do i = 1, n
              sum1 = 0.0D+00
              sum2 = 0.0D+00
              do j = 1, n
                sum1 = sum1 + a(i,j) * x(j,k-1)
                sum2 = sum2 - a(i,j) * x(j,k)
              end do
              sum1 = sum1 - wr(k) * x(i,k-1) - wi(k) * x(i,k)
              sum2 = sum2 - wi(k) * x(i,k-1) + wr(k) * x(i,k)
              write ( *, '(2g14.6)' ) sum1, sum2
            end do

          end if

        end do

      end if

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests RGG.
c
c  Discussion:
c
c    RGG is for a real generalized general eigenvalue problem.
c
c    A generalized eigenvalue problem.  Given matrices A and B, find
c    N numbers LAMBDA, and for each LAMBDA a vector X, so that
c
c      A*x = lambda*B*x
c
c    The matrix A is
c
c    ( -7 7  6  6)
c    (-10 8 10  8)
c    ( -8 3 10 11)
c    ( -4 0  4 12)
c
c    The matrix B is
c
c    (2 1 0 0)
c    (1 2 1 0)
c    (0 1 2 1)
c    (0 0 1 2)
c
c    The correct eigenvalues LAMBDA are
c
c    (1,2,3,4)
c
c    The correct eigenvectors X are
c
c    (4 3 2 1)
c    (3 3 2 1)
c    (2 2 2 1)
c    (1 1 1 1)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 January 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      double precision a(n,n)
      double precision acopy(n,n)
      double precision alfi(n)
      double precision alfr(n)
      double precision b(n,n)
      double precision bcopy(n,n)
      double precision beta(n)
      integer i
      integer ierr
      integer j
      integer k
      integer matz
      double precision sum3
      double precision sum1
      double precision sum2
      double precision x(n,n)
c
c  Set the values in the A matrix.
c
      a(1,1) = -7.0D+00
      a(1,2) = 7.0D+00
      a(1,3) = 6.0D+00
      a(1,4) = 6.0D+00

      a(2,1) = -10.0D+00
      a(2,2) = 8.0D+00
      a(2,3) = 10.0D+00
      a(2,4) = 8.0D+00

      a(3,1) = -8.0D+00
      a(3,2) = 3.0D+00
      a(3,3) = 10.0D+00
      a(3,4) = 11.0D+00

      a(4,1) = -4.0D+00
      a(4,2) = 0.0D+00
      a(4,3) = 4.0D+00
      a(4,4) = 12.0D+00
c
c  Save a copy of A.
c
      do j = 1, n
        do i = 1, n
          acopy(i,j) = a(i,j)
        end do
      end do
c
c  Set the values in the B matrix.
c
      b(1,1) = 2.0D+00
      b(1,2) = 1.0D+00
      b(1,3) = 0.0D+00
      b(1,4) = 0.0D+00

      b(2,1) = 1.0D+00
      b(2,2) = 2.0D+00
      b(2,3) = 1.0D+00
      b(2,4) = 0.0D+00

      b(3,1) = 0.0D+00
      b(3,2) = 1.0D+00
      b(3,3) = 2.0D+00
      b(3,4) = 1.0D+00

      b(4,1) = 0.0D+00
      b(4,2) = 0.0D+00
      b(4,3) = 1.0D+00
      b(4,4) = 2.0D+00
c
c  Save a copy of B.
c
      do j = 1, n
        do i = 1, n
          bcopy(i,j) = b(i,j)
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05:'
      write ( *, '(a)' ) '  RGG for real generalized problem.'
      write ( *, '(a)' ) '  Find scalars LAMBDA and vectors X so that'
      write ( *, '(a)' ) '    A*X = LAMBDA * B * X'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix order = ', n

      call r8mat_print ( n, n, a, '  The matrix A:' )

      call r8mat_print ( n, n, b, '  The matrix B:' )

      matz = 1

      call rgg ( n, n, a, b, alfr, alfi, beta, matz, x, ierr )

      if ( ierr .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Warning!'
        write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
      end if

      do i = 1, n
        alfr(i) = alfr(i) / beta(i)
        alfi(i) = alfi(i) / beta(i)
      end do

      call r8vec2_print ( n, alfr, alfi,
     &  '  Real and imaginary parts of eigenvalues:' )

      if ( matz .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  The eigenvectors are:'
        do i = 1, n
          write ( *, '(a)' ) ' '
          write ( *, '(a,i8)' ) '  Eigenvector ', i
          write ( *, '(a)' ) ' '
          do j = 1, n
            write ( *, '(g14.6)' ) x(i,j)
          end do
        end do
      end if
c
c  Check.
c  First, restore the original values of A and B.
c
      if ( matz .ne. 0 ) then

        do j = 1, n
          do i = 1, n
            a(i,j) = acopy(i,j)
          end do
        end do

        do j = 1, n
          do i = 1, n
            b(i,j) = bcopy(i,j)
          end do
        end do

        do k = 1, n
          write ( *, '(a)' ) ' '
          write ( *, '(a,i8)' )
     &      '  Residuals (A*x-(Alfr+Alfi*I)*B*x) for eigenvalue ', k
          write ( *, '(a)' ) ' '

          if ( alfi(k) .eq. 0.0D+00 ) then

            do i = 1, n

              sum3 = 0.0D+00
              do j = 1, n
                sum3 = sum3 + a(i,j) * x(j,k)
              end do

              do j = 1, n
                sum3 = sum3 - alfr(k) * b(i,j) * x(j,k)
              end do

              write ( *, '(g14.6)' ) sum3
            end do

          else if ( 0.0D+00 .lt. alfi(k) ) then

            do i = 1, n

              sum1 = 0.0D+00
              sum2 = 0.0D+00
              do j = 1, n
                sum1 = sum1 + a(i,j) * x(j,k)
                sum2 = sum2 + a(i,j) * x(j,k+1)
              end do

              do j = 1, n
                sum1 = sum1 - alfr(k) * b(i,j) * x(j,k)
     &                      + alfi(k) * b(i,j) * x(j,k+1)

                sum2 = sum2 - alfi(k) * b(i,j) * x(j,k)
     &                      - alfr(k) * b(i,j) * x(j,k+1)

              end do

              write ( *, '(2g14.6)' ) sum1, sum2
            end do

          else if ( alfi(k) .lt. 0.0D+00 ) then

            do i = 1, n

              sum1 = 0.0D+00
              sum2 = 0.0D+00
              do j = 1, n
                sum1 = sum1 + a(i,j) * x(j,k-1)
                sum2 = sum2 - a(i,j) * x(j,k)
              end do

              do j = 1, n

                sum1 = sum1 - alfr(k) * b(i,j) * x(j,k-1)
     &                      - alfi(k) * b(i,j) * x(j,k)

                sum2 = sum2 - alfi(k) * b(i,j) * x(j,k-1)
     &                      + alfr(k) * b(i,j) * x(j,k)
              end do

              write ( *, '(2g14.6)' ) sum1, sum2
            end do

          end if

        end do

      end if

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 tests RS.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 January 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      double precision a(n,n)
      double precision a2(n,n)
      double precision fv1(n)
      double precision fv2(n)
      integer i
      integer ierr
      integer j
      integer k
      integer matz
      double precision r(n,n)
      double precision w(n)
      double precision x(n,n)
c
c  Set the values in the matrix.
c
      a(1,1) = 5.0D+00
      a(1,2) = 4.0D+00
      a(1,3) = 1.0D+00
      a(1,4) = 1.0D+00

      a(2,1) = 4.0D+00
      a(2,2) = 5.0D+00
      a(2,3) = 1.0D+00
      a(2,4) = 1.0D+00

      a(3,1) = 1.0D+00
      a(3,2) = 1.0D+00
      a(3,3) = 4.0D+00
      a(3,4) = 2.0D+00

      a(4,1) = 1.0D+00
      a(4,2) = 1.0D+00
      a(4,3) = 2.0D+00
      a(4,4) = 4.0D+00
c
c  Save a copy of the matrix.
c
      do j = 1, n
        do i = 1, n
          a2(i,j) = a(i,j)
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' )
     &  '  RS computes the eigenvalues and eigenvectors'
      write ( *, '(a)' ) '  of a real symmetric matrix.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix order = ', n

      call r8mat_print ( n, n, a, '  The matrix A:' )

      matz = 1

      call rs ( n, n, a, w, matz, x, fv1, fv2, ierr )

      if ( ierr .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Warning!'
        write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
      end if

      call r8vec_print ( n, w, '  The eigenvalues Lambda:' )

      if ( matz .ne. 0 ) then

        call r8mat_print ( n, n, x, '  The eigenvector matrix:' )

        do i = 1, n
          do j = 1, n
            r(i,j) = 0.0D+00
            do k = 1, n
              r(i,j) = r(i,j) + a2(i,k) * x(k,j)
            end do
          end do
        end do

        do j = 1, n
          do i = 1, n
            r(i,j) = r(i,j) - w(j) * x(i,j)
          end do
        end do

        call r8mat_print ( n, n, r, '  The residual (A-Lambda*I)*X:' )

      end if

      return
      end
      subroutine test065 ( )

c*********************************************************************72
c
cc TEST065 tests RS.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 November 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      double precision a(n,n)
      double precision a2(n,n)
      double precision fv1(n)
      double precision fv2(n)
      integer i
      integer ierr
      integer j
      integer k
      integer matz
      double precision r(n,n)
      integer seed
      double precision t
      double precision w(n)
      double precision x(n,n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST065'
      write ( *, '(a)' ) '  RS computes the eigenvalues and'
      write ( *, '(a)' ) '  eigenvectors of a real symmetric matrix.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix order = ', n

      seed = 123456789

      call r8mat_uniform_01 ( n, n, seed, a )

      do i = 1, n - 1
        do j = i + 1, n
          t = ( a(i,j) + a(j,i) ) / 2.0D+00
          a(i,j) = t
          a(j,i) = t
        end do
      end do
c
c  Save a copy of the matrix.
c
      do j = 1, n
        do i = 1, n
          a2(i,j) = a(i,j)
        end do
      end do

      call r8mat_print ( n, n, a, '  The matrix A:' )

      matz = 1

      call rs ( n, n, a, w, matz, x, fv1, fv2, ierr )

      if ( ierr .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Warning!'
        write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
      end if

      call r8vec_print ( n, w, '  The eigenvalues Lambda:' )

      if ( matz .ne. 0 ) then

        call r8mat_print ( n, n, x, '  The eigenvector matrix:' )

        do i = 1, n
          do j = 1, n
            r(i,j) = 0.0D+00
            do k = 1, n
              r(i,j) = r(i,j) + a2(i,k) * x(k,j)
            end do
          end do
        end do

        do j = 1, n
          do i = 1, n
            r(i,j) = r(i,j) - w(j) * x(i,j)
          end do
        end do

        call r8mat_print ( n, n, r, '  The residual (A-Lambda*I)*X:' )

      end if

      return
      end
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07 tests RSB.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )
      integer mb
      parameter ( mb = 2 )

      double precision a(n,mb)
      double precision a2(n,n)
      double precision fv1(n)
      double precision fv2(n)
      integer i
      integer ierr
      integer j
      integer k
      integer matz
      double precision r(n,n)
      double precision w(n)
      double precision x(n,n)

      do i = 1, n
        do j = 1, mb
          a(i,j) = 0.0D+00
        end do
      end do

      do i = 1, n
        a(i,mb) = 2.0D+00
      end do

      do i = 2, n
        a(i,1) = -1.0D+00
      end do

      do i = 1, n
        do j = 1, n
          if ( i .eq. j ) then
            a2(i,j) = 2.0D+00
          else if ( abs ( i - j ) .eq. 1 ) then
            a2(i,j) = - 1.0D+00
          else
            a2(i,j) = 0.0D+00
          end if
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07'
      write ( *, '(a)' )
     &  '  RSB computes the eigenvalues and eigenvectors'
      write ( *, '(a)' ) '  of a real symmetric band matrix.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix order = ', n

      call r8mat_print ( n, n, a2, '  The matrix A:' )

      matz = 1

      call rsb ( n, n, mb, a, w, matz, x, fv1, fv2, ierr )

      if ( ierr .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Warning!'
        write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
      end if

      call r8vec_print ( n, w, '  The eigenvalues Lambda:' )

      if ( matz .ne. 0 ) then

        call r8mat_print ( n, n, x, '  The eigenvector matrix X:' )

        do i = 1, n
          do j = 1, n
            r(i,j) = 0.0D+00
            do k = 1, n
              r(i,j) = r(i,j) + a2(i,k) * x(k,j)
            end do
          end do
        end do

        do j = 1, n
          do i = 1, n
            r(i,j) = r(i,j) - w(j) * x(i,j)
          end do
        end do

        call r8mat_print ( n, n, r, '  The residual (A-Lambda*I)*X:' )

      end if

      return
      end
      subroutine test08 ( )

c*********************************************************************72
c
cc TEST08 tests RSG.
c
c  Discussion:
c
c    RGG is for a real generalized eigenvalue problem of the form
c
c      A*x = lambda*B*x
c
c    with A symmetric and B positive definite symmetric.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      double precision a(n,n)
      double precision a2(n,n)
      double precision b(n,n)
      double precision b2(n,n)
      double precision fv1(n)
      double precision fv2(n)
      integer i
      integer ierr
      integer j
      integer k
      integer matz
      double precision sum3
      double precision sum1
      double precision sum2
      double precision w(n)
      double precision x(n,n)

      do i = 1, n
        do j = 1, n
          a(i,j) = abs ( i - j )
        end do
      end do

      do i = 1, n
        do j = 1, n
          a2(i,j) = a(i,j)
        end do
      end do

      do i = 1, n
        do j = 1, n
          if ( i .eq. j ) then
            b(i,j) = 2.0D+00
          else if ( abs ( i - j ) .eq. 1 ) then
            b(i,j) = - 1.0D+00
          else
            b(i,j) = 0.0D+00
          end if
        end do
      end do

      do j = 1, n
        do i = 1, n
          b2(i,j) = b(i,j)
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08:'
      write ( *, '(a)' ) '  RSG for real symmetric generalized problem.'
      write ( *, '(a)' ) '  Find scalars LAMBDA and vectors X so that'
      write ( *, '(a)' ) '    A*X = LAMBDA * B * X'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix order = ', n

      call r8mat_print ( n, n, a, '  The matrix A:' )

      call r8mat_print ( n, n, b, '  The matrix B:' )

      matz = 1

      call rsg ( n, n, a, b, w, matz, x, fv1, fv2, ierr )

      if ( ierr .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Warning!'
        write ( *, '(a,i8)' ) '  The error return flag IERR = ',ierr
      end if

      call r8vec_print ( n, w, '  The eigenvalues Lambda:' )

      if ( matz .ne. 0 ) then

        call r8mat_print ( n, n, x, '  The eigenvector matrix X:' )

        do j = 1, n
          do i = 1, n
            a(i,j) = a2(i,j)
          end do
        end do

        do j = 1, n
          do i = 1, n
            b(i,j) = b2(i,j)
          end do
        end do

        do k = 1, n

          write ( *, '(a)' ) ' '
          write ( *, '(a,i8)' )
     &      'Residuals (A*x-(w*I)*B*x) for eigenvalue ', k
          write ( *, '(a)' ) ' '

            do i = 1, n

              sum3 = 0.0D+00
              do j = 1, n
                sum3 = sum3 + a(i,j) * x(j,k)
              end do

              do j = 1, n
                sum3 = sum3 - w(k) * b(i,j) * x(j,k)
              end do

              write ( *, '(g14.6)' ) sum3
            end do

        end do

      end if

      return
      end
      subroutine test09 ( )

c*********************************************************************72
c
cc TEST09 tests RSGAB.
c
c  Discussion:
c
c    RGGAB is for a real generalized eigenvalue problem of the form
c
c      A*B*x = lambda*x
c
c    with A symmetric and B positive definite symmetric.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      double precision a(n,n)
      double precision a2(n,n)
      double precision b(n,n)
      double precision b2(n,n)
      double precision fv1(n)
      double precision fv2(n)
      integer i
      integer ierr
      integer j
      integer k
      integer matz
      double precision r(n,n)
      double precision w(n)
      double precision x(n,n)

      do i = 1, n
        do j = 1, n
          a(i,j) = abs ( i - j )
        end do
      end do

      do j = 1, n
        do i = 1, n
          a2(i,j) = a(i,j)
        end do
      end do

      do i = 1, n
        do j = 1, n
          if ( i .eq. j ) then
            b(i,j) = 2.0D+00
          else if ( abs ( i - j ) .eq. 1 ) then
            b(i,j) = - 1.0D+00
          else
            b(i,j) = 0.0D+00
          end if
        end do
      end do

      do j = 1, n
        do i = 1, n
          b2(i,j) = b(i,j)
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST09:'
      write ( *, '(a)' )
     &  '  RSGAB for real symmetric generalized problem.'
      write ( *, '(a)' ) '  Find scalars LAMBDA and vectors X so that'
      write ( *, '(a)' ) '    A*B*X = LAMBDA * X'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix order = ', n

      call r8mat_print ( n, n, a, '  The matrix A:' )

      call r8mat_print ( n, n, b, '  The matrix B:' )

      matz = 1

      call rsgab ( n, n, a, b, w, matz, x, fv1, fv2, ierr )

      if ( ierr .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Warning!'
        write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
      end if

      call r8vec_print ( n, w, '  The eigenvalues Lambda:' )

      if ( matz .ne. 0 ) then

        call r8mat_print ( n, n, x, '  The eigenvector matrix X:' )

        do i = 1, n
          do j = 1, n
            r(i,j) = 0.0D+00
            do k = 1, n
              r(i,j) = r(i,j) + b2(i,k) * x(k,j)
            end do
          end do
        end do

        do i = 1, n
          do j = 1, n
            b2(i,j) = r(i,j)
          end do
        end do

        do i = 1, n
          do j = 1, n
            r(i,j) = 0.0D+00
            do k = 1, n
              r(i,j) = r(i,j) + a2(i,k) * b2(k,j)
            end do
          end do
        end do

        do j = 1, n
          do i = 1, n
            r(i,j) = r(i,j) - w(j) * x(i,j)
          end do
        end do

        call r8mat_print ( n, n, r,
     &  '  The residual matrix (A*B-Lambda*I)*X:' )

      end if

      return
      end
      subroutine test10 ( )

c*********************************************************************72
c
cc TEST10 tests RSGBA.
c
c  Discussion:
c
c    RGGBA is for a real generalized eigenvalue problem of the form
c
c      B*A*x = lambda*x
c
c    with A symmetric and B positive definite symmetric.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 January 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      double precision a(n,n)
      double precision a2(n,n)
      double precision b(n,n)
      double precision b2(n,n)
      double precision fv1(n)
      double precision fv2(n)
      integer i
      integer ierr
      integer j
      integer k
      integer matz
      double precision r(n,n)
      double precision w(n)
      double precision x(n,n)

      do i = 1, n
        do j = 1, n
          a(i,j) = abs ( i - j )
        end do
      end do

      do j = 1, n
        do i = 1, n
          a2(i,j) = a(i,j)
        end do
      end do

      do i = 1, n
        do j = 1, n
          if ( i .eq. j ) then
            b(i,j) = 2.0D+00
          else if ( abs ( i - j ) .eq. 1 ) then
            b(i,j) = - 1.0D+00
          else
            b(i,j) = 0.0D+00
          end if
        end do
      end do

      do j = 1, n
        do i = 1, n
          b2(i,j) = b(i,j)
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST10:'
      write ( *, '(a)' )
     &  '  RSGBA for real symmetric generalized problem.'
      write ( *, '(a)' ) '  Find scalars LAMBDA and vectors X so that'
      write ( *, '(a)' ) '    B*A*X = LAMBDA * X'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix order = ', n

      call r8mat_print ( n, n, a, '  The matrix A:' )

      call r8mat_print ( n, n, b, '  The matrix B:' )

      matz = 1

      call rsgba ( n, n, a, b, w, matz, x, fv1, fv2, ierr )

      if ( ierr .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Warning!'
        write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
      end if

      call r8vec_print ( n, w, '  The eigenvalues Lambda:' )

      if ( matz .ne. 0 ) then

        call r8mat_print ( n, n, x, '  The eigenvector matrix X:' )

        do i = 1, n
          do j = 1, n
            r(i,j) = 0.0D+00
            do k = 1, n
              r(i,j) = r(i,j) + a2(i,k) * x(k,j)
            end do
          end do
        end do

        do i = 1, n
          do j = 1, n
            a2(i,j) = r(i,j)
          end do
        end do

        do i = 1, n
          do j = 1, n
            r(i,j) = 0.0D+00
            do k = 1, n
              r(i,j) = r(i,j) + b2(i,k) * a2(k,j)
            end do
          end do
        end do

        do j = 1, n
          do i = 1, n
            r(i,j) = r(i,j) - w(j) * x(i,j)
          end do
        end do

        call r8mat_print ( n, n, r,
     &  '  The residual matrix (B*A-Lambda*I)*X:' )

      end if

      return
      end
      subroutine test11 ( )

c*********************************************************************72
c
cc TEST11 tests RSM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )
      integer m
      parameter ( m = n )

      double precision a(n,n)
      double precision a2(n,n)
      double precision fwork(8*n)
      integer i
      integer ierr
      integer iwork(n)
      integer j
      integer k
      integer matz
      double precision r(n,m)
      double precision w(n)
      double precision x(n,m)

      a(1,1) = 5.0D+00
      a(1,2) = 4.0D+00
      a(1,3) = 1.0D+00
      a(1,4) = 1.0D+00

      a(2,1) = 4.0D+00
      a(2,2) = 5.0D+00
      a(2,3) = 1.0D+00
      a(2,4) = 1.0D+00

      a(3,1) = 1.0D+00
      a(3,2) = 1.0D+00
      a(3,3) = 4.0D+00
      a(3,4) = 2.0D+00

      a(4,1) = 1.0D+00
      a(4,2) = 1.0D+00
      a(4,3) = 2.0D+00
      a(4,4) = 4.0D+00

      do j = 1, n
        do i = 1, n
          a2(i,j) = a(i,j)
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST11'
      write ( *, '(a)' )
     &  '  RSM computes some eigenvalues and eigenvectors'
      write ( *, '(a)' ) '  of a real symmetric matrix.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix order = ', n
      write ( *, '(a,i8)' ) '  Number of eigenvectors desired = ', m

      call r8mat_print ( n, n, a, '  The matrix A:' )

      call rsm ( n, n, a, w, m, x, fwork, iwork, ierr )

      if ( ierr .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Warning!'
        write ( *, '(a,i8)' ) '  The error return flag IERR = ',ierr
      end if

      call r8vec_print ( n, w, '  The eigenvalues Lambda:' )

      if ( 0 .lt. m ) then

        call r8mat_print ( n, m, x, '  The eigenvector matrix X:' )

        do i = 1, n
          do j = 1, m
            r(i,j) = 0.0D+00
            do k = 1, n
              r(i,j) = r(i,j) + a2(i,k) * x(k,j)
            end do
          end do
        end do

        do j = 1, m
          do i = 1, n
            r(i,j) = r(i,j) - w(j) * x(i,j)
          end do
        end do

        call r8mat_print ( n, m, r, '  The residual (A-Lambda*I)*X:' )

      end if

      return
      end
      subroutine test12 ( )

c*********************************************************************72
c
cc TEST12 tests RSP.
c
c  Discussion:
c
c    RSP is for the eigenvalues of a real symmetric packed matrix.
c
c    A is symmetric.  Because of this, we know that the eigenvalues
c    of A must be real (rather than complex) numbers.
c
c
c    The entries of A are
c
c    (5 4 1 1)
c    (4 5 1 1)
c    (1 1 4 2)
c    (1 1 2 4)
c
c    The eigenvalues of A are (10, 5, 2, 1)
c
c    One set of eigenvectors of A is:
c
c    ( 2 -1  0 -1)
c    ( 2 -1  0  1)
c    ( 1  2 -1  0)
c    ( 1  2  1  0)
c
c    However, this set is not orthonormal, and EISPACK will compute
c    a different set of values.
c
c    Note that the I-th eigenvector corresponding to the I-th eigenvalue
c    consists of the I-th column of the above matrix of eigenvectors.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )
      integer nv
      parameter ( nv = ( n * ( n + 1 ) ) / 2 )

      double precision a(nv)
      double precision a2(n,n)
      double precision fv1(n)
      double precision fv2(n)
      integer i
      integer ierr
      integer j
      integer k
      integer matz
      double precision r(n,n)
      double precision w(n)
      double precision x(n,n)
c
c  Set the values in the matrix.
c
      a(1) = 5.0D+00

      a(2) = 4.0D+00
      a(3) = 5.0D+00

      a(4) = 1.0D+00
      a(5) = 1.0D+00
      a(6) = 4.0D+00

      a(7) = 1.0D+00
      a(8) = 1.0D+00
      a(9) = 2.0D+00
      a(10) = 4.0D+00

      a2(1,1) = 5.0D+00
      a2(1,2) = 4.0D+00
      a2(1,3) = 1.0D+00
      a2(1,4) = 1.0D+00

      a2(2,1) = 4.0D+00
      a2(2,2) = 5.0D+00
      a2(2,3) = 1.0D+00
      a2(2,4) = 1.0D+00

      a2(3,1) = 1.0D+00
      a2(3,2) = 1.0D+00
      a2(3,3) = 4.0D+00
      a2(3,4) = 2.0D+00

      a2(4,1) = 1.0D+00
      a2(4,2) = 1.0D+00
      a2(4,3) = 2.0D+00
      a2(4,4) = 4.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST12'
      write ( *, '(a)' )
     &  '  RSP computes the eigenvalues and eigenvectors'
      write ( *, '(a)' ) '  of a real symmetric packed matrix.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix order = ', n

      call r8mat_print ( n, n, a2, '  The matrix A:' )

      matz = 1

      call rsp ( n, n, nv, a, w, matz, x, fv1, fv2, ierr )

      if ( ierr .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Warning!'
        write ( *, '(a,i8)' )
     &  '  The error return flag was IERR = ', ierr
      end if

      call r8vec_print ( n, w, '  The eigenvalues Lambda:' )

      if ( matz .ne. 0 ) then

        call r8mat_print ( n, n, x, '  The eigenvector matrix X:' )

        do i = 1, n
          do j = 1, n
            r(i,j) = 0.0D+00
            do k = 1, n
              r(i,j) = r(i,j) + a2(i,k) * x(k,j)
            end do
          end do
        end do

        do j = 1, n
          do i = 1, n
            r(i,j) = r(i,j) - w(j) * x(i,j)
          end do
        end do

        call r8mat_print ( n, n, r,
     &  '  The residual matrix (A-Lambda*I)*X:' )

      end if

      return
      end
      subroutine test13 ( )

c*********************************************************************72
c
cc TEST13 tests RSPP.
c
c  Discussion:
c
c    RSPP is for some eigenvalues of a real symmetric packed matrix.
c
c    A is symmetric.  Because of this, we know that the eigenvalues
c    of A must be real (rather than complex) numbers.
c
c    The entries of A are
c
c    (5 4 1 1)
c    (4 5 1 1)
c    (1 1 4 2)
c    (1 1 2 4)
c
c    The eigenvalues of A are (10, 5, 2, 1)
c
c    One set of eigenvectors of A is:
c
c    ( 2 -1  0 -1)
c    ( 2 -1  0  1)
c    ( 1  2 -1  0)
c    ( 1  2  1  0)
c
c    However, this set is not orthonormal, and EISPACK will compute
c    a different set of values.
c
c    Note that the I-th eigenvector corresponding to the I-th eigenvalue
c    consists of the I-th column of the above matrix of eigenvectors.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )
      integer m
      parameter ( m = n )
      integer nv
      parameter ( nv = ( n * ( n + 1 ) ) / 2 )

      double precision a(nv)
      double precision a2(n,n)
      double precision bd(n)
      double precision d(n )
      double precision e(n)
      double precision e2(n)
      integer i
      integer ierr
      integer j
      integer k
      integer matz
      double precision r(n,m)
      logical type
      double precision w(m)
      double precision work1(n)
      double precision work2(n)
      double precision work3(n)
      double precision work4(n)
      double precision work6(n)
      double precision x(n,m)
c
c  Set the values in the matrix.
c
      a(1) = 5.0D+00

      a(2) = 4.0D+00
      a(3) = 5.0D+00

      a(4) = 1.0D+00
      a(5) = 1.0D+00
      a(6) = 4.0D+00

      a(7) = 1.0D+00
      a(8) = 1.0D+00
      a(9) = 2.0D+00
      a(10) = 4.0D+00

      a2(1,1) = 5.0D+00
      a2(1,2) = 4.0D+00
      a2(1,3) = 1.0D+00
      a2(1,4) = 1.0D+00

      a2(2,1) = 4.0D+00
      a2(2,2) = 5.0D+00
      a2(2,3) = 1.0D+00
      a2(2,4) = 1.0D+00

      a2(3,1) = 1.0D+00
      a2(3,2) = 1.0D+00
      a2(3,3) = 4.0D+00
      a2(3,4) = 2.0D+00

      a2(4,1) = 1.0D+00
      a2(4,2) = 1.0D+00
      a2(4,3) = 2.0D+00
      a2(4,4) = 4.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST13'
      write ( *, '(a)' )
     &  '  RSPP finds some eigenvalues and eigenvectors of'
      write ( *, '(a)' ) '  a real symmetric packed matrix.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix order = ', n

      call r8mat_print ( n, n, a2, '  The matrix A:' )
c
c  Set MATZ = 0 for no eigenvectors, 1 for eigenvectors.
c
      matz = 1
c
c  TYPE = TRUE to find smallest eigenvalues, FALSE for largest.
c
      type = .true.

      call rspp ( n, nv, a, w, matz, x, ierr, m, type, bd, d, e, e2,
     &  work1, work2, work3, work4, work6 )

      if ( ierr .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Warning!'
        write ( *, '(a,i8)' )
     &  '  The error return flag was IERR = ', ierr
      end if

      call r8vec_print ( m, w, '  The eigenvalues Lambda:' )

      if ( matz .ne. 0 ) then

        call r8mat_print ( n, m, x, '  The eigenvector matrix X:' )

        do i = 1, n
          do j = 1, m
            r(i,j) = 0.0D+00
            do k = 1, n
              r(i,j) = r(i,j) + a2(i,k) * x(k,j)
            end do
          end do
        end do

        do j = 1, m
          do i = 1, n
            r(i,j) = r(i,j) - w(j) * x(i,j)
          end do
        end do

        call r8mat_print ( n, m, r,
     &  '  The residual matrix (A-Lambda*I)*X:' )

      end if

      return
      end
      subroutine test14 ( )

c*********************************************************************72
c
cc TEST14 tests RST.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      double precision a(n,n)
      double precision e(n)
      integer i
      integer ierr
      integer j
      integer k
      integer matz
      double precision r(n,n)
      double precision w(n)
      double precision x(n,n)
c
c  Here is where the matrix is defined.
c
      do i = 1, n
        w(i) = 2.0D+00
      end do

      e(1) = 0.0D+00
      do i = 2, n
        e(i) = -1.0D+00
      end do
c
c  We only set up and store the matrix A this way in order to make it easy
c  to compute the residual.
c
      do i = 1, n
        do j = 1, n
          if ( i .eq. j ) then
            a(i,j) = 2.0D+00
          else if ( abs ( i - j ) .eq. 1 ) then
            a(i,j) = - 1.0D+00
          else
            a(i,j) = 0.0D+00
          end if
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST14'
      write ( *, '(a)' )
     &  '  RST computes the eigenvalues and eigenvectors'
      write ( *, '(a)' ) '  of a real symmetric tridiagonal matrix.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix order = ', n

      call r8mat_print ( n, n, a, '  The matrix A:' )

      matz = 1

      call rst ( n, n, w, e, matz, x, ierr )

      if ( ierr .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Warning!'
        write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
      end if

      call r8vec_print ( n, w, '  The eigenvalues Lambda:' )

      if ( matz .ne. 0 ) then

        call r8mat_print ( n, n, x, '  The eigenvector matrix X:' )

        do i = 1, n
          do j = 1, n
            r(i,j) = 0.0D+00
            do k = 1, n
              r(i,j) = r(i,j) + a(i,k) * x(k,j)
            end do
          end do
        end do

        do j = 1, n
          do i = 1, n
            r(i,j) = r(i,j) - w(j) * x(i,j)
          end do
        end do

        call r8mat_print ( n, n, r,
     &  '  The residual matrix (A-Lambda*I)*X:' )

      end if

      return
      end
      subroutine test15 ( )

c*********************************************************************72
c
cc TEST15 tests RT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )

      double precision a(n,3)
      double precision a2(n,n)
      double precision e(n)
      double precision fv1(n)
      integer i
      integer ierr
      integer j
      integer k
      integer matz
      double precision r(n,n)
      double precision w(n)
      double precision x(n,n)
c
c  Here is where the matrix is defined.
c
      do i = 2, n
        a(i,1) = - 1.0D+00
      end do

      do i = 1, n
        a(i,2) =   2.0D+00
      end do

      do i = 1, n - 1
        a(i,3) = - 1.0D+00
      end do
c
c  We only set up and store the matrix A this way in order to make it easy
c  to compute the residual.
c
      do i = 1, n
        do j = 1, n
          if ( i .eq. j ) then
            a2(i,j) = 2.0D+00
          else if ( abs ( i - j ) .eq. 1 ) then
            a2(i,j) = - 1.0D+00
          else
            a2(i,j) = 0.0D+00
          end if
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST15'
      write ( *, '(a)' )
     &  '  RT computes the eigenvalues and eigenvectors'
      write ( *, '(a)' )
     &  '  of a real sign-symmetric tridiagonal matrix.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix order = ', n

      call r8mat_print ( n, n, a2, '  The matrix A:' )

      matz = 1

      call rt ( n, n, a, w, matz, x, fv1, ierr )

      if ( ierr .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Warning!'
        write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
      end if

      call r8vec_print ( n, w, '  The eigenvalues Lambda:' )

      if ( matz .ne. 0 ) then

        call r8mat_print ( n, n, x, '  The eigenvector matrix X:' )

        do i = 1, n
          do j = 1, n
            r(i,j) = 0.0D+00
            do k = 1, n
              r(i,j) = r(i,j) + a2(i,k) * x(k,j)
            end do
          end do
        end do

        do j = 1, n
          do i = 1, n
            r(i,j) = r(i,j) - w(j) * x(i,j)
          end do
        end do

        call r8mat_print ( n, n, r,
     &  '  The residual matrix (A-Lambda*I)*X:' )

      end if

      return
      end
      subroutine test16 ( )

c*********************************************************************72
c
cc TEST16 tests SVD.
c
c  Discussion:
c
c    In our special example, the matrix is square and symmetric.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )
      integer m
      parameter ( m = n )

      double precision a(m,n)
      integer i
      integer ierr
      integer j
      integer k
      logical matu
      logical matv
      double precision r(m,n)
      double precision rv1(n)
      double precision u(m,n)
      double precision v(n,n)
      double precision w(n)
c
c  Set the values of the matrix.
c
      a(1,1) = 0.9900D+00
      a(1,2) = 0.0020D+00
      a(1,3) = 0.0060D+00
      a(1,4) = 0.0020D+00

      a(2,1) = 0.0020D+00
      a(2,2) = 0.9900D+00
      a(2,3) = 0.0020D+00
      a(2,4) = 0.0060D+00

      a(3,1) = 0.0060D+00
      a(3,2) = 0.0020D+00
      a(3,3) = 0.9900D+00
      a(3,4) = 0.0020D+00

      a(4,1) = 0.0020D+00
      a(4,2) = 0.0060D+00
      a(4,3) = 0.0020D+00
      a(4,4) = 0.9900D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST16'
      write ( *, '(a)' )
     &  '  SVD computes the singular value decomposition'
      write ( *, '(a)' ) '  of a real general matrix.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix order = ', n

      call r8mat_print ( m, n, a, '  The matrix A:' )

      matu = .true.
      matv = .true.

      call svd ( m, m, n, a, w, matu, u, matv, v, ierr, rv1 )

      if ( ierr .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Warning!'
        write ( *, '(a,i8)' ) '  The error return flag IERR = ', ierr
      end if

      call r8vec_print ( n, w, '  The singular values S' )

      call r8mat_print ( m, n, u, '  The U matrix:' )

      call r8mat_print ( n, n, v, '  The V matrix:' )

      do j = 1, n
        do i = 1, n
          v(i,j) = w(j) * v(i,j)
        end do
      end do

        do i = 1, m
          do j = 1, n
            r(i,j) = 0.0D+00
            do k = 1, n
              r(i,j) = r(i,j) + u(i,k) * v(j,k)
            end do
          end do
        end do

      call r8mat_print ( m, n, r,
     &  '  The product U * S * Transpose(V):' )

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
c    Input, character ( len = * ) TITLE, an optional title.
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
      integer title_length

      title_length = len_trim ( title )

      if ( 0 .lt. title_length ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) title(1:title_length)
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

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

        end do

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
c    Input, character * ( * ) TITLE, an optional title.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      integer s_len_trim
      character ( len = * ) title
      integer title_length

      title_length = s_len_trim ( title )
      if ( 0 .lt. title_length ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) title(1:title_length)
      end if

      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,2x,g16.8)' ) i, a(i)
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
c    Input, character ( len = * ) TITLE, a title to be printed first.
c    TITLE may be blank.
c
      implicit none

      integer n

      double precision a1(n)
      double precision a2(n)
      integer i
      integer s_len_trim
      character ( len = * ) title
      integer title_len

      title_len = s_len_trim ( title )

      if ( 0 .lt. title_len ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) title(1:title_len)
      end if

      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) i, a1(i), a2(i)
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
