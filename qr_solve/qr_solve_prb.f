      program main

c*********************************************************************72
c
cc MAIN is the main program for QR_SOLVE_PRB.
c
c  Discussion:
c
c    QR_SOLVE_PRB tests the QR_SOLVE library.
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
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'QR_SOLVE_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the QR_SOLVE library.'
      write ( *, '(a)' ) '  QR_SOLVE needs the R8LIB library.'
      write ( *, '(a)' ) '  This test also needs the TEST_LS library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'QR_SOLVE_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests NORMAL_SOLVE.
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
      implicit none

      integer m_max
      parameter ( m_max = 10 )
      integer n_max
      parameter ( n_max = 10 )

      double precision a(m_max,n_max)
      double precision b(m_max)
      double precision b_norm
      integer flag
      integer i
      integer m
      integer n
      integer prob
      integer prob_num
      double precision r1(m_max)
      double precision r1_norm
      double precision r2(m_max)
      double precision r2_norm
      double precision r8vec_norm
      double precision r8vec_norm_affine
      double precision x_diff_norm
      double precision x1(n_max)
      double precision x1_norm
      double precision x2(n_max)
      double precision x2_norm

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) 
     &  '  NORMAL_SOLVE is a function with a simple interface which'
      write ( *, '(a)' ) 
     &  '  solves a linear system A*x = b in the least squares sense.'
      write ( *, '(a)' ) 
     &  '  Compare a tabulated solution X1 to NORMAL_SOLVE result X2.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  NORMAL_SOLVE cannot be applied when N < M,'
      write ( *, '(a)' ) 
     &  '  or if the matrix does not have full column rank.'

      call p00_prob_num ( prob_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  Number of problems = ', prob_num
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Index     M     N     ||B||         ' //
     &  '||X1 - X2||   ||X1||       ||X2||        ||R1||        ||R2||'
      write ( *, '(a)' ) ' '

      do prob = 1, prob_num
c
c  Get problem size.
c
        call p00_m ( prob, m )
        call p00_n ( prob, n )
c
c  Retrieve problem data.
c
        call p00_a ( prob, m, n, a )
        call p00_b ( prob, m, b )
        call p00_x ( prob, n, x1 )

        b_norm = r8vec_norm ( m, b )
        x1_norm = r8vec_norm ( n, x1 )
        call r8mat_mv ( m, n, a, x1, r1 )
        do i = 1, m
          r1(i) = r1(i) - b(i)
        end do
        r1_norm = r8vec_norm ( m, r1 )
c
c  Use NORMAL_SOLVE on the problem.
c
        call normal_solve ( m, n, a, b, x2, flag )

        if ( flag .ne. 0 ) then

          write ( *, 
     &      '(2x,i5,2x,i4,2x,i4,2x,g12.4,2x,a,2x,g12.4,' //
     &      '2x,a,2x,g12.4,2x,a)' ) 
     &      prob, m, n, b_norm, '------------', x1_norm, 
     &      '------------', r1_norm, '------------'

        else

          x2_norm = r8vec_norm ( n, x2 )
          call r8mat_mv ( m, n, a, x2, r2 )
          do i = 1, m
            r2(i) = r2(i) - b(i)
          end do
          r2_norm = r8vec_norm ( m, r2 )
c
c  Compare tabulated and computed solutions.
c
          x_diff_norm = r8vec_norm_affine ( n, x1, x2 )
c
c  Report results for this problem.
c
          write ( *, 
     &      '(2x,i5,2x,i4,2x,i4,2x,g12.4,2x,g12.4,2x,g12.4,' // 
     &      '2x,g12.4,2x,g12.4,2x,g12.4)' ) 
     &      prob, m, n, b_norm, x_diff_norm, x1_norm, x2_norm, 
     &      r1_norm, r2_norm

        end if

      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests QR_SOLVE.
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
      implicit none

      integer m_max
      parameter ( m_max = 10 )
      integer n_max
      parameter ( n_max = 10 )

      double precision a(m_max,n_max)
      double precision b(m_max)
      double precision b_norm
      integer i
      integer m
      integer n
      integer prob
      integer prob_num
      double precision r1(m_max)
      double precision r1_norm
      double precision r2(m_max)
      double precision r2_norm
      double precision r8vec_norm
      double precision r8vec_norm_affine
      double precision x_diff_norm
      double precision x1(n_max)
      double precision x1_norm
      double precision x2(n_max)
      double precision x2_norm

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) 
     &  '  QR_SOLVE is a function with a simple interface which'
      write ( *, '(a)' ) 
     &  '  solves a linear system A*x = b in the least squares sense.'
      write ( *, '(a)' ) 
     &  '  Compare a tabulated solution X1 to the QR_SOLVE result X2.'

      call p00_prob_num ( prob_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  Number of problems = ', prob_num
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Index     M     N     ||B||         ' //
     &  '||X1 - X2||   ||X1||       ||X2||        ||R1||        ||R2||'
      write ( *, '(a)' ) ' '

      do prob = 1, prob_num
c
c  Get problem size.
c
        call p00_m ( prob, m )
        call p00_n ( prob, n )
c
c  Retrieve problem data.
c
        call p00_a ( prob, m, n, a )
        call p00_b ( prob, m, b )
        call p00_x ( prob, n, x1 )

        b_norm = r8vec_norm ( m, b )
        x1_norm = r8vec_norm ( n, x1 )
        call r8mat_mv ( m, n, a, x1, r1 )
        do i = 1, m
          r1(i) = r1(i) - b(i)
        end do
        r1_norm = r8vec_norm ( m, r1 )
c
c  Use QR_SOLVE on the problem.
c
        call qr_solve ( m, n, a, b, x2 )

        x2_norm = r8vec_norm ( n, x2 )
        call r8mat_mv ( m, n, a, x2, r2 )
        do i = 1, m
          r2(i) = r2(i) - b(i)
        end do
        r2_norm = r8vec_norm ( m, r2 )
c
c  Compare tabulated and computed solutions.
c
        x_diff_norm = r8vec_norm_affine ( n, x1, x2 )
c
c  Report results for this problem.
c
        write ( *, 
     &    '(2x,i5,2x,i4,2x,i4,2x,g12.4,2x,g12.4,2x,g12.4,' //
     &    '2x,g12.4,2x,g12.4,2x,g12.4)' ) 
     &    prob, m, n, b_norm, x_diff_norm, x1_norm, x2_norm, 
     &    r1_norm, r2_norm

      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests SVD_SOLVE.
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
      implicit none

      integer m_max
      parameter ( m_max = 10 )
      integer n_max
      parameter ( n_max = 10 )

      double precision a(m_max,n_max)
      double precision b(m_max)
      double precision b_norm
      integer i
      integer m
      integer n
      integer prob
      integer prob_num
      double precision r1(m_max)
      double precision r1_norm
      double precision r2(m_max)
      double precision r2_norm
      double precision r8vec_norm
      double precision r8vec_norm_affine
      double precision x_diff_norm
      double precision x1(n_max)
      double precision x1_norm
      double precision x2(n_max)
      double precision x2_norm

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) 
     &  '  SVD_SOLVE is a function with a simple interface which'
      write ( *, '(a)' ) 
     &  '  solves a linear system A*x = b in the least squares sense.'
      write ( *, '(a)' ) 
     &  '  Compare a tabulated solution X1 to the SVD_SOLVE result X2.'

      call p00_prob_num ( prob_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  Number of problems = ', prob_num
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Index     M     N     ||B||         ||X1 - X2||   ' //
     &  '||X1||       ||X2||        ||R1||        ||R2||'
      write ( *, '(a)' ) ' '

      do prob = 1, prob_num
c
c  Get problem size.
c
        call p00_m ( prob, m )
        call p00_n ( prob, n )
c
c  Retrieve problem data.
c
        call p00_a ( prob, m, n, a )
        call p00_b ( prob, m, b )
        call p00_x ( prob, n, x1 )

        b_norm = r8vec_norm ( m, b )
        x1_norm = r8vec_norm ( n, x1 )
        call r8mat_mv ( m, n, a, x1, r1 )
        do i = 1, m
          r1(i) = r1(i) - b(i)
        end do
        r1_norm = r8vec_norm ( m, r1 )
c
c  Use SVD_SOLVE on the problem.
c
        call svd_solve ( m, n, a, b, x2 )

        x2_norm = r8vec_norm ( n, x2 )
        call r8mat_mv ( m, n, a, x2, r2 )
        do i = 1, m
          r2(i) = r2(i) - b(i)
        end do
        r2_norm = r8vec_norm ( m, r2 )
c
c  Compare tabulated and computed solutions.
c
        x_diff_norm = r8vec_norm_affine ( n, x1, x2 )
c
c  Report results for this problem.
c
        write ( *, 
     &    '(2x,i5,2x,i4,2x,i4,2x,g12.4,2x,g12.4,2x,g12.4,2x,' //
     &    'g12.4,2x,g12.4,2x,g12.4)' ) 
     &    prob, m, n, b_norm, x_diff_norm, x1_norm, x2_norm, 
     &    r1_norm, r2_norm

      end do

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests DQRLS.
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
      implicit none

      integer m
      parameter ( m = 5 )
      integer n
      parameter ( n = 3 )

      double precision a(m,n)
      double precision b(m)
      integer i
      integer ind
      integer itask
      integer j
      integer jpvt(n)
      integer kr
      double precision qraux(n)
      double precision tol
      double precision work(n)
      double precision x(n)

      save b

      data b / 1.0D+00, 2.3D+00, 4.6D+00, 3.1D+00, 1.2D+00 /
c
c  Set up least-squares problem
c  quadratic model, equally-spaced points
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) 
     &  '  DQRLS solves a linear system A*x = b ' //
     &  'in the least squares sense.'

      do i = 1, m
        a(i,1) = 1.0D+00
        do j = 2, n
          a(i,j) = a(i,j-1) * dble ( i )
        end do
      end do

      tol = 1.0D-06

      call r8mat_print ( m, n, a, '  Coefficient matrix A:' )

      call r8vec_print ( m, b, '  Right hand side b:' )
c
c  Solve least-squares problem
c
      itask = 1
      call dqrls ( a, m, m, n, tol, kr, b, x, b, jpvt, qraux, work, 
     &  itask, ind )
c
c  Print results
c
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  Error code =', ind
      write ( *, '(a,i4)' ) '  Estimated matrix rank =', kr

      call r8vec_print ( n, x, '  Least squares solution x:' )

      call r8vec_print ( m, b, '  Residuals A*x-b' )

      return
      end
