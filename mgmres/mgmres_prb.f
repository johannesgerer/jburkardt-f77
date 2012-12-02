      program main

c*********************************************************************72
c
cc MAIN is the main program for MGMRES_PRB.
c
c  Discussion:
c
c    MGMRES_PRB runs the quick checks for the MGMRES code.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 August 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MGMRES_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the MGMRES library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MGMRES_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests MGMRES_ST on the simple -1,2-1 matrix.
c
c  Discussion:
c
c    This is a very weak test, since the matrix has such a simple
c    structure, is diagonally dominant (though not strictly), 
c    and is symmetric.
c
c    To make the matrix bigger, simply increase the value of N.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 August 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n 
      parameter ( n = 20 )

      integer nz_num
      parameter ( nz_num = 3 * n - 2 )

      double precision a(nz_num)
      integer i
      integer ia(nz_num)
      integer itr_max
      integer j
      integer ja(nz_num)
      integer k
      integer mr
      double precision rhs(n)
      integer test
      double precision tol_abs
      double precision tol_rel
      double precision x_error
      double precision x_estimate(n)
      double precision x_exact(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Test MGMRES_ST on the simple -1,2-1 matrix.'
c
c  Set the matrix.
c
      k = 0

      do i = 1, n

        if ( 1 .lt. i ) then
          k = k + 1
          ia(k) = i
          ja(k) = i - 1
          a(k) = -1.0D+00
        end if

        k = k + 1
        ia(k) = i
        ja(k) = i
        a(k) = 2.0D+00

        if ( i .lt. n ) then
          k = k + 1
          ia(k) = i
          ja(k) = i + 1
          a(k) = -1.0D+00
        end if

      end do
c
c  Set the right hand side:
c
      do i = 1, n - 1
        rhs(i) = 0.0D+00
      end do
      rhs(n) = dble ( n + 1 )
c
c  Set the exact solution.
c
      do i = 1, n
        x_exact(i) = dble ( i )
      end do

      do test = 1, 3
c
c  Set the initial solution estimate.
c
        do i = 1, n
          x_estimate(i) = 0.0D+00
        end do

        x_error = 0.0D+00
        do i = 1, n
          x_error = x_error + ( x_exact(i) - x_estimate(i) )**2
        end do
        x_error = sqrt ( x_error )

        if ( test .eq. 1 ) then
          itr_max = 1
          mr = 20
        else if ( test .eq. 2 ) then
          itr_max = 2
          mr = 10
        else if ( test .eq. 3 ) then
          itr_max = 5
          mr = 4
        end if

        tol_abs = 1.0D-08
        tol_rel = 1.0D-08

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Test ', test
        write ( *, '(a,i8)' ) '  Matrix order N = ', n
        write ( *, '(a,i8)' ) '  Inner iteration limit = ', mr
        write ( *, '(a,i8)' ) '  Outer iteration limit = ', itr_max
        write ( *, '(a,g14.6)' ) '  Initial X_ERROR = ', x_error

        call mgmres_st ( n, nz_num, ia, ja, a, x_estimate, rhs, 
     &    itr_max, mr, tol_abs, tol_rel )

        x_error = 0.0D+00
        do i = 1, n
          x_error = x_error + ( x_exact(i) - x_estimate(i) )**2
        end do
        x_error = sqrt ( x_error )

        write ( *, '(a,g14.6)' ) '  Final X_ERROR = ', x_error

      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests MGMRES_ST on a 9 by 9 matrix.
c
c  Discussion:
c
c    A = 
c      2  0  0 -1  0  0  0  0  0
c      0  2 -1  0  0  0  0  0  0
c      0 -1  2  0  0  0  0  0  0
c     -1  0  0  2 -1  0  0  0  0
c      0  0  0 -1  2 -1  0  0  0
c      0  0  0  0 -1  2 -1  0  0
c      0  0  0  0  0 -1  2 -1  0
c      0  0  0  0  0  0 -1  2 -1
c      0  0  0  0  0  0  0 -1  2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 August 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 9 )
      integer nz_num
      parameter ( nz_num = 23 )

      double precision a(nz_num)
      integer i
      integer ia(nz_num)
      integer itr_max
      integer j
      integer ja(nz_num)
      integer k
      integer mr
      double precision rhs(n)
      integer seed
      integer test
      double precision tol_abs
      double precision tol_rel
      double precision x_error
      double precision x_estimate(n)
      double precision x_exact(n)

      save a
      save ia
      save ja
      save rhs
      save x_exact

      data a /
     &   2.0D+00, -1.0D+00, 
     &   2.0D+00, -1.0D+00, 
     &  -1.0D+00,  2.0D+00, 
     &  -1.0D+00,  2.0D+00, -1.0D+00, 
     &  -1.0D+00,  2.0D+00, -1.0D+00, 
     &  -1.0D+00,  2.0D+00, -1.0D+00, 
     &  -1.0D+00,  2.0D+00, -1.0D+00, 
     &  -1.0D+00,  2.0D+00, -1.0D+00, 
     &  -1.0D+00,  2.0D+00 /
      data ia /
     &  1, 1, 
     &  2, 2, 
     &  3, 3, 
     &  4, 4, 4, 
     &  5, 5, 5, 
     &  6, 6, 6, 
     &  7, 7, 7, 
     &  8, 8, 8, 
     &  9, 9 /
      data ja /
     &  1, 4, 
     &  2, 3, 
     &  2, 3, 
     &  1, 4, 5, 
     &  4, 5, 6, 
     &  5, 6, 7, 
     &  6, 7, 8, 
     &  7, 8, 9, 
     &  8, 9 /
      data rhs /
     &  1.0D+00, 
     &  1.0D+00, 
     &  1.0D+00, 
     &  1.0D+00, 
     &  1.0D+00, 
     &  1.0D+00, 
     &  1.0D+00, 
     &  1.0D+00, 
     &  1.0D+00 /
      data x_exact /
     &  3.5D+00, 
     &  1.0D+00, 
     &  1.0D+00, 
     &  6.0D+00, 
     &  7.5D+00, 
     &  8.0D+00, 
     &  7.5D+00, 
     &  6.0D+00, 
     &  3.5D+00 /

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) 
     &  '  Test MGMRES_ST on a matrix that is not quite '
      write ( *, '(a,i8)' ) '  the -1,2,-1 matrix, of order N = ', n

      do test = 1, 2

        if ( test .eq. 1 ) then

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  First try, use zero initial vector:'
          write ( *, '(a)' ) ' '

          do i = 1, n
            x_estimate(i) = 0.0D+00
          end do

        else

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Second try, use random initial vector:'
          write ( *, '(a)' ) ' '

          call r8vec_uniform_01 ( n, seed, x_estimate )

        end if

        x_error = 0.0D+00
        do i = 1, n
          x_error = x_error + ( x_exact(i) - x_estimate(i) )**2
        end do
        x_error = sqrt ( x_error )

        write ( *, '(a,g14.6)' ) '  Before solving, X_ERROR = ', x_error

        itr_max = 20
        mr = n - 1
        tol_abs = 1.0D-08
        tol_rel = 1.0D-08

        call mgmres_st ( n, nz_num, ia, ja, a, x_estimate, rhs, 
     &    itr_max, mr, tol_abs, tol_rel )

        x_error = 0.0D+00
        do i = 1, n
          x_error = x_error + ( x_exact(i) - x_estimate(i) )**2
        end do
        x_error = sqrt ( x_error )

        write ( *, '(a,g14.6)' ) '  After solving, X_ERROR = ', x_error

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Final solution estimate:'
        write ( *, '(a)' ) ' '
        do i = 1, n
          write ( *, '(2x,i8,2x,g14.6)' ) i, x_estimate(i)
        end do

      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests PMGMRES_ILU_CR on the simple -1,2-1 matrix.
c
c  Discussion:
c
c    This is a very weak test, since the matrix has such a simple
c    structure, is diagonally dominant (though not strictly), 
c    and is symmetric.
c
c    To make the matrix bigger, simply increase the value of N.
c
c    Note that PGMRES_ILU_CR expects the matrix to be stored using the
c    sparse compressed row format.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 August 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 20 )
      integer nz_num
      parameter ( nz_num = ( 3 * n - 2 ) )

      double precision a(nz_num)
      integer i
      integer ia(n+1)
      integer itr_max
      integer j
      integer ja(nz_num)
      integer k
      integer mr
      double precision rhs(n)
      integer test
      double precision tol_abs
      double precision tol_rel
      double precision x_error
      double precision x_estimate(n)
      double precision x_exact(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) 
     &  '  Test PMGMRES_ILU_CR on the simple -1,2-1 matrix.'
c
c  Set the matrix.
c  Note that we use 1-based index values in IA and JA.
c
      k = 1
      ia(1) = 1

      write ( *, '(a)' ) ' '
      write ( *, '(a,i4,a,i4)' ) '  ia(', 1, ') = ', ia(1)

      do i = 1, n

        ia(i+1) = ia(i)

        if ( 1 .lt. i ) then
          ia(i+1) = ia(i+1) + 1
          ja(k) = i - 1
          a(k) = -1.0D+00
          k = k + 1
        end if

        ia(i+1) = ia(i+1) + 1
        ja(k) = i
        a(k) = 2.0D+00
        k = k + 1

        if ( i .lt. n ) then
          ia(i+1) = ia(i+1) + 1
          ja(k) = i + 1
          a(k) = -1.0D+00
          k = k + 1
        end if
        write ( *, '(a,i4,a,i4)' ) '  ia(', i + 1, ') = ', ia(i+1)
      end do
c
c  Set the right hand side:
c
      do i = 1, n - 1
        rhs(i) = 0.0D+00
      end do
      rhs(n) = dble ( n + 1 )
c
c  Set the exact solution.
c
      do i = 1, n
        x_exact(i) = dble ( i )
      end do

      do test = 1, 3
c
c  Set the initial solution estimate.
c
        do i = 1, n
          x_estimate(i) = 0.0D+00
        end do
        x_error = 0.0D+00
        do i = 1, n
          x_error = x_error + ( x_exact(i) - x_estimate(i) ) ** 2
        end do
        x_error = sqrt ( x_error )

        if ( test .eq. 1 ) then
          itr_max = 1
          mr = 20
        else if ( test .eq. 2 ) then
          itr_max = 2
          mr = 10
        else if ( test .eq. 3 ) then
          itr_max = 5
          mr = 4
        end if

        tol_abs = 1.0D-08
        tol_rel = 1.0D-08

        write ( *, '(a)' ) ' '
        write ( *, '(a,i4)' ) '  Test ', test
        write ( *, '(a,i4)' ) '  Matrix order N = ', n
        write ( *, '(a,i4)' ) '  Inner iteration limit = ', mr
        write ( *, '(a,i4)' ) '  Outer iteration limit = ', itr_max
        write ( *, '(a,g14.6)' ) '  Initial X_ERROR = ', x_error

        call pmgmres_ilu_cr ( n, nz_num, ia, ja, a, x_estimate, rhs, 
     &    itr_max, mr, tol_abs, tol_rel )

        x_error = 0.0D+00
        do i = 1, n
          x_error = x_error + ( x_exact(i) - x_estimate(i) ) ** 2
        end do
        x_error = sqrt ( x_error )

        write ( *, '(a,g14.6)' ) '  Final X_ERROR = ', x_error

      end do

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests PMGMRES_ILU_CR on a simple 5 by 5 matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 August 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 5 )
      integer nz_num
      parameter ( nz_num = 9 )

      double precision a(nz_num)
      integer i
      integer ia(n+1)
      integer itr_max
      integer j
      integer ja(nz_num)
      integer k
      integer mr
      double precision rhs(n)
      integer test
      double precision tol_abs
      double precision tol_rel
      double precision x_error
      double precision x_estimate(n)
      double precision x_exact(n)

      save a
      save ia
      save ja
      save rhs
      save x_exact

      data a /
     &   1.0, 2.0, 1.0, 
     &   2.0, 
     &   3.0, 3.0, 
     &   4.0, 
     &   1.0, 5.0 /
      data ia /
     &  1, 4, 5, 7, 8, 10 /
      data ja /
     &  1, 4, 5, 
     &  2, 
     &  1, 3, 
     &  4, 
     &  2, 5 /
      data rhs /
     &  14.0, 4.0, 12.0, 16.0, 27.0 /
      data x_exact /
     &  1.0, 2.0, 3.0, 4.0, 5.0 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) 
     &  '  Test PMGMRES_ILU_CR on a simple 5 x 5 matrix.'

      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(a,i2,a,i2)' ) '  ia(', i, ') = ', ia(i)
      end do
     
      do test = 1, 3
c
c  Set the initial solution estimate.
c
        do i = 1, n
          x_estimate(i) = 0.0D+00
        end do
        x_error = 0.0D+00
        do i = 1, n
          x_error = x_error + ( x_exact(i) - x_estimate(i) ) ** 2
        end do
        x_error = sqrt ( x_error )

        if ( test .eq. 1 ) then
          itr_max = 1
          mr = 20
        else if ( test .eq. 2 ) then
          itr_max = 2
          mr = 10
        else if ( test .eq. 3 ) then
          itr_max = 5
          mr = 4
        end if

        tol_abs = 1.0D-08
        tol_rel = 1.0D-08

        write ( *, '(a)' ) ' '
        write ( *, '(a,i4)' ) '  Test ', test
        write ( *, '(a,i4)' ) '  Matrix order N = ', n
        write ( *, '(a,i4)' ) '  Inner iteration limit = ', mr
        write ( *, '(a,i4)' ) '  Outer iteration limit = ', itr_max
        write ( *, '(a,g14.6)' ) '  Initial X_ERROR = ', x_error

        call pmgmres_ilu_cr ( n, nz_num, ia, ja, a, x_estimate, rhs, 
     &    itr_max, mr, tol_abs, tol_rel )

        x_error = 0.0D+00
        do i = 1, n
          x_error = x_error + ( x_exact(i) - x_estimate(i) ) ** 2
        end do
        x_error = sqrt ( x_error )

        write ( *, '(a,g14.6)' ) '  Final X_ERROR = ', x_error

      end do

      return
      end
