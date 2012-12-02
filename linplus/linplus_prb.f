      program main

c*********************************************************************72
c
cc MAIN is the main program for LINPLUS_PRB.
c
c  Discussion:
c
c    LINPLUS_PRB calls the LINPLUS test routines.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 May 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LINPLUS_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version:'
      write ( *, '(a)' ) '  Test the LINPLUS library.'

      call test01 ( )
      call test0196 ( )

      call test29 ( )
      call test295 ( )
      call test31  ( )
      call test315 ( )
      call test317 ( )
      call test34 ( )
      call test345 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LINPLUS_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests C83_CR_FA, C83_CR_SLS.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 May 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )
      integer nb
      parameter ( nb = 1 )

      double complex a(3,n)
      double complex a_cr(3,0:2*n)
      double complex b(n,nb)
      integer i
      integer j
      double complex x(n,nb)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) 
     &  '  C83_CR_FA factors a complex tridiagonal matrix;'
      write ( *, '(a)' ) 
     &  '  C83_CR_SLS solves 1 or more factored systems.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix order N = ', n
c
c  Set the matrix values.
c
      a(1,1) = 0.0D+00
      do j = 2, n
        a(1,j) = dcmplx ( - 1, - j + 1 )
      end do

      do j = 1, n
        a(2,j) = dcmplx ( 2, 2 * j )
      end do

      do j = 1, n - 1
        a(3,j) = dcmplx ( - 1, - j - 1 )
      end do
      a(3,n) = 0.0D+00
c
c  Set the desired solution.
c
      do i = 1, n
        x(i,1) = dcmplx ( i,  10 * i )
      end do
c
c  Compute the corresponding right hand side.
c
      i = 1
      b(1,1) = dcmplx (   2, 2 * i ) * x(i,1) 
     &       + dcmplx ( - 1, - i ) * x(i+1,1)

      do i = 2, n - 1
        b(i,1) = dcmplx ( - 1, - i ) * x(i-1,1) 
     &         + dcmplx (   2, 2 * i ) * x(i,1) 
     &         + dcmplx ( - 1, - i ) * x(i+1,1)
      end do

      b(n,1) = dcmplx ( - 1, - i ) * x(i-1,1) 
     &       + dcmplx (   2, 2 * i ) * x(i,1)
c
c  Factor the matrix.
c
      call c83_cr_fa ( n, a, a_cr )
c
c  Solve the linear system.
c
      call c83_cr_sls ( n, a_cr, nb, b, x )

      call c8mat_print_some ( n, nb, x, 1, 1, 10, nb, '  Solution:' )

      return
      end
      subroutine test0196 ( )

c*********************************************************************72
c
cc TEST0196 tests R8VEC_TO_R8GE, R8GE_TO_R8VEC.
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
      implicit none

      integer m
      parameter ( m = 4 )
      integer n
      parameter ( n = 6 )

      double precision a(m,n)
      integer i
      integer j
      integer k
      double precision x(m*n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0196'
      write ( *, '(a)' ) '  For a general matrix,'
      write ( *, '(a)' ) 
     &  '  R8VEC_TO_R8GE converts a real vector to an R8GE matrix.'
      write ( *, '(a)' ) 
     &  '  R8GE_TO_R8VEC converts an R8GE matrix to a real vector.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
      write ( *, '(a,i8)' ) '  Matrix columns N = ', n

      call r8ge_indicator ( m, n, a )

      call r8ge_print ( m, n, a, '  The R8GE indicator matrix:' )

      call r8ge_to_r8vec ( m, n, a, x )

      k = 0
      do j = 1, n
        do i = 1, m
          k = k + 1
          write ( *, '(3i8,g14.6)' ) i, j, k, x(k)
        end do
      end do

      call r8vec_to_r8ge ( m, n, x, a )

      call r8ge_print ( m, n, a, 
     &  '  The recovered R8GE indicator matrix:' )

      return
      end
      subroutine test29

c*********************************************************************72
c
cc TEST29 tests R8GE_DET.
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
      implicit none

      integer n
      parameter ( n = 4 )

      double precision a(n,n)
      double precision det
      double precision exact
      integer i
      integer info
      integer j
      integer pivot(n)
      double precision x
      double precision y

      x = 2.0D+00
      y = 3.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST29'
      write ( *, '(a)' ) '  For a matrix in general storage,'
      write ( *, '(a)' ) '  R8GE_DET computes the determinant.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix order N = ', n
c
c  Set the matrix.
c
      do i = 1, n
        do j = 1, n
          if ( i .eq. j ) then
            a(i,j) = x + y
          else
            a(i,j) = y
          end if
        end do
      end do
c
c  Factor the matrix.
c
      call r8ge_fa ( n, a, pivot, info )
c
c  Compute the determinant.
c
      call r8ge_det ( n, a, pivot, det )

      exact = x**(n-1) * ( x + real ( n ) * y )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) 
     &  '  R8GE_DET computes the determinant = ', det
      write ( *, '(a,g14.6)' ) 
     &  '  Exact determinant =                 ', exact

      return
      end
      subroutine test295 ( )

c*********************************************************************72
c
cc TEST295 tests R8GE_DILU.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 August 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer ncol
      parameter ( ncol = 3 )
      integer nrow
      parameter ( nrow = 3 )
      integer n
      parameter ( n = nrow * ncol )
      integer m
      parameter ( m = n )

      double precision a(m,n)
      double precision d(m)
      integer i
      integer j

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST295'
      write ( *, '(a)' ) '  For a matrix in general storage,'
      write ( *, '(a)' ) '  R8GE_DILU returns the DILU factors.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
      write ( *, '(a,i8)' ) '  Matrix columns N = ', n

      do i = 1, nrow * ncol
        do j = 1, nrow * ncol

          if ( i == j ) then
            a(i,j) = 4.0D+00
          else if ( i .eq. j + 1 .or. i .eq. j - 1 .or. 
     &              i .eq. j + nrow .or. i .eq. j - nrow ) then
            a(i,j) = -1.0D+00
          else
            a(i,j) = 0.0D+00
          end if

        end do
      end do

      call r8ge_print ( m, n, a, '  Matrix A:' )
!
!  Compute the incomplete LU factorization.
!
      call r8ge_dilu ( m, n, a, d )

      call r8vec_print ( m, d, '  DILU factor:' )

      return
      end
      subroutine test31 ( )

c*********************************************************************72
c
cc TEST31 tests R8GE_FA, R8GE_SL.
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
      implicit none

      integer n
      parameter ( n = 5 )

      double precision a(n,n)
      double precision b(n)
      integer info
      integer job
      integer pivot(n)
      integer seed
      double precision x(n)

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST31'
      write ( *, '(a)' ) '  For a matrix in general storage,'
      write ( *, '(a)' ) '  R8GE_FA computes the LU factors,'
      write ( *, '(a)' ) '  R8GE_SL solves a factored system.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix order N = ', n
!
!  Set the matrix.
!
      call r8ge_random ( n, n, seed, a )

      call r8ge_print ( n, n, a, '  The matrix:' )
!
!  Set the desired solution.
!
      call r8vec_indicator ( n, x )
!
!  Compute the corresponding right hand side.
!
      call r8ge_mxv ( n, n, a, x, b )
!
!  Factor the matrix.
!
      call r8ge_fa ( n, a, pivot, info )
 
      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST31 - Fatal error!'
        write ( *, '(a)' ) '  R8GE_FA declares the matrix is singular!'
        write ( *, '(a,i8)' ) '  The value of INFO is ', info
        return
      end if
!
!  Display the gory details.
!
      call r8mat_print ( n, n, a, '  The compressed LU factors:' )

      call i4vec_print ( n, pivot, '  The pivot vector P:' )
!
!  Solve the linear system.
!
      job = 0
      call r8ge_sl ( n, a, pivot, b, job )

      call r8vec_print ( n, b, '  Solution:' )

      return
      end
      subroutine test315

!*****************************************************************************80
!
!! TEST315 tests R8GE_ILU.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    26 August 2006
!
!  Author:
!
!    John Burkardt
!
      implicit none

      integer ncol
      parameter ( ncol = 3 )
      integer nrow
      parameter ( nrow = 3 )
      integer n
      parameter ( n = nrow * ncol )
      integer m
      parameter ( m = n )

      double precision a(m,n)
      integer i
      integer j
      double precision l(m,m)
      double precision lu(m,n)
      double precision u(m,n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST315'
      write ( *, '(a)' ) '  For a matrix in general storage,'
      write ( *, '(a)' ) '  R8GE_ILU returns the ILU factors.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
      write ( *, '(a,i8)' ) '  Matrix columns N = ', n

      do i = 1, nrow * ncol
        do j = 1, nrow * ncol

          if ( i .eq. j ) then
            a(i,j) = 4.0D+00
          else if ( i .eq. j + 1 .or. i .eq. j - 1 .or. 
     &              i .eq. j + nrow .or. i .eq. j - nrow ) then
            a(i,j) = -1.0D+00
          else
            a(i,j) = 0.0D+00
          end if

        end do
      end do

      call r8ge_print ( m, n, a, '  Matrix A:' )
!
!  Compute the incomplete LU factorization.
!
      call r8ge_ilu ( m, n, a, l, u )

      call r8ge_print ( m, m, l, '  Factor L:' )

      call r8ge_print ( m, n, u, '  Factor U:' )
 
      call r8ge_mxm ( m, m, n, l, u, lu )

      call r8ge_print ( m, n, lu, '  Product L*U:' )

      return
      end
      subroutine test317 ( )

c*********************************************************************72
c
cc TEST317 tests R8GE_INDICATOR.
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
      implicit none

      integer m
      parameter ( m = 7 )
      integer n
      parameter ( n = 5 )

      double precision a(m,n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST317'
      write ( *, '(a)' ) '  For a matrix in general storage,'
      write ( *, '(a)' ) 
     &  '  R8GE_INDICATOR sets up the indicator matrix.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix rows M =    ', m
      write ( *, '(a,i8)' ) '  Matrix columns N = ', n

      call r8ge_indicator ( m, n, a )

      call r8ge_print ( m, n, a, '  The R8GE indicator matrix:' )

      return
      end
      subroutine test34 ( )

c*********************************************************************72
c
cc TEST34 tests R8GE_FS.
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
      implicit none

      integer n
      parameter ( n = 10 )

      double precision a(n,n)
      double precision b(n)
      integer info
      integer seed
      double precision x(n)

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST34'
      write ( *, '(a)' ) '  For a matrix in general storage,'
      write ( *, '(a)' ) '  R8GE_FS factors and solves a linear system.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix order N = ', n
c
c  Set the matrix.
c
      call r8ge_random ( n, n, seed, a )
c
c  Set the desired solution.
c
      call r8vec_indicator ( n, x )
c
c  Compute the corresponding right hand side.
c
      call r8ge_mxv ( n, n, a, x, b )
c
c  Factor and solve the system.
c
      call r8ge_fs ( n, a, b, info )
      
      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST34 - Fatal errorc'
        write ( *, '(a)' ) '  R8GE_FS reports the matrix is singular.'
        return
      end if

      call r8vec_print ( n, b, '  Solution:' )

      return
      end
      subroutine test345 ( )

c*********************************************************************72
c
cc TEST345 tests R8GE_FSS.
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
      implicit none

      integer n
      parameter ( n = 10 )
      integer nb
      parameter ( nb = 3 )

      double precision a(n,n)
      double precision b(n,nb)
      integer i
      integer info
      integer seed
      double precision x(n)

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST345'
      write ( *, '(a)' ) '  For a matrix in general storage,'
      write ( *, '(a)' ) 
     &  '  R8GE_FSS factors and solves multiple linear system.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix order N = ', n
c
c  Set the matrix.
c
      call r8ge_random ( n, n, seed, a )
c
c  Set the desired solutions.
c
      do i = 1, n
        x(i) = 1.0D+00
      end do
      call r8ge_mxv ( n, n, a, x, b(1,1) )

      do i = 1, n
        x(i) = dble ( i )
      end do
      call r8ge_mxv ( n, n, a, x, b(1,2) )

      do i = 1, n
        x(i) = dble ( mod ( i - 1, 3 ) + 1 )
      end do
      call r8ge_mxv ( n, n, a, x, b(1,3) )
c
c  Factor and solve the system.
c
      call r8ge_fss ( n, a, nb, b, info )
      
      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST345 - Fatal errorc'
        write ( *, '(a)' ) '  R8GE_FSS reports the matrix is singular.'
        return
      end if

      call r8mat_print ( n, nb, b, '  Solutions:' )

      return
      end
