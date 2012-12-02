      program main

c*********************************************************************72
c
cc MAIN is the main program for CYCLIC_REDUCTION_PRB.
c
c  Discussion:
c
c   CYCLIC_REDUCTION_PRB calls the CYCLIC_REDUCTION test routines.
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
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CYCLIC_REDUCTION_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version:'
      write ( *, '(a)' ) '  Test the CYCLIC_REDUCTION library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CYCLIC_REDUCTION_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests C83_CR_FA, C83_CR_SL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 May 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      double complex a(3,n)
      double complex a_cr(3,0:2*n)
      double complex b(n)
      integer i
      integer j
      double complex x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) 
     &  '  C83_CR_FA factors a complex tridiagonal matrix;'
      write ( *, '(a)' ) '  C83_CR_SL solves a factored system.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix order N = ', n
c
c  Set the matrix values.
c
      a(1,1) = 0.0D+00
      do j = 2, n
        a(1,j) = dcmplx ( -1, - ( j - 1 ) )
      end do

      do j = 1, n
        a(2,j) = dcmplx ( 2, 2 * j )
      end do

      do j = 1, n - 1
        a(3,j) = dcmplx ( -1, - ( j + 1 ) )
      end do
      a(3,n) = 0.0D+00
c
c  Set the desired solution.
c
      do i = 1, n
        x(i) = dcmplx ( i, 10 * i )
      end do
c
c  Compute the corresponding right hand side.
c
      call c83_mxv ( n, a, x, b )
c
c  Factor the matrix.
c
      call c83_cr_fa ( n, a, a_cr )
c
c  Solve the linear system.
c
      call c83_cr_sl ( n, a_cr, b, x )

      call c8vec_print ( n, x, '  Solution:' )

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests R83_CR_FA, R83_CR_SLS.
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
      implicit none

      integer n
      parameter ( n = 5 )
      integer nb
      parameter ( nb = 2 )

      double precision a(3,n)
      double precision a_cr(3,0:2*n)
      double precision b(n,nb)
      logical debug
      parameter ( debug = .false. )
      integer i
      integer j
      integer n2
      double precision x(n,nb)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) 
     &  '  R83_CR_FA factors a real tridiagonal matrix;'
      write ( *, '(a)' ) '  R83_CR_SLS solves 1 or more systems.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix order N = ', n
      write ( *, '(a)' ) 
     &  '  Demonstrate multiple system solution method.'
c
c  Set the matrix.
c
      a(1,1) = 0.0D+00
      do j = 2, n
        a(1,j) = -1.0D+00
      end do

      do j = 1, n
        a(2,j) = 2.0D+00
      end do

      do j = 1, n - 1
        a(3,j) = -1.0D+00
      end do
      a(3,n) = 0.0D+00

      if ( debug ) then
        call r83_print ( n, a, '  Input matrix:' )
      end if
c
c  Factor the matrix once.
c
      call r83_cr_fa ( n, a, a_cr )

      if ( debug ) then
        call r83_print ( 2 * n + 1, a_cr, 
     &  '  Cyclic reduction factor information:' )
      end if
c
c  Solve 2 systems simultaneously.
c
      do i = 1, n - 1
        b(i,1) = 0.0D+00
      end do
      b(n,1) = dble ( n + 1 )

      b(1,2) = 1.0D+00
      do i = 2, n - 1
        b(i,2) = 0.0D+00
      end do
      b(n,2) = 1.0D+00
c
c  Solve the linear systems.
c
      call r83_cr_sls ( n, a_cr, nb, b, x )

      call r8mat_print_some ( n, nb, x, 1, 1, 10, nb, '  Solutions:' )

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests R83_CR_FA, R83_CR_SL.
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
      implicit none

      integer n
      parameter ( n = 10 )

      double precision a(3,n)
      double precision a_cr(3,0:2*n)
      double precision b(n)
      logical debug
      parameter ( debug = .false. )
      integer i
      integer j
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  For a real tridiagonal matrix,'
      write ( *, '(a)' ) '  using CYCLIC REDUCTION,'
      write ( *, '(a)' ) '  R83_CR_FA factors;'
      write ( *, '(a)' ) '  R83_CR_SL solves a system.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix order N = ', n
      write ( *, '(a)' ) '  The matrix is NOT symmetric.'
c
c  Set the matrix values.
c
      a(1,1) = 0.0D+00
      do j = 2, n
        a(1,j) = dble ( j )
      end do

      do j = 1, n
        a(2,j) = 4.0D+00 * dble ( j )
      end do

      do j = 1, n - 1
        a(3,j) = dble ( j )
      end do
      a(3,n) = 0.0D+00

      if ( debug ) then
        call r83_print ( n, a, '  The matrix:' )
      end if
c
c  Set the desired solution.
c
      call r8vec_indicator ( n, x )
c
c  Compute the corresponding right hand side.
c
      call r83_mxv ( n, a, x, b )

      do i = 1, n
        x(i) = 0.0D+00
      end do
c
c  Factor the matrix.
c
      call r83_cr_fa ( n, a, a_cr )
c
c  Solve the linear system.
c
      call r83_cr_sl ( n, a_cr, b, x )

      call r8vec_print  ( n, x, '  Solution:' )

      return
      end
