      program main

c*********************************************************************72
c
cc MAIN is a test program for exponential propagator using Arnoldi approach.
c
c  Discussion:
c
c    This main program is a very simple test using diagonal matrices
c    (Krylov subspace methods are blind to the structure of the matrix
c    except for symmetry). This provides a good way of testing the
c    accuracy of the method as well as the error estimates.
c
c  Modified:
c
c    02 July 2005
c
      implicit none

      integer nmax
      parameter ( nmax = 150 )
      integer ih0
      parameter ( ih0 = 60 )
      integer ndmx
      parameter ( ndmx = 20 )

      integer nzmax
      parameter ( nzmax = 7 * nmax )

      double precision a(nzmax)
      double precision :: a0 = 0.0
      double precision :: b0 = 1.0
      double precision ddot
      double precision eps
      double precision h
      integer ioff(10)
      integer j
      integer k
      integer m
      integer n
      integer ndiag
      double precision t
      double precision tn
      double precision u(ih0*nmax)
      double precision w(nmax)
      double precision w1(nmax)
      double precision x(nmax)
      double precision y(nmax)

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPARSEKIT_PRB12'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Test EXPPRO, which computes the matrix exponential.'
c
c  Set dimension of matrix.
c
      n = 100
c
c  Define matrix.
c  A is a single diagonal matrix (ndiag = 1 and ioff(1) = 0 )
c
      ndiag = 1
      ioff(1) = 0
c
c  entries in the diagonal are uniformly distributed.
c
      h = 1.0 / dble ( n + 1 )
      do j = 1, n
        a(j) = dble ( j + 1 ) / dble ( n + 1 )
      end do
c
c  Set a value for TN
c
      tn = 2.0

      eps = 0.0001

      m = 5
      write ( *, '(a,i6)' ) '  Dimension of Krylov subspace M = ', m
c
c  Define initial conditions: chosen so that solution = (1,1,1,1..1)^T
c
      do j = 1, n
        w(j) = exp ( a(j) * tn )
        w1(j) = w(j)
      end do

      call expprod ( n, m, eps, tn, u, w, x, y, a, ioff, ndiag )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  First 10 components of final answer:'
      write ( *, '(a)' ) ' '
      do k = 1, 10
        write ( *, * ) w(k)
      end do
      write ( *, '(a)' ) ' '

      do k = 1, n
        w1(k) = exp ( -a(k) * tn ) * w1(k)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  First 10 components of exact solution '
      write ( *, '(a)' ) ' '
      do k = 1, 10
        write ( *, * ) w1(k)
      end do
c
c  Compute actual 2-norm of error.
c
      t = 0.0
      do k = 1, n
        t = t + ( w1(k) - w(k) )**2
      end do
      t = sqrt ( t / ddot ( n, w, 1, w, 1 ) )

      write ( *, '(a)' ) ' '
      write ( *, * ) '  RMS error (approx-exact)=', t

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPARSEKIT_PRB12'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine oped ( n, x, y, diag, ioff, ndiag )

c*****************************************************************************80
c
cc OPED performs a matrix by vector multiplication.
c
c  Discussion:
c
c    The matrix is diagonally structured and stored in diagonal
c    format.
c
      implicit none

      integer n
      integer ndiag

      double precision diag(n,ndiag)
      integer i
      integer i1
      integer i2
      integer io
      integer ioff(ndiag)
      integer j
      integer k
      double precision x(n)
      double precision y(n)

      do i = 1, n
        y(i) = 0.0
      end do

      do j = 1, ndiag
        io = ioff(j)
        i1 = max ( 1, 1 - io )
        i2 = min ( n, n - io )
        do k = i1, i2
          y(k) = y(k) + diag(k,j) * x(k+io)
        end do
      end do

      return
      end
