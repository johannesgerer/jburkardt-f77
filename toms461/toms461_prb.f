      program main

c*********************************************************************72
c
cc TOMS461_PRB tests TOMS461,
c
c  Modified:
c
c    12 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS461_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS461 library.'

      call test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS461_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests SPNBVP.
c
      implicit none

      integer np
      integer nk

      parameter ( np = 81 )
      parameter ( nk = np - 2 )

      real a
      real b
      real ep
      real exact
      integer i
      integer kg(np)
      real mat(nk,nk)
      integer gt(np)
      real t
      real vg(np)
      real vm(np)
      real vp(np)
      real vq(np)
      real vr(np)
      real x(nk)
      real xdp(nk)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  SPNBVP computes a cubic spline'
      write ( *, '(a)' ) '  approximation to the solution of'
      write ( *, '(a)' ) '    x"(t)=p(t)*x(t)+q(t)*x(g(t))+r(t)'
      write ( *, '(a)' ) '  on the interval [a,b], with boundary'
      write ( *, '(a)' ) '  conditions:'
      write ( *, '(a)' ) '    x(t) = u(t) for t <= a,'
      write ( *, '(a)' ) '    x(t) = v(t) for b <= t.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Because of the argument g(t), this is'
      write ( *, '(a)' ) '  a more complicated problem than it seems.'

      a = 0.0E+00
      b = 3.141592653589793E+00

      ep = 1.0E-04

      write ( *, '(a)' ) ' '
      write ( *, '(a,f10.4)' ) '  A = ', a
      write ( *, '(a,f10.4)' ) '  B = ', b
      write ( *, '(a,g14.6)' ) '  Error tolerance EP = ', ep

      call spnbvp ( a, b, np, nk, x, xdp, ep, gt, kg, vp, vq, vr,
     &  vg, mat, vm )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '        T         X                Exact        Error'
      write ( *, '(a)' ) ' '

      do i = 1, nk

        t = ( real ( nk - i + 1 ) * a
     &       + real (      i     ) * b )
     &       / real ( nk     + 1 )

        exact = sin ( t )

        write ( *, '(2x,f10.4,2x,g14.6,2x,g14.6,2x,g10.2)' )
     &    t, x(i), exact, abs ( x(i) - exact )

      end do

      return
      end
      function g ( t )

c*********************************************************************72
c
cc G
c
      implicit none

      real g
      real t

      g = t - 1.0E+00

      return
      end
      subroutine ludcmp ( a, n )

c*********************************************************************72
c
cc LUDCMP factors a matrix by nonpivoting Gaussian elimination.
c
c  Discussion:
c
c    The storage format is used for a general M by N matrix.  A storage
c    space is made for each logical entry.
c
c    This routine will fail if the matrix is singular, or if any zero
c    pivot is encountered.
c
c    If this routine successfully factors the matrix, LUSUB may be called
c    to solve linear systems involving the matrix.
c
c  Modified:
c
c    13 January 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, real A(N,N).
c    On input, A contains the matrix to be factored.
c    On output, A contains information about the factorization,
c    which must be passed unchanged to LUSUB for solutions.
c
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
      implicit none

      integer n

      real a(n,n)
      integer i
      integer j
      integer k

      do k = 1, n-1

        if ( a(k,k) == 0.0E+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'LUDCMP - Fatal error!'
          write ( *, '(a,i6)' ) '  Zero pivot on step ', k
          stop
        end if

        do i = k+1, n
          a(i,k) = -a(i,k) / a(k,k)
        end do

        do j = k+1, n
          do i = k+1, n
            a(i,j) = a(i,j) + a(i,k) * a(k,j)
          end do
        end do

      end do

      if ( a(n,n) == 0.0E+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'LUDCMP - Fatal error!'
        write ( *, '(a,i6)' ) '  Zero pivot on step ', n
        stop
      end if

      return
      end
      subroutine lusub ( b, a_lu, x, n )

c*********************************************************************72
c
cc LUSUB solves a system factored by LUDCMP.
c
c  Modified:
c
c    13 January 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real B(N), the right hand side vector.
c
c    Input, real A_LU(N,N), the LU factors from LUDCMP.
c
c    Output, real X(N), the solution.
c
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
      implicit none

      integer n

      real a_lu(n,n)
      real b(n)
      integer i
      integer j
      real x(n)

      do i = 1, n
        x(i) = b(i)
      end do

      do j = 1, n-1
        do i = j+1, n
          x(i) = x(i) + a_lu(i,j) * x(j)
        end do
      end do

      do j = n, 1, -1
        x(j) = x(j) / a_lu(j,j)
        do i = 1, j-1
          x(i) = x(i) - a_lu(i,j) * x(j)
        end do
      end do

      return
      end
      function p ( t )

c*********************************************************************72
c
cc P evaluates the coefficient function P(T).
c
      implicit none

      real p
      real t

      p = 1.0E+00

      return
      end
      function q ( t )

c*********************************************************************72
c
cc Q evaluates the coefficient function Q(T).
c
      implicit none

      real q
      real t

      q = 1.0E+00

      return
      end
      function r ( t )

c*********************************************************************72
c
cc R evaluates the right hand side function R(T).
c
      implicit none

      real g
      real p
      real q
      real r
      real t

      r = - sin ( t ) - ( p ( t ) * sin ( t ) + q(t) * sin ( g(t) ) )

      return
      end
      function u ( t )

c*********************************************************************72
c
cc U evaluates the left boundary condition.
c
      implicit none

      real t
      real u

      u = sin ( t )

      return
      end
      function v ( t )

c*********************************************************************72
c
cc V evaluates the right boundary condition.
c
      implicit none

      real t
      real v

      v = sin ( t )

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
c  Modified:
c
c    16 September 2005
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

      character ( len = 8 ) date
      character ( len = 10 ) time

      call date_and_time ( date, time )

      write ( *, '(a8,2x,a10)' ) date, time

      return
      end
