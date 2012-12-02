      program main

c*********************************************************************72
c
cc SPARSE_INTERP_ND_PRB tests SPARSE_INTERP_ND.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 November 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      integer sparse_max

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPARSE_INTERP_ND_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version.'
      write ( *, '(a)' ) '  Test the SPARSE_INTERP_ND library.'
      write ( *, '(a)' ) '  The R8LIB library is also required.'

      m = 1
      sparse_max = 9
      call test01 ( m, sparse_max )

      m = 2
      sparse_max = 9
      call test01 ( m, sparse_max )

      m = 3
      sparse_max = 9
      call test01 ( m, sparse_max )

      m = 4
      sparse_max = 7
      call test01 ( m, sparse_max )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SPARSE_INTERP_ND_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      return
      end
      subroutine test01 ( m, sparse_max )

c*********************************************************************72
c
cc TEST01: sequence of sparse interpolants to an M-dimensional function.
c
c  Discussion:
c
c    We have functions that can generate a Lagrange interpolant to data
c    in M dimensions, with specified order or level in each dimension.
c
c    We use the Lagrange function as the inner evaluator for a sparse
c    grid procedure. 
c
c    The procedure computes sparse interpolants of levels 0 to SPARSE_MAX
c    to a given function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 October 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Local, integer M, the spatial dimension.
c
c    Input, integer SPARSE_MAX, the maximum sparse grid level to try.
c
c  Local Parameters:
c
c    Local, double precision A(M), B(M), the upper and lower variable limits 
c    in each dimension.
c
c    Local, double precision APP_ERROR, the averaged Euclidean norm of the 
c    difference between the sparse interpolant and the exact function at 
c    the interpolation points.
c
c    Local, integer C(0:L_MAX), the sparse grid coefficient vector.
c    Results at level L are weighted by C(L).
c
c    Local, integer IND(M), the 1D indices defining a Lagrange grid.
c    Each entry is a 1d "level" that specifies the order of a 
c    Clenshaw Curtis 1D grid.
c
c    Local, integer L, the current Lagrange grid level.
c
c    Local, integer L_MAX, the current sparse grid level.
c
c    Local, integer MORE, used to control the enumeration of all the
c    Lagrange grids at a current grid level.
c
c    Local, integer ND, the number of points used in a Lagrange grid.
c
c    Local, integer ND_TOTAL, the total number of points used in all the
c    Lagrange interpolants for a given sparse interpolant points that occur
c    in more than one Lagrange grid are counted multiple times.
c
c    Local, integer NI, the number of interpolant evaluation points.
c
c    Local, integer SPARSE_MIN, the minimum sparse grid level to try.
c
c    Local, double precision XD(M,ND), the data points for a Lagrange grid.
c
c    Local, double precision XI(M,NI), the interpolant evaluation points.
c
c    Local, double precision ZD(ND), the data values for a Lagrange grid.
c
c    Local, double precision ZE(NI), the exact function values at XI.
c
c    Local, double precision ZI(NI), the sparse interpolant values at XI.
c
c    Local, double precision ZPI(NI), one set of Lagrange interpolant values at XI.
c
      implicit none

      integer nd_max
      parameter ( nd_max = 2000 )

      integer ni
      parameter ( ni = 100 )

      double precision a(m)
      double precision app_error
      double precision b(m)
      integer c(0:20)
      integer h
      integer i
      integer ind(m)
      integer j
      integer l
      integer l_max
      integer l_min
      integer m
      logical more
      integer nd
      integer nd_total
      double precision r8vec_norm_affine
      integer seed
      integer sparse_max
      integer sparse_min
      integer t
      integer w(0:20)
      double precision xd(m,nd_max)
      double precision xi(m,ni)
      double precision zd(nd_max)
      double precision ze(ni)
      double precision zi(ni)
      double precision zpi(ni)

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) '  Sparse interpolation for a function f(x) of'
      write ( *, '(a)' ) '  M-dimensional argument.  Use a sequence of '
      write ( *, '(a)' ) '  sparse grids of levels 0 through '
      write ( *, '(a)' ) '  SPARSE_MAX.  Invoke a general Lagrange '
      write ( *, '(a)' ) '  interpolant function to do this.'
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) '  Compare the exact function and the '
      write ( *, '(a)' ) '  interpolants at a grid of points.'
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) '  The "order" is the sum of the orders of '
      write ( *, '(a)' ) '  all the product grids used to make a '
      write ( *, '(a)' ) '  particular sparse grid.'
c
c  User input.
c
      write ( *, '(a)' ) ''
      write ( *, '(a,i4)' ) '  Spatial dimension M = ', m
      write ( *, '(a,i4)' ) '  Maximum sparse grid level = ', sparse_max
c
c  Define the region.
c
      do i = 1, m
        a(i) = 0.0D+00
        b(i) = 1.0D+00
      end do
c
c  Define the interpolation evaluation information.
c
      seed = 123456789
      call r8mat_uniform_abvec ( m, ni, a, b, seed, xi )

      write ( *, '(a,i6)' ) 
     &  '  Number of interpolation points is NI = ', ni

      call f_sinr ( m, ni, xi, ze )
c
c  Compute a sequence of sparse grid interpolants of increasing level.
c
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) '   L     Order    ApproxError'
      write ( *, '(a)' ) ''

      sparse_min = 0

      do l_max = sparse_min, sparse_max

        call smolyak_coefficients ( l_max, m, c, w )

        do i = 1, ni
          zi(i) = 0.0D+00
        end do

        nd_total = 0

        l_min = max ( l_max + 1 - m, 0 )

        do l = l_min, l_max

          more = .false.

10        continue
c
c  Define the next product grid.
c
            call comp_next ( l, m, ind, more, h, t )
c
c  Count the grid, find its points, evaluate the data there.
c
            call lagrange_interp_nd_size2 ( m, ind, nd )

            if ( nd_max .lt. nd ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'TEST01 - Fatal error!'
              write ( *, '(a,i4)' ) '  ND = ', nd
              write ( *, '(a,i4)' ) '  ND_MAX = ', nd_max
              stop
            end if

            call lagrange_interp_nd_grid2 ( m, ind, a, b, nd, xd )
            call f_sinr ( m, nd, xd, zd )
c
c  Use the grid to evaluate the interpolant.
c
            call lagrange_interp_nd_value2 ( m, ind, a, b, nd, zd, ni, 
     &        xi, zpi )
c
c  Weighted the interpolant values and add to the sparse grid interpolant.
c
            nd_total = nd_total + nd

            do i = 1, ni
              zi(i) = zi(i) + c(l) * zpi(i)
            end do

            if ( .not. more ) then
              go to 20
            end if

          go to 10

20        continue

        end do
c
c  Compare sparse interpolant and exact function at interpolation points.
c
        app_error = r8vec_norm_affine ( ni, zi, ze ) / dble ( ni )

        write ( *, '(2x,i2,2x,i8,2x,e8.2)' ) l, nd_total, app_error

      end do

      return
      end
      subroutine f_sinr ( m, n, x, z )

c*********************************************************************72
c
cc F_SINR is a scalar function of an M-dimensional argument, to be interpolated.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the spatial dimension.
c
c    Input, integer N, the number of points.
c
c    Input, double precision X(M,N), the points.
c
c    Output, double precision Z(N), the value of the function at each point.
c
      implicit none

      integer m
      integer n

      integer i
      integer j
      double precision r
      double precision x(m,n)
      double precision z(n)
             
      do j = 1, n
        r = 0.0D+00
        do i = 1, m
          r = r + x(i,j) ** 2
        end do
        r = sqrt ( r )
        z(j) = sin ( r )
      end do

      return
      end
