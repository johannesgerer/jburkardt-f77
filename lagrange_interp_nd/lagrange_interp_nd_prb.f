      program main

c*********************************************************************72
c
cc TEST tests the LAGRANGE_INTERP_ND library.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the LAGRANGE_INTERP_ND library.'
c
c  Use the interface that passes in the orders directly.
c
      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
c
c  Repeat tests 1, 2, 3, and 4,
c  using the interface that passes in the orders indirectly,
c  based on the "level".
c
      call test05 ( )
      call test06 ( )
      call test07 ( )
      call test08 ( )
c
c  Experiment with anisotropic orders.
c
      call test09 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ''
      call timestamp ( )

      return
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 interpolates in 1D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 1 )
      integer ni
      parameter ( ni = 5 )

      double precision ab(m,2)
      external f_sinr
      integer i
      integer j
      integer n_1d(m)
      integer nd
      integer seed
      double precision xd(m,5)
      double precision xi(m,ni)
      double precision zd(5)
      double precision ze(ni)
      double precision zi(ni)

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) '  Interpolate in 1D, using orders.'
      write ( *, '(a)' ) 
     &  '  LAGRANGE_INTERP_ND_GRID sets the interpolant.'
      write ( *, '(a)' ) '  LAGRANGE_INTERP_ND_VALUE evaluates it.'

      do i = 1, m
        n_1d(i) = 5
      end do

      do i = 1, m
        ab(i,1) = 0.0D+00
        ab(i,2) = 1.0D+00
      end do

      call lagrange_interp_nd_size ( m, n_1d, nd )

      call lagrange_interp_nd_grid ( m, n_1d, ab, nd, xd )
      call f_sinr ( m, nd, xd, zd )
c
c  Evaluate.
c
      seed = 123456789
      call r8mat_uniform_01 ( m, ni, seed, xi )

      call f_sinr ( m, ni, xi, ze )

      call lagrange_interp_nd_value ( m, n_1d, ab, nd, zd, ni, xi, zi )

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) '         Zinterp          Zexact      Error'
      write ( *, '(a)' ) ''

      do j = 1,  ni
        write ( *, '(2x,g14.6,2x,g14.6,2x,e10.2)' ) 
     &    zi(j), ze(j), abs ( zi(j) - ze(j) )
      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 interpolates in 2D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 2 )
      integer ni
      parameter ( ni = 5 )

      double precision ab(m,2)
      external f_sinr
      integer i
      integer j
      integer n_1d(m)
      integer nd
      integer seed
      double precision xd(m,25)
      double precision xi(m,ni)
      double precision zd(25)
      double precision ze(ni)
      double precision zi(ni)

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST02:'
      write ( *, '(a)' ) '  Interpolate in 2D, using orders.'
      write ( *, '(a)' ) 
     &  '  LAGRANGE_INTERP_ND_GRID sets the interpolant.'
      write ( *, '(a)' ) '  LAGRANGE_INTERP_ND_VALUE evaluates it.'

      do i = 1, m
        n_1d(i) = 5
      end do

      do i = 1, m
        ab(i,1) = 0.0D+00
        ab(i,2) = 1.0D+00
      end do

      call lagrange_interp_nd_size ( m, n_1d, nd )

      call lagrange_interp_nd_grid ( m, n_1d, ab, nd, xd )
      call f_sinr ( m, nd, xd, zd )
c
c  Evaluate.
c
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) '         Zinterp          Zexact      Error'
      write ( *, '(a)' ) ''

      seed = 123456789
      call r8mat_uniform_01 ( m, ni, seed, xi )

      call f_sinr ( m, ni, xi, ze )

      call lagrange_interp_nd_value ( m, n_1d, ab, nd, zd, ni, xi, zi )

      do j = 1,  ni
        write ( *, '(2x,g14.6,2x,g14.6,2x,e10.2)' ) 
     &    zi(j), ze(j), abs ( zi(j) - ze(j) ) 
      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 interpolates in 3D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 3 )
      integer ni
      parameter ( ni = 5 )

      double precision ab(m,2)
      external f_sinr
      integer i
      integer j
      integer n_1d(m)
      integer nd
      integer seed
      double precision xd(m,125)
      double precision xi(m,ni)
      double precision zd(125)
      double precision ze(ni)
      double precision zi(ni)

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST03:'
      write ( *, '(a)' ) '  Interpolate in 3D, using orders.'
      write ( *, '(a)' ) 
     &  '  LAGRANGE_INTERP_ND_GRID sets the interpolant.'
      write ( *, '(a)' ) '  LAGRANGE_INTERP_ND_VALUE evaluates it.'

      do i = 1, m
        n_1d(i) = 5
      end do

      do i = 1, m
        ab(i,1) = 0.0D+00
        ab(i,2) = 1.0D+00
      end do

      call lagrange_interp_nd_size ( m, n_1d, nd )

      call lagrange_interp_nd_grid ( m, n_1d, ab, nd, xd )
      call f_sinr ( m, nd, xd, zd )
c
c  Evaluate.
c
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) '         Zinterp          Zexact      Error'
      write ( *, '(a)' ) ''

      seed = 123456789
      call r8mat_uniform_01 ( m, ni, seed, xi )

      call f_sinr ( m, ni, xi, ze )

      call lagrange_interp_nd_value ( m, n_1d, ab, nd, zd, ni, xi, zi )

      do j = 1,  ni
        write ( *, '(2x,g14.6,2x,g14.6,2x,e10.2)' ) 
     &    zi(j), ze(j), abs ( zi(j) - ze(j) ) 
      end do

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 interpolates in 3D, using increasing resolution.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 3 )
      integer ni
      parameter ( ni = 20 )

      double precision ab(m,2)
      double precision e
      external f_sinr
      integer i
      integer j
      integer l
      integer n_1d(m)
      integer nd
      integer order
      double precision r8vec_norm_affine
      integer seed
      double precision xd(m,35937)
      double precision xi(m,ni)
      double precision zd(35937)
      double precision ze(ni)
      double precision zi(ni)

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST04:'
      write ( *, '(a)' ) '  Interpolate in 3D, using orders.'
      write ( *, '(a)' ) '  Use a sequence of increasing orders.'

      do i = 1, m
        ab(i,1) = 0.0D+00
        ab(i,2) = 1.0D+00
      end do

      seed = 123456789
      call r8mat_uniform_01 ( m, ni, seed, xi )

      call f_sinr ( m, ni, xi, ze )

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) '  Level     Order   Average Error'
      write ( *, '(a)' ) ''

      do l = 0, 5

        if ( l .eq. 0 ) then
          order = 1
        else
          order = 2 ** l + 1
        end if

        do i = 1, m
          n_1d(i) = order
        end do
    
        call lagrange_interp_nd_size ( m, n_1d, nd )

        call lagrange_interp_nd_grid ( m, n_1d, ab, nd, xd )
        call f_sinr ( m, nd, xd, zd )
c
c  Evaluate.
c
        call lagrange_interp_nd_value ( m, n_1d, ab, nd, zd, ni, 
     &    xi, zi )

        e = r8vec_norm_affine ( ni, zi, ze ) / dble ( ni )

        write ( *, '(2x,i5,2x,i8,2x,e10.2)' ) l, nd, e
        
      end do

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 repeats test 1 using levels.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 1 )
      integer ni
      parameter ( ni = 5 )

      double precision ab(m,2)
      external f_sinr
      integer i
      integer ind(m)
      integer j
      integer nd
      integer seed
      double precision xd(m,5)
      double precision xi(m,ni)
      double precision zd(5)
      double precision ze(ni)
      double precision zi(ni)

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST05:'
      write ( *, '(a)' ) '  Repeat test #1, using levels.'
      write ( *, '(a)' ) 
     &  '  LAGRANGE_INTERP_ND_GRID2 sets the interpolant.'
      write ( *, '(a)' ) '  LAGRANGE_INTERP_ND_VALUE2 evaluates it.'

      do i = 1, m
        ind(i) = 2
      end do

      do i = 1, m
        ab(i,1) = 0.0D+00
        ab(i,2) = 1.0D+00
      end do

      call lagrange_interp_nd_size2 ( m, ind, nd )

      call lagrange_interp_nd_grid2 ( m, ind, ab, nd, xd )
      call f_sinr ( m, nd, xd, zd )
c
c  Evaluate.
c
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) '         Zinterp          Zexact      Error'
      write ( *, '(a)' ) ''

      seed = 123456789
      call r8mat_uniform_01 ( m, ni, seed, xi )

      call f_sinr ( m, ni, xi, ze )

      call lagrange_interp_nd_value2 ( m, ind, ab, nd, zd, ni, xi, zi )

      do j = 1,  ni
        write ( *, '(2x,g14.6,2x,g14.6,2x,e10.2)' ) 
     &    zi(j), ze(j), abs ( zi(j) - ze(j) ) 
      end do

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 repeats test 2 using levels.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 2 )
      integer ni
      parameter ( ni = 5 )

      double precision ab(m,2)
      external f_sinr
      integer i
      integer ind(m)
      integer j
      integer nd
      integer seed
      double precision xd(m,25)
      double precision xi(m,ni)
      double precision zd(25)
      double precision ze(ni)
      double precision zi(ni)

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST06:'
      write ( *, '(a)' ) '  Repeat test #2, using levels.'
      write ( *, '(a)' ) 
     &  '  LAGRANGE_INTERP_ND_GRID2 sets the interpolant.'
      write ( *, '(a)' ) '  LAGRANGE_INTERP_ND_VALUE2 evaluates it.'

      do i = 1, m
        ind(i) = 2
      end do

      do i = 1, m
        ab(i,1) = 0.0D+00
        ab(i,2) = 1.0D+00
      end do

      call lagrange_interp_nd_size2 ( m, ind, nd )

      call lagrange_interp_nd_grid2 ( m, ind, ab, nd, xd )
      call f_sinr ( m, nd, xd, zd )
c
c  Evaluate.
c
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) '         Zinterp          Zexact      Error'
      write ( *, '(a)' ) ''

      seed = 123456789
      call r8mat_uniform_01 ( m, ni, seed, xi )

      call f_sinr ( m, ni, xi, ze )

      call lagrange_interp_nd_value2 ( m, ind, ab, nd, zd, ni, xi, zi )

      do j = 1,  ni
        write ( *, '(2x,g14.6,2x,g14.6,2x,e10.2)' ) 
     &    zi(j), ze(j), abs ( zi(j) - ze(j) ) 
      end do

      return
      end
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07 repeats test 3 using levels.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 3 )
      integer ni
      parameter ( ni = 5 )

      double precision ab(m,2)
      external f_sinr
      integer i
      integer ind(m)
      integer j
      integer nd
      integer seed
      double precision xd(m,125)
      double precision xi(m,ni)
      double precision zd(125)
      double precision ze(ni)
      double precision zi(ni)

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST07:'
      write ( *, '(a)' ) '  Repeat test #3,  using levels.'
      write ( *, '(a)' ) 
     &  '  LAGRANGE_INTERP_ND_GRID2 sets the interpolant.'
      write ( *, '(a)' ) '  LAGRANGE_INTERP_ND_VALUE2 evaluates it.'

      do i = 1, m
        ind(i) = 2
      end do

      do i = 1, m
        ab(i,1) = 0.0D+00
        ab(i,2) = 1.0D+00
      end do

      call lagrange_interp_nd_size2 ( m, ind, nd )

      call lagrange_interp_nd_grid2 ( m, ind, ab, nd, xd )
      call f_sinr ( m, nd, xd, zd )
c
c  Evaluate.
c
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) '         Zinterp          Zexact      Error'
      write ( *, '(a)' ) ''

      seed = 123456789
      call r8mat_uniform_01 ( m, ni, seed, xi )

      call f_sinr ( m, ni, xi, ze )

      call lagrange_interp_nd_value2 ( m, ind, ab, nd, zd, ni, xi, zi )

      do j = 1,  ni
        write ( *, '(2x,g14.6,2x,g14.6,2x,e10.2)' ) 
     &    zi(j), ze(j), abs ( zi(j) - ze(j) ) 
      end do

      return
      end
      subroutine test08 ( )

c*********************************************************************72
c
cc TEST08 repeats test 4 using levels.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 3 )
      integer ni
      parameter ( ni = 20 )

      double precision ab(m,2)
      double precision e
      external f_sinr
      integer i
      integer ind(m)
      integer j
      integer l
      integer nd
      integer order
      double precision r8vec_norm_affine
      integer seed
      double precision xd(m,35937)
      double precision xi(m,ni)
      double precision zd(35937)
      double precision ze(ni)
      double precision zi(ni)

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST08:'
      write ( *, '(a)' ) '  Interpolate in 3D, using levels.'
      write ( *, '(a)' ) '  Use a sequence of increasing levels.'

      do i = 1, m
        ab(i,1) = 0.0D+00
        ab(i,2) = 1.0D+00
      end do

      seed = 123456789
      call r8mat_uniform_01 ( m, ni, seed, xi )

      call f_sinr ( m, ni, xi, ze )

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) '  Level     Order   Average Error'
      write ( *, '(a)' ) ''

      do l = 0, 5

        do i = 1, m
          ind(i) = l
        end do
    
        call lagrange_interp_nd_size2 ( m, ind, nd )

        call lagrange_interp_nd_grid2 ( m, ind, ab, nd, xd )
        call f_sinr ( m, nd, xd, zd )
c
c  Evaluate.
c
        call lagrange_interp_nd_value2 ( m, ind, ab, nd, zd, ni, 
     &    xi, zi )

        e = r8vec_norm_affine ( ni, zi, ze ) / dble ( ni )

        write ( *, '(2x,i5,2x,i8,2x,e10.2)' ) l, nd, e

      end do

      return
      end
      subroutine test09 ( )

c*********************************************************************72
c
cc TEST09 interpolates in 3D, using anisotropic resolution.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 3 )
      integer ni
      parameter ( ni = 20 )

      double precision ab(m,2)
      double precision e
      external f_poly352
      integer i
      integer j
      integer l
      integer n_1d(m)
      integer nd
      integer order
      double precision r8vec_norm_affine
      integer seed
      double precision xd(m,96)
      double precision xi(m,ni)
      double precision zd(96)
      double precision ze(ni)
      double precision zi(ni)

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'TEST09:'
      write ( *, '(a)' ) '  Interpolate in 3D, using orders.'
      write ( *, '(a)' ) '  Use a sequence of increasing orders.'
      write ( *, '(a)' ) '  Use anisotropic resolution.'
      write ( *, '(a)' ) 
     &  '  The interpoland is a polynomial of degrees 3, 5, 2'
      write ( *, '(a)' ) 
     &  '  so our orders need to be at least 4, 6, 3 to match it.'

      do i = 1, m
        ab(i,1) = 0.0D+00
        ab(i,2) = 1.0D+00
      end do

      seed = 123456789
      call r8mat_uniform_01 ( m, ni, seed, xi )

      call f_poly352 ( m, ni, xi, ze )

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) '  Level     Orders   Average Error'
      write ( *, '(a)' ) ''

      do l = 0, 10

        if ( l .eq. 0 ) then
          n_1d(1) = 1
          n_1d(2) = 1
          n_1d(3) = 1
        else if ( l .eq. 1 ) then
          n_1d(1) = 2
          n_1d(2) = 1
          n_1d(3) = 1
        else if ( l .eq. 2 ) then
          n_1d(1) = 1
          n_1d(2) = 2
          n_1d(3) = 1
        else if ( l .eq. 3 ) then
          n_1d(1) = 1
          n_1d(2) = 1
          n_1d(3) = 2
        else if ( l .eq. 4 ) then
          n_1d(1) = 4
          n_1d(2) = 2
          n_1d(3) = 2
        else if ( l .eq. 5 ) then
          n_1d(1) = 2
          n_1d(2) = 4
          n_1d(3) = 2
        else if ( l .eq. 6 ) then
          n_1d(1) = 2
          n_1d(2) = 2
          n_1d(3) = 4
        else if ( l .eq. 8 ) then
          n_1d(1) = 6
          n_1d(2) = 4
          n_1d(3) = 4
        else if ( l .eq. 9 ) then
          n_1d(1) = 4
          n_1d(2) = 6
          n_1d(3) = 4
        else if ( l .eq. 10 ) then
          n_1d(1) = 4
          n_1d(2) = 4
          n_1d(3) = 6
        end if
      
        call lagrange_interp_nd_size ( m, n_1d, nd )

        call lagrange_interp_nd_grid ( m, n_1d, ab, nd, xd )
        call f_poly352 ( m, nd, xd, zd )
c
c  Evaluate.
c
        call lagrange_interp_nd_value ( m, n_1d, ab, nd, zd, ni, 
     &    xi, zi )

        e = r8vec_norm_affine ( ni, zi, ze ) / dble ( ni )

        write ( *, '(2x,i5,2x,i5,2x,i5,2x,i5,2x,e10.2)' ) 
     &    l, n_1d(1:3), e

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
c    28 September 2012
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
          r = r + x(i,j)**2
        end do
        r = sqrt ( r ) 
        z(j) = sin ( r )
      end do

      return
      end
      subroutine f_poly352 ( m, n, x, z )

c*********************************************************************72
c
cc F_POLY253 is a scalar function of a 3-dimensional argument, to be interpolated.
c
c  Discussion:
c
c    The polynomial is of maximum degrees 3, 5, and 2, in the three variables.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 September 2012
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

      integer j
      double precision x(m,n)
      double precision z(n)

      do j = 1, n
        z(j) = 1.0 + x(1,j) ** 2 * x(2,j) ** 5 * x(3,j) ** 2 
     &    + x(1,j) * x(2,j) ** 2 * x(3,j) + x(1,j) ** 3
      end do

      return
      end
