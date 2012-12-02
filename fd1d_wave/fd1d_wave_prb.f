      program main

c*********************************************************************72
c
cc FD1D_WAVE_PRB tests the FD1D finite difference wave computation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FD1D_WAVE_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the FD1D_WAVE library.'

      call fd1d_wave_test01 ( )
      call fd1d_wave_test02 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FD1D_WAVE_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      return
      end
      subroutine fd1d_wave_test01 ( )

c*********************************************************************72
c
cc FD1D_WAVE_TEST_01 tests the FD1D finite difference wave computation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer t_num
      parameter ( t_num = 41 )
      integer x_num
      parameter ( x_num = 16 )

      double precision alpha
      double precision c
      integer i
      integer j
      double precision t
      double precision t_delta
      double precision t_vec(t_num)
      double precision t1
      double precision t2
      double precision u(t_num,x_num)
      external u_x1_01
      external u_x2_01
      external u_t1_01
      external ut_t1_01
      double precision u1(x_num)
      double precision u2(x_num)
      double precision u3(x_num)
      double precision x_vec(x_num)
      double precision x1
      double precision x2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FD1D_WAVE_TEST01'
      write ( *, '(a)' ) '  Try the "shark" wave.'

      x1 = 0.0D+00
      x2 = 1.5D+00
      call r8vec_linspace ( x_num, x1, x2, x_vec )

      t1 = 0.0D+00
      t2 = 4.0D+00
      call r8vec_linspace ( t_num, t1, t2, t_vec )
      t_delta = ( t2 - t1 ) / dble ( t_num - 1 )

      c = 1.0D+00
      call fd1d_wave_alpha ( x_num, x1, x2, t_num, t1, t2, c, alpha )
c
c  Load the initial condition.
c
      call u_t1_01 ( x_num, x_vec, u1 )
      do j = 1, x_num
        u(1,j) = u1(j)
      end do
c
c  Take the first step.
c
      t = t_vec(2)
      call fd1d_wave_start ( x_num, x_vec, t, t_delta, alpha, u_x1_01, 
     &  u_x2_01, ut_t1_01, u1, u2 )
      do j = 1, x_num
        u(2,j) = u2(j)
      end do
c
c  Take all the other steps.
c
      do i = 3, t_num
        t = t_vec(i)
        call fd1d_wave_step ( x_num, t, alpha, u_x1_01, u_x2_01, u1, 
     &    u2, u3 )
        do j = 1, x_num
          u(i,j) = u3(j)
          u1(j) = u2(j)
          u2(j) = u3(j)
        end do
      end do
c
c  Write the solution to a file.
c
      call r8mat_write ( 'test01_plot.txt', t_num, x_num, u )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Plot data written to "test01_plot.txt".'

      return
      end
      subroutine u_x1_01 ( t, u )

c*********************************************************************72
c
cc U_X1_01 evaluates U at the boundary X1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T, the time.
c
c    Output, double precision U, the value of U(T,X1).
c
      implicit none

      integer nd
      parameter ( nd = 6 )

      double precision t
      double precision td(nd)
      double precision u
      double precision ud(nd)

      save td
      save ud

      data td /
     &  0.0D+00, 0.10D+00, 0.20D+00, 0.30D+00, 0.40D+00, 0.50D+00 /
      data ud /
     &  0.0D+00, 2.0D+00, 10.0D+00, 8.0D+00, 5.0D+00, 0.0D+00 /

      call piecewise_linear ( nd, td, ud, 1, t, u )

      return
      end
      subroutine u_x2_01 ( t, u )

c*********************************************************************72
c
cc U_X2_01 evaluates U at the boundary X2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T, the time.
c
c    Output, double precision U, the value of U(T,X2).
c
      implicit none

      double precision t
      double precision u

      u = 0.0D+00

      return
      end
      subroutine u_t1_01 ( x_num, x_vec, u )

c*********************************************************************72
c
cc U_T1_01 evaluates U at the initial time T1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer X_NUM, the number of nodes.
c
c    Input, double precision X_VEC(X_NUM), the coordinates of the nodes.
c
c    Output, double precision U(X_NUM), the value of U at the initial time. 
c
      implicit none

      integer x_num

      integer j
      double precision u(x_num)
      double precision x_vec(x_num)

      do j = 1, x_num
        u(j) = 0.0D+00
      end do

      return
      end
      subroutine ut_t1_01 ( x_num, x_vec, ut )

c*********************************************************************72
c
cc UT_T1_01 evaluates dUdT at the initial time T1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer X_NUM, the number of nodes.
c
c    Input, double precision X_VEC(X_NUM), the coordinates of the nodes.
c
c    Output, double precision UT(X_NUM), the value of dUdT at the initial time. 
c
      implicit none

      integer x_num

      integer j
      double precision ut(x_num)
      double precision x_vec(x_num)

      do j = 1, x_num
        ut(j) = 0.0D+00
      end do

      return
      end
      subroutine fd1d_wave_test02 ( )

c*********************************************************************72
c
cc FD1D_WAVE_TEST_02 tests the FD1D finite difference wave computation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer t_num 
      parameter ( t_num = 41 )
      integer x_num
      parameter ( x_num = 16 )

      double precision alpha
      double precision c
      integer i
      integer j
      double precision t
      double precision t_delta
      double precision t_vec(t_num)
      double precision t1
      double precision t2
      double precision u(t_num,x_num)
      double precision u1(x_num)
      double precision u2(x_num)
      double precision u3(x_num)
      external u_x1_02
      external u_x2_02
      external u_t1_02
      external ut_t1_02
      double precision x_vec(x_num)
      double precision x1
      double precision x2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FD1D_WAVE_TEST02'
      write ( *, '(a)' ) '  Try a sine curve.'

      x1 = 0.0D+00
      x2 = 1.5D+00
      call r8vec_linspace ( x_num, x1, x2, x_vec )

      t1 = 0.0D+00
      t2 = 4.0D+00
      call r8vec_linspace ( t_num, t1, t2, t_vec )
      t_delta = ( t2 - t1 ) / dble ( t_num - 1 )
c
c  Changing T2 to 4.5 is enough to push the algorithm into instability.
c
c  t2 = 4.5D+00
c
      c = 1.0D+00
      call fd1d_wave_alpha ( x_num, x1, x2, t_num, t1, t2, c, alpha )
c
c  Load the initial condition.
c
      call u_t1_02 ( x_num, x_vec, u1 )
      do j = 1, x_num
        u(1,j) = u1(j)
      end do
c
c  Take the first step.
c
      t = t_vec(2)
      call fd1d_wave_start ( x_num, x_vec, t, t_delta, alpha, u_x1_02, 
     &  u_x2_02, ut_t1_02, u1, u2 )
      do j = 1, x_num
        u(2,j) = u2(j)
      end do
c
c  Take all the other steps.
c
      do i = 3, t_num
        t = t_vec(i)
        call fd1d_wave_step ( x_num, t, alpha, u_x1_02, u_x2_02, u1, 
     &    u2, u3 )
        do j = 1, x_num
          u(i,j) = u3(j)
          u1(j) = u2(j)
          u2(j) = u3(j)
        end do
      end do
c
c  Write the solution to a file.
c
      call r8mat_write ( 'test02_plot.txt', t_num, x_num, u )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Plot data written to "test02_plot.txt".'

      return
      end
      subroutine u_x1_02 ( t, u )

c*********************************************************************72
c
cc U_X1_02 evaluates U at the boundary X1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T, the time.
c
c    Output, double precision U, the value of U(T,X1).
c
      implicit none

      double precision t
      double precision u

      u = 0.0D+00

      return
      end
      subroutine u_x2_02 ( t, u )

c*********************************************************************72
c
cc U_X2_02 evaluates U at the boundary X2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T, the time.
c
c    Output, double precision U, the value of U(T,X2).
c
      implicit none

      double precision t
      double precision u

      u = 0.0D+00

      return
      end
      subroutine u_t1_02 ( x_num, x_vec, u )

c*********************************************************************72
c
cc U_T1_02 evaluates U at the initial time T1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer X_NUM, the number of nodes.
c
c    Input, double precision X_VEC(X_NUM), the spatial node coordinates.
c
c    Output, double precision U(X_NUM), the value of U at the initial time,
c    and every node.
c
      implicit none

      integer x_num

      integer j
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision u(x_num)
      double precision x_vec(x_num)

      do j = 1, x_num
        u(j) = sin ( 2.0D+00 * pi * x_vec(j) )
      end do

      return
      end
      subroutine ut_t1_02 ( x_num, x_vec, ut )

c*********************************************************************72
c
cc UT_T1_02 evaluates dUdT at the initial time T1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer X_NUM, the number of spatial intervals.
c
c    Input, double precision X_VEC(X_NUM), the spatial node coordinates.
c
c    Output, double precision UT(X_NUM), the value of dUdT at the initial time,
c    and every node.
c
      implicit none

      integer x_num

      integer j
      double precision ut(x_num)
      double precision x_vec(x_num)

      do j = 1, x_num
        ut(j) = 0.0D+00
      end do

      return
      end
