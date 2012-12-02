      subroutine p00_a ( prob, m, n, a )

c*********************************************************************72
c
cc P00_A returns the matrix A for any least squares problem.
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
c  Parameters:
c
c    Input, integer PROB, the problem index.
c
c    Input, integer M, the number of equations.
c
c    Input, integer N, the number of variables.
c
c    Output, double precision A(M,N), the matrix.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer prob

      if ( prob .eq. 1 ) then
        call p01_a ( m, n, a )
      else if ( prob .eq. 2 ) then
        call p02_a ( m, n, a )
      else if ( prob .eq. 3 ) then
        call p03_a ( m, n, a )
      else if ( prob .eq. 4 ) then
        call p04_a ( m, n, a )
      else if ( prob .eq. 5 ) then
        call p05_a ( m, n, a )
      else if ( prob .eq. 6 ) then
        call p06_a ( m, n, a )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_A - Fatal error!'
        write ( *, '(a)' ) '  Illegal value of PROB = ', prob
        stop
      end if

      return
      end
      subroutine p00_b ( prob, m, b )

c*********************************************************************72
c
cc P00_B returns the right hand side B for any least squares problem.
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
c  Parameters:
c
c    Input, integer PROB, the problem index.
c
c    Input, integer M, the number of equations.
c
c    Output, double precision B(M), the right hand side.
c
      implicit none

      integer m

      double precision b(m)
      integer prob

      if ( prob .eq. 1 ) then
        call p01_b ( m, b )
      else if ( prob .eq. 2 ) then
        call p02_b ( m, b )
      else if ( prob .eq. 3 ) then
        call p03_b ( m, b )
      else if ( prob .eq. 4 ) then
        call p04_b ( m, b )
      else if ( prob .eq. 5 ) then
        call p05_b ( m, b )
      else if ( prob .eq. 6 ) then
        call p06_b ( m, b )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_B - Fatal error!'
        write ( *, '(a)' ) '  Illegal value of PROB = ', prob
        stop
      end if

      return
      end
      subroutine p00_m ( prob, m )

c*********************************************************************72
c
cc P00_M returns the number of equations M for any least squares problem.
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
c  Parameters:
c
c    Input, integer PROB, the problem index.
c
c    Output, integer M, the number of equations.
c
      implicit none

      integer m
      integer prob

      if ( prob .eq. 1 ) then
        call p01_m ( m )
      else if ( prob .eq. 2 ) then
        call p02_m ( m )
      else if ( prob .eq. 3 ) then
        call p03_m ( m )
      else if ( prob .eq. 4 ) then
        call p04_m ( m )
      else if ( prob .eq. 5 ) then
        call p05_m ( m )
      else if ( prob .eq. 6 ) then
        call p06_m ( m )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_M - Fatal error!'
        write ( *, '(a)' ) '  Illegal value of PROB = ', prob
        stop
      end if

      return
      end
      subroutine p00_n ( prob, n )

c*********************************************************************72
c
cc P00_N returns the number of variables N for any least squares problem.
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
c  Parameters:
c
c    Input, integer PROB, the problem index.
c
c    Output, integer N, the number of variables.
c
      implicit none

      integer n
      integer prob

      if ( prob .eq. 1 ) then
        call p01_n ( n )
      else if ( prob .eq. 2 ) then
        call p02_n ( n )
      else if ( prob .eq. 3 ) then
        call p03_n ( n )
      else if ( prob .eq. 4 ) then
        call p04_n ( n )
      else if ( prob .eq. 5 ) then
        call p05_n ( n )
      else if ( prob .eq. 6 ) then
        call p06_n ( n )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_N - Fatal error!'
        write ( *, '(a)' ) '  Illegal value of PROB = ', prob
        stop
      end if

      return
      end
      subroutine p00_prob_num ( prob_num )

c*********************************************************************72
c
cc P00_PROB_NUM returns the number of least squares problems.
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
c  Parameters:
c
c    Output, double precision PROB_NUM, the number of problems.
c
      implicit none

      integer prob_num

      prob_num = 6

      return
      end
      subroutine p00_x ( prob, n, x )

c*********************************************************************72
c
cc P00_X returns the least squares solution X for any least squares problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 April 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROB, the problem index.
c
c    Input, integer N, the number of variables.
c
c    Output, double precision X(N), the least squares solution.
c
      implicit none

      integer n

      integer prob
      double precision x(n)

      if ( prob .eq. 1 ) then
        call p01_x ( n, x )
      else if ( prob .eq. 2 ) then
        call p02_x ( n, x )
      else if ( prob .eq. 3 ) then
        call p03_x ( n, x )
      else if ( prob .eq. 4 ) then
        call p04_x ( n, x )
      else if ( prob .eq. 5 ) then
        call p05_x ( n, x )
      else if ( prob .eq. 6 ) then
        call p06_x ( n, x )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_X - Fatal error!'
        write ( *, '(a)' ) '  Illegal value of PROB = ', prob
        stop
      end if

      return
      end
      subroutine p01_a ( m, n, a )

c*********************************************************************72
c
cc P01_A returns the matrix A for problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 April 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of equations.
c
c    Input, integer N, the number of variables.
c
c    Output, double precision A(M,N), the matrix.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer i
      integer j

      do i = 1, m
        a(i,1) = 1.0D+00
        do j = 2, n
          a(i,j) = a(i,j-1) * dble ( i )
        end do
      end do

      return
      end
      subroutine p01_b ( m, b )

c*********************************************************************72
c
cc P01_B returns the right hand side B for problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 April 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of equations.
c
c    Output, double precision B(M), the right hand side.
c
      implicit none

      integer m

      double precision b(m)
      double precision b_save(5)

      save b_save

      data b_save / 1.0D+00, 2.3D+00, 4.6D+00, 3.1D+00, 1.2D+00 /

      call r8vec_copy ( m, b_save, b )

      return
      end
      subroutine p01_m ( m )

c*********************************************************************72
c
cc P01_M returns the number of equations M for problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 April 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer M, the number of equations.
c
      implicit none

      integer m

      m = 5

      return
      end
      subroutine p01_n ( n )

c*********************************************************************72
c
cc P01_N returns the number of variables N for problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 April 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer N, the number of variables.
c
      implicit none

      integer n

      n = 3

      return
      end
      subroutine p01_x ( n, x )

c*********************************************************************72
c
cc P01_X returns the least squares solution X for problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 April 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of variables.
c
c    Output, double precision X(N), the least squares solution.
c
      implicit none

      integer n

      double precision x(n)
      double precision x_save(3)

      save x_save

      data x_save / -3.0200000D+00, 4.4914286D+00, -0.72857143D+00 /

      call r8vec_copy ( n, x_save, x )

      return
      end
      subroutine p02_a ( m, n, a )

c*********************************************************************72
c
cc P02_A returns the matrix A for problem 2.
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
c  Reference:
c
c    Cleve Moler,
c    Numerical Computing with MATLAB,
c    SIAM, 2004,
c    ISBN13: 978-0-898716-60-3,
c    LC: QA297.M625,
c    ebook: http://www.mathworks.com/moler/chapters.html
c
c  Parameters:
c
c    Input, integer M, the number of equations.
c
c    Input, integer N, the number of variables.
c
c    Output, double precision A(M,N), the matrix.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer i
      integer j

      do i = 1, m
        a(i,n) = 1.0D+00
        do j = n - 1, 1, -1
          a(i,j) = a(i,j+1) * dble ( i - 1 ) / 5.0D+00
        end do
      end do

      return
      end
      subroutine p02_b ( m, b )

c*********************************************************************72
c
cc P02_B returns the right hand side B for problem 2.
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
c  Parameters:
c
c    Input, integer M, the number of equations.
c
c    Output, double precision B(M), the right hand side.
c
      implicit none

      integer m

      double precision b(m)
      double precision b_save(6)

      save b_save

      data b_save / 
     &  150.697D+00, 179.323D+00, 203.212D+00, 226.505D+00, 249.633D+00, 
     &  281.422D+00 /

      call r8vec_copy ( m, b_save, b )

      return
      end
      subroutine p02_m ( m )

c*********************************************************************72
c
cc P02_M returns the number of equations M for problem 2.
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
c  Parameters:
c
c    Output, integer M, the number of equations.
c
      implicit none

      integer m

      m = 6

      return
      end
      subroutine p02_n ( n )

c*********************************************************************72
c
cc P02_N returns the number of variables N for problem 2.
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
c  Parameters:
c
c    Output, integer N, the number of variables.
c
      implicit none

      integer n

      n = 3

      return
      end
      subroutine p02_x ( n, x )

c*********************************************************************72
c
cc P02_X returns the least squares solution X for problem 2.
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
c  Parameters:
c
c    Input, integer N, the number of variables.
c
c    Output, double precision X(N), the least squares solution.
c
      implicit none

      integer n

      double precision x(n)
      double precision x_save(3)

      save x_save

      data x_save / 5.7013D+00, 121.1341D+00, 152.4745D+00 /

      call r8vec_copy ( n, x_save, x )

      return
      end
      subroutine p03_a ( m, n, a )

c*********************************************************************72
c
cc P03_A returns the matrix A for problem 3.
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
c  Reference:
c
c    Cleve Moler,
c    Numerical Computing with MATLAB,
c    SIAM, 2004,
c    ISBN13: 978-0-898716-60-3,
c    LC: QA297.M625,
c    ebook: http://www.mathworks.com/moler/chapters.html
c
c  Parameters:
c
c    Input, integer M, the number of equations.
c
c    Input, integer N, the number of variables.
c
c    Output, double precision A(M,N), the matrix.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision a_save(5,3)

      save a_save

      data a_save /
     &  1.0D+00, 4.0D+00, 7.0D+00, 10.0D+00, 13.0D+00, 
     &  2.0D+00, 5.0D+00, 8.0D+00, 11.0D+00, 14.0D+00, 
     &  3.0D+00, 6.0D+00, 9.0D+00, 12.0D+00, 15.0D+00 /

      call r8mat_copy ( m, n, a_save, a )

      return
      end
      subroutine p03_b ( m, b )

c*********************************************************************72
c
cc P03_B returns the right hand side B for problem 3.
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
c  Parameters:
c
c    Input, integer M, the number of equations.
c
c    Output, double precision B(M), the right hand side.
c
      implicit none

      integer m

      double precision b(m)
      double precision b_save(5)

      save b_save

      data b_save / 16.0D+00, 17.0D+00, 18.0D+00, 19.0D+00, 20.0D+00 /

      call r8vec_copy ( m, b_save, b )

      return
      end
      subroutine p03_m ( m )

c*********************************************************************72
c
cc P03_M returns the number of equations M for problem 3.
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
c  Parameters:
c
c    Output, integer M, the number of equations.
c
      implicit none

      integer m

      m = 5

      return
      end
      subroutine p03_n ( n )

c*********************************************************************72
c
cc P03_N returns the number of variables N for problem 3.
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
c  Parameters:
c
c    Output, integer N, the number of variables.
c
      implicit none

      integer n

      n = 3

      return
      end
      subroutine p03_x ( n, x )

c*********************************************************************72
c
cc P03_X returns the least squares solution X for problem 3.
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
c  Parameters:
c
c    Input, integer N, the number of variables.
c
c    Output, double precision X(N), the least squares solution.
c
      implicit none

      integer n

      double precision x(n)
      double precision x_save(3)

      save x_save

      data x_save / -7.5555556D+00, 0.1111111D+00, 7.7777778D+00 /

      call r8vec_copy ( n, x_save, x )

      return
      end
      subroutine p04_a ( m, n, a )

c*********************************************************************72
c
cc P04_A returns the matrix A for problem 4.
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
c  Parameters:
c
c    Input, integer M, the number of equations.
c
c    Input, integer N, the number of variables.
c
c    Output, double precision A(M,N), the matrix.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer i
      integer j

      do j = 1, n
        do i = 1, m
          a(i,j) = dble ( j ** ( i - 1 ) )
        end do
      end do

      return
      end
      subroutine p04_b ( m, b )

c*********************************************************************72
c
cc P04_B returns the right hand side B for problem 4.
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
c  Parameters:
c
c    Input, integer M, the number of equations.
c
c    Output, double precision B(M), the right hand side.
c
      implicit none

      integer m

      double precision b(m)
      double precision b_save(3)

      save b_save

      data b_save / 15.0D+00, 55.0D+00, 225.0D+00 /

      call r8vec_copy ( m, b_save, b )

      return
      end
      subroutine p04_m ( m )

c*********************************************************************72
c
cc P04_M returns the number of equations M for problem 4.
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
c  Parameters:
c
c    Output, integer M, the number of equations.
c
      implicit none

      integer m

      m = 3

      return
      end
      subroutine p04_n ( n )

c*********************************************************************72
c
cc P04_N returns the number of variables N for problem 4.
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
c  Parameters:
c
c    Output, integer N, the number of variables.
c
      implicit none

      integer n

      n = 5

      return
      end
      subroutine p04_x ( n, x )

c*********************************************************************72
c
cc P04_X returns the least squares solution X for problem 4.
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
c  Parameters:
c
c    Input, integer N, the number of variables.
c
c    Output, double precision X(N), the least squares solution.
c
      implicit none

      integer n

      double precision x(n)
      double precision x_save(5)

      save x_save

      data x_save / 1.0D+00, 2.0D+00, 3.0D+00, 4.0D+00, 5.0D+00 /

      call r8vec_copy ( n, x_save, x )

      return
      end
      subroutine p05_a ( m, n, a )

c*********************************************************************72
c
cc P05_A returns the matrix A for problem 5.
c
c  Discussion:
c
c    A is the Hilbert matrix.
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
c  Parameters:
c
c    Input, integer M, the number of equations.
c
c    Input, integer N, the number of variables.
c
c    Output, double precision A(M,N), the matrix.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer i
      integer j

      do j = 1, n
        do i = 1, m
          a(i,j) = 1.0D+00 / dble ( i + j - 1 )
        end do
      end do

      return
      end
      subroutine p05_b ( m, b )

c*********************************************************************72
c
cc P05_B returns the right hand side B for problem 5.
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
c  Parameters:
c
c    Input, integer M, the number of equations.
c
c    Output, double precision B(M), the right hand side.
c
      implicit none

      integer m

      double precision b(m)
      integer i

      b(1) = 1.0D+00
      do i = 2, m
        b(i) = 0.0D+00
      end do

      return
      end
      subroutine p05_m ( m )

c*********************************************************************72
c
cc P05_M returns the number of equations M for problem 5.
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
c  Parameters:
c
c    Output, integer M, the number of equations.
c
      implicit none

      integer m

      m = 10

      return
      end
      subroutine p05_n ( n )

c*********************************************************************72
c
cc P05_N returns the number of variables N for problem 5.
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
c  Parameters:
c
c    Output, integer N, the number of variables.
c
      implicit none

      integer n

      n = 10

      return
      end
      subroutine p05_x ( n, x )

c*********************************************************************72
c
cc P05_X returns the least squares solution X for problem 5.
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
c  Parameters:
c
c    Input, integer N, the number of variables.
c
c    Output, double precision X(N), the least squares solution.
c
      implicit none

      integer n

      integer i
      double precision r8_choose
      double precision r8_mop
      double precision x(n)

      do i = 1, n
        x(i) = r8_mop ( i + 1 ) * dble ( i ) 
     &    * r8_choose ( n + i - 1, n - 1 ) * r8_choose ( n, n - i )
      end do

      return
      end
      subroutine p06_a ( m, n, a )

c*********************************************************************72
c
cc P06_A returns the matrix A for problem 6.
c
c  Discussion:
c
c    A is a symmetric, orthogonal matrix.
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
c  Parameters:
c
c    Input, integer M, the number of equations.
c
c    Input, integer N, the number of variables.
c
c    Output, double precision A(M,N), the matrix.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision angle
      integer i
      integer j
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      do i = 1, m
        do j = 1, n
          angle = dble ( i * j ) * pi / dble ( n + 1 )
          a(i,j) = sin ( angle ) * sqrt ( 2.0D+00 / dble ( n + 1 ) )
        end do
      end do

      return
      end
      subroutine p06_b ( m, b )

c*********************************************************************72
c
cc P06_B returns the right hand side B for problem 6.
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
c  Parameters:
c
c    Input, integer M, the number of equations.
c
c    Output, double precision B(M), the right hand side.
c
      implicit none

      integer m

      double precision b(m)
      integer i

      b(1) = 1.0D+00
      do i = 2, m
        b(i) = 0.0D+00
      end do

      return
      end
      subroutine p06_m ( m )

c*********************************************************************72
c
cc P06_M returns the number of equations M for problem 6.
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
c  Parameters:
c
c    Output, integer M, the number of equations.
c
      implicit none

      integer m

      m = 10

      return
      end
      subroutine p06_n ( n )

c*********************************************************************72
c
cc P06_N returns the number of variables N for problem 6.
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
c  Parameters:
c
c    Output, integer N, the number of variables.
c
      implicit none

      integer n

      n = 10

      return
      end
      subroutine p06_x ( n, x )

c*********************************************************************72
c
cc P06_X returns the least squares solution X for problem 6.
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
c  Parameters:
c
c    Input, integer N, the number of variables.
c
c    Output, double precision X(N), the least squares solution.
c
      implicit none

      integer n

      double precision angle
      integer i
      double precision pi 
      parameter ( pi = 3.141592653589793D+00 )
      double precision x(n)

      do i = 1, n
        angle = dble ( i ) * pi / dble ( n + 1 )
        x(i) = sin ( angle ) * sqrt ( 2.0D+00 / dble ( n + 1 ) )
      end do

      return
      end
