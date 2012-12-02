      program main

c*********************************************************************72
c
cc MAIN is the main program for TEST_INT_2D_PRB.
c
c  Discussion:
c
c    TEST_INT_2D_PRB demonstrates the TEST_INT_2D integration test functions.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_INT_2D_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TEST_INT_2D library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_INT_2D_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 applies a Monte Carlo rule.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a(2)
      double precision b(2)
      integer dim
      double precision error
      double precision exact
      double precision fx(4194304)
      integer i
      integer j
      integer n
      integer problem
      integer problem_num
      double precision quad
      double precision r8vec_sum
      integer seed
      double precision volume
      double precision x(2,4194304)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Use a Monte Carlo rule.'
      write ( *, '(a)' ) 
     &  '  Repeatedly multiply the number of points by 4.'

      call p00_problem_num ( problem_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '   Problem      Points       Approx         Error'

      do problem = 1, problem_num

        write ( *, '(a)' ) ' '

        n = 1

        do i = 1, 12

          seed = 123456789

          call r8mat_uniform_01 ( 2, n, seed, x )

          call p00_lim ( problem, a, b )

          do j = 1, n
            do dim = 1, 2
              x(dim,j) = ( 1.0D+00 - x(dim,j) ) * a(dim) 
     &                   +           x(dim,j)   * b(dim)
            end do
          end do

          volume = ( b(1) - a(1) ) * ( b(2) - a(2) )

          call p00_fun ( problem, n, x, fx )

          quad = volume * r8vec_sum ( n, fx ) / dble ( n )

          call p00_exact ( problem, exact )

          error = abs ( quad - exact )

          write ( *, '(2x,i8,2x,i10,2x,g14.6,2x,g14.6)' ) 
     &      problem, n, quad, error

          n = n * 4

        end do

        write ( *, '(2x,i8,2x,a10,2x,g14.6)' ) 
     &    problem, '     Exact', exact

      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 applies a product of composite midpoint rules.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a(2)
      double precision b(2)
      integer dim
      double precision error
      double precision exact
      double precision fx(4194304)
      integer i
      integer ix
      integer iy
      integer k
      integer n
      integer nx
      integer ny
      integer problem
      integer problem_num
      double precision quad
      double precision r8vec_sum
      double precision volume
      double precision x(2,4194304)
      double precision xval
      double precision yval

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Use a product of composite midpoint rules.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Repeatedly multiply the number of points by 4.'

      call p00_problem_num ( problem_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '   Problem      Points       Approx         Error'

      do problem = 1, problem_num

        write ( *, '(a)' ) ' '

        nx = 1
        ny = 1

        do i = 1, 12

          n = nx * ny

          call p00_lim ( problem, a, b )

          k = 0

          do ix = 1, nx

            xval = ( dble ( 2 * nx - 2 * ix + 1 ) * a(1)   
     &             + dble (          2 * ix - 1 ) * b(1) ) 
     &             / dble ( 2 * nx              )

            do iy = 1, ny

              yval = ( dble ( 2 * ny - 2 * iy + 1 ) * a(2)   
     &               + dble (          2 * iy - 1 ) * b(2) ) 
     &               / dble ( 2 * ny              )

              k = k + 1
              x(1,k) = xval
              x(2,k) = yval

            end do

          end do

          volume = ( b(1) - a(1) ) * ( b(2) - a(2) )

          call p00_fun ( problem, n, x, fx )

          quad = volume * r8vec_sum ( n, fx ) / dble ( n )

          call p00_exact ( problem, exact )

          error = abs ( quad - exact )

          write ( *, '(2x,i8,2x,i10,2x,g14.6,2x,g14.6)' ) 
     &      problem, n, quad, error

          nx = nx * 2
          ny = ny * 2

        end do

        write ( *, '(2x,i8,2x,a10,2x,g14.6)' ) 
     &    problem, '     Exact', exact

      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 applies a product of Gauss-Legendre rules.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a(2)
      double precision b(2)
      integer dim
      double precision error
      double precision exact
      double precision fxy(4194304)
      integer i
      integer ix
      integer iy
      integer k
      integer nx
      integer nxy
      integer ny
      integer problem
      integer problem_num
      double precision quad
      double precision r8vec_dot_product
      double precision volume
      double precision w(2048)
      double precision wxy(4194304)
      double precision x(2048)
      double precision xy(2,4194304)
      double precision xval
      double precision yval

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  Use a product of Gauss-Legendre rules.'
      write ( *, '(a)' ) '  The 1D rules essentially double in order.'

      call p00_problem_num ( problem_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '   Problem      Points       Approx         Error'

      do problem = 1, problem_num

        write ( *, '(a)' ) ' '

        nx = 1
        ny = 1

        do i = 1, 8

          call legendre_set ( nx, x, w )

          nxy = nx * ny

          call p00_lim ( problem, a, b )

          k = 0

          do ix = 1, nx

            xval = ( ( 1.0D+00 + x(ix) ) * a(1)   
     &             + ( 1.0D+00 - x(ix) ) * b(1) ) 
     &             /   2.0D+00

            do iy = 1, ny

              yval = ( ( 1.0D+00 + x(iy) ) * a(2)   
     &               + ( 1.0D+00 - x(iy) ) * b(2) ) 
     &               /   2.0D+00

              k = k + 1
              xy(1,k) = xval
              xy(2,k) = yval
              wxy(k) = w(ix) * w(iy)

            end do

          end do

          volume = ( b(1) - a(1) ) * ( b(2) - a(2) )

          call p00_fun ( problem, nxy, xy, fxy )

          quad = volume * r8vec_dot_product ( nxy, wxy, fxy ) / 4.0D+00

          call p00_exact ( problem, exact )

          error = abs ( quad - exact )

          write ( *, '(2x,i8,2x,i10,2x,g14.6,2x,g14.6)' ) 
     &      problem, nxy, quad, error

          nx = 2 * nx + 1
          ny = nx

        end do

        write ( *, '(2x,i8,2x,a10,2x,g14.6)' ) 
     &    problem, '     Exact', exact

      end do

      return
      end
