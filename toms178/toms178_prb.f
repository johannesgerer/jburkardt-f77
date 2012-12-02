      program main

c*********************************************************************72
c
cc MAIN is the main program for TOMS178_PRB.
c
c  Discussion:
c
c    TOMS178_PRB calls sample problems for the TOMS178 library.
c
c  Modified:
c
c    12 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS178_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS178 library.'

      call test01 ( )
      call test02 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS178_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests HOOKE with the Rosenbrock function.
c
c  Modified:
c
c    12 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nvars
      parameter ( nvars = 2 )

      double precision endpt(nvars)
      double precision eps
      integer hooke
      integer i
      integer it
      integer itermax
      double precision rho
      double precision rosenbrock
      external rosenbrock
      double precision startpt(nvars)
      double precision value

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  HOOKE seeks a minimizer of F(X).'
      write ( *, '(a)' ) '  Here we use the Rosenbrock function.'
c
c  Starting guess for Rosenbrock.
c
      startpt(1) = -1.2D+00
      startpt(2) = 1.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Initial estimate X = '
      write ( *, '(a)' ) ' '
      do i = 1, nvars
        write ( *, '(2x,i8,2x,g14.6)' ) i, startpt(i)
      end do

      value = rosenbrock ( startpt, nvars )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X) = ', value
c
c  Call HOOKE.
c
      itermax = 5000
      rho = 0.5D+00
      eps = 1.0D-06

      it = hooke ( nvars, startpt, endpt, rho, eps, itermax,
     &  rosenbrock )
c
c  Results.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of iterations taken = ', it
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X* = '
      write ( *, '(a)' ) ' '
      do i = 1, nvars
        write ( *, '(2x,i8,2x,g14.6)' ) i, endpt(i)
      end do

      value = rosenbrock ( endpt, nvars )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X*) = ', value

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests HOOKE with the WOODS function.
c
c  Discussion:
c
c    The Hooke and Jeeves algorithm works well when RHO = 0.5, but
c    does poorly when RHO = 0.6, and better when RHO = 0.8
c
c  Modified:
c
c    12 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nvars
      parameter ( nvars = 4 )

      double precision endpt(nvars)
      double precision eps
      integer hooke
      integer i
      integer it
      integer itermax
      double precision rho
      double precision startpt(nvars)
      double precision value
      double precision woods
      external woods

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  HOOKE seeks a minimizer of F(X).'
      write ( *, '(a)' ) '  Here we use the Woods function.'
c
c  Starting guess.
c
      startpt(1) = -3.0D+00
      startpt(2) = -1.0D+00
      startpt(3) = -3.0D+00
      startpt(4) = -1.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Initial estimate X = '
      write ( *, '(a)' ) ' '
      do i = 1, nvars
        write ( *, '(2x,i8,2x,g14.6)' ) i, startpt(i)
      end do

      value = woods ( startpt, nvars )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X) = ', value
c
c  Call HOOKE.
c
      itermax = 5000
      rho = 0.5D+00
      eps = 1.0D-06

      it = hooke ( nvars, startpt, endpt, rho, eps, itermax,
     &  woods )
c
c  Results.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of iterations taken = ', it
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  X* = '
      write ( *, '(a)' ) ' '
      do i = 1, nvars
        write ( *, '(2x,i8,2x,g14.6)' ) i, endpt(i)
      end do

      value = woods ( endpt, nvars )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  F(X*) = ', value

      return
      end
      function rosenbrock ( x, n )

c*********************************************************************72
c
cc ROSENBROCK evaluates the Rosenbrock function.
c
c  Discussion:
c
c    The Hooke and Jeeves algorithm works reasonably well on
c    Rosenbrock's test function, depending on the value of RHO chosen.
c
c  Modified:
c
c    12 February 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X(N), the argument of the function.
c
c    Input, integer N, the spatial dimension.
c
c    Output, double precision ROSENBROCK, the value of the function.
c
      implicit none

      integer n

      double precision rosenbrock
      double precision x(n)

      rosenbrock = 100.0 * ( x(2) - x(1) * x(1) )**2
     &           +         ( 1.0D+00 - x(1) )**2

      return
      end
      function woods ( x, n )

c*********************************************************************72
c
cc WOODS evaluates the Woods function.
c
c  Modified:
c
c    12 February 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X(N), the argument of the function.
c
c    Input, integer N, the spatial dimension.
c
c    Output, double precision WOODS, the value of the function.
c
      implicit none

      integer n

      double precision s1
      double precision s2
      double precision s3
      double precision t1
      double precision t2
      double precision t3
      double precision t4
      double precision t5
      double precision woods
      double precision x(n)

      s1 = x(2) - x(1) * x(1)
      s2 = 1.0D+00 - x(1)
      s3 = x(2) - 1.0D+00
      t1 = x(4) - x(3) * x(3)
      t2 = 1.0D+00 - x(3)
      t3 = x(4) - 1.0D+00
      t4 = s3 + t3
      t5 = s3 - t3

      woods = 100.0D+00 * s1**2
     &      +             s2**2
     &      +  90.0D+00 * t1**2
     &      +             t2**2
     &      +  10.0D+00 * t4**2
     &      +   0.1D+00 * t5**2

      return
      end

