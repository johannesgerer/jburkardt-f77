      program main

c*********************************************************************72
c
cc MAIN is the main program for TOMS502_PRB.
c
c  Modified:
c
c    03 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milan Kubicek,
c    Algorithm 502:
c    Dependence of solution of nonlinear systems on a parameter,
c    ACM Transactions on Mathematical Software,
c    Volume 2, Number 1, March 1976, pages 98-107.
c
      implicit none

      integer n
      parameter ( n = 4 )

      real e
      real eps
      real hh
      real hmax
      integer i
      integer inital
      integer itin
      integer maxout
      integer mxadms
      integer ncorr
      integer ncrad
      integer ndir(n+1)
      integer nout
      integer nprnt
      real out(100,n+2)
      real pref(n+1)
      real w(n+1)
      real x(n+1)
      real xlow(n+1)
      real xupp(n+1)

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS502_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS502 library'

      do i = 1, n + 1
        x(i) = 0.0E+00
      end do

      do i = 1, n + 1
        xlow(i) = -1.0E+00
      end do

      do i = 1, n + 1
        xupp(i) = 1.0E+00
      end do

      eps = 0.0001E+00

      do i = 1, n + 1
        w(i) = 1.0E+00
      end do

      inital = 1

      itin = 10

      hh = 0.02E+00

      hmax = 0.10E+00

      pref(1) = 1.0E+00
      pref(2) = 0.1E+00
      pref(3) = 1.0E+00
      pref(4) = 1.0E+00
      pref(5) = 0.02E+00

      do i = 1, n + 1
        ndir(i) = +1
      end do

      e = 0.0001E+00

      mxadms = 4

      ncorr = 4

      ncrad = 0

      nprnt = 2

      call derpar ( n, x, xlow, xupp, eps, w, inital, itin, hh,
     &  hmax, pref, ndir, e, mxadms, ncorr, ncrad, nout, out,
     &  maxout, nprnt )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS502_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine fctn ( n, x, f, g )

c*********************************************************************72
c
cc FCTN evaluates the function and jacobian.
c
c  Modified:
c
c    03 September 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Milan Kubicek,
c    Algorithm 502:
c    Dependence of solution of nonlinear systems on a parameter,
c    ACM Transactions on Mathematical Software,
c    Volume 2, Number 1, March 1976, pages 98-107.
c
c  Parameters:
c
c    Input, integer N, the number of equations, and the number of
c    variables (not including the parameter).
c
c    Input, real X(N+1), the current values of the N variables, and the
c    parameter.  The parameter is stored in X(N+1).
c
c    Output, real F(N), the equations, evaluated at X.
c
c    Output, real G(N,N+1), the jacobian matrix.
c    G(I,J) = dF(I)/dX(J).
c
      implicit none

      integer n

      real alpha
      real arg1
      real arg2
      real b
      parameter ( b = 22.0E+00 )
      real beta1
      parameter ( beta1 = 2.0E+00 )
      real beta2
      parameter ( beta2 = 2.0E+00 )
      real darg1
      real darg2
      real f(n)
      real g(n,n+1)
      real gamma
      parameter ( gamma = 1000.0E+00 )
      real theta1
      parameter ( theta1 = 0.0E+00 )
      real theta2
      parameter ( theta2 = 0.0E+00 )
      real x(n+1)

      alpha = x(n+1)

      arg1 = 10.0E+00 * x(1) / ( 1.0E+00 + 10.0 * x(1) / gamma )
      arg2 = 10.0E+00 * x(2) / ( 1.0E+00 + 10.0 * x(2) / gamma )
c
c  Evaluate the function.
c
      f(1) = alpha * ( 1.0E+00 - x(3) ) * exp ( arg1 ) - x(3)

      f(2) = alpha * b * ( 1.0E+00 - x(3) ) * exp ( arg1 )
     &  + beta1 * theta1 - 10.0E+00 * ( 1.0 + beta1 ) * x(1)

      f(3) = x(3) - x(4) + alpha * ( 1.0E+00 - x(4) ) * exp ( arg2 )

      f(4) = 10.0E+00 * x(1) - 10.0E+00 * ( 1.0E+00 + beta2 ) * x(2)
     &  + alpha * b * ( 1.0E+00 - x(4) ) * exp ( arg2 )
     &  + beta2 * theta2
c
c  Evaluate the jacobian.
c
      darg1 = 10.0E+00 / ( 1.0E+00 + 10.0 * x(1) / gamma )**2
      darg2 = 10.0E+00 / ( 1.0E+00 + 10.0 * x(2) / gamma )**2

      g(1,1) = alpha * ( 1.0E+00 - x(3) ) * exp ( arg1 ) * darg1
      g(1,2) = 0.0E+00
      g(1,3) = - alpha * exp ( arg1 ) - 1.0E+00
      g(1,4) = 0.0E+00
      g(1,5) = ( 1.0E+00 - x(3) ) * exp ( arg1 )

      g(2,1) = alpha * b * ( 1.0E+00 - x(3) ) * exp ( arg1 ) * darg1
     &  - 10.0E+00 * ( 1.0 + beta1 )
      g(2,2) = 0.0E+00
      g(2,3) = - alpha * b * exp ( arg1 )
      g(2,4) = 0.0E+00
      g(2,5) = b * ( 1.0E+00 - x(3) ) * exp ( arg1 )

      g(3,1) = 0.0E+00
      g(3,2) = alpha * ( 1.0E+00 - x(4) ) * exp ( arg2 ) * darg2
      g(3,3) = 1.0E+00
      g(3,4) = - 1.0E+00 - alpha * exp ( arg2 )
      g(3,5) = ( 1.0E+00 - x(4) ) * exp ( arg2 )

      g(4,1) = 10.0E+00
      g(4,2) = - 10.0E+00 * ( 1.0E+00 + beta2 )
     &  + alpha * b * ( 1.0E+00 - x(4) ) * exp ( arg2 ) * darg2
      g(4,3) = 0.0E+00
      g(4,4) = - alpha * b * exp ( arg2 )
      g(4,5) = b * ( 1.0E+00 - x(4) ) * exp ( arg2 )

      return
      end
