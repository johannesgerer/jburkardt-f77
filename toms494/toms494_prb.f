      program main

c*********************************************************************72
c
cc MAIN is a main program for the Burgers equation.
c
c  Discussion:
c
c    Burgers's equation:
c
c      du/dt = - u * du/dx + nu * d2u/dx2
c
c    for 0.0 <= x <= 1.0 and 0.0 <= t <= 1.1,
c
c    with nu = 0.003, and initial conditions
c
c      u(0,x) = phi(0,x),
c
c    and boundary conditions
c
c      u(t,0) = phi(t,0),
c      u(t,1) = phi(t,1),
c
c    where the "true solution" function phi(t,x) has the form
c
c      phi(t,x) = ( 0.1 *exp(a) + 0.5 * exp(b) + exp(c) ) /
c                 (      exp(a) +       exp(b) + exp(c) )
c
c    with
c
c      a = - ( x - 0.5   + 4.95 * t ) / ( 20 * nu )
c      b = - ( x - 0.5   + 0.75 * t ) / (  4 * nu )
c      c = - ( x - 0.375            ) / (  2 * nu )
c
c  Modified:
c
c    04 February 2011
c
c  Author:
c
c    Richard Sincovec, Niel Madsen
c
c  Reference:
c
c    Richard Sincovec, Niel Madsen,
c    Software for Nonlinear Partial Differential Equations,
c    ACM Transactions on Mathematical Software,
c    Volume 1, Number 3, September 1975, pages 232-260.
c
c    Richard Sincovec, Niel Madsen,
c    Algorithm 494: PDEONE,
c    ACM Transactions on Mathematical Software,
c    Volume 1, Number 3, September 1975, pages 261-263.
c
      implicit none

      real dt
      real eps
      integer i
      integer icord
      integer iflag
      integer mf
      integer nband
      integer neqn
      integer npde
      integer npts
      real t
      real tout
      real trusol
      real u(201)
      real x(201)

      common / coord / icord
      common / mesh / x
      common / sizes / npde, npts
c
c  Define problem parameters.
c
      icord = 0
      npts = 201
      npde = 1
      iflag = 1
      neqn = npde * npts
      mf = 10
      nband = 2 * npde - 1
      t = 0.0E+00
      tout = 0.0E+00
      eps = 1.0E-05
      dt = 1.0E-08
c
c  Define the mesh and initial data.
c
      do i = 1, npts
        x(i) = real ( i - 1 ) * 0.005E+00
      end do

      do i = 1, npts
        u(i) = trusol ( t, x(i) )
      end do

20    continue

      write ( *, '(a)' ) ' '
      write ( *, '(a,e10.3)' ) '  T = ', t
      write ( *, '(5e20.6)' ) ( u(i), i = 1, npts )

      if ( 1.1E+00 .le. t ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TOMS494_PRB:'
        write ( *, '(a)' ) '  Normal end of execution.'
        stop
      end if

      tout = tout + 0.1E+00

      call driveb ( neqn, t, dt, u, tout, eps, mf, iflag, nband, nband )

      if ( iflag .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TOMS494_PRB - Fatal error!'
        write ( *, '(a)' ) '  IFLAG = ', iflag
        stop
      end if

      go to 20

      end
      subroutine d ( t, x, u, dval, npde )

c*********************************************************************72
c
cc D evaluates the diffusion coefficients.
c
c  Modified:
c
c    04 February 2011
c
c  Author:
c
c    Richard Sincovec, Niel Madsen
c
c  Reference:
c
c    Richard Sincovec, Niel Madsen,
c    Software for Nonlinear Partial Differential Equations,
c    ACM Transactions on Mathematical Software,
c    Volume 1, Number 3, September 1975, pages 232-260.
c
c    Richard Sincovec, Niel Madsen,
c    Algorithm 494: PDEONE,
c    ACM Transactions on Mathematical Software,
c    Volume 1, Number 3, September 1975, pages 261-263.
c
      implicit none

      integer npde

      real dval(npde,npde)
      real t
      real u(npde)
      real x

      dval(1,1) = 0.003E+00

      return
      end
      subroutine f ( t, x, u, ux, uxx, fval, npde )

c*********************************************************************72
c
cc
c
c  Modified:
c
c    04 February 2011
c
c  Author:
c
c    Richard Sincovec, Niel Madsen
c
c  Reference:
c
c    Richard Sincovec, Niel Madsen,
c    Software for Nonlinear Partial Differential Equations,
c    ACM Transactions on Mathematical Software,
c    Volume 1, Number 3, September 1975, pages 232-260.
c
c    Richard Sincovec, Niel Madsen,
c    Algorithm 494: PDEONE,
c    ACM Transactions on Mathematical Software,
c    Volume 1, Number 3, September 1975, pages 261-263.
c
      implicit none

      integer npde

      real fval(npde)
      real t
      real u(npde)
      real ux(npde)
      real uxx(npde,npde)
      real x

      fval(1) = uxx(1,1) - u(1) * ux(1)

      return
      end
      subroutine bndry ( t, x, u, alpha, beta, gamma, npde )

c*********************************************************************72
c
cc BNDRY
c
c  Modified:
c
c    04 February 2011
c
c  Author:
c
c    Richard Sincovec, Niel Madsen
c
c  Reference:
c
c    Richard Sincovec, Niel Madsen,
c    Software for Nonlinear Partial Differential Equations,
c    ACM Transactions on Mathematical Software,
c    Volume 1, Number 3, September 1975, pages 232-260.
c
c    Richard Sincovec, Niel Madsen,
c    Algorithm 494: PDEONE,
c    ACM Transactions on Mathematical Software,
c    Volume 1, Number 3, September 1975, pages 261-263.
c
      implicit none

      integer npde

      real alpha(npde)
      real beta(npde)
      real gamma(npde)
      real t
      real trusol
      real u(npde)
      real x

      alpha(1) = 1.0E+00
      beta(1) = 0.0E+00
      gamma(1) = trusol ( t, x )

      return
      end
      function trusol ( t, x )

c*********************************************************************72
c
cc BNDRY
c
c  Modified:
c
c    04 February 2011
c
c  Author:
c
c    Richard Sincovec, Niel Madsen
c
c  Reference:
c
c    Richard Sincovec, Niel Madsen,
c    Software for Nonlinear Partial Differential Equations,
c    ACM Transactions on Mathematical Software,
c    Volume 1, Number 3, September 1975, pages 232-260.
c
c    Richard Sincovec, Niel Madsen,
c    Algorithm 494: PDEONE,
c    ACM Transactions on Mathematical Software,
c    Volume 1, Number 3, September 1975, pages 261-263.
c
      implicit none

      real a
      real anu
      real b
      real c
      real t
      real trusol
      real x

      anu = 0.003E+00

      a = - ( x - 0.5E+00 + 4.95E+00 * t ) / ( 20.0E+00 * anu )
      b = - ( x - 0.5E+00 * 0.75E+00 * t ) / (  4.0E+00 * anu )
      c = - ( x - 0.275E+00              ) / (  2.0E+00 * anu )
      a = exp ( a )
      b = exp ( b )
      c = exp ( c )
      trusol = ( 0.1E+00 * a + 0.5E+00 * b + c ) / ( a + b + c )

      return
      end
