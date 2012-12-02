      program main

c*********************************************************************72
c
cc MAIN is the main program for PCE_ODE_HERMITE_TEST.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( );
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PCE_ODE_HERMITE_TEST:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test PCE_ODE_HERMITE.'

      call pce_ode_hermite_test01 ( )
      call pce_ode_hermite_test02 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PCE_ODE_HERMITE_TEST:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine pce_ode_hermite_test01 ( )

c*********************************************************************72
c
cc PCE_ODE_HERMITE_TEST01 runs a test problem with PCE_ODE_HERMITE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer np
      parameter ( np = 4 )
      integer nt
      parameter ( nt = 200 )

      double precision alpha_mu
      double precision alpha_sigma
      integer i
      double precision t(0:nt)
      double precision tf
      double precision ti
      double precision u(0:nt,0:np)
      double precision uex(0:nt)
      double precision ui

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PCE_ODE_HERMITE_TEST01:'
      write ( *, '(a)' ) 
     &  '  Call PCE_ODE_HERMITE to compute a polynomial'
      write ( *, '(a)' ) '  chaos expansion for the ODE:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    u'' = - alpha * u,'
      write ( *, '(a)' ) '    u(0) = 1.'

      ti = 0.0D+00
      tf = 2.0D+00
      ui = 1.0D+00
      alpha_mu = 0.0D+00
      alpha_sigma = 1.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Initial time         TI = ', ti
      write ( *, '(a,g14.6)' ) '  Final time           TF = ', tf
      write ( *, '(a,i6)' ) '  Number of time steps NT = ', nt
      write ( *, '(a,g14.6)' ) '  Initial condition    UI = ', ui
      write ( *, '(a,i6)' ) '  Expansion degree     NP = ', np
      write ( *, '(a,g14.6)' ) '  E(ALPHA)       ALPHA_MU = ', alpha_mu
      write ( *, '(a,g14.6)' ) 
     &  '  STD(ALPHA)  ALPHA_SIGMA = ', alpha_sigma

      call pce_ode_hermite ( ti, tf, nt, ui, np, alpha_mu, 
     &  alpha_sigma, t, u )
c
c  Evaluate the exact expected value function.
c
      do i = 0, nt
        uex(i) = ui * exp ( t(i)**2 / 2.0D+00 )
      end do
c
c  Compare the first computed component against the exact expected value.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     i    T(i)    E(U(T(i)))' // 
     &  '       U(T(i),0)         Error'
      write ( *, '(a)' ) ' '
      do i = 0, nt, 10
        write ( *, '(2x,i4,2x,f6.3,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &    i, t(i), uex(i), u(i,0), abs ( uex(i) - u(i,0) )
      end do

      return
      end
      subroutine pce_ode_hermite_test02 ( )

c*********************************************************************72
c
cc PCE_ODE_HERMITE_TEST02 looks at convergence behavior for a fixed time.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer np_max
      parameter ( np_max = 5 )
      integer nt
      parameter ( nt = 2000 )

      double precision alpha_mu
      double precision alpha_sigma
      double precision ep(0:np_max)
      integer i
      integer np
      double precision t(0:nt)
      double precision tf
      double precision ti
      double precision u(0:nt,0:np_max)
      double precision uex(0:nt)
      double precision uexf
      double precision ui

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PCE_ODE_HERMITE_TEST02:'
      write ( *, '(a)' ) 
     &  '  Examine convergence behavior as the PCE degree increases:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    u'' = - alpha * u,'
      write ( *, '(a)' ) '    u(0) = 1.'

      ti = 0.0D+00
      tf = 2.0D+00
      ui = 1.0D+00
      alpha_mu = 0.0D+00
      alpha_sigma = 1.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Initial time         TI = ', ti
      write ( *, '(a,g14.6)' ) '  Final time           TF = ', tf
      write ( *, '(a,i6)' ) '  Number of time steps NT = ', nt
      write ( *, '(a,g14.6)' ) '  Initial condition    UI = ', ui
      write ( *, '(a,g14.6)' ) '  E(ALPHA)       ALPHA_MU = ', alpha_mu
      write ( *, '(a,g14.6)' ) 
     &  '  STD(ALPHA)  ALPHA_SIGMA = ', alpha_sigma

      uexf = ui * exp ( tf**2 / 2.0D+00 )

      do np = 0, np_max

        call pce_ode_hermite ( ti, tf, nt, ui, np, alpha_mu, 
     &    alpha_sigma, t, u )

        ep(np) = abs ( uexf - u(nt,0) )

      end do
c
c  Print error in expected value as a function of the PCE degree.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    NP     Error(NP)     Log(Error(NP))'
      write ( *, '(a)' ) ' '
      do np = 0, np_max
        write ( *, '(2x,i4,2x,g14.6,2x,g14.6)' ) 
     &    np, ep(np), log ( ep(np) )
      end do

      return
      end
