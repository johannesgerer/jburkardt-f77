      subroutine rk4 ( t0, u0, dt, f, u1 )

c*********************************************************************72
c
cc RK4 takes one Runge-Kutta step.
c
c  Discussion:
c
c    It is assumed that an initial value problem, of the form
c
c      du/dt = f ( t, u )
c      u(t0) = u0
c
c    is being solved.
c
c    If the user can supply current values of t, u, a stepsize dt, and a
c    function to evaluate the derivative, this function can compute the
c    fourth-order Runge Kutta estimate to the solution at time t+dt.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    31 January 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T0, the current time.
c
c    Input, double precision U0, the solution estimate at the current time.
c
c    Input, double precision DT, the time step.
c
c    Input, external F, a subroutine of the form 
c      subroutine f ( t, u, uprime ) 
c    which evaluates the derivative uprime given the time T and
c    solution vector U.
c
c    Output, double precision U0, the fourth-order Runge-Kutta solution 
c    estimate at time T0+DT.
c
      implicit none

      double precision dt
      external f
      double precision f1
      double precision f2
      double precision f3
      double precision f4
      double precision t0
      double precision u0
      double precision u1
c
c  Get four sample values of the derivative.
c
      call f ( t0,                u0,                     f1 )
      call f ( t0 + dt / 2.0D+00, u0 + dt * f1 / 2.0D+00, f2 )
      call f ( t0 + dt / 2.0D+00, u0 + dt * f2 / 2.0D+00, f3 )
      call f ( t0 + dt,           u0 + dt * f3,           f4 )
c
c  Combine them to estimate the solution U1 at time T1 = T0 + DT.
c
      u1 = u0 + dt * ( f1 + 2.0D+00 * f2 + 2.0D+00 * f3 + f4 ) / 6.0D+00

      return
      end
