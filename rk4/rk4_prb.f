      program main

c*********************************************************************72
c
cc RK4_PRB demonstrates the use of the RK4 one-step ODE solver.
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
c    Local, double precision DT, the time step.
c
c    Local, double precision T0, the time at which the solution is known.
c
c    Local, double precision TMAX, the maximum time at which a solution is desired.
c
c    Local, double precision U0, the estimated solution at time T0.
c
      implicit none

      double precision dt
      parameter ( dt = 0.1D+00 )
      external f3
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision t0
      double precision t1
      double precision tmax
      parameter ( tmax = 12.0D+00 * pi )
      double precision u0
      double precision u1

      t0 = 0.0D+00
      u0 = 0.5D+00

10    continue
c
c  Print (T0,U0).
c
        write ( *, '(2x,g14.6,2x,g14.6)' ) t0, u0
c
c  Stop if we've exceeded TMAX.
c
        if ( tmax .le. t0 ) then
          go to 20
        end if
c
c  Otherwise, advance to time T1, and have RK4 estimate 
c  the solution U1 there.
c
        t1 = t0 + dt
        call rk4 ( t0, u0, dt, f3, u1 )
c
c  Shift the data to prepare for another step.
c
        t0 = t1
        u0 = u1

      go to 10

20    continue

      stop
      end
      subroutine f3 ( t, u, uprime )

c*********************************************************************72
c
cc F3 evaluates the right hand side of a particular ODE.
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
c    Input, double precision T, the current time.
c
c    Input, double precision U, the current solution value.
c
c    Output, double precision UPRIME, the value of the derivative, dU/dT.
c
      implicit none

      double precision t
      double precision u
      double precision uprime
      
      uprime = u * cos ( t )
     
      return
      end
