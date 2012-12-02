      subroutine rk1_ti_step ( x, t, h, q, fi, gi, seed, xstar )

c*********************************************************************72
c
cc RK1_TI_STEP takes one step of a stochastic Runge Kutta scheme.
c
c  Discussion:
c
c    The Runge-Kutta scheme is first-order, and suitable for time-invariant
c    systems in which F and G do not depend explicitly on time.
c
c    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Jeremy Kasdin,
c    Runge-Kutta algorithm for the numerical integration of
c    stochastic differential equations,
c    Journal of Guidance, Control, and Dynamics,
c    Volume 18, Number 1, January-February 1995, pages 114-120.
c
c    Jeremy Kasdin,
c    Discrete Simulation of Colored Noise and Stochastic Processes
c    and 1/f^a Power Law Noise Generation,
c    Proceedings of the IEEE,
c    Volume 83, Number 5, 1995, pages 802-827.
c
c  Parameters:
c
c    Input, double precision X, the value at the current time.
c
c    Input, double precision T, the current time.
c
c    Input, double precision H, the time step.
c
c    Input, double precision Q, the spectral density of the input white noise.
c
c    Input, external double precision FI, the name of the deterministic
c    right hand side function.
c
c    Input, external double precision GI, the name of the stochastic
c    right hand side function.
c
c    Input/output, integer SEED, a seed for the random 
c    number generator.
c
c    Output, double precision XSTAR, the value at time T+H.
c
      implicit none

      double precision a21
      double precision fi
      external fi
      double precision gi
      external gi
      double precision h
      double precision k1
      double precision q
      double precision q1
      double precision r8_normal_01
      integer seed
      double precision t
      double precision t1
      double precision w1
      double precision x
      double precision x1
      double precision xstar

      a21 =   1.0D+00

      q1 = 1.0D+00

      t1 = t
      x1 = x
      w1 = r8_normal_01 ( seed ) * sqrt ( q1 * q / h )
      k1 = h * fi ( x1 ) + h * gi ( x1 ) * w1 

      xstar = x1 + a21 * k1

      return
      end
      subroutine rk2_ti_step ( x, t, h, q, fi, gi, seed, xstar )

c*********************************************************************72
c
cc RK2_TI_STEP takes one step of a stochastic Runge Kutta scheme.
c
c  Discussion:
c
c    The Runge-Kutta scheme is second-order, and suitable for time-invariant
c    systems.
c
c    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Jeremy Kasdin,
c    Runge-Kutta algorithm for the numerical integration of
c    stochastic differential equations,
c    Journal of Guidance, Control, and Dynamics,
c    Volume 18, Number 1, January-February 1995, pages 114-120.
c
c    Jeremy Kasdin,
c    Discrete Simulation of Colored Noise and Stochastic Processes
c    and 1/f^a Power Law Noise Generation,
c    Proceedings of the IEEE,
c    Volume 83, Number 5, 1995, pages 802-827.
c
c  Parameters:
c
c    Input, double precision X, the value at the current time.
c
c    Input, double precision T, the current time.
c
c    Input, double precision H, the time step.
c
c    Input, double precision Q, the spectral density of the input white noise.
c
c    Input, external double precision FI, the name of the deterministic
c    right hand side function.
c
c    Input, external double precision GI, the name of the stochastic
c    right hand side function.
c
c    Input/output, integer SEED, a seed for the random 
c    number generator.
c
c    Output, double precision XSTAR, the value at time T+H.
c
      implicit none

      double precision a21
      double precision a31
      double precision a32
      double precision fi
      external fi
      double precision gi
      external gi
      double precision h
      double precision k1
      double precision k2
      double precision q
      double precision q1
      double precision q2
      double precision r8_normal_01
      integer seed
      double precision t
      double precision t1
      double precision t2
      double precision w1
      double precision w2
      double precision x
      double precision x1
      double precision x2
      double precision xstar

      a21 =   1.0D+00
      a31 =   0.5D+00
      a32 =   0.5D+00

      q1 = 2.0D+00
      q2 = 2.0D+00

      t1 = t
      x1 = x
      w1 = r8_normal_01 ( seed ) * sqrt ( q1 * q / h )
      k1 = h * fi ( x1 ) + h * gi ( x1 ) * w1 

      t2 = t1 + a21 * h
      x2 = x1 + a21 * k1
      w2 = r8_normal_01 ( seed ) * sqrt ( q2 * q / h )
      k2 = h * fi ( x2 ) + h * gi ( x2 ) * w2

      xstar = x1 + a31 * k1 + a32 * k2

      return
      end
      subroutine rk3_ti_step ( x, t, h, q, fi, gi, seed, xstar )

c*********************************************************************72
c
cc RK3_TI_STEP takes one step of a stochastic Runge Kutta scheme.
c
c  Discussion:
c
c    The Runge-Kutta scheme is third-order, and suitable for time-invariant
c    systems in which F and G do not depend explicitly on time.
c
c    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Jeremy Kasdin,
c    Runge-Kutta algorithm for the numerical integration of
c    stochastic differential equations,
c    Journal of Guidance, Control, and Dynamics,
c    Volume 18, Number 1, January-February 1995, pages 114-120.
c
c    Jeremy Kasdin,
c    Discrete Simulation of Colored Noise and Stochastic Processes
c    and 1/f^a Power Law Noise Generation,
c    Proceedings of the IEEE,
c    Volume 83, Number 5, 1995, pages 802-827.
c
c  Parameters:
c
c    Input, double precision X, the value at the current time.
c
c    Input, double precision T, the current time.
c
c    Input, double precision H, the time step.
c
c    Input, double precision Q, the spectral density of the input white noise.
c
c    Input, external double precision FI, the name of the deterministic
c    right hand side function.
c
c    Input, external double precision GI, the name of the stochastic
c    right hand side function.
c
c    Input/output, integer SEED, a seed for the random 
c    number generator.
c
c    Output, double precision XSTAR, the value at time T+H.
c
      implicit none

      double precision a21
      double precision a31
      double precision a32
      double precision a41
      double precision a42
      double precision a43
      double precision fi
      external fi
      double precision gi
      external gi
      double precision h
      double precision k1
      double precision k2
      double precision k3
      double precision q
      double precision q1
      double precision q2
      double precision q3
      double precision r8_normal_01
      integer seed
      double precision t
      double precision t1
      double precision t2
      double precision t3
      double precision w1
      double precision w2
      double precision w3
      double precision x
      double precision x1
      double precision x2
      double precision x3
      double precision xstar

      a21 =   1.52880952525675D+00
      a31 =   0.0D+00
      a32 =   0.51578733443615D+00
      a41 =   0.53289582961739D+00
      a42 =   0.25574324768195D+00
      a43 =   0.21136092270067D+00

      q1 = 1.87653936176981D+00
      q2 = 3.91017166264989D+00
      q3 = 4.73124353935667D+00

      t1 = t
      x1 = x
      w1 = r8_normal_01 ( seed ) * sqrt ( q1 * q / h )
      k1 = h * fi ( x1 ) + h * gi ( x1 ) * w1 

      t2 = t1 + a21 * h
      x2 = x1 + a21 * k1
      w2 = r8_normal_01 ( seed ) * sqrt ( q2 * q / h )
      k2 = h * fi ( x2 ) + h * gi ( x2 ) * w2

      t3 = t1 + a31 * h  + a32 * h
      x3 = x1 + a31 * k1 + a32 * k2
      w3 = r8_normal_01 ( seed ) * sqrt ( q3 * q / h )
      k3 = h * fi ( x3 ) + h * gi ( x3 ) * w3

      xstar = x1 + a41 * k1 + a42 * k2 + a43 * k3

      return
      end
      subroutine rk4_ti_step ( x, t, h, q, fi, gi, seed, xstar )

c*********************************************************************72
c
cc RK4_TI_STEP takes one step of a stochastic Runge Kutta scheme.
c
c  Discussion:
c
c    The Runge-Kutta scheme is fourth-order, and suitable for time-invariant
c    systems in which F and G do not depend explicitly on time.
c
c    d/dx X(t,xsi) = F ( X(t,xsi) ) + G ( X(t,xsi) ) * w(t,xsi)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Jeremy Kasdin,
c    Runge-Kutta algorithm for the numerical integration of
c    stochastic differential equations,
c    Journal of Guidance, Control, and Dynamics,
c    Volume 18, Number 1, January-February 1995, pages 114-120.
c
c    Jeremy Kasdin,
c    Discrete Simulation of Colored Noise and Stochastic Processes
c    and 1/f^a Power Law Noise Generation,
c    Proceedings of the IEEE,
c    Volume 83, Number 5, 1995, pages 802-827.
c
c  Parameters:
c
c    Input, double precision X, the value at the current time.
c
c    Input, double precision T, the current time.
c
c    Input, double precision H, the time step.
c
c    Input, double precision Q, the spectral density of the input white noise.
c
c    Input, external double precision FI, the name of the deterministic
c    right hand side function.
c
c    Input, external double precision GI, the name of the stochastic
c    right hand side function.
c
c    Input/output, integer SEED, a seed for the random 
c    number generator.
c
c    Output, double precision XSTAR, the value at time T+H.
c
      implicit none

      double precision a21
      double precision a31
      double precision a32
      double precision a41
      double precision a42
      double precision a43
      double precision a51
      double precision a52
      double precision a53
      double precision a54
      double precision fi
      external fi
      double precision gi
      external gi
      double precision h
      double precision k1
      double precision k2
      double precision k3
      double precision k4
      double precision q
      double precision q1
      double precision q2
      double precision q3
      double precision q4
      double precision r8_normal_01
      integer seed
      double precision t
      double precision t1
      double precision t2
      double precision t3
      double precision t4
      double precision w1
      double precision w2
      double precision w3
      double precision w4
      double precision x
      double precision x1
      double precision x2
      double precision x3
      double precision x4
      double precision xstar

      a21 =   2.71644396264860D+00
      a31 = - 6.95653259006152D+00
      a32 =   0.78313689457981D+00
      a41 =   0.0D+00
      a42 =   0.48257353309214D+00
      a43 =   0.26171080165848D+00
      a51 =   0.47012396888046D+00
      a52 =   0.36597075368373D+00
      a53 =   0.08906615686702D+00
      a54 =   0.07483912056879D+00

      q1 =   2.12709852335625D+00
      q2 =   2.73245878238737D+00
      q3 =  11.22760917474960D+00
      q4 =  13.36199560336697D+00

      t1 = t
      x1 = x
      w1 = r8_normal_01 ( seed ) * sqrt ( q1 * q / h )
      k1 = h * fi ( x1 ) + h * gi ( x1 ) * w1 

      t2 = t1 + a21 * h
      x2 = x1 + a21 * k1
      w2 = r8_normal_01 ( seed ) * sqrt ( q2 * q / h )
      k2 = h * fi ( x2 ) + h * gi ( x2 ) * w2

      t3 = t1 + a31 * h  + a32 * h
      x3 = x1 + a31 * k1 + a32 * k2
      w3 = r8_normal_01 ( seed ) * sqrt ( q3 * q / h )
      k3 = h * fi ( x3 ) + h * gi ( x3 ) * w3

      t4 = t1 + a41 * h  + a42 * h + a43 * h
      x4 = x1 + a41 * k1 + a42 * k2
      w4 = r8_normal_01 ( seed ) * sqrt ( q4 * q / h )
      k4 = h * fi ( x4 ) + h * gi ( x4 ) * w4

      xstar = x1 + a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4

      return
      end
      subroutine rk1_tv_step ( x, t, h, q, fv, gv, seed, xstar )

c*********************************************************************72
c
cc RK1_TV_STEP takes one step of a stochastic Runge Kutta scheme.
c
c  Discussion:
c
c    The Runge-Kutta scheme is first-order, and suitable for time-varying
c    systems.
c
c    d/dx X(t,xsi) = F ( X(t,xsi), t ) + G ( X(t,xsi), t ) * w(t,xsi)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Jeremy Kasdin,
c    Runge-Kutta algorithm for the numerical integration of
c    stochastic differential equations,
c    Journal of Guidance, Control, and Dynamics,
c    Volume 18, Number 1, January-February 1995, pages 114-120.
c
c    Jeremy Kasdin,
c    Discrete Simulation of Colored Noise and Stochastic Processes
c    and 1/f^a Power Law Noise Generation,
c    Proceedings of the IEEE,
c    Volume 83, Number 5, 1995, pages 802-827.
c
c  Parameters:
c
c    Input, double precision X, the value at the current time.
c
c    Input, double precision T, the current time.
c
c    Input, double precision H, the time step.
c
c    Input, double precision Q, the spectral density of the input white noise.
c
c    Input, external double precision FV, the name of the deterministic
c    right hand side function.
c
c    Input, external double precision GV, the name of the stochastic
c    right hand side function.
c
c    Input/output, integer SEED, a seed for the random 
c    number generator.
c
c    Output, double precision XSTAR, the value at time T+H.
c
      implicit none

      double precision a21
      double precision fv
      external fv
      double precision gv
      external gv
      double precision h
      double precision k1
      double precision q
      double precision q1
      double precision r8_normal_01
      integer seed
      double precision t
      double precision t1
      double precision w1
      double precision x
      double precision x1
      double precision xstar

      a21 =   1.0D+00

      q1 = 1.0D+00

      t1 = t
      x1 = x
      w1 = r8_normal_01 ( seed ) * sqrt ( q1 * q / h )
      k1 = h * fv ( t1, x1 ) + h * gv ( t1, x1 ) * w1 

      xstar = x1 + a21 * k1

      return
      end
      subroutine rk2_tv_step ( x, t, h, q, fv, gv, seed, xstar )

c*********************************************************************72
c
cc RK2_TV_STEP takes one step of a stochastic Runge Kutta scheme.
c
c  Discussion:
c
c    The Runge-Kutta scheme is second-order, and suitable for time-varying
c    systems.
c
c    d/dx X(t,xsi) = F ( X(t,xsi), t ) + G ( X(t,xsi), t ) * w(t,xsi)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Jeremy Kasdin,
c    Runge-Kutta algorithm for the numerical integration of
c    stochastic differential equations,
c    Journal of Guidance, Control, and Dynamics,
c    Volume 18, Number 1, January-February 1995, pages 114-120.
c
c    Jeremy Kasdin,
c    Discrete Simulation of Colored Noise and Stochastic Processes
c    and 1/f^a Power Law Noise Generation,
c    Proceedings of the IEEE,
c    Volume 83, Number 5, 1995, pages 802-827.
c
c  Parameters:
c
c    Input, double precision X, the value at the current time.
c
c    Input, double precision T, the current time.
c
c    Input, double precision H, the time step.
c
c    Input, double precision Q, the spectral density of the input white noise.
c
c    Input, external double precision FV, the name of the deterministic
c    right hand side function.
c
c    Input, external double precision GV, the name of the stochastic
c    right hand side function.
c
c    Input/output, integer SEED, a seed for the random 
c    number generator.
c
c    Output, double precision XSTAR, the value at time T+H.
c
      implicit none

      double precision a21
      double precision a31
      double precision a32
      double precision fv
      external fv
      double precision gv
      external gv
      double precision h
      double precision k1
      double precision k2
      double precision q
      double precision q1
      double precision q2
      double precision r8_normal_01
      integer seed
      double precision t
      double precision t1
      double precision t2
      double precision w1
      double precision w2
      double precision x
      double precision x1
      double precision x2
      double precision xstar

      a21 =   1.0D+00
      a31 =   0.5D+00
      a32 =   0.5D+00

      q1 = 2.0D+00
      q2 = 2.0D+00

      t1 = t
      x1 = x
      w1 = r8_normal_01 ( seed ) * sqrt ( q1 * q / h )
      k1 = h * fv ( t1, x1 ) + h * gv ( t1, x1 ) * w1 

      t2 = t1 + a21 * h
      x2 = x1 + a21 * k1
      w2 = r8_normal_01 ( seed ) * sqrt ( q2 * q / h )
      k2 = h * fv ( t2, x2 ) + h * gv ( t2, x2 ) * w2

      xstar = x1 + a31 * k1 + a32 * k2

      return
      end
      subroutine rk4_tv_step ( x, t, h, q, fv, gv, seed, xstar )

c*********************************************************************72
c
cc RK4_TV_STEP takes one step of a stochastic Runge Kutta scheme.
c
c  Discussion:
c
c    The Runge-Kutta scheme is fourth-order, and suitable for time-varying
c    systems.
c
c    d/dx X(t,xsi) = F ( X(t,xsi), t ) + G ( X(t,xsi), t ) * w(t,xsi)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Jeremy Kasdin,
c    Runge-Kutta algorithm for the numerical integration of
c    stochastic differential equations,
c    Journal of Guidance, Control, and Dynamics,
c    Volume 18, Number 1, January-February 1995, pages 114-120.
c
c    Jeremy Kasdin,
c    Discrete Simulation of Colored Noise and Stochastic Processes
c    and 1/f^a Power Law Noise Generation,
c    Proceedings of the IEEE,
c    Volume 83, Number 5, 1995, pages 802-827.
c
c  Parameters:
c
c    Input, double precision X, the value at the current time.
c
c    Input, double precision T, the current time.
c
c    Input, double precision H, the time step.
c
c    Input, double precision Q, the spectral density of the input white noise.
c
c    Input, external double precision FV, the name of the deterministic
c    right hand side function.
c
c    Input, external double precision GV, the name of the stochastic
c    right hand side function.
c
c    Input/output, integer SEED, a seed for the random 
c    number generator.
c
c    Output, double precision XSTAR, the value at time T+H.
c
      implicit none

      double precision a21
      double precision a31
      double precision a32
      double precision a41
      double precision a42
      double precision a43
      double precision a51
      double precision a52
      double precision a53
      double precision a54
      double precision fv
      external fv
      double precision gv
      external gv
      double precision h
      double precision k1
      double precision k2
      double precision k3
      double precision k4
      double precision q
      double precision q1
      double precision q2
      double precision q3
      double precision q4
      double precision r8_normal_01
      integer seed
      double precision t
      double precision t1
      double precision t2
      double precision t3
      double precision t4
      double precision w1
      double precision w2
      double precision w3
      double precision w4
      double precision x
      double precision x1
      double precision x2
      double precision x3
      double precision x4
      double precision xstar

      a21 =   0.66667754298442D+00
      a31 =   0.63493935027993D+00
      a32 =   0.00342761715422D+00
      a41 = - 2.32428921184321D+00
      a42 =   2.69723745129487D+00
      a43 =   0.29093673271592D+00
      a51 =   0.25001351164789D+00
      a52 =   0.67428574806272D+00
      a53 = - 0.00831795169360D+00
      a54 =   0.08401868181222D+00

      q1 = 3.99956364361748D+00
      q2 = 1.64524970733585D+00
      q3 = 1.59330355118722D+00
      q4 = 0.26330006501868D+00

      t1 = t
      x1 = x
      w1 = r8_normal_01 ( seed ) * sqrt ( q1 * q / h )
      k1 = h * fv ( t1, x1 ) + h * gv ( t1, x1 ) * w1 

      t2 = t1 + a21 * h
      x2 = x1 + a21 * k1
      w2 = r8_normal_01 ( seed ) * sqrt ( q2 * q / h )
      k2 = h * fv ( t2, x2 ) + h * gv ( t2, x2 ) * w2

      t3 = t1 + a31 * h  + a32 * h
      x3 = x1 + a31 * k1 + a32 * k2
      w3 = r8_normal_01 ( seed ) * sqrt ( q3 * q / h )
      k3 = h * fv ( t3, x3 ) + h * gv ( t3, x3 ) * w3

      t4 = t1 + a41 * h  + a42 * h  + a43 * h
      x4 = x1 + a41 * k1 + a42 * k2 + a43 * k3
      w4 = r8_normal_01 ( seed ) * sqrt ( q4 * q / h )
      k4 = h * fv ( t4, x4 ) + h * gv ( t4, x4 ) * w4

      xstar = x1 + a51 * k1 + a52 * k2 + a53 * k3 + a54 * k4

      return
      end
      function r8_normal_01 ( seed )

c*********************************************************************72
c
cc R8_NORMAL_01 returns a unit pseudonormal R8.
c
c  Discussion:
c
c    Because this routine uses the Box Muller method, it requires pairs
c    of uniform random values to generate a pair of normal random values.
c    This means that on every other call, the code can use the second
c    value that it calculated.
c
c    However, if the user has changed the SEED value between calls,
c    the routine automatically resets itself and discards the saved data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, double precision R8_NORMAL_01, a sample of the standard normal PDF.
c
      implicit none

      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r1
      double precision r2
      double precision r8_normal_01
      double precision r8_uniform_01
      integer seed
      integer seed1
      integer seed2
      integer seed3
      integer used
      double precision v1
      double precision v2

      save seed1
      save seed2
      save seed3
      save used
      save v2

      data seed2 / 0 /
      data used / 0 /
      data v2 / 0.0D+00 /
c
c  If USED is odd, but the input SEED does not match
c  the output SEED on the previous call, then the user has changed
c  the seed.  Wipe out internal memory.
c
      if ( mod ( used, 2 ) == 1 ) then

        if ( seed .ne. seed2 ) then
          used = 0
          seed1 = 0
          seed2 = 0
          seed3 = 0
          v2 = 0.0D+00
        end if

      end if
c
c  If USED is even, generate two uniforms, create two normals,
c  return the first normal and its corresponding seed.
c
      if ( mod ( used, 2 ) .eq. 0 ) then

        seed1 = seed

        r1 = r8_uniform_01 ( seed )

        if ( r1 .eq. 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8_NORMAL_01 - Fatal error!'
          write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
          stop
        end if

        seed2 = seed

        r2 = r8_uniform_01 ( seed )

        seed3 = seed

        v1 = sqrt ( -2.0D+00 * log ( r1 ) ) * cos ( 2.0D+00 * pi * r2 )
        v2 = sqrt ( -2.0D+00 * log ( r1 ) ) * sin ( 2.0D+00 * pi * r2 )

        r8_normal_01 = v1
        seed = seed2
c
c  If USED is odd (and the input SEED matched the output value from
c  the previous call), return the second normal and its corresponding seed.
c
      else

        r8_normal_01 = v2
        seed = seed3

      end if

      used = used + 1

      return
      end
      function r8_uniform_01 ( seed )

c*********************************************************************72
c
cc R8_UNIFORM_01 returns a unit pseudorandom R8.
c
c  Discussion:
c
c    This routine implements the recursion
c
c      seed = 16807 * seed mod ( 2**31 - 1 )
c      r8_uniform_01 = seed / ( 2**31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      R8_UNIFORM_01
c      SEED      SEED
c
c         12345   207482415  0.096616
c     207482415  1790989824  0.833995
c    1790989824  2035175616  0.947702
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 July 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Springer Verlag, pages 201-202, 1983.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley Interscience, page 95, 1998.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, pages 362-376, 1986.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, pages 136-143, 1969.
c
c  Parameters:
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R8_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      integer k
      double precision r8_uniform_01
      integer seed

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + 2147483647
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r8_uniform_01 = dble ( seed ) * 4.656612875D-10

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Discussion:
c
c    This FORTRAN77 version is made available for cases where the
c    FORTRAN90 version cannot be used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    None
c
      implicit none

      character * ( 8 ) ampm
      integer d
      character * ( 8 ) date
      integer h
      integer m
      integer mm
      character * ( 9 ) month(12)
      integer n
      integer s
      character * ( 10 ) time
      integer y

      save month

      data month /
     &  'January  ', 'February ', 'March    ', 'April    ', 
     &  'May      ', 'June     ', 'July     ', 'August   ', 
     &  'September', 'October  ', 'November ', 'December ' /

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) y, m, d
      read ( time, '(i2,i2,i2,1x,i3)' ) h, n, s, mm

      if ( h .lt. 12 ) then
        ampm = 'AM'
      else if ( h .eq. 12 ) then
        if ( n .eq. 0 .and. s .eq. 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h .lt. 12 ) then
          ampm = 'PM'
        else if ( h .eq. 12 ) then
          if ( n .eq. 0 .and. s .eq. 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if

      write ( *, 
     &  '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) 
     &  d, month(m), y, h, ':', n, ':', s, '.', mm, ampm

      return
      end
