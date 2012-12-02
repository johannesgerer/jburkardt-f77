      function he_double_product_integral ( i, j )

c*********************************************************************72
c
cc HE_DOUBLE_PRODUCT_INTEGRAL: integral of He(i,x)*He(j,x)*e^(-x^2/2).
c
c  Discussion:
c
c    VALUE = integral ( -oo < x < +oo ) He(i,x)*He(j,x) exp(-x^2/2) dx
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 March 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Dongbin Xiu,
c    Numerical Methods for Stochastic Computations: A Spectral Method Approach,
c    Princeton, 2010,
c    ISBN13: 978-0-691-14212-8,
c    LC: QA274.23.X58.
c
c  Parameters:
c
c    Input, integer I, J, the polynomial indices.
c
c    Output, double precision HE_DOUBLE_PRODUCT_INTEGRAL, the value of 
c    the integral.
c
      implicit none

      double precision he_double_product_integral
      integer i
      integer j
      double precision r8_factorial
      double precision value

      if ( i .eq. j ) then
        value = r8_factorial ( i )
      else
        value = 0.0D+00
      end if

      he_double_product_integral = value

      return
      end
      function he_triple_product_integral ( i, j, k )

c*********************************************************************72
c
cc HE_TRIPLE_PRODUCT_INTEGRAL: integral of He(i,x)*He(j,x)*He(k,x)*e^(-x^2/2).
c
c  Discussion:
c
c    VALUE = integral ( -oo < x < +oo ) He(i,x)*He(j,x)*He(k,x) exp(-x^2/2) dx
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 March 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Dongbin Xiu,
c    Numerical Methods for Stochastic Computations: A Spectral Method Approach,
c    Princeton, 2010,
c    ISBN13: 978-0-691-14212-8,
c    LC: QA274.23.X58.
c
c  Parameters:
c
c    Input, integer I, J, K, the polynomial indices.
c
c    Output, double precision HE_TRIPLE_PRODUCT_INTEGRAL, the value of the integral.
c
      implicit none

      double precision he_triple_product_integral
      integer i
      integer j
      integer k
      double precision r8_factorial
      integer s
      double precision value

      s = ( i + j + k ) / 2

      if ( s .lt. max ( i, j, k ) ) then
        value = 0.0D+00
      else if ( mod ( i + j + k, 2 ) .ne. 0 ) then
        value = 0.0D+00
      else
        value = r8_factorial ( i ) / r8_factorial ( s - i ) 
     &        * r8_factorial ( j ) / r8_factorial ( s - j ) 
     &        * r8_factorial ( k ) / r8_factorial ( s - k )
      end if

      he_triple_product_integral = value

      return
      end
      subroutine pce_ode_hermite ( ti, tf, nt, ui, np, alpha_mu, 
     &  alpha_sigma, t, u )

c*********************************************************************72
c
cc PCE_ODE_HERMITE applies the polynomial chaos expansion to a scalar ODE.
c
c  Discussion:
c
c    The deterministic equation is
c
c      du/dt = - alpha * u,
c      u(0) = u0
c
c    In the stochastic version, it is assumed that the decay coefficient
c    ALPHA is a Gaussian random variable with mean value ALPHA_MU and variance
c    ALPHA_SIGMA^2.
c
c    The exact expected value of the stochastic equation will be
c
c      u(t) = u0 * exp ( t^2/2)
c
c    This should be matched by the first component of the polynomial chaos
c    expansion.
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
c  Parameters:
c
c    Input, double precision TI, TF, the initial and final times.
c
c    Input, integer NT, the number of output points.
c
c    Input, double precision UI, the initial condition.
c
c    Input, integer NP, the degree of the expansion.  Polynomials 
c    of degree 0 through NP will be used.
c
c    Input, double precision ALPHA_MU, ALPHA_SIGMA, the mean and standard 
c    deviation of the decay coefficient.
c
c    Output, double precision T(0:NT), U(0:NT,0:NP), the times and the PCE 
c    coefficients at the successive time steps.
c
      implicit none

      integer np
      integer nt

      double precision alpha_mu
      double precision alpha_sigma
      double precision dp
      double precision dt
      double precision he_double_product_integral
      double precision he_triple_product_integral
      integer i
      integer it
      integer j
      integer k
      double precision t(0:nt)
      double precision t1
      double precision t2
      double precision term
      double precision tf
      double precision ti
      double precision tp
      double precision u(0:nt,0:np)
      double precision u1(0:np)
      double precision u2(0:np)
      double precision ui

      dt = ( tf - ti ) / dble ( nt )
c
c  Set the PCE coefficients for the initial time.
c
      t1 = ti
      u1(0) = ui
      do i = 1, np
        u1(i) = 0.0D+00
      end do
c
c  Copy into the output arrays.
c
      t(0) = t1
      do j = 0, np
        u(0,j) = u1(j)
      end do
c
c  Time integration.
c
      do it = 1, nt

        t2 = ( dble ( nt - it ) * ti   
     &       + dble (      it ) * tf ) 
     &       / dble ( nt )

        do k = 0, np

          dp = he_double_product_integral ( k, k )

          term = - alpha_mu * u1(k)

          i = 1
          do j = 0, np
            tp = he_triple_product_integral ( i, j, k )
            term = term - alpha_sigma * u1(j) * tp / dp
          end do

          u2(k) = u1(k) + dt * term

        end do
c
c  Prepare for next step.
c
        t1 = t2
        do i = 0, np
          u1(i) = u2(i)
        end do
c
c  Copy into the output arrays.
c
        t(it) = t1
        do j = 0, np
          u(it,j) = u1(j)
        end do

      end do

      return
      end
      function r8_factorial ( n )

c*********************************************************************72
c
cc R8_FACTORIAL computes the factorial of N.
c
c  Discussion:
c
c    factorial ( N ) = product ( 1 <= I <= N ) I
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the argument of the factorial function.
c    If N is less than 1, the function value is returned as 1.
c
c    Output, double precision R8_FACTORIAL, the factorial of N.
c
      implicit none

      integer i
      integer n
      double precision r8_factorial

      r8_factorial = 1.0D+00

      do i = 1, n
        r8_factorial = r8_factorial * dble ( i )
      end do

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
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
