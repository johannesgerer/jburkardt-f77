      program main

c*********************************************************************72
c
cc MAIN is the main program for PCE_BURGERS.
c
c  Discussion:
c
c    The time-dependent viscous Burgers equation to be solved is:
c
c      du/dt = - d ( u*(1/2-u)) /dx + nu d2u/dx2
c
c    with boundary conditions
c
c      u(-3.0) = 0.0, u(+3.0) = 1.0.
c
c    The viscosity nu is assumed to be an uncertain quantity with
c    normal distribution of known mean and variance.
c
c    A polynomial chaos expansion is to be used, with Hermite polynomial
c    basis functions h(i,x), 0 <= i <= n.
c
c    Because the first two Hermite polynomials are simply 1 and x, 
c    we have that 
c
c      nu = nu_mean * h(0,x) + nu_variance * h(1,x).
c
c    We replace the time derivative by an explicit Euler approximation,
c    so that the equation now describes the value of U(x,t+dt) in terms
c    of known data at time t.
c
c    Now assume that the solution U(x,t) can be approximated
c    by the truncated expansion:
c
c      U(x,t) = sum ( 0 <= i <= n ) c(i,t) * h(i,x)
c
c    In the equation, we replace U by its expansion, and then multiply
c    successively by each of the basis functions h(*,x) to get a set of
c    n+1 equations that can be used to determine the values of c(i,t+dt).
c
c    This process is repeated until the desired final time is reached.
c
c    At any time, the coefficients c(0,t) contain information definining
c    the expected value of u(x,t) at that time, while the higher order 
c    coefficients can be used to deterimine higher moments.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 September 2012
c
c  Author:
c
c    Original FORTRAN90 version by Gianluca Iaccarino.
c    FORTRAN77 version is by John Burkardt.
c
c  Local parameters:
c
c    Local, double precision DT, the timestep.
c
c    Local, double precision DX, the spacing between grid points.
c
c    Local, integer N, the number of intervals in the spatial domain.
c
c    Local, double precision NUMEAN, the mean of viscosity.
c
c    Local, double precision NUVARIANCE, the variance of viscosity.
c
c    Local, integer P, the order of the PC expansion.
c
c    Local, double precision T, the current time.
c
c    Local, double precision TF, the final integration time.
c
c    Local, double precision U1(N+1,P+1), the PCE representation at the 
c    current time.
c
c    Local, double precision U2(N+1,P+1), the PCE representation for the next time.
c
c    Local, double precision X(N+1,1), the grid points.
c
      implicit none

      integer n
      parameter ( n = 32 )
      integer p
      parameter ( p = 5 )

      double precision conv
      double precision dp
      double precision dt
      double precision dx
      double precision he_double_product_integral
      double precision he_triple_product_integral
      integer i
      integer it
      integer ix
      integer j
      integer k
      integer nt
      double precision numean
      double precision nuvariance
      character * ( 80 ) output_filename
      integer output_unit
      double precision t1
      double precision t2
      double precision term1
      double precision term2
      double precision tf
      double precision ti
      double precision tp
      double precision u1(n+1,p+1)
      double precision u2(n+1,p+1)
      double precision umean(n+1)
      double precision uvariance(n+1)
      double precision visc(2)
      double precision x(n+1)

      nt = 2000
      ti = 0.0D+00
      tf = 2.0D+00
      dt = ( tf - ti ) / dble ( nt )
      numean = 0.25D+00
      nuvariance = 0.08D+00

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PCE_BURGERS:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Polynomial Chaos Solution'
      write ( *, '(a)' ) '  1D Burgers equation'
      write ( *, '(a)' ) '  Original version by Gianluca Iaccarino'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  PCE order       = ', p
      write ( *, '(a,i8)' ) '  Number of cells = ', n
      write ( *, '(a,g14.6)' ) '  Time step       = ', dt
      write ( *, '(a,g14.6)' ) '  Initial time    = ', ti
      write ( *, '(a,g14.6)' ) '  Final time      = ', tf
      write ( *, '(a,g14.6)' ) '  Viscosity Mean  = ', numean
      write ( *, '(a,g14.6)' ) '  Viscosity Var   = ', nuvariance
      write ( *, '(a)' ) ' '
c
c  Define some numerical parameters
c
      dx = 6.0D+00 / dble ( n )
      conv = dt / ( 2.0D+00 * dx )
c
c  The expansion for viscosity stops at the linear term.
c
      visc(1) = numean * dt / ( dx * dx )
      visc(2) = nuvariance * dt / ( dx * dx )
c
c  Define a uniform grid
c
      do i = 1, n + 1
        x(i) = ( dble ( n - i + 1 ) * ( -3.0D+00 )   
     &         + dble (     i - 1 ) * ( +3.0D+00 ) ) 
     &         / dble ( n )
      end do
c
c  Set the initial conditions
c
      do j = 1, p + 1
        do i = 1, n + 1
          u1(i,j) = 0.0D+00
        end do
      end do

      do i = 1, n + 1
        u1(i,1) = 0.5D+00 + x(i) / 6.0D+00
      end do
c
c  Write the current solution.
c
      output_filename = 'burgers.history.txt'

      call get_unit ( output_unit )

      open ( unit = output_unit, file = output_filename, 
     &  status = 'replace' )

      write ( output_unit, '(a)' ) '----------'
      write ( output_unit, '(a,g14.6)' ) 'T = ', t1
      do i = 1, n + 1
        write ( output_unit, '(10(2x,g14.6))' ) 
     &    x(i), ( u1(i,j), j  = 1 , p + 1 )
      end do
c
c  Time integration
c
      t1 = ti

      do it = 1, nt

        t2 = ( dble ( nt - it ) * ti   
     &       + dble (      it ) * tf ) 
     &       / dble ( nt )
c
c  Boundary conditions.
c
        do j = 1, p + 1
          u2(1,j) = 0.0D+00
        end do

        u2(n+1,1) = 1.0D+00

        do j = 2, p + 1
          u2(n+1,j) = 0.0D+00
        end do

        do k = 1, p + 1

          dp = he_double_product_integral ( k - 1, k - 1 )

          do ix = 2, n
c
c  Viscous term.
c
            term1 = visc(1) 
     &        * ( u1(ix+1,k) - 2.0D+00 * u1(ix,k) + u1(ix-1,k) )
            i = 2
            do j = 1, p + 1
              tp = he_triple_product_integral ( i - 1, j - 1, k - 1 )
              term1 = term1 + visc(i) * tp *
     &          ( u1(ix+1,j) - 2.0D+00 * u1(ix,j) + u1(ix-1,j) ) / dp
            end do
c
c  Convective term.
c
            term2 = - conv * 0.5D+00 * ( u1(ix+1,k) - u1(ix-1,k) )
            do j = 1, p + 1
              do i = 1, p + 1
                tp = he_triple_product_integral ( i - 1, j - 1, k - 1 )
                term2 = term2 + ( conv * u1(ix,i) 
     &            * ( u1(ix+1,j) - u1(ix-1,j) ) * tp ) / dp
              end do
            end do

            u2(ix,k) = u1(ix,k) + term1 + term2

          end do
        end do

        t1 = t2
        do j = 1, p + 1
          do i = 1, n + 1
            u1(i,j) = u2(i,j)
          end do
        end do
c
c  Print solution every 100 time steps.
c
        if ( mod ( it, 100 ) == 0 ) then
          write ( output_unit, '(a)' ) '----------'
          write ( output_unit, '(a,g14.6)' ) 'T = ', t1
          do i = 1, n + 1
            write ( output_unit, '(20(2x,g14.6))' ) 
     &        x(i), ( u1(i,j), j = 1, p + 1 )
          end do
        end if

      end do

      close ( unit = output_unit )
      write ( *, '(a)' ) '  Time history in "' 
     &  // trim ( output_filename ) // '".'
c
c  Compute the mean and variance.
c
      do i = 1, n + 1
        umean(i) = u1(i,1)
      end do

      do i = 1, n + 1
        uvariance(i) = 0.0D+00
        do j = 2, p + 1
          dp = he_double_product_integral ( j - 1, j - 1 )
          uvariance(i) = uvariance(i) + u1(i,j)**2 * dp
        end do
      end do
c
c  Write data about the solution at the final time.
c
      output_filename = 'burgers.moments.txt'
      open ( unit = output_unit, file = output_filename, 
     &  status = 'replace' )
      write ( output_unit, '(a)' ) 'X E[U] Var[U]'
      do i = 1, n + 1
        write ( output_unit, '(20(g18.8,2x))' ) 
     &    x(i), umean(i), uvariance(i)
      end do
      close ( unit = output_unit )
      write ( *, '(a)' ) '  Moments written to "' // 
     &  trim ( output_filename ) // '".'

      output_filename = 'burgers.modes.txt'
      open ( unit = output_unit, file = output_filename, 
     &  status = 'replace' )
      write ( output_unit, '(a)' ) 'X U_0 ... U_P'
      do i = 1, n + 1
        write ( output_unit, '(20(2x,g14.6))' ) 
     &    x(i), ( u1(i,j), j = 1, p + 1 )
      end do
      close ( unit = output_unit )
      write ( *, '(a)' ) '  Final modes written to "' 
     &  // trim ( output_filename ) // '".'
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PCE_BURGERS:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine get_unit ( unit )

c*********************************************************************72
c
cc GET_UNIT returns a free FORTRAN unit number.
c
c  Discussion:
c
c    A "free" FORTRAN unit number is a value between 1 and 99 which
c    is not currently associated with an I/O device.  A free FORTRAN unit
c    number is needed in order to open a file with the OPEN command.
c
c    If UNIT = 0, then no free FORTRAN unit could be found, although
c    all 99 units were checked (except for units 5, 6 and 9, which
c    are commonly reserved for console I/O).
c
c    Otherwise, UNIT is a value between 1 and 99, representing a
c    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
c    are special, and will never return those values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 October 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer UNIT, the free unit number.
c
      implicit none

      integer i
      integer unit

      unit = 0

      do i = 1, 99

        if ( i .ne. 5 .and. i .ne. 6 .and. i .ne. 9 ) then

          open ( unit = i, err = 10, status = 'scratch' )
          close ( unit = i )

          unit = i

          return
        end if

10      continue

      end do

      return
      end
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
