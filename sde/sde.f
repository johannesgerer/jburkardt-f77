      subroutine bpath ( seed, n, w )

c*********************************************************************72
c
cc BPATH performs a Brownian path simulation.
c
c  Discussion:
c
c    This routine computes one simulation of discretized Brownian 
c    motion over the time interval [0,1] using N time steps.
c    The user specifies a random number seed.  Different values of
c    the seed will result in different realizations of the path.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 September 2012
c
c  Author:
c
c    Original MATLAB version by Desmond Higham.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Desmond Higham,
c    An Algorithmic Introduction to Numerical Simulation of
c    Stochastic Differential Equations,
c    SIAM Review,
c    Volume 43, Number 3, September 2001, pages 525-546.
c
c  Parameters:
c
c    Input/output, integer SEED, a seed for the random number 
c    generator.
c
c    Input, integer N, the number of steps.
c
c    Output, double precision W(0:N), the Brownian path.
c
      implicit none

      integer n

      double precision dt
      double precision dw(n)
      integer i
      integer j
      integer seed
      double precision tmax
      double precision w(0:n)

      tmax = 1.0D+00
      dt = tmax / dble ( n )
c
c  Define the increments dW.
c
      call r8vec_normal_01 ( n, seed, dw )

      do i = 1, n
        dw(i) = sqrt ( dt ) * dw(i)
      end do
c
c  W is the sum of the previous increments.
c
      w(0) = 0.0D+00
      do j = 1, n
        w(j) = w(j-1) + dw(j)
      end do

      return
      end
      subroutine bpath_gnuplot ( n, w )

c*********************************************************************72
c
cc BPATH_GNUPLOT writes a GNUPLOT input file to plot BPATH data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 September 2012
c
c  Author:
c
c    John Burkardt.
c
c  Reference:
c
c    Desmond Higham,
c    An Algorithmic Introduction to Numerical Simulation of
c    Stochastic Differential Equations,
c    SIAM Review,
c    Volume 43, Number 3, September 2001, pages 525-546.
c
c  Parameters:
c
c    Input, integer N, the number of steps.
c
c    Input, double precision W(0:N), the Brownian path.
c
      implicit none

      integer n

      character * ( 80 ) command_filename
      integer command_unit
      character * ( 80 ) data_filename
      integer data_unit
      integer i
      double precision t
      double precision w(0:n)
c
c  Create the data file.
c
      call get_unit ( data_unit )

      data_filename = 'bpath_data.txt'

      open ( unit = data_unit, file = data_filename, 
     &  status = 'replace' )

      do i = 0, n
        t = dble ( i ) / dble ( n )
        write ( data_unit, '(2x,g14.6,2x,g14.6)' ) t, w(i)
      end do
      close ( unit = data_unit )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  BPATH data stored in "' 
     &  // trim ( data_filename ) // '".'
c
c  Create the command file.
c
      call get_unit ( command_unit )

      command_filename = 'bpath_commands.txt'

      open ( unit = command_unit, file = command_filename, 
     &  status = 'replace' )

      write ( command_unit, '(a)' ) '# bpath_commands.txt'
      write ( command_unit, '(a)' ) '# created by sde::bpath_gnuplot.'
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) '# Usage:'
      write ( command_unit, '(a)' ) '#  gnuplot < bpath_commands.txt'
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) 'set term png'
      write ( command_unit, '(a)' ) 'set output "bpath.png"'
      write ( command_unit, '(a)' ) 'set xlabel "t"'
      write ( command_unit, '(a)' ) 'set ylabel "W(t)"'
      write ( command_unit, '(a)' ) 
     &  'set title "Brownian motion by BPATH"'
      write ( command_unit, '(a)' ) 'set grid'
      write ( command_unit, '(a)' ) 'set style data lines'
      write ( command_unit, '(a)' ) 'plot "bpath_data.txt" using 1:2'
      write ( command_unit, '(a)' ) 'quit'

      close ( unit = command_unit )

      write ( *, '(a)' ) '  BPATH plot commands stored in "' 
     &  // trim ( command_filename ) // '".'

      return
      end
      subroutine bpath_average ( seed, m, n, u, umean, error )

c*********************************************************************72
c
cc BPATH_AVERAGE: displays the average of 1000 Brownian paths.
c
c  Discussion:
c
c    This routine computes M simulations of discretized Brownian 
c    motion W(t) over the time interval [0,1] using N time steps.
c    The user specifies a random number seed.  Different values of
c    the seed will result in a different set of realizations of the path.
c
c    Actually, we are interested in a function u(W(t)):
c
c      u(W(t)) = exp ( t + W(t)/2 )
c
c    The routine plots 5 of the simulations, as well as the average
c    of all the simulations.  
c
c    The plot of the average should be quite smooth.  Its expected
c    value is exp ( 9 * t / 8 ), and we compute the 'error', that is,
c    the difference between the averaged value and this expected
c    value.  This 'error' should decrease as the number of simulation
c    is increased.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 September 2012
c
c  Author:
c
c    Original Matlab version by Desmond Higham.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Desmond Higham,
c    An Algorithmic Introduction to Numerical Simulation of
c    Stochastic Differential Equations,
c    SIAM Review,
c    Volume 43, Number 3, September 2001, pages 525-546.
c
c  Parameters:
c
c    Input/output, integer SEED, a seed for the random 
c    number generator.
c
c    Input, integer M, the number of simulations to compute 
c    and average.  A typical value is 1000.
c
c    Input, integer N, the number of steps.  A typical value
c    is 500.
c
c    Output, double precision U(M,0:N), the M paths.
c
c    Output, double precision UMEAN(0:N), the averaged path.
c
c    Output, double precision ERROR, the maximum difference between the
c    averaged path and the exact expected value.
c
      implicit none

      integer m
      integer n

      double precision dt
      double precision dw(n)
      double precision e
      double precision error
      integer i
      integer j
      double precision r8vec_sum
      integer seed
      double precision t(0:n)
      double precision tmax
      double precision u(1:m,0:n)
      double precision umean(0:n)
      double precision w(0:n)

      tmax = 1.0D+00
      dt = tmax / dble ( n )
      do j = 0, n
        t(j) = dble ( j ) * tmax / dble ( n )
      end do

      do i = 1, m
c
c  Define the increments dW.
c
        call r8vec_normal_01 ( n, seed, dw )

        do j = 1, n
          dw(j) = sqrt ( dt ) * dw(j)
        end do
c
c  W is the sum of the previous increments.
c
        w(0) = 0.0D+00
        do j = 1, n
          w(j) = w(j-1) + dw(j)
        end do

        do j = 0, n
          u(i,j) = exp ( t(j) + 0.5D+00 * w(j) )
        end do

      end do
c
c  Average the M estimates of the path.
c
      do j = 0, n
        umean(j) = r8vec_sum ( m, u(1:m,j) ) / dble ( m )
      end do

      error = 0.0D+00
      do j = 0, n
        e = umean(j) - exp ( 9.0D+00 * t(j) / 8.0D+00 )
        error = max ( error, abs ( e ) )
      end do

      return
      end
      subroutine bpath_average_gnuplot ( m, n, u, umean )

c*********************************************************************72
c
cc BPATH_AVERAGE_GNUPLOT writes a GNUPLOT input file to plot BPATH_AVERAGE data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 September 2012
c
c  Author:
c
c    John Burkardt.
c
c  Reference:
c
c    Desmond Higham,
c    An Algorithmic Introduction to Numerical Simulation of
c    Stochastic Differential Equations,
c    SIAM Review,
c    Volume 43, Number 3, September 2001, pages 525-546.
c
c  Parameters:
c
c    Input, integer M, the number of simulations.
c
c    Input, integer N, the number of steps. 
c
c    Input, double precision U(M,0:N), the M paths.
c
c    Input, double precision UMEAN(0:N), the averaged path.
c
      implicit none

      integer m
      integer n

      character * ( 80 ) command_filename
      integer command_unit
      character * ( 80 ) data_filename
      integer data_unit
      integer i
      double precision t
      double precision u(m,0:n)
      double precision umean(0:n)
c
c  Create the data file.
c
      call get_unit ( data_unit )

      data_filename = 'bpath_average_data.txt'

      open ( unit = data_unit, file = data_filename, 
     &  status = 'replace' )

      do i = 0, n
        t = dble ( i ) / dble ( n )
        write ( data_unit, '(7(2x,g14.6))' ) t, u(1:5,i), umean(i)
      end do
      close ( unit = data_unit )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  BPATH_AVERAGE data stored in "' 
     &  // trim ( data_filename ) // '".'
c
c  Create the command file.
c
      call get_unit ( command_unit )

      command_filename = 'bpath_average_commands.txt'

      open ( unit = command_unit, file = command_filename, 
     &  status = 'replace' )

      write ( command_unit, '(a)' ) '# bpath_average_commands.txt'
      write ( command_unit, '(a)' ) 
     &  '# created by sde::bpath_average_gnuplot.'
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) '# Usage:'
      write ( command_unit, '(a)' ) 
     &  '#  gnuplot < bpath_average_commands.txt'
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) 'set term png'
      write ( command_unit, '(a)' ) 'set output "bpath_average.png"'
      write ( command_unit, '(a)' ) 'set xlabel "t"'
      write ( command_unit, '(a)' ) 'set ylabel "W(t)"'
      write ( command_unit, '(a)' ) 'set title "Averaged Brownian paths"'
      write ( command_unit, '(a)' ) 'set grid'
      write ( command_unit, '(a)' ) 'set style data lines'
      write ( command_unit, '(a)' ) 
     &  'plot "bpath_average_data.txt" using 1:2 title "sample 1", \'
      write ( command_unit, '(a)' ) 
     &  '     "bpath_average_data.txt" using 1:3 title "sample 2", \'
      write ( command_unit, '(a)' ) 
     &  '     "bpath_average_data.txt" using 1:4 title "sample 3", \'
      write ( command_unit, '(a)' ) 
     &  '     "bpath_average_data.txt" using 1:5 title "sample 4", \'
      write ( command_unit, '(a)' ) 
     &  '     "bpath_average_data.txt" using 1:6 title "sample 5", \'
      write ( command_unit, '(a)' ) 
     &  '     "bpath_average_data.txt" using 1:7 title "average" lw 3'
      write ( command_unit, '(a)' ) 'quit'

      close ( unit = command_unit )

      write ( *, '(a)' ) '  BPATH_AVERAGE plot commands stored in "' 
     &  // trim ( command_filename ) // '".'

      return
      end
      function ch_is_digit ( c )

c*********************************************************************72
c
cc CH_IS_DIGIT returns TRUE if a character is a decimal digit.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character C, the character to be analyzed.
c
c    Output, logical CH_IS_DIGIT, TRUE if C is a digit, FALSE otherwise.
c
      implicit none

      character c
      logical ch_is_digit

      if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then
        ch_is_digit = .true.
      else
        ch_is_digit = .false.
      end if

      return
      end
      subroutine ch_to_digit ( c, digit )

c*********************************************************************72
c
cc CH_TO_DIGIT returns the integer value of a base 10 digit.
c
c  Example:
c
c     C   DIGIT
c    ---  -----
c    '0'    0
c    '1'    1
c    ...  ...
c    '9'    9
c    ' '    0
c    'X'   -1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character C, the decimal digit, '0' through '9' or blank
c    are legal.
c
c    Output, integer DIGIT, the corresponding integer value.  If C was
c    'illegal', then DIGIT is -1.
c
      implicit none

      character c
      integer digit

      if ( lge ( c, '0' ) .and. lle ( c, '9' ) ) then

        digit = ichar ( c ) - 48

      else if ( c .eq. ' ' ) then

        digit = 0

      else

        digit = -1

      end if

      return
      end
      subroutine chain ( seed, n, xem, vem, diff )

c*********************************************************************72
c
cc CHAIN tests the stochastic Chain Rule.
c
c  Discussion:
c
c    This function solves a stochastic differential equation for 
c
c      V = sqrt(X) 
c
c    where X satisfies the stochastic differential equation:
c 
c      dX = ( alpha - X ) * dt + beta * sqrt(X) dW,
c      X(0) = Xzero,
c
c    with 
c
c      alpha = 2,
c      beta = 1,
c      Xzero = 1.
c
c    From the stochastic Chain Rule, the SDE for V is therefore:
c
c      dV = ( ( 4 * alpha - beta^2 ) / ( 8 * V ) - 1/2 V ) dt + 1/2 beta dW
c      V(0) = sqrt ( Xzero ).
c
c    Xem is the Euler-Maruyama solution for X. 
c
c    Vem is the Euler-Maruyama solution of the SDE for V from
c    the stochastic Chain Rule.
c
c    Hence, we compare sqrt(Xem) and Vem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 September 2012
c
c  Author:
c
c    Original Matlab version by Desmond Higham.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Desmond Higham,
c    An Algorithmic Introduction to Numerical Simulation of
c    Stochastic Differential Equations,
c    SIAM Review,
c    Volume 43, Number 3, September 2001, pages 525-546.
c
c  Parameters:
c
c    Input, integer SEED, a seed for the random number generator.
c
c    Input, integer N, the number of time steps.
c
c    Output, double precision XEM(0:N), the computed value of X.
c
c    Output, double precision VEM(0:N), the computed value of V.
c
c    Output, double precision DIFF, the maximum value of |sqrt(XEM)-V|.
c
      implicit none

      integer n

      double precision alpha
      double precision beta
      double precision diff
      double precision dt
      double precision dt2
      double precision dw(n)
      integer i
      integer j
      integer seed
      double precision tmax
      double precision vem(0:n)
      double precision xem(0:n)
c
c  Set problem parameters.
c
      alpha = 2.0D+00
      beta = 1.0D+00
c
c  Stepping parameters.
c  dt2 is the size of the Euler-Maruyama steps.
c
      tmax = 1.0D+00
      dt = tmax / dble ( n )
      dt2 = dt
c
c  Define the increments dW.
c
      call r8vec_normal_01 ( n, seed, dw )

      do j = 1, n
        dw(j) = sqrt ( dt ) * dw(j)
      end do
c
c  Solve for X(t).
c
      xem(0) = 1.0D+00
      do j = 1, n
        xem(j) = xem(j-1) + ( alpha - xem(j-1) ) * dt2 
     &                    + beta * sqrt ( xem(j-1) ) * dw(j)
      end do
c
c  Solve for V(t).
c
      vem(0) = sqrt ( xem(0) )
      do j = 1, n
        vem(j) = vem(j-1) 
     &    + ( ( 4.0D+00 * alpha - beta ** 2 ) / ( 8.0D+00 * vem(j-1) ) 
     &    - 0.5D+00 * vem(j-1) ) * dt2 
     &    + 0.5D+00 * beta * dw(j)
      end do
c
c  Compare sqrt(X) and V.
c
      diff = 0.0D+00
      do i = 0, n
        diff = max ( diff, abs ( sqrt ( xem(i) ) - vem(i) ) )
      end do

      return
      end
      subroutine chain_gnuplot ( n, x, v )

c*********************************************************************72
c
cc CHAIN_GNUPLOT writes a GNUPLOT input file to plot CHAIN data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 September 2012
c
c  Author:
c
c    John Burkardt.
c
c  Reference:
c
c    Desmond Higham,
c    An Algorithmic Introduction to Numerical Simulation of
c    Stochastic Differential Equations,
c    SIAM Review,
c    Volume 43, Number 3, September 2001, pages 525-546.
c
c  Parameters:
c
c    Input, integer N, the number of steps.
c
c    Input, double precision X(0:N), the value of X.
c
c    Input, double precision V(0:N), the value of V.
c
      implicit none

      integer n

      character * ( 80 ) command_filename
      integer command_unit
      character * ( 80 ) data_filename
      integer data_unit
      integer i
      double precision t
      double precision v(0:n)
      double precision x(0:n)
c
c  Create the data file.
c
      call get_unit ( data_unit )

      data_filename = 'chain_data.txt'

      open ( unit = data_unit, file = data_filename, 
     &  status = 'replace' )

      do i = 0, n
        t = dble ( i ) / dble ( n )
        write ( data_unit, '(3(2x,g14.6))' ) t, sqrt ( x(i) ), v(i)
      end do
      close ( unit = data_unit )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  CHAIN data stored in "' 
     &  // trim ( data_filename ) // '".'
c
c  Create the command file.
c
      call get_unit ( command_unit )

      command_filename = 'chain_commands.txt'

      open ( unit = command_unit, file = command_filename, 
     &  status = 'replace' )

      write ( command_unit, '(a)' ) '# chain_commands.txt'
      write ( command_unit, '(a)' ) '# created by sde::chain_gnuplot.'
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) '# Usage:'
      write ( command_unit, '(a)' ) '#  gnuplot < chain_commands.txt'
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) 'set term png'
      write ( command_unit, '(a)' ) 'set output "chain.png"'
      write ( command_unit, '(a)' ) 'set xlabel "t"'
      write ( command_unit, '(a)' ) 'set ylabel "Sqrt(X(t)) vs V(X(t))"'
      write ( command_unit, '(a)' ) 
     &  'set title "V(X(t)) from X(t) and from Chain Rule"'
      write ( command_unit, '(a)' ) 'set grid'
      write ( command_unit, '(a)' ) 'set style data lines'
      write ( command_unit, '(a)' ) 
     &  'plot "chain_data.txt" using 1:2 title "Sqrt(X(t))", \'
      write ( command_unit, '(a)' ) 
     &  '     "chain_data.txt" using 1:3 title "V(X(t))"'
      write ( command_unit, '(a)' ) 'quit'

      close ( unit = command_unit )

      write ( *, '(a)' ) '  CHAIN plot commands stored in "' 
     &  // trim ( command_filename ) // '".'

      return
      end
      subroutine digit_inc ( c )

c*********************************************************************72
c
cc DIGIT_INC increments a decimal digit.
c
c  Example:
c
c    Input  Output
c    -----  ------
c    '0'    '1'
c    '1'    '2'
c    ...
c    '8'    '9'
c    '9'    '0'
c    'A'    'A'
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character C, a digit to be incremented.
c
      implicit none

      character c
      integer digit

      call ch_to_digit ( c, digit )

      if ( digit .eq. -1 ) then
        return
      end if

      digit = digit + 1

      if ( digit .eq. 10 ) then
        digit = 0
      end if

      call digit_to_ch ( digit, c )

      return
      end
      subroutine digit_to_ch ( digit, c )

c*********************************************************************72
c
cc DIGIT_TO_CH returns the character representation of a decimal digit.
c
c  Example:
c
c    DIGIT   C
c    -----  ---
c      0    '0'
c      1    '1'
c    ...    ...
c      9    '9'
c     17    '*'
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIGIT, the digit value between 0 and 9.
c
c    Output, character C, the corresponding character, or '*' if DIGIT
c    was illegal.
c
      implicit none

      character c
      integer digit

      if ( 0 .le. digit .and. digit .le. 9 ) then

        c = char ( digit + 48 )

      else

        c = '*'

      end if

      return
      end
      subroutine em ( seed, n, t, xtrue, t2, xem, error )

c*********************************************************************72
c
cc EM applies the Euler-Maruyama method to a linear SDE.
c
c  Discussion:
c
c    The SDE is 
c
c      dX = lambda * X dt + mu * X dW,   
c      X(0) = Xzero,
c
c    where 
c
c      lambda = 2,
c      mu = 1,
c      Xzero = 1.
c
c    The discretized Brownian path over [0,1] uses
c    a stepsize dt = 2^(-8).
c
c    The Euler-Maruyama method uses a larger timestep Dt = R*dt,
c    where R is an integer.  For an SDE of the form
c
c      dX = f(X(t)) dt + g(X(t)) dW(t)
c
c    it has the form
c
c      X(j) = X(j-1) + f(X(j-1)) * Dt + g(X(j-1)) * ( W(j*Dt) - W((j-1)*Dt) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 September 2012
c
c  Author:
c
c    Original Matlab version by Desmond Higham.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Desmond Higham,
c    An Algorithmic Introduction to Numerical Simulation of
c    Stochastic Differential Equations,
c    SIAM Review,
c    Volume 43, Number 3, September 2001, pages 525-546
c
c  Parameters:
c
c    Input/output, integer SEED, a seed for the random 
c    number generator.
c
c    Input, integer N, the number of time steps.  A typical
c    value is 2^8.  N should be a multiple of 4.
c
c    Output, double precision T(0:N), the time values for the exact solution.
c
c    Output, double precision XTRUE(0:N), the exact solution.
c
c    Output, double precision T2(0:N/4), the time values for the 
c    Euler-Maruyama solution.
c
c    Output, double precision XEM(0:N/4), the Euler-Maruyama solution.
c
c    Output, double precision ERROR, the value of | XEM(T) - XTRUE(T) |.
c
      implicit none

      integer n

      double precision dt
      double precision dt2
      double precision dw(n)
      double precision dw2
      double precision error
      integer i
      integer j
      integer l
      double precision lambda
      double precision mu
      integer r
      integer seed
      double precision t(0:n)
      double precision t2(0:n/4)
      double precision tmax
      double precision w(0:n)
      double precision xem(0:n/4)
      double precision xtrue(0:n)
      double precision xzero
c
c  Set problem parameters.
c
      lambda = 2.0D+00
      mu = 1.0D+00
      xzero = 1.0D+00
c
c  Set stepping parameters.
c
      tmax = 1.0D+00
      dt = tmax / dble ( n )
c
c  Define the increments dW.
c
      call r8vec_normal_01 ( n, seed, dw )

      do j = 1, n
        dw(j) = sqrt ( dt ) * dw(j)
      end do
c
c  Sum the Brownian increments.
c
      w(0) = 0.0D+00
      do j = 1, n
        w(j) = w(j-1) + dw(j)
      end do

      do j = 0, n
        t(j) = dble ( j ) * tmax / dble ( n )
      end do
c
c  Compute the discretized Brownian path.
c
      do j = 0, n
        xtrue(j) = xzero * exp ( ( lambda - 0.5D+00 * mu**2 ) 
     &    * ( t(j) + mu * w(j) ) )
      end do
c
c  Set:
c  R, the multiplier for the EM step, 
c  Dt, the EM stepsize,
c  L, the number of EM steps (we need N to be a multiple of R!)
c
      r = 4
      dt2 = dble ( r ) * dt
      l = n / r

      do j = 0, l
        t2(j) = dble ( j ) * tmax / dble ( l )
      end do
c
c  Preallocate Xem for efficiency.
c
      xem(0) = xzero
      do j = 1, l
        dw2 = 0.0D+00
        do i = r * ( j - 1 ) + 1, r * j
          dw2 = dw2 + dw(i)
        end do
        xem(j) = xem(j-1) + dt2 * lambda * xem(j-1) 
     &    + mu * xem(j-1) * dw2
      end do

      error = abs ( xem(l) - xtrue(n) )

      return
      end
      subroutine em_gnuplot ( n, t, xtrue, t2, xem )

c*********************************************************************72
c
cc EM_GNUPLOT writes a GNUPLOT input file to plot EM data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 September 2012
c
c  Author:
c
c    John Burkardt.
c
c  Reference:
c
c    Desmond Higham,
c    An Algorithmic Introduction to Numerical Simulation of
c    Stochastic Differential Equations,
c    SIAM Review,
c    Volume 43, Number 3, September 2001, pages 525-546.
c
c  Parameters:
c
c    Input, integer N, the number of steps.
c
c    Input, double precision T(0:N), the time values for the exact solution.
c
c    Input, double precision XTRUE(0:N), the exact solution.
c
c    Input, double precision T2(0:N/4), the time values for the 
c    Euler-Maruyama solution.
c
c    Input, double precision XEM(0:N/4), the Euler-Maruyama solution.
c
      implicit none

      integer n

      character * ( 80 ) command_filename
      integer command_unit
      character * ( 80 ) data_filename
      integer data_unit
      integer i
      double precision t(0:n)
      double precision t2(0:n/4)
      double precision xem(0:n/4)
      double precision xtrue(0:n)
c
c  Create data file #1.
c
      call get_unit ( data_unit )

      data_filename = 'em1_data.txt'

      open ( unit = data_unit, file = data_filename, 
     &  status = 'replace' )

      do i = 0, n
        write ( data_unit, '(3(2x,g14.6))' ) t(i), xtrue(i)
      end do
      close ( unit = data_unit )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  EM data #1 stored in "' 
     &  // trim ( data_filename ) // '".'
c
c  Create data file #2.
c
      call get_unit ( data_unit )

      data_filename = 'em2_data.txt'

      open ( unit = data_unit, file = data_filename, 
     &  status = 'replace' )

      do i = 0, n / 4
        write ( data_unit, '(3(2x,g14.6))' ) t2(i), xem(i)
      end do
      close ( unit = data_unit )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  EM data #2 stored in "' 
     &  // trim ( data_filename ) // '".'
c
c  Create the command file.
c
      call get_unit ( command_unit )

      command_filename = 'em_commands.txt'

      open ( unit = command_unit, file = command_filename, 
     &  status = 'replace' )

      write ( command_unit, '(a)' ) '# em_commands.txt'
      write ( command_unit, '(a)' ) '# created by sde::em_gnuplot.'
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) '# Usage:'
      write ( command_unit, '(a)' ) '#  gnuplot < em_commands.txt'
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) 'set term png'
      write ( command_unit, '(a)' ) 'set output "em.png"'
      write ( command_unit, '(a)' ) 'set xlabel "t"'
      write ( command_unit, '(a)' ) 'set ylabel "X(t)"'
      write ( command_unit, '(a)' ) 
     &  'set title "Exact X(t) and Euler-Maruyama Estimate"'
      write ( command_unit, '(a)' ) 'set grid'
      write ( command_unit, '(a)' ) 'set style data lines'
      write ( command_unit, '(a)' ) 
     &  'plot "em1_data.txt" using 1:2 title "Exact X(t))", \'
      write ( command_unit, '(a)' ) 
     &  '     "em2_data.txt" using 1:2 title "EM X(t)"'
      write ( command_unit, '(a)' ) 'quit'

      close ( unit = command_unit )

      write ( *, '(a)' ) '  EM plot commands stored in "' 
     &  // trim ( command_filename ) // '".'

      return
      end
      subroutine emstrong ( seed, m, n, p_max, dtvals, xerr )

c*********************************************************************72
c
cc EMSTRONG tests the strong convergence of the EM method.
c
c  Discussion:
c
c    The SDE is 
c
c      dX = lambda * X dt + mu * X dW,   
c      X(0) = Xzero,
c
c    where 
c
c      lambda = 2,
c      mu = 1,
c      Xzero = 1.
c
c    The discretized Brownian path over [0,1] has dt = 2^(-9).
c
c    The Euler-Maruyama method uses 5 different timesteps: 
c      16*dt, 8*dt, 4*dt, 2*dt, dt.
c
c    We are interested in examining strong convergence at T=1,
c    that is
c
c      E | X_L - X(T) |.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 September 2012
c
c  Author:
c
c    Original Matlab version by Desmond Higham.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Desmond Higham,
c    An Algorithmic Introduction to Numerical Simulation of
c    Stochastic Differential Equations,
c    SIAM Review,
c    Volume 43, Number 3, September 2001, pages 525-546.
c
c  Parameters:
c
c    Input/output, integer SEED, a seed for the random 
c    number generator.
c
c    Input, integer M, the number of simulations to perform.
c    A typical value is M = 1000.
c
c    Input, integer N, the number of time steps to take.
c    A typical value is N = 512.
c
c    Input, integer P_MAX, the number of time step sizes to use.
c    A typical value is 5.
c
c    Output, double precision DTVALS(P_MAX), the time steps used.
c
c    Output, double precision XERR(P_MAX), the averaged absolute error in the
c    solution estimate at the final time.
c
      implicit none

      integer m
      integer n
      integer p_max

      double precision a(p_max,2)
      double precision dt
      double precision dt2
      double precision dtvals(p_max)
      double precision dw(n)
      double precision e
      integer i
      integer j
      integer l
      double precision lambda
      double precision mu
      integer p
      double precision q
      integer r
      double precision resid
      double precision rhs(p_max)
      integer s
      integer seed
      double precision sol(2)
      double precision tmax
      double precision w(0:n)
      double precision winc
      double precision xerr(p_max)
      double precision xtemp
      double precision xtrue
      double precision xzero
c
c  Set problem parameters.
c
      lambda = 2.0D+00
      mu = 1.0D+00
      xzero = 1.0D+00
c
c  Set stepping parameters.
c
      tmax = 1.0D+00
      dt = tmax / dble ( n )

      do p = 1, p_max
        dtvals(p) = dt * 2.0D+00 ** ( p - 1 )
      end do
c
c  Sample over discrete Brownian paths.
c
      do p = 1, p_max
        xerr(p) = 0.0D+00
      end do

      do s = 1, m
c
c  Define the increments dW.
c
        call r8vec_normal_01 ( n, seed, dw )

        do j = 1, n
          dw(j) = sqrt ( dt ) * dw(j)
        end do
c
c  Sum the increments to get the Brownian path.
c
        w(0) = 0.0D+00
        do j = 1, n
          w(j) = w(j-1) + dw(j)
        end do
c
c  Determine the true solution.
c
        xtrue = xzero * exp ( ( lambda - 0.5 * mu ** 2 ) + mu * w(n) )
c
c  Use the Euler-Maruyama method with 5 different time steps dt2 = r * dt
c  to estimate the solution value at time TMAX.
c
        do p = 1, p_max                       
          dt2 = dtvals(p)
          r = 2 ** ( p - 1 )
          l = n / r
          xtemp = xzero
          do j = 1, l
            winc = 0.0D+00
            do i = r * ( j - 1 ) + 1, r * j
              winc = winc + dw(i)
            end do
            xtemp = xtemp + dt2 * lambda * xtemp + mu * xtemp * winc
          end do
          xerr(p) = xerr(p) + abs ( xtemp - xtrue )
        end do

      end do

      do p = 1, p_max
        xerr(p) = xerr(p) / dble ( m )
      end do
c
c  Least squares fit of error = c * dt^q.
c
      do p = 1, p_max
        a(p,1) = 1.0D+00
        a(p,2) = log ( dtvals(p) )
        rhs(p) = log ( xerr(p) )
      end do

      call qr_solve ( p_max, 2, a, rhs, sol )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EMSTRONG:'
      write ( *, '(a)' ) 
     &  '  Least squares solution to Error = c * dt ^ q'
      write ( *, '(a)' ) '  Expecting Q to be about 1/2.'
      write ( *, '(a,g14.6)' ) '  Computed Q = ', sol(2)

      resid = 0.0D+00
      do i = 1, p_max
        e = a(i,1) * sol(1) + a(i,2) * sol(2) - rhs(i)
        resid = resid + e * e
      end do
      resid = sqrt ( resid )
      write ( *, '(a,g14.6)' ) '  Residual is ', resid

      return
      end
      subroutine emstrong_gnuplot ( p_max, dtvals, xerr )

c*********************************************************************72
c
cc EMSTRONG_GNUPLOT writes a GNUPLOT input file to plot EMSTRONG data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 September 2012
c
c  Author:
c
c    John Burkardt.
c
c  Reference:
c
c    Desmond Higham,
c    An Algorithmic Introduction to Numerical Simulation of
c    Stochastic Differential Equations,
c    SIAM Review,
c    Volume 43, Number 3, September 2001, pages 525-546.
c
c  Parameters:
c
c    Input, integer P_MAX, the number of time step sizes to use.
c
c    Input, double precision DTVALS(P_MAX), the time steps used.
c
c    Input, double precision XERR(P_MAX), the averaged absolute error in the
c    solution estimate at the final time.
c
      implicit none

      integer p_max

      character * ( 80 ) command_filename
      integer command_unit
      character * ( 80 ) data_filename
      integer data_unit
      integer i
      double precision dtvals(p_max)
      double precision xerr(p_max)
c
c  Create data file.
c
      call get_unit ( data_unit )

      data_filename = 'emstrong_data.txt'

      open ( unit = data_unit, file = data_filename, 
     &  status = 'replace' )

      do i = 1, p_max
        write ( data_unit, '(3(2x,g14.6))' ) 
     &    dtvals(i), xerr(i), sqrt ( dtvals(i) )
      end do
      close ( unit = data_unit )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  EMSTRONG data stored in "' 
     &  // trim ( data_filename ) // '".'
c
c  Create the command file.
c
      call get_unit ( command_unit )

      command_filename = 'emstrong_commands.txt'

      open ( unit = command_unit, file = command_filename, 
     &  status = 'replace' )

      write ( command_unit, '(a)' ) '# emstrong_commands.txt'
      write ( command_unit, '(a)' ) 
     &  '# created by sde::emstrong_gnuplot.'
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) '# Usage:'
      write ( command_unit, '(a)' ) '#  gnuplot < emstrong_commands.txt'
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) 'set term png'
      write ( command_unit, '(a)' ) 'set output "emstrong.png"'
      write ( command_unit, '(a)' ) 'set xlabel "Log(dt)"'
      write ( command_unit, '(a)' ) 
     &  'set ylabel "Log(Averaged Error at final T)"'
      write ( command_unit, '(a)' ) 'set logscale xy 10'
      write ( command_unit, '(a)' ) 
     &  'set title "Euler-Maruyama Error as function of DT"'
      write ( command_unit, '(a)' ) 'set grid'
      write ( command_unit, '(a)' ) 'set style data linespoints'
      write ( command_unit, '(a)' ) 
     &  'plot "emstrong_data.txt" using 1:2 title "Error", \'
      write ( command_unit, '(a)' ) 
     &  '     "emstrong_data.txt" using 1:3 title "Slope = 1/2"'
      write ( command_unit, '(a)' ) 'quit'

      close ( unit = command_unit )

      write ( *, '(a)' ) '  EMSTRONG plot commands stored in "' 
     &  // trim ( command_filename ) // '".'

      return
      end
      subroutine emweak ( seed, method, m, p_max, dtvals, xerr )

c*********************************************************************72
c
cc EMWEAK tests the weak convergence of the Euler-Maruyama method.
c
c  Discussion:
c
c    The SDE is 
c
c      dX = lambda * X dt + mu * X dW,   
c      X(0) = Xzero,
c
c    where 
c
c      lambda = 2,
c      mu = 1,
c      Xzero = 1.
c
c    The discretized Brownian path over [0,1] has dt = 2^(-9).
c
c    The Euler-Maruyama method will use 5 different timesteps:
c
c      2^(p-10),  p = 1,2,3,4,5.
c
c    We examine weak convergence at T=1:
c
c      | E (X_L) - E (X(T)) |.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 September 2012
c
c  Author:
c
c    Original MATLAB version by Desmond Higham.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Desmond Higham,
c    An Algorithmic Introduction to Numerical Simulation of
c    Stochastic Differential Equations,
c    SIAM Review,
c    Volume 43, Number 3, September 2001, pages 525-546
c
c  Parameters:
c
c    Input, integer SEED, a seed for the random number generator.
c
c    Input, integer METHOD.
c    0, use the standard Euler-Maruyama method;
c    1, use the weak Euler-Maruyama method.
c
c    Input, integer M, the number of simulations to perform.
c    A typical value is M = 1000.
c
c    Input, integer P_MAX, the number of time step sizes to use.
c    A typical value is 5.
c
c    Output, double precision DTVALS(P_MAX), the time steps used.
c
c    Output, double precision XERR(P_MAX), the averaged absolute error in the
c    solution estimate at the final time.
c
      implicit none

      integer m
      integer p_max

      double precision a(p_max,2)
      double precision dt
      double precision dt2
      double precision dtvals(p_max)
      double precision e
      integer i
      integer j
      integer l
      double precision lambda
      integer method
      double precision mu
      integer p
      double precision q
      integer r
      double precision r8_sign
      double precision resid
      double precision rhs(p_max)
      integer s
      integer seed
      double precision sol(2)
      double precision tmax
      double precision winc(m)
      double precision xem(p_max)
      double precision xerr(p_max)
      double precision xtemp(m)
      double precision xtrue
      double precision xzero
c
c  Problem parameters;
c
      lambda = 2.0D+00
      mu = 0.1D+00
      xzero = 1.0D+00
c
c  Stepping parameters.
c
      tmax = 1.0D+00
      do p = 1, p_max
        dtvals(p) = 2.0D+00 ** ( p - 10 )
      end do
c
c  Take various Euler timesteps.
c  For stepsize dt, we will need to take L Euler steps to reach time TMAX.
c
      do p = 1, p_max

        l = 2 ** ( 10 - p )
        dt = dtvals(p)

        do i = 1, m
          xtemp(i) = xzero
        end do

        do j = 1, l
        
          if ( method == 0 ) then
            call r8vec_normal_01 ( m, seed, winc )
            do i = 1, m
              winc(i) = sqrt ( dt ) * winc (i)
            end do
          else
            call r8vec_normal_01 ( m, seed, winc )
            do i = 1, m
              winc(i) = sqrt ( dt ) * r8_sign ( winc(i) )
            end do
          end if

          do i = 1, m
            xtemp(i) = xtemp(i) + dt * lambda * xtemp(i) 
     &        + mu * xtemp(i) * winc(i)
          end do

        end do
c
c  Average the M results for this stepsize.
c
        call r8vec_mean ( m, xtemp, xem(p) )

      end do
c
c  Compute the error in the estimates for each stepsize.
c
      do p = 1, p_max
        xerr(p) = abs ( xem(p) - exp ( lambda ) )
      end do
c
c  Least squares fit of error = c * dt^q.
c
      do p = 1, p_max
        a(p,1) = 1.0D+00
        a(p,2) = log ( dtvals(p) )
        rhs(p) = log ( xerr(p) )
      end do

      call qr_solve ( p_max, 2, a, rhs, sol )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EMWEAK:'
      if ( method == 0 ) then
        write ( *, '(a)' ) '  Using standard Euler-Maruyama method.'
      else
        write ( *, '(a)' ) '  Using weak Euler-Maruyama method.'
      end if
      write ( *, '(a)' ) 
     &  '  Least squares solution to Error = c * dt ^ q'
      write ( *, '(a)' ) '  Expecting Q to be about 1.'
      write ( *, '(a,g14.6)' ) '  Computed Q = ', sol(2)

      resid = 0.0D+00
      do i = 1, p_max
        e = a(i,1) * sol(1) + a(i,2) * sol(2) - rhs(i)
        resid = resid + e * e
      end do
      resid = sqrt ( resid )
      write ( *, '(a,g14.6)' ) '  Residual is ', resid

      return
      end
      subroutine emweak_gnuplot ( p_max, dtvals, xerr, method )

c*********************************************************************72
c
cc EMWEAK_GNUPLOT writes a GNUPLOT input file to plot EMWEAK data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 September 2012
c
c  Author:
c
c    John Burkardt.
c
c  Reference:
c
c    Desmond Higham,
c    An Algorithmic Introduction to Numerical Simulation of
c    Stochastic Differential Equations,
c    SIAM Review,
c    Volume 43, Number 3, September 2001, pages 525-546.
c
c  Parameters:
c
c    Input, integer P_MAX, the number of time step sizes to use.
c
c    Input, double precision DTVALS(P_MAX), the time steps used.
c
c    Input, double precision XERR(P_MAX), the averaged absolute error in the
c    solution estimate at the final time.
c
c    Input, integer METHOD.
c    0, use the standard Euler-Maruyama method;
c    1, use the weak Euler-Maruyama method.
c
      implicit none

      integer p_max

      character * ( 80 ) command_filename
      integer command_unit
      character * ( 80 ) data_filename
      integer data_unit
      integer i
      double precision dtvals(p_max)
      integer method
      double precision xerr(p_max)
c
c  Create data file.
c
      call get_unit ( data_unit )

      if ( method == 0 ) then
        data_filename = 'emweak0_data.txt'
      else 
        data_filename = 'emweak1_data.txt'
      end if

      open ( unit = data_unit, file = data_filename, 
     &  status = 'replace' )

      do i = 1, p_max
        write ( data_unit, '(3(2x,g14.6))' ) dtvals(i), xerr(i)
      end do
      close ( unit = data_unit )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  EMWEAK data stored in "' 
     &  // trim ( data_filename ) // '".'
c
c  Create the command file.
c
      call get_unit ( command_unit )

      if ( method == 0 ) then
        command_filename = 'emweak0_commands.txt'
      else
        command_filename = 'emweak1_commands.txt'
      end if

      open ( unit = command_unit, file = command_filename, 
     &  status = 'replace' )

      if ( method == 0 ) then
        write ( command_unit, '(a)' ) '# emweak0_commands.txt'
      else
        write ( command_unit, '(a)' ) '# emweak1_commands.txt'
      end if

      write ( command_unit, '(a)' ) '# created by sde::emweak_gnuplot.'
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) '# Usage:'
      if ( method == 0 ) then
        write ( command_unit, '(a)' ) 
     &    '#  gnuplot < emweak0_commands.txt'
      else
        write ( command_unit, '(a)' ) 
     &    '#  gnuplot < emweak1_commands.txt'
      end if
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) 'set term png'
      if ( method == 0 ) then
        write ( command_unit, '(a)' ) 'set output "emweak0.png"'
      else
        write ( command_unit, '(a)' ) 'set output "emweak1.png"'
      end if
      write ( command_unit, '(a)' ) 'set xlabel "Log(dt)"'
      write ( command_unit, '(a)' ) 
     &    'set ylabel "Log(Averaged Error at final T)"'
      write ( command_unit, '(a)' ) 'set logscale xy 10'
      if ( method == 0 ) then
        write ( command_unit, '(a)' ) 
     &    'set title "Standard Euler-Maruyama Error as function of DT"'
      else
        write ( command_unit, '(a)' ) 
     &    'set title "Weak Euler-Maruyama Error as function of DT"'
      end if
      write ( command_unit, '(a)' ) 'set grid'
      write ( command_unit, '(a)' ) 'set style data linespoints'
      if ( method == 0 ) then
        write ( command_unit, '(a)' ) 
     &    'plot "emweak0_data.txt" using 1:2 title "Error", \'
        write ( command_unit, '(a)' ) 
     &    '     "emweak0_data.txt" using 1:1 title "Slope = 1"'
      else
        write ( command_unit, '(a)' ) 
     &    'plot "emweak1_data.txt" using 1:2 title "Error", \'
        write ( command_unit, '(a)' ) 
     &    '     "emweak1_data.txt" using 1:1 title "Slope = 1"'
      end if

      write ( command_unit, '(a)' ) 'quit'

      close ( unit = command_unit )

      write ( *, '(a)' ) '  EMWEAK plot commands stored in "' 
     &  // trim ( command_filename ) // '".'

      return
      end
      subroutine filename_inc ( filename )

c*********************************************************************72
c
cc FILENAME_INC increments a partially numeric filename.
c
c  Discussion:
c
c    It is assumed that the digits in the name, whether scattered or
c    connected, represent a number that is to be increased by 1 on
c    each call.  Non-numeric letters of the name are unaffected.
c
c    If the name is empty, then the routine stops.
c
c    If the name contains no digits, the empty string is returned.
c
c  Example:
c
c      Input          Output
c      -----          ------
c      a7to11.txt     a7to12.txt
c      a7to99.txt     a8to00.txt
c      a9to99.txt     a0to00.txt
c      cat.txt        ' '
c      ' '            STOP!
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character*(*) FILENAME.
c    On input, a character string to be incremented.
c    On output, the incremented string.
c
      implicit none

      character c
      logical ch_is_digit
      integer change
      integer digit
      character*(*) filename
      integer i
      integer lens

      lens = len_trim ( filename )

      if ( lens .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FILENAME_INC - Fatal error!'
        write ( *, '(a)' ) '  The input string is empty.'
        stop
      end if

      change = 0

      do i = lens, 1, -1

        c = filename(i:i)

        if ( ch_is_digit ( c ) ) then

          change = change + 1

          call digit_inc ( c )

          filename(i:i) = c

          if ( c .ne. '0' ) then
            return
          end if

        end if

      end do
c
c  No digits were found.  Return blank.
c
      if ( change .eq. 0 ) then
        filename = ' '
        return
      end if

      return
      end
      subroutine get_unit ( unit )

c*********************************************************************72
c
cc GET_UNIT returns a free FORTRAN unit number.
c
c  Discussion:
c
c    A "free" FORTRAN unit number is an integer between 1 and 99 which
c    is not currently associated with an I/O device.  A free FORTRAN unit
c    number is needed in order to open a file with the OPEN command.
c
c    If UNIT = 0, then no free FORTRAN unit could be found, although
c    all 99 units were checked (except for units 5, 6 and 9, which
c    are commonly reserved for console I/O).
c
c    Otherwise, UNIT is an integer between 1 and 99, representing a
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
      subroutine milstrong ( seed, p_max, dtvals, xerr )

c*********************************************************************72
c
cc MILSTRONG tests the strong convergence of the Milstein method.
c
c  Discussion:
c
c    This function solves the stochastic differential equation
c
c      dX = sigma * X * ( k - X ) dt + beta * X dW,  
c      X(0) = Xzero,
c
c    where 
c
c       sigma = 2, 
c       k = 1, 
c       beta = 1,
c       Xzero = 0.5.
c
c    The discretized Brownian path over [0,1] has dt = 2^(-11).
c
c    The Milstein method uses timesteps 128*dt, 64*dt, 32*dt, 16*dt 
c    (also dt for reference).
c
c    We examine strong convergence at T=1:  
c
c      E | X_L - X(T) |.
c
c    The code is vectorized: all paths computed simultaneously.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 September 2012
c
c  Author:
c
c    Original MATLAB version by Desmond Higham.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Desmond Higham,
c    An Algorithmic Introduction to Numerical Simulation of
c    Stochastic Differential Equations,
c    SIAM Review,
c    Volume 43, Number 3, September 2001, pages 525-546.
c
c  Parameters:
c
c    Input/output, integer SEED, a seed for the random 
c    number generator.
c
c    Input, integer P_MAX, the number of time step sizes to use.
c    A typical value is 4.
c
c    Output, double precision DTVALS(P_MAX), the time steps used.
c
c    Output, double precision XERR(P_MAX), the averaged absolute error in the
c    solution estimate at the final time.
c
      implicit none

      integer m
      parameter ( m = 500 )
      integer n
      parameter ( n = 2 ** 11 )
      integer p_max

      double precision a(p_max,2)
      double precision beta
      double precision dt
      double precision dtp
      double precision dtvals(1:p_max)
      double precision dw(m,n)
      double precision e
      integer i
      integer i2
      integer j
      double precision k
      integer l
      integer p
      integer r
      double precision resid
      double precision rhs(p_max)
      integer seed
      double precision sigma
      double precision sol(2)
      double precision tmax
      double precision winc(m)
      double precision xerr(p_max)
      double precision xref(m)
      double precision xtemp(m)
      double precision xzero
c
c  Set problem parameters.
c
      sigma = 2.0D+00
      k = 1.0D+00
      beta = 0.25D+00
      xzero = 0.5D+00
c
c  Set stepping parameters.
c
      tmax = 1.0D+00
      dt = tmax / dble ( n )
c
c  Define the increments dW.
c
      call r8mat_normal_01 ( m, n, seed, dw )

      do j = 1, n
        do i = 1, m
          dw(i,j) = sqrt ( dt ) * dw(i,j)
        end do
      end do
c
c  Estimate the reference solution at time T M times.
c
      do i = 1, m
        xref(i) = xzero
      end do

      do j = 1, n
        do i = 1, m
          xref(i) = xref(i) 
     &      + dt * sigma * xref(i) * ( k - xref(i) ) 
     &      + beta * xref(i) * dw(i,j) 
     &      + 0.5D+00 * beta ** 2 * xref(i) * ( dw(i,j) ** 2 - dt )
        end do
      end do
c
c  Now compute M Milstein approximations at each of 4 timesteps,
c  and record the average errors.
c
      do p = 1, p_max
        dtvals(p) = dt * 8.0D+00 * 2.0D+00 ** p
      end do

      do p = 1, p_max
        xerr(p) = 0.0D+00
      end do

      do p = 1, p_max

        r = 8 * 2 ** p
        dtp = dtvals(p)
        l = n / r
        do i = 1, m
          xtemp(i) = xzero
        end do

        do j = 1, l
          do i = 1, m
            winc(i) = 0.0D+00
            do i2 = r * ( j - 1 ) + 1, r * j
              winc(i) = winc(i) + dw(i,i2)
            end do
            xtemp(i) = xtemp(i) 
     &        + dtp * sigma * xtemp(i) * ( k - xtemp(i) ) 
     &        + beta * xtemp(i) * winc(i) 
     &        + 0.5D+00 * beta ** 2 * xtemp(i) * ( winc(i) ** 2 - dtp )
          end do
        end do

        xerr(p) = 0.0D+00
        do i = 1, m
          xerr(p) = xerr(p) + abs ( xtemp(i) - xref(i) )
        end do
        xerr(p) = xerr(p) / dble ( m )

      end do
c
c  Least squares fit of error = C * dt^q
c
      do p = 1, p_max
        a(p,1) = 1.0D+00
        a(p,2) = log ( dtvals(p) )
        rhs(p) = log ( xerr(p) )
      end do

      call qr_solve ( p_max, 2, a, rhs, sol )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MILSTEIN:'
      write ( *, '(a)' ) 
     &  '  Least squares solution to Error = c * dt ^ q'
      write ( *, '(a)' ) '  Expecting Q to be about 1.'
      write ( *, '(a,g14.6)' ) '  Computed Q = ', sol(2)

      resid = 0.0D+00
      do i = 1, p_max
        e = a(i,1) * sol(1) + a(i,2) * sol(2) - rhs(i)
        resid = resid + e * e
      end do
      resid = sqrt ( resid )
      write ( *, '(a,g14.6)' ) '  Residual is ', resid

      return
      end
      subroutine milstrong_gnuplot ( p_max, dtvals, xerr )

c*********************************************************************72
c
cc MILSTRONG_GNUPLOT writes a GNUPLOT input file to plot MILSTRONG data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 September 2012
c
c  Author:
c
c    John Burkardt.
c
c  Reference:
c
c    Desmond Higham,
c    An Algorithmic Introduction to Numerical Simulation of
c    Stochastic Differential Equations,
c    SIAM Review,
c    Volume 43, Number 3, September 2001, pages 525-546.
c
c  Parameters:
c
c    Input, integer P_MAX, the number of time step sizes to use.
c
c    Input, double precision DTVALS(P_MAX), the time steps used.
c
c    Input, double precision XERR(P_MAX), the averaged absolute error in the
c    solution estimate at the final time.
c
      implicit none

      integer p_max

      character * ( 80 ) command_filename
      integer command_unit
      character * ( 80 ) data_filename
      integer data_unit
      integer i
      double precision dtvals(p_max)
      double precision xerr(p_max)
c
c  Create data file.
c
      call get_unit ( data_unit )

      data_filename = 'milstrong_data.txt'

      open ( unit = data_unit, file = data_filename, 
     &  status = 'replace' )

      do i = 1, p_max
        write ( data_unit, '(2(2x,g14.6))' ) dtvals(i), xerr(i)
      end do
      close ( unit = data_unit )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  MILSTRONG data stored in "' 
     &  // trim ( data_filename ) // '".'
c
c  Create the command file.
c
      call get_unit ( command_unit )

      command_filename = 'milstrong_commands.txt'

      open ( unit = command_unit, file = command_filename, 
     &  status = 'replace' )

      write ( command_unit, '(a)' ) '# milstrong_commands.txt'
      write ( command_unit, '(a)' ) 
     &  '# created by sde::milstrong_gnuplot.'
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) '# Usage:'
      write ( command_unit, '(a)' ) 
     &  '#  gnuplot < milstrong_commands.txt'
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) 'set term png'
      write ( command_unit, '(a)' ) 'set output "milstrong.png"'
      write ( command_unit, '(a)' ) 'set xlabel "Log(dt)"'
      write ( command_unit, '(a)' ) 
     &  'set ylabel "Log(Averaged Error at final T)"'
      write ( command_unit, '(a)' ) 'set logscale xy 10'
      write ( command_unit, '(a)' ) 
     &  'set title "Milstein Error as function of DT"'
      write ( command_unit, '(a)' ) 'set grid'
      write ( command_unit, '(a)' ) 'set style data linespoints'
      write ( command_unit, '(a)' ) 
     &  'plot "milstrong_data.txt" using 1:2 title "Error", \'
      write ( command_unit, '(a)' ) 
     &  '     "milstrong_data.txt" using 1:1 title "Slope = 1"'
      write ( command_unit, '(a)' ) 'quit'

      close ( unit = command_unit )

      write ( *, '(a)' ) '  MILSTRONG plot commands stored in "' 
     &  // trim ( command_filename ) // '".'

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
      if ( mod ( used, 2 ) .eq. 1 ) then

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
      function r8_sign ( x )

c*********************************************************************72
c
cc R8_SIGN returns the sign of an R8.
c
c  Discussion:
c
c    value = -1 if X < 0;
c    value = +1 if X => 0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 August 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the number whose sign is desired.
c
c    Output, double precision R8_SIGN, the sign of X.
c
      implicit none

      double precision r8_sign
      double precision x

      if ( x .lt. 0.0D+00 ) then
        r8_sign = -1.0D+00
      else
        r8_sign = +1.0D+00
      end if

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
c      seed = 16807 * seed mod ( 2^31 - 1 )
c      r8_uniform_01 = seed / ( 2^31 - 1 )
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
c    11 August 2004
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

      double precision r8_uniform_01
      integer k
      integer seed

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + 2147483647
      end if

      r8_uniform_01 = dble ( seed ) * 4.656612875D-10

      return
      end
      subroutine r8mat_normal_01 ( m, n, seed, r )

c*********************************************************************72
c
cc R8MAT_NORMAL_01 returns a unit pseudonormal R8MAT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 November 2010
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
c    Input, integer M, N, the number of rows and columns in the array.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R(M,N), the array of pseudonormal values.
c
      implicit none

      integer m
      integer n

      integer seed
      double precision r(m,n)

      call r8vec_normal_01 ( m * n, seed, r )

      return
      end
      subroutine r8vec_mean ( n, a, mean )

c*********************************************************************72
c
cc R8VEC_MEAN returns the mean of an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, double precision A(N), the vector whose mean is desired.
c
c    Output, double precision MEAN, the mean of the vector entries.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      double precision mean

      mean = 0.0D+00
      do i = 1, n
        mean = mean + a(i)
      end do
      mean = mean / dble ( n )

      return
      end
      subroutine r8vec_normal_01 ( n, seed, x )

c*********************************************************************72
c
cc R8VEC_NORMAL_01 returns a unit pseudonormal R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    The standard normal probability distribution function (PDF) has
c    mean 0 and standard deviation 1.
c
c    This routine can generate a vector of values on one call.  It
c    has the feature that it should provide the same results
c    in the same order no matter how we break up the task.
c
c    The Box-Muller method is used, which is efficient, but
c    generates an even number of values each time.  On any call
c    to this routine, an even number of new values are generated.
c    Depending on the situation, one value may be left over.
c    In that case, it is saved for the next call.
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
c  Parameters:
c
c    Input, integer N, the number of values desired.  If N is negative,
c    then the code will flush its internal memory; in particular,
c    if there is a saved value to be used on the next call, it is
c    instead discarded.  This is useful if the user has reset the
c    random number seed, for instance.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, double precision X(N), a sample of the standard normal PDF.
c
c  Local parameters:
c
c    Local, integer MADE, records the number of values that have
c    been computed.  On input with negative N, this value overwrites
c    the return value of N, so the user can get an accounting of
c    how much work has been done.
c
c    Local, integer SAVED, is 0 or 1 depending on whether there is a
c    single saved value left over from the previous call.
c
c    Local, integer X_LO_INDEX, X_HI_INDEX, records the range of entries of
c    X that we need to compute.  This starts off as 1:N, but is adjusted
c    if we have a saved value that can be immediately stored in X(1),
c    and so on.
c
c    Local, double precision Y, the value saved from the previous call, if
c    SAVED is 1.
c
      implicit none

      integer n

      integer i
      integer m
      integer made
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r(2)
      double precision r8_uniform_01
      integer saved
      integer seed
      double precision x(n)
      integer x_hi_index
      integer x_lo_index
      double precision y

      save made
      save saved
      save y

      data made / 0 /
      data saved / 0 /
      data y / 0.0D+00 /
c
c  I'd like to allow the user to reset the internal data.
c  But this won't work properly if we have a saved value Y.
c  I'm making a crock option that allows the user to signal
c  explicitly that any internal memory should be flushed,
c  by passing in a negative value for N.
c
      if ( n .lt. 0 ) then
        n = made
        made = 0
        saved = 0
        y = 0.0D+00
        return
      else if ( n .eq. 0 ) then
        return
      end if
c
c  Record the range of X we need to fill in.
c
      x_lo_index = 1
      x_hi_index = n
c
c  Use up the old value, if we have it.
c
      if ( saved .eq. 1 ) then
        x(1) = y
        saved = 0
        x_lo_index = 2
      end if
c
c  Maybe we don't need any more values.
c
      if ( x_hi_index - x_lo_index + 1 .eq. 0 ) then
c
c  If we need just one new value, do that here to avoid null arrays.
c
      else if ( x_hi_index - x_lo_index + 1 .eq. 1 ) then

        r(1) = r8_uniform_01 ( seed )

        if ( r(1) .eq. 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8VEC_NORMAL_01 - Fatal errorc'
          write ( *, '(a)' ) '  R8_UNIFORM_01 returned a value of 0.'
          stop
        end if

        r(2) = r8_uniform_01 ( seed )

        x(x_hi_index) =
     &           sqrt ( -2.0D+00 * log ( r(1) ) )
     &           * cos ( 2.0D+00 * pi * r(2) )
        y =      sqrt ( -2.0D+00 * log ( r(1) ) )
     &           * sin ( 2.0D+00 * pi * r(2) )

        saved = 1

        made = made + 2
c
c  If we require an even number of values, that's easy.
c
      else if ( mod ( x_hi_index - x_lo_index + 1, 2 ) .eq. 0 ) then

        do i = x_lo_index, x_hi_index, 2

          call r8vec_uniform_01 ( 2, seed, r )

          x(i) =
     &      sqrt ( -2.0D+00 * log ( r(1) ) )
     &      * cos ( 2.0D+00 * pi * r(2) )

          x(i+1) =
     &      sqrt ( -2.0D+00 * log ( r(1) ) )
     &      * sin ( 2.0D+00 * pi * r(2) )

        end do

        made = made + x_hi_index - x_lo_index + 1
c
c  If we require an odd number of values, we generate an even number,
c  and handle the last pair specially, storing one in X(N), and
c  saving the other for later.
c
      else

        do i = x_lo_index, x_hi_index - 1, 2

          call r8vec_uniform_01 ( 2, seed, r )

          x(i) =
     &      sqrt ( -2.0D+00 * log ( r(1) ) )
     &      * cos ( 2.0D+00 * pi * r(2) )

          x(i+1) =
     &      sqrt ( -2.0D+00 * log ( r(1) ) )
     &      * sin ( 2.0D+00 * pi * r(2) )

        end do

        call r8vec_uniform_01 ( 2, seed, r )

        x(n) = sqrt ( -2.0D+00 * log ( r(1) ) )
     &    * cos ( 2.0D+00 * pi * r(1) )

        y = sqrt ( -2.0D+00 * log ( r(2) ) )
     &    * sin ( 2.0D+00 * pi * r(2) )

        saved = 1

        made = made + x_hi_index - x_lo_index + 2

      end if

      return
      end
      function r8vec_sum ( n, v1 )

c*********************************************************************72
c
cc R8VEC_SUM sums the entries of an R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    In FORTRAN90, the system routine SUM should be called
c    directly.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 July 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the vectors.
c
c    Input, double precision V1(N), the vector.
c
c    Output, double precision R8VEC_SUM, the sum of the entries.
c
      implicit none

      integer n

      integer i
      double precision r8vec_sum
      double precision v1(n)
      double precision value

      value = 0.0D+00
      do i = 1, n
        value = value + v1(i)
      end do

      r8vec_sum = value

      return
      end
      subroutine r8vec_uniform_01 ( n, seed, r )

c*********************************************************************72
c
cc R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
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
c    Input, integer N, the number of entries in the vector.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R(N), the vector of pseudorandom values.
c
      implicit none

      integer n

      integer i
      integer k
      integer seed
      double precision r(n)

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed .lt. 0 ) then
          seed = seed + 2147483647
        end if

        r(i) = dble ( seed ) * 4.656612875D-10

      end do

      return
      end
      subroutine stab_asymptotic ( seed, n, p_max )

c*********************************************************************72
c
cc STAB_ASYMPTOTIC examines asymptotic stability.
c
c  Discussion:
c
c    The function tests the asymptotic stability
c    of the Euler-Maruyama method applied to a stochastic differential
c    equation (SDE).
c
c    The SDE is
c
c      dX = lambda*X dt + mu*X dW,
c      X(0) = Xzero,
c
c    where 
c
c      lambda is a constant,
c      mu is a constant,
c      Xzero = 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 September 2012
c
c  Author:
c
c    Original Matlab version by Desmond Higham.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Desmond Higham,
c    An Algorithmic Introduction to Numerical Simulation of
c    Stochastic Differential Equations,
c    SIAM Review,
c    Volume 43, Number 3, September 2001, pages 525-546.
c
c  Parameters:
c
c    Input/output, integer SEED, a seed for the random 
c    number generator.
c
c    Input, integer N, the number of time steps for the
c    first solution.
c
c    Input, integer P_MAX, the number of time step sizes.
c
      implicit none

      integer n
      integer p_max

      character * ( 80 ) command_filename
      integer command_unit
      character * ( 80 ) data_filename
      integer data_unit
      double precision dt
      double precision dtvals(p_max)
      integer i
      integer j
      double precision lambda
      double precision mu
      integer nval
      integer p
      double precision r8_normal_01
      integer seed
      double precision t
      double precision test
      double precision tmax
      double precision u(1000)
      double precision winc
      double precision xemabs(n*2**(p_max-1))
      double precision xmin
      double precision xtemp
      double precision xzero

      data_filename = 'stab_asymptotic0_data.txt'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'STAB_ASYMPTOTIC:'
      write ( *, '(a)' ) 
     &  '  Investigate asymptotic stability of Euler-Maruyama'
      write ( *, '(a)' ) '  solution with stepsize DT and MU.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  SDE is asymptotically stable if'
      write ( *, '(a)' ) '    Real ( lambda - 1/2 mu^2 ) < 0.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  EM with DT is asymptotically stable if'
      write ( *, '(a)' ) 
     &  '    E log ( | 1 + lambda dt - sqrt(dt) mu n(0,1) | ) < 0.'
      write ( *, '(a)' ) '  where n(0,1) is a normal random value.'
c
c  Problem parameters.
c
      lambda = 0.5D+00
      mu = sqrt ( 6.0D+00 )
      xzero = 1.0D+00
c
c  Test the SDE.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Lambda = ', lambda
      write ( *, '(a,g14.6)' ) '  Mu =     ', mu
      test = lambda - 0.5D+00 * mu ** 2
      write ( *, '(a,g14.6)' ) 
     &  '  SDE asymptotic stability test = ', test 
c
c  Step parameters.
c
      tmax = 500.0D+00
c
c  For each stepsize, compute the Euler-Maruyama solution.
c
      do p = 1, p_max

        nval = n * 2 ** ( p - 1 )
        dt = tmax / dble ( nval )
        dtvals(p) = dt
c
c  Test the EM for this DT.
c 
        write ( *, '(a)' ) ' '           
        write ( *, '(a,g14.6)' ) '  dt = ', dt
        call r8vec_normal_01 ( 1000, seed, u )
        test = 0.0D+00
        do i = 1, 1000
          u(i) = log ( 
     &      abs ( 1.0D+00 + lambda * dt - sqrt ( dt ) * mu * u(i) ) )
          test = test + u(i)
        end do
        test = test / 1000.0D+00
        write ( *, '(a,g14.6)' ) '  EM asymptotic test = ', test

        xtemp = xzero
        xemabs(0) = xtemp

        do j = 1, nval
          winc = sqrt ( dt ) * r8_normal_01 ( seed )
          xtemp = xtemp + dt * lambda * xtemp + mu * xtemp * winc
          xemabs(j) = abs ( xtemp )
        end do
c
c  Write this data to a file.
c
        call get_unit ( data_unit )

        call filename_inc ( data_filename )

        open ( unit = data_unit, file = data_filename, 
     &    status = 'replace' )
c
c  We have to impose a tiny lower bound on the values because we
c  will end up plotting their logs.
c
        xmin = exp ( -200.0D+00 )
        do i = 0, nval
          t = tmax * dble ( i ) / dble ( nval )
          write ( data_unit, '(2x,g14.6,2x,g14.6)' ) 
     &      t, max ( xemabs(i), xmin )
        end do
        close ( unit = data_unit )

        write ( *, '(a)' ) ' '
        write ( *, '(a,g14.6,a)' ) '  Data for DT = ', dt, 
     &    ' stored in "' // trim ( data_filename ) // '".'

      end do
c
c  Create the command file.
c
      call get_unit ( command_unit )

      command_filename = 'stab_asymptotic_commands.txt'

      open ( unit = command_unit, file = command_filename, 
     &  status = 'replace' )

      write ( command_unit, '(a)' ) '# stab_asymptotic_commands.txt'
      write ( command_unit, '(a)' ) 
     &  '# created by sde::stab_asymptotic.'
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) '# Usage:'
      write ( command_unit, '(a)' ) 
     &  '#  gnuplot < stab_asymptotic_commands.txt'
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) 'set term png'
      write ( command_unit, '(a)' ) 'set output "stab_asymptotic.png"'
      write ( command_unit, '(a)' ) 'set xlabel "t"'
      write ( command_unit, '(a)' ) 'set ylabel "|X(t)|"'
      write ( command_unit, '(a)' ) 
     &  'set title "Absolute value of EM Solution"'
      write ( command_unit, '(a)' ) 'set grid'
      write ( command_unit, '(a)' ) 'set logscale y 10'
      write ( command_unit, '(a)' ) 'set style data lines'

      data_filename = 'stab_asymptotic0_data.txt'

      call filename_inc ( data_filename )
      write ( command_unit, '(a)' ) 
     &  'plot "' // trim ( data_filename ) // '" using 1:2, \'

      do p = 2, p_max - 1
        call filename_inc ( data_filename )
        write ( command_unit, '(a)' ) 
     &    '     "' // trim ( data_filename ) // '" using 1:2, \'
      end do
      call filename_inc ( data_filename )
      write ( command_unit, '(a)' ) 
     &  '     "' // trim ( data_filename ) // '" using 1:2'

      write ( command_unit, '(a)' ) 'quit'

      close ( unit = command_unit )

      write ( *, '(a)' ) '  STAB_ASYMPTOTIC plot commands stored in "' 
     &  // trim ( command_filename ) // '".'

      return
      end
      subroutine stab_meansquare ( seed )

c*********************************************************************72
c
cc STAB_MEANSQUARE examines mean-square stability.
c
c  Discussion:
c
c    The function tests the mean-square stability
c    of the Euler-Maruyama method applied to a stochastic differential
c    equation (SDE).
c
c    The SDE is
c
c      dX = lambda*X dt + mu*X dW,
c      X(0) = Xzero,
c
c    where 
c
c      lambda is a constant,
c      mu is a constant,
c      Xzero = 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 September 2012
c
c  Author:
c
c    Original Matlab version by Desmond Higham.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Desmond Higham,
c    An Algorithmic Introduction to Numerical Simulation of
c    Stochastic Differential Equations,
c    SIAM Review,
c    Volume 43, Number 3, September 2001, pages 525-546.
c
c  Parameters:
c
c    Input, integer SEED, a seed for the random number generator.
c    In the reference, this value is set to 100.
c
      implicit none

      integer m
      parameter ( m = 50000 )

      character * ( 80 ) command_filename
      integer command_unit
      character * ( 80 ) data_filename
      integer data_unit
      double precision dt
      integer i
      integer j
      integer k
      double precision lambda
      double precision mu
      integer n
      integer seed
      double precision t
      double precision test
      double precision tmax
      double precision winc(m)
      double precision xms(0:80)
      double precision xtemp(m)
      double precision xzero

      data_filename = 'stab_meansquare0_data.txt'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'STAB_MEANSQUARE:'
      write ( *, '(a)' ) 
     &  '  Investigate mean square stability of Euler-Maruyama'
      write ( *, '(a)' ) '  solution with stepsize DT and MU.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  SDE is mean square stable if'
      write ( *, '(a)' ) '    Real ( lambda + 1/2 |mu|^2 ) < 0.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  EM with DT is mean square stable if'
      write ( *, '(a)' ) '    |1+dt^2| + dt * |mu|^2 - 1.0 < 0.'
c
c  Set problem parameters.
c
      tmax = 20.0D+00
      xzero = 1.0D+00 
c
c  Problem parameters.
c
      lambda = -3.0D+00
      mu = sqrt ( 3.0D+00 )
c
c  Test the SDE.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Lambda = ', lambda
      write ( *, '(a,g14.6)' ) '  Mu =     ', mu
      test = lambda + 0.5D+00 * mu ** 2
      write ( *, '(a,g14.6)' ) 
     &  '  SDE mean square stability test = ', test
c
c  XMS is the mean square estimate of M paths.
c
      do k = 1, 3

        dt = 2.0D+00 ** ( 1 - k )                      
        n = 20 * 2 ** ( k - 1 )
c
c  Test the EM for this DT.
c 
        write ( *, '(a)' ) ' '           
        write ( *, '(a,g14.6)' ) '  dt = ', dt
        test = ( 1.0D+00 + dt * lambda ) ** 2 + dt * mu ** 2 - 1.0D+00
        write ( *, '(a,g14.6)' ) 
     &    '  EM mean square stability test = ', test

        do i = 1, m
          xtemp(i) = xzero
        end do

        xms(0) = xzero

        do j = 1, n
          call r8vec_normal_01 ( m, seed, winc )
          do i = 1, m
            winc(i) = sqrt ( dt ) * winc(i) 
            xtemp(i) = xtemp(i) 
     &        + dt * lambda * xtemp(i) 
     &        + mu * xtemp(i) * winc(i)
          end do
          xms(j) = 0.0D+00
          do i = 1, m
            xms(j) = xms(j) + xtemp(i) ** 2
          end do
          xms(j) = xms(j) / dble ( m )
        end do
c
c  Write this data to a file.
c
        call get_unit ( data_unit )

        call filename_inc ( data_filename )

        open ( unit = data_unit, file = data_filename, 
     &    status = 'replace' )

        do j = 0, n
          t = tmax * dble ( j ) / dble ( n )
          write ( data_unit, '(2x,g14.6,2x,g14.6)' ) t, xms(j)
        end do
        close ( unit = data_unit )

        write ( *, '(a)' ) ' '
        write ( *, '(a,g14.6,a)' ) '  Data for DT = ', dt, 
     &    ' stored in "' // trim ( data_filename ) // '".'

      end do
c
c  Create the command file.
c
      call get_unit ( command_unit )

      command_filename = 'stab_meansquare_commands.txt'

      open ( unit = command_unit, file = command_filename, 
     &  status = 'replace' )

      write ( command_unit, '(a)' ) '# stab_meansquare_commands.txt'
      write ( command_unit, '(a)' ) '# created by sde::stab_meansquare.'
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) '# Usage:'
      write ( command_unit, '(a)' ) 
     &  '#  gnuplot < stab_meansquare_commands.txt'
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) 'set term png'
      write ( command_unit, '(a)' ) 'set output "stab_meansquare.png"'
      write ( command_unit, '(a)' ) 'set xlabel "t"'
      write ( command_unit, '(a)' ) 'set ylabel "E|X^2(t)|"'
      write ( command_unit, '(a)' ) 
     &  'set title "Mean Square of EM Solution"'
      write ( command_unit, '(a)' ) 'set grid'
      write ( command_unit, '(a)' ) 'set logscale y 10'
      write ( command_unit, '(a)' ) 'set style data lines'

      data_filename = 'stab_meansquare0_data.txt'

      call filename_inc ( data_filename )
      write ( command_unit, '(a)' ) 
     &  'plot "' // trim ( data_filename ) // '" using 1:2, \'

      do k = 2, 2
        call filename_inc ( data_filename )
        write ( command_unit, '(a)' ) 
     &    '     "' // trim ( data_filename ) // '" using 1:2, \'
      end do
      call filename_inc ( data_filename )
      write ( command_unit, '(a)' ) 
     &  '     "' // trim ( data_filename ) // '" using 1:2'

      write ( command_unit, '(a)' ) 'quit'

      close ( unit = command_unit )

      write ( *, '(a)' ) '  STAB_MEANSQUARE plot commands stored in "' 
     &  // trim ( command_filename ) // '".'

      return
      end
      subroutine stochastic_integral_ito ( n, seed, estimate, exact, 
     &  error )

c*********************************************************************72
c
cc STOCHASTIC_INTEGRAL_ITO approximates the Ito integral of W(t) dW.
c
c  Discussion:
c
c    This function estimates the Ito integral of W(t) dW over 
c    the interval [0,1].
c
c    The estimates is made by taking N steps.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 September 2012
c
c  Author:
c
c    Original Matlab version by Desmond Higham.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Desmond Higham,
c    An Algorithmic Introduction to Numerical Simulation of
c    Stochastic Differential Equations,
c    SIAM Review,
c    Volume 43, Number 3, September 2001, pages 525-546.
c
c  Parameters:
c
c    Input, integer N, the number of steps to take.
c
c    Input, integer SEED, a seed for the random number generator.
c
c    Output, double precision ESTIMATE, the estimate of the integral.
c
c    Output, double precision EXACT, the exact value of the integral.
c
c    Output, double precision ERROR, the error in the integral estimate.
c
      implicit none

      integer n

      double precision dt
      double precision dw(n)
      double precision error
      double precision estimate
      double precision exact
      integer j
      integer seed
      double precision tmax
      double precision w(0:n)
c
c  Set step parameters.
c
      tmax = 1.0D+00
      dt = tmax / dble ( n )
c
c  Define the increments dW.
c
      call r8vec_normal_01 ( n, seed, dw )

      do j = 1, n
        dw(j) = sqrt ( dt ) * dw(j)
      end do
c
c  Sum the increments to get the Brownian path.
c
      w(0) = 0.0D+00
      do j = 1, n
        w(j) = w(j-1) + dw(j)
      end do
c
c  Approximate the Ito integral.
c
      estimate = dot_product ( w, dw )
c
c  Compare with the exact solution.
c
      exact = 0.5D+00 * ( w(n) ** 2 - tmax )
      error = abs ( estimate - exact )

      return
      end
      subroutine stochastic_integral_strat ( n, seed, estimate, 
     &  exact, error )

c*********************************************************************72
c
cc STOCHASTIC_INTEGRAL_STRAT approximates the Stratonovich integral of W(t) dW.
c
c  Discussion:
c
c    This function estimates the Stratonovich integral of W(t) dW over 
c    the interval [0,1].
c
c    The estimates is made by taking N steps.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 September 2012
c
c  Author:
c
c    Original Matlab version by Desmond Higham.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Desmond Higham,
c    An Algorithmic Introduction to Numerical Simulation of
c    Stochastic Differential Equations,
c    SIAM Review,
c    Volume 43, Number 3, September 2001, pages 525-546.
c
c  Parameters:
c
c    Input, integer N, the number of steps to take.
c
c    Input, integer SEED, a seed for the random number generator.
c
c    Output, double precision ESTIMATE, the estimate of the integral.
c
c    Output, double precision EXACT, the exact value of the integral.
c
c    Output, double precision ERROR, the error in the integral estimate.
c
      implicit none

      integer n

      double precision dt
      double precision dw(n)
      double precision error
      double precision estimate
      double precision exact
      integer j
      double precision r8vec_dot_product
      integer seed
      double precision tmax
      double precision u(n)
      double precision v(n)
      double precision w(0:n)
c
c  Set step parameters.
c
      tmax = 1.0D+00
      dt = tmax / dble ( n )
c
c  Define the increments dW.
c
      call r8vec_normal_01 ( n, seed, dw )

      do j = 1, n
        dw(j) = sqrt ( dt ) * dw(j)
      end do
c
c  Sum the increments to get the Brownian path.
c
      w(0) = 0.0D+00
      do j = 1, n
        w(j) = w(j-1) + dw(j)
      end do
c
c  Approximate the Stratonovich integral.
c
      call r8vec_normal_01 ( n, seed, u )
      do j = 1, n
        v(j) = 0.5D+00 * ( w(j-1) + w(j) ) 
     &    + 0.5D+00 * sqrt ( dt ) * u(j)
      end do

      estimate = r8vec_dot_product ( n, v, dw )
c
c  Compare with the exact solution.
c
      exact = 0.5D+00 * w(n) ** 2
      error = abs ( estimate - exact )

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
