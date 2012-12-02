      program main

c*********************************************************************72
c
cc MAIN is the main program for FEYNMAN_KAC_2D.
c
c  Discussion:
c
c    This program is derived from section 2.5, exercise 2.2 of Petersen 
c    and Arbenz.
c
c    The problem is to determine the solution U(X,Y) of the following 
c    partial differential equation:
c
c      (1/2) Laplacian U - V(X,Y) * U = 0,
c
c    inside the elliptic domain D:
c 
c      D = { (X,Y) | (X/A)^2+(Y/B)^2 <= 1 }
c   
c    with the boundary condition U(boundary(D)) = 1.
c
c    V(X,Y) is the potential function:
c
c      V = 2 * ( (X/A^2)^2 + (Y/B^2)^2 ) + 1/A^2 + 1/B^2.
c
c    The analytic solution of this problem is already known:
c
c      U(X,Y,Z) = exp ( (X/A)^2 + (Y/B)^2 - 1 ).
c
c    Our method is via the Feynman-Kac Formula.
c
c    The idea is to start from any (x,y) in D, and
c    compute (x+Wx(t),y+Wy(t)) where 2D Brownian motion
c    (Wx,Wy) is updated each step by sqrt(h)*(z1,z2),
c    each z1,z2 are independent approximately Gaussian 
c    random variables with zero mean and variance 1. 
c
c    Each (x1(t),x2(t) ) is advanced until (x1,x2 ) exits the domain D.  
c
c    Upon its first exit from D, the sample path (x1,x2) is stopped and a 
c    new sample path at (x,y) is started until N such paths are completed.
c 
c    The Feynman-Kac formula gives the solution here as
c
c      U(X,Y) = (1/N) sum(1 <= I <= N) Y(tau_i),
c
c    where
c
c      Y(tau) = exp( -int(s=0..tau) v(x1(s),x2(s)) ds),
c
c    and tau = first exit time for path (x1,x2). 
c
c    The integration procedure is a second order weak accurate method:
c
c      X(t+h)  = [ x1(t) + sqrt ( h ) * z1 ]
c                [ x2(t) + sqrt ( h ) * z2 ]
c
c    Here Z1 and Z2 are approximately normal univariate Gaussians. 
c
c    An Euler predictor approximates Y at the end of the step
c
c      Y_e     = (1 - h*v(X(t)) * Y(t), 
c
c    A trapezoidal rule completes the step:
c
c      Y(t+h)  = Y(t) - (h/2)*[v(X(t+h))*Y_e + v(X(t))*Y(t)].
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 May 2012
c
c  Author:
c
c    Original C 3D version by Wesley Petersen.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Peter Arbenz, Wesley Petersen,
c    Introduction to Parallel Computing:
c    A Practical Guide with Examples in C,
c    Oxford, 2004,
c    ISBN: 0-19-851577-4,
c    LC: QA76.59.P47.
c
      implicit none

      double precision a
      parameter ( a = 2.0D+00 )
      double precision chk
      double precision dx
      double precision err
      double precision h
      parameter ( h = 0.0001D+00 )
      integer i
      integer it
      integer j
      integer k
      integer n1
      integer n
      parameter ( n = 10000 )
      integer n_int
      integer ni
      double precision r8_uniform_01
      double precision rth
      integer seed
      integer steps
      integer steps_ave
      double precision sum
      double precision test
      double precision us
      double precision vh
      double precision vs
      double precision x
      double precision x1
      double precision w
      double precision w_exact
      double precision we
      double precision wt

      seed = 123456789

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEYNMAN_KAC_1D:'
      write ( *, '(a)' ) '  FORTRAN77 version.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Program parameters:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  The calculation takes place inside an interval.'
      write ( *, '(a)' ) '  The solution will be estimated at NG points'
      write ( *, '(a)' ) 
     &  '  on a regular spaced grid within the interval.'
      write ( *, '(a,i8,a)' ) 
     &  '  Each solution will be estimated by computing ', n, 
     &  ' trajectories'
      write ( *, '(a)' ) '  from the point to the boundary.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    (X/A)^2 = 1'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The interval parameter A is:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '    A = ', a
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Path stepsize H = ', h
c
c  Choose the spacing so we have about ni points on or in the interval.
c
      ni = 21

      write ( *, '(a)' ) ' '
      write ( *, '(a,i4,a)' ) '  X coordinate discretized by ', ni+2, 
     &  ' points'
c
c  RTH is the scaled stepsize.
c
      rth = sqrt ( h )

      err = 0.0D+00
c
c  Loop over the points.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '     I     K       X           W exact' // 
     &  '      W Approx        Error      Ave Steps  Test'
      write ( *, '(a)' ) ' '

      k = 0
      n_int = 0

      do i = 0, ni + 1

        x = ( dble ( ni - i     ) * ( - a ) 
     &      + dble (      i - 1 ) *     a ) 
     &      / dble ( ni     - 1 )

        k = k + 1

        test = a**2 - x**2

        if ( test .lt. 0.0D+00  ) then
          w_exact = 1.0D+00
          wt = 1.0D+00
          steps_ave = 0
          write ( *, 
     &      '(2x,i4,2x,i4,4(2x,g12.4),2x,i8,2x,f8.4)' ) 
     &      i, k, x, w_exact, wt, abs ( w_exact - wt ), 
     &      steps_ave, test
          go to 20
        end if

        n_int = n_int + 1
c
c  Compute the exact solution at this point (x,y,z).
c
        w_exact = exp ( ( x / a )**2 - 1.0D+00 )
c
c  Now try to estimate the solution at this point.
c
        wt = 0.0D+00
        steps = 0

        do it = 1, n

          x1 = x
c 
c  W = exp(-int(s=0..t) v(X)ds) 
c
          w = 1.0D+00
c
c  CHK is < 1.0 while the point is inside the ellipsoid.
c
          chk = 0.0D+00

10        continue

          if ( chk .lt. 1.0D+00 ) then
c
c  Determine DX.
c
            us = r8_uniform_01 ( seed ) - 0.5D+00
            if ( us .lt. 0.0D+00 ) then
              dx = - rth
            else
              dx = + rth
            end if

            call potential ( a, x1, vs )
c
c  Move to the new point.
c
            x1 = x1 + dx

            steps = steps + 1

            call potential ( a, x1, vh )

            we = ( 1.0D+00 - h * vs ) * w
            w = w - 0.5D+00 * h * ( vh * we + vs * w ) 

            chk = ( x1 / a )**2

            go to 10

          end if

          wt = wt + w

        end do
c
c  WT is the average of the sum of the different trials.
c
        wt = wt / dble ( n )
        steps_ave = steps / n
c
c  Add error in WT to the running L2 error in the solution.
c
        err = err + ( w_exact - wt )**2

        write ( *, 
     &    '(2x,i4,2x,i4,4(2x,g12.4),2x,i8,2x,f8.4)' ) 
     &    i, k, x, w_exact, wt, abs ( w_exact - wt ), steps_ave, test

20      continue

      end do
c
c  Compute the RMS error for all the points.
c
      err = sqrt ( err / dble ( n_int ) )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) 
     &  '  RMS absolute error in solution = ', err
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEYNMAN_KAC_1D:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '

      call timestamp ( )

      stop
      end
      subroutine potential ( a, x, v )

c*********************************************************************72
c
cc POTENTIAL evaluates the potential function V(X,).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A, the parameters that define the ellipse.
c
c    Input, double precision X, the coordinates of the point.
c
c    Output, double precision V, the value of the potential function at X.
c
      implicit none

      double precision a
      double precision v
      double precision x

      v = 2.0D+00 * ( x / a**2 )**2 + 1.0D+00 / a**2

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
c    Second Edition,
c    Springer, 1987,
c    ISBN: 0387964673,
c    LC: QA76.9.C65.B73.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, December 1986, pages 362-376.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley, 1998,
c    ISBN: 0471134031,
c    LC: T57.62.H37.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, Number 2, 1969, pages 136-143.
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

      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer k
      double precision r8_uniform_01
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
        seed = seed + i4_huge
      end if

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
