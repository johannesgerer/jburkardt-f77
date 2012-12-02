      subroutine correlation_besselj ( n, rho, rho0, c )

c*********************************************************************72
c
cc CORRELATION_BESSELJ evaluates the Bessel J correlation function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Petter Abrahamsen,
c    A Review of Gaussian Random Fields and Correlation Functions,
c    Norwegian Computing Center, 1997.
c
c  Parameters:
c
c    Input, integer N, the number of arguments.
c
c    Input, double precision RHO(N), the arguments.
c
c    Input, double precision RHO0, the correlation length.
c
c    Output, double precision C(N), the correlations.
c
      implicit none

      integer n

      double precision c(n)
      integer i
      double precision r8_besj0
      double precision rho(n)
      double precision rho0
      double precision rhohat

      do i = 1, n
        rhohat = abs ( rho(i) ) / rho0
        c(i) = r8_besj0 ( rhohat )
      end do

      return
      end
      subroutine correlation_besselk ( n, rho, rho0, c )

c*********************************************************************72
c
cc CORRELATION_BESSELK evaluates the Bessel K correlation function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Petter Abrahamsen,
c    A Review of Gaussian Random Fields and Correlation Functions,
c    Norwegian Computing Center, 1997.
c
c  Parameters:
c
c    Input, integer N, the number of arguments.
c
c    Input, double precision RHO(N), the arguments.
c
c    Input, double precision RHO0, the correlation length.
c
c    Output, double precision C(N), the correlations.
c
      implicit none

      integer n

      double precision c(n)
      integer i
      double precision r8_besk1
      double precision rho(n)
      double precision rho0
      double precision rhohat

      do i = 1, n

        if ( rho(i) .eq. 0.0D+00 ) then
          c(i) = 1.0D+00
        else
          rhohat = abs ( rho(i) ) / rho0
          c(i) = rhohat * r8_besk1 ( rhohat )
        end if

      end do

      return
      end
      subroutine correlation_brownian ( m, n, s, t, rho0, c )

c*********************************************************************72
c
cc CORRELATION_BROWNIAN computes the Brownian correlation function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of arguments.
c
c    Input, double precision S(M), T(N), two samples.
c    0 <= S(*), T(*).
c
c    Input, double precision RHO0, the correlation length.
c
c    Output, double precision C(M,N), the correlations.
c
      implicit none

      integer m
      integer n

      double precision c(m,n)
      integer i
      integer j
      double precision rho0
      double precision s(n)
      double precision t(n)

      do j = 1, n
        do i = 1, m
          if ( 0.0D+00 .lt. max ( s(i), t(j) ) ) then
            c(i,j) = sqrt ( min ( s(i), t(j) ) / max ( s(i), t(j) ) )
          else
            c(i,j) = 1.0D+00
          end if
        end do
      end do

      return
      end
      subroutine correlation_brownian_display ( )

c*********************************************************************72
c
cc CORRELATION_BROWNIAN_DISPLAY displays 4 slices of the Brownian Correlation.
c
c  Discussion:
c
c    The correlation function is C(S,T) = sqrt ( min ( s, t ) / max ( s, t ) ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 November 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 101 )
      integer n2
      parameter ( n2 = 4 )

      double precision c(n,n2)
      character * ( 80 ) command_filename
      integer command_unit
      character * ( 80 ) data_filename
      integer data_unit
      integer i
      integer j
      double precision s(n)
      double precision t(n2)

      call r8vec_linspace ( n, 0.0D+00, 5.0D+00, s )

      t(1) = 0.25D+00
      t(2) = 1.50D+00
      t(3) = 2.50D+00
      t(4) = 3.75D+00

      do i = 1, n
        do j = 1, n2
          c(i,j) = sqrt ( min ( s(i), t(j) ) / max ( s(i), t(j) ) )
        end do
      end do

      call get_unit ( data_unit )
      data_filename = 'brownian_plots_data.txt'
      open ( unit = data_unit, file = data_filename, 
     &  status = 'replace' )
      do i = 1, n
        write ( data_unit, '(5(2x,g14.6))' ) s(i), c(i,1:n2)
      end do
      close ( unit = data_unit )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Created data file "' 
     &  // trim ( data_filename ) // '".'

      call get_unit ( command_unit )
      command_filename = 'brownian_plots_commands.txt'
      open ( unit = command_unit, file = command_filename, 
     &  status = 'replace' )
      write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) '# Usage:'
      write ( command_unit, '(a)' ) '#  gnuplot < ' 
     &  // trim ( command_filename )
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) 'set term png'
      write ( command_unit, '(a)' ) 'set key off'
      write ( command_unit, '(a)' ) 
     &  'set output "brownian_plots.png"'
      write ( command_unit, '(a)' ) 
     &  'set title "Brownian correlation C(S,T), ' // 
     &  'S = 0.25, 1.5, 2.5, 3.75"'
      write ( command_unit, '(a)' ) 'set xlabel "S"'
      write ( command_unit, '(a)' ) 'set ylabel "C(s,t)"'
      write ( command_unit, '(a)' ) 'set grid'
      write ( command_unit, '(a)' ) 'set style data lines'
      write ( command_unit, '(a)' ) 'plot "' 
     &  // trim ( data_filename ) // 
     &  '" using 1:2 lw 3 linecolor rgb "blue", \'
      write ( command_unit, '(a)' ) '     "' 
     &  // trim ( data_filename ) // 
     &  '" using 1:3 lw 3 linecolor rgb "blue", \'
      write ( command_unit, '(a)' ) '     "' 
     &  // trim ( data_filename ) // 
     &  '" using 1:4 lw 3 linecolor rgb "blue", \'
      write ( command_unit, '(a)' ) '     "' 
     &  // trim ( data_filename ) // 
     &  '" using 1:5 lw 3 linecolor rgb "blue"'
      write ( command_unit, '(a)' ) 'quit'
      close ( unit = command_unit )
      write ( *, '(a)' ) 
     &  '  Created command file "' 
     &  // trim ( command_filename ) // '".'

      return
      end
      subroutine correlation_circular ( n, rho, rho0, c )

c*********************************************************************72
c
cc CORRELATION_CIRCULAR evaluates the circular correlation function.
c
c  Discussion:
c
c    This correlation is based on the area of overlap of two circles
c    of radius RHO0 and separation RHO.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Petter Abrahamsen,
c    A Review of Gaussian Random Fields and Correlation Functions,
c    Norwegian Computing Center, 1997.
c
c  Parameters:
c
c    Input, integer N, the number of arguments.
c
c    Input, double precision RHO(N), the arguments.
c
c    Input, double precision RHO0, the correlation length.
c
c    Output, double precision C(N), the correlations.
c
      implicit none

      integer n

      double precision c(n)
      integer i
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision rho(n)
      double precision rho0
      double precision rhohat

      do i = 1, n

        rhohat = min ( abs ( rho(i) ) / rho0, 1.0D+00 )

        c(i) = ( 1.0D+00 - ( 2.0D+00 / pi ) 
     &    * ( rhohat * sqrt ( 1.0D+00 - rhohat ** 2 ) 
     &    + asin ( rhohat ) ) )

      end do

      return
      end
      subroutine correlation_constant ( n, rho, rho0, c )

c*********************************************************************72
c
cc CORRELATION_CONSTANT evaluates the constant correlation function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Petter Abrahamsen,
c    A Review of Gaussian Random Fields and Correlation Functions,
c    Norwegian Computing Center, 1997.
c
c  Parameters:
c
c    Input, integer N, the number of arguments.
c
c    Input, double precision RHO(N), the arguments.
c
c    Input, double precision RHO0, the correlation length.
c
c    Output, double precision C(N), the correlations.
c
      implicit none

      integer n

      double precision c(n)
      integer i
      double precision rho(n)
      double precision rho0

      do i = 1, n
        c(i) = 1.0D+00
      end do

      return
      end
      subroutine correlation_cubic ( n, rho, rho0, c )

c*********************************************************************72
c
cc CORRELATION_CUBIC evaluates the cubic correlation function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Petter Abrahamsen,
c    A Review of Gaussian Random Fields and Correlation Functions,
c    Norwegian Computing Center, 1997.
c
c  Parameters:
c
c    Input, integer N, the number of arguments.
c
c    Input, double precision RHO(N), the arguments.
c
c    Input, double precision RHO0, the correlation length.
c
c    Output, double precision C(N), the correlations.
c
      implicit none

      integer n

      double precision c(n)
      integer i
      double precision rho(n)
      double precision rho0
      double precision rhohat

      do i = 1, n

        rhohat = min ( abs ( rho(i) ) / rho0, 1.0D+00 )

        c(i) = 1.0D+00 
     &       - 7.0D+00  * rhohat ** 2 
     &       + 8.75D+00 * rhohat ** 3 
     &       - 3.5D+00  * rhohat ** 5 
     &       + 0.75D+00 * rhohat ** 7

      end do

      return
      end
      subroutine correlation_damped_cosine ( n, rho, rho0, c )

c*********************************************************************72
c
cc CORRELATION_DAMPED_COSINE evaluates the damped cosine correlation function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Petter Abrahamsen,
c    A Review of Gaussian Random Fields and Correlation Functions,
c    Norwegian Computing Center, 1997.
c
c  Parameters:
c
c    Input, integer N, the number of arguments.
c
c    Input, double precision RHO(N), the arguments.
c
c    Input, double precision RHO0, the correlation length.
c
c    Output, double precision C(N), the correlations.
c
      implicit none

      integer n

      double precision c(n)
      integer i
      double precision rho(n)
      double precision rho0

      do i = 1, n
        c(i) = exp ( - abs ( rho(i) ) / rho0 ) 
     &    * cos ( abs ( rho(i) ) / rho0 )
      end do

      return
      end
      subroutine correlation_damped_sine ( n, rho, rho0, c )

c*********************************************************************72
c
cc CORRELATION_DAMPED_SINE evaluates the damped sine correlation function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Petter Abrahamsen,
c    A Review of Gaussian Random Fields and Correlation Functions,
c    Norwegian Computing Center, 1997.
c
c  Parameters:
c
c    Input, integer N, the number of arguments.
c
c    Input, double precision RHO(N), the arguments.
c
c    Input, double precision RHO0, the correlation length.
c
c    Output, double precision C(N), the correlations.
c
      implicit none

      integer n

      double precision c(n)
      integer i
      double precision rho(n)
      double precision rho0
      double precision rhohat

      do i = 1, n

        if ( rho(i) .eq. 0.0D+00 ) then
          c(i) = 1.0D+00
        else
          rhohat = abs ( rho(i) ) / rho0
          c(i) = sin ( rhohat ) / rhohat
        end if

      end do

      return
      end
      subroutine correlation_exponential ( n, rho, rho0, c )

c*********************************************************************72
c
cc CORRELATION_EXPONENTIAL evaluates the exponential correlation function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Petter Abrahamsen,
c    A Review of Gaussian Random Fields and Correlation Functions,
c    Norwegian Computing Center, 1997.
c
c  Parameters:
c
c    Input, integer N, the number of arguments.
c
c    Input, double precision RHO(N), the arguments.
c
c    Input, double precision RHO0, the correlation length.
c
c    Output, double precision C(N), the correlations.
c
      implicit none

      integer n

      double precision c(n)
      integer i
      double precision rho(n)
      double precision rho0

      do i = 1, n
        c(i) = exp ( - abs ( rho(i) ) / rho0 )
      end do

      return
      end
      subroutine correlation_gaussian ( n, rho, rho0, c )

c*********************************************************************72
c
cc CORRELATION_GAUSSIAN evaluates the Gaussian correlation function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Petter Abrahamsen,
c    A Review of Gaussian Random Fields and Correlation Functions,
c    Norwegian Computing Center, 1997.
c
c  Parameters:
c
c    Input, integer N, the number of arguments.
c
c    Input, double precision RHO(N), the arguments.
c
c    Input, double precision RHO0, the correlation length.
c
c    Output, double precision C(N), the correlations.
c
      implicit none

      integer n

      double precision c(n)
      integer i
      double precision rho(n)
      double precision rho0

      do i = 1, n
        c(i) = exp ( - ( ( rho(i) / rho0 ) ** 2 ) )
      end do

      return
      end
      subroutine correlation_hole ( n, rho, rho0, c )

c*********************************************************************72
c
cc CORRELATION_HOLE evaluates the hole correlation function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of arguments.
c
c    Input, double precision RHO(N), the arguments.
c
c    Input, double precision RHO0, the correlation length.
c
c    Output, double precision C(N), the correlations.
c
      implicit none

      integer n

      double precision c(n)
      integer i
      double precision rho(n)
      double precision rho0

      do i = 1, n
        c(i) = ( 1.0D+00 - abs ( rho(i) ) / rho0 ) 
     &    * exp ( - abs ( rho(i) ) / rho0 )
      end do

      return
      end
      subroutine correlation_linear ( n, rho, rho0, c )

c*********************************************************************72
c
cc CORRELATION_LINEAR evaluates the linear correlation function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Petter Abrahamsen,
c    A Review of Gaussian Random Fields and Correlation Functions,
c    Norwegian Computing Center, 1997.
c
c  Parameters:
c
c    Input, integer N, the number of arguments.
c
c    Input, double precision RHO(N), the arguments.
c
c    Input, double precision RHO0, the correlation length.
c
c    Output, double precision C(N), the correlations.
c
      integer n

      double precision c(n)
      integer i
      double precision rho(n)
      double precision rho0

      do i = 1, n
        if ( rho0 .lt. abs ( rho(i) ) ) then
          c(i) = 0.0D+00
        else
          c(i) = ( rho0 - abs ( rho(i) ) ) / rho0
        end if
      end do

      return
      end
      subroutine correlation_matern ( n, rho, rho0, c )

c*********************************************************************72
c
cc CORRELATION_MATERN evaluates the Matern correlation function.
c
c  Discussion:
c
c    In order to call this routine under a dummy name, I had to drop NU from
c    the parameter list.
c
c    The Matern correlation is
c
c      rho1 = 2 * sqrt ( nu ) * rho / rho0
c
c      c(rho) = ( rho1 )^nu * BesselK ( nu, rho1 ) 
c               / gamma ( nu ) / 2 ^ ( nu - 1 )
c
c    The Matern covariance has the form:
c
c      K(rho) = sigma^2 * c(rho)
c
c    A Gaussian process with Matern covariance has sample paths that are
c    differentiable (nu - 1) times.
c
c    When nu = 0.5, the Matern covariance is the exponential covariance.
c
c    As nu goes to +oo, the correlation converges to exp ( - (rho/rho0)^2 ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of arguments.
c
c    Input, double precision RHO(N), the arguments.
c    0.0 <= RHO.
c
c    Input, double precision RHO0, the correlation length.
c    0.0 < RHO0.
c
c    Input, double precision NU, the smoothness parameter.
c    NU has a default value of 2.5;
c    0 < NU.
c
c    Output, double precision C(N), the correlations.
c
      implicit none

      integer n

      double precision c(n)
      integer i
      double precision nu
      double precision r8_besk
      double precision r8_gamma
      double precision rho(n)
      double precision rho0
      double precision rho1

      nu = 2.5D+00

      do i = 1, n

        rho1 = 2.0D+00 * sqrt ( nu ) * abs ( rho(i) ) / rho0

        if ( rho1 .eq. 0.0D+00 ) then
          c(i) = 1.0D+00
        else
          c(i) = rho1 ** nu * r8_besk ( nu, rho1 ) / r8_gamma ( nu ) 
     &      / 2.0 ** ( nu - 1.0D+00 )
        end if

      end do

      return
      end
      subroutine correlation_pentaspherical ( n, rho, rho0, c )

c*********************************************************************72
c
cc CORRELATION_PENTASPHERICAL evaluates the pentaspherical correlation function.
c
c  Discussion:
c
c    This correlation is based on the volume of overlap of two spheres
c    of radius RHO0 and separation RHO.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Petter Abrahamsen,
c    A Review of Gaussian Random Fields and Correlation Functions,
c    Norwegian Computing Center, 1997.
c
c  Parameters:
c
c    Input, integer N, the number of arguments.
c
c    Input, double precision RHO(N), the arguments.
c
c    Input, double precision RHO0, the correlation length.
c
c    Output, double precision C(N), the correlations.
c
      implicit none

      integer n

      double precision c(n)
      integer i
      double precision rho(n)
      double precision rho0
      double precision rhohat

      do i = 1, n

        rhohat = min ( abs ( rho(i) ) / rho0, 1.0D+00 )

        c(i) = 1.0D+00 - 1.875D+00 * rhohat + 1.25D+00 * rhohat ** 3 
     &    - 0.375D+00 * rhohat ** 5

      end do

      return
      end
      subroutine correlation_plot ( n, rho, c, header, title )

c*********************************************************************72
c
cc CORRELATION_PLOT makes a plot of a correlation function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of arguments.
c
c    Input, double precision RHO(N), the arguments.
c
c    Input, double precision C(N), the correlations.
c
c    Input, character * ( * ) HEADER, an identifier for the files.
c
c    Input, character * ( * ) TITLE, a title for the plot.
c
      implicit none

      integer n

      double precision c(n)
      character * ( 80 ) command_filename
      integer command_unit
      character * ( 80 ) data_filename
      integer data_unit
      character * ( * ) header
      integer i
      double precision rho(n)
      double precision rho0
      character * ( * ) title

      call get_unit ( data_unit )
      data_filename = trim ( header ) // '_data.txt'
      open ( unit = data_unit, file = data_filename, 
     &  status = 'replace' )
      do i = 1, n
        write ( data_unit, '(2x,g14.6,2x,g14.6)' ) rho(i), c(i)
      end do
      close ( unit = data_unit )
      write ( *, '(a)' ) '  Created data file "' 
     &  // trim ( data_filename ) // '".'

      call get_unit ( command_unit )
      command_filename = trim ( header ) // '_commands.txt'
      open ( unit = command_unit, file = command_filename, 
     &  status = 'replace' )
      write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) '# Usage:'
      write ( command_unit, '(a)' ) '#  gnuplot < ' 
     &  // trim ( command_filename )
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) 'set term png'
      write ( command_unit, '(a)' ) 
     &  'set output "' // trim ( header ) // '_plot.png"'
      write ( command_unit, '(a)' ) 'set xlabel "Distance Rho"'
      write ( command_unit, '(a)' ) 'set ylabel "Correlation C(Rho)"'
      write ( command_unit, '(a)' ) 'set title "' 
     &  // trim ( title ) // '"'
      write ( command_unit, '(a)' ) 'set grid'
      write ( command_unit, '(a)' ) 'set style data lines'
      write ( command_unit, '(a)' ) 'plot "' // trim ( data_filename ) 
     &  // '" using 1:2 lw 3 linecolor rgb "blue"'
      write ( command_unit, '(a)' ) 'quit'
      close ( unit = command_unit )
      write ( *, '(a)' ) 
     &  '  Created command file "' // trim ( command_filename ) // '".'

      return
      end
      subroutine correlation_plots ( n, n2, rho, rho0, c, header, 
     &  title )

c*********************************************************************72
c
cc CORRELATION_PLOTS plots correlations for a range of correlation lengths.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of values of RHO.
c
c    Input, integer N2, the number of values of RHO0.
c
c    Input, double precision RHO(N), the independent value.
c
c    Input, double precision RHO0(N2), the correlation lengths.
c
c    Input, double precision C(N,N2), the correlations.
c
c    Input, character * ( * ) HEADER, an identifier for the files.
c
c    Input, character * ( * ) TITLE, a title for the plot.
c
      implicit none

      integer n
      integer n2

      double precision c(n,n2)
      character * ( 80 ) command_filename
      integer command_unit
      character * ( 80 ) data_filename
      integer data_unit
      character * ( 40 ) format_string
      character * ( * ) header
      integer i
      double precision rho(n)
      double precision rho0(n2)
      character * ( * ) title

      write ( format_string, '(a1,i8,a1,i8,a1,i8,a1)' ) 
     &  '(', n2 + 1, 'g', 14, '.', 6, ')'

      call get_unit ( data_unit )
      data_filename = trim ( header ) // '_plots_data.txt'
      open ( unit = data_unit, file = data_filename, 
     &  status = 'replace' )
      do i = 1, n
        write ( data_unit, format_string ) rho(i), c(i,1:n2)
      end do
      close ( unit = data_unit )
      write ( *, '(a)' ) '  Created data file "' 
     &  // trim ( data_filename ) // '".'

      call get_unit ( command_unit )
      command_filename = trim ( header ) // '_plots_commands.txt'
      open ( unit = command_unit, file = command_filename, 
     &  status = 'replace' )
      write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) '# Usage:'
      write ( command_unit, '(a)' ) '#  gnuplot < ' 
     &  // trim ( command_filename )
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) 'set term png'
      write ( command_unit, '(a)' ) 
     &  'set output "' // trim ( header ) // '_plots.png"'
      write ( command_unit, '(a)' ) 'set xlabel "Rho"'
      write ( command_unit, '(a)' ) 'set ylabel "Correlation(Rho)"'
      write ( command_unit, '(a)' ) 'set title "' 
     &  // trim ( title ) // '"'
      write ( command_unit, '(a)' ) 'set grid'
      write ( command_unit, '(a)' ) 'set style data lines'
      write ( command_unit, '(a)' ) 'set key off'
      if ( n2 .eq. 1 ) then
        write ( command_unit, '(a)' ) 'plot "' 
     &  // trim ( data_filename ) // '" using 1:2 lw 3'
      else
        write ( command_unit, '(a)' ) 'plot "' 
     &    // trim ( data_filename ) // '" using 1:2 lw 3, \'
        do i = 2, n2 - 1
          write ( command_unit, '(a,i3,a)' ) '     "' 
     &      // trim ( data_filename ) // '" using 1:', i + 1, ' lw 3, \'
        end do
        write ( command_unit, '(a,i3,a)' ) '     "' 
     &    // trim ( data_filename ) // '" using 1:', n2 + 1, ' lw 3'
      end if
      write ( command_unit, '(a)' ) 'quit'
      close ( unit = command_unit )
      write ( *, '(a)' ) 
     &  '  Created command file "' // trim ( command_filename ) // '".'

      return
      end
      subroutine correlation_power ( n, rho, rho0, c )

c*********************************************************************72
c
cc CORRELATION_POWER evaluates the power correlation function.
c
c  Discussion:
c
c    In order to be able to call this routine under a dummy name, I had
c    to drop E from the argument list.
c
c    The power correlation is
c
c      C(rho) = ( 1 - |rho| )^e  if 0 <= |rho| <= 1
c             = 0                otherwise
c
c      The constraint on the exponent is 2 <= e.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of arguments.
c
c    Input, double precision RHO(N), the arguments.
c    0.0 <= RHO.
c
c    Input, double precision RHO0, the correlation length.
c    0.0 < RHO0.
c
c    Input, double precision E, the exponent.
c    E has a default value of 2.0;
c    2.0 <= E.
c
c    Output, double precision C(N), the correlations.
c
      implicit none

      integer n

      double precision c(n)
      double precision e
      integer i
      double precision rho(n)
      double precision rho0
      double precision rhohat

      e = 2.0D+00

      do i = 1, n
        rhohat = abs ( rho(i) ) / rho0
        if ( rhohat .le. 1.0D+00 ) then
          c(i) = ( 1.0D+00 - rhohat ) ** e
        else
          c(i) = 0.0D+00
        end if
      end do

      return
      end
      subroutine correlation_rational_quadratic ( n, rho, rho0, c )

c*********************************************************************72
c
cc CORRELATION_RATIONAL_QUADRATIC: rational quadratic correlation function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Petter Abrahamsen,
c    A Review of Gaussian Random Fields and Correlation Functions,
c    Norwegian Computing Center, 1997.
c
c  Parameters:
c
c    Input, integer N, the number of arguments.
c
c    Input, double precision RHO(N), the arguments.
c
c    Input, double precision RHO0, the correlation length.
c
c    Output, double precision C(N), the correlations.
c
      implicit none

      integer n

      double precision c(n)
      integer i
      double precision rho(n)
      double precision rho0

      do i = 1, n
        c(i) = 1.0D+00 / ( 1.0D+00 + ( rho(i) / rho0 ) ** 2 )
      end do

      return
      end
      subroutine correlation_spherical ( n, rho, rho0, c )

c*********************************************************************72
c
cc CORRELATION_SPHERICAL evaluates the spherical correlation function.
c
c  Discussion:
c
c    This correlation is based on the volume of overlap of two spheres
c    of radius RHO0 and separation RHO.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Petter Abrahamsen,
c    A Review of Gaussian Random Fields and Correlation Functions,
c    Norwegian Computing Center, 1997.
c
c  Parameters:
c
c    Input, integer N, the number of arguments.
c
c    Input, double precision RHO(N), the arguments.
c
c    Input, double precision RHO0, the correlation length.
c
c    Output, double precision C(N), the correlations.
c
      implicit none

      integer n

      double precision c(n)
      integer i
      double precision rho(n)
      double precision rho0
      double precision rhohat

      do i = 1, n
        rhohat = min ( abs ( rho(i) ) / rho0, 1.0D+00 )
        c(i) = 1.0D+00 - 1.5D+00 * rhohat + 0.5D+00 * rhohat ** 3
      end do

      return
      end
      subroutine correlation_to_covariance ( n, c, sigma, k )

c*********************************************************************72
c
cc CORRELATION_TO_COVARIANCE: covariance matrix from a correlation matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision C(N,N), the correlation matrix.
c
c    Input, double precision SIGMA(N), the standard deviations.
c
c    Output, double precision K(N,N), the covariance matrix.
c
      implicit none

      integer n

      double precision c(n,n)
      double precision c_max
      double precision c_min
      double precision e
      double precision error_frobenius
      integer i
      integer j
      double precision k(n,n)
      double precision r8_epsilon
      double precision r8mat_max
      double precision r8mat_min
      double precision sigma(n)
      double precision tol

      tol = sqrt ( r8_epsilon ( ) )
c
c  C must be symmetric.
c
      call r8mat_is_symmetric ( n, n, c, error_frobenius )

      if ( tol .lt. error_frobenius ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'CORRELATION_TO_COVARIANCE - Fatal error!'
        write ( *, '(a,g14.6)' ) 
     &    '  Input matrix C fails symmetry test with error ', 
     &    error_frobenius
        stop
      end if
c
c  The diagonal must be 1.
c
      do i = 1, n
        e = abs ( c(i,i) - 1.0D+00 )
        if ( tol .lt. e ) then
          write ( *, '(a)' ) ''
          write ( *, '(a)' ) 'CORRELATION_TO_COVARIANCE - Fatal error!'
          write ( *, '(a)' ) 
     &      '  Input matrix C has non-unit diagonal entries.'
          write ( *, '(a,i4,a,g14.6)' ) '  Error on row ', i, ' is ', e
          stop
        end if
      end do
c
c  Off-diagonals must be between -1 and 1.
c
      c_min = r8mat_min ( n, n, c )

      if ( c_min .lt. - 1.0D+00 - tol ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'CORRELATION_TO_COVARIANCE - Fatal error!'
        write ( *, '(a)' ) '  Input matrix C has entries less than -1.0'
        stop
      end if

      c_max = r8mat_max ( n, n, c )

      if ( 1.0D+00 + tol .lt. c_max ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'CORRELATION_TO_COVARIANCE - Fatal error!'
        write ( *, '(a)' ) 
     &    '  Input matrix C has entries greater than +1.0'
        stop
      end if
c
c  Form K.
c
      do j = 1, n
        do i = 1, n
          k(i,j) = sigma(i) * c(i,j) * sigma(j)
        end do
      end do

      return
      end
      subroutine correlation_white_noise ( n, rho, rho0, c )

c*********************************************************************72
c
cc CORRELATION_WHITE_NOISE evaluates the white noise correlation function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Petter Abrahamsen,
c    A Review of Gaussian Random Fields and Correlation Functions,
c    Norwegian Computing Center, 1997.
c
c  Parameters:
c
c    Input, integer N, the number of arguments.
c
c    Input, double precision RHO(N), the arguments.
c
c    Input, double precision RHO0, the correlation length.
c
c    Output, double precision C(N), the correlations.
c
      implicit none

      integer n

      double precision c(n)
      integer i
      double precision rho(n)
      double precision rho0

      do i = 1, n
        if ( rho(i) .eq. 0.0D+00 ) then
          c(i) = 1.0D+00
        else
          c(i) = 0.0D+00
        end if
      end do

      return
      end
      subroutine covariance_to_correlation ( n, k, c, sigma )

c*********************************************************************72
c
cc COVARIANCE_TO_CORRELATION: correlation matrix from a covariance matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c
c    Input, double precision K(N,N), the covariance matrix.
c
c    Output, double precision C(N,N), the correlation matrix.
c
c    Output, double precision SIGMA(N), the standard deviations.
c
      implicit none

      integer n

      double precision c(n,n)
      double precision e
      double precision error_frobenius
      integer i
      integer j
      double precision k(n,n)
      double precision r8_epsilon
      double precision r8vec_min
      double precision sigma(n)
      double precision sigma_min
      double precision tol

      tol = sqrt ( r8_epsilon ( ) )
c
c  K must be symmetric.
c
      call r8mat_is_symmetric ( n, n, k, error_frobenius )

      if ( tol .lt. error_frobenius ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'COVARIANCE_TO_CORRELATION - Fatal error!'
        write ( *, '(a,g14.6)' ) 
     &    '  Input matrix K fails symmetry test with error ', 
     &    error_frobenius
        stop
      end if
c
c  It must be the case that K(I,J)^2 <= K(I,I) * K(J,J).
c
      e = 0.0D+00
      do i = 1, n
        do j = i + 1, n
          e = max ( e, k(i,j) ** 2 - k(i,i) * k(j,j) )
        end do
      end do

      if ( tol .lt. e ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'COVARIANCE_TO_CORRELATION - Fatal error!'
        write ( *, '(a)' ) 
     &    '  Input matrix K fails K(I,J)^2 <= K(I,I)*K(J,J)'
        stop
      end if
c
c  Get the diagonal.
c
      do i = 1, n
        sigma(i) = k(i,i)
      end do
c
c  Ensure the diagonal is positive.
c
      sigma_min = r8vec_min ( n, sigma )

      if ( sigma_min .le. 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'COVARIANCE_TO_CORRELATION - Fatal error!'
        write ( *, '(a,g14.6)' ) 
     &    '  Input matrix K has nonpositive diagonal entry = ', 
     &    sigma_min
        stop
      end if
c
c  Convert from variance to standard deviation.
c
      do i = 1, n
        sigma(i) = sqrt ( sigma(i) )
      end do
c
c  Form C.
c
      do j = 1, n
        do i = 1, n
          c(i,j) = k(i,j) / sigma(i) / sigma(j)
        end do
      end do

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
      function i4_modp ( i, j )

c*********************************************************************72
c
cc I4_MODP returns the nonnegative remainder of integer division.
c
c  Discussion:
c
c    If
c      NREM = I4_MODP ( I, J )
c      NMULT = ( I - NREM ) / J
c    then
c      I = J * NMULT + NREM
c    where NREM is always nonnegative.
c
c    The MOD function computes a result with the same sign as the
c    quantity being divided.  Thus, suppose you had an angle A,
c    and you wanted to ensure that it was between 0 and 360.
c    Then mod(A,360) would do, if A was positive, but if A
c    was negative, your result would be between -360 and 0.
c
c    On the other hand, I4_MODP(A,360) is between 0 and 360, always.
c
c  Example:
c
c        I     J     MOD I4_MODP    Factorization
c
c      107    50       7       7    107 =  2 *  50 + 7
c      107   -50       7       7    107 = -2 * -50 + 7
c     -107    50      -7      43   -107 = -3 *  50 + 43
c     -107   -50      -7      43   -107 =  3 * -50 + 43
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 December 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the number to be divided.
c
c    Input, integer J, the number that divides I.
c
c    Output, integer I4_MODP, the nonnegative remainder when I is
c    divided by J.
c
      implicit none

      integer i
      integer i4_modp
      integer j
      integer value

      if ( j .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_MODP - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal divisor J = ', j
        stop
      end if

      value = mod ( i, j )

      if ( value .lt. 0 ) then
        value = value + abs ( j )
      end if

      i4_modp = value

      return
      end
      function i4_wrap ( ival, ilo, ihi )

c*********************************************************************72
c
cc I4_WRAP forces an I4 to lie between given limits by wrapping.
c
c  Example:
c
c    ILO = 4, IHI = 8
c
c    I  Value
c
c    -2     8
c    -1     4
c     0     5
c     1     6
c     2     7
c     3     8
c     4     4
c     5     5
c     6     6
c     7     7
c     8     8
c     9     4
c    10     5
c    11     6
c    12     7
c    13     8
c    14     4
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 December 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer IVAL, an integer value.
c
c    Input, integer ILO, IHI, the desired bounds for the integer value.
c
c    Output, integer I4_WRAP, a "wrapped" version of IVAL.
c
      implicit none

      integer i4_modp
      integer i4_wrap
      integer ihi
      integer ilo
      integer ival
      integer jhi
      integer jlo
      integer value
      integer wide

      jlo = min ( ilo, ihi )
      jhi = max ( ilo, ihi )

      wide = jhi - jlo + 1

      if ( wide .eq. 1 ) then
        value = jlo
      else
        value = jlo + i4_modp ( ival - jlo, wide )
      end if

      i4_wrap = value

      return
      end
      subroutine minij ( m, n, a )

c*********************************************************************72
c
cc MINIJ returns the MINIJ matrix.
c
c  Formula:
c
c    A(I,J) = min ( I, J )
c
c  Example:
c
c    N = 5
c
c    1 1 1 1 1
c    1 2 2 2 2
c    1 2 3 3 3
c    1 2 3 4 4
c    1 2 3 4 5
c
c  Properties:
c
c    A is integral, therefore det ( A ) is integral, and 
c    det ( A ) * inverse ( A ) is integral.
c
c    A is positive definite.
c
c    A is symmetric: A' = A.
c
c    Because A is symmetric, it is normal.
c
c    Because A is normal, it is diagonalizable.
c
c    The inverse of A is tridiagonal.
c
c    The eigenvalues of A are
c
c      LAMBDA(I) = 0.5 / ( 1 - cos ( ( 2 * I - 1 ) * pi / ( 2 * N + 1 ) ) ),
c
c    (N+1)*ONES(N) - A also has a tridiagonal inverse.
c
c    Gregory and Karney consider the matrix defined by
c
c      B(I,J) = N + 1 - MAX(I,J)
c
c    which is equal to the MINIJ matrix, but with the rows and
c    columns reversed.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Robert Gregory, David Karney,
c    Example 3.12, Example 4.14,
c    A Collection of Matrices for Testing Computational Algorithms,
c    Wiley, 1969, page 41, page 74, 
c    LC: QA263.G68.
c
c    Daniel Rutherford,
c    Some continuant determinants arising in physics and chemistry II,
c    Proceedings of the Royal Society Edinburgh,
c    Volume 63, A, 1952, pages 232-241.
c
c    John Todd,
c    Basic Numerical Mathematics, Vol. 2: Numerical Algebra,
c    Academic Press, 1977, page 158.
c
c    Joan Westlake,
c    A Handbook of Numerical Matrix Inversion and Solution of 
c    Linear Equations,
c    John Wiley, 1968,
c    ISBN13: 978-0471936756,
c    LC: QA263.W47.
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns 
c    of the matrix.
c
c    Output, double precision A(M,N), the matrix.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer i
      integer j

      do j = 1, n
        do i = 1, m
          a(i,j) = dble ( min ( i, j ) )
        end do
      end do

      return
      end
      subroutine paths_plot ( n, n2, rho, x, header, title )

c*********************************************************************72
c
cc PATHS_PLOT plots a sequence of paths or simulations.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of points in each path.
c
c    Input, integer N2, the number of paths.
c
c    Input, double precision RHO(N), the independent value.
c
c    Input, double precision X(N,N2), the path values.
c
c    Input, character * ( * ) HEADER, an identifier for the files.
c
c    Input, character * ( * ) TITLE, a title for the plot.
c
      implicit none

      integer n
      integer n2

      character * ( 80 ) command_filename
      integer command_unit
      character * ( 80 ) data_filename
      integer data_unit
      character * ( 40 ) format_string
      character * ( * ) header
      integer i
      double precision rho(n)
      double precision rho0
      character * ( * ) title
      double precision x(n,n2)

      write ( format_string, '(a1,i8,a1,i8,a1,i8,a1)' ) 
     &  '(', n2 + 1, 'g', 14, '.', 6, ')'

      call get_unit ( data_unit )
      data_filename = trim ( header ) // '_path_data.txt'
      open ( unit = data_unit, file = data_filename, 
     &  status = 'replace' )
      do i = 1, n
        write ( data_unit, format_string ) rho(i), x(i,1:n2)
      end do
      close ( unit = data_unit )
      write ( *, '(a)' ) '  Created data file "' 
     &  // trim ( data_filename ) // '".'

      call get_unit ( command_unit )
      command_filename = trim ( header ) // '_path_commands.txt'
      open ( unit = command_unit, file = command_filename, 
     &  status = 'replace' )
      write ( command_unit, '(a)' ) '# ' // trim ( command_filename )
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) '# Usage:'
      write ( command_unit, '(a)' ) '#  gnuplot < ' 
     &  // trim ( command_filename )
      write ( command_unit, '(a)' ) '#'
      write ( command_unit, '(a)' ) 'set term png'
      write ( command_unit, '(a)' ) 
     &  'set output "' // trim ( header ) // '_paths.png"'
      write ( command_unit, '(a)' ) 'set xlabel "Rho"'
      write ( command_unit, '(a)' ) 'set ylabel "X(Rho)"'
      write ( command_unit, '(a)' ) 'set title "' 
     &  // trim ( title ) // '"'
      write ( command_unit, '(a)' ) 'set grid'
      write ( command_unit, '(a)' ) 'set style data lines'
      write ( command_unit, '(a)' ) 'set key off'
      if ( n2 .eq. 1 ) then
        write ( command_unit, '(a)' ) 'plot "' 
     &    // trim ( data_filename ) // '" using 1:2 lw 3'
      else
        write ( command_unit, '(a)' ) 'plot "' 
     &    // trim ( data_filename ) // '" using 1:2, \'
        do i = 2, n2 - 1
          write ( command_unit, '(a,i3,a)' ) '     "' 
     &      // trim ( data_filename ) // '" using 1:', i + 1, ', \'
        end do
        write ( command_unit, '(a,i3)' ) '     "' 
     &    // trim ( data_filename ) // '" using 1:', n2 + 1
      end if
      write ( command_unit, '(a)' ) 'quit'
      close ( unit = command_unit )
      write ( *, '(a)' ) 
     &  '  Created command file "' // trim ( command_filename ) // '".'

      return
      end
      function pythag ( a, b )

c*********************************************************************72
c
cc PYTHAG computes sqrt ( A * A + B * B ) carefully.
c
c  Discussion:
c
c    The formula
c
c      PYTHAG = sqrt ( A * A + B * B )
c
c    is reasonably accurate, but can fail if, for example, A^2 is larger
c    than the machine overflow.  The formula can lose most of its accuracy
c    if the sum of the squares is very large or very small.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    18 October 2009
c
c  Author:
c
c    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
c    Klema, Moler.
c    FORTRAN90 version by John Burkardt.
c
c  Reference:
c
c    James Wilkinson, Christian Reinsch,
c    Handbook for Automatic Computation,
c    Volume II, Linear Algebra, Part 2,
c    Springer, 1971,
c    ISBN: 0387054146,
c    LC: QA251.W67.
c
c    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
c    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
c    Matrix Eigensystem Routines, EISPACK Guide,
c    Lecture Notes in Computer Science, Volume 6,
c    Springer Verlag, 1976,
c    ISBN13: 978-3540075462,
c    LC: QA193.M37.
c
c  Modified:
c
c    04 February 2003
c
c  Parameters:
c
c    Input, double precision A, B, the two legs of a right triangle.
c
c    Output, double precision PYTHAG, the length of the hypotenuse.
c
      implicit none

      double precision a
      double precision b
      double precision p
      double precision pythag
      double precision r
      double precision s
      double precision t
      double precision u

      p = max ( abs ( a ), abs ( b ) )

      if ( p .ne. 0.0D+00 ) then

        r = ( min ( abs ( a ), abs ( b ) ) / p ) ** 2

10      continue

          t = 4.0D+00 + r

          if ( t .eq. 4.0D+00 ) then
            go to 20
          end if

          s = r / t
          u = 1.0D+00 + 2.0D+00 * s
          p = u * p
          r = ( s / u ) ** 2 * r

        go to 10

20      continue

      end if

      pythag = p

      return
      end
      subroutine r8_b0mp ( x, ampl, theta )

c*********************************************************************72
c
cc R8_B0MP evaluates the modulus and phase for the Bessel J0 and Y0 functions.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 September 2011
c
c  Author:
c
c    Original FORTRAN77 version by Wayne Fullerton.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Wayne Fullerton,
c    Portable Special Function Routines,
c    in Portability of Numerical Software,
c    edited by Wayne Cowell,
c    Lecture Notes in Computer Science, Volume 57,
c    Springer 1977,
c    ISBN: 978-3-540-08446-4,
c    LC: QA297.W65.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision AMPL, THETA, the modulus and phase.
c
      implicit none

      double precision ampl
      double precision bm0cs(37)
      double precision bm02cs(40)
      double precision bt02cs(39)
      double precision bth0cs(44)
      double precision eta
      integer nbm0
      integer nbm02
      integer nbt02
      integer nbth0
      double precision pi4
      double precision r8_csevl
      integer r8_inits
      double precision r8_mach
      double precision theta
      double precision x
      double precision xmax
      double precision z

      save bm0cs
      save bm02cs
      save bt02cs
      save bth0cs
      save nbm0
      save nbm02
      save nbt02
      save nbth0
      save pi4
      save xmax

      data bm0cs(  1) / +0.9211656246827742712573767730182D-01/
      data bm0cs(  2) / -0.1050590997271905102480716371755D-02/
      data bm0cs(  3) / +0.1470159840768759754056392850952D-04/
      data bm0cs(  4) / -0.5058557606038554223347929327702D-06/
      data bm0cs(  5) / +0.2787254538632444176630356137881D-07/
      data bm0cs(  6) / -0.2062363611780914802618841018973D-08/
      data bm0cs(  7) / +0.1870214313138879675138172596261D-09/
      data bm0cs(  8) / -0.1969330971135636200241730777825D-10/
      data bm0cs(  9) / +0.2325973793999275444012508818052D-11/
      data bm0cs( 10) / -0.3009520344938250272851224734482D-12/
      data bm0cs( 11) / +0.4194521333850669181471206768646D-13/
      data bm0cs( 12) / -0.6219449312188445825973267429564D-14/
      data bm0cs( 13) / +0.9718260411336068469601765885269D-15/
      data bm0cs( 14) / -0.1588478585701075207366635966937D-15/
      data bm0cs( 15) / +0.2700072193671308890086217324458D-16/
      data bm0cs( 16) / -0.4750092365234008992477504786773D-17/
      data bm0cs( 17) / +0.8615128162604370873191703746560D-18/
      data bm0cs( 18) / -0.1605608686956144815745602703359D-18/
      data bm0cs( 19) / +0.3066513987314482975188539801599D-19/
      data bm0cs( 20) / -0.5987764223193956430696505617066D-20/
      data bm0cs( 21) / +0.1192971253748248306489069841066D-20/
      data bm0cs( 22) / -0.2420969142044805489484682581333D-21/
      data bm0cs( 23) / +0.4996751760510616453371002879999D-22/
      data bm0cs( 24) / -0.1047493639351158510095040511999D-22/
      data bm0cs( 25) / +0.2227786843797468101048183466666D-23/
      data bm0cs( 26) / -0.4801813239398162862370542933333D-24/
      data bm0cs( 27) / +0.1047962723470959956476996266666D-24/
      data bm0cs( 28) / -0.2313858165678615325101260800000D-25/
      data bm0cs( 29) / +0.5164823088462674211635199999999D-26/
      data bm0cs( 30) / -0.1164691191850065389525401599999D-26/
      data bm0cs( 31) / +0.2651788486043319282958336000000D-27/
      data bm0cs( 32) / -0.6092559503825728497691306666666D-28/
      data bm0cs( 33) / +0.1411804686144259308038826666666D-28/
      data bm0cs( 34) / -0.3298094961231737245750613333333D-29/
      data bm0cs( 35) / +0.7763931143074065031714133333333D-30/
      data bm0cs( 36) / -0.1841031343661458478421333333333D-30/
      data bm0cs( 37) / +0.4395880138594310737100799999999D-31/

      data bth0cs(  1) / -0.24901780862128936717709793789967D+00/
      data bth0cs(  2) / +0.48550299609623749241048615535485D-03/
      data bth0cs(  3) / -0.54511837345017204950656273563505D-05/
      data bth0cs(  4) / +0.13558673059405964054377445929903D-06/
      data bth0cs(  5) / -0.55691398902227626227583218414920D-08/
      data bth0cs(  6) / +0.32609031824994335304004205719468D-09/
      data bth0cs(  7) / -0.24918807862461341125237903877993D-10/
      data bth0cs(  8) / +0.23449377420882520554352413564891D-11/
      data bth0cs(  9) / -0.26096534444310387762177574766136D-12/
      data bth0cs( 10) / +0.33353140420097395105869955014923D-13/
      data bth0cs( 11) / -0.47890000440572684646750770557409D-14/
      data bth0cs( 12) / +0.75956178436192215972642568545248D-15/
      data bth0cs( 13) / -0.13131556016891440382773397487633D-15/
      data bth0cs( 14) / +0.24483618345240857495426820738355D-16/
      data bth0cs( 15) / -0.48805729810618777683256761918331D-17/
      data bth0cs( 16) / +0.10327285029786316149223756361204D-17/
      data bth0cs( 17) / -0.23057633815057217157004744527025D-18/
      data bth0cs( 18) / +0.54044443001892693993017108483765D-19/
      data bth0cs( 19) / -0.13240695194366572724155032882385D-19/
      data bth0cs( 20) / +0.33780795621371970203424792124722D-20/
      data bth0cs( 21) / -0.89457629157111779003026926292299D-21/
      data bth0cs( 22) / +0.24519906889219317090899908651405D-21/
      data bth0cs( 23) / -0.69388422876866318680139933157657D-22/
      data bth0cs( 24) / +0.20228278714890138392946303337791D-22/
      data bth0cs( 25) / -0.60628500002335483105794195371764D-23/
      data bth0cs( 26) / +0.18649748964037635381823788396270D-23/
      data bth0cs( 27) / -0.58783732384849894560245036530867D-24/
      data bth0cs( 28) / +0.18958591447999563485531179503513D-24/
      data bth0cs( 29) / -0.62481979372258858959291620728565D-25/
      data bth0cs( 30) / +0.21017901684551024686638633529074D-25/
      data bth0cs( 31) / -0.72084300935209253690813933992446D-26/
      data bth0cs( 32) / +0.25181363892474240867156405976746D-26/
      data bth0cs( 33) / -0.89518042258785778806143945953643D-27/
      data bth0cs( 34) / +0.32357237479762298533256235868587D-27/
      data bth0cs( 35) / -0.11883010519855353657047144113796D-27/
      data bth0cs( 36) / +0.44306286907358104820579231941731D-28/
      data bth0cs( 37) / -0.16761009648834829495792010135681D-28/
      data bth0cs( 38) / +0.64292946921207466972532393966088D-29/
      data bth0cs( 39) / -0.24992261166978652421207213682763D-29/
      data bth0cs( 40) / +0.98399794299521955672828260355318D-30/
      data bth0cs( 41) / -0.39220375242408016397989131626158D-30/
      data bth0cs( 42) / +0.15818107030056522138590618845692D-30/
      data bth0cs( 43) / -0.64525506144890715944344098365426D-31/
      data bth0cs( 44) / +0.26611111369199356137177018346367D-31/

      data bm02cs(  1) / +0.9500415145228381369330861335560D-01/
      data bm02cs(  2) / -0.3801864682365670991748081566851D-03/
      data bm02cs(  3) / +0.2258339301031481192951829927224D-05/
      data bm02cs(  4) / -0.3895725802372228764730621412605D-07/
      data bm02cs(  5) / +0.1246886416512081697930990529725D-08/
      data bm02cs(  6) / -0.6065949022102503779803835058387D-10/
      data bm02cs(  7) / +0.4008461651421746991015275971045D-11/
      data bm02cs(  8) / -0.3350998183398094218467298794574D-12/
      data bm02cs(  9) / +0.3377119716517417367063264341996D-13/
      data bm02cs( 10) / -0.3964585901635012700569356295823D-14/
      data bm02cs( 11) / +0.5286111503883857217387939744735D-15/
      data bm02cs( 12) / -0.7852519083450852313654640243493D-16/
      data bm02cs( 13) / +0.1280300573386682201011634073449D-16/
      data bm02cs( 14) / -0.2263996296391429776287099244884D-17/
      data bm02cs( 15) / +0.4300496929656790388646410290477D-18/
      data bm02cs( 16) / -0.8705749805132587079747535451455D-19/
      data bm02cs( 17) / +0.1865862713962095141181442772050D-19/
      data bm02cs( 18) / -0.4210482486093065457345086972301D-20/
      data bm02cs( 19) / +0.9956676964228400991581627417842D-21/
      data bm02cs( 20) / -0.2457357442805313359605921478547D-21/
      data bm02cs( 21) / +0.6307692160762031568087353707059D-22/
      data bm02cs( 22) / -0.1678773691440740142693331172388D-22/
      data bm02cs( 23) / +0.4620259064673904433770878136087D-23/
      data bm02cs( 24) / -0.1311782266860308732237693402496D-23/
      data bm02cs( 25) / +0.3834087564116302827747922440276D-24/
      data bm02cs( 26) / -0.1151459324077741271072613293576D-24/
      data bm02cs( 27) / +0.3547210007523338523076971345213D-25/
      data bm02cs( 28) / -0.1119218385815004646264355942176D-25/
      data bm02cs( 29) / +0.3611879427629837831698404994257D-26/
      data bm02cs( 30) / -0.1190687765913333150092641762463D-26/
      data bm02cs( 31) / +0.4005094059403968131802476449536D-27/
      data bm02cs( 32) / -0.1373169422452212390595193916017D-27/
      data bm02cs( 33) / +0.4794199088742531585996491526437D-28/
      data bm02cs( 34) / -0.1702965627624109584006994476452D-28/
      data bm02cs( 35) / +0.6149512428936330071503575161324D-29/
      data bm02cs( 36) / -0.2255766896581828349944300237242D-29/
      data bm02cs( 37) / +0.8399707509294299486061658353200D-30/
      data bm02cs( 38) / -0.3172997595562602355567423936152D-30/
      data bm02cs( 39) / +0.1215205298881298554583333026514D-30/
      data bm02cs( 40) / -0.4715852749754438693013210568045D-31/

      data bt02cs(  1) / -0.24548295213424597462050467249324D+00/
      data bt02cs(  2) / +0.12544121039084615780785331778299D-02/
      data bt02cs(  3) / -0.31253950414871522854973446709571D-04/
      data bt02cs(  4) / +0.14709778249940831164453426969314D-05/
      data bt02cs(  5) / -0.99543488937950033643468850351158D-07/
      data bt02cs(  6) / +0.85493166733203041247578711397751D-08/
      data bt02cs(  7) / -0.86989759526554334557985512179192D-09/
      data bt02cs(  8) / +0.10052099533559791084540101082153D-09/
      data bt02cs(  9) / -0.12828230601708892903483623685544D-10/
      data bt02cs( 10) / +0.17731700781805131705655750451023D-11/
      data bt02cs( 11) / -0.26174574569485577488636284180925D-12/
      data bt02cs( 12) / +0.40828351389972059621966481221103D-13/
      data bt02cs( 13) / -0.66751668239742720054606749554261D-14/
      data bt02cs( 14) / +0.11365761393071629448392469549951D-14/
      data bt02cs( 15) / -0.20051189620647160250559266412117D-15/
      data bt02cs( 16) / +0.36497978794766269635720591464106D-16/
      data bt02cs( 17) / -0.68309637564582303169355843788800D-17/
      data bt02cs( 18) / +0.13107583145670756620057104267946D-17/
      data bt02cs( 19) / -0.25723363101850607778757130649599D-18/
      data bt02cs( 20) / +0.51521657441863959925267780949333D-19/
      data bt02cs( 21) / -0.10513017563758802637940741461333D-19/
      data bt02cs( 22) / +0.21820381991194813847301084501333D-20/
      data bt02cs( 23) / -0.46004701210362160577225905493333D-21/
      data bt02cs( 24) / +0.98407006925466818520953651199999D-22/
      data bt02cs( 25) / -0.21334038035728375844735986346666D-22/
      data bt02cs( 26) / +0.46831036423973365296066286933333D-23/
      data bt02cs( 27) / -0.10400213691985747236513382399999D-23/
      data bt02cs( 28) / +0.23349105677301510051777740800000D-24/
      data bt02cs( 29) / -0.52956825323318615788049749333333D-25/
      data bt02cs( 30) / +0.12126341952959756829196287999999D-25/
      data bt02cs( 31) / -0.28018897082289428760275626666666D-26/
      data bt02cs( 32) / +0.65292678987012873342593706666666D-27/
      data bt02cs( 33) / -0.15337980061873346427835733333333D-27/
      data bt02cs( 34) / +0.36305884306364536682359466666666D-28/
      data bt02cs( 35) / -0.86560755713629122479172266666666D-29/
      data bt02cs( 36) / +0.20779909972536284571238399999999D-29/
      data bt02cs( 37) / -0.50211170221417221674325333333333D-30/
      data bt02cs( 38) / +0.12208360279441714184191999999999D-30/
      data bt02cs( 39) / -0.29860056267039913454250666666666D-31/

      data nbm0 / 0 /
      data nbm02 / 0 /
      data nbt02 / 0 /
      data nbth0 / 0 /
      data pi4 / 0.785398163397448309615660845819876D+00 /
      data xmax / 0.0D+00 /

      if ( nbm0 .eq. 0 ) then
        eta = 0.1D+00 * r8_mach ( 3 )
        nbm0 = r8_inits ( bm0cs, 37, eta )
        nbt02 = r8_inits ( bt02cs, 39, eta )
        nbm02 = r8_inits ( bm02cs, 40, eta )
        nbth0 = r8_inits ( bth0cs, 44, eta )
        xmax = 1.0D+00 / r8_mach ( 4 )
      end if

      if ( x .lt. 4.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_B0MP - Fatal error!'
        write ( *, '(a)' ) '  X < 4.'
        stop
      else if ( x .le. 8.0D+00 ) then
        z = ( 128.0D+00 / x / x - 5.0D+00 ) / 3.0D+00
        ampl = ( 0.75D+00 + r8_csevl ( z, bm0cs, nbm0 ) ) / dsqrt ( x )
        theta = x - pi4 + r8_csevl ( z, bt02cs, nbt02 ) / x
      else
        z = 128.0D+00 / x / x - 1.0D+00
        ampl = ( 0.75D+00 + r8_csevl ( z, bm02cs, nbm02) ) 
     &    / dsqrt ( x )
        theta = x - pi4 + r8_csevl ( z, bth0cs, nbth0 ) / x
      end if

      return
      end
      function r8_besi1 ( x )

c*********************************************************************72
c
cc R8_BESI1 evaluates the Bessel function I of order 1 of an R8 argument.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 March 2010
c
c  Author:
c
c    Original FORTRAN77 version by Wayne Fullerton.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Wayne Fullerton,
c    Portable Special Function Routines,
c    in Portability of Numerical Software,
c    edited by Wayne Cowell,
c    Lecture Notes in Computer Science, Volume 57,
c    Springer 1977,
c    ISBN: 978-3-540-08446-4,
c    LC: QA297.W65.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision R8_BESI1, the Bessel function I of order 1 of X.
c
      implicit none

      double precision bi1cs(17)
      integer nti1
      double precision r8_besi1
      double precision r8_besi1e
      double precision r8_csevl
      integer r8_inits
      double precision r8_mach
      double precision x
      double precision xmax
      double precision xmin
      double precision xsml
      double precision y

      save bi1cs
      save nti1
      save xmax
      save xmin
      save xsml

      data bi1cs(  1) / -0.19717132610998597316138503218149D-02 /
      data bi1cs(  2) / +0.40734887667546480608155393652014D+00 /
      data bi1cs(  3) / +0.34838994299959455866245037783787D-01 /
      data bi1cs(  4) / +0.15453945563001236038598401058489D-02 /
      data bi1cs(  5) / +0.41888521098377784129458832004120D-04 /
      data bi1cs(  6) / +0.76490267648362114741959703966069D-06 /
      data bi1cs(  7) / +0.10042493924741178689179808037238D-07 /
      data bi1cs(  8) / +0.99322077919238106481371298054863D-10 /
      data bi1cs(  9) / +0.76638017918447637275200171681349D-12 /
      data bi1cs( 10) / +0.47414189238167394980388091948160D-14 /
      data bi1cs( 11) / +0.24041144040745181799863172032000D-16 /
      data bi1cs( 12) / +0.10171505007093713649121100799999D-18 /
      data bi1cs( 13) / +0.36450935657866949458491733333333D-21 /
      data bi1cs( 14) / +0.11205749502562039344810666666666D-23 /
      data bi1cs( 15) / +0.29875441934468088832000000000000D-26 /
      data bi1cs( 16) / +0.69732310939194709333333333333333D-29 /
      data bi1cs( 17) / +0.14367948220620800000000000000000D-31 /

      data nti1 / 0 /
      data xmax / 0.0D+00 /
      data xmin / 0.0D+00 /
      data xsml / 0.0D+00 /

      if ( nti1 .eq. 0 ) then
        nti1 = r8_inits ( bi1cs, 17, 0.1D+00 * r8_mach ( 3 ) )
        xmin = 2.0D+00 * r8_mach ( 1 )
        xsml = dsqrt ( 8.0D+00 * r8_mach ( 3 ) )
        xmax = dlog ( r8_mach ( 2 ) )
      end if

      y = dabs ( x )

      if ( y .le. xmin ) then
        r8_besi1 = 0.0D+00
      else if ( y .le. xsml ) then
        r8_besi1 = 0.5D+00 * x
      else if ( y .le. 3.0D+00 ) then
        r8_besi1 = x * ( 0.875D+00 
     &    + r8_csevl ( y * y / 4.5D+00 - 1.0D+00, bi1cs, nti1 ) )
      else if ( y .le. xmax ) then
        r8_besi1 = dexp ( y ) * r8_besi1e ( x )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_BESI1 - Fatal error!'
        write ( *, '(a)' ) '  Result overflows.'
        stop
      end if

      return
      end
      function r8_besi1e ( x )

c*********************************************************************72
c
cc R8_BESI1E evaluates the exponentially scaled Bessel function I1(X).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 March 2010
c
c  Author:
c
c    Original FORTRAN77 version by Wayne Fullerton.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Wayne Fullerton,
c    Portable Special Function Routines,
c    in Portability of Numerical Software,
c    edited by Wayne Cowell,
c    Lecture Notes in Computer Science, Volume 57,
c    Springer 1977,
c    ISBN: 978-3-540-08446-4,
c    LC: QA297.W65.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision R8_BESI1E, the exponentially scaled Bessel 
c    function I1(X).
c
      implicit none

      double precision ai12cs(69)
      double precision ai1cs(46)
      double precision bi1cs(17)
      double precision eta
      integer ntai1
      integer ntai12
      integer nti1
      double precision r8_besi1e
      double precision r8_csevl
      integer r8_inits
      double precision r8_mach
      double precision x
      double precision xmin
      double precision xsml
      double precision y

      save ai12cs
      save ai1cs
      save bi1cs
      save ntai1
      save ntai12
      save nti1
      save xmin
      save xsml

      data bi1cs(  1) / -0.19717132610998597316138503218149D-02 /
      data bi1cs(  2) / +0.40734887667546480608155393652014D+00 /
      data bi1cs(  3) / +0.34838994299959455866245037783787D-01 /
      data bi1cs(  4) / +0.15453945563001236038598401058489D-02 /
      data bi1cs(  5) / +0.41888521098377784129458832004120D-04 /
      data bi1cs(  6) / +0.76490267648362114741959703966069D-06 /
      data bi1cs(  7) / +0.10042493924741178689179808037238D-07 /
      data bi1cs(  8) / +0.99322077919238106481371298054863D-10 /
      data bi1cs(  9) / +0.76638017918447637275200171681349D-12 /
      data bi1cs( 10) / +0.47414189238167394980388091948160D-14 /
      data bi1cs( 11) / +0.24041144040745181799863172032000D-16 /
      data bi1cs( 12) / +0.10171505007093713649121100799999D-18 /
      data bi1cs( 13) / +0.36450935657866949458491733333333D-21 /
      data bi1cs( 14) / +0.11205749502562039344810666666666D-23 /
      data bi1cs( 15) / +0.29875441934468088832000000000000D-26 /
      data bi1cs( 16) / +0.69732310939194709333333333333333D-29 /
      data bi1cs( 17) / +0.14367948220620800000000000000000D-31 /

      data ai1cs(  1) / -0.2846744181881478674100372468307D-01 /
      data ai1cs(  2) / -0.1922953231443220651044448774979D-01 /
      data ai1cs(  3) / -0.6115185857943788982256249917785D-03 /
      data ai1cs(  4) / -0.2069971253350227708882823777979D-04 /
      data ai1cs(  5) / +0.8585619145810725565536944673138D-05 /
      data ai1cs(  6) / +0.1049498246711590862517453997860D-05 /
      data ai1cs(  7) / -0.2918338918447902202093432326697D-06 /
      data ai1cs(  8) / -0.1559378146631739000160680969077D-07 /
      data ai1cs(  9) / +0.1318012367144944705525302873909D-07 /
      data ai1cs( 10) / -0.1448423418183078317639134467815D-08 /
      data ai1cs( 11) / -0.2908512243993142094825040993010D-09 /
      data ai1cs( 12) / +0.1266388917875382387311159690403D-09 /
      data ai1cs( 13) / -0.1664947772919220670624178398580D-10 /
      data ai1cs( 14) / -0.1666653644609432976095937154999D-11 /
      data ai1cs( 15) / +0.1242602414290768265232168472017D-11 /
      data ai1cs( 16) / -0.2731549379672432397251461428633D-12 /
      data ai1cs( 17) / +0.2023947881645803780700262688981D-13 /
      data ai1cs( 18) / +0.7307950018116883636198698126123D-14 /
      data ai1cs( 19) / -0.3332905634404674943813778617133D-14 /
      data ai1cs( 20) / +0.7175346558512953743542254665670D-15 /
      data ai1cs( 21) / -0.6982530324796256355850629223656D-16 /
      data ai1cs( 22) / -0.1299944201562760760060446080587D-16 /
      data ai1cs( 23) / +0.8120942864242798892054678342860D-17 /
      data ai1cs( 24) / -0.2194016207410736898156266643783D-17 /
      data ai1cs( 25) / +0.3630516170029654848279860932334D-18 /
      data ai1cs( 26) / -0.1695139772439104166306866790399D-19 /
      data ai1cs( 27) / -0.1288184829897907807116882538222D-19 /
      data ai1cs( 28) / +0.5694428604967052780109991073109D-20 /
      data ai1cs( 29) / -0.1459597009090480056545509900287D-20 /
      data ai1cs( 30) / +0.2514546010675717314084691334485D-21 /
      data ai1cs( 31) / -0.1844758883139124818160400029013D-22 /
      data ai1cs( 32) / -0.6339760596227948641928609791999D-23 /
      data ai1cs( 33) / +0.3461441102031011111108146626560D-23 /
      data ai1cs( 34) / -0.1017062335371393547596541023573D-23 /
      data ai1cs( 35) / +0.2149877147090431445962500778666D-24 /
      data ai1cs( 36) / -0.3045252425238676401746206173866D-25 /
      data ai1cs( 37) / +0.5238082144721285982177634986666D-27 /
      data ai1cs( 38) / +0.1443583107089382446416789503999D-26 /
      data ai1cs( 39) / -0.6121302074890042733200670719999D-27 /
      data ai1cs( 40) / +0.1700011117467818418349189802666D-27 /
      data ai1cs( 41) / -0.3596589107984244158535215786666D-28 /
      data ai1cs( 42) / +0.5448178578948418576650513066666D-29 /
      data ai1cs( 43) / -0.2731831789689084989162564266666D-30 /
      data ai1cs( 44) / -0.1858905021708600715771903999999D-30 /
      data ai1cs( 45) / +0.9212682974513933441127765333333D-31 /
      data ai1cs( 46) / -0.2813835155653561106370833066666D-31 /

      data ai12cs(  1) / +0.2857623501828012047449845948469D-01  /
      data ai12cs(  2) / -0.9761097491361468407765164457302D-02  /
      data ai12cs(  3) / -0.1105889387626237162912569212775D-03  /
      data ai12cs(  4) / -0.3882564808877690393456544776274D-05  /
      data ai12cs(  5) / -0.2512236237870208925294520022121D-06  /
      data ai12cs(  6) / -0.2631468846889519506837052365232D-07  /
      data ai12cs(  7) / -0.3835380385964237022045006787968D-08  /
      data ai12cs(  8) / -0.5589743462196583806868112522229D-09  /
      data ai12cs(  9) / -0.1897495812350541234498925033238D-10 /
      data ai12cs( 10) / +0.3252603583015488238555080679949D-10 /
      data ai12cs( 11) / +0.1412580743661378133163366332846D-10 /
      data ai12cs( 12) / +0.2035628544147089507224526136840D-11 /
      data ai12cs( 13) / -0.7198551776245908512092589890446D-12 /
      data ai12cs( 14) / -0.4083551111092197318228499639691D-12 /
      data ai12cs( 15) / -0.2101541842772664313019845727462D-13 /
      data ai12cs( 16) / +0.4272440016711951354297788336997D-13 /
      data ai12cs( 17) / +0.1042027698412880276417414499948D-13 /
      data ai12cs( 18) / -0.3814403072437007804767072535396D-14 /
      data ai12cs( 19) / -0.1880354775510782448512734533963D-14 /
      data ai12cs( 20) / +0.3308202310920928282731903352405D-15 /
      data ai12cs( 21) / +0.2962628997645950139068546542052D-15 /
      data ai12cs( 22) / -0.3209525921993423958778373532887D-16 /
      data ai12cs( 23) / -0.4650305368489358325571282818979D-16 /
      data ai12cs( 24) / +0.4414348323071707949946113759641D-17 /
      data ai12cs( 25) / +0.7517296310842104805425458080295D-17 /
      data ai12cs( 26) / -0.9314178867326883375684847845157D-18 /
      data ai12cs( 27) / -0.1242193275194890956116784488697D-17 /
      data ai12cs( 28) / +0.2414276719454848469005153902176D-18 /
      data ai12cs( 29) / +0.2026944384053285178971922860692D-18 /
      data ai12cs( 30) / -0.6394267188269097787043919886811D-19 /
      data ai12cs( 31) / -0.3049812452373095896084884503571D-19 /
      data ai12cs( 32) / +0.1612841851651480225134622307691D-19 /
      data ai12cs( 33) / +0.3560913964309925054510270904620D-20 /
      data ai12cs( 34) / -0.3752017947936439079666828003246D-20 /
      data ai12cs( 35) / -0.5787037427074799345951982310741D-22 /
      data ai12cs( 36) / +0.7759997511648161961982369632092D-21 /
      data ai12cs( 37) / -0.1452790897202233394064459874085D-21 /
      data ai12cs( 38) / -0.1318225286739036702121922753374D-21 /
      data ai12cs( 39) / +0.6116654862903070701879991331717D-22 /
      data ai12cs( 40) / +0.1376279762427126427730243383634D-22 /
      data ai12cs( 41) / -0.1690837689959347884919839382306D-22 /
      data ai12cs( 42) / +0.1430596088595433153987201085385D-23 /
      data ai12cs( 43) / +0.3409557828090594020405367729902D-23 /
      data ai12cs( 44) / -0.1309457666270760227845738726424D-23 /
      data ai12cs( 45) / -0.3940706411240257436093521417557D-24 /
      data ai12cs( 46) / +0.4277137426980876580806166797352D-24 /
      data ai12cs( 47) / -0.4424634830982606881900283123029D-25 /
      data ai12cs( 48) / -0.8734113196230714972115309788747D-25 /
      data ai12cs( 49) / +0.4045401335683533392143404142428D-25 /
      data ai12cs( 50) / +0.7067100658094689465651607717806D-26 /
      data ai12cs( 51) / -0.1249463344565105223002864518605D-25 /
      data ai12cs( 52) / +0.2867392244403437032979483391426D-26 /
      data ai12cs( 53) / +0.2044292892504292670281779574210D-26 /
      data ai12cs( 54) / -0.1518636633820462568371346802911D-26 /
      data ai12cs( 55) / +0.8110181098187575886132279107037D-28 /
      data ai12cs( 56) / +0.3580379354773586091127173703270D-27 /
      data ai12cs( 57) / -0.1692929018927902509593057175448D-27 /
      data ai12cs( 58) / -0.2222902499702427639067758527774D-28 /
      data ai12cs( 59) / +0.5424535127145969655048600401128D-28 /
      data ai12cs( 60) / -0.1787068401578018688764912993304D-28 /
      data ai12cs( 61) / -0.6565479068722814938823929437880D-29 /
      data ai12cs( 62) / +0.7807013165061145280922067706839D-29 /
      data ai12cs( 63) / -0.1816595260668979717379333152221D-29 /
      data ai12cs( 64) / -0.1287704952660084820376875598959D-29 /
      data ai12cs( 65) / +0.1114548172988164547413709273694D-29 /
      data ai12cs( 66) / -0.1808343145039336939159368876687D-30 /
      data ai12cs( 67) / -0.2231677718203771952232448228939D-30 /
      data ai12cs( 68) / +0.1619029596080341510617909803614D-30 /
      data ai12cs( 69) / -0.1834079908804941413901308439210D-31 /

      data ntai1 / 0 /
      data ntai12 / 0 /
      data nti1 / 0 /
      data xmin / 0.0D+00 /
      data xsml / 0.0D+00 /

      if ( nti1 .eq. 0 ) then
        eta = 0.1D+00 * r8_mach ( 3 )
        nti1 = r8_inits ( bi1cs, 17, eta )
        ntai1 = r8_inits ( ai1cs, 46, eta )
        ntai12 = r8_inits ( ai12cs, 69, eta )
        xmin = 2.0D+00 * r8_mach ( 1 )
        xsml = dsqrt ( 8.0D+00 * r8_mach ( 3 ) )
      end if

      y = dabs ( x )

      if ( y .le. xmin ) then
        r8_besi1e = 0.0D+00
      else if ( y .le. xsml ) then
        r8_besi1e = 0.5D+00 * x * dexp ( - y )
      else if ( y .le. 3.0D+00 ) then
        r8_besi1e = x * ( 0.875D+00 
     &    + r8_csevl ( y * y / 4.5D+00 - 1.0D+00, bi1cs, nti1 ) )
     &    * dexp ( - y ) 
      else if ( y .le. 8.0D+00 ) then
        r8_besi1e = ( 0.375D+00 
     &    + r8_csevl ( ( 48.0D+00 / y - 11.0D+00) / 5.0D+00, 
     &    ai1cs, ntai1 ) ) / dsqrt ( y )
        if ( x .lt. 0.0D+00 ) then
          r8_besi1e = - r8_besi1e
        end if
      else
        r8_besi1e = ( 0.375D+00 
     &    + r8_csevl ( 16.0D+00 / y - 1.0D+00, ai12cs, ntai12 ) ) 
     &    / dsqrt ( y )
        if ( x .lt. 0.0D+00 ) then
          r8_besi1e = - r8_besi1e
        end if
      end if

      return
      end
      function r8_besj0 ( x )

c*********************************************************************72
c
cc R8_BESJ0 evaluates the Bessel function J of order 0 of an R8 argument.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 March 2010
c
c  Author:
c
c    Original FORTRAN77 version by Wayne Fullerton.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Wayne Fullerton,
c    Portable Special Function Routines,
c    in Portability of Numerical Software,
c    edited by Wayne Cowell,
c    Lecture Notes in Computer Science, Volume 57,
c    Springer 1977,
c    ISBN: 978-3-540-08446-4,
c    LC: QA297.W65.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision R8_BESJ0, the Bessel function J of order 0 of X.
c
      implicit none

      double precision ampl
      double precision bj0cs(19)
      integer ntj0
      double precision r8_besj0
      double precision r8_csevl
      integer r8_inits
      double precision r8_mach
      double precision theta
      double precision x
      double precision xsml
      double precision y

      save bj0cs
      save ntj0
      save xsml

      data bj0cs(  1) / +0.10025416196893913701073127264074D+00 /
      data bj0cs(  2) / -0.66522300776440513177678757831124D+00 /
      data bj0cs(  3) / +0.24898370349828131370460468726680D+00 /
      data bj0cs(  4) / -0.33252723170035769653884341503854D-01 /
      data bj0cs(  5) / +0.23114179304694015462904924117729D-02 /
      data bj0cs(  6) / -0.99112774199508092339048519336549D-04 /
      data bj0cs(  7) / +0.28916708643998808884733903747078D-05 /
      data bj0cs(  8) / -0.61210858663032635057818407481516D-07 /
      data bj0cs(  9) / +0.98386507938567841324768748636415D-09 /
      data bj0cs( 10) / -0.12423551597301765145515897006836D-10 /
      data bj0cs( 11) / +0.12654336302559045797915827210363D-12 /
      data bj0cs( 12) / -0.10619456495287244546914817512959D-14 /
      data bj0cs( 13) / +0.74706210758024567437098915584000D-17 /
      data bj0cs( 14) / -0.44697032274412780547627007999999D-19 /
      data bj0cs( 15) / +0.23024281584337436200523093333333D-21 /
      data bj0cs( 16) / -0.10319144794166698148522666666666D-23 /
      data bj0cs( 17) / +0.40608178274873322700800000000000D-26 /
      data bj0cs( 18) / -0.14143836005240913919999999999999D-28 /
      data bj0cs( 19) / +0.43910905496698880000000000000000D-31 /

      data ntj0 / 0 /
      data xsml / 0.0D+00 /

      if ( ntj0 .eq. 0 ) then
        ntj0 = r8_inits ( bj0cs, 19, 0.1D+00 * r8_mach ( 3 ) )
        xsml = dsqrt ( 4.0D+00 * r8_mach ( 3 ) )
      end if

      y = dabs ( x )

      if ( y .le. xsml ) then
        r8_besj0 = 1.0D+00
      else if ( y .le. 4.0D+00 ) then
        r8_besj0 = r8_csevl ( 0.125D+00 * y * y - 1.0D+00, bj0cs, ntj0 )
      else
        call r8_b0mp ( y, ampl, theta )
        r8_besj0 = ampl * dcos ( theta )
      end if

      return
      end
      function r8_besk ( nu, x )

c*********************************************************************72
c
cc R8_BESK evaluates the Bessel function K of order NU of an R8 argument.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 November 2012
c
c  Author:
c
c    John Burkardt.
c
c  Reference:
c
c    Wayne Fullerton,
c    Portable Special Function Routines,
c    in Portability of Numerical Software,
c    edited by Wayne Cowell,
c    Lecture Notes in Computer Science, Volume 57,
c    Springer 1977,
c    ISBN: 978-3-540-08446-4,
c    LC: QA297.W65.
c
c  Parameters:
c
c    Input, double precision NU, the order.
c
c    Input, double precision X, the argument.
c
c    Output, double precision R8_BESK, the Bessel function K of order NU at X.
c
      implicit none

      double precision bke(101)
      integer nin
      double precision nu
      double precision r8_besk
      double precision x
      double precision xnu

      xnu = nu - int ( nu )
      nin = int ( nu ) + 1

      call r8_besks ( xnu, x, nin, bke )

      r8_besk = bke(nin)

      return
      end
      function r8_besk1 ( x )

c*********************************************************************72
c
cc R8_BESK1 evaluates the Bessel function K of order 1 of an R8 argument.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 March 2010
c
c  Author:
c
c    Original FORTRAN77 version by Wayne Fullerton.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Wayne Fullerton,
c    Portable Special Function Routines,
c    in Portability of Numerical Software,
c    edited by Wayne Cowell,
c    Lecture Notes in Computer Science, Volume 57,
c    Springer 1977,
c    ISBN: 978-3-540-08446-4,
c    LC: QA297.W65.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision R8_BESK1, the Bessel function K of order 1 of X.
c
      implicit none

      double precision bk1cs(16)
      integer ntk1
      double precision r8_besi1
      double precision r8_besk1
      double precision r8_besk1e
      double precision r8_csevl
      integer r8_inits
      double precision r8_mach
      double precision x
      double precision xmax
      double precision xmin
      double precision xsml
      double precision y

      save bk1cs
      save ntk1
      save xmax
      save xmin
      save xsml

      data bk1cs(  1) / +0.25300227338947770532531120868533D-01 /
      data bk1cs(  2) / -0.35315596077654487566723831691801D+00 /
      data bk1cs(  3) / -0.12261118082265714823479067930042D+00 /
      data bk1cs(  4) / -0.69757238596398643501812920296083D-02 /
      data bk1cs(  5) / -0.17302889575130520630176507368979D-03 /
      data bk1cs(  6) / -0.24334061415659682349600735030164D-05 /
      data bk1cs(  7) / -0.22133876307347258558315252545126D-07 /
      data bk1cs(  8) / -0.14114883926335277610958330212608D-09 /
      data bk1cs(  9) / -0.66669016941993290060853751264373D-12 /
      data bk1cs( 10) / -0.24274498505193659339263196864853D-14 /
      data bk1cs( 11) / -0.70238634793862875971783797120000D-17 /
      data bk1cs( 12) / -0.16543275155100994675491029333333D-19 /
      data bk1cs( 13) / -0.32338347459944491991893333333333D-22 /
      data bk1cs( 14) / -0.53312750529265274999466666666666D-25 /
      data bk1cs( 15) / -0.75130407162157226666666666666666D-28 /
      data bk1cs( 16) / -0.91550857176541866666666666666666D-31 /

      data ntk1 / 0 /
      data xmax / 0.0D+00 /
      data xmin / 0.0D+00 /
      data xsml / 0.0D+00 /

      if ( ntk1 .eq. 0 ) then
        ntk1 = r8_inits ( bk1cs, 16, 0.1D+00 * r8_mach ( 3 ) )
        xmin = dexp ( dmax1 ( dlog ( r8_mach ( 1 ) ), 
     &    - dlog ( r8_mach ( 2 ) ) ) + 0.01D+00 )
        xsml = dsqrt ( 4.0D+00 * r8_mach ( 3 ) )
        xmax = - dlog ( r8_mach ( 1 ) )
        xmax = xmax - 0.5D+00 * xmax * dlog ( xmax ) 
     &    / ( xmax + 0.5D+00 ) - 0.01D+00
      end if

      if ( x .le. 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_BESK1 = Fatal error!'
        write ( *, '(a)' ) '  X <= 0.'
        stop
      else if ( x .le. xsml ) then
        y = 0.0D+00
        r8_besk1 = dlog ( 0.5D+00 * x ) * r8_besi1 ( x ) + ( 0.75D+00 
     &    + r8_csevl ( 0.5D+00 * y - 1.0D+00, bk1cs, ntk1 ) ) / x
      else if ( x .le. 2.0D+00 ) then
        y = x * x
        r8_besk1 = dlog ( 0.5D+00 * x ) * r8_besi1 ( x ) + ( 0.75D+00 
     &    + r8_csevl ( 0.5D+00 * y - 1.0D+00, bk1cs, ntk1 ) ) / x
      else if ( x .le. xmax ) then
        r8_besk1 = dexp ( - x ) * r8_besk1e ( x )
      else
        r8_besk1 = 0.0D+00
      end if

      return
      end
      function r8_besk1e ( x )

c*********************************************************************72
c
cc R8_BESK1E evaluates the exponentially scaled Bessel function K1(X).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 March 2010
c
c  Author:
c
c    Original FORTRAN77 version by Wayne Fullerton.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Wayne Fullerton,
c    Portable Special Function Routines,
c    in Portability of Numerical Software,
c    edited by Wayne Cowell,
c    Lecture Notes in Computer Science, Volume 57,
c    Springer 1977,
c    ISBN: 978-3-540-08446-4,
c    LC: QA297.W65.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision R8_BESK1E, the exponentially scaled Bessel 
c    function K1(X).
c
      implicit none

      double precision ak12cs(33)
      double precision ak1cs(38)
      double precision bk1cs(16)
      double precision eta
      integer ntak1
      integer ntak12
      integer ntk1
      double precision r8_besi1
      double precision r8_besk1e
      double precision r8_csevl
      integer r8_inits
      double precision r8_mach
      double precision x
      double precision xmin
      double precision xsml
      double precision y

      save ak12cs
      save ak1cs
      save bk1cs
      save ntak1
      save ntak12
      save ntk1
      save xmin
      save xsml

      data bk1cs(  1) / +0.25300227338947770532531120868533D-01 /
      data bk1cs(  2) / -0.35315596077654487566723831691801D+00 /
      data bk1cs(  3) / -0.12261118082265714823479067930042D+00 /
      data bk1cs(  4) / -0.69757238596398643501812920296083D-02 /
      data bk1cs(  5) / -0.17302889575130520630176507368979D-03 /
      data bk1cs(  6) / -0.24334061415659682349600735030164D-05 /
      data bk1cs(  7) / -0.22133876307347258558315252545126D-07 /
      data bk1cs(  8) / -0.14114883926335277610958330212608D-09 /
      data bk1cs(  9) / -0.66669016941993290060853751264373D-12 /
      data bk1cs( 10) / -0.24274498505193659339263196864853D-14 /
      data bk1cs( 11) / -0.70238634793862875971783797120000D-17 /
      data bk1cs( 12) / -0.16543275155100994675491029333333D-19 /
      data bk1cs( 13) / -0.32338347459944491991893333333333D-22 /
      data bk1cs( 14) / -0.53312750529265274999466666666666D-25 /
      data bk1cs( 15) / -0.75130407162157226666666666666666D-28 /
      data bk1cs( 16) / -0.91550857176541866666666666666666D-31 /

      data ak1cs(  1) / +0.27443134069738829695257666227266D+00 /
      data ak1cs(  2) / +0.75719899531993678170892378149290D-01 /
      data ak1cs(  3) / -0.14410515564754061229853116175625D-02 /
      data ak1cs(  4) / +0.66501169551257479394251385477036D-04 /
      data ak1cs(  5) / -0.43699847095201407660580845089167D-05 /
      data ak1cs(  6) / +0.35402774997630526799417139008534D-06 /
      data ak1cs(  7) / -0.33111637792932920208982688245704D-07 /
      data ak1cs(  8) / +0.34459775819010534532311499770992D-08 /
      data ak1cs(  9) / -0.38989323474754271048981937492758D-09 /
      data ak1cs( 10) / +0.47208197504658356400947449339005D-10 /
      data ak1cs( 11) / -0.60478356628753562345373591562890D-11 /
      data ak1cs( 12) / +0.81284948748658747888193837985663D-12 /
      data ak1cs( 13) / -0.11386945747147891428923915951042D-12 /
      data ak1cs( 14) / +0.16540358408462282325972948205090D-13 /
      data ak1cs( 15) / -0.24809025677068848221516010440533D-14 /
      data ak1cs( 16) / +0.38292378907024096948429227299157D-15 /
      data ak1cs( 17) / -0.60647341040012418187768210377386D-16 /
      data ak1cs( 18) / +0.98324256232648616038194004650666D-17 /
      data ak1cs( 19) / -0.16284168738284380035666620115626D-17 /
      data ak1cs( 20) / +0.27501536496752623718284120337066D-18 /
      data ak1cs( 21) / -0.47289666463953250924281069568000D-19 /
      data ak1cs( 22) / +0.82681500028109932722392050346666D-20 /
      data ak1cs( 23) / -0.14681405136624956337193964885333D-20 /
      data ak1cs( 24) / +0.26447639269208245978085894826666D-21 /
      data ak1cs( 25) / -0.48290157564856387897969868800000D-22 /
      data ak1cs( 26) / +0.89293020743610130180656332799999D-23 /
      data ak1cs( 27) / -0.16708397168972517176997751466666D-23 /
      data ak1cs( 28) / +0.31616456034040694931368618666666D-24 /
      data ak1cs( 29) / -0.60462055312274989106506410666666D-25 /
      data ak1cs( 30) / +0.11678798942042732700718421333333D-25 /
      data ak1cs( 31) / -0.22773741582653996232867840000000D-26 /
      data ak1cs( 32) / +0.44811097300773675795305813333333D-27 /
      data ak1cs( 33) / -0.88932884769020194062336000000000D-28 /
      data ak1cs( 34) / +0.17794680018850275131392000000000D-28 /
      data ak1cs( 35) / -0.35884555967329095821994666666666D-29 /
      data ak1cs( 36) / +0.72906290492694257991679999999999D-30 /
      data ak1cs( 37) / -0.14918449845546227073024000000000D-30 /
      data ak1cs( 38) / +0.30736573872934276300799999999999D-31 /

      data ak12cs(  1) / +0.6379308343739001036600488534102D-01 /
      data ak12cs(  2) / +0.2832887813049720935835030284708D-01 /
      data ak12cs(  3) / -0.2475370673905250345414545566732D-03 /
      data ak12cs(  4) / +0.5771972451607248820470976625763D-05 /
      data ak12cs(  5) / -0.2068939219536548302745533196552D-06 /
      data ak12cs(  6) / +0.9739983441381804180309213097887D-08 /
      data ak12cs(  7) / -0.5585336140380624984688895511129D-09 /
      data ak12cs(  8) / +0.3732996634046185240221212854731D-10 /
      data ak12cs(  9) / -0.2825051961023225445135065754928D-11 /
      data ak12cs( 10) / +0.2372019002484144173643496955486D-12 /
      data ak12cs( 11) / -0.2176677387991753979268301667938D-13 /
      data ak12cs( 12) / +0.2157914161616032453939562689706D-14 /
      data ak12cs( 13) / -0.2290196930718269275991551338154D-15 /
      data ak12cs( 14) / +0.2582885729823274961919939565226D-16 /
      data ak12cs( 15) / -0.3076752641268463187621098173440D-17 /
      data ak12cs( 16) / +0.3851487721280491597094896844799D-18 /
      data ak12cs( 17) / -0.5044794897641528977117282508800D-19 /
      data ak12cs( 18) / +0.6888673850418544237018292223999D-20 /
      data ak12cs( 19) / -0.9775041541950118303002132480000D-21 /
      data ak12cs( 20) / +0.1437416218523836461001659733333D-21 /
      data ak12cs( 21) / -0.2185059497344347373499733333333D-22 /
      data ak12cs( 22) / +0.3426245621809220631645388800000D-23 /
      data ak12cs( 23) / -0.5531064394246408232501248000000D-24 /
      data ak12cs( 24) / +0.9176601505685995403782826666666D-25 /
      data ak12cs( 25) / -0.1562287203618024911448746666666D-25 /
      data ak12cs( 26) / +0.2725419375484333132349439999999D-26 /
      data ak12cs( 27) / -0.4865674910074827992378026666666D-27 /
      data ak12cs( 28) / +0.8879388552723502587357866666666D-28 /
      data ak12cs( 29) / -0.1654585918039257548936533333333D-28 /
      data ak12cs( 30) / +0.3145111321357848674303999999999D-29 /
      data ak12cs( 31) / -0.6092998312193127612416000000000D-30 /
      data ak12cs( 32) / +0.1202021939369815834623999999999D-30 /
      data ak12cs( 33) / -0.2412930801459408841386666666666D-31 /

      data ntak1 / 0 /
      data ntak12 / 0 /
      data ntk1 / 0 /
      data xmin / 0.0D+00 /
      data xsml / 0.0D+00 /

      if ( ntk1 .eq. 0 ) then
        eta = 0.1D+00 * r8_mach ( 3 )
        ntk1 = r8_inits ( bk1cs, 16, eta )
        ntak1 = r8_inits ( ak1cs, 38, eta )
        ntak12 = r8_inits ( ak12cs, 33, eta )
        xmin = dexp ( dmax1 ( dlog ( r8_mach ( 1 ) ), 
     &    - dlog ( r8_mach ( 2 ) ) ) + 0.01D+00 )
        xsml = dsqrt ( 4.0D+00 * r8_mach ( 3 ) )
      end if

      if ( x .le. 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_BESK1E = Fatal error!'
        write ( *, '(a)' ) '  X <= 0.'
        stop
      else if ( x .le. xsml ) then
        y = 0.0D+00
        r8_besk1e = dexp ( x ) * ( dlog ( 0.5D+00 * x ) * r8_besi1 ( x )
     &    + ( 0.75D+00 
     &    + r8_csevl ( 0.5D+00 * y - 1.0D+00, bk1cs, ntk1 ) ) / x )
      else if ( x .le. 2.0D+00 ) then
        y = x * x
        r8_besk1e = dexp ( x ) * ( dlog ( 0.5D+00 * x ) * r8_besi1 ( x )
     &    + ( 0.75D+00 
     &    + r8_csevl ( 0.5D+00 * y - 1.0D+00, bk1cs, ntk1 ) ) / x )
      else if ( x .le. 8.0D+00 ) then
        r8_besk1e = ( 1.25D+00 
     &    + r8_csevl ( ( 16.0D+00 / x - 5.0D+00 ) / 3.0D+00, ak1cs, 
     &    ntak1 ) ) / dsqrt ( x )
      else
        r8_besk1e = ( 1.25D+00 +
     &    r8_csevl ( 16.0D+00 / x - 1.0D+00, ak12cs, ntak12 ) ) 
     &    / dsqrt ( x )
      end if

      return
      end
      subroutine r8_beskes ( xnu, x, nin, bke )

c*********************************************************************72
c
cc R8_BESKES evaluates a sequence of exponentially scaled K Bessel functions at X.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 April 2010
c
c  Author:
c
c    Original FORTRAN77 version by Wayne Fullerton.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Wayne Fullerton,
c    Portable Special Function Routines,
c    in Portability of Numerical Software,
c    edited by Wayne Cowell,
c    Lecture Notes in Computer Science, Volume 57,
c    Springer 1977,
c    ISBN: 978-3-540-08446-4,
c    LC: QA297.W65.
c
c  Parameters:
c
c    Input, double precision XNU, ?
c    |XNU| < 1.
c
c    Input, double precision X, the argument.
c
c    Input, integer NIN, indicates the number of terms to compute.
c
c    Output, double precision BKE(abs(NIN)), the exponentially scaled 
c    K Bessel functions.
c
      implicit none

      double precision bke(*)
      double precision bknu1
      double precision direct
      integer i
      integer iswtch
      integer n
      integer nin
      double precision r8_mach
      double precision v
      double precision vend
      double precision vincr
      double precision x
      double precision xnu

      v = dabs ( xnu )
      n = iabs ( nin )

      if ( 1.0D+00 .le. v ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_BESKES - Fatal error!'
        write ( *, '(a)' ) '  |XNU| must be less than 1.'
        stop
      end if

      if ( x .le. 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_BESKES - Fatal error!'
        write ( *, '(a)' ) '  X <= 0.'
        stop
      end if

      if ( n .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_BESKES - Fatal error!'
        write ( *, '(a)' ) '  N = 0.'
        stop
      end if

      call r8_knus ( v, x, bke(1), bknu1, iswtch )

      if ( n .eq. 1 ) then
        return
      end if

      if ( nin .lt. 0 ) then
        vincr = - 1.0D+00
      else
        vincr = + 1.0D+00
      end if

      if ( xnu .lt. 0.0D+00 ) then
        direct = - vincr
      else
        direct = vincr
      end if

      bke(2) = bknu1

      if ( direct .lt. 0.0D+00 ) then
        call r8_knus ( dabs ( xnu + vincr ), x, bke(2), bknu1, iswtch )
      end if

      if ( n .eq. 2 ) then
        return
      end if

      vend = dabs ( xnu + dble ( nin ) ) - 1.0D+00

      v = xnu
      do i = 3, n
        v = v + vincr
        bke(i) = 2.0D+00 * v * bke(i-1) / x + bke(i-2)
      end do

      return
      end
      subroutine r8_besks ( xnu, x, nin, bk )

c*********************************************************************72
c
cc R8_BESKS evaluates a sequence of K Bessel functions at X.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 April 2010
c
c  Author:
c
c    Original FORTRAN77 version by Wayne Fullerton.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Wayne Fullerton,
c    Portable Special Function Routines,
c    in Portability of Numerical Software,
c    edited by Wayne Cowell,
c    Lecture Notes in Computer Science, Volume 57,
c    Springer 1977,
c    ISBN: 978-3-540-08446-4,
c    LC: QA297.W65.
c
c  Parameters:
c
c    Input, double precision XNU, ?
c    |XNU| < 1.
c
c    Input, double precision X, the argument.
c
c    Input, integer NIN, indicates the number of terms to compute.
c
c    Output, double precision BK(abs(NIN)), the K Bessel functions.
c
      implicit none

      integer nin

      double precision bk(abs(nin))
      double precision expxi
      integer i
      integer n
      double precision r8_mach
      double precision x
      double precision xmax
      double precision xnu

      save xmax

      data xmax / 0.0D+00 /

      if ( xmax .eq. 0.0D+00 ) then
        xmax = - dlog ( r8_mach ( 1 ) )
        xmax = xmax + 0.5D+00 * dlog ( 3.14D+00 * 0.5D+00 / xmax )
      end if

      call r8_beskes ( xnu, x, nin, bk )

      expxi = dexp ( - x )
      n = iabs ( nin )

      do i = 1, n
        bk(i) = expxi * bk(i)
      end do

      return
      end
      function r8_csevl ( x, a, n )

c*********************************************************************72
c
cc R8_CSEVL evaluates a Chebyshev series.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 March 2010
c
c  Author:
c
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Roger Broucke,
c    Algorithm 446:
c    Ten Subroutines for the Manipulation of Chebyshev Series,
c    Communications of the ACM,
c    Volume 16, Number 4, April 1973, pages 254-256.
c
c  Parameters:
c
c    Input, double precision X, the evaluation point.
c
c    Input, double precision CS(N), the Chebyshev coefficients.
c
c    Input, integer N, the number of Chebyshev coefficients.
c
c    Output, double precision R8_CSEVL, the Chebyshev series evaluated at X.
c
      implicit none

      integer n

      double precision a(n)
      double precision b0
      double precision b1
      double precision b2
      integer i
      double precision r8_csevl
      double precision twox
      double precision x

      if ( n .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_CSEVL - Fatal error!'
        write ( *, '(a)' ) '  Number of terms <= 0.'
        stop
      end if

      if ( 1000 .lt. n ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_CSEVL - Fatal error!'
        write ( *, '(a)' ) '  Number of terms > 1000.'
        stop
      end if

      if ( x .lt. -1.1D+00 .or. 1.1D+00 .lt. x ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_CSEVL - Fatal error!'
        write ( *, '(a)' ) '  X outside (-1,+1)'
        write ( *, '(a,g14.6)' ) '  X = ', x
        stop
      end if

      twox = 2.0D+00 * x
      b1 = 0.0D+00
      b0 = 0.0D+00

      do i = n, 1, -1
        b2 = b1
        b1 = b0
        b0 = twox * b1 - b2 + a(i)
      end do

      r8_csevl = 0.5D+00 * ( b0 - b2 )

      return
      end
      function r8_epsilon ( )

c*********************************************************************72
c
cc R8_EPSILON returns the R8 roundoff unit.
c
c  Discussion:
c
c    The roundoff unit is a number R which is a power of 2 with the
c    property that, to the precision of the computer's arithmetic,
c      1 .lt. 1 + R
c    but
c      1 = ( 1 + R / 2 )
c
c    FORTRAN90 provides the superior library routine
c
c      EPSILON ( X )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision R8_EPSILON, the R8 roundoff unit.
c
      implicit none

      double precision r8_epsilon

      r8_epsilon = 2.220446049250313D-016

      return
      end
      subroutine r8_gaml ( xmin, xmax )

c*********************************************************************72
c
cc R8_GAML evaluates bounds for an R8 argument of the gamma function.
c
c  Discussion:
c
c    This function calculates the minimum and maximum legal bounds 
c    for X in the evaluation of GAMMA ( X ).
c
c    XMIN and XMAX are not the only bounds, but they are the only 
c    non-trivial ones to calculate.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 March 2010
c
c  Author:
c
c    Original FORTRAN77 version by Wayne Fullerton.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Wayne Fullerton,
c    Portable Special Function Routines,
c    in Portability of Numerical Software,
c    edited by Wayne Cowell,
c    Lecture Notes in Computer Science, Volume 57,
c    Springer 1977,
c    ISBN: 978-3-540-08446-4,
c    LC: QA297.W65.
c
c  Parameters:
c
c    Output, double precision XMIN, XMAX, the bounds.
c
      implicit none

      double precision alnbig
      double precision alnsml
      integer i
      integer j
      double precision r8_mach
      double precision xln
      double precision xmax
      double precision xmin
      double precision xold

      alnsml = dlog ( r8_mach ( 1 ) )
      xmin = - alnsml

      do i = 1, 10

        xold = xmin
        xln = dlog ( xmin )
        xmin = xmin - xmin * ( ( xmin + 0.5D+00 ) * xln - xmin 
     &    - 0.2258D+00 + alnsml ) / ( xmin * xln + 0.5D+00 )

        if ( dabs ( xmin - xold ) .lt. 0.005D+00 ) then

          xmin = - xmin + 0.01D+00

          alnbig = dlog ( r8_mach ( 2 ) )
          xmax = alnbig

          do j = 1, 10

            xold = xmax
            xln = dlog ( xmax )
            xmax = xmax - xmax * ( ( xmax - 0.5D+00 ) * xln - xmax 
     &        + 0.9189D+00 - alnbig ) / ( xmax * xln - 0.5D+00 )

            if ( dabs ( xmax - xold ) .lt. 0.005D+00 ) then
              xmax = xmax - 0.01D+00
              xmin = dmax1 ( xmin, - xmax + 1.0D+00 )
              return
            end if

          end do

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8_GAML - Fatal error!'
          write ( *, '(a)' ) '  Unable to find XMAX.'
          stop

        end if

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_GAML - Fatal error!'
      write ( *, '(a)' ) '  Unable to find XMIN.'

      stop
      end
      function r8_gamma ( x )

c*********************************************************************72
c
cc R8_GAMMA evaluates the gamma function of an R8 argument.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2010
c
c  Author:
c
c    Original FORTRAN77 version by Wayne Fullerton.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Wayne Fullerton,
c    Portable Special Function Routines,
c    in Portability of Numerical Software,
c    edited by Wayne Cowell,
c    Lecture Notes in Computer Science, Volume 57,
c    Springer 1977,
c    ISBN: 978-3-540-08446-4,
c    LC: QA297.W65.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision R8_GAMMA, the gamma function of X.
c
      implicit none

      double precision dxrel
      double precision gcs(42)
      integer i
      integer n
      integer ngcs
      double precision pi
      double precision r8_csevl
      double precision r8_gamma
      integer r8_inits
      double precision r8_lgmc
      double precision r8_mach
      double precision sinpiy
      double precision sq2pil
      double precision x
      double precision xmax
      double precision xmin
      double precision xsml
      double precision y

      save dxrel
      save gcs
      save ngcs
      save pi
      save sq2pil
      save xmax
      save xmin
      save xsml

      data gcs(  1) / +0.8571195590989331421920062399942D-02 /
      data gcs(  2) / +0.4415381324841006757191315771652D-02 /
      data gcs(  3) / +0.5685043681599363378632664588789D-01 /
      data gcs(  4) / -0.4219835396418560501012500186624D-02 /
      data gcs(  5) / +0.1326808181212460220584006796352D-02 /
      data gcs(  6) / -0.1893024529798880432523947023886D-03 /
      data gcs(  7) / +0.3606925327441245256578082217225D-04 /
      data gcs(  8) / -0.6056761904460864218485548290365D-05 /
      data gcs(  9) / +0.1055829546302283344731823509093D-05 /
      data gcs( 10) / -0.1811967365542384048291855891166D-06 /
      data gcs( 11) / +0.3117724964715322277790254593169D-07 /
      data gcs( 12) / -0.5354219639019687140874081024347D-08 /
      data gcs( 13) / +0.9193275519859588946887786825940D-09 /
      data gcs( 14) / -0.1577941280288339761767423273953D-09 /
      data gcs( 15) / +0.2707980622934954543266540433089D-10 /
      data gcs( 16) / -0.4646818653825730144081661058933D-11 /
      data gcs( 17) / +0.7973350192007419656460767175359D-12 /
      data gcs( 18) / -0.1368078209830916025799499172309D-12 /
      data gcs( 19) / +0.2347319486563800657233471771688D-13 /
      data gcs( 20) / -0.4027432614949066932766570534699D-14 /
      data gcs( 21) / +0.6910051747372100912138336975257D-15 /
      data gcs( 22) / -0.1185584500221992907052387126192D-15 /
      data gcs( 23) / +0.2034148542496373955201026051932D-16 /
      data gcs( 24) / -0.3490054341717405849274012949108D-17 /
      data gcs( 25) / +0.5987993856485305567135051066026D-18 /
      data gcs( 26) / -0.1027378057872228074490069778431D-18 /
      data gcs( 27) / +0.1762702816060529824942759660748D-19 /
      data gcs( 28) / -0.3024320653735306260958772112042D-20 /
      data gcs( 29) / +0.5188914660218397839717833550506D-21 /
      data gcs( 30) / -0.8902770842456576692449251601066D-22 /
      data gcs( 31) / +0.1527474068493342602274596891306D-22 /
      data gcs( 32) / -0.2620731256187362900257328332799D-23 /
      data gcs( 33) / +0.4496464047830538670331046570666D-24 /
      data gcs( 34) / -0.7714712731336877911703901525333D-25 /
      data gcs( 35) / +0.1323635453126044036486572714666D-25 /
      data gcs( 36) / -0.2270999412942928816702313813333D-26 /
      data gcs( 37) / +0.3896418998003991449320816639999D-27 /
      data gcs( 38) / -0.6685198115125953327792127999999D-28 /
      data gcs( 39) / +0.1146998663140024384347613866666D-28 /
      data gcs( 40) / -0.1967938586345134677295103999999D-29 /
      data gcs( 41) / +0.3376448816585338090334890666666D-30 /
      data gcs( 42) / -0.5793070335782135784625493333333D-31 /

      data dxrel / 0.0D+00 /
      data ngcs / 0 /
      data pi / 3.14159265358979323846264338327950D+00 /
      data sq2pil / 0.91893853320467274178032973640562D+00 /
      data xmax / 0.0D+00 /
      data xmin / 0.0D+00 /
      data xsml / 0.0D+00 /

      if ( ngcs .eq. 0 ) then
        ngcs = r8_inits ( gcs, 42, 0.1D+00 * r8_mach ( 3 ) )
        call r8_gaml ( xmin, xmax )
        xsml = dexp ( dmax1 ( dlog ( r8_mach ( 1 ) ),
     &    - dlog ( r8_mach ( 2 ) ) ) + 0.01D+00 )
        dxrel = dsqrt ( r8_mach ( 4 ) )
      end if

      y = dabs ( x )

      if ( y .le. 10.0D+00 ) then

        n = int ( x )
        if ( x .lt. 0.0D+00 ) then
          n = n - 1
        end if
        y = x - dble ( n )
        n = n - 1
        r8_gamma = 0.9375D+00 + r8_csevl ( 2.0D+00 * y - 1.0D+00, 
     &    gcs, ngcs )

        if ( n .eq. 0 ) then

          return

        else if ( n .lt. 0 ) then

          n = - n

          if ( x .eq. 0.0D+00 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'R8_GAMMA - Fatal error!'
            write ( *, '(a)' ) '  X is 0.'
            stop
          end if

          if ( x .lt. 0.0D+00 .and. 
     &      x + dble ( n - 2 ) .eq. 0.0D+00 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'R8_GAMMA - Fatal error!'
            write ( *, '(a)' ) '  X is a negative integer.'
            stop
          end if

          if ( x .lt. - 0.5D+00 .and. 
     &      dabs ( ( x - dint ( x - 0.5D+00 ) ) / x ) .lt. dxrel ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'R8_GAMMA - Warning!'
            write ( *, '(a)' ) '  X too near a negative integer,'
            write ( *, '(a)' ) '  answer is half precision.'
          end if

          if ( y .lt. xsml ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'R8_GAMMA - Fatal error!'
            write ( *, '(a)' ) 
     &        '  X is so close to zero that Gamma overflows.'
            stop
          end if

          do i = 1, n
            r8_gamma = r8_gamma / ( x + dble ( i - 1 ) )
          end do

        else if ( n .eq. 0 ) then

        else

          do i = 1, n
            r8_gamma = ( y + dble ( i ) ) * r8_gamma
          end do

        end if

      else

        if ( xmax .lt. x ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8_GAMMA - Fatal error!'
          write ( *, '(a)' ) '  X so big that Gamma overflows.'
          stop
        end if
c
c  Underflow.
c
        if ( x .lt. xmin ) then
          r8_gamma = 0.0D+00
          return
        end if

        r8_gamma = dexp ( ( y - 0.5D+00 ) * dlog ( y ) - y + sq2pil 
     &    + r8_lgmc ( y ) )

        if ( 0.0D+00 .lt. x ) then
          return
        end if

        if ( dabs ( ( x - dint ( x - 0.5D+00 ) ) / x ) .lt. dxrel ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8_GAMMA - Warning!'
          write ( *, '(a)' ) '  X too near a negative integer,'
          write ( *, '(a)' ) '  answer is half precision.'
        end if

        sinpiy = dsin ( pi * y )

        if ( sinpiy .eq. 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8_GAMMA - Fatal error!'
          write ( *, '(a)' ) '  X is a negative integer.'
          stop
        end if

        r8_gamma = - pi / ( y * sinpiy * r8_gamma )

      end if

      return
      end
      function r8_huge ( )

c*********************************************************************72
c
cc R8_HUGE returns a "huge" R8.
c
c  Discussion:
c
c    The value returned by this function is NOT required to be the
c    maximum representable R8.  This value varies from machine to machine,
c    from compiler to compiler, and may cause problems when being printed.
c    We simply want a "very large" but non-infinite number.
c
c    FORTRAN90 provides a built-in routine HUGE ( X ) that
c    can return the maximum representable number of the same datatype
c    as X, if that is what is really desired.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 April 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision R8_HUGE, a huge number.
c
      implicit none

      double precision r8_huge

      r8_huge = 1.0D+30

      return
      end
      function r8_inits ( dos, nos, eta )

c*********************************************************************72
c
cc R8_INITS initializes a Chebyshev series.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 March 2010
c
c  Author:
c
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Roger Broucke,
c    Algorithm 446:
c    Ten Subroutines for the Manipulation of Chebyshev Series,
c    Communications of the ACM,
c    Volume 16, Number 4, April 1973, pages 254-256.
c
c  Parameters:
c
c    Input, double precision DOS(NOS), the Chebyshev coefficients.
c
c    Input, integer NOS, the number of coefficients.
c
c    Input, double precision ETA, the desired accuracy.
c
c    Output, integer R8_INITS, the number of terms of the series needed
c    to ensure the requested accuracy.
c
      implicit none

      integer nos

      double precision dos(nos)
      double precision err 
      double precision eta
      integer i
      integer r8_inits

      if ( nos .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_INITS - Fatal error!'
        write ( *, '(a)' ) '  Number of coefficients < 1.'
        stop
      end if

      err = 0.0D+00

      do i = nos, 1, -1
        err = err + dabs ( dos(i) )
        if ( eta .lt. err ) then
          r8_inits = i
          return
        end if
      end do

      r8_inits = nos
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_INITS - Warning!'
      write ( *, '(a)' ) '  ETA may be too small.'

      return
      end
      subroutine r8_knus ( xnu, x, bknu, bknu1, iswtch )

c*********************************************************************72
c
cc R8_KNUS computes a sequence of K Bessel functions.
c
c  Discussion:
c
c    This routine computes Bessel functions 
c      exp(x) * k-sub-xnu (x)  
c    and
c      exp(x) * k-sub-xnu+1 (x) 
c    for 0.0 <= xnu < 1.0.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 April 2010
c
c  Author:
c
c    Original FORTRAN77 version by Wayne Fullerton.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Wayne Fullerton,
c    Portable Special Function Routines,
c    in Portability of Numerical Software,
c    edited by Wayne Cowell,
c    Lecture Notes in Computer Science, Volume 57,
c    Springer 1977,
c    ISBN: 978-3-540-08446-4,
c    LC: QA297.W65.
c
c  Parameters:
c
c    Input, double precision XNU, the order parameter.
c
c    Input, double precision X, the argument.
c
c    Output, double precision BKNU, BKNU1, the two K Bessel functions.
c
c    Output, integer ISWTCH, ?
c
      implicit none

      double precision a(32)
      double precision a0
      double precision aln2
      double precision alnbig
      double precision alneps
      double precision alnsml
      double precision alnz
      double precision alpha(32)
      double precision an
      double precision b0
      double precision beta(32)
      double precision bknu
      double precision bknu0
      double precision bknu1
      double precision bknud
      double precision bn
      double precision c0
      double precision c0kcs(29)
      double precision eta
      double precision euler
      double precision expx
      integer i
      integer ii
      integer inu
      integer iswtch
      integer n
      integer ntc0k
      integer nterms
      integer ntznu1
      double precision p1
      double precision p2
      double precision p3
      double precision qq
      double precision r8_csevl
      double precision r8_gamma
      integer r8_inits
      double precision r8_mach
      double precision result
      double precision sqpi2
      double precision sqrtx
      double precision v
      double precision vlnz
      double precision x
      double precision x2n
      double precision x2tov
      double precision xi
      double precision xmu
      double precision xnu
      double precision xnusml
      double precision xsml
      double precision z
      double precision znu1cs(20)
      double precision ztov

      save aln2
      save alnbig
      save alneps
      save alnsml
      save c0kcs
      save euler
      save ntc0k
      save ntznu1
      save sqpi2
      save xnusml
      save xsml
      save znu1cs

      data c0kcs(  1) / +0.60183057242626108387577445180329D-01     /
      data c0kcs(  2) / -0.15364871433017286092959755943124D+00     /
      data c0kcs(  3) / -0.11751176008210492040068229226213D-01     /
      data c0kcs(  4) / -0.85248788891979509827048401550987D-03     /
      data c0kcs(  5) / -0.61329838767496791874098176922111D-04     /
      data c0kcs(  6) / -0.44052281245510444562679889548505D-05     /
      data c0kcs(  7) / -0.31631246728384488192915445892199D-06     /
      data c0kcs(  8) / -0.22710719382899588330673771793396D-07     /
      data c0kcs(  9) / -0.16305644608077609552274620515360D-08     /
      data c0kcs( 10) / -0.11706939299414776568756044043130D-09     /
      data c0kcs( 11) / -0.84052063786464437174546593413792D-11    /
      data c0kcs( 12) / -0.60346670118979991487096050737198D-12    /
      data c0kcs( 13) / -0.43326960335681371952045997366903D-13    /
      data c0kcs( 14) / -0.31107358030203546214634697772237D-14    /
      data c0kcs( 15) / -0.22334078226736982254486133409840D-15    /
      data c0kcs( 16) / -0.16035146716864226300635791528610D-16    /
      data c0kcs( 17) / -0.11512717363666556196035697705305D-17    /
      data c0kcs( 18) / -0.82657591746836959105169479089258D-19    /
      data c0kcs( 19) / -0.59345480806383948172333436695984D-20    /
      data c0kcs( 20) / -0.42608138196467143926499613023976D-21    /
      data c0kcs( 21) / -0.30591266864812876299263698370542D-22    /
      data c0kcs( 22) / -0.21963541426734575224975501815516D-23    /
      data c0kcs( 23) / -0.15769113261495836071105750684760D-24    /
      data c0kcs( 24) / -0.11321713935950320948757731048056D-25    /
      data c0kcs( 25) / -0.81286248834598404082792349714433D-27    /
      data c0kcs( 26) / -0.58360900893453226552829349315949D-28    /
      data c0kcs( 27) / -0.41901241623610922519452337780905D-29    /
      data c0kcs( 28) / -0.30083737960206435069530504212862D-30    /
      data c0kcs( 29) / -0.21599152067808647728342168089832D-31    /

      data znu1cs(  1) / +0.203306756994191729674444001216911D+00    /
      data znu1cs(  2) / +0.140077933413219771062943670790563D+00    /
      data znu1cs(  3) / +0.791679696100161352840972241972320D-02    /
      data znu1cs(  4) / +0.339801182532104045352930092205750D-03    /
      data znu1cs(  5) / +0.117419756889893366664507228352690D-04    /
      data znu1cs(  6) / +0.339357570612261680333825865475121D-06    /
      data znu1cs(  7) / +0.842594176976219910194629891264803D-08    /
      data znu1cs(  8) / +0.183336677024850089184748150900090D-09    /
      data znu1cs(  9) / +0.354969844704416310863007064469557D-11   /
      data znu1cs( 10) / +0.619032496469887332205244342078407D-13   /
      data znu1cs( 11) / +0.981964535680439424960346115456527D-15   /
      data znu1cs( 12) / +0.142851314396490474211473563005985D-16   /
      data znu1cs( 13) / +0.191894921887825298966162467488436D-18   /
      data znu1cs( 14) / +0.239430979739498914162313140597128D-20   /
      data znu1cs( 15) / +0.278890246815347354835870465474995D-22   /
      data znu1cs( 16) / +0.304606650633033442582845214092865D-24   /
      data znu1cs( 17) / +0.313173237042191815771564260932089D-26   /
      data znu1cs( 18) / +0.304133098987854951645174908005034D-28   /
      data znu1cs( 19) / +0.279840384636833084343185097659733D-30   /
      data znu1cs( 20) / +0.244637186274497596485238794922666D-32   /

      data aln2 / 0.69314718055994530941723212145818D+00 /
      data alnbig / 0.0D+00 /
      data alneps / 0.0D+00 /
      data alnsml / 0.0D+00 /
      data euler / 0.57721566490153286060651209008240D+00 /
      data ntc0k / 0 /
      data ntznu1 / 0 /
      data sqpi2 / +1.2533141373155002512078826424055D+00 /
      data xnusml / 0.0D+00 /
      data xsml / 0.0D+00 /

      if ( ntc0k .eq. 0 ) then
        eta = 0.1D+00 * r8_mach ( 3 )
        ntc0k = r8_inits ( c0kcs, 29, eta )
        ntznu1 = r8_inits ( znu1cs, 20, eta )
        xnusml = dsqrt ( r8_mach ( 3 ) / 8.0D+00 )
        xsml = 0.1D+00 * r8_mach ( 3 )
        alnsml = dlog ( r8_mach ( 1 ) )
        alnbig = dlog ( r8_mach ( 2 ) )
        alneps = dlog ( 0.1D+00 * r8_mach ( 3 ) )
      end if

      if ( xnu .lt. 0.0D+00 .or. 1.0D+00 .le. xnu ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_KNUS - Fatal error!'
        write ( *, '(a)' ) '  XNU < 0 or. 1 <= XNU.'
        stop
      end if

      if ( x .le. 0.0D+00 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_KNUS - Fatal error!'
        write ( *, '(a)' ) '  X <= 0.'
        stop
      end if

      iswtch = 0
c
c  X is small.  Compute k-sub-xnu (x) and the derivative of k-sub-xnu (x)
c  then find k-sub-xnu+1 (x).  xnu is reduced to the interval (-0.5,+0.5)
c  then to (0., .5), because k of negative order (-nu) = k of positive
c  order (+nu).
c
      if ( x .le. 2.0D+00 ) then

        if ( xnu .le. 0.5D+00 ) then
          v = xnu
        else
          v = 1.0D+00 - xnu
        end if
c
c  carefully find (x/2)**xnu and z**xnu where z = x*x/4.
c
        alnz = 2.0D+00 * ( dlog ( x ) - aln2 )

        if ( x .le. xnu ) then

          if ( alnbig .lt. 
     &      - 0.5D+00 * xnu * alnz - aln2 - dlog ( xnu ) ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'R8_KNUS - Fatal error!'
            write ( *, '(a)' ) '  Small X causing overflow.'
            stop
          end if

        end if

        vlnz = v * alnz
        x2tov = dexp ( 0.5D+00 * vlnz )

        if ( vlnz .le. alnsml ) then
          ztov = 0.0D+00
        else
          ztov = x2tov * x2tov
        end if

        a0 = 0.5D+00 * r8_gamma ( 1.0D+00 + v )
        b0 = 0.5D+00 * r8_gamma ( 1.0D+00 - v )
        c0 = - euler
        if ( 0.5D+00 .le. ztov .and. xnusml .lt. v ) then
          c0 = - 0.75D+00 +
     &      r8_csevl ( ( 8.0D+00 * v ) * v - 1.0D+00, c0kcs, ntc0k )
        end if

        if ( ztov .le. 0.5D+00 ) then
          alpha(1) = ( a0 - ztov * b0 ) / v
        else
          alpha(1) = c0 - alnz * ( 0.75D+00 +
     &      r8_csevl ( vlnz / 0.35D+00 + 1.0D+00, znu1cs, ntznu1 ) ) 
     &      * b0
        end if

        beta(1) = - 0.5D+00 * ( a0 + ztov * b0 )

        if ( x .le. xsml ) then
          z = 0.0D+00
        else
          z = 0.25D+00 * x * x
        end if

        nterms = max ( 2, int ( 11.0D+00 
     &    + ( 8.0D+00 * alnz - 25.19D+00 - alneps ) 
     &    / ( 4.28D+00 - alnz ) ) )

        do i = 2, nterms
          xi = dble ( i - 1 )
          a0 = a0 / ( xi * ( xi - v ) )
          b0 = b0 / ( xi * ( xi + v ) )
          alpha(i) = ( alpha(i-1) + 2.0D+00 * xi * a0 ) 
     &      / ( xi * ( xi + v ) )
          beta(i) = ( xi - 0.5D+00 * v ) * alpha(i) - ztov * b0
        end do

        bknu = alpha(nterms)
        bknud = beta(nterms)
        do ii = 2, nterms
          i = nterms + 1 - ii
          bknu = alpha(i) + bknu * z
          bknud = beta(i) + bknud * z
        end do

        expx = dexp ( x )
        bknu = expx * bknu / x2tov

        if ( alnbig .lt. 
     &    - 0.5D+00 * ( xnu + 1.0D+00 ) * alnz - 2.0D+00 * aln2 ) then
          iswtch = 1
          return
        end if

        bknud = expx * bknud * 2.0D+00 / ( x2tov * x )

        if ( xnu .le. 0.5D+00 ) then
          bknu1 = v * bknu / x - bknud
          return
        end if

        bknu0 = bknu
        bknu = - v * bknu / x - bknud
        bknu1 = 2.0D+00 * xnu * bknu / x + bknu0
c
c  x is large.  find k-sub-xnu (x) and k-sub-xnu+1 (x) with y. l. luke-s
c  rational expansion.
c
      else

        sqrtx = dsqrt ( x )

        if ( 1.0D+00 / xsml .lt. x ) then
          bknu = sqpi2 / sqrtx
          bknu1 = bknu
          return
        end if

        an = - 0.60D+00 - 1.02D+00 / x
        bn = - 0.27D+00 - 0.53D+00 / x
        nterms = min ( 32, max ( 3, int ( an + bn * alneps ) ) )

        do inu = 1, 2

          if ( inu .eq. 1 ) then
            if ( xnu .le. xnusml ) then
              xmu = 0.0D+00
            else
              xmu = ( 4.0D+00 * xnu ) * xnu
            end if
          else
            xmu = 4.0D+00 * ( dabs ( xnu ) + 1.0D+00 )**2
          end if

          a(1) = 1.0D+00 - xmu
          a(2) = 9.0D+00 - xmu
          a(3) = 25.0D+00 - xmu

          if ( a(2) .eq. 0.0D+00 ) then

            result = sqpi2 * ( 16.0D+00 * x + xmu + 7.0D+00 ) 
     &        / ( 16.0D+00 * x * sqrtx )

          else

            alpha(1) = 1.0D+00
            alpha(2) = ( 16.0D+00 * x + a(2) ) / a(2)
            alpha(3) = ( ( 768.0D+00 * x + 48.0D+00 * a(3) ) * x 
     &        + a(2) * a(3) ) / ( a(2) * a(3) )

            beta(1) = 1.0D+00
            beta(2) = ( 16.0D+00 * x + ( xmu + 7.0D+00 ) ) / a(2)
            beta(3) = ( ( 768.0D+00 * x 
     &        + 48.0D+00 * ( xmu + 23.0D+00 ) ) * x +
     &        ( ( xmu + 62.0D+00 ) * xmu + 129.0D+00 ) ) 
     &        / ( a(2) * a(3) )

            do i = 4, nterms

              n = i - 1
              x2n = dble ( 2 * n - 1 )

              a(i) = ( x2n + 2.0D+00 )**2 - xmu
              qq = 16.0D+00 * x2n / a(i)
              p1 = - x2n * ( dble ( 12 * n * n - 20 * n ) - a(1) ) 
     &          / ( ( x2n - 2.0D+00 ) * a(i) ) - qq * x
              p2 = ( dble ( 12 * n * n - 28 * n + 8 ) - a(1) ) / a(i) 
     &          - qq * x
              p3 = - x2n * a(i-3) / ( ( x2n - 2.0D+00 ) * a(i))

              alpha(i) = - p1 * alpha(i-1) 
     &                   - p2 * alpha(i-2) 
     &                   - p3 * alpha(i-3)

              beta(i) =  - p1 * beta(i-1)  
     &                   - p2 * beta(i-2) 
     &                   - p3 * beta(i-3)

            end do

            result = sqpi2 * beta(nterms) / ( sqrtx * alpha(nterms) )

          end if

          if ( inu .eq. 1 ) then
            bknu = result
          else
            bknu1 = result
          end if

        end do

      end if

      return
      end
      function r8_lgmc ( x )

c*********************************************************************72
c
cc R8_LGMC evaluates the log gamma correction factor for an R8 argument.
c
c  Discussion:
c
c    For 10 <= X, compute the log gamma correction factor so that
c
c      log ( gamma ( x ) ) = log ( sqrt ( 2 * pi ) ) 
c                          + ( x - 0.5 ) * log ( x ) - x 
c                          + r8_lgmc ( x )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 March 2010
c
c  Author:
c
c    Original FORTRAN77 version by Wayne Fullerton.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Wayne Fullerton,
c    Portable Special Function Routines,
c    in Portability of Numerical Software,
c    edited by Wayne Cowell,
c    Lecture Notes in Computer Science, Volume 57,
c    Springer 1977,
c    ISBN: 978-3-540-08446-4,
c    LC: QA297.W65.
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision R8_LGMC, the correction factor.
c
      implicit none

      double precision algmcs(15)
      integer nalgm
      double precision r8_csevl
      integer r8_inits
      double precision r8_lgmc
      double precision r8_mach
      double precision x
      double precision xbig
      double precision xmax

      save algmcs
      save nalgm
      save xbig
      save xmax

      data algmcs(  1) / +0.1666389480451863247205729650822D+00 /
      data algmcs(  2) / -0.1384948176067563840732986059135D-04 /
      data algmcs(  3) / +0.9810825646924729426157171547487D-08 /
      data algmcs(  4) / -0.1809129475572494194263306266719D-10 /
      data algmcs(  5) / +0.6221098041892605227126015543416D-13 /
      data algmcs(  6) / -0.3399615005417721944303330599666D-15 /
      data algmcs(  7) / +0.2683181998482698748957538846666D-17 /
      data algmcs(  8) / -0.2868042435334643284144622399999D-19 /
      data algmcs(  9) / +0.3962837061046434803679306666666D-21 /
      data algmcs( 10) / -0.6831888753985766870111999999999D-23 /
      data algmcs( 11) / +0.1429227355942498147573333333333D-24 /
      data algmcs( 12) / -0.3547598158101070547199999999999D-26 /
      data algmcs( 13) / +0.1025680058010470912000000000000D-27 /
      data algmcs( 14) / -0.3401102254316748799999999999999D-29 /
      data algmcs( 15) / +0.1276642195630062933333333333333D-30 /

      data nalgm / 0 /
      data xbig / 0.0D+00 /
      data xmax / 0.0D+00 /

      if ( nalgm .eq. 0 ) then
        nalgm = r8_inits ( algmcs, 15, r8_mach ( 3 ) )
        xbig = 1.0D+00 / dsqrt ( r8_mach ( 3 ) )
        xmax = dexp ( dmin1 ( dlog ( r8_mach ( 2 ) / 12.0D+00 ), 
     &    - dlog ( 12.0D+00 * r8_mach ( 1 ) ) ) )
      end if

      if ( x .lt. 10.0D+00 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_LGMC - Fatal error!'
        write ( *, '(a)' ) '  X must be at least 10.'
        stop

      else if ( x .lt. xbig ) then

        r8_lgmc = r8_csevl ( 2.0D+00 * ( 10.0D+00 / x ) 
     &    * ( 10.0D+00 / x ) - 1.0D+00, algmcs, nalgm ) / x

      else if ( x .lt. xmax ) then

        r8_lgmc = 1.0D+00 / ( 12.0D+00 * x )

      else

        r8_lgmc = 0.0D+00

      end if

      return
      end
      function r8_mach ( i )

c*********************************************************************72
c
cc R8_MACH returns double precision real machine-dependent constants.
c
c  Discussion:
c
c    R8_MACH can be used to obtain machine-dependent parameters
c    for the local machine environment.  It is a function
c    with one input argument, and can be called as follows:
c
c      D = R8_MACH ( I )
c
c    where I=1,...,5.  The output value of D above is
c    determined by the input value of I:.
c
c    R8_MACH ( 1) = B**(EMIN-1), the smallest positive magnitude.
c    R8_MACH ( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
c    R8_MACH ( 3) = B**(-T), the smallest relative spacing.
c    R8_MACH ( 4) = B**(1-T), the largest relative spacing.
c    R8_MACH ( 5) = LOG10(B)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 April 2007
c
c  Author:
c
c    Original FORTRAN77 version by Phyllis Fox, Andrew Hall, Norman Schryer.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Phyllis Fox, Andrew Hall, Norman Schryer,
c    Algorithm 528:
c    Framework for a Portable Library,
c    ACM Transactions on Mathematical Software,
c    Volume 4, Number 2, June 1978, page 176-188.
c
c  Parameters:
c
c    Input, integer I, the index of the desired constant.
c
c    Output, double precision R8_MACH, the value of the constant.
c
      implicit none

      double precision r8_mach
      integer i

      if ( i .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
        write ( *, '(a,i12)' ) '  I = ', i
        r8_mach = 0.0D+00
        stop
      else if ( i .eq. 1 ) then
        r8_mach = 4.450147717014403D-308
      else if ( i .eq. 2 ) then
        r8_mach = 8.988465674311579D+307
      else if ( i .eq. 3 ) then
        r8_mach = 1.110223024625157D-016
      else if ( i .eq. 4 ) then
        r8_mach = 2.220446049250313D-016
      else if ( i .eq. 5 ) then
        r8_mach = 0.301029995663981D+000
      else if ( 5 .lt. i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_MACH - Fatal error!'
        write ( *, '(a)' ) '  The input argument I is out of bounds.'
        write ( *, '(a)' ) '  Legal values satisfy 1 <= I <= 5.'
        write ( *, '(a,i12)' ) '  I = ', i
        r8_mach = 0.0D+00
        stop
      end if

      return
      end
      function r8_uniform_01 ( seed )

c*********************************************************************72
c
cc R8_UNIFORM_01 returns a pseudorandom R8 scaled to [0,1].
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
      subroutine r8mat_cholesky_factor ( n, a, c, flag )

c*********************************************************************72
c
cc R8MAT_CHOLESKY_FACTOR computes the Cholesky factor of a symmetric matrix.
c
c  Discussion:
c
c    An R8MAT is an array of R8's.
c
c    The matrix must be symmetric and positive semidefinite.
c
c    For a positive semidefinite symmetric matrix A, the Cholesky factorization
c    is a lower triangular matrix L such that:
c
c      A = L * L'
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 April 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of rows and columns of
c    the matrix A.
c
c    Input, double precision A(N,N), the N by N matrix.
c
c    Output, double precision C(N,N), the N by N lower triangular
c    Cholesky factor.
c
c    Output, integer FLAG:
c    0, no error occurred.
c    1, the matrix is not positive definite.
c
      implicit none

      integer n

      double precision a(n,n)
      double precision c(n,n)
      integer flag
      integer i
      integer j
      integer k
      double precision sum2

      flag = 0

      do j = 1, n
        do i = 1, n
          c(i,j) = a(i,j)
        end do
      end do

      do j = 1, n

        do i = 1, j - 1
          c(i,j) = 0.0D+00
        end do

        do i = j, n

          sum2 = 0.0D+00
          do k = 1, j - 1
            sum2 = sum2 + c(j,k) * c(i,k)
          end do
          sum2 = c(j,i) - sum2

          if ( i .eq. j ) then
            if ( sum2 .le. 0.0D+00 ) then
              flag = 1
              return
            else
              c(i,j) = sqrt ( sum2 )
            end if
          else
            if ( c(j,j) .ne. 0.0D+00 ) then
              c(i,j) = sum2 / c(j,j)
            else
              c(i,j) = 0.0D+00
            end if
          end if

        end do

      end do

      return
      end
      subroutine r8mat_is_symmetric ( m, n, a, error_frobenius )

c*********************************************************************72
c
cc R8MAT_IS_SYMMETRIC checks an R8MAT for symmetry.
c
c  Discussion:
c
c    An R8MAT is a matrix of double precision values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    16 October 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the order of the matrix.
c
c    Input, double precision A(M,N), the matrix.
c
c    Output, double precision ERROR_FROBENIUS, measures the 
c    Frobenius norm of ( A - A' ), which would be zero if the matrix
c    were exactly symmetric.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision error_frobenius
      integer i
      integer j
      double precision r8_huge

      if ( m .ne. n ) then
        error_frobenius = r8_huge ( )
        return
      end if
 
      error_frobenius = 0.0D+00
      do j = 1, n
        do i = 1, m
          error_frobenius = error_frobenius + ( a(i,j) - a(j,i) )**2
        end do
      end do
      error_frobenius = sqrt ( error_frobenius )

      return
      end
      function r8mat_max ( m, n, a )

c*********************************************************************72
c
cc R8MAT_MAX returns the maximum entry of an R8MAT.
c
c  Discussion:
c
c    An R8MAT is an array of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows in A.
c
c    Input, integer N, the number of columns in A.
c
c    Input, double precision A(M,N), the matrix.
c
c    Output, double precision R8MAT_MAX, the maximum entry of A.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer i
      integer j
      double precision r8mat_max
      double precision value

      value = a(1,1)
      do j = 1, n
        do i = 1, m
          value = max ( value, a(i,j) )
        end do
      end do

      r8mat_max = value

      return
      end
      function r8mat_min ( m, n, a )

c*********************************************************************72
c
cc R8MAT_MIN returns the minimum entry of an R8MAT.
c
c  Discussion:
c
c    An R8MAT is an array of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows in A.
c
c    Input, integer N, the number of columns in A.
c
c    Input, double precision A(M,N), the matrix.
c
c    Output, double precision R8MAT_MIN, the minimum entry of A.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer i
      integer j
      double precision r8mat_min
      double precision value

      value = a(1,1)
      do j = 1, n
        do i = 1, m
          value = min ( value, a(i,j) )
        end do
      end do

      r8mat_min = value

      return
      end
      subroutine r8mat_mm ( n1, n2, n3, a, b, c )

c*********************************************************************72
c
cc R8MAT_MM multiplies two R8MAT's.
c
c  Discussion:
c
c    An R8MAT is an array of R8 values.
c
c    In FORTRAN90, this operation is more efficiently done by the
c    command:
c
c      C(1:N1,1:N3) = MATMUL ( A(1:N1,1;N2), B(1:N2,1:N3) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 December 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N1, N2, N3, the order of the matrices.
c
c    Input, double precision A(N1,N2), B(N2,N3), the matrices to multiply.
c
c    Output, double precision C(N1,N3), the product matrix C = A * B.
c
      implicit none

      integer n1
      integer n2
      integer n3

      double precision a(n1,n2)
      double precision b(n2,n3)
      double precision c(n1,n3)
      integer i
      integer j
      integer k

      do i = 1, n1
        do j = 1, n3
          c(i,j) = 0.0D+00
          do k = 1, n2
            c(i,j) = c(i,j) + a(i,k) * b(k,j)
          end do
        end do
      end do

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
      subroutine r8mat_print ( m, n, a, title )

c*********************************************************************72
c
cc R8MAT_PRINT prints an R8MAT.
c
c  Discussion:
c
c    An R8MAT is an array of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 May 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows in A.
c
c    Input, integer N, the number of columns in A.
c
c    Input, double precision A(M,N), the matrix.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      character * ( * ) title

      call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

      return
      end
      subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi,
     &  title )

c*********************************************************************72
c
cc R8MAT_PRINT_SOME prints some of an R8MAT.
c
c  Discussion:
c
c    An R8MAT is an array of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, double precision A(M,N), an M by N matrix to be printed.
c
c    Input, integer ILO, JLO, the first row and column to print.
c
c    Input, integer IHI, JHI, the last row and column to print.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 5 )
      integer m
      integer n

      double precision a(m,n)
      character * ( 14 ) ctemp(incx)
      integer i
      integer i2hi
      integer i2lo
      integer ihi
      integer ilo
      integer inc
      integer j
      integer j2
      integer j2hi
      integer j2lo
      integer jhi
      integer jlo
      character * ( * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )

      if ( m .le. 0 .or. n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  (None)'
        return
      end if

      do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

        j2hi = j2lo + incx - 1
        j2hi = min ( j2hi, n )
        j2hi = min ( j2hi, jhi )

        inc = j2hi + 1 - j2lo

        write ( *, '(a)' ) ' '

        do j = j2lo, j2hi
          j2 = j + 1 - j2lo
          write ( ctemp(j2), '(i7,7x)') j
        end do

        write ( *, '(''  Col   '',5a14)' ) ( ctemp(j), j = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) ' '

        i2lo = max ( ilo, 1 )
        i2hi = min ( ihi, m )

        do i = i2lo, i2hi

          do j2 = 1, inc

            j = j2lo - 1 + j2

            write ( ctemp(j2), '(g14.6)' ) a(i,j)

          end do

          write ( *, '(i5,a,5a14)' ) i, ':', ( ctemp(j), j = 1, inc )

        end do

      end do

      return
      end
      subroutine r8vec_linspace ( n, a_first, a_last, a )

c*********************************************************************72
c
cc R8VEC_LINSPACE creates a vector of linearly spaced values.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    4 points evenly spaced between 0 and 12 will yield 0, 4, 8, 12.
c
c    In other words, the interval is divided into N-1 even subintervals,
c    and the endpoints of intervals are used as the points.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    14 March 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, double precision A_FIRST, A_LAST, the first and last entries.
c
c    Output, double precision A(N), a vector of linearly spaced data.
c
      implicit none

      integer n

      double precision a(n)
      double precision a_first
      double precision a_last
      integer i

      if ( n .eq. 1 ) then

        a(1) = ( a_first + a_last ) / 2.0D+00

      else

        do i = 1, n
          a(i) = ( dble ( n - i     ) * a_first 
     &           + dble (     i - 1 ) * a_last )
     &           / dble ( n     - 1 )
        end do

      end if

      return
      end
      subroutine r8vec_min ( n, a, amin )

c*********************************************************************72
c
cc R8VEC_MIN returns the minimum value in an R8VEC.
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
c    31 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in the array.
c
c    Input, double precision A(N), the array.
c
c    Output, double precision AMIN, the value of the smallest entry.
c
      implicit none

      integer n

      double precision a(n)
      double precision amin
      integer i

      amin = a(1)
      do i = 2, n
        amin = min ( amin, a(i) )
      end do

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
          write ( *, '(a)' ) 'R8VEC_NORMAL_01 - Fatal error!'
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
      subroutine r8vec_print ( n, a, title )

c*********************************************************************72
c
cc R8VEC_PRINT prints an R8VEC.
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
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of components of the vector.
c
c    Input, double precision A(N), the vector to be printed.
c
c    Input, character * ( * ) TITLE, a title.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      character * ( * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i8,a,1x,g16.8)' ) i, ':', a(i)
      end do

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
      subroutine sample_paths_cholesky ( n, n2, rhomax, rho0, 
     &  correlation, seed, x )

c*********************************************************************72
c
cc SAMPLE_PATHS_CHOLESKY: sample paths for stationary correlation functions.
c
c  Discussion:
c
c    This method uses the Cholesky factorization of the correlation matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of points on each path.
c
c    Input, integer N2, the number of paths.
c
c    Input, double precision RHOMAX, the maximum value of RHO.
c
c    Input, double precision RHO0, the correlation length.
c
c    Input, external CORRELATION, the name of the subroutine which evaluates
c    the correlation.
c
c    Input/output, integer SEED, a seed for the random number
c    generator.
c
c    Output, double precision X(N,N2), the sample paths.
c
      implicit none

      integer n
      integer n2

      double precision cor(n,n)
      double precision cor_vec(n)
      external correlation
      integer flag
      integer i
      integer i4_wrap
      integer j
      integer k
      double precision l(n,n)
      double precision r(n,n2)
      double precision rho_vec(n)
      double precision rho0
      double precision rhomax
      double precision rhomin
      integer seed
      double precision x(n,n2)
c
c  Choose N equally spaced sample points from 0 to RHOMAX.
c
      rhomin = 0.0D+00
      call r8vec_linspace ( n, rhomin, rhomax, rho_vec )
c
c  Evaluate the correlation function.
c
      call correlation ( n, rho_vec, rho0, cor_vec )
c
c  Construct the correlation matrix;
c
c  From the vector 
c    [ C(0), C(1), C(2), ... C(N-1) ]
c  construct the vector
c    [ C(N-1), ..., C(2), C(1), C(0), C(1), C(2), ...  C(N-1) ]
c  Every row of the correlation matrix can be constructed by a subvector
c  of this vector.
c
      do j = 1, n
        do i = 1, n
          k = i4_wrap ( abs ( j - i ) + 1, 1, n )
          cor(i,j) = cor_vec(k)
        end do
      end do
c
c  Get the Cholesky factorization of COR:
c
c    COR = L * L'.
c
      call r8mat_cholesky_factor ( n, cor, l, flag )
c
c  The matrix might not be nonnegative definite.
c
      if ( flag .eq. 2 ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'SAMPLE_PATHS_CHOLESKY - Fatal error!'
        write ( *, '(a)' ) '  The correlation matrix is not'
        write ( *, '(a)' ) '  symmetric nonnegative definite.'
        stop
      end if
c
c  Compute a matrix of N by N2 normally distributed values.
c
      call r8mat_normal_01 ( n, n2, seed, r )
c
c  Compute the sample path.
c
      call r8mat_mm ( n, n, n2, l, r, x )

      return
      end
      subroutine sample_paths_eigen ( n, n2, rhomax, rho0, 
     &  correlation, seed, x )

c*********************************************************************72
c
cc SAMPLE_PATHS_EIGEN: sample paths for stationary correlation functions.
c
c  Discussion:
c
c    This method uses the eigen-decomposition of the correlation matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of points on each path.
c
c    Input, integer N2, the number of paths.
c
c    Input, double precision RHOMAX, the maximum value of RHO.
c
c    Input, double precision RHO0, the correlation length.
c
c    Input, external CORRELATION, the name of the subroutine which evaluates
c    the correlation.
c
c    Input/output, integer SEED, a seed for the random number
c    generator.
c
c    Output, double precision X(N,N2), the sample paths.
c
      implicit none

      integer n
      integer n2

      double precision c(n,n)
      double precision cor(n,n)
      double precision cor_vec(n)
      external correlation
      double precision d(n)
      double precision dmin
      integer i
      integer i4_wrap
      integer ierr
      integer j
      integer k
      double precision r(n,n2)
      double precision r8_epsilon
      double precision r8vec_min
      double precision rho_vec(n)
      double precision rho0
      double precision rhomax
      double precision rhomin
      integer seed
      double precision v(n,n)
      double precision w(n)
      double precision x(n,n2)
c
c  Choose N equally spaced sample points from 0 to RHOMAX.
c
      rhomin = 0.0D+00
      call r8vec_linspace ( n, rhomin, rhomax, rho_vec )
c
c  Evaluate the correlation function.
c
      call correlation ( n, rho_vec, rho0, cor_vec )
c
c  Construct the correlation matrix;
c
c  From the vector 
c    [ C(0), C(1), C(2), ... C(N-1) ]
c  construct the vector
c    [ C(N-1), ..., C(2), C(1), C(0), C(1), C(2), ...  C(N-1) ]
c  Every row of the correlation matrix can be constructed by a subvector
c  of this vector.
c
      do j = 1, n
        do i = 1, n
          k = i4_wrap ( abs ( i - j ), 0, n - 1 ) + 1
          cor(i,j) = cor_vec(k)
        end do
      end do
c
c  Get the eigendecomposition of COR:
c
c    COR = V * D * V'.
c
c  Because COR is symmetric, V is orthogonal.
c
      call tred2 ( n, n, cor, d, w, v )

      call tql2 ( n, n, d, w, v, ierr )
c
c  We assume COR is non-negative definite, and hence that there
c  are no negative eigenvalues.
c
      dmin = r8vec_min ( n, d )

      if ( dmin .lt. - sqrt ( r8_epsilon ( ) ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SAMPLE_PATHS_EIGEN - Warning!'
        write ( *, '(a,g14.6)' ) 
     &    '  Negative eigenvalues observed as low as ', dmin
      end if

      do i = 1, n
        d(i) = max ( d(i), 0.0D+00 )
      end do
c
c  Compute the eigenvalues of the factor C.
c
      do i = 1, n
        d(i) = sqrt ( d(i) )
      end do
c
c  Compute C, such that C' * C = COR.
c
      do j = 1, n
        do i = 1, n
          c(i,j) = 0.0D+00
          do k = 1, n
            c(i,j) = c(i,j) + d(k) * v(i,k) * v(j,k)
          end do
        end do
      end do
c
c  Compute N by N2 independent random normal values.
c
      call r8mat_normal_01 ( n, n2, seed, r )
c
c  Multiply to get the variables X which have correlation COR.
c
      call r8mat_mm ( n, n, n2, c, r, x )

      return
      end
      subroutine sample_paths2_cholesky ( n, n2, rhomin, rhomax, rho0, 
     &  correlation2, seed, x )

c*********************************************************************72
c
cc SAMPLE_PATHS2_CHOLESKY: sample paths for stationary correlation functions.
c
c  Discussion:
c
c    This method uses the Cholesky factorization of the correlation matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of points on each path.
c
c    Input, integer N2, the number of paths.
c
c    Input, double precision RHOMIN, RHOMAX, the range of RHO.
c
c    Input, double precision RHO0, the correlation length.
c
c    Input, external CORRELATION2, the name of the subroutine which evaluates
c    the correlation.
c
c    Input/output, integer SEED, a seed for the random number
c    generator.
c
c    Output, double precision X(N,N2), the sample paths.
c
      implicit none

      integer n
      integer n2

      double precision cor(n,n)
      external correlation2
      integer flag
      integer i
      integer j
      integer k
      double precision l(n,n)
      double precision r(n,n2) 
      double precision rho0
      double precision rhomax
      double precision rhomin
      double precision s(n)
      integer seed
      double precision x(n,n2)
c
c  Choose N equally spaced sample points from RHOMIN to RHOMAX.
c
      call r8vec_linspace ( n, rhomin, rhomax, s )
c
c  Evaluate the correlation function.
c
      call correlation2 ( n, n, s, s, rho0, cor )
c
c  Get the Cholesky factorization of COR:
c
c    COR = L * L'.
c
      call r8mat_cholesky_factor ( n, cor, l, flag )
c
c  The matrix might not be nonnegative definite.
c
      if ( flag .eq. 2 ) then
        write ( *, '(a)' ) ''
        write ( *, '(a)' ) 'SAMPLE_PATHS2_CHOLESKY - Fatal error!'
        write ( *, '(a)' ) '  The correlation matrix is not'
        write ( *, '(a)' ) '  symmetric nonnegative definite.'
        stop
      end if
c
c  Compute a matrix of N by N2 normally distributed values.
c
      call r8mat_normal_01 ( n, n2, seed, r )
c
c  Compute the sample path.
c
      call r8mat_mm ( n, n, n2, l, r, x )

      return
      end
      subroutine sample_paths2_eigen ( n, n2, rhomin, rhomax, rho0, 
     &  correlation2, seed, x )

c*********************************************************************72
c
cc SAMPLE_PATHS2_EIGEN: sample paths for stationary correlation functions.
c
c  Discussion:
c
c    This method uses the eigen-decomposition of the correlation matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 November 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of points on each path.
c
c    Input, integer N2, the number of paths.
c
c    Input, double precision RHOMIN, RHOMAX, the range of RHO.
c
c    Input, double precision RHO0, the correlation length.
c
c    Input, external CORRELATION2, the name of the subroutine which evaluates
c    the correlation.
c
c    Input/output, integer SEED, a seed for the random number
c    generator.
c
c    Output, double precision X(N,N2), the sample paths.
c
      implicit none

      integer n
      integer n2

      double precision c(n,n)
      double precision cor(n,n)
      external correlation2
      double precision d(n)
      double precision dmin
      integer i
      integer ierr
      integer j
      integer k
      double precision r(n,n2)
      double precision r8_epsilon
      double precision r8vec_min
      double precision rho0
      double precision rhomax
      double precision rhomin
      double precision s(n)
      integer seed
      double precision v(n,n)
      double precision w(n)
      double precision x(n,n2)
c
c  Choose N equally spaced sample points from RHOMIN to RHOMAX.
c
      call r8vec_linspace ( n, rhomin, rhomax, s )
c
c  Evaluate the correlation function.
c
      call correlation2 ( n, n, s, s, rho0, cor )
c
c  Get the eigendecomposition of COR:
c
c    COR = V * D * V'.
c
c  Because COR is symmetric, V is orthogonal.
c
      call tred2 ( n, n, cor, d, w, v )

      call tql2 ( n, n, d, w, v, ierr )
c
c  We assume COR is non-negative definite, and hence that there
c  are no negative eigenvalues.
c
      dmin = r8vec_min ( n, d )

      if ( dmin .lt. - sqrt ( r8_epsilon ( ) ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SAMPLE_PATHS2_EIGEN - Warning!'
        write ( *, '(a,g14.6)' ) 
     &    '  Negative eigenvalues observed as low as ', dmin
      end if

      do i = 1, n
        d(i) = max ( d(i), 0.0D+00 )
      end do
c
c  Compute the eigenvalues of the factor C.
c
      do i = 1, n
        d(i) = sqrt ( d(i) )
      end do
c
c  Compute C, such that C' * C = COR.
c
      do j = 1, n
        do i = 1, n
          c(i,j) = 0.0D+00
          do k = 1, n
            c(i,j) = c(i,j) + d(k) * v(i,k) * v(j,k)
          end do
        end do
      end do
c
c  Compute N by N2 independent random normal values.
c
      call r8mat_normal_01 ( n, n2, seed, r )
c
c  Multiply to get the variables X which have correlation COR.
c
      call r8mat_mm ( n, n, n2, c, r, x )

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
      subroutine tql2(nm,n,d,e,z,ierr)

c*********************************************************************72
c
cc TQL2 computes all eigenvalues/vectors, real symmetric tridiagonal matrix.
c
c     this subroutine is a translation of the algol procedure tql2,
c     num. math. 11, 293-306(1968) by bowdler, martin, reinsch, and
c     wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 227-240(1971).
c
c     this subroutine finds the eigenvalues and eigenvectors
c     of a symmetric tridiagonal matrix by the ql method.
c     the eigenvectors of a full symmetric matrix can also
c     be found if  tred2  has been used to reduce this
c     full matrix to tridiagonal form.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        d contains the diagonal elements of the input matrix.
c
c        e contains the subdiagonal elements of the input matrix
c          in its last n-1 positions.  e(1) is arbitrary.
c
c        z contains the transformation matrix produced in the
c          reduction by  tred2, if performed.  if the eigenvectors
c          of the tridiagonal matrix are desired, z must contain
c          the identity matrix.
c
c      on output
c
c        d contains the eigenvalues in ascending order.  if an
c          error exit is made, the eigenvalues are correct but
c          unordered for indices 1,2,...,ierr-1.
c
c        e has been destroyed.
c
c        z contains orthonormal eigenvectors of the symmetric
c          tridiagonal (or full) matrix.  if an error exit is made,
c          z contains the eigenvectors associated with the stored
c          eigenvalues.
c
c        ierr is set to
c          zero       for normal return,
c          j          if the j-th eigenvalue has not been
c                     determined after 30 iterations.
c
c     calls pythag for  dsqrt(a*a + b*b) .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,l,m,n,ii,l1,l2,nm,mml,ierr
      double precision d(n),e(n),z(nm,n)
      double precision c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag

      ierr = 0
      if (n .eq. 1) go to 1001
c
      do 100 i = 2, n
  100 e(i-1) = e(i)
c
      f = 0.0d0
      tst1 = 0.0d0
      e(n) = 0.0d0
c
      do 240 l = 1, n
         j = 0
         h = dabs(d(l)) + dabs(e(l))
         if (tst1 .lt. h) tst1 = h
c     .......... look for small sub-diagonal element ..........
         do 110 m = l, n
            tst2 = tst1 + dabs(e(m))
            if (tst2 .eq. tst1) go to 120
c     .......... e(n) is always zero, so there is no exit
c                through the bottom of the loop ..........
  110    continue
c
  120    if (m .eq. l) go to 220
  130    if (j .eq. 30) go to 1000
         j = j + 1
c     .......... form shift ..........
         l1 = l + 1
         l2 = l1 + 1
         g = d(l)
         p = (d(l1) - g) / (2.0d0 * e(l))
         r = pythag(p,1.0d0)
         d(l) = e(l) / (p + dsign(r,p))
         d(l1) = e(l) * (p + dsign(r,p))
         dl1 = d(l1)
         h = g - d(l)
         if (l2 .gt. n) go to 145
c
         do 140 i = l2, n
  140    d(i) = d(i) - h
c
  145    f = f + h
c     .......... ql transformation ..........
         p = d(m)
         c = 1.0d0
         c2 = c
         el1 = e(l1)
         s = 0.0d0
         mml = m - l
c     .......... for i=m-1 step -1 until l do -- ..........
         do 200 ii = 1, mml
            c3 = c2
            c2 = c
            s2 = s
            i = m - ii
            g = c * e(i)
            h = c * p
            r = pythag(p,e(i))
            e(i+1) = s * r
            s = e(i) / r
            c = p / r
            p = c * d(i) - s * g
            d(i+1) = h + s * (c * g + s * d(i))
c     .......... form vector ..........
            do 180 k = 1, n
               h = z(k,i+1)
               z(k,i+1) = s * z(k,i) + c * h
               z(k,i) = c * z(k,i) - s * h
  180       continue
c
  200    continue
c
         p = -s * s2 * c3 * el1 * e(l) / dl1
         e(l) = s * p
         d(l) = c * p
         tst2 = tst1 + dabs(e(l))
         if (tst2 .gt. tst1) go to 130
  220    d(l) = d(l) + f
  240 continue

c     .......... order eigenvalues and eigenvectors ..........
      do 300 ii = 2, n
         i = ii - 1
         k = i
         p = d(i)
c
         do 260 j = ii, n
            if (d(j) .ge. p) go to 260
            k = j
            p = d(j)
  260    continue
c
         if (k .eq. i) go to 300
         d(k) = d(i)
         d(i) = p
c
         do 280 j = 1, n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
  280    continue
c
  300 continue
c
      go to 1001
c     .......... set error -- no convergence to an
c                eigenvalue after 30 iterations ..........
 1000 ierr = l
 1001 return
      end
      subroutine tred2(nm,n,a,d,e,z)

c*********************************************************************72
c
cc TRED2 transforms a real symmetric matrix to symmetric tridiagonal form.
c
c     this subroutine is a translation of the algol procedure tred2,
c     num. math. 11, 181-195(1968) by martin, reinsch, and wilkinson.
c     handbook for auto. comp., vol.ii-linear algebra, 212-226(1971).
c
c     this subroutine reduces a real symmetric matrix to a
c     symmetric tridiagonal matrix using and accumulating
c     orthogonal similarity transformations.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrix.
c
c        a contains the real symmetric input matrix.  only the
c          lower triangle of the matrix need be supplied.
c
c     on output
c
c        d contains the diagonal elements of the tridiagonal matrix.
c
c        e contains the subdiagonal elements of the tridiagonal
c          matrix in its last n-1 positions.  e(1) is set to zero.
c
c        z contains the orthogonal transformation matrix
c          produced in the reduction.
c
c        a and z may coincide.  if distinct, a is unaltered.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
      integer i,j,k,l,n,ii,nm,jp1
      double precision a(nm,n),d(n),e(n),z(nm,n)
      double precision f,g,h,hh,scale

      do 100 i = 1, n
c
         do 80 j = i, n
   80    z(j,i) = a(j,i)
c
         d(i) = a(n,i)
  100 continue
c
      if (n .eq. 1) go to 510
c     .......... for i=n step -1 until 2 do -- ..........
      do 300 ii = 2, n
         i = n + 2 - ii
         l = i - 1
         h = 0.0d0
         scale = 0.0d0
         if (l .lt. 2) go to 130
c     .......... scale row (algol tol then not needed) ..........
         do 120 k = 1, l
  120    scale = scale + dabs(d(k))
c
         if (scale .ne. 0.0d0) go to 140
  130    e(i) = d(l)
c
         do 135 j = 1, l
            d(j) = z(l,j)
            z(i,j) = 0.0d0
            z(j,i) = 0.0d0
  135    continue
c
         go to 290
c
  140    do 150 k = 1, l
            d(k) = d(k) / scale
            h = h + d(k) * d(k)
  150    continue
c
         f = d(l)
         g = -dsign(dsqrt(h),f)
         e(i) = scale * g
         h = h - f * g
         d(l) = f - g
c     .......... form a*u ..........
         do 170 j = 1, l
  170    e(j) = 0.0d0
c
         do 240 j = 1, l
            f = d(j)
            z(j,i) = f
            g = e(j) + z(j,j) * f
            jp1 = j + 1
            if (l .lt. jp1) go to 220
c
            do 200 k = jp1, l
               g = g + z(k,j) * d(k)
               e(k) = e(k) + z(k,j) * f
  200       continue
c
  220       e(j) = g
  240    continue
c     .......... form p ..........
         f = 0.0d0
c
         do 245 j = 1, l
            e(j) = e(j) / h
            f = f + e(j) * d(j)
  245    continue
c
         hh = f / (h + h)
c     .......... form q ..........
         do 250 j = 1, l
  250    e(j) = e(j) - hh * d(j)
c     .......... form reduced a ..........
         do 280 j = 1, l
            f = d(j)
            g = e(j)
c
            do 260 k = j, l
  260       z(k,j) = z(k,j) - f * e(k) - g * d(k)
c
            d(j) = z(l,j)
            z(i,j) = 0.0d0
  280    continue
c
  290    d(i) = h
  300 continue
c     .......... accumulation of transformation matrices ..........
      do 500 i = 2, n
         l = i - 1
         z(n,l) = z(l,l)
         z(l,l) = 1.0d0
         h = d(i)
         if (h .eq. 0.0d0) go to 380
c
         do 330 k = 1, l
  330    d(k) = z(k,i) / h
c
         do 360 j = 1, l
            g = 0.0d0
c
            do 340 k = 1, l
  340       g = g + z(k,i) * z(k,j)
c
            do 360 k = 1, l
               z(k,j) = z(k,j) - g * d(k)
  360    continue
c
  380    do 400 k = 1, l
  400    z(k,i) = 0.0d0
c
  500 continue
c
  510 do 520 i = 1, n
         d(i) = z(n,i)
         z(n,i) = 0.0d0
  520 continue
c
      z(n,n) = 1.0d0
      e(1) = 0.0d0
      return
      end
