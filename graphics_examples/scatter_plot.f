      program main

c*****************************************************************************80
c
cc SCATTER_PLOT uses DISLIN to draw a scatter plot of X Y data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 April 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Helmut Michels,
c    The Data Plotting Software DISLIN - version 10.4,
c    Shaker Media GmbH, January 2010,
c    ISBN13: 978-3-86858-517-9.
c
      implicit none

      integer n
      parameter ( n = 500 )

      integer i
      integer j
      integer nr
      integer nx
      integer ny
      integer pat
      real r
      real r4_uniform_01
      real s
      integer seed
      real x
      real xvec(n)
      real y
      real yvec(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCATTER_PLOT:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Use DISLIN routines to make a scatterplot.'
c
c  Generate the data.  
c  We average 4 random values to get data that tends to cluster
c  near (0.5,0.5).
c
      seed = 123456789

      do i = 1, n
        s = 0.0E+00
        do j = 1, 4
          s = s + r4_uniform_01 ( seed )
        end do
        xvec(i) = s / 4.0E+00
      end do

      do i = 1, n
        s = 0.0E+00
        do j = 1, 4
          s = s + r4_uniform_01 ( seed )
        end do
        yvec(i) = s / 4.0E+00
      end do
c
c  Specify the format of the output file.
c
      call metafl ( 'png' )
c
c  Indicate that new data overwrites old data.
c
      call filmod ( 'delete' )
c
c  Specify the name of the output graphics file.
c
      call setfil ( 'scatter_plot.png' )
c
c  Choose the page size and orientation.
c  'USA' is 2160 plot units wide and 2790 plot units high.
c  'P' requests PROFILE rather than LANDSCAPE orientation.
c
      call setpag ( 'usap' )
c
c  For PNG output, use reverse the default black background to white.
c
      call scrmod ( 'reverse' )
c
c  Open DISLIN.
c
      call disini ( )
c
c  Plot a border around the page.
c
      call pagera ( )
c
c  Use the COMPLEX font.
c
      call complx ( )
c
c  Define the X and Y sizes of the axis system in plot units.
c
      call axslen ( 1800, 1800 )
c
c  Specify how the lower X, left Y, upper X and right Y axes are labeled.
c
      call setgrf ( 'line', 'line', 'line', 'line' )
c
c  Set the axis origin 180 plot units to the right, and 2610 plot units DOWN.
c
      call axspos ( 180, 2610 )
c
c  Relate the physical coordinates to the axes.
c
      call graf ( 0.0, 1.0, 0.0, 0.1, 0.0, 1.0, 0.0, 0.1 )
c
c  Add a grid, with one grid line for every tick mark in the X and Y axes.
c
      call grid ( 1, 1 )
c
c  Select the shading pattern.
c
      pat = 16
      call shdpat ( pat )
c
c  Set the color to blue.
c
      call color ( 'blue' )
c
c  At every data point, draw a circle of radius 0.01.
c
      do i = 1, n
        call rlcirc ( xvec(i), yvec(i), 0.01 )
      end do
c
c  Select character height in plot units.
c
      call height ( 50 )
c
c  We want the title to be black, but we say "white" because they've
c  been reversed (remember?).
c
      call color ( 'white' )
c
c  Define axis system titles.
c
      call titlin ( 'Scatter Plot', 1 )
c
c  Draw the title.
c
      call title ( )
c
c  End this plot.
c
      call endgrf ( )
c
c  Close DISLIN.
c
      call disfin ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SCATTER_PLOT:'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
      function r4_uniform_01 ( seed )

c*********************************************************************72
c
cc R4_UNIFORM_01 returns a unit pseudorandom R4.
c
c  Discussion:
c
c    This routine implements the recursion
c
c      seed = 16807 * seed mod ( 2**31 - 1 )
c      r4_uniform_01 = seed / ( 2**31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      R4_UNIFORM_01
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
c    Output, real R4_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      integer k
      integer seed
      real r4_uniform_01

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R4_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + 2147483647
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r4_uniform_01 = real ( dble ( seed ) * 4.656612875E-10 )

      return
      end
