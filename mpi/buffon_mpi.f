      program main

c*********************************************************************72
c
cc MAIN is the main program for BUFFON_LAPLACE.
c
c  Discussion:
c
c    This program uses MPI to do a Buffon-Laplace simulation in parallel.
c
c    This is an example of an "embarassingly parallel" computation.  Each
c    processor does the same randomized computation, counting "successes",
c    and returning the number of successes to the master process.
c
c    The particular point made by this example is that each process must
c    use a different stream of random numbers, and that this can be done
c    if each process uses a different seed to initialize the random number
c    generator, and that this, in turn, can be done by simply adding the
c    process's rank to a common seed value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 February 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    William Gropp, Ewing Lusk, Anthony Skjellum,
c    Using MPI: Portable Parallel Programming with the
c    Message-Passing Interface,
c    Second Edition,
c    MIT Press, 1999,
c    ISBN: 0262571323,
c    LC: QA76.642.G76.
c
c    Marc Snir, Steve Otto, Steven Huss-Lederman, David Walker, 
c    Jack Dongarra,
c    MPI: The Complete Reference,
c    Volume I: The MPI Core,
c    Second Edition,
c    MIT Press, 1998,
c    ISBN: 0-262-69216-3,
c     LC: QA76.642.M65.
c
c    Sudarshan Raghunathan,
c    Making a Supercomputer Do What You Want: High Level Tools for
c    Parallel Programming,
c    Computing in Science and Engineering,
c    Volume 8, Number 5, September/October 2006, pages 70-80.
c
      include 'mpif.h'

      double precision a
      double precision b
      integer buffon_laplace_simulate
      integer hit_num
      integer hit_total
      integer ierr
      double precision l
      integer master
      double precision pdf_estimate
      double precision pi
      double precision pi_error
      double precision pi_estimate
      integer process_num
      integer process_rank
      double precision r8_uniform_01
      double precision random_value
      integer seed
      integer seed_init
      integer trial_num
      integer trial_total

      a = 1.0D+00
      b = 1.0D+00
      l = 1.0D+00
      master = 0
      pi = 3.141592653589793238462643D+00
      trial_num = 100000
c
c  Initialize MPI.
c
      call MPI_Init ( ierr )
c
c  Get the number of processes.
c
      call MPI_Comm_size ( MPI_COMM_WORLD, process_num, ierr )
c
c  Get the rank of this processor.
c
      call MPI_Comm_rank ( MPI_COMM_WORLD, process_rank, ierr )
c
c  The master process sets the value of the parameters,
c  and broadcasts them.
c
      if ( process_rank .eq. master ) then
        call timestamp ( )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BUFFON - Master process:'
        write ( *, '(a)' ) '  FORTRAN77 version'
        write ( *, '(a)' ) '  An MPI example program to estimate PI'
        write ( *, '(a)' ) '  in the Buffon-Laplace needle experiment.'
        write ( *, '(a)' ) 
     &  '  On a grid of cells of width A and height B,'
        write ( *, '(a)' ) '  a needle of length L is dropped.'
        write ( *, '(a)' ) '  We count the number of times it crosses'
        write ( *, '(a)' ) '  at least one grid line, and use this to'
        write ( *, '(a)' ) '  estimate the value of PI.'
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) 
     &  '  The number of processes is ', process_num
        write ( *, '(a)' ) ' '
        write ( *, '(a,f14.6)' ) '  Cell width A =    ', a
        write ( *, '(a,f14.6)' ) '  Cell height B =   ', b
        write ( *, '(a,f14.6)' ) '  Needle length L = ', l
      end if
c
c  Each process sets a random number seed.
c
      seed_init = 123456789 + process_rank * 100
      seed = seed_init
c
c  Just to make sure that we're all doing different things, have each
c  process print out its rank, seed value, and a first test random value.
c
      random_value = r8_uniform_01 ( seed )

      write ( *, '(2x,i8,2x,i12,2x,g14.6)' ) 
     &  process_rank, seed_init, random_value
c
c  Each process now carries out TRIAL_NUM trials, and then
c  sends the value back to the master process.
c
      hit_num = buffon_laplace_simulate ( a, b, l, trial_num, seed ) 

      call MPI_Reduce ( hit_num, hit_total, 1, MPI_INTEGER, 
     &  MPI_SUM, master, MPI_COMM_WORLD, ierr )
c
c  The master process can now estimate PI.
c
      if ( process_rank .eq. master ) then

        trial_total = trial_num * process_num

        pdf_estimate = dble ( hit_total ) / dble ( trial_total )

        if ( hit_total .eq. 0 ) then

          pi_estimate = r8_huge ( )

        else

          pi_estimate = l * ( 2.0D+00 * ( a + b ) - l ) 
     &      / ( a * b * pdf_estimate )

        end if

        pi_error = abs ( pi - pi_estimate )

        write ( *, '(a)' ) ' '
        write ( *, '(a,a)' ) 
     &    '    Trials      Hits    Estimated PDF  ',
     &    '     Estimated Pi         Error'
        write ( *, '(a)' ) ' '
        write ( *, '(2x,i8,2x,i8,2x,g20.12,2x,g20.12,2x,g20.12)' ) 
     &    trial_total, hit_total, pdf_estimate, pi_estimate, pi_error

      end if
c
c  Terminate MPI.
c
      call MPI_Finalize ( ierr )
c
c  Terminate.
c
      if ( process_rank .eq. master ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BUFFON - Master process:'
        write ( *, '(a)' ) '  Normal end of execution.'
        write ( *, '(a)' ) ' '
        call timestamp ( )
      end if

      stop
      end
      function buffon_laplace_simulate ( a, b, l, trial_num, seed )

c*********************************************************************72
c
cc BUFFON_LAPLACE_SIMULATE simulates a Buffon-Laplace needle experiment.
c
c  Discussion:
c
c    In the Buffon-Laplace needle experiment, we suppose that the plane has been
c    tiled into a grid of rectangles of width A and height B, and that a
c    needle of length L is dropped "at random" onto this grid.
c
c    We may assume that one end, the "eye" of the needle falls at the point
c    (X1,Y1), taken uniformly at random in the cell [0,A]x[0,B].
c
c    ANGLE, the angle that the needle makes is taken to be uniformly random.
c    The point of the needle, (X2,Y2), therefore lies at
c
c      (X2,Y2) = ( X1+L*cos(ANGLE), Y1+L*sin(ANGLE) )
c
c    The needle will have crossed at least one grid line if any of the
c    following are true:
c
c      X2 <= 0, A <= X2, Y2 <= 0, B <= Y2.
c
c    This routine simulates the tossing of the needle, and returns the number
c    of times that the needle crossed at least one grid line.
c
c    If L is larger than sqrt ( A*A + B*B ), then the needle will
c    cross every time, and the computation is uninteresting.  However, if
c    L is smaller than this limit, then the probability of a crossing on
c    a single trial is
c
c      P(L,A,B) = ( 2 * L * ( A + B ) - L * L ) / ( PI * A * B )
c
c    and therefore, a record of the number of hits for a given number of
c    trials can be used as a very roundabout way of estimating PI.
c    (Particularly roundabout, since we actually will use a good value of
c    PI in order to pick the random anglesc)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 February 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Sudarshan Raghunathan,
c    Making a Supercomputer Do What You Want: High Level Tools for
c    Parallel Programming,
c    Computing in Science and Engineering,
c    Volume 8, Number 5, September/October 2006, pages 70-80.
c
c  Parameters:
c
c    Input, double precision A, B, the horizontal and vertical dimensions
c    of each cell of the grid.  0 <= A, 0 <= B.
c
c    Input, double precision L, the length of the needle.
c    0 <= L <= min ( A, B ).
c
c    Input, integer TRIAL_NUM, the number of times the needle is
c    to be dropped onto the grid.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, integer BUFFON_LAPLACE_SIMULATE, the number of times the needle
c    crossed at least one line of the grid of cells.
c
      implicit none

      double precision a
      double precision angle
      double precision b
      integer buffon_laplace_simulate
      integer hits
      double precision l
      double precision pi
      parameter ( pi = 3.141592653589793238462643D+00 )
      double precision r8_uniform_01
      integer seed
      integer trial
      integer trial_num
      double precision x1
      double precision x2
      double precision y1
      double precision y2

      hits = 0

      do trial = 1, trial_num
c
c  Randomly choose the location of the eye of the needle in [0,0]x[A,B],
c  and the angle the needle makes.
c
        x1 = a * r8_uniform_01 ( seed )
        y1 = b * r8_uniform_01 ( seed )
	angle = 2.0D+00 * pi * r8_uniform_01 ( seed )
c
c  Compute the location of the point of the needle.
c
        x2 = x1 + l * cos ( angle )
        y2 = y1 + l * sin ( angle )
c
c  Count the end locations that lie outside the cell.
c
        if ( x2 .le. 0.0D+00 .or. a .le. x2 .or.
     &       y2 .le. 0.0D+00 .or. b .le. y2 ) then
          hits = hits + 1
        end if

      end do

      buffon_laplace_simulate = hits

      return
      end
      function i4_huge ( )

c*********************************************************************72
c
cc I4_HUGE returns a "huge" I4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer I4_HUGE, a huge number.
c
      implicit none

      integer i4_huge

      i4_huge = 2147483647

      return
      end
      function r8_huge ( )

c*********************************************************************72
c
cc R8_HUGE returns a "huge" R8.
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

      integer i4_huge
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
        seed = seed + i4_huge ( )
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
