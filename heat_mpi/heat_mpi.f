      program main

c*********************************************************************72
c
cc  MAIN is the main program for HEAT_MPI.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 October 2011
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
      include 'mpif.h'

      integer id
      integer ierr
      integer p
      double precision wtime

      call MPI_Init ( ierr )

      call MPI_Comm_rank ( MPI_COMM_WORLD, id, ierr )

      call MPI_Comm_size ( MPI_COMM_WORLD, p, ierr )

      if ( id .eq. 0 ) then
        call timestamp ( )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'HEAT_MPI:'
        write ( *, '(a)' ) '  FORTRAN77/MPI version'
        write ( *, '(a)' ) 
     &    '  Solve the 1D time-dependent heat equation.'
      end if
c
c  Record the starting time.
c
      if ( id .eq. 0 ) then
        wtime = MPI_Wtime ( )
      end if

      call update ( id, p )
c
c  Record the final time.
c
      if ( id .eq. 0 ) then
        wtime = MPI_Wtime ( ) - wtime
        write ( *, '(a)' ) ' '
        write ( *, '(a,g14.6)' ) '  Wall clock elapsed seconds = ', wtime
      end if
c
c  Terminate MPI.
c
      call MPI_Finalize ( )
c
c  Terminate.
c
      if ( id .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'HEAT_MPI:'
        write ( *, '(a)' ) '  Normal end of execution.'
        write ( *, '(a)' ) ' '
        call timestamp ( )
      end if

      stop
      end
      subroutine update ( id, p )

c*********************************************************************72
c
cc UPDATE computes the solution of the heat equation.
c
c  Discussion:
c
c    If there is only one processor (P .eq. 1 ), then the program writes the
c    values of X and H to files.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 April 2008
c
c  Author:
c
c    John Burkardtc
c
c  Parameters:
c
c    Input, integer ID, the id of this processor.
c
c    Input, integer P, the number of processors.
c
      include 'mpif.h'

      integer n
      parameter ( n = 11 )

      double precision boundary_condition
      double precision cfl
      double precision h(0:n+1)
      integer h_file
      double precision h_new(0:n+1)
      integer i
      integer id
      double precision initial_condition
      integer j
      integer j_max
      integer j_min
      double precision k
      integer p
      double precision rhs
      integer status(MPI_STATUS_SIZE)
      integer tag
      double precision time
      double precision time_delta
      double precision time_max
      double precision time_min
      double precision time_new
      double precision x(0:n+1)
      double precision x_delta
      integer x_file
      double precision x_max
      double precision x_min

      h_file = 11
      j_max = 400
      j_min = 0
      k = 0.002
      x_file = 12
      time_max = 10.0
      time_min = 0.0
      x_max = 1.0
      x_min = 0.0
c
c  Have process 0 print out some information.
c
      if ( id .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &  '  Compute an approximate solution to the time dependent'
        write ( *, '(a)' ) '  one dimensional heat equation:'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '    dH/dt - K * d2H/dx2 = f(x,t)'
        write ( *, '(a)' ) ' '
        write ( *, '(a,g14.6,a,g14.6)' ) 
     &  '  for ', x_min, ' = x_min < x < x_max = ', x_max
        write ( *, '(a)' ) ' '
        write ( *, '(a,g14.6,a,g14.6)' ) 
     &  '  and ', time_min, ' = time_min < t <= t_max = ', time_max
        write ( *, '(a)' ) ' '
        write ( *, '(a)' )
     &  '  Boundary conditions are specified at x_min and x_max.'
        write ( *, '(a)' ) 
     &  '  Initial conditions are specified at time_min.'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' )
     &  '  The finite difference method is used to discretize the'
        write ( *, '(a)' ) '  differential equation.'
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8,a)' ) 
     &  '  This uses ', p * n, ' equally spaced points in X'
        write ( *, '(a,i8,a)' ) 
     &  '  and ', j_max, ' equally spaced points in time.'
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8,a)' ) 
     &  '  Parallel execution is done using ', p, ' processors.'
        write ( *, '(a)' ) '  Domain decomposition is used.'
        write ( *, '(a,i8,a)' ) 
     &  '  Each processor works on ', n, ' nodes,'
        write ( *, '(a)' ) 
     &  '  and shares some information with its immediate neighbors.'
      end if
c
c  Set the X coordinates of the N nodes.
c  We don't actually need ghost values of X but we'll throw them in
c  as X(0) and X(N+1).
c
      do i = 0, n + 1
        x(i) = ( dble (         id * n + i - 1 ) * x_max
     &         + dble ( p * n - id * n - i     ) * x_min )
     &         / dble ( p * n              - 1 )
      end do
c
c  In single processor mode, write out the X coordinates for display.
c
      if ( p .eq. 1 ) then

        open ( unit = x_file, file = 'x_data.txt', status = 'unknown' )

        write ( x_file, '(11f14.6)' ) x(1:n)

        close ( unit = x_file )

      end if
c
c  Set the values of H at the initial time.
c
      time = time_min

      h(0) = 0.0
      do i = 1, n
        h(i) = initial_condition ( x(i), time )
      end do
      h(n+1) = 0.0
  
      time_delta = ( time_max - time_min ) / dble ( j_max - j_min )
      x_delta = ( x_max - x_min ) / dble ( p * n - 1 )
c
c  Check the CFL condition, have processor 0 print out its value,
c  and quit if it is too large.
c
      cfl = k * time_delta / x_delta / x_delta

      if ( id .eq. 0 ) then 
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'UPDATE'
        write ( *, '(a,g14.6)' ) 
     &  '  CFL stability criterion value = ', cfl 
      end if

      if ( 0.5 .le. cfl ) then
        if ( id .eq. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'UPDATE - Warning!'
          write ( *, '(a)' ) '  Computation cancelled!'
          write ( *, '(a)' ) '  CFL condition failed.'
          write ( *, '(a,g14.6)' ) '  0.5 <= K * dT / dX / dX = ', cfl
        end if
        return
      end if
c
c  In single processor mode, write out the values of H.
c
      if ( p .eq. 1 ) then

        open ( unit = h_file, file = 'h_data.txt', status = 'unknown' )

        write ( h_file, '(11f14.6)' ) h(1:n)

      end if
c
c  Compute the values of H at the next time, based on current data.
c
      do j = 1, j_max

        time_new = ( dble (         j - j_min ) * time_max
     &             + dble ( j_max - j         ) * time_min )
     &             / dble ( j_max     - j_min )
c
c  Send H(1) to ID-1.
c
        if ( 0 .lt. id ) then
          tag = 1
          call MPI_Send ( h(1), 1, MPI_DOUBLE_PRECISION, id-1, tag, 
     &      MPI_COMM_WORLD, ierr )
        end if
c
c  Receive H(N+1) from ID+1.
c
        if ( id .lt. p - 1 ) then
          tag = 1
          call MPI_Recv ( h(n+1), 1,  MPI_DOUBLE_PRECISION, id+1, tag, 
     &      MPI_COMM_WORLD, status, ierr )
        end if
c
c  Send H(N) to ID+1.
c
        if ( id .lt. p - 1 ) then
          tag = 2
          call MPI_Send ( h(n), 1, MPI_DOUBLE_PRECISION, id+1, tag, 
     &      MPI_COMM_WORLD, ierr )
        end if
c
c  Receive H(0) from ID-1.
c
        if ( 0 .lt. id ) then
          tag = 2
          call MPI_Recv ( h(0), 1, MPI_DOUBLE_PRECISION, id-1, tag, 
     &      MPI_COMM_WORLD, status, ierr )
        end if
c
c  Update the temperature based on the four point stencil.
c
        do i = 1, n
          h_new(i) = h(i) 
     &      + ( time_delta * k / x_delta / x_delta ) 
     &      * ( h(i-1) - 2.0 * h(i) + h(i+1) ) 
     &      + time_delta * rhs ( x(i), time )
        end do
c
c  H at the extreme left and right boundaries was incorrectly computed
c  using the differential equation.  Replace that calculation by
c  the boundary conditions.
c
        if ( 0 .eq. id ) then
          h_new(1) = boundary_condition ( x(1), time_new )
        end if

        if ( id .eq. p - 1 ) then
          h_new(n) = boundary_condition ( x(n), time_new )
        end if
c
c  Update time and temperature.
c
        time = time_new

        do i = 1, n
          h(i) = h_new(i)
        end do
c
c  In single processor mode, add current solution data to output file.
c
        if ( p .eq. 1 ) then
          write ( h_file, '(11f14.6)' ) h(1:n)
        end if

      end do

      if ( p .eq. 1 ) then
        close ( unit = h_file )
      end if

      return
      end
      function boundary_condition ( x, time )

c*********************************************************************72
c
cc BOUNDARY_CONDITION evaluates the boundary condition of the differential equation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, TIME, the position and time.
c
c    Output, double precision BOUNDARY_CONDITION, the value of the boundary condition.
c
      implicit none

      double precision boundary_condition
      double precision time
      double precision x
c
c  Left condition:
c
      if ( x .lt. 0.5 ) then
        boundary_condition = 100.0 + 10.0 * sin ( time )
      else
        boundary_condition = 75.0
      end if

      return
      end
      function initial_condition ( x, time )

c*********************************************************************72
c
cc INITIAL_CONDITION evaluates the initial condition of the differential equation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, TIME, the position and time.
c
c    Output, double precision INITIAL_CONDITION, the value of the initial condition.
c
      implicit none

      double precision initial_condition
      double precision time
      double precision x

      initial_condition = 95.0

      return
      end
      function rhs ( x, time )

c*********************************************************************72
c
cc RHS evaluates the right hand side of the differential equation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, TIME, the position and time.
c
c    Output, double precision RHS, the value of the right hand side.
c
      double precision rhs
      double precision time
      double precision x

      rhs = 0.0

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
