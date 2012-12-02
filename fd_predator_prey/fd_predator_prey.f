      program main

c*********************************************************************72
c
cc FD_PREDATOR_PREY solves a pair of predator-prey ODE's.
c
c  Discussion:
c
c    The physical system under consideration is a pair of animal populations.
c
c    The PREY reproduce rapidly; for each animal alive at the beginning of the
c    year, two more will be born by the end of the year.  The prey do not have
c    a natural death rate; instead, they only die by being eaten by the predator.
c    Every prey animal has 1 chance in 1000 of being eaten in a given year by
c    a given predator.
c
c    The PREDATORS only die of starvation, but this happens very quickly.
c    If unfed, a predator will tend to starve in about 1/10 of a year.
c    On the other hand, the predator reproduction rate is dependent on
c    eating prey, and the chances of this depend on the number of available prey.
c
c    The resulting differential equations can be written:
c
c      PREY(0) = 5000
c      PRED(0) =  100
c
c      d PREY / dT =    2 * PREY(T) - 0.001 * PREY(T) * PRED(T)
c      d PRED / dT = - 10 * PRED(T) + 0.002 * PREY(T) * PRED(T)
c
c    Here, the initial values (5000,100) are a somewhat arbitrary starting point.
c
c    The pair of ordinary differential equations that result have an interesting
c    behavior.  For certain choices of the interaction coefficients (such as
c    those given here), the populations of predator and prey will tend to
c    a periodic oscillation.  The two populations will be out of phase; the number
c    of prey will rise, then after a delay, the predators will rise as the prey
c    begins to fall, causing the predator population to crash again.
c
c    In this program, the pair of ODE's is solved with a simple finite difference
c    approximation using a fixed step size.  In general, this is NOT an efficient
c    or reliable way of solving differential equations.  However, this program is
c    intended to illustrate the ideas of finite difference approximation.
c
c    In particular, if we choose a fixed time step size DT, then a derivative
c    such as dPREY/dT is approximated by:
c
c      d PREY / dT = approximately ( PREY(T+DT) - PREY(T) ) / DT
c
c    which means that the first differential equation can be written as
c
c      PREY(T+DT) = PREY(T) + DT * ( 2 * PREY(T) - 0.001 * PREY(T) * PRED(T) ).
c
c    We can rewrite the equation for PRED as well.  Then, since we know the
c    values of PRED and PREY at time 0, we can use these finite difference
c    equations to estimate the values of PRED and PREY at time DT.  These values
c    can be used to get estimates at time 2*DT, and so on.  To get from time
c    T_START = 0 to time T_STOP = 5, we simply take STEP_NUM steps each of size
c    DT = ( T_STOP - T_START ) / STEP_NUM.
c
c    Because finite differences are only an approximation to derivatives, this
c    process only produces estimates of the solution.  And these estimates tend
c    to become more inaccurate for large values of time.  Usually, we can reduce
c    this error by decreasing DT and taking more, smaller time steps.
c
c    In this example, for instance, taking just 100 steps gives nonsensical
c    answers.  Using STEP_NUM = 1000 gives an approximate solution that seems
c    to have the right kind of oscillatory behavior, except that the amplitude
c    of the waves increases with each repetition.  Using 10000 steps, the
c    approximation begins to become accurate enough that we can see that the
c    waves seem to have a fixed period and amplitude.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    George Lindfield, John Penny,
c    Numerical Methods Using MATLAB,
c    Second Edition,
c    Prentice Hall, 1999,
c    ISBN: 0-13-012641-1,
c    LC: QA297.P45.
c
c  Parameters:
c
c    Input, integer STEP_NUM, the number of steps.
c
      implicit none

      integer step_max
      parameter ( step_max = 20000 )

      integer arg_num
      double precision dt
      character * ( 80 ) filename
      integer i
      integer iarg
      integer iargc
      integer ierror
      integer length
      double precision pred_init
      double precision prey_init
      integer step_num
      character ( len = 80 ) string
      double precision trf(3,step_max+1)
      double precision t_start
      double precision t_stop

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FD_PREDATOR_PREY'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  A finite difference approximate solution of a pair'
      write ( *, '(a)' )
     &  '  of ordinary differential equations for a population'
      write ( *, '(a)' ) '  of predators and prey.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  The exact solution shows wave behavior, with a fixed'
      write ( *, '(a)' )
     &  '  period and amplitude.  The finite difference approximation'
      write ( *, '(a)' )
     &  '  can provide a good estimate for this behavior if the'
      write ( *, '(a)' ) '  stepsize DT is small enough.'
c
c  STEP_NUM is an input argument or else read from the user interactively.
c
      arg_num = iargc ( )

      if ( 1 .le. arg_num ) then

        iarg = 1
        call getarg ( iarg, string )
        call s_to_i4 ( string, step_num, ierror, length )

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FD_PREDATOR_PREY:'
        write ( *, '(a)' ) '  Please enter the number of time steps:'

        read ( *, '(a)' ) step_num

      end if

      t_start = 0.0D+00
      t_stop =  5.0D+00
      dt = ( t_stop - t_start ) / dble ( step_num )
      prey_init = 5000.0D+00
      pred_init = 100.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Initial time    = ', t_start
      write ( *, '(a,g14.6)' ) '  Final time      = ', t_stop
      write ( *, '(a,g14.6)' ) '  Initial prey    = ', prey_init
      write ( *, '(a,g14.6)' ) '  Initial pred    = ', pred_init
      write ( *, '(a,i8)' ) '  Number of steps = ', step_num
c
c  TRF(1:3,1:STEP_NUM+1) contains TIME, PREY, and PREDATOR values for each step.
c
      trf(1,1) = t_start
      trf(2,1) = prey_init
      trf(3,1) = pred_init

      do i = 1, step_num
        trf(1,i+1) = trf(1,i) + dt
        trf(2,i+1) = trf(2,i) + dt
     &    * (    2.0D+00 * trf(2,i) - 0.001D+00 * trf(2,i) * trf(3,i) )
        trf(3,i+1) = trf(3,i) + dt
     &    * ( - 10.0D+00 * trf(3,i) + 0.002D+00 * trf(2,i) * trf(3,i) )
      end do
c
c  Write data to a file.
c
      write ( filename, '(a,i8,a)' ) 'trf_', step_num, '.txt'
      call s_blank_delete ( filename )

      call r8mat_write ( filename, 3, step_num + 1, trf )

      write ( *, '(a)' )
     &  '  T, R, F values written to "' // trim ( filename ) // '".'
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FD_PREDATOR_PREY'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

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

        if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

          open ( unit = i, err = 10, status = 'scratch' )
          close ( unit = i )

          unit = i

          return
        end if

10      continue

      end do

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
c    06 March 2006
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

      double precision r8
      double precision r8_epsilon
      double precision r8_test

      r8 = 1.0D+00
      r8_test = 1.0D+00 + ( r8 / 2.0D+00 )

10    continue

      if ( 1.0D+00 .lt. r8_test ) then
        r8 = r8 / 2.0D+00
        r8_test = 1.0D+00 + ( r8 / 2.0D+00 )
        go to 10
      end if

      r8_epsilon = r8

      return
      end
      subroutine r8mat_write ( output_filename, m, n, table )

c*********************************************************************72
c
cc R8MAT_WRITE writes a R8MAT file.
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
c    22 October 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) OUTPUT_FILENAME, the output file name.
c
c    Input, integer M, the spatial dimension.
c
c    Input, integer N, the number of points.
c
c    Input, double precision TABLE(M,N), the data.
c
      implicit none

      integer m
      integer n

      integer j
      character * ( * ) output_filename
      integer output_unit
      character * ( 30 ) string
      double precision table(m,n)
c
c  Open the file.
c
      call get_unit ( output_unit )

      open ( unit = output_unit, file = output_filename,
     &  status = 'replace' )
c
c  Create the format string.
c
      if ( 0 .lt. m .and. 0 .lt. n ) then

        write ( string, '(a1,i8,a1,i8,a1,i8,a1)' )
     &    '(', m, 'g', 24, '.', 16, ')'
c
c  Write the data.
c
        do j = 1, n
          write ( output_unit, string ) table(1:m,j)
        end do

      end if
c
c  Close the file.
c
      close ( unit = output_unit )

      return
      end
      subroutine s_blank_delete ( s )

c*********************************************************************72
c
cc S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
c
c  Discussion:
c
c    All TAB characters are also removed.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 July 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character*(*) S, the string to be transformed.
c
      implicit none

      character ch
      integer get
      integer put
      character*(*) s
      integer s_len_trim
      integer s_length
      character tab

      tab = char ( 9 )

      put = 0
      s_length = s_len_trim ( s )

      do get = 1, s_length

        ch = s(get:get)

        if ( ch .ne. ' ' .and. ch .ne. tab ) then
          put = put + 1
          s(put:put) = ch
        end if

      end do

      s(put+1:s_length) = ' '

      return
      end
      function s_len_trim ( s )

c*********************************************************************72
c
cc S_LEN_TRIM returns the length of a string to the last nonblank.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 March 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) S, a string.
c
c    Output, integer S_LEN_TRIM, the length of the string to the last nonblank.
c
      implicit none

      integer i
      character*(*) s
      integer s_len_trim

      do i = len ( s ), 1, -1

        if ( s(i:i) .ne. ' ' ) then
          s_len_trim = i
          return
        end if

      end do

      s_len_trim = 0

      return
      end
      subroutine s_to_i4 ( s, ival, ierror, length )

c*********************************************************************72
c
cc S_TO_I4 reads an I4 from a string.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) S, a string to be examined.
c
c    Output, integer IVAL, the integer value read from the string.
c    If the string is blank, then IVAL will be returned 0.
c
c    Output, integer IERROR, an error flag.
c    0, no error.
c    1, an error occurred.
c
c    Output, integer LENGTH, the number of characters of S
c    used to make IVAL.
c
      implicit none

      character c
      integer i
      integer ierror
      integer isgn
      integer istate
      integer ival
      integer length
      character * ( * ) s
      integer s_len_trim

      ierror = 0
      istate = 0
      isgn = 1
      ival = 0

      do i = 1, s_len_trim ( s )

        c = s(i:i)
c
c  Haven't read anything.
c
        if ( istate .eq. 0 ) then

          if ( c .eq. ' ' ) then

          else if ( c .eq. '-' ) then
            istate = 1
            isgn = -1
          else if ( c .eq. '+' ) then
            istate = 1
            isgn = + 1
          else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
            istate = 2
            ival = ichar ( c ) - ichar ( '0' )
          else
            ierror = 1
            return
          end if
c
c  Have read the sign, expecting digits.
c
        else if ( istate .eq. 1 ) then

          if ( c .eq. ' ' ) then

          else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
            istate = 2
            ival = ichar ( c ) - ichar ( '0' )
          else
            ierror = 1
            return
          end if
c
c  Have read at least one digit, expecting more.
c
        else if ( istate .eq. 2 ) then

          if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
            ival = 10 * ival + ichar ( c ) - ichar ( '0' )
          else
            ival = isgn * ival
            length = i - 1
            return
          end if

        end if

      end do
c
c  If we read all the characters in the string, see if we're OK.
c
      if ( istate .eq. 2 ) then
        ival = isgn * ival
        length = s_len_trim ( s )
      else
        ierror = 1
        length = 0
      end if

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
