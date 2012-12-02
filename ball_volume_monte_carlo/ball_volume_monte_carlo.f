      program main

c*********************************************************************72
c
cc MAIN is the main program for BALL_VOLUME_MONTE_CARLO.
c
c  Discussion:
c
c    DIM_NUM = 6 is a reasonable test.
c
c    N_LOG2_MAX = 25 puts a strain on the system, since we generate that
c    many temporary points at once.  To solve bigger problems, it would
c    be better to compute the new points in batches whose maximum size
c    is limited.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 August 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_max
      parameter ( dim_max = 50 )
      integer n_max
      parameter ( n_max = 1000 )

      integer arg_num
      character * 80 arg_string
      integer dim_num
      double precision estimate
      double precision error
      double precision exact
      double precision fx(n_max)
      integer i
      integer iarg
      integer iargc
      integer ierror
      integer j
      integer last
      integer n
      integer n_done
      integer n_more
      integer n_log2
      integer n_log2_max
      parameter ( n_log2_max = 25 )
      integer n2
      double precision quad
      double precision quad_more
      double precision r8vec_sum
      integer seed
      double precision volume
      double precision x(dim_max,n_max)

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BALL_VOLUME_MONTE_CARLO:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Estimate the volume of the unit sphere'
      write ( *, '(a)' ) '  using a Monte Carlo procedure.'
c
c  Get the number of command line arguments.
c
      arg_num = iargc ( )
c
c  Get the spatial dimension.
c
      if ( 1 .le. arg_num ) then

        iarg = 1
        call getarg ( iarg, arg_string )
        call s_to_i4 ( arg_string, dim_num, ierror, last )

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BALL_VOLUME_MONTE_CARLO:'
        write ( *, '(a)' ) '  Enter the spatial dimension of the sphere'

        read ( *, '(a)' ) dim_num

      end if

      if ( dim_max < dim_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BALL_VOLUME_MONTE_CARLO - Fatal error!'
        write ( *, '(a)' ) '  DIM_NUM is too large.'
        write ( *, '(a,i4)' ) '  Maximum value allowed is ', dim_max
        stop
      end if
c
c  Get the random number seed if it was supplied on the command line.
c
      if ( 2 .le. arg_num ) then

        iarg = 2
        call getarg ( iarg, arg_string )
        call s_to_i4 ( arg_string, seed, ierror, last )

      else

        seed = 123456789
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BALL_VOLUME_MONTE_CARLO:'
        write ( *, '(a)' ) 
     &    '  Using default seed for random number generator.'

      end if
c
c  Report user input.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' )  '  The spatial dimension is  ', dim_num
      write ( *, '(a,i12)' ) '  The random number seed is ', seed
c
c  Begin computation.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &    '    Log(N)         N      Estimate         Error'
      write ( *, '(a)' ) ' '

      quad = 0.0D+00
      volume = 2.0D+00**dim_num

      do n_log2 = 0, n_log2_max

        if ( n_log2 .eq. 0 ) then
          quad = 0.0D+00
          n_more = 1
          n = 0
        else if ( n_log2 .eq. 1 ) then
          n_more = 1
        else
          n_more = 2 * n_more
        end if
c
c  Evaluate N_MORE random points, in batches of 1000.
c
        n_done = 0
        quad_more = 0.0D+00

10      continue

        if ( n_done .lt. n_more ) then
          n2 = min ( 1000, n_more - n_done )
          call r8mat_uniform_01 ( dim_num, n2, seed, x )
          do j = 1, n2
            do i = 1, dim_num
              x(i,j) = 2.0D+00 * x(i,j) - 1.0D+00
            end do
          end do
          call sphere_indicator ( dim_num, n2, x, fx )
          quad_more = quad_more + r8vec_sum ( n2, fx )
          n_done = n_done + n2
          go to 10
        end if
c
c  Incorporate the new data into the totals.
c
        n = n + n_more
        quad = quad + quad_more

        estimate = volume * quad / real ( n, kind = 8 )
        call sphere_unit_volume_nd ( dim_num, exact )
        error = abs ( exact - estimate )
        write ( *, '(2x,i8,2x,i8,2x,g16.8,2x,g10.2)' ) 
     &    n_log2, n, estimate, error

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(8x,a2,8x,a2,2x,g16.8,2x,g10.2)' ) 
     &    'oo', 'oo', exact, 0.0D+00
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BALL_VOLUME_MONTE_CARLO:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine r8mat_uniform_01 ( m, n, seed, r )

c*********************************************************************72
c
cc R8MAT_UNIFORM_01 returns a unit pseudorandom R8MAT.
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
c    Output, double precision R(M,N), the array of pseudorandom values.
c
      implicit none

      integer m
      integer n

      integer i
      integer j
      integer k
      integer seed
      double precision r(m,n)

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8MAT_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do j = 1, n

        do i = 1, m

          k = seed / 127773

          seed = 16807 * ( seed - k * 127773 ) - k * 2836

          if ( seed .lt. 0 ) then
            seed = seed + 2147483647
          end if

          r(i,j) = dble ( seed ) * 4.656612875D-10

        end do
      end do

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
      subroutine sphere_indicator ( dim_num, point_num, x, fx )

c*********************************************************************72
c
cc SPHERE_INDICATOR evaluates the unit sphere indicator function.
c
c  Discussion:
c
c    F(X) = 1 if X is on or inside the unit sphere, and 0 elsewhere.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    19 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the spatial dimension.
c
c    Input, integer POINT_NUM, the number of points to evaluate.
c
c    Input, double precision X(DIM_NUM,POINT_NUM), the points.
c
c    Output, double precision FX(POINT_NUM), the unit sphere indicator 
c    function value.
c
      implicit none

      integer dim_num
      integer point_num

      double precision fx(point_num)
      integer j
      double precision x(dim_num,point_num)

      do j = 1, point_num
        if ( sum ( x(1:dim_num,j)**2 ) .le. 1.0D+00 ) then
          fx(j) = 1.0D+00
        else
          fx(j) = 0.0D+00
        end if
      end do

      return
      end
      subroutine sphere_unit_volume_nd ( dim_num, volume )

c*********************************************************************72
c
cc SPHERE_UNIT_VOLUME_ND computes the volume of a unit sphere in M-dimensions.
c
c  Discussion:
c
c    DIM_NUM  Volume
c
c    2             PI
c    3  (4/3)    * PI
c    4  (1/2)    * PI^2
c    5  (8/15)   * PI^2
c    6  (1/6)    * PI^3
c    7  (16/105) * PI^3
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 August 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer DIM_NUM, the dimension of the space.
c
c    Output, double precision VOLUME, the volume of the sphere.
c
      implicit none

      integer dim_num
      integer i
      integer m
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision volume

      volume = 1.0D+00

      if ( mod ( dim_num, 2 ) .eq. 0 ) then
        m = dim_num / 2
        do i = 1, m
          volume = volume * pi / dble ( i )
        end do
      else
        m = ( dim_num - 1 ) / 2
        do i = 1, m
          volume = volume * pi * 2.0D+00
        end do
        do i = m + 1, 2 * m + 1
          volume = volume * 2.0D+00 / dble ( i )
        end do
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
