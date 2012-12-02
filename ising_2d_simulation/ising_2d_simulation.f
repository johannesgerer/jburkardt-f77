      program main

c*********************************************************************72
c
cc MAIN is the main program for ISING_2D_SIMULATION.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      integer iterations
      integer m
      integer n
      double precision prob(5)
      integer seed
      integer step
      double precision thresh

      save prob

      data prob /
     &  0.98D+00, 0.85D+00, 0.50D+00, 0.15D+00, 0.02D+00 /

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ISING_2D_SIMULATION'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Monte Carlo simulation of a 2D Ising model.'
c
c  Get input from commandline or user.
c
      call get_input ( m, n, iterations, thresh, seed )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The number of rows is M = ', m
      write ( *, '(a,i8)' ) '  The number of columns is N = ', n
      write ( *, '(a,i8)' ) 
     &  '  The number of iterations taken is ITERATIONS = ', iterations
      write ( *, '(a,g14.6)' ) '  The threshhold THRESH = ', thresh
      write ( *, '(a,i12)' ) '  The seed SEED = ', seed
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  The transition probability table, based on the number of'
      write ( *, '(a)' ) '  neighbors with the same spin.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '      1         2         3         4         5'
      write ( *, '(a)' ) ' '
      write ( *, '(7f10.4)' ) ( prob(i), i = 1, 5 )
c
c  Do the simulation.
c  If I move the C1 array into TRANSITION, I think the current
c  FORTRAN77 will essentially let me dimension it on the fly.
c
      call transition ( m, n, iterations, prob, thresh, seed )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ISING_2D_SIMULATION'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine ch_cap ( ch )

c*********************************************************************72
c
cc CH_CAP capitalizes a single character.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character CH, the character to capitalize.
c
      implicit none

      character ch
      integer itemp

      itemp = ichar ( ch )

      if ( 97 .le. itemp .and. itemp .le. 122 ) then
        ch = char ( itemp - 32 )
      end if

      return
      end
      function ch_eqi ( c1, c2 )

c*********************************************************************72
c
cc CH_EQI is a case insensitive comparison of two characters for equality.
c
c  Example:
c
c    CH_EQI ( 'A', 'a' ) is TRUE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character C1, C2, the characters to compare.
c
c    Output, logical CH_EQI, the result of the comparison.
c
      implicit none

      character c1
      character c1_cap
      character c2
      character c2_cap
      logical ch_eqi

      c1_cap = c1
      c2_cap = c2

      call ch_cap ( c1_cap )
      call ch_cap ( c2_cap )

      if ( c1_cap == c2_cap ) then
        ch_eqi = .true.
      else
        ch_eqi = .false.
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
      subroutine get_input ( m, n, iterations, thresh, seed )

c*********************************************************************72
c
cc GET_INPUT gets input parameters.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 November 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer M, N, the number of rows and columns.
c
c    Output, integer ITERATIONS, the number of iterations.
c
c    Output, double precision THRESH, the threshhold.
c
c    Output, integer SEED, a seed for the random 
c    number generator.
c
      implicit none

      integer arg_num
      integer iarg
      integer iargc
      integer ierror
      integer iterations
      integer last
      integer m
      integer n
      integer seed
      character * ( 80 ) string
      double precision thresh

      arg_num = iargc ( )

      if ( 1 .le. arg_num ) then
        iarg = 1
        call getarg ( iarg, string )
        call s_to_i4 ( string, m, ierror, last )
      else if ( .true. ) then
        m = 10
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Enter M, the number of rows:'
        read ( *, * ) m
      end if

      if ( 2 .le. arg_num ) then
        iarg = 2
        call getarg ( iarg, string )
        call s_to_i4 ( string, n, ierror, last )
      else if ( .true. ) then
        n = 10
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Enter N, the number of columns:'
        read ( *, * ) n
      end if

      if ( 3 .le. arg_num ) then
        iarg = 3
        call getarg ( iarg, string )
        call s_to_i4 ( string, iterations, ierror, last )
      else if ( .true. ) then
        iterations = 15
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Enter the number of iterations to take.'
        read ( *, * ) iterations
      end if

      if ( 4 .le. arg_num ) then
        iarg = 4
        call getarg ( iarg, string )
        call s_to_r8 ( string, thresh, ierror, last )
      else if ( .true. ) then
        thresh = 0.50D+00
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Enter the threshhold.'
        read ( *, * ) thresh
      end if

      if ( 5 .le. arg_num ) then
        iarg = 5
        call getarg ( iarg, string )
        call s_to_i4 ( string, seed, ierror, last )
      else if ( .true. ) then
        seed = 123456789
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Enter the random number seed.'
        read ( *, * ) seed
      end if

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
      subroutine ising_2d_agree ( m, n, c1, c5 )

c*********************************************************************72
c
cc ISING_2D_AGREE returns the number of neighbors agreeing with each cell.
c
c  Discussion:
c
c    The count includes the cell itself, so it is between 1 and 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 November 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of cells in each 
c    spatial dimension.
c
c    Input, integer C1(M,N), an array of 1's and -1's.
c
c    Output, integer C5(M,N), the number of neighbors 
c    that agree.  1, 2, 3, 4, or 5.
c
      implicit none

      integer m
      integer n

      integer c1(m,n)
      integer c5(m,n)
      integer i
      integer i4_wrap
      integer im
      integer ip
      integer j
      integer jm
      integer jp

      do j = 1, n
        jp = i4_wrap ( j + 1, 1, n )
        jm = i4_wrap ( j - 1, 1, n )
        do i = 1, m
          ip = i4_wrap ( i + 1, 1, m )
          im = i4_wrap ( i - 1, 1, m )
          c5(i,j) = c1(i,j) + c1(ip,j) + c1(im,j) + c1(i,jm) + c1(i,jp)
        end do
      end do

      do j = 1, n
        do i = 1, m
          if ( 0 .lt. c1(i,j) ) then
            c5(i,j) = ( 5 + c5(i,j) ) / 2
          else
            c5(i,j) = ( 5 - c5(i,j) ) / 2
          end if
        end do
      end do

      return
      end
      subroutine ising_2d_initialize ( m, n, thresh, seed, c1 )

c*********************************************************************72
c
cc ISING_2D_INITIALIZE initializes the Ising array.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 November 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, double precision THRESH, the threshhold.
c
c    Input/output, integer SEED, a seed for the random 
c    number generator.
c
c    Output, integer C1(M,N), the initial Ising array.
c
      implicit none

      integer m
      integer n

      integer c1(m,n)
      integer i
      integer j
      double precision r(m,n)
      double precision thresh
      integer seed

      call r8mat_uniform_01 ( m, n, seed, r )

      do j = 1, n
        do i = 1, m
          if ( r(i,j) .le. thresh ) then
            c1(i,j) = -1
          else
            c1(i,j) = +1
          end if
        end do
      end do

      return
      end
      subroutine ising_2d_stats ( step, m, n, c1 )

c*********************************************************************72
c
cc ISING_2D_STATS prints information about the current step.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 November 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer STEP, the step number.
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, integer C1(M,N), the current state of the system.
c
      implicit none

      integer m
      integer n

      integer c1(m,n)
      integer i
      integer j
      integer pos_count
      double precision pos_percent
      integer step
      integer neg_count
      double precision neg_percent

      if ( step == 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Step     Positives       Negatives'
        write ( *, '(a)' ) '             #    %          #    %'
        write ( *, '(a)' ) ' '
      end if

      pos_count = 0
      do j = 1, n
        do i = 1, m
          if ( 0 .lt. c1(i,j) ) then
            pos_count = pos_count + 1
          end if
        end do
      end do

      neg_count = m * n - pos_count
      pos_percent = dble ( 100 * pos_count ) / dble ( m * n )
      neg_percent = dble ( 100 * neg_count ) / dble ( m * n )

      write ( *, '(2x,i4,2x,i6,2x,f6.2,2x,i6,2x,f6.2)' ) 
     &  step, pos_count, pos_percent, neg_count, neg_percent 

      return
      end
      subroutine neighbor_2d_stats ( step, m, n, c1, c5 )

c*********************************************************************72
c
cc NEIGHBOR_2D_STATS prints neighbor statistics about the current step.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 November 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer STEP, the step number.
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, integer C1(M,N), the current state of the system.
c
c    Input, integer C5(M,N), the number of agreeable neighbors.
c
      implicit none

      integer m
      integer n

      integer c1(m,n)
      integer c5(m,n)
      integer i
      integer j
      integer k
      integer stats(-5:5)
      integer step

      if ( step .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Step     Neighborhood Charge:'
        write ( *, '(a)' ) 
     &    '           -5    -4    -3    -2    -1' // 
     &    '    +1    +2    +3    +4    +5'
        write ( *, '(a)' ) ' '
      end if

      do i = -5, 5
        stats(i) = 0
      end do

      do j = 1, n
        do i = 1, m
          stats(c5(i,j)) = stats(c5(i,j)) + 1
        end do
      end do
      write ( *, '(2x,i4,10(2x,i4))' ) 
     &  step, ( stats(i), i = -5,-1), ( stats(i), i = 1, 5)
      
      return
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
      subroutine s_to_r8 ( s, dval, ierror, length )

c*********************************************************************72
c
cc S_TO_R8 reads an R8 from a string.
c
c  Discussion:
c
c    The routine will read as many characters as possible until it reaches
c    the end of the string, or encounters a character which cannot be
c    part of the number.
c
c    Legal input is:
c
c       1 blanks,
c       2 '+' or '-' sign,
c       2.5 blanks
c       3 integer part,
c       4 decimal point,
c       5 fraction part,
c       6 'E' or 'e' or 'D' or 'd', exponent marker,
c       7 exponent sign,
c       8 exponent integer part,
c       9 exponent decimal point,
c      10 exponent fraction part,
c      11 blanks,
c      12 final comma or semicolon,
c
c    with most quantities optional.
c
c  Example:
c
c    S                 DVAL
c
c    '1'               1.0
c    '     1   '       1.0
c    '1A'              1.0
c    '12,34,56'        12.0
c    '  34 7'          34.0
c    '-1E2ABCD'        -100.0
c    '-1X2ABCD'        -1.0
c    ' 2E-1'           0.2
c    '23.45'           23.45
c    '-4.2E+2'         -420.0
c    '17d2'            1700.0
c    '-14e-2'         -0.14
c    'e2'              100.0
c    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
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
c    Input, character * ( * ) S, the string containing the
c    data to be read.  Reading will begin at position 1 and
c    terminate at the end of the string, or when no more
c    characters can be read to form a legal real.  Blanks,
c    commas, or other nonnumeric data will, in particular,
c    cause the conversion to halt.
c
c    Output, double precision DVAL, the value read from the string.
c
c    Output, integer IERROR, error flag.
c    0, no errors occurred.
c    1, 2, 6 or 7, the input number was garbled.  The
c    value of IERROR is the last type of input successfully
c    read.  For instance, 1 means initial blanks, 2 means
c    a plus or minus sign, and so on.
c
c    Output, integer LENGTH, the number of characters read
c    to form the number, including any terminating
c    characters such as a trailing comma or blanks.
c
      implicit none

      logical ch_eqi
      character c
      double precision dval
      integer ierror
      integer ihave
      integer isgn
      integer iterm
      integer jbot
      integer jsgn
      integer jtop
      integer length
      integer nchar
      integer ndig
      double precision rbot
      double precision rexp
      double precision rtop
      character * ( * ) s
      integer s_len_trim

      nchar = s_len_trim ( s )

      ierror = 0
      dval = 0.0D+00
      length = -1
      isgn = 1
      rtop = 0
      rbot = 1
      jsgn = 1
      jtop = 0
      jbot = 1
      ihave = 1
      iterm = 0

10    continue

        length = length + 1

        if ( nchar .lt. length+1 ) then
          go to 20
        end if

        c = s(length+1:length+1)
c
c  Blank character.
c
        if ( c .eq. ' ' ) then

          if ( ihave .eq. 2 ) then

          else if ( ihave .eq. 6 .or. ihave .eq. 7 ) then
            iterm = 1
          else if ( 1 .lt. ihave ) then
            ihave = 11
          end if
c
c  Comma.
c
        else if ( c .eq. ',' .or. c .eq. ';' ) then

          if ( ihave .ne. 1 ) then
            iterm = 1
            ihave = 12
            length = length + 1
          end if
c
c  Minus sign.
c
        else if ( c .eq. '-' ) then

          if ( ihave .eq. 1 ) then
            ihave = 2
            isgn = -1
          else if ( ihave .eq. 6 ) then
            ihave = 7
            jsgn = -1
          else
            iterm = 1
          end if
c
c  Plus sign.
c
        else if ( c .eq. '+' ) then

          if ( ihave .eq. 1 ) then
            ihave = 2
          else if ( ihave .eq. 6 ) then
            ihave = 7
          else
            iterm = 1
          end if
c
c  Decimal point.
c
        else if ( c .eq. '.' ) then

          if ( ihave .lt. 4 ) then
            ihave = 4
          else if ( 6 .le. ihave .and. ihave .le. 8 ) then
            ihave = 9
          else
            iterm = 1
          end if
c
c  Scientific notation exponent marker.
c
        else if ( ch_eqi ( c, 'E' ) .or. ch_eqi ( c, 'D' ) ) then

          if ( ihave .lt. 6 ) then
            ihave = 6
          else
            iterm = 1
          end if
c
c  Digit.
c
        else if ( ihave .lt. 11 .and. lle ( '0', c )
     &    .and. lle ( c, '9' ) ) then

          if ( ihave .le. 2 ) then
            ihave = 3
          else if ( ihave .eq. 4 ) then
            ihave = 5
          else if ( ihave .eq. 6 .or. ihave .eq. 7 ) then
            ihave = 8
          else if ( ihave .eq. 9 ) then
            ihave = 10
          end if

          call ch_to_digit ( c, ndig )

          if ( ihave .eq. 3 ) then
            rtop = 10.0D+00 * rtop + dble ( ndig )
          else if ( ihave .eq. 5 ) then
            rtop = 10.0D+00 * rtop + dble ( ndig )
            rbot = 10.0D+00 * rbot
          else if ( ihave .eq. 8 ) then
            jtop = 10 * jtop + ndig
          else if ( ihave .eq. 10 ) then
            jtop = 10 * jtop + ndig
            jbot = 10 * jbot
          end if
c
c  Anything else is regarded as a terminator.
c
        else
          iterm = 1
        end if
c
c  If we haven't seen a terminator, and we haven't examined the
c  entire string, go get the next character.
c
        if ( iterm .eq. 1 ) then
          go to 20
        end if

        go to 10

20    continue
c
c  If we haven't seen a terminator, and we have examined the
c  entire string, then we're done, and LENGTH is equal to NCHAR.
c
      if ( iterm .ne. 1 .and. length+1 .eq. nchar ) then
        length = nchar
      end if
c
c  Number seems to have terminated.  Have we got a legal number?
c  Not if we terminated in states 1, 2, 6 or 7.
c
      if ( ihave .eq. 1 .or. ihave .eq. 2 .or.
     &     ihave .eq. 6 .or. ihave .eq. 7 ) then
        ierror = ihave
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'S_TO_R8 - Serious error!'
        write ( *, '(a)' ) '  Illegal or nonnumeric input:'
        write ( *, '(a,a)' ) '    ', s
        return
      end if
c
c  Number seems OK.  Form it.
c
      if ( jtop .eq. 0 ) then
        rexp = 1.0D+00
      else
        if ( jbot .eq. 1 ) then
          rexp = 10.0D+00 ** ( jsgn * jtop )
        else
          rexp = 10.0D+00 ** ( dble ( jsgn * jtop ) / dble ( jbot ) )
        end if
      end if

      dval = dble ( isgn ) * rexp * rtop / rbot

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
      subroutine transition ( m, n, iterations, prob, thresh, seed )

c*********************************************************************72
c
cc TRANSITION carries out a Monte Carlo simulation of a 3D Ising model.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 November 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, integer ITERATIONS, the number of iterations.
c
c    Input, double precision PROB(1:5).  PROB(I) represents the probability 
c    that the spin of a given cell will be reversed, given that it has I 
c    immediate neighbors (including itself) with spin the same as its own.
c
c    Input, double precision THRESH, the threshhold.
c
c    Input/output, integer SEED, a seed for the random number 
c    generator.
c
      implicit none

      integer m
      integer n

      integer c1(m,n)
      integer c5(m,n)
      integer i
      integer iterations
      integer j
      double precision prob(1:5)
      integer step
      double precision r(m,n)
      integer seed
      double precision thresh

      call ising_2d_initialize ( m, n, thresh, seed, c1 )

      step = 0
      call ising_2d_stats ( step, m, n, c1 )

      do step = 1, iterations
c
c  C5 contains 1 through 5, the number of cells that agree with the center cell.
c
        call ising_2d_agree ( m, n, c1, c5 )

        if ( .false. ) then
          call neighbor_2d_stats ( step, m, n, c1, c5 )
        end if
c
c  Determine the chances of flipping cell (I,J).
c
        call r8mat_uniform_01 ( m, n, seed, r )

        do j = 1, n
          do i = 1, m
            if ( r(i,j) < prob(c5(i,j)) ) then
              c1(i,j) = - c1(i,j)
            end if
          end do
        end do

        call ising_2d_stats ( step, m, n, c1 )

      end do

      return
      end
