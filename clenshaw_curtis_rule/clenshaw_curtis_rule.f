      program main

c*********************************************************************72
c
cc MAIN is the main program for CLENSHAW_CURTIS_RULE.
c
c  Discussion:
c
c    This program computes a Clenshaw Curtis quadrature rule
c    and writes it to a file.
c
c    The user specifies:
c    * the ORDER (number of points) in the rule
c    * the root name of the output files.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      integer arg_num
      double precision b
      character * ( 255 ) filename
      integer iarg
      integer iargc
      integer ierror
      integer last
      integer order
      character ( len = 255 ) string

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CLENSHAW_CURTIS_RULE'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Compute a Clenshaw Curtis rule for approximating'
      write ( *, '(a)' ) '    Integral ( A <= x <= B ) f(x) dx'
      write ( *, '(a)' ) '  of order ORDER.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  The user specifies ORDER, A, B and FILENAME.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  ORDER is the number of points;'
      write ( *, '(a)' ) '  A is the left endpoint;'
      write ( *, '(a)' ) '  B is the right endpoint;'
      write ( *, '(a)' ) '  FILENAME is used to generate 3 files:'
      write ( *, '(a)' ) '    filename_w.txt - the weight file'
      write ( *, '(a)' ) '    filename_x.txt - the abscissa file.'
      write ( *, '(a)' ) '    filename_r.txt - the region file.'
c
c  Get the number of command line arguments.
c
      arg_num = iargc ( )
c
c  Get ORDER.
c
      if ( 1 .le. arg_num ) then
        iarg = 1
        call getarg ( iarg, string )
        call s_to_i4 ( string, order, ierror, last )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Enter the rule order ORDER:'
        read ( *, * ) order
      end if
c
c  Get A.
c
      if ( 2 .le. arg_num ) then
        iarg = 2
        call getarg ( iarg, string )
        call s_to_r8 ( string, a, ierror, last )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Enter the left endpoint A:'
        read ( *, * ) a
      end if
c
c  Get B.
c
      if ( 3 .le. arg_num ) then
        iarg = 3
        call getarg ( iarg, string )
        call s_to_r8 ( string, b, ierror, last )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Enter the right endpoint B:'
        read ( *, * ) b
      end if
c
c  Get FILENAME.
c
      if ( 4 .le. arg_num ) then
        iarg = 4
        call getarg ( iarg, filename )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    '  Enter FILENAME, the "root name" of the quadrature files).'
        read ( *, '(a)' ) filename
      end if
c
c  Input summary.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  ORDER = ', order
      write ( *, '(a,g14.6)' ) '  A = ', a
      write ( *, '(a,g14.6)' ) '  B = ', b
      write ( *, '(a)' ) '  FILENAME = "' // trim ( filename ) // '".'
c
c  Create the rule and write it out.
c
      call clenshaw_curtis_handle ( order, a, b, filename )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CLENSHAW_CURTIS_RULE:'
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
      subroutine clenshaw_curtis_compute ( order, x, w )

c*********************************************************************72
c
cc CLENSHAW_CURTIS_COMPUTE computes a Clenshaw Curtis quadrature rule.
c
c  Discussion:
c
c    Our convention is that the abscissas are numbered from left to right.
c
c    The rule is defined on [-1,1].
c
c    The integral to approximate:
c
c      Integral ( -1 <= X <= 1 ) F(X) dX
c
c    The quadrature rule:
c
c      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 February 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer ORDER, the order of the rule.
c    1 <= ORDER.
c
c    Output, double precision X(ORDER), the abscissas.
c
c    Output, double precision W(ORDER), the weights.
c
      implicit none

      integer order

      double precision b
      integer i
      integer j
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision theta
      double precision w(order)
      double precision x(order)

      if ( order .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CLENSHAW_CURTIS_COMPUTE - Fatal error!'
        write ( *, '(a,i8)' ) '  Illegal value of ORDER = ', order
        stop
      end if

      if ( order .eq. 1 ) then
        x(1) = 0.0D+00
        w(1) = 2.0D+00
        return
      end if

      do i = 1, order
        x(i) = cos ( dble ( order - i ) * pi 
     &             / dble ( order - 1 ) )
      end do

      x(1) = -1.0D+00
      if ( mod ( order, 2 ) .eq. 1 ) then
        x((order+1)/2) = 0.0D+00
      end if
      x(order) = +1.0D+00

      do i = 1, order

        theta = dble ( i - 1 ) * pi 
     &        / dble ( order - 1 )

        w(i) = 1.0D+00

        do j = 1, ( order - 1 ) / 2

          if ( 2 * j .eq. ( order - 1 ) ) then
            b = 1.0D+00
          else
            b = 2.0D+00
          end if

          w(i) = w(i) - b * cos ( 2.0D+00 * dble ( j ) * theta ) 
     &         / dble ( 4 * j * j - 1 )

        end do

      end do

      w(1)         =           w(1)         / dble ( order - 1 )
      do i = 2, order - 1
        w(i) = 2.0D+00 * w(i) / dble ( order - 1 )
      end do
      w(order)     =           w(order)     / dble ( order - 1 )

      return
      end
      subroutine clenshaw_curtis_handle ( order, a, b, filename )

c*********************************************************************72
c
cc CLENSHAW_CURTIS_HANDLE constructs a Clenshaw Curtis quadrature rule.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    30 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer ORDER, the order of the rule.
c
c    Input, double precision A, B, the endpoints of the interval.
c
c    Input, character * ( * ) FILENAME, the root name for the files.
c
      implicit none

      double precision a
      double precision b
      character * ( * ) filename
      integer order
      double precision r(2)
      double precision w(order)
      double precision x(order)
c
c  Construct the rule.
c
      r(1) = a
      r(2) = b

      call clenshaw_curtis_compute ( order, x, w )
c
c  Rescale the rule.
c
      call rescale ( a, b, order, x, w )
c
c  Write the rule.
c
      call rule_write ( order, x, w, r, filename )

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
      subroutine rescale ( a, b, n, x, w )

c*********************************************************************72
c
cc RESCALE rescales a quadrature rule from [-1,+1] to [A,B].
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
c    John Burkardt.
c
c  Parameters:
c
c    Input, double precision A, B, the endpoints of the new interval.
c
c    Input, integer N, the order.
c
c    Input/output, double precision X(N), on input, the abscissas for [-1,+1].
c    On output, the abscissas for [A,B].
c
c    Input/output, double precision W(N), on input, the weights for [-1,+1].
c    On output, the weights for [A,B].
c
      implicit none

      integer n

      double precision a
      double precision b
      integer i
      double precision w(n)
      double precision x(n)

      do i = 1, n
        x(i) = ( ( a + b ) + ( b - a ) * x(i) ) / 2.0D+00
      end do

      do i = 1, n
        w(i) = ( b - a ) * w(i) / 2.0D+00
      end do

      return
      end
      subroutine rule_write ( order, x, w, r, filename )

c*********************************************************************72
c
cc RULE_WRITE writes a quadrature rule to a file.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    18 February 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer ORDER, the order of the rule.
c
c    Input, double precision X(ORDER), the abscissas.
c
c    Input, double precision W(ORDER), the weights.
c
c    Input, double precision R(2), defines the region.
c
c    Input, character ( len = * ) FILENAME, specifies the output.
c    'filename_w.txt', 'filename_x.txt', 'filename_r.txt' defining weights,
c    abscissas, and region.
c 
      implicit none

      integer order

      character * ( * ) filename
      character * ( 255 ) filename_r
      character * ( 255 ) filename_w
      character * ( 255 ) filename_x
      integer i
      integer kind
      double precision r(2)
      double precision w(order)
      double precision x(order)

      filename_w = trim ( filename ) // '_w.txt'
      filename_x = trim ( filename ) // '_x.txt'
      filename_r = trim ( filename ) // '_r.txt'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Creating quadrature files.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  "Root" file name is   "' // trim ( filename ) // '".'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Weight file will be   "' // trim ( filename_w ) // '".'
      write ( *, '(a)' ) 
     &  '  Abscissa file will be "' // trim ( filename_x ) // '".'
      write ( *, '(a)' ) 
     &  '  Region file will be   "' // trim ( filename_r ) // '".'
                
      call r8mat_write ( filename_w, 1, order, w )
      call r8mat_write ( filename_x, 1, order, x )
      call r8mat_write ( filename_r, 1, 2,     r )

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
