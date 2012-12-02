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
      subroutine file_name_inc ( file_name )

c*********************************************************************72
c
cc FILE_NAME_INC generates the next filename in a series.
c
c  Discussion:
c
c    It is assumed that the digits in the name, whether scattered or
c    connected, represent a number that is to be increased by 1 on
c    each call.  If this number is all 9's on input, the output number
c    is all 0's.  Non-numeric letters of the name are unaffected, and
c    if the name contains no digits, then nothing is done.
c
c  Example:
c
c      Input          Output
c      -----          ------
c      a7to11.txt     a7to12.txt
c      a7to99.txt     a8to00.txt
c      a9to99.txt     a0to00.txt
c      cat.txt        cat.txt
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    09 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character ( len = * ) FILE_NAME.
c    On input, a character string to be incremented.
c    On output, the incremented string.
c
      implicit none

      character c
      logical ch_is_digit
      character*(*) file_name
      integer i
      integer lens

      lens = len ( file_name )

      do i = lens, 1, -1

        c = file_name(i:i)

        if ( ch_is_digit ( c ) ) then

          call digit_inc ( c )

          file_name(i:i) = c

          if ( c .ne. '0' ) then
            return
          end if

        end if

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
c    06 March 2006
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

          open ( unit = i, err = 10 )
          close ( unit = i )

          unit = i

          return
        end if

10      continue

      end do

      return
      end
      subroutine number ( xpage, ypage, height, fpn, angle, ndec )

c*********************************************************************72
c
cc NUMBER draws a number at a given location.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real XPAGE, YPAGE, the coordinates of the location at
c    which the number is to be drawn.
c
c    Input, real HEIGHT, the height of the text.
c
c    Input, real FPN, the number to be drawn.
c
c    Input, real ANGLE, the angle at which the number is to be plotted.
c
c    Input, integer NDEC, controls the precision of the number.
c    If NDEC is greater than 0, it specifies the number of decimal
c    digits to be plotted.
c    If NDEC is 0, only the integer part is plotted.
c    If NDEC is negative, only the rounded integer part is plotted,
c    and the magnitude of NDEC determines the strength of the rounding.
c
      implicit none

      real angle
      real fpn
      real height
      integer ndec
      character ( len = 14 ) string
      integer x_ps
      real xpage
      integer y_ps
      real ypage

      x_ps = nint ( xpage * 72.0E+00 )
      y_ps = nint ( ypage * 72.0E+00 )

      write ( 99, '(a)' )        '  /Times-Roman findfont'
      write ( 99, '(2x,f8.4,2x,a)' ) height, 'inch scalefont setfont'
      write ( 99, '(2x,i8,2x,i8,2x,a)' ) x_ps, y_ps, 'moveto'
      if ( 0 < ndec ) then
        write ( string, '(f14.6)' ) fpn
        write ( 99, '(2x,a,a,a)' )  '(', string, ') show' 
      else
        write ( string, '(i6)' ) int ( fpn )
        write ( 99, '(2x,a,a,a)' )  '(', string(1:6), ') show' 
      end if

      return
      end
      subroutine plot ( xpage, ypage, ipen )

c*********************************************************************72
c
cc PLOT moves the plotting pen to a new position.
c
c  Discussion:
c
c    The PLOT command changes the location of the pen.
c
c    If the pen is down, a line will be drawn as the pen moves.
c    If the pen is up, no line will be drawn.
c
c    In the special case where IPEN is 999, then the call to
c    PLOT does not move the pen; instead, it is interpreted as
c    a request to terminate the current plot.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 January 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real XPAGE, YPAGE, the coordinates of the new pen
c    location.
c
c    Input, integer IPEN, determines the pen status.
c    2, the pen is down, and a line is drawn.
c    3, the pen is up, and no line is drawn.
c    -2, the pen is down, a line is drawn, and the location (XPAGE,YPAGE)
c    becomes the new plot origin.
c    -3, the pen is up, no line is drawn, and the location (XPAGE,YPAGE)
c    becomes the new plot origin.
c    999, the current plot is to be terminated.
c
      implicit none

      integer ipen
      integer unit
      integer x_ps
      integer x_ps_save
      real xpage
      integer y_ps
      integer y_ps_save
      real ypage

      save x_ps_save
      save y_ps_save

      data x_ps_save / 0 /
      data y_ps_save / 0 /

      x_ps = nint ( xpage * 72.0E+00 )
      y_ps = nint ( ypage * 72.0E+00 )

      if ( abs ( ipen ) == 2 ) then

        write ( 99, '(2x,i8,2x,i8,2x,a)' ) 
     &    x_ps_save, y_ps_save, 'moveto'
        write ( 99, '(2x,i8,2x,i8,2x,a)' ) 
     &    x_ps, y_ps, 'lineto'
        write ( 99, '(2x,a)' ) 'stroke'

        x_ps_save = x_ps
        y_ps_save = y_ps

      else if ( abs ( ipen ) == 3 ) then

        write ( 99, '(2x,i8,2x,i8,2x,a)' ) 
     &    x_ps, y_ps, 'moveto'

        x_ps_save = x_ps
        y_ps_save = y_ps

      else if ( ipen == 999 ) then

        write ( 99, '(a)' ) '  restore'
        write ( 99, '(a)' ) '  showpage'
        write ( 99, '(a)' ) '%%Trailer'
        write ( 99, '(a)' ) '%%Pages: 1'
        write ( 99, '(a)' ) '%%EOF'
        close ( unit = 99 )

      else 

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PLOT - Fatal error!'
        write ( *, '(a,i8)' ) '  Input value of IPEN = ', ipen
        stop

      end if

      return
      end
      subroutine plots ( i, j, ldev )

c*********************************************************************72
c
cc PLOTS is called to initialize a plot.
c
c  Discussion:
c
c    PLOTS must be the first CALCOMP routine called when creating
c    a plot.
c
c    The corresponding call to finalize a plot is through the
c    "PLOT" command, using a third argument with the flag value
c    of 999:
c
c      call plot ( x, y, 999 )
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
c    Input, integer I, J, arguments that are no longer used.
c
c    Input, integer LDEV, the FORTRAN logical unit number for the plot.
c
      implicit none

      character*(8) date
      character*(20) file_name
      integer i
      integer ierror
      integer j
      integer ldev
      integer length
      integer s_len_trim
      character*(10) time
      integer unit

      save file_name

      data file_name / 'calcomp_000.ps' /

      call file_name_inc ( file_name )

      open ( unit = 99, file = file_name, status = 'replace' )

      write ( 99, '(a)' )        '%!PS-Adobe-1.0'

      write ( 99, '(a)' )        '%%Creator: calcomp.f'
      length = s_len_trim ( file_name )
      write ( 99, '(a,a)' )      '%%Title: ', file_name(1:length)

      call date_and_time ( date, time )

      write ( 99, '(a,a,2x,a)' ) '%%CreationDate: ', date, time
      write ( 99, '(a)' )        '%%Pages: (atend)'
      write ( 99, '(a,4i8)' )    '%%BoundingBox:', 36, 36, 576, 756
      write ( 99, '(a)' )        '%%Document-Fonts: Times-Roman'
      write ( 99, '(a)' )        '%%LanguageLevel: 1'
      write ( 99, '(a)' )        '%%EndComments'
      write ( 99, '(a)' )        '%%BeginProlog'
      write ( 99, '(a)' )        '/inch {72 mul} def'
      write ( 99, '(a)' )        '%%EndProlog'
      write ( 99, '(a)' )        '  /Times-Roman findfont'
      write ( 99, '(a)' )        '  1.00 inch scalefont'
      write ( 99, '(a)' )        '  setfont'
      write ( 99, '(a)' )        '  1 setlinewidth'
      write ( 99, '(a)' )        '%%Page: 1 1'
      write ( 99, '(a)' )        '  save'

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
      subroutine symbol ( xpage, ypage, height, indx, angle, nchar )

c*********************************************************************72
c
cc SYMBOL plots a symbol.
c
c  Discussion:
c
c    This routine is intended to be used as a substitute for the
c    CALCOMP plotter routine SYMBOL, if that routines was being
c    called to display a single plotter symbol.
c
c    Such calls to SYMBOL will have a final argument that is nonpositive.
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
c    Input, real XPAGE, YPAGE, the coordinates of the location at
c    which the number is to be drawn.
c
c    Input, real HEIGHT, the height of the text.
c
c    Input, integer INDX, the index of the symbol.
c
c    Input, real ANGLE, the angle at which the number is to be plotted.
c
c    Input, integer NCHAR, intended to control how the symbol is drawn.
c    (WE IGNORE THIS ARGUMENT).
c    If NCHAR is zero, one character will be plotted.
c    If NCHAR is -1, the pen is up during the move, after which
c    a single character is plotted.
c    If NCHAR is less than -1, the pen is down during the move,
c    after which a single character is plotted.
c
      implicit none

      real angle
      real height
      integer indx
      integer marker_size
      parameter ( marker_size = 5 )
      integer nchar
      real xpage
      integer x_ps
      real ypage
      integer y_ps

      write ( 99, '(a)' ) '  newpath'

      x_ps = nint ( xpage * 72.0E+00 )
      y_ps = nint ( ypage * 72.0E+00 )

      write ( 99, '(2x,i8,2x,i8,2x,i8.2x,a)' ) 
     &  x_ps, y_ps, marker_size, '0 360 arc closepath fill'

      return
      end
      subroutine text ( xpage, ypage, height, string, angle, nchar )

c*********************************************************************72
c
cc TEXT plots text.
c
c  Discussion:
c
c    This routine is intended to be used as a substitute for the
c    CALCOMP plotter routine SYMBOL, if that routines was being
c    called to display text.
c
c    Such calls to SYMBOL will have a positive final argument.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    10 January 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real XPAGE, YPAGE, the coordinates of the location at
c    which the number is to be drawn.
c
c    Input, real HEIGHT, the height of the text.
c
c    Input, character*(*) STRING, contains the characters to plot.
c
c    Input, real ANGLE, the angle at which the number is to be plotted.
c
c    Input, integer NCHAR, the number of characters to display.
c    This argument is actually ignored, and the nonblank length
c    of STRING is used to determine what to display.
c
      implicit none

      real angle
      real height
      integer nchar
      integer s_len_trim
      character*(*) string
      integer string_length
      integer x_ps
      real xpage
      integer y_ps
      real ypage

      x_ps = nint ( xpage * 72.0E+00 )
      y_ps = nint ( ypage * 72.0E+00 )

      string_length = s_len_trim ( string )

      write ( 99, '(a)' )        '  /Times-Roman findfont'
      write ( 99, '(2x,f8.4,2x,a)' ) height, 'inch scalefont setfont'
      write ( 99, '(2x,i8,2x,i8,2x,a)' ) x_ps, y_ps, 'moveto'
      write ( 99, '(2x,a,a,a)' )  '(', string(1:string_length), ') show'  

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
