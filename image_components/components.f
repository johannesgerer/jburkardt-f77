      subroutine file_column_count ( input_filename, column_num )

c*********************************************************************72
c
cc FILE_COLUMN_COUNT counts the number of columns in the first line of a file.
c
c  Discussion:
c
c    The file is assumed to be a simple text file.
c
c    Most lines of the file is presumed to consist of COLUMN_NUM words,
c    separated by spaces.  There may also be some blank lines, and some
c    comment lines,
c    which have a "#" in column 1.
c
c    The routine tries to find the first non-comment non-blank line and
c    counts the number of words in that line.
c
c    If all lines are blanks or comments, it goes back and tries to analyze
c    a comment line.
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
c    Input, character * ( * ) INPUT_FILENAME, the name of the file.
c
c    Output, integer COLUMN_NUM, the number of columns in the file.
c
      implicit none

      integer column_num
      logical got_one
      character * ( * ) input_filename
      integer input_unit
      character * ( 255 ) line
      integer s_len_trim
c
c  Open the file.
c
      call get_unit ( input_unit )

      open ( unit = input_unit, file = input_filename,
     &  status = 'old', form = 'formatted', access = 'sequential' )
c
c  Read one line, but skip blank lines and comment lines.
c
      got_one = .false.

10    continue

        read ( input_unit, '(a)', err = 20 ) line

        if ( s_len_trim ( line ) .eq. 0 ) then
          go to 10
        end if

        if ( line(1:1) .eq. '#' ) then
          go to 10
        end if

        got_one = .true.
        go to 20

      go to 10

20    continue

      if ( .not. got_one ) then

        rewind ( input_unit )

30      continue

          read ( input_unit, '(a)', err = 40 ) line

          if ( s_len_trim ( line ) .eq. 0 ) then
            go to 30
          end if

          got_one = .true.
          go to 40

        go to 30

40    continue

      end if

      close ( unit = input_unit )

      if ( .not. got_one ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FILE_COLUMN_COUNT - Warning.'
        write ( *, '(a)' ) '  The file does not contain any data.'
        column_num = -1
        return
      end if

      call s_word_count ( line, column_num )

      return
      end
      subroutine file_row_count ( input_filename, row_num )

c*********************************************************************72
c
cc FILE_ROW_COUNT counts the number of row records in a file.
c
c  Discussion:
c
c    It does not count lines that are blank, or that begin with a
c    comment symbol '#'.
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
c    Input, character * ( * ) INPUT_FILENAME, the name of the input file.
c
c    Output, integer ROW_NUM, the number of rows found.
c
      implicit none

      integer bad_num
      integer comment_num
      integer ierror
      character * ( * ) input_filename
      integer input_status
      integer input_unit
      character * ( 255 ) line
      integer record_num
      integer row_num
      integer s_len_trim

      call get_unit ( input_unit )

      open ( unit = input_unit, file = input_filename,
     &  status = 'old' )

      comment_num = 0
      row_num = 0
      record_num = 0
      bad_num = 0

10    continue

        read ( input_unit, '(a)', err = 20, end = 20 ) line

        record_num = record_num + 1

        if ( line(1:1) .eq. '#' ) then
          comment_num = comment_num + 1
          go to 10
        end if

        if ( s_len_trim ( line ) .eq. 0 ) then
          comment_num = comment_num + 1
          go to 10
        end if

        row_num = row_num + 1

      go to 10

20    continue

      close ( unit = input_unit )

      return
      end
      subroutine get_unit ( unit )

c*********************************************************************72
c
cc GET_UNIT returns a free FORTRAN unit number.
c
c  Discussion:
c
c    A "free" FORTRAN unit number is a value between 1 and 99 which
c    is not currently associated with an I/O device.  A free FORTRAN unit
c    number is needed in order to open a file with the OPEN command.
c
c    If UNIT = 0, then no free FORTRAN unit could be found, although
c    all 99 units were checked (except for units 5, 6 and 9, which
c    are commonly reserved for console I/O).
c
c    Otherwise, UNIT is a value between 1 and 99, representing a
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
      subroutine i4block_components ( l, m, n, a, component_num, c )

c*********************************************************************72
c
cc I4BLOCK_COMPONENTS assigns contiguous nonzero pixels to a common component.
c
c  Discussion:
c
c    On input, the A array contains values of 0 or 1.
c
c    The 0 pixels are to be ignored.  The 1 pixels are to be grouped
c    into connected components.
c
c    The pixel A(I,J,K) is "connected" to the pixels:
c
c      A(I-1,J,  K  ),  A(I+1,J,  K  ),
c      A(I,  J-1,K  ),  A(I,  J+1,K  ),
c      A(I,  J,  K-1),  A(I,  J,  K+1),
c
c    so most pixels have 6 neighbors.
c
c    On output, COMPONENT_NUM reports the number of components of nonzero
c    data, and the array C contains the component assignment for
c    each nonzero pixel, and is 0 for zero pixels.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 February 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer L, M, N, the order of the array.
c
c    Input, integer A(L,M,N), the pixel array.
c
c    Output, integer COMPONENT_NUM, the number of components
c    of nonzero data.
c
c    Output, integer C(L,M,N), the component array.
c
      implicit none

      integer l
      integer m
      integer n

      integer a(l,m,n)
      integer b
      integer c(l,m,n)
      integer c1
      integer c2
      integer component
      integer component_num
      integer i
      integer j
      integer k
      integer north
      integer p(0:l*m*n)
      integer q(0:l*m*n)
      integer up
      integer west
c
c  Initialization.
c
      do k = 1, n
        do j = 1, m
          do i = 1, l
            c(i,j,k) = 0
          end do
        end do
      end do

      component_num = 0
c
c  P is simply used to store the component labels.  The dimension used
c  here is, of course, usually an absurd overestimate.
c
      do i = 0, l * m * n
        p(i) = i
      end do
c
c  "Read" the array one pixel at a time.  If a (nonzero) pixel's north or
c  west neighbor already has a label, the current pixel inherits it.
c  In case the labels disagree, we need to adjust the P array so we can
c  later deal with the fact that the two labels need to be merged.
c
      do i = 1, l

        do j = 1, m

          do k = 1, n

            if ( i .eq. 1 ) then
              north = 0
            else
              north = c(i-1,j,k)
            end if

            if ( j .eq. 1 ) then
              west = 0
            else
              west = c(i,j-1,k)
            end if

            if ( k .eq. 1 ) then
              up = 0
            else
              up = c(i,j,k-1)
            end if

            if ( a(i,j,k) .ne. 0 ) then
c
c  New component?
c
              if ( north .eq. 0 .and. west .eq. 0 .and. up .eq. 0 ) then
                component_num = component_num + 1
                c(i,j,k) = component_num
c
c  One predecessor is labeled.
c
              else if ( north .ne. 0 .and. west .eq. 0 .and. 
     &                  up .eq. 0 ) then
                c(i,j,k) = north
              else if ( north .eq. 0 .and. west .ne. 0 .and. 
     &                  up .eq. 0 ) then
                c(i,j,k) = west
              else if ( north .eq. 0 .and. west .eq. 0 .and. 
     &                  up .ne. 0 ) then
                c(i,j,k) = up
c
c  Two predecessors are labeled.
c
              else if ( north .eq. 0 .and. west .ne. 0 .and. 
     &          up .ne. 0 ) then
                c(i,j,k) = min ( west, up )
                c1 = min ( p(west), p(up) )
                p(west) = c1
                p(up) = c1
              else if ( north .ne. 0 .and. west .eq. 0 .and. 
     &          up .ne. 0 ) then
                c(i,j,k) = min ( north, up )
                c1 = min ( p(north), p(up) )
                p(north) = c1
                p(up) = c1
              else if ( north .ne. 0 .and. west .ne. 0 .and. 
     &          up .eq. 0 ) then
                c(i,j,k) = min ( north, west )
                c1 = min ( p(north), p(west) )
                p(north) = c1
                p(west) = c1
c
c  Three predecessors are labeled.
c
              else if ( north .ne. 0 .and. west .ne. 0 .and. 
     &          up .ne. 0 ) then
                c(i,j,k) = min ( north, west, up )
                c1 = min ( p(north), p(west), p(up) )
                p(north) = c1
                p(west) = c1
                p(up) = c1
              end if

            end if

          end do

        end do

      end do
c
c  When a component has multiple labels, have the higher labels
c  point to the lowest one.
c
      do component = component_num, 1, -1
        b = component
10      continue
        if ( p(b) .ne. b ) then
          b = p(b)
          go to 10
        end if
        p(component) = b
      end do
c
c  Locate the minimum label for each component.
c  Assign these mininum labels new consecutive indices.
c
      do i = 0, component_num
        q(i) = 0
      end do

      i = 0
      do component = 1, component_num
        if ( p(component) .eq. component ) then
          i = i + 1
          q(component) = i
        end if
      end do

      component_num = i
c
c  Replace the labels by consecutive labels.
c
      do i = 1, l
        do j = 1, m
          do k = 1, n
            c(i,j,k) = q ( p ( c(i,j,k) ) )
          end do
        end do
      end do

      return
      end
      subroutine i4mat_components ( m, n, a, component_num, c )

c*********************************************************************72
c
cc I4MAT_COMPONENTS assigns contiguous nonzero pixels to a common component.
c
c  Discussion:
c
c    On input, the A array contains values of 0 or 1.
c
c    The 0 pixels are to be ignored.  The 1 pixels are to be grouped
c    into connected components.
c
c    The pixel A(I,J) is "connected" to the pixels A(I-1,J), A(I+1,J),
c    A(I,J-1) and A(I,J+1), so most pixels have 4 neighbors.
c
c    (Another choice would be to assume that a pixel was connected
c    to the other 8 pixels in the 3x3 block containing it.)
c
c    On output, COMPONENT_NUM reports the number of components of nonzero
c    data, and the array C contains the component assignment for
c    each nonzero pixel, and is 0 for zero pixels.
c
c  Picture:
c
c    Input A:
c
c      0  2  0  0 17  0  3
c      0  0  3  0  1  0  4
c      1  0  4  8  8  0  7
c      3  0  6 45  0  0  0
c      3 17  0  5  9  2  5
c
c    Output:
c
c      COMPONENT_NUM = 4
c
c      C:
c
c      0  1  0  0  2  0  3
c      0  0  2  0  2  0  3
c      4  0  2  2  2  0  3
c      4  0  2  2  0  0  0
c      4  4  0  2  2  2  2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 February 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the order of the array.
c
c    Input, integer A(M,N), the pixel array.
c
c    Output, integer COMPONENT_NUM, the number of components
c    of nonzero data.
c
c    Output, integer C(M,N), the component array.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer b
      integer c(m,n)
      integer component
      integer component_num
      integer i
      integer j
      integer north
      integer p(0:m*n)
      integer q(0:m*n)
      integer west
c
c  Initialization.
c
      do j = 1, n
        do i = 1, m
          c(i,j) = 0
        end do
      end do

      component_num = 0
c
c  P is simply used to store the component labels.  The dimension used
c  here is, of course, usually an absurd overestimate.
c
      do i = 0, m * n
        p(i) = i
      end do
c
c  "Read" the array one pixel at a time.  If a (nonzero) pixel's north or
c  west neighbor already has a label, the current pixel inherits it.
c  In case the labels disagree, we need to adjust the P array so we can
c  later deal with the fact that the two labels need to be merged.
c
      do i = 1, m

        do j = 1, n

          if ( i .eq. 1 ) then
            north = 0
          else
            north = c(i-1,j)
          end if

          if ( j .eq. 1 ) then
            west = 0
          else
            west = c(i,j-1)
          end if

          if ( a(i,j) .ne. 0 ) then

            if ( north .eq. 0 ) then

              if ( west .eq. 0 ) then
                component_num = component_num + 1
                c(i,j) = component_num
              else
                c(i,j) = west
              end if

            else if ( north .ne. 0 ) then

              if ( west .eq. 0 .or. west .eq. north ) then
                c(i,j) = north
              else
                c(i,j) = min ( north, west )
                if ( north .lt. west ) then
                  p(west) = north
                else
                  p(north) = west
                end if
              end if

            end if

          end if

        end do

      end do
c
c  When a component has multiple labels, have the higher labels
c  point to the lowest one.
c
      do component = component_num, 1, -1
        b = component
10      continue
        if ( p(b) .ne. b ) then
          b = p(b)
          go to 10
        end if
        p(component) = b
      end do
c
c  Locate the minimum label for each component.
c  Assign these mininum labels new consecutive indices.
c
      do i = 0, component_num
        q(i) = 0
      end do

      i = 0
      do component = 1, component_num
        if ( p(component) .eq. component ) then
          i = i + 1
          q(component) = i
        end if
      end do

      component_num = i
c
c  Replace the labels by consecutive labels.
c
      do i = 1, m
        do j = 1, n
          c(i,j) = q ( p ( c(i,j) ) )
        end do
      end do

      return
      end
      subroutine i4mat_data_read ( input_filename, m, n, table )

c*********************************************************************72
c
cc I4MAT_DATA_READ reads data from an I4MAT file.
c
c  Discussion:
c
c    An I4MAT is an array of I4's.
c
c    The file may contain more than N points, but this routine
c    will return after reading N points.
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
c    Input, character * ( * ) INPUT_FILENAME, the name of the input file.
c
c    Input, integer M, the spatial dimension.
c
c    Input, integer N, the number of points.
c
c    Output, integer TABLE(M,N), the data.
c
      implicit none

      integer m
      integer n

      integer i
      integer ierror
      character * ( * ) input_filename
      integer input_status
      integer input_unit
      integer j
      character * ( 255 ) line
      integer table(m,n)
      integer x(m)

      ierror = 0

      call get_unit ( input_unit )

      open ( unit = input_unit, file = input_filename,
     &  status = 'old' )

      j = 0

10    continue

      if ( j .lt. n ) then

        read ( input_unit, '(a)' ) line

        if ( line(1:1) .eq. '#' .or. len_trim ( line ) .eq. 0 ) then
          go to 10
        end if

        call s_to_i4vec ( line, m, x, ierror )

        if ( ierror .ne. 0 ) then
          go to 10
        end if

        j = j + 1

        do i = 1, m
          table(i,j) = x(i)
        end do

        go to 10

      end if

      close ( unit = input_unit )

      return
      end
      subroutine i4mat_header_read ( input_filename, m, n )

c*********************************************************************72
c
cc I4MAT_HEADER_READ reads the header from an integer table file.
c
c  Discussion:
c
c    An I4MAT is an array of I4's.
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
c    Input, character * ( * ) INPUT_FILENAME, the name of the input file.
c
c    Output, integer M, spatial dimension.
c
c    Output, integer N, the number of points.
c
      implicit none

      character * ( * ) input_filename
      integer m
      integer n

      call file_column_count ( input_filename, m )

      if ( m .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4MAT_HEADER_READ - Fatal error!'
        write ( *, '(a)' ) '  There was an I/O problem while'
        write ( *, '(a)' ) '  trying to count the number of data'
        write ( *, '(a,a,a)' ) '  columns in "', input_filename, '".'
        stop
      end if

      call file_row_count ( input_filename, n )

      if ( n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4MAT_HEADER_READ - Fatal error!'
        write ( *, '(a)' ) '  There was an I/O problem while'
        write ( *, '(a)' ) '  trying to count the number of data rows'
        write ( *, '(a,a,a)' ) '  in "', input_filename, '".'
        stop
      end if

      return
      end
      subroutine i4vec_components ( n, a, component_num, c )

c*********************************************************************72
c
cc I4VEC_COMPONENTS assigns contiguous nonzero pixels to a common component.
c
c  Discussion:
c
c    This calculation is trivial compared to the 2D problem, and is included
c    primarily for comparison.
c
c    On input, the A array contains values of 0 or 1.
c
c    The 0 pixels are to be ignored.  The 1 pixels are to be grouped
c    into connected components.
c
c    The pixel A(I) is "connected" to the pixels A(I-1) and A(I+1).
c
c    On output, COMPONENT_NUM reports the number of components of nonzero
c    data, and the array C contains the component assignment for
c    each nonzero pixel, and is 0 for zero pixels.
c
c  Picture:
c
c    Input A:
c
c      0 0 1 2 4 0 0 4 0 0 0 8 9 9 1 2 3 0 0 5 0 1 6 0 0 0 4 0
c
c    Output:
c
c      COMPONENT_NUM = 6
c
c      C:
c
c      0 0 1 1 1 0 0 2 0 0 0 3 3 3 3 3 3 0 0 4 0 5 5 0 0 0 6 0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 February 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the vector.
c
c    Input, integer A(N), the pixel array.
c
c    Output, integer COMPONENT_NUM, the number of components
c    of nonzero data.
c
c    Output, integer C(N), the component array.
c
      implicit none

      integer n

      integer a(n)
      integer c(n)
      integer component_num
      integer j
      integer west
c
c  Initialization.
c
      do j = 1, n
        c(j) = 0
      end do
      component_num = 0
c
c  "Read" the array one pixel at a time.  If a (nonzero) pixel's west neighbor
c  already has a label, the current pixel inherits it.  Otherwise, we have
c  begun a new component.
c
      west = 0

      do j = 1, n

        if ( a(j) .ne. 0 ) then
          if ( west .eq. 0 ) then
            component_num = component_num + 1
          end if
          c(j) = component_num
        end if

        west = c(j)

      end do

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
      subroutine s_to_i4vec ( s, n, ivec, ierror )

c*********************************************************************72
c
cc S_TO_I4VEC reads an I4VEC from a string.
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
c    Input, character * ( * ) S, the string to be read.
c
c    Input, integer N, the number of values expected.
c
c    Output, integer IVEC(N), the values read from the string.
c
c    Output, integer IERROR, error flag.
c    0, no errors occurred.
c    -K, could not read data for entries -K through N.
c
      implicit none

      integer n

      integer i
      integer ierror
      integer ilo
      integer ivec(n)
      integer length
      character * ( * ) s

      i = 0
      ierror = 0
      ilo = 1

10    continue

      if ( i .lt. n ) then

        i = i + 1

        call s_to_i4 ( s(ilo:), ivec(i), ierror, length )

        if ( ierror .ne. 0 ) then
          ierror = -i
          go to 20
        end if

        ilo = ilo + length

      go to 10

      end if

20    continue

      return
      end
      subroutine s_word_count ( s, nword )

c*********************************************************************72
c
cc S_WORD_COUNT counts the number of "words" in a string.
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
c    Input, character * ( * ) S, the string to be examined.
c
c    Output, integer NWORD, the number of "words" in the string.
c    Words are presumed to be separated by one or more blanks.
c
      implicit none

      logical blank
      integer i
      integer lens
      integer nword
      character * ( * ) s

      nword = 0
      lens = len ( s )

      if ( lens .le. 0 ) then
        return
      end if

      blank = .true.

      do i = 1, lens

        if ( s(i:i) .eq. ' ' ) then
          blank = .true.
        else if ( blank ) then
          nword = nword + 1
          blank = .false.
        end if

      end do

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
