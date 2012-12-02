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
      subroutine getint ( done, ierror, inunit, ival, string )

c*********************************************************************72
c
cc GETINT reads an integer from a file.
c
c  Discussion:
c
c    The file, or at least the part read by GETINT, is assumed to
c    contain nothing but integers.  These integers may be separated
c    by spaces, or appear on separate lines.  Comments, which begin
c    with "#" and extend to the end of the line, may appear anywhere.
c
c    Each time GETINT is called, it tries to read the next integer
c    it can find.  It remembers where it was in the current line
c    of text.
c
c    The user should open a text file on FORTRAN unit INUNIT,
c    set STRING = ' ' and DONE = TRUE.  The GETINT routine will take
c    care of reading in a new STRING as necessary, and extracting
c    as many integers as possible from the line of text before 
c    reading in the next line.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, logical DONE.
c    On input, if this is the first call, or the user has changed
c    STRING, then set DONE = TRUE.
c    On output, if there is no more data to be read from STRING,
c    then DONE is TRUE.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred while trying to read the integer.
c
c    Input, integer INUNIT, the FORTRAN unit from which to read.
c
c    Output, integer IVAL, the integer that was read.
c
c    Input/output, character * ( * ) STRING, the text of the most recently 
c    read line of the file.
c
      implicit none

      logical done
      integer i
      integer ierror
      integer inunit
      integer ios
      integer ival
      integer last
      character * ( * )  string
      character * ( 80 ) word

10    continue

        call word_next_rd ( string, word, done )

        if ( .not. done ) then
          go to 20
        end if

        read ( inunit, '(a)', err = 30 ) string

        i = index ( string, '#' )

        if ( i .ne. 0 ) then
          string(i:) = ' '
        end if

      go to 10

20    continue

      call s_to_i4 ( word, ival, ierror, last )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'GETINT - Fatal error!'
        write ( *, '(a)' ) '  Error converting string to integer.'
        stop
      end if

      return

30    continue

      ierror = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'GETINT - Fatal error!'
      write ( *, '(a)' ) '  Error reading string from file.'
      stop
      end
      subroutine i4mat_histogram ( m, n, a, histo_num, histo_gram )

c*********************************************************************72
c
cc I4MAT_HISTOGRAM computes a histogram of the elements of an I4MAT.
c
c  Discussion:
c
c    An I4MAT is an array of I4's.
c
c    It is assumed that the entries in the vector A are nonnegative.
c    Only values between 0 and HISTO_NUM will be histogrammed.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the order of A.
c
c    Input, integer A(M,N), the array to examine.  
c
c    Input, integer HISTO_NUM, the maximum value for which a
c    histogram entry will be computed.
c
c    Output, integer HISTO_GRAM(0:HISTO_NUM), contains the 
c    number of entries of A with the values of 0 through HISTO_NUM.
c
      implicit none

      integer histo_num
      integer m
      integer n

      integer a(m,n)
      integer histo_gram(0:histo_num)
      integer i
      integer j

      do i = 0, histo_num
        histo_gram(i) = 0
      end do

      do j = 1, n
        do i = 1, m

          if ( 0 .le. a(i,j) .and. a(i,j) .le. histo_num ) then
            histo_gram(a(i,j)) = histo_gram(a(i,j)) + 1
          end if

        end do
      end do

      return
      end
      subroutine news ( m, n, a, e )

c*********************************************************************72
c
cc NEWS demonstrates the NEWS stencil for edge detection.
c
c  Discussion:
c
c    Given a black and white image A, which we regard as an M by N array
c    of pixels, we want to produce an array E of the same shape, which
c    contains information describing the location of edges.
c
c    A simple algorithm for trying to detect edges in an array that
c    represents an image is the NEWS scheme.  For each pixel A(C),
c    we consider its North, East, West, and South pixel neighbors.  The
c    indexing of arrays and images do not correspond, so we will use
c    these directions instead:
c
c             A(N)
c              |
c              |
c      A(W)---A(C)---A(E)
c              |
c              |
c             A(S)
c
c    Entry E(C) of the edge array will be computed by
c
c      E(C) = abs ( A(N) - A(S) ) + abs ( A(E) - A(W) )
c
c    Pixels of A that represent edges will tend to have high values
c    of E, while pixels that are interior to a region of roughly the
c    same shade will tend to have low values.
c
c    Thus, an edge detection scheme would use the NEWS stencil to
c    compute the E array, determine E_MAX, the maximum entry in E, 
c    choose some threshold value E_THRESH, and declare pixel A(I,J) 
c    to be associated with an edge whenever E(I,J) is greater than E_THRESH.
c
c    In this program, we demonstrate the NEWS stencil using a PGM
c    grayscale image of coins.  At the end, we use the edge information
c    to produce a color image in which the edges of the coins have been
c    outlined in red.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns in the image.
c
c    Input, integer A(M,N), the gray scale image data, presumably
c    integers between 0 and 255.
c
c    Output, integer E(M,N), is 1 for each pixel that is part of an
c    edge, and 0 otherwise.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer b(m+2,n+2)
      integer e(m,n)
      integer e_max
      integer i
      integer j
      integer thresh
c
c  For neatness, we add a border of zeros to the image,
c  then fill in the border by copying the nearby original values.
c  This will be our M+2 by N+2 data array B.
c
      do j = 1, n
        do i = 1, m
          b(i+1,j+1) = a(i,j)
        end do
      end do

      do j = 2, n + 1
        b(1,j) = b(2,j)
        b(m+2,j) = b(m+1,j)
      end do

      do i = 2, m + 1
        b(i,  1) = b(i,  2)
        b(i,n+2) = b(i,n+1)
      end do

      b(  1,  1) = ( b(  1,  2) + b(  2,  1) ) / 2
      b(m+2,  1) = ( b(m+2,  2) + b(m+1,  1) ) / 2
      b(  1,n+2) = ( b(  1,n+1) + b(  2,n+2) ) / 2
      b(m+2,n+2) = ( b(m+2,n+1) + b(m+1,n+2) ) / 2
c
c  Apply the NEWS Operator.  We don't process the boundary pixels.
c
c  The picture is:
c
c   |  0 +1  0 |     |  0  0   0 |
c   |  0  0  0 |  +  | -1  0  +1 |
c   |  0 -1  0 |     |  0  0   0 |
c
      do j = 1, n
        do i = 1, m
          e(i,j) = abs ( - b(i,  j+1) + b(i+2,j+1) ) 
     &           + abs ( - b(i+1,j  ) + b(i+1,j+2) )
        end do
      end do
c
c  Remap E so the largest value is 255.
c
      e_max = 0

      do j = 1, n
        do i = 1, m
          e_max = max ( e_max, e(i,j) )
        end do
      end do
c
c  Threshold the data.  Set the threshold to give enough detail 
c  to guess the coin denominations.
c
      thresh = e_max / 5

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'NEWS:'
      write ( *, '(a,i4)' ) '  E_MAX = ', e_max
      write ( *, '(a,i4)' ) '  Using threshold value THRESH = ', thresh

      do j = 1, n
        do i = 1, m
          if ( e(i,j) .lt. thresh ) then
            e(i,j) = 0
          else
            e(i,j) = 1
          end if
        end do
      end do

      return
      end
      subroutine pbma_write ( file_out_name, row_num, col_num, b )

c*********************************************************************72
c
cc PBMA_WRITE writes an ASCII PBM file.
c
c  Example:
c
c    P1
c    # feep.pbma created by PBMA_IO(PBMA_WRITE).
c    24 7
c    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
c    0 1 1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1 1 1 1 0
c    0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 1 0
c    0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 0 0 0 1 1 1 1 0
c    0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0 0 1 0 0 0 0
c    0 1 0 0 0 0 0 1 1 1 1 0 0 1 1 1 1 0 0 1 0 0 0 0
c    0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c 
c  Modified:
c
c    04 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) FILE_OUT_NAME, the name of the file.
c
c    Input, integer ROW_NUM, COL_NUM, the number of rows 
c    and columns of data.
c
c    Input, integer B(ROW_NUM,COL_NUM), the bit value of each 
c    pixel.  These should be 0 or 1.
c
      implicit none

      integer col_num
      integer row_num

      integer b(row_num,col_num)
      character * ( * ) file_out_name
      integer file_out_unit
      integer ios
c
c  Open the file.
c
      call get_unit ( file_out_unit )

      open ( unit = file_out_unit, file = file_out_name, 
     &  status = 'replace', err = 10 )
c
c  Write the header.
c
      call pbma_write_header ( file_out_name, file_out_unit, row_num, 
     &  col_num )
c
c  Write the data.
c
      call pbma_write_data ( file_out_unit, row_num, col_num, b )
c
c  Close the file.
c
      close ( unit = file_out_unit )

      return

10    continue

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PBMA_WRITE - Fatal error!'
      write ( *, '(a)' ) '  Could not open the file.'

      stop
      end
      subroutine pbma_write_data ( file_out_unit, row_num, col_num, b )

c*********************************************************************72
c
cc PBMA_WRITE_DATA writes the data of an ASCII PBM file.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer FILE_OUT_UNIT, the output file unit number.
c
c    Input, integer ROW_NUM, COL_NUM, the number of rows and 
c    columns of data.
c
c    Input, integer B(ROW_NUM,COL_NUM), the bit value of each 
c    pixel.  These should be 0 or 1.
c
      implicit none

      integer col_num
      integer row_num

      integer b(row_num,col_num)
      integer file_out_unit
      integer i
      integer jhi
      integer jlo
c
c  Write the header.
c
      do i = 1, row_num
        do jlo = 1, col_num, 60
          jhi = min ( jlo + 59, col_num )
          write ( file_out_unit, '(60i1)' ) b(i,jlo:jhi)
        end do
      end do

      return
      end
      subroutine pbma_write_header ( file_out_name, file_out_unit, 
     &  row_num, col_num )

c*********************************************************************72
c
cc PBMA_WRITE_HEADER writes the header of an ASCII PBM file.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) FILE_OUT_NAME, the name of the output file.
c
c    Input, integer FILE_OUT_UNIT, the output file unit number.
c
c    Input, integer ROW_NUM, COL_NUM, the number of rows and 
c    columns of data.
c
      implicit none

      character * ( * )  file_out_name
      integer file_out_unit
      character * ( 2 ) magic
      parameter ( magic = 'P1' )
      integer col_num
      integer row_num
c
c  Write the header.
c
      write ( file_out_unit, '(a2)' ) magic
      write ( file_out_unit, '(a)' ) '# ' // trim ( file_out_name ) 
     &  // ' created by PBMA_IO::PBMA_WRITE.F90.'
      write ( file_out_unit, '(i8,2x,i8)' ) col_num, row_num

      return
      end
      subroutine pgma_read_data ( file_in_unit, row_num, col_num, g )

c*********************************************************************72
c
cc PGMA_READ_DATA reads the data in an ASCII PGM file.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer FILE_IN_UNIT, the unit number of the file.
c
c    Input, integer ROW_NUM, COL_NUM, the number of rows and 
c    columns of data.
c
c    Output, integer G(ROW_NUM,COL_NUM), the gray data.
c
      implicit none

      integer col_num
      integer row_num

      logical done
      integer file_in_unit
      integer g(row_num,col_num)
      integer i
      integer ierror
      integer j
      character ( len = 80 ) string

      ierror = 0
      done = .true.
      string = ' '

      do i = 1, row_num
        do j = 1, col_num

          call getint ( done, ierror, file_in_unit, g(i,j), string )

          if ( ierror .ne. 0 ) then
            ierror = 5
            close ( unit = file_in_unit )
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'PGMA_READ_DATA - Fatal error!'
            write ( *, '(a)' ) '  Problem reading G data.'
            stop
          end if

        end do
      end do

      return
      end
      subroutine pgma_read_header ( file_in_unit, row_num, col_num, 
     &  g_max )

c*********************************************************************72
c
cc PGMA_READ_HEADER reads the header of an ASCII PGM file.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer FILE_IN_UNIT, the unit number of the file.
c
c    Output, integer ROW_NUM, COL_NUM, the number of rows and 
c    columns of data.
c
c    Output, integer G_MAX, the maximum gray value.
c
      implicit none

      logical done
      integer file_in_unit
      integer ierror
      integer ios
      character ( 2 )  magic
      integer g_max
      integer col_num
      integer row_num
      logical s_eqi
      character ( len = 80 ) string
c
c  Read the first line of data, which must begin with the magic number.
c
      read ( file_in_unit, '(a)', err = 10 ) magic

      if ( .not. s_eqi ( magic, 'P2' ) ) then
        ierror = 3
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PGMA_READ_HEADER - Fatal error.'
        write ( *, '(a)' ) '  First two bytes not magic number "P2".'
        write ( *, '(a)' ) '  First two bytes are: "' // magic // '".'
        stop
      end if
c
c  Now search for COL_NUM, ROW_NUM, and G_MAX.
c
      done = .true.
      string = ' '

      call getint ( done, ierror, file_in_unit, col_num, string )

      if ( ierror .ne. 0 ) then
        close ( unit = file_in_unit )
        ierror = 4
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PGMA_READ_HEADER - Fatal error!'
        write ( *, '(a)' ) '  Problem reading COL_NUM.'
        stop
      end if

      call getint ( done, ierror, file_in_unit, row_num, string )

      if ( ierror .ne. 0 ) then
        ierror = 4
        close ( unit = file_in_unit )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PGMA_READ_HEADER - Fatal error!'
        write ( *, '(a)' ) '  Problem reading ROW_NUM.'
        stop
      end if

      call getint ( done, ierror, file_in_unit, g_max, string )

      if ( ierror .ne. 0 ) then
        ierror = 4
        close ( unit = file_in_unit )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PGMA_READ_HEADER - Fatal error!'
        write ( *, '(a)' ) '  Problem reading G_MAX.'
        stop
      end if

      return

10    continue

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PGMA_READ_HEADER - Fatal error!'
      write ( *, '(a)' ) '  End or error while reading file.'
      ierror = 2
      stop
      end
      function s_eqi ( s1, s2 )

c*********************************************************************72
c
cc S_EQI is a case insensitive comparison of two strings for equality.
c
c  Example:
c
c    S_EQI ( 'Anjana', 'ANJANA' ) is TRUE.
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
c    Input, character*(*) S1, S2, the strings to compare.
c
c    Output, logical S_EQI, the result of the comparison.
c
      implicit none

      character c1
      character c2
      integer i
      integer lenc
      logical s_eqi
      character*(*) s1
      integer s1_length
      character*(*) s2
      integer s2_length

      s1_length = len ( s1 )
      s2_length = len ( s2 )
      lenc = min ( s1_length, s2_length )

      s_eqi = .false.

      do i = 1, lenc

        c1 = s1(i:i)
        c2 = s2(i:i)
        call ch_cap ( c1 )
        call ch_cap ( c2 )

        if ( c1 .ne. c2 ) then
          return
        end if

      end do

      do i = lenc + 1, s1_length
        if ( s1(i:i) .ne. ' ' ) then
          return
        end if
      end do

      do i = lenc + 1, s2_length
        if ( s2(i:i) .ne. ' ' ) then
          return
        end if
      end do

      s_eqi = .true.

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
      subroutine word_next_rd ( line, word, done )

c*********************************************************************72
c
cc WORD_NEXT_RD "reads" words from a string, one at a time.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) LINE, a string, presumably containing
c    words separated by spaces.
c
c    Output, character * ( * ) WORD.
c    If DONE is FALSE,
c      WORD contains the "next" word read from LINE.
c    Else
c      WORD is blank.
c
c    Input/output, logical DONE.
c    On input, on the first call, or with a fresh value of LINE,
c      set DONE to TRUE.
c    Else
c      leave it at the output value of the previous call.
c    On output, if a new nonblank word was extracted from LINE
c      DONE is FALSE
c    ELSE
c      DONE is TRUE.
c    If DONE is TRUE, then you need to provide a new LINE of data.
c
c  Local Parameters:
c
c    NEXT is the next location in LINE that should be searched.
c
      implicit none

      logical done
      integer ilo
      integer lenl
      character * ( * )  line
      integer next
      character TAB
      parameter ( TAB = char ( 9 ) )
      character * ( * )  word

      save next

      data next / 1 /

      lenl = len_trim ( line )

      if ( done ) then
        next = 1
        done = .false.
      end if
c
c  Beginning at index NEXT, search LINE for the next nonblank.
c
      ilo = next

10    continue
c
c  ...LINE(NEXT:LENL) is blank.  Return with WORD=' ', and DONE=TRUE.
c
        if ( lenl < ilo ) then
          word = ' '
          done = .true.
          next = lenl + 1
          return
        end if
c
c  ...If the current character is blank, skip to the next one.
c
        if ( line(ilo:ilo) .ne. ' ' .and. line(ilo:ilo) .ne. TAB ) then
          go to 20
        end if

        ilo = ilo + 1

      go to 10

20    continue
c
c  To get here, ILO must be the index of the nonblank starting
c  character of the next word.
c
c  Now search for the LAST nonblank character.
c
      next = ilo + 1

30    continue

        if ( lenl .lt. next ) then
          word = line(ilo:next-1)
          return
        end if

        if ( line(next:next) .eq. ' ' .or. 
     &       line(next:next) .eq. TAB ) then
          go to 40
        end if

        next = next + 1

      go to 30

40    continue

      word = line(ilo:next-1)

      return
      end
