      program main

c*********************************************************************72
c
cc MAIN is the main program for FEM1D_SAMPLE.
c
c  Discussion:
c
c    FEM1D_SAMPLE reads files defining a 1D FEM representation of data,
c    and writes out two files containing the arguments and values of the
c    finite element function at sample points.
c
c  Usage:
c
c    fem1d_sample fem_prefix sample_prefix
c
c    where 'fem_prefix' is the common prefix for the FEM files:
c
c    * fem_prefix_nodes.txt,    the node coordinates.
c    * fem_prefix_elements.txt, the nodes that make up each element;
c    * fem_prefix_values.txt,   the values defined at each node.
c
c    and 'sample_prefix' is the common prefix for the SAMPLE files.
c    (the node file is input, and the values file is created by the program.)
c
c    * sample_prefix_nodes.txt,  the node coordinates where samples are desired.
c    * sample_prefix_values.txt, the values computed at each sample node.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer fem_element_num_max
      integer fem_element_order_max
      integer fem_node_dim_max
      integer fem_node_num_max
      integer fem_value_dim_max
      integer fem_value_num_max
      integer sample_node_dim_max
      integer sample_node_num_max
      integer sample_value_dim_max
      integer sample_value_num_max

      parameter ( fem_element_num_max = 1000 )
      parameter ( fem_element_order_max = 2 )
      parameter ( fem_node_dim_max = 1 )
      parameter ( fem_node_num_max = 1000 )
      parameter ( fem_value_dim_max = 3 )
      parameter ( fem_value_num_max = 1000 )
      parameter ( sample_node_dim_max = 1 )
      parameter ( sample_node_num_max = 1000 )
      parameter ( sample_value_dim_max = 3 )
      parameter ( sample_value_num_max = 1000 )

      integer 
     &  fem_element_node(fem_element_order_max*fem_element_num_max)
      integer fem_element_num
      integer fem_element_order
      character * ( 255 ) fem_element_filename
      integer fem_node_dim
      character * ( 255 ) fem_node_filename
      integer fem_node_num
      double precision fem_node_x(fem_node_dim_max*fem_node_num_max)
      character * ( 255 ) fem_prefix
      double precision 
     &  fem_value(fem_value_dim_max*fem_value_num_max)
      integer fem_value_dim
      character * ( 255 ) fem_value_filename
      integer fem_value_num
      integer iarg
      integer iargc
      integer ios
      integer num_arg
      character * ( 255 ) sample_prefix
      integer sample_node_dim
      character * ( 255 ) sample_node_filename
      integer sample_node_num
      double precision 
     &  sample_node_x(sample_node_dim_max*sample_node_num_max)
      integer sample_value_dim
      integer sample_value_num
      double precision 
     &  sample_value(sample_value_dim_max*sample_value_num_max)
      character * ( 255 ) sample_value_filename

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM1D_SAMPLE'
      write ( *, '(a)' ) '  FORTRAN77 version.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Read files defining an FEM function of 1 argument.'
      write ( *, '(a)' ) '  Read a file of sample arguments.'
      write ( *, '(a)' ) 
     &  '  Write a file of function values at the arguments.'
c
c  Get the number of command line arguments.
c
      num_arg = iargc ( )

      if ( num_arg .lt. 1 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Enter the FEM file prefix:'
        read ( *, '(a)' ) fem_prefix

      else

        iarg = 1

        call getarg ( iarg, fem_prefix )

      end if

      if ( num_arg .lt. 2 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Enter the sample file prefix:'
        read ( *, '(a)', iostat = ios ) sample_prefix

      else

        iarg = 2

        call getarg ( iarg, sample_prefix )

      end if
c
c  Create the filenames.
c
      fem_node_filename = trim ( fem_prefix ) // '_nodes.txt'
      fem_element_filename = trim ( fem_prefix ) // '_elements.txt'
      fem_value_filename = trim ( fem_prefix ) // '_values.txt'

      sample_node_filename = trim ( sample_prefix ) // '_nodes.txt'
      sample_value_filename = trim ( sample_prefix ) // '_values.txt'
c
c  Read the FEM data.
c
      call dtable_header_read ( fem_node_filename, fem_node_dim, 
     &  fem_node_num )

      if ( fem_node_dim_max * fem_node_num_max .lt. 
     &  fem_node_dim * fem_node_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FEM1D_SAMPLE - Fatal error!'
        write ( *, '(a)' ) '  Not enough storage for FEM_NODE_X.'
        write ( *, '(a,i8)' ) 
     &    '  FEM_NODE_DIM_MAX * FEM_NODE_NUM_MAX = ', 
     &    fem_node_dim_max * fem_node_num_max
        write ( *, '(a,i8)' ) 
     &    '  FEM_NODE_DIM     * FEM_NODE_NUM     = ', 
     &    fem_node_dim * fem_node_num
        stop
      end if

      call dtable_data_read ( fem_node_filename, fem_node_dim, 
     &  fem_node_num, fem_node_x )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  The FEM node dimension is        ', fem_node_dim
      write ( *, '(a,i8)' ) 
     &  '  The FEM node number is           ', fem_node_num

      if ( fem_node_dim .ne. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FEM1D_SAMPLE - Fatal error!'
        write ( *, '(a)' ) '  Spatial dimension of the nodes is not 1.'
        stop
      end if

      call itable_header_read ( fem_element_filename, fem_element_order, 
     &  fem_element_num )

      if ( fem_element_order_max * fem_element_num_max .lt. 
     &  fem_element_order * fem_element_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FEM1D_SAMPLE - Fatal error!'
        write ( *, '(a)' ) '  Not enough storage for FEM_ELEMENT_NODE.'
        write ( *, '(a,i8)' ) 
     &    '  FEM_ELEMENT_ORDER_MAX * FEM_ELEMENT_NUM_MAX = ', 
     &    fem_element_order_max * fem_element_num_max
        write ( *, '(a,i8)' ) 
     &    '  FEM_ELEMENT_DIM     * FEM_ELEMENT_NUM     = ', 
     &    fem_element_order * fem_element_num
        stop
      end if

      call itable_data_read ( fem_element_filename, fem_element_order, 
     &  fem_element_num, fem_element_node )

      write ( *, '(a,i8)' ) 
     &  '  The FEM element order is         ', fem_element_order
      write ( *, '(a,i8)' ) 
     &  '  The FEM element number is        ', fem_element_num

      call dtable_header_read ( fem_value_filename, fem_value_dim, 
     &  fem_value_num )

      write ( *, '(a,i8)' ) 
     &  '  The FEM value order is           ', fem_value_dim
      write ( *, '(a,i8)' ) 
     &  '  the FEM value number is          ', fem_value_num

      if ( fem_value_num .ne. fem_node_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FEM1D_SAMPLE - Fatal error!'
        write ( *, '(a)' ) '  Number of values and nodes differ.'
        stop
      end if

      if ( fem_value_dim_max * fem_value_num_max .lt. 
     &  fem_value_dim * fem_value_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FEM1D_SAMPLE - Fatal error!'
        write ( *, '(a)' ) '  Not enough storage for FEM_VALUE.'
        write ( *, '(a,i8)' ) 
     &    '  FEM_VALUE_DIM_MAX * FEM_VALUE_NUM_MAX = ', 
     &    fem_value_dim_max * fem_value_num_max
        write ( *, '(a,i8)' ) 
     &    '  FEM_VALUE_DIM     * FEM_VALUE_NUM     = ', 
     &    fem_value_dim * fem_value_num
        stop
      end if

      call dtable_data_read ( fem_value_filename, fem_value_dim, 
     &  fem_value_num, fem_value )
c
c  Read the SAMPLE node data.
c
      call dtable_header_read ( sample_node_filename, sample_node_dim, 
     &  sample_node_num )

      if ( sample_node_dim_max * sample_node_num_max .lt. 
     &  sample_node_dim * sample_node_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FEM1D_SAMPLE - Fatal error!'
        write ( *, '(a)' ) '  Not enough storage for SAMPLE_NODE_X.'
        write ( *, '(a,i8)' ) 
     &  '  SAMPLE_NODE_DIM_MAX * SAMPLE_NODE_NUM_MAX = ', 
     &    sample_node_dim_max * sample_node_num_max
        write ( *, '(a,i8)' ) 
     &  '  SAMPLE_NODE_DIM     * SAMPLE_NODE_NUM     = ', 
     &    sample_node_dim * sample_node_num
        stop
      end if

      call dtable_data_read ( sample_node_filename, sample_node_dim, 
     &  sample_node_num, sample_node_x )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  Sample node spatial dimension is ', sample_node_dim
      write ( *, '(a,i8)' ) 
     &  '  Sample node number is            ', sample_node_num

      if ( sample_node_dim .ne. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FEM1D_SAMPLE - Fatal error!'
        write ( *, '(a)' ) 
     &    '  Spatial dimension of the sample nodes is not 1.'
        stop
      end if
c
c  Compute the SAMPLE values.
c
      sample_value_dim = fem_value_dim
      sample_value_num = sample_node_num

      if ( sample_value_dim_max * sample_value_num_max .lt. 
     &  sample_value_dim * sample_value_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FEM1D_SAMPLE - Fatal error!'
        write ( *, '(a)' ) '  Not enough storage for SAMPLE_VALUE.'
        write ( *, '(a,i8)' ) 
     &  '  SAMPLE_VALUE_DIM_MAX * SAMPLE_VALUE_NUM_MAX = ', 
     &    sample_value_dim_max * sample_value_num_max
        write ( *, '(a,i8)' ) 
     &    '  SAMPLE_VALUE_DIM     * SAMPLE_VALUE_NUM     = ', 
     &    sample_value_dim * sample_value_num
        stop
      end if

      call fem1d_evaluate ( fem_node_num, fem_node_x, 
     &  fem_element_order, fem_element_num, fem_value_dim, 
     &  fem_value, sample_node_num, sample_node_x, sample_value )
c
c  Write the sample values.
c
      call dtable_write0 ( sample_value_filename, sample_value_dim, 
     &  sample_value_num, sample_value )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Interpolated FEM data written to "' 
     &  // trim ( sample_value_filename ) // '".'
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM1D_SAMPLE'
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
      subroutine dtable_data_read ( input_file_name, m, n, table )

c*********************************************************************72
c
cc DTABLE_DATA_READ reads data from a DTABLE file.
c
c  Discussion:
c
c    The file may contain more than N points, but this routine will
c    return after reading N of them.
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
c    Input, character * ( * ) INPUT_FILE_NAME, the name of the input file.
c
c    Input, integer M, the spatial dimension.
c
c    Input, integer N, the number of points.
c
c    Output, double precision TABLE(M,N), the table data.
c
      implicit none

      integer m
      integer  n

      integer i
      integer ierror
      character * ( * ) input_file_name
      integer input_status
      integer input_unit
      integer j
      character * ( 255 ) line
      integer s_len_trim
      double precision table(m,n)
      double precision x(m)

      ierror = 0

      call get_unit ( input_unit )

      open ( unit = input_unit, file = input_file_name, 
     &  status = 'old' )

      j = 0

10    continue

      if ( j .lt. n ) then

        read ( input_unit, '(a)' ) line

        if ( line(1:1) == '#' .or. s_len_trim ( line ) .eq. 0 ) then
          go to 10
        end if

        call s_to_r8vec ( line, m, x, ierror )

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
      subroutine dtable_header_read ( input_file_name, m, n )

c*********************************************************************72
c
cc DTABLE_HEADER_READ reads the header from a DTABLE file.
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
c    Input, character * ( * ) INPUT_FILE_NAME, the name of the input file.
c
c    Output, integer M, spatial dimension.
c
c    Output, integer N, the number of points.
c
      implicit none

      character * ( * ) input_file_name
      integer m
      integer n

      call file_column_count ( input_file_name, m )

      if ( m .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DTABLE_HEADER_READ - Fatal errorc'
        write ( *, '(a)' ) '  There was an I/O problem while trying'
        write ( *, '(a)' ) '  to count the number of data columns in'
        write ( *, '(a,a,a)' ) '  the file "', input_file_name, '".'
        stop
      end if

      call file_row_count ( input_file_name, n )

      if ( n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'DTABLE_HEADER_READ - Fatal errorc'
        write ( *, '(a)' ) '  There was an I/O problem while trying'
        write ( *, '(a)' ) '  to count the number of data rows in'
        write ( *, '(a,a,a)' ) '  the file "', input_file_name, '".'
        stop
      end if

      return
      end
      subroutine dtable_write0 ( output_file_name, m, n, table )

c*********************************************************************72
c
cc DTABLE_WRITE0 writes a DTABLE file with no header.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) OUTPUT_FILE_NAME, the output file name.
c
c    Input, integer M, the spatial dimension.
c
c    Input, integer N, the number of points.
c
c    Input, double precision TABLE(M,N), the table data.
c
      implicit none

      integer m
      integer n

      integer j
      character * ( * ) output_file_name
      integer output_unit
      character * ( 30 ) string
      double precision table(m,n)
c
c  Open the file.
c
      call get_unit ( output_unit )

      open ( unit = output_unit, file = output_file_name, 
     &  status = 'replace' )
c
c  Create the format string.
c
      write ( string, '(a1,i8,a1,i8,a1,i8,a1)' ) 
     &  '(', m, 'g', 14, '.', 6, ')'
c
c  Write the data.
c
      do j = 1, n
        write ( output_unit, string ) table(1:m,j)
      end do
c
c  Close the file.
c
      close ( unit = output_unit )

      return
      end
      subroutine fem1d_evaluate ( node_num, node_x, element_order, 
     &  element_num, value_dim, value, sample_node_num, sample_node_x, 
     &  sample_value )

c*********************************************************************72
c
cc FEM1D_EVALUATE evaluates a 1D FEM function at sample points.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NODE_NUM, the number of FEM nodes.
c
c    Input, double precision NODE_X(NODE_NUM), the nodes.  
c
c    Input, integer ELEMENT_ORDER, the element order.
c
c    Input, integer ELEMENT_NUM, the number of elements.
c
c    Input, integer VALUE_DIM, the value dimension.
c
c    Input, double precision VALUE(VALUE_DIM,NODE_NUM), the FEM values.
c
c    Input, integer SAMPLE_NODE_NUM, the number of sample points.
c
c    INput, double precision SAMPLE_NODE_X(SAMPLE_NODE_NUM), the sample nodes.
c
c    Output, double precision SAMPLE_VALUE(VALUE_DIM,SAMPLE_NODE_NUM),
c    the interpolated FEM values at sample nodes.
c
      implicit none

      integer node_num
      integer sample_node_num
      integer value_dim
      integer value_num

      integer element_num
      integer element_order
      integer i
      integer l
      double precision node_x(node_num)
      integer r
      integer sample
      integer sample_left(sample_node_num)
      double precision sample_node_x(sample_node_num)
      double precision sample_value(value_dim,sample_node_num)
      double precision value(value_dim,node_num)
c
c  For each sample point, find NODE_LEFT and NODE_RIGHT that bracket it.
c
      call r8vec_bracket4 ( node_num, node_x, sample_node_num, 
     &  sample_node_x, sample_left )

      if ( element_order == 1 ) then

        do sample = 1, sample_node_num
          do i = 1, value_dim
            sample_value(i,sample) = value(i,sample_left(sample)) 
          end do
        end do
      
      else if ( element_order == 2 ) then

        do sample = 1, sample_node_num

          l = sample_left(sample)
          r = sample_left(sample)+1

          do i = 1, value_dim
            sample_value(i,sample) = 
     &        ( value(i,l) * ( node_x(r) - sample_node_x(sample) )
     &        + value(i,r) * ( sample_node_x(sample) - node_x(l) ) ) 
     &        /              ( node_x(r) - node_x(l) )

          end do

        end do

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FEM_EVALUATE - Fatal errorc'
        write ( *, '(a)' ) '  Cannot handle elements of this order.'
        stop

      end if

      return
      end
      subroutine file_column_count ( input_file_name, column_num )

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
c    Input, character * ( * ) INPUT_FILE_NAME, the name of the file.
c
c    Output, integer COLUMN_NUM, the number of columns in the file.
c
      implicit none

      integer column_num
      logical got_one
      character * ( * ) input_file_name
      integer input_unit
      character * ( 256 ) line
      integer s_len_trim
c
c  Open the file.
c
      call get_unit ( input_unit )

      open ( unit = input_unit, file = input_file_name, 
     &  status = 'old',  form = 'formatted', access = 'sequential' )
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
      subroutine file_row_count ( input_file_name, row_num )

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
c    Input, character * ( * ) INPUT_FILE_NAME, the name of the input file.
c
c    Output, integer ROW_NUM, the number of rows found.
c
      implicit none

      integer bad_num
      integer comment_num
      integer ierror
      character * ( * ) input_file_name
      integer input_status
      integer input_unit
      character * ( 100 ) line
      integer record_num
      integer row_num
      integer s_len_trim

      call get_unit ( input_unit )

      open ( unit = input_unit, file = input_file_name, 
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
      subroutine itable_data_read ( input_file_name, m, n, table )

c*********************************************************************72
c
cc ITABLE_DATA_READ reads data from an ITABLE file.
c
c  Discussion:
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
c    Input, character * ( * ) INPUT_FILE_NAME, the name of the input file.
c
c    Input, integer M, the spatial dimension.
c
c    Input, integer N, the number of points.
c
c    Output, integer TABLE(M,N), the table data.
c
      implicit none

      integer m
      integer n

      integer i
      integer ierror
      character * ( * ) input_file_name
      integer input_status
      integer input_unit
      integer j
      character * ( 255 ) line
      integer table(m,n)
      integer x(m)

      ierror = 0

      call get_unit ( input_unit )

      open ( unit = input_unit, file = input_file_name, 
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
      subroutine itable_header_read ( input_file_name, m, n )

c*********************************************************************72
c
cc ITABLE_HEADER_READ reads the header from an integer table file.
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
c    Input, character * ( * ) INPUT_FILE_NAME, the name of the input file.
c
c    Output, integer M, spatial dimension.
c
c    Output, integer N, the number of points.
c
      implicit none

      character * ( * ) input_file_name
      integer m
      integer n

      call file_column_count ( input_file_name, m )

      if ( m .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'ITABLE_HEADER_READ - Fatal error!'
        write ( *, '(a)' ) '  There was an I/O problem while'
        write ( *, '(a)' ) '  trying to count the number of data'
        write ( *, '(a,a,a)' ) '  columns in "', input_file_name, '".'
        stop
      end if

      call file_row_count ( input_file_name, n )

      if ( n .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'ITABLE_HEADER_READ - Fatal error!'
        write ( *, '(a)' ) '  There was an I/O problem while'
        write ( *, '(a)' ) '  trying to count the number of data rows'
        write ( *, '(a,a,a)' ) '  in "', input_file_name, '".'
        stop
      end if

      return
      end
      subroutine r8vec_bracket4 ( nt, t, ns, s, left )

c*********************************************************************72
c
cc R8VEC_BRACKET4 finds the nearest interval to each of a vector of values.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    The routine always returns the index LEFT of the sorted array
c    T with the property that either
c    *  T is contained in the interval [ T(LEFT), T(LEFT+1) ], or
c    *  T .lt. T(LEFT) = T(1), or
c    *  T > T(LEFT+1) = T(NT).
c
c    The routine is useful for interpolation problems, where
c    the abscissa must be located within an interval of data
c    abscissas for interpolation, or the "nearest" interval
c    to the (extreme) abscissa must be found so that extrapolation
c    can be carried out.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NT, length of the input array.
c
c    Input, double precision T(NT), an array that has been sorted
c    into ascending order.
c
c    Input, integer NS, the number of points to be bracketed.
c
c    Input, double precision S(NS), values to be bracketed by entries of T.
c
c    Output, integer LEFT(NS).
c    LEFT(I) is set so that the interval [ T(LEFT(I)), T(LEFT(I)+1) ]
c    is the closest to S(I); it either contains S(I), or else S(I)
c    lies outside the interval [ T(1), T(NT) ].
c
      implicit none

      integer ns
      integer nt

      integer high
      integer i
      integer left(ns)
      integer low
      integer mid
      double precision s(ns)
      double precision t(nt)
c
c  Check the input data.
c
      if ( nt .lt. 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8VEC_BRACKET4 - Fatal error!'
        write ( *, '(a)' ) '  NT must be at least 2.'
        stop
      end if

      do i = 1, ns

        left(i) = ( nt + 1 ) / 2
c
c  CASE 1: S .lt. T(LEFT):
c  Search for S in [T(I), T(I+1)] for intervals I = 1 to LEFT-1.
c
        if ( s(i) .lt. t(left(i)) ) then

          if ( left(i) .eq. 1 ) then
            go to 100
          else if ( left(i) .eq. 2 ) then
            left(i) = 1
            go to 100
          else if ( t(left(i)-1) .le. s(i) ) then
            left(i) = left(i) - 1
            go to 100
          else if ( s(i) .le. t(2) ) then
            left(i) = 1
            go to 100
          end if
c
c  ...Binary search for S in [T(I), T(I+1)] for intervals I = 2 to LEFT-2.
c
          low = 2
          high = left(i) - 2

10        continue
      
            if ( low .eq. high ) then
              left(i) = low
              go to 20
            end if

            mid = ( low + high + 1 ) / 2

            if ( t(mid) .le. s(i) ) then
              low = mid
            else
              high = mid - 1
            end if

          go to 10

20        continue
c
c  CASE2: T(LEFT+1) .lt. S:
c  Search for S in [T(I),T(I+1)] for intervals I = LEFT+1 to N-1.
c
        else if ( t(left(i)+1) .lt. s(i) ) then

          if ( left(i) .eq. nt - 1 ) then
            go to 100
          else if ( left(i) .eq. nt - 2 ) then
            left(i) = left(i) + 1
            go to 100
          else if ( s(i) .le. t(left(i)+2) ) then
            left(i) = left(i) + 1
            go to 100
          else if ( t(nt-1) .le. s(i) ) then
            left(i) = nt - 1
            go to 100
          end if
c
c  ...Binary search for S in [T(I), T(I+1)] for intervals I = LEFT+2 to NT-2.
c
          low = left(i) + 2
          high = nt - 2

30        continue

            if ( low .eq. high ) then
              left(i) = low
              go to 40
            end if

            mid = ( low + high + 1 ) / 2

            if ( t(mid) .le. s(i) ) then
              low = mid
            else
              high = mid - 1
            end if

          go to 30

40        continue
c
c  CASE3: T(LEFT) .le. S .le. T(LEFT+1):
c  S is in [T(LEFT), T(LEFT+1)].
c
        else

        end if

100     continue

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
      subroutine s_to_r8vec ( s, n, rvec, ierror )

c*********************************************************************72
c
cc S_TO_R8VEC reads an R8VEC from a string.
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
c    Output, double precision RVEC(N), the values read from the string.
c
c    Output, integer IERROR, error flag.
c    0, no errors occurred.
c    -K, could not read data for entries -K through N.
c
      implicit none

      integer  n

      integer i
      integer ierror
      integer ilo
      integer lchar
      double precision rvec(n)
      character * ( * ) s

      i = 0
      ierror = 0
      ilo = 1

10    continue

      if ( i .lt. n ) then

        i = i + 1

        call s_to_r8 ( s(ilo:), rvec(i), ierror, lchar )

        if ( ierror .ne. 0 ) then
          ierror = -i
          go to 20
        end if

        ilo = ilo + lchar

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
