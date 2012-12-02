      program main

c*********************************************************************72
c
cc MAIN is the main program for FEM1D_PROJECT.
c
c  Discussion:
c
c    FEM1D_PROJECT reads files defining a sampling of a (scalar or vector)
c    function of 1 argument, and a list of nodes and elements to use for
c    a finite element representation of the data.
c
c    It computes a set of finite element coefficients to be associated with
c    the given finite element mesh, and writes that information to a file
c    so that an FEM representation is formed by the node, element and value 
c    files.
c
c  Usage:
c
c    fem1d_project sample_prefix fem_prefix
c
c    where 'sample_prefix' is the common prefix for the SAMPLE files:
c
c    * sample_prefix_nodes.txt,  the node coordinates where samples were taken,
c    * sample_prefix_values.txt, the sample values.
c
c    and 'fem_prefix' is the common prefix for the FEM files:
c
c    * fem_prefix_nodes.txt,    the node coordinates.
c    * fem_prefix_elements.txt, the nodes that make up each element;
c    * fem_prefix_values.txt,   the values defined at each node.
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
      parameter ( fem_node_dim_max = 3 )
      parameter ( fem_node_num_max = 1000 )
      parameter ( fem_value_dim_max = 3 )
      parameter ( fem_value_num_max = 1000 )
      parameter ( sample_node_dim_max = 3 )
      parameter ( sample_node_num_max = 1000 )
      parameter ( sample_value_dim_max = 3 )
      parameter ( sample_value_num_max = 1000 )

      character * ( 255 ) fem_element_filename
      integer 
     &  fem_element_node(fem_element_order_max*fem_element_num_max)
      integer fem_element_num
      integer fem_element_order
      character * ( 255 ) fem_prefix
      integer fem_node_dim
      character * ( 255 ) fem_node_filename
      integer fem_node_num
      double precision fem_node_x(fem_node_dim_max*fem_node_num_max)
      double precision fem_value(fem_value_dim_max*fem_value_num_max)
      integer fem_value_dim
      integer fem_value_num
      character * ( 255 ) fem_value_filename
      integer iarg
      integer iargc
      integer num_arg
      integer sample_node_dim
      character ( len = 255 ) sample_node_filename
      integer sample_node_num
      double precision 
     &  sample_node_x(sample_node_dim_max*sample_node_num_max)
      character * ( 255 ) sample_prefix
      integer sample_value_dim
      integer sample_value_num
      double precision 
     &  sample_value(sample_value_dim_max*sample_value_num_max)
      character * ( 255 ) sample_value_filename

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM1D_PROJECT'
      write ( *, '(a)' ) '  FORTRAN77 version.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Read files defining a sampling of a function of 1 argument.'
      write ( *, '(a)' ) '  Read files defining a finite element mesh.'
      write ( *, '(a)' ) '  Project the sample data onto the mesh, and'
      write ( *, '(a)' ) '  write a file of FEM coefficient values.'
c
c  Get the number of command line arguments.
c
      num_arg = iargc ( )

      if ( num_arg .lt. 1 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Enter the sample file prefix:'
        read ( *, '(a)' ) sample_prefix

      else

        iarg = 1

        call getarg ( iarg, sample_prefix )

      end if

      if ( num_arg .lt. 2 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Enter the FEM file prefix:'
        read ( *, '(a)' ) fem_prefix

      else

        iarg = 2

        call getarg ( iarg, fem_prefix )

      end if
c
c  Create the filenames.
c
      sample_node_filename = trim ( sample_prefix ) // '_nodes.txt'
      sample_value_filename = trim ( sample_prefix ) // '_values.txt'

      fem_node_filename = trim ( fem_prefix ) // '_nodes.txt'
      fem_element_filename = trim ( fem_prefix ) // '_elements.txt'
      fem_value_filename = trim ( fem_prefix ) // '_values.txt'
c
c  Read the SAMPLE data.
c
      call dtable_header_read ( sample_node_filename, sample_node_dim, 
     &  sample_node_num )

      if ( sample_node_dim_max * sample_node_num_max .lt. 
     &  sample_node_dim * sample_node_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FEM1D_PROJECT - Fatal error!'
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
        write ( *, '(a)' ) 'FEM1D_PROJECT - Fatal error!'
        write ( *, '(a)' ) 
     &    '  Spatial dimension of the sample nodes is not 1.'
        stop
      end if

      call dtable_header_read ( sample_value_filename, 
     &  sample_value_dim, sample_value_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Sample value dimension is        ', 
     &  sample_value_dim
      write ( *, '(a,i8)' ) '  Sample value number is           ', 
     &  sample_value_num

      if ( sample_value_num /= sample_node_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FEM1D_PROJECT - Fatal errorc'
        write ( *, '(a)' ) 
     &    '  Number of sample nodes and values are not equal.'
        stop
      end if

      if ( sample_value_dim_max * sample_value_num_max .lt. 
     &  sample_value_dim * sample_value_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FEM1D_PROJECT - Fatal error!'
        write ( *, '(a)' ) '  Not enough storage for SAMPLE_VALUE_X.'
        write ( *, '(a,i8)' ) 
     &  '  SAMPLE_VALUE_DIM_MAX * SAMPLE_VALUE_NUM_MAX = ', 
     &    sample_value_dim_max * sample_value_num_max
        write ( *, '(a,i8)' ) 
     &    '  SAMPLE_VALUE_DIM     * SAMPLE_VALUE_NUM     = ', 
     &    sample_value_dim * sample_value_num
        stop
      end if

      call dtable_data_read ( sample_value_filename, sample_value_dim, 
     &  sample_value_num, sample_value )
c
c  Read the FEM data.
c
      call dtable_header_read ( fem_node_filename, fem_node_dim, 
     &  fem_node_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The FEM node dimension is        ', 
     &  fem_node_dim
      write ( *, '(a,i8)' ) '  The FEM node number is           ', 
     &  fem_node_num

      if ( fem_node_dim .ne. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FEM1D_PROJECT - Fatal error!'
        write ( *, '(a)' ) '  Spatial dimension of the nodes is not 1.'
        stop
      end if

      if ( fem_node_dim_max * fem_node_num_max .lt. 
     &  fem_node_dim * fem_node_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FEM1D_PROJECT - Fatal error!'
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

      call itable_header_read ( fem_element_filename, fem_element_order, 
     &  fem_element_num )

      write ( *, '(a,i8)' ) '  The FEM element order is         ', 
     &  fem_element_order
      write ( *, '(a,i8)' ) '  The FEM element number is        ', 
     &  fem_element_num

      if ( fem_element_order_max * fem_element_num_max .lt. 
     &  fem_element_order * fem_element_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FEM1D_PROJECT - Fatal error!'
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
c
c  Compute the FEM values.
c
      fem_value_dim = sample_value_dim
      fem_value_num = fem_node_num

      if ( fem_value_dim_max * fem_value_num_max .lt. 
     &  fem_value_dim * fem_value_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FEM1D_PROJECT - Fatal error!'
        write ( *, '(a)' ) '  Not enough storage for FEM_VALUE.'
        write ( *, '(a,i8)' ) 
     &    '  FEM_VALUE_DIM_MAX * FEM_VALUE_NUM_MAX = ', 
     &    fem_value_dim_max * fem_value_num_max
        write ( *, '(a,i8)' ) 
     &    '  FEM_VALUE_DIM     * FEM_VALUE_NUM     = ', 
     &    fem_value_dim * fem_value_num
        stop
      end if

      call fem1d_approximate ( sample_node_num, sample_value_dim, 
     &  sample_node_x, sample_value, fem_node_num, fem_node_x, 
     &  fem_element_order, fem_element_num, fem_value_dim, 
     &  fem_value_num, fem_value )
c
c  Write the FEM values.
c
      call dtable_write0 ( fem_value_filename, fem_value_dim, 
     &  fem_value_num, fem_value )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  FEM value data written to "' 
     &  // trim ( fem_value_filename ) // '".'
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM1D_PROJECT'
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
      subroutine fem1d_approximate ( sample_node_num, sample_value_dim, 
     &  sample_node_x, sample_value, fem_node_num, fem_node_x, 
     &  fem_element_order, fem_element_num, fem_value_dim, 
     &  fem_value_num, fem_value )

c*********************************************************************72
c
cc FEM1D_APPROXIMATE approximates data at sample points with an FEM function.
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
c    Input, integer SAMPLE_NODE_NUM, the number of sample points.
c
c    Input, integer SAMPLE_VALUE_DIM, the value dimension.
c
c    Input, double precision SAMPLE_NODE_X(SAMPLE_NODE_NUM), the sample nodes.
c
c    Input, double precision SAMPLE_VALUE(VALUE_DIM,SAMPLE_NODE_NUM),
c    the values at sample nodes.
c
c    Input, integer FEM_NODE_NUM, the number of FEM nodes.
c
c    Input, double precision FEM_NODE_X(FEM_NODE_NUM), the FEM nodes.  
c
c    Input, integer FEM_ELEMENT_ORDER, the element order.
c
c    Input, integer FEM_ELEMENT_NUM, the number of elements.
c
c    Input, integer FEM_VALUE_DIM, the FEM value dimension.
c
c    Input, integer FEM_VALUE_NUM, the number of FEM values.
c
c    Output, double precision FEM_VALUE(FEM_VALUE_DIM,FEM_VALUE_NUM), the FEM values.
c
      implicit none

      integer fem_node_num
      integer fem_value_dim
      integer fem_value_num
      integer quad_num
      parameter ( quad_num = 2 )
      integer sample_node_num
      integer sample_value_dim
      integer sample_value_num

      double precision a(3,fem_node_num)
      double precision a1
      double precision b(fem_node_num,fem_value_dim)
      double precision b1
      integer dim
      integer fem_element_num
      integer fem_element_order
      double precision fem_node_x(fem_node_num)
      double precision fem_value(fem_value_dim,fem_value_num)
      integer i
      double precision integral
      integer j
      integer l
      integer phi_num
      double precision phi_v(3)
      double precision phi_x(3)
      double precision phil
      double precision phir
      double precision phis
      integer quad
      double precision quad_x(quad_num)
      double precision quad_w(quad_num)
      integer r
      integer sample
      double precision sample_node_x(sample_node_num)
      double precision sample_value(sample_value_dim,sample_node_num)
      double precision svec(sample_node_num)
      double precision wq
      double precision x(fem_value_num,fem_value_dim)
      double precision xl
      double precision xq
      double precision xr

      save quad_x
      save quad_w

      data quad_x /
     &  -0.577350269189625764509148780502D+00, 
     &   0.577350269189625764509148780502D+00 /

      data quad_w / 1.0D+00, 1.0D+00 /
c
c  Set up the matrix A.
c
      do j = 1, fem_node_num
        do i = 1, 3
          a(i,j) = 0.0D+00
        end do
      end do

      do l = 1, fem_node_num - 1

        r = l + 1
        xl = fem_node_x(l)
        xr = fem_node_x(r)

        do quad = 1, quad_num

          xq = ( ( 1.0D+00 - quad_x(quad) ) * xl   
     &         + ( 1.0D+00 + quad_x(quad) ) * xr ) 
     &         /   2.0D+00

          wq = quad_w(quad) * ( xr - xl ) / 2.0D+00

          phil = (      xq - xr ) 
     &         / ( xl      - xr )

          phir = ( xl - xq      )
     &         / ( xl      - xr )

          a(2,l) = a(2,l) + wq * phil * phil
          a(3,l) = a(3,l) + wq * phil * phir

          a(1,r) = a(1,r) + wq * phir * phil
          a(2,r) = a(2,r) + wq * phir * phir

        end do

      end do
c
c  Set up the right hand side b.
c
      do j = 1, fem_value_dim
        do i = 1,fem_node_num
          b(i,j) = 0.0D+00
        end do
      end do

      do i = 1, fem_node_num
        
        if ( i .eq. 1 ) then
          phi_num = 2
          phi_x(1) = fem_node_x(1)
          phi_x(2) = fem_node_x(2)
          phi_v(1) = 1.0D+00
          phi_v(2) = 0.0D+00
        else if ( i .lt. fem_node_num ) then
          phi_num = 3
          phi_x(1) = fem_node_x(i-1)
          phi_x(2) = fem_node_x(i)
          phi_x(3) = fem_node_x(i+1)
          phi_v(1) = 0.0D+00
          phi_v(2) = 1.0D+00
          phi_v(3) = 0.0D+00
        else if ( i .eq. fem_node_num ) then
          phi_num = 2
          phi_x(1) = fem_node_x(fem_node_num-1)
          phi_x(2) = fem_node_x(fem_node_num)
          phi_v(1) = 0.0D+00
          phi_v(2) = 1.0D+00
        end if

        a1 = phi_x(1)
        b1 = phi_x(phi_num)

        do dim = 1, fem_value_dim

          do j = 1, sample_node_num
            svec(j) = sample_value(dim,j)
          end do

          call piecewise_linear_product_quad ( a1, b1, phi_num, 
     &      phi_x, phi_v, sample_node_num, sample_node_x, 
     &      svec, integral )

          b(i,dim) = integral

        end do

      end do
c
c  Solve A * X = B.
c
      call r83_np_fss ( fem_node_num, a, fem_value_dim, b, x )

      do j = 1, fem_value_num
        do i = 1, fem_value_dim
          fem_value(i,j) = x(j,i)
        end do
      end do

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
      subroutine piecewise_linear_product_quad ( a, b, f_num, f_x, f_v, 
     &  g_num, g_x, g_v, quad )

c*********************************************************************72
c
cc PIECEWISE_LINEAR_PRODUCT_QUAD: piecewise linear product integral.
c
c  Discussion:
c
c    We are given two piecewise linear functions F(X) and G(X) and we wish
c    to compute the exact value of the integral
c
c      INTEGRAL = Integral ( A <= X <= B ) F(X) * G(X) dx
c
c    The functions F(X) and G(X) are defined as tables of coordinates X and
c    values V.  A piecewise linear function is evaluated at a point X by 
c    evaluating the interpolant to the data at the endpoints of the interval 
c    containing X.  
c
c    It must be the case that A <= B.
c
c    It must be the case that the node coordinates F_X(*) and G_X(*) are
c    given in ascending order.
c
c    It must be the case that:
c
c      F_X(1) <= A and B <= F_X(F_NUM)
c      G_X(1) <= A and B <= G_X(G_NUM)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 April 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision A, B, the limits of integration.
c
c    Input, integer F_NUM, the number of nodes for F.
c
c    Input, double precision F_X(F_NUM), the node coordinates for F.
c
c    Input, double precision F_V(F_NUM), the nodal values for F.
c
c    Input, integer G_NUM, the number of nodes for G.
c
c    Input, double precision G_X(G_NUM), the node coordinates for G.
c
c    Input, double precision G_V(G_NUM), the nodal values for G.
c
c    Output, double precision QUAD, the integral of F(X) * G(X)
c    from A to B.
c
      implicit none

      integer f_num
      integer g_num

      double precision a
      double precision b
      double precision bit
      integer f_left
      double precision f_v(f_num)
      double precision f_x(f_num)
      double precision f0
      double precision f1
      double precision fl
      double precision fr
      integer g_left
      double precision g_v(g_num)
      double precision g_x(g_num)
      double precision g0
      double precision g1
      double precision gl
      double precision gr
      double precision h0
      double precision h1
      double precision h2
      integer i
      double precision quad
      double precision r8_epsilon
      double precision xl
      double precision xr
      double precision xr_max

      quad = 0.0D+00

      if ( f_x(f_num) .le. a .or. g_x(g_num) .le. a ) then
        return
      end if

      if ( f_num .lt. 2 .or. g_num .lt. 2 ) then
        return
      end if

      xr = a

      f_left = 1
      call r8vec_bracket3 ( f_num, f_x, xr, f_left )
      fr = f_v(f_left) + ( xr - f_x(f_left) ) 
     &  * ( f_v(f_left+1) - f_v(f_left) ) 
     &  / ( f_x(f_left+1) - f_x(f_left) )

      g_left = 1
      call r8vec_bracket3 ( g_num, g_x, xr, g_left )
      gr = g_v(g_left) + ( xr - g_x(g_left) ) 
     &  * ( g_v(g_left+1) - g_v(g_left) ) 
     &  / ( g_x(g_left+1) - g_x(g_left) )

      xr_max = b
      xr_max = min ( xr_max, f_x(f_num) )
      xr_max = min ( xr_max, g_x(g_num) )

10    continue

      if ( xr < xr_max ) then
c
c  Shift right values to left.
c
        xl = xr
        fl = fr
        gl = gr
c
c  Determine the new right values.
c  The hard part is figuring out how to advance XR some, but not too much.
c
        xr = xr_max

        do i = 1, 2
          if ( f_left + i .le. f_num ) then
            if ( xl .lt. f_x(f_left+i) .and. 
     &           f_x(f_left+i) .lt. xr ) then
              xr = f_x(f_left+i)
             go to 20
            end if
          end if
        end do

20      continue

        do i = 1, 2
          if ( g_left + i .le. g_num ) then
            if ( xl .lt. g_x(g_left+i) .and. 
     &           g_x(g_left+i) .lt. xr ) then
              xr = g_x(g_left+i)
              go to 30
            end if
          end if
        end do

30      continue

        call r8vec_bracket3 ( f_num, f_x, xr, f_left )
        fr = f_v(f_left) + ( xr - f_x(f_left) ) 
     &    * ( f_v(f_left+1) - f_v(f_left) ) 
     &    / ( f_x(f_left+1) - f_x(f_left) )

        call r8vec_bracket3 ( g_num, g_x, xr, g_left )
        gr = g_v(g_left) + ( xr - g_x(g_left) ) 
     &    * ( g_v(g_left+1) - g_v(g_left) ) 
     &    / ( g_x(g_left+1) - g_x(g_left) )
c
c  Form the linear polynomials for F(X) and G(X) over [XL,XR],
c  then the product H(X), integrate H(X) and add to the running total.
c
        if ( r8_epsilon ( ) .le. abs ( xr - xl ) ) then

          f1 = fl - fr
          f0 = fr * xl - fl * xr

          g1 = gl - gr
          g0 = gr * xl - gl * xr

          h2 = f1 * g1
          h1 = f1 * g0 + f0 * g1
          h0 = f0 * g0

          h2 = h2 / 3.0D+00
          h1 = h1 / 2.0D+00

          bit = ( ( h2 * xr + h1 ) * xr + h0 ) * xr 
     &        - ( ( h2 * xl + h1 ) * xl + h0 ) * xl

          quad = quad + bit / ( xr - xl ) / ( xr - xl )

        end if

        go to 10

      end if

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
      subroutine r83_np_fss ( n, a, nb, b, x )

c*********************************************************************72
c
cc R83_NP_FSS factors and solves multiple R83 systems.
c
c  Discussion:
c
c    The R83 storage format is used for a tridiagonal matrix.
c    The superdiagonal is stored in entries (1,2:N), the diagonal in
c    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
c    original matrix is "collapsed" vertically into the array.
c
c    This algorithm requires that each diagonal entry be nonzero.
c    It does not use pivoting, and so can fail on systems that
c    are actually nonsingular.
c
c  Example:
c
c    Here is how an R83 matrix of order 5 would be stored:
c
c       *  A12 A23 A34 A45
c      A11 A22 A33 A44 A55
c      A21 A32 A43 A54  *
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the linear system.
c
c    Input/output, double precision A(3,N).
c    On input, the tridiagonal matrix.
c    On output, the data in these vectors has been overwritten
c    by factorization information.
c
!    Input, integer ( kind = 4 ) NB, the number of right hand sides.
!
c    Input, double precision B(N,NB), the right hand sides of the linear system.
c
c    Output, double precision X(N,NB), the solutions of the linear system.
c
      implicit none

      integer n
      integer nb

      double precision a(3,n)
      double precision b(n,nb)
      integer i
      integer j
      double precision x(n,nb)
      double precision xmult
c
c  The diagonal entries can't be zero.
c
      do i = 1, n
        if ( a(2,i) .eq. 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R83_NP_FSS - Fatal error!'
          write ( *, '(a,i8,a)' ) '  A(2,', i, ') = 0.'
          return
        end if
      end do

      do j = 1, nb
        do i = 1, n
          x(i,j) = b(i,j)
        end do
      end do

      do i = 2, n
        xmult = a(3,i-1) / a(2,i-1)
        a(2,i) = a(2,i) - xmult * a(1,i)
        do j = 1, nb
          x(i,j)   = x(i,j) - xmult * x(i-1,j)
        end do
      end do

      do j = 1, nb
        x(n,j) = x(n,j) / a(2,n)
      end do

      do i = n-1, 1, -1
        do j = 1, nb
          x(i,j) = ( x(i,j) - a(1,i+1) * x(i+1,j) ) / a(2,i)
        end do
      end do

      return
      end
      subroutine r8vec_bracket3 ( n, t, tval, left )

c*********************************************************************72
c
cc R8VEC_BRACKET3 finds the interval containing or nearest a given value.
c
c  Discussion:
c
c    An R8VEC is a vector of R8's.
c
c    The routine always returns the index LEFT of the sorted array
c    T with the property that either
c    *  T is contained in the interval [ T(LEFT), T(LEFT+1) ], or
c    *  T .lt. T(LEFT) = T(1), or
c    *  T > T(LEFT+1) = T(N).
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
c    Input, integer N, length of the input array.
c
c    Input, double precision T(N), an array that has been sorted
c    into ascending order.
c
c    Input, double precision TVAL, a value to be bracketed by entries of T.
c
c    Input/output, integer LEFT.
c    On input, if 1 .le. LEFT .le. N-1, LEFT is taken as a suggestion for the
c    interval [ T(LEFT), T(LEFT+1) ] in which TVAL lies.  This interval
c    is searched first, followed by the appropriate interval to the left
c    or right.  After that, a binary search is used.
c    On output, LEFT is set so that the interval [ T(LEFT), T(LEFT+1) ]
c    is the closest to TVAL; it either contains TVAL, or else TVAL
c    lies outside the interval [ T(1), T(N) ].
c
      implicit none

      integer n

      integer high
      integer left
      integer low
      integer mid
      double precision t(n)
      double precision tval
c
c  Check the input data.
c
      if ( n .lt. 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8VEC_BRACKET3 - Fatal error!'
        write ( *, '(a)' ) '  N must be at least 2.'
        stop
      end if
c
c  If LEFT is not between 1 and N-1, set it to the middle value.
c
      if ( left .lt. 1 .or. n - 1 .lt. left ) then
        left = ( n + 1 ) / 2
      end if
c
c  CASE 1: TVAL .lt. T(LEFT):
c  Search for TVAL in [T(I), T(I+1)] for intervals I = 1 to LEFT-1.
c
      if ( tval .lt. t(left) ) then

        if ( left .eq. 1 ) then
          return
        else if ( left .eq. 2 ) then
          left = 1
          return
        else if ( t(left-1) .le. tval ) then
          left = left - 1
          return
        else if ( tval .le. t(2) ) then
          left = 1
          return
        end if
c
c  ...Binary search for TVAL in [T(I), T(I+1)] for intervals I = 2 to LEFT-2.
c
        low = 2
        high = left - 2

10      continue

          if ( low .eq. high ) then
            left = low
            return
          end if

          mid = ( low + high + 1 ) / 2

          if ( t(mid) .le. tval ) then
            low = mid
          else
            high = mid - 1
          end if

        go to 10
c
c  CASE2: T(LEFT+1) .lt. TVAL:
c  Search for TVAL in [T(I),T(I+1)] for intervals I = LEFT+1 to N-1.
c
      else if ( t(left+1) .lt. tval ) then

        if ( left .eq. n - 1 ) then
          return
        else if ( left .eq. n - 2 ) then
          left = left + 1
          return
        else if ( tval .le. t(left+2) ) then
          left = left + 1
          return
        else if ( t(n-1) .le. tval ) then
          left = n - 1
          return
        end if
c
c  ...Binary search for TVAL in [T(I), T(I+1)] for intervals I = LEFT+2 to N-2.
c
        low = left + 2
        high = n - 2

20      continue

          if ( low .eq. high ) then
            left = low
            return
          end if

          mid = ( low + high + 1 ) / 2

          if ( t(mid) .le. tval ) then
            low = mid
          else
            high = mid - 1
          end if

        go to 20
c
c  CASE3: T(LEFT) .le. TVAL .le. T(LEFT+1):
c  T is in [T(LEFT), T(LEFT+1)], as the user said it might be.
c
      else

      end if

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
