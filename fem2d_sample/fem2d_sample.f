      program main

c*****************************************************************************80
c
cc MAIN is the main program for FEM2D_SAMPLE.
c
c  Discussion:
c
c    FEM2D_SAMPLE reads files defining a 2D FEM representation of data,
c    and a set of sample points, and writes out a file containing the 
c    value of the finite element function at the sample points.
c
c  Usage:
c
c    fem2d_sample fem_prefix sample_prefix
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
      integer fem_node_num_max
      integer fem_value_dim_max
      integer fem_value_num_max
      integer sample_node_num_max
      integer sample_value_dim_max
      integer sample_value_num_max

      parameter ( fem_element_num_max = 1000 )
      parameter ( fem_element_order_max = 6 )
      parameter ( fem_node_num_max = 1000 )
      parameter ( fem_value_dim_max = 3 )
      parameter ( fem_value_num_max = 1000 )
      parameter ( sample_node_num_max = 1000 )
      parameter ( sample_value_dim_max = 3 )
      parameter ( sample_value_num_max = 1000 )

      integer fem_element
      character ( len = 255 ) fem_element_filename
      integer fem_element_neighbor(3*fem_element_num_max)
      integer 
     &  fem_element_node(fem_element_order_max*fem_element_num_max)
      integer fem_element_num
      integer fem_element_order
      integer fem_node_dim
      integer fem_node_num
      character ( len = 255 ) fem_node_filename
      double precision fem_node_xy(2*fem_node_num_max)
      character ( len = 255 ) fem_prefix
      double precision fem_value(fem_value_dim_max*fem_value_num_max)
      integer fem_value_dim
      character ( len = 255 ) fem_value_filename
      integer fem_value_num
      integer i
      integer iarg
      integer iargc
      integer ios
      integer j
      integer num_arg
      integer sample_node_dim
      character ( len = 255 ) sample_node_filename
      integer sample_node_num
      double precision sample_node_xy(2*sample_node_num_max)
      character ( len = 255 ) sample_prefix
      double precision 
     &  sample_value(sample_value_dim_max*sample_value_num_max)
      integer sample_value_dim
      character ( len = 255 ) sample_value_filename
      integer sample_value_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM2D_SAMPLE'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Read files defining an FEM function of 2 arguments.'
      write ( *, '(a)' ) '  Read a file of sample arguments.'
      write ( *, '(a)' ) 
     &  '  Write a file of function values at the arguments.'
c
c  Get the number of command line arguments.
c
      num_arg = iargc ( )

      if ( num_arg < 1 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Enter the FEM file prefix:'
        read ( *, '(a)' ) fem_prefix

      else

        iarg = 1

        call getarg ( iarg, fem_prefix )

      end if

      if ( num_arg < 2 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Enter the sample file prefix:'
        read ( *, '(a)' ) sample_prefix

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

      if ( 2 * fem_node_num_max .lt. fem_node_dim * fem_node_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FEM2D_SAMPLE - Fatal error!'
        write ( *, '(a)' ) '  Not enough storage for FEM_NODE_XY.'
        write ( *, '(a,i8)' ) '  2 * FEM_NODE_NUM_MAX = ', 
     &    2 * fem_node_num_max
        write ( *, '(a,i8)' ) '  2 * FEM_NODE_NUM     = ', 
     &    2 * fem_node_num
        stop
      end if

      call dtable_data_read ( fem_node_filename, fem_node_dim, 
     &  fem_node_num, fem_node_xy )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The FEM node dimension is        ',
     &  fem_node_dim
      write ( *, '(a,i8)' ) '  The FEM node number is           ',
     &  fem_node_num

      if ( fem_node_dim .ne. 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FEM2D_SAMPLE - Fatal error!'
        write ( *, '(a)' ) '  Spatial dimension of the nodes is not 2.'
        stop
      end if

      call itable_header_read ( fem_element_filename, fem_element_order,
     &  fem_element_num )

      if ( fem_element_order_max * fem_element_num_max .lt. 
     &  fem_element_order * fem_element_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FEM2D_SAMPLE - Fatal error!'
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

      write ( *, '(a,i8)' ) '  The FEM element order is         ',
     &  fem_element_order
      write ( *, '(a,i8)' ) '  The FEM element number is        ',
     &  fem_element_num

      call dtable_header_read ( fem_value_filename, fem_value_dim, 
     &  fem_value_num )

      write ( *, '(a,i8)' ) '  The FEM value order is           ',
     &  fem_value_dim
      write ( *, '(a,i8)' ) '  the FEM value number is          ',
     &  fem_value_num

      if ( fem_value_num .ne. fem_node_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FEM2D_SAMPLE - Fatal error!'
        write ( *, '(a)' ) '  Number of values and nodes differ.'
        stop
      end if

      if ( fem_value_dim_max * fem_value_num_max .lt. 
     &  fem_value_dim * fem_value_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FEM2D_SAMPLE - Fatal error!'
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
c  Create the element neighbor array.
c
      if ( fem_element_order .eq. 3 ) then

        call triangulation_order3_neighbor_triangles ( fem_element_num, 
     &    fem_element_node, fem_element_neighbor )

      else if ( fem_element_order .eq. 6 ) then

        call triangulation_order6_neighbor_triangles ( fem_element_num, 
     &    fem_element_node, fem_element_neighbor )

      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FEM2D_SAMPLE - Fatal error!'
        write ( *, '(a)' ) '  The element order must be 3 or 6.'
        write ( *, '(a,i8)' ) '  But this data has element order = ', 
     &    fem_element_order
        return

      end if

      write ( *, '(a)' ) 
     &  '  The element neighbor array has been computed.'
c
c  Read the SAMPLE node data.
c
      call dtable_header_read ( sample_node_filename, sample_node_dim, 
     &  sample_node_num )

      if ( sample_node_dim .ne. 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FEM2D_SAMPLE - Fatal error!'
        write ( *, '(a)' ) 
     &    '  Spatial dimension of the sample nodes is not 2.'
        stop
      end if

      if ( 2 * sample_node_num_max .lt. 
     &  2 * sample_node_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'FEM2D_SAMPLE - Fatal error!'
        write ( *, '(a)' ) '  Not enough storage for SAMPLE_NODE_XY.'
        write ( *, '(a,i8)' ) 
     &  '  2 * SAMPLE_NODE_NUM_MAX = ', 
     &    2 * sample_node_num_max
        write ( *, '(a,i8)' ) 
     &  '  2     * SAMPLE_NODE_NUM     = ', 
     &    2 * sample_node_num
        stop
      end if

      call dtable_data_read ( sample_node_filename, sample_node_dim, 
     &  sample_node_num, sample_node_xy )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Sample node spatial dimension is ', 
     &  sample_node_dim
      write ( *, '(a,i8)' ) '  Sample node number is            ', 
     &  sample_node_num
c
c  Compute the sample values.
c
      sample_value_dim = fem_value_dim
      sample_value_num = sample_node_num

      call fem2d_evaluate ( fem_node_num, fem_node_xy, 
     &  fem_element_order, fem_element_num, fem_element_node, 
     &  fem_element_neighbor, fem_value_dim, fem_value, sample_node_num, 
     &  sample_node_xy, sample_value )
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
      write ( *, '(a)' ) 'FEM2D_SAMPLE'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine basis_mn_t3 ( t, n, p, phi, dphidx, dphidy )

c*********************************************************************72
c
cc BASIS_MN_T3: all bases at N points for a T3 element.
c
c  Discussion:
c
c    The routine is given the coordinates of the vertices of a triangle.
c    It works directly with these coordinates, and does not refer to a
c    reference element.
c
c    The sides of the triangle DO NOT have to lie along a coordinate
c    axis.
c
c    The routine evaluates the basis functions associated with each vertex,
c    and their derivatives with respect to X and Y.
c
c  Physical Element T3:
c
c            3
c           / \
c          /   \
c         /     \
c        /       \
c       1---------2
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
c    Input, double precision T(2,3), the coordinates of the vertices
c    of the triangle.  It is common to list these points in counter clockwise
c    order.
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision P(2,N), the points where the basis functions
c    are to be evaluated.
c
c    Output, double precision PHI(3,N), the value of the basis functions
c    at the evaluation points.
c
c    Output, double precision DPHIDX(3,N), DPHIDY(3,N), the value of the
c    derivatives at the evaluation points.
c
c  Local parameters:
c
c    Local, double precision AREA, is (twice) the area of the triangle.
c
      implicit none

      integer n

      double precision area
      double precision dphidx(3,n)
      double precision dphidy(3,n)
      integer i
      integer j
      double precision p(2,n)
      double precision phi(3,n)
      double precision t(2,3)

      area = t(1,1) * ( t(2,2) - t(2,3) ) 
     &     + t(1,2) * ( t(2,3) - t(2,1) ) 
     &     + t(1,3) * ( t(2,1) - t(2,2) )

      if ( area .eq. 0.0D+00 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BASIS_MN_T3 - Fatal error!'
        write ( *, '(a)' ) '  Element has zero area.'
        stop

      end if

      do j = 1, n

        phi(1,j) =     ( ( t(1,3) - t(1,2) ) * ( p(2,j) - t(2,2) )     
     &                 - ( t(2,3) - t(2,2) ) * ( p(1,j) - t(1,2) ) )
        dphidx(1,j) =  - ( t(2,3) - t(2,2) )
        dphidy(1,j) =    ( t(1,3) - t(1,2) )

        phi(2,j) =     ( ( t(1,1) - t(1,3) ) * ( p(2,j) - t(2,3) )     
     &                 - ( t(2,1) - t(2,3) ) * ( p(1,j) - t(1,3) ) )
        dphidx(2,j) =  - ( t(2,1) - t(2,3) )
        dphidy(2,j) =    ( t(1,1) - t(1,3) )

        phi(3,j) =     ( ( t(1,2) - t(1,1) ) * ( p(2,j) - t(2,1) )     
     &                 - ( t(2,2) - t(2,1) ) * ( p(1,j) - t(1,1) ) )
        dphidx(3,j) =  - ( t(2,2) - t(2,1) )
        dphidy(3,j) =    ( t(1,2) - t(1,1) )

      end do
c
c  Normalize.
c
      do j = 1, n
        do i = 1, 3
          phi(1:3,j) = phi(1:3,j) / area
          dphidx(1:3,j) = dphidx(1:3,j) / area
          dphidy(1:3,j) = dphidy(1:3,j) / area
        end do
      end do

      return
      end
      subroutine basis_mn_t6 ( t, n, p, phi, dphidx, dphidy )

c*********************************************************************72
c
cc BASIS_MN_T6: all bases at N points for a T6 element.
c
c  Discussion:
c
c    The routine is given the coordinates of the vertices and midside
c    nodes of a triangle.  It works directly with these coordinates, and does
c    not refer to a reference element.
c
c    This routine requires that the midside nodes be "in line"
c    with the vertices, that is, that the sides of the triangle be
c    straight.  However, the midside nodes do not actually have to
c    be halfway along the side of the triangle.
c
c  Physical element T6:
c
c    This picture indicates the assumed ordering of the six nodes
c    of the triangle.
c
c             3
c            / \
c           /   \
c          6     5
c         /       \
c        /         \
c       1-----4-----2
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
c    Input, double precision T(2,6), the nodal oordinates of the element.
c    It is common to list these points in counter clockwise order.
c
c    Input, integer N, the number of evaluation points.
c
c    Input, double precision P(2,N), the coordinates of the point where
c    the basis functions are to be evaluated.
c
c    Output, double precision PHI(6,N), the basis functions at the
c    evaluation points.
c
c    Output, double precision DPHIDX(6,N), DPHIDY(6,N), the derivatives
c    of the basis functions at the evaluation points.
c
c  Local Parameters:
c
c    Local, double precision AREA, is (twice) the area of the triangle.
c
      implicit none

      integer n

      double precision dphidx(6,n)
      double precision dphidy(6,n)
      double precision gn(n)
      double precision gx(n)
      double precision hn(n)
      double precision hx(n)
      integer j
      double precision p(2,n)
      double precision phi(6,n)
      double precision t(2,6)
c
c  Basis function 1: PHI(X,Y) = G(3,2) * H(6,4) / normalization.
c
      do j = 1, n

        gx(j) = ( p(1,j) - t(1,2) ) * ( t(2,3) - t(2,2) ) 
     &        - ( t(1,3) - t(1,2) ) * ( p(2,j) - t(2,2) )

        gn(j) = ( t(1,1) - t(1,2) ) * ( t(2,3) - t(2,2) ) 
     &        - ( t(1,3) - t(1,2) ) * ( t(2,1) - t(2,2) )

        hx(j) = ( p(1,j) - t(1,4) ) * ( t(2,6) - t(2,4) ) 
     &        - ( t(1,6) - t(1,4) ) * ( p(2,j) - t(2,4) )

        hn(j) = ( t(1,1) - t(1,4) ) * ( t(2,6) - t(2,4) ) 
     &        - ( t(1,6) - t(1,4) ) * ( t(2,1) - t(2,4) )

        phi(1,j) =     ( gx(j) * hx(j) ) / ( gn(j) * hn(j) )
        dphidx(1,j) =  (      ( t(2,3) - t(2,2) ) * hx(j) 
     &              + gx(j) * ( t(2,6) - t(2,4) ) ) / ( gn(j) * hn(j) )
        dphidy(1,j) = -(      ( t(1,3) - t(1,2) ) * hx(j) 
     &              + gx(j) * ( t(1,6) - t(1,4) ) ) / ( gn(j) * hn(j) )
c
c  Basis function 2: PHI(X,Y) = G(3,1) * H(4,5) / normalization.
c
        gx(j) = ( p(1,j) - t(1,1) ) * ( t(2,3) - t(2,1) ) 
     &        - ( t(1,3) - t(1,1) ) * ( p(2,j) - t(2,1) )

        gn(j) = ( t(1,2) - t(1,1) ) * ( t(2,3) - t(2,1) ) 
     &        - ( t(1,3) - t(1,1) ) * ( t(2,2) - t(2,1) )

        hx(j) = ( p(1,j) - t(1,5) ) * ( t(2,4) - t(2,5) ) 
     &        - ( t(1,4) - t(1,5) ) * ( p(2,j) - t(2,5) )

        hn(j) = ( t(1,2) - t(1,5) ) * ( t(2,4) - t(2,5) ) 
     &        - ( t(1,4) - t(1,5) ) * ( t(2,2) - t(2,5) )

        phi(2,j) = ( gx(j) * hx(j) ) / ( gn(j) * hn(j) )
        dphidx(2,j) =  (      ( t(2,3) - t(2,1) ) * hx(j) 
     &              + gx(j) * ( t(2,4) - t(2,5) ) ) / ( gn(j) * hn(j) )
        dphidy(2,j) = -(      ( t(1,3) - t(1,1) ) * hx(j) 
     &              + gx(j) * ( t(1,4) - t(1,5) ) ) / ( gn(j) * hn(j) )
c
c  Basis function 3: PHI(X,Y) = G(1,2) * H(5,6) / normalization.
c
        gx(j) = ( p(1,j) - t(1,2) ) * ( t(2,1) - t(2,2) ) 
     &        - ( t(1,1) - t(1,2) ) * ( p(2,j) - t(2,2) )

        gn(j) = ( t(1,3) - t(1,2) ) * ( t(2,1) - t(2,2) ) 
     &        - ( t(1,1) - t(1,2) ) * ( t(2,3) - t(2,2) )

        hx(j) = ( p(1,j) - t(1,6) ) * ( t(2,5) - t(2,6) ) 
     &        - ( t(1,5) - t(1,6) ) * ( p(2,j) - t(2,6) )

        hn(j) = ( t(1,3) - t(1,6) ) * ( t(2,5) - t(2,6) ) 
     &        - ( t(1,5) - t(1,6) ) * ( t(2,3) - t(2,6) )

        phi(3,j) = ( gx(j) * hx(j) ) / ( gn(j) * hn(j) )
        dphidx(3,j) =  (      ( t(2,1) - t(2,2) ) * hx(j)
     &              + gx(j) * ( t(2,5) - t(2,6) ) ) / ( gn(j) * hn(j) )
        dphidy(3,j) = -(      ( t(1,1) - t(1,2) ) * hx(j)
     &              + gx(j) * ( t(1,5) - t(1,6) ) ) / ( gn(j) * hn(j) )
c
c  Basis function 4: PHI(X,Y) = G(1,3) * H(2,3) / normalization.
c
        gx(j) = ( p(1,j) - t(1,3) ) * ( t(2,1) - t(2,3) )
     &        - ( t(1,1) - t(1,3) ) * ( p(2,j) - t(2,3) )

        gn(j) = ( t(1,4) - t(1,3) ) * ( t(2,1) - t(2,3) )
     &        - ( t(1,1) - t(1,3) ) * ( t(2,4) - t(2,3) )

        hx(j) = ( p(1,j) - t(1,3) ) * ( t(2,2) - t(2,3) )
     &        - ( t(1,2) - t(1,3) ) * ( p(2,j) - t(2,3) )

        hn(j) = ( t(1,4) - t(1,3) ) * ( t(2,2) - t(2,3) )
     &        - ( t(1,2) - t(1,3) ) * ( t(2,4) - t(2,3) )

        phi(4,j) = ( gx(j) * hx(j) ) / ( gn(j) * hn(j) )
        dphidx(4,j) =  (      ( t(2,1) - t(2,3) ) * hx(j)
     &             + gx(j) * ( t(2,2) - t(2,3) ) ) / ( gn(j) * hn(j) )
        dphidy(4,j) = -(      ( t(1,1) - t(1,3) ) * hx(j)
     &             + gx(j) * ( t(1,2) - t(1,3) ) ) / ( gn(j) * hn(j) )
c
c  Basis function 5: PHI(X,Y) = G(2,1) * H(3,1) / normalization.
c
        gx(j) = ( p(1,j) - t(1,1) ) * ( t(2,2) - t(2,1) )
     &        - ( t(1,2) - t(1,1) ) * ( p(2,j) - t(2,1) )

        gn(j) = ( t(1,5) - t(1,1) ) * ( t(2,2) - t(2,1) )
     &        - ( t(1,2) - t(1,1) ) * ( t(2,5) - t(2,1) )

        hx(j) = ( p(1,j) - t(1,1) ) * ( t(2,3) - t(2,1) )
     &        - ( t(1,3) - t(1,1) ) * ( p(2,j) - t(2,1) )

        hn(j) = ( t(1,5) - t(1,1) ) * ( t(2,3) - t(2,1) )
     &        - ( t(1,3) - t(1,1) ) * ( t(2,5) - t(2,1) )

        phi(5,j) = ( gx(j) * hx(j) ) / ( gn(j) * hn(j) )
        dphidx(5,j) =  (      ( t(2,2) - t(2,1) ) * hx(j)
     &             + gx(j) * ( t(2,3) - t(2,1) ) ) / ( gn(j) * hn(j) )
        dphidy(5,j) = -(      ( t(1,2) - t(1,1) ) * hx(j)
     &             + gx(j) * ( t(1,3) - t(1,1) ) ) / ( gn(j) * hn(j) )
c
c  Basis function 6: PHI(X,Y) = G(1,2) * H(3,2) / normalization.
c
        gx(j) = ( p(1,j) - t(1,2) ) * ( t(2,1) - t(2,2) )
     &        - ( t(1,1) - t(1,2) ) * ( p(2,j) - t(2,2) )

        gn(j) = ( t(1,6) - t(1,2) ) * ( t(2,1) - t(2,2) )
     &        - ( t(1,1) - t(1,2) ) * ( t(2,6) - t(2,2) )

        hx(j) = ( p(1,j) - t(1,2) ) * ( t(2,3) - t(2,2) )
     &        - ( t(1,3) - t(1,2) ) * ( p(2,j) - t(2,2) )

        hn(j) = ( t(1,6) - t(1,2) ) * ( t(2,3) - t(2,2) )
     &        - ( t(1,3) - t(1,2) ) * ( t(2,6) - t(2,2) )

        phi(6,j) = ( gx(j) * hx(j) ) / ( gn(j) * hn(j) )
        dphidx(6,j) =  (      ( t(2,1) - t(2,2) ) * hx(j)
     &             + gx(j) * ( t(2,3) - t(2,2) ) ) / ( gn(j) * hn(j) )
        dphidy(6,j) = -(      ( t(1,1) - t(1,2) ) * hx(j)
     &             + gx(j) * ( t(1,3) - t(1,2) ) ) / ( gn(j) * hn(j) )

      end do

      return
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
      subroutine fem2d_evaluate ( fem_node_num, fem_node_xy, 
     &  fem_element_order, fem_element_num, fem_element_node, 
     &  fem_element_neighbor, fem_value_dim, fem_value, 
     &  sample_node_num, sample_node_xy, sample_value )

c*********************************************************************72
c
cc FEM2D_EVALUATE samples an FEM function on a T3 or T6 triangulation.
c
c  Discussion:
c
c    Note that the sample values returned are true values of the underlying
c    finite element function.  They are NOT produced by constructing some
c    other function that interpolates the data at the finite element nodes
c    (something which MATLAB's griddata function can easily do.)  Instead, 
c    each sampling node is located within one of the associated finite
c    element triangles, and the finite element function is developed and 
c    evaluated there.  
c
c    MATLAB's scattered data interpolation is wonderful, but it cannot
c    be guaranteed to reproduce the finite element function corresponding
c    to nodal data.  This routine can (or at least tries toc).
c 
c    So if you are using finite elements, then using THIS routine
c    (but not MATLAB's griddata function), what you see is what you havec
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
c    Input, integer FEM_NODE_NUM, the number of nodes.
c
c    Input, double precision FEM_NODE_XY(2,FEM_NODE_NUM), the coordinates 
c    of the nodes.
c
c    Input, integer FEM_ELEMENT_ORDER, the order of the elements, 
c    either 3 or 6.
c
c    Input, integer FEM_ELEMENT_NUM, the number of triangles.
c
c    Input, integer 
c    FEM_ELEMENT_NODE(FEM_ELEMENT_ORDER,FEM_ELEMENT_NUM), the
c    nodes that make up each triangle.
c
c    Input, integer FEM_ELEMENT_NEIGHBOR(3,FEM_ELEMENT_NUM), the 
c    index of the neighboring triangle on each side, or -1 if no neighbor there.
c
c    Input, integer FEM_VALUE_DIM, the "dimension" of the values.
c
c    Input, double precision FEM_VALUE(FEM_VALUE_DIM,FEM_NODE_NUM), the 
c    finite element coefficient values at each node.
c
c    Input, integer SAMPLE_NODE_NUM, the number of sample nodes.
c
c    Input, double precision SAMPLE_NODE_XY(2,SAMPLE_NODE_NUM), the sample nodes.
c
c    Output, double precision SAMPLE_VALUE(FEM_VALUE_DIM,SAMPLE_NODE_NUM),
c    the sampled values.
c
      implicit none

      integer fem_element_num
      integer fem_element_order
      integer fem_node_num
      integer fem_value_dim
      integer sample_node_num

      double precision b(fem_element_order)
      double precision dbdx(fem_element_order)
      double precision dbdy(fem_element_order)
      double precision dot
      integer edge
      integer fem_element_neighbor(3,fem_element_num)
      integer fem_element_node(fem_element_order,fem_element_num)
      double precision fem_node_xy(2,fem_node_num)
      double precision fem_value(fem_value_dim,fem_node_num)
      integer i
      integer j
      integer k
      double precision p_xy(2)
      double precision sample_node_xy(2,sample_node_num)
      double precision sample_value(fem_value_dim,sample_node_num)
      integer t
      integer t_node(fem_element_order)
      double precision t_xy(2,fem_element_order)
c
c  For each sample point: find the triangle T that contains it,
c  and evaluate the finite element function there.
c
      do j = 1, sample_node_num

        p_xy(1:2) = sample_node_xy(1:2,j)
c
c  Find the triangle T that contains the point.
c
        call triangulation_search ( fem_node_num, fem_node_xy, 
     &    fem_element_order, fem_element_num, fem_element_node, 
     &    fem_element_neighbor, p_xy, t, edge )
c
c  Evaluate the finite element basis functions at the point in T.
c
        t_node(1:fem_element_order) = 
     &    fem_element_node(1:fem_element_order,t)

        t_xy(1:2,1:fem_element_order) = fem_node_xy(1:2,t_node)

        if ( fem_element_order .eq. 3 ) then
          call basis_mn_t3 ( t_xy, 1, p_xy, b, dbdx, dbdy )
        else if ( fem_element_order .eq. 6 ) then
          call basis_mn_t6 ( t_xy, 1, p_xy, b, dbdx, dbdy )
        end if
c
c  Multiply by the finite element values to get the sample values.
c
        do i = 1, fem_value_dim
          dot = 0.0D+00
          do k = 1, fem_element_order
            dot = dot + fem_value(i,t_node(k)) * b(k)
          end do
          sample_value(i,j) = dot
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
      character * ( 255 ) line
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
      character * ( 255 ) line
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
      subroutine get_seed ( seed )

c*********************************************************************72
c
cc GET_SEED returns a seed for the random number generator.
c
c  Discussion:
c
c    The seed depends on the current time, and ought to be (slightly)
c    different every millisecond.  Thus, calling this routine several
c    times in succession will probably return the SAME seed, but
c    calling it a few minutes or days apart will turn a suitably
c    "random" seed.
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
c    Output, integer SEED, a pseudorandom seed value.
c
      implicit none

      integer day
      integer hour
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer milli
      integer minute
      integer month
      integer second
      integer seed
      double precision temp
      character * ( 10 ) time
      character * ( 8 ) date
      integer year

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) year, month, day
      read ( time, '(i2,i2,i2,1x,i3)' ) hour, minute, second, milli

      temp = 0.0D+00
      temp = temp + dble ( month - 1 ) / 11.0D+00
      temp = temp + dble ( day   - 1 ) / 30.0D+00
      temp = temp + dble ( hour      ) / 23.0D+00
      temp = temp + dble ( minute    ) / 59.0D+00
      temp = temp + dble ( second    ) / 59.0D+00
      temp = temp + dble ( milli     ) / 999.0D+00

      temp = temp / 6.0D+00
c
c  Force 0 < TEMP <= 1.
c
10    continue

      if ( temp .le. 0.0D+00 ) then
        temp = temp + 1.0D+00
        go to 10
      end if

20    continue

      if ( 1.0D+00 .lt. temp ) then
        temp = temp - 1.0D+00
        go to 20
      end if

      seed = int ( dble ( i4_huge ) * temp )
c
c  Never use a seed of 0 or maximum integer.
c
      if ( seed .eq. 0 ) then
        seed = 1
      end if

      if ( seed .eq. i4_huge ) then
        seed = seed - 1
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
      function i4_uniform ( a, b, seed )

c*********************************************************************72
c
cc I4_UNIFORM returns a scaled pseudorandom I4.
c
c  Discussion:
c
c    An I4 is an integer value.
c
c    The pseudorandom number should be uniformly distributed
c    between A and B.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Second Edition,
c    Springer, 1987,
c    ISBN: 0387964673,
c    LC: QA76.9.C65.B73.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, December 1986, pages 362-376.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley, 1998,
c    ISBN: 0471134031,
c    LC: T57.62.H37.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, Number 2, 1969, pages 136-143.
c
c  Parameters:
c
c    Input, integer A, B, the limits of the interval.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, integer I4_UNIFORM, a number between A and B.
c
      implicit none

      integer a
      integer b
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer i4_uniform
      integer k
      real r
      integer seed
      integer value

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'I4_UNIFORM - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + i4_huge
      end if

      r = real ( seed ) * 4.656612875E-10
c
c  Scale R to lie between A-0.5 and B+0.5.
c
      r = ( 1.0E+00 - r ) * ( real ( min ( a, b ) ) - 0.5E+00 )
     &  +             r   * ( real ( max ( a, b ) ) + 0.5E+00 )
c
c  Use rounding to convert R to an integer between A and B.
c
      value = nint ( r )

      value = max ( value, min ( a, b ) )
      value = min ( value, max ( a, b ) )

      i4_uniform = value

      return
      end
      subroutine i4row_compare ( m, n, a, i, j, isgn )

c*********************************************************************72
c
cc I4ROW_COMPARE compares two rows of an I4ROW.
c
c  Discussion:
c
c    An I4ROW is an M by N array of integer values, regarded
c    as an array of M rows of length N.
c
c  Example:
c
c    Input:
c
c    M = 3, N = 4, I = 2, J = 3
c
c    A = (
c      1  2  3  4
c      5  6  7  8
c      9 10 11 12 )
c
c    Output:
c
c      ISGN = -1
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
c    Input, integer M, N, the number of rows and columns.
c
c    Input, integer A(M,N), an array of M rows of vectors 
c    of length N.
c
c    Input, integer I, J, the rows to be compared.
c    I and J must be between 1 and M.
c
c    Output, integer ISGN, the results of the comparison:
c    -1, row I .lt. row J,
c     0, row I = row J,
c    +1, row J .lt. row I.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer i
      integer isgn
      integer j
      integer k
c
c  Check that I and J are legal.
c
      if ( i .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IROW_COMPARE - Fatal error!'
        write ( *, '(a)' ) '  Row index I is less than 1.'
        write ( *, '(a,i8)' ) '  I = ', i
        stop
      else if ( m .lt. i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IROW_COMPARE - Fatal error!'
        write ( *, '(a)' ) '  Row index I is out of bounds.'
        write ( *, '(a,i8)' ) '  I = ', i
        write ( *, '(a,i8)' ) '  Maximum legal value is M = ', m
        stop
      end if

      if ( j .lt. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IROW_COMPARE - Fatal error!'
        write ( *, '(a)' ) '  Row index J is less than 1.'
        write ( *, '(a,i8)' ) '  J = ', j
        stop
      else if ( m .lt. j ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IROW_COMPARE - Fatal error!'
        write ( *, '(a)' ) '  Row index J is out of bounds.'
        write ( *, '(a,i8)' ) '  J = ', j
        write ( *, '(a,i8)' ) '  Maximum legal value is M = ', m
        stop
      end if

      isgn = 0

      if ( i .eq. j ) then
        return
      end if

      k = 1

10    continue

      if ( k .le. n ) then

        if ( a(i,k) .lt. a(j,k) ) then
          isgn = -1
          return
        else if ( a(j,k) .lt. a(i,k) ) then
          isgn = +1
          return
        end if

        k = k + 1

        go to 10

      end if

      return
      end
      subroutine i4row_sort_a ( m, n, a )

c*********************************************************************72
c
cc I4ROW_SORT_A ascending sorts the rows of an I4ROW.
c
c  Discussion:
c
c    An I4ROW is an M by N array of integer values, regarded
c    as an array of M rows of length N.
c
c    In lexicographic order, the statement "X .lt. Y", applied to two
c    vectors X and Y of length M, means that there is some index I, with
c    1 .le. I .le. M, with the property that
c
c      X(J) = Y(J) for J .lt. I,
c    and
c      X(I) .lt. Y(I).
c
c    In other words, X is less than Y if, at the first index where they
c    differ, the X value is less than the Y value.
c
c  Example:
c
c    Input:
c
c      M = 5, N = 3
c
c      A =
c        3  2  1
c        2  4  3
c        3  1  8
c        2  4  2
c        1  9  9
c
c    Output:
c
c      A =
c        1  9  9
c        2  4  2
c        2  4  3
c        3  1  8
c        3  2  1
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
c    Input, integer M, the number of rows of A.
c
c    Input, integer N, the number of columns of A.
c
c    Input/output, integer A(M,N).
c    On input, the array of M rows of N-vectors.
c    On output, the rows of A have been sorted in ascending
c    lexicographic order.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer i
      integer indx
      integer isgn
      integer j

      if ( m .le. 1 ) then
        return
      end if

      if ( n .le. 0 ) then
        return
      end if
c
c  Initialize.
c
      i = 0
      indx = 0
      isgn = 0
      j = 0
c
c  Call the external heap sorter.
c
10    continue

        call sort_heap_external ( m, indx, i, j, isgn )
c
c  Interchange the I and J objects.
c
        if ( 0 .lt. indx ) then

          call i4row_swap ( m, n, a, i, j )
c
c  Compare the I and J objects.
c
        else if ( indx .lt. 0 ) then

          call i4row_compare ( m, n, a, i, j, isgn )

        else if ( indx .eq. 0 ) then

          go to 20

        end if

      go to 10

20    continue

      return
      end
      subroutine i4row_swap ( m, n, a, i1, i2 )

c*********************************************************************72
c
cc I4ROW_SWAP swaps two rows of an I4ROW.
c
c  Discussion:
c
c    An I4ROW is an M by N array of integer values, regarded
c    as an array of M rows of length N.
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
c    Input, integer M, N, the number of rows and columns.
c
c    Input/output, integer A(M,N), an array of data.
c
c    Input, integer I1, I2, the two rows to swap.
c
      implicit none

      integer m
      integer n

      integer a(m,n)
      integer i1
      integer i2
      integer row(n)
c
c  Check.
c
      if ( i1 .lt. 1 .or. m .lt. i1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IROW_SWAP - Fatal error!'
        write ( *, '(a)' ) '  I1 is out of range.'
        stop
      end if

      if ( i2 .lt. 1 .or. m .lt. i2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'IROW_SWAP - Fatal error!'
        write ( *, '(a)' ) '  I2 is out of range.'
        stop
      end if

      if ( i1 .eq. i2 ) then
        return
      end if

      row(1:n)  = a(i1,1:n)
      a(i1,1:n) = a(i2,1:n)
      a(i2,1:n) = row(1:n)

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
      subroutine sort_heap_external ( n, indx, i, j, isgn )

c*********************************************************************72
c
cc SORT_HEAP_EXTERNAL externally sorts a list of items into ascending order.
c
c  Discussion:
c
c    The actual list of data is not passed to the routine.  Hence this
c    routine may be used to sort integers, reals, numbers, names,
c    dates, shoe sizes, and so on.  After each call, the routine asks
c    the user to compare or interchange two items, until a special
c    return value signals that the sorting is completed.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 January 2007
c
c  Author:
c
c    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Albert Nijenhuis, Herbert Wilf,
c    Combinatorial Algorithms for Computers and Calculators,
c    Second Edition,
c    Academic Press, 1978,
c    ISBN: 0-12-519260-6,
c    LC: QA164.N54.
c
c  Parameters:
c
c    Input, integer N, the number of items to be sorted.
c
c    Input/output, integer INDX, the main communication signal.
c
c    The user must set INDX to 0 before the first call.
c    Thereafter, the user should not change the value of INDX until
c    the sorting is done.
c
c    On return, if INDX is
c
c      greater than 0,
c      * interchange items I and J;
c      * call again.
c
c      less than 0,
c      * compare items I and J;
c      * set ISGN = -1 if I .lt. J, ISGN = +1 if J .lt. I;
c      * call again.
c
c      equal to 0, the sorting is done.
c
c    Output, integer I, J, the indices of two items.
c    On return with INDX positive, elements I and J should be interchanged.
c    On return with INDX negative, elements I and J should be compared, and
c    the result reported in ISGN on the next call.
c
c    Input, integer ISGN, results of comparison of elements I and J.
c    (Used only when the previous call returned INDX less than 0).
c    ISGN .le. 0 means I is less than or equal to J;
c    0 .le. ISGN means I is greater than or equal to J.
c
      implicit none

      integer i
      integer i_save
      integer indx
      integer isgn
      integer j
      integer j_save
      integer k
      integer k1
      integer n
      integer n1

      save i_save
      save j_save
      save k
      save k1
      save n1

      data i_save / 0 /
      data j_save / 0 /
      data k / 0 /
      data k1 / 0 /
      data n1 / 0 /
c
c  INDX = 0: This is the first call.
c
      if ( indx .eq. 0 ) then

        i_save = 0
        j_save = 0
        k = n / 2
        k1 = k
        n1 = n
c
c  INDX .lt. 0: The user is returning the results of a comparison.
c
      else if ( indx .lt. 0 ) then

        if ( indx .eq. -2 ) then

          if ( isgn .lt. 0 ) then
            i_save = i_save + 1
          end if

          j_save = k1
          k1 = i_save
          indx = -1
          i = i_save
          j = j_save
          return

        end if

        if ( 0 .lt. isgn ) then
          indx = 2
          i = i_save
          j = j_save
          return
        end if

        if ( k .le. 1 ) then

          if ( n1 .eq. 1 ) then
            i_save = 0
            j_save = 0
            indx = 0
          else
            i_save = n1
            n1 = n1 - 1
            j_save = 1
            indx = 1
          end if

          i = i_save
          j = j_save
          return

        end if

        k = k - 1
        k1 = k
c
c  0 .lt. INDX, the user was asked to make an interchange.
c
      else if ( indx .eq. 1 ) then

        k1 = k

      end if

10    continue

        i_save = 2 * k1

        if ( i_save .eq. n1 ) then
          j_save = k1
          k1 = i_save
          indx = -1
          i = i_save
          j = j_save
          return
        else if ( i_save .le. n1 ) then
          j_save = i_save + 1
          indx = -2
          i = i_save
          j = j_save
          return
        end if

        if ( k .le. 1 ) then
          go to 20
        end if

        k = k - 1
        k1 = k

      go to 10

20    continue

      if ( n1 .eq. 1 ) then
        i_save = 0
        j_save = 0
        indx = 0
        i = i_save
        j = j_save
      else
        i_save = n1
        n1 = n1 - 1
        j_save = 1
        indx = 1
        i = i_save
        j = j_save
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
      subroutine triangulation_order3_neighbor_triangles ( 
     &  triangle_num, triangle_node, triangle_neighbor )

c*********************************************************************72
c
cc TRIANGULATION_ORDER3_NEIGHBOR_TRIANGLES determines triangle neighbors.
c
c  Discussion:
c
c    A triangulation of a set of nodes can be completely described by
c    the coordinates of the nodes, and the list of nodes that make up
c    each triangle.  However, in some cases, it is necessary to know
c    triangle adjacency information, that is, which triangle, if any,
c    is adjacent to a given triangle on a particular side.
c
c    This routine creates a data structure recording this information.
c
c    The primary amount of work occurs in sorting a list of 3 * TRIANGLE_NUM
c    data items.
c
c    Note that ROW is a work array allocated dynamically inside this
c    routine.  It is possible, for very large values of TRIANGLE_NUM,
c    that the necessary amount of memory will not be accessible, and the
c    routine will fail.  This is a limitation of the implementation of
c    dynamic arrays in FORTRAN90.  One way to get around this would be
c    to require the user to declare ROW in the calling routine
c    as an allocatable array, get the necessary memory explicitly with
c    an ALLOCATE statement, and then pass ROW into this routine.
c
c    Of course, the point of dynamic arrays was to make it easy to
c    hide these sorts of temporary work arrays from the poor user!
c
c  Example:
c
c    The input information from TRIANGLE_NODE:
c
c    Triangle   Nodes
c    --------   ---------------
c     1         3      4      1
c     2         3      1      2
c     3         3      2      8
c     4         2      1      5
c     5         8      2     13
c     6         8     13      9
c     7         3      8      9
c     8        13      2      5
c     9         9     13      7
c    10         7     13      5
c    11         6      7      5
c    12         9      7      6
c    13        10      9      6
c    14         6      5     12
c    15        11      6     12
c    16        10      6     11
c
c    The output information in TRIANGLE_NEIGHBOR:
c
c    Triangle  Neighboring Triangles
c    --------  ---------------------
c
c     1        -1     -1      2
c     2         1      4      3
c     3         2      5      7
c     4         2     -1      8
c     5         3      8      6
c     6         5      9      7
c     7         3      6     -1
c     8         5      4     10
c     9         6     10     12
c    10         9      8     11
c    11        12     10     14
c    12         9     11     13
c    13        -1     12     16
c    14        11     -1     15
c    15        16     14     -1
c    16        13     15     -1
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
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer TRIANGLE_NODE(3,TRIANGLE_NUM), the nodes 
c    that make up each triangle.
c
c    Output, integer TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), the 
c    triangles that are direct neighbors of a given triangle.  
c    TRIANGLE_NEIGHBOR(1,I) is the index of the triangle which touches side 1, 
c    defined by nodes 2 and 3, and so on.  TRIANGLE_NEIGHBOR(1,I) is negative
c    if there is no neighbor on that side.  In this case, that side of the 
c    triangle lies on the boundary of the triangulation.
c
      implicit none

      integer triangle_num
      integer triangle_order
      parameter ( triangle_order = 3 )

      integer i
      integer irow
      integer j
      integer k
      integer row(3*triangle_num,4)
      integer side1
      integer side2
      integer triangle_neighbor(3,triangle_num)
      integer tri
      integer triangle_node(triangle_order,triangle_num)
      integer tri1
      integer tri2
c
c  Step 1.
c  From the list of nodes for triangle T, of the form: (I,J,K)
c  construct the three neighbor relations:
c
c    (I,J,1,T) or (J,I,1,T),
c    (J,K,2,T) or (K,J,2,T),
c    (K,I,3,T) or (I,K,3,T)
c
c  where we choose (I,J,1,T) if I .lt. J, or else (J,I,1,T)
c
      do tri = 1, triangle_num

        i = triangle_node(1,tri)
        j = triangle_node(2,tri)
        k = triangle_node(3,tri)

        if ( i .lt. j ) then
          row(3*(tri-1)+1,1) = i
          row(3*(tri-1)+1,2) = j
          row(3*(tri-1)+1,3) = 1
          row(3*(tri-1)+1,4) = tri
        else
          row(3*(tri-1)+1,1) = j
          row(3*(tri-1)+1,2) = i
          row(3*(tri-1)+1,3) = 1
          row(3*(tri-1)+1,4) = tri
        end if

        if ( j .lt. k ) then
          row(3*(tri-1)+2,1) = j
          row(3*(tri-1)+2,2) = k
          row(3*(tri-1)+2,3) = 2
          row(3*(tri-1)+2,4) = tri
        else
          row(3*(tri-1)+2,1) = k
          row(3*(tri-1)+2,2) = j
          row(3*(tri-1)+2,3) = 2
          row(3*(tri-1)+2,4) = tri
        end if

        if ( k .lt. i ) then
          row(3*(tri-1)+3,1) = k
          row(3*(tri-1)+3,2) = i
          row(3*(tri-1)+3,3) = 3
          row(3*(tri-1)+3,4) = tri
        else
          row(3*(tri-1)+3,1) = i
          row(3*(tri-1)+3,2) = k
          row(3*(tri-1)+3,3) = 3
          row(3*(tri-1)+3,4) = tri
        end if

      end do
c
c  Step 2. Perform an ascending dictionary sort on the neighbor relations.
c  We only intend to sort on columns 1 and 2; the routine we call here
c  sorts on columns 1 through 4 but that won't hurt us.
c
c  What we need is to find cases where two triangles share an edge.
c  Say they share an edge defined by the nodes I and J.  Then there are
c  two rows of ROW that start out ( I, J, ?, ? ).  By sorting ROW,
c  we make sure that these two rows occur consecutively.  That will
c  make it easy to notice that the triangles are neighbors.
c
      call i4row_sort_a ( 3*triangle_num, 4, row )
c
c  Step 3. Neighboring triangles show up as consecutive rows with
c  identical first two entries.  Whenever you spot this happening,
c  make the appropriate entries in TRIANGLE_NEIGHBOR.
c
      triangle_neighbor(1:3,1:triangle_num) = -1

      irow = 1

10    continue

        if ( 3 * triangle_num .le. irow ) then
          go to 20
        end if

        if ( row(irow,1) .ne. row(irow+1,1) .or. 
     &       row(irow,2) .ne. row(irow+1,2) ) then
          irow = irow + 1
          go to 10
        end if

        side1 = row(irow,3)
        tri1 = row(irow,4)
        side2 = row(irow+1,3)
        tri2 = row(irow+1,4)

        triangle_neighbor(side1,tri1) = tri2
        triangle_neighbor(side2,tri2) = tri1

        irow = irow + 2

      go to 10

20    continue

      return
      end
      subroutine triangulation_order6_neighbor_triangles ( triangle_num, 
     &  triangle_node, triangle_neighbor )

c*********************************************************************72
c
cc TRIANGULATION_ORDER6_NEIGHBOR_TRIANGLES determines triangle neighbors.
c
c  Discussion:
c
c    A triangulation of a set of nodes can be completely described by
c    the coordinates of the nodes, and the list of nodes that make up
c    each triangle.  However, in some cases, it is necessary to know
c    triangle adjacency information, that is, which triangle, if any,
c    is adjacent to a given triangle on a particular side.
c
c    This routine creates a data structure recording this information.
c
c    The primary amount of work occurs in sorting a list of 3 * TRIANGLE_NUM
c    data items.
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
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer TRIANGLE_ORDER(6,TRIANGLE_NUM), the nodes 
c    that make up each triangle.
c
c    Output, integer TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), the 
c    three triangles that are direct neighbors of a given triangle. 
c    TRIANGLE_NEIGHBOR(1,I) is the index of the triangle which touches side 1, 
c    defined by nodes 2 and 3, and so on.  TRIANGLE_NEIGHBOR(1,I) is negative 
c    if there is no neighbor on that side.  In this case, that side of the 
c    triangle lies on the boundary of the triangulation.
c
      implicit none

      integer triangle_num
      integer triangle_order
      parameter ( triangle_order = 6 )

      integer i
      integer irow
      integer j
      integer k
      integer row(3*triangle_num,4)
      integer side1
      integer side2
      integer tri
      integer triangle_node(triangle_order,triangle_num)
      integer triangle_neighbor(3,triangle_num)
      integer tri1
      integer tri2
c
c  Step 1.
c  From the list of vertices for triangle T, of the form: (I,J,K)
c  construct the three neighbor relations:
c
c    (I,J,1,T) or (J,I,1,T),
c    (J,K,2,T) or (K,J,2,T),
c    (K,I,3,T) or (I,K,3,T)
c
c  where we choose (I,J,1,T) if I .lt. J, or else (J,I,1,T)
c
      do tri = 1, triangle_num

        i = triangle_node(1,tri)
        j = triangle_node(2,tri)
        k = triangle_node(3,tri)

        if ( i .lt. j ) then
          row(3*(tri-1)+1,1) = i
          row(3*(tri-1)+1,2) = j
          row(3*(tri-1)+1,3) = 1
          row(3*(tri-1)+1,4) = tri
        else
          row(3*(tri-1)+1,1) = j
          row(3*(tri-1)+1,2) = i
          row(3*(tri-1)+1,3) = 1
          row(3*(tri-1)+1,4) = tri
        end if

        if ( j .lt. k ) then
          row(3*(tri-1)+2,1) = j
          row(3*(tri-1)+2,2) = k
          row(3*(tri-1)+2,3) = 2
          row(3*(tri-1)+2,4) = tri
        else
          row(3*(tri-1)+2,1) = k
          row(3*(tri-1)+2,2) = j
          row(3*(tri-1)+2,3) = 2
          row(3*(tri-1)+2,4) = tri
        end if

        if ( k .lt. i ) then
          row(3*(tri-1)+3,1) = k
          row(3*(tri-1)+3,2) = i
          row(3*(tri-1)+3,3) = 3
          row(3*(tri-1)+3,4) = tri
        else
          row(3*(tri-1)+3,1) = i
          row(3*(tri-1)+3,2) = k
          row(3*(tri-1)+3,3) = 3
          row(3*(tri-1)+3,4) = tri
        end if

      end do
c
c  Step 2. Perform an ascending dictionary sort on the neighbor relations.
c  We only intend to sort on columns 1 and 2; the routine we call here
c  sorts on columns 1 through 4 but that won't hurt us.
c
c  What we need is to find cases where two triangles share an edge.
c  Say they share an edge defined by the nodes I and J.  Then there are
c  two rows of ROW that start out ( I, J, ?, ? ).  By sorting ROW,
c  we make sure that these two rows occur consecutively.  That will
c  make it easy to notice that the triangles are neighbors.
c
      call i4row_sort_a ( 3*triangle_num, 4, row )
c
c  Step 3. Neighboring triangles show up as consecutive rows with
c  identical first two entries.  Whenever you spot this happening,
c  make the appropriate entries in TRIANGLE_NEIGHBOR.
c
      triangle_neighbor(1:3,1:triangle_num) = -1

      irow = 1

10    continue

        if ( 3 * triangle_num .le. irow ) then
          go to 20
        end if

        if ( row(irow,1) .ne. row(irow+1,1) .or. 
     &       row(irow,2) .ne. row(irow+1,2) ) then
          irow = irow + 1
          go to 10
        end if

        side1 = row(irow,3)
        tri1 = row(irow,4)
        side2 = row(irow+1,3)
        tri2 = row(irow+1,4)

        triangle_neighbor(side1,tri1) = tri2
        triangle_neighbor(side2,tri2) = tri1

        irow = irow + 2

      go to 10

20    continue


      return
      end
      subroutine triangulation_search ( node_num, node_xy, 
     & triangle_order, triangle_num, triangle_node, triangle_neighbor, 
     &  p, triangle_index, edge )

c*********************************************************************72
c
cc TRIANGULATION_SEARCH searches a triangulation for a point.
c
c  Discussion:
c
c    The algorithm "walks" from one triangle to its neighboring triangle,
c    and so on, until a triangle is found containing point P, or P is found
c    to be outside the convex hull.
c
c    The algorithm computes the barycentric coordinates of the point with
c    respect to the current triangle.  If all three quantities are positive,
c    the point is contained in the triangle.  If the I-th coordinate is
c    negative, then P lies on the far side of edge I, which is opposite
c    from vertex I.  This gives a hint as to where to search next.
c
c    For a Delaunay triangulation, the search is guaranteed to terminate.
c    For other triangulations, a cycle may occur.
c
c    Note the surprising fact that, even for a Delaunay triangulation of
c    a set of nodes, the nearest node to P need not be one of the
c    vertices of the triangle containing P.
c
c    The code can be called for triangulations of any order, but only
c    the first three nodes in each triangle are considered.  Thus, if
c    higher order triangles are used, and the extra nodes are intended
c    to give the triangle a polygonal shape, these will have no effect,
c    and the results obtained here might be misleading.
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
c    John Burkardt.
c
c  Reference:
c
c    Barry Joe,
c    GEOMPACK - a software package for the generation of meshes
c    using geometric algorithms,
c    Advances in Engineering Software,
c    Volume 13, pages 325-331, 1991.
c
c  Parameters:
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Input, double precision NODE_XY(2,NODE_NUM), the coordinates of the nodes.
c
c    Input, integer TRIANGLE_ORDER, the order of the triangles.
c
c    Input, integer TRIANGLE_NUM, the number of triangles.
c
c    Input, integer TRIANGLE_NODE(TRIANGLE_ORDER,TRIANGLE_NUM),
c    the nodes that make up each triangle.
c
c    Input, integer TRIANGLE_NEIGHBOR(3,TRIANGLE_NUM), the 
c    triangle neighbor list.
c
c    Input, double precision P(2), the coordinates of a point.
c
c    Output, integer TRIANGLE_INDEX, the index of the triangle 
c    where the search ended.  If a cycle occurred, then TRIANGLE_INDEX = -1.
c
c    Output, integer EDGE, indicates the position of the point P in
c    triangle TRIANGLE_INDEX:
c    0, the interior or boundary of the triangle;
c    -1, outside the convex hull of the triangulation, past edge 1;
c    -2, outside the convex hull of the triangulation, past edge 2;
c    -3, outside the convex hull of the triangulation, past edge 3.
c
      implicit none

      integer dim_num
      parameter ( dim_num = 2 )
      integer node_num
      integer triangle_num
      integer triangle_order

      integer a
      double precision alpha
      integer b
      double precision beta
      integer c
      integer count
      double precision det
      double precision dxp
      double precision dxa
      double precision dxb
      double precision dyp
      double precision dya
      double precision dyb
      integer edge
      double precision gamma
      integer i4_uniform
      double precision node_xy(dim_num,node_num)
      double precision p(dim_num)
      integer seed
      integer triangle_node(triangle_order,triangle_num)
      integer triangle_index
      integer triangle_neighbor(3,triangle_num)

      count = 0
      edge = 0

      call get_seed ( seed )

      triangle_index = i4_uniform ( 1, triangle_num, seed )

10    continue

        count = count + 1

        if ( triangle_num .lt. count ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'TRIANGULATION_SEARCH - Fatal error!'
          write ( *, '(a)' ) '  The algorithm seems to be cycling.'
          triangle_index = -1
          edge = -1
          stop
        end if
c
c  Get the nodes of triangle TRIANGLE_INDEX.
c
        a = triangle_node(1,triangle_index)
        b = triangle_node(2,triangle_index)
        c = triangle_node(3,triangle_index)
c
c  Using vertex C as a base, compute the distances to vertices A and B,
c  and the point P.
c
        dxa = node_xy(1,a) - node_xy(1,c)
        dya = node_xy(2,a) - node_xy(2,c)

        dxb = node_xy(1,b) - node_xy(1,c)
        dyb = node_xy(2,b) - node_xy(2,c)

        dxp = p(1)         - node_xy(1,c)
        dyp = p(2)         - node_xy(2,c)

        det = dxa * dyb - dya * dxb
c
c  Compute the barycentric coordinates of the point P with respect
c  to this triangle.
c
        alpha = ( dxp * dyb - dyp * dxb ) / det
        beta =  ( dxa * dyp - dya * dxp ) / det
        gamma = 1.0D+00 - alpha - beta
c
c  If the barycentric coordinates are all positive, then the point
c  is inside the triangle and we're done.
c
        if ( 0.0D+00 .le. alpha .and. 
     &       0.0D+00 .le. beta  .and. 
     &       0.0D+00 .le. gamma ) then
          go to 20
        end if
c
c  At least one barycentric coordinate is negative.
c
c  If there is a negative barycentric coordinate for which there exists
c  an opposing triangle neighbor closer to the point, move to that triangle.
c
c  (Two coordinates could be negative, in which case we could go for the
c  most negative one, or the most negative one normalized by the actual
c  distance it represents).
c
        if ( alpha .lt. 0.0D+00 .and. 
     &    0 .lt. triangle_neighbor(2,triangle_index) ) then
          triangle_index = triangle_neighbor(2,triangle_index)
          go to 10
        else if ( beta .lt. 0.0D+00 .and.
     &    0 .lt. triangle_neighbor(3,triangle_index) ) then
          triangle_index = triangle_neighbor(3,triangle_index)
          go to 10
        else if ( gamma .lt. 0.0D+00 .and.
     &    0 .lt. triangle_neighbor(1,triangle_index) ) then
          triangle_index = triangle_neighbor(1,triangle_index)
          go to 10
        end if
c
c  All negative barycentric coordinates correspond to vertices opposite
c  sides on the convex hull.
c
c  Note the edge and exit.
c
        if ( alpha .lt. 0.0D+00 ) then
          edge = -2
          go to 20
        else if ( beta .lt. 0.0D+00 ) then
          edge = -3
          go to 20
        else if ( gamma .lt. 0.0D+00 ) then
          edge = -1
          go to 20
        end if

      go to 10

20    continue

      return
      end
