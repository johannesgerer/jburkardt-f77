      subroutine bandwidth_mesh ( element_order, element_num, 
     &  element_node, ml, mu, m )

c*********************************************************************72
c
cc BANDWIDTH_MESH: bandwidth of finite element mesh.
c
c  Discussion:
c
c    The quantity computed here is the "geometric" bandwidth determined
c    by the finite element mesh alone.
c
c    If a single finite element variable is associated with each node
c    of the mesh, and if the nodes and variables are numbered in the
c    same way, then the geometric bandwidth is the same as the bandwidth
c    of a typical finite element matrix.
c
c    The bandwidth M is defined in terms of the lower and upper bandwidths:
c
c      M = ML + 1 + MU
c
c    where 
c
c      ML = maximum distance from any diagonal entry to a nonzero
c      entry in the same row, but earlier column,
c
c      MU = maximum distance from any diagonal entry to a nonzero
c      entry in the same row, but later column.
c
c    Because the finite element node adjacency relationship is symmetric,
c    we are guaranteed that ML = MU.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer ELEMENT_ORDER, the order of the elements.
c
c    Input, integer ELEMENT_NUM, the number of elements.
c
c    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM); 
c    ELEMENT_NODE(I,J) is the global index of local node I in element J.
c
c    Output, integer ML, MU, the lower and upper bandwidths of 
c    the matrix.
c
c    Output, integer M, the bandwidth of the matrix.
c
      implicit none

      integer element_num
      integer element_order

      integer element
      integer element_node(element_order,element_num)
      integer global_i
      integer global_j
      integer local_i
      integer local_j
      integer m
      integer ml
      integer mu

      ml = 0
      mu = 0
      
      do element = 1, element_num

        do local_i = 1, element_order
          global_i = element_node(local_i,element)

          do local_j = 1, element_order
            global_j = element_node(local_j,element)

            mu = max ( mu, global_j - global_i )
            ml = max ( ml, global_i - global_j )

          end do
        end do
      end do

      m = ml + 1 + mu

      return
      end
      subroutine bandwidth_var ( element_order, element_num, 
     &  element_node, node_num, var_node, var_num, var, ml, mu, m )

c*********************************************************************72
c
cc BANDWIDTH_VAR determines the bandwidth for finite element variables.
c
c  Discussion:
c
c    We assume that, attached to each node in the finite element mesh
c    there are a (possibly zero) number of finite element variables.
c    We wish to determine the bandwidth necessary to store the stiffness
c    matrix associated with these variables.
c
c    An entry K(I,J) of the stiffness matrix must be zero unless the
c    variables I and J correspond to nodes N(I) and N(J) which are
c    common to some element.
c
c    In order to determine the bandwidth of the stiffness matrix, we
c    essentially seek a nonzero entry K(I,J) for which abs ( I - J )
c    is maximized.
c
c    The bandwidth M is defined in terms of the lower and upper bandwidths:
c
c      M = ML + 1 + MU
c
c    where
c
c      ML = maximum distance from any diagonal entry to a nonzero
c      entry in the same row, but earlier column,
c
c      MU = maximum distance from any diagonal entry to a nonzero
c      entry in the same row, but later column.
c
c    We assume the finite element variable adjacency relationship is 
c    symmetric, so we are guaranteed that ML = MU.
c
c    Note that the user is free to number the variables in any way
c    whatsoever, and to associate variables to nodes in any way,
c    so that some nodes have no variables, some have one, and some
c    have several.  
c
c    The storage of the indices of the variables is fairly simple.
c    In VAR, simply list all the variables associated with node 1, 
c    then all those associated with node 2, and so on.  Then set up
c    the pointer array VAR_NODE so that we can jump to the section of
c    VAR where the list begins for any particular node.
c
c    The routine does not check that each variable is only associated
c    with a single node.  This would normally be the case in a finite
c    element setting.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer ELEMENT_ORDER, the order of the elements.
c
c    Input, integer ELEMENT_NUM, the number of elements.
c
c    Input, integer ELEMENT_NODE(ELEMENT_ORDER,ELEMENT_NUM);
c    ELEMENT_NODE(I,J) is the global index of local node I in element J.
c
c    Input, integer NODE_NUM, the number of nodes.
c
c    Input, integer VAR_NODE(NODE_NUM+1), used to find the 
c    variables associated with a given node, which are in VAR in locations 
c    VAR_NODE(NODE) to VAR_NODE(NODE+1)-1.  Note that the last entry of
c    this array points to the location just after the last location in VAR.
c
c    Input, integer VAR_NUM, the number of variables.
c
c    Input, integer VAR(VAR_NUM), the indexes of the variables, 
c    which are presumably (but not necessarily) 1, 2, 3, ..., VAR_NUM.
c
c    Output, integer ML, MU, the lower and upper bandwidths of the 
c    matrix.
c
c    Output, integer M, the bandwidth of the matrix.
c
      implicit none

      integer element_num
      integer element_order
      integer node_num
      integer var_num

      integer element
      integer element_node(element_order,element_num)
      integer m
      integer ml
      integer mu
      integer node_global_i
      integer node_global_j
      integer node_local_i
      integer node_local_j
      integer var(var_num)
      integer var_global_i
      integer var_global_j
      integer var_local_i
      integer var_local_j
      integer var_node(node_num+1)

      ml = 0
      mu = 0

      do element = 1, element_num

        do node_local_i = 1, element_order
          node_global_i = element_node(node_local_i,element)

          do var_local_i = var_node(node_global_i), 
     &      var_node(node_global_i+1)-1
            var_global_i = var(var_local_i)

            do node_local_j = 1, element_order
              node_global_j = element_node(node_local_j,element)

              do var_local_j = var_node(node_global_j), 
     &          var_node(node_global_j+1)-1
                var_global_j = var(var_local_j)

                mu = max ( mu, var_global_j - var_global_i )
                ml = max ( ml, var_global_i - var_global_j )

              end do
            end do
          end do
        end do
      end do

      m = ml + 1 + mu

      return
      end
      subroutine basis_11_q4 ( q, i, p, phi, dphidx, dphidy )

c*********************************************************************72
c
cc BASIS_11_Q4: one basis at one point for a Q4 element.
c
c  Discussion:
c
c    The routine is given the coordinates of the vertices of a quadrilateral.
c    It works directly with these coordinates, and does not refer to a 
c    reference element.
c
c    The sides of the element are presumed to lie along coordinate axes.
c
c    The routine evaluates the basis functions, and their X and Y derivatives.
c
c  Physical Element Q4:
c
c    |
c    |  4-----3
c    |  |     |
c    |  |     |
c    Y  |     |
c    |  |     |
c    |  |     |
c    |  1-----2
c    |
c    +-----X------>
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision Q(2,4), the coordinates of the vertices.
c    It is common to list these points in counter clockwise order.
c
c    Input, integer I, the index of the basis function.
c
c    Input, double precision P(2), the evaluation point.
c
c    Output, double precision PHI(4), the basis functions 
c    at the evaluation points.
c
c    Output, double precision DPHIDX(4), DPHIDY(4), the basis
c    derivatives at the evaluation points.
c
c  Local Parameter:
c
c    Local, double precision AREA, the area of the rectangle.
c
      implicit none

      double precision area
      double precision dphidx(4)
      double precision dphidy(4)
      integer i
      double precision p(2)
      double precision phi
      double precision q(2,4)

      area = ( q(1,3) - q(1,1) ) * ( q(2,3) - q(2,1) )

      if ( i .eq. 1 ) then
        phi    =   ( q(1,3) - p(1) ) * ( q(2,3) - p(2) ) / area
        dphidx = -                     ( q(2,3) - p(2) ) / area
        dphidy = - ( q(1,3) - p(1) )                     / area
      else if ( i .eq. 2 ) then
        phi    =   ( p(1) - q(1,1) ) * ( q(2,3) - p(2) ) / area
        dphidx =                       ( q(2,3) - p(2) ) / area
        dphidy = - ( p(1) - q(1,1) )                     / area
      else if ( i .eq. 3 ) then
        phi    =   ( p(1) - q(1,1) ) * ( p(2) - q(2,1) ) / area
        dphidx =                       ( p(2) - q(2,1) ) / area
        dphidy =   ( p(1) - q(1,1) )                     / area
      else if ( i .eq. 4 ) then
        phi    =   ( q(1,3) - p(1) ) * ( p(2) - q(2,1) ) / area
        dphidx = -                     ( p(2) - q(2,1) ) / area
        dphidy =   ( q(1,3) - p(1) )                     / area
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BASIS_11_Q4 - Fatal error!'
        write ( *, '(a)' ) '  Illegal basis function index.'
        stop
      end if

      return
      end
      subroutine basis_11_q4_test ( )

c*********************************************************************72
c
cc BASIS_11_Q4_TEST verifies BASIS_11_Q4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 June 2009
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

      integer node_num
      parameter ( node_num = 4 )

      double precision dphidx(node_num,node_num)
      double precision dphidy(node_num,node_num)
      integer i
      integer j
      double precision phi(node_num,node_num)
      double precision q(2,node_num)
      double precision sum_x
      double precision sum_y

      save q

      data q /
     &  2.0D+00, 0.0D+00, 
     &  4.0D+00, 4.0D+00, 
     &  0.0D+00, 3.0D+00, 
     &  0.0D+00, 1.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BASIS_11_Q4_TEST:'
      write ( *, '(a)' ) '  Verify basis functions for element Q4.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of nodes = ', node_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Physical Nodes:'
      write ( *, '(a)' ) ' '
      do j = 1, node_num
        write ( *, '(2x,i8,2x,f7.3,2x,f7.3)' ) j, q(1:2,j)
      end do
     
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The basis function values at basis nodes'
      write ( *, '(a)' ) '  should form the identity matrix.'
      write ( *, '(a)' ) ' '

      do i = 1, node_num
        do j = 1, node_num
          call basis_11_q4 ( q, i, q(1,j), phi(i,j), dphidx(i,j), 
     &      dphidy(i,j) )
        end do
      end do

      do i = 1, node_num
        write ( *, '(2x,10f7.3)' ) phi(i,1:node_num)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The X and Y derivatives should sum to 0.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      dPhidX sum      dPhidY sum'
      write ( *, '(a)' ) ' '
      do j = 1, node_num
        sum_x = 0.0D+00
        sum_y = 0.0D+00
        do i = 1, node_num
          sum_x = sum_x + dphidx(i,j)
          sum_y = sum_y + dphidy(i,j)
        end do
        write ( *, '(2x,f14.8,2x,f14.8)' ) sum_x, sum_y
      end do

      return
      end
      subroutine basis_11_t3 ( t, i, p, qi, dqidx, dqidy )

c*********************************************************************72
c
cc BASIS_11_T3: one basis at one point for the T3 element.
c
c  Discussion:
c
c    The routine is given the coordinates of the nodes of a triangle. 
c        
c           3
c          / \
c         /   \
c        /     \
c       1-------2
c
c    It evaluates the linear basis function Q(I)(X,Y) associated with
c    node I, which has the property that it is a linear function
c    which is 1 at node I and zero at the other two nodes.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T(2,3), the coordinates of the nodes.
c
c    Input, integer I, the index of the desired basis function.
c    I should be between 1 and 3.
c
c    Input, double precision P(2), the coordinates of a point at which the
c    basis function is to be evaluated.
c
c    Output, double precision QI, DQIDX, DQIDY, the values of the basis
c    function and its X and Y derivatives.
c
      implicit none

      double precision area
      double precision dqidx
      double precision dqidy
      integer i
      integer i4_wrap
      integer ip1
      integer ip2
      double precision p(2)
      double precision qi
      double precision t(2,3)

      area = abs ( t(1,1) * ( t(2,2) - t(2,3) ) 
     &           + t(1,2) * ( t(2,3) - t(2,1) ) 
     &           + t(1,3) * ( t(2,1) - t(2,2) ) )

      if ( area .eq. 0.0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BASIS_11_T3 - Fatal error!'
        write ( *, '(a)' ) '  Element has zero area.'
        stop
      end if

      if ( i .lt. 1 .or. 3 .lt. i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BASIS_11_T3 - Fatal error!'
        write ( *, '(a)' ) '  Basis index I is not between 1 and 3.'
        write ( *, '(a,i8)' ) '  I = ', i
        stop
      end if

      ip1 = i4_wrap ( i + 1, 1, 3 )
      ip2 = i4_wrap ( i + 2, 1, 3 )

      qi = ( ( t(1,ip2) - t(1,ip1) ) * ( p(2) - t(2,ip1) ) 
     &     - ( t(2,ip2) - t(2,ip1) ) * ( p(1) - t(1,ip1) ) ) / area

      dqidx = - ( t(2,ip2) - t(2,ip1) ) / area
      dqidy =   ( t(1,ip2) - t(1,ip1) ) / area

      return
      end
      subroutine basis_11_t3_test ( )

c*********************************************************************72
c
cc BASIS_11_T3_TEST verifies BASIS_11_T3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 June 2009
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

      integer node_num
      parameter ( node_num = 3 )

      double precision dphidx(node_num,node_num)
      double precision dphidy(node_num,node_num)
      integer i
      integer j
      double precision phi(node_num,node_num)
      double precision sum_x
      double precision sum_y
      double precision t(2,node_num)

      save t

      data t /
     &  2.0D+00, 0.0D+00, 
     &  4.0D+00, 3.0D+00, 
     &  0.0D+00, 4.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BASIS_11_T3_TEST:'
      write ( *, '(a)' ) '  Verify basis functions for element T3.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of nodes = ', node_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Physical Nodes:'
      write ( *, '(a)' ) ' '
      do j = 1, node_num
        write ( *, '(2x,i8,2x,f7.3,2x,f7.3)' ) j, t(1,j), t(2,j)
      end do
     
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The basis function values at basis nodes'
      write ( *, '(a)' ) '  should form the identity matrix.'
      write ( *, '(a)' ) ' '

      do i = 1, node_num
        do j = 1, node_num
          call basis_11_t3 ( t, i, t(1,j), phi(i,j), dphidx(i,j), 
     &      dphidy(i,j) )
        end do
      end do

      do i = 1, node_num
        write ( *, '(2x,10f7.3)' ) ( phi(i,j), j = 1, node_num )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The X and Y derivatives should sum to 0.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      dPhidX sum    dPhidY sum'
      write ( *, '(a)' ) ' '

      do j = 1, node_num

        sum_x = 0.0D+00
        sum_y = 0.0D+00
        do i = 1, node_num
          sum_x = sum_x + dphidx(i,j)
          sum_y = sum_y + dphidy(i,j)
        end do
        write ( *, '(2x,f14.8,f14.8)' ) sum_x, sum_y

      end do

      return
      end
      subroutine basis_11_t4 ( t, i, p, phi, dphidx, dphidy )

c*********************************************************************72
c
cc BASIS_11_T4: one basis at one point for a T4 element.
c
c  Discussion:
c
c    The T4 element is the cubic bubble triangle.
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
c  Physical Element T4: 
c        
c            3
c           / \
c          /   \
c         /  4  \
c        /       \
c       1---------2
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T(2,4), the coordinates of the vertices
c    of the triangle, and the coordinates of the centroid.  
c    It is common to list the first three points in counter clockwise
c    order.
c
c    Input, integer I, the index of the basis function.
c
c    Input, double precision P(2), the points where the basis function
c    is to be evaluated.
c
c    Output, double precision PHI, the value of the basis function
c    at the evaluation point.
c
c    Output, double precision DPHIDX, DPHIDY, the value of the 
c    derivatives at the evaluation point.
c
c  Local parameters:
c
c    Local, double precision AREA, is (twice) the area of the triangle.
c
      implicit none

      double precision area
      double precision dphidx
      double precision dphidy
      double precision dpsidx(4)
      double precision dpsidy(4)
      integer i
      integer j
      double precision p(2)
      double precision phi
      double precision psi(4)
      double precision t(2,4)

      area = t(1,1) * ( t(2,2) - t(2,3) ) 
     &     + t(1,2) * ( t(2,3) - t(2,1) ) 
     &     + t(1,3) * ( t(2,1) - t(2,2) )

      psi(1) =     (   ( t(1,3) - t(1,2) ) * ( p(2) - t(2,2) )     
     &               - ( t(2,3) - t(2,2) ) * ( p(1) - t(1,2) ) ) 
      dpsidx(1) =    - ( t(2,3) - t(2,2) )
      dpsidy(1) =      ( t(1,3) - t(1,2) )

      psi(2) =     (   ( t(1,1) - t(1,3) ) * ( p(2) - t(2,3) )     
     &               - ( t(2,1) - t(2,3) ) * ( p(1) - t(1,3) ) )
      dpsidx(2) =    - ( t(2,1) - t(2,3) )
      dpsidy(2) =      ( t(1,1) - t(1,3) )

      psi(3) =     (   ( t(1,2) - t(1,1) ) * ( p(2) - t(2,1) )     
     &               - ( t(2,2) - t(2,1) ) * ( p(1) - t(1,1) ) )
      dpsidx(3) =    - ( t(2,2) - t(2,1) )
      dpsidy(3) =      ( t(1,2) - t(1,1) )
c
c  Normalize the first three functions.
c
      psi(1:3)    =    psi(1:3) / area
      dpsidx(1:3) = dpsidx(1:3) / area
      dpsidy(1:3) = dpsidy(1:3) / area
c
c  Compute the cubic bubble function.
c
      psi(4) = 27.0D+00 * psi(1) * psi(2) * psi(3)

      dpsidx(4) = 27.0D+00 * ( 
     &              dpsidx(1) *    psi(2) *    psi(3) 
     &              +  psi(1) * dpsidx(2) *    psi(3) 
     &              +  psi(1) *    psi(2) * dpsidx(3) )

      dpsidy(4) = 27.0D+00 * ( 
     &              dpsidy(1) *    psi(2) *    psi(3) 
     &              +  psi(1) * dpsidy(2) *    psi(3) 
     &              +  psi(1) *    psi(2) * dpsidy(3) )
c
c  Subtract 1/3 of the cubic bubble function from each of the three linears.
c
      do j = 1, 3
        psi(j)    =    psi(j) -    psi(4) / 3.0D+00
        dpsidx(j) = dpsidx(j) - dpsidx(4) / 3.0D+00
        dpsidy(j) = dpsidy(j) - dpsidy(4) / 3.0D+00
      end do

      phi    = psi(i)
      dphidx = dpsidx(i)
      dphidy = dpsidy(i)

      return
      end
      subroutine basis_11_t4_test ( )

c*********************************************************************72
c
cc BASIS_11_T4_TEST verifies BASIS_11_T4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 June 2009
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

      integer node_num
      parameter ( node_num = 4 )

      double precision dphidx(node_num,node_num)
      double precision dphidy(node_num,node_num)
      integer i
      integer j
      double precision phi(node_num,node_num)
      double precision sum_x
      double precision sum_y
      double precision t(2,node_num)

      save t

      data t /
     &  2.0D+00, 0.0D+00, 
     &  4.0D+00, 3.0D+00, 
     &  0.0D+00, 4.0D+00, 
     &  0.0D+00, 0.0D+00 /
c
c  The node associated with the fourth basis function is the centroid.
c
      t(1,4) = ( t(1,1) + t(1,2) + t(1,3) ) / 3.0D+00
      t(2,4) = ( t(2,1) + t(2,2) + t(2,3) ) / 3.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BASIS_11_T4_TEST:'
      write ( *, '(a)' ) '  Verify basis functions for element T4.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of nodes = ', node_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Physical Nodes:'
      write ( *, '(a)' ) ' '
      do j = 1, node_num
        write ( *, '(2x,i8,2x,f7.3,2x,f7.3)' ) j, t(1,j), t(2,j)
      end do
     
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The basis function values at basis nodes'
      write ( *, '(a)' ) '  should form the identity matrix.'
      write ( *, '(a)' ) ' '

      do i = 1, node_num
        do j = 1, node_num
          call basis_11_t4 ( t, i, t(1,j), phi(i,j), dphidx(i,j), 
     &      dphidy(i,j) )
        end do
      end do

      do i = 1, node_num
        write ( *, '(2x,10f7.3)' ) ( phi(i,j), j = 1, node_num )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The X and Y derivatives should sum to 0.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      dPhidX sum    dPhidY sum'
      write ( *, '(a)' ) ' '

      do j = 1, node_num

        sum_x = 0.0D+00
        sum_y = 0.0D+00
        do i = 1, node_num
          sum_x = sum_x + dphidx(i,j)
          sum_y = sum_y + dphidy(i,j)
        end do
        write ( *, '(2x,f14.8,f14.8)' ) sum_x, sum_y

      end do

      return
      end
      subroutine basis_11_t6 ( t, i, p, bi, dbidx, dbidy )

c*********************************************************************72
c
cc BASIS_11_T6: one basis at one point for the T6 element.
c
c  Discussion:
c
c    The routine is given the coordinates of the nodes of a triangle. 
c        
c           3
c          / \
c         6   5
c        /     \
c       1---4---2
c
c    It evaluates the quadratic basis function B(I)(X,Y) associated with
c    node I, which has the property that it is a quadratic function
c    which is 1 at node I and zero at the other five nodes.
c
c    This routine assumes that the sides of the triangle are straight,
c    so that the midside nodes fall on the line between two vertices.
c
c    This routine relies on the fact that each basis function can be
c    written as the product of two linear factors, which are easily
c    computed and normalized.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T(2,6), the coordinates of the nodes.
c
c    Input, integer I, the index of the desired basis function.
c    I should be between 1 and 6.
c
c    Input, double precision P(2), the coordinates of a point at which the
c    basis function is to be evaluated.
c
c    Output, double precision BI, DBIDX, DBIDY, the values of the basis
c    function and its X and Y derivatives.
c
      implicit none

      double precision bi
      double precision dbidx
      double precision dbidy
      double precision gf
      double precision gn
      double precision hf
      double precision hn
      integer i
      integer i4_wrap
      integer j1
      integer j2
      integer k1
      integer k2
      double precision p(2)
      double precision t(2,6)

      if ( i .lt. 1 .or. 6 .lt. i ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'BASIS_11_T6 - Fatal error!'
        write ( *, '(a)' ) '  Basis index I is not between 1 and 6.'
        write ( *, '(a,i8)' ) '  I = ', i
        stop
      end if
c
c  Determine the pairs of nodes.
c
      if ( i .le. 3 ) then
        j1 = i4_wrap ( i + 1, 1, 3 )
        j2 = i4_wrap ( i + 2, 1, 3 )
        k1 = i + 3
        k2 = i4_wrap ( i + 5, 4, 6 )
      else
        j1 = i - 3
        j2 = i4_wrap ( i - 3 + 2, 1, 3 )
        k1 = i4_wrap ( i - 3 + 1, 1, 3 )
        k2 = i4_wrap ( i - 3 + 2, 1, 3 )
      end if
c
c  Evaluate the two linear factors GF and HF, 
c  and their normalizers GN and HN.
c
      gf = ( p(1)    - t(1,j1) ) * ( t(2,j2) - t(2,j1) ) 
     &   - ( t(1,j2) - t(1,j1) ) * ( p(2)    - t(2,j1) ) 

      gn = ( t(1,i)  - t(1,j1) ) * ( t(2,j2) - t(2,j1) ) 
     &   - ( t(1,j2) - t(1,j1) ) * ( t(2,i)  - t(2,j1) )   

      hf = ( p(1)    - t(1,k1) ) * ( t(2,k2) - t(2,k1) ) 
     &   - ( t(1,k2) - t(1,k1) ) * ( p(2)    - t(2,k1) ) 

      hn = ( t(1,i)  - t(1,k1) ) * ( t(2,k2) - t(2,k1) ) 
     &   - ( t(1,k2) - t(1,k1) ) * ( t(2,i)  - t(2,k1) ) 
c
c  Construct the basis function and its derivatives.
c
      bi =        ( gf                  / gn ) 
     &          * ( hf                  / hn )

      dbidx =   ( ( t(2,j2) - t(2,j1) ) / gn ) 
     &          * ( hf                  / hn ) 
     &          + ( gf                  / gn ) 
     &        * ( ( t(2,k2) - t(2,k1) ) / hn )

      dbidy = - ( ( t(1,j2) - t(1,j1) ) / gn ) 
     &        * (   hf                  / hn ) 
     &        - (   gf                  / gn ) 
     &        * ( ( t(1,k2) - t(1,k1) ) / hn )

      return
      end
      subroutine basis_11_t6_test ( )

c*********************************************************************72
c
cc BASIS_11_T6_TEST verifies BASIS_11_T6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 June 2009
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

      integer node_num
      parameter ( node_num = 6 )

      double precision dphidx(node_num,node_num)
      double precision dphidy(node_num,node_num)
      integer i
      integer j
      double precision phi(node_num,node_num)
      double precision sum_x
      double precision sum_y
      double precision t(2,node_num)

      save t

      data t /
     &  2.0D+00, 0.0D+00, 
     &  4.0D+00, 3.0D+00, 
     &  0.0D+00, 4.0D+00, 
     &  3.0D+00, 1.5D+00, 
     &  2.0D+00, 3.5D+00, 
     &  1.0D+00, 2.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BASIS_11_T6_TEST:'
      write ( *, '(a)' ) '  Verify basis functions for element T6.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of nodes = ', node_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Physical Nodes:'
      write ( *, '(a)' ) ' '
      do j = 1, node_num
        write ( *, '(2x,i8,2x,f7.3,2x,f7.3)' ) j, t(1:2,j)
      end do
     
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The basis function values at basis nodes'
      write ( *, '(a)' ) '  should form the identity matrix.'
      write ( *, '(a)' ) ' '

      do i = 1, node_num
        do j = 1, node_num
          call basis_11_t6 ( t, i, t(1,j), phi(i,j), dphidx(i,j), 
     &      dphidy(i,j) )
        end do
      end do

      do i = 1, node_num
        write ( *, '(2x,10f7.3)' ) ( phi(i,j), j = 1, node_num)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The X and Y derivatives should sum to 0.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      dPhidX sum      dPhidY sum'
      write ( *, '(a)' ) ' '
      do j = 1, node_num
        sum_x = 0.0D+00
        sum_y = 0.0D+00
        do i = 1, node_num
          sum_x = sum_x + dphidx(i,j)
          sum_y = sum_y + dphidy(i,j)
        end do
        write ( *, '(2x,f14.8,2x,f14.8)' ) sum_x, sum_y
      end do

      return
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
      function degrees_to_radians ( degrees )

c*********************************************************************72
c
cc DEGREES_TO_RADIANS converts an angle measure from degrees to radians.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision DEGREES, the angle measure in degrees.
c
c    Output, double precision DEGREES_TO_RADIANS, the angle measure in radians.
c
      implicit none

      double precision degrees
      double precision degrees_to_radians
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )

      degrees_to_radians = ( degrees / 180.0D+00 ) * pi

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
      function element_code ( i )

c*********************************************************************72
c
cc ELEMENT_CODE returns the code for each element.
c
c  Discussion:
c
c     I  ELEMENT_CODE   Definition
c     -  ------------   ----------
c     1  Q4             4 node linear Lagrange/serendipity quadrilateral;
c     2  Q8             8 node quadratic serendipity quadrilateral;
c     3  Q9             9 node quadratic Lagrange quadrilateral;
c     4  Q12            12 node cubic serendipity quadrilateral;
c     5  Q16            16 node cubic Lagrange quadrilateral;
c     6  QL             6 node linear/quadratic quadrilateral;
c     7  T3             3 node linear triangle;
c     8  T4             4 node cubic bubble triangle
c     9  T6             6 node quadratic triangle;
c    10  T10            10 node cubic triangle.
c 
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    12 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, the index of the element type.  
c
c    Output, character * ( 3 ) ELEMENT_CODE, the code for the element type.
c
      implicit none

      character * ( 3 ) element_code
      integer i

      if ( i .eq. 1 ) then
        element_code = 'Q4'
      else if ( i .eq. 2 ) then
        element_code = 'Q8'
      else if ( i .eq. 3 ) then
        element_code = 'Q9'
      else if ( i .eq. 4 ) then
        element_code = 'Q12'
      else if ( i .eq. 5 ) then
        element_code = 'Q16'
      else if ( i .eq. 6 ) then
        element_code = 'QL'
      else if ( i .eq. 7 ) then
        element_code = 'T3'
      else if ( i .eq. 8 ) then
        element_code = 'T4'
      else if ( i .eq. 9 ) then
        element_code = 'T6'
      else if ( i .eq. 10 ) then
        element_code = 'T10'
      else
        element_code = '???'
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
      subroutine itable_write0 ( output_file_name, m, n, table )

c*********************************************************************72
c
cc ITABLE_WRITE0 writes an ITABLE file with no header.
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
c    Input, integer TABLE(M,N), the table data.
c
      implicit none

      integer m
      integer n

      integer j
      character * ( * ) output_file_name
      integer output_unit
      character * ( 30 ) string
      integer table(m,n)
c
c  Open the file.
c
      call get_unit ( output_unit )

      open ( unit = output_unit, file = output_file_name, 
     &  status = 'replace' )
c
c  Create the format string.
c
      write ( string, '(a1,i8,a4)' ) '(', m, 'i10)'
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
      function r8_huge ( )

c*********************************************************************72
c
cc R8_HUGE returns a "huge" R8.
c
c  Discussion:
c
c    The value returned by this function is NOT required to be the
c    maximum representable R8.  This value varies from machine to machine,
c    from compiler to compiler, and may cause problems when being printed.
c    We simply want a "very large" but non-infinite number.
c
c    FORTRAN90 provides a built-in routine HUGE ( X ) that
c    can return the maximum representable number of the same datatype
c    as X, if that is what is really desired.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 April 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision R8_HUGE, a huge number.
c
      implicit none

      double precision r8_huge

      r8_huge = 1.0D+30

      return
      end
      function r8_power ( r, p )

c*********************************************************************72
c
cc R8_POWER computes the P-th power of an R8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision R, the base.
c
c    Input, integer P, the power, which may be negative.
c
c    Output, double precision R8_POWER, the value of the P-th power of R.
c
      implicit none

      integer p
      double precision r
      double precision r8_power
      double precision value
c
c  Special case.  R^0 = 1.
c
      if ( p .eq. 0 ) then

        value = 1.0D+00
c
c  Special case.  Positive powers of 0 are 0.
c  For negative powers of 0, we go ahead and compute R**P,
c  relying on the software to complain.
c
      else if ( r .eq. 0.0D+00 ) then

        if ( 0 .lt. p ) then
          value = 0.0D+00
        else
          value = r**p
        end if

      else if ( 1 .le. p ) then
        value = r**p
      else
        value = 1.0D+00 / r**(-p)
      end if

      r8_power = value

      return
      end
      subroutine r8_swap ( x, y )

c*********************************************************************72
c
cc R8_SWAP switches two R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 November 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, double precision X, Y.  On output, the values of X and
c    Y have been interchanged.
c
      implicit none

      double precision x
      double precision y
      double precision z

      z = x
      x = y
      y = z

      return
      end
      function r8_uniform_01 ( seed )

c*********************************************************************72
c
cc R8_UNIFORM_01 returns a unit pseudorandom R8.
c
c  Discussion:
c
c    This routine implements the recursion
c
c      seed = 16807 * seed mod ( 2**31 - 1 )
c      r8_uniform_01 = seed / ( 2**31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      R8_UNIFORM_01
c      SEED      SEED
c
c         12345   207482415  0.096616
c     207482415  1790989824  0.833995
c    1790989824  2035175616  0.947702
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
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley Interscience, page 95, 1998.
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
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R8_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      double precision r8_uniform_01
      integer k
      integer seed

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + 2147483647
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r8_uniform_01 = dble ( seed ) * 4.656612875D-10

      return
      end
      subroutine r8mat_print ( m, n, a, title )

c*********************************************************************72
c
cc R8MAT_PRINT prints an R8MAT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 May 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of rows in A.
c
c    Input, integer N, the number of columns in A.
c
c    Input, double precision A(M,N), the matrix.
c
c    Input, character ( len = * ) TITLE, a title to be printed.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      character ( len = * ) title

      call r8mat_print_some ( m, n, a, 1, 1, m, n, title )

      return
      end
      subroutine r8mat_print_some ( m, n, a, ilo, jlo, ihi, jhi,
     &  title )

c*********************************************************************72
c
cc R8MAT_PRINT_SOME prints some of an R8MAT.
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
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, double precision A(M,N), an M by N matrix to be printed.
c
c    Input, integer ILO, JLO, the first row and column to print.
c
c    Input, integer IHI, JHI, the last row and column to print.
c
c    Input, character ( len = * ) TITLE, an optional title.
c
      implicit none

      integer incx
      parameter ( incx = 5 )
      integer m
      integer n

      double precision a(m,n)
      character * ( 14 ) ctemp(incx)
      integer i
      integer i2hi
      integer i2lo
      integer ihi
      integer ilo
      integer inc
      integer j
      integer j2
      integer j2hi
      integer j2lo
      integer jhi
      integer jlo
      character * ( * ) title
      integer title_length

      title_length = len_trim ( title )

      if ( 0 .lt. title_length ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) title(1:title_length)
      end if

      do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

        j2hi = j2lo + incx - 1
        j2hi = min ( j2hi, n )
        j2hi = min ( j2hi, jhi )

        inc = j2hi + 1 - j2lo

        write ( *, '(a)' ) ' '

        do j = j2lo, j2hi
          j2 = j + 1 - j2lo
          write ( ctemp(j2), '(i7,7x)') j
        end do

        write ( *, '(''  Col   '',5a14)' ) ( ctemp(j), j = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) ' '

        i2lo = max ( ilo, 1 )
        i2hi = min ( ihi, m )

        do i = i2lo, i2hi

          do j2 = 1, inc

            j = j2lo - 1 + j2

            write ( ctemp(j2), '(g14.6)' ) a(i,j)

          end do

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

        end do

      end do

      write ( *, '(a)' ) ' '

      return
      end
      subroutine reference_to_physical_t3 ( t, n, ref, phy )

c*********************************************************************72
c
cc REFERENCE_TO_PHYSICAL_T3 maps T3 reference points to physical points.
c
c  Discussion:
c
c    Given the vertices of an order 3 physical triangle and a point
c    (XSI,ETA) in the reference triangle, the routine computes the value
c    of the corresponding image point (X,Y) in physical space.
c
c    This routine is also appropriate for an order 4 triangle,
c    as long as the fourth node is the centroid of the triangle.
c
c    This routine may also be appropriate for an order 6
c    triangle, if the mapping between reference and physical space
c    is linear.  This implies, in particular, that the sides of the
c    image triangle are straight and that the "midside" nodes in the
c    physical triangle are halfway along the sides of
c    the physical triangle.
c
c  Reference Element T3:
c
c    |
c    1  3
c    |  |\
c    |  | \
c    S  |  \
c    |  |   \
c    |  |    \
c    0  1-----2
c    |
c    +--0--R--1-->
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T(2,3), the coordinates of the vertices.
c    The vertices are assumed to be the images of (0,0), (1,0) and
c    (0,1) respectively.
c
c    Input, integer N, the number of points to transform.
c
c    Input, double precision REF(2,N), points in the reference triangle.
c
c    Output, double precision PHY(2,N), corresponding points in the
c    physical triangle.
c
      implicit none

      integer n

      integer i
      integer j
      double precision phy(2,n)
      double precision ref(2,n)
      double precision t(2,3)

      do i = 1, 2
        do j = 1, n
          phy(i,j) = t(i,1) * ( 1.0D+00 - ref(1,j) - ref(2,j) )
     &             + t(i,2) *             ref(1,j)
     &             + t(i,3) *                        ref(2,j)
        end do
      end do

      return
      end
      subroutine reference_to_physical_t6 ( t, n, ref, phy )

c*********************************************************************72
c
cc REFERENCE_TO_PHYSICAL_T6 maps T6 reference points to physical points.
c
c  Discussion:
c
c    Given the vertices of an order 6 physical triangle and a point
c    (XSI,ETA) in the reference triangle, the routine computes the value
c    of the corresponding image point (X,Y) in physical space.
c
c    The mapping from (XSI,ETA) to (X,Y) has the form:
c
c      X(ETA,XSI) = A1 * XSI**2 + B1 * XSI*ETA + C1 * ETA**2
c                 + D1 * XSI    + E1 * ETA     + F1
c
c      Y(ETA,XSI) = A2 * XSI**2 + B2 * XSI*ETA + C2 * ETA**2
c                 + D2 * XSI    + E2 * ETA     + F2
c
c  Reference Element T6:
c
c    |
c    1  3
c    |  |\
c    |  | \
c    S  6  5
c    |  |   \
c    |  |    \
c    0  1--4--2
c    |
c    +--0--R--1-->
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 June 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision T(2,6), the coordinates of the vertices.
c    The vertices are assumed to be the images of (0,0), (1,0),
c    (0,1),(1/2,0), (1/2,1/2) and (0,1/2) respectively.
c
c    Input, integer N, the number of points to transform.
c
c    Input, double precision REF(2,N), points in the reference triangle.
c
c    Output, double precision PHY(2,N), corresponding points in the
c    physical triangle.
c
      implicit none

      integer n

      double precision a(2)
      double precision b(2)
      double precision c(2)
      double precision d(2)
      double precision e(2)
      double precision f(2)
      integer i
      integer j
      double precision phy(2,n)
      double precision ref(2,n)
      double precision t(2,6)

      do i = 1, 2

        a(i) =   2.0D+00 * t(i,1) + 2.0D+00 * t(i,2)                    
     &         - 4.0D+00 * t(i,4)

        b(i) =   4.0D+00 * t(i,1)                                       
     &         - 4.0D+00 * t(i,4) + 4.0D+00 * t(i,5) - 4.0D+00 * t(i,6)

        c(i) =   2.0D+00 * t(i,1)                    + 2.0D+00 * t(i,3) 
     &                                               - 4.0D+00 * t(i,6)

        d(i) = - 3.0D+00 * t(i,1) -           t(i,2)                    
     &         + 4.0D+00 * t(i,4)

        e(i) = - 3.0D+00 * t(i,1)                    -           t(i,3) 
     &                                               + 4.0D+00 * t(i,6)
        f(i) =             t(i,1)

      end do

      do i = 1, 2
        do j = 1, n
          phy(i,j) = a(i) * ref(1,j) * ref(1,j) 
     &             + b(i) * ref(1,j) * ref(2,j) 
     &             + c(i) * ref(2,j) * ref(2,j) 
     &             + d(i) * ref(1,j)
     &             + e(i) * ref(2,j)
     &             + f(i)
        end do
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
