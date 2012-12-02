      program main

c*********************************************************************72
c
cc MAIN is the main routine of the finite element program FEM2D_POISSON_RECTANGLE_LINEAR.
c
c  Discussion:
c
c    This program solves
c
c      - d2U(X,Y)/dx2 - d2U(X,Y)/dy2 = F(X,Y)
c
c    in a rectangular region in the plane.
c
c    Along the boundary of the region, Dirichlet conditions
c    are imposed:
c
c      U(X,Y) = G(X,Y)
c
c    The code uses continuous piecewise linear basis functions on
c    triangles determined by a uniform grid of NX by NY points.
c
c    u    =      sin ( pi * x ) * sin ( pi * y ) + x
c
c    dudx = pi * cos ( pi * x ) * sin ( pi * y ) + 1
c    dudy = pi * sin ( pi * x ) * cos ( pi * y )
c
c    d2udx2 = - pi * pi * sin ( pi * x ) * sin ( pi * y )
c    d2udy2 = - pi * pi * sin ( pi * x ) * sin ( pi * y )
c
c    rhs  = 2 * pi * pi * sin ( pi * x ) * sin ( pi * y )
c
c  THINGS YOU CAN EASILY CHANGE:
c
c    1) Change NX or NY, the number of nodes in the X and Y directions.
c    2) Change XL, XR, YB, YT, the left, right, bottom and top limits of the rectangle.
c    3) Change the exact solution in the EXACT routine, but make sure you also
c       modify the formula for RHS in the assembly portion of the programc
c
c  HARDER TO CHANGE:
c
c    4) Change from "linear" to "quadratic" triangles;
c    5) Change the region from a rectangle to a general triangulated region;
c    6) Store the matrix as a sparse matrix so you can solve bigger systems.
c    7) Handle Neumann boundary conditions.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 November 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nx
      parameter ( nx = 17 )
      integer ny
      parameter ( ny = 17 )

      integer node_num
      parameter ( node_num = nx * ny )

      integer element_num
      parameter ( element_num = 2 * ( nx - 1 ) * ( ny - 1 ) )

      double precision a(node_num,node_num)
      double precision area
      double precision b(node_num)
      double precision dqidx
      double precision dqidy
      double precision dqjdx
      double precision dqjdy
      double precision dudx
      double precision dudy
      integer e
      integer element_node(3,element_num)
      integer i
      integer i1
      integer i2
      integer i3
      integer info
      integer j
      integer k
      integer k2
      integer nq1
      integer nq2
      integer nti1
      integer nti2
      integer nti3
      integer ntj1
      integer ntj2
      integer ntj3
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      integer q1
      integer q2
      double precision qi
      double precision qj
      integer ti1
      integer ti2
      integer ti3
      integer tj1
      integer tj2
      integer tj3
      double precision rhs
      double precision u
      double precision wq
      double precision x(node_num)
      double precision xl
      parameter ( xl = 0.0D+00 )
      double precision xq
      double precision xr
      parameter ( xr = 1.0D+00 )
      double precision y(node_num)
      double precision yb
      parameter ( yb = 0.0D+00 )
      double precision yq
      double precision yt
      parameter ( yt = 1.0D+00 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM2D_POISSON_RECTANGLE_LINEAR'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Solution of the Poisson equation:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  - Uxx - Uyy = F(x,y) inside the region,'
      write ( *, '(a)' ) '       U(x,y) = G(x,y) on the boundary.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The region is a rectangle, defined by:'
      write ( *, '(a)' ) ' '
      write ( *, '(g14.6,a,g14.6)' ) xl,' = XL<= X <= XR = ', xr
      write ( *, '(g14.6,a,g14.6)' ) yb,' = YB<= Y <= YT = ', yt
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The finite element method is used, with'
      write ( *, '(a)' ) '  piecewise linear basis functions on 3 node'
      write ( *, '(a)' ) '  triangular elements.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The vertices of the triangles are generated'
      write ( *, '(a)' ) '  by an underlying grid whose dimensions are'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  NX =                       ', nx
      write ( *, '(a,i8)' ) '  NY =                       ', ny
c
c  NODE COORDINATES
c
c  Numbering of nodes is suggested by the following 5x10 example:
c
c    J=5 | K=41  K=42 ... K=50
c    ... |
c    J=2 | K=11  K=12 ... K=20
c    J=1 | K= 1  K= 2     K=10
c        +--------------------
c          I= 1  I= 2 ... I=10
c
      write ( *, '(a,i8)' ) '  Number of nodes =          ', node_num

      k = 0
      do j = 1, ny
        do i = 1, nx

          k = k + 1

          x(k) = ( dble ( nx - i     ) * xl
     &           + dble (      i - 1 ) * xr )
     &           / dble ( nx     - 1 )

          y(k) = ( dble ( ny - j     ) * yb
     &           + dble (      j - 1 ) * yt )
     &           / dble ( ny     - 1 )

        end do
      end do
c
c  ELEMENT array
c
c  Organize the nodes into a grid of 3-node triangles.
c  Here is part of the diagram for a 5x10 example:
c
c    |  \ |  \ |  \ |
c    |   \|   \|   \|
c   21---22---23---24--
c    |\ 8 |\10 |\12 |
c    | \  | \  | \  |
c    |  \ |  \ |  \ |  \ |
c    |  7\|  9\| 11\|   \|
c   11---12---13---14---15---16---17---18---19---20
c    |\ 2 |\ 4 |\ 6 |\  8|                   |\ 18|
c    | \  | \  | \  | \  |                   | \  |
c    |  \ |  \ |  \ |  \ |      ...          |  \ |
c    |  1\|  3\|  5\| 7 \|                   |17 \|
c    1----2----3----4----5----6----7----8----9---10
c
      write ( *, '(a,i8)' ) '  Number of elements =       ', element_num

      k = 0

      do j = 1, ny - 1
        do i = 1, nx - 1

          k = k + 1
          element_node(1,k) = i     + ( j - 1 ) * nx
          element_node(2,k) = i + 1 + ( j - 1 ) * nx
          element_node(3,k) = i     +   j       * nx

          k = k + 1
          element_node(1,k) = i + 1 +   j       * nx
          element_node(2,k) = i     +   j       * nx
          element_node(3,k) = i + 1 + ( j - 1 ) * nx

        end do
      end do
c
c  ASSEMBLE THE SYSTEM
c
c  Assemble the coefficient matrix A and the right-hand side B of the
c  finite element equations, ignoring boundary conditions.
c
      do i = 1, node_num
        b(i) = 0.0D+00
      end do

      do j = 1, node_num
        do i = 1, node_num
          a(i,j) = 0.0D+00
        end do
      end do

      do e = 1, element_num

        i1 = element_node(1,e)
        i2 = element_node(2,e)
        i3 = element_node(3,e)

        area = 0.5D+00 *
     &    ( x(i1) * ( y(i2) - y(i3) )
     &    + x(i2) * ( y(i3) - y(i1) )
     &    + x(i3) * ( y(i1) - y(i2) ) )
c
c  Consider each quadrature point.
c  Here, we use the midside nodes as quadrature points.
c
        do q1 = 1, 3

          q2 = mod ( q1, 3 ) + 1

          nq1 = element_node(q1,e)
          nq2 = element_node(q2,e)

          xq = 0.5D+00 * ( x(nq1) + x(nq2) )
          yq = 0.5D+00 * ( y(nq1) + y(nq2) )
          wq = 1.0D+00 / 3.0D+00
c
c  Consider each test function in the element.
c
          do ti1 = 1, 3

            ti2 = mod ( ti1,     3 ) + 1
            ti3 = mod ( ti1 + 1, 3 ) + 1

            nti1 = element_node(ti1,e)
            nti2 = element_node(ti2,e)
            nti3 = element_node(ti3,e)

            qi = 0.5D+00 * (
     &          ( x(nti3) - x(nti2) ) * ( yq - y(nti2) )
     &        - ( y(nti3) - y(nti2) ) * ( xq - x(nti2) ) ) / area
            dqidx = - 0.5D+00 * ( y(nti3) - y(nti2) ) / area
            dqidy =   0.5D+00 * ( x(nti3) - x(nti2) ) / area

            rhs = 2.0D+00 * pi * pi * sin ( pi * xq ) * sin ( pi * yq )

            b(nti1) = b(nti1) + area * wq * rhs * qi
c
c  Consider each basis function in the element.
c
            do tj1 = 1, 3

              tj2 = mod ( tj1,     3 ) + 1
              tj3 = mod ( tj1 + 1, 3 ) + 1

              ntj1 = element_node(tj1,e)
              ntj2 = element_node(tj2,e)
              ntj3 = element_node(tj3,e)

              qj = 0.5D+00 * (
     &            ( x(ntj3) - x(ntj2) ) * ( yq - y(ntj2) )
     &          - ( y(ntj3) - y(ntj2) ) * ( xq - x(ntj2) ) ) / area
              dqjdx = - 0.5D+00 * ( y(ntj3) - y(ntj2) ) / area
              dqjdy =   0.5D+00 * ( x(ntj3) - x(ntj2) ) / area

              a(nti1,ntj1) = a(nti1,ntj1)
     &          + area * wq * ( dqidx * dqjdx + dqidy * dqjdy )

            end do

          end do

        end do

      end do
c
c  BOUNDARY CONDITIONS
c
c  If the K-th variable is at a boundary node, replace the K-th finite
c  element equation by a boundary condition that sets the variable to U(K).
c
      k = 0

      do j = 1, ny

        do i = 1, nx

          k = k + 1

          if ( i .eq. 1 .or.
     &         i .eq. nx .or.
     &         j .eq. 1 .or.
     &         j .eq. ny ) then

            call exact ( x(k), y(k), u, dudx, dudy )

            do k2 = 1, node_num
              a(k,k2) = 0.0D+00
            end do
            a(k,k)          = 1.0D+00
            b(k)            = u

          end if
        end do
      end do
c
c  SOLVE the linear system A * X = B.
c
c  The solution X is actually returned in the space occupied by B.
c
      call r8ge_fs ( node_num, a, b, info )

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' )
     &    'FEM2D_POISSON_RECTANGLE_LINEAR - Fatal error!'
        write ( *, '(a)' ) '  R8GEFS returned an error condition.'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  The linear system was not solved, and the'
        write ( *, '(a)' ) '  algorithm cannot proceed.'
        stop
      end if
c
c  COMPARE computed and exact solutions.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     K     I     J          X           ' //
     &  'Y        U               U               Error'
      write ( *, '(a)' ) '                                        ' //
     &  '         exact           computed '
      write ( *, '(a)' ) ' '

      k = 0

      do j = 1, ny
        do i = 1, nx

          k = k + 1

          call exact ( x(k), y(k), u, dudx, dudy )
          write ( *, '(i6,i6,i6,f12.2,f12.2,g16.6,g16.6,g16.6)' )
     &      k, i, j, x(k), y(k), u, b(k), abs ( u - b(k) )

        end do
        write ( *, '(a)' ) ' '
      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM2D_POISSON_RECTANGLE_LINEAR:'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
      subroutine r8ge_fs ( n, a, b, info )

c*********************************************************************72
c
cc R8GE_FS factors and solves a R8GE system.
c
c  Discussion:
c
c    The R8GE storage format is used for a general M by N matrix.  A storage
c    space is made for each entry.  The two dimensional logical
c    array can be thought of as a vector of M*N entries, starting with
c    the M entries in the column 1, then the M entries in column 2
c    and so on.  Considered as a vector, the entry A(I,J) is then stored
c    in vector location I+(J-1)*M.
c
c    R8GE storage is used by LINPACK and LAPACK.
c
c    R8GE_FS does not save the LU factors of the matrix, and hence cannot
c    be used to efficiently solve multiple linear systems, or even to
c    factor A at one time, and solve a single linear system at a later time.
c
c    R8GE_FS uses partial pivoting, but no pivot vector is required.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 November 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input/output, double precision A(N,N).
c    On input, A is the coefficient matrix of the linear system.
c    On output, A is in unit upper triangular form, and
c    represents the U factor of an LU factorization of the
c    original coefficient matrix.
c
c    Input/output, double precision B(N).
c    On input, B is the right hand side of the linear system.
c    On output, B is the solution of the linear system.
c
c    Output, integer INFO, singularity flag.
c    0, no singularity detected.
c    nonzero, the factorization failed on the INFO-th step.
c
      implicit none

      integer n

      double precision a(n,n)
      double precision b(n)
      integer i
      integer info
      integer ipiv
      integer j
      integer jcol
      double precision piv
      double precision temp

      info = 0

      do jcol = 1, n
c
c  Find the maximum element in column I.
c
        piv = abs ( a(jcol,jcol) )
        ipiv = jcol
        do i = jcol + 1, n
          if ( piv .lt. abs ( a(i,jcol) ) ) then
            piv = abs ( a(i,jcol) )
            ipiv = i
          end if
        end do

        if ( piv .eq. 0.0D+00 ) then
          info = jcol
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8GE_FS - Fatal errorc'
          write ( *, '(a,i8)' ) '  Zero pivot on step ', info
          return
        end if
c
c  Switch rows JCOL and IPIV, and B.
c
        if ( jcol .ne. ipiv ) then

          do j = 1, n
            temp = a(jcol,j)
            a(jcol,j) = a(ipiv,j)
            a(ipiv,j) = temp
          end do

          temp = b(jcol)
          b(jcol) = b(ipiv)
          b(ipiv) = temp

        end if
c
c  Scale the pivot row.
c
        do j = jcol + 1, n
          a(jcol,j) = a(jcol,j) / a(jcol,jcol)
        end do
        b(jcol) = b(jcol) / a(jcol,jcol)
        a(jcol,jcol) = 1.0D+00
c
c  Use the pivot row to eliminate lower entries in that column.
c
        do i = jcol + 1, n
          if ( a(i,jcol) .ne. 0.0D+00 ) then
            temp = - a(i,jcol)
            a(i,jcol) = 0.0D+00
            do j = jcol + 1, n
              a(i,j) = a(i,j) + temp * a(jcol,j)
            end do
            b(i) = b(i) + temp * b(jcol)
          end if
        end do

      end do
c
c  Back solve.
c
      do j = n, 2, -1
        do i = 1, j - 1
          b(i) = b(i) - a(i,j) * b(j)
        end do
      end do

      return
      end
      subroutine exact ( x, y, u, dudx, dudy )

c*********************************************************************72
c
cc EXACT calculates the exact solution and its first derivatives.
c
c  Discussion:
c
c    The function specified here depends on the problem being
c    solved.  The user must be sure to change both EXACT and RHS
c    or the program will have inconsistent data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 November 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, Y, the coordinates of a point
c    in the region, at which the exact solution is to be evaluated.
c
c    Output, double precision U, DUDX, DUDY, the value of
c    the exact solution U and its derivatives dUdX
c    and dUdY at the point (X,Y).
c
      implicit none

      double precision dudx
      double precision dudy
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision u
      double precision x
      double precision y

      u    =      sin ( pi * x ) * sin ( pi * y ) + x
      dudx = pi * cos ( pi * x ) * sin ( pi * y ) + 1.0D+00
      dudy = pi * sin ( pi * x ) * cos ( pi * y )

      return
      end
