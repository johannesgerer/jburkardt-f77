      program main

c*********************************************************************72
c
cc MAIN is the main program for FEM1D.
c
c  Discussion:
c
c    The differential equation is:
c
c    - d/dx ( P(X) * dU(X)/dx ) + Q(X) * U(X) = F(X)
c
c    An approximate solution to this equation is determined
c    by the finite-element method using piecewise linear basis
c    functions.
c
c    Here U is an unknown scalar function of X defined on the
c    interval [XL,XR], and P, Q and F are given functions of X.
c
c    The values of U or U' at XL and XR are also specified.
c
c
c    The interval [XL,XR] is "meshed" with NSUB+1 points,
c
c    XN(0)=XL, XN(1)=XL+H, XN(2)=XL+2*H, ..., XN(NSUB)=XR.
c
c    This creates NSUB subintervals, with interval number 1
c    having endpoints XN(0) and XN(1), and so on up to interval
c    NSUB, which has endpoints XN(NSUB-1) and XN(NSUB).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    double precision ADIAG(NU).
c    ADIAG(I) is the "diagonal" coefficient of the I-th
c    equation in the linear system.  That is, ADIAG(I) is
c    the coefficient of the I-th unknown in the I-th equation.
c
c    double precision ALEFT(NU).
c    ALEFT(I) is the "left hand" coefficient of the I-th
c    equation in the linear system.  That is, ALEFT(I) is the
c    coefficient of the (I-1)-th unknown in the I-th equation.
c    There is no value in ALEFT(1), since the first equation
c    does not refer to a "0-th" unknown.
c
c    double precision ARITE(NU).
c    ARITE(I) is the "right hand" coefficient of the I-th
c    equation in the linear system.  ARITE(I) is the coefficient
c    of the (I+1)-th unknown in the I-th equation.  There is
c    no value in ARITE(NU) because the NU-th equation does not
c    refer to an "NU+1"-th unknown.
c
c    double precision F(NSUB+1) or F(NU).
c    ASSEMB stores into F the right hand side of the linear
c    equations.
c    SOLVE replaces those values of F by the solution of the
c    linear equations.
c
c    double precision H(NSUB)
c    H(I) is the length of subinterval I.  This code uses
c    equal spacing for all the subintervals.
c
c    integer IBC.
c    IBC declares what the boundary conditions are.
c    1, at the left endpoint, U has the value UL,
c    at the right endpoint, U' has the value UR.
c    2, at the left endpoint, U' has the value UL,
c    at the right endpoint, U has the value UR.
c    3, at the left endpoint, U has the value UL,
c    and at the right endpoint, U has the value UR.
c    4, at the left endpoint, U' has the value UL,
c    at the right endpoint U' has the value UR.
c
c    INDX   integer INDX(0:NSUB).
c    For a node I, INDX(I) is the index of the unknown
c    associated with node I.
c    If INDX(I) is equal to -1, then no unknown is associated
c    with the node, because a boundary condition fixing the
c    value of U has been applied at the node instead.
c    Unknowns are numbered beginning with 1.
c    If IBC is 2 or 4, then there is an unknown value of U
c    at node 0, which will be unknown number 1.  Otherwise,
c    unknown number 1 will be associated with node 1.
c    If IBC is 1 or 4, then there is an unknown value of U
c    at node NSUB, which will be unknown NSUB or NSUB+1,
c    depending on whether there was an unknown at node 0.
c
c    integer NL.
c    The number of basis functions used in a single
c    subinterval.  (NL-1) is the degree of the polynomials
c    used.  For this code, NL is fixed at 2, meaning that
c    piecewise linear functions are used as the basis.
c
c    integer NODE(NL,NSUB).
c    For each subinterval I:
c    NODE(1,I) is the number of the left node, and
c    NODE(2,I) is the number of the right node.
c
c    integer NQUAD.
c    The number of quadrature points used in a subinterval.
c    This code uses NQUAD=1.
c
c    integer NSUB.
c    The number of subintervals into which the interval
c    [XL,XR] is broken.
c
c    integer NU.
c    NU is the number of unknowns in the linear system.
c    Depending on the value of IBC, there will be NSUB-1,
c    NSUB, or NSUB+1 unknown values, which are the coefficients
c    of basis functions.
c
c    double precision UL.
c    If IBC is 1 or 3, UL is the value that U is required
c    to have at X=XL.
c    If IBC is 2 or 4, UL is the value that U' is required
c    to have at X=XL.
c
c    double precision UR.
c    If IBC is 2 or 3, UR is the value that U is required
c    to have at X=XR.
c    If IBC is 1 or 4, UR is the value that U' is required
c    to have at X=XR.
c
c    double precision XL.
c    XL is the left endpoint of the interval over which the
c    differential equation is being solved.
c
c    double precision XN(0:NSUB).
c    XN(I) is the location of the I-th node.  XN(0) is XL,
c    and XN(NSUB) is XR.
c
c    double precision XQUAD(NSUB)
c    XQUAD(I) is the location of the single quadrature point
c    in interval I.
c
c    double precision XR.
c    XR is the right endpoint of the interval over which the
c    differential equation is being solved.
c
      implicit none

      integer nsub
      parameter (nsub=5)
c
c  Set the number of basis functions per interval.
c  2 means linear functions are used.
c
      integer nl
      parameter (nl=2)

      double precision adiag(nsub+1)
      double precision aleft(nsub+1)
      double precision arite(nsub+1)
      double precision f(nsub+1)
      double precision h(nsub)
      integer ibc
      integer indx(0:nsub)
      integer node(nl,nsub)
      integer nquad
      integer nu
      double precision ul
      double precision ur
      double precision xl
      double precision xn(0:nsub)
      double precision xquad(nsub)
      double precision xr

      call timestamp ( )
      write ( *, * ) ' '
      write ( *, * ) 'FEM1D'
      write ( *, * ) '  FORTRAN77 version.'
      write ( *, * ) ' '
      write ( *, * ) '  A finite element solver for a 1D problem.'
      write ( *, * ) ' '
      write ( *, * ) '  Solve the two-point boundary value problem'
      write ( *, * ) ' '
      write ( *, * ) '  -d/dx ( P(x) dU(x)/dx ) + Q(x) U(x) = F(x)'
      write ( *, * ) ' '
      write ( *, * ) '  on the interval [XL,XR], specifying'
      write ( *, * ) '  the value of U or U'' at each end.'
      write ( *, * ) ' '
      write ( *, * ) 
     &  '  The interval [XL,XR] is broken into NSUB = ', nsub,
     &  ' subintervals'
      write ( *, * ) 
     &  '  Piecewise linear finite element functions are used.'
      write ( *, * ) '  The number of basis functions that are nonzero'
      write ( *, * ) '  in any element is at most NL = ', nl
c
c  Initialize the data that defines the problem.
c
      call init ( ibc, nquad, ul, ur, xl, xr )
c
c  Compute the quantities which define the geometry of the
c  problem.
c
      call geometry ( h, ibc, indx, nl, node, nsub, nu, xl, xn, 
     &  xquad, xr )
c
c  Assemble the linear system.
c
      call assemble(adiag,aleft,arite,f,h,indx,nl,node,nquad,
     &  nsub,nu,ul,ur,xn,xquad)
c
c  Print out the linear system.
c
      call prsys(adiag,aleft,arite,f,nu)
c
c  Solve the linear system.
c
      call solve(adiag,aleft,arite,f,nu)
c
c  Print out the solution.
c
      call output(f,ibc,indx,nsub,nu,ul,ur,xn)
c
c  Terminate.
c
      write ( *, * ) ' '
      write ( *, * ) 'FEM1D'
      write ( *, * ) '  Normal end of execution.'

      write ( *, * ) ' '
      call timestamp ( )

      stop
      end
      subroutine assemble ( adiag, aleft, arite, f, h, indx, nl, node, 
     &  nquad, nsub, nu, ul, ur, xn, xquad )

c*********************************************************************72
c
cc ASSEMBLE assembles the matrix and the right-hand-side.
c
c  Discussion:
c
c    The differential equation is:
c
c    - d/dx ( P(X) * dU(X)/dx ) + Q(X) * U(X) = F(X)
c
c    The linear system is
c
c      K * C = F
c
c    that is to be solved for the coefficients C.
c
c    Numerical integration is used to compute the entries of K
c    and F.
c
c    Note that a 1 point quadrature rule, which is sometimes used to
c    assemble the matrix and right hand side, is just barely accurate
c    enough for simple problems.  If you want better results, you
c    should use a quadrature rule that is more accurate.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    29 April 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    double precision ADIAG(NU).
c    ADIAG(I) is the "diagonal" coefficient of the I-th
c    equation in the linear system.  That is, ADIAG(I) is
c    the coefficient of the I-th unknown in the I-th equation.
c
c    double precision ALEFT(NU).
c    ALEFT(I) is the "left hand" coefficient of the I-th
c    equation in the linear system.  That is, ALEFT(I) is the
c    coefficient of the (I-1)-th unknown in the I-th equation.
c    There is no value in ALEFT(1), since the first equation
c    does not refer to a "0-th" unknown.
c
c    double precision ARITE(NU).
c    ARITE(I) is the "right hand" coefficient of the I-th
c    equation in the linear system.  ARITE(I) is the coefficient
c    of the (I+1)-th unknown in the I-th equation.  There is
c    no value in ARITE(NU) because the NU-th equation does not
c    refer to an "NU+1"-th unknown.
c
c    double precision F(NSUB+1) or F(NU).
c    ASSEMB stores into F the right hand side of the linear
c    equations.
c    SOLVE replaces those values of F by the solution of the
c    linear equations.
c
c    double precision H(NSUB)
c    H(I) is the length of subinterval I.  This code uses
c    equal spacing for all the subintervals.
c
c    INDX   integer INDX(0:NSUB).
c    For a node I, INDX(I) is the index of the unknown
c    associated with node I.
c    If INDX(I) is equal to -1, then no unknown is associated
c    with the node, because a boundary condition fixing the
c    value of U has been applied at the node instead.
c    Unknowns are numbered beginning with 1.
c    If IBC is 2 or 4, then there is an unknown value of U
c    at node 0, which will be unknown number 1.  Otherwise,
c    unknown number 1 will be associated with node 1.
c    If IBC is 1 or 4, then there is an unknown value of U
c    at node NSUB, which will be unknown NSUB or NSUB+1,
c    depending on whether there was an unknown at node 0.
c
c    integer NL.
c    The number of basis functions used in a single
c    subinterval.  (NL-1) is the degree of the polynomials
c    used.  For this code, NL is fixed at 2, meaning that
c    piecewise linear functions are used as the basis.
c
c    integer NODE(NL,NSUB).
c    For each subinterval I:
c    NODE(1,I) is the number of the left node, and
c    NODE(2,I) is the number of the right node.
c
c    integer NQUAD.
c    The number of quadrature points used in a subinterval.
c    This code uses NQUAD=1.
c
c    integer NSUB.
c    The number of subintervals into which the interval
c    [XL,XR] is broken.
c
c    integer NU.
c    NU is the number of unknowns in the linear system.
c    Depending on the value of IBC, there will be NSUB-1,
c    NSUB, or NSUB+1 unknown values, which are the coefficients
c    of basis functions.
c
c    double precision UL.
c    If IBC is 1 or 3, UL is the value that U is required
c    to have at X=XL.
c    If IBC is 2 or 4, UL is the value that U' is required
c    to have at X=XL.
c
c    double precision UR.
c    If IBC is 2 or 3, UR is the value that U is required
c    to have at X=XR.
c    If IBC is 1 or 4, UR is the value that U' is required
c    to have at X=XR.
c
c    double precision XN(0:NSUB).
c    XN(I) is the location of the I-th node.  XN(0) is XL,
c    and XN(NSUB) is XR.
c
c    double precision XQUAD(NSUB)
c    XQUAD(I) is the location of the single quadrature point
c    in interval I.
c
      implicit none

      integer nl
      integer nsub
      integer nu

      double precision aij
      double precision adiag(nu)
      double precision aleft(nu)
      double precision arite(nu)
      double precision f(nu)
      double precision ff
      double precision h(nsub)
      double precision he
      integer i
      integer ie
      integer ig
      integer il
      integer indx(0:nsub)
      integer iq
      integer iu
      integer jg
      integer jl
      integer ju
      integer node(nl,nsub)
      integer nquad
      double precision phii
      double precision phiix
      double precision phij
      double precision phijx
      double precision pp
      double precision qq
      double precision ul
      double precision ur
      double precision x
      double precision xleft
      double precision xn(0:nsub)
      double precision xquad(nsub)
      double precision xquade
      double precision xrite
c
c  Zero out the arrays that hold the coefficients of the matrix
c  and the right hand side.
c
      do i = 1, nu
        f(i)=0.0D+00
        adiag(i)=0.0D+00
        aleft(i)=0.0D+00
        arite(i)=0.0D+00
      end do
c
c  For interval number IE,
c
      do ie = 1, nsub

        he=h(ie)
        xleft=xn(node(1,ie))
        xrite=xn(node(2,ie))
c
c  consider each quadrature point IQ,
c
        do iq = 1, nquad

          xquade=xquad(ie)
c
c  and evaluate the integrals associated with the basis functions
c  for the left, and for the right nodes.
c
          do 30 il = 1, nl

            ig=node(il,ie)
            iu=indx(ig)

            if(iu.le.0)go to 30

            call phi(il,xquade,phii,phiix,xleft,xrite)

            f(iu) = f(iu) + he * ff ( xquade ) * phii
c
c  Take care of boundary nodes at which U' was specified.
c
            if ( ig .eq. 0 ) then

              x = 0.0D+00
              f(iu) = f(iu) - pp ( x ) * ul

            else if (ig.eq.nsub ) then

              x = 1.0D+00
              f(iu) = f(iu) + pp ( x ) * ur

            end if
c
c  Evaluate the integrals that take a product of the basis
c  function times itself, or times the other basis function
c  that is nonzero in this interval.
c
            do jl = 1, nl

              jg = node(jl,ie)
              ju = indx(jg)

              call phi ( jl, xquade, phij, phijx, xleft, xrite )

              aij = he * ( pp ( xquade ) * phiix * phijx
     &                   + qq ( xquade ) * phii  * phij   )
c
c  If there is no variable associated with the node, then it's
c  a specified boundary value, so we multiply the coefficient
c  times the specified boundary value and subtract it from the
c  right hand side.
c
              if ( ju .le. 0 ) then

                if ( jg .eq. 0 ) then

                  f(iu) = f(iu) - aij * ul

                else if (jg .eq. nsub ) then

                  f(iu) = f(iu) - aij * ur

                end if
c
c  Otherwise, we add the coefficient we've just computed to the
c  diagonal, or left or right entries of row IU of the matrix.
c
              else

                if ( iu .eq. ju ) then

                  adiag(iu) = adiag(iu) + aij

                else if ( ju .lt. iu ) then

                  aleft(iu) = aleft(iu) + aij

                else

                  arite(iu) = arite(iu) + aij

                end if

              end if

            end do

30        continue

        end do

      end do

      return
      end
      function ff ( x )

c*********************************************************************72
c
cc FF evaluates the function F in the differential equation.
c
c    The differential equation is:
c
c    - d/dx ( P(X) * dU(X)/dx ) + Q(X) * U(X) = F(X)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the point where F is to be evaluated.
c
c    Output, double precision FF, the value of F(X).
c
      implicit none

      double precision ff
      double precision x

      ff = 0.0D+00

      return
      end
      subroutine geometry ( h, ibc, indx, nl, node, nsub, nu, xl, 
     &  xn, xquad, xr )

c*********************************************************************72
c
cc GEOMETRY sets up the geometry for the interval [XL,XR].
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision H(NSUB)
c    H(I) is the length of subinterval I.  This code uses
c    equal spacing for all the subintervals.
c
c    Input, integer IBC.
c    IBC declares what the boundary conditions are.
c    1, at the left endpoint, U has the value UL,
c    at the right endpoint, U' has the value UR.
c    2, at the left endpoint, U' has the value UL,
c    at the right endpoint, U has the value UR.
c    3, at the left endpoint, U has the value UL,
c    and at the right endpoint, U has the value UR.
c    4, at the left endpoint, U' has the value UL,
c    at the right endpoint U' has the value UR.
c
c    Output, INDX   integer INDX(0:NSUB).
c    For a node I, INDX(I) is the index of the unknown
c    associated with node I.
c    If INDX(I) is equal to -1, then no unknown is associated
c    with the node, because a boundary condition fixing the
c    value of U has been applied at the node instead.
c    Unknowns are numbered beginning with 1.
c    If IBC is 2 or 4, then there is an unknown value of U
c    at node 0, which will be unknown number 1.  Otherwise,
c    unknown number 1 will be associated with node 1.
c    If IBC is 1 or 4, then there is an unknown value of U
c    at node NSUB, which will be unknown NSUB or NSUB+1,
c    depending on whether there was an unknown at node 0.
c
c    Input, integer NL.
c    The number of basis functions used in a single
c    subinterval.  (NL-1) is the degree of the polynomials
c    used.  For this code, NL is fixed at 2, meaning that
c    piecewise linear functions are used as the basis.
c
c    Output, integer NODE(NL,NSUB).
c    For each subinterval I:
c    NODE(1,I) is the number of the left node, and
c    NODE(2,I) is the number of the right node.
c
c    Input, integer NSUB.
c    The number of subintervals into which the interval
c    [XL,XR] is broken.
c
c    Output, integer NU.
c    NU is the number of unknowns in the linear system.
c    Depending on the value of IBC, there will be NSUB-1,
c    NSUB, or NSUB+1 unknown values, which are the coefficients
c    of basis functions.
c
c    Input, double precision XL.
c    XL is the left endpoint of the interval over which the
c    differential equation is being solved.
c
c    Output, double precision XN(0:NSUB).
c    XN(I) is the location of the I-th node.  XN(0) is XL,
c    and XN(NSUB) is XR.
c
c    Output, double precision XQUAD(NSUB)
c    XQUAD(I) is the location of the single quadrature point
c    in interval I.
c
c    Input, double precision XR.
c    XR is the right endpoint of the interval over which the
c    differential equation is being solved.
c
      implicit none

      integer nl
      integer nsub

      double precision h(nsub)
      integer i
      integer ibc
      integer indx(0:nsub)
      integer node(nl,nsub)
      integer nu
      double precision xl
      double precision xn(0:nsub)
      double precision xquad(nsub)
      double precision xr
c
c  Set the value of XN, the locations of the nodes.
c
      write(*,*)' '
      write(*,*)'  Node      Location'
      write(*,*)' '
      do i = 0, nsub
        xn(i) = ( dble ( nsub - i ) * xl   
     &          + dble (        i ) * xr ) 
     &          / dble ( nsub )
        write(*,*) i, xn(i)
      end do
c
c  Set the lengths of each subinterval.
c
      write(*,*)' '
      write(*,*)'Subint    Length'
      write(*,*)' '
      do i = 1, nsub
        h(i) = xn(i) - xn(i-1)
        write(*,*) i, h(i)
      end do
c
c  Set the quadrature points, each of which is the midpoint
c  of its subinterval.
c
      write(*,*)' '
      write(*,*)'Subint    Quadrature point'
      write(*,*)' '
      do i = 1, nsub
        xquad(i) = 0.5D+00 * ( xn(i-1) + xn(i) )
        write(*,*)i, xquad(i)
      end do
c
c  Set the value of NODE, which records, for each interval,
c  the node numbers at the left and right.
c
      write(*,*)' '
      write(*,*)'Subint  Left Node  Right Node'
      write(*,*)' '
      do i = 1, nsub
        node(1,i) = i-1
        node(2,i) = i
        write(*,*) i, node(1,i), node(2,i)
      end do
c
c  Starting with node 0, see if an unknown is associated with
c  the node.  If so, give it an index.
c
      nu = 0
c
c  Handle first node.
c
      i = 0
      if ( ibc .eq. 1 .or. ibc .eq. 3 ) then
        indx(i) = -1
      else
        nu = nu + 1
        indx(i) = nu
      end if
c
c  Handle nodes 1 through nsub-1
c
      do i = 1, nsub-1
        nu = nu + 1
        indx(i) = nu
      end do
c
c  Handle the last node.
c
      i = nsub

      if ( ibc .eq. 2 .or. ibc .eq. 3 ) then
        indx(i) = -1
       else
        nu = nu + 1
        indx(i) = nu
       end if

      write(*,*)' '
      write(*,*)'  Number of unknowns NU = ', nu
      write(*,*)' '
      write(*,*)'  Node  Unknown'
      write(*,*)' '
      do i = 0, nsub
        write(*,*) i, indx(i)
      end do

      return
      end
      subroutine init ( ibc, nquad, ul, ur, xl, xr )

c*********************************************************************72
c
cc INIT assigns values to variables which define the problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
      implicit none

      integer ibc
      integer nquad
      double precision ul
      double precision ur
      double precision xl
      double precision xr
c
c  IBC declares what the boundary conditions are.
c
      ibc = 1
c
c  NQUAD is the number of quadrature points per subinterval.
c  The program as currently written cannot handle any value for
c  NQUAD except 1!
c
      nquad = 1
c
c  Set the values of U or U' at the endpoints.
c
      ul = 0.0D+00
      ur = 1.0D+00
c
c  Define the location of the endpoints of the interval.
c
      xl = 0.0D+00
      xr = 1.0D+00
c
c  Print out the values that have been set.
c
      write(*,*)' '
      write(*,*)'The equation is to be solved for'
      write(*,*)'X greater than XL = ',xl,' and less than XR = ',xr
      write(*,*)' '
      write(*,*)'The boundary conditions are:'
      write(*,*)' '

      if ( ibc .eq. 1 .or. ibc .eq. 3 ) then
        write(*,*)'  At X=XL, U=',ul
      else
        write(*,*)'  At X=XL, U''=',ul
      end if

      if ( ibc .eq. 2 .or. ibc .eq. 3 ) then
        write(*,*)'  At X=XR, U=',ur
      else
        write(*,*)'  At X=XR, U''=',ur
      end if

      write(*,*)' '
      write(*,*)'Number of quadrature points per element is ', nquad

      return
      end
      subroutine output ( f, ibc, indx, nsub, nu, ul, ur, xn )

c*********************************************************************72
c
cc OUTPUT prints out the computed solution.
c
c  Discussion:
c
c    We simply print out the solution vector F, except that, for
c    certain boundary conditions, we are going to have to get the
c    value of the solution at XL or XR by using the specified
c    boundary value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
      implicit none

      integer nsub
      integer nu

      double precision f(nu)
      integer i
      integer ibc
      integer indx(0:nsub)
      double precision u
      double precision ul
      double precision ur
      double precision xn(0:nsub)

      write(*,*)' '
      write(*,*)'Computed solution:'
      write(*,*)' '
      write(*,*)'Node    X(I)        U(X(I))'
      write(*,*)' '

      do i = 0, nsub
c
c  If we're at the first node, check the boundary condition.
c
        if ( i .eq. 0 ) then

          if ( ibc .eq. 1 .or. ibc .eq. 3 ) then
            u = ul
          else
            u = f(indx(i))
          end if
c
c  If we're at the last node, check the boundary condition.
c
        else if ( i .eq. nsub ) then

          if ( ibc .eq. 2 .or. ibc .eq. 3 ) then
            u = ur
          else
            u = f(indx(i))
          end if
c
c  Any other node, we're sure the value is stored in F.
c
        else

          u=f(indx(i))

        end if

        write(*,*) i, xn(i), u

      end do

      return
      end
      subroutine phi ( il, x, phii, phiix, xleft, xrite )

c*********************************************************************72
c
cc PHI evaluates a linear basis function and its derivative.
c
c  Discussion:
c
c    The evaluation is done at a point X in an interval [XLEFT,XRITE]. 
c 
c    In this interval, there are just two nonzero basis functions.  
c    The first basis function is a line which is 1 at the left 
c    endpoint and 0 at the right.  The second basis function is 0 at
c    the left endpoint and 1 at the right.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
      implicit none

      integer il
      double precision phii
      double precision phiix
      double precision x
      double precision xleft
      double precision xrite

      if ( xleft .le. x .and. x .le. xrite ) then

        if ( il .eq. 1 ) then
          phii = ( xrite - x ) / ( xrite - xleft )
          phiix = -1.0D+00 / ( xrite - xleft )
        else
          phii = ( x - xleft )/ ( xrite - xleft )
          phiix = 1.0D+00 / ( xrite - xleft )
        end if
c
c  If X is outside of the interval, just set everything to 0.
c
      else

        phii = 0.0D+00
        phiix = 0.0D+00

      end if

      return
      end
      function pp ( x )

c*********************************************************************72
c
cc PP evaluates the function P(X) in the differential equation:
c
c    The differential equation is:
c
c    - d/dx ( P(X) * dU(X)/dx ) + Q(X) * U(X) = F(X)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the point where P is to be evaluated.
c
c    Output, double precision PP, the value of P(X).
c
      implicit none

      double precision pp
      double precision x

      pp = 1.0D+00

      return
      end
      subroutine prsys ( adiag, aleft, arite, f, nu )

c*********************************************************************72
c
cc PRSYS prints out the tridiagonal linear system.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
      implicit none

      integer nu

      double precision adiag(nu)
      double precision aleft(nu)
      double precision arite(nu-1)
      double precision f(nu)
      integer i

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'The tridiagonal linear system:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Equation  ALEFT  ADIAG  ARITE  RHS'
      write ( *, '(a)' ) ' '

      do i = 1, nu

        write ( *, * ) i, aleft(i), adiag(i), arite(i), f(i)

      end do

      return
      end
      function qq ( x )

c*********************************************************************72
c
cc QQ evaluates the function Q(X) in the differential equation.
c
c  Discussion:
c
c    The differential equation is:
c
c    - d/dx ( P(X) * dU(X)/dx ) + Q(X) * U(X) = F(X)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the point where F is to be evaluated.
c
c    Output, double precision QQ, the value of Q(X).
c
      implicit none

      double precision qq
      double precision x

      qq = 0.0D+00

      return
      end
      subroutine solve ( adiag, aleft, arite, f, nu )

c*********************************************************************72
c
cc SOLVE solves a tridiagonal matrix system of the form A*x=b.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    01 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, double precision ADIAG(NU), ALEFT(NU), ARITE(NU).
c    On input, ADIAG, ALEFT, and ARITE contain the diagonal,
c    left and right entries of the equations.
c    On output, ADIAG and ARITE have been changed in order
c    to compute the solution.
c    Note that for the first equation, there is no ALEFT
c    coefficient, and for the last, there is no ARITE.
c    So there is no need to store a value in ALEFT(1), nor
c    in ARITE(NU).
c
c    Input/output, double precision F(NU).
c    On input, F contains the right hand side of the linear
c    system to be solved.
c    On output, F contains the solution of the linear system.
c
c    Input, INTEGER NU.
c    NU is the number of equations to be solved.
c
      implicit none

      integer nu

      double precision adiag(nu)
      double precision aleft(nu)
      double precision arite(nu-1)
      double precision f(nu)
      integer i
c
c  Carry out Gauss elimination on the matrix, saving information
c  needed for the backsolve.
c
      arite(1) = arite(1) / adiag(1)
      do i = 2, nu-1
        adiag(i) = adiag(i) - aleft(i) * arite(i-1)
        arite(i) = arite(i) / adiag(i)
      end do
      adiag(nu) = adiag(nu) - aleft(nu) * arite(nu-1)
c
c  Carry out the same elimination steps on F that were done to the
c  matrix.
c
      f(1)=f(1) / adiag(1)
      do i = 2, nu
        f(i) = ( f(i) - aleft(i) * f(i-1) ) / adiag(i)
      end do
c
c  And now carry out the steps of "back substitution".
c
      do i = nu - 1, 1, -1
        f(i) = f(i) - arite(i) * f(i+1)
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
