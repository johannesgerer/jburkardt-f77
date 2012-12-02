      program main

c*********************************************************************72
c
cc MAIN is the main program for FEM1D_ADAPTIVE.
c
c  Discussion:
c
c    FEM1D_ADAPTIVE solves the one dimensional problem:
c
c    - d/dX (P dU/dX) + Q U = F
c
c    by the finite-element method using piecewise linear basis
c    functions.
c
c    An adaptive method is used to try to reduce the maximum
c    error by refining the mesh in certain places.
c
c    Here U is an unknown scalar function of X defined on the
c    interval [XL,XR], and P, Q and F are given functions of X.
c
c    The values of U at XL and XR are also specified.
c
c
c    The interval [XL,XR] is "meshed" with N+1 points,
c
c    XN(0) = XL, XN(1) = XL+H, XN(2) = XL+2*H, ..., XN(N) = XR.
c
c    This creates N subintervals, with interval number 1
c    having endpoints XN(0) and XN(1), and so on up to interval
c    N, which has endpoints XN(N-1) and XN(N).
c
c
c    The algorithm employed by ADAPT tries to guarantee a certain amount
c    of accuracy by examining the current solution, estimating the error
c    in each subinterval, and, if necessary, subdividing one or more
c    subintervals and repeating the calculation.
c
c    We can think of the adaptive part of the algorithm as a refined
c    problem.  The program re-solves the problem on the pair of
c    intervals J and J+1, which extend from node J-1 to node J+1.
c    The values of U that were just computed at nodes J-1 and J+1
c    will be used as the boundary values for this refined problem.
c    The intervals J and J+1 will each be evenly divided into NY
c    smaller subintervals.  This boundary value problem is solved,
c    and the derivatives of the original and refined solutions are
c    then compared to get an estimate of the error.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 November 2006
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
c    double precision ETA(N).
c    ETA(I) is the error estimate for interval I.  It is computed
c    as the sum of two quantities, one associated with the left
c    and one with the right node of the interval.
c
c    double precision F(NU).
c    ASSEMBLE stores into F the right hand side of the linear
c    equations.
c    SOLVE replaces those values of F by the solution of the
c    linear equations.
c
c    double precision FY(M).
c    FY is the right hand side of the linear system of the refined
c    problem.
c
c    double precision H(N)
c    H(I) is the length of subinterval I.  This code uses
c    equal spacing for all the subintervals.
c
c    double precision HY(M).
c    HY(I) is the length of subinterval I in the refined problem.
c
c    integer IBC.
c    IBC declares what the boundary conditions are.
c    1, at the left endpoint, U has the value UL,
c       at the right endpoint, U' has the value UR.
c    2, at the left endpoint, U' has the value UL,
c       at the right endpoint, U has the value UR.
c    3, at the left endpoint, U has the value UL,
c       and at the right endpoint, U has the value UR.
c    4, at the left endpoint, U' has the value UL,
c       at the right endpoint U' has the value UR.
c
c    integer IBCY.
c    IBCY declares the boundary conditions for the refined problem
c    which should always be that the value of U is specified at
c    both the left and right endpoints.  This corresponds to a
c    value of IBCY = 3.
c
c    integer INDX(0:N).
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
c    at node N, which will be unknown N or N+1,
c    depending on whether there was an unknown at node 0.
c
c    integer INDY(0:M).
c    INDY(I) records the index of the unknown associated with
c    node I for the refined problem.
c
c    integer JADD(N).
c    JADD(I) is 1 if the error estimates show that interval I
c    should be subdivided.
c
c    integer KOUNT, the number of adaptive steps that have been taken.
c
c    integer M.
c    M is the number of subintervals used in the refined problem.
c    M is equal to NY for computations centered at node 0 or node N,
c    and otherwise, M is equal to 2*NY.
c
c    integer N.
c    The number of subintervals into which the interval
c    [XL,XR] is broken.
c
c    integer NL.
c    The number of basis functions used in a single
c    subinterval.  (NL-1) is the degree of the polynomials
c    used.  For this code, NL is fixed at 2, meaning that
c    piecewise linear functions are used as the basis.
c
c    integer NMAX, the maximum number of unknowns that can be handled.
c
c    integer NODE(NL,N).
c    For each subinterval I:
c    NODE(1,I) is the number of the left node, and
c    NODE(2,I) is the number of the right node.
c
c    integer NODEY(NL,M).
c    NODEY performs the same function for the refined problem that
c    NODE performs for the full problem, recording the node numbers
c    associated with a particular subinterval.
c
c    integer NQUAD
c    The number of quadrature points used in a subinterval.
c
c    integer NU.
c    NU is the number of unknowns in the linear system.
c    Depending on the value of IBC, there will be N-1,
c    N, or N+1 unknown values, which are the coefficients
c    of basis functions.
c
c    integer NUY.
c    The number of unknowns in the refined problem.
c
c    integer NY.
c    NY is the number of subintervals into which a given interval
c    will be subdivided, before solving the refined probelm.
c
c    integer PROBLEM, chooses the problem to be solved.
c    The user must choose this value by setting it in routine get_problem.
c    * 1, u = x, p = 1, q = 0, f = 0, ibc = 3, ul = 0, ur = 1.
c    The program should find the solution exactly, and the
c    adaptive code should find that there is no reason to
c    subdivide any interval.
c    * 2, u = x*x, p = 1, q = 0, f = -2, ibc = 3, ul = 0, ur = 1.
c    This problem should find the solution exactly, and
c    the adaptive code should again find there is nothing
c    to do.
c    *3, u = sin(pi*x/2), p = 1, q = 0, ibc = 3, f = 0.25*pi*pi*sin(pi*x/2), 
c    ul = 0, ur = 1.
c    *4, u = cos(pi*x/2), p = 1, q = 0, ibc = 3, f = 0.25*pi*pi*cos(pi*x/2), 
c    ul = 1, ur = 0.
c    *5: u = x**(beta+2)/((beta+2)*(beta+1)), p = 1, q = 1, ibc = 3, 
c    f = -x**beta + (x**(beta+2))/((beta+2)*(beta+1)),
c    ul = 0, ur = 1/((beta+2)*(beta+1))
c    (beta must be greater than -2, and not equal to -1)
c    *6: u = atan((x-0.5)/alpha), p = 1, q = 0, ibc = 3, 
c    f =  2*alpha*(x-0.5) / (alpha**2 + (x-0.5)**2) **2,
c    ul = u(0), ur = u(1)
c
c    integer STATUS, reports status of subdivision.
c    0, a new subdivision was carried out.
c    1, no more subdivisions are needed.
c    -1, no more subdivisions can be carried out.
c
c    double precision TOL.
c    A tolerance that is used to determine whether the estimated
c    error in an interval is so large that it should be subdivided
c    and the problem solved again.
c
c    double precision UL.
c    If IBC is 1 or 3, UL is the value that U is required
c    to have at X = XL.
c    If IBC is 2 or 4, UL is the value that U' is required
c    to have at X = XL.
c
c    double precision UR.
c    If IBC is 2 or 3, UR is the value that U is required
c    to have at X = XR.
c    If IBC is 1 or 4, UR is the value that U' is required
c    to have at X = XR.
c
c    double precision WQUAD(NQUAD).
c    WQUAD(I) is the weight associated with the I-th point
c    of an NQUAD point Gaussian quadrature rule.
c
c    double precision XL.
c    XL is the left endpoint of the interval over which the
c    differential equation is being solved.
c
c    double precision XN(0:N).
c    XN(I) is the location of the I-th node.  XN(0) is XL,
c    and XN(N) is XR.
c
c    double precision XQUAD(NQUAD,NMAX), the I-th quadrature point
c    in interval J.
c
c    double precision XQUADY(NQUAD,NMAY), the I-th quadrature point
c    in subinterval J of the refined problem.
c
c    double precision XR.
c    XR is the right endpoint of the interval over which the
c    differential equation is being solved.
c
c    Workspace, double precision XT(0:NMAX), used to compute a new
c    set of nodes.
c
c    double precision YN(0:M).
c    YN(I) is the location of the I-th node in the refined
c    problem.
c
      implicit none

      integer nl
      parameter (nl = 2)

      integer nmax
      parameter (nmax = 30)

      integer nquad
      parameter (nquad = 2)

      double precision adiag(nmax)
      double precision aleft(nmax)
      double precision alpha
      double precision arite(nmax)
      double precision beta
      double precision eta(nmax)
      double precision f(nmax)
      double precision h(nmax)
      integer ibc
      integer indx(0:nmax)
      integer problem
      integer jadd(nmax)
      integer kount
      integer n
      integer node(nl,nmax)
      integer nu
      integer status
      double precision tol
      double precision ul
      double precision ur
      double precision wquad(nquad)
      double precision xl
      double precision xn(0:nmax)
      double precision xquad(nquad,nmax)
      double precision xr
      double precision xt(0:nmax)

      call timestamp ( )
      write(*,*)' '
      write(*,*)'FEM1D_ADAPTIVE'
      write(*,*)'  FORTRAN77 version'
      write(*,*)' '
      write(*,*)'Solve the two-point boundary value problem:'
      write(*,*)' '
      write(*,*)'  -d/dx P du/dx + Q U = F'
      write(*,*)' '
      write(*,*)'on the interval [0,1], specifying the value'
      write(*,*)'of U or U'' at each endpoint.'
      write(*,*)' '
      write(*,*)'The number of basis functions per element is ', nl
      write(*,*)' '
      write(*,*)' '
      write(*,*)'The number of quadrature points per element is ', nquad

      call get_problem ( problem )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Problem index = ', problem
      write ( *, '(a)' ) ' '

      if ( problem == 1 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  "Linear" problem:'
        write ( *, '(a)' ) '  (No refinement needed)'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  U(X) =  X'
        write ( *, '(a)' ) '  P(X) =  1.0'
        write ( *, '(a)' ) '  Q(X) =  0.0'
        write ( *, '(a)' ) '  F(X) =  0.0'
        write ( *, '(a)' ) '  IBC  =  3'
        write ( *, '(a)' ) '  UL   =  0.0'
        write ( *, '(a)' ) '  UR   =  1.0'

      else if ( problem == 2 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  "Quadratic" problem:'
        write ( *, '(a)' ) '  (No refinement needed)'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  U(X) =  X*X'
        write ( *, '(a)' ) '  P(X) =  1.0'
        write ( *, '(a)' ) '  Q(X) =  0.0'
        write ( *, '(a)' ) '  F(X) = -2.0'
        write ( *, '(a)' ) '  IBC  =  3'
        write ( *, '(a)' ) '  UL   =  0.0'
        write ( *, '(a)' ) '  UR   =  1.0'

      else if ( problem == 3 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  "SINE" problem:'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  U(X) =  SIN(PI*X/2)'
        write ( *, '(a)' ) '  P(X) =  1.0'
        write ( *, '(a)' ) '  Q(X) =  0.0'
        write ( *, '(a)' ) '  F(X) =  PI*PI*SIN(PI*X/2)/4'
        write ( *, '(a)' ) '  IBC  =  3'
        write ( *, '(a)' ) '  UL   =  0.0'
        write ( *, '(a)' ) '  UR   =  1.0'

      else if ( problem == 4 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  "COSINE" problem:'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  U(X) =  COS(PI*X/2)'
        write ( *, '(a)' ) '  P(X) =  1.0'
        write ( *, '(a)' ) '  Q(X) =  0.0'
        write ( *, '(a)' ) '  F(X) =  PI*PI*COS(PI*X/2)/4'
        write ( *, '(a)' ) '  IBC  =  3'
        write ( *, '(a)' ) '  UL   =  0.0'
        write ( *, '(a)' ) '  UR   =  1.0'

      else if ( problem == 5 ) then

        call get_beta ( beta )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  "RHEINBOLDT" problem:'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  U(X) =  X**(B+2)/((B+2)*(B+1))'
        write ( *, '(a)' ) '  P(X) =  1.0'
        write ( *, '(a)' ) '  Q(X) =  1.0'
        write ( *, '(a)' ) '  F(X) =  -X**B+(X**B+2))/((B+2)*(B+1))'
        write ( *, '(a)' ) '  IBC  =  3'
        write ( *, '(a)' ) '  UL   =  0.0'
        write ( *, '(a)' ) '  UR   =  1/((B+2)*(B+1))'
        write ( *, '(a,g14.6)' ) '  B    = ', beta

      else if ( problem == 6 ) then

        call get_alpha ( alpha )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  "ARCTAN" problem:'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  U(X) =  ATAN((X-0.5)/A)'
        write ( *, '(a)' ) '  P(X) =  1.0'
        write ( *, '(a)' ) '  Q(X) =  0.0'
        write ( *, '(a)' ) 
     &  '  F(X) =  2*A*(X-0.5)/(A**2+(X-0.5)**2)**2'
        write ( *, '(a)' ) '  IBC  =  3'
        write ( *, '(a)' ) '  UL   =  ATAN(-0.5/A)'
        write ( *, '(a)' ) '  UR   =  ATAN( 0.5/A)'
        write ( *, '(a,g14.6)' ) '  A    = ', alpha

      end if
c
c  Start out with just 4 subintervals.
c
      n = 4
c
c  Initialize values that define the problem.
c
      call init ( ibc, n, tol, ul, ur, xl, xn, xr )
c
c  Start the iteration counter off at 0.
c
      kount = 0
c
c  Begin the next iteration.
c
10    continue

      kount = kount + 1

      write(*,*)' '
      write(*,*)'Begin new iteration with ', n, ' nodes.'
      write(*,*)' '
c
c  Solve the regular problem.
c
      call solvex ( adiag, aleft, arite, f, h, ibc, indx, kount, 
     &  n, nl, nmax, node, nquad, nu, ul, ur, wquad, xn, xquad )
c
c  Solve N subproblems to get the error estimators.
c
      call solvey ( eta, f, h, n, nu, ul, ur, xn )
c
c  Examine the error estimators, and see how many intervals should
c  be subdivided.
c
      call subdiv ( eta, jadd, kount, n, nmax, tol, xn, xt, status )
c
c  Solve the problem again, with the new nodes.
c
      if ( status = =  0 ) then
        go to 10
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM1D_ADAPTIVE:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine assemble ( adiag, aleft, arite, f, h, n, indx, node,
     &  nu, nl, nquad, nmax, ul, ur, wquad, xn, xquad )

c*********************************************************************72
c
cc ASSEMBLE assembles the global matrix.
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
c    Output, double precision ADIAG(NU).
c    ADIAG(I) is the "diagonal" coefficient of the I-th
c    equation in the linear system.  That is, ADIAG(I) is
c    the coefficient of the I-th unknown in the I-th equation.
c
c    Output, double precision ALEFT(NU).
c    ALEFT(I) is the "left hand" coefficient of the I-th
c    equation in the linear system.  That is, ALEFT(I) is the
c    coefficient of the (I-1)-th unknown in the I-th equation.
c    There is no value in ALEFT(1), since the first equation
c    does not refer to a "0-th" unknown.
c
c    Output, double precision ARITE(NU).
c    ARITE(I) is the "right hand" coefficient of the I-th
c    equation in the linear system.  ARITE(I) is the coefficient
c    of the (I+1)-th unknown in the I-th equation.  There is
c    no value in ARITE(NU) because the NU-th equation does not
c    refer to an "NU+1"-th unknown.
c
c    Output, double precision F(NU).
c    ASSEMBLE stores into F the right hand side of the linear
c    equations.
c    SOLVE replaces those values of F by the solution of the
c    linear equations.
c
c    Input, double precision H(N)
c    H(I) is the length of subinterval I.  This code uses
c    equal spacing for all the subintervals.
c
c    integer N.
c    The number of subintervals into which the interval
c    [XL,XR] is broken.
c
c    Input, integer INDX(0:N).
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
c    at node N, which will be unknown N or N+1,
c    depending on whether there was an unknown at node 0.
c
c    Input, integer NODE(NL,N).
c    For each subinterval I:
c    NODE(1,I) is the number of the left node, and
c    NODE(2,I) is the number of the right node.
c
c    Input, integer NU.
c    NU is the number of unknowns in the linear system.
c    Depending on the value of IBC, there will be N-1,
c    N, or N+1 unknown values, which are the coefficients
c    of basis functions.
c
c    Input, integer NL.
c    The number of basis functions used in a single
c    subinterval.  (NL-1) is the degree of the polynomials
c    used.  For this code, NL is fixed at 2, meaning that
c    piecewise linear functions are used as the basis.
c
c    Input, integer NQUAD
c    The number of quadrature points used in a subinterval.
c
c    Input, integer NMAX, the maximum number of unknowns that can be handled.
c
c    Input, double precision UL.
c    If IBC is 1 or 3, UL is the value that U is required
c    to have at X = XL.
c    If IBC is 2 or 4, UL is the value that U' is required
c    to have at X = XL.
c
c    Input, double precision UR.
c    If IBC is 2 or 3, UR is the value that U is required
c    to have at X = XR.
c    If IBC is 1 or 4, UR is the value that U' is required
c    to have at X = XR.
c
c    Input, double precision WQUAD(NQUAD).
c    WQUAD(I) is the weight associated with the I-th point
c    of an NQUAD point Gaussian quadrature rule.
c
c    Input, double precision XN(0:N).
c    XN(I) is the location of the I-th node.  XN(0) is XL,
c    and XN(N) is XR.
c
c    double precision XQUAD(NQUAD,NMAX), the I-th quadrature point
c    in interval J.
c
      implicit none

      integer n
      integer nl
      integer nmax
      integer nquad
      integer nu

      double precision adiag(nu)
      double precision aleft(nu)
      double precision arite(nu)
      double precision aij
      double precision f(nu)
      double precision ff
      external ff
      double precision h(n)
      double precision he
      integer i
      integer ie
      integer ig
      integer il
      integer indx(0:n)
      integer iq
      integer iu
      integer jg
      integer jl
      integer ju
      integer node(nl,nmax)
      double precision phii
      double precision phiix
      double precision phij
      double precision phijx
      double precision pp
      external pp
      double precision qq
      external qq
      double precision ul
      double precision ur
      double precision wquad(nquad)
      double precision wquade
      double precision x
      double precision xleft
      double precision xn(0:n)
      double precision xquad(nquad,nmax)
      double precision xquade
      double precision xrite
c
c  Zero out the entries.
c
      do i = 1, nu
        f(i) = 0.0D+00
        aleft(i) = 0.0D+00
        arite(i) = 0.0D+00
        adiag(i) = 0.0D+00
      end do
c
c  For each interval,
c
      do ie = 1, n

        he = h(ie)
        xleft = xn(node(1,ie))
        xrite = xn(node(2,ie))
c
c  For each quadrature point in the interval,
c
        do iq = 1, nquad

          xquade = xquad(iq,ie)
          wquade = wquad(iq)
c
c  Pick a basis function which defines the equation,
c
          do il = 1, nl

            ig = node(il,ie)
            iu = indx(ig)

            if ( 0 .lt. iu ) then

              call phi ( il, xquade, phii, phiix, xleft, xrite )

              f(iu) = f(iu) + he * wquade * ff ( xquade ) * phii
c
c  Take care of boundary conditions specifying the value of U'.
c
              if ( ig .eq. 0 ) then

                x = 0.0D+00
                f(iu) = f(iu) - pp ( x ) * ul

              else if ( ig .eq. n ) then

                x = 1.0D+00
                f(iu) = f(iu) + pp ( x ) * ur

              end if
c
c  Pick a basis function which defines the coefficient being computed.
c
              do jl = 1, nl

                jg = node(jl,ie)
                ju = indx(jg)
                call phi ( jl, xquade, phij, phijx, xleft, xrite )

                aij = he * wquade *
     &            ( pp ( xquade ) * phiix * phijx 
     &            + qq ( xquade ) * phii * phij )
c
c  Decide where the coefficient is to be added.
c
                if ( ju .le. 0 ) then

                  if ( jg .eq. 0 ) then
                    f(iu) = f(iu) - aij * ul
                  else if ( jg .eq. n ) then
                    f(iu) = f(iu)-aij*ur
                  end if

                else if ( iu .eq. ju ) then
                  adiag(iu) = adiag(iu) + aij
                else if ( ju .lt. iu ) then
                  aleft(iu) = aleft(iu) + aij
                else
                  arite(iu) = arite(iu) + aij
                end if

              end do

            end if

          end do

        end do

      end do

      return
      end
      function ff ( x )

c*********************************************************************72
c
cc FF evaluates the right hand side function F(X).
c
c  Discussion:
c
c    The function F(X) occurs on the right hand side of the equation:
c
c    - d/dx (p du/dx) + q u = f
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the evaluation point.
c
c    Output, double precision FF, the value of FF(X).
c
      implicit none

      double precision alpha
      double precision beta
      double precision ff
      integer problem
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision x
c
c  Find out which problem we're working on.
c
      call get_problem ( problem )

      if ( problem .eq. 1 ) then

        ff = 0.0D+00

      else if ( problem .eq. 2 ) then

        ff = -2.0D+00 * x

      else if ( problem .eq. 3 ) then

        ff = 0.25D+00 * pi * pi * sin ( 0.5D+00 * pi * x )

      else if ( problem .eq. 4 ) then

        ff = 0.25D+00 * pi * pi * cos ( 0.5D+00 * pi * x )

      else if ( problem .eq. 5 ) then

        call get_beta ( beta )
        ff = - ( x**beta ) + ( x**( beta + 2.0D+00 ) ) 
     &    / ( ( beta + 2.0D+00 ) * ( beta + 1.0D+00 ) )

      else if ( problem .eq. 6 ) then

        call get_alpha ( alpha )
        ff = 2.0D+00 * alpha * ( x - 0.5D+00 ) 
     &    / ( alpha**2 + ( x - 0.5D+00 )**2 )**2

      end if

      return
      end
      subroutine geometry ( h, ibc, indx, n, nl, nmax, node, nquad,
     &  nu, wquad, xn, xquad )

c*********************************************************************72
c
cc GEOMETRY sets geometric information for the problem.  
c
c  Discussion:
c
c    Note, however, that the location of the nodes
c    is done outside of this routine, and, in fact, before this
c    routine is called.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision H(N)
c    H(I) is the length of subinterval I.  This code uses
c    equal spacing for all the subintervals.
c
c    Input, integer IBC.
c    IBC declares what the boundary conditions are.
c    1, at the left endpoint, U has the value UL,
c       at the right endpoint, U' has the value UR.
c    2, at the left endpoint, U' has the value UL,
c       at the right endpoint, U has the value UR.
c    3, at the left endpoint, U has the value UL,
c       and at the right endpoint, U has the value UR.
c    4, at the left endpoint, U' has the value UL,
c       at the right endpoint U' has the value UR.
c
c    Output, integer INDX(0:N).
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
c    at node N, which will be unknown N or N+1,
c    depending on whether there was an unknown at node 0.
c
c    Input, integer N.
c    The number of subintervals into which the interval
c    [XL,XR] is broken.
c
c    Input, integer NL.
c    The number of basis functions used in a single
c    subinterval.  (NL-1) is the degree of the polynomials
c    used.  For this code, NL is fixed at 2, meaning that
c    piecewise linear functions are used as the basis.
c
c    Input, integer NMAX, the maximum number of unknowns that can be handled.
c
c    Output, integer NODE(NL,N).
c    For each subinterval I:
c    NODE(1,I) is the number of the left node, and
c    NODE(2,I) is the number of the right node.
c
c    Input, integer NQUAD
c    The number of quadrature points used in a subinterval.
c
c    Output, integer NU.
c    NU is the number of unknowns in the linear system.
c    Depending on the value of IBC, there will be N-1,
c    N, or N+1 unknown values, which are the coefficients
c    of basis functions.
c
c    Output, double precision WQUAD(NQUAD).
c    WQUAD(I) is the weight associated with the I-th point
c    of an NQUAD point Gaussian quadrature rule.
c
c    Input, double precision XN(0:N).
c    XN(I) is the location of the I-th node.  XN(0) is XL,
c    and XN(N) is XR.
c
c    Output, double precision XQUAD(NQUAD,NMAX), the I-th quadrature point
c    in interval J.
c
      implicit none

      integer n
      integer nl
      integer nmax
      integer nquad

      double precision alfa
      double precision h(n)
      integer i
      integer ibc
      integer igl
      integer igr
      integer indx(0:n)
      integer node(nl,nmax)
      integer nu
      double precision wquad(nquad)
      double precision xl
      double precision xn(0:n)
      double precision xquad(nquad,nmax)
      double precision xr
c
c  Store in NODE the fact that interval I has node I-1
c  as its left endpoint, and node I as its right endpoint.
c
      do i = 1, n
        node(1,i) = i-1
        node(2,i) = i
      end do
c
c  For every node that is associated with an unknown, we
c  record the number of the unknown in INDX.
c
      nu = 0

      do i = 0, n

        if ( i .eq. 0 .and. ( ibc .eq. 1 .or. ibc .eq. 3 ) ) then
          indx(i) = -1
        else if ( i .eq. n .and. ( ibc .eq. 2 .or. ibc .eq. 3 ) ) then
          indx(i) = -1
        else
          nu = nu + 1
          indx(i) = nu
        end if

      end do
c
c  We compute the width of each interval.
c
      do i = 1, n
        igl = node(1,i)
        igr = node(2,i)
        h(i) = xn(igr) - xn(igl)
      end do
c
c  We compute the location of the quadrature points in each
c  interval.
c
      do i = 1, n

        xl = xn(node(1,i))
        xr = xn(node(2,i))

        if ( nquad .eq. 1 ) then

          xquad(1,i) = 0.5D+00 * ( xl + xr )

        else if ( nquad .eq. 2 ) then

          alfa = -0.577350D+00

          xquad(1,i) = ( ( 1.0D+00 - alfa ) * xl 
     &                 + ( 1.0D+00 + alfa ) * xr ) 
     &                 /   2.0D+00

          alfa = +0.577350D+00

          xquad(2,i) = ( ( 1.0D+00 - alfa ) * xl 
     &                 + ( 1.0D+00 + alfa ) * xr ) 
     &                 /   2.0D+00

        else if ( nquad .eq. 3 ) then

          alfa = -0.774597D+00
          xquad(1,i) = ( ( 1.0D+00 - alfa ) * xl 
     &                 + ( 1.0D+00 + alfa ) * xr ) 
     &                 /   2.0D+00

          xquad(2,i) = 0.5D+00 * ( xl + xr )

          alfa = +0.774597D+00
          xquad(3,i) = ( ( 1.0D+00 - alfa ) * xl 
     &                 + ( 1.0D+00 + alfa ) * xr ) 
     &                 /   2.0D+00

        end if

      end do
c
c  We store the weights for the quadrature rule.
c
      if ( nquad .eq. 1 ) then
        wquad(1) = 1.0D+00
      else if ( nquad .eq. 2 ) then
        wquad(1) = 0.5D+00
        wquad(2) = 0.5D+00
      else if ( nquad .eq. 3 ) then
        wquad(1) = 4.0D+00 / 9.0D+00
        wquad(2) = 5.0D+00 / 18.0D+00
        wquad(3) = 4.0D+00 / 9.0D+00
      end if

      return
      end
      subroutine get_alpha ( alpha )

c*********************************************************************72
c
cc GET_ALPHA returns the value of ALPHA, for use by problem 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision ALPHA, the value of ALPHA.
c
      implicit none

      double precision alpha

      alpha = 0.01D+00

      return
      end
      subroutine get_beta ( beta )

c*********************************************************************72
c
cc GET_BETA returns the value of BETA, for use by problem 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 November 2006
c
c  Parameters:
c
c    Output, double precision BETA, the value of Beta.
c
      implicit none

      double precision beta

      beta = -0.9D+00

      return
      end
      subroutine get_problem ( problem )

c*********************************************************************72
c
cc GET_PROBLEM returns the value of the current problem number.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer PROBLEM, chooses the problem to be solved.
c    The user must choose this value by setting it in routine get_problem.
c    * 1, u = x, p = 1, q = 0, f = 0, ibc = 3, ul = 0, ur = 1.
c    The program should find the solution exactly, and the
c    adaptive code should find that there is no reason to
c    subdivide any interval.
c    * 2, u = x*x, p = 1, q = 0, f = -2, ibc = 3, ul = 0, ur = 1.
c    This problem should find the solution exactly, and
c    the adaptive code should again find there is nothing
c    to do.
c    *3, u = sin(pi*x/2), p = 1, q = 0, f = 0.25*pi*pi*sin(pi*x/2), 
c    ul = 0, ur = 1.
c    *4, u = cos(pi*x/2), p = 1, q = 0, f = 0.25*pi*pi*cos(pi*x/2), 
c    ul = 1, ur = 0.
c    *5: u = x**(beta+2)/((beta+2)*(beta+1)), p = 1, q = 1,
c    f = -x**beta + (x**(beta+2))/((beta+2)*(beta+1)),
c    ul = 0, ur = 1/((beta+2)*(beta+1))
c    (beta must be greater than -2, and not equal to -1)
c    *6: u = atan((x-0.5)/alpha), p = 1, q = 0,
c    f =  2*alpha*(x-0.5) / (alpha**2 + (x-0.5)**2) **2,
c    ul = u(0), ur = u(1)
c
      implicit none

      integer problem

      problem = 6

      return
      end
      subroutine init ( ibc, n, tol, ul, ur, xl, xn, xr )

c*********************************************************************72
c
cc INIT initializes some parameters that define the problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IBC.
c    IBC declares what the boundary conditions are.
c    1, at the left endpoint, U has the value UL,
c       at the right endpoint, U' has the value UR.
c    2, at the left endpoint, U' has the value UL,
c       at the right endpoint, U has the value UR.
c    3, at the left endpoint, U has the value UL,
c       and at the right endpoint, U has the value UR.
c    4, at the left endpoint, U' has the value UL,
c       at the right endpoint U' has the value UR.
c
c    Input, integer N.
c    The number of subintervals into which the interval
c    [XL,XR] is broken.
c
c    Output, double precision TOL.
c    A tolerance that is used to determine whether the estimated
c    error in an interval is so large that it should be subdivided
c    and the problem solved again.
c
c    Output, double precision UL.
c    If IBC is 1 or 3, UL is the value that U is required
c    to have at X = XL.
c    If IBC is 2 or 4, UL is the value that U' is required
c    to have at X = XL.
c
c    Output, double precision UR.
c    If IBC is 2 or 3, UR is the value that U is required
c    to have at X = XR.
c    If IBC is 1 or 4, UR is the value that U' is required
c    to have at X = XR.
c
c    Output, double precision XL.
c    XL is the left endpoint of the interval over which the
c    differential equation is being solved.
c
c    Output, double precision XN(0:N).
c    XN(I) is the location of the I-th node.  XN(0) is XL,
c    and XN(N) is XR.
c
c    Output, double precision XR.
c    XR is the right endpoint of the interval over which the
c    differential equation is being solved.
c
      implicit none

      integer n

      double precision alpha
      double precision beta
      integer i
      integer ibc
      integer problem
      double precision tol
      double precision u_exact
      external u_exact
      double precision ul
      double precision ur
      double precision xl
      double precision xn(0:n)
      double precision xr

      tol = 0.01D+00
c
c  Find out which problem we're working on.
c
      call get_problem ( problem )
c
c  Set the boundary conditions for the problem, and
c  print out its title.
c
      if ( problem .eq. 1 ) then

        ibc = 3
        ul = 0.0D+00
        ur = 1.0D+00
        xl = 0.0D+00
        xr = 1.0D+00
        write(*,*)' '
        write(*,*)'Exact solution is U = X'

      else if ( problem .eq. 2 ) then

        ibc = 3
        ul = 0.0D+00
        ur = 1.0D+00
        xl = 0.0D+00
        xr = 1.0D+00
        write(*,*)' '
        write(*,*)'Exact solution is U = X*X'

      else if ( problem .eq. 3 ) then

        ibc = 3
        ul = 0.0D+00
        ur = 1.0D+00
        xl = 0.0D+00
        xr = 1.0D+00
        write(*,*)' '
        write(*,*)'Exact solution is U = SIN(PI*X/2)'

      else if ( problem .eq. 4 ) then

        ibc = 3
        ul = 1.0D+00
        ur = 0.0D+00
        xl = 0.0D+00
        xr = 1.0D+00
        write(*,*)' '
        write(*,*)'Exact solution is U = COS(PI*X/2)'

      else if ( problem .eq. 5 ) then

        ibc = 3
        call get_beta ( beta )
        ul = 0.0D+00
        ur = 1.0D+00 / ( ( beta + 2.0D+00 ) * ( beta + 1.0D+00 ) )
        xl = 0.0D+00
        xr = 1.0D+00
        write(*,*)' '
        write(*,*)'Rheinboldt problem'

      else if ( problem .eq. 6 ) then

        ibc = 3
        call get_alpha ( alpha )
        xl = 0.0D+00
        xr = 1.0D+00
        ul = u_exact ( xl )
        ur = u_exact ( xr )
        write(*,*)' '
        write(*,*)'Arctangent problem'

      end if
c
c  The nodes are defined here, and not in the GEOMETRY routine.
c  This is because each new iteration chooses the location
c  of the new nodes in a special way.
c
      do i = 0, n
        xn(i) = ( dble ( n - i ) * xl 
     &          + dble (     i ) * xr )
     &          / dble ( n     )
      end do

      write(*,*)'The equation is to be solved for '
      write(*,*)'X greater than ', xl
      write(*,*)' and less than ', xr
      write(*,*)' '
      write(*,*)'The boundary conditions are:'
      write(*,*)' '

      if ( ibc .eq. 1 .or. ibc .eq. 3 ) then
        write(*,*)'  At X = XL, U = ', ul
      else
        write(*,*)'  At X = XL, U'' = ', ul
      end if

      if ( ibc .eq. 2 .or. ibc .eq. 3 ) then
        write(*,*)'  At X = XR, U = ', ur
      else
        write(*,*)'  At X = XR, U'' = ', ur
      end if

      return
      end
      subroutine output ( f, ibc, indx, n, nu, ul, ur, xn )

c*********************************************************************72
c
cc OUTPUT prints out the computed solution.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision F(NU).
c    ASSEMBLE stores into F the right hand side of the linear
c    equations.
c    SOLVE replaces those values of F by the solution of the
c    linear equations.
c
c    Input, integer IBC.
c    IBC declares what the boundary conditions are.
c    1, at the left endpoint, U has the value UL,
c       at the right endpoint, U' has the value UR.
c    2, at the left endpoint, U' has the value UL,
c       at the right endpoint, U has the value UR.
c    3, at the left endpoint, U has the value UL,
c       and at the right endpoint, U has the value UR.
c    4, at the left endpoint, U' has the value UL,
c       at the right endpoint U' has the value UR.
c
c    Input, integer INDX(0:N).
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
c    at node N, which will be unknown N or N+1,
c    depending on whether there was an unknown at node 0.
c
c    Input, integer N.
c    The number of subintervals into which the interval
c    [XL,XR] is broken.
c
c    Input, integer NU.
c    NU is the number of unknowns in the linear system.
c    Depending on the value of IBC, there will be N-1,
c    N, or N+1 unknown values, which are the coefficients
c    of basis functions.
c
c    Input, double precision UL.
c    If IBC is 1 or 3, UL is the value that U is required
c    to have at X = XL.
c    If IBC is 2 or 4, UL is the value that U' is required
c    to have at X = XL.
c
c    Input, double precision UR.
c    If IBC is 2 or 3, UR is the value that U is required
c    to have at X = XR.
c    If IBC is 1 or 4, UR is the value that U' is required
c    to have at X = XR.
c
c    Input, double precision XN(0:N).
c    XN(I) is the location of the I-th node.  XN(0) is XL,
c    and XN(N) is XR.
c
      implicit none

      integer n
      integer nu

      double precision error
      double precision f(nu)
      integer i
      integer ibc
      integer indx(0:n)
      double precision u
      double precision uex
      double precision u_exact
      external u_exact
      double precision ul
      double precision ur
      double precision xn(0:n)

      write(*,*)' '
      write(*,*)'Node    X(I)        U(X(I))        U exact   '
     &  // '     Error'
      write(*,*)' '

      do i = 0, n

        if ( i .eq. 0 ) then

          if ( ibc .eq. 1 .or. ibc .eq. 3 ) then
            u = ul
          else
            u = f(indx(i))
          end if

        else if ( i .eq. n ) then

          if ( ibc .eq. 2 .or. ibc .eq. 3 ) then
            u = ur
          else
            u = f(indx(i))
          end if

        else

          u = f(indx(i))

        end if

        uex = u_exact ( xn(i) )
        error = u - uex

        write(*,'(i4,4g14.6)') i, xn(i), u, uex, error

      end do

      return
      end
      subroutine phi ( il, x, phii, phiix, xleft, xrite )

c*********************************************************************72
c
cc PHI evaluates a particular finite element basis function.
c
c  Discussion:
c
c    This routine evaluates a linear basis function and its 
c    derivative at a point X in an interval.  In any
c    interval, there are just two basis functions.  The first
c    basis function is a line which is 1 at the left endpoint
c    and 0 at the right.  The second basis function is 0 at
c    the left endpoint and 1 at the right.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer IL, the local index of the basis function.
c
c    Input, double precision X, the evaluation point.
c
c    Output, double precision PHII, PHIIX, the value of the basis
c    function and its derivative.
c
c    Input, double precision XLEFT, XRITE, the left and right
c    endpoints of the interval.
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
          phii = ( x - xleft ) / ( xrite - xleft )
          phiix = 1.0D+00 / ( xrite - xleft )
        end if
c
c  If X is outside of the interval, then the basis function
c  is always zero.
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
cc PP evaluates the function P in the differential equation.
c
c  Discussion:
c
c    P(X) is the coefficient function that occurs in the equation:
c
c    - d/dx ( p du/dx ) + q u = f
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the evaluation point.
c
c    Output, double precision PP, the value of P(X).
c
      implicit none

      integer problem
      double precision pp
      double precision x
c
c  Find out which problem we're working on.
c
      call get_problem ( problem )

      if ( problem .eq. 1 ) then
        pp = 1.0D+00
      else if ( problem .eq. 2 ) then
        pp = 1.0D+00
      else if ( problem .eq. 3 ) then
        pp = 1.0D+00
      else if ( problem .eq. 4 ) then
        pp = 1.0D+00
      else if ( problem .eq. 5 ) then
        pp = 1.0D+00
      else if ( problem .eq. 6 ) then
        pp = 1.0D+00
      end if

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
c    03 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ADIAG(NU).
c    ADIAG(I) is the "diagonal" coefficient of the I-th
c    equation in the linear system.  That is, ADIAG(I) is
c    the coefficient of the I-th unknown in the I-th equation.
c
c    Input, double precision ALEFT(NU).
c    ALEFT(I) is the "left hand" coefficient of the I-th
c    equation in the linear system.  That is, ALEFT(I) is the
c    coefficient of the (I-1)-th unknown in the I-th equation.
c    There is no value in ALEFT(1), since the first equation
c    does not refer to a "0-th" unknown.
c
c    Input, double precision ARITE(NU).
c    ARITE(I) is the "right hand" coefficient of the I-th
c    equation in the linear system.  ARITE(I) is the coefficient
c    of the (I+1)-th unknown in the I-th equation.  There is
c    no value in ARITE(NU) because the NU-th equation does not
c    refer to an "NU+1"-th unknown.
c
c    Input, double precision F(NU).
c    ASSEMBLE stores into F the right hand side of the linear
c    equations.
c    SOLVE replaces those values of F by the solution of the
c    linear equations.
c
c    Input, integer NU.
c    NU is the number of unknowns in the linear system.
c    Depending on the value of IBC, there will be N-1,
c    N, or N+1 unknown values, which are the coefficients
c    of basis functions.
c
      implicit none

      integer nu

      double precision adiag(nu)
      double precision aleft(nu)
      double precision arite(nu-1)
      double precision f(nu)
      integer i

      write(*,*)' '
      write(*,*)'Printout of tridiagonal linear system:'
      write(*,*)' '
      write(*,*)'Equation  A-Left  A-Diag  A-Rite  RHS'
      write(*,*)' '
      do i = 1, nu

        if ( i .eq. 1 ) then
          write(*,'(i3,14x,3g14.6)') i, adiag(i), arite(i), f(i)
        else if ( i .lt. nu ) then
          write(*,'(i3,4g14.6)') i, aleft(i), adiag(i), arite(i), f(i)
        else
          write(*,'(i3,2g14.6,14x,g14.6)') i, aleft(i), adiag(i), f(i)
        end if

      end do

      return
      end
      function qq ( x )

c***********************************************************************
c
cc QQ evaluates the function Q in the differential equation.
c
c  Discussion:
c
c    Q(X) is the coefficient that occurs in the equation:
c
c    - d/dx ( p du/dx ) + q u = f
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the evaluation point.
c
c    Output, double precision QQ, the value of Q(X).
c
      implicit none

      integer problem
      double precision qq
      double precision x
c
c  Find out which problem we're working on.
c
      call get_problem(problem)

      if ( problem .eq. 1 ) then
        qq = 0.0D+00
      else if ( problem .eq. 2 ) then
        qq = 0.0D+00
      else if ( problem .eq. 3 ) then
        qq = 0.0D+00
      else if ( problem .eq. 4 ) then
        qq = 0.0D+00
      else if ( problem .eq. 5 ) then
        qq = 1.0D+00
      else if ( problem .eq. 6 ) then
        qq = 0.0D+00
      end if

      return
      end
      subroutine solve ( adiag, aleft, arite, f, nu )

c*********************************************************************72
c
cc SOLVE solves a tridiagonal matrix system of the form A*x = b.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 November 2006
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
c    system to be solve
c    On output, F contains the solution of the linear system.
c
c    Input, INTEGER NU.
c    NU is the number of equations to be solved.
c
      implicit none

      integer nu

      double precision adiag(nu)
      double precision aleft(nu)
      double precision arite(*)
      double precision f(nu)
      integer i
c
c  Handle the special case of a single equation.
c
      if ( nu .eq. 1 ) then

        f(1) = f(1) / adiag(1)
c
c  The general case, when NU is greater than 1.
c
      else

        arite(1) = arite(1) / adiag(1)
        do i = 2, nu-1
          adiag(i) = adiag(i) - aleft(i) * arite(i-1)
          arite(i) = arite(i) / adiag(i)
        end do
        adiag(nu) = adiag(nu) - aleft(nu) * arite(nu-1)

        f(1) = f(1) / adiag(1)
        do i = 2, nu
          f(i) = ( f(i) - aleft(i) * f(i-1) ) / adiag(i)
        end do

        do i = nu-1, 1, -1
          f(i) = f(i) - arite(i) * f(i+1)
        end do

      end if
 
      return
      end
      subroutine solvex ( adiag, aleft, arite, f, h, ibc, indx,
     &  kount, n, nl, nmax, node, nquad, nu, ul, ur, wquad, xn,
     &  xquad )

c*********************************************************************72
c
cc SOLVEX solves the finite element system.
c
c  Discussion:
c
c    For the given set of nodes XN, determined externally, this 
c    routine sets up the associated finite element linear system 
c    and solves it for the coefficients.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Workspace, double precision ADIAG(NU).
c    ADIAG(I) is the "diagonal" coefficient of the I-th
c    equation in the linear system.  That is, ADIAG(I) is
c    the coefficient of the I-th unknown in the I-th equation.
c
c    Workspace, double precision ALEFT(NU).
c    ALEFT(I) is the "left hand" coefficient of the I-th
c    equation in the linear system.  That is, ALEFT(I) is the
c    coefficient of the (I-1)-th unknown in the I-th equation.
c    There is no value in ALEFT(1), since the first equation
c    does not refer to a "0-th" unknown.
c
c    Workspace, double precision ARITE(NU).
c    ARITE(I) is the "right hand" coefficient of the I-th
c    equation in the linear system.  ARITE(I) is the coefficient
c    of the (I+1)-th unknown in the I-th equation.  There is
c    no value in ARITE(NU) because the NU-th equation does not
c    refer to an "NU+1"-th unknown.
c
c    Output, double precision F(NU), used for the right hand side,
c    and then solution, of a linear system.
c
c    Output, double precision H(N)
c    H(I) is the length of subinterval I.  This code uses
c    equal spacing for all the subintervals.
c
c    Input, integer IBC.
c    IBC declares what the boundary conditions are.
c    1, at the left endpoint, U has the value UL,
c       at the right endpoint, U' has the value UR.
c    2, at the left endpoint, U' has the value UL,
c       at the right endpoint, U has the value UR.
c    3, at the left endpoint, U has the value UL,
c       and at the right endpoint, U has the value UR.
c    4, at the left endpoint, U' has the value UL,
c       at the right endpoint U' has the value UR.
c
c    Output, integer INDX(0:N).
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
c    at node N, which will be unknown N or N+1,
c    depending on whether there was an unknown at node 0.
c
c    Input, integer KOUNT, the number of adaptive steps that have been taken.
c
c    Input, integer N.
c    The number of subintervals into which the interval
c    [XL,XR] is broken.
c
c    Input, integer NL.
c    The number of basis functions used in a single
c    subinterval.  (NL-1) is the degree of the polynomials
c    used.  For this code, NL is fixed at 2, meaning that
c    piecewise linear functions are used as the basis.
c
c    Input, integer NMAX, the maximum number of unknowns that can be handled.
c
c    Output, integer NODE(NL,N).
c    For each subinterval I:
c    NODE(1,I) is the number of the left node, and
c    NODE(2,I) is the number of the right node.
c
c    Input, integer NQUAD
c    The number of quadrature points used in a subinterval.
c
c    Output, integer NU.
c    NU is the number of unknowns in the linear system.
c    Depending on the value of IBC, there will be N-1,
c    N, or N+1 unknown values, which are the coefficients
c    of basis functions.
c
c    Input, double precision UL.
c    If IBC is 1 or 3, UL is the value that U is required
c    to have at X = XL.
c    If IBC is 2 or 4, UL is the value that U' is required
c    to have at X = XL.
c
c    Input, double precision UR.
c    If IBC is 2 or 3, UR is the value that U is required
c    to have at X = XR.
c    If IBC is 1 or 4, UR is the value that U' is required
c    to have at X = XR.
c
c    Output, double precision WQUAD(NQUAD).
c    WQUAD(I) is the weight associated with the I-th point
c    of an NQUAD point Gaussian quadrature rule.
c
c    Input, double precision XN(0:N).
c    XN(I) is the location of the I-th node.  XN(0) is XL,
c    and XN(N) is XR.
c
c    Output, double precision XQUAD(NQUAD,NMAX), the I-th quadrature point
c    in interval J.
c
      implicit none

      integer n
      integer nl
      integer nmax
      integer nquad

      double precision adiag(nmax)
      double precision aleft(nmax)
      double precision arite(nmax)
      double precision f(nmax)
      double precision h(n)
      integer ibc
      integer indx(0:n)
      integer kount
      integer node(nl,nmax)
      integer nu
      double precision ul
      double precision ur
      double precision wquad(nquad)
      double precision xn(0:n)
      double precision xquad(nquad,nmax)
c
c  Given a set of N nodes (where N increases on each iteration),
c  compute the other geometric information.
c
      call geometry ( h, ibc, indx, n, nl, nmax, node, nquad, 
     &  nu, wquad, xn, xquad )
c
c  Assemble the linear system.
c
      call assemble ( adiag, aleft, arite, f, h, n, indx, node, 
     &  nu, nl, nquad, nmax, ul, ur, wquad, xn, xquad )
c
c  Print out the linear system, just once.
c
      if ( kount .eq. 1 ) then
        call prsys ( adiag, aleft, arite, f, nu )
      end if
c
c  Solve the linear system.
c
      call solve ( adiag, aleft, arite, f, nu )
c
c  Print out the solution.
c
      write(*,*)' '
      write(*,*)'Basic solution'

      call output ( f, ibc, indx, n, nu, ul, ur, xn )

      return
      end
      subroutine solvey ( eta, f, h, n, nu, ul, ur, xn )

c*********************************************************************72
c
cc SOLVEY computes error estimators for a computed solution.
c
c  Discussion:
c
c    SOLVEY accepts information about the solution of a finite element
c    problem on a grid of nodes with coordinates XN.  It then starts
c    at node 0, and for each node, computes two "error estimators",
c    one for the left, and one for the right interval associated with the
c    node.  These estimators are found by solving a finite element problem
c    over the two intervals, using the known values of the original
c    solution as boundary data, and using a mesh that is "slightly"
c    refined over the original one.
c
c    The computations at the first and last nodes only involve
c    a single interval.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision ETA(N).
c    ETA(I) is the error estimate for interval I.  It is computed
c    as the sum of two quantities, one associated with the left
c    and one with the right node of the interval.
c
c    Input, double precision F(NU).
c    ASSEMBLE stores into F the right hand side of the linear
c    equations.
c    SOLVE replaces those values of F by the solution of the
c    linear equations.
c
c    Input, double precision H(N)
c    H(I) is the length of subinterval I.  This code uses
c    equal spacing for all the subintervals.
c
c    Input, integer N.
c    The number of subintervals into which the interval
c    [XL,XR] is broken.
c
c    Input, integer NU.
c    NU is the number of unknowns in the linear system.
c    Depending on the value of IBC, there will be N-1,
c    N, or N+1 unknown values, which are the coefficients
c    of basis functions.
c
c    Input, double precision UL.
c    If IBC is 1 or 3, UL is the value that U is required
c    to have at X = XL.
c    If IBC is 2 or 4, UL is the value that U' is required
c    to have at X = XL.
c
c    Input, double precision UR.
c    If IBC is 2 or 3, UR is the value that U is required
c    to have at X = XR.
c    If IBC is 1 or 4, UR is the value that U' is required
c    to have at X = XR.
c
c    Input, double precision XN(0:N).
c    XN(I) is the location of the I-th node.  XN(0) is XL,
c    and XN(N) is XR.
c
      implicit none

      integer nl
      parameter (nl = 2)

      integer ny
      parameter (ny = 2)

      integer nquad
      parameter (nquad = 2)

      integer nmay
      parameter (nmay = 2*ny)

      integer n
      integer nu

      double precision adiag(nmay)
      double precision aleft(nmay)
      double precision arite(nmay)
      double precision eta(n)
      double precision f(nu)
      double precision fy(nmay)
      double precision h(n)
      double precision hy(nmay)
      integer i
      integer ibcy
      integer indy(0:nmay)
      integer j
      integer jhi
      integer jlo
      integer jmid
      integer k
      integer m
      integer nodey(nl,nmay)
      integer nuy
      double precision pp
      external pp
      double precision qq
      external qq
      double precision total
      double precision ul
      double precision uleft
      double precision ulval
      double precision uly
      double precision uprime
      double precision ur
      double precision urite
      double precision urval
      double precision ury
      double precision uval
      double precision vlval
      double precision vprime
      double precision vrval
      double precision vval
      double precision wquad(nquad)
      double precision xn(0:n)
      double precision xquady(nquad,nmay)
      double precision y
      double precision yl
      double precision ym
      double precision yn(0:nmay)
      double precision yr
c
c  Initialize the error estimators to zero.
c
      do j = 1, n
        eta(j) = 0.0D+00
      end do
c
c  Set the boundary conditions for each subproblem to be
c  known values of U at the left and right.
c
c
c  For each node, subdivide its left and right hand intervals
c  into NY subintervals.
c
c  Set up and solve the differential equation again on this
c  smaller region.
c
c  The 0-th and N-th nodes are special cases.
c
      ibcy = 3

      do j = 0, n

        if ( j .eq. 0 ) then
          m = ny
          jlo = j
          jmid = j + 1
          jhi = j + 1
        else if ( j .eq. n ) then
          m = ny
          jlo = j - 1
          jmid = j
          jhi = j
        else
          m = 2 * ny
          jlo = j - 1
          jmid = j
          jhi = j+1
        end if
c
c  Set the location of the nodes in the subintervals.
c
        yl = xn(jlo)
        ym = xn(jmid)
        yr = xn(jhi)

        do i = 0, ny
          yn(i) = ( dble ( ny - i ) * yl
     &            + dble (      i ) * ym )
     &            / dble ( ny     )
        end do

        do i = ny+1, m
          yn(i) = ( dble ( m - i      ) * ym
     &            + dble (     i - ny ) * yr )
     &            / dble ( m -     ny )
        end do
c
c  Set up the geometry of the sub-problem.
c
        call geometry ( hy, ibcy, indy, m, nl, nmay, nodey, 
     &    nquad, nuy, wquad, yn, xquady )
c
c  Set the boundary values for the sub-problem.
c
        if ( j .le. 1 ) then
          uly = ul
        else
          uly = f(j-1)
        end if

        if ( n - 1 .le. j ) then
          ury = ur
        else
          ury = f(j+1)
        end if
c
c  Assemble the matrix for the sub-problem.
c
        call assemble ( adiag, aleft, arite, fy, hy, m, indy,
     &    nodey, nuy, nl, nquad, nmay, uly, ury, wquad, yn,
     &    xquady )
c
c  Solve the system.
c
        call solve ( adiag, aleft, arite, fy, nuy )
c
c  Compute the weighted sum of the squares of the differences
c  of the original computed slope and the refined computed slopes.
c
c  Calculation for left interval.
c
        if ( 1 .le. j ) then

          if ( j .le. 1 ) then
            uleft = ul
            urite = f(1)
          else if ( j .eq. n ) then
            uleft = f(j-1)
            urite = ur
          else
            uleft = f(j-1)
            urite = f(j)
          end if

          uprime = ( urite - uleft ) / h(j)

          total = 0.0D+00

          do i = 1, ny

            yl = yn(i-1)
            yr = yn(i)

            if ( i .eq. 1 ) then
              vlval = uly
              vrval = fy(i)
            else if ( i .eq. m ) then
              vlval = fy(i-1)
              vrval = ury
            else
              vlval = fy(i-1)
              vrval = fy(i)
            end if

            vprime = ( vrval - vlval ) / hy(i)

            ulval = ( dble ( ny - i + 1 ) * uleft
     &              + dble (      i - 1 ) * urite )
     &              / dble ( ny         )

            urval = ( dble ( ny - i ) * uleft
     &              + dble (      i ) * urite )
     &              / dble ( ny     )
c
c  Compute the integral of
c
c    p(x)*(u'(x)-v'(x))**2 + q(x)*(u(x)-v(x))**2
c
            do k = 1, nquad

              y = xquady(k,i)

              uval = ( ( yl - y      ) * urval
     &               + (      y - yr ) * ulval )
     &               / ( yl     - yr )

              vval = ( ( yl - y      ) * vrval
     &               + (      y - yr ) * vlval )
     &               / ( yl     - yr )

              total = total + 0.5D+00 * wquad(k) * hy(i) *
     &          ( pp ( y ) * ( uprime - vprime )**2
     &          + qq ( y ) * ( uval - vval )**2 )

            end do

          end do

          eta(j) = eta(j) + 0.5D+00 * sqrt ( total )

        end if
c
c  Calculation for right interval.
c
        if ( j .le. n - 1 ) then

          if ( j .eq. 0 ) then
            uleft = ul
            urite = f(j+1)
          else if ( n-1 .le. j ) then
            uleft = f(j)
            urite = ur
          else
            uleft = f(j)
            urite = f(j+1)
          end if

          uprime = ( urite - uleft ) / h(j+1)

          total = 0.0D+00

          do i = m+1-ny, m

            yl = yn(i-1)
            yr = yn(i)

            if ( i .eq. 1 ) then
              vlval = uly
              vrval = fy(i)
            else if ( i .eq. m ) then
              vlval = fy(i-1)
              vrval = ury
            else
              vlval = fy(i-1)
              vrval = fy(i)
            end if

            vprime = ( vrval - vlval ) / hy(i)

            ulval = ( dble (      m - i + 1 ) * uleft
     &              + dble ( ny - m + i - 1 ) * urite )
     &              / dble ( ny )

            urval = ( dble (      m - i ) * uleft
     &              + dble ( ny - m + i ) * urite )
     &              / dble ( ny )
c
c  Compute the integral of
c
c    p(x)*(u'(x)-v'(x))**2 + q(x)*(u(x)-v(x))**2
c
            do k = 1, nquad

              y = xquady(k,i)

              uval = ( ( yl - y      ) * urval
     &               + (      y - yr ) * ulval ) 
     &               / ( yl     - yr )

              vval = ( ( yl - y      ) * vrval
     &               + (      y - yr ) * vlval )
     &               / ( yl     - yr )

              total = total + 0.5D+00 * wquad(k) * hy(i) *
     &          ( pp ( y ) * ( uprime - vprime )**2
     &          + qq ( y ) * ( uval - vval )**2 )

            end do

          end do

          eta(j+1) = eta(j+1) + 0.5D+00 * sqrt ( total )

        end if

      end do
c
c  Print out the error estimators.
c
      write(*,*)' '
      write(*,*)'ETA'
      write(*,*)' '
      do j = 1, n
        write(*,*) eta(j)
      end do

      return
      end
      subroutine subdiv ( eta, jadd, kount, n, nmax, tol, xn, xt, 
     &  status )

c*********************************************************************72
c
cc SUBDIV decides which intervals should be subdivided.
c
c  Discussion:
c
c    The ETA vector contains the information that is used to decide
c    where to subdivide.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision ETA(N).
c    ETA(I) is the error estimate for interval I.  It is computed
c    as the sum of two quantities, one associated with the left
c    and one with the right node of the interval.
c
c    Output, integer JADD(N).
c    JADD(I) is 1 if the error estimates show that interval I
c    should be subdivided.
c
c    Input, integer KOUNT, the number of adaptive steps that have been taken.
c
c    Input/output, integer N.
c    The number of subintervals into which the interval
c    [XL,XR] is broken.  On output, N has been increased.
c
c    Input, integer NMAX, the maximum number of unknowns that can be handled.
c
c    Input, double precision TOL.
c    A tolerance that is used to determine whether the estimated
c    error in an interval is so large that it should be subdivided
c    and the problem solved again.
c
c    Input/output, double precision XN(0:N).
c    XN(I) is the location of the I-th node.  XN(0) is XL,
c    and XN(N) is XR.  On output, the nodes have been updated.
c
c    Workspace, double precision XT(0:NMAX).
c
c    Output, integer STATUS, reports status of subdivision.
c    0, a new subdivision was carried out.
c    1, no more subdivisions are needed.
c    -1, no more subdivisions can be carried out.
c
      implicit none

      integer n

      double precision ave
      double precision eta(n)
      integer j
      integer jadd(n)
      integer k
      integer kount
      integer nmax
      integer status
      double precision temp
      double precision tol
      double precision total
      double precision xn(0:nmax)
      double precision xt(0:nmax)

      status = 0
c
c  Add up the ETA's, and get their average.
c
      total = 0.0D+00
      do j = 1, n
        total = total + eta(j)
      end do
      ave = total / dble ( n )
c
c  Look for intervals whose ETA value is relatively large,
c  and note in JADD that these intervals should be subdivided.
c
      k = 0
      temp = max ( 1.2D+00 * ave + 0.00001D+00, tol**2 / dble ( n ))

      write(*,*)' '
      write(*,*)'Tolerance = ', temp
      write(*,*)' '

      do j = 1, n

        if ( temp .lt. eta(j) ) then
          k = k+1
          jadd(j) = 1
          write(*,*)'Subdivide interval ', j
        else
          jadd(j) = 0
        end if

      end do
c
c  If no subdivisions needed, we're done.
c
      if ( k .le. 0 ) then
        write(*,*)'Success on step ', kount
        status = 1
        return
      end if
c
c  See if we're about to go over our limit.
c
      if ( nmax .lt. n + k ) then
        write(*,*)' '
        write(*,*)'The iterations did not reach their goal.'
        write(*,*)'The next value of N is ', n + k
        write(*,*)'which exceeds NMAX = ', nmax
        status = -1
        return
      end if
c
c  Insert new nodes where needed.
c
      k = 0
      xt(0) = xn(0)

      do j = 1, n

        if ( 0 .lt. jadd(j) ) then
          xt(j+k) = 0.5D+00 * ( xn(j) + xn(j-1) )
          k = k + 1
        end if

        xt(j+k) = xn(j)

      end do
c
c  Update the value of N, and copy the new nodes into XN.
c
      n = n + k
      do j = 0, n
        xn(j) = xt(j)
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
      function u_exact ( x )

c*********************************************************************72
c
cc U_EXACT returns the value of the exact solution at any point X.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 November 2006
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the evaluation point.
c
c    Output, double precision U_EXACT, the value of the exact solution.
c
      implicit none

      double precision alpha
      double precision beta
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      integer problem
      double precision u_exact
      double precision x
c
c  Find out which problem we're working on.
c
      call get_problem ( problem )

      if ( problem .eq. 1 ) then
        u_exact = x
      else if ( problem .eq. 2 ) then
        u_exact = x * x
      else if ( problem .eq. 3 ) then
        u_exact = sin ( pi * x / 2.0D+00 )
      else if ( problem .eq. 4 ) then
        u_exact = cos ( pi * x / 2.0D+00 )
      else if ( problem .eq. 5 ) then
        call get_beta ( beta )
        u_exact = ( x**( beta + 2.0D+00 ) ) 
     &    / ( ( beta + 2.0D+00 ) * ( beta + 1.0D+00 ) )
      else if ( problem .eq. 6 ) then
        call get_alpha ( alpha )
        u_exact = atan ( ( x - 0.5D+00 ) / alpha )
      end if

      return
      end
