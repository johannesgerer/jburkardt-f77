      program main

c*********************************************************************72
c
cc MAIN is the main program for FEM1D_NONLINEAR.
c
c  Discussion:
c
c    FEM1D_NONLINEAR solves the nonlinear one dimensional boundary value
c    problem:
c
c      -d/dx (p(x) du/dx) + q(x)*u + u*u' = f(x)
c
c    by the finite-element method using piecewise linear basis
c    functions.
c
c    Here U is an unknown scalar function of X defined on the
c    interval [XL,XR], and P, Q and F are given functions of X.
c
c    The values of U or U' at XL and XR are also specified.
c
c    Sample problem #1:
c
c      u(x)  = x, 
c      p(x)  = 1, 
c      q(x)  = 0, 
c      f(x)  = x,
c      u(0)  = 0,
c      u'(1) = 1.
c
c      The code should solve this problem exactly.
c
c    Sample problem #2:
c
c      u(x)  = 2*(1-cos(0.5*pi*x))/pi, 
c      p(x)  = 1,
c      q(x)  = 0,
c      f(x)  = -0.5*pi*cos(0.5*pi*x) + 2*sin(0.5*pi*x)*(1-cos(0.5*pi*x)/pi
c      u(0)  = 0,
c      u'(1) = 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 April 2007
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
c    double precision F(N+1) or F(NU).
c    ASSEMB stores into F the right hand side of the linear
c    equations.
c    SOLVE replaces those values of F by the solution of the
c    linear equations.
c
c    double precision FOLD(N+1) or FOLD(NU).
c    FOLD contains the value of F from the previous iteration,
c    and is used in ASSEMB to add correction terms to the
c    matrix and right hand side.
c
c    double precision H(N)
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
c    integer IMAX.
c    The number of iterations to carry out.
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
c    integer NL.
c    The number of basis functions used in a single
c    subinterval.  (NL-1) is the degree of the polynomials
c    used.  For this code, NL is fixed at 2, meaning that
c    piecewise linear functions are used as the basis.
c
c    integer NODE(NL,N).
c    For each subinterval I:
c    NODE(1,I) is the number of the left node, and
c    NODE(2,I) is the number of the right node.
c
c    integer NPRINT.
c    The number of points at which the computed solution
c    should be printed out when compared to the exact
c    solution.
c
c    integer NQUAD.
c    The number of quadrature points used in a subinterval.
c    This code uses NQUAD=1.
c
c    integer N.
c    The number of subintervals into which the interval
c    [XL,XR] is broken.
c
c    integer NU.
c    NU is the number of unknowns in the linear system.
c    Depending on the value of IBC, there will be N-1,
c    N, or N+1 unknown values, which are the coefficients
c    of basis functions.
c
c    integer PROBLEM, indicates which problem to be solved.
c    * 1, u = x, p=1, q=0, f=x, u(0)=0, u'(1)=1.
c    * 2, u = 2*(1-cos(0.5*pi*x))/pi, p=1, q=0,
c    f = -0.5*pi*cos(0.5*pi*x) + 2*sin(0.5*pi*x)*(1-cos(0.5*pi*x)/pi
c    u(0) = 0, u'(1)=1.
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
c    double precision XN(0:N).
c    XN(I) is the location of the I-th node.  XN(0) is XL,
c    and XN(N) is XR.
c
c    double precision XQUAD(N)
c    XQUAD(I) is the location of the single quadrature point
c    in interval I.
c
c    double precision XR.
c    XR is the right endpoint of the interval over which the
c    differential equation is being solved.
c
      implicit none

      integer n
      parameter ( n = 10 )

      integer nl
      parameter ( nl = 2 )

      double precision adiag(n+1)
      double precision aleft(n+1)
      double precision arite(n+1)
      double precision f(n+1)
      double precision fold(n+1)
      double precision h(n)
      integer i
      integer ibc
      integer imax
      integer indx(0:n)
      integer j
      integer node(nl,n)
      integer nprint
      integer nquad
      integer nu
      integer problem
      double precision ul
      double precision ur
      double precision xl
      double precision xn(0:n)
      double precision xquad(n)
      double precision xr

      write ( *, '(a)' ) ' '
      call timestamp ( )
      write ( *, '(a)' )' '
      write ( *, '(a)' )'FEM1D_NONLINEAR'
      write ( *, '(a)' )'  FORTRAN77 version.'
      write ( *, '(a)' )' '
      write ( *, '(a)' )'  Solve a nonlinear boundary value problem:'
      write ( *, '(a)' )' '
      write ( *, '(a)' )'  -d/dx (p(x) du/dx) + q(x)*u + u*u'' = f(x)'
      write ( *, '(a)' )' '
      write ( *, '(a)' )'  on an interval [xl,xr], with the values of'
      write ( *, '(a)' )'  u or u'' specified at xl and xr.'
c
c  Initialize variables that define the problem.
c
      call init ( ibc, imax, nprint, nquad, problem, ul, ur, xl, xr )
c
c  Compute the quantities that describe the geometry of the problem.
c
      call geometry ( h, ibc, indx, nl, node, n, nu, xl, xn, xquad, xr )
c
c  Initialize the "previous" solution to 0.
c
      do i = 1, nu
        fold(i) = 0.0D+00
      end do
c
c  Begin the iteration.
c
      do i = 1, imax
c
c  Assemble the matrix and right hand side.
c
        if ( i .le. 3 ) then
          call assemble_picard ( adiag, aleft, arite, f, fold, h, 
     &      indx, n, nl, node, nquad, nu, problem, ul, ur, xn, xquad )
        else
          call assemble_newton ( adiag, aleft, arite, f, fold, h, 
     &      indx, n, nl, node, nquad, nu, problem, ul, ur, xn, xquad )
        end if
c
c  Print out the linear system, just once.
c
        if ( i .eq. 1 ) then
          call prsys ( adiag, aleft, arite, f, nu )
        end if
c
c  Solve the linear system.
c
        call solve ( adiag, aleft, arite, f, nu )
c
c  Print the current solution.
c
        call output ( f, ibc, indx, n, nu, ul, ur, xn )
c
c  Save a copy of the current solution in FOLD.
c
        do j = 1, nu
          fold(j) = f(j)
        end do

      end do
c
c  Compare the solution to the exact solution.
c
      call compare ( f, indx, n, nl, node, nprint, nu, problem, ul, ur, 
     &  xl, xn, xr )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FEM1D_NONLINEAR:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine assemble_newton ( adiag, aleft, arite, f, fold, h, 
     &  indx, n, nl, node, nquad, nu, problem, ul, ur, xn, xquad )

c*********************************************************************72
c
cc ASSEMBLE_NEWTON assembles the Newton linear system.
c
c  Discussion:
c
c    The linear system being solved here is for the Newton correction
c    to an approximate solution of a nonlinear system.
c
c    Thus, we suppose that we have a nonlinear function F(X),
c    and an approximate solution X0.  If we can suppose there is an
c    exact solution X* that is "nearby", and in fact close enough
c    that Taylor's theorem gives us a useful estimate, then we
c    may write:
c
c      F(X*) = F(X0) + F'(X0) * ( X* - X0 ) + Order ( X* - X0 )^2 
c
c    and by rearranging, we get the Newton system (which is only
c    approximately correct):
c
c      F'(X0) * ( X* - X0 ) = - F(X0)
c
c    We solve this system and add the solution to X0 to get a
c    new approximate solution that, we hope, is much improved.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 April 2007
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
c    Output, double precision F(NU), the right hand side of the linear
c    equations.
c
c    Input, double precision FOLD(NU), the solution value 
c    from the previous iteration,
c
c    Input, double precision H(N), the length of the subintervals.  
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
c    Input, integer NL.
c    The number of basis functions used in a single
c    subinterval.  (NL-1) is the degree of the polynomials
c    used.  For this code, NL is fixed at 2, meaning that
c    piecewise linear functions are used as the basis.
c
c    Input, integer NODE(NL,N).
c    For each subinterval I:
c    NODE(1,I) is the number of the left node, and
c    NODE(2,I) is the number of the right node.
c
c    Input, integer NQUAD.
c    The number of quadrature points used in a subinterval.
c    This code uses NQUAD = 1.
c
c    Input, integer NU.
c    NU is the number of unknowns in the linear system.
c    Depending on the value of IBC, there will be N-1,
c    N, or N+1 unknown values, which are the coefficients
c    of basis functions.
c
c    Input, integer PROBLEM, indicates which problem to be solved.
c    * 1, u = x, p=1, q=0, f=x, u(0)=0, u'(1)=1.
c    * 2, u = 2*(1-cos(0.5*pi*x))/pi, p=1, q=0,
c    f = -0.5*pi*cos(0.5*pi*x) + 2*sin(0.5*pi*x)*(1-cos(0.5*pi*x)/pi
c    u(0) = 0, u'(1)=1.
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
c    Input, double precision XQUAD(N)
c    XQUAD(I) is the location of the single quadrature point
c    in interval I.
c
      implicit none

      integer n
      integer nu

      double precision aij
      double precision adiag(nu)
      double precision aleft(nu)
      double precision arite(nu)
      double precision f(nu)
      double precision ff
      double precision fold(nu)
      double precision h(n)
      double precision he
      integer i
      integer ie
      integer ig
      integer il
      integer indx(0:n)
      integer iu
      integer iul
      integer iur
      integer iq
      integer jg
      integer jl
      integer jr
      integer ju
      integer nl
      integer node(nl,n)
      integer nquad
      double precision phii
      double precision phiix
      double precision phij
      double precision phijx
      external pp
      double precision pp
      integer problem
      double precision qq
      external qq
      double precision total
      double precision ul
      double precision uold
      double precision uoldx
      double precision ur
      double precision x
      double precision xleft
      double precision xn(0:n)
      double precision xqe
      double precision xquad(n)
      double precision xrite

      do i = 1, nu
        f(i) = 0.0D+00
        adiag(i) = 0.0D+00
        aleft(i) = 0.0D+00
        arite(i) = 0.0D+00
      end do
c
c  For element IE...
c
      do ie = 1, n

        he = h(ie)
        xleft = xn(node(1,ie))
        xrite = xn(node(2,ie))
c
c  For quadrature point IQ...
c
        do iq = 1, nquad

          xqe = xquad(ie)
c
c  Compute value of U for previous solution.
c
          total = 0.0

          do il = 1, nl

            ig = node(il,ie)
            iu = indx(ig)

            if ( iu .le. 0 ) then

              if ( il .eq. 1 ) then
                total = total + ul
              else
                total = total + ur
              end if

            else
              total = total + fold(iu)
            end if

          end do

          uold = total / dble ( nl )
c
c  Compute value of U' for previous solution.
c
          jl = node(1,ie)
          jr = node(2,ie)
          iul = indx(jl)
          iur = indx(jr)

          if ( iul .le. 0 ) then
            uoldx = ( fold(iur) - ul ) / he
          elseif ( iur .le. 0 ) then
            uoldx = ( ur - fold(iul) ) / he
          else
            uoldx = ( fold(iur) - fold(iul) ) / he
          end if
c
c  For basis function IL...
c
          do il = 1, nl

            ig = node(il,ie)
            iu = indx(ig)

            if ( 0 < iu ) then

              call phi ( il, xqe, phii, phiix, xleft, xrite )

              f(iu) = f(iu) + he * phii * ( ff ( xqe, problem ) 
     &              + uold * uoldx )
c
c  Handle boundary conditions that prescribe the value of U'.
c
              if ( ig .eq. 0 ) then

                x = 0.0D+00
                f(iu) = f(iu) - pp ( x, problem ) * ul

              else if ( ig .eq. n ) then

                x = 1.0D+00
                f(iu) = f(iu) + pp ( x, problem ) * ur

              end if
c
c  For basis function JL...
c
              do jl = 1, nl

                jg = node(jl,ie)
                ju = indx(jg)

                call phi ( jl, xqe, phij, phijx, xleft, xrite )

                aij = he * ( pp ( xqe, problem ) * phiix * phijx
     &                  + qq ( xqe, problem ) * phii * phij
     &                  + uold * phii * phijx
     &                  + uoldx * phij * phii )

                if ( ju .le. 0 ) then

                  if ( jg .eq. 0 ) then
                    f(iu) = f(iu) - aij * ul
                  else if ( jg .eq. n ) then
                    f(iu) = f(iu) - aij * ur
                  end if

                else if ( iu .eq. ju ) then
                  adiag(iu) = adiag(iu) + aij
                else if( ju .lt. iu )then
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
      subroutine assemble_picard ( adiag, aleft, arite, f, fold, h, 
     &  indx, n, nl, node, nquad, nu, problem, ul, ur, xn, xquad )

c*********************************************************************72
c
cc ASSEMBLE_PICARD assembles the Picard linear system.
c
c  Discussion:
c
c    The equation we are trying to solve has the form:
c
c      -d/dx ( p(x) du/dx ) + q(x) * u + u * u' = f(x)
c
c    For the Picard iteration, we need to modify the nonlinear term u * u'
c    so that it is linear in the unknown u, and any other factors of u are
c    lagged.  One way to do this gives us the following equation:
c
c      -d/dx ( p(x) du/dx ) + q(x) * u + u * uold' = f(x)
c
c    where uold is the previous iterate.
c
c    Now we can formulate this system as a (linear) finite element problem
c
c      A * u = rhs
c
c    to be solved for the new approximate solution u.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 April 2007
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
c    Output, double precision F(NU), the right hand side of the linear
c    equations.
c
c    Input, double precision FOLD(NU), the solution value 
c    from the previous iteration,
c
c    Input, double precision H(N), the length of the subintervals.  
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
c    Input, integer NL.
c    The number of basis functions used in a single
c    subinterval.  (NL-1) is the degree of the polynomials
c    used.  For this code, NL is fixed at 2, meaning that
c    piecewise linear functions are used as the basis.
c
c    Input, integer NODE(NL,N).
c    For each subinterval I:
c    NODE(1,I) is the number of the left node, and
c    NODE(2,I) is the number of the right node.
c
c    Input, integer NQUAD.
c    The number of quadrature points used in a subinterval.
c    This code uses NQUAD = 1.
c
c    Input, integer NU.
c    NU is the number of unknowns in the linear system.
c    Depending on the value of IBC, there will be N-1,
c    N, or N+1 unknown values, which are the coefficients
c    of basis functions.
c
c    Input, integer PROBLEM, indicates which problem to be solved.
c    * 1, u = x, p=1, q=0, f=x, u(0)=0, u'(1)=1.
c    * 2, u = 2*(1-cos(0.5*pi*x))/pi, p=1, q=0,
c    f = -0.5*pi*cos(0.5*pi*x) + 2*sin(0.5*pi*x)*(1-cos(0.5*pi*x)/pi
c    u(0) = 0, u'(1)=1.
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
c    Input, double precision XQUAD(N)
c    XQUAD(I) is the location of the single quadrature point
c    in interval I.
c
      implicit none

      integer n
      integer nu

      double precision aij
      double precision adiag(nu)
      double precision aleft(nu)
      double precision arite(nu)
      double precision f(nu)
      double precision ff
      double precision fold(nu)
      double precision h(n)
      double precision he
      integer i
      integer ie
      integer ig
      integer il
      integer indx(0:n)
      integer iu
      integer iul
      integer iur
      integer iq
      integer jg
      integer jl
      integer jr
      integer ju
      integer nl
      integer node(nl,n)
      integer nquad
      double precision phii
      double precision phiix
      double precision phij
      double precision phijx
      external pp
      double precision pp
      integer problem
      double precision qq
      external qq
      double precision total
      double precision ul
      double precision uold
      double precision uoldx
      double precision ur
      double precision x
      double precision xleft
      double precision xn(0:n)
      double precision xqe
      double precision xquad(n)
      double precision xrite

      do i = 1, nu
        f(i) = 0.0D+00
        adiag(i) = 0.0D+00
        aleft(i) = 0.0D+00
        arite(i) = 0.0D+00
      end do
c
c  For element IE...
c
      do ie = 1, n

        he = h(ie)
        xleft = xn(node(1,ie))
        xrite = xn(node(2,ie))
c
c  For quadrature point IQ...
c
        do iq = 1, nquad

          xqe = xquad(ie)
c
c  Compute value of U for previous solution.
c
          total = 0.0

          do il = 1, nl

            ig = node(il,ie)
            iu = indx(ig)

            if ( iu .le. 0 ) then

              if ( il .eq. 1 ) then
                total = total + ul
              else
                total = total + ur
              end if

            else
              total = total + fold(iu)
            end if

          end do

          uold = total / dble ( nl )
c
c  Compute value of U' for previous solution.
c
          jl = node(1,ie)
          jr = node(2,ie)
          iul = indx(jl)
          iur = indx(jr)

          if ( iul .le. 0 ) then
            uoldx = ( fold(iur) - ul ) / he
          elseif ( iur .le. 0 ) then
            uoldx = ( ur - fold(iul) ) / he
          else
            uoldx = ( fold(iur) - fold(iul) ) / he
          end if
c
c  For basis function IL...
c
          do il = 1, nl

            ig = node(il,ie)
            iu = indx(ig)

            if ( 0 < iu ) then

              call phi ( il, xqe, phii, phiix, xleft, xrite )

              f(iu) = f(iu) + he * phii * ( ff ( xqe, problem ) )
c
c  Handle boundary conditions that prescribe the value of U'.
c
              if ( ig .eq. 0 ) then

                x = 0.0D+00
                f(iu) = f(iu) - pp ( x, problem ) * ul

              else if ( ig .eq. n ) then

                x = 1.0D+00
                f(iu) = f(iu) + pp ( x, problem ) * ur

              end if
c
c  For basis function JL...
c
              do jl = 1, nl

                jg = node(jl,ie)
                ju = indx(jg)

                call phi ( jl, xqe, phij, phijx, xleft, xrite )

                aij = he * ( pp ( xqe, problem ) * phiix * phijx
     &                  + qq ( xqe, problem ) * phii * phij
     &                  + uold * phii * phijx )

                if ( ju .le. 0 ) then

                  if ( jg .eq. 0 ) then
                    f(iu) = f(iu) - aij * ul
                  else if ( jg .eq. n ) then
                    f(iu) = f(iu) - aij * ur
                  end if

                else if ( iu .eq. ju ) then
                  adiag(iu) = adiag(iu) + aij
                else if( ju .lt. iu )then
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
      subroutine compare ( f, indx, n, nl, node, nprint, nu, problem, 
     &  ul, ur, xl, xn, xr )

c*********************************************************************72
c
cc COMPARE compares the computed and exact solutions.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 November 2006
c
c  Parameters:
c
c    Input, double precision F(NU), the solution of the linear equations.
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
c    Input, integer NL.
c    The number of basis functions used in a single
c    subinterval.  (NL-1) is the degree of the polynomials
c    used.  For this code, NL is fixed at 2, meaning that
c    piecewise linear functions are used as the basis.
c
c    Input, integer NODE(NL,N).
c    For each subinterval I:
c    NODE(1,I) is the number of the left node, and
c    NODE(2,I) is the number of the right node.
c
c    Input, integer NPRINT.
c    The number of points at which the computed solution
c    should be printed out when compared to the exact solution.
c
c    Input, integer NU.
c    NU is the number of unknowns in the linear system.
c    Depending on the value of IBC, there will be N-1,
c    N, or N+1 unknown values, which are the coefficients
c    of basis functions.
c
c    Input, integer PROBLEM, indicates which problem to be solved.
c    * 1, u = x, p=1, q=0, f=x, u(0)=0, u'(1)=1.
c    * 2, u = 2*(1-cos(0.5*pi*x))/pi, p=1, q=0,
c    f = -0.5*pi*cos(0.5*pi*x) + 2*sin(0.5*pi*x)*(1-cos(0.5*pi*x)/pi
c    u(0) = 0, u'(1)=1.
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
c    Input, double precision XL.
c    XL is the left endpoint of the interval over which the
c    differential equation is being solved.
c
c    Input, double precision XN(0:N).
c    XN(I) is the location of the I-th node.  XN(0) is XL,
c    and XN(N) is XR.
c
c    Input, double precision XR.
c    XR is the right endpoint of the interval over which the
c    differential equation is being solved.
c
      implicit none

      integer n
      integer nl
      integer nu

      double precision f(nu)
      integer i
      integer ig
      integer indx(0:n)
      integer iu
      integer j
      integer k
      integer node(nl,n)
      integer nprint
      double precision phii
      double precision phiix
      integer problem
      double precision u
      double precision u_exact
      double precision ul
      double precision ur
      double precision ux
      double precision x
      double precision xl
      double precision xleft
      double precision xn(0:n)
      double precision xr
      double precision xrite

      write(*,*)' '
      write(*,*)'Compare computed and exact solutions:'
      write(*,*)' '
      write(*,*)'      X      Computed U      Exact U'
      write(*,*)' '

      do i = 1, nprint

        x = ( dble ( nprint - i     ) * xl 
     &      + dble (          i - 1 ) * xr ) 
     &      / dble ( nprint     - 1 )

        ux = u_exact ( x, problem )

        do j = 1, n

          xleft = xn(j-1)
          xrite = xn(j)
c
c  Search for the interval that X lies in.
c
          if ( xleft .le. x .and. x .le. xrite ) then
            u = 0.0D+00
            do k = 1, nl
              ig = node(k,j)
              iu = indx(ig)
              call phi ( k, x, phii, phiix, xleft, xrite )
              if ( iu .le. 0 ) then
                if ( j .eq. 1 .and. k .eq. 1 ) then
                  u = u + ul * phii
                else if ( j .eq. n .and. k .eq. nl ) then
                  u = u + ur * phii
                end if
              else
                u = u + f(iu) * phii
              end if

            end do

            go to 25

          end if

        end do

25      continue

        write(*,'(3g14.6)') x, u, ux

      end do

      return
      end
      function ff ( x, problem )

c*********************************************************************72
c
cc FF returns the right hand side of the differential equation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 November 2006
c
c  Parameters:
c
c    Input, double precision X, the evaluation point.
c
c    Input, integer PROBLEM, indicates which problem to be solved.
c    * 1, u = x, p=1, q=0, f=x, u(0)=0, u'(1)=1.
c    * 2, u = 2*(1-cos(0.5*pi*x))/pi, p=1, q=0,
c    f = -0.5*pi*cos(0.5*pi*x) + 2*sin(0.5*pi*x)*(1-cos(0.5*pi*x)/pi
c    u(0) = 0, u'(1)=1.
c
c    Output, double precision FF, the value of F(X).
c
      implicit none

      double precision ff
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      integer problem
      double precision x
c
c  Test problem 1
c
      if ( problem == 1 ) then

        ff = x
c
c  Test problem 2
c
      else if ( problem == 2 ) then

        ff = -0.5D+00 * pi * cos ( 0.5D+00 * pi * x )
     &    + 2.0D+00 * sin ( 0.5D+00 * pi * x ) 
     &    * ( 1.0D+00 - cos ( 0.5D+00 * pi * x ) ) / pi

      end if

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
c  Parameters:
c
c    Output, double precision H(N), the length of the subintervals.  
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
c    Input, integer NL.
c    The number of basis functions used in a single
c    subinterval.  (NL-1) is the degree of the polynomials
c    used.  For this code, NL is fixed at 2, meaning that
c    piecewise linear functions are used as the basis.
c
c    Output, integer NODE(NL,N).
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
c    Depending on the value of IBC, there will be N-1,
c    N, or N+1 unknown values, which are the coefficients
c    of basis functions.
c
c    Input, double precision XL.
c    XL is the left endpoint of the interval over which the
c    differential equation is being solved.
c
c    Output, double precision XN(0:N).
c    XN(I) is the location of the I-th node.  XN(0) is XL,
c    and XN(N) is XR.
c
c    Output, double precision XQUAD(N)
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
        xn(i)= ( dble ( nsub - i ) * xl 
     &         + dble (        i ) * xr ) 
     &         / dble ( nsub     )
        write(*,'(i6,g14.6)') i, xn(i)
      end do
c
c  Set the lengths of each subinterval.
c
      write(*,*)' '
      write(*,*)'Subint    Length'
      write(*,*)' '
      do i = 1, nsub
        h(i) = xn(i) - xn(i-1)
        write(*,'(i6,g14.6)') i, h(i)
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
        write(*,'(i6,g14.6)')i,xquad(i)
      end do
c
c  Set the value of NODE, which records, for each interval,
c  the node numbers at the left and right.
c
      write(*,*)' '
      write(*,*)'Subint  Left Node  Right Node'
      write(*,*)' '
      do i = 1, nsub
        node(1,i) = i - 1
        node(2,i) = i
        write(*,'(3i6)') i, node(1,i), node(2,i)
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
      write(*,*)'  Node  Unknown'
      write(*,*)' '
      do i = 0, nsub
        write(*,'(2i6)') i, indx(i)
      end do

      return
      end
      subroutine init ( ibc, imax, nprint, nquad, problem, ul, ur, 
     &  xl, xr )

c*********************************************************************72
c
cc INIT initializes variables that define the problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 November 2006
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
c    Output, integer IMAX.
c    The number of Newton iterations to carry out.
c
c    Output, integer NPRINT.
c    The number of points at which the computed solution
c    should be printed out when compared to the exact solution.
c
c    Output, integer NQUAD.
c    The number of quadrature points used in a subinterval.
c    This code uses NQUAD = 1.
c
c    Output, integer PROBLEM, indicates which problem to be solved.
c    * 1, u = x, p=1, q=0, f=x, u(0)=0, u'(1)=1.
c    * 2, u = 2*(1-cos(0.5*pi*x))/pi, p=1, q=0,
c    f = -0.5*pi*cos(0.5*pi*x) + 2*sin(0.5*pi*x)*(1-cos(0.5*pi*x)/pi
c    u(0) = 0, u'(1)=1.
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
c    Output, double precision XR.
c    XR is the right endpoint of the interval over which the
c    differential equation is being solved.
c
      implicit none

      integer ibc
      integer imax
      integer nprint
      integer nquad
      integer problem
      double precision ul
      double precision ur
      double precision xl
      double precision xr

      ibc = 1

      imax = 10

      nprint = 9

      nquad = 1
      problem = 2
      ul = 0.0D+00
      ur = 1.0D+00

      xl = 0.0D+00
      xr = 1.0D+00
c
c  Print out the values that have been set.
c
      write(*,*)' '
      write(*,*)'The equation is to be solved for'
      write(*,*)'X greater than XL = ',xl
      write(*,*)' and less than XR = ', xr
      write(*,*)' '
      write(*,*)'The boundary conditions are:'
      write(*,*)' '

      if ( ibc .eq. 1 .or. ibc .eq. 3 ) then
        write(*,*)'  At X=XL, U=', ul
      else
        write(*,*)'  At X=XL, U''=', ul
      end if

      if(ibc.eq.2.or.ibc.eq.3)then
        write(*,*)'  At X=XR, U=', ur
      else
        write(*,*)'  At X=XR, U''=', ur
      end if

      if ( problem == 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  This is test problem #1:'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  P(X) = 1, Q(X) = 0, F(X) = X.'
        write ( *, '(a)' ) 
     &  '  Boundary conditions: U(0) = 0, U''(1) = 1.'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  The exact solution is U(X) = X'
      else if ( problem == 2 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  This is test problem #2:'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  P(X) = 1, Q(X) = 0, '
        write ( *, '(a)' ) '  F(X) = -0.5*pi*cos(0.5*pi*X)'
        write ( *, '(a)' ) 
     &  '        + 2*sin(0.5*pi*X)*(1-cos(0.5*pi*X)/pi.'
        write ( *, '(a)' ) 
     &  '  Boundary conditions: U(0) = 0, U''(1) = 1.'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    '  The exact solution is U(X) = 2*(1-cos(pi*x/2))/pi'
      end if

      write(*,*)' '
      write(*,*)'Number of quadrature points per element is ',nquad
      write(*,*)'Number of iterations is ',imax

      return
      end
      subroutine output ( f, ibc, indx, nsub, nu, ul, ur, xn )

c*********************************************************************72
c
cc OUTPUT prints out the computed solution at the nodes.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 November 2006
c
c  Parameters:
c
c    Input, double precision F(NU), the solution of the linear equations.
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
c    integer NSUB.
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
        if ( i .eq. 0 ) then
          if ( ibc .eq. 1 .or. ibc .eq. 3 ) then
            u = ul
          else
            u = f(indx(i))
          end if
        else if ( i .eq. nsub ) then
          if ( ibc .eq. 2 .or. ibc .eq. 3 ) then
            u = ur
          else
            u = f(indx(i))
          end if
        else
          u = f(indx(i))
        end if
        write ( *, '(i4,2g14.6)') i, xn(i), u
      end do

      return
      end
      subroutine phi ( il, x, phii, phiix, xleft, xrite )

c*********************************************************************72
c
cc PHI evaluates a linear basis function.
c
c  Discussion:
c
c    This routine evaluates a linear basis function and its derivative 
c    at a point X in an interval.  In any interval, there are just two 
c    basis functions.  The first basis function is a line which is 1 
c    at the left endpoint and 0 at the right.  The second basis 
c    function is 0 at the left endpoint and 1 at the right.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 November 2006
c
c  Parameters:
c
c    Input, integer IL, the index of the basis function.
c    1, the function which is 1 at XLEFT and 0 at XRITE.
c    2, the function which is 0 at XLEFT and 1 at XRITE.
c
c    Input, double precision X, the evaluation point.
c
c    Output, double precision PHII, PHIIX, the value of the
c    basis function and its derivative at X.
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
      function pp ( x, problem )

c*********************************************************************72
c
cc PP evaluates the coefficient function P(X).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 November 2006
c
c  Parameters:
c
c    Input, double precision X, the evaluation point.
c
c    Input, integer PROBLEM, indicates which problem to be solved.
c    * 1, u = x, p=1, q=0, f=x, u(0)=0, u'(1)=1.
c    * 2, u = 2*(1-cos(0.5*pi*x))/pi, p=1, q=0,
c    f = -0.5*pi*cos(0.5*pi*x) + 2*sin(0.5*pi*x)*(1-cos(0.5*pi*x)/pi
c    u(0) = 0, u'(1)=1.
c
c    Output, double precision PP, the value of P(X).
c
      implicit none

      double precision pp
      integer problem
      double precision x
c
c  Test problem 1
c
      if ( problem == 1 ) then

        pp = 1.0D+00
c
c  Test problem 2
c
      else if ( problem == 2 ) then

        pp = 1.0D+00

      end if

      return
      end
      subroutine prsys ( adiag, aleft, arite, f, nu )

c*********************************************************************72
c
cc PRSYS prints out the linear system.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 November 2006
c
c  Parameters:
c
c    Input, double precision ADIAG(NU), ALEFT(NU), ARITE(NU),
c    the diagonal, left and right entries of the equations.
c
c    Input, double precision F(NU), the right hand side of the linear
c    system to be solved.
c
c    Input, integer NU.
c    NU is the number of equations to be solved.
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
      write(*,*)'Equation  ALEFT  ADIAG  ARITE  RHS'
      write(*,*)' '
      do i = 1, nu
        if ( i .eq. 1 ) then
          write(*,'(i3,14x,3g14.6)')i,adiag(i),arite(i),f(i)
        else if ( i .lt. nu ) then
          write(*,'(i3,4g14.6)')i,aleft(i),adiag(i),arite(i),f(i)
        else
          write(*,'(i3,2g14.6,14x,g14.6)')i,aleft(i),adiag(i),f(i)
        end if
      end do

      return
      end
      function qq ( x, problem )

c*********************************************************************72
c
cc QQ evaluates the coefficient function Q(X).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 November 2006
c
c  Parameters:
c
c    Input, double precision X, the evaluation point.
c
c    Input, integer PROBLEM, indicates which problem to be solved.
c    * 1, u = x, p=1, q=0, f=x, u(0)=0, u'(1)=1.
c    * 2, u = 2*(1-cos(0.5*pi*x))/pi, p=1, q=0,
c    f = -0.5*pi*cos(0.5*pi*x) + 2*sin(0.5*pi*x)*(1-cos(0.5*pi*x)/pi
c    u(0) = 0, u'(1)=1.
c
c    Output, double precision QQ, the value of Q(X).
c
      implicit none

      integer problem
      double precision qq
      double precision x
c
c  Test problem 1
c
      if ( problem == 1 ) then

        qq = 0.0D+00
c
c  Test problem 2
c
      else if ( problem == 2 ) then

        qq = 0.0D+00

      end if

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
c  Parameters:
c
c    Input/output, double precision ADIAG(NU), ALEFT(NU), ARITE(NU).
c    On input, ADIAG, ALEFT, and ARITE contain the diagonal,
c    left and right entries of the equations.
c    On output, ADIAG and ARITE have been changed in order
c    to compute the solution.
c
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
c    Input, integer NU.
c    NU is the number of equations to be solved.
c
      implicit none

      integer nu

      double precision adiag(nu)
      double precision aleft(nu)
      double precision arite(nu-1)
      double precision f(nu)
      integer i

      arite(1) = arite(1) / adiag(1)
      do i = 2, nu-1
        adiag(i) = adiag(i) - aleft(i) * arite(i-1)
        arite(i) = arite(i) / adiag(i)
      end do
      adiag(nu) = adiag(nu) - aleft(nu) * arite(nu-1)

      f(1)=f(1) / adiag(1)
      do i = 2, nu
        f(i) = ( f(i) - aleft(i) * f(i-1) ) / adiag(i)
      end do

      do i = nu-1, 1, -1
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
      function u_exact ( x, problem )

c*********************************************************************72
c
cc U_EXACT returns the value of the exact solution at a point X.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 November 2006
c
c  Parameters:
c
c    Input, double precision X, the evaluation point.
c
c    Input, integer PROBLEM, indicates which problem to be solved.
c    * 1, u = x, p=1, q=0, f=x, u(0)=0, u'(1)=1.
c    * 2, u = 2*(1-cos(0.5*pi*x))/pi, p=1, q=0,
c    f = -0.5*pi*cos(0.5*pi*x) + 2*sin(0.5*pi*x)*(1-cos(0.5*pi*x)/pi
c    u(0) = 0, u'(1)=1.
c
c    Output, double precision U_EXACT, the value of the exact solution at X.
c
      implicit none

      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      integer problem
      double precision u_exact
      double precision x
c
c  Test problem 1
c
      if ( problem == 1 ) then

        u_exact = x
c
c  Test problem 2
c
      else if ( problem == 2 ) then

        u_exact = 2.0D+00 * ( 1.0D+00 - cos ( 0.5D+00 * pi * x ) ) / pi

      end if

      return
      end
