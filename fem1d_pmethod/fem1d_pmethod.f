      program main

c*********************************************************************72
c
cc MAIN is the main program for FEM1D_PMETHOD.
c
c  Discussion:
c
c    FEM1D_PMETHOD is the P-method of the finite element method for 
c    one dimensional elliptic problems.
c
c    Program to solve the one dimensional problem:
c
c     - d/dX (P dU/dX) + Q U = F
c
c    by the finite-element method using a sequence of polynomials
c    which satisfy the boundary conditions and are orthogonal
c    with respect to the inner product:
c
c      (U,V) = Integral (-1 to 1) P U' V' + Q U V dx
c
c    Here U is an unknown scalar function of X defined on the
c    interval [-1,1], and P, Q and F are given functions of X.
c
c    The boundary values are U(-1)=U(1)=0.
c
c    Sample problem #1:
c
c      U=1-x**4,        P=1, Q=1, F=1.0+12.0*x**2-x**4
c
c    Sample problem #2:
c
c      U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x)
c
c    The program should be able to get the exact solution for
c    the first problem, using NP=2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 November 2006
c
c  Author:
c
c    Max Gunzburger
c    Teresa Hodge.
c
c  Parameters:
c
c    double precision A(0:NP), the squares of the norms of the 
c    basis functions.
c
c    double precision ALPHA(NP), coefficients in the recurrence
c    relationship that defines the basis functions.
c
c    double precision BETA(NP), coefficients in the recurrence
c    relationship that defines the basis functions.
c
c    double precision F(0:NP).
c    F contains the basis function coefficients that form the
c    representation of the solution U.  That is,
c      U(X) = SUM (I=0 to NP) F(I) * BASIS(I)(X)
c    where "BASIS(I)(X)" means the I-th basis function
c    evaluated at the point X.
c
c    integer NP, the highest degree polynomial to use.
c
c    integer NPRINT, the number of printout points.
c
c    integer NQUAD, the number of points in the quadrature rule.
c
c    integer PROBLEM, indicates the problem being solved.
c    1, U=1-x**4, P=1, Q=1, F=1.0+12.0*x**2-x**4.
c    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
c
c    double precision WQUAD(NQUAD), the quadrature weights.
c
c    double precision XQUAD(NQUAD), the quadrature abscissas.
c
      implicit none

      integer np
      integer nquad

      parameter ( np = 2 )
      parameter ( nquad = 10 )

      double precision a(0:np)
      double precision alpha(np)
      double precision beta(np)
      double precision f(0:np)
      integer nprint
      parameter ( nprint = 10 )
      integer problem
      double precision wquad(nquad)
      double precision xquad(nquad)

      problem = 2

      call timestamp ( )
      write(*,*)' '
      write(*,*)'FEM1D_PMETHOD'
      write(*,*)'  FORTRAN77 version'
      write(*,*)' '
      write(*,*)'Solve the two-point boundary value problem'
      write(*,*)' '
      write(*,*)'  - d/dX (P dU/dX) + Q U = F'
      write(*,*)' '
      write(*,*)'on the interval [-1,1], with'
      write(*,*)'U(-1)=U(1)=0.'
      write(*,*)' '
      write(*,*)'The P method is used, which represents U as'
      write(*,*)'a weighted sum of orthogonal polynomials.'
      write(*,*)' '
      write(*,*)' '
      write(*,*)'Highest degree polynomial to use is ', np
      write(*,*)'Number of points to be used for output=', nprint

      if ( problem == 1 ) then
        write(*,*)' '
        write(*,*)'  Problem #1:'
        write(*,*)'  U=1-x**4,'
        write(*,*)'  P=1,'
        write(*,*)'  Q=1,'
        write(*,*)'  F=1 + 12 * x**2 - x**4'
      else if ( problem == 2 ) then
        write(*,*)' '
        write(*,*)'  Problem #2:'
        write(*,*)'  U=cos(0.5*pi*x),'
        write(*,*)'  P=1,'
        write(*,*)'  Q=0,'
        write(*,*)'  F=0.25*pi*pi*cos(0.5*pi*x)'
      end if
c
c  Get quadrature abscissas and weights for interval [-1,1].
c
      call quad ( nquad, wquad, xquad )
c
c  Compute the constants for the recurrence relationship
c  that defines the basis functions.
c
      call alpbet ( a, alpha, beta, np, nquad, problem, 
     &  wquad, xquad )
c
c  Test the orthogonality of the basis functions.
c
      call ortho ( a, alpha, beta, np, nquad, problem,
     &  wquad, xquad )
c
c  Solve for the solution of the problem, in terms of coefficients
c  of the basis functions.
c
      call sol ( a, alpha, beta, f, np, nquad, problem,
     &  wquad, xquad )
c
c  Print out the solution, evaluated at each of the NPRINT points.
c
      call output ( alpha, beta, f, np, nprint )
c
c  Compare the computed and exact solutions.
c
      call compare ( alpha, beta, f, np, nprint, problem )

      write(*,*)' '
      write(*,*)'FEM1D_PMETHOD'
      write(*,*)'  Normal end of execution.'

      write(*,*)' '
      call timestamp ( )

      stop
      end
      subroutine alpbet ( a, alpha, beta, np, nquad, problem, 
     &  wquad, xquad )

c*********************************************************************72
c
cc ALPBET calculates the recurrence coefficients.
c
c  Discussion:
c
c    This routine computes ALPHA and BETA, the coefficients in the 
c    three-term recurrence relation for the orthogonal basis functions
c    on [-1,1].
c
c    It also calculates A, the square of the norm of each basis
c    function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 November 2006
c
c  Author:
c
c    Max Gunzburger
c    Teresa Hodge.
c
c  Parameters:
c
c    Output, double precision A(0:NP), the squares of the norms of the 
c    basis functions.
c
c    Output, double precision ALPHA(NP), coefficients in the recurrence
c    relationship that defines the basis functions.
c
c    Output, double precision BETA(NP), coefficients in the recurrence
c    relationship that defines the basis functions.
c
c    Input, integer NP, the highest degree polynomial to use.
c
c    Input, integer NQUAD, the number of points in the quadrature rule.
c
c    Input, integer PROBLEM, indicates the problem being solved.
c    1, U=1-x**4, P=1, Q=1, F=1.0+12.0*x**2-x**4.
c    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
c
c    Input, double precision WQUAD(NQUAD), the quadrature weights.
c
c    Input, double precision XQUAD(NQUAD), the quadrature abscissas.
c
      implicit none

      integer np
      integer nquad

      double precision a(0:np)
      double precision alpha(np)
      double precision beta(np)
      integer i
      integer iq
      integer k
      double precision pp
      external pp
      integer problem
      double precision q
      double precision qm1
      double precision qm1x
      double precision qm2
      double precision qm2x
      double precision qq
      external qq
      double precision qx
      double precision s
      double precision ss
      double precision su
      double precision sv
      double precision t
      double precision u
      double precision v
      double precision wquad(nquad)
      double precision x
      double precision xquad(nquad)

      ss = 0.0D+00
      su = 0.0D+00

      do iq = 1, nquad

        x = xquad(iq)

        s = 4.0D+00 * pp ( x, problem ) * x**2
     &    + qq ( x, problem ) * ( 1.0D+00 - x * x )**2

        u = 2.0D+00 * pp ( x, problem ) 
     &    * x * ( 3.0D+00 * x**2 - 1.0D+00 )
     &    + x * qq ( x, problem ) * ( 1.0D+00 - x * x )**2

        ss = ss + s * wquad(iq)

        su = su + u * wquad(iq)

      end do

      beta(1) = 0.0D+00
      a(0) = ss
      alpha(1) = su / ss

      do i = 2, np+1

        ss = 0.0D+00
        su = 0.0D+00
        sv = 0.0D+00

        do iq = 1, nquad

          x = xquad(iq)
          q = 1.0D+00
          qm1 = 0.0D+00
          qx = 0.0D+00
          qm1x = 0.0D+00

          do k = 1, i-1
            qm2 = qm1
            qm1 = q
            qm2x = qm1x
            qm1x = qx
            q = ( x - alpha(k) ) * qm1 - beta(k) * qm2
            qx = qm1 + ( x - alpha(k) ) * qm1x - beta(k) * qm2x
          end do

          t = 1.0D+00 - x**2

          s = pp ( x, problem ) * ( t * qx - 2.0D+00 * x * q )**2 
     &      + qq ( x, problem ) * ( t * q )**2

          u = pp ( x, problem ) 
     &      * ( x * t * qx + ( 1.0D+00 - 3.0D+00 * x**2 ) * q ) 
     &      * ( t * qx - 2.0D+00 * x * q ) 
     &      + x * qq ( x, problem ) * ( t * q )**2

          v = pp ( x, problem ) 
     &      * ( x * t * qx + ( 1.0D+00 - 3.0D+00 * x**2 ) * q ) 
     &      * ( t * qm1x - 2.0D+00 * x * qm1 )
     &      + x * qq ( x, problem ) * t**2 * q * qm1

          ss = ss + s * wquad(iq)
          su = su + u * wquad(iq)
          sv = sv + v * wquad(iq)

        end do

        a(i-1) = ss

        if ( i .le. np ) then
          alpha(i) = su / ss
          beta(i) = sv / a(i-2)
        end if

      end do

      return
      end
      subroutine compare ( alpha, beta, f, np, nprint, problem )

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
c    02 November 2006
c
c  Author:
c
c    Max Gunzburger
c    Teresa Hodge.
c
c  Parameters:
c
c    Input, double precision ALPHA(NP), coefficients in the recurrence
c    relationship that defines the basis functions.
c
c    Input, double precision BETA(NP), coefficients in the recurrence
c    relationship that defines the basis functions.
c
c    Input, double precision F(0:NP).
c    F contains the basis function coefficients that form the
c    representation of the solution U.  That is,
c      U(X) = SUM (I=0 to NP) F(I) * BASIS(I)(X)
c    where "BASIS(I)(X)" means the I-th basis function
c    evaluated at the point X.
c
c    Input, integer NP, the highest degree polynomial to use.
c
c    Input, integer NPRINT, the number of printout points.
c
c    Input, integer PROBLEM, indicates the problem being solved.
c    1, U=1-x**4, P=1, Q=1, F=1.0+12.0*x**2-x**4.
c    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
c
      implicit none

      integer np

      double precision alpha(np)
      double precision beta(np)
      double precision el2
      double precision error
      double precision f(0:np)
      integer i
      integer ip
      integer nprint
      double precision phii
      double precision phiix
      integer problem
      double precision ue
      double precision u_exact
      double precision up
      double precision x

      write(*,*)' '
      write(*,*)'Comparison of computed and exact solutions:'
      write(*,*)' '
      write(*,*)'    X        U computed    U exact     Error'
      write(*,*)' '
      el2 = 0.0D+00

      do ip = 0, nprint

        x = dble ( 2 * ip - nprint ) / dble ( nprint )
        ue = u_exact ( x, problem )
        up = 0.0D+00

        do i = 0, np
          call phi ( alpha, beta, i, np, phii, phiix, x )
          up = up + phii * f(i)
        end do

        error = ue - up
        write(*,'(4e12.4)') x, up, ue, error
        el2 = el2 + error**2

      end do

      el2 = sqrt ( el2 )
      write(*,*)' '
      write(*,*)'Little L2 error =', el2
      write(*,*)' '

      return
      end
      function ff ( x, problem )

c*********************************************************************72
c
cc FF evaluates the right hand side function F(X).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 November 2006
c
c  Author:
c
c    Max Gunzburger
c    Teresa Hodge.
c
c  Parameters:
c
c    Input, double precision X, the evaluation point.
c
c    Input, integer PROBLEM, indicates the problem being solved.
c    1, U=1-x**4, P=1, Q=1, F=1.0+12.0*x**2-x**4.
c    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
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

        ff = 1.0D+00 + 12.0D+00 * x**2 - x**4
c
c  Test problem 2
c
      else if ( problem == 2 ) then

        ff = 0.25D+00 * pi * pi * cos ( 0.5D+00 * pi * x )

      end if

      return
      end
      subroutine ortho ( a, alpha, beta, np, nquad, problem,
     &  wquad, xquad )

c*********************************************************************72
c
cc ORTHO tests the basis functions for orthogonality.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 November 2006
c
c  Author:
c
c    Max Gunzburger
c    Teresa Hodge.
c
c  Parameters:
c
c    Input, double precision A(0:NP), the squares of the norms of the 
c    basis functions.
c
c    Input, double precision ALPHA(NP), coefficients in the recurrence
c    relationship that defines the basis functions.
c
c    Input, double precision BETA(NP), coefficients in the recurrence
c    relationship that defines the basis functions.
c
c    Input, integer NP, the highest degree polynomial to use.
c
c    Input, integer NQUAD, the number of points in the quadrature rule.
c
c    Input, integer PROBLEM, indicates the problem being solved.
c    1, U=1-x**4, P=1, Q=1, F=1.0+12.0*x**2-x**4.
c    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
c
c    Input, double precision WQUAD(NQUAD), the quadrature weights.
c
c    Input, double precision XQUAD(NQUAD), the quadrature abscissas.
c
      implicit none

      integer np
      integer nquad

      double precision a(0:np)
      double precision alpha(np)
      double precision b(0:np,0:np)
      double precision beta(np)
      double precision bij
      integer i
      integer iq
      integer j
      double precision phii
      double precision phiix
      double precision phij
      double precision phijx
      double precision pp
      external pp
      integer problem
      double precision qq
      external qq
      double precision wquad(nquad)
      double precision x
      double precision xquad(nquad)
c
c  Zero out the B array, so we can start summing up the
c  dot products.
c
      do i = 0, np
        do j = 0, np
          b(i,j) = 0.0D+00
        end do
      end do
c
c  Approximate the integral of the product of basis function
c  I and basis function J over the interval [-1,1].
c
c  We expect to get zero, except when I and J are equal,
c  when we should get A(I).
c
      do iq = 1, nquad
        x = xquad(iq)
        do i = 0, np
          call phi ( alpha, beta, i, np, phii, phiix, x )
          do j = 0, np
            call phi ( alpha, beta, j, np, phij, phijx, x )
            bij = pp ( x, problem ) * phiix * phijx 
     &          + qq ( x, problem ) * phii * phij
            b(i,j) = b(i,j) + bij * wquad(iq)
          end do
        end do
      end do
c
c  Print out the results of the test.
c
      write(*,*)' '
      write(*,*)'Basis function orthogonality test:'
      write(*,*)' '
      write(*,*)'   i   j     b(i,j)/a(i)'
      write(*,*)' '
      do i = 0, np
        write(*,*)' '
        do j = 0, np
          write(*,'(2i4,g14.6)') i, j, b(i,j) / a(i)
        end do
      end do

      return
      end
      subroutine output ( alpha, beta, f, np, nprint )

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
c    02 November 2006
c
c  Author:
c
c    Max Gunzburger
c    Teresa Hodge.
c
c  Parameters:
c
c    Input, double precision ALPHA(NP), coefficients in the recurrence
c    relationship that defines the basis functions.
c
c    Input, double precision BETA(NP), coefficients in the recurrence
c    relationship that defines the basis functions.
c
c    Input, double precision F(0:NP).
c    F contains the basis function coefficients that form the
c    representation of the solution U.  That is,
c      U(X) = SUM (I=0 to NP) F(I) * BASIS(I)(X)
c    where "BASIS(I)(X)" means the I-th basis function
c    evaluated at the point X.
c
c    Input, integer NP, the highest degree polynomial to use.
c
c    Input, integer NPRINT, the number of printout points.
c
      implicit none

      integer np

      double precision alpha(np)
      double precision beta(np)
      double precision f(0:np)
      integer i
      integer ip
      integer nprint
      double precision phii
      double precision phiix
      double precision up
      double precision x

      write(*,*)' '
      write(*,*)'Representation of solution:'
      write(*,*)' '
      write(*,*)'Basis function coefficients:'
      write(*,*)' '

      do i = 0, np
        write ( *, '(i4,g14.6)' ) i, f(i)
      end do

      write(*,*)' '
      write(*,*)' '
      write(*,*)'       X     Approximate Solution'
      write(*,*)' '

      do ip = 0, nprint

        x = dble ( 2 * ip - nprint ) / dble ( nprint )
        up = 0.0D+00

        do i = 0, np
          call phi ( alpha, beta, i, np, phii, phiix, x )
          up = up + phii * f(i)
        end do

        write ( *, '(2g14.6)') x, up

      end do

      write ( *, * ) ' '

      return
      end
      subroutine phi ( alpha, beta, i, np, phii, phiix, x )

c*********************************************************************72
c
cc PHI evaluate the I-th basis function at the point X.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 November 2006
c
c  Author:
c
c    Max Gunzburger
c    Teresa Hodge.
c
c  Parameters:
c
c    Input, double precision ALPHA(NP), coefficients in the recurrence
c    relationship that defines the basis functions.
c
c    Input, double precision BETA(NP), coefficients in the recurrence
c    relationship that defines the basis functions.
c
c    Input, integer I, the index of the basis function.
c
c    Input, integer NP, the highest degree polynomial to use.
c
c    Output, double precision PHII, PHIIX, the value of the basis function
c    and its derivative, evaluated at X.
c
c    Input, double precision X, the evaluation point.
c
      implicit none

      integer np

      double precision alpha(np)
      double precision beta(np)
      integer i
      integer j
      double precision phii
      double precision phiix
      double precision q
      double precision qm1
      double precision qm1x
      double precision qm2
      double precision qm2x
      double precision qx
      double precision t
      double precision x

      qm1 = 0.0D+00
      q = 1.0D+00
      qm1x = 0.0D+00
      qx = 0.0D+00

      do j = 1, i
        qm2 = qm1
        qm1 = q
        qm2x = qm1x
        qm1x = qx
        t = x - alpha(j)
        q = t * qm1 - beta(j) * qm2
        qx = qm1 + t * qm1x - beta(j) * qm2x
      end do

      t = 1.0D+00 - x**2
      phii = t * q
      phiix = t * qx - 2.0D+00 * x * q

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
c    02 November 2006
c
c  Author:
c
c    Max Gunzburger
c    Teresa Hodge.
c
c  Parameters:
c
c    Input, double precision X, the evaluation point.
c
c    Input, integer PROBLEM, indicates the problem being solved.
c    1, U=1-x**4, P=1, Q=1, F=1.0+12.0*x**2-x**4.
c    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
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
      else

        pp = 1.0D+00

      end if

      return
      end
      subroutine quad ( nquad, wquad, xquad )

c*********************************************************************72
c
cc QUAD returns the quadrature rule.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 November 2006
c
c  Author:
c
c    Max Gunzburger
c    Teresa Hodge.
c
c  Parameters:
c
c    Input, integer NQUAD, the number of points in the quadrature rule.
c
c    Output, double precision WQUAD(NQUAD), the quadrature weights.
c
c    Output, double precision XQUAD(NQUAD), the quadrature abscissas.
c
      implicit none

      integer nquad

      double precision wquad(nquad)
      double precision xquad(nquad)
c
c  Quadrature points on [-1,1]
c
      xquad(1)  = -0.973906528517172D+00
      xquad(2)  = -0.865063366688985D+00
      xquad(3)  = -0.679409568299024D+00
      xquad(4)  = -0.433395394129247D+00
      xquad(5)  = -0.148874338981631D+00
      xquad(6)  =  0.148874338981631D+00
      xquad(7)  =  0.433395394129247D+00
      xquad(8)  =  0.679409568299024D+00
      xquad(9)  =  0.865063366688985D+00
      xquad(10) =  0.973906528517172D+00
c
c  Weight factors
c
      wquad(1)  = 0.066671344308688D+00
      wquad(2)  = 0.149451349150581D+00
      wquad(3)  = 0.219086362515982D+00
      wquad(4)  = 0.269266719309996D+00
      wquad(5)  = 0.295524224714753D+00
      wquad(6)  = 0.295524224714753D+00
      wquad(7)  = 0.269266719309996D+00
      wquad(8)  = 0.219086362515982D+00
      wquad(9)  = 0.149451349150581D+00
      wquad(10) = 0.066671344308688D+00

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
c    02 November 2006
c
c  Author:
c
c    Max Gunzburger
c    Teresa Hodge.
c
c  Parameters:
c
c    Input, double precision X, the evaluation point.
c
c    Input, integer PROBLEM, indicates the problem being solved.
c    1, U=1-x**4, P=1, Q=1, F=1.0+12.0*x**2-x**4.
c    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
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

        qq = 1.0D+00
c
c  Test problem 2
c
      else if ( problem == 2 ) then

        qq = 0.0D+00

      end if

      return
      end
      subroutine sol ( a, alpha, beta, f, np, nquad, problem, 
     &  wquad, xquad )

c*********************************************************************72
c
cc SOL solves for the coefficients.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    02 November 2006
c
c  Author:
c
c    Max Gunzburger
c    Teresa Hodge.
c
c  Parameters:
c
c    Input, double precision A(0:NP), the squares of the norms of the 
c    basis functions.
c
c    Input, double precision ALPHA(NP), coefficients in the recurrence
c    relationship that defines the basis functions.
c
c    Input, double precision BETA(NP), coefficients in the recurrence
c    relationship that defines the basis functions.
c
c    Output, double precision F(0:NP).
c    F contains the basis function coefficients that form the
c    representation of the solution U.  That is,
c      U(X) = SUM (I=0 to NP) F(I) * BASIS(I)(X)
c    where "BASIS(I)(X)" means the I-th basis function
c    evaluated at the point X.
c
c    Input, integer NP, the highest degree polynomial to use.
c
c    Input, integer NQUAD, the number of points in the quadrature rule.
c
c    Input, integer PROBLEM, indicates the problem being solved.
c    1, U=1-x**4, P=1, Q=1, F=1.0+12.0*x**2-x**4.
c    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
c
c    Input, double precision WQUAD(NQUAD), the quadrature weights.
c
c    Input, double precision XQUAD(NQUAD), the quadrature abscissas.
c
      implicit none

      integer np
      integer nquad

      double precision a(0:np)
      double precision alpha(np)
      double precision beta(np)
      double precision f(0:np)
      double precision ff
      integer i
      integer iq
      double precision phii
      double precision phiix
      integer problem
      double precision t
      double precision wquad(nquad)
      double precision x
      double precision xquad(nquad)

      do i = 0, np
        f(i) = 0.0D+00
      end do

      do iq = 1, nquad
        x = xquad(iq)
        t = ff ( x, problem ) * wquad(iq)
        do i = 0, np
          call phi ( alpha, beta, i, np, phii, phiix, x )
          f(i) = f(i) + phii * t
        end do
      end do

      do i = 0, np
        f(i) = f(i) / a(i)
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
c    02 November 2006
c
c  Author:
c
c    Max Gunzburger
c    Teresa Hodge.
c
c  Parameters:
c
c    Input, double precision X, the evaluation point.
c
c    Input, integer PROBLEM, indicates the problem being solved.
c    1, U=1-x**4, P=1, Q=1, F=1.0+12.0*x**2-x**4.
c    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
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

        u_exact = 1.0D+00 - x**4
c
c  Test problem 2
c
      else if ( problem == 2 ) then

        u_exact = cos ( 0.5D+00 * pi * x )

      end if

      return
      end
