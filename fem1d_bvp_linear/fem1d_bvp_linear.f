      subroutine compute_l2_error ( n, x, u, exact, l2_error )

c*********************************************************************72
c
cc COMPUTE_L2_ERROR estimates the L2 error norm of a finite element solution.
c
c  Discussion:
c
c    We assume the finite element method has been used, over an interval [A,B]
c    involving N nodes, with piecewise linear elements used for the basis.
c    The coefficients U(1:N) have been computed, and a formula for the
c    exact solution is known.
c
c    This function estimates the L2 norm of the error:
c
c      L2_NORM = Integral ( A <= X <= B ) ( U(X) - EXACT(X) )^2 dX
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of nodes.
c
c    Input, double precision X(N), the mesh points.
c
c    Input, double precision U(N), the finite element coefficients.
c
c    Input, function EQ = EXACT ( X ), returns the value of the exact
c    solution at the point X.
c
c    Output, double precision L2_ERROR, the estimated L2 norm of the error.
c
      implicit none

      integer n

      double precision abscissa(2)
      double precision eq
      double precision exact
      external exact
      integer i
      double precision l2_error
      integer q
      integer quad_num
      double precision u(n)
      double precision ul
      double precision ur
      double precision uq
      double precision weight(2)
      double precision wq
      double precision x(n)
      double precision xl
      double precision xq
      double precision xr

      l2_error = 0.0D+00
c
c  Quadrature definitions.
c
      quad_num = 2
      abscissa(1) = -0.577350269189625764509148780502D+00
      abscissa(2) = +0.577350269189625764509148780502D+00
      weight(1) = 1.0D+00
      weight(2) = 1.0D+00
c
c  Integrate over each interval.
c
      do i = 1, n - 1

        xl = x(i)
        xr = x(i+1)
        ul = u(i)
        ur = u(i+1)

        do q = 1, quad_num

          xq = ( ( 1.0D+00 - abscissa(q) ) * xl   
     &         + ( 1.0D+00 + abscissa(q) ) * xr ) 
     &         /   2.0D+00

          wq = weight(q) * ( xr - xl ) / 2.0D+00
c
c  Use the fact that U is a linear combination of piecewise linears.
c
          uq = ( ( xr - xq      ) * ul 
     &         + (      xq - xl ) * ur ) 
     &         / ( xr      - xl )

          eq = exact ( xq )

          l2_error = l2_error + wq * ( uq - eq )**2

        end do

      end do

      l2_error = sqrt ( l2_error )

      return
      end
      subroutine compute_seminorm_error ( n, x, u, exact_ux, 
     &  seminorm_error )

c*********************************************************************72
c
cc COMPUTE_SEMINORM_ERROR estimates the seminorm error of a finite element solution.
c
c  Discussion:
c
c    We assume the finite element method has been used, over an interval [A,B]
c    involving N nodes, with piecewise linear elements used for the basis.
c    The coefficients U(1:N) have been computed, and a formula for the
c    exact derivative is known.
c
c    This function estimates the seminorm of the error:
c
c      SEMINORM = Integral ( A <= X <= B ) ( dU(X)/dx - EXACT_UX(X) )^2 dX
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 February 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of nodes.
c
c    Input, double precision X(N), the mesh points.
c
c    Input, double precision U(N), the finite element coefficients.
c
c    Input, function EQ = EXACT_UX ( X ), returns the value of the exact
c    derivative at the point X.
c
c    Output, double precision SEMINORM_ERROR, the estimated seminorm of 
c    the error.
c
      implicit none

      integer n

      double precision abscissa(2)
      double precision exact_ux
      external exact_ux
      double precision exq
      integer i
      integer q
      integer quad_num
      double precision seminorm_error
      double precision u(n)
      double precision ul
      double precision ur
      double precision uxq
      double precision weight(2)
      double precision wq
      double precision x(n)
      double precision xl
      double precision xq
      double precision xr

      seminorm_error = 0.0D+00
c
c  Quadrature definitions.
c
      quad_num = 2
      abscissa(1) = -0.577350269189625764509148780502D+00
      abscissa(2) = +0.577350269189625764509148780502D+00
      weight(1) = 1.0D+00
      weight(2) = 1.0D+00
c
c  Integrate over each interval.
c
      do i = 1, n - 1

        xl = x(i)
        xr = x(i+1)
        ul = u(i)
        ur = u(i+1)

        do q = 1, quad_num

          xq = ( ( 1.0D+00 - abscissa(q) ) * xl   
     &         + ( 1.0D+00 + abscissa(q) ) * xr ) 
     &         /   2.0D+00

          wq = weight(q) * ( xr - xl ) / 2.0D+00
c
c  The piecewise linear derivative is a constant in the interval.
c
          uxq = ( ur - ul ) / ( xr - xl )

          exq = exact_ux ( xq )
     
          seminorm_error = seminorm_error + wq * ( uxq - exq )**2

        end do

      end do

      seminorm_error = sqrt ( seminorm_error )

      return
      end
      subroutine fem1d_bvp_linear ( n, a, c, f, x, u )

c*********************************************************************72
c
cc FEM1D_BVP_LINEAR solves a two point boundary value problem.
c
c  Discussion:
c
c    The program uses the finite element method, with piecewise linear basis
c    functions to solve a boundary value problem in one dimension.
c
c    The problem is defined on the region 0 <= x <= 1.
c
c    The following differential equation is imposed between 0 and 1:
c
c      - d/dx a(x) du/dx + c(x) * u(x) = f(x)
c
c    where a(x), c(x), and f(x) are given functions.
c
c    At the boundaries, the following conditions are applied:
c
c      u(0.0) = 0.0
c      u(1.0) = 0.0
c
c    A set of N equally spaced nodes is defined on this
c    interval, with 0 = X(1) < X(2) < ... < X(N) = 1.0.
c
c    At each node I, we associate a piecewise linear basis function V(I,X),
c    which is 0 at all nodes except node I.  This implies that V(I,X) is
c    everywhere 0 except that
c
c    for X(I-1) <= X <= X(I):
c
c      V(I,X) = ( X - X(I-1) ) / ( X(I) - X(I-1) ) 
c
c    for X(I) <= X <= X(I+1):
c
c      V(I,X) = ( X(I+1) - X ) / ( X(I+1) - X(I) )
c
c    We now assume that the solution U(X) can be written as a linear
c    sum of these basis functions:
c
c      U(X) = sum ( 1 <= J <= N ) U(J) * V(J,X)
c
c    where U(X) on the left is the function of X, but on the right,
c    is meant to indicate the coefficients of the basis functions.
c
c    To determine the coefficient U(J), we multiply the original
c    differential equation by the basis function V(J,X), and use
c    integration by parts, to arrive at the I-th finite element equation:
c
c        Integral A(X) * U'(X) * V'(I,X) + C(X) * U(X) * V(I,X) dx 
c      = Integral F(X) * V(I,X) dx
c
c    We note that the functions U(X) and U'(X) can be replaced by
c    the finite element form involving the linear sum of basis functions,
c    but we also note that the resulting integrand will only be nonzero
c    for terms where J = I - 1, I, or I + 1.
c
c    By writing this equation for basis functions I = 2 through N - 1,
c    and using the boundary conditions, we have N linear equations
c    for the N unknown coefficients U(1) through U(N), which can
c    be easily solved.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of nodes.
c
c    Input, function A ( X ), evaluates a(x);
c
c    Input, function C ( X ), evaluates c(x);
c
c    Input, function F ( X ), evaluates f(x);
c
c    Input, double precision X(N), the mesh points.
c
c    Output, double precision U(N), the finite element coefficients, which are also
c    the value of the computed solution at the mesh points.
c
      implicit none

      integer n
      integer quad_num
      parameter ( quad_num = 2 )

      double precision a
      external a
      double precision abscissa(quad_num)
      double precision al
      double precision am
      double precision ar
      double precision amat(n,n)
      double precision axq
      double precision b(n)
      double precision bm
      double precision c
      external c
      double precision cxq
      double precision f
      external f
      double precision fxq
      integer i
      integer ierror
      integer j
      integer q
      double precision weight(quad_num)
      double precision wq
      double precision u(n)
      double precision vl
      double precision vlp
      double precision vm
      double precision vmp
      double precision vr
      double precision vrp
      double precision x(n)
      double precision xl
      double precision xm
      double precision xq
      double precision xr
c
c  Quadrature definitions.
c
      abscissa(1) = -0.577350269189625764509148780502D+00
      abscissa(2) = +0.577350269189625764509148780502D+00
      weight(1) = 1.0D+00
      weight(2) = 1.0D+00
c
c  Zero out the matrix and right hand side.
c
      do j = 1, n
        do i = 1, n
          amat(1:n,1:n) = 0.0D+00
        end do
      end do

      do i = 1, n
        b(i) = 0.0D+00
      end do
c
c  Equation 1 is the left boundary condition, U(0.0) = 0.0;
c
      amat(1,1) = 1.0D+00
      b(1) = 0.0D+00
c
c  Equation I involves the basis function at node I.
c  This basis function is nonzero from X(I-1) to X(I+1).
c  Equation I looks like this:
c
c    Integral A(X) U'(X) V'(I,X) 
c           + C(X) * U(X) V(I,X) dx 
c  = Integral F(X) V(I,X) dx
c
c  Then, we realize that U(X) = sum ( 1 <= J <= N ) U(J) * V(J,X), 
c  (U(X) means the function; U(J) is the coefficient of V(J,X) ).
c
c  The only V functions that are nonzero when V(I,X) is nonzero are
c  V(I-1,X) and V(I+1,X). 
c
c  Let's use the shorthand 
c
c    VL(X) = V(I-1,X)
c    VM(X) = V(I,X)
c    VR(X) = V(I+1,X)
c
c  So our equation becomes
c
c    Integral A(X) [ VL'(X) U(I-1) + VM'(X) U(I) + VR'(X) U(I+1) ] * VM'(X)
c           + C(X) [ VL(X)  U(I-1) + VM(X)  U(I) + VR(X)  U(I+1) ] * VM(X) dx
c  = Integral F(X) VM(X) dx.
c
c  
c
c  This is actually a set of N-2 linear equations for the N coefficients U.
c
c  Now gather the multipliers of U(I-1) to get the matrix entry A(I,I-1), 
c  and so on.
c
      do i = 2, n - 1
c
c  Get the left, right and middle coordinates.
c
        xl = x(i-1)
        xm = x(i)
        xr = x(i+1)
c
c  Make temporary variables for A(I,I-1), A(I,I), A(I,I+1) and B(I).
c
        al = 0.0D+00
        am = 0.0D+00
        ar = 0.0D+00
        bm = 0.0D+00
c
c  We approximate the integrals by using a weighted sum of
c  the integrand values at quadrature points.
c
        do q = 1, quad_num
c
c  Integrate over the LEFT interval, between XL and XM, where:
c
c  VL(X) = ( XM - X       ) / ( XM - XL )
c  VM(X) = (      X  - XL ) / ( XM - XL )
c  VR(X) = 0
c
c  VL'(X) =             - 1 / ( XM - XL )
c  VM'(X) =             + 1 / ( XM - XL ) 
c  VR'(X) = 0
c
          xq = ( ( 1.0D+00 - abscissa(q) ) * xl   
     &         + ( 1.0D+00 + abscissa(q) ) * xm ) 
     &         /   2.0D+00

          wq = weight(q) * ( xm - xl ) / 2.0D+00

          vl =  ( xm - xq ) / ( xm - xl )
          vlp =  - 1.0D+00  / ( xm - xl )

          vm =  ( xq - xl ) / ( xm - xl )
          vmp =  + 1.0D+00  / ( xm - xl )

          vr =  0.0D+00
          vrp = 0.0D+00

          axq = a ( xq )
          cxq = c ( xq )
          fxq = f ( xq )

          al = al + wq * ( axq * vlp * vmp + cxq * vl * vm )
          am = am + wq * ( axq * vmp * vmp + cxq * vm * vm )
          ar = ar + wq * ( axq * vrp * vmp + cxq * vr * vm )
          bm = bm + wq * ( fxq * vm )
c
c  Integrate over the RIGHT interval, between XM and XR, where:
c
c  VL(X) = 0
c  VM(X) = ( XR - X       ) / ( XR - XM )
c  VR(X) = (      X  - XM ) / ( XR - XM )
c
c  VL'(X) = 0
c  VM'(X) =             - 1 / ( XR - XM )
c  VR'(X) =             + 1 / ( XR - XM ) 
c
          xq = ( ( 1.0D+00 - abscissa(q) ) * xm   
     &         + ( 1.0D+00 + abscissa(q) ) * xr ) 
     &         /   2.0D+00

          wq = weight(q) * ( xr - xm ) / 2.0D+00

          vl = 0.0D+00
          vlp = 0.0D+00

          vm = ( xr - xq ) / ( xr - xm )
          vmp = - 1.0D+00  / ( xr - xm )

          vr = ( xq - xm ) / ( xr - xm )
          vrp =  1.0D+00   / ( xr - xm )

          axq = a ( xq )
          cxq = c ( xq )
          fxq = f ( xq )

          al = al + wq * ( axq * vlp * vmp + cxq * vl * vm )
          am = am + wq * ( axq * vmp * vmp + cxq * vm * vm )
          ar = ar + wq * ( axq * vrp * vmp + cxq * vr * vm )
          bm = bm + wq * ( fxq * vm )

        end do

        amat(i,i-1) = al
        amat(i,i)   = am
        amat(i,i+1) = ar
        b(i)        = bm

      end do
c
c  Equation N is the right boundary condition, U(1.0) = 0.0;
c
      amat(n,n) = 1.0D+00
      b(n) = 0.0D+00
c
c  Solve the linear system.
c
      call r8mat_solve2 ( n, amat, b, u, ierror )

      return
      end
      subroutine r8mat_solve2 ( n, a, b, x, ierror )

c*********************************************************************72
c
cc R8MAT_SOLVE2 computes the solution of an N by N linear system.
c
c  Discussion:
c
c    An R8MAT is an array of R8 values.
c
c    The linear system may be represented as
c
c      A*X = B
c
c    If the linear system is singular, but consistent, then the routine will
c    still produce a solution.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 August 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of equations.
c
c    Input/output, double precision A(N,N).
c    On input, A is the coefficient matrix to be inverted.
c    On output, A has been overwritten.
c
c    Input/output, double precision B(N).
c    On input, B is the right hand side of the system.
c    On output, B has been overwritten.
c
c    Output, double precision X(N), the solution of the linear system.
c
c    Output, integer IERROR.
c    0, no error detected.
c    1, consistent singularity.
c    2, inconsistent singularity.
c
      implicit none

      integer n

      double precision a(n,n)
      double precision amax
      double precision b(n)
      integer i
      integer ierror
      integer imax
      integer ipiv(n)
      integer j
      integer k
      double precision x(n)

      ierror = 0

      do i = 1, n
        ipiv(i) = 0
      end do

      do i = 1, n
        x(i) = 0.0D+00
      end do
c
c  Process the matrix.
c
      do k = 1, n
c
c  In column K:
c    Seek the row IMAX with the properties that:
c      IMAX has not already been used as a pivot;
c      A(IMAX,K) is larger in magnitude than any other candidate.
c
        amax = 0.0D+00
        imax = 0
        do i = 1, n
          if ( ipiv(i) .eq. 0 ) then
            if ( amax .lt. abs ( a(i,k) ) ) then
              imax = i
              amax = abs ( a(i,k) )
            end if
          end if
        end do
c
c  If you found a pivot row IMAX, then,
c    eliminate the K-th entry in all rows that have not been used for pivoting.
c
        if ( imax .ne. 0 ) then

          ipiv(imax) = k
          do j = k + 1, n
            a(imax,j) = a(imax,j) / a(imax,k)
          end do
          b(imax) = b(imax) / a(imax,k)
          a(imax,k) = 1.0D+00

          do i = 1, n

            if ( ipiv(i) .eq. 0 ) then
              do j = k + 1, n
                a(i,j) = a(i,j) - a(i,k) * a(imax,j)
              end do
              b(i) = b(i) - a(i,k) * b(imax)
              a(i,k) = 0.0D+00
            end if

          end do

        end if

      end do
c
c  Now, every row with nonzero IPIV begins with a 1, and
c  all other rows are all zero.  Begin solution.
c
      do j = n, 1, -1

        imax = 0
        do k = 1, n
          if ( ipiv(k) .eq. j ) then
            imax = k
          end if
        end do

        if ( imax .eq. 0 ) then

          x(j) = 0.0D+00

          if ( b(j) .eq. 0.0D+00 ) then
            ierror = 1
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'R8MAT_SOLVE2 - Warning:'
            write ( *, '(a,i8)' ) 
     &        '  Consistent singularity, equation = ', j
          else
            ierror = 2
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'R8MAT_SOLVE2 - Error:'
            write ( *, '(a,i8)' ) 
     &        '  Inconsistent singularity, equation = ', j
          end if

        else

          x(j) = b(imax)

          do i = 1, n
            if ( i .ne. imax ) then
              b(i) = b(i) - a(i,j) * x(j)
            end if
          end do

        end if

      end do

      return
      end
      subroutine r8vec_even ( n, alo, ahi, a )

c*********************************************************************72
c
cc R8VEC_EVEN returns an R8VEC of evenly spaced values.
c
c  Discussion:
c
c    An R8VEC is a vector of R8 values.
c
c    If N is 1, then the midpoint is returned.
c
c    Otherwise, the two endpoints are returned, and N-2 evenly
c    spaced points between them.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    09 December 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of values.
c
c    Input, double precision ALO, AHI, the low and high values.
c
c    Output, double precision A(N), N evenly spaced values.
c    Normally, A(1) = ALO and A(N) = AHI.
c    However, if N = 1, then A(1) = 0.5*(ALO+AHI).
c
      implicit none

      integer n

      double precision a(n)
      double precision ahi
      double precision alo
      integer i

      if ( n .eq. 1 ) then

        a(1) = 0.5D+00 * ( alo + ahi )

      else

        do i = 1, n
          a(i) = ( dble ( n - i     ) * alo   
     &           + dble (     i - 1 ) * ahi ) 
     &           / dble ( n     - 1 )
        end do

      end if

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
