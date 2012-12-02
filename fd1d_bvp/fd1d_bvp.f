      subroutine fd1d_bvp  ( n, a, aprime, c, f, x, u )

c*********************************************************************72
c
cc FD1D_BVP solves a two point boundary value problem.
c
c  Discussion:
c
c    The program uses the finite difference method to solve a BVP
c    (boundary value problem) in one dimension.
c
c    The problem is defined on the region x(1) <= x <= x(N).
c
c    The following differential equation is imposed in the region:
c
c      - d/dx a(x) du/dx + c(x) * u(x) = f(x)
c
c    where a(x), c(x), and f(x) are given functions.  We write out
c    the equation in full as
c
c      - a(x) * u''(x) - a'(x) * u'(x) + c(x) * u(x) = f(x)
c
c    At the boundaries, the following conditions are applied:
c
c      u(X(1)) = 0.0
c      u(X(N)) = 0.0
c
c    We replace the function U(X) by a vector of N values U, associated
c    with the nodes.
c
c    The first and last values of U are determined by the boundary conditions.
c
c    At each interior node I, we write an equation to help us determine
c    U(I).  We do this by approximating the derivatives of U(X) by
c    finite differences.  Let us write XL, XM, and XR for X(I-1), X(I) and X(I+1).
c    Similarly we have UL, UM, and UR.  Other quantities to be evaluated at
c    X(I) = XM will also be labeled with an M:
c
c      - AM * ( UL - 2 UM + UR ) / DX^2 - A'M * ( UL - UR ) / ( 2 * DX ) = FM
c
c    These N-2 linear equations for the unknown coefficients complete the
c    linear system and allow us to compute the finite difference approximation
c    to the solution of the BVP.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 August 2010
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
c    Input, function APRIME ( X ), evaluates a'(x);
c
c    Input, function C ( X ), evaluates c(x);
c
c    Input, function F ( X ), evaluates f(x);
c
c    Input, double precision X(N), the mesh points, which may be
c    nonuniformly spaced.
c
c    Output, double precision U(N), the value of the finite difference
c    approximation to the solution.
c
      implicit none

      integer n

      double precision a
      external a
      double precision am
      double precision apm
      double precision aprime
      external aprime
      double precision c
      external c
      double precision cm
      double precision f
      external f
      double precision fm
      integer i
      double precision rhs(n)
      double precision tri(3,n)
      double precision u(n)
      double precision x(n)
      double precision x1
      double precision x2
      double precision xm
c
c  Equation 1 is the left boundary condition, U(X(1)) = 0.0;
c
      tri(1,1) = 0.0D+00
      tri(2,1) = 1.0D+00
      tri(3,1) = 0.0D+00
      rhs(1) = 0.0D+00
c
c  Now gather the multipliers of U(I-1) to get the matrix entry A(I,I-1),
c  and so on.
c
      do i = 2, n - 1

        xm  = x(i)
        am  = a ( xm )
        apm = aprime ( xm )
        cm  = c ( xm )
        fm  = f ( xm )

        tri(1,i) =
     &    - 2.0D+00 * am / ( x(i) - x(i-1) ) / ( x(i+1) - x(i-1) )
     &    + apm / ( x(i+1) - x(i-1) )

        tri(2,i) =
     &    + 2.0D+00 * am / ( x(i) - x(i-1) ) / ( x(i+1) - x(i) )
     &    + cm

        tri(3,i) =
     &    - 2.0D+00 * am / ( x(i+1) - x(i) ) / ( x(i+1) - x(i-1) )
     &    - apm / ( x(i+1) - x(i-1) )

        rhs(i)   = fm

      end do
c
c  Equation N is the right boundary condition, U(X(N)) = 0.0;
c
      tri(1,n) = 0.0D+00
      tri(2,n) = 1.0D+00
      tri(3,n) = 0.0D+00
      rhs(n) = 0.0D+00
c
c  Solve the linear system.
c
      call r83np_fs ( n, tri, rhs, u )

      return
      end
      subroutine get_unit ( unit )

c*********************************************************************72
c
cc GET_UNIT returns a free FORTRAN unit number.
c
c  Discussion:
c
c    A "free" FORTRAN unit number is a value between 1 and 99 which
c    is not currently associated with an I/O device.  A free FORTRAN unit
c    number is needed in order to open a file with the OPEN command.
c
c    If UNIT = 0, then no free FORTRAN unit could be found, although
c    all 99 units were checked (except for units 5, 6 and 9, which
c    are commonly reserved for console I/O).
c
c    Otherwise, UNIT is a value between 1 and 99, representing a
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
      subroutine r83np_fs ( n, a, b, x )

c*********************************************************************72
c
cc R83NP_FS factors and solves an R83NP system.
c
c  Discussion:
c
c    The R83NP storage format is used for a tridiagonal matrix.
c    The subdiagonal   is in entries (1,2:N),
c    the diagonal      is in entries (2,1:N),
c    the superdiagonal is in entries (3,1:N-1).
c
c    This algorithm requires that each diagonal entry be nonzero.
c    It does not use pivoting, and so can fail on systems that
c    are actually nonsingular.
c
c    The "R83NP" format used for this routine is different from the R83 format.
c    Here, we insist that the nonzero entries
c    for a given row now appear in the corresponding column of the
c    packed array.
c
c  Example:
c
c    Here is how an R83NP matrix of order 5 would be stored:
c
c       *  A21 A32 A43 A54
c      A11 A22 A33 A44 A55
c      A12 A23 A34 A45  *
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 May 2009
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
c    Input, double precision B(N), the right hand side of the linear system.
c
c    Output, double precision X(N), the solution of the linear system.
c
      implicit none

      integer n

      double precision a(3,n)
      double precision b(n)
      integer i
      double precision x(n)
      double precision xmult
c
c  The diagonal entries can't be zero.
c
      do i = 1, n
        if ( a(2,i) .eq. 0.0D+00 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R83NP_FS - Fatal error!'
          write ( *, '(a,i8,a)' ) '  A(2,', i, ') = 0.'
          return
        end if
      end do

      do i = 1, n
        x(i) = b(i)
      end do

      do i = 2, n
        a(2,i) = a(2,i) - a(3,i-1) * a(1,i) / a(2,i-1)
        x(i)   = x(i)   - x(i-1)   * a(1,i) / a(2,i-1)
      end do

      x(n) = x(n) / a(2,n)
      do i = n-1, 1, -1
        x(i) = ( x(i) - a(3,i) * x(i+1) ) / a(2,i)
      end do

      return
      end
      subroutine r8mat_write ( output_filename, m, n, table )

c*********************************************************************72
c
cc R8MAT_WRITE writes a R8MAT file.
c
c  Discussion:
c
c    An R8MAT is an array of R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 October 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) OUTPUT_FILENAME, the output file name.
c
c    Input, integer M, the spatial dimension.
c
c    Input, integer N, the number of points.
c
c    Input, double precision TABLE(M,N), the data.
c
      implicit none

      integer m
      integer n

      integer j
      character * ( * ) output_filename
      integer output_unit
      character * ( 30 ) string
      double precision table(m,n)
c
c  Open the file.
c
      call get_unit ( output_unit )

      open ( unit = output_unit, file = output_filename,
     &  status = 'replace' )
c
c  Create the format string.
c
      if ( 0 .lt. m .and. 0 .lt. n ) then

        write ( string, '(a1,i8,a1,i8,a1,i8,a1)' )
     &    '(', m, 'g', 24, '.', 16, ')'
c
c  Write the data.
c
        do j = 1, n
          write ( output_unit, string ) table(1:m,j)
        end do

      end if
c
c  Close the file.
c
      close ( unit = output_unit )

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
