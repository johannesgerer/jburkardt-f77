      subroutine fd1d_heat_steady ( n, a, b, ua, ub, k, f, x, u )

c*********************************************************************72
c
cc FD1D_HEAT_STEADY solves the steady 1D heat equation.
c
c  Discussion:
c
c    This program seeks a solution of the steady heat equation:
c
c      - d/dx ( K(X) dUdx ) = F(X)
c
c    over the interval [A,B] with boundary conditions
c
c      U(A) = UA,
c      U(B) = UB.
c
c    The code uses the finite difference method to approximate the
c    second derivative in space.  This results in a sparse linear system.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 May 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of grid points.
c
c    Input, double precision A, B, the interval endpoints.
c
c    Input, double precision UA, UB, the values prescribed for U
c    at the endpoints.
c
c    Input, function K(X), evaluates the thermal conductance at the N
c    points X.  Set K(X) = 1 if you don't care about this coefficient.
c
c    Input, function F(X), evaluates the heat source term at the N
c    points X.  Set F(X) = 0 if you don't want any heat sources.
c
c    Output, double precision X(N), the grid points.
c
c    Output, double precision U(N), the approximation to the solution
c    at the grid points.
c
      implicit none

      integer n

      double precision a
      double precision b
      double precision dx
      double precision f
      external f
      integer i
      double precision k
      external k
      double precision rhs(n)
      double precision tri(3,n)
      double precision u(n)
      double precision ua
      double precision ub
      double precision x(n)
      double precision xm
      double precision xp

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FD1D_HEAT_STEADY'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Finite difference solution of'
      write ( *, '(a)' ) '  the steady 1D heat equation'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    - d/dx ( k(x) dUdx ) = F(x)'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  for space interval A <= X <= B with boundary conditions'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    U(A) = UA'
      write ( *, '(a)' ) '    U(B) = UB'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  A second order difference approximation is used.'
c
c  Set the X values.
c
      dx = ( b - a ) / dble ( n - 1 )

      do i = 1, n
        x(i) = ( dble ( n - i     ) * a
     &         + dble (     i - 1 ) * b )
     &         / dble ( n -     1 )
      end do
c
c  Set up the tridiagonal matrix.
c
      tri(1,1) = 0.0D+00
      tri(2,1) = 1.0D+00
      tri(3,1) = 0.0D+00
      rhs(1) = ua

      do i = 2, n - 1

        xm = ( x(i-1) + x(i) ) / 2.0D+00
        xp = ( x(i) + x(i+1) ) / 2.0D+00

        tri(1,i) = - k ( xm )              / dx / dx
        tri(2,i) = ( k ( xm ) + k ( xp ) ) / dx / dx
        tri(3,i) =            - k ( xp )   / dx / dx

        rhs(i) = f ( x(i) )

      end do

      tri(1,n) = 0.0D+00
      tri(2,n) = 1.0D+00
      tri(3,n) = 0.0D+00
      rhs(n) = ub
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
      subroutine s_blank_delete ( s )

c*********************************************************************72
c
cc S_BLANK_DELETE removes blanks from a string, left justifying the remainder.
c
c  Discussion:
c
c    All TAB characters are also removed.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 July 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character*(*) S, the string to be transformed.
c
      implicit none

      character ch
      integer get
      integer put
      character*(*) s
      integer s_len_trim
      integer s_length
      character tab

      tab = char ( 9 )

      put = 0
      s_length = s_len_trim ( s )

      do get = 1, s_length

        ch = s(get:get)

        if ( ch .ne. ' ' .and. ch .ne. tab ) then
          put = put + 1
          s(put:put) = ch
        end if

      end do

      s(put+1:s_length) = ' '

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
