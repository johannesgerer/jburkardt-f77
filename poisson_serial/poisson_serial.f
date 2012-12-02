      program main

c*********************************************************************72
c
cc MAIN is the main program for POISSON_SERIAL.
c
c  Discussion:
c
c    POISSON_SERIAL is a program for solving the Poisson problem.
c
c    This program runs serially.  Its output is used as a benchmark for
c    comparison with similar programs run in a parallel environment.
c
c    The Poisson equation
c
c      - DEL^2 U(X,Y) = F(X,Y)
c
c    is solved on the unit square [0,1] x [0,1] using a grid of [0,NX+1] by
c    [0,NY+1] evenly spaced points.  
c
c    The boundary conditions and F are set so that the exact solution is
c
c      U(x,y) = sin ( pi * x * y )
c
c    so that
c
c      - DEL^2 U(x,y) = pi^2 * ( x^2 + y^2 ) * sin ( pi * x * y )
c
c    The Jacobi iteration is repeatedly applied until convergence is detected.
c
c    For convenience in writing the discretized equations, we assume that NX = NY.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 October 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nx
      integer ny

      parameter ( nx = 11 )
      parameter ( ny = 11 )

      logical converged
      double precision diff
      double precision dx
      double precision dy
      double precision error
      double precision f(nx,ny)
      integer i
      integer it
      integer it_max
      integer j
      double precision r8mat_rms
      double precision tolerance
      double precision u(nx,ny)
      double precision u_exact
      double precision u_norm
      double precision udiff(nx,ny)
      double precision uexact(nx,ny)
      double precision unew(nx,ny)
      double precision unew_norm
      double precision x
      double precision y

      it_max = 1000
      tolerance = 0.000001D+00
      dx = 1.0D+00 / dble ( nx - 1 )
      dy = dx
c
c  Print a message.
c
      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POISSON_SERIAL:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  A program for solving the Poisson equation.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  -DEL^2 U = F(X,Y)'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  on the rectangle 0 <= X <= 1, 0 <= Y <= 1.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  F(X,Y) = pi^2 * ( x^2 + y^2 ) * sin ( pi * x * y )'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) 
     &  '  The number of interior X grid points is ', nx
      write ( *, '(a,i6)' ) 
     &  '  The number of interior Y grid points is ', ny
      write ( *, '(a,g14.6)' ) '  The X grid spacing is ', dx
      write ( *, '(a,g14.6)' ) '  The Y grid spacing is ', dy
c
c  Initialize the data.
c
      call rhs ( nx, ny, f )
c
c  Set the initial solution estimate.
c  We are "allowed" to pick up the boundary conditions exactly.
c
      do j = 1, ny
        do i = 1, nx
          if ( i == 1 .or. i == nx .or.
     &         j == 1 .or. j == ny ) then
            unew(i,j) = f(i,j)
          else
            unew(i,j) = 0.0D+00
          end if
        end do
      end do

      unew_norm = r8mat_rms ( nx, ny, unew )
c
c  Set up the exact solution.
c
      do j = 1, ny 
        y = dble ( j - 1 ) / dble ( ny - 1 )
        do i = 1, nx
          x = dble ( i - 1 ) / dble ( nx - 1 )
          uexact(i,j) = u_exact ( x, y )
        end do
      end do
      u_norm = r8mat_rms ( nx, ny, uexact )
      write ( *, '(a,g14.6)' ) '  L2 norm of exact solution = ', u_norm
c
c  Do the iteration.
c
      converged = .false.

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  'Step    ||Unew||     ||Unew-U||     ||Unew-Exact||'
      write ( *, '(a)' ) ' '

      do j = 1, ny
        do i = 1, nx
          udiff(i,j) = unew(i,j) - uexact(i,j)
        end do
      end do

      error = r8mat_rms ( nx, ny, udiff )
      write ( *, '(2x,i4,2x,g14.6,2x,14x,2x,g14.6)' ) 
     &  0, unew_norm, error

      do it = 1, it_max

        do j = 1, ny
          do i = 1, nx
            u(i,j) = unew(i,j)
          end do
        end do
c
c  UNEW is derived from U by one Jacobi step.
c
        call sweep ( nx, ny, dx, dy, f, u, unew )
c
c  Check for convergence.
c
        u_norm = unew_norm
        unew_norm = r8mat_rms ( nx, ny, unew ) 

        do j = 1, ny
          do i = 1, nx
            udiff(i,j) = unew(i,j) - u(i,j)
          end do
        end do
        diff = r8mat_rms ( nx, ny, udiff )

        do j = 1, ny
          do i = 1, nx
            udiff(i,j) = unew(i,j) - uexact(i,j)
          end do
        end do
        error = r8mat_rms ( nx, ny, udiff )

        write ( *, '(2x,i4,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &    it, unew_norm, diff, error

        if ( diff .le. tolerance ) then
          converged = .true.
          exit
        end if

      end do

      if ( converged ) then
        write ( *, '(a)' ) '  The iteration has converged.'
      else
        write ( *, '(a)' ) '  The iteration has NOT converged.'
      end if
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POISSON_SERIAL:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      function r8mat_rms ( m, n, a )

c*********************************************************************72
c
cc R8MAT_RMS returns the RMS of a vector stored as an R8MAT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns in A.
c
c    Input, double precision A(M,N), the data whose RMS is desired.
c
c    Output, double precision R8MAT_RMS, the RMS of A.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer i
      integer j
      double precision r8mat_rms
      double precision value

      value = 0.0D+00
      do j = 1, n
        do i = 1, m
          value = value + a(i,j) * a(i,j)
        end do
      end do
      value = sqrt ( value / dble ( m * n ) )

      r8mat_rms = value

      return
      end
      subroutine rhs ( nx, ny, f )

c*********************************************************************72
c
cc RHS initializes the arrays.
c
c  Discussion:
c
c    It is convenient for us to set up RHS as a 2D array.  However, each
c    entry of RHS is really the right hand side of a linear system of the
c    form
c
c      A * U = F
c
c    In cases where U(I,J) is a boundary value, then the equation is simply
c
c      U(I,J) = F(i,j)
c
c    and F(I,J) holds the boundary data.
c
c    Otherwise, the equation has the form
c
c      (1/DX^2) * ( U(I+1,J)+U(I-1,J)+U(I,J-1)+U(I,J+1)-4*U(I,J) ) = F(I,J)
c
c    where DX is the spacing and F(I,J) is the value at X(I), Y(J) of
c
c      pi^2 * ( x^2 + y^2 ) * sin ( pi * x * y )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NX, NY, the X and Y grid dimensions.
c
c    Output, double precision F(0:NX+1,0:NY+1), the right hand side data.
c
      implicit none

      integer nx
      integer ny

      double precision f(nx,ny)
      double precision fnorm
      integer i
      integer j
      double precision r8mat_rms
      double precision u_exact
      double precision uxxyy_exact
      double precision x
      double precision y
c
c  The "boundary" entries of F store the boundary values of the solution.
c  The "interior" entries of F store the right hand sides of the Poisson equation.
c
      do j = 1, ny
        y = dble ( j - 1 ) / dble ( ny - 1 )
        do i = 1, nx
          x = dble ( i - 1 ) / dble ( nx - 1 )
          if ( i .eq. 1 .or. i .eq. nx .or. 
     &         j .eq. 1 .or. j .eq. ny ) then
            f(i,j) = u_exact ( x, y )
          else
            f(i,j) = - uxxyy_exact ( x, y )
          end if
        end do
      end do

      fnorm = r8mat_rms ( nx, ny, f )

      write ( *, '(a,g14.6)' ) '  L2 norm of F = ', fnorm

      return
      end
      subroutine sweep ( nx, ny, dx, dy, f, u, unew )

c*********************************************************************72
c
cc SWEEP carries out one step of the Jacobi iteration.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NX, NY, the X and Y grid dimensions.
c
c    Input, double precision DX, DY, the spacing between grid points.
c
c    Input, double precision F(NX,NY), the right hand side data.
c
c    Input, double precision U(NX,NY), the previous solution estimate.
c
c    Output, double precision UNEW(NX,NY), the updated solution estimate.
c
      implicit none

      integer nx
      integer ny

      double precision dx
      double precision dy
      double precision f(nx,ny)
      integer i
      integer j
      double precision u(nx,ny)
      double precision unew(nx,ny)
!
!  The boundary values are stored in F.
!  The interior values are set by a Jacobi iteration.
!
      do j = 1, ny
        do i = 1, nx

          if ( j .eq. 1 .or. j .eq. ny .or. 
     &         i .eq. 1 .or. i .eq. nx ) then
            unew(i,j) = f(i,j)
          else
            unew(i,j) = 0.25D+00 * (
     &        ( u(i-1,j) + u(i,j+1) + u(i,j-1) + u(i+1,j) ) 
     &        + f(i,j) * dx * dy )
          end if

        end do
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
      function u_exact ( x, y )

c*********************************************************************72
c
cc U_EXACT evaluates the exact solution.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, Y, the coordinates of a point.
c
c    Output, double precision U_EXACT, the value of the exact solution 
c    at (X,Y).
c
      implicit none

      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision u_exact
      double precision x
      double precision y

      u_exact = sin ( pi * x * y )

      return
      end
      function uxxyy_exact ( x, y )

c*********************************************************************72
c
cc UXXYY_EXACT evaluates ( d/dx d/dx + d/dy d/dy ) of the exact solution.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 October 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, Y, the coordinates of a point.
c
c    Output, double precision UXXYY_EXACT, the value of ( d/dx d/dx + d/dy d/dy ) 
c    of the exact solution at (X,Y).
c
      implicit none

      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision uxxyy_exact
      double precision x
      double precision y

      uxxyy_exact = - pi * pi * ( x * x + y * y ) * sin ( pi * x * y )

      return
      end
