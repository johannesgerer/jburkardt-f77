      program main

c*********************************************************************72
c
cc MAIN is the main program for POISSON_OPENMP.
c
c  Discussion:
c
c    POISSON_OPENMP is a program for solving the Poisson problem.
c
c    This program runs in parallel, using OpenMP in the SWEEP function.
c
c    The Poisson equation
c
c      - DEL^2 U(x,y) = F(x,y)
c
c    is solved on the unit square [0,1] x [0,1] using a grid of NX by
c    NX evenly spaced points.  The first and last points in each direction
c    are boundary points.
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
c    13 December 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      include 'omp_lib.h'

      integer nx
      parameter ( nx = 161 )
      integer ny 
      parameter ( ny = 161 )

      logical converged
      double precision diff
      double precision dx
      double precision dy
      double precision error
      double precision f(nx,ny)
      integer i
      integer id
      integer itnew
      integer itold
      integer jt
      integer jt_max 
      parameter ( jt_max = 20 )
      integer j
      double precision r8mat_rms
      double precision tolerance
      parameter ( tolerance = 0.000001D+00 )
      double precision u(nx,ny)
      double precision u_exact
      double precision u_norm
      double precision udiff(nx,ny)
      double precision uexact(nx,ny)
      double precision unew(nx,ny)
      double precision unew_norm
      double precision wtime
      double precision x
      double precision y

      dx = 1.0D+00 / dble ( nx - 1 )
      dy = dx

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POISSON_OPENMP:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  A program for solving the Poisson equation.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Use OpenMP for parallel execution.'
      write ( *, '(a,i4)' ) '  The number of processors available is ',
     &  omp_get_num_procs ( )
c$omp parallel
      id = omp_get_thread_num ( )
      if ( id .eq. 0 ) then
        write ( *, '(a,i4)' ) '  The number of threads available is ',
     &  omp_get_num_threads ( )
      end if
c$omp end parallel
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  -DEL^2 U = F(X,Y)'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  on the rectangle 0 <= X <= 1, 0 <= Y <= 1.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  F(X,Y) = pi^2 * ( x^2 + y^2 ) * sin ( pi * x * y )'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  The number of interior X grid points is ', nx
      write ( *, '(a,i8)' ) 
     &  '  The number of interior Y grid points is ', ny
      write ( *, '(a,f10.4)' ) '  The X grid spacing is ', dx
      write ( *, '(a,f10.4)' ) '  The Y grid spacing is ', dy
c
c  Set the right hand side array F.
c
      call rhs ( nx, ny, f )
c
c  Set the initial solution estimate UNEW.
c  We are "allowed" to pick up the boundary conditions exactly.
c
      do j = 1, ny
        do i = 1, nx
          unew(i,j) = 0.0D+00
        end do
      end do

      do j = 1, ny
        unew(1,j) = f(1,j) 
      end do 

      do j = 1, ny
        unew(nx,j) = f(nx,j)
      end do

      do i = 1, nx
        unew(i,1) = f(i,1)
      end do

      do i = 1, nx
        unew(i,ny) = f(i,ny)
      end do

      unew_norm = r8mat_rms ( nx, ny, unew )
c
c  Set up the exact solution UEXACT.
c
      do j = 1, ny 
        y = dble ( j - 1 ) / dble ( ny - 1 )
        do i = 1, nx
          x = dble ( i - 1 ) / dble ( nx - 1 )
          uexact(i,j) = u_exact ( x, y )
        end do
      end do
      u_norm = r8mat_rms ( nx, ny, uexact )
      write ( *, '(a,g14.6)' ) '  RMS of exact solution = ', u_norm
c
c  Do the iteration.
c
      converged = .false.

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Step    ||Unew||     ||Unew-U||     ||Unew-Exact||'
      write ( *, '(a)' ) ' '

      do j = 1, ny
        do i = 1, nx
          udiff(i,j) = unew(i,j) - uexact(i,j)
        end do
      end do

      error = r8mat_rms ( nx, ny, udiff )

      write ( *, '(2x,i6,2x,g14.6,2x,14x,2x,g14.6)' ) 
     &  0, unew_norm, error

      wtime = omp_get_wtime ( )

      itnew = 0

      do

        itold = itnew
        itnew = itold + 500
c
c  SWEEP carries out 500 Jacobi steps in parallel before we come
c  back to check for convergence.
c
        call sweep ( nx, ny, dx, dy, f, itold, itnew, u, unew )
c
c  We declare convergence if successive iterates change very little.
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

        write ( *, '(2x,i6,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &    itnew, unew_norm, diff, error

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

      wtime = omp_get_wtime ( ) - wtime
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Elapsed seconds = ', wtime
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POISSON_OPENMP:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      function r8mat_rms ( m, n, a )

c*********************************************************************72
c
cc R8MAT_RMS returns the root mean square of data stored as an R8MAT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 December 2011
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
c    Output, double precision R8MAT_RMS, the root mean square of A.
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
          value = value + a(i,j)**2
        end do
      end do
      value = sqrt ( value / dble ( m * n ) ) 

      r8mat_rms = value

      return
      end
      subroutine rhs ( nx, ny, f )

c*********************************************************************72
c
cc RHS initializes the right hand side "vector".
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
c    13 December 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NX, NY, the X and Y grid dimensions.
c
c    Output, double precision F(NX,NY), the right hand side data.
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
c  The "boundary" entries of F will store the boundary values of the solution.
c
      x = 0.0D+00
      do j = 1, ny
        y = dble ( j - 1 ) / dble ( ny - 1 )
        f(1,j) = u_exact ( x, y )
      end do

      x = 1.0D+00
      do j = 1, ny
        y = dble ( j - 1 ) / dble ( ny - 1 )
        f(nx,j) = u_exact ( x, y )
      end do

      y = 0.0D+00
      do i = 1, nx
        x = dble ( i - 1 ) / dble ( nx - 1 )
        f(i,1) = u_exact ( x, y )
      end do

      y = 1.0D+00
      do i = 1, nx
        x = dble ( i - 1 ) / dble ( nx - 1 )
        f(i,ny) = u_exact ( x, y )
      end do
c
c  The "interior" entries of F store the right hand sides of the Poisson equation.
c
      do j = 2, ny - 1
        y = dble ( j - 1 ) / dble ( ny - 1 )
        do i = 2, nx - 1
          x = dble ( i - 1 ) / dble ( nx - 1 )
          f(i,j) = - uxxyy_exact ( x, y )
        end do
      end do

      fnorm = r8mat_rms ( nx, ny, f ) 

      write ( *, '(a,g14.6)' ) '  RMS of F = ', fnorm

      return
      end
      subroutine sweep ( nx, ny, dx, dy, f, itold, itnew, u, unew )

c*********************************************************************72
c
cc SWEEP carries out several steps of the Jacobi iteration.
c
c  Discussion:
c
c    Assuming DX = DY, we can approximate
c
c      - ( d/dx d/dx + d/dy d/dy ) U(X,Y) 
c
c    by
c
c      ( U(i-1,j) + U(i+1,j) + U(i,j-1) + U(i,j+1) - 4*U(i,j) ) / dx / dy
c
c    The discretization employed below will not be correct in the general
c    case where DX and DY are not equal.  It's only a little more complicated
c    to allow DX and DY to be different, but we're not going to worry about 
c    that right now.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    13 December 2011
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
c    Input, integer ITOLD, the iteration index on input.
c
c    Input, integer ITNEW, the desired iteration index
c    on output.
c
c    Output, double precision U(NX,NY), the solution estimate on 
c    iteration ITNEW-1.
c
c    Input/output, double precision UNEW(NX,NY); on input, the solution 
c    estimate on iteration ITOLD.  On output, the solution estimate on 
c    iteration ITNEW.
c
      implicit none

      integer nx
      integer ny

      double precision dx
      double precision dy
      double precision f(nx,ny)
      integer i
      integer it
      integer itnew
      integer itold
      integer j
      double precision u(nx,ny)
      double precision unew(nx,ny)
c
c$omp parallel shared ( dx, dy, f, itnew, itold, nx, ny, u, unew ) 
c$omp&         private ( i, it, j )
c
c  Carry out iterations ITOLD+1 through ITNEW.
c
      do it = itold + 1, itnew
c
c  Save the current estimate.
c
c$omp do
        do j = 1, ny
          do i = 1, nx
            u(i,j) = unew(i,j)
          end do
        end do
c$omp end do
c
c  Compute a new estimate.
c
c$omp do
        do j = 1, ny
          do i = 1, nx

            if ( i .eq. 1 .or.
     &           i .eq. nx .or.
     &           j .eq. 1 .or. 
     &           j .eq. ny ) then
              unew(i,j) = f(i,j)
            else
              unew(i,j) = 0.25D+00 * ( 
     &          u(i-1,j) + u(i,j+1) + u(i,j-1) + u(i+1,j) 
     &          + f(i,j) * dx * dy )
            end if

          end do
        end do
c$omp end do

      end do

c$omp end parallel

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
c    13 December 2011
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
c    13 December 2011
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, Y, the coordinates of a point.
c
c    Output, double precision UXXYY_EXACT, the value of 
c    ( d/dx d/dx + d/dy d/dy ) of the exact solution at (X,Y).
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
