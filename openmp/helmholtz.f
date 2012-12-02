      program main

c*********************************************************************72
c
cc MAIN is the main program for HELMHOLTZ.
c
c  Discussion:
c
c    HELMHOLTZ solves a discretized Helmholtz equation.
c
c    The two dimensional region given is:
c
c      -1 <= X <= +1
c      -1 <= Y <= +1
c
c    The region is discretized by a set of M by N nodes:
c
c      P(I,J) = ( X(I), Y(J) )
c
c      X(I) = ( - M + 2*I - 1 ) / ( M - 1 )
c      Y(J) = ( - N + 2*J - 1 ) / ( N - 1 )
c
c    The Helmholtz equation for the scalar function U(X,Y) is
c
c      - Uxx(X,Y) -Uyy(X,Y) + ALPHA * U(X,Y) = F(X,Y)
c
c    where ALPHA is a positive constant.  We suppose that Dirichlet
c    boundary conditions are specified, that is, that the value of
c    U(X,Y) is given for all points along the boundary.
c
c    We suppose that the right hand side function F(X,Y) is specified in 
c    such a way that the exact solution is
c
c      U(X,Y) = ( 1 - X^2 ) * ( 1 - Y^2 )
c
c    Using standard finite difference techniques, the second derivatives
c    of U can be approximated by linear combinations of the values
c    of U at neighboring points.  Using this fact, the discretized
c    differential equation becomes a set of linear equations of the form:
c
c      A * U = F
c
c    These linear equations are then solved using a form of the Jacobi 
c    iterative method with a relaxation factor.
c
c    Directives are used in this code to achieve parallelism.
c    All do loops are parallized with default 'static' scheduling.
c
c  Modified:
c
c    18 February 2008
c
c  Author:
c
c    Joseph Robicheaux, Sanjiv Shah.
c    Modifications by John Burkardt
c
      implicit none

      include 'omp_lib.h'

      double precision alpha
      integer id
      integer it_max
      integer m
      integer n
      double precision omega
      double precision tol
      double precision wtime

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HELMHOLTZ'
      write ( *, '(a)' ) '  FORTRAN77/OpenMP version'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  A program for the 2D Helmholtz equation.'

      write ( *, '(a,i8)' ) 
     &  '  The number of processors = ', omp_get_num_procs ( )
      write ( *, '(a,i8)' ) 
     &  '  The number of threads    = ', omp_get_max_threads ( )
c
c  Set the parameters.
c
      m = 500
      n = 500
      alpha = 0.25D+00
      omega = 1.1D+00
      tol = 1.0D-08
      it_max = 100

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The region is [-1,1] x [-1,1].'
      write ( *, '(a,i8)' ) 
     &  '  The number of nodes in the X direction is M = ', m
      write ( *, '(a,i8)' ) 
     &  '  The number of nodes in the Y direction is N = ', n
      write ( *, '(a,g14.6)' ) 
     &  '  The scalar coefficient is ALPHA = ', alpha
      write ( *, '(a,g14.6)' ) 
     &  '  The relaxation value is OMEGA = ', omega
      write ( *, '(a,g14.6)' ) 
     &  '  The error tolerance is TOL = ', tol
      write ( *, '(a,i8)' ) 
     &  '  The maximum number of Jacobi iterations is IT_MAX = ', it_max
c
c  Call the driver routine.
c
      wtime = omp_get_wtime ( )

      call driver ( m, n, it_max, alpha, omega, tol )

      wtime = omp_get_wtime ( ) - wtime

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Elapsed wall clock time = ', wtime
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HELMHOLTZ'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
      subroutine driver ( m, n, it_max, alpha, omega, tol )

c*********************************************************************72
c
cc DRIVER allocates arrays and solves the problem.
c
c  Modified:
c
c    13 November 2007
c
c  Author:
c
c    Joseph Robicheaux, Sanjiv Shah.
c
c  Parameters:
c
c    Input, integer M, N, the number of grid points in the X and Y directions.
c
c    Input, integer IT_MAX, the maximum number of Jacobi iterations allowed.
c
c    Input, double precision ALPHA, the scalar coefficient in the
c    Helmholtz equation.
c
c    Input, double precision OMEGA, the relaxation parameter, which
c    should be strictly between 0 and 2.  For a pure Jacobi method,
c    use OMEGA = 1.
c
c    Input, double precision TOL, an error tolerance for the linear
c    equation solver.
c
      implicit none

      integer m
      integer n

      double precision alpha
      double precision f(m,n)
      integer it_max
      double precision omega
      double precision tol
      double precision u(m,n)
c
c  Initialize the data.
c
      call initialize ( m, n, alpha, u, f )
c
c  Solve the Helmholtz equation.
c
      call jacobi ( m, n, alpha, omega, u, f, tol, it_max )
c
c  Determine the error.
c
      call error_check ( m, n, alpha, u, f )

      return
      end
      subroutine initialize ( m, n, alpha, u, f )

c*********************************************************************72
c
cc INITIALIZE initializes the data.
c
c  Discussion:
c
c    The routine assumes that the exact solution and its second
c    derivatives are given by the routine EXACT.
c
c    The appropriate Dirichlet boundary conditions are determined
c    by getting the value of U returned by EXACT.
c
c    The appropriate right hand side function is determined by
c    having EXACT return the values of U, UXX and UYY, and setting
c
c      F = -UXX - UYY + ALPHA * U
c
c  Modified:
c
c    13 November 2007
c
c  Author:
c
c    Joseph Robicheaux, Sanjiv Shah.
c
c  Parameters:
c
c    Input, integer M, N, the number of grid points in the X and Y directions.
c
c    Input, double precision ALPHA, the scalar coefficient in the
c    Helmholtz equation.  ALPHA should be positive.
c
c    Output, double precision U(M,N), (the initial estimate of) the solution
c    of the Helmholtz equation at the grid points.
c
c    Output, double precision F(M,N), values of the right hand side function 
c    for the Helmholtz equation at the grid points.
c
      implicit none

      integer m
      integer n

      double precision alpha
      double precision dx
      double precision dy
      double precision f(m,n)
      double precision f_norm
      integer i
      integer j
      double precision u(m,n)
      double precision u_exact
      double precision uxx
      double precision uyy
      double precision x
      double precision y

c$omp parallel
c$omp& shared ( f, m, n, u )
c$omp& private ( i, j )

c$omp do
      do j = 1, n
        do i = 1, m
          u(i,j) = 0.0D+00
          f(i,j) = 0.0D+00
        end do
      end do
c$omp end do

c$omp end parallel

      dx = 2.0D+00 / dble ( m - 1 )
      dy = 2.0D+00 / dble ( n - 1 )
c
c  Set the boundary conditions.
c
      j = 1
      y = -1.0D+00 + dy * dble ( j - 1 )

      do i = 1, m
        x = -1.0D+00 + dx * dble ( i - 1 )
        call exact ( x, y, u_exact, uxx, uyy )
        f(i,j) = u_exact
        u(i,j) = u_exact
      end do

      j = n
      y = -1.0D+00 + dy * dble ( j - 1 )
      do i = 1, m
        x = -1.0D+00 + dx * dble ( i - 1 )
        call exact ( x, y, u_exact, uxx, uyy )
        f(i,j) = u_exact
        u(i,j) = u_exact
      end do

      i = 1
      x = -1.0D+00 + dx * dble ( i - 1 )
      do j = 1, n
        y = -1.0D+00 + dy * dble ( j - 1 )
        call exact ( x, y, u_exact, uxx, uyy )
        f(i,j) = u_exact
        u(i,j) = u_exact
      end do

      i = m
      x = -1.0D+00 + dx * dble ( i - 1 )
      do j = 1, n
        y = -1.0D+00 + dy * dble ( j - 1 )
        call exact ( x, y, u_exact, uxx, uyy )
        f(i,j) = u_exact
        u(i,j) = u_exact
      end do
c
c  Set the right hand side F.
c
      f_norm = 0.0D+00

c$omp parallel
c$omp& shared ( alpha, dx, dy, f, m, n ) 
c$omp& private ( i, j, u_exact, uxx, uyy, x, y )

c$omp do reduction ( + : f_norm )
      do j = 2, n - 1
        do i = 2, m - 1
          x = -1.0D+00 + dx * dble ( i - 1 )
          y = -1.0D+00 + dy * dble ( j - 1 )
          call exact ( x, y, u_exact, uxx, uyy )
          f(i,j) = - uxx - uyy + alpha * u_exact
          f_norm = f_norm + f(i,j)**2
        end do
      end do
c$omp end do

c$omp end parallel

      f_norm = sqrt ( f_norm )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Right hand side RMS = ', f_norm
      write ( *, '(a)' ) ' '

      return
      end
      subroutine jacobi ( m, n, alpha, omega, u, f, tol, it_max )

c*********************************************************************72
c
cc JACOBI applies the Jacobi iterative method to solve the linear system.
c
c  Modified:
c
c    13 November 2007
c
c  Author:
c
c    Joseph Robicheaux, Sanjiv Shah.
c
c  Parameters:
c
c    Input, integer M, N, the number of grid points in the X and Y directions.
c
c    Input, double precision ALPHA, the scalar coefficient in the
c    Helmholtz equation.  ALPHA should be positive.
c
c    Input, double precision OMEGA, the relaxation parameter, which
c    should be strictly between 0 and 2.  For a pure Jacobi method,
c    use OMEGA = 1.
c
c    Input/output, double precision U(M,N), the solution of the Helmholtz
c    equation at the grid points.
c
c    Input, double precision F(M,N), values of the right hand side function 
c    for the Helmholtz equation at the grid points.
c
c    Input, double precision TOL, an error tolerance for the linear
c    equation solver.
c
c    Input, integer IT_MAX, the maximum number of Jacobi iterations allowed.
c
      implicit none

      integer m
      integer n

      double precision alpha
      double precision ax
      double precision ay
      double precision b
      double precision dx
      double precision dy
      double precision error
      double precision f(m,n)
      integer i
      integer it
      integer it_max
      integer j
      double precision omega
      double precision resid
      double precision tol
      double precision u(m,n)
      double precision uold(m,n)
c
c  Initialize the coefficients.
c
      dx = 2.0D+00 / dble ( m - 1 )
      dy = 2.0D+00 / dble ( n - 1 )

      ax = -1.0D+00 / dx**2
      ay = -1.0D+00 / dy**2
      b  = +2.0D+00 / dx**2 + 2.0D+00 / dy**2 + alpha

      do it = 1, it_max

        error = 0.0D+00
c
c  Copy new solution into old.
c

c$omp parallel 
c$omp& shared ( m, n, u, uold )
c$omp& private ( i, j )

c$omp do
        do j = 1, n
          do i = 1, m
            uold(i,j) = u(i,j)
          end do
        end do
c$omp end do

c$omp end parallel
c
c  Compute stencil, residual, and update.
c

c$omp parallel
c$omp& shared ( ax, ay, b, f, m, n, omega, u, uold )
c$omp& private ( i, j, resid )

c$omp do reduction ( + : error )
        do j = 1, n
          do i = 1, m
c
c  Evaluate the residual.
c
            if ( i .eq. 1 .or. 
     &           i .eq. m .or. 
     &           j .eq. 1 .or. 
     &           j .eq. n ) then

              resid = uold(i,j) - f(i,j)

            else

              resid = ( ax * ( uold(i-1,j) + uold(i+1,j) ) 
     &                + ay * ( uold(i,j-1) + uold(i,j+1) ) 
     &                + b  *   uold(i,j) - f(i,j) ) / b

            end if
c
c  Update the solution.
c
            u(i,j) = uold(i,j) - omega * resid
c
c  Accumulate the residual error.
c
            error = error + resid * resid

          end do
        end do
c$omp end do

c$omp end parallel
c
c  Error check.
c
        error = sqrt ( error )

        write ( *, '(a,g14.6)' ) 
     &  '  Jacobi residual RMS            ', error

        if ( error .le. tol * dble ( m * n ) .and. 1 < it ) then
          write ( *, *) '  Test = ', tol * dble ( m * n )
          exit
        end if

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 'Total number of iterations ', it

      return
      end
      subroutine error_check ( m, n, alpha, u, f )

c*********************************************************************72
c
cc ERROR_CHECK determines the error in the numerical solution.
c
c  Modified:
c
c    13 November 2007
c
c  Author:
c
c    Joseph Robicheaux, Sanjiv Shah.
c
c  Parameters:
c
c    Input, integer M, N, the number of grid points in the X and Y directions.
c
c    Input, double precision ALPHA, the scalar coefficient in the
c    Helmholtz equation.  ALPHA should be positive.
c
c    Input, double precision U(M,N), the solution of the Helmholtz equation 
c    at the grid points.
c
c    Input, double precision F(M,N), values of the right hand side function 
c    for the Helmholtz equation at the grid points.
c
      implicit none

      integer m
      integer n

      double precision alpha
      double precision dx
      double precision dy
      double precision error
      double precision f(m,n)
      integer i
      integer j
      double precision u(m,n)
      double precision u_exact
      double precision u_exact_norm
      double precision u_norm
      double precision uxx
      double precision uyy
      double precision x
      double precision y

      dx = 2.0D+00 / dble ( m - 1 )
      dy = 2.0D+00 / dble ( n - 1 )

      u_norm = 0.0D+00

c$omp parallel 
c$omp& shared ( m, n, u ) 
c$omp& private ( i, j )

c$omp do reduction ( + : u_norm )
      do j = 1, n
        do i = 1, m
          u_norm = u_norm + u(i,j)**2
        end do
      end do
c$omp end do

c$omp end parallel

      u_norm = sqrt ( u_norm )

      error = 0.0D+00
      u_exact_norm = 0.0D+00

c$omp parallel 
c$omp& shared ( dx, dy, m, n, u ) 
c$omp& private ( i, j, u_exact, uxx, uyy, x, y ) 

c$omp do reduction ( + : error, u_exact_norm )
      do j = 1, n
        do i = 1, m
          x = -1.0D+00 + dx * dble ( i - 1 )
          y = -1.0D+00 + dy * dble ( j - 1 )
          call exact ( x, y, u_exact, uxx, uyy )
          error = error + ( u(i,j) - u_exact )**2
          u_exact_norm = u_exact_norm + u_exact**2
        end do
      end do
c$omp end do

c$omp end parallel

      error = sqrt ( error )
      u_exact_norm = sqrt ( u_exact_norm )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) 
     &  '  Computed U RMS norm :       ', u_norm
      write ( *, '(a,g14.6)' ) 
     &  '  Computed U_EXACT RMS norm : ', u_exact_norm
      write ( *, '(a,g14.6)' ) 
     &  '  RMS error :                 ', error

      return
      end
      subroutine exact ( x, y, u, uxx, uyy )

c*********************************************************************72
c
cc EXACT returns the exact value of the solution and derivatives.
c
c  Discussion:
c
c    EXACT can be used to compute the error, and to formulate the
c    appropriate boundary conditions and right hand side for any
c    desired solution.
c
c  Modified:
c
c    13 November 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, Y, the point at which the values are needed.
c
c    Output, double precision U, UXX, UYY, the values of the function,
c    and its second derivatives with respect to X and Y.
c
      implicit none

      double precision u
      double precision uxx
      double precision uyy
      double precision x
      double precision y

      u = ( 1.0D+00 - x * x ) * ( 1.0D+00 - y * y )
      uxx = -2.0D+00 * ( 1.0D+00 + y ) * ( 1.0D+00 - y )
      uyy = -2.0D+00 * ( 1.0D+00 + x ) * ( 1.0D+00 - x )

      return
      end
