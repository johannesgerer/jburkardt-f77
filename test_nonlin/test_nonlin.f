      subroutine p00_dif ( problem, n, fjac, x )

c*********************************************************************72
c
cc P00_DIF approximates the jacobian via finite differences.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    07 May 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROBLEM, the problem number.
c
c    Input, integer N, the number of equations.
c
c    Output, double precision FJAC(N,N), the approximante N by N
c    jacobian matrix.
c
c    Input, double precision X(N), the point at which the jacobian
c    is to be approximated.
c
      implicit none

      integer n

      double precision dxj
      double precision eps
      parameter ( eps = 0.0001D+00 )
      double precision f(n)
      double precision fjac(n,n)
      double precision fplus(n)
      integer i
      integer j
      integer problem
      double precision x(n)
      double precision xsave

      call p00_fx ( problem, f, n, x )

      do j = 1, n

        if ( 0.0D+00 <= x(j) ) then
          dxj = eps * ( x(j) + 1.0D+00 )
        else
          dxj = eps * ( x(j) - 1.0D+00 )
        end if

        xsave = x(j)
        x(j) = xsave + dxj

        call p00_fx ( problem, fplus, n, x )

        do i = 1, n
          fjac(i,j) = ( fplus(i) - f(i) ) / dxj
        end do

        x(j) = xsave

      end do

      return
      end
      subroutine p00_fx ( problem, f, n, x )

c*********************************************************************72
c
cc P00_FX evaluates the function for any problem.
c
c  Discussion:
c
c    Most of the problems were originally part of ACM algorithm 566, 
c    and were used to test the MINPACK package.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    John Dennis, David Gay, Phuong Vu,
c    A new nonlinear equations test problem,
c    Technical Report 83-16,
c    Mathematical Sciences Department,
c    Rice University, 1983.
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c    Noel deVilliers, David Glasser,
c    A continuation method for nonlinear regression,
c    SIAM Journal on Numerical Analysis,
c    Volume 18, Number 6, December 1981, pages 1139-1154.
c
c    Chris Fraley, 
c    Solution of nonlinear least-squares problems, 
c    Technical Report STAN-CS-1165, 
c    Computer Science Department,
c    Stanford University, 1987.
c
c    Chris Fraley, 
c    Software performance on nonlinear least-squares problems,
c    Technical Report SOL 88-17, 
c    Systems Optimization Laboratory,
c    Department of Operations Research, 
c    Stanford University, 1988.
c
c    JJ McKeown,
c    Specialized versus general-purpose algorithms for functions 
c    that are sums of squared terms,
c    Mathematical Programming,
c    Volume 9, 1975, pages 57-68.
c
c    JJ McKeown,
c    On algorithms for sums of squares problems,
c    in Towards Global Optimisation,
c    edited by Laurence Dixon, Gabor Szego,
c    North-Holland, 1975, pages 229-257,
c    ISBN: 0444109552,
c    LC: QA402.5.T7.
c
c    Jorge More, Burton Garbow, Kenneth Hillstrom,
c    Algorithm 566:
c    Testing unconstrained optimization software,
c    ACM Transactions on Mathematical Software,
c    Volume 7, Number 1, March 1981, pages 17-41.
c 
c    Douglas Salane, 
c    A continuation approach for solving large residual nonlinear
c    least squares  problems, 
c    SIAM Journal of Scientific and Statistical Computing, 
c    Volume 8, Number 4, July 1987, pages 655-671.
c
c  Parameters:
c
c    Input, integer PROBLEM, the number of the problem.
c
c    Output, double precision F(N), the value of the function at X.
c
c    Input, integer N, the number of equations.
c
c    Input, double precision X(N), the point at which F is to be evaluated.
c
      implicit none

      integer n

      double precision f(n)
      integer problem
      double precision x(n)

      if ( problem == 1 ) then
        call p01_fx ( n, x, f )
      else if ( problem == 2 ) then
        call p02_fx ( n, x, f )
      else if ( problem  ==  3 ) then
        call p03_fx ( n, x, f )
      else if ( problem == 4 ) then
        call p04_fx ( n, x, f )
      else if ( problem == 5 ) then
        call p05_fx ( n, x, f )
      else if ( problem == 6 ) then
        call p06_fx ( n, x, f )
      else if ( problem == 7 ) then
        call p07_fx ( n, x, f )
      else if ( problem == 8 ) then
        call p08_fx ( n, x, f )
      else if ( problem == 9 ) then
        call p09_fx ( n, x, f )
      else if ( problem == 10 ) then
        call p10_fx ( n, x, f )
      else if ( problem == 11 ) then
        call p11_fx ( n, x, f )
      else if ( problem == 12 ) then
        call p12_fx ( n, x, f )
      else if ( problem == 13 ) then
        call p13_fx ( n, x, f )
      else if ( problem == 14 ) then
        call p14_fx ( n, x, f )
      else if ( problem == 15 ) then
        call p15_fx ( n, x, f )
      else if ( problem == 16 ) then
        call p16_fx ( n, x, f )
      else if ( problem == 17 ) then
        call p17_fx ( n, x, f )
      else if ( problem == 18 ) then
        call p18_fx ( n, x, f )
      else if ( problem == 19 ) then
        call p19_fx ( n, x, f )
      else if ( problem == 20 ) then
        call p20_fx ( n, x, f )
      else if ( problem == 21 ) then
        call p21_fx ( n, x, f )
      else if ( problem == 22 ) then
        call p22_fx ( n, x, f )
      else if ( problem == 23 ) then
        call p23_fx ( n, x, f )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_FX - Fatal error!'
        write ( *, '(a,i6)' ) '  Illegal problem number = ', problem
        stop
      end if
     
      return
      end
      subroutine p00_jac ( problem, n, fjac, x )

c*********************************************************************72
c
cc P00_JAC evaluates the jacobian for any problem.
c
c  Discussion:
c
c    P00_JAC evaluates the matrix FJAC(I,J)  =  D F(I) / D X(J)
c    given the problem chosen, the number of equations, and the value
c    of the point X.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROBLEM, the problem number.
c
c    Input, integer N, the number of equations.
c
c    Output, double precision FJAC(N,N), the N by N jacobian matrix.
c
c    Input, double precision X(N), the point at which the jacobian
c    is to be evaluated.
c
      implicit none

      integer n

      double precision fjac(n,n)
      integer problem
      double precision x(n)

      if ( problem == 1 ) then
        call p01_jac ( fjac, n, x )
      else if ( problem == 2 ) then
        call p02_jac ( fjac, n, x )
      else if ( problem  ==  3 ) then
        call p03_jac ( fjac, n, x )
      else if ( problem == 4 ) then
        call p04_jac ( fjac, n, x )
      else if ( problem == 5 ) then
        call p05_jac ( fjac, n, x )
      else if ( problem == 6 ) then
        call p06_jac ( fjac, n, x )
      else if ( problem == 7 ) then
        call p07_jac ( fjac, n, x )
      else if ( problem == 8 ) then
        call p08_jac ( fjac, n, x )
      else if ( problem == 9 ) then
        call p09_jac ( fjac, n, x )
      else if ( problem == 10 ) then
        call p10_jac ( fjac, n, x )
      else if ( problem == 11 ) then
        call p11_jac ( fjac, n, x )
      else if ( problem == 12 ) then
        call p12_jac ( fjac, n, x )
      else if ( problem == 13 ) then
        call p13_jac ( fjac, n, x )
      else if ( problem == 14 ) then
        call p14_jac ( fjac, n, x )
      else if ( problem == 15 ) then
        call p15_jac ( fjac, n, x )
      else if ( problem == 16 ) then
        call p16_jac ( fjac, n, x )
      else if ( problem == 17 ) then
        call p17_jac ( fjac, n, x )
      else if ( problem == 18 ) then
        call p18_jac ( fjac, n, x )
      else if ( problem == 19 ) then
        call p19_jac ( fjac, n, x )
      else if ( problem == 20 ) then
        call p20_jac ( fjac, n, x )
      else if ( problem == 21 ) then
        call p21_jac ( fjac, n, x )
      else if ( problem == 22 ) then
        call p22_jac ( fjac, n, x )
      else if ( problem == 23 ) then
        call p23_jac ( fjac, n, x )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_JAC - Fatal error!'
        write ( *, '(a,i6)' ) '  Illegal problem number = ', problem
        stop
      end if

      return
      end
      subroutine p00_n ( problem, n )

c*********************************************************************72
c
cc P00_N returns the number of equations for a problem.
c
c  Discussion:
c
c    If N is positive, then N is the only value for the number
c    of equations for this problem.
c
c    If N is negative, then the absolute value of N is the
c    MINIMUM possible value for the number of equations for
c    this problem, and all larger values may also be used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROBLEM, the problem number for which N is desired.
c
c    Output, integer N, the number of equations.
c
      implicit none

      integer n
      integer problem

      if ( problem == 1 ) then
        call p01_n ( n )
      else if ( problem == 2 ) then
        call p02_n ( n )
      else if ( problem == 3 ) then
        call p03_n ( n )
      else if ( problem == 4 ) then
        call p04_n ( n )
      else if ( problem == 5 ) then
        call p05_n ( n )
      else if ( problem == 6 ) then
        call p06_n ( n )
      else if ( problem == 7 ) then
        call p07_n ( n )
      else if ( problem == 8 ) then
        call p08_n ( n )
      else if ( problem == 9 ) then
        call p09_n ( n )
      else if ( problem == 10 ) then
        call p10_n ( n )
      else if ( problem == 11 ) then
        call p11_n ( n )
      else if ( problem == 12 ) then
        call p12_n ( n )
      else if ( problem == 13 ) then
        call p13_n ( n )
      else if ( problem == 14 ) then
        call p14_n ( n )
      else if ( problem == 15 ) then
        call p15_n ( n )
      else if ( problem == 16 ) then
        call p16_n ( n )
      else if ( problem == 17 ) then
        call p17_n ( n )
      else if ( problem == 18 ) then
        call p18_n ( n )
      else if ( problem == 19 ) then
        call p19_n ( n )
      else if ( problem == 20 ) then
        call p20_n ( n )
      else if ( problem == 21 ) then
        call p21_n ( n )
      else if ( problem == 22 ) then
        call p22_n ( n )
      else if ( problem == 23 ) then
        call p23_n ( n )
      else
        n = 0
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_N - Fatal error!'
        write ( *, '(a,i6)' ) '  Illegal problem number = ', problem
        stop
      end if
     
      return
      end
      subroutine p00_problem_num ( problem_num )

c*********************************************************************72
c
cc P00_PROBLEM_NUM returns the number of problems available.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer PROBLEM_NUM, the number of problems available.
c
      implicit none

      integer problem_num

      problem_num = 23

      return
      end
      subroutine p00_sol ( problem, iknow, n, x )

c*********************************************************************72
c
cc P00_SOL returns the solution of any problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROBLEM, the problem number.
c
c    Output, integer IKNOW.
c    -1, there is no solution to the system.
c    0, the solution of the system is not known.
c    positive, IKNOW solutions are known.
c
c    Input, integer N, the order of the system.  N is primarily
c    needed for problems where N may vary.
c
c    Output, double precision X(N), the value of a solution of the system,
c   if known.
c
      implicit none

      integer n

      integer iknow
      integer problem
      double precision x(n)
c
c  Call the appropriate function-specific routine.
c
      if ( problem == 1 ) then
        call p01_sol ( iknow, n, x )
      else if ( problem == 2 ) then
        call p02_sol ( iknow, n, x )
      else if ( problem == 3 ) then
        call p03_sol ( iknow, n, x )
      else if ( problem == 4 ) then
        call p04_sol ( iknow, n, x )
      else if ( problem == 5 ) then
        call p05_sol ( iknow, n, x )
      else if ( problem == 6 ) then
        call p06_sol ( iknow, n, x )
      else if ( problem == 7 ) then
        call p07_sol ( iknow, n, x )
      else if ( problem == 8 ) then
        call p08_sol ( iknow, n, x )
      else if ( problem == 9 ) then
        call p09_sol ( iknow, n, x )
      else if ( problem == 10 ) then
        call p10_sol ( iknow, n, x )
      else if ( problem == 11 ) then
        call p11_sol ( iknow, n, x )
      else if ( problem == 12 ) then
        call p12_sol ( iknow, n, x )
      else if ( problem == 13 ) then
        call p13_sol ( iknow, n, x )
      else if ( problem == 14 ) then
        call p14_sol ( iknow, n, x )
      else if ( problem == 15 ) then
        call p15_sol ( iknow, n, x )
      else if ( problem == 16 ) then
        call p16_sol ( iknow, n, x )
      else if ( problem == 17 ) then
        call p17_sol ( iknow, n, x )
      else if ( problem == 18 ) then
        call p18_sol ( iknow, n, x )
      else if ( problem == 19 ) then
        call p19_sol ( iknow, n, x )
      else if ( problem == 20 ) then
        call p20_sol ( iknow, n, x )
      else if ( problem == 21 ) then
        call p21_sol ( iknow, n, x )
      else if ( problem == 22 ) then
        call p22_sol ( iknow, n, x )
      else if ( problem == 23 ) then
        call p23_sol ( iknow, n, x )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_SOL - Fatal error!'
        write ( *, '(a,i6)' ) '  Illegal problem number = ', problem
        stop
      end if

      return
      end
      subroutine p00_start ( problem, n, x )

c*********************************************************************72
c
cc P00_START specifies a standard approximate solution.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROBLEM, is the problem number.
c
c    Input, integer N, is the number of equations for the problem.
c
c    Output, double precision X(N), the starting point.
c
      implicit none

      integer n

      integer problem
      double precision x(n)

      if ( problem == 1 ) then
        call p01_start ( n, x )
      else if ( problem == 2 ) then
        call p02_start ( n, x )
      else if ( problem == 3 ) then
        call p03_start ( n, x )
      else if ( problem == 4 ) then
        call p04_start ( n, x )
      else if ( problem == 5 ) then
        call p05_start ( n, x )
      else if ( problem == 6 ) then
        call p06_start ( n, x )
      else if ( problem == 7 ) then
        call p07_start ( n, x )
      else if ( problem == 8 ) then
        call p08_start ( n, x )
      else if ( problem == 9 ) then
        call p09_start ( n, x )
       else if ( problem == 10 ) then
        call p10_start ( n, x )
      else if ( problem == 11 ) then
        call p11_start ( n, x )
      else if ( problem == 12 ) then
        call p12_start ( n, x )
      else if ( problem == 13 ) then
        call p13_start ( n, x )
      else if ( problem == 14 ) then
        call p14_start ( n, x )
      else if ( problem == 15 ) then
        call p15_start ( n, x )
      else if ( problem == 16 ) then
        call p16_start ( n, x )
      else if ( problem == 17 ) then
        call p17_start ( n, x )
      else if ( problem == 18 ) then
        call p18_start ( n, x )
      else if ( problem == 19 ) then
        call p19_start ( n, x )
      else if ( problem == 20 ) then
        call p20_start ( n, x )
      else if ( problem == 21 ) then
        call p21_start ( n, x )
      else if ( problem == 22 ) then
        call p22_start ( n, x )
      else if ( problem == 23 ) then
        call p23_start ( n, x )
      else

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_START - Fatal error!'
        write ( *, '(a,i6)' ) '  Illegal problem number = ', problem
        stop

      end if

      return
      end
      subroutine p00_title ( problem, title )

c*********************************************************************72
c
cc P00_TITLE returns the title of the problem.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer PROBLEM, the problem number.
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      integer problem
      character ( len = * ) title

      if ( problem == 1 ) then
        call p01_title ( title )
      else if ( problem == 2 ) then
        call p02_title ( title )
      else if ( problem == 3 ) then
        call p03_title ( title )
      else if ( problem == 4 ) then
        call p04_title ( title )
      else if ( problem == 5 ) then
        call p05_title ( title )
      else if ( problem == 6 ) then
        call p06_title ( title )
      else if ( problem == 7 ) then
        call p07_title ( title )
      else if ( problem == 8 ) then
        call p08_title ( title )
      else if ( problem == 9 ) then
        call p09_title ( title )
      else if ( problem == 10 ) then
        call p10_title ( title )
      else if ( problem == 11 ) then
        call p11_title ( title )
      else if ( problem == 12 ) then
        call p12_title ( title )
      else if ( problem == 13 ) then
        call p13_title ( title )
      else if ( problem == 14 ) then
        call p14_title ( title )
      else if ( problem == 15 ) then
        call p15_title ( title )
      else if ( problem == 16 ) then
        call p16_title ( title )
      else if ( problem == 17 ) then
        call p17_title ( title )
      else if ( problem == 18 ) then
        call p18_title ( title )
      else if ( problem == 19 ) then
        call p19_title ( title )
      else if ( problem == 20 ) then
        call p20_title ( title )
      else if ( problem == 21 ) then
        call p21_title ( title )
      else if ( problem == 22 ) then
        call p22_title ( title )
      else if ( problem == 23 ) then
        call p23_title ( title )
      else
        title = ' '
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'P00_TITLE - Fatal error!'
        write ( *, '(a,i6)' ) '  Illegal problem number = ', problem
        stop
      end if
     
      return
      end
      subroutine p01_fx ( n, x, f )

c*********************************************************************72
c
cc P01_FX evaluates the function for problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of variables.
c
c    Input, double precision X(N), the point where the function
c    is to be evaluated.
c
c    Output, double precision F(N), the value of the function.
c
      implicit none

      integer n

      double precision f(n)
      integer i
      double precision x(n)

      f(1) = 1.0D+00 - x(1)
      do i = 2, n
        f(i) = 10.0D+00 * ( x(i) - x(i-1)**2 )
      end do
     
      return
      end
      subroutine p01_n ( n )

c*********************************************************************72
c
cc P01_N returns the number of equations for problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer N.
c
c    If N is positive, then N is the only value for the number
c    of equations for this problem.
c
c    If N is negative, then the absolute value of N is the
c    MINIMUM possible value for the number of equations for
c    this problem, and all larger values may also be used.
c
      implicit none

      integer n

      n = -2

      return
      end
      subroutine p01_jac ( fjac, n, x )

c*********************************************************************72
c
cc P01_JAC sets the jacobian for problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision FJAC(N,N), the N by N jacobian matrix.
c
c    Input, integer N, the number of equations.
c
c    Input, double precision X(N), the point where the jacobian 
c    is to be evaluated.
c
      implicit none

      integer n

      double precision fjac(n,n)
      integer i
      integer j
      double precision x(n)

      do j = 1, n
        do i = 1, n
          fjac(i,j) = 0.0D+00
        end do
      end do

      fjac(1,1) = - 1.0D+00

      do i = 2, n
        fjac(i,i-1) = - 20.0D+00 * x(i-1)
        fjac(i,i) = 10.0D+00
      end do

      return
      end
      subroutine p01_sol ( iknow, n, x )

c*********************************************************************72
c
cc P01_SOL returns the solution of problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IKNOW.
c    -1, there is no solution to the system.
c    0, the solution of the system is not known.
c    positive, IKNOW solutions are known.
c
c    Input, integer N, the order of the system.  N is primarily
c    needed for problems where N may vary.
c
c    Output, double precision X(N), the value of a solution of the system,
c    if known.
c
      implicit none

      integer n

      integer i
      integer iknow
      double precision x(n)

      iknow = 1

      do i = 1, n
        x(i) = 1.0D+00
      end do

      return
      end
      subroutine p01_start ( n, x )

c*********************************************************************72
c
cc P01_START specifies a standard approximate solution for problem 1.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of equations for the problem.
c
c    Output, double precision X(N), the starting point.
c
      implicit none

      integer n

      double precision x(n)

      x(1) = - 1.2D+00
      x(2:n) = 1.0D+00

      return
      end
      subroutine p01_title ( title )

c*********************************************************************72
c
cc P01_TITLE returns the title of problem 1.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'Generalized Rosenbrock function.'
     
      return
      end
      subroutine p02_fx ( n, x, f )

c*********************************************************************72
c
cc P02_FX evaluates the function for problem 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of variables.
c
c    Input, double precision X(N), the point where the function
c    is to be evaluated.
c
c    Output, double precision F(N), the value of the function.
c
      implicit none

      integer n

      double precision f(n)
      double precision x(n)

      f(1) = x(1) + 10.0D+00 * x(2)
      f(2) = sqrt ( 5.0D+00 ) * ( x(3) - x(4) )
      f(3) = ( x(2) - 2.0D+00 * x(3) )**2
      f(4) = sqrt ( 10.0D+00 ) * ( x(1) - x(4) )**2
     
      return
      end
      subroutine p02_jac ( fjac, n, x )

c*********************************************************************72
c
cc P02_JAC sets the jacobian for problem 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision FJAC(N,N), the N by N jacobian matrix.
c
c    Input, integer N, the number of equations.
c
c    Input, double precision X(N), the point at which the jacobian
c    is to be evaluated.
c
      implicit none

      integer n

      double precision fjac(n,n)
      double precision x(n)

      fjac(1,1) = 1.0D+00
      fjac(1,2) = 10.0D+00
      fjac(1,3) = 0.0D+00
      fjac(1,4) = 0.0D+00

      fjac(2,1) = 0.0D+00
      fjac(2,2) = 0.0D+00
      fjac(2,3) = sqrt ( 5.0D+00 )
      fjac(2,4) = - sqrt ( 5.0D+00 )

      fjac(3,1) = 0.0D+00
      fjac(3,2) = 2.0D+00 * ( x(2) - 2.0D+00 * x(3) )
      fjac(3,3) = - 4.0D+00 * ( x(2) - 2.0D+00 * x(3) )
      fjac(3,4) = 0.0D+00

      fjac(4,1) = 2.0D+00 * sqrt ( 10.0D+00 ) * ( x(1) - x(4) )
      fjac(4,2) = 0.0D+00
      fjac(4,3) = 0.0D+00
      fjac(4,4) = - 2.0D+00 * sqrt ( 10.0D+00 ) * ( x(1) - x(4) )

      return
      end
      subroutine p02_n ( n )

c*********************************************************************72
c
cc P02_N returns the number of equations for problem 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer N.
c
c    If N is positive, then N is the only value for the number
c    of equations for this problem.
c
c    If N is negative, then the absolute value of N is the
c    MINIMUM possible value for the number of equations for
c    this problem, and all larger values may also be used.
c
      implicit none

      integer n

      n = 4
     
      return
      end
      subroutine p02_sol ( iknow, n, x )

c*********************************************************************72
c
cc P02_SOL returns the solution of problem 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IKNOW.
c    -1, there is no solution to the system.
c    0, the solution of the system is not known.
c    positive, IKNOW solutions are known.
c
c    Input, integer N, the order of the system.  N is primarily
c    needed for problems where N may vary.
c
c    Output, double precision X(N), the value of a solution of the system,
c    if known.
c
      implicit none

      integer n

      integer i
      integer iknow
      double precision x(n)

      iknow = 1

      do i = 1, n
        x(i) = 0.0D+00
      end do

      return
      end
      subroutine p02_start ( n, x )

c*********************************************************************72
c
cc P02_START specifies a standard approximate solution for problem 2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of equations for the problem.
c
c    Output, double precision X(N), the starting point.
c
      implicit none

      integer n

      double precision x(n)

      x(1) =  3.0D+00
      x(2) = -1.0D+00
      x(3) =  0.0D+00
      x(4) =  1.0D+00

      return
      end
      subroutine p02_title ( title )

c*********************************************************************72
c
cc P02_TITLE returns the title of problem 2.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'Powell singular function.'

      return
      end
      subroutine p03_fx ( n, x, f )

c*********************************************************************72
c
cc P03_FX evaluates the function for problem 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of variables.
c
c    Input, double precision X(N), the point where the function
c    is to be evaluated.
c
c    Output, double precision F(N), the value of the function.
c
      implicit none

      integer n

      double precision f(n)
      double precision x(n)

      f(1) = 10000.0D+00 * x(1) * x(2) - 1.0D+00
      f(2) = exp ( - x(1) ) + exp ( - x(2) ) - 1.0001D+00
     
      return
      end
      subroutine p03_jac ( fjac, n, x )

c*********************************************************************72
c
cc P03_JAC sets the jacobian for problem 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision FJAC(N,N), the N by N jacobian matrix.
c
c    Input, integer N, the number of equations.
c
c    Input, double precision X(N), the point at which the jacobian
c    is to be evaluated.
c
      implicit none

      integer n

      double precision fjac(n,n)
      double precision x(n)

      fjac(1,1) = 10000.0D+00 * x(2)
      fjac(1,2) = 10000.0D+00 * x(1)

      fjac(2,1) = - exp ( - x(1) )
      fjac(2,2) = - exp ( - x(2) )

      return
      end
      subroutine p03_n ( n )

c*********************************************************************72
c
cc P03_N returns the number of equations for problem 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer N.
c
c    If N is positive, then N is the only value for the number
c    of equations for this problem.
c
c    If N is negative, then the absolute value of N is the
c    MINIMUM possible value for the number of equations for
c    this problem, and all larger values may also be used.
c
      implicit none

      integer n

      n = 2
     
      return
      end
      subroutine p03_sol ( iknow, n, x )

c*********************************************************************72
c
cc P03_SOL returns the solution of problem 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IKNOW.
c    -1, there is no solution to the system.
c    0, the solution of the system is not known.
c    positive, IKNOW solutions are known.
c
c    Input, integer N, the order of the system.  N is primarily
c    needed for problems where N may vary.
c
c    Output, double precision X(N), the value of a solution of the system,
c    if known.
c
      implicit none

      integer n

      integer iknow
      double precision x(n)

      iknow = 1

      x(1:2) = (/ 1.098159D-05, 9.106146D+00 /)

      return
      end
      subroutine p03_start ( n, x )

c*********************************************************************72
c
cc P03_START specifies a standard approximate solution for problem 3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of equations for the problem.
c
c    Output, double precision X(N), the starting point.
c
      implicit none

      integer n

      double precision x(n)

      x(1:2) = (/ 0.0D+00, 1.0D+00 /)

      return
      end
      subroutine p03_title ( title )

c*********************************************************************72
c
cc P03_TITLE returns the title of problem 3.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'Powell badly scaled function.'

      return
      end
      subroutine p04_fx ( n, x, f )

c*********************************************************************72
c
cc P04_FX evaluates the function for problem 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of variables.
c
c    Input, double precision X(N), the point where the function
c    is to be evaluated.
c
c    Output, double precision F(N), the value of the function.
c
      implicit none

      integer n

      double precision f(n)
      double precision temp1
      double precision temp2
      double precision x(n)

      temp1 = x(2) - x(1)**2
      temp2 = x(4) - x(3)**2

      f(1) = - 200.0D+00 * x(1) * temp1 - ( 1.0D+00 - x(1) )

      f(2) = 200.0D+00 * temp1 + 20.2D+00 * ( x(2) - 1.0D+00 ) 
     &  + 19.8D+00 * ( x(4) - 1.0D+00 )

      f(3) = - 180.0D+00 * x(3) * temp2 - ( 1.0D+00 - x(3) )

      f(4) = 180.0D+00 * temp2 + 20.2D+00 * ( x(4) - 1.0D+00 )
     &  + 19.8D+00 * ( x(2) - 1.0D+00 )
     
      return
      end
      subroutine p04_jac ( fjac, n, x )

c*********************************************************************72
c
cc P04_JAC sets the jacobian for problem 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision FJAC(N,N), the N by N jacobian matrix.
c
c    Input, integer N, the number of equations.
c
c    Input, double precision X(N), the point at which the jacobian
c    is to be evaluated.
c
      implicit none

      integer n

      double precision fjac(n,n)
      double precision x(n)

      fjac(1,1) = - 200.0D+00 * ( x(2) - 3.0D+00 * x(1)**2 ) + 1.0D+00
      fjac(1,2) = - 200.0D+00 * x(1)
      fjac(1,3) = 0.0D+00
      fjac(1,4) = 0.0D+00

      fjac(2,1) = - 400.0D+00 * x(1)
      fjac(2,2) =   220.2D+00
      fjac(2,3) =     0.0D+00
      fjac(2,4) =    19.8D+00

      fjac(3,1) = 0.0D+00
      fjac(3,2) = 0.0D+00
      fjac(3,3) = - 180.0D+00 * ( x(4) - 3.0D+00 * x(3)**2 ) + 1.0D+00
      fjac(3,4) = - 180.0D+00 * x(3)

      fjac(4,1) =     0.0D+00
      fjac(4,2) =    19.8D+00
      fjac(4,3) = - 360.0D+00 * x(3)
      fjac(4,4) =   300.2D+00

      return
      end
      subroutine p04_n ( n )

c*********************************************************************72
c
cc P04_N returns the number of equations for problem 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer N.
c
c    If N is positive, then N is the only value for the number
c    of equations for this problem.
c
c    If N is negative, then the absolute value of N is the
c    MINIMUM possible value for the number of equations for
c    this problem, and all larger values may also be used.
c
      implicit none

      integer n

      n = 4

      return
      end
      subroutine p04_sol ( iknow, n, x )

c*********************************************************************72
c
cc P04_SOL returns the solution of problem 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IKNOW.
c    -1, there is no solution to the system.
c    0, the solution of the system is not known.
c    positive, IKNOW solutions are known.
c
c    Input, integer N, the order of the system.  N is primarily
c    needed for problems where N may vary.
c
c    Output, double precision X(N), the value of a solution of the system,
c    if known.
c
      implicit none

      integer n

      integer iknow
      double precision x(n)

      iknow = 1

      x(1:4) = 1.0D+00

      return
      end
      subroutine p04_start ( n, x )

c*********************************************************************72
c
cc P04_START specifies a standard approximate solution for problem 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of equations for the problem.
c
c    Output, double precision X(N), the starting point.
c
      implicit none

      integer n

      double precision x(n)

      x(1:4) = (/ -3.0D+00, -1.0D+00, -3.0D+00, -1.0D+00 /)

      return
      end
      subroutine p04_title ( title )

c*********************************************************************72
c
cc P04_TITLE returns the title of problem 4.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'Wood function.'

      return
      end
      subroutine p05_fx ( n, x, f )

c*********************************************************************72
c
cc P05_FX evaluates the function for problem 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of variables.
c
c    Input, double precision X(N), the point where the function
c    is to be evaluated.
c
c    Output, double precision F(N), the value of the function.
c
      implicit none

      integer n

      double precision f(n)
      double precision :: pi = 3.141592653589793D+00
      double precision temp
      double precision x(n)

      if ( 0.0D+00 .lt. x(1) ) then
        temp = atan ( x(2) / x(1) ) / ( 2.0D+00 * pi )
      else if (x(1) .lt. 0.0D+00 ) then
        temp = atan ( x(2) / x(1) ) / ( 2.0D+00 * pi ) + 0.5D+00
      else
        temp = sign ( 0.25D+00, x(2) )
      end if
     
      f(1) = 10.0D+00 * ( x(3) - 10.0D+00 * temp )

      f(2) = 10.0D+00 * ( sqrt ( x(1)**2 + x(2)**2 ) - 1.0D+00 )

      f(3) = x(3)
     
      return
      end
      subroutine p05_jac ( fjac, n, x )

c*********************************************************************72
c
cc P05_JAC sets the jacobian for problem 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision FJAC(N,N), the N by N jacobian matrix.
c
c    Input, integer N, the number of equations.
c
c    Input, double precision X(N), the point at which the jacobian
c    is to be evaluated.
c
      implicit none

      integer n

      double precision fjac(n,n)
      double precision :: pi = 3.141592653589793D+00
      double precision x(n)

      fjac(1,1) =  100.0D+00 * x(2) 
     &  / ( 2.0D+00 * pi * ( x(1)**2 + x(2)**2 ) )

      fjac(1,2) = -100.0D+00 * x(1) 
     &  / ( 2.0D+00 * pi * ( x(1)**2 + x(2)**2 ) )

      fjac(1,3) =   10.0D+00

      fjac(2,1) = 10.0D+00 * x(1) / sqrt ( x(1)**2 + x(2)**2 )
      fjac(2,2) = 10.0D+00 * x(2) / sqrt ( x(1)**2 + x(2)**2 )
      fjac(2,3) = 0.0D+00

      fjac(3,1) = 0.0D+00
      fjac(3,2) = 0.0D+00
      fjac(3,3) = 1.0D+00

      return
      end
      subroutine p05_n ( n )

c*********************************************************************72
c
cc P05_N returns the number of equations for problem 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer N.
c
c    If N is positive, then N is the only value for the number
c    of equations for this problem.
c
c    If N is negative, then the absolute value of N is the
c    MINIMUM possible value for the number of equations for
c    this problem, and all larger values may also be used.
c
      implicit none

      integer n

      n = 3

      return
      end
      subroutine p05_sol ( iknow, n, x )

c*********************************************************************72
c
cc P05_SOL returns the solution of problem 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IKNOW.
c    -1, there is no solution to the system.
c    0, the solution of the system is not known.
c    positive, IKNOW solutions are known.
c
c    Input, integer N, the order of the system.  N is primarily
c    needed for problems where N may vary.
c
c    Output, double precision X(N), the value of a solution of the system,
c    if known.
c
      implicit none

      integer n

      integer iknow
      double precision x(n)

      iknow = 1

      x(1:3) = (/ 1.0D+00, 0.0D+00, 0.0D+00 /)

      return
      end
      subroutine p05_start ( n, x )

c*********************************************************************72
c
cc P05_START specifies a standard approximate solution for problem 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, is the number of equations for the problem.
c
c    Output, double precision X(N), the starting point.
c
      implicit none

      integer n

      double precision x(n)

      x(1:3) = (/ -1.0D+00, 0.0D+00, 0.0D+00 /)

      return
      end
      subroutine p05_title ( title )

c*********************************************************************72
c
cc P05_TITLE returns the title of problem 5.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'Helical valley function.'

      return
      end
      subroutine p06_fx ( n, x, f )

c*********************************************************************72
c
cc P06_FX evaluates the function for problem 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of variables.
c
c    Input, double precision X(N), the point where the function
c    is to be evaluated.
c
c    Output, double precision F(N), the value of the function.
c
      implicit none

      integer n

      double precision f(n)
      integer i
      integer j
      integer k
      double precision sum1
      double precision sum2
      double precision temp
      double precision ti
      double precision x(n)

      do i = 1, n
        f(i) = 0.0D+00
      end do
   
      do i = 1, 29

        ti = dble ( i ) / 29.0D+00

        sum1 = 0.0D+00
        temp = 1.0D+00
        do j = 2, n
          sum1 = sum1 + dble ( j - 1 ) * temp * x(j)
          temp = ti * temp
        end do

        sum2 = 0.0D+00
        temp = 1.0D+00
        do j = 1, n
          sum2 = sum2 + temp * x(j)
          temp = ti * temp
        end do

        temp = 1.0D+00 / ti

        do k = 1, n
          f(k) = f(k) + temp * ( sum1 - sum2**2 - 1.0D+00 ) 
     &      * ( dble ( k - 1 ) - 2.0D+00 * ti * sum2 )
          temp = ti * temp
        end do

      end do
     
      f(1) = f(1) + 3.0D+00 * x(1) - 2.0D+00 * x(1) 
     &  * x(2) + 2.0D+00 * x(1)**3 
      f(2) = f(2) + x(2) - x(1)**2 - 1.0D+00
     
      return
      end
      subroutine p06_jac ( fjac, n, x )

c*********************************************************************72
c
cc P06_JAC sets the jacobian for problem 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 June 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision FJAC(N,N), the N by N jacobian matrix.
c
c    Input, integer N, the number of equations.
c
c    Input, double precision X(N), the point at which the jacobian
c    is to be evaluated.
c
      implicit none

      integer n

      double precision fjac(n,n)
      integer i
      integer j
      integer k
      double precision sum1
      double precision sum2
      double precision temp
      double precision temp1
      double precision ti
      double precision tj
      double precision tk
      double precision x(n)

      do j = 1, n
        do i = 1, n
          fjac(i,j) = 0.0D+00
        end do
      end do

      do i = 1, 29

        ti = dble ( i ) / 29.0D+00

        sum1 = 0.0D+00
        temp = 1.0D+00
        do j = 2, n
          sum1 = sum1 + dble ( j - 1 ) * temp * x(j)
          temp = ti * temp
        end do

        sum2 = 0.0D+00
        temp = 1.0D+00
        do j = 1, n
          sum2 = sum2 + temp * x(j)
          temp = ti * temp
        end do

        temp1 = 2.0D+00 * ( sum1 - sum2**2 - 1.0D+00 )
        tk = 1.0D+00

        do k = 1, n
          tj = tk
          do j = k, n
            fjac(k,j) = fjac(k,j) 
     &        + tj * ( ( dble ( k - 1 ) / ti - 2.0D+00 * sum2 ) 
     &        * ( dble ( j - 1 ) / ti - 2.0D+00 * sum2 ) - temp1 )
            tj = ti * tj
          end do
          tk = ti**2 * tk
        end do

      end do

      fjac(1,1) = fjac(1,1) + 3.0D+00 
     &  - 2.0D+00 * x(2) + 6.0D+00 * x(1)**2  
      fjac(1,2) = fjac(1,2) - 2.0D+00 * x(1)

      fjac(2,1) = fjac(2,1) - 2.0D+00 * x(1)
      fjac(2,2) = fjac(2,2) + 1.0D+00

      do k = 1, n
        fjac(k:n,k) = fjac(k,k:n)
      end do

      return
      end
      subroutine p06_n ( n )

c*********************************************************************72
c
cc P06_N returns the number of equations for problem 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer N.
c
c    If N is positive, then N is the only value for the number
c    of equations for this problem.
c
c    If N is negative, then the absolute value of N is the
c    MINIMUM possible value for the number of equations for
c    this problem, and all larger values may also be used.
c
      implicit none

      integer n

      n = - 2

      return
      end
      subroutine p06_sol ( iknow, n, x )

c*********************************************************************72
c
cc P06_SOL returns the solution of problem 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IKNOW.
c    -1, there is no solution to the system.
c    0, the solution of the system is not known.
c    positive, IKNOW solutions are known.
c
c    Input, integer N, the order of the system.  N is primarily
c    needed for problems where N may vary.
c
c    Output, double precision X(N), the value of a solution of the system,
c    if known.
c
      implicit none

      integer n

      integer i
      integer iknow
      double precision x(n)

      iknow = 0

      do i = 1, n
        x(i) = 0.0D+00
      end do

      return
      end
      subroutine p06_start ( n, x )

c*********************************************************************72
c
cc P06_START specifies a standard approximate solution for problem 6.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, is the number of equations for the problem.
c
c    Output, double precision X(N), the starting point.
c
      implicit none

      integer n

      integer i
      double precision x(n)

      do i = 1, n
        x(i) = 0.0D+00
      end do

      return
      end
      subroutine p06_title ( title )

c*********************************************************************72
c
cc P06_TITLE returns the title of problem 6.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'Watson function.'

      return
      end
      subroutine p07_fx ( n, x, f )

c*********************************************************************72
c
cc P07_FX evaluates the function for problem 7.
c
c  Discussion:
c
c    The Chebyquad function is related to the computation of the 
c    abscissas for Chebyshev quadrature.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of variables.
c
c    Input, double precision X(N), the point where the function
c    is to be evaluated.
c
c    Output, double precision F(N), the value of the function.
c
      implicit none

      integer n

      double precision f(n)
      integer i
      integer, parameter :: i4_2 = 2
      integer j
      double precision t1
      double precision t2
      double precision t3
      double precision x(n)

      do i = 1, n
        f(i) = 0.0D+00
      end do
   
      do j = 1, n

        t1 = 1.0D+00
        t2 = x(j)

        do i = 1, n

          f(i) = f(i) + t2

          t3 = 2.0D+00 * x(j) * t2 - t1
          t1 = t2
          t2 = t3

        end do

      end do

      do i = 1, n
        f(i) = f(i) / dble ( n )
      end do

      do i = 1, n

        if ( mod ( i, i4_2 ) == 0 ) then
          f(i) = f(i) + 1.0D+00 / dble ( i * i - 1 )
        end if

      end do
     
      return
      end
      subroutine p07_jac ( fjac, n, x )

c*********************************************************************72
c
cc P07_JAC sets the jacobian for problem 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision FJAC(N,N), the N by N jacobian matrix.
c
c    Input, integer N, the number of equations.
c
c    Input, double precision X(N), the point at which the jacobian
c    is to be evaluated.
c
      implicit none

      integer n

      double precision fjac(n,n)
      integer i
      integer j
      double precision t1
      double precision t2
      double precision t3
      double precision t4
      double precision t5
      double precision t6
      double precision x(n)

      do j = 1, n

        t1 = 1.0D+00
        t2 = x(j)

        t4 = 0.0D+00
        t5 = 1.0D+00

        do i = 1, n

          fjac(i,j) = t5 

          t6 = 2.0D+00 * t2 + 2.0D+00 * t5 * x(j) - t4
          t4 = t5
          t5 = t6

          t3 = 2.0D+00 * x(j) * t2 - t1
          t1 = t2
          t2 = t3

        end do

      end do

      do j = 1, n
        do i = 1, n
          fjac(i,j) = fjac(i,j) / dble ( n )
        end do
      end do

      return
      end
      subroutine p07_n ( n )

c*********************************************************************72
c
cc P07_N returns the number of equations for problem 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer N.
c
c    If N is positive, then N is the only value for the number
c    of equations for this problem.
c
c    If N is negative, then the absolute value of N is the
c    MINIMUM possible value for the number of equations for
c    this problem, and all larger values may also be used.
c
      implicit none

      integer n
c
c  Only the values N = 1, 2, 3, 4, 5, 6, 7 and 9 may be used.
c
      n = 9

      return
      end
      subroutine p07_sol ( iknow, n, x )

c*********************************************************************72
c
cc P07_SOL returns the solution of problem 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IKNOW.
c    -1, there is no solution to the system.
c    0, the solution of the system is not known.
c    positive, IKNOW solutions are known.
c
c    Input, integer N, the order of the system.  N is primarily
c    needed for problems where N may vary.
c
c    Output, double precision X(N), the value of a solution of the system,
c    if known.
c
      implicit none

      integer n

      integer i
      integer iknow
      double precision x(n)

      if ( n == 1 ) then
        iknow = 1
        x(1) = 0.5000000000000000D+00   
      else if ( n == 2 ) then
        iknow = 1
        x(1) = 0.2113248654051871D+00   
        x(2) = 0.7886751345948129D+00   
      else if ( n == 3 ) then
        iknow = 1
        x(1) = 0.1464466094067263D+00   
        x(2) = 0.5000000000000000D+00  
        x(3) = 0.8535533905932737D+00   
      else if ( n == 4 ) then
        iknow = 1
        x(1) = 0.1026727638500000D+00 
        x(2) = 0.4062037629500000D+00   
        x(3) = 0.5937962370500000D+00  
        x(4) = 0.8973272361500000D+00  
      else if ( n == 5 ) then
        iknow = 1
        x(1) = 8.3751256499999982D-02
        x(2) = 0.3127292952000000D+00  
        x(3) = 0.5000000000000000D+00 
        x(4) = 0.6872707048000000D+00  
        x(5) = 0.9162487435000000D+00  
      else if ( n == 6 ) then
        iknow = 1
        x(1) = 6.6876590949999981D-02
        x(2) = 0.2887406731000000D+00  
        x(3) = 0.3666822992500000D+00  
        x(4) = 0.6333177007499999D+00  
        x(5) = 0.7112593269000000D+00  
        x(6) = 0.9331234090500000D+00  
      else if ( n == 7 ) then
        iknow = 1
        x(1) = 5.8069149599999981D-02
        x(2) = 0.2351716123500000D+00 
        x(3) = 0.3380440947500000D+00  
        x(4) = 0.5000000000000000D+00  
        x(5) = 0.6619559052500000D+00  
        x(6) = 0.7648283876499999D+00  
        x(7) = 0.9419308504000000D+00 
      else if ( n == 9 ) then
        iknow = 1
        x(1) = 4.4205346149999991D-02
        x(2) = 0.1994906723000000D+00  
        x(3) = 0.2356191084500000D+00  
        x(4) = 0.4160469079000000D+00  
        x(5) = 0.5000000000000000D+00   
        x(6) = 0.5839530921000000D+00
        x(7) = 0.7643808915500000D+00
        x(8) = 0.8005093276999999D+00
        x(9) = 0.9557946538500000D+00  
      else
        iknow = - 1
        do i = 1, n
          x(i) = 0.5D+00
        end do
      end if

      do i = 1, n
        x(i) = 2.0D+00 * x(i) - 1.0D+00
      end do

      return
      end
      subroutine p07_start ( n, x )

c*********************************************************************72
c
cc P07_START specifies a standard approximate solution for problem 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, is the number of equations for the problem.
c
c    Output, double precision X(N), the starting point.
c
      implicit none

      integer n

      integer i
      double precision x(n)

      do i = 1, n
        x(i) = dble ( 2 * i - 1 - n ) / dble ( n + 1 )
      end do

      return
      end
      subroutine p07_title ( title )

c*********************************************************************72
c
cc P07_TITLE returns the title of problem 7.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'Chebyquad function.'
     
      return
      end
      subroutine p08_fx ( n, x, f )

c*********************************************************************72
c
cc P08_FX evaluates the function for problem 8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of variables.
c
c    Input, double precision X(N), the point where the function
c    is to be evaluated.
c
c    Output, double precision F(N), the value of the function.
c
      implicit none

      integer n

      double precision f(n)
      double precision x(n)

      f(1:n-1) = x(1:n-1) + sum ( x(1:n) ) - dble ( n + 1 )
     
      f(n) = product ( x(1:n) ) - 1.0D+00
     
      return
      end
      subroutine p08_jac ( fjac, n, x )

c*********************************************************************72
c
cc P08_JAC sets the jacobian for problem 8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision FJAC(N,N), the N by N jacobian matrix.
c
c    Input, integer N, the number of equations.
c
c    Input, double precision X(N), the point at which the jacobian
c    is to be evaluated.
c
      implicit none

      integer n

      double precision fjac(n,n)
      integer i
      integer j
      integer k
      double precision prod
      double precision x(n)

      do j = 1, n
        do i = 1, n - 1
          fjac(i,j) = 1.0D+00
        end do
      end do

      do i = 1, n-1
        fjac(i,i) = 2.0D+00
      end do
c
c  Last row:
c
      do j = 1, n
        prod = 1.0D+00
        do k = 1, n
          if ( k .ne. j ) then
            prod = x(k) * prod
          end if
        end do
        fjac(n,j) = prod
      end do

      return
      end
      subroutine p08_n ( n )

c*********************************************************************72
c
cc P08_N returns the number of equations for problem 8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer N.
c
c    If N is positive, then N is the only value for the number
c    of equations for this problem.
c
c    If N is negative, then the absolute value of N is the
c    MINIMUM possible value for the number of equations for
c    this problem, and all larger values may also be used.
c
      implicit none

      integer n

      n = - 1

      return
      end
      subroutine p08_sol ( iknow, n, x )

c*********************************************************************72
c
cc P08_SOL returns the solution of problem 8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IKNOW.
c    -1, there is no solution to the system.
c    0, the solution of the system is not known.
c    positive, IKNOW solutions are known.
c
c    Input, integer N, the order of the system.  N is primarily
c    needed for problems where N may vary.
c
c    Output, double precision X(N), the value of a solution of the system,
c    if known.
c
      implicit none

      integer n

      integer iknow
      double precision x(n)

      iknow = 1

      x(1:n) = 1.0D+00

      return
      end
      subroutine p08_start ( n, x )

c*********************************************************************72
c
cc P08_START specifies a standard approximate solution for problem 8.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, is the number of equations for the problem.
c
c    Output, double precision X(N), the starting point.
c
      implicit none

      integer n

      double precision x(n)

      x(1:n) = 0.5D+00

      return
      end
      subroutine p08_title ( title )

c*********************************************************************72
c
cc P08_TITLE returns the title of problem 8.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'Brown almost linear function.'
     
      return
      end
      subroutine p09_fx ( n, x, f )

c*********************************************************************72
c
cc P09_FX evaluates the function for problem 9.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of variables.
c
c    Input, double precision X(N), the point where the function
c    is to be evaluated.
c
c    Output, double precision F(N), the value of the function.
c
      implicit none

      integer n

      double precision f(n)
      double precision h
      integer k
      double precision x(n)

      h = 1.0D+00 / dble ( n + 1 )
     
      do k = 1, n
     
        f(k) = 2.0D+00 * x(k) + 0.5D+00 * h**2 
     &    * ( x(k) + dble ( k ) * h + 1.0D+00 )**3
     
        if ( 1 .lt. k ) then
          f(k) = f(k) - x(k-1)
        end if
     
        if ( k .lt. n ) then
          f(k) = f(k) - x(k+1)
        end if
     
      end do
     
      return
      end
      subroutine p09_jac ( fjac, n, x )

c*********************************************************************72
c
cc P09_JAC sets the jacobian for problem 9.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision FJAC(N,N), the N by N jacobian matrix.
c
c    Input, integer N, the number of equations.
c
c    Input, double precision X(N), the point at which the jacobian
c    is to be evaluated.
c
      implicit none

      integer n

      double precision fjac(n,n)
      integer i
      integer j
      double precision x(n)

      do j = 1, n
        do i = 1, n
          fjac(i,j) = 0.0D+00
        end do
      end do

      do i = 1, n
     
        fjac(i,i) = 2.0D+00 + 1.5D+00 * ( x(i) + 1.0D+00 + dble ( i ) 
     &    / dble ( n + 1 )  )**2 / dble ( n + 1 )**2

        if ( 1 .lt. i ) then
          fjac(i,i-1) = - 1.0D+00
        end if

        if ( i .lt. n ) then
          fjac(i,i+1) = - 1.0D+00
       end if

      end do

      return
      end
      subroutine p09_n ( n )

c*********************************************************************72
c
cc P09_N returns the number of equations for problem 9.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer N.
c
c    If N is positive, then N is the only value for the number
c    of equations for this problem.
c
c    If N is negative, then the absolute value of N is the
c    MINIMUM possible value for the number of equations for
c    this problem, and all larger values may also be used.
c
      implicit none

      integer n

      n = - 1

      return
      end
      subroutine p09_sol ( iknow, n, x )

c*********************************************************************72
c
cc P09_SOL returns the solution of problem 9.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IKNOW.
c    -1, there is no solution to the system.
c    0, the solution of the system is not known.
c    positive, IKNOW solutions are known.
c
c    Input, integer N, the order of the system.  N is primarily
c    needed for problems where N may vary.
c
c    Output, double precision X(N), the value of a solution of the system,
c    if known.
c
      implicit none

      integer n

      integer i
      integer iknow
      double precision x(n)

      iknow = 0

      do i = 1, n
        x(i) = 0.0D+00
      end do

      return
      end
      subroutine p09_start ( n, x )

c*********************************************************************72
c
cc P09_START specifies a standard approximate solution for problem 9.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, is the number of equations for the problem.
c
c    Output, double precision X(N), the starting point.
c
      implicit none

      integer n

      integer i
      double precision x(n)

      do i = 1, n
        x(i) = dble ( i * ( i - n - 1 ) ) / dble ( n + 1 )**2
      end do

      return
      end
      subroutine p09_title ( title )

c*********************************************************************72
c
cc P09_TITLE returns the title of problem 9.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'Discrete boundary value function.'

      return
      end
      subroutine p10_fx ( n, x, f )

c*********************************************************************72
c
cc P10_FX evaluates the function for problem 10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of variables.
c
c    Input, double precision X(N), the point where the function 
c    is to be evaluated.
c
c    Output, double precision F(N), the value of the function.
c
      implicit none

      integer n

      double precision f(n)
      double precision h
      integer j
      integer k
      double precision sum1
      double precision sum2
      double precision tj
      double precision tk
      double precision x(n)

      h = 1.0D+00 / dble ( n + 1 )
     
      do k = 1, n
     
        tk = dble ( k ) / dble ( n + 1 )
     
        sum1 = 0.0D+00
        do j = 1, k
          tj = dble ( j ) * h
          sum1 = sum1 + tj * ( x(j) + tj + 1.0D+00 )**3
        end do
     
        sum2 = 0.0D+00
        do j = k+1, n
          tj = dble ( j ) * h
          sum2 = sum2 + ( 1.0D+00 - tj ) * ( x(j) + tj + 1.0D+00 )**3
        end do
     
        f(k) = x(k) + h * ( ( 1.0D+00 - tk ) * sum1 + tk * sum2 ) 
     &    / 2.0D+00
     
      end do
     
      return
      end
      subroutine p10_jac ( fjac, n, x )

c*********************************************************************72
c
cc P10_JAC sets the jacobian for problem 10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision FJAC(N,N), the N by N jacobian matrix.
c
c    Input, integer N, the number of equations.
c
c    Input, double precision X(N), the point at which the jacobian
c    is to be evaluated.
c
      implicit none

      integer n

      double precision fjac(n,n)
      integer i
      integer j
      double precision temp1
      double precision temp2
      double precision ti
      double precision tj
      double precision x(n)

      do i = 1, n

        ti = dble ( i ) / dble ( n + 1 )

        do j = 1, n
          tj = dble ( j ) / dble ( n + 1 )
          temp1 = ( x(j) + tj + 1.0D+00 )**2
          temp2 = min ( ti, tj ) - ti * tj
          fjac(i,j) = 1.5D+00 * temp2 * temp1 / dble ( n + 1 )
        end do

        fjac(i,i) = fjac(i,i) + 1.0D+00

      end do

      return
      end
      subroutine p10_n ( n )

c*********************************************************************72
c
cc P10_N returns the number of equations for problem 10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer N.
c
c    If N is positive, then N is the only value for the number
c    of equations for this problem.
c
c    If N is negative, then the absolute value of N is the
c    MINIMUM possible value for the number of equations for
c    this problem, and all larger values may also be used.
c
      implicit none

      integer n

      n = - 1

      return
      end
      subroutine p10_sol ( iknow, n, x )

c*********************************************************************72
c
cc P10_SOL returns the solution of problem 10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IKNOW.
c    -1, there is no solution to the system.
c    0, the solution of the system is not known.
c    positive, IKNOW solutions are known.
c
c    Input, integer N, the order of the system.  N is primarily
c    needed for problems where N may vary.
c
c    Output, double precision X(N), the value of a solution of the system,
c    if known.
c
      implicit none

      integer n

      integer i
      integer iknow
      double precision x(n)

      iknow = 0

      do i = 1, n
        x(i) = 0.0D+00
      end do

      return
      end
      subroutine p10_start ( n, x )

c*********************************************************************72
c
cc P10_START specifies a standard approximate solution for problem 10.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, is the number of equations for the problem.
c
c    Output, double precision X(N), the starting point.
c
      implicit none

      integer n

      integer i
      double precision x(n)

      do i = 1, n
        x(i) = dble ( i * ( i - n - 1 ) ) / dble ( n + 1 )**2
      end do

      return
      end
      subroutine p10_title ( title )

c*********************************************************************72
c
cc P10_TITLE returns the title of problem 10.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'Discrete integral equation function.'
     
      return
      end
      subroutine p11_fx ( n, x, f )

c*********************************************************************72
c
cc P11_FX evaluates the function for problem 11.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of variables.
c
c    Input, double precision X(N), the point where the function
c    is to be evaluated.
c
c    Output, double precision F(N), the value of the function.
c
      implicit none

      integer n

      double precision c_sum
      double precision f(n)
      integer k
      double precision x(n)

      c_sum = sum ( cos ( x(1:n) ) )
     
      do k = 1, n
        f(k) =  dble ( n ) 
     &    - c_sum + dble ( k ) * ( 1.0D+00 - cos ( x(k) ) ) 
     &    - sin ( x(k) )
      end do
     
      return
      end
      subroutine p11_jac ( fjac, n, x )

c*********************************************************************72
c
cc P11_JAC sets the jacobian for problem 11.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision FJAC(N,N), the N by N jacobian matrix.
c
c    Input, integer N, the number of equations.
c
c    Input, double precision X(N), the point at which the jacobian
c    is to be evaluated.
c
      implicit none

      integer n

      double precision fjac(n,n)
      integer i
      integer j
      double precision x(n)

      do i = 1, n

        do j = 1, n

          if ( i .ne. j ) then
            fjac(i,j) = sin ( x(j) )
          else
            fjac(i,j) = dble ( j + 1 ) * sin ( x(j) ) - cos ( x(j) )
          end if

        end do

      end do

      return
      end
      subroutine p11_n ( n )

c*********************************************************************72
c
cc P11_N returns the number of equations for problem 11.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer N.
c
c    If N is positive, then N is the only value for the number
c    of equations for this problem.
c
c    If N is negative, then the absolute value of N is the
c    MINIMUM possible value for the number of equations for
c    this problem, and all larger values may also be used.
c
      implicit none

      integer n

      n = - 1

      return
      end
      subroutine p11_sol ( iknow, n, x )

c*********************************************************************72
c
cc P11_SOL returns the solution of problem 11.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IKNOW.
c    -1, there is no solution to the system.
c    0, the solution of the system is not known.
c    positive, IKNOW solutions are known.
c
c    Input, integer N, the order of the system.  N is primarily
c    needed for problems where N may vary.
c
c    Output, double precision X(N), the value of a solution of the system,
c    if known.
c
      implicit none

      integer n

      integer i
      integer iknow
      double precision x(n)

      iknow = 0

      do i = 1, n
        x(i) = 0.0D+00
      end do

      return
      end
      subroutine p11_start ( n, x )

c*********************************************************************72
c
cc P11_START specifies a standard approximate solution for problem 11.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, is the number of equations for the problem.
c
c    Output, double precision X(N), the starting point.
c
      implicit none

      integer n

      integer i
      double precision x(n)

      do i = 1, n
        x(i) = 1.0D+00 / dble ( n )
      end do

      return
      end
      subroutine p11_title ( title )

c*********************************************************************72
c
cc P11_TITLE returns the title of problem 11.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'Trigonometric function.'
     
      return
      end
      subroutine p12_fx ( n, x, f )

c*********************************************************************72
c
cc P12_FX evaluates the function for problem 12.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of variables.
c
c    Input, double precision X(N), the point where the function
c    is to be evaluated.
c
c    Output, double precision F(N), the value of the function.
c
      implicit none

      integer n

      double precision f(n)
      integer j
      integer k
      double precision sum1
      double precision x(n)

      sum1 = 0.0D+00
      do j = 1, n
        sum1 = sum1 + dble ( j ) * ( x(j) - 1.0D+00 )
      end do
     
      do k = 1, n
        f(k) = x(k) - 1.0D+00 + dble ( k ) * sum1 
     &    * ( 1.0D+00 + 2.0D+00 * sum1**2 )
      end do
     
      return
      end
      subroutine p12_jac ( fjac, n, x )

c*********************************************************************72
c
cc P12_JAC sets the jacobian for problem 12.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision FJAC(N,N), the N by N jacobian matrix.
c
c    Input, integer N, the number of equations.
c
c    Input, double precision X(N), the point at which the jacobian
c    is to be evaluated.
c
      implicit none

      integer n

      double precision fjac(n,n)
      integer i
      integer j
      double precision sum1
      double precision x(n)

      sum1 = 0.0D+00
      do i = 1, n
        sum1 = sum1 + dble ( i ) * ( x(i) - 1.0D+00 )
      end do

      do i = 1, n
        do j = 1, n

          fjac(i,j) = dble ( i * j ) * ( 1.0D+00 + 6.0D+00 * sum1**2 )

          if ( i == j ) then
            fjac(i,j) = fjac(i,j) + 1.0D+00
          end if

        end do
      end do

      return
      end
      subroutine p12_n ( n )

c*********************************************************************72
c
cc P12_N returns the number of equations for problem 12.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer N.
c
c    If N is positive, then N is the only value for the number
c    of equations for this problem.
c
c    If N is negative, then the absolute value of N is the
c    MINIMUM possible value for the number of equations for
c    this problem, and all larger values may also be used.
c
      implicit none

      integer n

      n = - 1

      return
      end
      subroutine p12_sol ( iknow, n, x )

c*********************************************************************72
c
cc P12_SOL returns the solution of problem 12.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IKNOW.
c    -1, there is no solution to the system.
c    0, the solution of the system is not known.
c    positive, IKNOW solutions are known.
c
c    Input, integer N, the order of the system.  N is primarily
c    needed for problems where N may vary.
c
c    Output, double precision X(N), the value of a solution of the system,
c    if known.
c
      implicit none

      integer n

      integer i
      integer iknow
      double precision x(n)

      iknow = 1

      do i = 1, n
        x(i) = 1.0D+00
      end do

      return
      end
      subroutine p12_start ( n, x )

c*********************************************************************72
c
cc P12_START specifies a standard approximate solution for problem 12.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, is the number of equations for the problem.
c
c    Output, double precision X(N), the starting point.
c
      implicit none

      integer n

      integer i
      double precision x(n)

      do i = 1, n
        x(i) = 1.0D+00 - dble ( i ) / dble ( n )
      end do

      return
      end
      subroutine p12_title ( title )

c*********************************************************************72
c
cc P12_TITLE returns the title of problem 12.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'Variably dimensioned function.'
     
      return
      end
      subroutine p13_fx ( n, x, f )

c*********************************************************************72
c
cc P13_FX evaluates the function for problem 13.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of variables.
c
c    Input, double precision X(N), the point where the function
c    is to be evaluated.
c
c    Output, double precision F(N), the value of the function.
c
      implicit none

      integer n

      double precision f(n)
      integer k
      double precision x(n)

      do k = 1, n
     
        f(k) = ( 3.0D+00 - 2.0D+00 * x(k) ) * x(k) + 1.0D+00
     
        if ( 1 .lt. k ) then
          f(k) = f(k) - x(k-1)
        end if
     
        if ( k .lt. n ) then
          f(k) = f(k) - 2.0D+00 * x(k+1)
        end if
     
      end do
     
      return
      end
      subroutine p13_jac ( fjac, n, x )

c*********************************************************************72
c
cc P13_JAC sets the jacobian for problem 13.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision FJAC(N,N), the N by N jacobian matrix.
c
c    Input, integer N, the number of equations.
c
c    Input, double precision X(N), the point at which the jacobian
c    is to be evaluated.
c
      implicit none

      integer n

      double precision fjac(n,n)
      integer i
      integer j
      integer k
      double precision x(n)

      do j = 1, n
        do i = 1, n
          fjac(i,j) = 0.0D+00
        end do
      end do

      do k = 1, n

        fjac(k,k) = 3.0D+00 - 4.0D+00 * x(k)

        if ( 1 .lt. k ) then
          fjac(k,k-1) = - 1.0D+00
        end if

        if ( k .lt. n ) then
          fjac(k,k+1) = - 2.0D+00
        end if

      end do

      return
      end
      subroutine p13_n ( n )

c*********************************************************************72
c
cc P13_N returns the number of equations for problem 13.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer N.
c
c    If N is positive, then N is the only value for the number
c    of equations for this problem.
c
c    If N is negative, then the absolute value of N is the
c    MINIMUM possible value for the number of equations for
c    this problem, and all larger values may also be used.
c
      integer n

      n = - 1

      return
      end
      subroutine p13_sol ( iknow, n, x )

c*********************************************************************72
c
cc P13_SOL returns the solution of problem 13.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IKNOW.
c    -1, there is no solution to the system.
c    0, the solution of the system is not known.
c    positive, IKNOW solutions are known.
c
c    Input, integer N, the order of the system.  N is primarily
c    needed for problems where N may vary.
c
c    Output, double precision X(N), the value of a solution of the system,
c    if known.
c
      implicit none

      integer n

      integer i
      integer iknow
      double precision x(n)

      iknow = 0

      do i = 1, n
        x(i) = 0.0D+00
      end do

      return
      end
      subroutine p13_start ( n, x )

c*********************************************************************72
c
cc P13_START specifies a standard approximate solution for problem 13.
c
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, is the number of equations for the problem.
c
c    Output, double precision X(N), the starting point.
c
      implicit none

      integer n

      integer i
      double precision x(n)

      do i = 1, n
        x(i) = -1.0D+00
      end do

      return
      end
      subroutine p13_title ( title )

c*********************************************************************72
c
cc P13_TITLE returns the title of problem 13.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'Broyden tridiagonal function.'
     
      return
      end
      subroutine p14_fx ( n, x, f )

c*********************************************************************72
c
cc P14_FX evaluates the function for problem 14.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of variables.
c
c    Input, double precision X(N), the point where the function
c    is to be evaluated.
c
c    Output, double precision F(N), the value of the function.
c
      implicit none

      integer n

      double precision f(n)
      integer j
      integer k
      integer k1
      integer k2
      integer ml
      integer mu
      double precision temp
      double precision x(n)

      ml = 5
      mu = 1

      do k = 1, n

        k1 = max ( 1, k - ml )
        k2 = min ( n, k + mu )

        temp = 0.0D+00
        do j = k1, k2
          if ( j .ne. k ) then
            temp = temp + x(j) * ( 1.0D+00 + x(j) )
          end if
        end do

        f(k) = x(k) * ( 2.0D+00 + 5.0D+00 * x(k)**2 ) + 1.0D+00 - temp

      end do
     
      return
      end
      subroutine p14_jac ( fjac, n, x )

c*********************************************************************72
c
cc P14_JAC sets the jacobian for problem 14.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision FJAC(N,N), the N by N jacobian matrix.
c
c    Input, integer N, the number of equations.
c
c    Input, double precision X(N), the point at which the jacobian
c    is to be evaluated.
c
      implicit none

      integer n

      double precision fjac(n,n)
      integer i
      integer j
      integer k
      integer k1
      integer k2
      integer ml
      integer mu
      double precision x(n)

      do j = 1, n
        do i = 1, n
          fjac(i,j) = 0.0D+00
        end do
      end do

      ml = 5
      mu = 1

      do k = 1, n

        k1 = max ( 1, k - ml )
        k2 = min ( n, k + mu )

        do j = k1, k2
          if ( j .ne. k ) then
            fjac(k,j) = - ( 1.0D+00 + 2.0D+00 * x(j) )
          else
            fjac(k,j) = 2.0D+00 + 15.0D+00 * x(k)**2
          end if
        end do

      end do

      return
      end
      subroutine p14_n ( n )

c*********************************************************************72
c
cc P14_N returns the number of equations for problem 14.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer N.
c
c    If N is positive, then N is the only value for the number
c    of equations for this problem.
c
c    If N is negative, then the absolute value of N is the
c    MINIMUM possible value for the number of equations for
c    this problem, and all larger values may also be used.
c
      implicit none

      integer n

      n = - 1
     
      return
      end
      subroutine p14_sol ( iknow, n, x )

c*********************************************************************72
c
cc P14_SOL returns the solution of problem 14.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IKNOW.
c    -1, there is no solution to the system.
c    0, the solution of the system is not known.
c    positive, IKNOW solutions are known.
c
c    Input, integer N, the order of the system.  N is primarily
c    needed for problems where N may vary.
c
c    Output, double precision X(N), the value of a solution of the system,
c    if known.
c
      implicit none

      integer n

      integer i
      integer iknow
      double precision x(n)

      iknow = 0

      do i = 1, n
        x(i) = 0.0D+00
      end do

      return
      end
      subroutine p14_start ( n, x )

c*********************************************************************72
c
cc P14_START specifies a standard approximate solution for problem 14.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, is the number of equations for the problem.
c
c    Output, double precision X(N), the starting point.
c
      implicit none

      integer n

      integer i
      double precision x(n)

      do i = 1, n
        x(i) = -1.0D+00
      end do

      return
      end
      subroutine p14_title ( title )

c*********************************************************************72
c
cc P14_TITLE returns the title of problem 14.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'Broyden banded function.'
     
      return
      end
      subroutine p15_fx ( n, x, f )

c*********************************************************************72
c
cc P15_FX evaluates the function for problem 15.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of variables.
c
c    Input, double precision X(N), the point where the function
c    is to be evaluated.
c
c    Output, double precision F(N), the value of the function.
c
      implicit none

      integer n

      double precision f(n)
      double precision x(n)

      f(1) = ( x(1) * x(1) + x(2) * x(3) ) - 0.0001D+00
      f(2) = ( x(1) * x(2) + x(2) * x(4) ) - 1.0D+00
      f(3) = ( x(3) * x(1) + x(4) * x(3) ) - 0.0D+00
      f(4) = ( x(3) * x(2) + x(4) * x(4) ) - 0.0001D+00
     
      return
      end
      subroutine p15_jac ( fjac, n, x )

c*********************************************************************72
c
cc P15_JAC sets the jacobian for problem 15.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision FJAC(N,N), the N by N jacobian matrix.
c
c    Input, integer N, the number of equations.
c
c    Input, double precision X(N), the point at which the jacobian
c    is to be evaluated.
c
      implicit none

      integer n

      double precision fjac(n,n)
      double precision x(n)

      fjac(1,1) = 2.0D+00 * x(1)
      fjac(1,2) = x(3)
      fjac(1,3) = x(2)
      fjac(1,4) = 0.0D+00

      fjac(2,1) = x(2)
      fjac(2,2) = x(1) + x(4)
      fjac(2,3) = 0.0D+00
      fjac(2,4) = x(2)

      fjac(3,1) = x(3)
      fjac(3,2) = 0.0D+00
      fjac(3,3) = x(1) + x(4)
      fjac(3,4) = x(3)

      fjac(4,1) = 0.0D+00
      fjac(4,2) = x(3)
      fjac(4,3) = x(2)
      fjac(4,4) = 2.0D+00 * x(4)

      return
      end
      subroutine p15_n ( n )

c*********************************************************************72
c
cc P15_N returns the number of equations for problem 15.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer N.
c
c    If N is positive, then N is the only value for the number
c    of equations for this problem.
c
c    If N is negative, then the absolute value of N is the
c    MINIMUM possible value for the number of equations for
c    this problem, and all larger values may also be used.
c
      implicit none

      integer n

      n = 4

      return
      end
      subroutine p15_sol ( iknow, n, x )

c*********************************************************************72
c
cc P15_SOL returns the solution of problem 15.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IKNOW.
c    -1, there is no solution to the system.
c    0, the solution of the system is not known.
c    positive, IKNOW solutions are known.
c
c    Input, integer N, the order of the system.  N is primarily
c    needed for problems where N may vary.
c
c    Output, double precision X(N), the value of a solution of the system,
c    if known.
c
      implicit none

      integer n

      integer iknow
      double precision x(n)

      iknow = 1

      x(1) = 0.01D+00
      x(2) = 50.0D+00
      x(3) = 0.0D+00
      x(4) = 0.01D+00

      return
      end
      subroutine p15_start ( n, x )

c*********************************************************************72
c
cc P15_START specifies a standard approximate solution for problem 15.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, is the number of equations for the problem.
c
c    Output, double precision X(N), the starting point.
c
      implicit none

      integer n

      double precision x(n)

      x(1) = 1.0D+00
      x(2) = 0.0D+00
      x(3) = 0.0D+00
      x(4) = 1.0D+00

      return
      end
      subroutine p15_title ( title )

c*********************************************************************72
c
cc P15_TITLE returns the title of problem 15.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'Hammarling 2 by 2 matrix square root problem.'
     
      return
      end
      subroutine p16_fx ( n, x, f )

c*********************************************************************72
c
cc P16_FX evaluates the function for problem 16.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of variables.
c
c    Input, double precision X(N), the point where the function
c    is to be evaluated.
c
c    Output, double precision F(N), the value of the function.
c
      implicit none

      integer n

      double precision f(n)
      double precision x(n)

      f(1) = ( x(1) * x(1) + x(2) * x(4) + x(3) * x(7) ) - 0.0001D+00
      f(2) = ( x(1) * x(2) + x(2) * x(5) + x(3) * x(8) ) - 1.0D+00
      f(3) = ( x(1) * x(3) + x(2) * x(6) + x(3) * x(9) )

      f(4) = ( x(4) * x(1) + x(5) * x(4) + x(6) * x(7) )
      f(5) = ( x(4) * x(2) + x(5) * x(5) + x(6) * x(8) ) - 0.0001D+00
      f(6) = ( x(4) * x(3) + x(5) * x(6) + x(6) * x(9) )

      f(7) = ( x(7) * x(1) + x(8) * x(4) + x(9) * x(7) )
      f(8) = ( x(7) * x(2) + x(8) * x(5) + x(9) * x(8) )
      f(9) = ( x(7) * x(3) + x(8) * x(6) + x(9) * x(9) ) - 0.0001D+00
     
      return
      end
      subroutine p16_jac ( fjac, n, x )

c*********************************************************************72
c
cc P16_JAC sets the jacobian for problem 16.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision FJAC(N,N), the N by N jacobian matrix.
c
c    Input, integer N, the number of equations.
c
c    Input, double precision X(N), the point at which the jacobian
c    is to be evaluated.
c
      implicit none

      integer n

      double precision fjac(n,n)
      integer i
      integer j
      double precision x(n)

      do j = 1, n
        do i = 1, n
          fjac(i,j) = 0.0D+00
        end do
      end do

      fjac(1,1) = 2.0D+00 * x(1)
      fjac(1,2) = x(4)
      fjac(1,3) = x(7)
      fjac(1,4) = x(2)
      fjac(1,7) = x(3)

      fjac(2,1) = x(2)
      fjac(2,2) = x(1) + x(5)
      fjac(2,3) = x(8)
      fjac(2,5) = x(2)
      fjac(2,8) = x(3)

      fjac(3,1) = x(3)
      fjac(3,2) = x(6)
      fjac(3,3) = x(1) + x(9)
      fjac(3,6) = x(2)
      fjac(3,9) = x(3)

      fjac(4,1) = x(4)
      fjac(4,4) = x(1) + x(5)
      fjac(4,5) = x(4)
      fjac(4,6) = x(7)
      fjac(4,7) = x(6)

      fjac(5,2) = x(4)
      fjac(5,4) = x(2)
      fjac(5,5) = 2.0D+00 * x(5)
      fjac(5,6) = x(8)
      fjac(5,8) = x(6)

      fjac(6,3) = x(4)
      fjac(6,4) = x(3)
      fjac(6,5) = x(6)
      fjac(6,6) = x(5) + x(9)
      fjac(6,9) = x(6)

      fjac(7,1) = x(7)
      fjac(7,4) = x(8)
      fjac(7,7) = x(1) + x(9)
      fjac(7,8) = x(4)
      fjac(7,9) = x(7)

      fjac(8,2) = x(7)
      fjac(8,5) = x(8)
      fjac(8,7) = x(2)
      fjac(8,8) = x(5) + x(9)
      fjac(8,9) = x(8)

      fjac(9,3) = x(7)
      fjac(9,6) = x(8)
      fjac(9,7) = x(3)
      fjac(9,8) = x(6)
      fjac(9,9) = 2.0D+00 * x(9)

      return
      end
      subroutine p16_n ( n )

c*********************************************************************72
c
cc P16_N returns the number of equations for problem 16.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer N.
c
c    If N is positive, then N is the only value for the number
c    of equations for this problem.
c
c    If N is negative, then the absolute value of N is the
c    MINIMUM possible value for the number of equations for
c    this problem, and all larger values may also be used.
c
      implicit none

      integer n

      n = 9

      return
      end
      subroutine p16_sol ( iknow, n, x )

c*********************************************************************72
c
cc P16_SOL returns the solution of problem 16.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IKNOW.
c    -1, there is no solution to the system.
c    0, the solution of the system is not known.
c    positive, IKNOW solutions are known.
c
c    Input, integer N, the order of the system.  N is primarily
c    needed for problems where N may vary.
c
c    Output, double precision X(N), the value of a solution of the system,
c    if known.
c
      implicit none

      integer n

      integer iknow
      double precision x(n)

      iknow = 1

      x(1) = 0.01D+00
      x(2) = 50.0D+00
      x(3) = 0.0D+00

      x(4) = 0.0D+00
      x(5) = 0.01D+00
      x(6) = 0.0D+00

      x(7) = 0.0D+00
      x(8) = 0.0D+00
      x(9) = 0.01D+00

      return
      end
      subroutine p16_start ( n, x )

c*********************************************************************72
c
cc P16_START specifies a standard approximate solution for problem 16.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, is the number of equations for the problem.
c
c    Output, double precision X(N), the starting point.
c
      implicit none

      integer n

      double precision x(n)

      x(1) = 1.0D+00
      x(2) = 0.0D+00
      x(3) = 0.0D+00
      x(4) = 0.0D+00
      x(5) = 1.0D+00
      x(6) = 0.0D+00
      x(7) = 0.0D+00
      x(8) = 0.0D+00
      x(9) = 1.0D+00

      return
      end
      subroutine p16_title ( title )

c*********************************************************************72
c
cc P16_TITLE returns the title of problem 16.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'Hammarling 3 by 3 matrix square root problem.'
     
      return
      end
      subroutine p17_fx ( n, x, f )

c*********************************************************************72
c
cc P17_FX evaluates the function for problem 17.
c
c  Discussion:
c
c    The equations are:
c
c      F1(X) = X(1) + X(2) - 3
c      F2(X) = X(1)**2 + X(2)**2 - 9
c
c    with roots (3,0) and (0,3).
c
c    Using a starting point of (1,5), here are the iterates for Broyden's
c    method and Newton's method:
c
c    Broyden's Method
c 
c    0   1.0              5.0
c    1  -0.625            3.625
c    2  -0.0757575757575  3.0757575757575
c    3  -0.0127942681679  3.0127942681679
c    4  -0.0003138243387  3.0003138243387
c    5  -0.0000013325618  3.0000013325618
c    6  -0.0000000001394  3.0000000001394
c    7   0.0              3.0
c
c    Newton's Method
c 
c    0   1.0              5.0
c    1  -0.625            3.625
c    2  -0.0919117647059  3.0919117647059
c    3  -0.0026533419372  3.0026533419372
c    4  -0.0000023425973  3.0000023425973
c    5  -0.0000000000018  3.0000000000018
c    6   0.0              3.0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    John Dennis, Robert Schnabel,
c    Numerical Methods for Unconstrained Optimization 
c    and Nonlinear Equations,
c    SIAM, 1996,
c    ISBN13: 978-0-898713-64-0,
c    LC: QA402.5.D44.
c
c  Parameters:
c
c    Input, integer N, the number of variables.
c
c    Input, double precision X(N), the point where the function 
c    is to be evaluated.
c
c    Output, double precision F(N), the value of the function.
c
      implicit none

      integer n

      double precision f(n)
      double precision x(n)

      f(1) = x(1) + x(2) - 3.0D+00
      f(2) = x(1)**2 + x(2)**2 - 9.0D+00
     
      return
      end
      subroutine p17_jac ( fjac, n, x )

c*********************************************************************72
c
cc P17_JAC sets the jacobian for problem 17.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision FJAC(N,N), the N by N jacobian matrix.
c
c    Input, integer N, the number of equations.
c
c    Input, double precision X(N), the point at which the jacobian
c    is to be evaluated.
c
      implicit none

      integer n

      double precision fjac(n,n)
      double precision x(n)

      fjac(1,1) = 1.0D+00
      fjac(1,2) = 1.0D+00

      fjac(2,1) = 2.0D+00 * x(1)
      fjac(2,2) = 2.0D+00 * x(2)

      return
      end
      subroutine p17_n ( n )

c*********************************************************************72
c
cc P17_N returns the number of equations for problem 17.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer N.
c
c    If N is positive, then N is the only value for the number
c    of equations for this problem.
c
c    If N is negative, then the absolute value of N is the
c    MINIMUM possible value for the number of equations for
c    this problem, and all larger values may also be used.
c
      implicit none

      integer n

      n = 2
     
      return
      end
      subroutine p17_sol ( iknow, n, x )

c*********************************************************************72
c
cc P17_SOL returns the solution of problem 17.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IKNOW.
c    -1, there is no solution to the system.
c    0, the solution of the system is not known.
c    positive, IKNOW solutions are known.
c
c    Input, integer N, the order of the system.  N is primarily
c    needed for problems where N may vary.
c
c    Output, double precision X(N), the value of a solution of the system,
c    if known.
c
      implicit none

      integer n

      integer iknow
      double precision x(n)

      iknow = 1

      x(1) = 0.0D+00
      x(2) = 3.0D+00

      return
      end
      subroutine p17_start ( n, x )

c*********************************************************************72
c
cc P17_START specifies a standard approximate solution for problem 17.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, is the number of equations for the problem.
c
c    Output, double precision X(N), the starting point.
c
      implicit none

      integer n

      double precision x(n)

      x(1) = 1.0D+00
      x(2) = 5.0D+00

      return
      end
      subroutine p17_title ( title )

c*********************************************************************72
c
cc P17_TITLE returns the title of problem 17.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'Dennis and Schnabel 2 by 2 example.'
     
      return
      end
      subroutine p18_fx ( n, x, f )

c*********************************************************************72
c
cc P18_FX evaluates the function for problem 18.
c
c  Discussion:
c
c    This problem has as roots any point (x,y) with x or y equal to
c    zero.  The jacobian is of rank one at all roots, except at the
c    origin, where it completely vanishes.  
c
c    Newton iterates can converge to the origin, or to points on
c    the X or Y axes, or can even "converge" to a point at infinity,
c    generating a sequence of ever larger points with ever smaller
c    function values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of variables.
c
c    Input, double precision X(N), the point where the function
c    is to be evaluated.
c
c    Output, double precision F(N), the value of the function.
c
      implicit none

      integer n

      double precision f(n)
      double precision x(n)

      if ( x(1) .ne. 0.0D+00 ) then
        f(1) = x(2)**2 * ( 1.0D+00 - exp ( - x(1)**2 ) ) / x(1)
      else
        f(1) = 0.0D+00
      end if
     
      if ( x(2) .ne. 0.0D+00 ) then
        f(2) = x(1) * ( 1.0D+00 - exp ( - x(2)**2 ) ) / x(2)
      else
        f(2) = 0.0D+00
      end if
     
      return
      end
      subroutine p18_jac ( fjac, n, x )

c*********************************************************************72
c
cc P18_JAC sets the jacobian for problem 18.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision FJAC(N,N), the N by N jacobian matrix.
c
c    Input, integer N, the number of equations.
c
c    Input, double precision X(N), the point at which the jacobian
c    is to be evaluated.
c
      implicit none

      integer n

      double precision fjac(n,n)
      double precision x(n)

      if ( x(1) .ne. 0.0D+00 ) then
        fjac(1,1) = x(2)**2 * ( 2.0D+00 * exp ( - x(1)**2 ) - 
     &    ( 1.0D+00 - exp ( - x(1)**2 ) ) / x(1)**2)

        fjac(1,2) = 2.0D+00 * x(2) * ( 1.0D+00 - exp ( - x(1)**2 ) ) 
     &    / x(1)

      else

        fjac(1,1) = x(2)**2

        fjac(1,2) = 0.0D+00

      end if

      if ( x(2) .ne. 0.0D+00 ) then

        fjac(2,1) = ( 1.0D+00 - exp ( - x(2)**2 ) ) / x(2)

        fjac(2,2) = x(1) * ( 2.0D+00 * exp ( - x(2)**2 ) - 
     &    ( 1.0D+00 - exp ( - x(2)**2 ) ) / x(2)**2 )

      else

        fjac(2,1) = 0.0D+00

        fjac(2,2) = x(1)

      end if

      return
      end
      subroutine p18_n ( n )

c*********************************************************************72
c
cc P18_N returns the number of equations for problem 18.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer N.
c
c    If N is positive, then N is the only value for the number
c    of equations for this problem.
c
c    If N is negative, then the absolute value of N is the
c    MINIMUM possible value for the number of equations for
c    this problem, and all larger values may also be used.
c
      implicit none

      integer n

      n = 2

      return
      end
      subroutine p18_sol ( iknow, n, x )

c*********************************************************************72
c
cc P18_SOL returns the solution of problem 18.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IKNOW.
c    -1, there is no solution to the system.
c    0, the solution of the system is not known.
c    positive, IKNOW solutions are known.
c
c    Input, integer N, the order of the system.  N is primarily
c    needed for problems where N may vary.
c
c    Output, double precision X(N), the value of a solution of the system,
c    if known.
c
      implicit none

      integer n

      integer iknow
      double precision x(n)

      iknow = 1

      x(1) = 0.0D+00
      x(2) = 0.0D+00

      return
      end
      subroutine p18_start ( n, x )

c*********************************************************************72
c
cc P18_START specifies a standard approximate solution for problem 18.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, is the number of equations for the problem.
c
c    Output, double precision X(N), the starting point.
c
      implicit none

      integer n

      double precision x(n)

      x(1) = 2.0D+00
      x(2) = 2.0D+00

      return
      end
      subroutine p18_title ( title )

c*********************************************************************72
c
cc P18_TITLE returns the title of problem 18.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'Sample problem 18.'
     
      return
      end
      subroutine p19_fx ( n, x, f )

c*********************************************************************72
c
cc P19_FX evaluates the function for problem 19.
c
c  Discussion:
c
c    This problem has a single root at the origin.  Convergence of the
c    Newton iterates should be monotonic, but only linear in rate,
c    since the jacobian is singular at the origin.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of variables.
c
c    Input, double precision X(N), the point where the function 
c    is to be evaluated.
c
c    Output, double precision F(N), the value of the function.
c
      implicit none

      integer n

      double precision f(n)
      double precision x(n)

      f(1) = x(1) * ( x(1)**2 + x(2)**2 )
      f(2) = x(2) * ( x(1)**2 + x(2)**2 )
     
      return
      end
      subroutine p19_jac ( fjac, n, x )

c*********************************************************************72
c
cc P19_JAC sets the jacobian for problem 19.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision FJAC(N,N), the N by N jacobian matrix.
c
c    Input, integer N, the number of equations.
c
c    Input, double precision X(N), the point at which the jacobian
c    is to be evaluated.
c
      implicit none

      integer n

      double precision fjac(n,n)
      double precision x(n)

      fjac(1,1) = 3.0D+00 * x(1)**2 + x(2)**2
      fjac(1,2) = 2.0D+00 * x(1) * x(2)

      fjac(2,1) = 2.0D+00 * x(1) * x(2)
      fjac(2,2) = x(1)**2 + 3.0D+00 * x(2)**2

      return
      end
      subroutine p19_n ( n )

c*********************************************************************72
c
cc P19_N returns the number of equations for problem 19.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer N.
c
c    If N is positive, then N is the only value for the number
c    of equations for this problem.
c
c    If N is negative, then the absolute value of N is the
c    MINIMUM possible value for the number of equations for
c    this problem, and all larger values may also be used.
c
      implicit none

      integer n

      n = 2
     
      return
      end
      subroutine p19_sol ( iknow, n, x )

c*********************************************************************72
c
cc P19_SOL returns the solution of problem 19.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IKNOW.
c    -1, there is no solution to the system.
c    0, the solution of the system is not known.
c    positive, IKNOW solutions are known.
c
c    Input, integer N, the order of the system.  N is primarily
c    needed for problems where N may vary.
c
c    Output, double precision X(N), the value of a solution of the system,
c    if known.
c
      implicit none

      integer n

      integer iknow
      double precision x(n)

      iknow = 1

      x(1) = 0.0D+00
      x(2) = 0.0D+00

      return
      end
      subroutine p19_start ( n, x )

c*********************************************************************72
c
cc P19_START specifies a standard approximate solution for problem 19.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, is the number of equations for the problem.
c
c    Output, double precision X(N), the starting point.
c
      implicit none

      integer n

      double precision x(n)

      x(1:2) = (/ 3.0D+00, 3.0D+00 /)

      return
      end
      subroutine p19_title ( title )

c*********************************************************************72
c
cc P19_TITLE returns the title of problem 19.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'Sample problem 19.'
     
      return
      end
      subroutine p20_fx ( n, x, f )

c*********************************************************************72
c
cc P20_FX evaluates the function for problem 20.
c
c  Discussion:
c
c    This problem has a single root at the origin, and a multiple root
c    at x = 5.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of variables.
c
c    Input, double precision X(N), the point where the function
c    is to be evaluated.
c
c    Output, double precision F(N), the value of the function.
c
      implicit none

      integer n

      double precision f(n)
      double precision x(n)

      f(1) = x(1) * ( x(1) - 5.0D+00 )**2
     
      return
      end
      subroutine p20_jac ( fjac, n, x )

c*********************************************************************72
c
cc P20_JAC sets the jacobian for problem 20.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision FJAC(N,N), the N by N jacobian matrix.
c
c    Input, integer N, the number of equations.
c
c    Input, double precision X(N), the point at which the jacobian
c    is to be evaluated.
c
      implicit none

      integer n

      double precision fjac(n,n)
      double precision x(n)

      fjac(1,1) = ( 3.0D+00 * x(1) - 5.0D+00 ) * ( x(1) - 5.0D+00 )

      return
      end
      subroutine p20_n ( n )

c*********************************************************************72
c
cc P20_N returns the number of equations for problem 20.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer N.
c
c    If N is positive, then N is the only value for the number
c    of equations for this problem.
c
c    If N is negative, then the absolute value of N is the
c    MINIMUM possible value for the number of equations for
c    this problem, and all larger values may also be used.
c
      implicit none

      integer n

      n = 1
     
      return
      end
      subroutine p20_sol ( iknow, n, x )

c*********************************************************************72
c
cc P20_SOL returns the solution of problem 20.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IKNOW.
c    -1, there is no solution to the system.
c    0, the solution of the system is not known.
c    positive, IKNOW solutions are known.
c
c    Input, integer N, the order of the system.  N is primarily
c    needed for problems where N may vary.
c
c    Output, double precision X(N), the value of a solution of the system,
c    if known.
c
      implicit none

      integer n

      integer, save :: icall = -1
      integer iknow
      double precision x(n)

      icall = icall + 1
      icall = mod ( icall, iknow )

      iknow = 2

      if ( icall == 0 ) then
        x(1) = 0.0D+00
      else if ( icall == 1 ) then
        x(1) = 5.0D+00
      end if

      return
      end
      subroutine p20_start ( n, x )

c*********************************************************************72
c
cc P20_START specifies a standard approximate solution for problem 20.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, is the number of equations for the problem.
c
c    Output, double precision X(N), the starting point.
c
      implicit none

      integer n

      double precision x(n)

      x(1) = 1.0D+00

      return
      end
      subroutine p20_title ( title )

c*********************************************************************72
c
cc P20_TITLE returns the title of problem 20.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'Scalar problem f(x) = x * ( x - 5 )**2.'
     
      return
      end
      subroutine p21_fx ( n, x, f )

c*********************************************************************72
c
cc P21_FX evaluates the function for problem 21.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of variables.
c
c    Input, double precision X(N), the point where the function
c    is to be evaluated.
c
c    Output, double precision F(N), the value of the function.
c
      implicit none

      integer n

      double precision f(n)
      double precision x(n)

      f(1) = x(1) - x(2)**3 + 5.0D+00 * x(2)**2 -  2.0D+00 * x(2) - 13.0D+00
      f(2) = x(1) + x(2)**3 +           x(2)**2 - 14.0D+00 * x(2) - 29.0D+00
     
      return
      end
      subroutine p21_jac ( fjac, n, x )

c*********************************************************************72
c
cc P21_JAC sets the jacobian for problem 21
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision FJAC(N,N), the N by N jacobian matrix.
c
c    Input, integer N, the number of equations.
c
c    Input, double precision X(N), the point at which the jacobian
c    is to be evaluated.
c
      implicit none

      integer n

      double precision fjac(n,n)
      double precision x(n)

      fjac(1,1) = 1.0D+00
      fjac(1,2) = - 3.0D+00 * x(2)**2 + 10.0D+00 * x(2) -  2.0D+00

      fjac(2,1) = 1.0D+00
      fjac(2,2) =   3.0D+00 * x(2)**2 +  2.0D+00 * x(2) - 14.0D+00

      return
      end
      subroutine p21_n ( n )

c*********************************************************************72
c
cc P21_N returns the number of equations for problem 21.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer N.
c
c    If N is positive, then N is the only value for the number
c    of equations for this problem.
c
c    If N is negative, then the absolute value of N is the
c    MINIMUM possible value for the number of equations for
c    this problem, and all larger values may also be used.
c
      implicit none

      integer n

      n = 2

      return
      end
      subroutine p21_sol ( iknow, n, x )

c*********************************************************************72
c
cc P21_SOL returns the solution of problem 21.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IKNOW.
c    -1, there is no solution to the system.
c    0, the solution of the system is not known.
c    positive, IKNOW solutions are known.
c
c    Input, integer N, the order of the system.  N is primarily
c    needed for problems where N may vary.
c
c    Output, double precision X(N), the value of a solution of the system,
c    if known.
c
      implicit none

      integer n

      integer iknow
      double precision x(n)

      iknow = 1

      x(1) = 5.0D+00
      x(2) = 4.0D+00

      return
      end
      subroutine p21_start ( n, x )

c*********************************************************************72
c
cc P21_START specifies a standard approximate solution for problem 21.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, is the number of equations for the problem.
c
c    Output, double precision X(N), the starting point.
c
      implicit none

      integer n

      double precision x(n)

      x(1) =   0.5D+00
      x(2) = - 2.0D+00

      return
      end
      subroutine p21_title ( title )

c*********************************************************************72
c
cc P21_TITLE returns the title of problem 21.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'Freudenstein-Roth function.'
     
      return
      end
      subroutine p22_fx ( n, x, f )

c*********************************************************************72
c
cc P22_FX evaluates the function for problem 22.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 December 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of variables.
c
c    Input, double precision X(N), the point where the function
c    is to be evaluated.
c
c    Output, double precision F(N), the value of the function.
c
      implicit none

      integer n

      double precision f(n)
      double precision :: pi = 3.141592653589793D+00
      double precision x(n)

      f(1) = x(1)**2 - x(2) + 1.0D+00
      f(2) = x(1) - cos ( 0.5D+00 * pi * x(2) )
     
      return
      end
      subroutine p22_jac ( fjac, n, x )

c*********************************************************************72
c
cc P22_JAC sets the jacobian for problem 22.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 June 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision FJAC(N,N), the N by N jacobian matrix.
c
c    Input, integer N, the number of equations.
c
c    Input, double precision X(N), the point at which the jacobian
c    is to be evaluated.
c
      implicit none

      integer n

      double precision fjac(n,n)
      double precision :: pi = 3.141592653589793D+00
      double precision x(n)

      fjac(1,1) = 2.0D+00 * x(1)
      fjac(1,2) = - 1.0D+00

      fjac(2,1) = 1.0D+00
      fjac(2,2) = 0.5D+00 * pi * sin ( 0.5D+00 * pi * x(2) )

      return
      end
      subroutine p22_n ( n )

c*********************************************************************72
c
cc P22_N returns the number of equations for problem 22.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 January 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer N.
c
c    If N is positive, then N is the only value for the number
c    of equations for this problem.
c
c    If N is negative, then the absolute value of N is the
c    MINIMUM possible value for the number of equations for
c    this problem, and all larger values may also be used.
c
      implicit none

      integer n

      n = 2
     
      return
      end
      subroutine p22_sol ( iknow, n, x )

c*********************************************************************72
c
cc P22_SOL returns the solution of problem 22.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IKNOW.
c    -1, there is no solution to the system.
c    0, the solution of the system is not known.
c    positive, IKNOW solutions are known.
c
c    Input, integer N, the order of the system.  N is primarily
c    needed for problems where N may vary.
c
c    Output, double precision X(N), the value of a solution of the system,
c    if known.
c
      implicit none

      integer n

      integer iknow
      double precision x(n)

      iknow = 1

      x(1) = 0.0D+00
      x(2) = 1.0D+00

      return
      end
      subroutine p22_start ( n, x )

c*********************************************************************72
c
cc P22_START specifies a standard approximate solution for problem 22.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, is the number of equations for the problem.
c
c    Output, double precision X(N), the starting point.
c
      implicit none

      integer n

      double precision x(n)

      x(1:2) = (/ 1.0D+00, 0.0D+00 /)

      return
      end
      subroutine p22_title ( title )

c*********************************************************************72
c
cc P22_TITLE returns the title of problem 22.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'Boggs function.'

      return
      end
      subroutine p23_fx ( n, x, f )

c*********************************************************************72
c
cc P23_FX evaluates the function for problem 23.
c
c    The Chandrasekhar H function is
c
c      F(H)(mu) = H(mu) 
c        - 1 / ( 1 - c/2 
c        * Integral ( 0 <= nu <= 1 ) [ mu * H(mu) / ( mu + nu ) ] dnu )
c
c    Applying the composite midpoint rule, and setting, for i = 1 to N,
c
c      mu(i) = ( i - 0.5 ) / N
c
c    we have the discrete function
c
c      F(h)(i) = h(i) 
c        - 1 / ( 1 - c/(2N) 
c        * sum ( 1 <= j <= N ) ( mu(i) * h(j) / ( mu(i) + mu(j) ) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 June 2007
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Subramanyan Chandrasekhar,
c    Radiative Transfer,
c    Dover, 1960,
c    ISBN13: 978-0486605906,
c    LC: QB461.C46.
c
c  Parameters:
c
c    Input, integer N, the number of variables.
c
c    Input, double precision X(N), the point where the function
c    is to be evaluated.
c
c    Output, double precision F(N), the value of the function.
c
c  Local Parameters:
c
c    Local, double precision C, the constant, which must be between
c    0 and 1.
c
      implicit none

      integer n

      double precision, parameter :: c = 0.9D+00
      double precision f(n)
      integer i
      integer j
      double precision mu(n)
      double precision temp
      double precision term
      double precision x(n)

      do i = 1, n
        f(i) = x(i)
      end do

      do i = 1, n
        mu(i) = dble ( 2 * i - 1 ) / dble ( 2 * n )
      end do

      do i = 1, n

        temp = 0.0D+00
        do j = 1, n
          temp = temp + mu(i) * x(j) / ( mu(i) + mu(j) )
        end do

        term = 1.0D+00 - c * temp / dble ( 2 * n )

        f(i) = f(i) - 1.0D+00 / term

      end do
     
      return
      end
      subroutine p23_n ( n )

c*********************************************************************72
c
cc P23_N returns the number of equations for problem 23.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 June 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer N.
c    If N is positive, then N is the only value for the number
c    of equations for this problem.
c    If N is negative, then the absolute value of N is the
c    MINIMUM possible value for the number of equations for
c    this problem, and all larger values may also be used.
c
      implicit none

      integer n

      n = -1

      return
      end
      subroutine p23_jac ( fjac, n, x )

c*********************************************************************72
c
cc P23_JAC sets the jacobian for problem 23.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    25 June 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision FJAC(N,N), the N by N jacobian matrix.
c
c    Input, integer N, the number of equations.
c
c    Input, double precision X(N), the point where the jacobian 
c    is to be evaluated.
c
      implicit none

      integer n

      double precision, parameter :: c = 0.9D+00
      double precision fjac(n,n)
      integer i
      integer j
      double precision mu(n)
      double precision temp
      double precision term
      double precision termp
      double precision x(n)

      do j = 1, n
        do i = 1, n
          fjac(i,j) = 0.0D+00
        end do
      end do

      do i = 1, n
        mu(i) = dble ( 2 * i - 1 ) / dble ( 2 * n )
      end do

      do i = 1, n

        temp = 0.0
        do j = 1, n
          temp = temp + mu(i) * x(j) / ( mu(i) + mu(j) )
        end do

        term = 1.0D+00 - c * temp / dble ( 2 * n )

        do j = 1, n

          termp = c * mu(i) / ( mu(i) + mu(j) ) / dble ( 2 * n )

          fjac(i,j) = termp / term**2

        end do

      end do

      return
      end
      subroutine p23_sol ( iknow, n, x )

c*********************************************************************72
c
cc P23_SOL returns the solution of problem 23.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 June 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IKNOW.
c    -1, there is no solution to the system.
c    0, the solution of the system is not known.
c    positive, IKNOW solutions are known.
c
c    Input, integer N, the order of the system.  N is primarily
c    needed for problems where N may vary.
c
c    Output, double precision X(N), the value of a solution of the system,
c    if known.
c
      implicit none

      integer n

      integer i
      integer iknow
      double precision x(n)

      iknow = -1

      do i = 1, n
        x(i) = 0.0D+00
      end do

      return
      end
      subroutine p23_start ( n, x )

c*********************************************************************72
c
cc P23_START specifies a standard approximate solution for problem 23.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 June 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, is the number of equations for the 
c    problem.
c
c    Output, double precision X(N), the starting point.
c
      implicit none

      integer n

      integer i
      double precision x(n)

      do i = 1, n
        x(i) = 1.0D+00
      end do

      return
      end
      subroutine p23_title ( title )

c*********************************************************************72
c
cc P23_TITLE returns the title of problem 23.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 June 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character ( len = * ) TITLE, the title of the problem.
c
      implicit none

      character ( len = * ) title

      title = 'Chandrasekhar function.'
     
      return
      end
      subroutine r8_swap ( x, y )

c*********************************************************************72
c
cc R8_SWAP switches two R8's.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 May 2000
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, double precision X, Y.  On output, the values of X and
c    Y have been interchanged.
c
      implicit none

      double precision x
      double precision y
      double precision z

      z = x
      x = y
      y = z

      return
      end
      function r8vec_norm2 ( n, a )

c*********************************************************************72
c
cc R8VEC_NORM2 returns the 2-norm of a vector.
c
c  Definition:
c
c    The vector 2-norm is defined as:
c
c      value = sqrt ( sum ( 1 <= I <= N ) A(I)**2 ).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    29 November 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries in A.
c
c    Input, double precision A(N), the vector whose 2-norm is desired.
c
c    Output, double precision R8VEC_NORM2, the 2-norm of A.
c
      implicit none

      integer n

      double precision a(n)
      integer i
      double precision r8vec_norm2

      r8vec_norm2 = 0.0D+00
      do i = 1, n
        r8vec_norm2 = r8vec_norm2 + a(i)**2
      end do

      r8vec_norm2 = sqrt ( r8vec_norm2 )

      return
      end
      subroutine r8ge_fa ( n, a, pivot, info )

c*********************************************************************72
c
cc R8GE_FA factors a general matrix.
c
c  Discussion:
c
c    R8GE_FA is a simplified version of the LINPACK routine R8GEFA.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 February 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input/output, double precision A(N,N), the matrix to be factored.
c    On output, A contains an upper triangular matrix and the multipliers
c    which were used to obtain it.  The factorization can be written
c    A = L * U, where L is a product of permutation and unit lower
c    triangular matrices and U is upper triangular.
c
c    Output, integer PIVOT(N), a vector of pivot indices.
c
c    Output, integer INFO, singularity flag.
c    0, no singularity detected.
c    nonzero, the factorization failed on the INFO-th step.
c
      implicit none

      integer n

      double precision a(n,n)
      integer i
      integer info
      integer pivot(n)
      integer j
      integer k
      integer l

      info = 0

      do k = 1, n-1
c
c  Find L, the index of the pivot row.
c
        l = k
        do i = k+1, n
          if ( abs ( a(l,k) ) .lt. abs ( a(i,k) ) ) then
            l = i
          end if
        end do

        pivot(k) = l
c
c  If the pivot index is zero, the algorithm has failed.
c
        if ( a(l,k) .eq. 0.0D+00 ) then
          info = k
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'R8GE_FA - Fatal error!'
          write ( *, '(a,i6)' ) '  Zero pivot on step ', info
          return
        end if
c
c  Interchange rows L and K if necessary.
c
        if ( l .ne. k ) then
          call r8_swap ( a(l,k), a(k,k) )
        end if
c
c  Normalize the values that lie below the pivot entry A(K,K).
c
        do i = k + 1, n
          a(i,k) = -a(i,k) / a(k,k)
        end do
c
c  Row elimination with column indexing.
c
        do j = k+1, n

          if ( l .ne. k ) then
            call r8_swap ( a(l,j), a(k,j) )
          end if

          do i = k + 1, n
            a(i,j) = a(i,j) + a(i,k) * a(k,j)
          end do

        end do

      end do

      pivot(n) = n

      if ( a(n,n) .eq. 0.0D+00 ) then
        info = n
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8GE_FA - Fatal error!'
        write ( *, '(a,i6)' ) '  Zero pivot on step ', info
      end if

      return
      end
      subroutine r8ge_sl ( n, a, pivot, b, job )

c*********************************************************************72
c
cc R8GE_SL solves a system factored by SGE_FA.
c
c  Discussion:
c
c    R8GE_SL is a simplified version of the LINPACK routine R8GESL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 March 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the matrix.
c    N must be positive.
c
c    Input, double precision A(N,N), the LU factors from R8GE_FA.
c
c    Input, integer PIVOT(N), the pivot vector from R8GE_FA.
c
c    Input/output, double precision B(N).
c    On input, the right hand side vector.
c    On output, the solution vector.
c
c    Input, integer JOB, specifies the operation.
c    0, solve A * x = b.
c    nonzero, solve transpose ( A ) * x = b.
c
      implicit none

      integer n

      double precision a(n,n)
      double precision b(n)
      integer i
      integer job
      integer k
      integer l
      integer pivot(n)
      double precision t
c
c  Solve A * x = b.
c
      if ( job == 0 ) then
c
c  Solve PL * Y = B.
c
        do k = 1, n-1

          l = pivot(k)

          if ( l .ne. k ) then
            call r8_swap ( b(l), b(k) )
          end if

          do i = k + 1, n
            b(i) = b(i) + a(i,k) * b(k)
          end do

        end do
c
c  Solve U * X = Y.
c
        do k = n, 1, -1
          b(k) = b(k) / a(k,k)
          do i = 1, k - 1
            b(i) = b(i) - a(i,k) * b(k)
          end do
        end do
c
c  Solve transpose ( A ) * X = B.
c
      else
c
c  Solve transpose ( U ) * Y = B.
c
        do k = 1, n
          t = 0.0D+00
          do i = 1, k - 1
            t = t + b(i) * a(i,k)
          end do
          b(k) = ( b(k) - t ) / a(k,k)
        end do
c
c  Solve transpose ( PL ) * X = Y.
c
        do k = n-1, 1, -1

          do i = k + 1, n
            b(k) = b(k) + b(i) * a(i,k)
          end do

          l = pivot(k)

          if ( l .ne. k ) then
            call r8_swap ( b(l), b(k) )
          end if

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
