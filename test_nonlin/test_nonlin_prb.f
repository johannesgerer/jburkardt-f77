      program main

c*********************************************************************72
c
cc MAIN is the main program for TEST_NONLIN_PRB.
c
c  Discussion:
c
c    TEST_NONLIN_PRB demonstrates the use of the test functions in TEST_NONLIN.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 May 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_NONLIN_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TEST_NONLIN library.'

      call test01 ( )
      call test02 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_NONLIN_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 prints the title of each problem.
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
      implicit none

      integer problem
      integer problem_num
      character * ( 80 ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Print the title of each problem.'
      write ( *, '(a)' ) ' '
c
c  Get the number of problems.
c
      call p00_problem_num ( problem_num )

      do problem = 1, problem_num
        call p00_title ( problem, title )
        write ( *, '(2x,i2,2x,a)' )
     &    problem, '"' // trim ( title ) // '"'
      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 uses a simple Newton method on all the problems.
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
      implicit none

      integer n_max
      parameter ( n_max = 20 )

      double precision abstol
      logical converged
      double precision damp
      double precision damptol
      double precision f(n_max)
      double precision fnrm
      double precision fnrm_old
      double precision fnrm_zero
      double precision fprime(n_max,n_max)
      integer i
      integer ijac
      integer iknow
      integer info
      integer it
      integer maxit
      integer n
      integer pivot(n_max)
      integer problem
      integer problem_num
      double precision reltol
      double precision r8vec_norm2
      double precision snrm
      double precision step(n_max)
      character * ( 80 ) title
      double precision x(n_max)
      double precision xnew(n_max)
      double precision xnrm
      double precision xold

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Seek roots of the nonlinear functions'
      write ( *, '(a)' ) '  in TEST_NONLIN using Newton''s method.'
      write ( *, '(a)' ) ' '
c
c  Initialization.
c
      abstol = 1.0D-06
      damptol = 1.0D-06
      reltol = 1.0D-06
      maxit = 40
c
c  Get the number of problems.
c
      call p00_problem_num ( problem_num )
c
c  Solve each problem.
c
      do problem = 1, problem_num
c
c  Get the default problem size.
c  If none, use N = 10.
c
        call p00_n ( problem, n )

        if ( n .lt. 0 ) then
          n = 10
        end if
c
c  Get the problem title.
c
        call p00_title ( problem, title )

        do ijac = 0, 1

          write ( *, '(a)' ) ' '
          write ( *, '(a,i6)' ) '  Problem index ', problem
          write ( *, '(3x,a)' ) trim ( title )
          write ( *, '(a,i6)' ) '  Problem size = ', n
          if ( ijac .eq. 0 ) then
            write ( *, '(a)' ) '  Use analytic jacobian'
          else
            write ( *, '(a)' ) '  Use finite difference jacobian'
          end if
c
c  Get the problem starting point.
c
          call p00_start ( problem, n, x )

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Initial X:'
          write ( *, '(a)' ) ' '
          write ( *, '(5g14.6)' ) ( x(i), i = 1, n )
          xnrm = r8vec_norm2 ( n, x )
c
c  Evaluate function at starting point.
c
          call p00_fx ( problem, f, n, x )

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Function value at initial X:'
          write ( *, '(a)' ) ' '
          write ( *, '(5g14.6)' ) ( f(i), i = 1, n )
          fnrm = r8vec_norm2 ( n, f )
          fnrm_old = fnrm
          fnrm_zero = fnrm

          write ( *, '(a)' ) ' '
          write ( *, '(a)' )
     &      'Iteration   ||X||         ||F(X)||   Damping   ||dX||'
          write ( *, '(a)' ) ' '
          it = 0
          write ( *, '(7x,i3,g14.6,g14.6)' ) it, xnrm, fnrm

          converged = .false.

          do it = 1, maxit
c
c  Get the Jacobian matrix.
c
            if ( ijac == 0 ) then
              call p00_jac ( problem, n, fprime, x )
            else
              call p00_dif ( problem, n, fprime, x )
            end if
c
c  Factor the jacobian matrix.
c
            call r8ge_fa ( n, fprime, pivot, info )

            if ( info .ne. 0 ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'TEST02 - Warning:'
              write ( *, '(a)' ) '  The iteration must be halted.'
              write ( *, '(a)' ) '  The jacobian matrix is singularc'
              exit
            end if
c
c  Solve Fprime * DeltaX = F.
c
            do i = 1, n
              step(i) = f(i)
            end do

            call r8ge_sl ( n, fprime, pivot, step, 0 )

            snrm = r8vec_norm2 ( n, step )

            if ( 10.0D+00 * ( xnrm + 10.0D+00 ) .lt. snrm .and.
     &        1 .lt. it ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) 'TEST02'
              write ( *, '(a,g14.6)' )
     &          '  Excessive correction norm = ', snrm
              write ( *, '(a,g14.6)' )
     &          '  Solution norm was           ', xnrm
              write ( *, '(a)' ) '  The iteration is terminated.'
              exit
            end if
c
c  Update XNEW := X - DeltaX.
c
            if ( snrm .lt. 2.0D+00 * ( xnrm + 1.0D+00 ) ) then
              damp = 1.0D+00
            else
              damp = 2.0D+00 * ( xnrm + 1.0D+00 ) / snrm
            end if

            damp = max ( damp, damptol )

10          continue

            do i = 1, n
              xnew(i) = x(i) - damp * step(i)
            end do

            xold = xnrm
            xnrm = r8vec_norm2 ( n, xnew )
            call p00_fx ( problem, f, n, xnew )
            fnrm = r8vec_norm2 ( n, f )
c
c  Should we reject this step?
c
            if ( 1.2D+00 * ( fnrm_old + abstol ) .lt. fnrm ) then

              if ( damptol .lt. damp ) then
                damp = 0.50D+00 * damp
                go to 10
              else
                write ( *, '(7x,3x,14x,14x,14x,g14.6)' ) damp * snrm
                if ( damp .eq. 1.0D+00 ) then
                  write ( *, '(7x,i3,g14.6,g14.6)' ) it, xnrm, fnrm
                else
                  write ( *, '(7x,i3,g14.6,g14.6,g14.6)' )
     &              it, xnrm, fnrm, damp
                end if
                write ( *, '(a)' ) ' '
                write ( *, '(a)' ) '  Norm of function increased!'
                write ( *, '(a)' ) '  Iteration halted.'
                exit
              end if

            end if

            if ( 1000.0D+00 * ( xold + 1.0D+00 ) .lt. xnrm ) then
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) '  Norm of X blew up!'
              write ( *, '(a)' ) '  Iteration halted.'

              if ( damptol .lt. damp ) then
                damp = 0.25D+00 * damp
                go to 10
              else
                write ( *, '(7x,3x,14x,14x,14x,g14.6)' ) damp * snrm
                if ( damp == 1.0D+00 ) then
                  write ( *, '(7x,i3,g14.6,g14.6)' ) it, xnrm, fnrm
                else
                  write ( *, '(7x,i3,g14.6,g14.6,g14.6)' )
     &              it, xnrm, fnrm, damp
                end if
                exit
              end if

            end if

            write ( *, '(7x,3x,14x,14x,14x,g14.6)' ) damp * snrm
            if ( damp .eq. 1.0D+00 ) then
              write ( *, '(7x,i3,g14.6,g14.6)' ) it, xnrm, fnrm
            else
              write ( *, '(7x,i3,g14.6,g14.6,g14.6)' )
     &          it, xnrm, fnrm, damp
            end if
c
c  Accept the point as the next iterate.
c
            do i = 1, n
              x(i) = xnew(i)
            end do

            fnrm_old = fnrm
c
c  Should we stop at this step?
c
            if ( snrm .lt. reltol * ( xnrm + 1.0D+00 ) .and.
     &          fnrm .lt. abstol ) then
              converged = .true.
              write ( *, '(a)' ) ' '
              write ( *, '(a)' ) '  Convergence criteria satisfied.'
              exit
            end if

          end do

          if ( .not. converged ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' )
     &        '  The Newton iteration did not converge.'
          end if

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Final X:'
          write ( *, '(a)' ) ' '
          write ( *, '(5g14.6)' ) x(1:n)
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '  Final F(X):'
          write ( *, '(a)' ) ' '
          write ( *, '(5g14.6)' ) f(1:n)
c
c  Check against the exact solution, if it is known.
c
          call p00_sol ( problem, iknow, n, x )

          if ( 0 .lt. iknow ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) '  The exact solution X:'
            write ( *, '(a)' ) ' '
            write ( *, '(5g14.6)' ) ( x(i), i = 1, n )
            write ( *, '(a)' ) ' '
            call p00_fx ( problem, f, n, x )
            write ( *, '(a)' ) '  F(X):'
            write ( *, '(a)' ) ' '
            write ( *, '(5g14.6)' ) ( f(i), i = 1, n )
          else if ( iknow .eq. 0 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) '  The exact solution is not known.'
          else if ( iknow .eq. -1 ) then
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) '  An exact solution is not given.'
          end if

        end do

      end do

      return
      end
