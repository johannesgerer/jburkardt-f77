      program main

c*********************************************************************72
c
cc MAIN is the main program for SDE_PRB.
c
c  Discussion:
c
c    SDE_PRB demonstrates the use of SDE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SDE_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the SDE library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
      call test06 ( )
      call test07 ( )
      call test08 ( )
      call test09 ( )
      call test10 ( )
      call test11 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SDE_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests BPATH.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 500 )

      integer seed
      double precision w(0:n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) 
     &  '  BPATH generates a sample Brownian motion path'

      seed = 123456789

      call bpath ( seed, n, w )

      call bpath_gnuplot ( n, w )

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests BPATH_AVERAGE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 1000 )
      integer n
      parameter ( n = 500 )

      double precision error
      integer seed
      double precision u(m,0:n)
      double precision umean(0:n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02:'
      write ( *, '(a)' ) '  BPATH_AVERAGE generates many Brownian paths'
      write ( *, '(a)' ) '  and averages them.'

      seed = 123456789

      call bpath_average ( seed, m, n, u, umean, error )

      call bpath_average_gnuplot ( m, n, u, umean )

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests CHAIN.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 200 )

      double precision diff
      integer seed
      double precision vem(0:n)
      double precision xem(0:n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03:'
      write ( *, '(a)' ) 
     &  '  CHAIN solves a stochastic differential equation for'
      write ( *, '(a)' ) '  a function of a stochastic variable X.'
      write ( *, '(a)' ) 
     &  '  We can solve for X(t), and then evaluate V(X(t)).'
      write ( *, '(a)' ) 
     &  '  Or, we can apply the stochastic chain rule to derive an'
      write ( *, '(a)' ) '  an SDE for V, and solve that.'

      seed = 123456789

      call chain ( seed, n, xem, vem, diff )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Maximum | Sqrt(X) - V | = ', diff

      call chain_gnuplot ( n, xem, vem )

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests EM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 256 )

      double precision diff
      integer seed
      double precision t(0:n)
      double precision t2(0:n/4)
      double precision xtrue(0:n)
      double precision xem(0:n/4)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04:'
      write ( *, '(a)' ) 
     &  '  EM solves a stochastic differential equation'
      write ( *, '(a)' ) '  using the Euler-Maruyama method.'

      seed = 123456789

      call em ( seed, n, t, xtrue, t2, xem, diff )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  | Exact X(T) - EM X(T) | = ', diff

      call em_gnuplot ( n, t, xtrue, t2, xem )

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests EMSTRONG.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 100 )
      integer n
      parameter ( n = 512 )
      integer p_max
      parameter ( p_max = 6 )

      double precision dtvals(p_max)
      integer seed
      double precision xerr(p_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05:'
      write ( *, '(a)' ) 
     &  '  EMSTRONG investigates the strong convergence'
      write ( *, '(a)' ) '  of the Euler-Maruyama method.'

      seed = 123456789

      call emstrong ( seed, m, n, p_max, dtvals, xerr )

      call emstrong_gnuplot ( p_max, dtvals, xerr )

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 tests EMWEAK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 50000 )
      integer p_max
      parameter ( p_max = 5 )

      double precision dtvals(p_max)
      integer method
      integer seed
      double precision xerr(p_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06:'
      write ( *, '(a)' ) '  EMWEAK investigates the weak convergence'
      write ( *, '(a)' ) '  of the Euler-Maruyama method.'

      seed = 123456789
      method = 0

      call emweak ( seed, method, m, p_max, dtvals, xerr )

      call emweak_gnuplot ( p_max, dtvals, xerr, method )

      seed = 123456789
      method = 1

      call emweak ( seed, method, m, p_max, dtvals, xerr )

      call emweak_gnuplot ( p_max, dtvals, xerr, method )

      return
      end
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07 tests MILSTRONG.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer p_max
      parameter ( p_max = 4 )

      double precision dtvals(p_max)
      integer seed
      double precision xerr(p_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07:'
      write ( *, '(a)' ) 
     &  '  MILSTRONG investigates the strong convergence'
      write ( *, '(a)' ) '  of the Milstein method.'

      seed = 123456789

      call milstrong ( seed, p_max, dtvals, xerr )

      call milstrong_gnuplot ( p_max, dtvals, xerr )

      return
      end
      subroutine test08 ( )

c*********************************************************************72
c
cc TEST08 tests STAB_ASYMPTOTIC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 1000 )
      integer p_max
      parameter ( p_max = 3 )
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08:'
      write ( *, '(a)' ) '  STAB_ASYMPTOTIC investigates the asymptotic'
      write ( *, '(a)' ) '  stability of the Euler-Maruyama method.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  For technical reasons, the plotting is done'
      write ( *, '(a)' ) '  in the same routine as the computations.'

      seed = 123456789

      call stab_asymptotic ( seed, n, p_max )

      return
      end
      subroutine test09 ( )

c*********************************************************************72
c
cc TEST09 tests STAB_MEANSQUARE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST09:'
      write ( *, '(a)' ) 
     &  '  STAB_MEANSQUARE investigates the mean square'
      write ( *, '(a)' ) '  stability of the Euler-Maruyama method.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  For technical reasons, the plotting is done'
      write ( *, '(a)' ) '  in the same routine as the computations.'

      seed = 123456789

      call stab_meansquare ( seed )

      return
      end
      subroutine test10 ( )

c*********************************************************************72
c
cc TEST10 tests STOCHASTIC_INTEGRAL_ITO.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    21 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision error
      double precision estimate
      double precision exact
      integer i
      integer n
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST10:'
      write ( *, '(a)' ) 
     &  '  Estimate the Ito integral of W(t) dW over [0,1].'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '                                                 ' //
     &  'Abs          Rel'
      write ( *, '(a)' ) 
     &  '         N        Exact        Estimate          ' //
     &  'Error        Error'
      write ( *, '(a)' ) ' '

      n = 100
      seed = 123456789

      do i = 1, 7

        call stochastic_integral_ito ( n, seed, estimate, exact, error )

        write ( *, '(2x,i8,2x,g16.8,2x,g16.8,2x,g10.2,2x,g10.2)' ) 
     &    n, exact, estimate, error, error / exact

        n = n * 4

      end do

      return
      end
      subroutine test11 ( )

c*********************************************************************72
c
cc TEST11 tests STOCHASTIC_INTEGRAL_STRAT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    22 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision error
      double precision estimate
      double precision exact
      integer i
      integer n
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST11:'
      write ( *, '(a)' ) 
     &  '  Estimate the Stratonovich integral of W(t) dW over [0,1].'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '                                                 ' //
     &  'Abs          Rel'
      write ( *, '(a)' ) 
     &  '         N        Exact        Estimate          ' //
     &  'Error        Error'
      write ( *, '(a)' ) ' '

      n = 100
      seed = 123456789

      do i = 1, 7

        call stochastic_integral_strat ( n, seed, estimate, 
     &   exact, error )

        write ( *, '(2x,i8,2x,g16.8,2x,g16.8,2x,g10.2,2x,g10.2)' ) 
     &    n, exact, estimate, error, error / exact

        n = n * 4

      end do

      return
      end


