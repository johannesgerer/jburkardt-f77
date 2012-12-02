      program main

c*********************************************************************72
c
cc MAIN is the main program for STOCHASTIC_RK_PRB.
c
c  Discussion:
c
c    STOCHASTIC_RK_PRB calls a set of problems for STOCHASTIC_RK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 June 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'STOCHASTIC_RK_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the STOCHASTIC_RK library.'
     
      call test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'STOCHASTIC_RK_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'
       
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests RK1_TI_STEP.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 June 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      double precision fi
      external fi
      double precision gi
      external gi
      double precision h
      integer i
      double precision q
      integer seed
      double precision t
      double precision t0
      parameter ( t0 = 0.0D+00 )
      double precision tn
      parameter ( tn = 1.0D+00 )
      double precision x(0:n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  RK1_TI_STEP uses a first order RK method'
      write ( *, '(a)' ) '  for a problem whose right hand side does'
      write ( *, '(a)' ) '  not depend explicitly on time.'

      h = ( tn - t0 ) / dble ( n )
      q = 1.0D+00
      seed = 123456789

      i = 0
      t = t0
      x(i) = 0.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         I           T             X'
      write ( *, '(a)' ) ' '
      write ( *, '(2x,i8,2x,f14.6,2x,g14.6)' ) i, t, x(i)

      do i = 1, n

        t = ( dble ( n - i ) * t0   
     &      + dble (     i ) * tn ) 
     &      / dble ( n     )

        call rk1_ti_step ( x(i-1), t, h, q, fi, gi, seed, x(i) )

        write ( *, '(2x,i8,2x,f14.6,2x,g14.6)' ) i, t, x(i)

      end do

      return
      end
      function fi ( x )

c*********************************************************************72
c
cc FI is a time invariant deterministic right hand side.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision FI, the value.
c
      implicit none

      double precision fi
      double precision x

      fi = 1.0D+00

      return
      end
      function gi ( x )

c*********************************************************************72
c
cc GI is a time invariant stochastic right hand side.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    20 June 2010
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision X, the argument.
c
c    Output, double precision GI, the value.
c
      implicit none

      double precision gi
      double precision x

      gi = 1.0D+00

      return
      end
