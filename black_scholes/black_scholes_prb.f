      program main

c*********************************************************************72
c
cc BLACK_SCHOLES_PRB tests BLACK_SCHOLES.
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
      implicit none

      call timestamp ( );
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BLACK_SCHOLES_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the BLACK_SCHOLES library.'

      call asset_path_test ( )
      call binomial_test ( )
      call bsf_test ( )
      call forward_test ( )
      call mc_test ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BLACK_SCHOLES_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      return
      end
      subroutine asset_path_test ( )

c*********************************************************************72
c
cc ASSET_PATH_TEST tests ASSET_PATH.
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
      implicit none

      integer n
      parameter ( n = 100 )

      double precision mu
      character * ( 100 ) output_filename
      double precision s(0:n)
      double precision s0
      integer seed
      double precision sigma
      double precision t1

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASSET_PATH_TEST:'
      write ( *, '(a)' ) 
     &  '  Demonstrate the simulated of an asset price path.'

      s0 = 2.0D+00
      mu = 0.1D+00
      sigma = 0.3D+00
      t1 = 1.0D+00
      seed = 123456789

      write ( *, '(a,g14.6)' ) ' '
      write ( *, '(a,g14.6)' ) 
     &  '  The asset price at time 0      S0    = ', s0
      write ( *, '(a,g14.6)' ) 
     &  '  The asset expected growth rate MU    = ', mu
      write ( *, '(a,g14.6)' ) 
     &  '  The asset volatility           SIGMA = ', sigma
      write ( *, '(a,g14.6)' ) 
     &  '  The expiry date                T1    = ', t1
      write ( *, '(a,i6)' )    
     &  '  The number of time steps       N     = ', n
      write ( *, '(a,i12)' )   
     &  '  The random number seed was     SEED  = ', seed

      call asset_path ( s0, mu, sigma, t1, n, seed, s )

      call r8vec_print_part ( n + 1, s, 10, '  Partial results:' )

      output_filename = 'asset_path.txt'
      call r8vec_write ( output_filename, n + 1, s )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Full results written to "' 
     &  // trim ( output_filename ) // '".'

      return
      end
      subroutine binomial_test ( )

c*********************************************************************72
c
cc BINOMIAL_TEST tests BINOMIAL.
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
      implicit none

      double precision c
      double precision e
      integer m
      double precision r
      double precision s0
      double precision sigma
      double precision t1

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BINOMIAL_TEST:'
      write ( *, '(a)' ) '  A demonstration of the binomial method'
      write ( *, '(a)' ) '  for option valuation.'

      s0 = 2.0D+00
      e = 1.0D+00
      r = 0.05D+00
      sigma = 0.25D+00
      t1 = 3.0D+00
      m = 256

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) 
     &  '  The asset price at time 0 S0    = ', s0
      write ( *, '(a,g14.6)' ) 
     &  '  The exercise price        E     = ', e
      write ( *, '(a,g14.6)' ) 
     &  '  The interest rate         R     = ', r
      write ( *, '(a,g14.6)' ) 
     &  '  The asset volatility      SIGMA = ', sigma
      write ( *, '(a,g14.6)' ) 
     &  '  The expiry date           T1    = ', t1
      write ( *, '(a,i8)' )    
     &  '  The number of intervals   M     = ', m

      call binomial ( s0, e, r, sigma, t1, m, c )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  The option value is ', c

      return
      end
      subroutine bsf_test ( )

c*********************************************************************72
c
cc BSF_TEST tests BSF.
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
      implicit none

      double precision c
      double precision e
      double precision r
      double precision s0
      double precision sigma
      double precision t0
      double precision t1

      write ( *, '(a)' )
      write ( *, '(a)' ) 'BSF_TEST:'
      write ( *, '(a)' ) 
     &  '  A demonstration of the Black-Scholes formula'
      write ( *, '(a)' ) '  for option valuation.'

      s0 = 2.0D+00
      t0 = 0.0D+00
      e = 1.0D+00
      r = 0.05D+00
      sigma = 0.25D+00
      t1 = 3.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) 
     &  '  The asset price at time T0 S0    = ', s0
      write ( *, '(a,g14.6)' ) 
     &  '  The time                   T0    = ', t0
      write ( *, '(a,g14.6)' ) 
     &  '  The exercise price         E     = ', e
      write ( *, '(a,g14.6)' ) 
     &  '  The interest rate          R     = ', r
      write ( *, '(a,g14.6)' ) 
     &  '  The asset volatility       SIGMA = ', sigma
      write ( *, '(a,g14.6)' ) 
     &  '  The expiry date            T1    = ', t1

      call bsf ( s0, t0, e, r, sigma, t1, c )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  The option value C = ', c

      return
      end
      subroutine forward_test ( )

c*********************************************************************72
c
cc FORWARD_TEST tests FORWARD.
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
      implicit none

      integer nt
      parameter ( nt = 29 )
      integer nx
      parameter ( nx = 11 )

      double precision e
      integer i
      double precision r
      double precision s
      double precision sigma
      double precision smax
      double precision smin
      double precision t1
      double precision u(nx-1,nt+1)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FORWARD_TEST:'
      write ( *, '(a)' ) 
     &  '  A demonstration of the forward difference method'
      write ( *, '(a)' ) '  for option valuation.'

      e = 4.0D+00
      r = 0.03D+00
      sigma = 0.50D+00
      t1 = 1.0D+00
      smax = 10.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) 
     &  '  The exercise price        E =     ', e
      write ( *, '(a,g14.6)' ) 
     &  '  The interest rate         R =     ', r
      write ( *, '(a,g14.6)' ) 
     &  '  The asset volatility      SIGMA = ', sigma;
      write ( *, '(a,g14.6)' ) 
     &  '  The expiry date           T1 =    ', t1
      write ( *, '(a,i8)' ) 
     &  '  The number of space steps NX =    ', nx
      write ( *, '(a,i8)' ) 
     &  '  The number of time steps  NT =    ', nt
      write ( *, '(a,g14.6)' ) 
     &  '  The value of              SMAX =  ', smax

      call forward ( e, r, sigma, t1, nx, nt, smax, u )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Initial   Option'
      write ( *, '(a)' ) '  Value     Value' 
      write ( *, '(a)' ) ' '

      smin = 0.0D+00
      do i = 1, nx - 1 
        s = ( ( nx - i - 1 ) * smin + i * smax ) / dble ( nx - 1 )
        write ( *, '(2x,g14.6,2x,g14.6)' ) s, u(i,nt+1)
      end do

      return
      end
      subroutine mc_test ( )

c*********************************************************************72
c
cc MC_TEST tests MC.
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
      implicit none

      double precision conf(2)
      double precision e
      integer m
      double precision r
      double precision s0
      integer seed
      double precision sigma
      double precision t1

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MC_TEST:'
      write ( *, '(a)' ) '  A demonstration of the Monte Carlo method'
      write ( *, '(a)' ) '  for option valuation.'

      s0 = 2.0D+00
      e = 1.0D+00
      r = 0.05D+00
      sigma = 0.25D+00
      t1 = 3.0D+00
      m = 1000000
      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a, g14.6)' ) 
     &  '  The asset price at time 0, S0    = ', s0
      write ( *, '(a, g14.6)' ) 
     &  '  The exercise price         E     = ', e
      write ( *, '(a, g14.6)' ) 
     &  '  The interest rate          R     = ', r
      write ( *, '(a, g14.6)' ) 
     &  '  The asset volatility       SIGMA = ', sigma
      write ( *, '(a, g14.6)' ) 
     &  '  The expiry date            T1    = ', t1
      write ( *, '(a, i8)' )    
     &  '  The number of simulations  M     = ', m

      call mc ( s0, e, r, sigma, t1, m, seed, conf )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,g14.6,a)' ) 
     &  '  The confidence interval is [', conf(1), ',', conf(2), '].'

      return
      end
