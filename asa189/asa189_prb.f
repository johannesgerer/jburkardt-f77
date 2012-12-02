      program main

c*********************************************************************72
c
cc MAIN is the main program for ASA189_PRB.
c
c  Discussion:
c
c    ASA189_PRB controls the test of ASA189.
c
c  Modified:
c
c    08 December 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer sample_max
      parameter ( sample_max = 100000 )

      double precision a
      double precision a_est
      double precision b
      double precision b_est
      integer c
      integer in(sample_max)
      double precision mu
      double precision mu_est
      double precision p(sample_max)
      integer sample_log
      integer sample_num
      integer seed
      double precision theta
      double precision theta_est
      double precision w(sample_max)
      double precision x(sample_max)

      call timestamp ( )

      a = 2.0D+00
      b = 3.0D+00
      c = 4

      call beta_binomial_check ( a, b, c )

      mu = a / ( a + b )
      theta = 1.0D+00 / ( a + b )

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA189_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the ASA189 library.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Applied Statistics Algorithm 189'
      write ( *, '(a)' ) '  Estimate the parameters A and B of a '
      write ( *, '(a)' ) '  beta binomial distribution.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Exact values:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '          A             B                MU         THETA'
      write ( *, '(a)' ) ' '
      write ( *, '(6x,4g14.6)' ) a, b, mu, theta
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Estimated values:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  'Samples   A_est         B_est        MU_est     THETA_est'
      write ( *, '(a)' ) ' '

      do sample_log = 2, 5

        sample_num = 10**sample_log

        call analyze ( sample_num, a, b, c, in, w, p, x, seed,
     &  mu_est, theta_est )
c
c  Convert the ASA189 "THETA" and "MU" parameters to "A" and "B".
c
        a_est = mu_est / theta_est
        b_est = ( 1.0D+00 - mu_est ) / theta_est

        write ( *, '(i6,4g14.6)' )
     &  sample_num, a_est, b_est, mu_est, theta_est

      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA189_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine analyze ( sample_num, a, b, c, in, w, p, x, seed,
     &  mu_est, theta_est )

c*********************************************************************72
c
cc ANALYZE generates data and analyzes it with ASA189.
c
c  Discussion:
c
c    The routine to generate the samples uses parameter A, B and C,
c    while ASA189 prefers a related form MU, THETA.  The calling routine
c    has to figure out how these are related.
c
c  Modified:
c
c    09 December 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer SAMPLE_NUM, the number of samples to generate.
c
c    Input, double precision A, double precision B, integer C, the parameters to use in
c    the beta binomial distribution.
c
c    Workspace, integer IN(SAMPLE_NUM).
c
c    Workspace, double precision W(SAMPLE_NUM), P(SAMPLE_NUM).
c
c    Workspace, double precision X(SAMPLE_NUM).
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, double precision MU_EST, THETA_EST, the estimates of MU and THETA
c    produced by ASA189.
c
      implicit none

      integer mrl
      parameter ( mrl = 4 )
      integer sample_num

      double precision a
      double precision b
      integer c
      double precision ccrit
      integer i
      integer ifault
      integer in(sample_num)
      integer iter
      double precision mu_est
      double precision mu_se
      double precision p(sample_num)
      integer rl(mrl,3)
      double precision rnl
      integer sample_i
      integer seed
      double precision theta_est
      double precision theta_se
      double precision w(sample_num)
      integer x(sample_num)
c
c  Generate the sample data using the exact parameters A, B.
c
      do sample_i = 1, sample_num
        call beta_binomial_sample ( a, b, c, seed, x(sample_i) )
      end do
c
c  Analyze the sample data, trying to estimate the parameters
c  in the form "MU", "THETA".  Note that ASA189 expects to be told
c  the value of C (although C could vary from one sample to the next).
c
      do i = 1, sample_num
        in(i) = c
      end do

      iter = 10
      ccrit = 0.001D+00

      call bbml ( sample_num, x, in, w, p, rl, mrl, iter, ccrit,
     &  mu_est, theta_est, mu_se, theta_se, rnl, ifault )

      return
      end
