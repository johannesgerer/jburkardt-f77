      program main

c*********************************************************************72
c
cc POLYGON_MOMENTS_PRB tests POLYGON_MOMENTS.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 October 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POLYGON_MOMENTS_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test POLYGON_MOMENTS library.'

      call test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'POLYGON_MOMENTS_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 carries out a test on a rectangle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 October 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      double precision alpha_exact(6)
      double precision alpha_pq
      integer k
      double precision mu_exact(6)
      double precision mu_pq
      double precision nu_exact(6)
      double precision nu_pq
      integer p
      integer q
      integer s
      double precision x(n)
      double precision y(n)

      save alpha_exact
      save mu_exact
      save nu_exact
      save x
      save y

      data alpha_exact / 
     &  1.0D+00, 
     &  5.0D+00, 4.0D+00, 
     &  30.66666666666667D+00, 22.0D+00, 18.66666666666666D+00 /
      data mu_exact /
     &  1.0D+00, 
     &  0.0D+00, 0.0D+00, 
     &  5.666666666666667D+00, 2.0D+00, 2.666666666666667D+00 /
      data nu_exact /
     &  40.0D+00, 
     &  200.0D+00, 160.0D+00, 
     &  1226.66666666666667D+00, 880.0D+00, 746.66666666666666D+00 /
      data x/
     &  2.0D+00, 10.0D+00, 8.0D+00, 0.0D+00 /
      data y /
     &  0.0D+00,  4.0D+00, 8.0D+00, 4.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Check normalized moments of a rectangle.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   P   Q             Nu(P,Q)'
      write ( *, '(a)' ) '            Computed         Exact'
      write ( *, '(a)' ) ' '
      k = 0
      do s = 0, 2
        do p = s, 0, -1
          q = s - p
          k = k + 1
          call moment ( n, x, y, p, q, nu_pq )
          write ( *, '(2x,i2,2x,i2,2x,g14.6,2x,g14.6)' ) 
     &      p, q, nu_pq, nu_exact(k)
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   P   Q           Alpha(P,Q)'
      write ( *, '(a)' ) '            Computed         Exact'
      write ( *, '(a)' ) ' '
      k = 0
      do s = 0, 2
        do p = s, 0, -1
          q = s - p
          k = k + 1
          call moment_normalized ( n, x,y, p, q, alpha_pq )
          write ( *, '(2x,i2,2x,i2,2x,g14.6,2x,g14.6)' ) 
     &      p, q, alpha_pq, alpha_exact(k)
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   P   Q             Mu(P,Q)'
      write ( *, '(a)' ) '            Computed         Exact'
      write ( *, '(a)' ) ' '
      k = 0
      do s = 0, 2
        do p = s, 0, -1
          q = s - p
          k = k + 1
          call moment_central ( n, x, y , p, q, mu_pq )
          write ( *, '(2x,i2,2x,i2,2x,g14.6,2x,g14.6)' ) 
     &      p, q, mu_pq, mu_exact(k)
        end do
      end do

      return
      end
