      program main

c*********************************************************************72
c
cc BROWNIAN_MOTION_SIMULATION_PRB tests the BROWNIAN_MOTION_SIMULATION library.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 September 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer k
      parameter ( k = 40 )
      integer m_max
      parameter ( m_max = 3 )
      integer n
      parameter ( n = 1001 )

      double precision d
      double precision dsq(k,n)
      character * ( 80 ) header
      integer m
      integer seed
      double precision t
      double precision x(m_max,n)

      call timestamp ( )
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'BROWNIAN_MOTION_SIMULATION_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) 
     &  '  Test the BROWNIAN_MOTION_SIMULATION library.'
c
c  Compute the path of a particle undergoing Brownian motion.
c
      do m = 1, 2

        d = 10.0D+00
        t = 1.0D+00
        seed = 123456789
        call brownian_motion_simulation ( m, n, d, t, seed, x )
        if ( m .eq. 1 ) then
          header = 'motion_1d'
        else if ( m .eq. 2 ) then
          header = 'motion_2d'
        end if
        call brownian_motion_display ( m, n, x, header )

      end do
c
c  Estimate the average displacement of the particle from the origin
c  as a function of time.
c
      do m = 1, 3

        d = 10.0D+00
        t = 1.0D+00
        seed = 123456789

        call brownian_displacement_simulation ( k, n, m, d, t, seed, 
     &    dsq )
        if ( m .eq. 1 ) then
          header = 'displacement_1d'
        else if ( m .eq. 2 ) then
          header = 'displacement_2d'
        else if ( m .eq. 3 ) then
          header = 'displacement_3d'
        end if
        call brownian_displacement_display ( k, n, d, t, dsq, header )

      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'BROWNIAN_MOTION_SIMULATION_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ''
      call timestamp ( );

      return
      end
