      program main

c*********************************************************************72
c
c MAIN is the main program for the reactor shielding simulation.
c
c  Discussion:
c
c    This is a Monte Carlo simulation, using
c    uniform random numbers, which investigates the
c    effectiveness of a shield intended to absorb the
c    neutrons emitted from a nuclear reactor.
c   
c    The reactor is modeled as a point source,
c    located at (0,0,0).
c   
c    A particle emitted from the reactor has a random
c    initial direction, and an energy selected from
c    [Emin,Emax] with a 1/Sqrt(E) distribution.
c   
c    The shield is modeled as a wall of thickness THICK,
c    extending from 0 to THICK in the X direction, and
c    extending forever in the Y and Z directions.
c   
c    Based on the particle energy, a distance D is computed
c    which measures how far the particle could travel through
c    the shield before colliding.
c   
c    Based on the particle direction, the position is updated
c    by D units.
c   
c    If the particle is now to the left of the shield, it is
c    counted as being REFLECTED.
c   
c    If the particle is to the right of the shield, it is 
c    counted as being ABSORBED.
c   
c    If the particle is inside the shield, it has COLLIDED.
c    A particle that collides is either absorbed (end of story)
c    or SCATTERED with a new random direction and a new (lower)
c    energy.
c   
c    Every particle is followed from origin to its final fate,
c    which is reflection, transmission, or absorption.
c    At the end, a summary is printed, giving the number of
c    particles with each fate, and the average energy of each
c    group of particles.
c   
c    Increasing NTOT, the number of particles used, will improve the
c    expected reliability of the results.
c   
c    Increasing THICK, the thickness of the shield, should 
c    result in more absorptions and reflections.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 September 2012
c
c  Author:
c
c    Original FORTRAN77 version by Kahaner, Moler, Nash.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Local Parameters:
c
c    Local, double precision AZM, the azimuthal angle of the particle's
c    direction.
c
c    Local, double precision D, the distance that the particle can
c    travel through the slab, given its current energy.
c
c    Local, double precision E, the energy of the particle.
c
c    Local, double precision EA, energy absorbed by the slab.
c
c    Local, double precision ER, energy reflected by the slab.
c
c    Local, double precision ET, energy transmitted through the slab.
c
c    Local, double precision MU, the cosine of the angle between the
c    particle's direction and the X axis.
c
c    Local, integer NA, number of particles absorbed by the slab.
c
c    Local, integer NPART, the index of the current particle.
c
c    Local, integer NR, number of particles reflected by the slab.
c
c    Local, integer NT, number of particles transmitted by the slab.
c
c    Local, integer NTOT, the total number of particles to be
c    emitted from the neutron source.
c
c    Local, double precision SA, standard deviation of absorbed energy.
c
c    Local, double precision SR, standard deviation of reflected energy.
c
c    Local, double precision ST, standard deviation of transmitted energy.
c
c    Local, double precision THICK, the thickness of the slab that is
c    intended to absorb most of the particles.
c
c    Local, double precision X, Y, Z, the current position of the particle.
c
      implicit none

      logical absorb
      double precision azm
      double precision d
      double precision dist2c
      double precision e
      double precision ea
      double precision er
      double precision et
      integer i
      double precision mu
      integer na
      integer npart
      integer nr
      integer nt
      integer ntot
      parameter ( ntot = 100000 )
      double precision sa
      integer seed
      double precision sr
      double precision st
      integer test
      integer test_num
      parameter ( test_num = 5 )
      double precision thick
      parameter ( thick = 2.0D+00 )
      double precision x
      double precision y
      double precision z

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'REACTOR_SIMULATION'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  The reactor shielding simulation.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Shield thickness is THICK = ', thick
      write ( *, '(a,i8)' ) 
     &  '  Number of simulated particles is NTOT = ', ntot
      write ( *, '(a,i8)' ) '  Number of tests TEST_NUM = ', test_num

      seed = 123456789

      do test = 1, test_num

        write ( *, '(a)' ) ' '
        write ( *, '(a,i2)' ) '  Test # ', test
        write ( *, '(a,i12)' ) '  SEED = ', seed
c
c  Initialize.
c
        ea = 0.0D+00
        er = 0.0D+00
        et = 0.0D+00
        na = 0
        nr = 0
        nt = 0
        sa = 0.0D+00
        sr = 0.0D+00
        st = 0.0D+00
c
c  Main loop over NTOT particles.
c
        do npart = 1, ntot
c
c  Generate a new particle.
c
          call source ( seed, e, mu, azm, x, y, z )

10        continue
c
c  Compute the distance that the particle can travel through the slab,
c  based on its current energy.
c
            d = dist2c ( e, seed )
c
c  Update the particle's position by D units.
c
            call update ( mu, azm, d, x, y, z )
c
c  The particle was reflected by the shield, and this path is complete.
c
            if ( x .lt. 0.0D+00 ) then

              nr = nr + 1
              er = er + e
              sr = sr + e * e
              go to 20
c
c  The particle was transmitted through the shield, and this path is complete.
c
            else if ( thick < x ) then

              nt = nt + 1
              et = et + e
              st = st + e * e
              go to 20
c
c  The particle collided with the shield, and was absorbed.  This path is done.
c
            else if ( absorb ( seed ) ) then

              na = na + 1
              ea = ea + e
              sa = sa + e * e
              go to 20
c
c  The particle collided with the shield and was scattered.
c  Find the scattering angle and energy, and continue along the new path.
c
            else

              call scatter ( seed, e, mu, azm )

            end if

          go to 10

20        continue

        end do
c
c  Print the results of the simulation.
c
        call output ( na, ea, sa, nr, er, sr, nt, et, st, ntot )

      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'REACTOR_SIMULATION:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      function absorb ( seed )

c*********************************************************************72
c
c ABSORB determines if a colliding particle is absorbed.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 September 2012
c
c  Author:
c
c    Original FORTRAN77 version by Kahaner, Moler, Nash.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input/output, integer SEED, a seed for the random
c    number generator.
c
c    Output, logical ABSORB, is TRUE if the particle is absorbed.
c
c  Local parameters:
c
c    Local, double precision PA, the probability of absorption.
c
      implicit none

      logical absorb
      double precision pa
      parameter ( pa = 0.1D+00 )
      double precision r8_uniform_01
      integer seed
      double precision u

      u = r8_uniform_01 ( seed )

      if ( u .le. pa ) then
        absorb = .true.
      else
        absorb = .false.
      end if

      return
      end
      function cross ( e )

c*********************************************************************72
c
c CROSS returns the "cross section" of a particle based on its energy.
c
c  Discussion:
c
c    The particle's cross section is a measure of its likelihood to collide
c    with the material of the slab.  This quantity typically depends on both
c    the particle's energy and the kind of medium through which it is traveling.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 September 2012
c
c  Author:
c
c    Original FORTRAN77 version by Kahaner, Moler, Nash.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, double precision E, the energy of the particle.
c
c    Output, double precision CROSS, the cross section.
c
      implicit none

      double precision cross
      double precision e
      double precision s
      double precision y

      s = abs ( sin ( 100.0D+00 * ( exp ( e ) - 1.0D+00 ) ) 
     &  + sin ( 18.81D+00 * ( exp ( e ) - 1.0D+00 ) ) )

      y = max ( 0.02D+00, s )

      cross = 10.0D+00 * exp ( -0.1D+00 / y )

      return
      end
      function dist2c ( e, seed )

c*********************************************************************72
c
c DIST2C returns the distance to collision.
c
c  Discussion:
c
c    Assuming the particle has a given energy, and assuming it is currently
c    somewhere inside the shield, it is possible to determine a typical distance
c    which the particle can travel before it collides with the material of
c    the shield.
c
c    The computation of the collision distance is made by estimating a
c    "cross section" (as though having more energy made the particle "bigger"
c    and hence more likely to collide) and then randomly selecting a distance
c    that is logarithmically distributed.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 September 2012
c
c  Author:
c
c    Original FORTRAN77 version by Kahaner, Moler, Nash.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, double precision E, the energy of the particle.
c
c    Input/output, integer SEED, a seed for the random
c    number generator.
c
c    Output, double precision DIST2C, the distance the particle can travel
c    through the slab before colliding.
c
      implicit none

      double precision cross
      double precision dist2c
      double precision e
      double precision r8_uniform_01
      integer seed
      double precision u

      u = r8_uniform_01 ( seed )

      dist2c = - log ( u ) / cross ( e )

      return
      end
      function energy ( seed )

c*********************************************************************72
c
c ENERGY assigns an energy to an emitted particle.
c
c  Discussion:
c
c    The energy E is in the range [EMIN,EMAX], with distribution
c    const/sqrt(energy).
c
c    An inverse function approach is used to compute this.
c
c    The energies are measured in MeV.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 September 2012
c
c  Author:
c
c    Original FORTRAN77 version by Kahaner, Moler, Nash.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input/output, integer SEED, a seed for the random
c    number generator.
c
c    Output, double precision ENERGY, a randomly chosen energy that is
c    distributed as described above.
c
c  Local parameters:
c
c    Local, double precision EMIN, EMAX, the minimum and maximum
c    energies.
c
      implicit none

      double precision c
      double precision emax
      parameter ( emax = 2.5D+00 )
      double precision emin
      parameter ( emin = 1.0D-03 )
      double precision energy
      double precision r8_uniform_01
      integer seed
      double precision u

      u = r8_uniform_01 ( seed )

      c = 1.0D+00 / ( 2.0D+00 * ( sqrt ( emax ) - sqrt ( emin ) ) )

      energy = ( u / ( 2.0D+00 * c ) + sqrt ( emin ) )
      energy = energy * energy

      return
      end
      subroutine output ( na, ea, sa, nr, er, sr, nt, et, st, ntot )

c*********************************************************************72
c
c OUTPUT prints the results of the reactor shielding simulation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 September 2012
c
c  Author:
c
c    Original FORTRAN77 version by Kahaner, Moler, Nash.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, integer NA, number of particles absorbed by the slab.
c
c    Input, double precision EA, energy absorbed by the slab.
c
c    Input, double precision SA, the sum of the squares of the 
c    absorbed energies.
c
c    Input, integer NR, number of particles reflected by the slab.
c
c    Input, double precision ER, energy reflected by the slab.
c
c    Input, double precision SR, the sum of the squares of the 
c    reflected energies.
c
c    Input, integer NT, number of particles transmitted by the slab.
c
c    Input, double precision ET, energy transmitted through the slab.
c
c    Input, double precision ST, the sum of the squares of the 
c    transmitted energies.
c
c    Input, integer NTOT, the total number of particles.
c
      implicit none

      double precision ea
      double precision ea_ave
      double precision er
      double precision er_ave
      double precision et
      double precision et_ave
      double precision etot
      integer na
      integer nr
      integer nt
      integer ntot
      double precision pa
      double precision pr
      double precision pt
      double precision ptot
      double precision sa
      double precision sr
      double precision st

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The Reactor Shielding Problem:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '                           Total                   Average'
      write ( *, '(a)' ) '                    #      Energy      ' 
     &  // 'Percent     Energy         StDev'
      write ( *, '(a)' ) ' '

      etot = ea + er + et

      if ( 0 .lt. na ) then
        ea_ave = ea / dble ( na )
        sa = sqrt ( sa / dble ( na ) - ea_ave * ea_ave )
      else
        ea_ave = 0.0D+00
      end if

      pa = dble ( na * 100 ) / dble ( ntot )

      write ( *, '(a,2x,i8,2x,g14.6,2x,f6.2,2x,g14.6,2x,g14.6)' ) 
     &  'Absorbed   ', na, ea, pa, ea_ave, sa

      if ( 0 .lt. nr ) then
        er_ave = er / dble ( nr )
        sr = sqrt ( sr / dble ( nr ) - er_ave * er_ave )
      else
        er_ave = 0.0D+00
      end if

      pr = dble ( nr * 100 ) / dble ( ntot )

      write ( *, '(a,2x,i8,2x,g14.6,2x,f6.2,2x,g14.6,2x,g14.6)' )  
     &  'Reflected  ', nr, er, pr, er_ave, sr

      if ( 0 .lt. nt ) then
        et_ave = et / dble ( nt )
        st = sqrt ( st / dble ( nt ) - et_ave * et_ave )
      else
        et_ave = 0.0D+00
      end if

      pt = dble ( nt * 100 ) / dble ( ntot )

      write ( *, '(a,2x,i8,2x,g14.6,2x,f6.2,2x,g14.6,2x,g14.6)' )  
     &  'Transmitted', nt, et, pt, et_ave, st

      ptot = 100.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,2x,i8,2x,g14.6,2x,f6.2,2x,g14.6,2x,g14.6)' )  
     &  'Total      ', ntot, etot, ptot

      return
      end
      function r8_uniform_01 ( seed )

c*********************************************************************72
c
cc R8_UNIFORM_01 returns a unit pseudorandom R8.
c
c  Discussion:
c
c    This routine implements the recursion
c
c      seed = 16807 * seed mod ( 2^31 - 1 )
c      r8_uniform_01 = seed / ( 2^31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      R8_UNIFORM_01
c      SEED      SEED
c
c         12345   207482415  0.096616
c     207482415  1790989824  0.833995
c    1790989824  2035175616  0.947702
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 August 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Springer Verlag, pages 201-202, 1983.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley Interscience, page 95, 1998.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, pages 362-376, 1986.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, pages 136-143, 1969.
c
c  Parameters:
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R8_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      double precision r8_uniform_01
      integer k
      integer seed

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + 2147483647
      end if

      r8_uniform_01 = dble ( seed ) * 4.656612875D-10

      return
      end
      subroutine scatter ( seed, e, mu, azm )

c*********************************************************************72
c
c SCATTER returns the new direction and energy of a particle that is scattered.
c
c  Discussion:
c
c    The scattering direction is chosen uniformly on the sphere.
c
c    The energy of the scattered particle is chosen uniformly in
c    [ 0.3*E, E ].
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 September 2012
c
c  Author:
c
c    Original FORTRAN77 version by Kahaner, Moler, Nash.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input/output, integer SEED, a seed for the random
c    number generator.
c
c    Input/output, double precision E.  On input, the particle energy
c    before collision.  On output, the particle energy after collision
c    and scattering.
c
c    Output, double precision MU, the cosine of the angle between the
c    particle's direction and the X axis.
c
c    Output, double precision AZM, the azimuthal angle of the particle's
c    direction.
c
      implicit none

      double precision azm
      double precision e
      double precision mu
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r8_uniform_01
      integer seed
      double precision u

      u = r8_uniform_01 ( seed )
      mu = - 1.0D+00 + 2.0D+00 * u

      u = r8_uniform_01 ( seed )
      azm = u * 2.0D+00 * pi

      u = r8_uniform_01 ( seed )
      e = ( u * 0.7D+00 + 0.3D+00 ) * e

      return
      end
      subroutine source ( seed, e, mu, azm, x, y, z )

c*********************************************************************72
c
c SOURCE generates a new particle from the neutron source.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 September 2012
c
c  Author:
c
c    Original FORTRAN77 version by Kahaner, Moler, Nash.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input/output, integer SEED, a seed for the random
c    number generator.
c
c    Output, double precision E, the initial energy of the particle.
c
c    Output, double precision MU, the cosine of the angle between the
c    particle's direction and the X axis.
c
c    Output, double precision AZM, the azimuthal angle of the particle's
c    direction.
c
c    Output, double precision X, Y, Z, the initial coordinates of the particle.
c
      implicit none

      double precision azm
      double precision e
      double precision energy
      double precision mu
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision r8_uniform_01
      integer seed
      double precision u
      double precision x
      double precision y
      double precision z

      u = r8_uniform_01 ( seed )
      mu = u

      u = r8_uniform_01 ( seed )
      azm = u * 2.0D+00 * pi

      x = 0.0D+00
      y = 0.0D+00
      z = 0.0D+00

      e = energy ( seed )

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
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
      subroutine update ( mu, azm, d, x, y, z )

c*********************************************************************72
c
c UPDATE determines the position of the particle after it has traveled D units.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    02 September 2012
c
c  Author:
c
c    Original FORTRAN77 version by Kahaner, Moler, Nash.
c    This FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    David Kahaner, Cleve Moler, Steven Nash,
c    Numerical Methods and Software,
c    Prentice Hall, 1989,
c    ISBN: 0-13-627258-4,
c    LC: TA345.K34.
c
c  Parameters:
c
c    Input, double precision MU, the cosine of the angle between the
c    particle's direction and the X axis.
c
c    Input, double precision AZM, the azimuthal angle of the particle's
c    direction.
c
c    Input, double precision D, the distance the particle traveled.
c
c    Input/output, double precision X, Y, Z.  On input, the previous
c    coordinates of the particle.  On output, the updated coordinates of the
c    particle.
c
      implicit none

      double precision azm
      double precision d
      double precision mu
      double precision s
      double precision x
      double precision y
      double precision z

      s = sqrt ( 1.0D+00 - mu * mu )

      x = x + d * mu
      y = y + d * s * cos ( azm )
      z = z + d * s * sin ( azm )

      return
      end
