      program main

c*********************************************************************72
c
cc DUEL_SIMULATION simulates a duel.
c
c  Discussion:
c
c    Player 1 fires at player 2, and hits with a probability of P(1).
c    If Player 2 misses, then Player 2 fires at Player 1, hitting with
c    a probability of P(2).
c
c    The duel continues with alternating shots until only one player 
c    survives.
c
c    The simulation is intended to estimate the probabilities that a
c    player will survive, and the number of turns required.
c
c    The exact probability that player 1 will survive is
c
c      P(1) / ( P(1) + P(2) - P(1) * P(2) )
c
c    Player 2's chance is
c
c     P(2) * ( 1 - P(1) ) / ( P(1) + P(2) - P(1) * P(2) )
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    17 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Nahin,
c    Duelling Idiots and Other Probability Puzzlers,
c    Princeton University Press, 2000,
c    ISBN13: 978-0691009797,
c    LC: QA273.N29.
c
c    Martin Shubik,
c    "Does the Fittest Necessarily Survive?",
c    in Readings in Game Theory and Political Behavior,
c    edited by Martin Shubik,
c    Doubleday, 1954,
c    LC: H61.S53.
c
c  Parameters:
c
c    Input, double precision A_ACCURACY, B_ACCURACY, the probabilities that 
c    players A and B will hit their opponent in a single shot.
c
c    Input, integer DUEL_NUM, the number of duels to run.
c
c    Output, double precision A_PROB, B_PROB, the estimated probablities that 
c    players A and B will survive.
c
c    Output, double precision TURN_AVERAGE, the average number of turns 
c    required to complete the duel.
c
      implicit none

      double precision a_accuracy
      double precision a_prob
      integer  a_wins
      double precision b_accuracy
      double precision b_prob
      integer  b_wins
      integer duel
      integer duel_num
      integer seed
      integer winner

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'DUEL_SIMULATION:'
      write ( *, '(a)' ) '  FORTRAN77 version'

      write ( *, '(a)' ) '  Enter number of duels to run:'
      read ( *, * ) duel_num

      write ( *, '(a)' ) 
     &  '  Enter player A''s accuracy between 0.0 and 1.0:'
      read ( *, * ) a_accuracy

      write ( *, '(a)' ) 
     &  '  Enter player B''s accuracy between 0.0 and 1.0:'
      read ( *, * ) b_accuracy

      a_wins = 0
      b_wins = 0

      do duel = 1, duel_num

        call duel_result ( a_accuracy, b_accuracy, seed, winner )

        if ( winner .eq. 1 ) then
          a_wins = a_wins + 1
        else
          b_wins = b_wins + 1
        end if

      end do

      a_prob = dble ( a_wins ) / dble ( duel_num )
      b_prob = dble ( b_wins ) / dble ( duel_num )

      write ( *, '(a)' ) ''
      write ( *, '(a,g14.6)' ) 
     &  '  Player A wins with probability ', a_prob
      write ( *, '(a,g14.6)' ) 
     &  '  Player B wins with probability ', b_prob

      write ( *, '(a)' ) ''
      write ( *, '(a)' ) 'DUEL_SIMULATION:'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
      subroutine duel_result ( a_accuracy, b_accuracy, seed, winner )

c*********************************************************************72
c
cc DUEL_RESULT returns the outcome of a single duel.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 September 2012
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Martin Shubik,
c    “Does the Fittest Necessarily Survive?”,
c    in Readings in Game Theory and Political Behavior,
c    edited by Martin Shubik,
c    Doubleday, 1954,
c    LC: H61.S53.
c
c  Parameters:
c
c    Input, double precision A_ACCURACY, B_ACCURACY, the probabilities that 
c    player A and B will hit their opponent in a single shot.
c
c    Input/output, integer SEED, a seed for the random number generator.
c
c    Output, integer WINNER, the survivor of the duel.
c
      implicit none

      double precision a_accuracy
      double precision b_accuracy
      double precision r
      double precision r8_uniform_01
      integer seed
      integer winner

10    continue

        r = r8_uniform_01 ( seed )

        if ( r .le. a_accuracy ) then
          winner = 1
          go to 20
        end if

        r = r8_uniform_01 ( seed )

        if ( r .le. b_accuracy ) then
          winner = 2
          go to 20
        end if

      go to 10

20    continue

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
