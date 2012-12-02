      program main

c*********************************************************************72
c
cc MAIN is the main program for UNIFORM_PRB.
c
c  Discussion:
c
c    UNIFORM_PRB calls sample problems for the UNIFORM library.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    04 April 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'UNIFORM_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the UNIFORM library.'

      call test01
      call test02
      call test03
      call test04
      call test05
      call test06
      call test065
      call test07
      call test08
      call test09

      call test10
      call test11
      call test111
      call test118

      if ( .true. ) then
        call test119
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST119 is being SKIPPED!'
      end if

      call test12
      call test13
      call test14
      call test15
      call test16
      call test17
      call test18
      call test19

      call test20

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'UNIFORM_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests C4_UNIFORM_01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex c4_uniform_01
      integer i
      integer seed
      integer seed_init

      seed_init = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  C4_UNIFORM_01 computes pseudorandom'
      write ( *, '(a)' ) '  complex values uniformly distributed '
      write ( *, '(a)' ) '  in the unit circle.'

      seed = seed_init

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
      write ( *, '(a)' ) ' '

      do i = 1, 10
        write ( *, '(2x,i8,2x,2g14.6)' ) i, c4_uniform_01 ( seed )
      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests C4VEC_UNIFORM_01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      complex c(n)
      integer i
      integer seed
      integer seed_init

      seed_init = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  C4VEC_UNIFORM_01 computes pseudorandom '
      write ( *, '(a)' ) '  complex values uniformly distributed in '
      write ( *, '(a)' ) '  the unit circle.'

      seed = seed_init

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  The initial seed is ', seed_init

      call c4vec_uniform_01 ( n, seed, c )

      write ( *, '(a)' ) ' '

      do i = 1, 10
        write ( *, '(2x,i8,2x,2g14.6)' ) i, c(i)
      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests C8_UNIFORM_01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      integer seed
      integer seed_init
      parameter ( seed_init = 123456789 )
      double complex c8_uniform_01

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  C8_UNIFORM_01 computes pseudorandom '
      write ( *, '(a)' ) '  double precision complex values '
      write ( *, '(a)' ) '  uniformly distributed in the unit circle.'

      seed = seed_init

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
      write ( *, '(a)' ) ' '

      do i = 1, 10
        write ( *, '(2x,i8,2x,2g14.6)' ) i, c8_uniform_01 ( seed )
      end do

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests C8VEC_UNIFORM_01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      double complex c(n)
      integer i
      integer seed
      integer seed_init
      parameter ( seed_init = 123456789 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  C8VEC_UNIFORM_01 computes pseudorandom '
      write ( *, '(a)' ) '  double complex values uniformly '
      write ( *, '(a)' ) '  distributed in the unit circle.'

      seed = seed_init

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  The initial seed is ', seed_init

      call c8vec_uniform_01 ( n, seed, c )

      write ( *, '(a)' ) ' '

      do i = 1, 10
        write ( *, '(2x,i8,2x,2g14.6)' ) i, c(i)
      end do

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests CH_UNIFORM_AB.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      character ch_uniform_ab
      character chi
      character clo
      integer i
      integer seed
      integer seed_init

      seed_init = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  CH_UNIFORM_AB computes pseudorandom '
      write ( *, '(a)' ) '  characters in an interval [CLO,CHI].'

      clo = 'A'
      chi = 'J'
      seed = seed_init

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The lower endpoint CLO = "' // clo // '".'
      write ( *, '(a)' ) '  The upper endpoint CHI = "' // chi // '".'
      write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
      write ( *, '(a)' ) ' '

      do i = 1, 10
        write ( *, '(2x,i8,2x,a1)' ) i, ch_uniform_ab ( clo, chi, seed )
      end do

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 tests GET_SEED.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision r8_uniform_01
      integer i
      integer j
      integer seed
      integer seed_old

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) 
     &  '  GET_SEED picks an initial seed value for UNIFORM.'
      write ( *, '(a)' ) 
     &  '  The value chosen should vary over time, because'
      write ( *, '(a)' ) 
     &  '  the seed is based on reading the clock.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  This is just the "calendar" clock, which does'
      write ( *, '(a)' ) 
     &  '  not change very fast, so calling GET_SEED several'
      write ( *, '(a)' ) 
     &  '  times in a row may result in the same value.'

      seed = 12345678
      seed_old = seed

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  Initial SEED is ', seed
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Next 3 values of R8_UNIFORM:'
      write ( *, '(a)' ) ' '

      do j = 1, 3
        write ( *, '(g14.6)' ) r8_uniform_01 ( seed )
      end do

      do i = 1, 4

10      continue

          call get_seed ( seed )

          if ( seed .ne. seed_old ) then
            seed_old = seed
            go to 20
          end if

        go to 10

20      continue

        write ( *, '(a)' ) ' '
        write ( *, '(a,i12)' ) '  New seed from GET_SEED is ', seed
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Next 3 values of R8_UNIFORM_01:'
        write ( *, '(a)' ) ' '

        do j = 1, 3
          write ( *, '(g14.6)' ) r8_uniform_01 ( seed )
        end do

      end do

      return
      end
      subroutine test065 ( )

c*********************************************************************72
c
cc TEST065 tests I4_SEED_ADVANCE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 May 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i4_seed_advance
      integer seed
      integer seed_new
      integer step

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST065'
      write ( *, '(a)' ) '  I4_SEED_ADVANCE advances the seed.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Step        SEED input       SEED output'
      write ( *, '(a)' ) ' '

      seed_new = 12345

      do step = 1, 10

        seed = seed_new
        seed_new = i4_seed_advance ( seed )

        write ( *, '(2x,i4,2x,i16,2x,i16)' ) step, seed, seed_new

      end do

      return
      end
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07 tests I4_UNIFORM_AB.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 November 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer a
      parameter ( a = 6 )
      integer b
      parameter ( b = 10 )
      integer freq(6:10)
      integer i
      integer i4_uniform_ab
      integer j
      integer seed
      integer seed_init

      seed_init = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07'
      write ( *, '(a)' ) '  I4_UNIFORM_AB computes pseudorandom values '
      write ( *, '(a)' ) '  in an interval [A,B].'

      seed = seed_init

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The lower endpoint A = ', a
      write ( *, '(a,i8)' ) '  The upper endpoint B = ', b
      write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
      write ( *, '(a)' ) ' '

      do i = a, b
        freq(i) = 0
      end do

      do i = 1, 10000

        j = i4_uniform_ab ( a, b, seed )

        if ( j .lt. a ) then
          write ( *, '(a,i8)' ) '  Illegal value J = ', j
        else if ( j .le. b ) then
          freq(j) = freq(j) + 1
        else
          write ( *, '(a,i8)' ) '  Illegal value J = ', j
        end if

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         I    Frequency'
      write ( *, '(a)' ) ' '
      do i = a, b
        write ( *, '(2x,i8,2x,i8)' ) i, freq(i)
      end do

      return
      end
      subroutine test08 ( )

c*********************************************************************72
c
cc TEST08 tests I4_UNIFORM_AB.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 May 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer a
      parameter ( a = -100 )
      integer b
      parameter ( b = 200 )
      integer i
      integer i4_uniform_ab
      integer j
      integer seed
      integer seed_init
      parameter ( seed_init = 123456789 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08'
      write ( *, '(a)' ) '  I4_UNIFORM_AB computes pseudorandom values '
      write ( *, '(a)' ) '  in an interval [A,B].'

      seed = seed_init

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  The lower endpoint A = ', a
      write ( *, '(a,i12)' ) '  The upper endpoint B = ', b
      write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
      write ( *, '(a)' ) ' '

      do i = 1, 20

        j = i4_uniform_ab ( a, b, seed )

        write ( *, '(2x,i8,2x,i8)' ) i, j

      end do

      return
      end
      subroutine test09 ( )

c*********************************************************************72
c
cc TEST09 tests I4_UNIFORM_0I
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 1000 )

      integer i
      integer i4_uniform_0i
      real mean
      integer seed
      real variance
      integer x(n)
      integer x_max
      integer x_min

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST09'
      write ( *, '(a)' ) '  I4_UNIFORM_0I samples a uniform random'
      write ( *, '(a)' ) '  integer distribution in [0,2**31-1].'

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  Starting with seed = ', seed

      do i = 1, n
        x(i) = i4_uniform_0i ( seed )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  First few values:'
      write ( *, '(a)' ) ' '
      do i = 1, 5
        write ( *, '(2x,i8,2x,i12)' ) i, x(i)
      end do

      mean = 0.0E+00
      do i = 1, n
        mean = mean + real ( x(i) )
      end do
      mean = mean / dble ( n )

      variance = 0.0E+00
      do i = 1, n
        variance = variance + ( real ( x(i) ) - mean )**2
      end do
      variance = variance / real ( n - 1 )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of values computed was N = ', n
      write ( *, '(a,g14.6)' ) '  Average value was ', mean

      x_min = x(1)
      x_max = x(1)
      do i = 2, n
        x_min = min ( x_min, x(i) )
        x_max = max ( x_max, x(i) )
      end do
      write ( *, '(a,i12)' ) '  Minimum value was ', x_min
      write ( *, '(a,i12)' ) '  Maximum value was ', x_max
      write ( *, '(a,g14.6)' ) '  Variance was ', variance

      return
      end
      subroutine test10 ( )

c*********************************************************************72
c
cc TEST10 tests I4VEC_UNIFORM_AB.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 November 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10000 )

      integer a
      parameter ( a = 6 )
      integer b
      parameter ( b = 10 )
      integer freq(6:10)
      integer i
      integer i4vec(n)
      integer j
      integer seed
      integer seed_init

      seed_init = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST10'
      write ( *, '(a)' ) '  I4VEC_UNIFORM_AB computes a vector of'
      write ( *, '(a)' ) '  pseudorandom values in an interval [A,B].'

      seed = seed_init

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The lower endpoint A = ', a
      write ( *, '(a,i8)' ) '  The upper endpoint B = ', b
      write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
      write ( *, '(a)' ) ' '

      do i = a, b
        freq(i) = 0
      end do

      call i4vec_uniform_ab ( n, a, b, seed, i4vec )

      do i = 1, n

        if ( i4vec(i) .lt. a ) then
          write ( *, '(a,i8)' ) '  Illegal value J = ', i4vec(i)
        else if ( i4vec(i) .le. b ) then
          freq(i4vec(i)) = freq(i4vec(i)) + 1
        else
          write ( *, '(a,i8)' ) '  Illegal value J = ', i4vec(i)
        end if

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         I    Frequency'
      write ( *, '(a)' ) ' '
      do i = a, b
        write ( *, '(2x,i8,2x,i8)' ) i, freq(i)
      end do

      return
      end
      subroutine test11 ( )

c*********************************************************************72
c
cc TEST11 tests I8_UNIFORM_AB.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 May 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer*8 a
      integer*8 b
      integer*8 i
      integer*8 i8_uniform_ab
      integer*8 seed
      integer*8 seed_init 
c     parameter ( seed_init = 123456789 )
      parameter ( seed_init = 123456789987654321_8 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST11'
      write ( *, '(a)' ) '  I8_UNIFORM_AB computes pseudorandom values '
      write ( *, '(a)' ) '  in an interval [A,B].'

c     a = 100000
      a = 1000000000_8
c     b = 800000
      b = 8000000000_8

      seed = seed_init

      write ( *, '(a)' ) ' '
      write ( *, '(a,i24)' ) '  The lower endpoint A = ', a
      write ( *, '(a,i24)' ) '  The upper endpoint B = ', b
      write ( *, '(a,i24)' ) '  The initial seed is    ', seed_init
      write ( *, '(a)' ) ' '

      do i = 1, 10
        write ( *, '(2x,i8,2x,i24)' ) i, i8_uniform_ab ( a, b, seed )
      end do

      return
      end
      subroutine test111 ( )

c*********************************************************************72
c
cc TEST111 tests L_UNIFORM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 May 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      logical l_uniform
      integer seed
      integer seed_init
      parameter ( seed_init = 123456789 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST111'
      write ( *, '(a)' ) 
     &  '  L_UNIFORM computes pseudorandom logical values.'

      seed = seed_init

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
      write ( *, '(a)' ) ' '

      do i = 1, 10
        write ( *, '(2x,i8,2x,l1)' ) i, l_uniform ( seed )
      end do

      return
      end
      subroutine test118 ( )

c*********************************************************************72
c
cc TEST118 tests LCRG_ANBN.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 April 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer a
      integer an
      integer b
      integer bn
      integer c
      integer n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST118'
      write ( *, '(a)' ) 
     &  '  LCRG_ANBN determines a linear congruential random'
      write ( *, '(a)' ) 
     &  '  number generator equivalent to N steps of a given one.'
c
c  These parameters define the old (1969) IBM 360 random number generator:
c
      a = 16807
      b = 0
      c = 2147483647

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  LCRG parameters:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  A = ', a
      write ( *, '(a,i12)' ) '  B = ', b
      write ( *, '(a,i12)' ) '  C = ', c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             N             A             B'
      write ( *, '(a)' ) ' '

      do n = 0, 10
        call lcrg_anbn ( a, b, c, n, an, bn )
        write ( *, '(2x,i12,2x,i12,2x,i12)' ) n, an, bn
      end do

      return
      end
      subroutine test119 ( )

c*********************************************************************72
c
cc TEST119 tests LCRG_ANBN.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 April 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4 )

      integer a
      integer an
      integer b
      integer bn
      integer c
      integer j
      integer k
      integer u
      integer v
      integer x(n)
      integer y(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST119'
      write ( *, '(a)' ) 
     &  '  LCRG_ANBN determines a linear congruential random'
      write ( *, '(a)' ) 
     &  '  number generator equivalent to N steps of a given one.'
c
c  These parameters define the old (1969) IBM 360 random number generator:
c
      a = 16807
      b = 0
      c = 2147483647

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  LCRG parameters:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  A  = ', a
      write ( *, '(a,i12)' ) '  B  = ', b
      write ( *, '(a,i12)' ) '  C  = ', c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '                           N            In           Out'
      write ( *, '(a)' ) ' '

      k = 0
      u = 12345
      write ( *, '(2x,12x,2x,i12,2x,12x,2x,i12)' ) k, u
      do k = 1, 11
        call lcrg_evaluate ( a, b, c, u, v )
        write ( *, '(2x,12x,2x,i12,2x,i12,2x,i12)' ) k, u, v
        u = v
      end do
c
c  Now try to replicate these results using N procesors.
c
      call lcrg_anbn ( a, b, c, n, an, bn )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  LCRG parameters:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  AN = ', an
      write ( *, '(a,i12)' ) '  BN = ', bn
      write ( *, '(a,i12)' ) '  C  = ', c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '             J             N            In           Out'
      write ( *, '(a)' ) ' '

      x(1) = 12345
      do j = 2, n
        call lcrg_evaluate ( a, b, c, x(j-1), x(j) )
      end do

      do j = 1, n
        write ( *, '(2x,i12,2x,i12,2x,12x,2x,i12)' ) j, j-1, x(j)
      end do

      do k = n + 1, 12, n
        do j = 1, n
          call lcrg_evaluate ( an, bn, c, x(j), y(j) )
          write ( *, '(2x,i12,2x,i12,2x,i12,2x,i12)' ) 
     &      j, k+j-2, x(j), y(j)
          x(j) = y(j)
        end do
      end do

      return
      end
      subroutine test12 ( )

c*********************************************************************72
c
cc TEST12 tests LCRG_SEED.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 November 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer a
      integer b
      integer c
      integer i
      real r4_uniform_01
      integer seed
      integer seed_in
      integer seed_lcrg
      integer seed_out
      integer seed_start
      real u

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST12'
      write ( *, '(a)' ) '  LCRG_SEED directly computes the updated'
      write ( *, '(a)' ) '  value of a seed used by an linear'
      write ( *, '(a)' ) '  congruential random number generator.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '       I          SEED          SEED          SEED    U'
      write ( *, '(a)' ) 
     &  '                 Input        Output          LCRG'
      write ( *, '(a)' ) ' '
c
c  These parameters define the old (1969) IBM 360 random number generator:
c
      a = 16807
      b = 0
      c = 2147483647
c
c  This seed value was used in Pierre L'Ecuyer's article.
c
      seed_start = 12345

      seed = seed_start
c
c  Compute 1000 random numbers "the hard way", that is, sequentially.
c  Every now and then, call LCRG_SEED to compute SEED directly.
c
      do i = 1, 1000

        seed_in = seed
        u = r4_uniform_01 ( seed )
        seed_out = seed

        if ( i .le. 10 .or. i .eq. 100 .or. i .eq. 1000 ) then

          call lcrg_seed ( a, b, c, i, seed_start, seed_lcrg )

          write ( *, '(2x,i8,2x,i12,2x,i12,2x,i12,2x,g14.6)' ) 
     &      i, seed_in, seed_out, seed_lcrg, u

        end if

      end do

      return
      end
      subroutine test13 ( )

c*********************************************************************72
c
cc TEST13 tests R4_UNIFORM_AB.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      real a
      real b
      integer i
      real r4_uniform_ab
      integer seed
      integer seed_init
      parameter ( seed_init = 123456789 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST13'
      write ( *, '(a)' ) '  R4_UNIFORM_AB computes pseudorandom values '
      write ( *, '(a)' ) '  in an interval [A,B].'

      a = 5.0E+00
      b = 10.0E+00
      seed = seed_init

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  The lower endpoint A = ', a
      write ( *, '(a,g14.6)' ) '  The upper endpoint B = ', b
      write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
      write ( *, '(a)' ) ' '

      do i = 1, 10
        write ( *, '(2x,i8,2x,g14.6)' ) i, r4_uniform_ab ( a, b, seed )
      end do

      return
      end
      subroutine test14 ( )

c*********************************************************************72
c
cc TEST14 tests R4_UNIFORM_01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      real r4_uniform_01
      integer seed
      integer seed_init

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST14'
      write ( *, '(a)' ) '  R4_UNIFORM_01 computes pseudorandom values '
      write ( *, '(a)' ) '  in the interval [0,1].'

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  The initial seed is ', seed
      write ( *, '(a)' ) ' '

      do i = 1, 10
        write ( *, '(2x,i8,2x,g14.6)' ) i, r4_uniform_01 ( seed )
      end do

      return
      end
      subroutine test15 ( )

c*********************************************************************72
c
cc TEST15 tests R8_UNIFORM_AB.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision b
      integer i
      double precision r8_uniform_ab
      integer seed
      integer seed_init

      seed_init = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST15'
      write ( *, '(a)' ) '  R8_UNIFORM_AB computes pseudorandom values '
      write ( *, '(a)' ) '  in an interval [A,B].'

      a = 5.0D+00
      b = 10.0D+00
      seed = seed_init

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  The lower endpoint A = ', a
      write ( *, '(a,g14.6)' ) '  The upper endpoint B = ', b
      write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
      write ( *, '(a)' ) ' '

      do i = 1, 10
        write ( *, '(2x,i8,2x,g14.6)' ) i, r8_uniform_ab ( a, b, seed )
      end do

      return
      end
      subroutine test16 ( )

c*********************************************************************72
c
cc TEST16 tests R8_UNIFORM_01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      double precision r8_uniform_01
      integer seed
      integer seed_init

      seed_init = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST16'
      write ( *, '(a)' ) '  R8_UNIFORM_01 computes pseudorandom values '
      write ( *, '(a)' ) '  in the interval [0,1].'

      seed = seed_init

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
      write ( *, '(a)' ) ' '

      do i = 1, 10
        write ( *, '(2x,i8,2x,g14.6)' ) i, r8_uniform_01 ( seed )
      end do

      return
      end
      subroutine test17 ( )

c*********************************************************************72
c
cc TEST17 tests R8_UNIFORM_01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n

      parameter ( n = 1000 )

      integer i
      double precision r8_uniform_01
      integer seed
      integer seed_in
      integer seed_out
      double precision u(n)
      double precision u_avg
      double precision u_var

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST17'
      write ( *, '(a)' ) '  R8_UNIFORM_01 computes a sequence of '
      write ( *, '(a)' ) '  uniformly distributed pseudorandom numbers.'
c
c  Start with a known seed.
c
      seed = 12345

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  Initial SEED = ', seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  First 10 values:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '     I         Input        Output    R8_UNIFORM_01'
      write ( *, '(a)' ) '                SEED          SEED'
      write ( *, '(a)' ) ' '

      do i = 1, 10

        seed_in = seed
        u(i) = r8_uniform_01 ( seed )
        seed_out = seed

        write ( *, '(i8,2x,i12,2x,i12,2x,g14.6)' ) 
     &    i, seed_in, seed_out, u(i)

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,i10,a)' ) '  Now compute ', n, ' elements.'
      write ( *, '(a)' ) ' '

      u_avg = 0.0D+00
      do i = 1, n
        u(i) = r8_uniform_01 ( seed )
        u_avg = u_avg + u(i)
      end do

      u_avg = u_avg / dble ( n )

      u_var = 0.0D+00
      do i = 1, n
        u_var = u_var + ( u(i) - u_avg )**2
      end do

      u_var = u_var / dble ( n - 1 )

      write ( *, '(a)' ) ' '
      write ( *, '(a,f10.6)' ) '  Average value = ', u_avg
      write ( *, '(a,f10.6)' ) '  Expecting       ', 0.5D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,f10.6)' ) '  Variance =      ', u_var
      write ( *, '(a,f10.6)' ) '  Expecting       ', 1.0D+00 / 12.0D+00

      return
      end
      subroutine test18 ( )

c*********************************************************************72
c
cc TEST18 tests R8_UNIFORM_01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      double precision r8_uniform_01
      integer seed
      integer seed_in
      integer seed_out
      integer seed_save
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST18'
      write ( *, '(a)' ) '  R8_UNIFORM_01 computes a sequence of '
      write ( *, '(a)' ) '  pseudorandom numbers but all computations'
      write ( *, '(a)' ) '   depend on the seed value.'
      write ( *, '(a)' ) '  In this test, we show how a sequence of '
      write ( *, '(a)' ) '  "random" values can be manipulated by '
      write ( *, '(a)' ) '  accessing the seed.'

      seed = 1066

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  Initial SEED is ', seed
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Call R8_UNIFORM_01 10 times, and watch SEED.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '     I         Input        Output    R8_UNIFORM_01'
      write ( *, '(a)' ) '                SEED          SEED'
      write ( *, '(a)' ) ' '

      do i = 1, 10

        seed_in = seed

        if ( i == 5 ) then
          seed_save = seed
        end if

        x = r8_uniform_01 ( seed )
        seed_out = seed
        write ( *, '(i8,2x,i12,2x,i12,2x,g14.6)' ) 
     &    i, seed_in, seed_out,x

      end do

      seed = seed_save

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) 
     &  '  Reset SEED to its value at step 5, = ', seed
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Now call R8_UNIFORM_01 10 times, and watch'
      write ( *, '(a)' ) 
     &  '  how SEED and R8_UNIFORM_01 restart themselves.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '     I         Input        Output    R8_UNIFORM_01'
      write ( *, '(a)' ) '                SEED          SEED'
      write ( *, '(a)' ) ' '

      do i = 1, 10
        seed_in = seed
        x = r8_uniform_01 ( seed )
        seed_out = seed
        write ( *, '(i8,2x,i12,2x,i12,2x,g14.6)' ) 
     &    i, seed_in, seed_out, x
      end do

      seed = -12345678

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  What happens with a negative SEED?'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '     I         Input        Output    R8_UNIFORM_01'
      write ( *, '(a)' ) '                SEED          SEED'
      write ( *, '(a)' ) ' '

      do i = 1, 10
        seed_in = seed
        x = r8_uniform_01 ( seed )
        seed_out = seed
        write ( *, '(i8,2x,i12,2x,i12,2x,g14.6)' ) 
     &    i, seed_in, seed_out, x
      end do

      return
      end
      subroutine test19 ( )

c*********************************************************************72
c
cc TEST19 tests R8_UNIFORM_01 and R8MAT_UNIFORM_01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      integer n

      parameter ( m = 100 )
      parameter ( n = 10 )

      double precision a(m,n)
      double precision b(m,n)
      integer i
      integer j
      integer k
      double precision r8_uniform_01
      integer seed
      integer seed_init

      seed_init = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST19'
      write ( *, '(a)' ) '  R8_UNIFORM_01 computes pseudorandom values '
      write ( *, '(a)' ) '    one at a time.'
      write ( *, '(a)' ) 
     &  '  R8MAT_UNIFORM_01 computes a matrix of values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  For the same initial seed, the results'
      write ( *, '(a)' ) '  should be identical, but R8MAT_UNIFORM_01 '
      write ( *, '(a)' ) '  might be faster.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  The initial seed is ', seed_init

      seed = seed_init

      do j = 1, n
        do i = 1, m
          a(i,j) = r8_uniform_01 ( seed )
        end do
      end do

      seed = seed_init;
      call r8mat_uniform_01 ( m, n, seed, b )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '       I       J      A(I,J)        B(I,J)'
      write ( *, '(a)' ) 
     &  '                   (R8_UNIFORM_01)  (R8MAT_UNIFORM_01)'
      write ( *, '(a)' ) ' '

      do k = 0, 10
        i = ( k * m + ( 10 - k ) * 1 ) / 10
        j = ( k * n + ( 10 - k ) * 1 ) / 10
        write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g14.6)' ) 
     &    i, j, a(i,j), b(i,j)
      end do
  
      return
      end
      subroutine test20 ( )

c*********************************************************************72
c
cc TEST20 tests R8_UNIFORM_01 and R8VEC_UNIFORM_01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      double precision a(n)
      double precision b(n)
      integer i
      double precision r8_uniform_01
      integer seed
      integer seed_init
      parameter ( seed_init = 123456789 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST20'
      write ( *, '(a)' ) '  R8_UNIFORM_01 computes pseudorandom values '
      write ( *, '(a)' ) '  one at a time.'
      write ( *, '(a)' ) 
     &  '  R8VEC_UNIFORM_01 computes a vector of values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  For the same initial seed, the results'
      write ( *, '(a)' ) '  should be identical, but R8VEC_UNIFORM_01'
      write ( *, '(a)') '   might be faster.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  The initial seed is ', seed_init

      seed = seed_init

      do i = 1, n
        a(i) = r8_uniform_01 ( seed )
      end do

      seed = seed_init;
      call r8vec_uniform_01 ( n, seed, b )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       I      A(I)            B(I)'
      write ( *, '(a)' ) 
     &  '           (R8_UNIFORM_01)  (R8VEC_UNIFORM_01)'
      write ( *, '(a)' ) ' '

      do i = 1, 10
        write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) i, a(i), b(i)
      end do
  
      return
      end
