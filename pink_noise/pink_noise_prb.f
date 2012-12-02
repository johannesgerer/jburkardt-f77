      program main

c*********************************************************************72
c
cc MAIN tests the PINK_NOISE routines.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'PINK_NOISE_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )

      write ( *, '(a)' ) ' ' 
      write ( *, '(a)' ) 'PINK_NOISE_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      return
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests WRAP2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      integer m
      integer q
      integer q_in

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  WRAP2 performs a circular wrap.'
      write ( *, '(a)' ) '  Q is expected to range between 0 and M.'
      write ( *, '(a)' ) '  WRAP2 takes an input value of Q, and either'
      write ( *, '(a)' ) '  increments it by M+1 until in the range, or'
      write ( *, '(a)' ) '  decrements it by M+1 until in the range,'
      write ( *, '(a)' ) 
     &  '  and returns the result as the function value.'

      do m = 2, 4
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '   M  Qin  Qout'
        write ( *, '(a)' ) ' '
        do i = -5, 3 * m - 1
          q = i
          q_in = q
          call wrap2 ( m, q )
          write ( *, '(2x,i2,2x,i2,2x,i2)' ) m, q_in, q
        end do
      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests CDELAY2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      integer m
      integer q
      integer q_in

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  CDELAY2 is a circular buffer implementation'
      write ( *, '(a)' ) '  of an M-fold delay.  Q is a counter'
      write ( *, '(a)' ) '  which is decremented by CDELAY2, but reset'
      write ( *, '(a)' ) '  to M after it reaches 0.'

      do m = 2, 4
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '   I   M  Qin  Qout'
        write ( *, '(a)' ) ' '
        q = m
        do i = 1, 3 * ( m + 1 )
          q_in = q
          call cdelay2 ( m, q )
          write ( *, '(2x,i2,2x,i2,2x,i2,2x,i2)' ) i, m, q_in, q
        end do
      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests RANH.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer d
      integer i
      integer q
      double precision ranh
      integer seed
      double precision u
      double precision y

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  RANH is a random hold function.'
      write ( *, '(a)' ) '  Given a value U and a delay D, it returns'
      write ( *, '(a)' ) '  the value U for D calls, then resets U.'

      do d = 5, 1, -1
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '   I   D   Q      U           Y'
        write ( *, '(a)' ) ' '
        u = 0.5D+00
        q = 3
        do i = 1, 20
          y = ranh ( d, u, q, seed )
          write ( *, '(2x,i2,2x,i2,2x,i2,2x,f10.6,2x,f10.6)' ) 
     &      i, d, q, u, y
        end do
      end do

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests RAN1F.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    01 June 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer b_max
      parameter ( b_max = 31 )

      integer b
      integer i
      integer q(31)
      double precision r8_uniform_01
      double precision ran1f
      integer rep
      integer seed
      double precision u(31)
      double precision y

      seed = 12345689

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  RAN1F generates random values with an'
      write ( *, '(a)' ) '  approximate 1/F distribution.'

      b = 1

10    continue

      if ( b .lt. 32 ) then

        do rep = 1, 4

          do i = 1, b
            u(i) = r8_uniform_01 ( seed )
          end do

          do i = 1, b
            q(i) = 0
          end do

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) '   B   I      Y'
          write ( *, '(a)' ) ' '

          do i = 1, 20
            y = ran1f ( b, u, q, seed )
            write ( *, '(2x,i2,2x,i2,2x,f10.6)' ) b, i, y 
          end do

        end do

        b = b * 2

        go to 10

      end if

      return
      end
