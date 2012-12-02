      program main

c*********************************************************************72
c
cc TOMS448_PRB tests COUNT.
c
c  Modified:
c
c    12 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS448_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS448 library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS448_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests COUNT.
c
c  Modified:
c
c    12 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer k_max
      integer n_max

      parameter ( k_max = 10 )
      parameter ( n_max = 10 )

      integer c(k_max)
      external count
      integer i
      integer k
      integer n
      integer p(n_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Let C = 1, 2, 3, ..., N'
      write ( *, '(a)' ) ' '

      n = 10
      k = 10

      do i = 1, k
        c(i) = i
      end do

      call count ( c, k, p, n )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       I    P(I)'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,i6)' ) i, p(i)
      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests COUNT.
c
c  Modified:
c
c    12 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer k_max
      integer n_max

      parameter ( k_max = 10 )
      parameter ( n_max = 10 )

      integer c(k_max)
      external count
      integer i
      integer k
      integer n
      integer p(n_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02:'
      write ( *, '(a)' ) '  Let C = 1, 3, 5, 7, ( N + 1 ) / 2'
      write ( *, '(a)' ) ' '

      n = 10
      k = ( n + 1 ) / 2

      do i = 1, k
        c(i) = 2 * i - 1
      end do

      call count ( c, k, p, n )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       I    P(I)'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,i6)' ) i, p(i)
      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests COUNT.
c
c  Modified:
c
c    12 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer k_max
      integer n_max

      parameter ( k_max = 10 )
      parameter ( n_max = 10 )

      integer c(k_max)
      external count
      integer i
      integer k
      integer n
      integer p(n_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03:'
      write ( *, '(a)' ) '  Let C = 1, 2, 2'
      write ( *, '(a)' ) ' '

      n = 10
      k = 3

      c(1) = 1
      c(2) = 2
      c(3) = 2

      call count ( c, k, p, n )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       I    P(I)'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,i6)' ) i, p(i)
      end do

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests COUNT.
c
c  Modified:
c
c    12 January 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer k_max
      integer n_max

      parameter ( k_max = 10 )
      parameter ( n_max = 100 )

      integer c(k_max)
      external count
      integer i
      integer k
      integer n
      integer p(n_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04:'
      write ( *, '(a)' ) '  Let C = 1, 5, 10, 25, 50, 100'
      write ( *, '(a)' ) ' '

      n = 100
      k = 6

      c(1) = 1
      c(2) = 5
      c(3) = 10
      c(4) = 25
      c(5) = 50
      c(6) = 100

      call count ( c, k, p, n )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       I    P(I)'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i6,2x,i6)' ) i, p(i)
      end do

      return
      end

