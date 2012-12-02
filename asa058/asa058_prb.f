      program main

c*********************************************************************72
c
cc MAIN is the main program for ASA058_PRB.
c
c  Discussion:
c
c    ASA058_PRB tests the ASA058 clustering algorithm.
c
c  Modified:
c
c    22 February 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA058_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the ASA058 library.'

      call test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA058_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tries out the ASA058 routine.
c
c  Modified:
c
c    04 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer k
      parameter ( k = 5 )
      integer m
      parameter ( m = 2 )
      integer n
      parameter ( n = 100 )

      integer b(n)
      double precision d(k,m)
      double precision dev(k)
      double precision dev_sum
      integer e(k)
      integer e_sum
      double precision f(n)
      integer i
      integer j
      integer k2
      integer nz
      double precision x(n,m)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Test the CLUSTR algorithm.'
      write ( *, '(a)' ) '  Applied Statistics Algorithm 58'
c
c  Read the data.
c
      write (  *, '(a)' ) ' '
      write (  *, '(a)' ) '  Reading the data.'

      open ( unit = 1, file = 'points_100.txt', status = 'old' )

      do i = 1, n
        read ( 1, * ) ( x(i,j), j = 1, m )
      end do

      close ( unit = 1 )
c
c  Print a few data values.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  First 5 data values:'
      write ( *, '(a)' ) ' '

      do i = 1, 5
        write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) i, ( x(i,j), j = 1, m )
      end do
c
c  Initialize the cluster centers arbitrarily.
c
      do i = 1, k
        do j = 1, m
          d(i,j) = x(i,j)
        end do
      end do
c
c  Compute the clusters.
c
      nz = 1
      k2 = k

      call clustr ( x, d, dev, b, f, e, n, m, k, nz, k2 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Cluster  Population  Energy'
      write ( *, '(a)' ) ' '

      do i = 1, k
        write ( *, '(2x,i8,2x,i8,2x,g14.6)' ) i, e(i), dev(i)
      end do

      e_sum = 0
      dev_sum = 0.0D+00

      do i = 1, k
        e_sum = e_sum + e(i)
        dev_sum = dev_sum + dev(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(2x,a8,2x,i8,2x,g14.6)' ) '   Total', e_sum, dev_sum

      return
      end
