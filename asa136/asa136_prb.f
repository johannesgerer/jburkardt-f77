      program main

c*********************************************************************72
c
cc MAIN is the main program for ASA136_PRB.
c
c  Discussion:
c
c    ASA136_PRB tests the ASA136 clustering algorithm.
c
c  Modified:
c
c    21 February 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA136_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the ASA136 library.'

      call test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA136_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tries out the ASA136 routine.
c
c  Modified:
c
c    14 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer k
      parameter ( k = 5 )
      integer m
      parameter ( m = 100 )
      integer n
      parameter ( n = 2 )

      double precision a(m,n)
      double precision an1(k)
      double precision an2(k)
      double precision c(k,n)
      double precision d(m)
      integer i
      integer ic1(m)
      integer ic2(m)
      integer ifault
      integer iter
      integer itran(k)
      integer j
      integer live(k)
      integer nc(k)
      integer nc_sum
      integer ncp(k)
      double precision wss(k)
      double precision wss_sum

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Test the KMNS algorithm,'
      write ( *, '(a)' ) '  Applied Statistics Algorithm #136.'
c
c  Read the data.
c
      open ( unit = 1, file = 'points_100.txt', status = 'old' )

      do i = 1, m
        read ( 1, * ) ( a(i,j), j = 1, n )
      end do

      close ( unit = 1 )
c
c  Print a few data values.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  First 5 data values:'
      write ( *, '(a)' ) ' '

      do i = 1, 5
        write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) i, ( a(i,j), j = 1, n )
      end do
c
c  Initialize the cluster centers.
c   Here, we arbitrarily make the first K data points cluster centers.
c
      do i = 1, k
        do j = 1, n
          c(i,j) = a(i,j)
        end do
      end do

      iter = 50
c
c  Compute the clusters.
c
      call kmns ( a, m, n, c, k, ic1, ic2, nc, an1, an2, ncp, d,
     &  itran, live, iter, wss, ifault )

      if ( ifault .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST01 - Fatal error!'
        write ( *, '(a,i8)' ) '  KMNS returned IFAULT = ', ifault
        return
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Cluster  Population  Energy'
      write ( *, '(a)' ) ' '

      nc_sum = 0
      wss_sum = 0.0D+00

      do i = 1, k
        write ( *, '(2x,i8,2x,i8,2x,g14.6)' ) i, nc(i), wss(i)
        nc_sum = nc_sum + nc(i)
        wss_sum = wss_sum + wss(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(2x,a8,2x,i8,2x,g14.6)' ) '   Total', nc_sum, wss_sum

      return
      end
