      program main

c*********************************************************************72
c
cc MAIN is the main program for ASA113_PRB.
c
c  Discussion:
c
c    ASA113_PRB tests the ASA113 clustering algorithm.
c
c  Modified:
c
c    16 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA113_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the ASA113 library.'

      call test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA113_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tries out the ASA113 routine.
c
c  Modified:
c
c   16 February 2008
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
      integer c(m)
      double precision c_center(k,n)
      integer c_size(k)
      integer ci
      double precision critvl
      integer i
      integer ifault
      integer j
      integer ntrans1
      integer ntrans2
      double precision wss(k)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Test the ASA113 classification algorithm.'
c
c  Read the data
c
      open ( unit = 1, file = 'points_100.txt', status = 'old' )

      do i = 1, m
        read ( 1, * ) ( a(i,j), j = 1, n )
      end do

      close ( unit = 1 )
c
c  Print first five points.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  First five points:'
      write ( *, '(a)' ) ' '
      do i = 1, 5
        write ( *, '(2x,i8,2x,g14.6,2x,g14.6)' ) i, ( a(i,j), j = 1, n )
      end do
c
c  Assign points randomly to classes.
c
      do i = 1, m
        c(i) = mod ( i, k ) + 1
      end do
c
c  Define the critical value as the sum of the squares of the distances
c  of the points to their cluster center.
c
      do i = 1, k
        c_size(i) = 0
        do j = 1, n
          c_center(i,j) = 0.0D+00
        end do
      end do

      do i = 1, m
        ci = c(i)
        c_size(ci) = c_size(ci) + 1
        do j = 1, n
          c_center(ci,j) = c_center(ci,j) + a(i,j)
        end do
      end do

      do i = 1, k
        do j = 1, n
          c_center(i,j) = c_center(i,j) / dble ( c_size(i) )
        end do
      end do

      do i = 1, k
        wss(i) = 0.0D+00
      end do

      do i = 1, m
        ci = c(i)
        do j = 1, n
          wss(ci) = wss(ci) + ( a(i,j) - c_center(ci,j) )**2
        end do
      end do

      critvl = 0.0D+00
      do i = 1, k
        critvl = critvl + wss(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '        Initial CRITVL = ', critvl
c
c  Compute the clusters.
c
      ntrans1 = -1
      ntrans2 = -1

10    continue

        call trnsfr ( a, c, c_size, m, k, n, critvl, ntrans1, ifault )

        if ( ntrans1 .eq. 0 .and. ntrans2 .eq. 0 ) then
          go to 20
        end if

        write ( *, '(a,g14.6)' ) '  After TRNSFR, CRITVL = ', critvl

        call swap ( a, c, c_size, m, k, n, critvl, ntrans2, ifault )

        if ( ntrans1 .eq. 0 .and. ntrans2 .eq. 0 ) then
          go to 20
        end if

        write ( *, '(a,g14.6)' ) '    After SWAP, CRITVL = ', critvl

      go to 10

20    continue
c
c  Define the critical value as the sum of the squares of the distances
c  of the points to their cluster center.
c
      do i = 1, k
        do j = 1, n
          c_center(i,j) = 0.0D+00
        end do
      end do

      do i = 1, m
        ci = c(i)
        do j = 1, n
          c_center(ci,j) = c_center(ci,j) + a(i,j)
        end do
      end do

      do i = 1, k
        do j = 1, n
          c_center(i,j) = c_center(i,j) / dble ( c_size(i) )
        end do
      end do

      do i = 1, k
        wss(i) = 0.0D+00
      end do

      do i = 1, m
        ci = c(i)
        do j = 1, n
          wss(ci) = wss(ci) + ( a(i,j) - c_center(ci,j) )**2
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Cluster  Population  Energy'
      write ( *, '(a)' ) ' '

      do i = 1, k
        write ( *, '(2x,i8,2x,i8,2x,g14.6)' ) i, c_size(i), wss(i)
      end do
      write ( *, '(a)' ) ' '
      write ( *, '(2x,a8,2x,i8,2x,g14.6)' ) '   Total', m, critvl

      return
      end
      subroutine crswap ( a, c, c_size, m, k, n, critvl, i1, i2,
     &  c1, c2, iswitch, inc )

c*********************************************************************72
c
cc CRSWAP determines the effect of swapping two objects.
c
c  Discussion:
c
c    This computation is very inefficient.  It is only set up so that we
c    can compare algorithm ASA 113 to the K-means algorithms ASA 058 and
c    ASA 136.
c
c  Modified:
c
c    15 February 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Colin Banfield, LC Bassill,
c    Algorithm AS 113:
c    A transfer for non-hierarchichal classification,
c    Applied Statistics,
c    Volume 26, Number 2, 1977, pages 206-210.
c
c  Parameters:
c
c    Input, double precision A(M,N), the data values.  There are M objects,
c    each having spatial dimension N.
c
c    Input, integer C(M), the classification of each object.
c
c    Input, integer C_SIZE(K), the number of objects in each class.
c
c    Input, integer M, the number of objects.
c
c    Input, integer K, the number of classes.
c
c    Input, integer N, the number of spatial dimensions, or variates,
c    of the objects.
c
c    Input, double precision CRITVL, the current value of the criterion.
c
c    Input, integer I1, I2, the objects to be swapped.
c
c    Input, integer C1, C2, the current classes of objects I1 and I2.
c
c    Input, integer ISWITCH:
c    1, indicates that I1 and I2 should be temporarily swapped, the
c       change in CRITVL should be computed, and then I1 and I2 restored.
c    2, indicates that I1 and I2 will be swapped.
c
c    Output, double precision INC, the change to CRITVL that would occur if I1 and
c    I2 were swapped.  This is only computed for ISWITCH = 1.
c
      implicit none

      integer k
      integer m
      integer n

      integer k_copy
      parameter ( k_copy = 5 )
      integer m_copy
      parameter ( m_copy = 100 )
      integer n_copy
      parameter ( n_copy = 2 )

      double precision a(m,n)
      integer c(m)
      double precision c_center(k_copy,n_copy)
      integer c_size(k)
      integer c1
      integer c2
      integer ci
      double precision critvl
      double precision critvl_new
      integer i
      integer i1
      integer i2
      double precision inc
      integer iswitch
      integer j
      double precision wss(k_copy)

      if ( iswitch .eq. 2 ) then
        return
      end if
c
c  Move object I1 from class C1 to class C2.
c  Move object I2 from class C2 to class C1.
c
      c(i1) = c2
      c(i2) = c1
c
c  Define the critical value as the sum of the squares of the distances
c  of the points to their cluster center.
c
      do i = 1, k
        c_size(i) = 0
        do j = 1, n
          c_center(i,j) = 0.0D+00
        end do
      end do

      do i = 1, m
        ci = c(i)
        c_size(ci) = c_size(ci) + 1
        do j = 1, n
          c_center(ci,j) = c_center(ci,j) + a(i,j)
        end do
      end do

      do i = 1, k
        do j = 1, n
          c_center(i,j) = c_center(i,j) / dble ( c_size(i) )
        end do
      end do

      do i = 1, k
        wss(i) = 0.0D+00
      end do

      do i = 1, m
        ci = c(i)
        do j = 1, n
          wss(ci) = wss(ci) + ( a(i,j) - c_center(ci,j) )**2
        end do
      end do

      critvl_new = 0.0D+00
      do i = 1, k
        critvl_new = critvl_new + wss(i)
      end do

      inc = critvl_new - critvl
!
!  Move object I1 from class C2 to class C1.
!  Move object I2 from class C1 to class C2.
!
      c(i1) = c1
      c(i2) = c2

      return
      end
      subroutine crtran ( a, c, c_size, m, k, n, critvl, i1, c1, c2,
     &  iswitch, inc )

c*********************************************************************72
c
cc CRTRAN determines the effect of moving an object to another class.
c
c  Discussion:
c
c    This computation is very inefficient.  It is only set up so that we
c    can compare algorithm ASA 113 to the K-means algorithms ASA 058 and
c    ASA 136.
c
c  Modified:
c
c    15 February 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Colin Banfield, LC Bassill,
c    Algorithm AS 113:
c    A transfer for non-hierarchichal classification,
c    Applied Statistics,
c    Volume 26, Number 2, 1977, pages 206-210.
c
c  Parameters:
c
c    Input, double precision A(M,N), the data values.  There are M objects,
c    each having spatial dimension N.
c
c    Input, integer C(M), the classification of each object.
c
c    Input, integer C_SIZE(K), the number of objects in each class.
c
c    Input, integer M, the number of objects.
c
c    Input, integer K, the number of classes.
c
c    Input, integer N, the number of spatial dimensions, or variates,
c    of the objects.
c
c    Input, double precision CRITVL, the current value of the criterion.
c
c    Input, integer I1, the object to be transferred.
c
c    Input, integer C1, C2, the current class of object I1, and the
c    class to which it may be transferred.
c
c    Input, integer ISWITCH:
c    1, indicates that I1 should be temporarily transferred, the
c       change in CRITVL should be computed, and then I1 restored.
c    2, indicates that I1 will be permanently transferred.
c
c    Output, double precision INC, the change to CRITVL that would occur if I1 were
c    transferred from class C1 to C2.  This is only computed for ISWITCH = 1.
c
      implicit none

      integer k
      integer m
      integer n

      integer k_copy
      parameter ( k_copy = 5 )
      integer m_copy
      parameter ( m_copy = 100 )
      integer n_copy
      parameter ( n_copy = 2 )

      double precision a(m,n)
      integer c(m)
      double precision c_center(k_copy,n_copy)
      integer c_size(k)
      integer c1
      integer c2
      integer ci
      double precision critvl
      double precision critvl_new
      integer i
      integer i1
      double precision inc
      integer iswitch
      integer j
      double precision wss(k_copy)

      if ( iswitch .eq. 2 ) then
        return
      end if
c
c  Move object I from class C1 to class C2.
c
      c(i1) = c2
      c_size(c1) = c_size(c1) - 1
      c_size(c2) = c_size(c2) + 1
c
c  Define the critical value as the sum of the squares of the distances
c  of the points to their cluster center.
c
      do i = 1, k
        c_size(i) = 0
        do j = 1, n
          c_center(i,j) = 0.0D+00
        end do
      end do

      do i = 1, m
        ci = c(i)
        c_size(ci) = c_size(ci) + 1
        do j = 1, n
          c_center(ci,j) = c_center(ci,j) + a(i,j)
        end do
      end do

      do i = 1, k
        do j = 1, n
          c_center(i,j) = c_center(i,j) / dble ( c_size(i) )
        end do
      end do

      do i = 1, k
        wss(i) = 0.0D+00
      end do

      do i = 1, m
        ci = c(i)
        do j = 1, n
          wss(ci) = wss(ci) + ( a(i,j) - c_center(ci,j) )**2
        end do
      end do

      critvl_new = 0.0D+00
      do i = 1, k
        critvl_new = critvl_new + wss(i)
      end do

      inc = critvl_new - critvl
!
!  Move object I1 from class C2 to class C1.
!
      c(i1) = c1
      c_size(c1) = c_size(c1) + 1
      c_size(c2) = c_size(c2) - 1

      return
      end
