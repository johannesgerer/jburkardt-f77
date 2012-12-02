      subroutine clustr ( x, d, dev, b, f, e, i, j, n, nz, k )
 
c*********************************************************************72
c
cc CLUSTR uses the K-means algorithm to cluster data.
c
c  Discussion:
c
c    Given a matrix of I observations on J variables, the
c    observations are allocated to N clusters in such a way that the
c    within-cluster sum of squares is minimised.
c
c  Modified:
c
c    04 February 2008
c
c  Author:
c
c    David Sparks
c    Modifications by John Burkardt
c
c  Reference:
c
c    David Sparks,
c    Algorithm AS 58:
c    Euclidean Cluster Analysis,
c    Applied Statistics,
c    Volume 22, Number 1, 1973, pages 126-130.
c
c  Parameters:
c
c    Input, double precision X(I,J), the observed data.
c
c    Input/output, double precision D(K,J), the cluster centers.  On input,
c    the user has chosen these.  On output, they have been updated.
c
c    Output, double precision DEV(K), the sums of squared deviations of 
c    observations from their cluster centers.
c
c    Output, integer B(I), indicates the cluster to which each observation
c    has been assigned.
c
c    Workspace, double precision F(I).
c
c    Output, integer E(K), the number of observations assigned to each cluster.
c
c    Input, integer I, the number of observations.
c
c    Input, integer J, the number of variables.
c
c    Input, integer N, the number of clusters.
c
c    Input, integer NZ, the minimum number of observations which any cluster
c    is allowed to have.
c
c    Input, integer K, the maximum number of clusters.
c
      implicit none

      integer i
      integer k

      integer b(i)
      double precision big
      parameter ( big = 1.0D+10 )
      double precision d(k,j)
      double precision da
      double precision db
      double precision dc
      double precision de
      double precision dev(k)
      integer e(k)
      double precision f(i)
      double precision fl
      double precision fm
      double precision fq
      integer ia
      integer ic
      integer id
      integer ie
      integer ig
      integer ih
      integer ii
      integer ij
      integer ik
      integer il
      integer in
      integer ip
      integer ir
      integer is
      integer it
      integer iu
      integer iw
      integer ix
      integer iy
      integer j
      integer n
      integer nz
      double precision one
      parameter ( one = 1.0D+00 )
      double precision two
      parameter ( two = 2.0D+00 )
      double precision x(i,j)

      do ia = 1, n
        e(ia) = 0
      end do
c
c  For each observation, calculate the distance from each cluster
c  center, and assign to the nearest.
c
      do ic = 1, i
        f(ic) = 0.0D+00
        da = big
        do id = 1, n
          db = 0.0D+00
          do ie=1, j
            dc = x(ic,ie) - d(id,ie)
            db = db + dc * dc
            if ( da .le. db ) then
              go to 30
            end if
          end do
          da = db
          b(ic) = id
   30     continue
        end do
        ig = b(ic)
        e(ig) = e(ig) + 1
      end do
c
c  Calculate the mean and sum of squares for each cluster.
c
      do ix = 1, n
        dev(ix) = 0.0D+00
        do iy = 1, j
          d(ix,iy) = 0.0D+00
        end do
      end do

      do ic = 1, i
        ig = b(ic)
        do ih = 1, j
          d(ig,ih) = d(ig,ih) + x(ic,ih)
        end do
      end do

      do ij = 1, j
        do ii = 1, n
          if ( 0 .lt. e(ii) ) then
            d(ii,ij) = d(ii,ij) / real ( e(ii) )
          end if
        end do
      end do

      do ij = 1, j
        do ik = 1, i
          il = b(ik)
          da = x(ik,ij) - d(il,ij)
          db = da * da
          f(ik) = f(ik) + db
          dev(il) = dev(il) + db
        end do
      end do

      do ik = 1, i
        il = b(ik)
        fl = e(il)
        if ( 2.0D+00 .le. fl ) then
          f(ik) = f(ik) * fl / ( fl - one )
        end if
      end do
c
c  Examine each observation in turn to see if it should be
c  reassigned to a different cluster.
c
      if ( nz .le. 0 ) then
        nz = 1
      end if

   90 continue

      iw = 0

      do ik = 1, i

        il = b(ik)
        ir = il
c
c  If the number of cluster points is less than or equal to the
c  specified minimum, NZ, then bypass this iteration.
c
        if ( nz .lt. e(il) ) then

          fl = e(il)
          dc = f(ik)

          do in = 1, n
            if ( in .ne. il ) then
              fm = e(in)
              fm = fm / ( fm + one )
              de = 0.0D+00
              do ip = 1, j
                da = x(ik,ip) - d(in,ip)
                de = de + da * da * fm
                if ( dc .le. de ) then
                  go to 100
                end if
              end do
              dc = de
              ir = in
  100         continue
            end if
          end do
c
c  Reassignment is made here if necessary.
c
          if ( ir .ne. il ) then

            fq = e(ir)
            dev(il) = dev(il) - f(ik)
            dev(ir) = dev(ir) + dc
            e(ir) = e(ir) + 1
            e(il) = e(il) - 1

            do is = 1, j
              d(il,is) = ( d(il,is) * fl - x(ik,is) ) / ( fl - one )
              d(ir,is) = ( d(ir,is) * fq + x(ik,is) ) / ( fq + one )
            end do

            b(ik) = ir

            do it = 1, i

              ij = b(it)

              if ( ij .eq. il .or. ij .eq. ir ) then
                f(it) = 0.0D+00
                do iu = 1, j
                  da = x(it,iu) - d(ij,iu)
                  f(it) = f(it) + da * da
                end do
                fl = e(ij)
                f(it) = f(it) * fl / ( fl - one )
              end if

            end do

            iw = iw + 1

          end if

        end if

      end do
c
c  If any reassignments were made on this pass, then do another pass.
c
      if ( iw .ne. 0 ) then
        go to 90
      end if

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Discussion:
c
c    This FORTRAN77 version is made available for cases where the
c    FORTRAN90 version cannot be used.
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
