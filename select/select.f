cccc from n&w's book, p.108
cccc this is typed by lei zhao on 07/12/95.

      subroutine select ( family, task, n, k, mu, nu, edge, m, 
     &  newone, rank, b )

c*********************************************************************72
c
cc SELECT ...
c
c  Discussion:
c
c    The SELECT program seems to write data into entries (N+1) through (N+K)
c    of EDGE, MU and NU, beyond the original formal dimension of N.  This extra data
c    is of no interest to the user, but it must be allocated in the
c    calling program to avoid memory overwriting!
c
c  Modified:
c
c    02 December 2008
c
      implicit none

      integer k
      integer n

      integer b(10,10)
      integer b1
      integer b2
      integer bnew
      integer edge(n+k)
      integer flast
      integer family
      integer j
      integer k1
      integer klast
      integer l
      integer m
      integer m1
      integer mu(n+k)
      logical newone
      integer nlast
      integer nu(n+k)
      integer phi
      integer psi
      real rand
      integer rank
      integer rankp
      integer t
      integer t1
      integer task
      integer xnew

      data flast / 0 /
      data klast / 0 /
      data nlast / 0 /

      if (n .le. nlast.and.k.le.klast.and.family.eq.flast) go to 205

      nlast = n
      klast = k
      flast = family

      do m1 = 1, n
        do k1 = 1, k
          b(m1,k1)  = phi(m1,k1,family)*bnew(m1,k1,family,1,b)
     &      +psi(m1,k1,family)*bnew(m1,k1,family,2,b)
        end do
      end do

205   go to (100, 400, 200, 300), task

300   continue

      rank = rand(1) * b(n,k)

200   continue

      j = 1
      rankp = rank
      mu(j) = n
      nu(j) = k
      m = j

510   continue

      t1 = phi(mu(m), nu(m), family)

      if ( t1 + psi ( mu(m), nu(m), family ) .eq. 0 ) then
        m = m - 1
        return
      end if

518   continue

      b1 = bnew(mu(m), nu(m), family,1,b)
      if ( rankp .ge. t1*b1) go to 520
512   edge(m) = rankp/b1
      rankp = rankp-edge(m)*b1
      mu(m+1) = xnew(mu(m), nu(m), family, 1)
      nu(m+1) = nu(m)
515   m = m+1
      go to 510

520   rankp = rankp-b1*t1
      b2 = bnew(mu(m), nu(m), family, 2, b)
      edge(m) = t1+rankp/b2
      rankp = rankp-(edge(m)-t1)*b2
      mu(m+1) = xnew(mu(m), nu(m), family,2)
      nu(m+1) = nu(m)-1
      go to 515

400   continue

      rank = 0
      m1 = m-1

      do j = 1, m1

        if ( nu(j+1) .eq. nu(j) ) then

          rank = rank+edge(j)*bnew(mu(j),nu(j),family,1,b)

        else

          rank = rank+phi(mu(j),nu(j),family)
     &      *bnew(mu(j),nu(j),family,1,b)
     *      +(edge(j)-phi(mu(j),nu(j),family))
     &      *bnew(mu(j),nu(j),family,2,b)

        end if

      end do

      return

100   continue

      if ( .not. newone ) then
        newone = .true.
        rank = 0
        go to 200
      end if

105   j = m
120   continue

      t = phi(mu(j),nu(j),family)

      if (edge(j).lt.t+psi(mu(j),nu(j),family)-1) go to 130
      j = j - 1

      if ( j .eq. 0 ) then
        newone = .false.
        return
      end if

      go to 120

130   edge(j) = edge(j)+1
      l = j + 1

      if ( edge(j) .eq. t ) then
        nu(l) = nu(l)-1
        mu(l) = xnew(mu(l-1), nu(l-1), family, 2)
      end if

140   continue

      rankp = 0
      m = l
      go to 510

      end
      function xnew ( m1, m2, mf, mgo )

c*********************************************************************72
c
cc XNEW ???
c
c  Modified:
c
c    02 December 2008
c
      implicit none

      integer m1
      integer m2
      integer mf
      integer mgo
      integer xnew

      xnew = m1-1
      if(mf.eq.6.and.mgo.eq.1) xnew = m1-m2
      if(mf.eq.7.and.mgo.eq.2) xnew = m1

      return
      end
      function bnew ( m1, m2, mf, mgo, b )

c*********************************************************************72
c
cc BNEW ???
c
c  Modified:
c
c    02 December 2008
c
      implicit none

      integer b(10,10)
      integer bnew
      integer lx
      integer ly
      integer m1
      integer m2
      integer mf
      integer mgo
      integer xnew

      lx = xnew ( m1, m2, mf, mgo )
      ly = m2 - mgo + 1

      if ( 0 .lt. lx * ly ) then
        bnew = b(lx,ly)
      else
        bnew = 1
      end if

      return
      end
      function phi ( mu, nu, family )

c*********************************************************************72
c
cc PHI ???
c
c  Modified:
c
c    02 December 2008
c
c  Parameters:
c
c    Input, integer MU, ?
c
c    Input, integer NU, ?
c
c    Input, integer FAMILY, ?
c
c    Output, integer PHI, ?
c
      implicit none

      integer family
      integer mu
      integer nu
      integer phi
      integer q

      save q

      data q / 2 /

      phi = 0

      if ( family .eq. 1 ) then
        if (mu.ge.nu+1.and.nu.ge.0) phi = 1
      else if ( family .eq. 2 ) then
        if(mu.ge.nu+1.and.nu.ge.1) phi = nu
      else if ( family .eq. 3 ) then
        if (mu.ge.nu+1.and.nu.ge.0) phi = mu-1
      else if ( family .eq. 4 ) then
        if (mu.ge.nu+1.and.nu.ge.0) phi = q**nu
      else if ( family .eq. 5 ) then
        if (mu.ge.nu+1.and.nu.ge.1) phi = nu
        if (mu.eq.1.and.nu.eq.1) phi = 1
      else if ( family .eq. 6 ) then
        if (mu.ge.2*nu.and.nu.ge.1) phi = 1
      else if ( family .eq. 7 ) then
        if (mu.ge.1.and.nu.ge.1) phi = 1
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'PHI - Fatal error!'
        write ( *, '(a)' ) '  Illegal value of FAMILY.'
        stop
      end if

      return
      end
      function psi ( mu, nu, family )

c*********************************************************************72
c
cc PSI ???
c
c  Modified:
c
c    02 December 2008
c
c  Parameters:
c
c    Input, integer MU, ?
c
c    Input, integer NU, ?
c
c    Input, integer FAMILY, ?
c
c    Output, integer PSI, ?
c
      implicit none

      integer family
      integer mu
      integer nu
      integer psi

      psi = 0
      go to (10,20,20,10,30,40,50), family
10    if (mu.ge.nu.and.nu.ge.1) psi = 1
      return
20    if ((mu.ge.nu.and.nu.ge.2) .or. (mu.eq.1.and.nu.eq.1)) psi = 1
      return
30    if (mu.ge.nu.and.nu.ge.2) psi = mu-nu+1
      return
40    if ((mu.ge.nu.and.nu.ge.2).or.(mu.eq.1.and.nu.eq.1)) psi = 1
      return
50    if (mu.ge.0.and.nu.ge.2) psi = 1

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
