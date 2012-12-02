      program main

c*********************************************************************72
c
cc MAIN is the main program for TOMS571_PRB.
c
c  Modified:
c
c    30 December 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      write ( *, '(a)' ) ' '
      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS571_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS571 library.'

      call test01 ( )
      call test02 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS571_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tabulates VKAPPA and BESRAT.
c
c  Discussion:
c
c    This test code tabulates VKAPPA(R) and BESRAT(V) to 8d.
c
c    comparison with tables b,c to 5s in batschelet(1965), amer. inst.
c    biol. sciences monograph, or with tables to 5-6s in appendices 2.3,
c    2.2 of mardia (1972), statistics of directional data, academic
c    press, london and new york, or, for r up to 0.87, in table 2 of
c    gumbel et al. (1953) j. a. s. a. 48, reveals incorrect table
c    values for r = 0.94(0.01)0.99 in batschelet's
c    extension of gumbel's table and thence copied by mardia (app. 2.3).
c    incorrect values, r(12) = 0.95730, r(24) = 0.97937 and r(40) =
c    0.98739 occur in mardia's (appendix 2.2) extension of batschelet's
c    table c.
c
c    Using a variable precision version of besrat the actual number of
c    correct decimal digits for vkappa(r) is also tabulated to show
c    achievment of at least 8s precision.
c
c    The actual precision of besrat versions for a range of target
c    precisions is also tabulated to demonstrate validity of formulae
c    for forming data constants c1 and c2.
c
c  Modified:
c
c    30 December 2006
c
      implicit none

      real besra2
      real besrat
      real d
      real e
      real eps
      integer j
      integer k
      integer l
      integer n
      real r4_precis
      real t
      real u
      real vkappa
      real x(16)

      save eps

      data eps / 1.0E-06 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Validation of VKAPPA and BESRAT routines'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Table of VKAPPA(R) and BESRAT(V)'
      write ( *, '(a)' ) ' '

      n = 14
      k = n + 1
      write  ( *,99999) (j,j=5,n)

      do j = 1, 101

        x(1) = real ( j - 1 ) * 0.01E+00
        x(2) = vkappa ( x(1) )
        t = besra2 ( x(2), 16 )
        d = eps + x(2)
        u = besra2 ( d, 16 )
c
c  Linearly inverse interpolate to more precise VKAPPA value = E.
c
        if ( t - u .ne. 0.0E+00 ) then
          e = x(2) + ( t - x(1) ) / ( t - u ) * eps
        else
          e = x(2)
        end if

        x(3) = r4_precis ( x(2), e )

        if ( x(2) .ne. -1.0E+00 ) then
          x(2) = e
        end if
c
c  Generate arguments in x(4) to match table entries.
c
        if ( j .le. 80 ) then
          x(4) = real ( j - 1 ) * 0.1E+00
        else if ( j .lt. 92 ) then
          x(4) = real ( j - 41 ) * 0.2E+00
        else if ( j .eq. 92 ) then
          x(4) = 12.0E+00
        else if ( j .eq. 93 ) then
          x(4) = 15.0E+00
        else if ( j .gt. 93 ) then
          x(4) = 120.0E+00 / ( real ( 101 - j ) + eps * eps )
        end if

        t = besra2 ( x(4), 16 )
        x(5) = besrat ( x(4) )

        do l = 5, n
          x(l+1) = r4_precis ( besra2 ( x(4), l ), t )
        end do

        write  ( *,99998) (x(l),l=1,k)

      end do

      write  ( *,99997)
      return
99999 format (
     & 'compare',4x, 'table b,', 4x, 'batschelet(1965)', 4x, 'table c.'/
     & 8x,'appendix2.3', 6x, 'mardia(1972)', 6x, 'appendix 2.2',10x,
     & 'tabulated precision of besrat should exceed       '/
     & '       r',
     & '       vkappa(r)  precis', '       v       besrat(v)', 2x,
     & 11(i5, 's') )
99998 format (f8.3, f16.8, 2f8.2, f16.8, 4x, 11f6.2)
99997 format (14x, 'error flag'/
     &  '   actual precision .ge.', '    8.1',
     & '1', 30x, '5.07  6.08  7.08  8.07  9.10 10.12 11.17 12.20 13.4',
     & '5   -')
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tabulates SPHERR and CAPPA3.
c
c  Discussion:
c
c    validation code for functions SPHERR(cappa) and cappa3(r).
c
c    first demonstrate for precisions 5s - 20s that upper bound of
c    level of recursion for continued fraction in SPHERR achieves
c    at least target precision for argument up to 1.0 and SPHERR
c    precision over 14 decimals for 48b precision cdc7600.
c
c    secondly, demonstrate 48b precision performance of cappa3 and
c    the 8d, 16d, and 25d versions; checking error conditions and
c    precision for exact arguments near unity.
c
c  Modified:
c
c    30 December 2006
c
      implicit none

      double precision aa
      real b
      double precision bb
      real c
      double precision cc
      double precision d
      double precision dd
      double precision e
      integer j
      integer k
      double precision p(16)
      real r4_cappa3
      real r4_spherr
      double precision r8_cappa3
      double precision r8_precis
      double precision r8_spherr
      real x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Validation of SPHERR and CAPPA3 routines'
      write ( *, '(a)' ) ' '

      write ( *,99999) (k,k=5,20)
99999 format('test precision of continued fraction to level l =',
     * ' s/4(x+0.88)+2 for d =         '/7x,16i7,9x,'SPHERR')
      write ( *, '(a)' ) ' '

      do j = 1, 100

        x = real ( j ) * 0.01E+00
        d = 0.0E+00
        bb = x
        aa = r8_spherr ( bb, d )

        do k = 1, 16
          d = dble ( int ( ( x + 0.88D+00 ) * dble ( k + 4 )
     &      * 0.25D+00 + 2.0D+00 ) )
          cc = r8_spherr ( bb, d )
          p(k) = r8_precis ( cc, aa )
        end do

        c = r4_spherr ( x )
        cc = c
        d = r8_precis ( cc, aa )
        write ( *,99998) x, p, d

      end do

      write ( *, '(a)' ) ' '
      write ( *,99997) (j,j=1,4)
      write ( *, '(a)' ) ' '

      do j = 1, 105

        if ( j .lt. 102 ) then
          x = real ( j - 1 ) * 0.01E+00
        else if ( j .eq. 102 ) then
          x = 1023.0E+00 / 1024.0E+00
        else if ( j .eq. 103 ) then
          x = 0.9999E+00
        else if ( j .eq. 104 ) then
          x = 0.0001E+00
        else if ( j .eq. 105 ) then
          x = -0.0001E+00
        end if

        c = r4_cappa3 ( x )

        if ( c .ne. -1.0E+00 ) then

          bb = x
          aa = r8_cappa3 ( bb, 30.0D+00 )
          dd = c
          p(1) = r8_precis ( aa, dd )
          b = r4_spherr ( c)
          cc = r8_spherr ( dd, 0.0D+00 )
          dd = b
          e = r8_precis ( dd, cc )

          do k = 1, 4
            d = dble ( k - 1 ) * 8.0D+00
            dd = r8_cappa3 ( bb, d )
            p(k+1) = r8_precis ( aa, dd )
          end do

        end if

        write ( *,99996) x, c, (p(k),k=1,5), b, e

      end do

      return
99998 format (f5.2, 5x, 16f7.2, f12.2)
99997 format ('test precision of cappa3 and SPHERR procedures'
     & /'0   x', 16x, 'cappa3  its precis', 4x, 'approx', i2, 3i7,
     & '   SPHERR(cappa3(x))', '  its precis'/)
99996 format (f7.4, f20.10, 2f12.2, 3f7.2, f20.10, f12.2)
      end
