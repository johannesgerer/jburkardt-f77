      program main

c*********************************************************************72
c
cc MAIN is the main program for testing TESTPACK.
c
c  Discussion:
c
c     MULTST Test Package with example test program and early ADAPT
c
c  Modified:
c
c    13 March 2007
c
c  Author:
c
c    Alan Genz
c
      implicit none

      integer ndim1
      parameter ( ndim1 = 5 )
      integer tstlim
      parameter ( tstlim = 6 )
      integer tstmax
      parameter ( tstmax = 6 )

      external adapt
      character*6 sbname
      double precision difclt(tstmax)
      double precision expnts(tstmax)
      integer i
      integer maxcls
      integer ndims(ndim1)
      integer nsamp
      double precision rel_tol
      integer tstfns(tstlim)

      data difclt / 
     &  110.0D+00, 600.0D+00, 600.0D+00, 100.0D+00, 150.0D+00, 
     &  100.0D+00 /
      data expnts / 
     &  1.5D+00, 2.0D+00, 2.0D+00, 1.0D+00, 2.0D+00, 2.0D+00 /
      data tstfns / 1, 2, 3, 4, 5, 6 /
      data ndims / 2, 3, 4, 6, 8 /

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TESTPACK'
      write ( *, '(a)' ) '  FORTRAN77 version'

      nsamp = 20
      sbname = ' ADAPT'
      rel_tol = 1.0D-06
      maxcls = 10000

      call multst ( nsamp, tstlim, tstfns, tstmax, difclt,
     &  expnts, ndim1, ndims, sbname, adapt, rel_tol, maxcls )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TESTPACK'
      write ( *, '(a)' ) '  Normal end of execution'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine adapt ( ndim, a, b, minpts, maxpts, functn, rel_tol,
     &  relerr, lenwrk, wrkstr, finest, ifail )

c*********************************************************************72
c
cc ADAPT carries out adaptive multidimensional quadrature.
c
c  Modified:
c
c    13 March 2007
c
c  Author:
c
c    Alan Genz
c
c  Parameters:
c
c    Input, integer NDIM, the number of variables.
c    2 <= NDIM <= 100.
c
c    Input, double precision A(NDIM), the lower limits of integration.
c
c    Input, double precision B(NDIM), the upper limits of integration.
c
c    Input/output, integer MINPTS.  On input, the minimum number of function 
c    evaluations to be allowed,  MINPTS must not exceed MAXPTS.  If MINPTS < 0
c    then the routine assumes a previous call has been made with the same 
c    integrand and continues that calculation.
c    On output, the actual number of function evaluations used.
c
c    Input, integer MAXPTS.  On input, the maximum number of function
c    evaluations allowed, which must be at least RULCLS, where
c    RULCLS = 2**NDIM + 2 * NDIM**2 + 2 * NDIM + 1, when NDIM < 16 and
c    RULCLS = ( NDIM * ( 14 - NDIM * ( 6 - 4 * NDIM ) ) / 3 + 1, when NDIM > 15.
c    for NDIM  =  2   3   4   5   6   7   8   9   10   11   12
c    RULCLS   =  17  33  57  93 149 241 401 693 1245 2313 4409
c    A suggested starting value for MAXPTS is 100*RULCLS.  If
c    this is not large enough for the required accuracy, then
c    MAXPTS (and LENWRK) should be increased accordingly.
c
c    Input, external FUNCTN, the user defined function to be integrated.
c    It must have the form
c      subroutine functn ( indx, ndim, z, alpha, beta, f )
c    where
c      INDX is the index of the test function,
c      NDIM is the spatial dimension,
c      Z is the evaluation point,
c      ALPHA is a parameter vector,
c      BETA is a parameter vector,
c      F is the function value.
c
c    Input, double precision REL_TOL, the user's requested relative accuracy.
c
c    Input, integer LENWRK, the length of the array WRKSTR.
c    The routine needs (2*NDIM+3)*(1+MAXPTS/RULCLS)/2 for LENWRK if
c    MAXPTS function calls are used.
c
c    Workspace, double precision WRKSTR(LENWRK).
c
c    Output, double precision RELERR, the estimated relative accuracy of FINEST
c
c    Output, double precision FINEST, the estimated value of integral.
c
c    Output, integer IFAIL
c    * 0, for normal exit, when estimated relative error RELERR is less
c    than REL_TOL, and with MAXPTS or less function calls made.
c    * 1, if MAXPTS was too small for ADAPT to obtain the required relative
c    error REL_TOL.  In this case ADAPT returns a value of FINEST with 
c    estimated relative error RELERR.
c    * 2, if LENWRK was too small for MAXPTS function calls.  In
c    this case ADAPT returns a value of FINEST with estimated error
c    RELERR using the working storage available, but RELERR is likely to
c    be greater than REL_TOL.
c    * 3, if NDIM < 2 or 100 < NDIM or MAXPTS < MINPTS or MAXPTS < RULCLS.
c
      implicit none

      integer lenwrk
      integer ndim

      double precision a(ndim)
      double precision b(ndim)
      double precision center(100)
      double precision df1
      double precision df2
      double precision dif
      double precision difmax
      integer divaxn
      integer divaxo
      integer divflg
      double precision f1
      double precision f2
      double precision f3
      double precision f4
      double precision finest
      integer funcls
      external functn
      integer i
      integer ifail
      integer index1
      integer index2
      integer j
      integer k
      integer l
      double precision lamda2
      double precision lamda4
      double precision lamda5
      integer m
      integer maxpts
      integer minpts
      integer n
      double precision ratio
      double precision rel_tol
      double precision relerr
      double precision rgncmp
      double precision rgnerr
      integer rgnstr
      double precision rgnval
      double precision rgnvol
      integer rulcls
      integer sbrgns
      integer sbtmpp
      integer subrgn
      integer subtmp
      double precision sum1
      double precision sum2
      double precision sum3
      double precision sum4
      double precision sum5
      double precision twondm
      double precision weit1
      double precision weit2
      double precision weit3
      double precision weit4
      double precision weit5
      double precision weitp1
      double precision weitp2
      double precision weitp3
      double precision weitp4
      double precision width(100)
      double precision widthl(100)
      double precision wrkstr(lenwrk)
      double precision z(100)

      save rgnstr
      save sbrgns

      ifail = 3
      relerr = 1.0D+00
      funcls = 0

      if ( ndim .lt. 2 .or. 100 .lt. ndim ) then
        minpts = 0
        wrkstr(lenwrk-1) = sbrgns
        return
      end if

      if ( maxpts .lt. minpts ) then
        minpts = 0
        wrkstr(lenwrk-1) = sbrgns
        return
      end if
c
c  Initialisation.
c
      twondm = 2.0D+00**ndim
      rgnstr = 2 * ndim + 3
      divaxo = 0
c
c  Basic rule initialisation.
c
      lamda5 = 9.0D+00 / 19.0D+00

      if ( ndim .le. 15 ) then

        rulcls = 2**ndim + 2 * ndim * ndim + 2 * ndim + 1
        lamda4 = 9.0D+00 / 10.0D+00
        lamda2 = 9.0D+00 / 70.0D+00
        weit5 = 1.0D+00 / ( 3.0D+00 * lamda5 )**3 / twondm

      else

        rulcls = 1 + ( ndim * ( 12 + ( ndim - 1 ) 
     &    * ( 6 + ( ndim - 2 ) * 4 ) ) ) / 3

        ratio = dble ( ndim - 2 ) / 9.0D+00

        lamda4 = ( 1.0D+00 / 5.0D+00 - ratio ) 
     &    / ( 1.0D+00 / 3.0D+00 - ratio / lamda5 )

        ratio = ( 1.0D+00 - lamda4 / lamda5 ) 
     &    * dble ( ndim - 1 ) * ratio / 6.0D+00

        lamda2 = ( 1.0D+00 / 7.0D+00 - lamda4 / 5.0D+00 - ratio )
     &    / ( 1.0D+00 / 5.0D+00 - lamda4 / 3.0D+00 - ratio / lamda5 )

        weit5 = 1.0D+00 / ( 6.0D+00 * lamda5 )**3

      end if

      weit4 = ( 1.0D+00 / 15.0D+00 - lamda5 / 9.0D+00 ) 
     &  / ( 4.0D+00 * ( lamda4 - lamda5 ) * lamda4**2 )

      weit3 = ( 1.0D+00 / 7.0D+00 - ( lamda5 + lamda2 ) / 5.0D+00
     &  + lamda5 * lamda2 / 3.0D+00 ) / ( 2.0D+00 * lamda4 
     &  * ( lamda4 - lamda5 ) * ( lamda4 - lamda2 ) ) 
     &  - 2.0D+00 * dble ( ndim - 1 ) * weit4

      weit2 = ( 1.0D+00 / 7.0D+00 - ( lamda5 + lamda4 ) / 5.0D+00
     &  + lamda5 * lamda4 / 3.0D+00 ) / ( 2.0D+00 * lamda2 
     &  * ( lamda2 - lamda5 ) * ( lamda2 - lamda4 ) )

      if ( ndim .le. 15 ) then
        weit1 = 1.0D+00 - 2.0D+00 * dble ( ndim )
     &    * ( weit2 + weit3 + dble ( ndim - 1 ) * weit4 ) 
     &    - twondm * weit5
      else
        weit1 = 1.0D+00 - 2.0D+00 * dble ( ndim ) 
     &    * ( weit2 + weit3 + dble ( ndim - 1 ) *
     &    ( weit4 + 2.0D+00 * dble ( ndim - 2 ) * weit5 / 3.0D+00 ) )
      end if

      weitp4 = 1.0D+00 / ( 6.0D+00 * lamda4 )**2

      weitp3 = ( 1.0D+00 / 5.0D+00 - lamda2 / 3.0D+00 ) / 
     &  ( 2.0D+00 * lamda4 * ( lamda4 - lamda2 ) ) 
     &  - 2.0D+00 * dble ( ndim - 1 ) * weitp4

      weitp2 = ( 1.0D+00 / 5.0D+00 - lamda4 / 3.0D+00 ) 
     &  / ( 2.0D+00 * lamda2 * ( lamda2 - lamda4 ) )

      weitp1 = 1.0D+00 - 2.0D+00 * dble ( ndim ) * 
     &  ( weitp2 + weitp3 + dble ( ndim - 1 ) * weitp4 )

      ratio = lamda2 / lamda4

      lamda5 = sqrt ( lamda5 )
      lamda4 = sqrt ( lamda4 )
      lamda2 = sqrt ( lamda2 )

      if ( maxpts .lt. rulcls ) then
        return
      end if
c
c  End basic rule initialisation.
c
      if ( minpts .lt. 0 ) then

        sbrgns = int ( wrkstr(lenwrk-1) )
        divflg = 0
        subrgn = rgnstr
        wrkstr(lenwrk) = wrkstr(lenwrk) - wrkstr(subrgn)
        finest = finest - wrkstr(subrgn-1)
        divaxo = int ( wrkstr(subrgn-2) )

        do j = 1, ndim
          subtmp = subrgn - 2 * ( j + 1 )
          center(j) = wrkstr(subtmp+1)
          width(j) = wrkstr(subtmp)
        end do

        width(divaxo) = width(divaxo) / 2.0D+00
        center(divaxo) = center(divaxo) - width(divaxo)

      else

        do j = 1, ndim
          width(j) = ( b(j) - a(j) ) / 2.0D+00
          center(j) = a(j) + width(j)
        end do

        finest = 0.0D+00
        wrkstr(lenwrk) = 0.0D+00
        divflg = 1
        subrgn = rgnstr
        sbrgns = rgnstr

      end if
c
c  Begin basic rule.
c
10    continue

      rgnvol = twondm
      do j = 1, ndim
        rgnvol = rgnvol * width(j)
      end do

      do j = 1, ndim
        z(j) = center(j)
      end do

      call functn ( ndim, z, sum1 )
c
c  Compute symmetric sums of functn(lamda2,0,0,...,0) and
c  functn(lamda4,0,0,...,0), and maximum fourth difference.
c
      difmax = -1.0D+00
      sum2 = 0.0D+00
      sum3 = 0.0D+00

      do j = 1, ndim

        z(j) = center(j) - lamda2 * width(j)
        call functn ( ndim, z, f1 )
        z(j) = center(j) + lamda2 * width(j)
        call functn ( ndim, z, f2 )
        widthl(j) = lamda4 * width(j)
        z(j) = center(j) - widthl(j)
        call functn ( ndim, z, f3 )
        z(j) = center(j) + widthl(j)
        call functn ( ndim, z, f4 )
        sum2 = sum2 + f1 + f2
        sum3 = sum3 + f3 + f4
        df1 = f1 + f2 - 2.0D+00 * sum1
        df2 = f3 + f4 - 2.0D+00 * sum1
        dif = abs ( df1 - ratio * df2 )

        if ( difmax .lt. dif ) then
          difmax = dif
          divaxn = j
        end if

        z(j) = center(j)

      end do

      if ( sum1 .eq. sum1 + difmax / 8.0D+00 ) then
        divaxn = mod ( divaxo, ndim ) + 1
      end if
c
c  Compute symmetric sum of functn(lamda4,lamda4,0,0,...,0).
c
      sum4 = 0.0D+00

      do j = 2, ndim

        do k = j, ndim

          do l = 1, 2

            widthl(j-1) = -widthl(j-1)
            z(j-1) = center(j-1) + widthl(j-1)

            do m = 1, 2
              widthl(k) = -widthl(k)
              z(k) = center(k) + widthl(k)
              call functn ( ndim, z, f1 )
              sum4 = sum4 + f1
            end do

          end do

          z(k) = center(k)

        end do

        z(j-1) = center(j-1)

      end do
c
c  If NDIM < 16 compute symmetric sum of functn(lamda5,lamda5,...,lamda5).
c
      sum5 = 0.0D+00

      if ( ndim .le. 15 ) then

        do j = 1, ndim
          widthl(j) = -lamda5 * width(j)
        end do

        do j = 1, ndim
          z(j) = center(j) + widthl(j)
        end do

20      continue

        call functn ( ndim, z, f1 )
        sum5 = sum5 + f1

        do j = 1, ndim

          widthl(j) = -widthl(j)
          z(j) = center(j) + widthl(j)

          if ( widthl(j) .gt. 0.0D+00 ) then
            go to 20
          end if

        end do
c
c  If 15 < NDIM, compute symmetric sum of
c  FUNCTN(LAMDA5,LAMDA5,LAMDA5,0,0,...,0).
c
      else

        do j = 1, ndim
          widthl(j) = lamda5 * width(j)
        end do

        do i = 3, ndim
          do j = i, ndim
            do k = j, ndim

              do l = 1, 2
                widthl(i-2) = -widthl(i-2)
                z(i-2) = center(i-2) + widthl(i-2)
                do m = 1, 2
                  widthl(j-1) = -widthl(j-1)
                  z(j-1) = center(j-1) + widthl(j-1)
                  do n = 1, 2
                    widthl(k) = -widthl(k)
                    z(k) = center(k) + widthl(k)
                    call functn ( ndim, z, f1 )
                    sum5 = sum5 + f1
                  end do
                end do
              end do

              z(k) = center(k)

            end do

            z(j-1) = center(j-1)

          end do

          z(i-2) = center(i-2)

        end do

      end if
c
c  Compute fifth and seventh degree rules and error.
c
      rgncmp = rgnvol * ( weitp1 * sum1 
     &                  + weitp2 * sum2
     &                  + weitp3 * sum3
     &                  + weitp4 * sum4 )

      rgnval = rgnvol * ( weit1 * sum1
     &                  + weit2 * sum2
     &                  + weit3 * sum3
     &                  + weit4 * sum4
     &                  + weit5 * sum5 )

      rgnerr = abs ( rgnval - rgncmp )
c
c  End basic rule.
c
      finest = finest + rgnval
      wrkstr(lenwrk) = wrkstr(lenwrk) + rgnerr
      funcls = funcls + rulcls
c
c  Place results of basic rule into partially ordered list
c  according to subregion error.
c
c  When DIVFLG = 0 start at top of list and move down list tree to
c  find correct position for results from first half of recently
c  divided subregion.
c
      if ( divflg .ne. 1 ) then

30      continue

        subtmp = 2 * subrgn

        if ( subtmp .gt. sbrgns ) then
          go to 50
        end if

        if ( subtmp .ne. sbrgns ) then
          sbtmpp = subtmp + rgnstr
          if ( wrkstr(subtmp) .lt. wrkstr(sbtmpp) ) then
            subtmp = sbtmpp
          end if
        end if

        if ( rgnerr .ge. wrkstr(subtmp) ) then
          go to 50
        end if

        do k = 1, rgnstr
          wrkstr(subrgn-k+1) = wrkstr(subtmp-k+1)
        end do

        subrgn = subtmp
        go to 30

      end if
c
c  When DIVFLG = 1 start at bottom right branch and move up list
c  tree to find correct position for results from second half of
c  recently divided subregion.
c
40    continue

      subtmp = ( subrgn / ( 2 * rgnstr ) ) * rgnstr

      if ( subtmp .lt. rgnstr ) then
        go to 50
      end if

      if ( rgnerr .le. wrkstr(subtmp) ) then
        go to 50
      end if

      do k = 1, rgnstr
        index1 = subrgn - k + 1
        index2 = subtmp - k + 1
        wrkstr(index1) = wrkstr(index2)
      end do

      subrgn = subtmp

      go to 40
c
c  Store results of basic rule in correct position in list.
c
50    continue

      wrkstr(subrgn) = rgnerr
      wrkstr(subrgn-1) = rgnval
      wrkstr(subrgn-2) = divaxn

      do j = 1, ndim
        subtmp = subrgn - 2 * ( j + 1 )
        wrkstr(subtmp+1) = center(j)
        wrkstr(subtmp) = width(j)
      end do
c
c  When DIVFLG = 0 prepare for second application of basic rule.
c
      if ( divflg .ne. 1 ) then
        center(divaxo) = center(divaxo) + 2.0D+00 * width(divaxo)
        sbrgns = sbrgns + rgnstr
        subrgn = sbrgns
        divflg = 1
        go to 10
      end if
c
c  End ordering and storage of basic rule results.
c  Make checks for possible termination of routine.
c
      relerr = 1.0D+00

      if ( wrkstr(lenwrk) .le. 0.0D+00 ) then
        wrkstr(lenwrk) = 0.0D+00
      end if

      if ( abs ( finest ) .ne. 0.0D+00 ) then
        relerr = wrkstr(lenwrk) / abs ( finest )
      end if

      if ( 1.0D+00 .lt. relerr ) then
        relerr = 1.0D+00
      end if

      if ( lenwrk .lt. sbrgns + rgnstr + 2 ) then
        ifail = 2
      end if

      if ( maxpts .lt. funcls + 2 * rulcls ) then
        ifail = 1
      end if

      if ( relerr .lt. rel_tol .and. minpts .le. funcls ) then
        ifail = 0
      end if

      if ( ifail .lt. 3 ) then
        minpts = funcls
        wrkstr(lenwrk-1) = sbrgns
        return
      end if
c
c  Prepare to use basic rule on each half of subregion with largest
c  error.
c
      divflg = 0
      subrgn = rgnstr
      wrkstr(lenwrk) = wrkstr(lenwrk) - wrkstr(subrgn)
      finest = finest - wrkstr(subrgn-1)
      divaxo = int ( wrkstr(subrgn-2) )

      do j = 1, ndim
        subtmp = subrgn - 2 * ( j + 1 )
        center(j) = wrkstr(subtmp+1)
        width(j) = wrkstr(subtmp)
      end do

      width(divaxo) = width(divaxo) / 2.0D+00
      center(divaxo) = center(divaxo) - width(divaxo)
c
c  Loop back to apply basic rule.
c
      go to 10

      end
      subroutine genz_function ( ndim, z, f )

c*********************************************************************72
c
cc GENZ_FUNCTION evaluates one of the test integrand functions.
c
c  Modified:
c
c    31 May 2007
c
c  Author:
c
c    Alan Genz
c
c  Reference:
c
c    Alan Genz,
c    A Package for Testing Multiple Integration Subroutines,
c    in Numerical Integration: 
c    Recent Developments, Software and Applications,
c    edited by Patrick Keast, Graeme Fairweather,
c    D Reidel, 1987, pages 337-340,
c    LC: QA299.3.N38.
c
c  Parameters:
c
c    Input, integer NDIM, the spatial dimension.
c
c    Input, double precision Z(NDIM), the point at which the integrand 
c    is to be evaluated.
c
c    Output, double precision F, the value of the test function .
c
      implicit none

      integer maxdim
      parameter ( maxdim = 20 )
      integer ndim

      double precision alpha(maxdim)
      double precision beta(maxdim)
      double precision dfclt
      double precision exn
      double precision expmax
      parameter ( expmax = 100.0D+00 )
      double precision f
      integer itst
      integer j
      double precision pi
      parameter ( pi = 3.14159265358979323844D+00 )
      logical test
      double precision total
      double precision value
      double precision z(ndim)

      common /tstblk/ itst, alpha, beta, dfclt, exn

      value = 0.0D+00
c
c  Oscillatory.
c
      if ( itst .eq. 1 ) then

        total = 2.0D+00 * pi * beta(1)
        do j = 1, ndim
          total = total + z(j)
        end do
        value = cos ( total )
c
c  Product Peak.
c
      else if ( itst .eq. 2 ) then

        total = 1.0D+00
        do j = 1, ndim
          total = total 
     &      / ( 1.0D+00 / alpha(j)**2 + ( z(j) - beta(j) )**2 )
        end do
        value = total
c
c  Corner Peak.
c
      else if ( itst .eq. 3 ) then
c
c  For this case, the BETA's are used to randomly select
c  a corner for the peak.
c
        total = 1.0D+00
        do j = 1, ndim
          if ( beta(j) .lt. 0.5D+00 ) then
            total = total + z(j)
          else
            total = total + alpha(j) - z(j)
          end if
        end do
        value = 1.0D+00 / total**( ndim + 1 )
c
c  Gaussian.
c
      else if ( itst .eq. 4 ) then

        total = 0.0D+00
        do j = 1, ndim
          total = total + ( alpha(j) * ( z(j) - beta(j) ) )**2
        end do

        if ( total .lt. expmax ) then
          value = exp ( - total )
        end if
c
c  C0 Function.
c
      else if ( itst .eq. 5 ) then

        total = 0.0D+00
        do j = 1, ndim
          total = total + alpha(j) * abs ( z(j) - beta(j) )
        end do

        if ( total .lt. expmax ) then
          value = exp ( - total )
        end if
c
c  Discontinuous.
c
      else if ( itst .eq. 6 ) then

        test = .false.
        do j = 1, ndim

          if ( beta(j) .lt. z(j) ) then
            test = .true.
          end if

        end do

        if ( test ) then

          value = 0.0D+00

        else

          total = 0.0D+00
          do j = 1, ndim
            total = total + alpha(j) * z(j)
          end do

          value = exp ( total )

        end if

      end if

      f = value

      return
      end
      function genz_integral ( indx, ndim, a, b )

c*********************************************************************72
c
cc GENZ_INTEGRAL computes the exact integrals of the test functions.
c
c  Modified:
c
c    26 May 2007
c
c  Author:
c
c    Alan Genz
c
c  Reference:
c
c    Alan Genz,
c    A Package for Testing Multiple Integration Subroutines,
c    in Numerical Integration: 
c    Recent Developments, Software and Applications,
c    edited by Patrick Keast, Graeme Fairweather,
c    D Reidel, 1987, pages 337-340,
c    LC: QA299.3.N38.
c
c  Parameters:
c
c    Input, integer ITEST, the index of the test.
c
c    Input, integer NDIM, the spatial dimension.
c
c    Input, double precision A(NDIM), B(NDIM), the lower and upper limits 
c    of integration.
c
c    Output, double precision GENZ_INTEGRAL, the exact value of the integral.
c
      implicit none

      integer maxdim
      parameter ( maxdim = 20 )
      integer ndim

      double precision a(ndim)
      double precision ab
      double precision alpha(maxdim)
      double precision b(ndim)
      double precision beta(maxdim)
      double precision dfact
      double precision dfclt
      double precision exn
      double precision genz_integral
      double precision genz_phi
      integer ic(maxdim)
      integer indx
      integer isum
      integer itst
      integer j
      double precision pi
      parameter ( pi = 3.14159265358979323844D+00 )
      double precision sgndm
      double precision sign
      double precision total
      double precision value

      common /tstblk/ itst, alpha, beta, dfclt, exn

      total = 0.0D+00
      do j = 1, ndim
        total = total + alpha(j)
      end do

      dfact = total * dble ( ndim )**exn / dfclt
      do j = 1, ndim
        alpha(j) = alpha(j) / dfact
      end do
c
c  Oscillatory.
c
      if ( indx .eq. 1 ) then

        do j = 1, ndim
          b(j) = alpha(j)
        end do

        do j = 1, ndim
          ic(j) = 0
        end do

        value = 0.0D+00

 10     continue

        isum = 0
        total = 2.0D+00 * pi * beta(1)
        do j = 1, ndim
          if ( ic(j) .ne. 1 ) then
            total = total + alpha(j)
          end if
          isum = isum + ic(j)
        end do

        sign = 1 + 2 * ( ( isum / 2 ) * 2 - isum )

        if ( mod ( ndim, 2 ) .eq. 0 ) then
          value = value + sign * cos ( total )
        else
          value = value + sign * sin ( total ) 
        end if

        do j = 1, ndim
          ic(j) = ic(j) + 1
          if ( ic(j) .lt. 2 ) then
            go to 10
          end if
          ic(j) = 0
        end do

        value = value * dble ( (-1)**(ndim/2) )
c
c  Product Peak.
c
      else if ( indx .eq. 2 ) then

        value = 1.0D+00
        do j = 1, ndim
          value = value * alpha(j) * ( 
     &        atan ( ( 1.0D+00 - beta(j) ) * alpha(j) ) 
     &      + atan (           + beta(j)   * alpha(j) ) )
        end do
c
c  Corner Peak.
c
      else if ( indx .eq. 3 ) then

        value = 0.0D+00
        sgndm = 1.0D+00
        do j = 1, ndim
          sgndm = - sgndm / dble ( j )
          b(j) = alpha(j)
          ic(j) = 0
        end do

 20     continue

        total = 1.0D+00
        isum = 0

        do j = 1, ndim
          if ( ic(j) .ne. 1 ) then
            total = total + alpha(j)
          end if
          isum = isum + ic(j)
        end do

        sign = 1 + 2 * ( ( isum / 2 ) * 2 - isum )
        value = value + sign / total
        do j = 1, ndim
          ic(j) = ic(j) + 1
          if ( ic(j) .lt. 2 ) then
            go to 20
          end if
          ic(j) = 0
        end do
        value = value * sgndm
c
c  Gaussian.
c
      else if ( indx .eq. 4 ) then

        value = 1.0D+00
        ab = sqrt ( 2.0D+00 )
        do j = 1, ndim
          value = value * ( sqrt ( pi ) / alpha(j) ) *
     &      (   genz_phi ( ( 1.0D+00 - beta(j) ) * ab * alpha(j) )
     &        - genz_phi (           - beta(j)   * ab * alpha(j) ) )
        end do
c
c  C0 Function.
c
      else if ( indx .eq. 5 ) then

        value = 1.0D+00
        do j = 1, ndim
          ab = alpha(j) * beta(j)
          value = value * 
     &      ( 2.0D+00 - exp ( - ab ) - exp ( ab - alpha(j) ) ) 
     &      / alpha(j)
        end do
c
c  Discontinuous.
c
      else if ( indx .eq. 6 ) then

        value = 1.0D+00
        do j = 1, ndim
          if ( 2 .lt. j ) then
            beta(j) = 1
          end if
          value = value * ( exp ( alpha(j) * beta(j) ) - 1.0D+00 ) 
     &      / alpha(j)
        end do

      end if

      genz_integral = value

      return
      end
      function genz_name ( indx )

c*********************************************************************72
c
cc GENZ_NAME returns the name of a Genz test integrand.
c
c  Modified:
c
c    26 May 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer INDX, the index of the test integrand.
c
c    Output, character*13 GENZ_NAME, the name of the test integrand.
c
      implicit none

      integer indx
      character*13 genz_name

      if ( indx == 1 ) then
        genz_name = 'Oscillatory  '
      else if ( indx == 2 ) then
        genz_name = 'Product Peak '
      else if ( indx == 3 ) then
        genz_name = 'Corner Peak  '
      else if ( indx == 4 ) then
        genz_name = 'Gaussian     '
      else if ( indx == 5 ) then
        genz_name = 'C0 Function  '
      else if ( indx == 6 ) then
        genz_name = 'Discontinuous'
      else
        genz_name = '?????????????'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  GENZ_NAME - Fatal error!'
        write ( *, '(a)' ) '  1 <= INDX <= 6 is required.'
        stop
      end if

      return
      end
      function genz_phi ( z )

c*********************************************************************72
c
cc GENZ_PHI estimates the normal cumulative density function.
c
c  Discussion:
c
c    The approximation is accurate to 1.0D-07.
c
c    This routine is based upon algorithm 5666 for the error function,
c    from Hart et al.
c
c  Modified:
c
c    13 March 2007
c
c  Author:
c
c    Alan Miller
c
c  Reference:
c
c    John Hart, Ward Cheney, Charles Lawson, Hans Maehly, 
c    Charles Mesztenyi, John Rice, Henry Thatcher, 
c    Christoph Witzgall,
c    Computer Approximations,
c    Wiley, 1968,
c    LC: QA297.C64.
c
c  Parameters:
c
c    Input, double precision Z, a value which can be regarded as the distance,
c    in standard deviations, from the mean.
c
c    Output, double precision GENZ_PHI, the integral of the normal PDF 
c    from negative infinity to Z.
c
      implicit none

      double precision expntl
      double precision genz_phi
      double precision p
      double precision p0
      parameter ( p0 = 220.2068679123761D+00 )
      double precision p1
      parameter ( p1 = 221.2135961699311D+00 )
      double precision p2
      parameter ( p2 = 112.0792914978709D+00 )
      double precision p3
      parameter ( p3 = 33.91286607838300D+00 )
      double precision p4
      parameter ( p4 = 6.373962203531650D+00 )
      double precision p5
      parameter ( p5 = 0.7003830644436881D+00 )
      double precision p6
      parameter ( p6 = 0.03526249659989109D+00 )
      double precision q0
      parameter ( q0 = 440.4137358247522D+00 )
      double precision q1
      parameter ( q1 = 793.8265125199484D+00 )
      double precision q2
      parameter ( q2 = 637.3336333788311D+00 )
      double precision q3
      parameter ( q3 = 296.5642487796737D+00 )
      double precision q4
      parameter ( q4 = 86.78073220294608D+00 )
      double precision q5
      parameter ( q5 = 16.06417757920695D+00 )
      double precision q6
      parameter ( q6 = 1.755667163182642D+00 )
      double precision q7
      parameter ( q7 = 0.08838834764831844D+00 )
c
c  This is square root of TWO pi.
c
      double precision rootpi
      parameter ( rootpi = 2.506628274631001D+00 )
      double precision z
      double precision zabs

      zabs = abs ( z )
c
c  |Z| > 12
c
      if ( 12.0D+00 .lt. zabs ) then

        p = 0.0D+00

      else
c
c  |Z| <= 12
c
        expntl = exp ( - zabs * zabs / 2.0D+00 )
c
c  |Z| < 7
c
        if ( zabs .lt. 7.0D+00 ) then

          p = expntl * ((((((
     &                p6
     &       * zabs + p5 )
     &       * zabs + p4 )
     &       * zabs + p3 )
     &       * zabs + p2 )
     &       * zabs + p1 )
     &       * zabs + p0 ) / (((((((
     &                q7
     &       * zabs + q6 )
     &       * zabs + q5 )
     &       * zabs + q4 )
     &       * zabs + q3 )
     &       * zabs + q2 )
     &       * zabs + q1 )
     &       * zabs + q0 )
c
c  CUTOFF <= |Z|
c
        else

          p = expntl / (
     &      zabs + 1.0D+00 / (
     &      zabs + 2.0D+00 / (
     &      zabs + 3.0D+00 / (
     &      zabs + 4.0D+00 / (
     &      zabs + 0.65D+00 ))))) / rootpi

        end if

      end if

      if ( 0.0D+00 .lt. z ) then
        p = 1.0D+00 - p
      end if

      genz_phi = p

      return
      end
      function genz_random ( )

c*********************************************************************72
c
cc GENZ_RANDOM is a portable random number generator
c
c  Modified:
c
c    12 March 2007
c
c  Author:
c
c    Linus Schrage
c
c  Reference:
c
c    Linus Schrage,
c    A More Portable Fortran Random Number Generator,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 2, June 1979, pages 132-138.
c
c  Parameters:
c
c    Output, double precision GENZ_RANDOM, a pseudorandom value.
c
      implicit none

      integer a
      parameter ( a = 16807 )
      integer b15
      parameter ( b15 = 32768 )
      integer b16
      parameter ( b16 = 65536 )
      integer fhi
      double precision genz_random
      integer k
      integer leftlo
      integer p
      parameter ( p = 2147483647 )
      integer seed
      integer xalo
      integer xhi

      save seed

      data seed / 123456 /

      xhi = seed / b16
      xalo = ( seed - xhi * b16 ) * a
      leftlo = xalo / b16
      fhi = xhi * a + leftlo
      k = fhi / b15

      seed = (
     &         (
     &           ( xalo - leftlo * b16 ) - p
     &         )
     &       + ( fhi - k * b15 ) * b16
     &       ) + k

      if ( seed .lt. 0 ) then
        seed = seed + p
      end if

      genz_random = dble ( seed ) / dble ( p )

      return
      end
      subroutine median ( n, r )

c*********************************************************************72
c
cc MEDIAN computes the median of an array of reals.
c
c  Discussion:
c
c    On exit R(1) contains the median, and (R(2),R(3)) specifies
c    the confidence interval.
c
c  Modified:
c
c    12 March 2007
c
c  Author:
c
c    Alan Genz
c
c  Parameters:
c
c    Input, integer N, the dimension of the array.
c
c    Input/output, double precision R(N).  On input, the array to be examined.
c    On output, the information in R has been overwritten, and only
c    R(1) through R(3) are meaningful.  R(1) contains the median,
c    R(2) and R(3) specify the confidence interval.
c
      implicit none

      integer n

      integer j
      integer k
      integer kmax
      integer nconf
      integer nd
      double precision r(n)
      double precision rmax
      double precision rmed

      if ( n .lt. 3 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MEDIAN - Fatal error!'
        write ( *, '(a)' ) '  This routine is not valid for N < 3.'
        stop
      end if

      do j = 1, n

        kmax = j

        do k = j, n
          if ( r(kmax) .lt. r(k) ) then
            kmax = k
          end if
        end do

        rmax = r(kmax)
        r(kmax) = r(j)
        r(j) = rmax

      end do

      nd = n / 2

      if ( mod ( n, 2 ) .eq. 0 ) then
        rmed = ( r(nd) + r(nd+1) ) / 2.0D+00
      else
        rmed = r(nd+1)
      end if

      nconf = max ( 1, ( 2 * n ) / 5 - 2 )
      rmax = r(n-nconf+1)

      r(3) = r(nconf)
      r(2) = rmax

      r(1) = rmed

      return
      end
      subroutine multst ( nsamp, tstlim, tstfns, tstmax, difclt,
     &  expnts, ndiml, ndims, sbname, subrtn, rel_tol, maxcls )

c*********************************************************************72
c
cc MULTST tests a multidimensional integration routine.
c
c  Discussion:
c
c    The routine uses the Genz test integrand functions, with
c    the user selecting the particular subset of test integrands, 
c    the set of difficulty factors, and the spatial dimensions.
c
c  Modified:
c
c    27 May 2007
c
c  Author:
c
c    Alan Genz
c
c  Reference:
c
c    Alan Genz,
c    A Package for Testing Multiple Integration Subroutines,
c    in Numerical Integration: 
c    Recent Developments, Software and Applications,
c    edited by Patrick Keast, Graeme Fairweather,
c    D Reidel, 1987, pages 337-340,
c    LC: QA299.3.N38.
c
c  Parameters:
c
c    Input, integer NSAMP, the number of samples.
c    1 <= NSAMP <= 50.
c
c    Input, integer TSTLIM, the number of test integrands.
c
c    Input, integer TSTFNS(TSTLIM), the indices of the test integrands.
c    Each index is between 1 and 6.
c
c    Input, integer TSTMAX, the number of difficulty levels to be tried.
c
c    Input, double precision DIFCLT(TSTMAX) of difficulty levels.
c
c    Input, double precision EXPNTS(TSTMAX), the difficulty exponents.
c
c    Input, integer NDIML, the number of sets of variable sizes.
c
c    Input, integer NDIMS(NDIML), the number of variables for the integrals
c    in each test.  The entries of NDIMS must not exceed 20.
c
c    Input, character*6 SBNAME, the name of the integration subroutine 
c    to be tested.
c
c    Input, exernal SUBRTN, the integration subroutine to be tested.
c
c    Input, double precision REL_TOL, the relative error tolerance.
c
c    Input, integer MAXCLS, the maximum number of integrand calls 
c    for all tests.
c
      implicit none

      integer lenwrk
      parameter ( lenwrk = 50000 )
      integer maxdim
      parameter ( maxdim = 20 )
      integer maxsmp
      parameter ( maxsmp = 50 )
      integer mxtsfn
      parameter ( mxtsfn = 6 )
      integer tstlim
      integer tstmax

      double precision a(maxdim)
      double precision alpha(maxdim)
      double precision b(maxdim)
      double precision beta(maxdim)
      double precision callsa(mxtsfn,mxtsfn)
      double precision callsb(mxtsfn,mxtsfn)
      double precision concof
      double precision dfclt
      double precision difclt(tstmax)
      integer digits
      double precision rel_tol
      double precision errest
      double precision errlog
      double precision ersacb(mxtsfn,mxtsfn)
      double precision ersact(mxtsfn,mxtsfn)
      double precision ersdsb(mxtsfn,mxtsfn)
      double precision ersdsc(mxtsfn,mxtsfn)
      double precision ersesb(mxtsfn,mxtsfn)
      double precision ersest(mxtsfn,mxtsfn)
      double precision ersrel(mxtsfn,mxtsfn)
      double precision estlog
      double precision exn
      double precision expnts(tstmax)
      double precision expons(mxtsfn)
      double precision finest
      external genz_function
      double precision genz_integral
      character*13 genz_name
      double precision genz_random
      integer i
      integer idfclt(mxtsfn)
      integer ifail
      integer ifails
      integer it
      integer itest
      integer itst
      integer j
      integer k
      integer maxcls
      double precision medacb(mxtsfn)
      double precision medact(maxsmp)
      double precision medcla(mxtsfn)
      double precision medclb(mxtsfn)
      double precision medcls(maxsmp)
      double precision meddsb(mxtsfn)
      double precision meddsc(maxsmp)
      double precision medesb(mxtsfn)
      double precision medest(maxsmp)
      double precision medrel
      double precision medrll(maxsmp)
      integer mincls
      integer n
      character*13 name
      integer nconf
      integer ndim
      integer ndiml
      integer ndims(ndiml)
      integer ndimv
      integer nsamp
      double precision qality
      double precision qallty(maxsmp)
      double precision qualty(mxtsfn,mxtsfn)
      integer rcalsa
      integer rcalsb
      double precision relerr
      character*6 sbname
      double precision small
      external subrtn
      double precision tactrb(mxtsfn)
      double precision tactrs(mxtsfn)
      double precision tcalsa(mxtsfn)
      double precision tcalsb(mxtsfn)
      double precision terdsb(mxtsfn)
      double precision terdsc(mxtsfn)
      double precision testrb(mxtsfn)
      double precision testrs(mxtsfn)
      double precision tqualt(mxtsfn)
      double precision trelib(mxtsfn)
      integer tstfns(tstlim)
      double precision value
      double precision wrkstr(lenwrk)

      common /tstblk/ itst, alpha, beta, dfclt, exn
c
c  Initialize and compute confidence coefficient.
c
      concof = 0.0D+00
      nconf = max ( 1, ( 2 * nsamp ) / 5 - 2 )

      do i = 1, nconf
        concof = 1.0D+00 + dble ( nsamp - nconf + i ) * concof 
     &    / dble ( nconf - i + 1 )
      end do

      concof = 1.0D+00 - concof / dble ( 2**(nsamp-1) )

      small = 1.0D+00

10    continue

      small = small / 2.0D+00
      if ( 1.0D+00 < 1.0D+00 + small ) then
        go to 10
      end if

      small = 2.0D+00 * small

      do it = 1, tstlim
        itest = tstfns(it)
        idfclt(it) = difclt(itest)
        expons(it) = expnts(itest)
      end do
c
c  Begin main loop for different numbers of variables.
c
      do ndimv = 1, ndiml

        ndim = ndims(ndimv)

        if ( mod ( ndimv - 1, 6 ) .eq. 0 ) then

          write ( *, '(a)' ) ' '
          write ( *, '(a,i4,a)' ) '         Test results with', 
     &      nsamp, ' samples per test'
          write ( *, '(a)' ) ' '

          write ( *, '(a,10i6)' ) 
     &      '  Difficulty levels', (idfclt(j),j=1,tstlim)

          write ( *, '(a,10f6.1)' ) 
     &      '      Exponents    ', (expons(j),j=1,tstlim)

          digits = int ( - log10 ( rel_tol ) )

          write ( *, '(a)' ) ' '
          write ( *, '(a,i3,a,i8)') '   Requested digits =', digits, 
     &      ', Maximum values =', maxcls
          write ( *, '(a)' ) ' '
          write ( *, '(2x,a,a,f5.2)' ) sbname, 
     &      ' tests, variable results with confidence', concof
          write ( *, '(a,a)' )
     &      ' Vari-  integrand     Correct digits   Relia-  Wrong',
     &      '   Integrand   Quality Total'
          write ( *, '(a,a)' )
     &      ' ables              Estimated   Actual bility Digits',
     &      '    Values             Fails'
          write ( *, '(a)' ) ' '

        end if
c
c  Begin loop for different test integrands.
c
        do it = 1, tstlim

          itest = tstfns(it)
          itst = itest
          exn = expnts(itest)
          dfclt = difclt(itest)

          do j = 1, ndim
            a(j) = 0.0D+00
            b(j) = 1.0D+00
          end do

          ifails = 0
          medrel = 0
c
c  Begin loop for different samples.
c
          do k = 1, nsamp

            ifail = 1

            do n = 1, ndim
              alpha(n) = genz_random ( )
              beta(n) = genz_random ( )
            end do

            value = genz_integral ( itest, ndim, a, b )
c
c  Call integration subroutine.
c
            mincls = 4 * 2**ndim

            call subrtn ( ndim, a, b, mincls, maxcls, genz_function, 
     &        rel_tol, errest, lenwrk, wrkstr, finest, ifail )

            relerr = abs ( ( finest - value ) / value )
            ifails = ifails + min ( ifail, 1 )
            relerr = max ( min ( 1.0D+00, relerr ), small )
            errlog = max ( 0.0D+00, -log10 ( relerr ) )
            errest = max ( min ( 1.0D+00, errest ), small )
            estlog = max ( 0.0D+00, -log10 ( errest ) )
            meddsc(k) = max ( 0.0D+00, estlog - errlog )
            medest(k) = estlog
            medact(k) = errlog
            medcls(k) = mincls

            if ( errest .ge. relerr ) then
              medrel = medrel + 1
            end if

          end do
c
c  End loop for different samples and compute medians.
c
          call median ( nsamp, medest )
          call median ( nsamp, medact )
          call median ( nsamp, medcls )
          call median ( nsamp, meddsc )

          medrel = medrel / dble ( nsamp )

          trelib(it) = medrel
          tactrs(it) = medact(2)
          testrs(it) = medest(2)
          terdsc(it) = meddsc(2)
          tcalsa(it) = medcls(2)
          tcalsb(it) = medcls(3)
          tactrb(it) = medact(3)
          testrb(it) = medest(3)
          terdsb(it) = meddsc(3)

          ersrel(itest,ndimv) = medrel
          ersest(itest,ndimv) = medest(2)
          ersact(itest,ndimv) = medact(2)
          ersdsc(itest,ndimv) = meddsc(2)
          ersesb(itest,ndimv) = medest(3)
          ersacb(itest,ndimv) = medact(3)
          ersdsb(itest,ndimv) = meddsc(3)
          callsa(itest,ndimv) = medcls(2)
          callsb(itest,ndimv) = medcls(3)

          qality = 0.0D+00

          if ( medcls(1) .ne. 0.0D+00 ) then
            qality = ( medact(1) + 1.0D+00 )*
     &        ( medest(1) + 1.0D+00 - meddsc(1) ) / log ( medcls(1) )
          end if

          tqualt(it) = qality
          qualty(itest,ndimv) = qality
          rcalsa = medcls(2)
          rcalsb = medcls(3)
          name = genz_name ( itest )

          write ( *,99995) ndim, name, medest(2),
     &           medest(3), medact(2), medact(3), medrel, meddsc(2),
     &           meddsc(3), rcalsa, rcalsb, qality, ifails
99995     format ( 2x, i2, 2x, a14, 2(f4.1, ',', f4.1, 1x), f4.2,
     &           f4.1, ',', f3.1, i7, ',', i7, f6.2, i5)

        end do
c
c  End loop for different test integrands.
c
        call median ( tstlim, tactrs )
        call median ( tstlim, trelib )
        call median ( tstlim, testrs )
        call median ( tstlim, terdsc )
        call median ( tstlim, tactrb )
        call median ( tstlim, testrb )
        call median ( tstlim, terdsb )
        call median ( tstlim, tqualt )
        call median ( tstlim, tcalsa )
        call median ( tstlim, tcalsb )

        rcalsa = tcalsa(1)
        rcalsb = tcalsb(1)

        write (*,99994) ndim, testrs(1), testrb(1), tactrs(1),
     &        tactrb( 1), trelib(1), terdsc(1), terdsb(1),
     &        rcalsa, rcalsb, tqualt(1)
        write ( *, '(a)' ) ' '

99994   format (2x, i2, '   Medians', 6x, 2(f4.1, ',', f4.1, 1x),
     &        f4.2, f4.1, ',', f3.1, i7, ',', i7, f6.2 )

      end do
c
c  End loop for different numbers of variables.
c
      if ( ndiml .gt. 1 ) then

        write ( *, '(a)' ) ' '
        write ( *, '(6x,a,a,12i3)' ) sbname, 
     &    ' Test integrand medians for variables',
     &    (ndims(ndimv),ndimv=1,ndiml)

        write ( *, '(a,a)' )
     &    '        Integrand     Correct digits   Relia-  Wrong',
     &    '   Integrand   Quality'
        write ( *, '(a,a)' )
     &    '                    Estimated   Actual bility digits',
     &    '     Values'
        write ( *, '(a)' ) ' '

        do it = 1, tstlim

          itest = tstfns(it)

          do ndimv = 1, ndiml
            medact(ndimv) = ersact(itest,ndimv)
            medest(ndimv) = ersest(itest,ndimv)
            meddsc(ndimv) = ersdsc(itest,ndimv)
            medacb(ndimv) = ersacb(itest,ndimv)
            medesb(ndimv) = ersesb(itest,ndimv)
            meddsb(ndimv) = ersdsb(itest,ndimv)
            medrll(ndimv) = ersrel(itest,ndimv)
            qallty(ndimv) = qualty(itest,ndimv)
            medcla(ndimv) = callsa(itest,ndimv)
            medclb(ndimv) = callsb(itest,ndimv)
          end do

          call median ( ndiml, medrll )
          call median ( ndiml, medact )
          call median ( ndiml, medest )
          call median ( ndiml, meddsc )
          call median ( ndiml, medacb )
          call median ( ndiml, medesb )
          call median ( ndiml, meddsb )
          call median ( ndiml, qallty )
          call median ( ndiml, medcla )
          call median ( ndiml, medclb )

          rcalsa = medcla(1)
          rcalsb = medclb(1)
          name = genz_name ( itest )

          write (*,99991) name, medest(1), medesb(1),
     &           medact(1), medacb(1),medrll(1), meddsc(1), meddsb(1),
     &           rcalsa, rcalsb, qallty(1)
99991     format (6x, a14, 2(f4.1, ',', f4.1, 1x), f4.2, f4.1, ',',
     &           f3.1, i7, ',', i7, f6.2)

          tactrs(it) = medact(1)
          testrs(it) = medest(1)
          terdsc(it) = meddsc(1)
          tactrb(it) = medacb(1)
          testrb(it) = medesb(1)
          terdsb(it) = meddsb(1)
          tcalsa(it) = medcla(1)
          tcalsb(it) = medclb(1)
          trelib(it) = medrll(1)
          tqualt(it) = qallty(1)

        end do

        call median ( tstlim, tactrs)
        call median ( tstlim, testrs)
        call median ( tstlim, terdsc)
        call median ( tstlim, tactrb)
        call median ( tstlim, testrb)
        call median ( tstlim, terdsb)
        call median ( tstlim, trelib)
        call median ( tstlim, tqualt)
        call median ( tstlim, tcalsa)
        call median ( tstlim, tcalsb)

        rcalsa = tcalsa(1)
        rcalsb = tcalsb(1)

        write (*,99990) testrs(1), testrb(1), tactrs(1), tactrb(1),
     &        trelib(1), terdsc(1), terdsb(1), rcalsa, rcalsb, tqualt(1)
        write ( *, '(a)' ) ' '
99990   format ( '     Global medians ', 2(f4.1, ',', f4.1, 1x), f4.2,
     &        f4.1, ',', f3.1, i7, ',', i7, f6.2)

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
