      subroutine nelmin ( fn, n, start, xmin, ynewlo, reqmin, step,
     &  konvge, kcount, icount, numres, ifault )

c*********************************************************************72
c
cc NELMIN minimizes a function using the Nelder-Mead algorithm.
c
c  Discussion:
c
c    This routine seeks the minimum value of a user-specified function.
c
c     Simplex function minimisation procedure due to Nelder+Mead(1965),
c     as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
c     subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
c     25, 97) and Hill(1978, 27, 380-2)
c
c    The function to be minimized must be defined by a function of
c    the form
c
c      function fn ( x, f )
c      double precision fn
c      double precision x(*)
c
c    and the name of this subroutine must be declared EXTERNAL in the
c    calling routine and passed as the argument FN.
c
c    This routine does not include a termination test using the
c    fitting of a quadratic surface.
c
c  Modified:
c
c    27 February 2008
c
c  Author:
c
c    FORTRAN77 version by R ONeill
c    Modifications by John Burkardt
c
c  Reference:
c
c    John Nelder, Roger Mead,
c    A simplex method for function minimization,
c    Computer Journal,
c    Volume 7, 1965, pages 308-313.
c
c    R ONeill,
c    Algorithm AS 47:
c    Function Minimization Using a Simplex Procedure,
c    Applied Statistics,
c    Volume 20, Number 3, 1971, pages 338-345.
c
c  Parameters:
c
c    Input, external FN, the name of the function which evaluates
c    the function to be minimized.
c
c    Input, integer N, the number of variables.
c
c    Input/output, double precision START(N).  On input, a starting point
c    for the iteration.  On output, this data may have been overwritten.
c
c    Output, double precision XMIN(N), the coordinates of the point which
c    is estimated to minimize the function.
c
c    Output, double precision YNEWLO, the minimum value of the function.
c
c    Input, double precision REQMIN, the terminating limit for the variance
c    of function values.
c
c    Input, double precision STEP(N), determines the size and shape of the
c    initial simplex.  The relative magnitudes of its elements should reflect
c    the units of the variables.
c
c    Input, integer KONVGE, the convergence check is carried out every
c    KONVGE iterations.
c
c    Input, integer KCOUNT, the maximum number of function evaluations.
c
c    Output, integer ICOUNT, the number of function evaluations used.
c
c    Output, integer NUMRES, the number of restarts.
c
c    Output, integer IFAULT, error indicator.
c    0, no errors detected.
c    1, REQMIN, N, or KONVGE has an illegal value.
c    2, iteration terminated because KCOUNT was exceeded without convergence.
c
      implicit none

      integer n
      integer n_max
      parameter ( n_max = 20 )

      double precision ccoeff
      parameter ( ccoeff = 0.5D+00 )
      double precision del
      double precision dn
      double precision dnn
      double precision ecoeff
      parameter ( ecoeff = 2.0D+00 )
      double precision eps
      parameter ( eps = 0.001D+00 )
      double precision fn
      external fn
      integer i
      integer icount
      integer ifault
      integer ihi
      integer ilo
      integer j
      integer jcount
      integer kcount
      integer konvge
      integer l
      integer nn
      integer numres
      double precision p(n_max,n_max+1)
      double precision pstar(n_max)
      double precision p2star(n_max)
      double precision pbar(n_max)
      double precision rcoeff
      parameter ( rcoeff = 1.0D+00 )
      double precision reqmin
      double precision rq
      double precision start(n)
      double precision step(n)
      double precision x
      double precision xmin(n)
      double precision y(n_max+1)
      double precision y2star
      double precision ylo
      double precision ynewlo
      double precision ystar
      double precision z
c
c  Check the input parameters.
c
      if ( reqmin .le. 0.0D+00 ) then
        ifault = 1
        return
      end if

      if ( n .lt. 1 ) then
        ifault = 1
        return
      end if

      if ( n_max .lt. n ) then
        ifault = 1
        return
      end if

      if ( konvge .lt. 1 ) then
        ifault = 1
        return
      end if

      icount = 0
      numres = 0

      jcount = konvge  
      dn = dble ( n )   
      nn = n + 1         
      dnn = dble ( nn ) 
      del = 1.0D+00
      rq = reqmin * dn
c
c  Construction of initial simplex.
c
   10 continue

      do i = 1, n     
        p(i,nn) = start(i)
      end do

      y(nn) = fn(start)

      do j = 1, n     
        x = start(j)
        start(j) = start(j) + step(j) * del
        do i = 1, n
          p(i,j) = start(i)
        end do
        y(j) = fn ( start )
        start(j) = x
      end do

      icount = icount + nn
c                    
c  The simplex construction is complete.
c                    
c  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
c  the vertex of the simplex to be replaced.
c                    
      ylo = y(1)
      ilo = 1

      do i = 2, nn
        if ( y(i) .lt. ylo ) then
          ylo = y(i) 
          ilo = i
        end if
      end do

   50 continue

      ynewlo = y(1)
      ihi = 1

      do i = 2, nn
        if ( ynewlo .lt. y(i) ) then
          ynewlo = y(i)
          ihi = i
        end if
      end do
c
c  Calculate PBAR, the centroid of the simplex vertices
c  excepting the vertex with Y value YNEWLO.
c
      do i = 1, n
        z = 0.0D+00
        do j = 1, nn    
          z = z + p(i,j)
        end do
        z = z - p(i,ihi)   
        pbar(i) = z / dn   
      end do
c
c  Reflection through the centroid.
c
      do i = 1, n
        pstar(i) = pbar(i) + rcoeff * ( pbar(i) - p(i,ihi) )
      end do

      ystar = fn ( pstar )
      icount = icount + 1
c
c  Successful reflection, so extension.
c
      if ( ystar .lt. ylo ) then

        do i = 1, n
          p2star(i) = pbar(i) + ecoeff * ( pstar(i) - pbar(i) )
        end do

        y2star = fn ( p2star )
        icount = icount + 1
c
c  Check extension.
c
        if ( ystar .lt. y2star ) then

          do i = 1, n
            p(i,ihi) = pstar(i)
          end do

          y(ihi) = ystar
c
c  Retain extension or contraction.
c
        else

          do i = 1, n
            p(i,ihi) = p2star(i)
          end do

          y(ihi) = y2star

        end if
c
c  No extension.
c
      else

        l = 0
        do i = 1, nn
          if ( ystar .lt. y(i) ) then
            l = l + 1
          end if
        end do

        if ( 1 .lt. l ) then

          do i = 1, n
            p(i,ihi) = pstar(i)
          end do

          y(ihi) = ystar
c
c  Contraction on the  Y(IHI) side of the centroid.
c
        else if ( l .eq. 0 ) then

          do i = 1, n
            p2star(i) = pbar(i) + ccoeff * ( p(i,ihi) - pbar(i) )
          end do
          y2star = fn ( p2star )
          icount = icount + 1
c
c  Contract the whole simplex.
c
          if ( y(ihi) .lt. y2star ) then

            do j = 1, nn
              do i = 1, n
                p(i,j) = ( p(i,j) + p(i,ilo) ) * 0.5D+00
                xmin(i) = p(i,j)
              end do
              y(j) = fn ( xmin )
            end do

            icount = icount + nn
            if ( kcount .lt. icount ) then
               go to 260
            end if

            ylo = y(1)
            ilo = 1

            do i = 2, nn
              if ( y(i) .lt. ylo ) then
                ylo = y(i) 
                ilo = i
              end if
            end do

            go to 50
c
c  Retain contraction.
c
          else

            do i = 1, n
              p(i,ihi) = p2star(i)
            end do
            y(ihi) = y2star

          end if
c
c  Contraction on the reflection side of the centroid.
c
        else if ( l .eq. 1 ) then

          do i = 1, n
            p2star(i) = pbar(i) + ccoeff * ( pstar(i) - pbar(i) )
          end do

          y2star = fn ( p2star )
          icount = icount + 1
c
c  Retain reflection?
c
          if ( y2star .le. ystar ) then

            do i = 1, n
              p(i,ihi) = p2star(i)
            end do
            y(ihi) = y2star

          else

            do i = 1, n
              p(i,ihi) = pstar(i)
            end do
            y(ihi) = ystar  

          end if
 
        end if

      end if
c
c  Check if YLO improved.
c
      if ( y(ihi) .lt. ylo ) then
        ylo = y(ihi)
        ilo = ihi
      end if

      jcount = jcount - 1

      if ( jcount .ne. 0 ) then
        go to 50
      end if
c
c  Check to see if minimum reached.
c
      if ( icount .le. kcount ) then

        jcount = konvge

        z = 0.0D+00
        do i = 1, nn
          z = z + y(i)
        end do
        x = z / dnn

        z = 0.0D+00
        do i = 1, nn
          z = z + ( y(i) - x )**2
        end do

        if ( rq .lt. z ) then
          go to 50
        end if

      end if
c
c  Factorial tests to check that YNEWLO is a local minimum.
c
  260 continue

      do i = 1, n
        xmin(i) = p(i,ilo)
      end do

      ynewlo = y(ilo)

      if ( kcount .lt. icount ) then
        ifault = 2
        return
      end if

      ifault = 0

      do i = 1, n
        del = step(i) * eps
        xmin(i) = xmin(i) + del
        z = fn ( xmin )
        icount = icount + 1
        if ( z .lt. ynewlo ) then
          ifault = 2
          go to 290
        end if
        xmin(i) = xmin(i) - del - del
        z = fn ( xmin )
        icount = icount + 1
        if ( z .lt. ynewlo ) then
          ifault = 2
          go to 290
        end if
        xmin(i) = xmin(i) + del
      end do

290   continue

      if ( ifault == 0 ) then
        return
      end if
c
c  Restart the procedure.
c
      do i = 1, n
        start(i) = xmin(i)
      end do

      del = eps
      numres = numres + 1
      go to 10

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

