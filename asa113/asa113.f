      subroutine swap ( varval, class, clsize, in, ik, iv, critvl,
     &  ntrans, ifault )

c*********************************************************************72
c
cc SWAP interchanges objects between different classes to improve a criterion.
c
c  Discussion:
c
c    This routine is given a classification of objects, including the
c    number of objects in each class, and the current value of some criterion
c    which is desired to be minimized.
c
c    The routine calculates the change in criterion for all possible swaps,
c    that is, operations in which two objects in different classes exchange 
c    places. Each swap that would result in a lowering of the criterion is 
c    executed, and the related quantities are updated.
c
c    When no more advantageous swaps can be found, the routine returns.
c
c    The routine relies on a user-supplied routine, CRSWAP, to report the
c    expected change in the criterion for a given swap, and to carry
c    out that transfer if requested.
c
c    The variables CLASS and CRITVL have been added to the argument list
c    of CRSWAP.
c
c    Also, the order of the two classes "L" and "M" was interchanged in
c    the call to CRSWAP.  The original order was counterintuitive.
c
c  Modified:
c
c    16 February 2008
c
c  Author:
c
c    FORTRAN77 original by Banfield, Bassill
c    Modifications by John Burkardt
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
c    Input, double precision VARVAL(IN,IV), the data values.  There are IN objects,
c    each having spatial dimension IV.
c
c    Input/output, integer CLASS(IN), the classification of each object.
c
c    Input/output, integer CLSIZE(IK), the number of objects in each class.
c
c    Input, integer IN, the number of objects.
c
c    Input, integer IK, the number of classes.
c
c    Input, integer IV, the number of spatial dimensions, or variates, 
c    of the objects.
c
c    Input/output, double precision CRITVL, the current value of the criterion.
c
c    Output, integer NTRANS, the number of transfers executed.
c
c    Output, integer IFAULT, error indicator.
c    0, no error detected.
c    1, the number of classes was less than 2.
c    2, the number of objects was less than the number of classes.
c
      implicit none

      integer ik
      integer in
      integer iv

      integer class(in)
      integer clsize(ik)
      double precision critvl
      double precision eps
      parameter ( eps = 1.0D-38 )
      integer i
      integer icount
      integer ifault
      double precision inc
      integer iswitch
      integer it
      integer itop
      integer j
      integer k
      integer l
      integer m
      integer ntrans
      double precision varval(in,iv)

      if ( ik .le. 1 ) then
        ifault = 1
        return
      end if

      if ( in .le. ik ) then
        ifault = 2
        return
      end if

      ifault = 0
      icount = 0
      ntrans = 0
      itop = ( in * ( in - 1 ) ) / 2
      i = 1

2     continue

      i = i + 1

      if ( itop .le. icount ) then
        return
      end if

      if ( in .lt. i ) then
        i = 1
        go to 2
      end if

      l = class(i)
      k = l
      it = i - 1
c
c  Test the swap of object I from class M to L, 
c  and object J from class L to M.
c
      do j = 1, it

        icount = icount + 1
        m = class(j)

        if ( l .ne. j ) then

          if ( clsize(l) .ne. 1 .or. clsize(m) .ne. 1 ) then

            iswitch = 1
            call crswap ( varval, class, clsize, in, ik, iv, critvl,
     &        i, j, l, m, iswitch, inc )
            
            if ( inc .lt. - eps ) then

              critvl = critvl + inc
              icount = 0

              iswitch = 2
              call crswap ( varval, class, clsize, in, ik, iv, critvl,
     &          i, j, l, m, iswitch, inc )

              ntrans = ntrans + 1
              class(i) = m
              class(j) = l
              l = m

            end if

          end if
        end if

      end do

      go to 2
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
      subroutine trnsfr ( varval, class, clsize, in, ik, iv, critvl,
     &  ntrans, ifault )

c*********************************************************************72
c
cc TRNSFR transfers objects between classes to improve a criterion.
c
c  Discussion:
c
c    This routine is given a classification of objects, including the
c    number of objects in each class, and the current value of some criterion
c    which is desired to be minimized.
c
c    The routine calculates the change in criterion for all possible transfers
c    of any object from its current class to a different class.  Each transfer
c    that would result in a lowering of the criterion is executed, and the
c    related quantities are updated.
c
c    When no more advantageous transfers can be found, the routine returns.
c
c    The routine relies on a user-supplied routine, CRTRAN, to report the
c    expected change in the criterion for a given transfer, and to carry
c    out that transfer if requested.
c
c    The variables CLASS and CRITVL have been added to the argument list
c    of CRTRAN.
c
c    Also, the order of the two classes "L" and "M" was interchanged in
c    the call to CRTRAN.  The original order was counterintuitive.
c
c  Modified:
c
c    15 February 2008
c
c  Author:
c
c    FORTRAN77 original by Banfield, Bassill
c    Modifications by John Burkardt
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
c    Input, double precision VARVAL(IN,IV), the data values.  There are IN objects,
c    each having spatial dimension IV.
c
c    Input/output, integer CLASS(IN), the classification of each object.
c
c    Input/output, integer CLSIZE(IK), the number of objects in each class.
c
c    Input, integer IN, the number of objects.
c
c    Input, integer IK, the number of classes.
c
c    Input, integer IV, the number of spatial dimensions, or variates, 
c    of the objects.
c
c    Input/output, double precision CRITVL, the current value of the criterion.
c
c    Output, integer NTRANS, the number of transfers executed.
c
c    Output, integer IFAULT, error indicator.
c    0, no error detected.
c    1, the number of classes was less than 2.
c    2, the number of objects was less than the number of classes.
c
      implicit none

      integer ik
      integer in
      integer iv

      integer class(in)
      integer clsize(ik)
      double precision critvl
      double precision eps
      parameter ( eps = 1.0D-38 )
      integer i
      integer icount
      integer ifault
      double precision inc
      double precision inco
      integer iswitch
      integer l
      integer lo
      integer m
      integer ntrans
      double precision varval(in,iv)

      if ( ik .le. 1 ) then
        ifault = 1
        return
      end if

      if ( in .le. ik ) then
        ifault = 2
        return
      end if

      ifault = 0
      ntrans = 0
      i = 0
      icount = 0

2     continue

      i = i + 1

      if ( in .le. icount ) then
        return
      end if

      if ( in .lt. i ) then
        i = 0
        icount = 0
        go to 2
      end if

      m = class(i)
      if ( clsize(m) .le. 1 ) then
        icount = icount + 1
        go to 2
      end if

      inco = - eps
      lo = m
c
c  Test the transfer of object I from class M to class L.
c
      do l = 1, ik

        if ( l .ne. m ) then

          iswitch = 1
          call crtran ( varval, class, clsize, in, ik, iv, critvl, 
     &      i, m, l, iswitch, inc )
c
c  Remember the values of L and INC.
c
          if ( inc .lt. inco ) then
            lo = l
            inco = inc
          end if

        end if

      end do

      icount = icount + 1
c
c  Execute the transfer of object I from class M to class LO.
c
      if ( lo .ne. m ) then

        l = lo
        critvl = critvl + inco
        icount = 0

        iswitch = 2
        call crtran ( varval, class, clsize, in, ik, iv, critvl,
     &    i, m, l, iswitch, inc )

        ntrans = ntrans + 1
        class(i) = l
        clsize(l) = clsize(l) + 1
        clsize(m) = clsize(m) - 1

      end if

      go to 2
      end
