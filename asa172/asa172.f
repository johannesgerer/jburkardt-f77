      subroutine revers ( ivec, kdim )

c*********************************************************************72
c
cc REVERS reorders the subscript vector, if required.
c
c  Modified:
c
c    27 July 2008
c
c  Author:
c
c    M O'Flaherty, G MacKenzie
c
c  Reference:
c
c    M O'Flaherty, G MacKenzie,
c    Algorithm AS 172:
c    Direct Simulation of Nested Fortran DO-LOOPS,
c    Applied Statistics,
c    Volume 31, Number 1, 1982, pages 71-74.
c
c  Parameters:
c
c    Input/output, integer IVEC(KDIM), the subscript vector.
c
c    Input, integer KDIM, the dimension of the subscript vector.
c
      implicit none

      integer kdim

      integer i
      integer itemp
      integer ivec(kdim)

      do i = 1, kdim / 2
        itemp          = ivec(i)
        ivec(i)        = ivec(kdim+1-i)
        ivec(kdim+1-i) = itemp
      end do

      return
      end
      subroutine simdo ( qind, qfor, iprod, kdim, jsub, ivec, ifault )

c*********************************************************************72
c
cc SIMDO generates multi-indices, simulating nested DO-loops.
c
c  Discussion:
c
c    The loops are assumed to be nested to a depth of K.
c
c    The R-th loop is assumed to have upper limit N(R) and increment Inc(R).
c
c    The total number of executions of the innermost loop is 
c
c      N = product ( 1 <= R <= K ) N(R).
c
c    Let these executions be indexed by the single integer J, which
c    we call the index subscript.
c
c    Each value of J corresponds to a particular set of loop indices,
c    which we call the subscript vector I(J).
c
c    This routine can start with J and find I(J), or determine
c    J from I(J).
c    
c  Modified:
c
c    26 July 2008
c
c  Author:
c
c    M O'Flaherty, G MacKenzie
c
c  Reference:
c
c    M O'Flaherty, G MacKenzie,
c    Algorithm AS 172:
c    Direct Simulation of Nested Fortran DO-LOOPS,
c    Applied Statistics,
c    Volume 31, Number 1, 1982, pages 71-74.
c
c  Parameters:
c
c    Input, logical QIND.
c    TRUE to convert an index subscript J to the subscript vector I(J).
c    FALSE to convert the subscript vector I(J) to the index subscript J.
c
c    Input, logical QFOR,
c    TRUE if conversion is required in standard Fortran subscripting order,
c    FALSE otherwise.
c
c    Input, integer IPROD(KDIM), contains the partial products.
c    If QFOR is FALSE, then
c      IPROD(S) = product ( 1 <= R <= S ) N(R).
c    If QFOR is TRUE, then
c      IPROD(S) = product ( 1 <= R <= S ) N(KDIM+1-R).
c
c    Input, integer KDIM, the nesting depth of the loops.
c
c    Input/output, integer JSUB.
c    If QIND is TRUE, then JSUB is an input quantity, an index subscript J
c    to be converted into the subscript vector I(J).
c    If QIND is FALSE, then JSUB is an output quantity, the index subscript J
c    corresponding to the subscript vector I(J).
c
c    Input/output, integer IVEC(KDIM).
c    if QIND is TRUE, then IVEC is an output quantity, the subscript vector I(J)
c    corresponding to the index subscript J.
c    If QIND is FALSE, then IVEC is an input quantity, a subscript vector I(J)
c    for which the corresponding index subscript J is to be computed.
c
c    Output, integer IFAULT, error flag.
c    0, no error was detected.
c    1, if QIND is TRUE, and the input value of JSUB exceeds IPROD(KDIM).
c    2, if QIND is FALSE, and IVEC contains an illegal component.
c
      implicit none

      integer kdim

      integer i
      integer ifault
      integer ij
      integer ik
      integer iprod(kdim)
      integer itempv
      integer ivec(kdim)
      integer jsub
      logical qfor
      logical qind

      ifault = 0
c
c  Index subscript to subscript vector conversion.
c
      if ( qind ) then

        if ( iprod(kdim) .lt. jsub ) then
          ifault = 1
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SIMDO - Fatal error!'
          write ( *, '(a)' ) '  JSUB is out of bounds.'
          stop
        end if

        itempv = jsub - 1
        ij = kdim - 1

        do i = 1, ij
          ik = kdim - i
          ivec(i) = itempv / iprod(ik)
          itempv = itempv - iprod(ik) * ivec(i)
          ivec(i) = ivec(i) + 1
        end do

        ivec(kdim) = itempv + 1

        if ( qfor ) then
          call revers ( ivec, kdim )
        end if
c
c  Subscript vector to index subscript conversion.
c
      else

        if ( .not. qfor ) then
          call revers ( ivec, kdim )
        end if

        if ( iprod(1) .lt. ivec(1) ) then
          ifault = 2
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'SIMDO - Fatal error!'
          write ( *, '(a)' ) '  An entry of IVEC is illegal.'
          stop
        end if

        do i = 2, kdim
          if ( iprod(i) / iprod(i-1) .lt. ivec(i) ) then
            ifault = 2
            write ( *, '(a)' ) ' '
            write ( *, '(a)' ) 'SIMDO - Fatal error!'
            write ( *, '(a)' ) '  An entry of IVEC is illegal.'
            stop
          end if
        end do

        jsub = ivec(1)

        do i = 2, kdim
          jsub = jsub + ( ivec(i) - 1 ) * iprod(i-1)
        end do
c
c  As a courtesy to the caller, UNREVERSE the IVEC vector
c  if you reversed it.
c
        if ( .not. qfor ) then
          call revers ( ivec, kdim )
        end if

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
