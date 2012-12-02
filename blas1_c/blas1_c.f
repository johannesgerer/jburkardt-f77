      function cabs1 ( z )

c*********************************************************************72
c
cc CABS1 returns the L1 norm of a single precision complex number.
c
c  Discussion:
c
c    The L1 norm of a complex number is the sum of the absolute values
c    of the real and imaginary components.
c
c    CABS1 ( Z ) = abs ( real ( Z ) ) + abs ( imaginary ( Z ) )
c
c  Modified:
c
c    22 May 2002
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, complex Z, the number whose norm is desired.
c
c    Output, real CABS1, the L1 norm of Z.
c
      implicit none

      real cabs1
      complex z

      cabs1 = abs ( real ( z ) ) + abs ( aimag ( z ) )

      return
      end
      function cabs2 ( z )

c*********************************************************************72
c
cc CABS2 returns the L2 norm of a single precision complex number.
c
c  Discussion:
c
c    The L2 norm of a complex number is the square root of the sum of the 
c    squares of the real and imaginary components.
c
c    CABS2 ( Z ) = sqrt ( ( real ( Z ) )**2 + ( imaginary ( Z ) )**2 )
c
c  Modified:
c
c    14 April 2006
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, complex Z, the number whose norm is desired.
c
c    Output, real CABS2, the L2 norm of Z.
c
      implicit none

      real cabs2
      complex z

      cabs2 = sqrt ( ( real ( z ) )**2 + ( aimag ( z ) )**2 )

      return
      end
      subroutine caxpy ( n, ca, cx, incx, cy, incy )

c*********************************************************************72
c
cc CAXPY computes constant times a vector plus a vector.
c
c  Discussion:
c
c    This routine uses single precision complex arithmetic.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Jack Dongarra
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of elements in CX and CY.
c
c    Input, complex CA, the multiplier of CX.
c
c    Input, complex CX(*), the first vector.
c
c    Input, integer INCX, the increment between successive entries of CX.
c
c    Input/output, complex CY(*), the second vector.
c    On output, CY(*) has been replaced by CY(*) + CA * CX(*).
c
c    Input, integer INCY, the increment between successive entries of CY.
c
      implicit none

      complex ca
      complex cx(*)
      complex cy(*)
      integer i
      integer incx
      integer incy
      integer ix
      integer iy
      integer n

      if ( n .le. 0 ) then
        return
      end if

      if ( abs ( real ( ca ) ) + abs ( aimag ( ca ) ) .eq. 0.0 ) then
        return
      end if
c
c  Code for both increments equal to 1.
c
      if ( incx .eq. 1 .and. incy .eq. 1 ) then

        do i = 1, n
          cy(i) = cy(i) + ca * cx(i)
        end do
c
c  Code for unequal increments or equal increments
c  not equal to 1.
c
      else

        if ( incx .lt. 0 ) then
          ix = ( - n + 1 ) * incx + 1
        else
          ix = 1
        end if

        if ( incy .lt. 0 ) then
          iy = ( - n + 1 ) * incy + 1
        else
          iy = 1
        end if

        do i = 1, n
          cy(iy) = cy(iy) + ca * cx(ix)
          ix = ix + incx
          iy = iy + incy
        end do

      end if

      return
      end
      subroutine ccopy ( n, cx, incx, cy, incy )

c*********************************************************************72
c
cc CCOPY copies a vector, x, to a vector, y.
c
c  Discussion:
c
c    This routine uses single precision complex arithmetic.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Jack Dongarra
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of elements in CX and CY.
c
c    Input, complex CX(*), the first vector.
c
c    Input, integer INCX, the increment between successive entries of CX.
c
c    Output, complex CY(*), the second vector.
c
c    Input, integer INCY, the increment between successive entries of CY.
c
      implicit none

      complex cx(*),cy(*)
      integer i,incx,incy,ix,iy,n

      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c  code for unequal increments or equal increments
c  not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        cy(iy) = cx(ix)
        ix = ix + incx
        iy = iy + incy
      end do
      return
c
c  code for both increments equal to 1
c
   20 do i = 1,n
        cy(i) = cx(i)
      end do

      return
      end
      function cdotc ( n, cx, incx, cy, incy )

c*********************************************************************72
c
cc CDOTC forms the conjugated  dot product of two vectors.
c
c  Discussion:
c
c    This routine uses single precision complex arithmetic.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Jack Dongarra
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vectors.
c
c    Input, complex CX(*), the first vector.
c
c    Input, integer INCX, the increment between successive entries in CX.
c
c    Input, complex CY(*), the second vector.
c
c    Input, integer INCY, the increment between successive entries in CY.
c
c    Output, complex CDOTC, the conjugated dot product of 
c    the corresponding entries of CX and CY.
c
      implicit none

      complex cdotc
      complex cx(*),cy(*),ctemp
      integer i,incx,incy,ix,iy,n

      ctemp = (0.0,0.0)
      cdotc = (0.0,0.0)
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c  code for unequal increments or equal increments
c  not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        ctemp = ctemp + conjg(cx(ix))*cy(iy)
        ix = ix + incx
        iy = iy + incy
      end do
      cdotc = ctemp
      return
c
c  code for both increments equal to 1
c
   20 do i = 1,n
        ctemp = ctemp + conjg(cx(i))*cy(i)
      end do

      cdotc = ctemp
      return
      end
      function cdotu ( n, cx, incx, cy, incy )

c*********************************************************************72
c
cc CDOTU forms the dot product of two vectors.
c
c  Discussion:
c
c    This routine uses single precision complex arithmetic.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Jack Dongarra
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vectors.
c
c    Input, complex CX(*), the first vector.
c
c    Input, integer INCX, the increment between successive entries in CX.
c
c    Input, complex CY(*), the second vector.
c
c    Input, integer INCY, the increment between successive entries in CY.
c
c    Output, complex CDOTU, the unconjugated dot product of 
c    the corresponding entries of CX and CY.
c
      implicit none

      complex cdotu
      complex cx(*),cy(*),ctemp
      integer i,incx,incy,ix,iy,n

      ctemp = (0.0,0.0)
      cdotu = (0.0,0.0)
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c  code for unequal increments or equal increments
c  not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        ctemp = ctemp + cx(ix)*cy(iy)
        ix = ix + incx
        iy = iy + incy
      end do
      cdotu = ctemp
      return
c
c  code for both increments equal to 1
c
   20 do i = 1,n
        ctemp = ctemp + cx(i)*cy(i)
      end do

      cdotu = ctemp

      return
      end
      function cmach ( job )

c*********************************************************************72
c
cc CMACH computes single precision complex machine parameters.
c
c  Discussion:
c
c    Assume the computer has
c
c      B = base of arithmetic;
c      T = number of base B digits;
c      L = smallest possible exponent;
c      U = largest possible exponent;
c
c    then
c
c      EPS = B**(1-T)
c      TINY = 100.0 * B**(-L+T)
c      HUGE = 0.01 * B**(U-T)
c
c    If complex division is done by
c
c      1 / (X+i*Y) = (X-i*Y) / (X**2+Y**2)
c
c    then
c
c      TINY = sqrt ( TINY )
c      HUGE = sqrt ( HUGE )
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer JOB:
c    1, EPS is desired;
c    2, TINY is desired;
c    3, HUGE is desired.
c
c    Output, real CMACH, the requested value.
c
      implicit none

      real cmach
      integer job
      real eps,tiny,huge,s

      eps = 1.0
   10 eps = eps/2.0
      s = 1.0 + eps
      if (s .gt. 1.0) go to 10
      eps = 2.0*eps
      cmach =eps
      if( job .eq. 1) return

      s = 1.0
   20 tiny = s
      s = s/16.0
      if (s*1.0 .ne. 0.0) go to 20
      tiny = (tiny/eps)*100.
      s = real((1.0,0.0)/cmplx(tiny,0.0))
      if (s .ne. 1.0/tiny) tiny = sqrt(tiny)
      huge = 1.0/tiny

      if ( job .eq. 1 ) then
        cmach = eps
      else if ( job .eq. 2 ) then
        cmach = tiny
      else if ( job .eq. 3 ) then
        cmach = huge
      end if

      return
      end
      subroutine crotg ( ca, cb, c, s )

c*********************************************************************72
c
cc CROTG generates a Givens rotation.
c
c  Discussion:
c
c    This routine uses single precision complex arithmetic.
c
c    Given values A and B, this routine computes:
c
c    If A = 0:
c
c      R = B
c      C = 0
c      S = (1,0).
c
c    If A /= 0:
c
c      ALPHA = A / abs ( A )
c      NORM  = sqrt ( ( abs ( A ) )**2 + ( abs ( B ) )**2 )
c      R     = ALPHA * NORM
c      C     = abs ( A ) / NORM
c      S     = ALPHA * conj ( B ) / NORM
c
c    In either case, the computed numbers satisfy the equation:
c
c    (         C    S ) * ( A ) = ( R )
c    ( -conj ( S )  C )   ( B ) = ( 0 )
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Jack Dongarra
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input/output, complex CA; on input, the value A.  On output,
c    the value R.
c
c    Input, complex CB, the value B.
c
c    Output, real C, the cosine of the Givens rotation.
c
c    Output, complex S, the sine of the Givens rotation.
c
      implicit none

      complex ca,cb,s
      real c
      real norm,scale
      complex alpha

      if ( cabs(ca) .eq. 0.0 ) then
         c = 0.
         s = (1.,0.)
         ca = cb
      else
         scale = cabs(ca) + cabs(cb)
         norm = scale * sqrt((cabs(ca/scale))**2 + (cabs(cb/scale))**2)
         alpha = ca /cabs(ca)
         c = cabs(ca) / norm
         s = alpha * conjg(cb) / norm
         ca = alpha * norm
      end if

      return
      end
      subroutine cscal ( n, ca, cx, incx )

c*********************************************************************72
c
cc CSCAL scales a vector by a constant.
c
c  Discussion:
c
c    This routine uses single precision complex arithmetic.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Jack Dongarra
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, complex CA, the multiplier.
c
c    Input/output, complex CX(*), the vector to be scaled.
c
c    Input, integer INCX, the increment between successive entries of CX.
c
      implicit none

      complex ca,cx(*)
      integer i,incx,n,nincx

      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c  code for increment not equal to 1
c
      nincx = n*incx
      do i = 1,nincx,incx
        cx(i) = ca*cx(i)
      end do
      return
c
c  code for increment equal to 1
c
   20 do i = 1,n
        cx(i) = ca*cx(i)
      end do

      return
      end
      function csign1 ( z1, z2 )

c*********************************************************************72
c
cc CSIGN1 is a single precision complex transfer-of-sign function.
c
c  Discussion:
c
c    The L1 norm is used.
c
c  Modified:
c
c    14 May 2004
c
c  Parameters:
c
c    Input, complex Z1, Z2, the arguments.
c
c    Output, complex CSIGN1,  a complex value, with the magnitude of
c    Z1, and the argument of Z2.
c
      implicit none

      real cabs1
      complex csign1
      complex z1
      complex z2

      if ( cabs1 ( z2 ) == 0.0E+00 ) then
        csign1 = ( 0.0E+00, 0.0E+00 )
      else
        csign1 = cabs1 ( z1 ) * ( z2 / cabs1 ( z2 ) )
      end if

      return
      end
      function csign2 ( z1, z2 )

c*********************************************************************72
c
cc CSIGN2 is a single precision complex transfer-of-sign function.
c
c  Discussion:
c
c    The L2 norm is used.
c
c  Modified:
c
c    19 March 2006
c
c  Parameters:
c
c    Input, complex Z1, Z2, the arguments.
c
c    Output, complex CSIGN2,  a complex value, with the magnitude of
c    Z1, and the argument of Z2.
c
      implicit none

      real cabs2
      complex csign2
      complex z1
      complex z2

      if ( cabs2 ( z2 ) == 0.0E+00 ) then
        csign2 = ( 0.0E+00, 0.0E+00 )
      else
        csign2 = cabs2 ( z1 ) * ( z2 / cabs2 ( z2 ) )
      end if

      return
      end
      subroutine csrot ( n, cx, incx, cy, incy, c, s )

c*********************************************************************72
c
cc CSROT applies a plane rotation.
c
c  Discussion:
c
c    The cos and sin (c and s) are real and the vectors cx and cy are complex.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Jack Dongarra
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vectors.
c
c    Input/output, complex CX(*), one of the vectors to be rotated.
c
c    Input, integer INCX, the increment between successive entries of CX.
c
c    Input/output, complex CY(*), one of the vectors to be rotated.
c
c    Input, integer INCY, the increment between successive elements of CY.
c
c    Input, real C, S, parameters (presumably the cosine and sine of
c    some angle) that define a plane rotation.
c
      implicit none

      complex cx(1),cy(1),ctemp
      real c,s
      integer i,incx,incy,ix,iy,n

      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        ctemp = c*cx(ix) + s*cy(iy)
        cy(iy) = c*cy(iy) - s*cx(ix)
        cx(ix) = ctemp
        ix = ix + incx
        iy = iy + incy
      end do

      return
c
c       code for both increments equal to 1
c
   20 do i = 1,n
        ctemp = c*cx(i) + s*cy(i)
        cy(i) = c*cy(i) - s*cx(i)
        cx(i) = ctemp
      end do
      return
      end
      subroutine csscal ( n, sa, cx, incx )

c*********************************************************************72
c
cc CSSCAL scales a complex vector by a real constant.
c
c  Discussion:
c
c    This routine uses single precision complex arithmetic.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Jack Dongarra
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, real SA, the multiplier.
c
c    Input/output, complex CX(*), the vector to be scaled.
c
c    Input, integer INCX, the increment between successive entries of 
c    the vector CX.
c
      implicit none

      complex cx(*)
      real sa
      integer i,incx,n,nincx

      if( n.le.0 .or. incx.le.0 )return

      if(incx.eq.1)go to 20
c
c  code for increment not equal to 1
c
      nincx = n*incx
      do i = 1,nincx,incx
        cx(i) = cmplx ( sa*real(cx(i)), sa*aimag(cx(i)) )
      end do
      return
c
c  code for increment equal to 1
c
   20 do i = 1,n
        cx(i) = cmplx ( sa*real(cx(i)), sa*aimag(cx(i)) )
      end do

      return
      end
      subroutine cswap ( n, cx, incx, cy, incy )

c*********************************************************************72
c
cc CSWAP interchanges two vectors.
c
c  Discussion:
c
c    This routine uses single precision complex arithmetic.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Jack Dongarra
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vectors.
c
c    Input/output, complex CX(*), one of the vectors to swap.
c
c    Input, integer INCX, the increment between successive entries of CX.
c
c    Input/output, complex CY(*), one of the vectors to swap.
c
c    Input, integer INCY, the increment between successive elements of CY.
c
      implicit none

      complex cx(*),cy(*),ctemp
      integer i,incx,incy,ix,iy,n

      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c  code for unequal increments or equal increments not equal
c  to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        ctemp = cx(ix)
        cx(ix) = cy(iy)
        cy(iy) = ctemp
        ix = ix + incx
        iy = iy + incy
      end do
      return
c
c       code for both increments equal to 1
c
   20 do i = 1,n
        ctemp = cx(i)
        cx(i) = cy(i)
        cy(i) = ctemp
      end do

      return
      end
      function icamax ( n, cx, incx )

c*********************************************************************72
c
cc ICAMAX finds the index of element having maximum absolute value.
c
c  Discussion:
c
c    This routine uses single precision complex arithmetic.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Jack Dongarra
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, complex X(*), the vector.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Output, integer ICAMAX, the index of the element of maximum
c    absolute value.
c
      implicit none

      complex cx(*)
      integer icamax
      real smax
      integer i,incx,ix,n
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))

      icamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      icamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      smax = cabs1(cx(1))
      ix = ix + incx
      do i = 2,n
         if(cabs1(cx(ix)).le.smax) go to 5
         icamax = i
         smax = cabs1(cx(ix))
    5    ix = ix + incx
      end do

      return
c
c        code for increment equal to 1
c
   20 smax = cabs1(cx(1))
      do i = 2,n
        if ( smax .lt. cabs1(cx(i)) ) then
          icamax = i
          smax = cabs1(cx(i))
        end if
      end do

      return
      end
      function lsame ( ca, cb )

c*********************************************************************72
c
cc LSAME returns TRUE if CA is the same letter as CB regardless of case.
c
c  Modified:
c
c    07 July 2007
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, character CA, CB, the character to compare.
c
c    Output, logical LSAME, is TRUE if the characters are equal,
c    disregarding case.
c
      implicit none

      character          ca, cb
      intrinsic          ichar
      logical lsame
      integer            inta, intb, zcode
c
c     Test if the characters are equal
c
      lsame = ca.eq.cb
      if( lsame ) then
        return
      end if
c
c     Now test for equivalence if both characters are alphabetic.
c
      zcode = ichar( 'Z' )
c
c     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
c     machines, on which ICHAR returns a value with bit 8 set.
c     ICHAR('A') on Prime machines returns 193 which is the same as
c     ICHAR('A') on an EBCDIC machine.
c
      inta = ichar( ca )
      intb = ichar( cb )
c
      if( zcode.eq.90 .or. zcode.eq.122 ) then
c
c        ASCII is assumed - ZCODE is the ASCII code of either lower or
c        upper case 'Z'.
c
         if( inta.ge.97 .and. inta.le.122 ) inta = inta - 32
         if( intb.ge.97 .and. intb.le.122 ) intb = intb - 32
c
      else if( zcode.eq.233 .or. zcode.eq.169 ) then
c
c        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
c        upper case 'Z'.
c
         if( inta.ge.129 .and. inta.le.137 .or.
     &       inta.ge.145 .and. inta.le.153 .or.
     &       inta.ge.162 .and. inta.le.169 ) inta = inta + 64
         if( intb.ge.129 .and. intb.le.137 .or.
     &       intb.ge.145 .and. intb.le.153 .or.
     &       intb.ge.162 .and. intb.le.169 ) intb = intb + 64
c
      else if( zcode.eq.218 .or. zcode.eq.250 ) then
c
c        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
c        plus 128 of either lower or upper case 'Z'.
c
         if( inta.ge.225 .and. inta.le.250 ) inta = inta - 32
         if( intb.ge.225 .and. intb.le.250 ) intb = intb - 32
      end if
      lsame = inta.eq.intb

      return
      end
      function scasum ( n, cx, incx )

c*********************************************************************72
c
cc SCASUM takes the sum of the absolute values of a complex vector.
c
c  Discussion:
c
c    This routine uses single precision complex arithmetic.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Jack Dongarra
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, complex X(*), the vector.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Output, real SCASUM, the sum of the absolute values.
c
      implicit none

      complex cx(*)
      real scasum
      real stemp
      integer i,incx,n,nincx

      scasum = 0.0e0
      stemp = 0.0e0
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do i = 1,nincx,incx
        stemp = stemp + abs(real(cx(i))) + abs(aimag(cx(i)))
      end do
      scasum = stemp
      return
c
c        code for increment equal to 1
c
   20 do i = 1,n
        stemp = stemp + abs(real(cx(i))) + abs(aimag(cx(i)))
      end do

      scasum = stemp

      return
      end
      function scnrm2 ( n, x, incx )

c*********************************************************************72
c
cc SCNRM2 returns the euclidean norm of a complex vector.
c
c  Discussion:
c
c    This routine uses single precision complex arithmetic.
c
c    SCNRM2 := sqrt ( sum ( conjg ( x(1:n) ) * x(1:n) ) )
c            = sqrt ( dot_product ( x(1:n), x(1:n) ) )
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Sven Hammarling
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input, complex X(*), the vector.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Output, real SCNRM2, the norm of the vector.
c
      implicit none

      integer incx
      integer ix
      integer n
      real one
      parameter ( one = 1.0e+0 )
      real scnrm2
      complex x( * )
      real zero
      parameter ( zero = 0.0e+0 )


      real                  norm, scale, ssq, temp
      intrinsic             abs, aimag, real, sqrt

      if( n.lt.1 .or. incx.lt.1 )then
         norm  = zero
      else
         scale = zero
         ssq   = one
c        The following loop is equivalent to this call to the LAPACK
c        auxiliary routine:
c        CALL CLASSQ( N, X, INCX, SCALE, SSQ )
c
         do ix = 1, 1 + ( n - 1 )*incx, incx
            if( real( x( ix ) ).ne.zero )then
               temp = abs( real( x( ix ) ) )
               if( scale.lt.temp )then
                  ssq   = one   + ssq*( scale/temp )**2
                  scale = temp
               else
                  ssq   = ssq   +     ( temp/scale )**2
               end if
            end if
            if( aimag( x( ix ) ).ne.zero )then
               temp = abs( aimag( x( ix ) ) )
               if( scale.lt.temp )then
                  ssq   = one   + ssq*( scale/temp )**2
                  scale = temp
               else
                  ssq   = ssq   +     ( temp/scale )**2
               end if
            end if
         end do
         norm  = scale * sqrt( ssq )
      end if

      scnrm2 = norm

      return
      end
      subroutine xerbla ( srname, info )

c*********************************************************************72
c
cc XERBLA is an error handler.
c
c  Discussion:
c
c    XERBLA is called if an input parameter has an invalid value.  
c    A message is printed and execution stops.
c
c  Modified:
c
c    07 July 2007
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User's Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh,
c    Basic Linear Algebra Subprograms for FORTRAN usage,
c    ACM Transactions on Mathematical Software,
c    Volume 5, Number 3, pages 308-323, 1979.
c
c  Parameters:
c
c    Input, character*6 SRNAME, the name of the routine which called XERBLA.
c
c    Input, integer INFO, the position of the invalid parameter in the 
c    parameter list of the calling routine.
c
      implicit none

      integer info
      character*6 srname

      write ( *, '(a,a6,a,i2,a)' ) ' ** On entry to ', srname, 
     &  ' parameter number ', info, ' had an illegal value.'

      stop
      end
