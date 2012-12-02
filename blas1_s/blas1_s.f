      function isamax ( n, sx, incx )

c*********************************************************************72
c
cc ISAMAX finds the index of element having maximum absolute value.
c
c  Discussion:
c
c    This routine uses single precision real arithmetic.
c
c  Modified:
c
c    25 October 2008
c
c  Author:
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh
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
c    Input, real X(*), the vector to be examined.
c
c    Input, integer INCX, the increment between successive entries of SX.
c
c    Output, integer ISAMAX, the index of the element of SX of maximum
c    absolute value.
c
      implicit none

      integer i
      integer incx
      integer isamax
      integer ix
      integer n
      real sx(*)
      real smax

      isamax = 0

      if ( n .lt. 1 .or. incx .le. 0 ) then
        return
      end if

      isamax = 1
      if ( n .eq. 1 ) then
        return
      end if

      if ( incx .ne. 1 ) then

        ix = 1
        smax = abs ( sx(1) )
        ix = ix + incx
        do i = 2, n
          if ( smax .lt. abs ( sx(ix) ) ) then
            isamax = i
            smax = abs ( sx(ix) )
          end if
          ix = ix + incx
        end do

      else

        smax = abs ( sx(1) )
        do i = 2, n
          if ( smax .lt. abs ( sx(i) ) ) then
            isamax = i
            smax = abs ( sx(i) )
          end if
        end do

      end if

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
c  Author:
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh
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
      logical lsame
      intrinsic          ichar

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
      ZCODE = ICHAR( 'Z' )
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
      function sasum ( n, sx, incx )

c*********************************************************************72
c
cc SASUM takes the sum of the absolute values.
c
c  Discussion:
c
c    This routine uses single precision real arithmetic.
c
c    This routine uses unrolled loops for increment equal to one.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh
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
c    Input, real X(*), the vector to be examined.
c
c    Input, integer INCX, the increment between successive entries of X.
c    INCX must not be negative.
c
c    Output, real SASUM, the sum of the absolute values of X.
c
      implicit none

      real sasum
      real sx(*),stemp
      integer i,incx,m,n,nincx

      sasum = 0.0e0
      stemp = 0.0e0
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do i = 1,nincx,incx
        stemp = stemp + abs(sx(i))
      end do
      sasum = stemp
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        stemp = stemp + abs(sx(i))
      end do
      if( n .lt. 6 ) go to 60

   40 continue

      do i = m+1, n, 6
        stemp = stemp + abs(sx(i)) + abs(sx(i + 1)) + abs(sx(i + 2))
     &  + abs(sx(i + 3)) + abs(sx(i + 4)) + abs(sx(i + 5))
      end do

   60 sasum = stemp

      return
      end
      subroutine saxpy ( n, sa, sx, incx, sy, incy )

c*********************************************************************72
c
cc SAXPY computes constant times a vector plus a vector.
c
c  Discussion:
c
c    This routine uses single precision real arithmetic.
c
c    This routine uses unrolled loop for increments equal to one.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh
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
c    Input, real X(*), the vector to be scaled and added to Y.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Input/output, real Y(*), the vector to which a multiple of X is to 
c    be added.
c
c    Input, integer INCY, the increment between successive entries of Y.
c
      implicit none

      real sx(*),sy(*),sa
      integer i,incx,incy,ix,iy,m,n

      if(n.le.0)return
      if (sa .eq. 0.0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        sy(iy) = sy(iy) + sa*sx(ix)
        ix = ix + incx
        iy = iy + incy
      end do
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        sy(i) = sy(i) + sa*sx(i)
      end do
      if( n .lt. 4 ) return
   40 continue

      do i = m+1, n, 4
        sy(i) = sy(i) + sa*sx(i)
        sy(i + 1) = sy(i + 1) + sa*sx(i + 1)
        sy(i + 2) = sy(i + 2) + sa*sx(i + 2)
        sy(i + 3) = sy(i + 3) + sa*sx(i + 3)
      end do

      return
      end
      subroutine scopy ( n, sx, incx, sy, incy )

c*********************************************************************72
c
cc SCOPY copies a vector, x, to a vector, y.
c
c  Discussion:
c
c    This routine uses single precision real arithmetic.
c
c    This routine uses unrolled loops for increments equal to 1.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh
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
c    Input, real X(*), the vector to be copied into Y.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Output, real Y(*), the copy of X.
c
c    Input, integer INCY, the increment between successive elements of Y.
c
      implicit none

      real sx(*),sy(*)
      integer i,incx,incy,ix,iy,m,n

      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        sy(iy) = sx(ix)
        ix = ix + incx
        iy = iy + incy
      end do

      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        sy(i) = sx(i)
      end do

      if( n .lt. 7 ) return
   40 continue

      do i = m+1, n, 7
        sy(i) = sx(i)
        sy(i + 1) = sx(i + 1)
        sy(i + 2) = sx(i + 2)
        sy(i + 3) = sx(i + 3)
        sy(i + 4) = sx(i + 4)
        sy(i + 5) = sx(i + 5)
        sy(i + 6) = sx(i + 6)
      end do

      return
      end
      function sdot ( n, sx, incx, sy, incy )

c*********************************************************************72
c
cc SDOT forms the dot product of two vectors.
c
c  Discussion:
c
c    This routine uses single precision real arithmetic.
c
c    This routine uses unrolled loops for increments equal to one.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh
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
c    Input, real X(*), one of the vectors to be multiplied.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Input, real Y(*), one of the vectors to be multiplied.
c
c    Input, integer INCY, the increment between successive elements of Y.
c
c    Output, real SDOT, the dot product of X and Y.
c
      implicit none

      real sdot
      real sx(*),sy(*),stemp
      integer i,incx,incy,ix,iy,m,n

      stemp = 0.0e0
      sdot = 0.0e0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do i = 1,n
        stemp = stemp + sx(ix)*sy(iy)
        ix = ix + incx
        iy = iy + incy
      end do
      sdot = stemp
      return
c
c  code for both increments equal to 1
c
c
c  clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        stemp = stemp + sx(i)*sy(i)
      end do

      if( n .lt. 5 ) go to 60
   40 continue

      do i = m+1, n, 5
        stemp = stemp + sx(i)*sy(i) + sx(i + 1)*sy(i + 1) +
     &   sx(i + 2)*sy(i + 2) + sx(i + 3)*sy(i + 3) + sx(i + 4)*sy(i + 4)
      end do

   60 sdot = stemp

      return
      end
      REAL FUNCTION SDSDOT(N,SB,SX,INCX,SY,INCY)

cc SDSDOT computes the inner product of two vectors with extended precision.
c
*     .. Scalar Arguments ..
      REAL SB
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      REAL SX(*),SY(*)
*     ..
*
*  PURPOSE
*  =======
*
*  Compute the inner product of two vectors with extended
*  precision accumulation.
*
*  Returns S.P. result with dot product accumulated in D.P.
*  SDSDOT = SB + sum for I = 0 to N-1 of SX(LX+I*INCX)*SY(LY+I*INCY),
*  where LX = 1 if INCX .GE. 0, else LX = 1+(1-N)*INCX, and LY is
*  defined in a similar way using INCY.
*
*  AUTHOR
*  ======
*  Lawson, C. L., (JPL), Hanson, R. J., (SNLA),
*  Kincaid, D. R., (U. of Texas), Krogh, F. T., (JPL)
*
*  ARGUMENTS 
*  =========
*
*  N      (input) INTEGER
*         number of elements in input vector(s)
*
*  SB     (input) REAL
*         single precision scalar to be added to inner product
*
*  SX     (input) REAL array, dimension (N)
*         single precision vector with N elements
*
*  INCX   (input) INTEGER
*         storage spacing between elements of SX
*
*  SY     (input) REAL array, dimension (N)
*         single precision vector with N elements
*
*  INCY   (input) INTEGER
*         storage spacing between elements of SY
*
*  SDSDOT (output) REAL
*         single precision dot product (SB if N .LE. 0)
*
*  Further Details
*  ===============
*
*  REFERENCES
*
*  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
*  Krogh, Basic linear algebra subprograms for Fortran
*  usage, Algorithm No. 539, Transactions on Mathematical
*  Software 5, 3 (September 1979), pp. 308-323.
*
*  REVISION HISTORY  (YYMMDD)
*      
*  791001  DATE WRITTEN
*  890531  Changed all specific intrinsics to generic.  (WRB)
*  890831  Modified array declarations.  (WRB)
*  890831  REVISION DATE from Version 3.2
*  891214  Prologue converted to Version 4.0 format.  (BAB)
*  920310  Corrected definition of LX in DESCRIPTION.  (WRB)
*  920501  Reformatted the REFERENCES section.  (WRB)
*  070118  Reformat to LAPACK coding style
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION DSDOT
      INTEGER I,KX,KY,NS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC DBLE
*     ..
      DSDOT = SB
      IF (N.LE.0) THEN
         SDSDOT = DSDOT
         RETURN
      END IF   
      IF (INCX.EQ.INCY .AND. INCX.GT.0) THEN
*
*     Code for equal and positive increments.
*
         NS = N*INCX
         DO I = 1,NS,INCX
            DSDOT = DSDOT + DBLE(SX(I))*DBLE(SY(I))
         END DO
      ELSE
*
*     Code for unequal or nonpositive increments.
*
         KX = 1
         KY = 1
         IF (INCX.LT.0) KX = 1 + (1-N)*INCX
         IF (INCY.LT.0) KY = 1 + (1-N)*INCY
         DO I = 1,N
            DSDOT = DSDOT + DBLE(SX(KX))*DBLE(SY(KY))
            KX = KX + INCX
            KY = KY + INCY
         END DO
      END IF
      SDSDOT = DSDOT
      RETURN
      END
      function smach ( job )

c*********************************************************************72
c
cc SMACH computes machine parameters of floating point arithmetic.
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
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh
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
c    Output, real SMACH, the requested value.
c
      implicit none

      integer job

      real eps,tiny,huge,s
      real smach

      eps = 1.0
   10 eps = eps/2.0
      s = 1.0 + eps
      if (s .gt. 1.0) go to 10
      eps = 2.0*eps

      s = 1.0
   20 tiny = s
      s = s/16.0
      if (s*100. .ne. 0.0) go to 20
      tiny = (tiny/eps)*100.0
      huge = 1.0/tiny

      if (job .eq. 1) smach = eps
      if (job .eq. 2) smach = tiny
      if (job .eq. 3) smach = huge

      return
      end
      function snrm2 ( n, x, incx )

c*********************************************************************72
c
cc SNRM2 returns the euclidean norm of a real vector.
c
c  Discussion:
c
c    This routine uses single precision real arithmetic.
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
c    Input, real X(*), the vector whose norm is to be computed.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Output, real SNRM2, the Euclidean norm of X.
c
      implicit none

      integer                           incx, n
      real snrm2
      real                              x( * )

      real                  one         , zero
      parameter           ( one = 1.0e+0, zero = 0.0e+0 )

      integer               ix
      real                  absxi, norm, scale, ssq

      intrinsic             abs, sqrt

      if( n.lt.1 .or. incx.lt.1 )then
         norm  = zero
      else if( n.eq.1 )then
         norm  = abs( x( 1 ) )
      else
         scale = zero
         ssq   = one
c        The following loop is equivalent to this call to the LAPACK
c        auxiliary routine:
c        CALL SLASSQ( N, X, INCX, SCALE, SSQ )
c
         do ix = 1, 1 + ( n - 1 )*incx, incx
            if( x( ix ).ne.zero )then
               absxi = abs( x( ix ) )
               if( scale.lt.absxi )then
                  ssq   = one   + ssq*( scale/absxi )**2
                  scale = absxi
               else
                  ssq   = ssq   +     ( absxi/scale )**2
               end if
            end if
        end do
        norm  = scale * sqrt( ssq )
      end if

      snrm2 = norm
      return
      end
      subroutine srot ( n, sx, incx, sy, incy, c, s )

c*********************************************************************72
c
cc SROT applies a plane rotation.
c
c  Discussion:
c
c    This routine uses single precision real arithmetic.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh
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
c    Input/output, real X(*), one of the vectors to be rotated.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Input/output, real Y(*), one of the vectors to be rotated.
c
c    Input, integer INCY, the increment between successive elements of Y.
c
c    Input, real C, S, parameters (presumably the cosine and sine of
c    some angle) that define a plane rotation.
c
      implicit none

      real sx(*),sy(*),stemp,c,s
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
        stemp = c*sx(ix) + s*sy(iy)
        sy(iy) = c*sy(iy) - s*sx(ix)
        sx(ix) = stemp
        ix = ix + incx
        iy = iy + incy
      end do
      return
c
c       code for both increments equal to 1
c
   20 do i = 1,n
        stemp = c*sx(i) + s*sy(i)
        sy(i) = c*sy(i) - s*sx(i)
        sx(i) = stemp
      end do

      return
      end
      subroutine srotg ( sa, sb, c, s )

c*********************************************************************72
c
cc SROTG constructs a Givens plane rotation.
c
c  Discussion:
c
c    This routine uses single precision real arithmetic.
c
c    Given values A and B, this routine computes
c
c    SIGMA = sign ( A ) if abs ( A ) >  abs ( B )
c          = sign ( B ) if abs ( A ) <= abs ( B );
c
c    R     = SIGMA * ( A * A + B * B );
c
c    C = A / R if R is not 0
c      = 1     if R is 0;
c
c    S = B / R if R is not 0,
c        0     if R is 0.
c
c    The computed numbers then satisfy the equation
c
c    (  C  S ) ( A ) = ( R )
c    ( -S  C ) ( B ) = ( 0 )
c
c    The routine also computes
c
c    Z = S     if abs ( A ) > abs ( B ),
c      = 1 / C if abs ( A ) <= abs ( B ) and C is not 0,
c      = 1     if C is 0.
c
c    The single value Z encodes C and S, and hence the rotation:
c
c    If Z = 1, set C = 0 and S = 1;
c    If abs ( Z ) < 1, set C = sqrt ( 1 - Z * Z ) and S = Z;
c    if abs ( Z ) > 1, set C = 1/ Z and S = sqrt ( 1 - C * C );
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh
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
c    Input/output, real SA, SB.  On input, SA and SB are the values
c    A and B.  On output, SA is overwritten with R, and SB is
c    overwritten with Z.
c
c    Output, real C, S, the cosine and sine of the Givens rotation.
c
      implicit none

      real sa,sb,c,s,roe,scale,r,z

      roe = sb
      if( abs(sa) .gt. abs(sb) ) roe = sa
      scale = abs(sa) + abs(sb)

      if( scale .eq. 0.0 ) then
         c = 1.0
         s = 0.0
         r = 0.0
         z = 0.0
      else
        r = scale*sqrt((sa/scale)**2 + (sb/scale)**2)
        r = sign(1.0,roe)*r
        c = sa/r
        s = sb/r
        z = 1.0
        if( abs(sa) .gt. abs(sb) ) z = s
        if( abs(sb) .ge. abs(sa) .and. c .ne. 0.0 ) z = 1.0/c
      end if

      sa = r
      sb = z

      return
      end
      SUBROUTINE SROTM(N,SX,INCX,SY,INCY,SPARAM)

cc SROTM applies a modified Givens transformation matrix.
c
*     .. Scalar Arguments ..
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      REAL SPARAM(5),SX(*),SY(*)
*     ..
*
*  Purpose
*  =======
*
*     APPLY THE MODIFIED GIVENS TRANSFORMATION, H, TO THE 2 BY N MATRIX
*
*     (SX**T) , WHERE **T INDICATES TRANSPOSE. THE ELEMENTS OF SX ARE IN
*     (SX**T)
*
*     SX(LX+I*INCX), I = 0 TO N-1, WHERE LX = 1 IF INCX .GE. 0, ELSE
*     LX = (-INCX)*N, AND SIMILARLY FOR SY USING USING LY AND INCY.
*     WITH SPARAM(1)=SFLAG, H HAS ONE OF THE FOLLOWING FORMS..
*
*     SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0
*
*       (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0)
*     H=(          )    (          )    (          )    (          )
*       (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0).
*     SEE  SROTMG FOR A DESCRIPTION OF DATA STORAGE IN SPARAM.
*
*
*  Arguments
*  =========
*
*  N      (input) INTEGER
*         number of elements in input vector(s)
*
*  SX     (input/output) REAL array, dimension N
*         double precision vector with N elements
*
*  INCX   (input) INTEGER
*         storage spacing between elements of SX
*
*  SY     (input/output) REAL array, dimension N
*         double precision vector with N elements
*
*  INCY   (input) INTEGER
*         storage spacing between elements of SY
*
*  SPARAM (input/output)  REAL array, dimension 5
*     SPARAM(1)=SFLAG
*     SPARAM(2)=SH11
*     SPARAM(3)=SH21
*     SPARAM(4)=SH12
*     SPARAM(5)=SH22
*
*  =====================================================================
*
*     .. Local Scalars ..
      REAL SFLAG,SH11,SH12,SH21,SH22,TWO,W,Z,ZERO
      INTEGER I,KX,KY,NSTEPS
*     ..
*     .. Data statements ..
      DATA ZERO,TWO/0.E0,2.E0/
*     ..
*
      SFLAG = SPARAM(1)
      IF (N.LE.0 .OR. (SFLAG+TWO.EQ.ZERO)) RETURN
      IF (INCX.EQ.INCY.AND.INCX.GT.0) THEN
*
         NSTEPS = N*INCX
         IF (SFLAG.LT.ZERO) THEN
            SH11 = SPARAM(2)
            SH12 = SPARAM(4)
            SH21 = SPARAM(3)
            SH22 = SPARAM(5)
            DO I = 1,NSTEPS,INCX
               W = SX(I)
               Z = SY(I)
               SX(I) = W*SH11 + Z*SH12
               SY(I) = W*SH21 + Z*SH22
            END DO
         ELSE IF (SFLAG.EQ.ZERO) THEN
            SH12 = SPARAM(4)
            SH21 = SPARAM(3)
            DO I = 1,NSTEPS,INCX
               W = SX(I)
               Z = SY(I)
               SX(I) = W + Z*SH12
               SY(I) = W*SH21 + Z
            END DO
         ELSE
            SH11 = SPARAM(2)
            SH22 = SPARAM(5)
            DO I = 1,NSTEPS,INCX
               W = SX(I)
               Z = SY(I)
               SX(I) = W*SH11 + Z
               SY(I) = -W + SH22*Z
            END DO
         END IF
      ELSE
         KX = 1
         KY = 1
         IF (INCX.LT.0) KX = 1 + (1-N)*INCX
         IF (INCY.LT.0) KY = 1 + (1-N)*INCY
*
         IF (SFLAG.LT.ZERO) THEN
            SH11 = SPARAM(2)
            SH12 = SPARAM(4)
            SH21 = SPARAM(3)
            SH22 = SPARAM(5)
            DO I = 1,N
               W = SX(KX)
               Z = SY(KY)
               SX(KX) = W*SH11 + Z*SH12
               SY(KY) = W*SH21 + Z*SH22
               KX = KX + INCX
               KY = KY + INCY
            END DO
         ELSE IF (SFLAG.EQ.ZERO) THEN
            SH12 = SPARAM(4)
            SH21 = SPARAM(3)
            DO I = 1,N
               W = SX(KX)
               Z = SY(KY)
               SX(KX) = W + Z*SH12
               SY(KY) = W*SH21 + Z
               KX = KX + INCX
               KY = KY + INCY
            END DO
         ELSE
             SH11 = SPARAM(2)
             SH22 = SPARAM(5)
             DO I = 1,N
                W = SX(KX)
                Z = SY(KY)
                SX(KX) = W*SH11 + Z
                SY(KY) = -W + SH22*Z
                KX = KX + INCX
                KY = KY + INCY
            END DO
         END IF
      END IF
      RETURN
      END
      SUBROUTINE SROTMG(SD1,SD2,SX1,SY1,SPARAM)

cc SROTMG constructs a modified Givens transformation matrix
c
*     .. Scalar Arguments ..
      REAL SD1,SD2,SX1,SY1
*     ..
*     .. Array Arguments ..
      REAL SPARAM(5)
*     ..
*
*  Purpose
*  =======
*
*     CONSTRUCT THE MODIFIED GIVENS TRANSFORMATION MATRIX H WHICH ZEROS
*     THE SECOND COMPONENT OF THE 2-VECTOR  (SQRT(SD1)*SX1,SQRT(SD2)*
*     SY2)**T.
*     WITH SPARAM(1)=SFLAG, H HAS ONE OF THE FOLLOWING FORMS..
*
*     SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0
*
*       (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0)
*     H=(          )    (          )    (          )    (          )
*       (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0).
*     LOCATIONS 2-4 OF SPARAM CONTAIN SH11,SH21,SH12, AND SH22
*     RESPECTIVELY. (VALUES OF 1.E0, -1.E0, OR 0.E0 IMPLIED BY THE
*     VALUE OF SPARAM(1) ARE NOT STORED IN SPARAM.)
*
*     THE VALUES OF GAMSQ AND RGAMSQ SET IN THE DATA STATEMENT MAY BE
*     INEXACT.  THIS IS OK AS THEY ARE ONLY USED FOR TESTING THE SIZE
*     OF SD1 AND SD2.  ALL ACTUAL SCALING OF DATA IS DONE USING GAM.
*
*
*  Arguments
*  =========
*
*
*  SD1    (input/output) REAL
*
*  SD2    (input/output) REAL
*
*  SX1    (input/output) REAL
*
*  SY1    (input) REAL
*
*
*  SPARAM (input/output)  REAL array, dimension 5
*     SPARAM(1)=SFLAG
*     SPARAM(2)=SH11
*     SPARAM(3)=SH21
*     SPARAM(4)=SH12
*     SPARAM(5)=SH22
*
*  =====================================================================
*
*     .. Local Scalars ..
      REAL GAM,GAMSQ,ONE,RGAMSQ,SFLAG,SH11,SH12,SH21,SH22,SP1,SP2,SQ1,
     $     SQ2,STEMP,SU,TWO,ZERO
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC ABS
*     ..
*     .. Data statements ..
*
      DATA ZERO,ONE,TWO/0.E0,1.E0,2.E0/
      DATA GAM,GAMSQ,RGAMSQ/4096.E0,1.67772E7,5.96046E-8/
*     ..

      IF (SD1.LT.ZERO) THEN
*        GO ZERO-H-D-AND-SX1..
         SFLAG = -ONE
         SH11 = ZERO
         SH12 = ZERO
         SH21 = ZERO
         SH22 = ZERO
*
         SD1 = ZERO
         SD2 = ZERO
         SX1 = ZERO
      ELSE
*        CASE-SD1-NONNEGATIVE
         SP2 = SD2*SY1
         IF (SP2.EQ.ZERO) THEN
            SFLAG = -TWO
            SPARAM(1) = SFLAG
            RETURN
         END IF 
*        REGULAR-CASE..
         SP1 = SD1*SX1
         SQ2 = SP2*SY1
         SQ1 = SP1*SX1
*
         IF (ABS(SQ1).GT.ABS(SQ2)) THEN
            SH21 = -SY1/SX1
            SH12 = SP2/SP1
*
            SU = ONE - SH12*SH21
*
           IF (SU.GT.ZERO) THEN
             SFLAG = ZERO
             SD1 = SD1/SU
             SD2 = SD2/SU
             SX1 = SX1*SU
           END IF
         ELSE

            IF (SQ2.LT.ZERO) THEN
*              GO ZERO-H-D-AND-SX1..
               SFLAG = -ONE
               SH11 = ZERO
               SH12 = ZERO
               SH21 = ZERO
               SH22 = ZERO
*
               SD1 = ZERO
               SD2 = ZERO
               SX1 = ZERO
            ELSE
               SFLAG = ONE
               SH11 = SP1/SP2
               SH22 = SX1/SY1
               SU = ONE + SH11*SH22
               STEMP = SD2/SU
               SD2 = SD1/SU
               SD1 = STEMP
               SX1 = SY1*SU
            END IF
         END IF

*     PROCESURE..SCALE-CHECK
         IF (SD1.NE.ZERO) THEN
            DO WHILE ((SD1.LE.RGAMSQ) .OR. (SD1.GE.GAMSQ))
               IF (SFLAG.EQ.ZERO) THEN
                  SH11 = ONE
                  SH22 = ONE
                  SFLAG = -ONE
               ELSE
                  SH21 = -ONE
                  SH12 = ONE
                  SFLAG = -ONE
               END IF
               IF (SD1.LE.RGAMSQ) THEN
                  SD1 = SD1*GAM**2
                  SX1 = SX1/GAM
                  SH11 = SH11/GAM
                  SH12 = SH12/GAM
               ELSE
                  SD1 = SD1/GAM**2
                  SX1 = SX1*GAM
                  SH11 = SH11*GAM
                  SH12 = SH12*GAM
               END IF
            ENDDO
         END IF
  
         IF (SD2.NE.ZERO) THEN
            DO WHILE ( (ABS(SD2).LE.RGAMSQ) .OR. (ABS(SD2).GE.GAMSQ) )
               IF (SFLAG.EQ.ZERO) THEN
                  SH11 = ONE
                  SH22 = ONE
                  SFLAG = -ONE
               ELSE
                  SH21 = -ONE
                  SH12 = ONE
                  SFLAG = -ONE
               END IF
               IF (ABS(SD2).LE.RGAMSQ) THEN
                  SD2 = SD2*GAM**2
                  SH21 = SH21/GAM
                  SH22 = SH22/GAM
               ELSE
                  SD2 = SD2/GAM**2
                  SH21 = SH21*GAM
                  SH22 = SH22*GAM
               END IF      
            END DO
         END IF
     
      END IF

      IF (SFLAG.LT.ZERO) THEN
         SPARAM(2) = SH11
         SPARAM(3) = SH21
         SPARAM(4) = SH12
         SPARAM(5) = SH22
      ELSE IF (SFLAG.EQ.ZERO) THEN
         SPARAM(3) = SH21
         SPARAM(4) = SH12 
      ELSE
         SPARAM(2) = SH11
         SPARAM(5) = SH22
      END IF

  260 CONTINUE
      SPARAM(1) = SFLAG
      RETURN
      END
      subroutine sscal ( n, sa, sx, incx )

c*********************************************************************72
c
cc SSCAL scales a vector by a constant.
c
c  Discussion:
c
c    This routine uses single precision real arithmetic.
c
c    This routine uses unrolled loops for increment equal to 1.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh
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
c    Input/output, real X(*), the vector to be scaled.
c
c    Input, integer INCX, the increment between successive entries of X.
c
      implicit none

      real sa,sx(*)
      integer i,incx,m,n,nincx

      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c  code for increment not equal to 1
c
      nincx = n*incx
      do i = 1,nincx,incx
        sx(i) = sa*sx(i)
      end do
      return
c
c  code for increment equal to 1
c
c
c  clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        sx(i) = sa*sx(i)
      end do
      if( n .lt. 5 ) return
   40 continue

      do i = m+1, n, 5
        sx(i) = sa*sx(i)
        sx(i + 1) = sa*sx(i + 1)
        sx(i + 2) = sa*sx(i + 2)
        sx(i + 3) = sa*sx(i + 3)
        sx(i + 4) = sa*sx(i + 4)
      end do

      return
      end
      subroutine sswap ( n, sx, incx, sy, incy )

c*********************************************************************72
c
cc SSWAP interchanges two vectors.
c
c  Discussion:
c
c    This routine uses single precision real arithmetic.
c
c    This routine uses unrolled loops for increments equal to 1.
c
c  Modified:
c
c    07 July 2007
c
c  Author:
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh
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
c    Input/output, real X(*), one of the vectors to swap.
c
c    Input, integer INCX, the increment between successive entries of X.
c
c    Input/output, real Y(*), one of the vectors to swap.
c
c    Input, integer INCY, the increment between successive elements of Y.
c
      implicit none

      real sx(*),sy(*),stemp
      integer i,incx,incy,ix,iy,m,n

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
        stemp = sx(ix)
        sx(ix) = sy(iy)
        sy(iy) = stemp
        ix = ix + incx
        iy = iy + incy
      end do
      return
c
c       code for both increments equal to 1
c
c
c       clean-up loop
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do i = 1,m
        stemp = sx(i)
        sx(i) = sy(i)
        sy(i) = stemp
      end do
      if( n .lt. 3 ) return
   40 continue

      do i = m+1, n, 3
        stemp = sx(i)
        sx(i) = sy(i)
        sy(i) = stemp
        stemp = sx(i + 1)
        sx(i + 1) = sy(i + 1)
        sy(i + 1) = stemp
        stemp = sx(i + 2)
        sx(i + 2) = sy(i + 2)
        sy(i + 2) = stemp
      end do

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
c  Author:
c
c    Charles Lawson, Richard Hanson, David Kincaid, Fred Krogh
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
