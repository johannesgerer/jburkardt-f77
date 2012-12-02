      program main

c*********************************************************************72
c
cc MAIN is the main program for TOMS370_PRB.
c
c  Modified:
c
c    05 February 2011
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Edgar Butler,
c    Algorithm 370:
c    General Random Number Generator,
c    Communications of the ACM,
c    Volume 13, Number 1, January 1970, pages 49-52.
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS370_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TOMS370 library.'

      call test01 ( )
      call test02 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TOMS370_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests RANDGI and RANDG for the normal PDF.
c
c  Modified:
c
c    05 February 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 21 )

      integer i
      integer ier
      integer k
      integer l
      integer l_in
      real mu
      real p(n)
      real pi
      parameter ( pi = 3.141592653589793E+00 )
      real q(257)
      real r(256)
      real sigma
      real xdata(n)
      real y
      real ydata(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) '  RANDGI can generate data from a PDF'
      write ( *, '(a)' ) '  if K is not equal to 1.'
      write ( *, '(a)' ) '  Here we work with the normal PDF.'
c
c  Generate N samples of normal PDF between -1 and +1
c
      mu = 0.0E+00
      sigma = 1.0E+00

      do i = 1, n
        xdata(i) = ( real ( n - i ) * ( -1.0E+00 )
     &             + real (     i ) * ( +1.0E+00 ) )
     &             / real ( n - 1 )
        ydata(i) =
     &    exp ( - ( xdata(i) - mu )**2 / ( 2.0E+00 * sigma**2 ) )
     &    / ( sqrt ( 2.0E+00 * pi ) * sigma )
      end do

      k = 0

      call randgi ( n, xdata, ydata, p, q, r, k, ier )

      if ( ier .eq. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST01 - Warning!'
        write ( *, '(a,i8)' ) '  RANDGI returned IER = ', ier
        write ( *, '(a)' ) '  The PDF data does not integrate to 1.'
        write ( *, '(a)' ) '  A scaling factor will be used.'
      else if ( ier .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST01 - Error!'
        write ( *, '(a,i8)' ) '  RANDGI returned IER = ', ier
        return
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I      Seed(in)     Seed(out)      Y'
      write ( *, '(a)' ) ' '

      l = 123456789
      do i = 1, 100
        l_in = l
        call randg ( l, q, r, y )
        write ( *, '(2x,i4,2x,i12,2x,i12,2x,g14.6)' ) i, l_in, l, y
      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests RANDGI and RANDG for the uniform PDF.
c
c  Modified:
c
c    05 February 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 3 )

      integer i
      integer ier
      integer k
      integer l
      integer l_in
      real p(n)
      real q(257)
      real r(256)
      real xdata(n)
      real y
      real ydata(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02:'
      write ( *, '(a)' ) '  RANDGI can generate data from a PDF'
      write ( *, '(a)' ) '  if K is not equal to 1.'
      write ( *, '(a)' ) '  Here we work with the uniform PDF.'
c
c  Generate 3 samples of uniform PDF between -1 and +1
c
      xdata(1) = -1.0E+00
      ydata(1) =  0.5E+00
      xdata(2) =  0.0E+00
      ydata(2) =  0.5E+00
      xdata(3) = +1.0E+00
      ydata(3) =  0.5E+00

      k = 0

      call randgi ( n, xdata, ydata, p, q, r, k, ier )

      if ( ier .eq. 1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST02 - Warning!'
        write ( *, '(a,i8)' ) '  RANDGI returned IER = ', ier
        write ( *, '(a)' ) '  The PDF data does not integrate to 1.'
        write ( *, '(a)' ) '  A scaling factor will be used.'
      else if ( ier .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST02 - Error!'
        write ( *, '(a,i8)' ) '  RANDGI returned IER = ', ier
        return
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I      Seed(in)     Seed(out)      Y'
      write ( *, '(a)' ) ' '

      l = 123456789
      do i = 1, 100
        l_in = l
        call randg ( l, q, r, y )
        write ( *, '(2x,i4,2x,i12,2x,i12,2x,g14.6)' ) i, l_in, l, y
      end do

      return
      end
