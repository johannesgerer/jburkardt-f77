      program main

c*********************************************************************72
c
cc COMPLEX_NUMBERS is a program which demonstrates the use of complex numbers.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 November 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'COMPLEX_NUMBERS:'
      write ( *, '(a)' ) '  FORTRAN77 version.'
      write ( *, '(a)' ) '  Demonstrate complex number usage.'
c
c  Single precision complex.
c
      call test01 ( )
      call test02 ( )
      call test03 ( )
c
c  Double precision complex.
c
      call test04 ( )
      call test05 ( )
      call test06 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'COMPLEX_NUMBERS:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 declaration and assignment.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 November 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none
c
c  Declare a complex number A.
c  Declare a complex vector B.
c  Declare a complex array C.
c
      complex a
      complex b(3)
      complex c(2,2)
      integer i

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Declare a COMPLEX variable.'
      write ( *, '(a)' ) '  Assign value with an = statement.'
c
c  Assign values to A, B, and C.
c
      a = ( 1.0, 2.0 )

      b(1) = ( 1.0, 2.0 )
      b(2) = ( 3.0, 4.0 )
      b(3) = ( 5.0, 6.0 )

      c(1,1) = ( 1.0, 0.1 )
      c(2,1) = ( 2.0, 0.1 )
      c(1,2) = ( 1.0, 0.2 )
      c(2,2) = ( 2.0, 0.2 )
c
c  Print them.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Scalar A:'
      write ( *, '(a)' ) ' '

      print *, a
      write ( *, * ) a
      write ( *, '(2g14.6)' ) a

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Vector B:'
      write ( *, '(a)' ) ' '

      do i = 1, 3
        write ( *, '(2g14.6)' ) b(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Array C:'
      write ( *, '(a)' ) ' '

      do i = 1, 2
        write ( *, '(2(2g14.6))' ) c(i,1:2)
      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02: declaration with initialization.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 November 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none
c
c  Declare and initialize a complex number A.
c  Declare and initialize a complex vector B.
c  Declare and initialize a complex array C.
c
      complex a
      complex b(3)
      complex c(2,2)
      integer i

      save a
      save b
      save c

      data a / ( 1.0, 2.0 ) /
      data b / ( 1.0, 2.0 ), ( 3.0, 4.0 ), ( 5.0, 6.0 ) /
      data c / ( 1.0, 0.1 ), ( 2.0, 0.1 ), ( 1.0, 0.2 ), ( 2.0, 0.2 ) /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Declare a COMPLEX variable.'
      write ( *, '(a)' ) '  Initialize with a data statement.'
c
c  Print them.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Scalar A:'
      write ( *, '(a)' ) ' '

      print *, a
      write ( *, * ) a
      write ( *, '(2g14.6)' ) a

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Vector B:'
      write ( *, '(a)' ) ' '

      do i = 1, 3
        write ( *, '(2g14.6)' ) b(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Array C:'
      write ( *, '(a)' ) ' '

      do i = 1, 2
        write ( *, '(2(2g14.6))' ) c(i,1:2)
      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03: intrinsic functions.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 November 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex a

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) 
     &  '  Apply intrinsic functions to COMPLEX variable'

      a = ( 1.0, 2.0 )
c
c  Print them.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a,2g14.6)' ) '  a =              ', a
      write ( *, '(a,2g14.6)' ) '  - a =            ', - a
      write ( *, '(a,2g14.6)' ) '  a + 3 =          ', a + 3
      write ( *, '(a,2g14.6)' ) '  a + (0,5) =      ', a + ( 0, 5 )
      write ( *, '(a,2g14.6)' ) '  4 * a =          ', 4 * a
      write ( *, '(a,2g14.6)' ) '  a / 8 =          ', a / 8
      write ( *, '(a,2g14.6)' ) '  a * a =          ', a * a
      write ( *, '(a,2g14.6)' ) '  a**2 =           ', a**2
      write ( *, '(a,2g14.6)' ) '  1/a =            ', 1.0 / a
      write ( *, '(a)' ) ' '
      write ( *, '(a, g14.6)' ) '  cabs(a) =        ', cabs ( a )
      write ( *, '(a,2g14.6)' ) '  ccos(a) =        ', ccos ( a )
      write ( *, '(a,2g14.6)' ) '  cexp(a) =        ', cexp ( a )
      write ( *, '(a,2g14.6)' ) '  clog(a) =        ', clog ( a )
      write ( *, '(a,2g14.6)' ) '  cmplx(1) =       ', cmplx ( 1 )
      write ( *, '(a,2g14.6)' ) '  cmplx(2,3) =     ', cmplx ( 2, 3 )
      write ( *, '(a,2g14.6)' ) '  cmplx(4.0) =     ', cmplx ( 4.0 )
      write ( *, '(a,2g14.6)' ) 
     &  '  cmplx(5.0,6.0) = ', cmplx ( 5.0, 6.0 )
      write ( *, '(a,2g14.6)' ) '  conjg(a) =       ', conjg ( a )
      write ( *, '(a,2g14.6)' ) '  csin(a) =        ', csin ( a )
      write ( *, '(a,2g14.6)' ) '  csqrt(a) =       ', csqrt ( a )
      write ( *, '(a, g14.6)' ) '  imag(a) =        ', imag ( a )
      write ( *, '(a, i8)'    ) '  int(a) =         ', int ( a )
      write ( *, '(a, g14.6)' ) '  real(a) =        ', real ( a )

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 declaration and assignment for double precision complex variables.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 November 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none
c
c  Declare a double precision complex number A.
c  Declare a double precision complex vector B.
c  Declare a double precision complex array C.
c
      double complex a
      double complex b(3)
      double complex c(2,2)
      integer i

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  Declare a DOUBLE COMPLEX variable.'
      write ( *, '(a)' ) '  Assign value with an = statement.'
c
c  Assign values to A, B, and C.
c
      a = ( 1.0D+00, 2.0D+00 )

      b(1) = ( 1.0D+00, 2.0D+00 )
      b(2) = ( 3.0D+00, 4.0D+00 )
      b(3) = ( 5.0D+00, 6.0D+00 )

      c(1,1) = ( 1.0D+00, 0.1D+00 )
      c(2,1) = ( 2.0D+00, 0.1D+00 )
      c(1,2) = ( 1.0D+00, 0.2D+00 )
      c(2,2) = ( 2.0D+00, 0.2D+00 )
c
c  Print them.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Scalar A:'
      write ( *, '(a)' ) ' '

      print *, a
      write ( *, * ) a
      write ( *, '(2g14.6)' ) a

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Vector B:'
      write ( *, '(a)' ) ' '

      do i = 1, 3
        write ( *, '(2g14.6)' ) b(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Array C:'
      write ( *, '(a)' ) ' '

      do i = 1, 2
        write ( *, '(2(2g14.6))' ) c(i,1:2)
      end do

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05: declaration with initialization for double precision complex variables.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 November 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none
c
c  Declare and initialize a double complex number A.
c  Declare and initialize a double complex vector B.
c  Declare and initialize a double complex array C.
c
      double complex a
      double complex b(3)
      double complex c(2,2)
      integer i

      save a
      save b
      save c

      data a / ( 1.0D+00, 2.0D+00 ) /
      data b / ( 1.0D+00, 2.0D+00 ), 
     &         ( 3.0D+00, 4.0D+00 ), 
     &         ( 5.0D+00, 6.0D+00 ) /
      data c / ( 1.0D+00, 0.1D+00 ), 
     &         ( 2.0D+00, 0.1D+00 ), 
     &         ( 1.0D+00, 0.2D+00 ), 
     &         ( 2.0D+00, 0.2D+00 ) /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  Declare a DOUBLE COMPLEX variable.'
      write ( *, '(a)' ) '  Initialize with a data statement.'
c
c  Print them.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Scalar A:'
      write ( *, '(a)' ) ' '

      print *, a
      write ( *, * ) a
      write ( *, '(2g14.6)' ) a

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Vector B:'
      write ( *, '(a)' ) ' '

      do i = 1, 3
        write ( *, '(2g14.6)' ) b(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Array C:'
      write ( *, '(a)' ) ' '

      do i = 1, 2
        write ( *, '(2(2g14.6))' ) c(i,1:2)
      end do

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST03: intrinsic functions for double precision complex variables.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 November 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double complex a

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) 
     &  '  Apply intrinsic functions to DOUBLE COMPLEX variable'

      a = ( 1.0D+00, 2.0D+00 )
c
c  Print them.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a,2g14.6)' ) '  a =              ', a
      write ( *, '(a,2g14.6)' ) '  - a =            ', - a
      write ( *, '(a,2g14.6)' ) '  a + 3 =          ', a + 3
      write ( *, '(a,2g14.6)' ) '  a + (0,5) =      ', a + ( 0, 5 )
      write ( *, '(a,2g14.6)' ) '  4 * a =          ', 4 * a
      write ( *, '(a,2g14.6)' ) '  a / 3 =          ', a / 3
      write ( *, '(a,2g14.6)' ) '  a * a =          ', a * a
      write ( *, '(a,2g14.6)' ) '  a**2 =           ', a**2
      write ( *, '(a,2g14.6)' ) '  1/a =            ', 1.0 / a
      write ( *, '(a)' ) ' '
      write ( *, '(a, g14.6)' ) '  cdabs(a) =        ', cdabs ( a )
      write ( *, '(a,2g14.6)' ) '  cdcos(a) =        ', cdcos ( a )
      write ( *, '(a,2g14.6)' ) '  cdexp(a) =        ', cdexp ( a )
      write ( *, '(a,2g14.6)' ) '  cdlog(a) =        ', cdlog ( a )
      write ( *, '(a,2g14.6)' ) '  dcmplx(1) =       ', dcmplx ( 1 )
      write ( *, '(a,2g14.6)' ) '  dcmplx(2,3) =     ', dcmplx ( 2, 3 )
      write ( *, '(a,2g14.6)' ) '  dcmplx(4.0) =     ', dcmplx ( 4.0 )
      write ( *, '(a,2g14.6)' ) 
     &  '  dcmplx(5.0,6.0) = ', dcmplx ( 5.0, 6.0 )
      write ( *, '(a,2g14.6)' ) '  dconjg(a) =       ', dconjg ( a )
      write ( *, '(a,2g14.6)' ) '  cdsin(a) =        ', cdsin ( a )
      write ( *, '(a,2g14.6)' ) '  cdsqrt(a) =       ', cdsqrt ( a )
      write ( *, '(a, g14.6)' ) '  dimag(a) =        ', dimag ( a )
      write ( *, '(a, i8)'    ) '  int(a) =          ', int ( a )
      write ( *, '(a, g14.6)' ) '  dreal(a) =        ', dreal ( a )

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
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
