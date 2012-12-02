      program main

c*********************************************************************72
c
cc MAIN is the main program for MACHAR_PRB.
c
c  Discussion:
c
c    MACHAR_PRB runs the MACHAR tests.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 November 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MACHAR_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the MACHAR library.'

      call test01 ( )
      call test02 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MACHAR_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests R4_MACHAR.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 June 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      real eps
      real epsneg
      integer ibeta
      integer iexp
      integer irnd
      integer it
      integer machep
      integer maxexp
      integer minexp
      integer negep
      integer ngrd
      real xmax
      real xmin

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  R4_MACHAR computes single'
      write ( *, '(a)' ) '  precision machine constants.'

      call r4_machar ( ibeta, it, irnd, ngrd, machep, negep, iexp, 
     &  minexp, maxexp, eps, epsneg, xmin, xmax )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  IBETA is the internal base for machine arithmetic.'
      write ( *, '(a,i8)' ) '    IBETA =  ', ibeta
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  IT is the number of digits, base IBETA, in the'
      write ( *, '(a)' ) '  floating point significand.'
      write ( *, '(a,i8)' ) '    IT =     ', it
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  IRND reports on floating point addition rounding:'
      write ( *, '(a)' ) '  0, for chopping;'
      write ( *, '(a)' ) '  1, for non-IEEE rounding;'
      write ( *, '(a)' ) '  2, for IEEE rounding;'
      write ( *, '(a)' ) '  3, for chopping with partial underflow;'
      write ( *, '(a)' ) 
     &  '  4, for non-IEEE rounding with partial underflow.'
      write ( *, '(a)' ) 
     &  '  5, for IEEE rounding with partial underflow.'
      write ( *, '(a,i8)' ) '    IRND =   ', irnd
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  NGRD is the number of guard digits for floating point'
      write ( *, '(a)' ) '  multiplication with truncating arithmetic.'
      write ( *, '(a,i8)' ) '    NGRD =   ', ngrd
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  MACHEP is the largest negative integer such that'
      write ( *, '(a)' ) '  1.0 < 1.0 + BETA^MACHEP.'
      write ( *, '(a,i8)' ) '    MACHEP = ', machep
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  NEGEPS is the largest negative integer such that'
      write ( *, '(a)' ) '  1.0 - BETA^NEGEPS < 1.0:'
      write ( *, '(a,i8)' ) '    NEGEP =  ', negep
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  IEXP is the number of bits reserved for the exponent'
      write ( *, '(a)' ) '  of a floating point number:'
      write ( *, '(a,i8)' ) '    IEXP =   ', iexp
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  MINEXP is the most negative power of BETA such that'
      write ( *, '(a)' ) '  BETA^MINEXP is positive and normalized.'
      write ( *, '(a,i8)' ) '    MINEXP = ', minexp
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  MAXEXP is the smallest positive power of BETA that'
      write ( *, '(a)' ) '  overflows:'
      write ( *, '(a,i8)' ) '    MAXEXP = ', maxexp
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  EPS is a small positive floating point number'
      write ( *, '(a)' ) '  such that 1.0 < 1.0 + EPS.'
      write ( *, '(a,e25.13)' ) '    EPS    = ', eps
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  EPSNEG is a small positive floating point number'
      write ( *, '(a)' ) '  such that 1.0 - EPSNEG < 1.0.'
      write ( *, '(a,e25.13)' ) '    EPSNEG = ', epsneg
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  XMIN is the smallest positive normalized floating'
      write ( *, '(a)' ) '  point power of the radix:'
      write ( *, '(a,e25.13)' ) '    XMIN =   ', xmin
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  XMAX is the largest finite floating point number:'
      write ( *, '(a,e25.13)' ) '    XMAX   = ', xmax

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Repeat floating point data using * format:'
      write ( *, '(a)' ) ' '
      write ( *, * ) '    EPS    = ', eps
      write ( *, * ) '    EPSNEG = ', epsneg
      write ( *, * ) '    XMIN   = ', xmin
      write ( *, * ) '    XMAX   = ', xmax

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests R8_MACHAR.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 June 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision eps
      double precision epsneg
      integer ibeta
      integer iexp
      integer irnd
      integer it
      integer machep
      integer maxexp
      integer minexp
      integer negep
      integer ngrd
      double precision xmax
      double precision xmin

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  R8_MACHAR computes double'
      write ( *, '(a)' ) '  precision machine constants.'

      call r8_machar ( ibeta, it, irnd, ngrd, machep, negep, iexp, 
     &  minexp, maxexp, eps, epsneg, xmin, xmax )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  IBETA is the internal base for machine arithmetic.'
      write ( *, '(a,i8)' ) '    IBETA =  ', ibeta
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  IT is the number of digits, base IBETA, in the'
      write ( *, '(a)' ) '  floating point significand.'
      write ( *, '(a,i8)' ) '    IT =     ', it
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  IRND reports on floating point addition rounding:'
      write ( *, '(a)' ) '  0, for chopping;'
      write ( *, '(a)' ) '  1, for non-IEEE rounding;'
      write ( *, '(a)' ) '  2, for IEEE rounding;'
      write ( *, '(a)' ) '  3, for chopping with partial underflow;'
      write ( *, '(a)' ) 
     &  '  4, for non-IEEE rounding with partial underflow.'
      write ( *, '(a)' ) 
     &  '  5, for IEEE rounding with partial underflow.'
      write ( *, '(a,i8)' ) '    IRND =   ', irnd
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  NGRD is the number of guard digits for floating point'
      write ( *, '(a)' ) '  multiplication with truncating arithmetic.'
      write ( *, '(a,i8)' ) '    NGRD =   ', ngrd
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  MACHEP is the largest negative integer such that'
      write ( *, '(a)' ) '  1.0 < 1.0 + BETA^MACHEP.'
      write ( *, '(a,i8)' ) '    MACHEP = ', machep
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  NEGEPS is the largest negative integer such that'
      write ( *, '(a)' ) '  1.0 - BETA^NEGEPS < 1.0:'
      write ( *, '(a,i8)' ) '    NEGEP =  ', negep
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  IEXP is the number of bits reserved for the exponent'
      write ( *, '(a)' ) '  of a floating point number:'
      write ( *, '(a,i8)' ) '    IEXP =   ', iexp
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  MINEXP is the most negative power of BETA such that'
      write ( *, '(a)' ) '  BETA^MINEXP is positive and normalized.'
      write ( *, '(a,i8)' ) '    MINEXP = ', minexp
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  MAXEXP is the smallest positive power of BETA that'
      write ( *, '(a)' ) '  overflows:'
      write ( *, '(a,i8)' ) '    MAXEXP = ', maxexp
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  EPS is a small positive floating point number'
      write ( *, '(a)' ) '  such that 1.0 < 1.0 + EPS.'
      write ( *, '(a,e25.13)' ) '    EPS    = ', eps
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  EPSNEG is a small positive floating point number'
      write ( *, '(a)' ) '  such that 1.0 - EPSNEG < 1.0.'
      write ( *, '(a,e25.13)' ) '    EPSNEG = ', epsneg
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  XMIN is the smallest positive normalized floating'
      write ( *, '(a)' ) '  point power of the radix:'
      write ( *, '(a,e25.13)' ) '    XMIN =   ', xmin
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  XMAX is the largest finite floating point number:'
      write ( *, '(a,e25.13)' ) '    XMAX   = ', xmax

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Repeat floating point data using * format:'
      write ( *, '(a)' ) ' '
      write ( *, * ) '    EPS    = ', eps
      write ( *, * ) '    EPSNEG = ', epsneg
      write ( *, * ) '    XMIN   = ', xmin
      write ( *, * ) '    XMAX   = ', xmax

      return
      end
