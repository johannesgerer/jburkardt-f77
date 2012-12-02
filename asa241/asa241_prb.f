      program main

c*********************************************************************72
c
cc MAIN is the main program for ASA241_PRB.
c
c  Discussion:
c
c    ASA241_PRB tests ASA241.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 January 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA241_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Test the routines in ASA241.'

      call test01
      call test02

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA241_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01

c*********************************************************************72
c
cc TEST01 tests PPND7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 January 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      real fx2
      integer ifault
      integer n_data
      real ppnd7
      double precision x
      real x2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) '  Let FX = NormalCDF ( X ).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  NORMAL_01_CDF_VALUES returns some values of ( X, FX ).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  PPND7 takes the value of FX, and'
      write ( *, '(a)' ) 
     &  '    computes an estimate X2, of the corresponding input,'
      write ( *, '(a)' ) 
     &  '    argument, accurate to about 7 decimal places.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) 
     &  '      FX                         X                         X2',
     &  '                     DIFF'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call normal_01_cdf_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = real ( fx )
        x2 = ppnd7 ( fx2, ifault )

        write ( *, '(2x,g24.16,2x,g24.16,2x,g24.16,2x,g10.4)' ) 
     &  fx, x, x2, abs ( x - x2 )

      go to 10

20    continue

      return
      end
      subroutine test02

c*********************************************************************72
c
cc TEST02 tests R8_NORMAL_01_CDF_INVERSE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 January 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision fx
      integer ifault
      integer n_data
      double precision ppnd16
      double precision x
      double precision x2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02:'
      write ( *, '(a)' ) '  Let FX = NormalCDF ( X ).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  NORMAL_01_CDF_VALUES returns some values of ( X, FX ).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  PPND16 takes the value of FX, and'
      write ( *, '(a)' ) 
     &  '    computes an estimate X2, of the corresponding input, '
      write ( *, '(a)' ) 
     &  '    argument, accurate to about 16 decimal places.'
      write ( *, '(a)' ) ' ' 
      write ( *, '(a,a)' ) 
     &  '      FX                         X                         X2',
     &  '                     DIFF'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call normal_01_cdf_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        x2 = ppnd16 ( fx, ifault )

        write ( *, '(2x,g24.16,2x,g24.16,2x,g24.16,2x,g10.4)' ) 
     &  fx, x, x2, dabs ( x - x2 )

      go to 10

20    continue

      return
      end
