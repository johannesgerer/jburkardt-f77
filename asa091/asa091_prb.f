      program main

c*********************************************************************72
c
cc MAIN is the main program for ASA091_PRB.
c
c  Modified:
c
c    04 January 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      write ( *, '(a)' ) ' '
      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA091_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the ASA091 library.'

      call test01 ( )
      call test02 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA091_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 makes a single simple calculation with PPCHI2.
c
c  Modified:
c
c    02 January 2008
c
c  Author:
c
c    John Burkardt
c
      double precision alngam
      double precision g
      integer ifault
      double precision p
      double precision ppchi2
      double precision v
      double precision value
      double precision value_correct
      parameter ( value_correct = 0.4D+00 )

      p = 0.017523D+00
      v = 4.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) '  Perform a simple sample calculation using'
      write ( *, '(a)' ) '  PPCHI2 to invert the Chi-Squared CDF.'

      g = alngam ( v / 2.0D+00, ifault )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g24.16)' ) '  P =                  ', p
      write ( *, '(a,g24.16)' ) '  V =                  ', v
      write ( *, '(a,g24.16)' ) '  G Log(Gamma(V/2)) =  ', g

      value = ppchi2 ( p, v, g, ifault )

      write ( *, '(a,g24.16)' ) '  VALUE =              ', value
      write ( *, '(a,g24.16)' ) '  VALUE (correct) =    ', value_correct

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)'  ) '  Error flag IFAULT = ', ifault

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 compare PPCHI2 against tabulated values.
c
c  Modified:
c
c    21 January 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer a
      double precision alngam
      double precision fx
      integer n_data
      double precision x
      double precision x2
      double precision g
      double precision ppchi2
      integer ifault
      double precision v

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02:'
      write ( *, '(a)' ) '  Compare tabulated values of the Chi-Squared'
      write ( *, '(a)' ) '  Cumulative Density Function against values'
      write ( *, '(a)' ) '  computed by PPCHI2.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '         N        CDF       X                        ',
     &  ' X2                      DIFF'
      write ( *, '(a,a)' )
     &  '                          (tabulated)                ',
     &  '(PPCHI2)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call chi_square_cdf_values ( n_data, a, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        v = dble ( a )

        g = alngam ( dble ( a ) / 2.0D+00, ifault )

        x2 = ppchi2 ( fx, v, g, ifault )

        write ( *, '(2x,i8,2x,f10.4,2x,g24.16,2x,g24.16,2x,g10.4)' )
     &  a, fx, x, x2, dabs ( x - x2 )

      go to 10

20    continue

      return
      end
