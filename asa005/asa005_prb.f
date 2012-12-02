      program main

c*********************************************************************72
c
cc MAIN is the main program for ASA005_PRB.
c
c  Discussion:
c
c    ASA005_PRB calls the ASA005 routines.
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

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA005_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the ASA005 library.'

      call test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA005_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 demonstrates the use of PRNCST.
c
c  Modified:
c
c    12 January 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer df
      double precision fx
      double precision fx2
      integer ifault
      double precision lambda
      integer n_data
      double precision prncst
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) '  PRNCST computes the noncentral Student T '
      write ( *, '(a)' ) '  Cumulative Density Function.'
      write ( *, '(a)' ) '  Compare with tabulated values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) '        X         LAMBDA        DF     ',
     &  ' CDF             CDF           DIFF'
      write ( *, '(a,a)' ) '                                       ',
     &  ' Tabulated       PRNCST'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call student_noncentral_cdf_values ( n_data, df, lambda,
     &    x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = prncst ( x, df, lambda, ifault )

        write ( *,
     &  '(2x,f10.4,2x,f10.4,2x,i8,2x,g14.6,2x,g14.6,2x,g10.4)' )
     &  x, lambda, df, fx, fx2, dabs ( fx - fx2 )

      go to 10

20    continue

      return
      end
