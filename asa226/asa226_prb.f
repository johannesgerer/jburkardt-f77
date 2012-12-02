      program main

c*********************************************************************72
c
cc MAIN is the main program for ASA226_PRB.
c
c  Discussion:
c
c    ASA226_PRB calls the ASA226 routines.
c
c  Modified:
c
c    11 January 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA226_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the ASA226 library.'

      call test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA226_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests BETANC against tabulated values.
c
c  Modified:
c
c    11 January 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision b
      double precision betanc
      double precision fx
      double precision fx2
      integer ifault
      double precision lambda
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) '  Compare tabulated values of the noncentral'
      write ( *, '(a)' ) '  incomplete Beta Function against values'
      write ( *, '(a)' ) '  computed by BETANC.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '      A        B     LAMBDA        X       ',
     &  ' CDF               CDF          DIFF'
      write ( *, '(a,a)' )
     &  '                                           ',
     &  '(tabulated)       (BETANC)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call beta_noncentral_cdf_values ( n_data, a, b, lambda, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = betanc ( x, a, b, lambda, ifault )

        write ( *,
     &  '(2x,f7.1,2x,f7.1,2x,f7.1,2x,f10.4,
     &  2x,g14.6,2x,g14.6,2x,g10.4)' )
     &  a, b, lambda, x, fx, fx2, dabs ( fx - fx2 )

      go to 10

20    continue

      return
      end
