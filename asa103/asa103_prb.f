      program main

c*********************************************************************72
c
cc MAIN is the main program for ASA103_PRB.
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
      write ( *, '(a)' ) 'ASA103_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the ASA103 library.'

      call test01 ( * )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA103_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 compare DIGAMA against tabulated values.
c
c  Modified:
c
c    17 January 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision digama
      double precision fx
      double precision fx2
      integer ifault
      integer n_data
      double precision x

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) '  Compare tabulated values of the digamma'
      write ( *, '(a)' ) '  or Psi function against values computed'
      write ( *, '(a)' ) '  computed by DIGAMA.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '         X         Psi                     ',
     &  '  Psi                    DIFF'
      write ( *, '(a,a)' )
     &  '               (tabulated)                 ',
     &  '(DIGAMA)'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call psi_values ( n_data, x, fx )

        if ( n_data .eq. 0 ) then
          go to 20
        end if

        fx2 = digama ( x, ifault )

        write ( *, '(2x,f10.4,2x,g24.16,2x,g24.16,2x,g10.4)' )
     &  x, fx, fx2, dabs ( fx - fx2 )

      go to 10

20    continue

      return
      end
