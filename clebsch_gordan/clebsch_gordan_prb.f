      program main

c*********************************************************************72
c
cc MAIN is the main program for CLEBSCH_GORDAN_PRB.
c
c  Modified:
c
c    09 February 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CLEBSCH_GORDAN_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the CLEBSCH_GORDAN library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CLEBSCH_GORDAN_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests COF3J.
c
c  Modified:
c
c    07 February 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision cof3j
      double precision fx
      double precision j1
      double precision j2
      double precision j3
      double precision m1
      double precision m2
      double precision m3
      integer n_data
      double precision value

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' )
     &  '  COF3J evaluates the Wigner 3J coefficient.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '      J1      J2      J3      ',
     &  'M1      M2      M3        THREE_J'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call three_j_values ( n_data, j1, j2, j3, m1, m2, m3, fx )

        if ( n_data .le. 0 ) then
          go to 20
        end if

        write ( *,
     &  '(2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,g24.16)' )
     &  j1, j2, j3, m1, m2, m3, fx

        value = cof3j ( j1, j2, j3, m1, m2, m3 )

        write ( *,
     &  '(2x,6x,2x,6x,2x,6x,2x,6x,2x,6x,2x,6x,2x,g24.16)' )
     &  value

        write ( *, '(a)') ' '

      go to 10

20    continue

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests COF6J.
c
c  Modified:
c
c    07 February 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision cof6j
      double precision fx
      double precision j1
      double precision j2
      double precision j3
      double precision j4
      double precision j5
      double precision j6
      integer n_data
      double precision value

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' )
     &  '  COF6J evaluates the Wigner 6J coefficient.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '      J1      J2      J3      ',
     &  'J4      J5      J6        SIX_J'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call six_j_values ( n_data, j1, j2, j3, j4, j5, j6, fx )

        if ( n_data .le. 0 ) then
          go to 20
        end if

        write ( *,
     &  '(2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,g24.16)' )
     &  j1, j2, j3, j4, j5, j6, fx

        value = cof6j ( j1, j2, j3, j4, j5, j6 )

        write ( *,
     &  '(2x,6x,2x,6x,2x,6x,2x,6x,2x,6x,2x,6x,2x,g24.16)' )
     &  value

        write ( *, '(a)') ' '

      go to 10

20    continue

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests COF9J.
c
c  Modified:
c
c    09 February 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision cof9j
      double precision fx
      double precision j1
      double precision j2
      double precision j3
      double precision j4
      double precision j5
      double precision j6
      double precision j7
      double precision j8
      double precision j9
      integer n_data
      double precision value

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' )
     &  '  COF9J evaluates the Wigner 9J coefficient.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a,a)' )
     &  '      J1      J2      J3',
     &  '      J4      J5      J6',
     &  '      J7      J8      J9        NINE_J'
      write ( *, '(a)' ) ' '

      n_data = 0

10    continue

        call nine_j_values ( n_data, j1, j2, j3, j4, j5, j6,
     &  j7, j8, j9, fx )

        if ( n_data .le. 0 ) then
          go to 20
        end if

        write ( *,
     &  '(2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,f6.2,2x,
     &  f6.2,2x,f6.2,2x,f6.2,2x,g24.16)' )
     &  j1, j2, j3, j4, j5, j6, j7, j8, j9, fx

        value = cof9j ( j1, j2, j3, j4, j5, j6, j7, j8, j9 )

        write ( *,
     &  '(2x,6x,2x,6x,2x,6x,2x,6x,2x,6x,2x,6x,2x,6x,2x,6x,2x,
     &  6x,2x,g24.16)' )
     &  value

        write ( *, '(a)') ' '

      go to 10

20    continue

      return
      end
