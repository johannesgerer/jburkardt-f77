      program main

c*********************************************************************72
c
cc MAIN is the main program for COMBINATION_LOCK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 May 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of dials.
c
c    Input, integer N, the number of symbols on each dial.
c    We assume the symbols are the integers 0 to N-1.
c
c    Input, integer C(M), the combination.
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'COMBINATION_LOCK'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the COMBINATION_LOCK libary.'

      call test01 ( )
      call test02 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'COMBINATION_LOCK'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      return
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests BICYCLE_LOCK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 May 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer c
      integer i4_uniform
      integer m
      parameter ( m = 3 )
      integer n
      parameter ( n = 10 )
      integer seed
      integer step
      
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) 
     &  '  A bicycle combination lock consists of 3 dials,'
      write ( *, '(a)' ) '  each having 10 symbols, 0 through 9.'
      write ( *, '(a)' ) '  We seek to determine the combination C.'
c
c  Report on the problem data.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  The number of dials is M = ', m
      write ( *, '(a,i4)' ) '  The number of symbols is N = ', n
      write ( *, '(a,i8)' ) 
     &  '  The number of possible combinations is M^N = ', n ** m

      call get_seed ( seed )
      c = i4_uniform ( 0, 999, seed )

      write ( *, '(a,i3)' ) '  The "secret" combination is ', c

      call bicycle_lock ( c, step )

      if ( step .eq. -1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  The combination was not found!'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) 
     &    '  The combination was found on step ', step
      end if

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests COMBINATION_LOCK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 May 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, the number of dials.
c
c    Input, integer N, the number of symbols on each dial.
c    We assume the symbols are the integers 0 to N-1.
c
c    Input, integer C(M), the combination.
c
      implicit none

      integer m
      parameter ( m = 4 )

      integer c(m)
      integer n
      parameter ( n = 5 )
      integer step

      save c

      data c / 1, 2, 3, 4 /
      
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  A combination lock consists of M dials,'
      write ( *, '(a)' ) '  each having N symbols.'
      write ( *, '(a)' ) '  We seek to determine the combination C.'
c
c  Report on the problem data.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  The number of dials is M = ', m
      write ( *, '(a,i4)' ) '  The number of symbols is N = ', n
      write ( *, '(a,i8)' ) 
     &  '  The number of possible combinations is M^N = ', n ** m

      call i4vec_print ( m, c, '  The "secret" combination:' );

      call combination_lock ( m, n, c, step )

      if ( step .eq. -1 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  The combination was not found!'
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) 
     &    '  The combination was found on step ', step
      end if

      return
      end
