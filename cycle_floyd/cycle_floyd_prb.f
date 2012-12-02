      program main

c*********************************************************************72
c
cc MAIN is the main program for CYCLE_FLOYD_PRB.
c
c  Discussion:
c
c    CYCLE_FLOYD_PRB tests the CYCLE_FLOYD library.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 June 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CYCLE_FLOYD_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the CYCLE_FLOYD library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CYCLE_FLOYD_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests CYCLE_FLOYD for a tiny example.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 June 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer f1
      external f1
      integer lam
      integer mu
      integer x0

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Test CYCLE_FLOYD on F1().'
      write ( *, '(a)' ) '  f1(0) = 6.'
      write ( *, '(a)' ) '  f1(1) = 6.'
      write ( *, '(a)' ) '  f1(2) = 0.'
      write ( *, '(a)' ) '  f1(3) = 1.'
      write ( *, '(a)' ) '  f1(4) = 4.'
      write ( *, '(a)' ) '  f1(5) = 3.'
      write ( *, '(a)' ) '  f1(6) = 3.'
      write ( *, '(a)' ) '  f1(7) = 4.'
      write ( *, '(a)' ) '  f1(8) = 0.'

      x0 = 2
      write ( *, '(a)' ) ' '
      write ( *, '(a,i3)' ) '  Starting argument X0 = ', x0

      call cycle_floyd ( f1, x0, lam, mu )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i3)' ) '  Reported cycle length is ', lam
      write ( *, '(a)' ) '  Expected value is 3'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i3)' ) 
     &  '  Reported distance to first cycle element is ', mu
      write ( *, '(a)' ) '  Expected value is 2'

      return
      end
      function f1 ( i )

c*********************************************************************72
c
cc F1 is the iteration function for example 1.
c
c  Discussion:
c
c    This function has two cycles:
c
c    6, 3, 1, of length 3
c    4, of length 1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
      implicit none

      integer f1
      integer f_table(0:8)
      integer i

      save f_table

      data f_table / 6, 6, 0, 1, 4, 3, 3, 4, 0 /

      f1 = f_table ( i )

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests CYCLE_FLOYD for F2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 June 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer f2
      external f2
      integer lam
      integer mu
      integer x0

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Test CYCLE_FLOYD for F2().'
      write ( *, '(a)' ) '  f2(i) = mod ( 22 * i + 1, 72 ).'

      x0 = 0
      write ( *, '(a)' ) ' '
      write ( *, '(a,i3)' ) '  Starting argument X0 = ', x0

      call cycle_floyd ( f2, x0, lam, mu )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i3)' ) '  Reported cycle length is ', lam
      write ( *, '(a)' ) '  Expected value is 9'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i3)' ) 
     &  '  Reported distance to first cycle element is ', mu
      write ( *, '(a)' ) '  Expected value is 3'

      return
      end
      function f2 ( i )

c*********************************************************************72
c
cc F2 is the iteration function for example 2.
c
c  Discussion:
c
c    This function has a cycle
c
c    3, 67, 35, 51, 43, 11, 27, 19, 59, of length 9
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
      implicit none

      integer f2
      integer i

      f2 = mod ( 22 * i + 1, 72 )

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests CYCLE_FLOYD for F3.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 June 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer f3
      external f3
      integer lam
      integer mu
      integer x0

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  Test CYCLE_FLOYD for F3().'
      write ( *, '(a)' ) '  f3(i) = mod ( 123 * i + 456, 100000 ).'

      x0 = 789
      write ( *, '(a)' ) ' '
      write ( *, '(a,i3)' ) '  Starting argument X0 = ', x0

      call cycle_floyd ( f3, x0, lam, mu )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Reported cycle length is ', lam
      write ( *, '(a)' ) '  Expected value is 50000'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  Reported distance to first cycle element is ', mu
      write ( *, '(a)' ) '  Expected value is 0'

      return
      end
      function f3 ( i )

c*********************************************************************72
c
cc F3 is the iteration function for example 3.
c
c  Discussion:
c
c    This function has a cycle of length 50000
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
      implicit none

      integer f3
      integer i

      f3 = mod ( 123 * i + 456, 1000000 )

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests CYCLE_FLOYD for F4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 June 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer f4
      external f4
      integer lam
      integer mu
      integer x0

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  Test CYCLE_FLOYD for F4().'
      write ( *, '(a)' ) '  f4(i) = mod ( 31421 * i + 6927, 65536 ).'

      x0 = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a,i3)' ) '  Starting argument X0 = ', x0

      call cycle_floyd ( f4, x0, lam, mu )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Reported cycle length is ', lam
      write ( *, '(a)' ) '  Expected value is 65536'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  Reported distance to first cycle element is ', mu
      write ( *, '(a)' ) '  Expected value is 0'

      return
      end
      function f4 ( i )

c*********************************************************************72
c
cc F4 is the iteration function for example 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
      implicit none

      integer f4
      integer i

      f4 = mod ( 31421 * i + 6927, 65536 )

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests CYCLE_FLOYD for F5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 June 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer f5
      external f5
      integer i
      integer lam
      integer mu
      integer x0

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  Test CYCLE_FLOYD for F5().'
      write ( *, '(a)' ) '  f5(i) = mod ( 16383 * i + 1, 65536 ).'

      x0 = 1
      write ( *, '(a)' ) ' '
      write ( *, '(a,i3)' ) '  Starting argument X0 = ', x0

      call cycle_floyd ( f5, x0, lam, mu )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Reported cycle length is ', lam
      write ( *, '(a)' ) '  Expected value is 8'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  Reported distance to first cycle element is ', mu
      write ( *, '(a)' ) '  Expected value is 0'

      i = 0
      x0 = 1
      write ( *, * ) i, x0
      do i = 1, 10
        x0 = f5 ( x0 )
        write ( *, * ) i, x0
      end do

      return
      end
      function f5 ( i )

c*********************************************************************72
c
cc F5 is the iteration function for example 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 June 2012
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
      implicit none

      integer f5
      integer i

      f5 = mod ( 16383 * i + 1, 65536 )

      return
      end
