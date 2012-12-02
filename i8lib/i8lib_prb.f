      program main

c*********************************************************************72
c
cc MAIN is the main program for I8LIB_PRB.
c
c  Discussion:
c
c    I8LIB_PRB tests routines from the I8LIB library.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 June 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'I8LIB_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the routines in the I8LIB library.'

      call test015 ( )
      call test190 ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'I8LIB_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test015 ( )

c*********************************************************************72
c
cc TEST015 tests I8_CHOOSE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 June 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer*8 cnk
      integer*8 i8_choose
      integer*8 k
      integer*8 n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST015'
      write ( *, '(a)' ) '  I8_CHOOSE evaluates C(N,K).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       N       K      CNK'
      write ( *, '(a)' ) ' '

      do n = 0, 4
        do k = 0, n
          cnk = i8_choose ( n, k )
          write ( *, '(2x,i8,i8,i8)' ) n, k, cnk
        end do
      end do

      return
      end
      subroutine test190 ( )

c*********************************************************************72
c
cc TEST190 tests I8_XOR.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 June 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer*8 i
      integer*8 ihi
      integer*8 ilo
      integer*8 i8_uniform
      integer*8 i8_xor
      integer*8 j
      integer*8 k
      integer*8 l
      integer*8 seed
      integer*8 test
      integer*8 test_num

      ihi = 100
      ilo = 0
      test_num = 10

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST190'
      write ( *, '(a)' ) '  I8_XOR returns the bitwise exclusive OR of'
      write ( *, '(a)' ) '  two integers.'
      write ( *, '(a)' ) '  Compare the FORTRAN intrinsic IEOR.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '         I         J    I8_XOR      IEOR'
      write ( *, '(a)' ) ' '

      do test = 1, test_num
        i = i8_uniform ( ilo, ihi, seed )
        j = i8_uniform ( ilo, ihi, seed )
        k = i8_xor ( i, j )
        l = ieor ( i, j )
        write ( *, '(2x,i8,2x,i8,2x,i8,2x,i8)' ) i, j, k, l
      end do

      return
      end