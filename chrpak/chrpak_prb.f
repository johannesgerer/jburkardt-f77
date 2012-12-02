      program main

c*********************************************************************72
c
cc MAIN is the main program for CHRPAK_PRB.
c
c  Discussion:
c
c    CHRPAK_PRB tests routines from the CHRPAK library.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CHRPAK_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the CHRPAK library.'

      call test001 ( )
      call test021 ( )
      call test029 ( )
      call test0325 ( )
      call test108 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CHRPAK_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test001 ( )

c*********************************************************************72
c
cc TEST001 tests A_TO_I4 and I4_TO_A.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 March 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      character a
      integer a_to_i4
      integer i
      integer i2
      character i4_to_a

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST001'
      write ( *, '(a)' ) '  A_TO_I4: Alphabetic character => I'
      write ( *, '(a)' ) '  I4_TO_A: I => Alphabetic character'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  1:26 = A:Z'
      write ( *, '(a)' ) '  27:52 = a:z'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   I  ==>  A  ==>  I'
      write ( *, '(a)' ) ' '

      do i = 0, 55, 3

        a = i4_to_a ( i )
        i2 = a_to_i4 ( a )

        write ( *, '(i8,5x,a1,5x,i8)' ) i, a, i2

      end do

      return
      end
      subroutine test021 ( )

c*********************************************************************72
c
cc TEST021 tests CH_TO_DIGIT and DIGIT_TO_CH.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 November 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      character c
      integer i
      integer i2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST021'
      write ( *, '(a)' ) '  CH_TO_DIGIT: character -> decimal digit'
      write ( *, '(a)' ) '  DIGIT_TO_C: decimal digit -> character.'
      write ( *, '(a)' ) ' '

      do i = -2, 11
        call digit_to_ch ( i, c )
        call ch_to_digit ( c, i2 )
        write ( *, '(2x,i6,a6,i6)' ) i, c, i2
      end do

      return
      end
      subroutine test029 ( )

c*********************************************************************72
c
cc TEST029 tests CH_UNIFORM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 November 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer a_to_i4
      character ch
      character ch_uniform
      character chi
      character clo
      integer count(26)
      integer i
      character i4_to_a
      integer j
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST029'
      write ( *, '(a)' ) '  CH_UNIFORM returns a random character.'
      write ( *, '(a)' ) ' '

      do i = 1, 26
        count(i) = 0
      end do

      clo = 'D'
      chi = 'W'
      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '   I  A  Count'
      write ( *, '(a)' ) ' '

      do i = 1, 100000

        ch = ch_uniform ( clo, chi, seed )

        j = a_to_i4 ( ch )
        count(j) = count(j) + 1

      end do

      do i = 1, 26
        write ( *, '(2x,i2,2x,a1,2x,i5)' ) i, i4_to_a(i), count(i)
      end do

      return
      end
      subroutine test0325

c*********************************************************************72
c
cc TEST0325 tests FILE_NAME_INC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    15 November 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      character ( len = 30 ) s

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST0325'
      write ( *, '(a)' ) '  FILE_NAME_INC can "increment" the numeric'
      write ( *, '(a)' ) '  part of a file name.'

      s = 'data01.txt'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Starting string: "' // trim ( s ) // '"'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Incremented forms:'
      write ( *, '(a)' ) ' '

      do i = 1, 5
        call file_name_inc ( s )
        write ( *, '(2x,a)' ) '  "' // trim ( s ) // '"'
      end do

      s = 'mat09lab98.txt'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Starting string: "' // trim ( s ) // '"'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Incremented forms:'
      write ( *, '(a)' ) ' '

      do i = 1, 5
        call file_name_inc ( s )
        write ( *, '(2x,a)' ) '  "' // trim ( s ) // '"'
      end do

      return
      end
      subroutine test108 ( )

c*********************************************************************72
c
cc TEST108 tests S_INDEX_LAST and S_INDEXI.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    04 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i1
      integer i2
      integer i3
      integer i4
      integer s_indexi
      integer s_index_last
      character * ( 30 ) s
      character * ( 10 ) substring

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST108'
      write ( *, '(a)' ) '  S_INDEXI reports the first occurrence of a'
      write ( *, '(a)' ) '    substring, case and trailing space'
      write ( *, '(a)' ) '    insensitive.'
      write ( *, '(a)' ) '  S_INDEX_LAST reports the LAST occurrence'
      write ( *, '(a)' ) '    of a substring.'
      write ( *, '(a)' ) 
     &  '  INDEX is a case and trailing space sensitive'
      write ( *, '(a)' ) 
     &  '    routine which reports the first occurrence'
      write ( *, '(a)' ) '    of a substring.'

      s = 'Bob is debobbing the bobber!'
      substring = 'bob'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  String =     "' // trim ( s ) // '"'
      write ( *, '(a)' ) '  Substring is "' // trim ( substring ) // '"'

      i1 = index ( s, substring )
      i2 = index ( s, trim ( substring ) )
      i3 = s_indexi ( s, substring )
      i4 = s_index_last ( s, substring )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  INDEX =              ', i1
      write ( *, '(a,i8)' ) '  INDEX (restricted) = ', i2
      write ( *, '(a,i8)' ) '  INDEXI =             ', i3
      write ( *, '(a,i8)' ) '  S_INDEX_LAST =       ', i4

      return
      end