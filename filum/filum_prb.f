      program main

c*********************************************************************72
c
cc MAIN is the main program for FILUM_PRB.
c
c  Discussion:
c
c    FILUM_PRB tests routines from the FILUM library.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FILUM_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the FILUM library.'
     
      call test03 ( )
      call test06 ( )

      call test14 ( )
      call test145 ( )
      call test15 ( )
      call test16 ( )
      call test165 ( )

      call test22 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FILUM_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests FILE_COLUMN_COUNT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 June 2007
c
c  Author:
c
c    John Burkardt 
c
      implicit none

      integer column_num
      character * ( 100 ) file_name

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) 
     &  '  FILE_COLUMN_COUNT counts the columns in a file.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  It is assumed that the file contains a number of lines,'
      write ( *, '(a)' ) 
     &  '  with each line containing the same number of words.'
      write ( *, '(a)' ) 
     &  '  The task is to determine the number of words in a line,'
      write ( *, '(a)' ) '  that is, the number of "columns" of text.'

      file_name = 'filum_prb_4by5.txt'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Examining the file: ' 
      write ( *, '(a)' ) '    "' // trim ( file_name ) // '".'

      call file_column_count ( file_name, column_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of columns: ', column_num

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 tests FILE_EXIST.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 November 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      logical file_exist
      character * ( 255 ) file_name

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) '  FILE_EXIST reports whether a file "exists".'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Exist?  File_name'
      write ( *, '(a)' ) ' '

      file_name = 'filum_prb.f'
      write ( *, '(7x,l1,2x,a)' ) 
     &  file_exist ( file_name ), trim ( file_name )
      file_name = 'filum.f'
      write ( *, '(7x,l1,2x,a)' ) 
     &  file_exist ( file_name ), trim ( file_name )
      file_name = 'raisin.txt'
      write ( *, '(7x,l1,2x,a)' ) 
     &  file_exist ( file_name ), trim ( file_name )
      file_name = 'make.money.fast'
      write ( *, '(7x,l1,2x,a)' ) 
     &  file_exist ( file_name ), trim ( file_name )

      return
      end
      subroutine test14 ( )

c*********************************************************************72
c
cc TEST14 tests FILENAME_INC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    14 September 2005
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer test_num
      parameter ( test_num = 4 )

      character * ( 20 ) filename
      character * ( 20 ) filename_old
      integer j
      character * ( 20 ) string(test_num)
      integer test

      string(1) = 'file???.dat'
      string(2) = 'file072.dat'
      string(3) = '2cat9.dat'
      string(4) = 'fred99.txt'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST14'
      write ( *, '(a)' ) '  FILENAME_INC increments a string'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     Input             Output'

      do test = 1, test_num
        filename = string(test)
        write ( *, '(a)' ) ' '
        do j = 1, 4
          filename_old = filename
          call filename_inc ( filename )
          write ( *, '(2x,a,2x,a)' ) filename_old, filename
          if ( len_trim ( filename ) .le. 0 ) then
            write ( *, '(a)' ) '  (Empty output string.  Quit loop!)'
            exit
          end if
        end do
      end do

      return
      end
      subroutine test145 ( )

c*********************************************************************72
c
cc TEST145 tests FILENAME_INC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      character * ( 20 ) filename
      character * ( 20 ) filename_old
      integer test

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST145'
      write ( *, '(a)' ) '  FILENAME_INC increments a string.'
      write ( *, '(a)' ) 
     &  '  This test checks that a file name is properly'
      write ( *, '(a)' ) '  incremented when carrying is required.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     Input             Output'

      filename = 'file???.dat'
      write ( *, '(a)' ) ' '
      do test = 1, 10
        filename_old = filename
        call filename_inc ( filename )
        write ( *, '(2x,a,2x,a)' ) filename_old, filename
        if ( len_trim ( filename ) .le. 0 ) then
          write ( *, '(a)' ) 
     &      '  (File name not incrementable.  Quit loop!)'
          go to 10
        end if
      end do

10    continue

      filename = 'file072.dat'
      write ( *, '(a)' ) ' '
      do test = 1, 10
        filename_old = filename
        call filename_inc ( filename )
        write ( *, '(2x,a,2x,a)' ) filename_old, filename
        if ( len_trim ( filename ) .le. 0 ) then
          write ( *, '(a)' ) 
     &      '  (File name not incrementable.  Quit loop!)'
          go to 20
        end if
      end do

20    continue

      filename = '2cat9.dat'
      write ( *, '(a)' ) ' '
      do test = 1, 10
        filename_old = filename
        call filename_inc ( filename )
        write ( *, '(2x,a,2x,a)' ) filename_old, filename
        if ( len_trim ( filename ) .le. 0 ) then
          write ( *, '(a)' ) 
     &      '  (File name not incrementable.  Quit loop!)'
          go to 30
        end if
      end do

30    continue

      filename = 'fred98.txt'
      write ( *, '(a)' ) ' '
      do test = 1, 10
        filename_old = filename
        call filename_inc ( filename )
        write ( *, '(2x,a,2x,a)' ) filename_old, filename
        if ( len_trim ( filename ) .le. 0 ) then
          write ( *, '(a)' ) 
     &      '  (File name not incrementable.  Quit loop!)'
          go to 40
        end if
      end do

40    continue

      return
      end
      subroutine test15 ( )

c*********************************************************************72
c
cc TEST15 tests FILENAME_INC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 June 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer ios
      character ( len = 20 ) s

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST15'
      write ( *, '(a)' ) '  FILENAME_INC "increments" a string, such '
      write ( *, '(a)' ) '  as a file name.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  In this example, the string is a file name'
      write ( *, '(a)' ) '  of the form'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     "file_00.txt".'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  We know we have a sequence of files named'
      write ( *, '(a)' ) '    file_001.txt, file_002.txt, ...'
      write ( *, '(a)' ) '  and we want to generate the name of the'
      write ( *, '(a)' ) '  next file and open it.  If it doesn''t '
      write ( *, '(a)' ) '  exist, exit.'

      write ( *, '(a)' ) ' '

      s = 'file_00.txt'

10    continue

        call filename_inc ( s )

        write ( *, '(a)' ) '  Looking for file "' // trim ( s ) // '".'

        open ( unit = 1, file = s, status = 'old', err = 20 )

        close ( unit = 1 )

      go to 10

20    continue

      return
      end
      subroutine test16 ( )

c*********************************************************************72
c
cc TEST16 tests FILENAME_INC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 June 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer sim
      character * ( 20 ) s
      character * ( 20 ) s1
      character * ( 20 ) s2
      integer time_step

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST16'
      write ( *, '(a)' ) '  FILENAME_INC "increments" a string, such'
      write ( *, '(a)' ) '  as a file name.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  In this example, the string is a file name'
      write ( *, '(a)' ) '  of the form'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     "file_s00_t000.txt'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  The user is going to carry out several simulations.'
      write ( *, '(a)' ) 
     &  '  For each simulation, a number of time steps are done.'
      write ( *, '(a)' ) 
     &  '  In the file name, the "s" file records the simulation,'
      write ( *, '(a)' ) '  and the "t" field records the time step.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  So a typical file name is "file_s05_t017.txt".'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Assuming we have 5 simulations, and 4 time steps,'
      write ( *, '(a)' ) 
     &  '  the following double loop will generate all the'
      write ( *, '(a)' ) '  file names, from'
      write ( *, '(a)' ) '    "file_s01_t001.txt"'
      write ( *, '(a)' ) '  to'
      write ( *, '(a)' ) '    "file_s05_t004.txt".'

      s1 = 'file_s00'
      s2 = '_t000.txt'

      do sim = 1, 5

        call filename_inc ( s1 )

        s = trim ( s1 ) // trim ( s2 )

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Begin simulation ', sim
        write ( *, '(a)' ) ' '

        do time_step = 1, 4

          call filename_inc ( s )

          write ( *, '(2x,a)' ) trim ( s )

        end do

      end do

      return
      end
      subroutine test165 ( )

c*********************************************************************72
c
cc TEST165 tests FILENAME_INC_NOWRAP.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 November 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      character * ( 20 ) filename
      character * ( 20 ) filename_old
      integer test

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST165'
      write ( *, '(a)' ) '  FILENAME_INC_NOWRAP increments a string.'
      write ( *, '(a)' ) 
     &  '  This test checks that a file name is properly'
      write ( *, '(a)' ) '  incremented when carrying is required.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     Input             Output'

      filename = 'file???.dat'
      write ( *, '(a)' ) ' '
      do test = 1, 10
        filename_old = filename
        call filename_inc_nowrap ( filename )
        write ( *, '(2x,a,2x,a)' ) filename_old, filename
        if ( len_trim ( filename ) .le. 0 ) then
          write ( *, '(a)' ) 
     &      '  (File name not incrementable.  Quit loop!)'
          go to 10
        end if
      end do

10    continue

      filename = 'file072.dat'
      write ( *, '(a)' ) ' '
      do test = 1, 10
        filename_old = filename
        call filename_inc_nowrap ( filename )
        write ( *, '(2x,a,2x,a)' ) filename_old, filename
        if ( len_trim ( filename ) .le. 0 ) then
          write ( *, '(a)' ) 
     &      '  (File name not incrementable.  Quit loop!)'
          go to 20
        end if
      end do

20    continue

      filename = '2cat9.dat'
      write ( *, '(a)' ) ' '
      do test = 1, 10
        filename_old = filename
        call filename_inc_nowrap ( filename )
        write ( *, '(2x,a,2x,a)' ) filename_old, filename
        if ( len_trim ( filename ) .le. 0 ) then
          write ( *, '(a)' ) 
     &      '  (File name not incrementable.  Quit loop!)'
          go to 30
        end if
      end do

30    continue

      filename = 'fred98.txt'
      write ( *, '(a)' ) ' '
      do test = 1, 10
        filename_old = filename
        call filename_inc_nowrap ( filename )
        write ( *, '(2x,a,2x,a)' ) filename_old, filename
        if ( len_trim ( filename ) .le. 0 ) then
          write ( *, '(a)' ) 
     &      '  (File name not incrementable.  Quit loop!)'
          go to 40
        end if
      end do

40    continue

      return
      end
      subroutine test22 ( )

c*********************************************************************72
c
cc TEST22 tests FILE_ROW_COUNT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 June 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      character * ( 255 ) file_name
      integer line_num

      file_name = 'filum_prb_test.txt'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST22'
      write ( *, '(a)' ) '  FILE_ROW_COUNT counts the "rows" in a file.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Examining file:'
      write ( *, '(a)' ) '    "' // trim ( file_name ) // '".'
      write ( *, '(a)' ) ' '

      call file_row_count ( file_name, line_num )
      write ( *, '(a,i8)' ) '  Number of lines:      ', line_num

      return
      end
