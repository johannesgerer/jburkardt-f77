      program main

c*********************************************************************72
c
cc MAIN is the main program for TEST_INTERP_2D_PRB.
c
c  Discussion:
c
c    TEST_INTERP_2D_PRB tests the TEST_INTERP_2D library.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 October 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_INTERP_2D_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TEST_INTERP_2D library.'
      write ( *, '(a)' ) 
     &  '  The test requires access to the R8LIB library.'

      call test01 ( )
      call test02 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_INTERP_2D_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 simply prints the title of each grid and function.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer f_num
      integer fi
      character ( len = 50 ) ft
      integer g_num
      integer gi
      character ( len = 50 ) gt

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) 
     &  '  For each grid and function, print the title.'

      call g00_num ( g_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  GRIDS:'
      write ( *, '(a)' ) '  Index  Title'
      write ( *, '(a)' ) ' '

      do gi = 1, g_num

        call g00_title ( gi, gt )

        write ( *, '(2x,i2,2x,a)' ) gi, gt

      end do

      call f00_num ( f_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  FUNCTIONS:'
      write ( *, '(a)' ) '  Index  Title'
      write ( *, '(a)' ) ' '

      do fi = 1, f_num

        call f00_title ( fi, ft )

        write ( *, '(2x,i2,2x,a)' ) fi, ft

      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 samples each function using each grid.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 January 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision f(100)
      double precision f_ave
      double precision f_max
      double precision f_min
      integer f_num
      integer fi
      character * ( 50 ) ft
      integer g_num
      integer gi
      integer gn
      character * ( 50 ) gt
      double precision gx(100)
      double precision gy(100)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Sample each function over each grid.'

      call g00_num ( g_num )
      call f00_num ( f_num )

      do fi = 1, f_num

        call f00_title ( fi, ft )
        write ( *, '(a)' ) ' '
        write ( *, '(2x,i2,2x,a)' ) fi, trim ( ft )
        write ( *, '(a)' ) '        Grid Title                     ' //
     &    'Min(F)          Ave(F)           Max(F)'
        write ( *, '(a)' ) ' '

        do gi = 1, g_num

          call g00_title ( gi, gt )
          call g00_size ( gi, gn )

          call g00_xy ( gi, gn, gx, gy )

          call f00_f0 ( fi, gn, gx, gy, f )

          call r8vec_max ( gn, f, f_max )
          call r8vec_min ( gn, f, f_min )
          call r8vec_sum ( gn, f, f_ave )
          f_ave = f_ave / dble ( gn )

          write ( *, '(2x,i4,2x,a25,2x,g14.6,2x,g14.6,2x,g14.6)' ) 
     &      gi, gt, f_min, f_ave, f_max

        end do

      end do

      return
      end
