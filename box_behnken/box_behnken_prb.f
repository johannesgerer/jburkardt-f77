      program main

c*********************************************************************72
c
cc MAIN calls a set of problems for BOX_BEHNKEN.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 October 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BOX_BEHNKEN_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the BOX_BEHNKEN library.'

      call test01 ( )
      call test02 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BOX_BEHNKEN_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests BOX_BEHNKEN.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 October 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 3 )
      integer x_max
      parameter ( x_max = 50 )

      double precision range(dim_num,2)
      integer x_num
      double precision x(dim_num,x_max)

      data range /
     &  0.0D+00, 10.0D+00,  5.0D+00,
     &  1.0D+00, 11.0D+00, 15.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  BOX_BEHNKEN computes a Box-Behnken dataset.'

      call r8mat_transpose_print ( dim_num, 2, range,
     &  '  The ranges:' )

      call box_behnken_size ( dim_num, x_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6,a,i6)' ) '  For dimension DIM_NUM = ', dim_num,
     &  ' the Box-Behnken design is of size ', x_num

      if ( x_max .lt. x_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  There is not enough storage space to'
        write ( *, '(a)' ) '  generate the requested data.'
        return
      end if

      call box_behnken ( dim_num, x_num, range, x )

      call r8mat_transpose_print ( dim_num, x_num, x,
     &  '  The Box-Behnken design:' )

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests R8MAT_WRITE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer dim_num
      parameter ( dim_num = 4 )
      integer x_max
      parameter ( x_max = 50 )

      character * ( 80 ) file_out_name
      double precision range(dim_num,2)

      integer x_num
      double precision x(dim_num,x_max)

      data range/
     &  0.0D+00, 0.0D+00, 0.0D+00, 0.0D+00,
     &  1.0D+00, 1.0D+00, 1.0D+00, 1.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  R8MAT_WRITE writes a Box-Behnken'
      write ( *, '(a)' ) '  dataset to a file.'

      call r8mat_transpose_print ( dim_num, 2, range,
     &  '  The ranges:' )

      call box_behnken_size ( dim_num, x_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6,a,i6)' ) '  For dimension DIM_NUM = ', dim_num,
     &  ' the Box-Behnken design is of size ', x_num

      if ( x_max .lt. x_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  There is not enough storage space to'
        write ( *, '(a)' ) '  generate the requested data.'
        return
      end if

      call box_behnken ( dim_num, x_num, range, x )

      file_out_name = 'box_behnken_04_33.txt'

      call r8mat_write ( file_out_name, dim_num, x_num, x )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The data was written to the file "' //
     &  trim ( file_out_name ) // '".'

      return
      end

