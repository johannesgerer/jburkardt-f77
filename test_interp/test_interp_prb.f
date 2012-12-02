      program main

c*********************************************************************72
c
cc MAIN is the main program for TEST_INTERP_PRB.
c
c  Discussion:
c
c    TEST_INTERP_PRB calls the TEST_INTERP tests.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp (  )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_INTERP_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TEST_INTERP library.'
      write ( *, '(a)' ) '  This test also requires the R8LIB library.'

      call test01 ( )
      call test02 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_INTERP_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 shows how P00_STORY can be called.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer prob
      integer prob_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  P00_STORY prints the problem "story".'

      call p00_prob_num ( prob_num )

      do prob = 1, prob_num

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Problem ', prob

        call p00_story ( prob )

      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 prints the data for each problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    06 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer data_max
      parameter ( data_max = 100 )

      integer data_num
      integer dim_num
      real p(2,data_max)
      integer prob
      integer prob_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) 
     &  '  P00_DATA_NUM returns N, the number of data points.'
      write ( *, '(a)' ) 
     &  '  P00_DIM_NUM returns M, the dimension of data.'
      write ( *, '(a)' ) 
     &  '  P00_DATA returns the actual (MxN) data.'

      call p00_prob_num ( prob_num )

      do prob = 1, prob_num

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Problem  ', prob

        call p00_data_num ( prob, data_num )
        write ( *, '(a,i8)' ) '  DATA_NUM ', data_num
        call p00_dim_num ( prob, dim_num )
        write ( *, '(a,i8)' ) '  DIM_NUM  ', dim_num

        call p00_data ( prob, dim_num, data_num, p )

        call r8mat_transpose_print ( dim_num, data_num, p, 
     &    '  Data array:' )

      end do

      return
      end
