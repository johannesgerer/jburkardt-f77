      program main

c*********************************************************************72
c
cc BINS_PRB tests routines from the BINS library.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 July 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BINS_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the BINS library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
c     call test05 ( )
c     call test07 ( )
c     call test08 ( )
      call test09 ( )

      call test10 ( )
c     call test11 ( )
c     call test12 ( )
c     call test13 ( )
c     call test14 ( )
c     call test15 ( )
c     call test16 ( )
c
c  Terminate
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'BINS_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests BIN_TO_R8_EVEN, R8_TO_BIN_EVEN.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 July 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      parameter ( a = 10.0D+00 )
      double precision b
      parameter ( b = 20.0D+00 )
      integer bin
      double precision c
      double precision cmax
      double precision cmin
      integer i
      integer nbin
      parameter ( nbin = 7 )
      double precision r8_uniform
      double precision rmax
      parameter ( rmax = 23.0D+00 )
      double precision rmin
      parameter ( rmin = 8.0D+00 )
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  R8_TO_BIN_EVEN puts a number into a bin.'
      write ( *, '(a)' ) '  BIN_TO_R8_EVEN returns the bin limits.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  The bins are equally spaced between A and B,'
      write ( *, '(a)' )
     &  '  with two extra bins, for things less than A,'
      write ( *, '(a)' ) '  or greater than B.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  A = ', a
      write ( *, '(a,g14.6)' ) '  B = ', b
      write ( *, '(a,i6)' ) '  Total number of bins = ', nbin

      call get_seed ( seed )
      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  Using random seed = ', seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  Generate some random values C and put them in bins.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '          C      Bin   Bin_Min       Bin_Max'
      write ( *, '(a)' ) ' '

      do i = 1, 30
        c = r8_uniform ( rmin, rmax, seed )
        call r8_to_bin_even ( nbin, a, b, c, bin )
        call bin_to_r8_even ( nbin, bin, a, b, cmin, cmax )
        write ( *, '(2x,g14.6,i4,2g14.6)' ) c, bin, cmin, cmax
      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests BIN_TO_R8_EVEN2, R8_TO_BIN_EVEN2.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 July 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      parameter ( a = 10.0D+00 )
      double precision b
      parameter ( b = 20.0D+00 )
      integer bin
      double precision c
      double precision cmax
      double precision cmin
      integer i
      integer nbin
      parameter ( nbin = 5 )
      double precision r8_uniform
      double precision rmax
      parameter ( rmax = 23.0D+00 )
      double precision rmin
      parameter ( rmin = 8.0D+00 )
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  BIN_TO_R8_EVEN2 returns the bin limits.'
      write ( *, '(a)' ) '  R8_TO_BIN_EVEN2 puts a number into a bin.'
      write ( *, '(a)' )
     &  '  The bins are equally spaced between A and B.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  A = ', a
      write ( *, '(a,g14.6)' ) '  B = ', b
      write ( *, '(a,i6)' ) '  Total number of bins = ', nbin

      call get_seed ( seed )
      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  Using random seed = ', seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  Generate some random values C and put them in bins.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C      Bin   Bin_Min  Bin_Max'
      write ( *, '(a)' ) ' '

      do i = 1, 30
        c = r8_uniform ( rmin, rmax, seed )
        call r8_to_bin_even2 ( nbin, a, b, c, bin )
        call bin_to_r8_even2 ( nbin, bin, a, b, cmin, cmax )
        write ( *, '(2x,g14.6,i4,2g14.6)' ) c, bin, cmin, cmax
      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests BIN_TO_R82_EVEN, R82_TO_BIN_EVEN.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    24 July 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a(2)
      double precision b(2)
      integer bin(2)
      double precision c(2)
      double precision cmax(2)
      double precision cmin(2)
      integer i
      integer nbin
      parameter ( nbin = 7 )
      double precision rmax(2)
      double precision rmin(2)
      integer seed

      save a
      save b
      save rmax
      save rmin

      data a /  5.0D+00,  0.0D+00 /
      data b / 15.0D+00, 20.0D+00 /
      data rmax / 23.0D+00, 21.0D+00 /
      data rmin / 3.0D+00, -2.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  BIN_TO_R82_EVEN returns the bin limits.'
      write ( *, '(a)' )
     &  '  R82_TO_BIN_EVEN puts a R82 number into a bin.'
      write ( *, '(a)' ) ''
      write ( *, '(a)' )
     &  '  The bins are equally spaced between A and B,'
      write ( *, '(a)' )
     &  '  with two extra bins, for things less than A,'
      write ( *, '(a)' ) '  or greater than B.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  A(1) = ', a(1)
      write ( *, '(a,g14.6)' ) '  B(1) = ', b(1)
      write ( *, '(a,g14.6)' ) '  A(2) = ', a(2)
      write ( *, '(a,g14.6)' ) '  B(2) = ', b(2)
      write ( *, '(a,i6)' ) '  Total number of bins = ', nbin
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  Generate some random values C and put them in bins.'
      write ( *, '(a)' )
     &  '  We list the X and Y components on separate lines.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       C      Bin   Bin_Min  Bin_Max'
      write ( *, '(a)' ) ' '

      call get_seed ( seed )

      do i = 1, 30
        call r82_uniform ( rmin, rmax, seed, c )
        call r82_to_bin_even ( nbin, a, b, c, bin )
        call bin_to_r82_even ( nbin, bin, a, b, cmin, cmax )
        write ( *, '(a)' ) ' '
        write ( *, '(2x,g14.6,i4,2g14.6)' )
     &    c(1), bin(1), cmin(1), cmax(1)
        write ( *, '(2x,g14.6,i4,2g14.6)' )
     &    c(2), bin(2), cmin(2), cmax(2)
      end do

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests R82VEC_BIN_EVEN, R82VEC_BINNED_REORDER, R82VEC_BINNED_SORT_A.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 July 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 30 )
      integer nbin
      parameter ( nbin = 4 )

      double precision a(2,n)
      double precision amax(2)
      double precision amin(2)
      double precision bin_max(2)
      double precision bin_min(2)
      integer bin_last(nbin,nbin)
      integer bin_next(n)
      integer bin_start(nbin,nbin)
      integer i1
      integer i2
      integer j
      integer k
      integer seed

      save amax
      save amin
      save bin_max
      save bin_min

      data amax / 23.0D+00, 12.0D+00 /
      data amin / 8.0D+00, 3.0D+00 /
      data bin_max / 20.0D+00, 10.0D+00 /
      data bin_min / 10.0D+00, 5.0D+00 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' )
     &  '  R82VEC_BIN_EVEN constructs evenly spaced bins and'
      write ( *, '(a)' )
     &  '    assigns each element of a R82VEC to a bin.'
      write ( *, '(a)' )
     &  '  R82VEC_BINNED_REORDER can reorder the array'
      write ( *, '(a)' ) '    to correspond to the bin ordering.'
      write ( *, '(a)' )
     &  '  R82VEC_BINNED_SORT_A can sort the individual bins'
      write ( *, '(a)' ) '    after the array has been reordered.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The bins are equally spaced between '
      write ( *, '(a)' ) '  BIN_MIN and BIN_MAX,'
      write ( *, '(a)' )
     &  '  with two extra bins, for things less than BIN_MIN,'
      write ( *, '(a)' ) '  or greater than BIN_MAX.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,2g14.6)' )
     &  '  Component 1 range: ', bin_min(1), bin_max(1)
      write ( *, '(a,2g14.6)' )
     &  '  Component 2 range: ', bin_min(2), bin_max(2)
      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' )
     &  '  Number of bins per row and column = ', nbin
      write ( *, '(a)' ) ' '

      call get_seed ( seed )

      call r82vec_uniform ( n, amin, amax, seed, a )

      call r82vec_print ( n, a, '  The data vector A to be binned:' )

      call r82vec_bin_even ( n, a, nbin, bin_min, bin_max, bin_start,
     &  bin_last, bin_next )

      call i4mat_print ( nbin, nbin, bin_start,
     &  '  The BIN_START array:' )

      call i4mat_print ( nbin, nbin, bin_start,
     &  '  The BIN_LAST array:' )

      call i4vec_print ( n, bin_next, '  The BIN_NEXT array:' )

      do i1 = 1, nbin

        do i2 = 1, nbin

          write ( *, '(a)' ) ' '
          write ( *, '(a,2i6)' ) '  Contents of bin number ', i1, i2
          write ( *, '(a)' ) ' '

          j = bin_start(i1,i2)
          k = 0

10        continue

          if ( 0 .lt. j ) then
            k = k + 1
            write ( *, '(2x,2i4,2g14.6)' ) k, j, a(1,j), a(2,j)
            j = bin_next(j)
            go to 10
          end if

        end do

      end do
c
c  Now reorder the data to correspond to the bins.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  Call R82VEC_BINNED_REORDER to reorder the array.'
      write ( *, '(a)' ) ' '

      call r82vec_binned_reorder ( n, a, nbin, bin_start, bin_last,
     &  bin_next )

      call r82vec_print ( n, a, '  The data vector, sorted by bins:' )

      call i4mat_print ( nbin, nbin, bin_start,
     &  '  The BIN_START array:' )

      call i4mat_print ( nbin, nbin, bin_last,
     &  '  The BIN_LAST array:' )

      call i4vec_print ( n, bin_next, '  The BIN_NEXT array:' )
c
c  Now sort the bins.
c
      call r82vec_binned_sort_a ( n, a, nbin, bin_start, bin_last )

      call r82vec_print ( n, a,
     &  '  The data vector, with sorted bins:' )

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests R82VEC_PART_QUICK_A.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    30 January 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 12 )

      double precision a(2,n)
      double precision alo(2)
      double precision ahi(2)
      integer l
      integer r
      integer seed

      save alo
      save ahi

      data alo / 0.0D+00, 2.0D+00 /
      data ahi / 10.0D+00, 3.0D+00 /

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  R82VEC_PART_QUICK_A reorders an R82VEC'
      write ( *, '(a)' ) '    as part of a quick sort.'
      write ( *, '(a,i12)' )
     &  '  Using initial random number seed = ', seed

      call r82vec_uniform ( n, alo, ahi, seed, a )

      call r82vec_print ( n, a, '  Before rearrangment:' )

      call r82vec_part_quick_a ( n, a, l, r )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Rearranged array'
      write ( *, '(a,i6)' ) '  Left index =  ', l
      write ( *, '(a,i6)' ) '  Key index =   ', l + 1
      write ( *, '(a,i6)' ) '  Right index = ', r
      write ( *, '(a)' ) ' '

      call r82vec_print ( l,     a(1:2,1:l),   '  Left half:' )
      call r82vec_print ( 1,     a(1:2,l+1),   '  Key:' )
      call r82vec_print ( n-l-1, a(1:2,l+2:n), '  Right half:' )

      return
      end
      subroutine test09 ( )

c*********************************************************************72
c
cc TEST09 tests R8VEC_BIN.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    25 July 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 25 )
      integer nbin
      parameter ( nbin = 5 )

      integer bin(0:nbin+1)
      double precision bin_limit(0:nbin)
      double precision bin_max
      double precision bin_min
      integer i
      integer seed
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST09'
      write ( *, '(a)' ) '  R8VEC_BIN computes bins for an R8VEC.'

      call get_seed ( seed )

      call r8vec_uniform ( n, -2.0D+00, 11.0D+00, seed, x )

      call r8vec_print ( n, x, '  The vector to be binned:' )

      bin_min =  0.0D+00
      bin_max = 10.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,i6)' ) '  Number of bins is ', nbin
      write ( *, '(a,g14.6)' ) '  Bin minimum is ', bin_min
      write ( *, '(a,g14.6)' ) '  Bin maximum is ', bin_max

      call r8vec_bin ( n, x, nbin, bin_min, bin_max, bin, bin_limit )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Lower Limit    Upper Limit    Count'
      write ( *, '(a)' ) ' '

      write ( *, '(2x,2f8.4,i4)' ) bin_min, bin_limit(0), bin(0)
      do i = 1, nbin
        write ( *, '(2x,2f8.4,i4)' )
     &    bin_limit(i-1), bin_limit(i), bin(i)
      end do
      write ( *, '(2x,f8.4,8x,i4)' ) bin_limit(nbin), bin(nbin+1)

      return
      end
      subroutine test10 ( )

c*********************************************************************72
c
cc TEST10 tests R8VEC_BIN_EVEN, R8VEC_BINNED_REORDER, R8VEC_BINNED_SORT_A.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    23 July 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 30 )
      integer nbin
      parameter ( nbin = 7 )

      double precision a(n)
      double precision amax
      parameter ( amax = 23.0D+00 )
      double precision amin
      parameter ( amin = 8.0D+00 )
      double precision bin_max
      parameter ( bin_max = 20.0D+00 )
      double precision bin_min
      parameter ( bin_min = 10.0D+00 )
      integer bin_last(nbin)
      integer bin_next(n)
      integer bin_start(nbin)
      integer i
      integer j
      integer k
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST10'
      write ( *, '(a)' )
     &  '  R8VEC_BIN_EVEN constructs evenly spaced bins and'
      write ( *, '(a)' ) '    assigns each element of a DVEC to a bin.'
      write ( *, '(a)' ) '  R8VEC_BINNED_REORDER can reorder the array'
      write ( *, '(a)' ) '    to correspond to the bin ordering.'
      write ( *, '(a)' ) '  R8VEC_BINNED_SORT_A can sort the array'
      write ( *, '(a)' ) '    once it has been reordered.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The bins are equally spaced between '
      write ( *, '(a)' ) '  BIN_MIN and BIN_MAX,'
      write ( *, '(a)' )
     &  '  with two extra bins, for things less than BIN_MIN,'
      write ( *, '(a)' ) '  or greater than BIN_MAX.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  BIN_MIN = ', bin_min
      write ( *, '(a,g14.6)' ) '  BIN_MAX = ', bin_max
      write ( *, '(a,i6)' ) '  Total number of bins = ', nbin
      write ( *, '(a)' ) ' '

      call get_seed ( seed )

      call r8vec_uniform ( n, amin, amax, seed, a )

      call r8vec_print ( n, a, '  The data vector A to be binned:' )

      call r8vec_bin_even ( n, a, nbin, bin_min, bin_max, bin_start,
     &  bin_last, bin_next )

      call i4vec_print ( nbin, bin_start, '  The BIN_START array:' )

      call i4vec_print ( nbin, bin_last, '  The BIN_LAST array:' )

      call i4vec_print ( n, bin_next, '  The BIN_NEXT array:' )

      do i = 1, nbin

        write ( *, '(a)' ) ' '
        write ( *, '(a,i6)' ) '  Contents of bin number ', i
        write ( *, '(a)' ) ' '

        j = bin_start(i)
        k = 0

        do while ( 0 < j )
          k = k + 1
          write ( *, '(2x,2i4,g14.6)' ) k, j, a(j)
          j = bin_next(j)
        end do

      end do
c
c  Now reorder the data to correspond to the bins.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  Call R8VEC_BINNED_REORDER to reorder the array.'
      write ( *, '(a)' ) ' '

      call r8vec_binned_reorder ( n, a, nbin, bin_start, bin_last,
     &  bin_next )

      call r8vec_print ( n, a, '  The data vector A:' )

      call i4vec_print ( nbin, bin_start, '  The BIN_START array:' )

      call i4vec_print ( nbin, bin_last, '  The BIN_LAST array:' )

      call i4vec_print ( n, bin_next, '  The BIN_NEXT array:' )
c
c  Now sort the data, one bin at a time
c
      call r8vec_binned_sort_a ( n, a, nbin, bin_start, bin_last )

      call r8vec_print ( n, a, '  The sorted data vector A:' )

      return
      end
