      program main

c*********************************************************************72
c
cc MAIN is the main program for TEST_APPROX_PRB.
c
c  Discussion:
c
c    TEST_APPROX_PRB calls the TEST_APPROX tests.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp (  )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_APPROX_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the TEST_APPROX library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
      call test06 ( )
      call test07 ( )
      call test08 ( )
      call test09 ( )

      call test10 ( )
      call test11 ( )
      call test12 ( )
      call test13 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST_APPROX_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 shows how P00_TITLE can be called.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer prob
      integer prob_num
      character * ( 80 ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) 
     &  '  Demonstrate some of the bookkeeping routines.'
      write ( *, '(a)' ) 
     &  '  P00_PROB_NUM returns the number of problems.'
      write ( *, '(a)' ) '  P00_TITLE returns the problem title.'
      write ( *, '(a)' ) '  P00_LIMIT returns the problem limits.'

      call p00_prob_num ( prob_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Number of problems = ', prob_num
      write ( *, '(a)' ) ' '

      do prob = 1, prob_num

        call p00_title ( prob, title )
        write ( *, '(2x,i2,2x,a)' ) prob, '"' // trim ( title ) // '".'

      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 shows how P00_STORY can be called.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer prob
      integer prob_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  P00_STORY prints the problem "story".'

      call p00_prob_num ( prob_num )

      do prob = 1, prob_num

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Problem ', prob

        call p00_story ( prob )

      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 uses polynomial interpolation on data vector problems.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer max_tab
      parameter ( max_tab = 12 )

      double precision diftab(max_tab)
      integer i
      integer j
      integer jhi
      character mark
      integer ntab
      integer data_num
      integer prob
      integer prob_num
      character * ( 80 ) title
      double precision x
      double precision xdata(max_tab)
      double precision yapprox
      double precision ydata(max_tab)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) 
     &  '  Polynomial interpolation to a vector of data.'

      call p00_prob_num ( prob_num )

      do prob = 1, prob_num

        call p00_title ( prob, title )

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Problem ', prob
        write ( *, '(2x,a)' ) trim ( title )

        call p00_data_num ( prob, data_num )

        write ( *, '(2x,a,i8)' ) '  DATA_NUM = ', data_num

        if ( max_tab .lt. data_num ) then

          write ( *, '(a)' ) ' '
          write ( *, '(a,i8)' ) '  Skipped problem ', prob
          write ( *, '(a)' ) '  Too big.'

        else

          call p00_dat ( prob, data_num, xdata, ydata )

          ntab = data_num

          write ( *, '(a)' ) ' '
          write ( *, '(a,i8)' ) 
     &      '  Interpolating polynomial order = ', ntab
          write ( *, '(a)' ) ' '
c
c  Construct the interpolating polynomial via finite differences.
c
          call data_to_dif ( ntab, xdata, ydata, diftab )
c
c  Print out the approximation, including midpoints of the intervals.
c
          do i = 1, ntab

            if ( i .lt. ntab ) then
              jhi = 2
            else
              jhi = 1
            end if

            do j = 1, jhi

              if ( i .lt. ntab ) then
                x = ( dble ( jhi - j + 1 ) * xdata(i)     
     &              + dble (       j - 1 ) * xdata(i+1) ) 
     &              / dble ( jhi         )
              else
                x = xdata(ntab)
              end if

              if ( j .eq. 1 ) then
                mark = '*'
              else
                mark = ' '
              end if

              call dif_val ( ntab, xdata, diftab, x, yapprox )

              write ( *, '(2x,a,2g14.6)' ) mark, x, yapprox

            end do

          end do

        end if

      end do

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 uses linear spline interpolation on all problems.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer max_data
      parameter ( max_data = 300 )

      double precision a
      double precision b
      integer i
      integer imax
      character mark
      integer data_num
      integer prob
      integer prob_num
      character * ( 80 ) title
      double precision xdata(max_data)
      double precision xval
      double precision ydata(max_data)
      double precision ypval
      double precision yval

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  Linear spline interpolation.'

      call p00_prob_num ( prob_num )

      do prob = 1, prob_num

        call p00_title ( prob, title )

        call p00_data_num ( prob, data_num )

        if ( max_data < data_num ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'TEST04 - Fatal error!'
          write ( *, '(a)' ) '  MAX_DATA < DATA_NUM'
          stop
        end if

        call p00_dat ( prob, data_num, xdata, ydata )

        a = xdata(1)
        b = xdata(data_num)

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Problem ', prob
        write ( *, '(2x,a)' ) trim ( title )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '       X          Y          Y'' '
        write ( *, '(a)' ) ' '
c
c  Evaluate the interpolation function.
c
        imax = 2 * data_num - 1

        do i = 1, imax

          xval = ( dble ( imax - i     ) * a   
     &           + dble (        i - 1 ) * b ) 
     &           / dble ( imax     - 1 )

          call spline_linear_val ( data_num, xdata, ydata, xval, yval, 
     &      ypval )

          if ( mod ( i, 2 ) .eq. 1 ) then
            mark = '*'
          else
            mark = ' '
          end if

          write ( *, '(2x,a,3g14.6)' ) mark, xval, yval, ypval

        end do

      end do

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 uses Overhauser spline interpolation on all problems.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer max_data
      parameter ( max_data = 300 )

      integer num_dim
      parameter ( num_dim = 1 )

      double precision a
      double precision b
      integer i
      integer j
      integer jhi
      integer jmax
      character mark
      integer data_num
      integer prob
      integer prob_num
      character * ( 80 ) title
      double precision xdata(max_data)
      double precision xval
      double precision ydata(max_data)
      double precision yval

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  Overhauser spline interpolation.'

      call p00_prob_num ( prob_num )

      do prob = 1, prob_num

        call p00_title ( prob, title )

        call p00_data_num ( prob, data_num )

        if ( max_data < data_num ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'TEST05 - Fatal error!'
          write ( *, '(a)' ) '  MAX_DATA < DATA_NUM'
          stop
        end if

        call p00_dat ( prob, data_num, xdata, ydata )

        a = xdata(1)
        b = xdata(data_num)

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Problem ', prob
        write ( *, '(2x,a)' ) trim ( title )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  X   Y'
        write ( *, '(a)' ) ' '
c
c  Evaluate the interpolation function.
c
        do i = 1, data_num - 1

          jmax = 3

          if ( i .eq. data_num - 1 ) then
            jhi = jmax
          else
            jhi = jmax - 1
          end if

          do j = 1, jhi

            xval = ( dble ( jmax - j      ) * xdata(i)     
     &             + dble (        j - 1 ) * xdata(i+1) ) 
     &             / dble ( jmax     - 1 )

            call spline_overhauser_val ( num_dim, data_num, xdata, 
     &        ydata, xval, yval )

            if ( j .eq. 1 .or. j .eq. 3 ) then
              mark = '*'
            else
              mark = ' '
            end if

            write ( *, '(2x,a,2g14.6)' ) mark, xval, yval

          end do

        end do

      end do

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 uses cubic spline interpolation on all problems.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer max_data
      parameter ( max_data = 300 )

      double precision a
      double precision b
      integer i
      integer ibcbeg
      integer ibcend
      integer j
      integer jhi
      integer jmax
      character mark
      integer data_num
      integer prob
      integer prob_num
      character * ( 80 ) title
      double precision xdata(max_data)
      double precision xval
      double precision ybcbeg
      double precision ybcend
      double precision ydata(max_data)
      double precision ypp(max_data)
      double precision yppval
      double precision ypval
      double precision yval

      ibcbeg = 0
      ibcend = 0
      ybcbeg = 0.0D+00
      ybcend = 0.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) '  Cubic spline interpolation.'

      call p00_prob_num ( prob_num )

      do prob = 1, prob_num

        call p00_title ( prob, title )

        call p00_data_num ( prob, data_num )

        if ( max_data < data_num ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'TEST06 - Fatal error!'
          write ( *, '(a)' ) '  MAX_DATA < DATA_NUM'
          stop
        end if

        call p00_dat ( prob, data_num, xdata, ydata )

        a = xdata(1)
        b = xdata(data_num)
c
c  Set up the interpolation function.
c
        call spline_cubic_set ( data_num, xdata, ydata, ibcbeg, 
     &    ybcbeg, ibcend, ybcend, ypp )

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Problem ', prob
        write ( *, '(2x,a)' ) trim ( title )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '    X   Y'
        write ( *, '(a)' ) ' '
c
c  Evaluate the interpolation function.
c
        do i = 1, data_num - 1

          jmax = 3

          if ( i .eq. data_num - 1 ) then
            jhi = jmax
          else
            jhi = jmax - 1
          end if

          do j = 1, jhi

            xval = ( dble ( jmax - j     ) * xdata(i)     
     &             + dble (        j - 1 ) * xdata(i+1) ) 
     &             / dble ( jmax     - 1 )

            call spline_cubic_val ( data_num, xdata, ydata, ypp, 
     &        xval, yval, ypval, yppval )

            if ( j .eq. 1 .or. j .eq. 3 ) then
              mark = '*'
            else
              mark = ' '
            end if

            write ( *, '(2x,a,2g14.6)' ) mark, xval, yval

          end do

        end do

      end do

      return
      end
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07 plots an Overhauser spline interpolant for problem 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer max_data
      parameter ( max_data = 300 )
      integer jmax
      parameter ( jmax = 7 )
      integer nplot
      parameter ( nplot = ( jmax - 1 ) * ( max_data - 1 ) + 1 )

      character * ( 80 ) approx_filename
      character * ( 80 ) data_filename
      integer data_num
      integer i
      integer j
      integer jhi
      integer num_dim
      parameter ( num_dim = 1 )
      integer plot
      integer prob
      double precision xdata(max_data)
      double precision xplot(nplot)
      double precision xval
      double precision ydata(max_data)
      double precision yplot(nplot)
      double precision yval

      approx_filename = 'test07_approx.txt'
      data_filename = 'test07_data.txt'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07'
      write ( *, '(a)' ) 
     &  '  Plot an Overhauser spline interpolant for problem 7.'
c
c  Get the problem data.
c
      prob = 7

      call p00_data_num ( prob, data_num )

      if ( max_data < data_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST07 - Fatal error!'
        write ( *, '(a)' ) '  MAX_DATA < DATA_NUM'
        stop
      end if

      call p00_dat ( prob, data_num, xdata, ydata )

      call r8vec2_write ( data_filename, data_num, xdata, ydata )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Data values stored in "' 
     &  // trim ( data_filename ) // '".'
c
c  Evaluate the approximating function.
c
      plot = 0

      do i = 1, data_num - 1

        if ( i .eq. data_num - 1 ) then
          jhi = jmax
        else
          jhi = jmax - 1
        end if

        do j = 1, jhi

          xval = ( dble ( jmax - j     ) * xdata(i)     
     &           + dble (        j - 1 ) * xdata(i+1) ) 
     &           / dble ( jmax     - 1 )

          call spline_overhauser_val ( num_dim, data_num, xdata, 
     &      ydata, xval, yval )

          plot = plot + 1
          xplot(plot) = xval
          yplot(plot) = yval

        end do

      end do

      call r8vec2_write ( approx_filename, nplot, xplot, yplot )

      write ( *, '(a)' ) '  Approximant values stored in "' 
     &  // trim ( approx_filename ) // '".'

      return
      end
      subroutine test08 ( )

c*********************************************************************72
c
cc TEST08 plots a cubic spline interpolant for problem 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer max_data
      parameter ( max_data = 300 )
      integer jmax
      parameter ( jmax = 7 )
      integer nplot
      parameter ( nplot = ( jmax - 1 ) * ( max_data - 1 ) + 1 )

      character * ( 80 ) approx_filename
      character * ( 80 ) data_filename
      integer i
      integer ibcbeg
      integer ibcend
      integer j
      integer jhi
      integer data_num
      integer plot
      integer prob
      double precision xdata(max_data)
      double precision xplot(nplot)
      double precision xval
      double precision ybcbeg
      double precision ybcend
      double precision ydata(max_data)
      double precision yplot(nplot)
      double precision ypp(max_data)
      double precision yppval
      double precision ypval
      double precision yval

      approx_filename = 'test08_approx.txt'
      data_filename = 'test08_data.txt'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08'
      write ( *, '(a)' ) 
     &  '  Plot a cubic spline interpolant for problem 7.'

      prob = 7
c
c  Get the data.
c
      call p00_data_num ( prob, data_num )

      if ( max_data < data_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST08 - Fatal error!'
        write ( *, '(a)' ) '  MAX_DATA < DATA_NUM'
        stop
      end if

      call p00_dat ( prob, data_num, xdata, ydata )

      call r8vec2_write ( data_filename, data_num, xdata, ydata )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Data values stored in "' 
     &  // trim ( data_filename ) // '".'
c
c  Set up the interpolation function.
c
      ibcbeg = 0
      ibcend = 0
      ybcbeg = 0.0D+00
      ybcend = 0.0D+00

      call spline_cubic_set ( data_num, xdata, ydata, ibcbeg, ybcbeg, 
     &  ibcend, ybcend, ypp )
c
c  Evaluate the interpolation function.
c
      plot = 0

      do i = 1, data_num - 1

        if ( i .eq. data_num - 1 ) then
          jhi = jmax
        else
          jhi = jmax - 1
        end if

        do j = 1, jhi

          xval = ( dble ( jmax - j     ) * xdata(i)     
     &           + dble (        j - 1 ) * xdata(i+1) ) 
     &           / dble ( jmax     - 1 )

          call spline_cubic_val ( data_num, xdata, ydata, ypp, xval, 
     &      yval, ypval, yppval )

          plot = plot + 1
          xplot(plot) = xval
          yplot(plot) = yval

        end do

      end do

      call r8vec2_write ( approx_filename, nplot, xplot, yplot )

      write ( *, '(a)' ) '  Approximant values stored in "' 
     &  // trim ( approx_filename ) // '".'

      return
      end
      subroutine test09 ( )

c*********************************************************************72
c
cc TEST09 uses B spline approximation on all problems.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer max_data
      parameter ( max_data = 300 )

      double precision a
      double precision b
      integer i
      integer j
      integer jhi
      integer jmax
      character mark
      integer data_num
      integer prob
      integer prob_num
      character * ( 80 ) title
      double precision xdata(max_data)
      double precision xval
      double precision ydata(max_data)
      double precision yval

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST09'
      write ( *, '(a)' ) '  B spline approximation.'

      call p00_prob_num ( prob_num )

      do prob = 1, prob_num

        call p00_title ( prob, title )

        call p00_data_num ( prob, data_num )

        if ( max_data < data_num ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'TEST09 - Fatal error!'
          write ( *, '(a)' ) '  MAX_DATA < DATA_NUM'
          stop
        end if

        call p00_dat ( prob, data_num, xdata, ydata )

        a = xdata(1)
        b = xdata(data_num)

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Problem ', prob
        write ( *, '(2x,a)' ) trim ( title )
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '       X        Y'
        write ( *, '(a)' ) ' '
c
c  Evaluate the interpolation function.
c
        do i = 1, data_num - 1

          jmax = 3

          if ( i .eq. data_num - 1 ) then
            jhi = jmax
          else
            jhi = jmax - 1
          end if

          do j = 1, jhi

            xval = ( dble ( jmax - j     ) * xdata(i)     
     &             + dble (        j - 1 ) * xdata(i+1) ) 
     &             / dble ( jmax     - 1 )

            call spline_b_val ( data_num, xdata, ydata, xval, yval )

            if ( j .eq. 1 .or. j .eq. 3 ) then
              mark = '*'
            else
              mark = ' '
            end if

            write ( *, '(2x,a,2g14.6)' ) mark, xval, yval

          end do

        end do

      end do

      return
      end
      subroutine test10 ( )

c*********************************************************************72
c
cc TEST10 plots a B spline approximant for problem 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer max_data
      parameter ( max_data = 300 )
      integer jmax
      parameter ( jmax = 7 )
      integer nplot
      parameter ( nplot = ( jmax - 1 ) * ( max_data - 1 ) + 1 )

      character * ( 80 ) approx_filename
      character * ( 80 ) data_filename
      integer i
      integer j
      integer jhi
      integer data_num
      integer plot
      integer prob
      character * ( 80 ) title
      double precision xdata(max_data)
      double precision xplot(nplot)
      double precision xval
      double precision ydata(max_data)
      double precision yplot(nplot)
      double precision yval

      approx_filename = 'test10_approx.txt'
      data_filename = 'test10_data.txt'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST10'
      write ( *, '(a)' ) '  Plot a B spline approximant for problem 7'

      prob = 7

      call p00_title ( prob, title )
c
c  Get the data.
c
      call p00_data_num ( prob, data_num )

      if ( max_data < data_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST10 - Fatal error!'
        write ( *, '(a)' ) '  MAX_DATA < DATA_NUM'
        stop
      end if

      call p00_dat ( prob, data_num, xdata, ydata )

      call r8vec2_write ( data_filename, data_num, xdata, ydata )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Data values stored in "' 
     &  // trim ( data_filename ) // '".'
c
c  Evaluate the approximation function.
c
      plot = 0

      do i = 1, data_num - 1

        if ( i .eq. data_num - 1 ) then
          jhi = jmax
        else
          jhi = jmax - 1
        end if

        do j = 1, jhi

          xval = ( dble ( jmax - j     ) * xdata(i)     
     &           + dble (        j - 1 ) * xdata(i+1) ) 
     &           / dble ( jmax     - 1 )

          call spline_b_val ( data_num, xdata, ydata, xval, yval )

          plot = plot + 1
          xplot(plot) = xval
          yplot(plot) = yval

        end do

      end do

      call r8vec2_write ( approx_filename, nplot, xplot, yplot )

      write ( *, '(a)' ) '  Approximant values stored in "' 
     &  // trim ( approx_filename ) // '".'

      return
      end
      subroutine test11 ( )

c*********************************************************************72
c
cc TEST11 plots a beta spline approximant for problem 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer max_data
      parameter ( max_data = 300 )
      integer jmax
      parameter ( jmax = 7 )
      integer nplot
      parameter ( nplot = ( jmax - 1 ) * ( max_data - 1 ) + 1 )

      character * ( 80 ) approx_filename
      double precision beta1
      double precision beta2
      character * ( 80 ) data_filename
      integer i
      integer j
      integer jhi
      integer data_num
      integer plot
      integer prob
      character * ( 80 ) title
      double precision xdata(max_data)
      double precision xplot(nplot)
      double precision xval
      double precision ydata(max_data)
      double precision yplot(nplot)
      double precision yval

      approx_filename = 'test11_approx.txt'
      data_filename = 'test11_data.txt'

      beta1 = 100.0D+00
      beta2 = 0.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST11'
      write ( *, '(a)' ) 
     &  '  Plot a beta spline approximant for problem 7'
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  BETA1 = ', beta1
      write ( *, '(a,g14.6)' ) '  BETA2 = ', beta2

      prob = 7

      call p00_title ( prob, title )
c
c  Get the data.
c
      call p00_data_num ( prob, data_num )

      if ( max_data < data_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST11 - Fatal error!'
        write ( *, '(a)' ) '  MAX_DATA < DATA_NUM'
        stop
      end if

      call p00_dat ( prob, data_num, xdata, ydata )

      call r8vec2_write ( data_filename, data_num, xdata, ydata )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Data values stored in "' 
     &  // trim ( data_filename ) // '".'
c
c  Evaluate the interpolation function.
c
      plot = 0

      do i = 1, data_num - 1

        if ( i .eq. data_num - 1 ) then
          jhi = jmax
        else
          jhi = jmax - 1
        end if

        do j = 1, jhi

          xval = ( dble ( jmax - j     ) * xdata(i)     
     &           + dble (        j - 1 ) * xdata(i+1) ) 
     &           / dble ( jmax     - 1 )

          call spline_beta_val ( beta1, beta2, data_num, xdata, 
     &      ydata, xval, yval )

          plot = plot + 1
          xplot(plot) = xval
          yplot(plot) = yval

        end do

      end do

      call r8vec2_write ( approx_filename, nplot, xplot, yplot )

      write ( *, '(a)' ) '  Approximant values stored in "' 
     &  // trim ( approx_filename ) // '".'

      return
      end
      subroutine test12 ( )

c*********************************************************************72
c
cc TEST12 plots a Bernstein spline approximant for problem 7.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer max_data
      parameter ( max_data = 300 )
      integer nplot
      parameter ( nplot = 101 )

      double precision a
      character * ( 80 ) approx_filename
      double precision b
      character * ( 80 ) data_filename
      integer i
      integer data_num
      integer plot
      integer prob
      double precision xdata(max_data)
      double precision xplot(nplot)
      double precision xval
      double precision ydata(max_data)
      double precision yplot(nplot)
      double precision yval

      approx_filename = 'test12_approx.txt'
      data_filename = 'test12_data.txt'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST12'
      write ( *, '(a)' ) '  Plot a Bernstein approximant for problem 5.'
      write ( *, '(a)' ) '  The Bernstein approximant requires'
      write ( *, '(a)' ) '  equally spaced data!'

      prob = 5
c
c  Get the data.
c
      call p00_data_num ( prob, data_num )

      if ( max_data < data_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST12 - Fatal error!'
        write ( *, '(a)' ) '  MAX_DATA < DATA_NUM'
        stop
      end if

      call p00_dat ( prob, data_num, xdata, ydata )

      call r8vec2_write ( data_filename, data_num, xdata, ydata )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Data values stored in "' 
     &  // trim ( data_filename ) // '".'
c
c  Evaluate the approximant function.
c
      a = xdata(1)
      b = xdata(data_num)

      do plot = 1, nplot

        xval = ( dble ( nplot - plot     ) * a     
     &         + dble (         plot - 1 ) * b ) 
     &         / dble ( nplot        - 1 )

        call bpab_approx ( data_num - 1, a, b, ydata, xval, yval )

        xplot(plot) = xval
        yplot(plot) = yval

      end do

      call r8vec2_write ( approx_filename, nplot, xplot, yplot )

      write ( *, '(a)' ) '  Approximant values stored in "' 
     &  // trim ( approx_filename ) // '".'

      return
      end
      subroutine test13 ( )

c*********************************************************************72
c
cc TEST13 plots a cubic spline interpolant for problem 5.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    07 February 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer max_data
      parameter ( max_data = 300 )

      integer nplot
      parameter ( nplot = 101 )

      character * ( 80 ) approx_filename 
      character * ( 80 ) data_filename
      integer ibcbeg
      integer ibcend
      integer j
      integer data_num
      integer prob
      double precision xdata(max_data)
      double precision xplot(nplot)
      double precision xval
      double precision ybcbeg
      double precision ybcend
      double precision ydata(max_data)
      double precision yplot(nplot)
      double precision ypp(max_data)
      double precision yppval
      double precision ypval
      double precision yval

      approx_filename = 'test13_approx.txt'
      data_filename = 'test13_data.txt'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST13'
      write ( *, '(a)' ) 
     &  '  Plot a cubic spline interpolant for problem 5'

      prob = 5

      call p00_data_num ( prob, data_num )

      if ( max_data < data_num ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST13 - Fatal error!'
        write ( *, '(a)' ) '  MAX_DATA < DATA_NUM'
        stop
      end if

      call p00_dat ( prob, data_num, xdata, ydata )

      call r8vec2_write ( data_filename, data_num, xdata, ydata )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Data values stored in "' 
     &  // trim ( data_filename ) // '".'
c
c  Set up the interpolation function.
c
      ibcbeg = 0
      ibcend = 0
      ybcbeg = 0.0D+00
      ybcend = 0.0D+00

      call spline_cubic_set ( data_num, xdata, ydata, ibcbeg, ybcbeg, 
     &  ibcend, ybcend, ypp )
c
c  Evaluate the interpolation function.
c
      do j = 1, nplot

        xval = ( dble ( nplot - j     ) * xdata(1)    
     &         + dble (         j - 1 ) * xdata(data_num) ) 
     &         / dble ( nplot     - 1 )

        call spline_cubic_val ( data_num, xdata, ydata, ypp, xval, 
     &    yval, ypval, yppval )

        xplot(j) = xval
        yplot(j) = yval

      end do

      call r8vec2_write ( approx_filename, nplot, xplot, yplot )

      write ( *, '(a)' ) '  Interpolant values stored in "' 
     &  // trim ( approx_filename ) // '".'

      return
      end
