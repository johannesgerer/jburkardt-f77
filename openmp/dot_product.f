      program main

c*********************************************************************72
c
cc MAIN is the main program for DOT_PRODUCT.
c
c  Discussion:
c
c    This program illustrates how a vector dot product could be set up
c    in a FORTRAN77 program using OpenMP.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 April 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      include 'omp_lib.h'
c
c  I wanted to have N_MAX = 1,000,000, but that seems to be too
c  big for my PowerPC G5 system.
c
      integer n_max
      parameter ( n_max = 100000 )

      double precision factor
      integer i
      integer n
      double precision wtime
      double precision x(n_max)
      double precision xdoty 
      double precision y(n_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DOT_PRODUCT'
      write ( *, '(a)' ) '  FORTRAN77/OpenMP version'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  A program to compute a vector dot product.'

      write ( *, '(a,i8)' ) 
     &  '  The number of processors = ', omp_get_num_procs ( )
      write ( *, '(a,i8)' ) 
     &  '  The number of threads    = ', omp_get_max_threads ( )
c
c  Set up the vector data.
c  N may be increased to get better timing data.
c
c  The value FACTOR is chosen so that the correct value of the dot product 
c  of X and Y is N.
c
      n = 100

10    continue

      if ( n < n_max ) then

        n = n * 10

        factor = dble ( n )
        factor = 1.0D+00 
     &    / sqrt ( 2.0D+00 * factor * factor + 3 * factor + 1.0D+00 )

        do i = 1, n
          x(i) = i * factor
        end do

        do i = 1, n
          y(i) = i * 6 * factor
        end do

        write ( *, '(a)' ) ' '
c
c  Test #1
c
        wtime = omp_get_wtime ( )

        call test01 ( n, x, y, xdoty )

        wtime = omp_get_wtime ( ) - wtime

        write ( * , '(2x,a10,2x,i8,2x,g14.6,2x,f15.10)' ) 
     &  'Sequential', n, xdoty, wtime
c
c  Test #2
c
        wtime = omp_get_wtime ( )

        call test02 ( n, x, y, xdoty )

        wtime = omp_get_wtime ( ) - wtime

        write ( * , '(2x,a10,2x,i8,2x,g14.6,2x,f15.10)' ) 
     &  'Parallel  ', n, xdoty, wtime
      
        go to 10

      end if
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DOT_PRODUCT'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
      subroutine test01 ( n, x, y, xdoty )

c*********************************************************************72
c
cc  TEST01 computes the dot product with no parallel processing directives.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the vectors.
c
c    Input, real X(N), Y(N), the vectors.
c
c    Output, real XDOTY, the dot product of X and Y.
c
      implicit none

      integer n

      integer i
      double precision xdoty
      double precision x(n)
      double precision y(n)

      xdoty = 0.0D+00

      do i = 1, n
        xdoty = xdoty + x(i) * y(i)
      end do

      return
      end
      subroutine test02 ( n, x, y, xdoty )

c*********************************************************************72
c
cc  TEST02 computes the dot product with parallel processing directives.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the order of the vectors.
c
c    Input, real X(N), Y(N), the vectors.
c
c    Output, real XDOTY, the dot product of X and Y.
c
      implicit none

      integer n

      integer i
      double precision xdoty
      double precision x(n)
      double precision y(n)

      xdoty = 0.0D+00

c$omp parallel 
c$omp& shared ( n, x, y )
c$omp& private ( i )

c$omp do reduction ( + : xdoty )
      do i = 1, n
        xdoty = xdoty + x(i) * y(i)
      end do
c$omp end do

c$omp end parallel

      return
      end
