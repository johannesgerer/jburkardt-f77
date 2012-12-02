      program main

c*********************************************************************72
c
cc MAIN is the main program for MXV_OPENMP.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 May 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      include 'omp_lib.h'

      integer m_max
      parameter ( m_max = 10000 )
      integer n_max
      parameter ( n_max = 10000 )
      integer mn_max
      parameter ( mn_max = 1000000 )

      double precision a(mn_max)
      integer i
      integer m
      integer n
      integer proc_num
      integer thread_num
      double precision x(n_max)
      double precision y(m_max)

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MXV_OPENMP:'
      write ( *, '(a)' ) '  FORTRAN77/OpenMP version'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Compute matrix vector products y = A*x.'

      write ( *, '(a,i8)' ) 
     &  '  The number of processors = ', omp_get_num_procs ( )
      write ( *, '(a,i8)' ) 
     &  '  The number of threads    = ', omp_get_max_threads ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Compare various algorithms:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  MXV_PLAIN          - "plain vanilla" FORTRAN.'
      write ( *, '(a)' ) '  MXV_PLAIN_OPENMP  - PLAIN + OpenMP.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Algorithm                  M         N      Seconds'
c
c  N = M
c
      m = 10

      do i = 1, 3

        write ( *, '(a)' ) ' '

        n = m
        call test01 ( m, n, a, x, y )

        m = m * 10

      end do
c
c  N = 10 * M
c
      m = 1

      do i = 1, 3

        write ( *, '(a)' ) ' '

        n = 10 * m
        call test01 ( m, n, a, x, y )

        m = m * 10

      end do
c
c  M = 10 * N
c
      n = 1

      do i = 1, 3

        write ( *, '(a)' ) ' '

        m = 10 * n
        call test01 ( m, n, a, x, y )

        n = n * 10

      end do
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MXV_OPENMP:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( m, n, a, x, y )

c*********************************************************************72
c
cc TEST01 compares various algorithms for a given matrix size MxN.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns 
c    of the matrix.
c
c    Input, double precision A(M,N), X(N), Y(M), space to hold the
c    matrix and vectors.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      double precision omp_get_wtime
      double precision seconds
      double precision x(n)
      double precision y(m)

      call matgen ( m, n, a, x )

      seconds = omp_get_wtime ( )
      call mxv_plain ( m, n, a, x, y )
      seconds = omp_get_wtime ( ) - seconds
      write ( *, '(2x,a18,2x,i8,2x,i8,2x,g14.6)' ) 
     &  'MXV_PLAIN         ', m, n, seconds

      seconds = omp_get_wtime ( )
      call mxv_plain_openmp ( m, n, a, x, y )
      seconds = omp_get_wtime ( ) - seconds
      write ( *, '(2x,a18,2x,i8,2x,i8,2x,g14.6)' ) 
     &  'MXV_PLAIN_OPENMP ', m, n, seconds

      return
      end
      subroutine matgen ( m, n, a, x )

c*********************************************************************72
c
cc MATGEN generates a random matrix A and vector X.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns 
c    of the matrix.
c
c    Output, real ( kind = 8 ) A(M,N), the matrix.
c
c    Output, real ( kind = 8 ) X(N), the vector.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer i
      integer j
      integer seed
      double precision x(n)

      seed = 1325
c
c  Set the matrix A.
c
      do j = 1, n
        do i = 1, m
          seed = mod ( ( 3125 * seed ), 65536 )
          a(i,j) = ( seed - 32768.0 ) / 16384.0
        end do
      end do
c
c  Set X.
c
      do i = 1, n
        x(i) = i
      end do

      return
      end
      subroutine mxv_plain ( m, n, a, x, y )

c*********************************************************************72
c
cc MXV_PLAIN computes y = A * x, using "plain" code.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns
c    of the matrix.
c
c    Input, real ( kind = 8 ) A(M,N), the matrix.
c
c    Input, real ( kind = 8 ) X(N), the vector to be multiplied.
c
c    Output, real ( kind = 8 ) Y(M), the product vector.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer i
      integer j
      double precision x(n)
      double precision y(m)

      do i = 1, m
        y(i) = 0.0
        do j = 1, n
          y(i) = y(i) + a(i,j) * x(j)
        end do
      end do

      return
      end
      subroutine mxv_plain_openmp ( m, n, a, x, y )

c*********************************************************************72
c
cc MXV_PLAIN_OPENMP computes y = A * x, using OpenMP parallel directives.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Barbara Chapman, Gabriele Jost, Ruud vanderPas, David Kuck,
c    Using OpenMP: Portable Shared Memory Parallel Processing,
c    MIT Press, 2007,
c    ISBN13: 978-0262533027,
c    LC: QA76.642.C49.
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns 
c    of the matrix.
c
c    Input, real ( kind = 8 ) A(M,N), the matrix.
c
c    Input, real ( kind = 8 ) X(N), the vector to be multiplied.
c
c    Output, real ( kind = 8 ) Y(M), the product vector.
c
      implicit none

      integer m
      integer n

      double precision a(m,n)
      integer i
      integer j
      double precision x(n)
      double precision y(m)

c$omp parallel
c$omp& shared ( m, n, a, x, y ) 
c$omp& private ( i, j )

c$omp do
      do i = 1, m
        y(i) = 0.0
        do j = 1, n
          y(i) = y(i) + a(i,j) * x(j)
        end do
      end do
c$omp end do

c$omp end parallel

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    None
c
      implicit none

      character * ( 8 ) ampm
      integer d
      character * ( 8 ) date
      integer h
      integer m
      integer mm
      character * ( 9 ) month(12)
      integer n
      integer s
      character * ( 10 ) time
      integer y

      save month

      data month /
     &  'January  ', 'February ', 'March    ', 'April    ', 
     &  'May      ', 'June     ', 'July     ', 'August   ', 
     &  'September', 'October  ', 'November ', 'December ' /

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) y, m, d
      read ( time, '(i2,i2,i2,1x,i3)' ) h, n, s, mm

      if ( h .lt. 12 ) then
        ampm = 'AM'
      else if ( h .eq. 12 ) then
        if ( n .eq. 0 .and. s .eq. 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h .lt. 12 ) then
          ampm = 'PM'
        else if ( h .eq. 12 ) then
          if ( n .eq. 0 .and. s .eq. 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if

      write ( *, 
     &  '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) 
     &  d, month(m), y, h, ':', n, ':', s, '.', mm, ampm

      return
      end
