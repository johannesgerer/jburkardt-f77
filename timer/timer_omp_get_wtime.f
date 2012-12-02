      program main

c*********************************************************************72
c
cc MAIN is the main program for TIMER_OMP_GET_WTIME.
c
c  Discussion:
c
c    TIMER_OMP_GET_WTIME uses OMP_GET_WTIME as the timer.
c
c    OMP_GET_WTIME is a timing utility accessible to C codes that
c    support OpenMP.  It returns the elapsed wallclock time in seconds.
c
c    Here, we run on as many threads as there are processors.  We could
c    force the number of threads to be 1 to make a better comparison to
c    timers that run on a single processor.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 July 2008
c
c  Author:
c
c    John Burkardt
c
      include 'omp_lib.h'

      integer proc_num
      integer thread_num

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TIMER_OMP_GET_WTIME'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Demonstrate the OMP_GET_WTIME timer.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  omp_get_wtime ( ) is an OpenMP function.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  It returns elapsed wall clock time in seconds.'
c
c  How many processors are available?
c
      proc_num = omp_get_num_procs ( )

      thread_num = proc_num

      call omp_set_num_threads ( thread_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  Number of processors is       ', proc_num
      write ( *, '(a,i8)' ) 
     &  '  Number of threads requested = ', thread_num

      call test03
      call test04
      call test05

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TIMER_OMP_GET_WTIME'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test03

c*********************************************************************72
c
cc TEST03 times some unvectorized loops.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 July 2008
c
c  Author:
c
c    John Burkardt
c
      include 'omp_lib.h'

      integer n_log_min
      integer n_log_max
      integer n_min
      integer n_max
      integer n_rep

      parameter ( n_log_min = 12 )
      parameter ( n_log_max = 19 )
      parameter ( n_min = 2**n_log_min )
      parameter ( n_max = 2**n_log_max )
      parameter ( n_rep = 5 )

      real delta(n_log_max,n_rep)
      real etime
      integer func
      integer i
      integer i_rep
      integer j
      integer n
      integer n_log
      real pi
      parameter ( pi = 3.141592653589793E+00 )
      integer seed
      double precision wtime1
      double precision wtime2
      real x(n_max)
      real y(n_max)

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  Time the following operations:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    y(1:n) =        x(1:n)  '
      write ( *, '(a)' ) '    y(1:n) = PI *   x(1:n) )'
      write ( *, '(a)' ) '    y(1:n) = sqrt ( x(1:n) )'
      write ( *, '(a)' ) '    y(1:n) = exp  ( x(1:n) )'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) 
     &  '  Data vectors will be of minimum size ', n_min
      write ( *, '(a,i12)' ) 
     &  '  Data vectors will be of maximum size ', n_max
      write ( *, '(a,i12)' ) 
     &  '  Number of repetitions of the operation: ', n_rep

      do func = 1, 4

        do i_rep = 1, n_rep
      
          do n_log = n_log_min, n_log_max

            n = 2**( n_log )

            call r4vec_uniform_01 ( n, seed, x )

            wtime1 = omp_get_wtime ( )

            if ( func == 1 ) then
              do i = 1, n
                y(i) = x(i)
              end do
            else if ( func == 2 ) then
              do i = 1, n
                y(i) = pi * x(i)
              end do
            else if ( func == 3 ) then
              do i = 1, n
                y(i) = sqrt ( x(i) )
              end do
            else if ( func == 4 ) then
              do i = 1, n
                y(i) = exp ( x(i) )
              end do
            end if

            wtime2 = omp_get_wtime ( )

            delta(n_log,i_rep) = wtime2 - wtime1

          end do

        end do

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) '  Timing results:'
        write ( *, '(a)' ) ' '
        write ( *, '(a,a)' ) 
     &    '    Vector Size  Rep #1        Rep #2        ', 
     &    'Rep #3        Rep #4        Rep #5'
        write ( *, '(a)' ) ' '

        do n_log = n_log_min, n_log_max
          n = 2**( n_log )
          write ( *, '(i10,5f14.6)' ) 
     &      n, ( delta(n_log,j), j = 1, n_rep )
        end do

      end do

      return
      end
      subroutine test04

c*********************************************************************72
c
cc TEST04 times the 2D nearest neighbor problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 July 2008
c
c  Author:
c
c    John Burkardt
c
      include 'omp_lib.h'

      integer n_log_min
      integer n_log_max
      integer n_min
      integer n_max
      integer n_rep
      integer n_test

      parameter ( n_log_min = 10 )
      parameter ( n_log_max = 19 )
      parameter ( n_min = 2**n_log_min )
      parameter ( n_max = 2**n_log_max )
      parameter ( n_rep = 5 )
      parameter ( n_test = 1 )

      real delta(n_log_max,n_rep)
      real dist_i
      real dist_min
      real etime
      integer i
      integer i_min
      integer i_rep
      integer j
      integer n
      integer n_log
      integer seed
      double precision wtime1
      double precision wtime2
      real x(2,n_max)
      real y(2)

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  Time the 2D nearest neighbor problem.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Given X(2,N) and Y(2),' 
      write ( *, '(a)' ) '    find X(2,*) closest to Y(2).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    do i = 1, n'
      write ( *, '(a)' ) 
     &  '      if distance ( x(2,i), y ) < minimum so far'
      write ( *, '(a)' ) '        x_min = x(2,i)'
      write ( *, '(a)' ) '    end do'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) 
     &  '  Data vectors will be of minimum size ', n_min
      write ( *, '(a,i12)' ) 
     &  '  Data vectors will be of maximum size ', n_max
      write ( *, '(a,i12)' ) 
     &  '  Number of repetitions of the operation: ', n_rep

      call r4vec_uniform_01 ( 2*n_max, seed, x )
      call r4vec_uniform_01 ( 2, seed, y )

      do i_rep = 1, n_rep

        do n_log = n_log_min, n_log_max

          n = 2**( n_log )

          wtime1 = omp_get_wtime ( )

          dist_min = 1000000.0E+00
          i_min = 0
          do i = 1, n
            dist_i = ( x(1,i) - y(1) )**2 + ( x(2,i) - y(2) )**2
            if ( dist_i < dist_min ) then
              dist_min = dist_i
              i_min = i
            end if
          end do

          wtime2 = omp_get_wtime ( )
          delta(n_log,i_rep) = wtime2 - wtime1

        end do

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Timing results:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) 
     &  '    Vector Size  Rep #1        Rep #2        Rep #3        ' ,
     &  'Rep #4        Rep #5'
      write ( *, '(a)' ) ' '
      do n_log = n_log_min, n_log_max
        n = 2**( n_log )
        write ( *, '(i10,5f14.6)' ) n, ( delta(n_log,j), j = 1, n_rep )
      end do

      return
      end
      subroutine test05

c*********************************************************************72
c
cc TEST05 times the matrix multiplication problem problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 July 2008
c
c  Author:
c
c    John Burkardt
c
      include 'omp_lib.h'

      integer l_log_min
      parameter ( l_log_min = 1 )
      integer l_log_max
      parameter ( l_log_max = 4 )
      integer l_min
      parameter ( l_min = 4**l_log_min )
      integer l_max
      parameter ( l_max = 4**l_log_max )
      integer rep_num
      parameter ( rep_num = 5 )

      double precision a(l_max,l_max)
      double precision b(l_max,l_max)
      double precision c(l_max,l_max)
      double precision delta(l_log_max,rep_num)
      real etime
      integer i
      integer j
      integer k
      integer l
      integer l_log
      integer m
      integer n
      integer rep
      integer seed
      double precision wtime1
      double precision wtime2

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  Time the matrix multiplication problem.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Compute C = A * B'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  where'
      write ( *, '(a)' ) '    A is an L by M matrix,'
      write ( *, '(a)' ) '    B is an M by N matrix,'
      write ( *, '(a)' ) '  and so'
      write ( *, '(a)' ) '    C is an L by N matrix.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  Minimum value of L = M = N = ', l_min
      write ( *, '(a,i12)' ) '  Maximum value of L = M = N = ', l_max
      write ( *, '(a,i12)' ) 
     &  '  Number of repetitions of the operation: ', rep_num

      seed = 123456789

      call r8mat_uniform_01 ( l_max, l_max, seed, a )
      call r8mat_uniform_01 ( l_max, l_max, seed, b )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Use nested DO loops for matrix multiplication.'

      do rep = 1, rep_num

        do l_log = l_log_min, l_log_max - 1

          l = 4**( l_log )
          m = l
          n = l

          wtime1 = omp_get_wtime ( )

          do i = 1, l
            do j = 1, l
              c(i,j) = 0.0D+00
              do k = 1, l
                c(i,j) = c(i,j) + a(i,k) * b(k,j)
              end do
            end do
          end do

          wtime2 = omp_get_wtime ( )

          delta(l_log,rep) = wtime2 - wtime1

        end do

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Timing results using nested DO loops:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' ) 
     &  '    Vector Size  Rep #1        Rep #2        Rep #3        ',
     &  'Rep #4        Rep #5'
      write ( *, '(a)' ) ' '
      do l_log = l_log_min, l_log_max - 1
        l = 4**( l_log )
        write ( *, '(i10,5f14.6)' ) l, ( delta(l_log,j), j = 1, rep_num)
      end do

      return
      end
      subroutine r4vec_uniform_01 ( n, seed, r )

c*********************************************************************72
c
cc R4VEC_UNIFORM_01 sets a real vector to unit pseudorandom numbers.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    05 March 2006
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Springer Verlag, pages 201-202, 1983.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, pages 362-376, 1986.
c
c    P A Lewis, A S Goodman, J M Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, pages 136-143, 1969.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, real R(N), the vector of pseudorandom values.
c
      implicit none

      integer n

      integer i
      integer k
      integer seed
      real r(n)

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed < 0 ) then
          seed = seed + 2147483647
        end if

        r(i) = real ( seed ) * 4.656612875E-10

      end do

      return
      end
      subroutine r8mat_uniform_01 ( m, n, seed, r )

c*********************************************************************72
c
cc R8MAT_UNIFORM_01 returns a unit pseudorandom R8MAT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    11 August 2004
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, Linus Schrage,
c    A Guide to Simulation,
c    Second Edition,
c    Springer, 1987,
c    ISBN: 0387964673,
c    LC: QA76.9.C65.B73.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, December 1986, pages 362-376.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley, 1998,
c    ISBN: 0471134031,
c    LC: T57.62.H37.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, Number 2, 1969, pages 136-143.
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns in the array.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R(M,N), the array of pseudorandom values.
c
      implicit none

      integer m
      integer n

      integer i
      integer i4_huge
      parameter ( i4_huge = 2147483647 )
      integer j
      integer k
      integer seed
      double precision r(m,n)

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8MAT_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      do j = 1, n

        do i = 1, m

          k = seed / 127773

          seed = 16807 * ( seed - k * 127773 ) - k * 2836

          if ( seed .lt. 0 ) then
            seed = seed + i4_huge
          end if

          r(i,j) = dble ( seed ) * 4.656612875D-10

        end do
      end do

      return
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Discussion:
c
c    This FORTRAN77 version is made available for cases where the
c    FORTRAN90 version cannot be used.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 September 2005
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

      character ( len = 8 ) date
      character ( len = 10 ) time

      call date_and_time ( date, time )

      write ( *, '(a8,2x,a10)' ) date, time

      return
      end
