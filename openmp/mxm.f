      program main

c********************************************************************72
c
cc MAIN is the main program for MXM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 February 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      include 'omp_lib.h'

      integer l_max
      parameter ( l_max = 500 )
      integer m_max
      parameter ( m_max = 500 )
      integer n_max
      parameter ( n_max = 500 )

      double precision a(l_max,n_max)
      double precision b(l_max,m_max)
      double precision c(m_max,n_max)
      integer l
      integer m
      integer n

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MXM'
      write ( *, '(a)' ) '  FORTRAN77/OpenMP version.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Matrix multiplication tests.'

      write ( *, '(a,i8)' ) 
     &  '  The number of processors = ', omp_get_num_procs ( )
      write ( *, '(a,i8)' ) 
     &  '  The number of threads    = ', omp_get_max_threads ( )

      l = 500
      m = 500
      n = 500

      call r8_mxm ( l, m, n, a, b, c )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MXM:'
      write ( *, '(a)' ) '  Normal end of execution.'

      return
      end
      subroutine r8_mxm ( l, m, n, a, b, c )

c********************************************************************72
c
c  Purpose:
c
c    R8_MXM carries out a matrix-matrix multiplication in R8 arithmetic.
c
c  Discussion:
c
c    A(LxN) = B(LxM) * C(MxN).
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    13 February 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer L, M, N, the dimensions that specify the sizes of the
c    A, B, and C matrices.
c
c    Workspace, double precision A(L,N), B(L,M), C(M,N), space for the
c    matrices.
c
      implicit none

      include 'omp_lib.h'

      double precision a(l,n)
      double precision b(l,m)
      double precision c(m,n)
      integer i
      integer j
      integer k
      integer l
      integer m
      integer n
      integer ops
      double precision r8_uniform_01
      double precision rate
      integer seed
      double precision time_begin
      double precision time_elapsed
      double precision time_stop
c
c  Assign values to the B and C matrices.
c
      seed = 123456789

      do j = 1, m
        do i = 1, l
          b(i,j) = r8_uniform_01 ( seed )
        end do
      end do

      do j = 1, n
        do i = 1, m
          c(i,j) = r8_uniform_01 ( seed )
        end do
      end do
c
c  Compute A = B * C.
c
      time_begin = omp_get_wtime ( )

c$omp parallel
c$omp& shared ( a, b, c, l, m, n )
c$omp& private ( i, j, k )

c$omp do

      do j = 1, n
        do i = 1, l
          a(i,j) = 0.0D+00
          do k = 1, m
            a(i,j) = a(i,j) + b(i,k) * c(k,j)
          end do
        end do
      end do
c$omp end do

c$omp end parallel

      time_stop = omp_get_wtime ( )
c
c  Report.
c
      ops = l * n * ( 2 * m )
      time_elapsed = time_stop - time_begin
      rate = dble ( ops ) / time_elapsed / 1000000.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'R8_MXM matrix multiplication timing.'
      write ( *, '(a)' ) '  A(LxN) = B(LxM) * C(MxN).'
      write ( *, '(a,i8)' ) '  L = ', l
      write ( *, '(a,i8)' ) '  M = ', m
      write ( *, '(a,i8)' ) '  N = ', n
      write ( *, '(a,i12)' ) '  Floating point OPS roughly ', ops
      write ( *, '(a,g14.6)' ) '  Elapsed time dT = ', time_elapsed
      write ( *, '(a,g14.6)' ) '  Rate = MegaOPS/dT = ', rate

      return
      end
      function r8_uniform_01 ( seed )

c*********************************************************************72
c
cc R8_UNIFORM_01 returns a unit pseudorandom R8.
c
c  Discussion:
c
c    This routine implements the recursion
c
c      seed = 16807 * seed mod ( 2**31 - 1 )
c      r8_uniform_01 = seed / ( 2**31 - 1 )
c
c    The integer arithmetic never requires more than 32 bits,
c    including a sign bit.
c
c    If the initial seed is 12345, then the first three computations are
c
c      Input     Output      R8_UNIFORM_01
c      SEED      SEED
c
c         12345   207482415  0.096616
c     207482415  1790989824  0.833995
c    1790989824  2035175616  0.947702
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
c    Springer Verlag, pages 201-202, 1983.
c
c    Pierre L'Ecuyer,
c    Random Number Generation,
c    in Handbook of Simulation,
c    edited by Jerry Banks,
c    Wiley Interscience, page 95, 1998.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, pages 362-376, 1986.
c
c    Peter Lewis, Allen Goodman, James Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, pages 136-143, 1969.
c
c  Parameters:
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, double precision R8_UNIFORM_01, a new pseudorandom variate,
c    strictly between 0 and 1.
c
      implicit none

      double precision r8_uniform_01
      integer k
      integer seed

      if ( seed .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'R8_UNIFORM_01 - Fatal error!'
        write ( *, '(a)' ) '  Input value of SEED = 0.'
        stop
      end if

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed .lt. 0 ) then
        seed = seed + 2147483647
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r8_uniform_01 = dble ( seed ) * 4.656612875D-10

      return
      end
