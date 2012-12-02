      program main

c*********************************************************************72
c
cc FFTW3_PRB demonstrates the use of FFTW3.
c
c  Modified:
c
c    05 November 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FFTW3_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the FFTW3 library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FFTW3_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'
 
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01: complex 1D data.
c
c  Modified:
c
c    26 October 2005
c
      implicit none

      integer n
      parameter ( n = 100 )

      include "fftw3.f"

      double precision a
      double precision b
      double precision r8_uniform_01
      integer i
      double complex in(n)
      double complex in2(n)
      double complex out(n)
      integer*8 plan_backward
      integer*8 plan_forward
      integer seed

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  Demonstrate FFTW3 on a single vector '
      write ( *, '(a)' ) '  of complex data.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Transform data to FFT coefficients.'
      write ( *, '(a)' ) '  Backtransform FFT coefficients to recover'
      write ( *, '(a)' ) '  the data.'
      write ( *, '(a)' ) '  Compare recovered data to original data.'
c
c  Compute the data.
c
      do i = 1, n
        a = r8_uniform_01 ( seed )
        b = r8_uniform_01 ( seed )
        in(i) = dcmplx ( a, b )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Input Data:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,i4,2x,2g14.6)' ) i, in(i)
      end do
c
c  Make a plan for the FFT, and forward transform the data.
c
      call dfftw_plan_dft_1d_ ( plan_forward, n, in, out, 
     &  FFTW_FORWARD, FFTW_ESTIMATE )

      call dfftw_execute_ ( plan_forward )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Output FFT Coefficients:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,i4,2x,2g14.6)' ) i, out(i)
      end do
c
c  Make a plan for the backward FFT, and recover the original data.
c
      call dfftw_plan_dft_1d_ ( plan_backward, n, out, in2, 
     &  FFTW_BACKWARD, FFTW_ESTIMATE )

      call dfftw_execute_ ( plan_backward )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Recovered input data divided by N:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,i4,2x,2g14.6)' ) i, in2(i) / dble ( n )
      end do
c
c  Discard the information associated with the plans.
c
      call dfftw_destroy_plan_ ( plan_forward )
      call dfftw_destroy_plan_ ( plan_backward )

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02: real 1D data.
c
c  Modified:
c
c    23 October 2005
c
      implicit none

      integer n
      integer nc

      parameter ( n = 100 )
      parameter ( nc = 51 )

      include "fftw3.f"

      double precision r8_uniform_01
      integer i
      double precision in(n)
      double precision in2(n)
      double complex out(nc)
      integer*8 plan_backward
      integer*8 plan_forward
      integer seed

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  Demonstrate FFTW3 on a single vector'
      write ( *, '(a)' ) '  of real data.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Transform data to FFT coefficients.'
      write ( *, '(a)' ) '  Backtransform FFT coefficients to recover '
      write ( *, '(a)' ) '  data.'
      write ( *, '(a)' ) '  Compare recovered data to original data.'
c
c  Set up the input data, a real vector of length N.
c
      do i = 1, n
        in(i) = r8_uniform_01 ( seed )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Input Data:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,i4,2x,g14.6)' ) i, in(i)
      end do
c
c  Set up a plan, and execute the plan to transform the IN data to
c  the OUT FFT coefficients.
c
      call dfftw_plan_dft_r2c_1d_ ( plan_forward, n, in, out, 
     &  FFTW_ESTIMATE )

      call dfftw_execute_ ( plan_forward )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Output FFT Coefficients:'
      write ( *, '(a)' ) ' '

      do i = 1, nc
        write ( *, '(2x,i4,2x,2g14.6)' ) i, out(i)
      end do
c
c  Set up a plan, and execute the plan to backtransform the
c  complex FFT coefficients in OUT to real data.
c
      call dfftw_plan_dft_c2r_1d_ ( plan_backward, n, out, in2,
     &  FFTW_ESTIMATE )

      call dfftw_execute_ ( plan_backward )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Recovered input data divide by N:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,i4,2x,g14.6)' ) i, in2(i) / dble ( n )
      end do
c
c  Release the memory associated with the plans.
c
      call dfftw_destroy_plan_ ( plan_forward )
      call dfftw_destroy_plan_ ( plan_backward )

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03: complex 2D data.
c
c  Modified:
c
c    05 November 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nx
      parameter ( nx = 8 )
      integer ny
      parameter ( ny = 10 )

      include "fftw3.f"

      double precision a
      double precision b
      double precision r8_uniform_01
      integer i
      double complex in(nx,ny)
      double complex in2(nx,ny)
      integer j
      double complex out(nx,ny)
      integer*8 plan_backward
      integer*8 plan_forward
      integer seed

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  Demonstrate FFTW3 on a 2D complex array'
      write ( *, '(a,i8)' ) '  NX = ', nx
      write ( *, '(a,i8)' ) '  NY = ', ny
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Transform data to FFT coefficients.'
      write ( *, '(a)' ) '  Backtransform FFT coefficients to recover'
      write ( *, '(a)' ) '  the data.'
      write ( *, '(a)' ) '  Compare recovered data to original data.'
c
c  Compute the data.
c
      do j = 1, ny
        do i = 1, nx
          a = r8_uniform_01 ( seed )
          b = r8_uniform_01 ( seed )
          in(i,j) = dcmplx ( a, b )
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Input Data:'
      write ( *, '(a)' ) ' '

      do i = 1, nx
        do j = 1, ny
          write ( *, '(2x,i4,2x,i4,2x,2g14.6)' ) i, j, in(i,j)
        end do
      end do
c
c  Make a plan for the FFT, and forward transform the data.
c
      call dfftw_plan_dft_2d_ ( plan_forward, nx, ny, in, out, 
     &  FFTW_FORWARD, FFTW_ESTIMATE )

      call dfftw_execute_ ( plan_forward )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Output FFT Coefficients:'
      write ( *, '(a)' ) ' '

      do i = 1, nx
        do j = 1, ny
          write ( *, '(2x,i4,2x,i4,2x,2g14.6)' ) i, j, out(i,j)
        end do
      end do
c
c  Make a plan for the backward FFT, and recover the original data.
c
      call dfftw_plan_dft_2d_ ( plan_backward, nx, ny, out, in2, 
     &  FFTW_BACKWARD, FFTW_ESTIMATE )

      call dfftw_execute_ ( plan_backward )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Recovered input data divided by NX * NY:'
      write ( *, '(a)' ) ' '

      do i = 1, nx
        do j = 1, ny
          write ( *, '(2x,i4,2x,i4,2x,2g14.6)' )
     &      i, j, in2(i,j) / dble ( nx * ny )
        end do
      end do
c
c  Discard the information associated with the plans.
c
      call dfftw_destroy_plan_ ( plan_forward )
      call dfftw_destroy_plan_ ( plan_backward )

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04: real 2D data.
c
c  Discussion:
c
c    In contrast to the C example, in FORTRAN it is the FIRST dimension
c    of the complex coefficient array that is "half" the size.
c
c  Modified:
c
c    05 November 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nx
      parameter ( nx = 8 )
      integer ny
      parameter ( ny = 10 )

      integer nh
      parameter ( nh = ( nx / 2 ) + 1 )

      include "fftw3.f"

      double precision a
      double precision b
      double precision r8_uniform_01
      integer i
      double precision in(nx,ny)
      double precision in2(nx,ny)
      integer j
      double complex out(nh,ny)
      integer*8 plan_backward
      integer*8 plan_forward
      integer seed

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  Demonstrate FFTW3 on a 2D real array'
      write ( *, '(a,i8)' ) '  NX = ', nx
      write ( *, '(a,i8)' ) '  NY = ', ny
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Transform data to FFT coefficients.'
      write ( *, '(a)' ) '  Backtransform FFT coefficients to recover'
      write ( *, '(a)' ) '  the data.'
      write ( *, '(a)' ) '  Compare recovered data to original data.'
c
c  Compute the data.
c
      do j = 1, ny
        do i = 1, nx
          in(i,j) = r8_uniform_01 ( seed )
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Input Data:'
      write ( *, '(a)' ) ' '

      do i = 1, nx
        do j = 1, ny
          write ( *, '(2x,i4,2x,i4,2x,g14.6)' ) i, j, in(i,j)
        end do
      end do
c
c  Make a plan for the FFT, and forward transform the data.
c
      call dfftw_plan_dft_r2c_2d_ ( plan_forward, nx, ny, in, out, 
     &  FFTW_ESTIMATE )

      call dfftw_execute_ ( plan_forward )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Output FFT Coefficients:'
      write ( *, '(a)' ) ' '

      do i = 1, nh
        do j = 1, ny
          write ( *, '(2x,i4,2x,i4,2x,2g14.6)' ) i, j, out(i,j)
        end do
      end do
c
c  Make a plan for the backward FFT, and recover the original data.
c
      call dfftw_plan_dft_c2r_2d_ ( plan_backward, nx, ny, out, in2, 
     &  FFTW_ESTIMATE )

      call dfftw_execute_ ( plan_backward )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Recovered input data divided by NX * NY:'
      write ( *, '(a)' ) ' '

      do i = 1, nx
        do j = 1, ny
          write ( *, '(2x,i4,2x,i4,2x,g14.6)' )
     &      i, j, in2(i,j) / dble ( nx * ny )
        end do
      end do
c
c  Discard the information associated with the plans.
c
      call dfftw_destroy_plan_ ( plan_forward )
      call dfftw_destroy_plan_ ( plan_backward )

      return
      end
      function r8_uniform_01 ( seed )

c*********************************************************************72
c
cc R8_UNIFORM_01 returns a unit double precision pseudorandom number.
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
c  Modified:
c
c    23 October 2005
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, L E Schrage,
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
c    P A Lewis, A S Goodman, J M Miller,
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

      k = seed / 127773

      seed = 16807 * ( seed - k * 127773 ) - k * 2836

      if ( seed < 0 ) then
        seed = seed + 2147483647
      end if
c
c  Although SEED can be represented exactly as a 32 bit integer,
c  it generally cannot be represented exactly as a 32 bit real number!
c
      r8_uniform_01 = dble ( seed ) * 4.656612875D-10

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
