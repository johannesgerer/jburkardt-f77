      function wtime ( )

c*********************************************************************72
c
cc WTIME returns a reading of the wall clock time.
c
c  Discussion:
c
c    To get the elapsed wall clock time, call WTIME before and after a given
c    operation, and subtract the first reading from the second.
c
c    This function is meant to suggest the similar routines:
c
c      "omp_get_wtime ( )" in OpenMP,
c      "MPI_Wtime ( )" in MPI,
c      and "tic" and "toc" in MATLAB.
c
c    This function is only going to work if SYSTEM_CLOCK is available.
c    SYSTEM_CLOCK is not a FORTRAN77 standard library function.  It does
c    come with FORTRAN90, however.  So this function may work, but if so,
c    it is only by the grace of compiler writers!
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 April 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision WTIME, the wall clock reading, in seconds.
c
      implicit none

      integer clock_max
      integer clock_rate
      integer clock_reading
      double precision wtime

      call system_clock ( clock_reading, clock_rate, clock_max )

      wtime = dble ( clock_reading ) / dble ( clock_rate )

      return
      end
