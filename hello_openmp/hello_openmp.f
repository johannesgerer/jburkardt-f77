      program main

c*********************************************************************72
c
cc MAIN is the main program for HELLO.
c
c  Discussion:
c
c    HELLO is a "Hello, World" program for OpenMP.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 June 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      include 'omp_lib.h'

      integer id
      double precision wtime

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HELLO_OPENMP'
      write ( *, '(a)' ) '  FORTRAN77/OpenMP version'

      wtime = omp_get_wtime ( )

      write ( *, '(a,i8)' ) 
     &  '  The number of processors = ', omp_get_num_procs ( )
      write ( *, '(a,i8)' ) 
     &  '  The number of threads    = ', omp_get_max_threads ( )
c
c  OUTSIDE THE PARALLEL REGION, have each thread say hello (there's only 1).
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  OUTSIDE the parallel region.'
      write ( *, '(a)' ) ' '

      id = omp_get_thread_num ( )
      write ( *, '(a,i8,a,i8)' ) '  HELLO from process ', id

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Going INSIDE the parallel region:'
      write ( *, '(a)' ) ' '
c
c  INSIDE THE PARALLEL REGION, have each thread say hello.
c
c$omp parallel
c$omp& private ( id )
c
c  Have each thread say hello.
c
      id = omp_get_thread_num ( )
      write ( *, '(a,i8)' ) '  Hello from process ', id

c$omp end parallel

c
c  Finish up by measuring the elapsed time.
c
      wtime = omp_get_wtime ( ) - wtime

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Back OUTSIDE the parallel region.'
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'HELLO_OPENMP'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  Elapsed wall clock time = ', wtime

      stop
      end
