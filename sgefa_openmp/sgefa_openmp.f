      program main

c*********************************************************************72
c
cc MAIN is the main program for the SGEFA_OPENMP test program.
c
c  Discussion:
c
c    We want to compare methods of solving the linear system A*x=b.
c
c    The first way uses the standard sequential algorithm "SGEFA".
c
c    We also compare variant codes, "SGEFA_C" and "SGEFA_R", with
c    and without the use of OpenMP.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 April 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      include 'omp_lib.h'

      integer n_max
      parameter ( n_max = 1000 )

      real a(n_max,n_max)
      real b(n_max)
      integer ipvt(n_max)
      integer n
      integer proc_num
      integer thread_num
      real x(n_max)

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SGEFA_OPENMP'
      write ( *, '(a)' ) '  FORTRAN77/OpenMP version'

      write ( *, '(a,i8)' ) 
     &  '  The number of processors = ', omp_get_num_procs ( )
      write ( *, '(a,i8)' ) 
     &  '  The number of threads    = ', omp_get_max_threads ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  ' Algorithm        Mode          N    Error       Time'

      write ( *, '(a)' ) ' '
      n = 10
      call test_sgefa_c ( n, a, b, ipvt, x )
      call test_sgefa_c_omp ( n, a, b, ipvt, x )
      call test_sgefa_r ( n, a, b, ipvt, x )
      call test_sgefa_r_omp ( n, a, b, ipvt, x )

      write ( *, '(a)' ) ' '
      n = 100
      call test_sgefa_c ( n, a, b, ipvt, x )
      call test_sgefa_c_omp ( n, a, b, ipvt, x )
      call test_sgefa_r ( n, a, b, ipvt, x )
      call test_sgefa_r_omp ( n, a, b, ipvt, x )

      write ( *, '(a)' ) ' '
      n = 1000
      call test_sgefa_c ( n, a, b, ipvt, x )
      call test_sgefa_c_omp ( n, a, b, ipvt, x )
      call test_sgefa_r ( n, a, b, ipvt, x )
      call test_sgefa_r_omp ( n, a, b, ipvt, x )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SGEFA_OPENMP'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      return
      end
      subroutine test_sgefa_c ( n, a, b, ipvt, x )

c*********************************************************************72
c
cc TEST_SGEFA_C tests SGEFA_C.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 April 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      include 'omp_lib.h'

      real a(n,n)
      real b(n)
      real err
      integer i
      integer info
      integer ipvt(n)
      integer job
      integer lda
      integer n
      double precision wtime
      real x(n)
c
c  Generate the linear system A * x = b.
c
      lda = n

      call matgen ( lda, n, a, x, b )
c
c  Factor the linear system.
c
      wtime = omp_get_wtime ( )
      call sgefa_c ( a, lda, n, ipvt, info )
      wtime = omp_get_wtime ( ) - wtime

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST_SGEFA_C - Fatal error.'
        write ( *, '(a)' ) '  SGEFA_C reports singularity.'
        return
      end if
c
c  Solve the linear system.
c
      job = 0
      call sgesl ( a, lda, n, ipvt, b, job )

      err = 0.0
      do i = 1, n
        err = err + abs ( x(i) - b(i) )
      end do

      write ( *, '(a,i8,2x,g10.4,2x,f10.4)' ) 
     &  '  SGEFA_C   Sequential   ', n, err, wtime

      return
      end
      subroutine test_sgefa_c_omp ( n, a, b, ipvt, x )

c*********************************************************************72
c
cc TEST_SGEFA_C_OMP tests SGEFA_C_OMP.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 April 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      include 'omp_lib.h'

      real a(n,n)
      real b(n)
      real err
      integer i
      integer info
      integer ipvt(n)
      integer job
      integer lda
      integer n
      double precision wtime
      real x(n)
c
c  Generate the linear system A * x = b.
c
      lda = n

      call matgen ( lda, n, a, x, b )
c
c  Factor the linear system.
c
      wtime = omp_get_wtime ( )
      call sgefa_c_omp ( a, lda, n, ipvt, info )
      wtime = omp_get_wtime ( ) - wtime

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST_SGEFA_C_OMP - Fatal error.'
        write ( *, '(a)' ) '  SGEFA_C_OMP reports singularity.'
        return
      end if
c
c  Solve the linear system.
c
      job = 0
      call sgesl ( a, lda, n, ipvt, b, job )

      err = 0.0
      do i = 1, n
        err = err + abs ( x(i) - b(i) )
      end do

      write ( *, '(a,i8,2x,g10.4,2x,f10.4)' ) 
     &  '  SGEFA_C     Parallel   ', n, err, wtime

      return
      end
      subroutine test_sgefa_r ( n, a, b, ipvt, x )

c*********************************************************************72
c
cc TEST_SGEFA_R tests SGEFA_R.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 April 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      include 'omp_lib.h'

      real a(n,n)
      real b(n)
      real err
      integer i
      integer info
      integer ipvt(n)
      integer job
      integer lda
      integer n
      double precision wtime
      real x(n)
c
c  Generate the linear system A * x = b.
c
      lda = n

      call matgen ( lda, n, a, x, b )
c
c  Factor the linear system.
c
      wtime = omp_get_wtime ( )
      call sgefa_r ( a, lda, n, ipvt, info )
      wtime = omp_get_wtime ( ) - wtime

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST_SGEFA_R - Fatal error.'
        write ( *, '(a)' ) '  SGEFA_R reports singularity.'
        return
      end if
c
c  Solve the linear system.
c
      job = 0
      call sgesl ( a, lda, n, ipvt, b, job )

      err = 0.0
      do i = 1, n
        err = err + abs ( x(i) - b(i) )
      end do

      write ( *, '(a,i8,2x,g10.4,2x,f10.4)' ) 
     &  '  SGEFA_R   Sequential   ', n, err, wtime

      return
      end
      subroutine test_sgefa_r_omp ( n, a, b, ipvt, x )

c*********************************************************************72
c
cc TEST_SGEFA_R_OMP tests SGEFA_R_OMP.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 April 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      include 'omp_lib.h'

      real a(n,n)
      real b(n)
      real err
      integer i
      integer info
      integer ipvt(n)
      integer job
      integer lda
      integer n
      double precision wtime
      real x(n)
c
c  Generate the linear system A * x = b.
c
      lda = n

      call matgen ( lda, n, a, x, b )
c
c  Factor the linear system.
c
      wtime = omp_get_wtime ( )
      call sgefa_r_omp ( a, lda, n, ipvt, info )
      wtime = omp_get_wtime ( ) - wtime

      if ( info .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST_SGEFA_R_OMP - Fatal error.'
        write ( *, '(a)' ) '  SGEFA_R_OMP reports singularity.'
        return
      end if
c
c  Solve the linear system.
c
      job = 0
      call sgesl ( a, lda, n, ipvt, b, job )

      err = 0.0
      do i = 1, n
        err = err + abs ( x(i) - b(i) )
      end do

      write ( *, '(a,i8,2x,g10.4,2x,f10.4)' ) 
     &  '  SGEFA_R     Parallel   ', n, err, wtime

      return
      end
      subroutine matgen ( lda, n, a, x, b )

c*********************************************************************72
c 
cc MATGEN generates a "random" matrix for testing.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer LDA, the leading dimension of the matrix.
c
c    Input, integer N, the order of the matrix, and the length of the vector.
c
c    Output, real A(LDA,N), the matrix.
c
c    Output, real X(N), the solution vector.
c
c    Output, real B(N), the right hand side vector.
c
      implicit none

      integer lda
      integer n

      real a(lda,n)
      real b(n)
      integer i
      integer j
      integer seed
      real value
      real x(n)

      seed = 1325
c
c  Set the matrix A.
c
      do j = 1, n
        do i = 1, n
          seed = mod ( 3125 * seed, 65536 )
          value = ( real ( seed ) - 32768.0 ) / 16384.0
          a(i,j) = value
        end do
      end do
c
c  Set x.
c
      do i = 1, n
        x(i) = real ( i ) / real ( n )
      end do
c
c  Set b = A * x.
c
      do i = 1, n
        b(i) = 0.0
        do j = 1, n
          b(i) = b(i) + a(i,j) * x(j)
        end do
      end do

      return
      end
      subroutine sgefa_c ( a, lda, n, ipvt, info )

c*********************************************************************72
c 
cc SGEFA_C factors a matrix by gaussian elimination.
c
c  Discussion:
c
c    The step in which multiples of the pivot row are added to individual
c    rows has been replaced by a single call which updates the entire
c    matrix sub-block.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c  Parameters:
c
c    Input/output, real A(LDA,N).  On input, the matrix to be factored.
c    On output, an upper triangular matrix and the multipliers which were 
c    used to obtain it.  The factorization can be written A = L * U where
c    L is a product of permutation and unit lower triangular matrices and
c    U is upper triangular.
c
c    Input, integer LDA, the leading dimension of the matrix.
c
c    Input, integer N, the order of the matrix.
c
c    Output, integer IPVT(N), the pivot indices.
c
c    Output, integer INFO, indicates singularity.
c    If 0, this is the normal value, and the algorithm succeeded.
c    If K, then on the K-th elimination step, a zero pivot was encountered.
c    The matrix is numerically not invertible.
c
      implicit none

      integer lda
      integer n

      real a(lda,n)
      integer i
      integer info
      integer ipvt(n)
      integer j
      integer k
      integer l
      real t

      info = 0
      do k = 1, n - 1
c
c  Find the pivot index L.
c
        l = k
        t = abs ( a(k,k) )
        do i = k + 1, n
          if ( t .lt. abs ( a(i,k) ) ) then
            l = i
            t = abs ( a(i,k) )
          end if
        end do
        ipvt(k) = l

        if ( a(l,k) .eq. 0.0E+00 ) then
          info = k
          return
        end if
c
c  Interchange rows K and L.
c
        if ( l .ne. k ) then
          do j = k, n
            t      = a(l,j)
            a(l,j) = a(k,j)
            a(k,j) = t
          end do
        end if
c
c  Compute column K of the lower triangular factor.
c
        do i = k + 1, n
          a(i,k) = - a(i,k) / a(k,k)
        end do
c
c  Add multiples of the pivot row to the remaining rows.
c
        do j = k + 1, n
          do i = k + 1, n
            a(i,j) = a(i,j) + a(i,k) * a(k,j)
          end do
        end do

      end do

      ipvt(n) = n
      if ( a(n,n) .eq. 0.0E+00 ) then
        info = n
      end if

      return
      end
      subroutine sgefa_c_omp ( a, lda, n, ipvt, info )

c*********************************************************************72
c 
cc SGEFA_C_OMP factors a matrix by gaussian elimination.
c
c  Discussion:
c
c    This variant of SGEFA_C uses OpenMP for parallel execution.
c
c    The step in which multiples of the pivot row are added to individual
c    rows has been replaced by a single call which updates the entire
c    matrix sub-block.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c  Parameters:
c
c    Input/output, real A(LDA,N).  On input, the matrix to be factored.
c    On output, an upper triangular matrix and the multipliers which were 
c    used to obtain it.  The factorization can be written A = L * U where
c    L is a product of permutation and unit lower triangular matrices and
c    U is upper triangular.
c
c    Input, integer LDA, the leading dimension of the matrix.
c
c    Input, integer N, the order of the matrix.
c
c    Output, integer IPVT(N), the pivot indices.
c
c    Output, integer INFO, indicates singularity.
c    If 0, this is the normal value, and the algorithm succeeded.
c    If K, then on the K-th elimination step, a zero pivot was encountered.
c    The matrix is numerically not invertible.
c
      implicit none

      integer lda
      integer n

      real a(lda,n)
      integer i
      integer info
      integer ipvt(n)
      integer j
      integer k
      integer l
      real t
      info = 0
      do k = 1, n - 1
c
c  Find the pivot index L.
c
        l = k
        t = abs ( a(k,k) )
        do i = k + 1, n
          if ( t .lt. abs ( a(i,k) ) ) then
            l = i
            t = abs ( a(i,k) )
          end if
        end do
        ipvt(k) = l

        if ( a(l,k) .eq. 0.0E+00 ) then
          info = k
          return
        end if
c
c  Interchange rows K and L.
c
        if ( l .ne. k ) then
          do j = k, n
            t      = a(l,j)
            a(l,j) = a(k,j)
            a(k,j) = t
          end do
        end if
c
c  Compute column K of the lower triangular factor.
c  Add multiples of the pivot row to the remaining rows.
c
c$omp parallel 
c$omp& shared ( a, k, n ) 
c$omp& private ( i, j )

c$omp do
        do i = k + 1, n
          a(i,k) = - a(i,k) / a(k,k)
          do j = k + 1, n
            a(i,j) = a(i,j) + a(i,k) * a(k,j)
          end do
        end do
c$omp end do
c$omp end parallel

      end do

      ipvt(n) = n
      if ( a(n,n) .eq. 0.0E+00 ) then
        info = n
      end if

      return
      end
      subroutine sgefa_r ( a, lda, n, ipvt, info )

c*********************************************************************72
c 
cc SGEFA_R factors a matrix by gaussian elimination.
c
c  Discussion:
c
c    The step in which multiples of the pivot row are added to individual
c    rows has been replaced by a single call which updates the entire
c    matrix sub-block.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c  Parameters:
c
c    Input/output, real A(LDA,N).  On input, the matrix to be factored.
c    On output, an upper triangular matrix and the multipliers which were 
c    used to obtain it.  The factorization can be written A = L * U where
c    L is a product of permutation and unit lower triangular matrices and
c    U is upper triangular.
c
c    Input, integer LDA, the leading dimension of the matrix.
c
c    Input, integer N, the order of the matrix.
c
c    Output, integer IPVT(N), the pivot indices.
c
c    Output, integer INFO, indicates singularity.
c    If 0, this is the normal value, and the algorithm succeeded.
c    If K, then on the K-th elimination step, a zero pivot was encountered.
c    The matrix is numerically not invertible.
c
      implicit none

      integer lda
      integer n

      real a(lda,n)
      integer i
      integer info
      integer ipvt(n)
      integer j
      integer k
      integer l
      real t

      info = 0
      do k = 1, n - 1
c
c  Find the pivot index L.
c
        l = k
        t = abs ( a(k,k) )
        do i = k + 1, n
          if ( t .lt. abs ( a(i,k) ) ) then
            l = i
            t = abs ( a(i,k) )
          end if
        end do
        ipvt(k) = l

        if ( a(l,k) .eq. 0.0E+00 ) then
          info = k
          return
        end if
c
c  Interchange rows K and L.
c
        if ( l .ne. k ) then
          do j = k, n
            t      = a(l,j)
            a(l,j) = a(k,j)
            a(k,j) = t
          end do
        end if
c
c  Compute column K of the lower triangular factor.
c
        do i = k + 1, n
          a(i,k) = - a(i,k) / a(k,k)
        end do
c
c  Add multiples of the pivot row to the remaining rows.
c
        do i = k + 1, n
          do j = k + 1, n
            a(i,j) = a(i,j) + a(i,k) * a(k,j)
          end do
        end do

      end do

      ipvt(n) = n
      if ( a(n,n) .eq. 0.0E+00 ) then
        info = n
      end if

      return
      end
      subroutine sgefa_r_omp ( a, lda, n, ipvt, info )

c*********************************************************************72
c 
cc SGEFA_R_OMP factors a matrix by gaussian elimination.
c
c  Discussion:
c
c    This variant of SGEFA_R uses OpenMP for parallel execution.
c
c    The step in which multiples of the pivot row are added to individual
c    rows has been replaced by a single call which updates the entire
c    matrix sub-block.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    27 April 2008
c
c  Author:
c
c    C original version by Wesley Petersen.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
c    LINPACK User Guide,
c    SIAM, 1979,
c    ISBN13: 978-0-898711-72-1,
c    LC: QA214.L56.
c
c  Parameters:
c
c    Input/output, real A(LDA,N).  On input, the matrix to be factored.
c    On output, an upper triangular matrix and the multipliers which were 
c    used to obtain it.  The factorization can be written A = L * U where
c    L is a product of permutation and unit lower triangular matrices and
c    U is upper triangular.
c
c    Input, integer LDA, the leading dimension of the matrix.
c
c    Input, integer N, the order of the matrix.
c
c    Output, integer IPVT(N), the pivot indices.
c
c    Output, integer INFO, indicates singularity.
c    If 0, this is the normal value, and the algorithm succeeded.
c    If K, then on the K-th elimination step, a zero pivot was encountered.
c    The matrix is numerically not invertible.
c
      implicit none

      integer lda
      integer n

      real a(lda,n)
      integer i
      integer info
      integer ipvt(n)
      integer isamax
      integer j
      integer k
      integer l
      real t

      info = 0
      do k = 1, n - 1
c
c  Find the pivot index L.
c
        l = k
        t = abs ( a(k,k) )
        do i = k + 1, n
          if ( t .lt. abs ( a(i,k) ) ) then
            l = i
            t = abs ( a(i,k) )
          end if
        end do
        ipvt(k) = l

        if ( a(l,k) .eq. 0.0E+00 ) then
          info = k
          return
        end if
c
c  Interchange rows K and L.
c
        if ( l .ne. k ) then
          do j = k, n
            t      = a(l,j)
            a(l,j) = a(k,j)
            a(k,j) = t
          end do
        end if
c
c  Compute column K of the lower triangular factor.
c  Add multiples of the pivot row to the remaining rows.
c
c$omp parallel 
c$omp& shared ( a, k, n )
c$omp& private ( i, j )

c$omp do
        do i = k + 1, n
          a(i,k) = - a(i,k)  / a(k,k)
          do j = k + 1, n
            a(i,j) = a(i,j) + a(i,k) * a(k,j)
          end do
        end do
c$omp end do

c$omp end parallel

      end do

      ipvt(n) = n
      if ( a(n,n) .eq. 0.0E+00 ) then
        info = n
      end if

      return
      end
      subroutine sgesl(a,lda,n,ipvt,b,job)

c*********************************************************************72
c
cc SGESL solves a linear system factored by SGEFA.
c
c     sgesl solves the real system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by sgeco or sgefa.
c
c     on entry
c
c        a       real(lda, n)
c                the output from sgeco or sgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from sgeco or sgefa.
c
c        b       real(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if sgeco has set rcond .gt. 0.0
c        or sgefa has set info .eq. 0 .
c
      implicit none

      integer lda
      integer n

      real a(lda,n)
      real b(1)
      integer i
      integer ipvt(n)
      integer j
      integer job
      integer k
      integer l
      real t

      if ( job .eq. 0 ) then

        do k = 1, n - 1
          l = ipvt(k)
          t = b(l)
          if ( l .ne. k ) then
             b(l) = b(k)
             b(k) = t
          end if
          do i = k + 1, n
            b(i) = b(i) + t * a(i,k)
          end do
        end do

        do k = n, 1, -1
          b(k) = b(k) / a(k,k)
          do i = 1, k - 1
            b(i) = b(i) - a(i,k) * b(k)
          end do
        end do

      end if

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
