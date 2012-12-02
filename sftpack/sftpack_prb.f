      program main

c*********************************************************************72
c
cc MAIN is the main program for SFTPACK_PRB.
c
c  Discussion:
c
c    SFTPACK_PRB calls sample problems for the SFTPACK library.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 June 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      write ( *, '(a)' ) ' '
      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SFTPACK_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the SFTPACK library.'

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
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'SFTPACK_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests R8VEC_SCT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 February 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 256 )

      double precision ahi
      parameter ( ahi = 5.0D+00 )
      double precision alo
      parameter ( alo = 0.0D+00 )
      double precision c(n)
      double precision d(n)
      double precision e(n)
      integer i
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) 
     &  '  For double precision slow cosine transforms,'
      write ( *, '(a)' ) 
     &  '  R8VEC_SCT does a forward or backward transform.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  The number of data items is N = ', n
c
c  Set the data values.
c
      seed = 123456789

      call r8vec_uniform ( n, alo, ahi, seed, c )

      call r8vec_print_part ( n, c, 10, '  The original data:' )
c
c  Compute the coefficients.
c
      call r8vec_sct ( n, c, d )

      call r8vec_print_part ( n, d, 10, '  The cosine coefficients:' )
c
c  Now compute inverse transform of coefficients.  Should get back the
c  original data.

      call r8vec_sct ( n, d, e )

      do i = 1, n
        e(i) = e(i) / dble ( 2 * n )
      end do

      call r8vec_print_part ( n, e, 10, '  The retrieved data:' )

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests R8VEC_SFTB and R8VEC_SFTF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 February 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 36 )

      double precision a(n/2)
      double precision ahi
      parameter ( ahi = 5.0D+00 )
      double precision alo
      parameter ( alo = 0.0D+00 )
      double precision azero
      double precision b(n/2)
      integer i
      integer seed
      double precision x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  For real slow Fourier transforms,'
      write ( *, '(a)' ) '  R8VEC_SFTF computes the forward transform.'
      write ( *, '(a)' ) '  R8VEC_SFTB computes the backward transform.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The number of data values, N = ', n

      seed = 123456789

      call r8vec_uniform ( n, alo, ahi, seed, x )

      call r8vec_print_part ( n, x, 10, '  The original data:' )
c
c  Compute the slow Fourier transform of the data.
c
      call r8vec_sftf ( n, x, azero, a, b )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  A (cosine) coefficients:'
      write ( *, '(a)' ) ' '

      write ( *, '(2x,i3,g14.6)' ) 0, azero

      do i = 1, n/2
        write ( *, '(2x,i3,g14.6)' ) i, a(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  B (sine) coefficients:'
      write ( *, '(a)' ) ' '

      do i = 1, n/2
        write ( *, '(2x,i3,g14.6)' ) i, b(i)
      end do
c
c  Now try to retrieve the data from the coefficients.
c
      call r8vec_sftb ( n, azero, a, b, x )

      call r8vec_print_part ( n, x, 10, '  The retrieved data:' )

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests R8VEC_SHT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 February 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 17 )

      double precision ahi
      parameter ( ahi = 5.0D+00 )
      double precision alo
      parameter ( alo = 0.0D+00 )
      double precision c(n)
      double precision d(n)
      double precision e(n)
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  For real slow Hartley transforms,'
      write ( *, '(a)' ) 
     &  '  R8VEC_SHT does a forward or backward transform.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  The number of data items is N = ', n
c
c  Set the data values.
c
      seed = 123456789

      call r8vec_uniform ( n, alo, ahi, seed, c )

      call r8vec_print_part ( n, c, 10, '  The original data:' )
c
c  Compute the coefficients.
c
      call r8vec_sht ( n, c, d )

      call r8vec_print_part ( n, d, 10, '  The Hartley coefficients:' )
c
c  Now compute inverse transform of coefficients.  Should get back the
c  original data.

      call r8vec_sht ( n, d, e )

      call r8vec_print_part ( n, e, 10, '  The retrieved data:' )

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests R8VEC_SQCTB and R8VEC_SQCTF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 February 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 256 )

      double precision ahi
      parameter ( ahi = 5.0D+00 )
      double precision alo
      parameter ( alo = 0.0D+00 )
      integer seed
      double precision x(n)
      double precision y(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) 
     &  '  For real slow quarter wave cosine transforms,'
      write ( *, '(a)' ) '  R8VEC_SQCTF does a forward transform;'
      write ( *, '(a)' ) '  R8VEC_SQCTB does a backward transform.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  The number of data items is N = ', n
c
c  Set the data values.
c
      seed = 123456789

      call r8vec_uniform ( n, alo, ahi, seed, x )

      call r8vec_print_part ( n, x, 10, '  The original data:' )
c
c  Compute the coefficients.
c
      call r8vec_sqctf ( n, x, y )

      call r8vec_print_part ( n, y, 10, '  The cosine coefficients:' )
c
c  Now compute inverse transform of coefficients.  Should get back the
c  original data.

      call r8vec_sqctb ( n, y, x )

      call r8vec_print_part ( n, x, 10, '  The retrieved data:' )

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests R8VEC_SQSTB and R8VEC_SQSTF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 February 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 256 )

      double precision ahi
      parameter ( ahi = 5.0D+00 )
      double precision alo
      parameter ( alo = 0.0D+00 )
      integer seed
      double precision x(n)
      double precision y(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  For real slow quarter wave sine transforms,'
      write ( *, '(a)' ) '  R8VEC_SQSTF does a forward transform;'
      write ( *, '(a)' ) '  R8VEC_SQSTB does a backward transform.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  The number of data items is N = ', n
c
c  Set the data values.
c
      seed = 123456789

      call r8vec_uniform ( n, alo, ahi, seed, x )

      call r8vec_print_part ( n, x, 10, '  The original data:' )
c
c  Compute the coefficients.
c
      call r8vec_sqstf ( n, x, y )

      call r8vec_print_part ( n, y, 10, '  The sine coefficients:' )
c
c  Now compute inverse transform of coefficients.  Should get back the
c  original data.

      call r8vec_sqstb ( n, y, x )

      call r8vec_print_part ( n, x, 10, '  The retrieved data:' )

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 tests R8VEC_SST.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 February 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 256 )

      double precision ahi
      parameter ( ahi = 5.0D+00 )
      double precision alo
      parameter ( alo = 0.0D+00 )
      double precision c(n)
      double precision d(n)
      double precision e(n)
      integer i
      integer seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) '  For slow sine transforms,'
      write ( *, '(a)' ) 
     &  '  R8VEC_SST does a forward or backward transform.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  The number of data items is N = ', n
c
c  Set the data values.
c
      seed = 123456789

      call r8vec_uniform ( n, alo, ahi, seed, c )

      call r8vec_print_part ( n, c, 10, '  The original data:' )
c
c  Compute the coefficients.
c
      call r8vec_sst ( n, c, d )

      call r8vec_print_part ( n, d, 10, '  The sine coefficients:' )
c
c  Now compute inverse transform of coefficients.  Should get back the
c  original data.

      call r8vec_sst ( n, d, e )

      do i = 1, n
        e(i) = e(i) / dble ( 2 * ( n + 1 ) )
      end do

      call r8vec_print_part ( n, e, 10, '  The retrieved data:' )

      return
      end
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07 tests C4VEC_SFTB and C4VEC_SFTF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 June 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 36 )

      integer i
      integer seed
      complex x(n)
      complex x2(n)
      complex y(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07'
      write ( *, '(a)' ) '  For complex slow Fourier transforms,'
      write ( *, '(a)' ) '  C4VEC_SFTF computes the forward transform.'
      write ( *, '(a)' ) '  C4VEC_SFTB computes the backward transform.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The number of data values, N = ', n

      seed = 123456789

      call c4vec_uniform_01 ( n, seed, x )

      call c4vec_print_part ( n, x, 10, '  The data X:' )
!
!  Compute the slow Fourier transform of the data.
!
      call c4vec_sftf ( n, x, y )

      call c4vec_print_part ( n, y, 10, 
     &  '  The Fourier coefficients Y:' )

      call c4vec_sftb ( n, y, x2 )

      call c4vec_print_part ( n, x2, 10, '  The recovered data:' )

      return
      end
      subroutine test08 ( )

c*********************************************************************72
c
cc TEST08 tests C8VEC_SFTB and C8VEC_SFTF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 June 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 36 )

      integer i
      integer seed
      double complex x(n)
      double complex x2(n)
      double complex y(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08'
      write ( *, '(a)' ) '  For complex slow Fourier transforms,'
      write ( *, '(a)' ) '  C8VEC_SFTF computes the forward transform.'
      write ( *, '(a)' ) '  C8VEC_SFTB computes the backward transform.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The number of data values, N = ', n

      seed = 123456789

      call c8vec_uniform_01 ( n, seed, x )

      call c8vec_print_part ( n, x, 10, '  The data X:' )
!
!  Compute the slow Fourier transform of the data.
!
      call c8vec_sftf ( n, x, y )

      call c8vec_print_part ( n, y, 10, 
     &  '  The Fourier coefficients Y:' )

      call c8vec_sftb ( n, y, x2 )

      call c8vec_print_part ( n, x2, 10, '  The recovered data:' )

      return
      end
      subroutine test09 ( )

c*********************************************************************72
c
cc TEST09 tests R4VEC_SFTB and R4VEC_SFTF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    10 June 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 36 )

      real a(n/2)
      real ahi
      parameter ( ahi = 5.0E+00 )
      real alo
      parameter ( alo = 0.0E+00 )
      real azero
      real b(n/2)
      integer i
      integer seed
      real x(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST09'
      write ( *, '(a)' ) '  For real slow Fourier transforms,'
      write ( *, '(a)' ) '  R4VEC_SFTF computes the forward transform.'
      write ( *, '(a)' ) '  R4VEC_SFTB computes the backward transform.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The number of data values, N = ', n

      seed = 123456789

      call r4vec_uniform ( n, alo, ahi, seed, x )

      call r4vec_print_part ( n, x, 10, '  The original data:' )
c
c  Compute the slow Fourier transform of the data.
c
      call r4vec_sftf ( n, x, azero, a, b )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  A (cosine) coefficients:'
      write ( *, '(a)' ) ' '

      write ( *, '(2x,i3,g14.6)' ) 0, azero

      do i = 1, n/2
        write ( *, '(2x,i3,g14.6)' ) i, a(i)
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  B (sine) coefficients:'
      write ( *, '(a)' ) ' '

      do i = 1, n/2
        write ( *, '(2x,i3,g14.6)' ) i, b(i)
      end do
c
c  Now try to retrieve the data from the coefficients.
c
      call r4vec_sftb ( n, azero, a, b, x )

      call r4vec_print_part ( n, x, 10, '  The retrieved data:' )

      return
      end
      subroutine test10 ( )

c*********************************************************************72
c
cc TEST10 tests C4MAT_SFTB and C4MAT_SFTF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    21 June 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n1
      parameter ( n1 = 10 )
      integer n2
      parameter ( n2 = 4 )

      integer i
      integer seed
      complex x(n1,n2)
      complex x2(n1,n2)
      complex y(n1,n2)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST10'
      write ( *, '(a)' ) '  For complex slow Fourier transforms,'
      write ( *, '(a)' ) '  C4MAT_SFTF computes the forward transform.'
      write ( *, '(a)' ) '  C4MAT_SFTB computes the backward transform.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a,i8)' ) 
     &  '  The data has dimension N1 = ', n1, ' by N2 = ', n2

      seed = 123456789

      call c4mat_uniform_01 ( n1, n2, seed, x )

      call c4mat_print_some ( n1, n2, x, 1, 1, 10, 10, '  The data X:' )
c
c  Compute the slow Fourier transform of the data.
c
      call c4mat_sftf ( n1, n2, x, y )

      call c4mat_print_some ( n1, n2, y, 1, 1, 10, 10, 
     &  '  The Fourier coefficients Y:' )

      call c4mat_sftb ( n1, n2, y, x2 )

      call c4mat_print_some ( n1, n2, x2, 1, 1, 10, 10, 
     &  '  The recovered data:' )

      return
      end
      subroutine test11 ( )

c*********************************************************************72
c
cc TEST11 tests C8MAT_SFTB and C8MAT_SFTF.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 June 2010
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n1
      parameter ( n1 = 10 )
      integer n2
      parameter ( n2 = 4 )

      integer i
      integer seed
      double complex x(n1,n2)
      double complex x2(n1,n2)
      double complex y(n1,n2)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST11'
      write ( *, '(a)' ) '  For complex slow Fourier transforms,'
      write ( *, '(a)' ) '  C8MAT_SFTF computes the forward transform.'
      write ( *, '(a)' ) '  C8MAT_SFTB computes the backward transform.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8,a,i8)' ) 
     &  '  The data has dimension N1 = ', n1, ' by N2 = ', n2

      seed = 123456789

      call c8mat_uniform_01 ( n1, n2, seed, x )

      call c8mat_print_some ( n1, n2, x, 1, 1, 10, 10, '  The data X:' )
c
c  Compute the slow Fourier transform of the data.
c
      call c8mat_sftf ( n1, n2, x, y )

      call c8mat_print_some ( n1, n2, y, 1, 1, 10, 10, 
     &  '  The Fourier coefficients Y:' )

      call c8mat_sftb ( n1, n2, y, x2 )

      call c8mat_print_some ( n1, n2, x2, 1, 1, 10, 10, 
     &  '  The recovered data:' )

      return
      end
