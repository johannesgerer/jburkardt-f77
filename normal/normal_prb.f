      program main

c*********************************************************************72
c
cc MAIN is the main program for NORMAL_PRB.
c
c  Discussion:
c
c    NORMAL_PRB calls sample problems for the NORMAL library.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    08 January 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'NORMAL_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version.'
      write ( *, '(a)' ) '  Test the NORMAL library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      if ( .false. ) then
        call test04 ( )
      else
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'SKIPPING TEST04 FOR NOW!'
      end if
      call test05 ( )
      call test06 ( )
      call test07 ( )
      call test08 ( )
      call test09 ( )

      call test10 ( )
      call test11 ( )
      call test12 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'NORMAL_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests C4_NORMAL_01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      complex c4_normal_01
      integer i
      integer seed
      integer seed_init

      seed_init = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  C4_NORMAL_01 computes pseudorandom'
      write ( *, '(a)' ) '  complex values normally distributed '
      write ( *, '(a)' ) '  in the unit circle.'

      seed = seed_init

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
      write ( *, '(a)' ) ' '

      do i = 1, 10
        write ( *, '(2x,i8,2x,2g14.6)' ) i, c4_normal_01 ( seed )
      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests C8_NORMAL_01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double complex c8_normal_01
      integer i
      integer seed
      integer seed_init
      parameter ( seed_init = 123456789 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  C8_NORMAL_01 computes pseudorandom '
      write ( *, '(a)' ) '  double precision complex values '
      write ( *, '(a)' ) '  normally distributed in the unit circle.'

      seed = seed_init

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
      write ( *, '(a)' ) ' '

      do i = 1, 10
        write ( *, '(2x,i6,2x,2g14.6)' ) i, c8_normal_01 ( seed )
      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests I4_NORMAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      real a
      real b
      integer i
      integer i4_normal
      integer seed
      integer seed_init
      parameter ( seed_init = 123456789 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  I4_NORMAL computes integer pseudorandom'
      write ( *, '(a)' ) '  values in an interval [A,B].'

      a = 70.0E+00
      b = 10.0E+00
      seed = seed_init

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  The lower endpoint A = ', a
      write ( *, '(a,g14.6)' ) '  The upper endpoint B = ', b
      write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
      write ( *, '(a)' ) ' '

      do i = 1, 10
        write ( *, '(2x,i8,2x,i8)' ) i, i4_normal ( a, b, seed )
      end do

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests I8_NORMAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision b
      integer*8 i
      integer*8 i8_normal
      integer*8 seed
      integer*8 seed_init
      parameter ( seed_init = 123456789 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  I8_NORMAL computes integer pseudorandom'
      write ( *, '(a)' ) '  values in an interval [A,B].'

      a = 70.0D+00
      b = 10.0D+00
      seed = seed_init

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  The lower endpoint A = ', a
      write ( *, '(a,g14.6)' ) '  The upper endpoint B = ', b
      write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
      write ( *, '(a)' ) ' '

      do i = 1, 10
        write ( *, '(2x,i8,2x,i8)' ) i, i8_normal ( a, b, seed )
      end do

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests R4_NORMAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      real a
      real b
      integer i
      real r4_normal
      integer seed
      integer seed_init
      parameter ( seed_init = 123456789 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  R4_NORMAL computes real pseudorandom values'
      write ( *, '(a)' ) '  in an interval [A,B].'

      a = 5.0E+00
      b = 10.0E+00
      seed = seed_init

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  The lower endpoint A = ', a
      write ( *, '(a,g14.6)' ) '  The upper endpoint B = ', b
      write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
      write ( *, '(a)' ) ' '

      do i = 1, 10
        write ( *, '(2x,i6,2x,g14.6)' ) i, r4_normal ( a, b, seed )
      end do

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 tests R4_NORMAL_01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      real r4_normal_01
      integer seed
      integer seed_init

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) '  R4_NORMAL_01 computes normal pseudorandom'
      write ( *, '(a)' ) '  values in the interval [0,1].'

      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  The initial seed is ', seed
      write ( *, '(a)' ) ' '

      do i = 1, 10
        write ( *, '(2x,i6,2x,g14.6)' ) i, r4_normal_01 ( seed )
      end do

      return
      end
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07 tests R8_NORMAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      double precision a
      double precision b
      integer i
      double precision r8_normal
      integer seed
      integer seed_init

      seed_init = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07'
      write ( *, '(a)' ) '  R8_NORMAL computes pseudonormal values '
      write ( *, '(a)' ) '  with mean A and standard deviation B.'

      a = 10.0D+00
      b = 2.0D+00
      seed = seed_init

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) '  The mean A = ', a
      write ( *, '(a,g14.6)' ) '  The standard deviation B = ', b
      write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
      write ( *, '(a)' ) ' '

      do i = 1, 10
        write ( *, '(2x,i6,2x,g14.6)' ) i, r8_normal ( a, b, seed )
      end do

      return
      end
      subroutine test08 ( )

c*******************************************************************************
c
cc TEST08 tests R8_NORMAL_01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      double precision r8_normal_01
      integer seed
      integer seed_init

      seed_init = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08'
      write ( *, '(a)' ) '  R8_NORMAL_01 computes pseudonormal values '
      write ( *, '(a)' ) '  with mean 0 and standard deviation 1.'

      seed = seed_init

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  The initial seed is ', seed_init
      write ( *, '(a)' ) ' '

      do i = 1, 10
        write ( *, '(2x,i6,2x,g14.6)' ) i, r8_normal_01 ( seed )
      end do

      return
      end
      subroutine test09 ( )

c*********************************************************************72
c
cc TEST09 tests R8_NORMAL_01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n

      parameter ( n = 1000 )

      integer i
      double precision r8_normal_01
      integer seed
      integer seed_in
      integer seed_out
      double precision u(n)
      double precision u_avg
      double precision u_var

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST09'
      write ( *, '(a)' ) '  R8_NORMAL_01 computes a sequence of '
      write ( *, '(a)' ) '  normally distributed pseudorandom numbers.'
c
c  Start with a given seed.
c
      seed = 12345

      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  Initial SEED = ', seed

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  First 10 values:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '     I         Input        Output    R8_NORMAL_01'
      write ( *, '(a)' ) '                SEED          SEED'
      write ( *, '(a)' ) ' '

      do i = 1, 10

        seed_in = seed
        u(i) = r8_normal_01 ( seed )
        seed_out = seed

        write ( *, '(i6,2x,i12,2x,i12,2x,g14.6)' ) 
     &    i, seed_in, seed_out, u(i)

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a,i10,a)' ) '  Now compute ', n, ' elements.'
      write ( *, '(a)' ) ' '

      u_avg = 0.0D+00
      do i = 1, n
        u(i) = r8_normal_01 ( seed )
        u_avg = u_avg + u(i)
      end do

      u_avg = u_avg / dble ( n )

      u_var = 0.0D+00
      do i = 1, n
        u_var = u_var + ( u(i) - u_avg )**2
      end do

      u_var = u_var / dble ( n - 1 )

      write ( *, '(a)' ) ' '
      write ( *, '(a,f10.6)' ) '  Average value = ', u_avg
      write ( *, '(a,f10.6)' ) '  Expecting       ', 0.0D+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,f10.6)' ) '  Variance =      ', u_var
      write ( *, '(a,f10.6)' ) '  Expecting       ', 1.0D+00

      return
      end
      subroutine test10 ( )

c*********************************************************************72
c
cc TEST10 tests R8_NORMAL_01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer i
      double precision r8_normal_01
      integer seed
      integer seed_init
      integer seed_input
      integer seed_output
      double precision value

      seed_init = 123456789

      write ( *, '(a)' ) '  R8_NORMAL_01 computes pseudonormal values '
      write ( *, '(a)' ) '  with mean 0 and standard deviation 1.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Verify that we can change the seed'
      write ( *, '(a)' ) '  and get the desired results.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  The initial seed is ', seed_init

      seed = seed_init

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '         I    Seed(in)   Seed(out)   R8_NORMAL_01'
      write ( *, '(a)' ) ' '

      do i = 1, 5

        seed_input = seed
        value = r8_normal_01 ( seed )
        seed_output = seed

        write ( *, '(2x,i8,2x,i12,2x,i12,2x,g14.6)' ) 
     &    i, seed_input, seed_output, value

      end do

      seed = seed_init

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Resetting seed to repeat, after an ODD number of steps.'
      write ( *, '(a)' ) ' '

      do i = 6, 10

        seed_input = seed
        value = r8_normal_01 ( seed )
        seed_output = seed

        write ( *, '(2x,i8,2x,i12,2x,i12,2x,g14.6)' ) 
     &    i, seed_input, seed_output, value

      end do

      seed = seed_init

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Resetting seed to repeat, after an EVEN number of steps.'
      write ( *, '(a)' ) ' '

      do i = 11, 15

        seed_input = seed
        value = r8_normal_01 ( seed )
        seed_output = seed

        write ( *, '(2x,i8,2x,i12,2x,i12,2x,g14.6)' ) 
     &    i, seed_input, seed_output, value

      end do

      return
      end
      subroutine test11 ( )

c*********************************************************************72
c
cc TEST11 tests R8_NORMAL_01 and R8MAT_NORMAL_01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      integer n

      parameter ( m = 100 )
      parameter ( n = 10 )

      double precision a(m,n)
      double precision b(m,n)
      integer i
      integer j
      integer k
      double precision r8_normal_01
      integer seed
      integer seed_init

      seed_init = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST11'
      write ( *, '(a)' ) '  R8_NORMAL_01 computes pseudorandom normal '
      write ( *, '(a)' ) '    values one at a time.'
      write ( *, '(a)' ) 
     &  '  R8MAT_NORMAL_01 computes a matrix of values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  For the same initial seed, the results'
      write ( *, '(a)' ) '  should be identical, but R8MAT_NORMAL_01 '
      write ( *, '(a)' ) '  might be faster.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  The initial seed is ', seed_init

      seed = seed_init

      do j = 1, n
        do i = 1, m
          a(i,j) = r8_normal_01 ( seed )
        end do
      end do

      seed = seed_init;
      call r8mat_normal_01 ( m, n, seed, b )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '       I       J      A(I,J)        B(I,J)'
      write ( *, '(a)' ) 
     &  '                   (R8_NORMAL_01)  (R8MAT_NORMAL_01)'
      write ( *, '(a)' ) ' '

      do k = 0, 10
        i = ( k * m + ( 10 - k ) * 1 ) / 10
        j = ( k * n + ( 10 - k ) * 1 ) / 10
        write ( *, '(2x,i8,2x,i8,2x,g14.6,2x,g14.6)' ) 
     &    i, j, a(i,j), b(i,j)
      end do
  
      return
      end
      subroutine test12 ( )

c*********************************************************************72
c
cc TEST12 tests R8_NORMAL_01 and R8VEC_NORMAL_01.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 June 2006
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      double precision a(n)
      double precision b(n)
      integer i
      double precision r8_normal_01
      integer seed
      integer seed_init
      parameter ( seed_init = 123456789 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST12'
      write ( *, '(a)' ) '  R8_NORMAL_01 computes pseudorandom normal '
      write ( *, '(a)' ) '  values one at a time.'
      write ( *, '(a)' ) 
     &  '  R8VEC_NORMAL_01 computes a vector of values.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  For the same initial seed, the results'
      write ( *, '(a)' ) '  should be identical, but R8VEC_NORMAL_01'
      write ( *, '(a)') '   might be faster.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i12)' ) '  The initial seed is ', seed_init

      seed = seed_init

      do i = 1, n
        a(i) = r8_normal_01 ( seed )
      end do

      seed = seed_init;
      call r8vec_normal_01 ( n, seed, b )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '       I      A(I)            B(I)'
      write ( *, '(a)' ) '           (R8_NORMAL_01)  (R8VEC_NORMAL_01)'
      write ( *, '(a)' ) ' '

      do i = 1, 10
        write ( *, '(2x,i6,2x,g14.6,2x,g14.6)' ) i, a(i), b(i)
      end do
  
      return
      end
