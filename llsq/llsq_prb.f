      program main

c*********************************************************************72
c
cc LLSQ_PRB tests LLSQ.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 March 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LLSQ_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the LLSQ library.'

      call test01 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LLSQ_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 calls LLSQ to match 15 data values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    08 March 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 15 )

      double precision a
      double precision b
      double precision error
      integer i
      double precision x(n)
      double precision y(n)

      save x
      save y

      data x /
     &  1.47, 1.50, 1.52, 1.55, 1.57, 
     &  1.60, 1.63, 1.65, 1.68, 1.70, 
     &  1.73, 1.75, 1.78, 1.80, 1.83 /
      data y /
     &  52.21, 53.12, 54.48, 55.84, 57.20, 
     &  58.57, 59.93, 61.29, 63.11, 64.47, 
     &  66.28, 68.10, 69.92, 72.19, 74.46 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) 
     &  '  LLSQ can compute the formula for a line of the form'
      write ( *, '(a)' ) '    y = A * x + B'
      write ( *, '(a)' ) 
     &  '  which minimizes the RMS error to a set of N data values.'

      call llsq ( n, x, y, a, b )

      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6,a,g14.6)' ) 
     &  '  Estimated relationship is y = ', a, ' * x + ', b
      write ( *, '(a)' ) 
     &  '  Expected value is         y = 61.272 * x - 39.062'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I      X       Y      B+A*X    |error|'
      write ( *, '(a)' ) ' '
      error = 0.0D+00
      do i = 1, n
        write ( *, '(2x,i4,2x,f7.4,2x,f7.4,2x,f7.4,2x,f7.4)' ) 
     &    i, x(i), y(i), b + a * x(i), b + a * x(i) - y(i)
        error = error + ( b + a * x(i) - y(i) )**2
      end do
      error = sqrt ( error / dble ( n ) )
      write ( *, '(a)' ) ' '
      write ( *, '(a,g14.6)' ) 
     &  '  RMS error =                     ', error

      return
      end

