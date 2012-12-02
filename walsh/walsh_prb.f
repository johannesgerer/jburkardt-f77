      program main

c*********************************************************************72
c
cc MAIN is the main program for WALSH_PRB.
c
c  Discussion:
c
c    WALSH_PRB tests the WALSH library.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 March 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'WALSH_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the WALSH library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'WALSH_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests FWT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 March 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer, parameter :: n = 16

      integer i
      integer j
      integer seed
      double precision w(n)
      double precision work(n)
      double precision x(n)
      double precision y(n)
      double precision z(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  FWT computes a fast Walsh transform.'

      do j = 1, 2

        if ( j == 1 ) then
          seed = 123456789
          call r8vec_uniform_01 ( n, seed, w )
        else
          do i = 1, n
            w(i) = dble ( i )
          end do
        end if

        call r8vec_copy ( n, w, x )
        call fwt ( n, w, work )
        call r8vec_copy ( n, w, y )
        do i = 1, n
           y(i) = y(i) / dble ( n )
        end do
        call fwt ( n, w, work )
        call r8vec_copy ( n, w, z )
        do i = 1, n
           z(i) = z(i) / dble ( n )
        end do

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    '     I        X(I)    Y=FWT(X)/N   Z=FWT(Y)/N'
        write ( *, '(a)' ) ' '
        do i = 1, n
          write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4)' ) 
     &      i, x(i), y(i), z(i)
        end do

      end do
        
      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests WALSH.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 March 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer, parameter :: n = 16

      integer i
      integer j
      integer seed
      double precision w(n)
      double precision work(n)
      double precision x(n)
      double precision y(n)
      double precision z(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  WALSH computes a fast Walsh transform.'

      do j = 1, 2

        if ( j == 1 ) then
          seed = 123456789
          call r8vec_uniform_01 ( n, seed, w )
        else
          do i = 1, n
            w(i) = dble ( i )
          end do
        end if

        call r8vec_copy ( n, w, x )
        call walsh ( n, w, work )
        call r8vec_copy ( n, w, y )
        do i = 1, n
           y(i) = y(i) / dble ( n )
        end do
        call walsh ( n, w, work )
        call r8vec_copy ( n, w, z )
        do i = 1, n
           z(i) = z(i) / dble ( n )
        end do

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    '     I        X(I)    Y=FWT(X)/N   Z=FWT(Y)/N'
        write ( *, '(a)' ) ' '
        do i = 1, n
          write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4)' ) 
     &      i, x(i), y(i), z(i)
        end do

      end do
        
      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests HAAR, HAARIN and HNORM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 March 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer, parameter :: n = 16

      integer i
      integer j
      integer seed
      double precision w(n)
      double precision work(n)
      double precision x(n)
      double precision y(n)
      double precision z(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  HAAR computes a Haar transform.'
      write ( *, '(a)' ) '  HNORM normalizes the transformed data.'
      write ( *, '(a)' ) '  HAARIN computes an inverse Haar transform.'

      do j = 1, 2

        if ( j == 1 ) then
          seed = 123456789
          call r8vec_uniform_01 ( n, seed, w )
        else
          do i = 1, n
            w(i) = dble ( i )
          end do
        end if

        call r8vec_copy ( n, w, x )

        call haar ( n, w, work )

        call r8vec_copy ( n, w, y )

        call hnorm ( n, w )

        call r8vec_copy ( n, w, z )

        call haarin ( n, w, work )

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    '     I        X(I)    Y=HAAR(X)  Z=HNORM(Y)  W=HAARIN(Z)'
        write ( *, '(a)' ) ' '
        do i = 1, n
          write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4,2x,f10.4)' ) 
     &      i, x(i), y(i), z(i), w(i)
        end do

      end do
        
      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests FFWT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    11 March 2011
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer, parameter :: n = 16

      integer i
      integer j
      integer seed
      double precision w(n)
      double precision x(n)
      double precision y(n)
      double precision z(n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  FFWT computes a fast Walsh transform.'

      do j = 1, 2

        if ( j == 1 ) then
          seed = 123456789
          call r8vec_uniform_01 ( n, seed, w )
        else
          do i = 1, n
            w(i) = dble ( i )
          end do
        end if

        call r8vec_copy ( n, w, x )
        call ffwt ( n, w )
        call r8vec_copy ( n, w, y )
        do i = 1, n
           y(i) = y(i) / dble ( n )
        end do
        call ffwt ( n, w )
        call r8vec_copy ( n, w, z )
        do i = 1, n
           z(i) = z(i) / dble ( n )
        end do

        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    '     I        X(I)   Y=FFWT(X)/N  Z=FFWT(Y)/N'
        write ( *, '(a)' ) ' '
        do i = 1, n
          write ( *, '(2x,i4,2x,f10.4,2x,f10.4,2x,f10.4)' ) 
     &      i, x(i), y(i), z(i)
        end do

      end do
        
      return
      end
