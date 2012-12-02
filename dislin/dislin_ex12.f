      program main

c*********************************************************************72
c
cc DISLIN_EX12 demonstrates a shaded contour plot.
c
c  Modified:
c
c    09 April 2011
c
c  Reference:
c
c    Helmut Michels,
c    The Data Plotting Software DISLIN - version 10.4,
c    Shaker Media GmbH, January 2010,
c    ISBN13: 978-3-86858-517-9.
c
      implicit none

      integer n
      parameter ( n = 50 )

      integer i
      integer j
      real step
      real x
      real xray(n)
      real y
      real yray(n)
      real zlev(12)
      real zmat(n,n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DISLIN_EX12:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Demonstrate the creation of '
      write ( *, '(a)' ) '  a shaded contour plot.'

      step = 1.6 / real ( n - 1 )

      do i = 1, n
        x = real ( i - 1 ) * step
        xray(i) = x
        do j = 1, n
          y = real ( j - 1 ) * step
          yray(j) = y
          zmat(i,j) = ( x**2 - 1.0 )**2 + ( y**2 - 1.0 )**2
        end do
      end do
c
c  Specify the format of the output file.
c
      call metafl ( 'png' )
c
c  Indicate that new data overwrites old data.
c
      call filmod ( 'delete' )
c
c  Specify the name of the output graphics file.
c
      call setfil ( 'dislin_ex12.png' )
c
c  Choose the page size and orientation.
c
      call setpag ( 'usap' )
c
c  For PNG output, reverse the default black background to white.
c
      call scrmod ( 'reverse' )
c
c  Open DISLIN.
c
      call disini ( )
c
c  Plot a border around the page.
c
      call pagera ( )
c
c  Use the COMPLEX font.
c
      call complx ( )

      call mixalf ( )
      call titlin ( 'Shaded Contour Plot', 1 )
      call titlin ( 'F(X,Y) = (X[2$ - 1)[2$ + (Y[2$ - 1)[2$', 3 )
      call name ( 'X-axis', 'X' )
      call name ( 'Y-axis', 'Y' )

      call shdmod ( 'POLY', 'CONTUR' )
      call axspos ( 450, 2670 )
      call graf ( 0.0, 1.6, 0.0, 0.2, 0.0, 1.6, 0.0, 0.2 )

      do i = 1, 12
        zlev(13-i) = 0.1 + real ( i - 1 ) * 0.1
      end do

      call conshd ( xray, n, yray, n, zmat, zlev, 12 )

      call height ( 50 )
      call title ( )
c
c  Close DISLIN.
c
      call disfin ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DISLIN_EX12:'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
