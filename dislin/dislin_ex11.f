      program main

c*********************************************************************72
c
cc DISLIN_EX11 demonstrates a contour plot.
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

      real fpi
      integer i
      integer j
      real step
      real xray(n)
      real yray(n)
      real zlev
      real zmat(n,n)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DISLIN_EX11:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Demonstrate the creation of '
      write ( *, '(a)' ) '  a contour plot.'

      fpi = 3.14159 / 180.0
      step = 360.0 / real ( n - 1 )

      do i = 1, n
        xray(i) = real ( i - 1 ) * step
        yray(i) = real ( i - 1 ) * step
      end do

      do i = 1, n
        do j = 1, n
          zmat(i,j) = 2.0 * sin ( xray(i) * fpi ) 
     &                     * sin ( yray(j) * fpi )
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
      call setfil ( 'dislin_ex11.png' )
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

      call titlin ( 'Contour Plot', 1 )
      call titlin ( 'F(X,Y) = 2 * SIN(X) * SIN(Y)', 3 )

      call name ( 'X-axis', 'X' )
      call name ( 'Y-axis', 'Y' )

      call intax ( )
      call axspos ( 450, 2670 )
      call graf ( 0.0, 360.0, 0.0, 90.0, 0.0, 360.0, 0.0, 90.0 )

      call height ( 30 )

      do i = 1, 9
        call setclr ( i * 25 )
        zlev = -2.0 + ( i - 1 ) * 0.5
        if ( i .eq. 5 ) then
          call labels ( 'NONE', 'CONTUR' )
        else
          call labels ( 'FLOAT', 'CONTUR' )
        end if
        call contur ( xray, n, yray, n, zmat, zlev )
      end do

      call height ( 50 )
      call color ( 'FORE' )
      call title ( )
c
c  Close DISLIN.
c
      call disfin ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DISLIN_EX11:'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
