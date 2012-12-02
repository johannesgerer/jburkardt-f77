      program main

c*********************************************************************72
c
cc DISLIN_EX07B demonstrates 3D bar graphs.
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
      parameter ( n = 18 )

      character*80 cbuf
      integer i
      integer icray(n)
      real xray(n)
      real xwray(n)
      real yray(n)
      real ywray(n)
      real z1ray(n)
      real z2ray(n)

      data icray / 30, 30, 30, 30, 30, 30, 100, 100, 100, 100,
     &  100, 100, 170, 170, 170, 170, 170, 170 /
      data xray / 1.0, 3.0, 8.0, 1.5, 9.0, 6.3, 5.8, 2.3, 8.1, 3.5,
     &  2.2, 8.7, 9.2, 4.8, 3.4, 6.9, 7.5, 3.8 /
      data yray / 5.0, 8.0, 3.5, 2.0, 7.0, 1.0, 4.3, 7.2, 6.0, 8.5,
     &  4.1, 5.0, 7.3, 2.8, 1.6, 8.9, 9.5, 3.2 /
      data z1ray / 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
     &  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /
      data z2ray / 4.0, 5.0, 3.0, 2.0, 3.5, 4.5, 2.0, 1.6, 3.8, 4.7,
     &  2.1, 3.5, 1.9, 4.2, 4.9, 2.8, 3.6, 4.3 / 

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DISLIN_EX07B:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Demonstrate the creation of 3D bar graphs.'

      do i = 1, n
        xwray(i) = 0.5
        ywray(i) = 0.5
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
      call setfil ( 'dislin_ex07b.png' )
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
c  Use the HARDWARE font.
c
      call hwfont ( )
      call axspos ( 200, 2600 )
      call axslen ( 1800, 1800 ) 

      call name ( 'X-axis', 'X' )
      call name ( 'Y-axis', 'Y' )
      call name ( 'Z-axis', 'Z' )

      call titlin('3-D Bars / BARS3D', 3 )

      call labl3d ( 'HORI' )
      call graf3d ( 0.0, 10.0, 0.0, 2.0, 0.0, 10.0, 0.0, 2.0, 0.0, 
     &  5.0, 0.0, 1.0 )
      call grid3d ( 1, 1, 'BOTTOM' )

      call bars3d ( xray, yray, z1ray, z2ray, xwray, ywray, icray, n ) 

      call legini ( cbuf, 3, 20 )
      call legtit ( ' ' )
      call legpos ( 1300, 1100 )
      call leglin ( cbuf, 'First', 1 )
      call leglin ( cbuf, 'Second', 2 )
      call leglin ( cbuf, 'Third', 3 )
      call legend ( cbuf, 3 )

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
      write ( *, '(a)' ) 'DISLIN_EX07B:'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
