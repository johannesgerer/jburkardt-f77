      program main

c*********************************************************************72
c
cc DISLIN_EX07 demonstrates 3D bar graphs and pie charts.
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

      character*80 cbuf
      integer ic1ray(5)
      integer ic2ray(5)
      real xray(5)
      real y1ray(5)
      real y2ray(5)

      data xray / 2.0, 4.0, 6.0, 8.0, 10.0 /
      data y1ray / 0.0, 0.0, 0.0, 0.0, 0.0 /
      data y2ray / 3.2, 1.5, 2.0, 1.0, 3.0 /
      data ic1ray / 50, 150, 100, 200, 175 /
      data ic2ray / 50, 150, 100, 200, 175 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DISLIN_EX07:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Demonstrate the creation of 3D bar '
      write ( *, '(a)' ) '  and pie charts.'
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
      call setfil ( 'dislin_ex07.png' )
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

      call titlin ( '3-D Bar Graph / 3-D Pie Chart', 2 )
      call htitle ( 40 )

      call shdpat ( 16 )
      call axslen ( 1500, 1000 )
      call axspos ( 300, 1400 )

      call barwth ( 0.5 )
      call bartyp ( '3DVERT' )
      call labels ( 'SECOND', 'BARS' )
      call labpos ( 'OUTSIDE', 'BARS' )
      call labclr ( 255, 'BARS' )
      call graf ( 0.0, 12.0, 0.0, 2.0, 0.0, 5.0, 0.0, 1.0 )
      call title ( )
      call color ( 'RED' )
      call bars ( xray, y1ray, y2ray, 5 )
      call endgrf ( )

      call shdpat ( 16 )
      call labels ( 'DATA', 'PIE' )
      call labclr ( 255, 'PIE' )
      call chnpie ( 'NONE' )
      call pieclr ( ic1ray, ic2ray, 5 )
      call pietyp ( '3D' )
      call axspos ( 300, 2700 )
      call piegrf ( cbuf, 0, y2ray, 5 )
c
c  Close DISLIN.
c     
      call disfin ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DISLIN_EX07:'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
