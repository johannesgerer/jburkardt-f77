      program main

c*********************************************************************72
c
cc DISLIN_EX13 demonstrates a map plot.
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

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DISLIN_EX13:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Demonstrate the creation of a map plot.'
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
      call setfil ( 'dislin_ex13.png' )
c
c  Choose the page size and orientation.
c
      call setpag ( 'usal' )
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

      call frame ( 3 )
      call axspos ( 400, 1850 )
      call axslen ( 2400, 1400 )

      call name ( 'Longitude', 'X' )
      call name ( 'Latitude', 'Y' )
      call titlin( 'World Coastlines and Lakes', 3 )

      call labels ( 'MAP', 'XY' )
      call grafmp ( -180.0, 180.0, -180.0, 90.0, -90.0, 90.0, 
     &   -90.0, 30.0 )

      call gridmp ( 1, 1 )
      call color ( 'GREEN' )
      call world ( )
      call color ( 'FORE' )

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
      write ( *, '(a)' ) 'DISLIN_EX13:'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
