      program main

c*********************************************************************72
c
cc DISLIN_EX10 demonstrates a 3D surface plot.
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

      character*60 ctit1
      character*60 ctit2
      external zfun

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DISLIN_EX10:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Demonstrate the creation of '
      write ( *, '(a)' ) '  a surface plot.'

      ctit1 = 'Surface Plot (SURFUN)'
      ctit2 = 'F(X,Y) = 2*SIN(X)*SIN(Y)' 
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
      call setfil ( 'dislin_ex10.png' )
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

      call axspos ( 200, 2600 )
      call axslen ( 1800, 1800 )

      call name ( 'X-axis', 'X' )
      call name ( 'Y-axis', 'Y' )
      call name ( 'Z-axis', 'Z' )

      call titlin ( ctit1, 2 )
      call titlin ( ctit2, 4 )

      call view3d ( -5.0, -5.0, 4.0, 'ABS' )

      call graf3d ( 0.0, 360.0, 0.0, 90.0, 0.0, 360.0, 0.0, 90.0, -3.0, 
     &  3.0, -3.0, 1.0 )
      call height ( 50 )
      call title ( )

      call surfun ( zfun, 1, 10.0, 1, 10.0 )
c
c  Close DISLIN.
c
      call disfin ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DISLIN_EX10:'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
      function zfun ( x, y )

c*********************************************************************72
c
cc ZFUN evaluates the function Z(X,Y).
c
c  Modified:
c
c    09 April 2011
c
c  Parameters:
c
c    Input, real X, Y, the evaluation point.
c
c    Output, real ZFUN, the value of Z(X,Y).
c
      real fpi
      real x
      real y
      real zfun

      fpi = 3.14159 / 180.0
      zfun = 2.0 * sin ( x * fpi ) * sin ( y * fpi )

      return
      end
