      program main

c*********************************************************************72
c
cc DISLIN_EX04 demonstrates various interpolation methods for data.
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

      character*8 cpol(6)
      character*60 ctit
      integer i
      integer nx
      integer nxposn
      integer ny
      integer nya
      integer nyposn
      real x(16)
      real y(16)

      data x / 0.0, 1.0, 3.0, 4.5, 6.0, 8.0, 9.0, 11.0, 12.0, 12.5,
     &  13.0, 15.0, 16.0, 17.0, 19.0, 20.0 /

      data y / 2.0, 4.0, 4.5, 3.0, 1.0, 7.0, 2.0, 3.0, 5.0, 2.0,
     &  2.5, 2.0, 4.0, 6.0, 5.5, 4.0 /

      data cpol / 'SPLINE', 'STEM', 'BARS', 'STAIRS', 'STEP', 'LINEAR' /

      data nya / 2700 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DISLIN_EX04:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Demonstrate the various interpolation'
      write ( *, '(a)' ) '  methods available for (X,Y) data.'

      ctit = 'Interpolation Methods'
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
      call setfil ( 'dislin_ex04.png' )
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
      call incmrk ( 1 )
      call hsymbl ( 25 )
      call titlin ( ctit, 1 )
      call axslen ( 1500, 350 )
      call setgrf ( 'LINE', 'LINE', 'LINE', 'LINE' )

      do i = 1, 6

        call axspos ( 350, nya-(i-1)*350 )
        call polcrv ( cpol(i) )
        call marker ( 0 )

        call graf ( 0.0, 20.0, 0.0, 5.0, 0.0, 10.0, 0.0, 5.0 )
        nx = nxposn ( 1.0 )
        ny = nyposn ( 8.0 )
        call messag ( cpol(i), nx, ny )
        call curve ( x, y, 16 )

        if ( i .eq. 6 ) then
          call height ( 50 )
          call title ( )
        end if

        call endgrf ( )

      end do
c
c  Close DISLIN.
c
      call disfin ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DISLIN_EX04:'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
