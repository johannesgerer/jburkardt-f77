      program main

c*********************************************************************72
c
cc DISLIN_EX05 demonstrates the creation of bar graphs.
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

      character*24 cbuf
      character*60 ctit
      integer i
      integer nya
      real x(9)
      real y(9)
      real y1(9)
      real y2(9)
      real y3(9)

      data x  / 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0 / 
      data y  / 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 /
      data y1 / 1.0, 1.5, 2.5, 1.3, 2.0, 1.2, 0.7, 1.4, 1.1 /
      data y2 / 2.0, 2.7, 3.5, 2.1, 3.2, 1.9, 2.0, 2.3, 1.8 /
      data y3 / 4.0, 3.5, 4.5, 3.7, 4.0, 2.9, 3.0, 3.2, 2.6 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DISLIN_EX05:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Demonstrate the creation of bar graphs.'

      nya = 2700
      ctit = 'Bar Graphs (BARS)'
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
      call setfil ( 'dislin_ex05.png' )
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
      call ticks ( 1, 'X' )
      call intax ( )
      call axslen ( 1600, 700 )
      call titlin ( ctit, 3 )

      call legini ( cbuf, 3, 8 )
      call leglin ( cbuf, 'FIRST', 1 )
      call leglin ( cbuf, 'SECOND', 2 )
      call leglin ( cbuf, 'THIRD', 3 )
      call legtit ( ' ' )

      call shdpat ( 5 ) 

      do i = 1, 3

        if ( 1 .lt. i ) then
          call labels ( 'NONE', 'X' )
        end if

        call axspos ( 300, nya-(i-1)*800 )
  
        call graf ( 0.0, 10.0, 0.0, 1.0, 0.0, 5.0, 0.0, 1.0 )
  
        if ( i .eq. 1 ) then
          call bargrp ( 3, 0.15 )
          call bars ( X, Y, Y1, 9 )
          call bars ( X, Y, Y2, 9 )
          call bars ( X, Y, Y3, 9 )
          call reset ( 'BARGRP' )
        else if ( i .eq. 2 ) then
          call height ( 30 )
          call labels ( 'DELTA', 'BARS' )
          call labpos ( 'CENTER', 'BARS' )
          call bars ( x, y, y1, 9 )
          call bars ( x, y1, y2, 9 )
          call bars ( x, y2, y3, 9 )
          call height ( 36 )
        else if ( i .eq. 3 ) then
          call labels ( 'SECOND', 'BARS' )
          call labpos ( 'OUTSIDE', 'BARS' )
          call bars ( x, y, y1, 9 )
        end if

        if ( i .lt. 3 ) then
          call legend ( cbuf, 7 )
        else
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
      write ( *, '(a)' ) 'DISLIN_EX05:'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
