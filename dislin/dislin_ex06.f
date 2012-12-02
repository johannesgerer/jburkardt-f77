      program main

c*********************************************************************72
c
cc DISLIN_EX06 demonstrates the creation of pie charts.
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

      character*40 cbuf
      character*60 ctit
      integer i
      integer nya
      real xray(5)

      data xray / 1.0, 2.5, 2.0, 2.7, 1.8 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DISLIN_EX06:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Demonstrate the creation of pie charts.'

      ctit = 'Pie Charts (PIEGRF)'
      nya = 2800
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
      call setfil ( 'dislin_ex06.png' )
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
      call axslen ( 1600, 1000 )
      call titlin ( ctit, 2 )

      call legini ( cbuf, 5, 8 )
      call leglin ( cbuf, 'FIRST', 1 )
      call leglin ( cbuf, 'SECOND', 2 )
      call leglin ( cbuf, 'THIRD', 3 )
      call leglin ( cbuf, 'FOURTH', 4 )
      call leglin ( cbuf, 'FIFTH', 5 )
c
c     Selecting shading patterns.
c
      call patcyc ( 1, 7 )
      call patcyc ( 2, 4 )
      call patcyc ( 3, 13 )
      call patcyc ( 4, 3 )
      call patcyc ( 5, 5 )

      do i = 1, 2

        call axspos ( 250, nya-(i-1)*1200 )

        if ( i .eq. 2 ) then
          call labels ( 'DATA', 'PIE' )
          call labpos ( 'EXTERNAL', 'PIE' )
        end if

        call piegrf ( cbuf, 1, xray, 5 )

        if ( i .eq. 2 ) then
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
      write ( *, '(a)' ) 'DISLIN_EX06:'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
