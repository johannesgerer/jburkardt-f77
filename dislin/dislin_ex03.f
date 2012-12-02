      program main

c*********************************************************************72
c
cc DISLIN_EX03 demonstrates the creation and display of special symbols.
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

      character*2 cstr
      character*20 ctit
      integer i
      integer nl
      integer nlmess
      integer nxp
      integer ny

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DISLIN_EX03:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Demonstrate the use of the SYMBOL routine'
      write ( *, '(a)' ) '  to create and display special symbols.'

      ctit = 'Symbols'
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
      call setfil ( 'dislin_ex03.png' )
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

      call height ( 60 )

      nl = nlmess ( ctit )
      call messag ( ctit, ( 2100 - nl ) / 2, 200 )

      call height ( 50 )
      call hsymbl ( 120 )

      ny = 150

      do i = 0, 23
        if ( mod ( i, 4 ) .eq. 0 ) then
          ny = ny + 400
          nxp = 550
        else
          nxp = nxp + 350
        end if

        if ( i .lt. 10 ) then
          write ( cstr, '(i1)' ) i
        else
          write ( cstr, '(i2)' ) i
        end if

        nl = nlmess ( cstr ) / 2
        call messag ( cstr, nxp-nl, ny+150 )
        call symbol ( i, nxp, ny )

      end do
c
c  Close DISLIN.
c
      call disfin ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DISLIN_EX03:'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
