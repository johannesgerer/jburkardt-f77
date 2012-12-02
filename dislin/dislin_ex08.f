      program main

c*********************************************************************72
c
cc DISLIN_EX08 demonstrates various shading patterns.
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
      character*60 ctit
      integer i
      integer iclr
      integer ii
      integer ix(4)
      integer ixp(4)
      integer iy(4)
      integer iyp(4)
      integer j
      integer k
      integer nl
      integer nlmess
      integer nx
      integer nx0
      integer ny
      integer ny0

      data ix / 0, 300, 300, 0 /
      data iy / 0, 0, 400, 400 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DISLIN_EX08:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Demonstrate the use of shading patterns.'

      ctit = 'Shading Patterns (AREAF)'
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
      call setfil ( 'dislin_ex08.png' )
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
      call setvlt ( 'SMALL' )

      call height ( 50 )
      nl = nlmess ( ctit )
      nx = ( 2970 - nl ) / 2
      call messag ( ctit, nx, 200 )

      nx0 = 335
      ny0 = 350

      do i = 1, 3

        ny = ny0 + ( i - 1 ) * 600

        do j = 1, 6

          iclr = ( i - 1 ) * 6 + j - 1
          iclr = mod ( iclr, 15 )

          if ( iclr .eq. 0 ) then
            iclr = 15
          end if

          call setclr ( iclr )

          nx = nx0 + ( j - 1 ) * 400
          ii = ( i - 1 ) * 6 + j - 1
          call shdpat ( ii )
          write ( cstr, '(i2)' ) ii

          do k = 1, 4
            ixp(k) = ix(k) + nx
            iyp(k) = iy(k) + ny
          end do

          call areaf ( ixp, iyp, 4 )

          nl = nlmess ( cstr )
          nx = nx + ( 300 - nl ) / 2
          call messag ( cstr, nx, ny+460 )

        end do

      end do
c
c  Close DISLIN.
c
      call disfin ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'DISLIN_EX08:'
      write ( *, '(a)' ) '  Normal end of execution.'

      stop
      end
