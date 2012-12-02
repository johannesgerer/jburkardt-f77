c  ANYBUG.F   11 October 1994
c
c  ANYPLT/BUGPLT interface
c
c  ANYBUG was written by
c
c  John Burkardt
c  Staff Programmer
c  Mathematics Department
c  University of Pittsburgh
c  Pittsburgh, PA
c
      subroutine anyplt(icom)

c***********************************************************************
c
cc ANYPLT is a generic graphics interface.
c
c  This version is an interface to a dummy debugging device.
c
c  BUGPLT graphics is actually a debugging device which writes
c  a file with extension BUG containing a list of the calls
c  to ANYPLT
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    24 January 2009
c
c  Author:
c
c    John Burkardt
c
      integer ldunit
      parameter (ldunit=3)
c
      character*80 carray
      integer icom
      integer iplt1
      integer iplt2
      integer ixplt1
      integer ixplt2
      integer iyplt1
      integer iyplt2
      integer marray
      real xplt1
      real xplt2
      real yplt1
      real yplt2
c
      common /anycom/ iplt1,iplt2,ixplt1,ixplt2,iyplt1,
     &                iyplt2,marray,xplt1,xplt2,yplt1,yplt2
      common /anychr/ carray
c
c
c  ICOM=0  Enable graphics
c
      if(icom.eq.0)then

        open(unit=ldunit,file='anyplt.bug',status='old',err=10)
        close(unit=ldunit,status='delete')

10      continue

        open(unit=ldunit,file='anyplt.bug',status='new')

        write(ldunit,*)'BugPlt - 00 - Enable graphics.'
        write(ldunit,*)'  Xmin=',xplt1,' Xmax=',xplt2
        write(ldunit,*)'  Ymin=',yplt1,' Ymax=',yplt2
c
c  ICOM=1  Disable graphics
c
      elseif(icom.eq.1)then

        write(ldunit,*)'BugPlt - 01 - Disable graphics.'
        close(unit=ldunit)
c
c  ICOM=2  Begin plot
c
      elseif(icom.eq.2)then
        write(ldunit,*)'BugPlt - 02 - Begin plot'
c
c  ICOM=3  Define plot size
c
      elseif(icom.eq.3)then
        write(ldunit,*)'BugPlt - 03 - Define plot limits:'
        write(ldunit,*)'  Xmin=',xplt1,' Xmax=',xplt2
        write(ldunit,*)'  Ymin=',yplt1,' Ymax=',yplt2
c
c  ICOM=4  Move to point
c
      elseif(icom.eq.4)then
        write(ldunit,*)'BugPlt - 04 - Move to the point:'
        write(ldunit,'(2x,2g14.6)')xplt1,yplt1
c
c  ICOM=5  Draw to point
c
      elseif(icom.eq.5)then
        write(ldunit,*)'BugPlt - 05 - Draw to the point:'
        write(ldunit,'(2x,2g14.6)')xplt1,yplt1
c
c  ICOM=6  Clear screen
c
      elseif(icom.eq.6)then
        write(ldunit,*)'BugPlt - 06 - Clear the screen.'
c
c  ICOM=7,  Write string at position
c
      elseif(icom.eq.7)then

        write(ldunit,*)'BugPlt - 07 - Write a string.'

        write(ldunit,*)'  At the point ',xplt1,yplt1
        write(ldunit,*)'  Write the ',marray,' characters:'
        write(ldunit,'(2x,a)')carray(1:marray)
c
c  ICOM=8  Use virtual cursor
c
      elseif(icom.eq.8)then
        write(ldunit,*)'BugPlt - 08 - Use the virtual cursor.'
c
c  ICOM=9  End plot
c
      elseif(icom.eq.9)then
        write(ldunit,*)'BugPlt - 09 - End this plot.'
c
c  ICOM=10  Ring bell
c
      elseif(icom.eq.10)then
        write(ldunit,*)'BugPlt - 10 - Ring the bell.'
c
c  ICOM=11  Mark data
c
      elseif(icom.eq.11)then

        write(ldunit,*)'BugPlt - 11 - Mark data point.'
        write(ldunit,*)'  At the point ',xplt1,yplt1
        write(ldunit,*)'  Mark the data with ',carray(1:1)
c
c  ICOM=12  Return screen data
c
      elseif(icom.eq.12)then
        write(ldunit,*)'BugPlt - 12 - Return maximum screen data.'
c
c  ICOM=13  Return version
c
      elseif(icom.eq.13)then
        write(ldunit,*)'Bugplt - 13 - Return version.'
        carray='AnyPlt - Version 1.02  11 October 1994  BugPlt'
        write(ldunit,'(a)')carray
c
c  Unknown value of ICOM.
c
      else
        write(*,*)'AnyPlt - Fatal error!'
        write(*,*)'  Unknown value of ICOM=',icom
        stop
      endif

      return
      end
 
