c  ANYNUL.FOR  20 September 1990  Version 1.01
c  ANYPLT dummy graphics interface
c
c  ANYPLT is a subroutine which provides a simple, standard interface
c  between FORTRAN programs and various output devices.  To run a
c  program which calls ANYPLT on a different machine, the program
c  is not modified in any way, but a different version of the ANYPLT
c  program is provided.  Currently, the following versions are available:
c
c  ANYATT - AT&T PC6300 graphics (640 by 400).  Requires ATTPLT.ASM.
c  ANYBUG - Simple debugging output to a file.
c  ANYCAL - CALCOMP file output.  Available on many mainframes.
c  ANYIBM - IBM PC hi resolution (640 by 200).  Requires IBMPLT.ASM.
c  ANYMAC - Macintosh graphics.  Requires auxilliary routine TOOLBX.SUB.
c  ANYNCR - NCAR graphics package.
c  ANYNUL - Does nothing.
c  ANYP10 - PLOT10 interactive graphics. (1024 by 768)
c  ANYTTY - Simple 'typewriter' graphics (80 by 24)
c
c  The personal computer routines generally require an auxilliary
c  set of routines written in assembly language.
c
c  ANYATT was written by
c
c  John Burkardt
c  Staff Programmer
c  Mathematics Department
c  University of Pittsburgh
c  Pittsburgh, PA
c
      subroutine anyplt(icom)
c
c***********************************************************************
c
cc ANYPLT
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
      character carray*80
      common /anycom/ iplt1,iplt2,ixplt1,ixplt2,iyplt1,
     &                iyplt2,marray,xplt1,xplt2,yplt1,yplt2
      common /anychr/ carray

      if(icom.eq.13)then
        carray='anynul - version 1.01  20 september 1990  null'
      end if

      return
      end
