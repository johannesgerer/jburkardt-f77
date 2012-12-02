      program main

c*********************************************************************72
c
cc MAIN is the main program for ASA172_PRB.
c
c  Discussion:
c
c    ASA172_PRB calls the ASA172 routines.
c
c  Modified:
c
c    26 July 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA172_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Test routines in the ASA172 library.'

      call test01
      call test02

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA172_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01

c*********************************************************************72
c
cc TEST01 compares indices computed by a triple loop.
c
c  Modified:
c
c    26 July 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer kdim
      parameter ( kdim = 3 )

      integer i
      integer i1
      integer i2
      integer i3
      integer ifault
      integer iprod(kdim)
      integer ivec(kdim)
      integer j
      integer jsub
      integer n
      integer nr(kdim)
      logical qfor
      logical qind

      save nr

      data nr / 3, 2, 4 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) '  SIMDO can convert between compressed and'
      write ( *, '(a)' ) '  vector indices representing a nested loop.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Here, we set QFOR = FALSE, meaning we do'
      write ( *, '(a)' ) '  NOT want to convert from FORTRAN ordering'
      write ( *, '(a)' ) '  to lexical ordering.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Here, we actually carry out a triple loop'
      write ( *, '(a)' ) '  list the indices, and then compare.'

      qfor = .false.
c
c  If QFOR is FALSE, then the definition of IPROD is reversed...
c
      iprod(1) = nr(kdim)
      do i = 2, kdim
        iprod(i) = iprod(i-1) * nr(kdim+1-i)
      end do

      n = iprod(kdim)
c
c  Carry out the nested loops, and use JSUB to count each iteration.
c  In the inmost loop, print JSUB and the corresponding (I1,I2,I3) vector.
c
      jsub = 0

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  #1: Generate JSUB by counting as we DO the loops:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  DO I1 = 1, N1'
      write ( *, '(a)' ) '    DO I2 = 1, N2'
      write ( *, '(a)' ) '      DO I3 = 1, N3'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      JSUB            I1        I2        I3'
      write ( *, '(a)' ) ' '
      do i1 = 1, nr(1)
        ivec(1) = i1
        do i2 = 1, nr(2)
          ivec(2) = i2
          do i3 = 1, nr(3)
            ivec(3) = i3
            jsub = jsub + 1
            write ( *, '(2x,i8,6x,i8,2x,i8,2x,i8)' ) 
     &        jsub, i1, i2, i3
          end do
        end do
      end do
c
c  Now for each value of JSUB, retrieve the corresponding index subscript.
c  In order to use the QFOR = .FALSE. switch, I apparently have to reverse
c  the sense of the NR vector!
c
      qind = .true.

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  #2: Loop on JSUB, retrieve loop indices'
      write ( *, '(a)' ) '      QIND = TRUE J ->I(J)'
      write ( *, '(a)' ) '      QFOR = FALSE'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      JSUB            I1        I2        I3'
      write ( *, '(a)' ) ' '

      do j = 1, n
        jsub = j
        call simdo ( qind, qfor, iprod, kdim, jsub, ivec, ifault )
        write ( *, '(2x,i8,6x,i8,2x,i8,2x,i8)' ) 
     &    jsub, ( ivec(i), i = 1, kdim )
      end do
c
c  Carry out the nested loops, and DO NOT compute JSUB.
c  Have SIMDO determine JSUB.
c
      qind = .false.

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  #3: For any set of loop indices, retrieve JSUB'
      write ( *, '(a)' ) '      QIND = FALSE I(J) -> J'
      write ( *, '(a)' ) '      QFOR = FALSE'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      JSUB            I1        I2        I3'
      write ( *, '(a)' ) ' '
      do i1 = 1, nr(1)
        ivec(1) = i1
        do i2 = 1, nr(2)
          ivec(2) = i2
          do i3 = 1, nr(3)
            ivec(3) = i3
            call simdo ( qind, qfor, iprod, kdim, jsub, ivec, ifault )
            write ( *, '(2x,i8,6x,i8,2x,i8,2x,i8)' ) 
     &        jsub, i1, i2, i3
          end do
        end do
      end do

      return
      end
      subroutine test02

c*********************************************************************72
c
cc TEST02 compares indices computed by a triple loop.
c
c  Modified:
c
c    26 July 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer kdim
      parameter ( kdim = 3 )

      integer i
      integer i1
      integer i2
      integer i3
      integer ifault
      integer iprod(kdim)
      integer ivec(kdim)
      integer j
      integer jsub
      integer n
      integer nr(kdim)
      logical qfor
      logical qind

      save nr

      data nr / 3, 2, 4 /

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02:'
      write ( *, '(a)' ) '  SIMDO can convert between compressed and'
      write ( *, '(a)' ) '  vector indices representing a nested loop.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Here, we set QFOR = TRUE, meaning we DO'
      write ( *, '(a)' ) '  want to convert from the FORTRAN '
      write ( *, '(a)' ) '  ordering to lexical convention.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Here, we actually carry out a triple loop'
      write ( *, '(a)' ) '  list the indices, and then compare.'

      qfor = .true.

      iprod(1) = nr(1)
      do i = 2, kdim
        iprod(i) = iprod(i-1) * nr(i)
      end do

      n = iprod(kdim)
c
c  Carry out the nested loops, and use JSUB to count each iteration.
c  In the inmost loop, print JSUB and the corresponding (I1,I2,I3) vector.
c
      jsub = 0

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  #1: Generate JSUB by counting as we do the loops.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  DO I3 = 1, N3'
      write ( *, '(a)' ) '    DO I2 = 1, N2'
      write ( *, '(a)' ) '      DO I1 = 1, N1'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      JSUB            I1        I2        I3'
      write ( *, '(a)' ) ' '
      do i3 = 1, nr(3)
        ivec(3) = i3
        do i2 = 1, nr(2)
          ivec(2) = i2
          do i1 = 1, nr(1)
            ivec(1) = i1
            jsub = jsub + 1
            write ( *, '(2x,i8,6x,i8,2x,i8,2x,i8)' ) 
     &        jsub, i1, i2, i3
          end do
        end do
      end do
c
c  Reverse the order, so that the loop indices are generated in lexical order.
c
      qind = .true.

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  #2: Setting QFOR false means loop indices'
      write ( *, '(a)' ) '  are generated in lexical order.'
      write ( *, '(a)' ) '      QIND = TRUE J -> I(J)'
      write ( *, '(a)' ) '      QFOR = TRUE'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      JSUB            I1        I2        I3'
      write ( *, '(a)' ) ' '

      do j = 1, n
        jsub = j
        call simdo ( qind, qfor, iprod, kdim, jsub, ivec, ifault )
        write ( *, '(2x,i8,6x,i8,2x,i8,2x,i8)' ) 
     &    jsub, ( ivec(i), i = 1, kdim )
      end do
c
c  Carry out the nested loops, and DO NOT compute JSUB.
c  Have SIMDO determine JSUB.
c
      qind = .false.

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  #3: For any set of loop indices, retrieve JSUB'
      write ( *, '(a)' ) '      QIND = FALSE I(J) -> J'
      write ( *, '(a)' ) '      QFOR = TRUE'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '      JSUB            I1        I2        I3'
      write ( *, '(a)' ) ' '
      do i3 = 1, nr(3)
        ivec(3) = i3
        do i2 = 1, nr(2)
          ivec(2) = i2
          do i1 = 1, nr(1)
            ivec(1) = i1
            call simdo ( qind, qfor, iprod, kdim, jsub, ivec, ifault )
            write ( *, '(2x,i8,6x,i8,2x,i8,2x,i8)' ) 
     &        jsub, i1, i2, i3
          end do
        end do
      end do

      return
      end
