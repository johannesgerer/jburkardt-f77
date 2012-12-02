      program main

c*********************************************************************72
c
cc MAIN is the main program for HB_TO_ST.
c
c  Discussion:
c
c    HB_TO_ST converts a Harwell-Boeing file to sparse triplet format.
c
c    Read a sparse matrix in the Harwell/Boeing format and output a matrix
c    in zero-based sparse triplet format.  Only the lower triangular part of a
c    symmetric matrix is provided.  Does not handle skew-symmetric
c    matrices.
c
c  Usage:
c
c    hb_to_st < hb_file > st_file
c
c  Modified:
c
c    09 March 2006
c
c  Reference:
c
c    Timothy Davis,
c    Direct Methods for Sparse Linear Systems,
c    SIAM, Philadelphia, 2006.
c
      implicit none

      integer nmax
      integer nzmax
c
c  The values of NMAX and NZMAX specified here determine the maximum size
c  of the matrices that can be read by the program.  If the program fails
c  to read a matrix, and the error is caused by insufficient space, a message
c  will be printed suggesting the correct larger value of NMAX or NZMAX
c  that should be used here.
c
      parameter ( nmax =  10000 )
      parameter ( nzmax = 30000 )

      integer col
      integer indcrd
      integer index(nzmax)
      character*16 indfmt
      character*30 key
      integer n
      integer ncol
      integer nel
      integer nrhs
      integer nrow
      integer nz
      integer nzrhs
      integer p
      integer ptr(nmax)
      integer ptrcrd
      character*16 ptrfmt
      integer rhscrd
      character*20 rhsfmt
      character*3 rhstyp
      integer row
      double precision skew
      integer stype
      logical sym
      character*72 title
      integer totcrd
      character*3 type
      integer valcrd
      character*20 valfmt
      double precision value(nzmax)

      write ( 0, '(a)' ) ' '
      call timestamp ( )

      write ( 0, '(a)' ) ' '
      write ( 0, '(a)' ) 'HB_TO_ST:'
      write ( 0, '(a)' ) '  FORTRAN77 version'
      write ( 0, '(a)' ) '  Convert sparse matrix from HB to ST.'
c
c  Read header information.
c
      read ( 5, '(a72,a8)', err = 998 ) title, key

      read ( 5, '(5i14)', err = 998 ) 
     &  totcrd, ptrcrd, indcrd, valcrd, rhscrd

      read ( 5, '(a3,11x,4i14)', err = 998 ) 
     &  type, nrow, ncol, nz, nel

      read ( 5, '(2a16,2a20)', err = 998 ) 
     &  ptrfmt, indfmt, valfmt, rhsfmt

      if ( 0 .lt. rhscrd ) then
        read ( 5, '(a3,11x,i14,i14)', err = 998 ) 
     &    rhstyp, nrhs, nzrhs
      end if

      if ( type(2:2) .eq. 'Z' .or. type(2:2) .eq. 'z' ) then
        skew = -1.0D+00
      else if ( type(2:2) .eq. 'S' .or. type(2:2) .eq. 's' ) then
        skew =  1.0D+00
      else
        skew = 0.0D+00
      end if

      sym = skew .ne. 0.0D+00

      write ( 0, '(a)' ) ' '
      write ( 0, '(a,a)'   ) 'Title:  ', title
      write ( 0, '(a,a)'   ) 'Key:    ', key
      write ( 0, '(a,a)'   ) 'Type:   ', type
      write ( 0, '(a,i14)' ) 'NROW:   ', nrow
      write ( 0, '(a,i14)' ) 'NCOL:   ', ncol
      write ( 0, '(a,i14)' ) 'NZ:     ', nz

      if ( 0 .lt. rhscrd ) then
        write ( 0, '(a,a)'   ) 'RHSTYP: ', rhstyp
        write ( 0, '(a,i14)' ) 'NRHS:   ', nrhs
        write ( 0, '(a,i14)' ) 'NZRHS:  ', nzrhs
      end if

      write ( 0, '(a,l1)'  ) 'SYM:   ', sym
      write ( 0, '(a,g14.6)' ) 'SKEW:  ', skew

      if ( skew .eq. -1.0D+00 ) then
       write ( 0, '(a)' ) ' '
       write ( 0, '(a)' ) 'HB_TO_ST - Fatal error!'
       write ( 0, '(a)' ) '  Cannot handle skew-symmetric matrices.'
       stop
      end if

      n = max ( nrow, ncol )

      if ( nmax .le. ncol ) then
        write ( 0, '(a)' ) ' '
        write ( 0, '(a)' ) 'HB_TO_ST - Fatal error!'
        write ( 0, '(a)' ) '  The internal arrays are too small to'
        write ( 0, '(a)' ) '  handle this matrix.  To proceed,'
        write ( 0, '(a,i14)' ) '  increase NMAX to more than ', ncol
        stop
      end if

      if ( nzmax .lt. nz ) then
        write ( 0, '(a)' ) ' '
        write ( 0, '(a)' ) 'HB_TO_ST - Fatal error!'
        write ( 0, '(a)' ) '  The internal arrays are too small to'
        write ( 0, '(a)' ) '  handle this matrix.  To proceed,'
        write ( 0, '(a,i14)' ) '  increase NZMAX to at least ', nz
        stop
      end if
c
c  Read the pointers, indices, and possibly data values.
c
      read ( 5, ptrfmt, err = 998 ) ( ptr(p), p = 1, ncol+1 )
      read ( 5, indfmt, err = 998 ) ( index(p), p = 1, nz )

      if ( 0 .lt. valcrd ) then
        read ( 5, valfmt, err = 998 ) ( value(p), p = 1, nz )
      end if
c
c  Create the triplet form of the input matrix
c  STYPE = 0: unsymmetric
c  STYPE = -1: symmetric, lower triangular part present
c
      stype = -skew

      do col = 1, ncol

        do p = ptr(col), ptr(col+1) - 1

          row = index(p)

          if ( 0 .lt. valcrd ) then
            write ( 6, '(i8,i8,e30.18)' ) row-1, col-1, value(p)
            if ( sym .and. row .ne. col ) then
      	      write ( 6, '(i8,i8,e30.18)' ) 
     &           col-1, row-1, skew * value (p)
      	    end if
          else
            write ( 6, '(i8,i8,e30.18)' ) row-1, col-1, 1.0D+00
          end if

        end do

      end do

      write ( 0, '(a)' ) ' '
      write ( 0, '(a)' ) 'HB_TO_ST:'
      write ( 0, '(a)' ) '  Normal end of execution.'
      write ( 0, '(a)' ) ' '
      call timestamp ( )

      stop

998   continue
      write ( 0, '(a)' ) ' '
      write ( 0, '(a)' ) 'HB_TO_ST:'
      write ( 0, '(a)' ) '  I/O error reading the matrix.'
      write ( 0, '(a)' ) '  ABNORMAL end of execution.'
      stop
      end
      subroutine timestamp ( )

c*********************************************************************72
c
cc TIMESTAMP prints out the current YMDHMS date as a timestamp.
c
c  Discussion:
c
c    This FORTRAN77 version is made available for cases where the
c    FORTRAN90 version cannot be used.
c
c  Modified:
c
c    12 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    None
c
      implicit none

      character * ( 8 ) ampm
      integer d
      character * ( 8 ) date
      integer h
      integer m
      integer mm
      character * ( 9 ) month(12)
      integer n
      integer s
      character * ( 10 ) time
      integer y

      save month

      data month /
     &  'January  ', 'February ', 'March    ', 'April    ', 
     &  'May      ', 'June     ', 'July     ', 'August   ', 
     &  'September', 'October  ', 'November ', 'December ' /

      call date_and_time ( date, time )

      read ( date, '(i4,i2,i2)' ) y, m, d
      read ( time, '(i2,i2,i2,1x,i3)' ) h, n, s, mm

      if ( h .lt. 12 ) then
        ampm = 'AM'
      else if ( h .eq. 12 ) then
        if ( n .eq. 0 .and. s .eq. 0 ) then
          ampm = 'Noon'
        else
          ampm = 'PM'
        end if
      else
        h = h - 12
        if ( h .lt. 12 ) then
          ampm = 'PM'
        else if ( h .eq. 12 ) then
          if ( n .eq. 0 .and. s .eq. 0 ) then
            ampm = 'Midnight'
          else
            ampm = 'AM'
          end if
        end if
      end if

      write ( 0, 
     &  '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) 
     &  d, month(m), y, h, ':', n, ':', s, '.', mm, ampm

      return
      end
