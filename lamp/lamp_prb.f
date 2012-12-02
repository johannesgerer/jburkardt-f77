      program main

c*********************************************************************72
c
cc MAIN is the main program for LAMP_PRB.
c
c  Modified:
c
c    15 December 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LAMP_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the LAMP library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
      call test06 ( )
      call test07 ( )
      if ( .true. ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST08 is being skipped for now.'
        write ( *, '(a)' ) '  The QAP routine is being debugged.'
      else
        call test08 ( )
      end if
      call test09 ( )
      call test10 ( )
      call test11 ( )
      call test12 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LAMP_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests BMP.
c
c  Modified:
c
c    20 November 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 8 )

      integer basis(n)
      integer c(n,n)
      integer cc((n*(n-1))/2)
      integer i
      integer j
      integer mem(n)
      integer nmatch(n)
      integer nmatch_correct(n)
      integer p(n)
      integer sm(n)
      integer stack(n)
      integer sup
      integer tma(n)
      integer tmb(n)
      integer z
      integer z_correct

      data c /
     &   0, 33, 55, 46, 29, 68, 38, 37,
     &  33,  0, 57, 95, 46, 30, 38, 28,
     &  55, 57,  0, 43, 71, 60, 51, 42,
     &  46, 95, 43,  0, 20, 37, 14, 57,
     &  29, 46, 71, 20,  0, 48, 46, 93,
     &  68, 30, 60, 37, 48,  0, 68, 77,
     &  38, 38, 51, 14, 46, 68,  0, 61,
     &  37, 28, 42, 57, 93, 77, 61,  0 /
c
c  Here are the entries of the strict upper triangle,
c  listed as a vector, columnwise.
c
      data cc /
     & 33,
     & 55, 57,
     & 46, 95, 43,
     & 29, 46, 71, 20,
     & 68, 30, 60, 37, 48,
     & 38, 38, 51, 14, 46, 68,
     & 37, 28, 42, 57, 93, 77, 61 /

      data nmatch_correct /
     &  5, 6, 8, 7, 1, 2, 4, 3 /

      z_correct = 42

      sup = 100000000

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  BMP solves the'
      write ( *, '(a)' ) '  Bottleneck Matching Problem.'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Cost matrix:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(8i4)' ) ( c(i,j), j = 1, n )
      end do

      call bmp ( n, cc, sup, nmatch, z, basis, mem, stack, sm, tma,
     &  tmb, p )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Optimal Matching:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I      K           K'
      write ( *, '(a)' ) '        (computed)  (correct)'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i4,4x,i8,4x,i8)' )
     &  i, nmatch(i), nmatch_correct(I)
      end do
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Optimal cost (computed) = ', z
      write ( *, '(a,i8)' ) '  Optimal cost (correct)  = ', z_correct

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests CMP.
c
c  Modified:
c
c    22 November 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 46 )
      integer n
      parameter ( n = 27 )

      integer basis(n)
      integer i
      integer index(n+1)
      integer mem(n)
      integer nbl(2*m)
      integer ncard
      integer ncard_correct
      integer nmatch(n)
      integer nmatch_correct(n)
      integer sm(n)
      integer stack(n)
      integer tma(n)
      integer tmb(n)

      data index /
     &  1,  3, 10, 14, 18, 21, 28, 30, 33, 36,
     & 39, 41, 43, 46, 51, 55, 57, 60, 63, 64,
     & 67, 73, 76, 80, 83, 87, 90, 93 /

      data nbl /
     &   2, 12,
     &   1,  3, 15, 17, 19, 24, 27,
     &   2,  4,  6,  8,
     &   3,  5,  6, 10,
     &   4,  6,  9,
     &   3,  4,  5,  7,  8,  9, 10,
     &   6, 11,
     &   3,  6,  9,
     &   5,  6,  8,
     &   4,  6, 11,
     &   7, 10,
     &   1, 13,
     &  12, 14, 15,
     &  13, 16, 18, 20, 21,
     &   2, 13, 17, 20,
     &  14, 17,
     &   2, 15, 16,
     &  14, 21, 25,
     &   2,
     &  14, 15, 21,
     &  14, 18, 20, 22, 23, 25,
     &  21, 23, 25,
     &  21, 22, 24, 26,
     &   2, 23, 27,
     &  18, 21, 22, 26,
     &  23, 25, 27,
     &   2, 24, 26 /

      data nmatch_correct /
     &  12, 19,  8, 10,  9,  0, 11, 3,  5,  4,
     &   7,  1, 15, 20, 13, 17, 16, 21, 2, 14,
     &  18, 25, 24, 23, 22, 27, 26 /

      ncard_correct = 13

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  CMP solves the'
      write ( *, '(a)' ) '  Cardinality matching problem.'

      call cmp ( n, m, nbl, index, nmatch, ncard, basis, mem, stack,
     &  sm, tma, tmb )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Optimal Matching:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I      K           K'
      write ( *, '(a)' ) '        (computed)  (correct)'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i4,4x,i8,4x,i8)' )
     &  i, nmatch(i), nmatch_correct(I)
      end do
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Cardinality (computed) = ', ncard
      write ( *, '(a,i8)' ) '  Cardinality (correct)  = ', ncard_correct

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests CONNECT.
c
c  Modified:
c
c    06 December 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 78 )
      integer n
      parameter ( n = 52 )

      integer basis(n)
      integer icon
      integer icon_correct
      integer index(n+1)
      integer mem(n)
      integer nb(2*m)

      data index /
     &   1,   4,   8,  12,  16,  19,  22,  26,  29,  32,
     &  36,  40,  43,  46,  49,  52,  55,  58,  61,  65,
     &  66,  69,  72,  75,  78,  81,  84,  87,  90,  91,
     &  94,  96, 100, 103, 107, 110, 113, 117, 118, 121,
     & 125, 126, 129, 132, 135, 138, 141, 144, 147, 150,
     & 153, 156, 157 /

      data nb /
     &    2,  3,  4,
     &    1, 17, 28, 29,
     &    1,  5, 21, 52,
     &    1,  5,  7, 17,
     &    3,  4,  6,
     &    5,  8,  9,
     &    4,  8, 17, 18,
     &    6,  7, 18,
     &    6, 10, 19,
     &    9, 11, 18, 20,
     &   10, 12, 19, 42,
     &   11, 13, 23,
     &   12, 14, 37,
     &   13, 15, 39,
     &   14, 16, 24,
     &   15, 25, 26,
     &    2,  4,  7,
     &    7,  8, 10,
     &    9, 11, 22, 23,
     &   10,
     &    3, 22, 26,
     &   19, 21, 27,
     &   12, 19, 24,
     &   15, 23, 25,
     &   16, 24, 27,
     &   16, 21, 27,
     &   22, 25, 26,
     &    2, 47, 50,
     &    2,
     &   31, 49, 50,
     &   30, 50,
     &   33, 34, 40, 43,
     &   32, 34, 51,
     &   32, 33, 35, 37,
     &   34, 36, 39,
     &   35, 39, 51,
     &   13, 34, 38, 40,
     &   37,
     &   14, 35, 36,
     &   32, 37, 41, 42,
     &   40,
     &   11, 40, 43,
     &   32, 42, 44,
     &   43, 45, 48,
     &   44, 46, 48,
     &   45, 47, 51,
     &   28, 46, 49,
     &   44, 45, 49,
     &   30, 47, 48,
     &   28, 30, 31,
     &   33, 36, 46,
     &   3 /

      icon_correct = 1

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  CONNECT checks that a graph is connected.'

      call connect ( n, nb, index, basis, mem, icon )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  ICON (computed) = ', icon
      write ( *, '(a,i8)' ) '  ICON (correct)  = ', icon_correct

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests CPP.
c
c  Modified:
c
c    10 December 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer m
      parameter ( m = 78 )
      integer n
      parameter ( n = 52 )

      integer basis(n)
      real dminus(n)
      real dplus(n)
      real eps
      integer i
      integer icon
      integer icon_correct
      integer index(n+1)
      integer j
      integer jhi
      integer jlo
      integer ka(n)
      integer kb(n)
      integer kost(2*m)
      integer kosten
      integer kosten_correct
      integer kst
      integer kst_correct
      integer kurs
      integer mem(n)
      integer nb(2*m)
      integer sm(n)
      integer tma(n)
      integer tmb(n)
      integer top
      real y1(n)
      real y2(n)

      data index /
     &   1,   4,   8,  12,  16,  19,  22,  26,  29,  32,
     &  36,  40,  43,  46,  49,  52,  55,  58,  61,  65,
     &  66,  69,  72,  75,  78,  81,  84,  87,  90,  91,
     &  94,  96, 100, 103, 107, 110, 113, 117, 118, 121,
     & 125, 126, 129, 132, 135, 138, 141, 144, 147, 150,
     & 153, 156, 157 /

      data kost /
     &      26,   52,   18,
     &      26,   17,  141,   14,
     &      52,   20,  161, 1605,
     &      18,   41,   32,   26,
     &      20,   41,   37,
     &      37,   39,   43,
     &      32,   13,   43,   82,
     &      39,   13,    7,
     &      43,   47,   47,
     &      47,   44,   37,   27,
     &      44,   23,   36,   35,
     &      23,   14,   35,
     &      14,   57,   37,
     &      57,   34,   54,
     &      34,  101,   55,
     &     101,   25,   59,
     &      17,   26,   43,
     &      82,    7,   37,
     &      47,   36,   16,   26,
     &      27,
     &     161,   47,   12,
     &      16,   47,   16,
     &      35,   26,   22,
     &      55,   22,   24,
     &      25,   24,   40,
     &      59,   12,   25,
     &      16,   40,   25,
     &     141,   19,   12,
     &      14,
     &      42,    8,   18,
     &      42,   25,
     &      62,   18,   31,   17,
     &      62,   48,   31,
     &      18,   48,   23,   44,
     &      23,   45,   46,
     &      45,   72,   24,
     &      37,   44,   17,   18,
     &      17,
     &      54,   46,   72,
     &      31,   18,   28,   19,
     &      28,
     &      35,   19,   42,
     &      17,   42,   18,
     &      18,   27,   29,
     &      27,   48,    9,
     &      48,  113,   32,
     &      19,  113,   38,
     &      29,    9,   17,
     &       8,   38,   17,
     &      12,   18,   25,
     &      31,   24,   32,
     &    1605 /

      data nb /
     &    2,  3,  4,
     &    1, 17, 28, 29,
     &    1,  5, 21, 52,
     &    1,  5,  7, 17,
     &    3,  4,  6,
     &    5,  8,  9,
     &    4,  8, 17, 18,
     &    6,  7, 18,
     &    6, 10, 19,
     &    9, 11, 18, 20,
     &   10, 12, 19, 42,
     &   11, 13, 23,
     &   12, 14, 37,
     &   13, 15, 39,
     &   14, 16, 24,
     &   15, 25, 26,
     &    2,  4,  7,
     &    7,  8, 10,
     &    9, 11, 22, 23,
     &   10,
     &    3, 22, 26,
     &   19, 21, 27,
     &   12, 19, 24,
     &   15, 23, 25,
     &   16, 24, 27,
     &   16, 21, 27,
     &   22, 25, 26,
     &    2, 47, 50,
     &    2,
     &   31, 49, 50,
     &   30, 50,
     &   33, 34, 40, 43,
     &   32, 34, 51,
     &   32, 33, 35, 37,
     &   34, 36, 39,
     &   35, 39, 51,
     &   13, 34, 38, 40,
     &   37,
     &   14, 35, 36,
     &   32, 37, 41, 42,
     &   40,
     &   11, 40, 43,
     &   32, 42, 44,
     &   43, 45, 48,
     &   44, 46, 48,
     &   45, 47, 51,
     &   28, 46, 49,
     &   44, 45, 49,
     &   30, 47, 48,
     &   28, 30, 31,
     &   33, 36, 46,
     &   3 /

      icon_correct = 1
      kst_correct = 2248
      kosten_correct = 6700

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  CPP solves the Chinese Postman Problem.'
!
!  Verify that the graph is connected.
!
      call connect ( n, nb, index, basis, mem, icon )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Call CONNECT to verify that the graph'
      write ( *, '(a)' ) '  is connected.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  ICON (computed) = ', icon
      write ( *, '(a,i8)' ) '  ICON (correct)  = ', icon_correct

      if ( icon .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TEST04'
        write ( *, '(a)' ) '  The graph is not connected.'
        write ( *, '(a)' ) '  We cannot call CPP.'
        return
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  List of edges and associated costs:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        jlo = index(i)
        jhi = index(i+1) - 1
        do j = jlo, jhi
          write ( *, '(2x,i4,2x,i4,2x,i4)' ) i, nb(j), kost(j)
        end do
      end do

      kosten = 0
      do i = 1, 2 * m
        kosten = kosten + kost(i)
      end do

      top = 100000000
      kurs = 52
      eps = 0.000001E+00

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  The tour starts at node ', kurs

      call cpp ( n, m, kst, nb, kost, index, top, basis, mem, ka, kb,
     &  sm, tma, tmb, y1, y2, dplus, dminus, kurs, eps )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Next Node characterization of the tour:'
      write ( *, '(a)' ) ' '

      do i = 1, n
        write ( *, '(2x,i4,10(2x,i4))' )
     &  i, ( kost(j), j = index(i), kb(i) )
      end do

      kosten = ( kosten / 2 ) + kst

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' )
     &  '  Cost of the duplicated edges (computed) = ', kst
      write ( *, '(a,i8)' )
     &  '  Cost of the duplicated edges (correct)  = ', kst_correct
      write ( *, '(a,i8)' )
     &  '  Cost of the postman tour     (computed) = ', kosten
      write ( *, '(a,i8)' )
     &  '  Cost of the postman tour     (correct)  = ', kosten_correct

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests LBAP.
c
c  Modified:
c
c    20 November 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 8 )

      integer c(n,n)
      integer dminus(n)
      integer dplus(n)
      integer i
      integer j
      logical label(n)
      integer spalte(n)
      integer spalte_correct(n)
      integer sup
      integer vor(n)
      integer vos(n)
      integer z
      integer z_correct
      integer zeile(n)

      data c /
     &  13, 12, 35, 34, 21, 42, 16, 26,
     &  21, 36, 32, 54,  6, 19, 34, 20,
     &  20, 25, 13,  7, 45, 39, 38,  5,
     &  12, 41, 36,  8, 18, 15,  3, 17,
     &   8, 40, 26, 12, 24, 14, 34, 45,
     &  26, 11, 21, 22, 34, 16, 40, 31,
     &  22,  4, 13, 11, 12, 28, 22, 37,
     &  11,  8, 37, 40, 48, 46, 24, 43 /

      data spalte_correct /
     & 1, 8, 7, 5, 2, 6, 4, 3 /

      z_correct = 16

      sup = 100000000
!
!  The matrix is the transpose of what I think!
!
      do i = 1, n
        do j = 1, i - 1
          z = c(i,j)
          c(i,j) = c(j,i)
          c(j,i) = z
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  LBAP solves the'
      write ( *, '(a)' ) '  Linear bottleneck assignment problem'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix order N = ', n
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Cost matrix C:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(8(2x,i4))' ) ( c(i,j), j = 1, n )
      end do

      call lbap ( n, sup, c, z, zeile, spalte, dminus, dplus, vor,
     &  vos, label )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Optimal assignment:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I      K           K'
      write ( *, '(a)' ) '        (computed)  (correct)'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i4,4x,i8,4x,i8)' )
     &  i, spalte(i), spalte_correct(I)
      end do
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Optimal cost (computed) = ', z
      write ( *, '(a,i8)' ) '  Optimal cost (correct)  = ', z_correct

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 tests LSAPI.
c
c  Modified:
c
c    20 November 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 8 )

      integer c(n,n)
      integer dminus(n)
      integer dplus(n)
      integer i
      integer j
      logical label(n)
      integer spalte(n)
      integer spalte_correct(n)
      integer sup
      integer vor(n)
      integer ys(n)
      integer yt(n)
      integer z
      integer z_correct
      integer zeile(n)

      data c /
     &  13, 12, 35, 34, 21, 42, 16, 26,
     &  21, 36, 32, 54,  6, 19, 34, 20,
     &  20, 25, 13,  7, 45, 39, 38,  5,
     &  12, 41, 36,  8, 18, 15,  3, 17,
     &   8, 40, 26, 12, 24, 14, 34, 45,
     &  26, 11, 21, 22, 34, 16, 40, 31,
     &  22,  4, 13, 11, 12, 28, 22, 37,
     &  11,  8, 37, 40, 48, 46, 24, 43 /

      data spalte_correct /
     & 1, 8, 7, 5, 2, 6, 4, 3 /

      z_correct = 76

      sup = 100000000
!
!  The matrix is the transpose of what I think!
!
      do i = 1, n
        do j = 1, i - 1
          z = c(i,j)
          c(i,j) = c(j,i)
          c(j,i) = z
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) '  LSAPI solves the '
      write ( *, '(a)' ) '  linear bottleneck assignment problem'
      write ( *, '(a)' ) '  for an INTEGER cost matrix.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix order N = ', n
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Cost matrix C:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(8(2x,i4))' ) ( c(i,j), j = 1, n )
      end do

      call lsapi ( n, sup, c, z, zeile, spalte, dminus, dplus, ys,
     &  yt, vor, label )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Optimal assignment:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I      K           K'
      write ( *, '(a)' ) '        (computed)  (correct)'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i4,4x,i8,4x,i8)' )
     &  i, spalte(i), spalte_correct(I)
      end do
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Optimal cost (computed) = ', z
      write ( *, '(a,i8)' ) '  Optimal cost (correct)  = ', z_correct

      return
      end
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07 tests LSAPR.
c
c  Modified:
c
c    22 November 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 8 )

      real c(n,n)
      real dminus(n)
      real dplus(n)
      real eps
      integer i
      integer j
      logical label(n)
      integer spalte(n)
      integer spalte_correct(n)
      integer sup
      integer vor(n)
      real ys(n)
      real yt(n)
      real z
      real z_correct
      integer zeile(n)

      data c /
     &  13.0, 12.0, 35.0, 34.0, 21.0, 42.0, 16.0, 26.0,
     &  21.0, 36.0, 32.0, 54.0,  6.0, 19.0, 34.0, 20.0,
     &  20.0, 25.0, 13.0,  7.0, 45.0, 39.0, 38.0,  5.0,
     &  12.0, 41.0, 36.0,  8.0, 18.0, 15.0,  3.0, 17.0,
     &   8.0, 40.0, 26.0, 12.0, 24.0, 14.0, 34.0, 45.0,
     &  26.0, 11.0, 21.0, 22.0, 34.0, 16.0, 40.0, 31.0,
     &  22.0,  4.0, 13.0, 11.0, 12.0, 28.0, 22.0, 37.0,
     &  11.0,  8.0, 37.0, 40.0, 48.0, 46.0, 24.0, 43.0 /

      data spalte_correct /
     & 1, 8, 7, 5, 2, 6, 4, 3 /

      eps = 0.000001E+00
      z_correct = 76.0E+00

      sup = 100000000
!
!  The matrix is the transpose of what I think!
!
      do i = 1, n
        do j = 1, i - 1
          z = c(i,j)
          c(i,j) = c(j,i)
          c(j,i) = z
        end do
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07'
      write ( *, '(a)' ) '  LSAPR solves the '
      write ( *, '(a)' ) '  linear bottleneck assignment problem'
      write ( *, '(a)' ) '  for a REAL cost matrix.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix order N = ', n
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Cost matrix C:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(8(2x,f6.2))' ) ( c(i,j), j = 1, n )
      end do

      call lsapr ( n, sup, c, z, zeile, spalte, dminus, dplus, ys,
     &  yt, vor, label, eps )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Optimal assignment:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I      K           K'
      write ( *, '(a)' ) '        (computed)  (correct)'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i4,4x,i8,4x,i8)' )
     &  i, spalte(i), spalte_correct(I)
      end do
      write ( *, '(a)' ) ' '
      write ( *, '(a,f8.2)' ) '  Optimal cost (computed) = ', z
      write ( *, '(a,f8.2)' ) '  Optimal cost (correct)  = ', z_correct

      return
      end
      subroutine test08 ( )

c*********************************************************************72
c
cc TEST08 tests QAP.
c
c  Modified:
c
c    14 December 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 8 )

      integer a(n,n)
      integer alter(n-2)
      integer aspei((n*(n+1)*(n-1))/3)
      integer b(n,n)
      logical bool(n)
      logical bool1(n)
      integer bspei((n*(n+1)*(n-1))/3)
      integer c(n,n)
      integer cspei((n*(n+1)*(2*n+1))/6-1)
      integer dd(n)
      integer h1(n)
      logical hl1(n)
      integer i
      integer j
      integer lab(n)
      integer loesg(n)
      integer loesg_correct(n)
      integer meng(n-2)
      integer menge(n-2)
      integer olwert
      integer olwert_correct
      integer partpe(n)
      integer phiofm(n-2)
      integer u(n)
      integer umspei(n,n)
      integer unendl
      integer v(n)
      integer vekqu(n-2)
      integer veksum(n-2)
      integer y(n)
      integer z1(n)
      integer zul(n,n)

      data a /
     &  0, 1, 2, 3, 1, 2, 3, 4,
     &  1, 0, 1, 2, 2, 1, 2, 3,
     &  2, 1, 0, 1, 3, 2, 1, 2,
     &  3, 2, 1, 0, 4, 3, 2, 1,
     &  1, 2, 3, 4, 0, 1, 2, 3,
     &  2, 1, 2, 3, 1, 0, 1, 2,
     &  3, 2, 1, 2, 2, 1, 0, 1,
     &  4, 3, 2, 1, 3, 2, 1, 0 /

      data b /
     &  0,  5,  2,  4,  1,  0,  0,  6,
     &  5,  0,  3,  0,  2,  2,  2,  0,
     &  2,  3,  0,  0,  0,  0,  0,  5,
     &  4,  0,  0,  0,  5,  2,  2, 10,
     &  1,  2,  0,  5,  0, 10,  0,  0,
     &  0,  2,  0,  2, 10,  0,  5,  1,
     &  0,  2,  0,  2,  0,  5,  0, 10,
     &  6,  0,  5, 10,  0,  1, 10,  0 /

      data loesg_correct /
     & 3,  8,  7,  6,  2,  1,  4,  5 /

      unendl = 1000000000
      olwert_correct = 214

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08'
      write ( *, '(a)' )
     &  '  QAP solves the Quadratic Assignment Problem.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix order N = ', n
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Distance matrix A:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(12(2x,i4))' ) ( a(i,j), j = 1, n )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Connection matrix B:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(12(2x,i4))' ) ( b(i,j), j = 1, n )
      end do

      call qap ( n, a, b, unendl, loesg, olwert, c, umspei, zul,
     &  u, v, dd, partpe, y, lab, z1, h1, meng, phiofm, menge, veksum,
     &  vekqu, alter, aspei, bspei, cspei, bool, bool1, hl1 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Optimal permutation:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I           P          P'
      write ( *, '(a)' ) '        (computed)  (example)'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i4,4x,i8,4x,i8)' )
     &  i, loesg(i), loesg_correct(I)
      end do
      write ( *, '(a)' ) ' '
      write ( *, '(a,I8)' )
     &  '  Objective functional value (computed) = ', olwert
      write ( *, '(a,I8)' )
     &  '  Objective functional value (example)  = ', olwert_correct

      return
      end
      subroutine test09 ( )

c*********************************************************************72
c
cc TEST09 tests QAPH1.
c
c  Modified:
c
c    27 November 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 12 )

      integer a(n,n)
      integer b(n,n)
      logical bool(n)
      logical bool1(n)
      integer i
      integer iperm(n)
      integer j
      integer menge(n)
      integer olwert
      integer olwert_correct
      integer perm(n)
      integer perm_correct(n)
      integer phiofm(n)
      integer rep
      integer seed
      integer unendl

      data a /
     &  0, 1, 2, 3, 1, 2, 3, 4, 2, 3, 4, 5,
     &  1, 0, 1, 2, 2, 1, 2, 3, 3, 2, 3, 4,
     &  2, 1, 0, 1, 3, 2, 1, 2, 4, 3, 2, 3,
     &  3, 2, 1, 0, 4, 3, 2, 1, 5, 4, 3, 2,
     &  1, 2, 3, 4, 0, 1, 2, 3, 1, 2, 3, 4,
     &  2, 1, 2, 3, 1, 0, 1, 2, 2, 1, 2, 3,
     &  3, 2, 1, 2, 2, 1, 0, 1, 3, 2, 1, 2,
     &  4, 3, 2, 1, 3, 2, 1, 0, 4, 3, 2, 1,
     &  2, 3, 4, 5, 1, 2, 3, 4, 0, 1, 2, 3,
     &  3, 2, 3, 4, 2, 1, 2, 3, 1, 0, 1, 2,
     &  4, 3, 2, 3, 3, 2, 1, 2, 2, 1, 0, 1,
     &  5, 4, 3, 2, 4, 3, 2, 1, 3, 2, 1, 0 /

      data b /
     &  0,  5,  2,  4,  1,  0,  0,  6,  2,  1,  1,  1,
     &  5,  0,  3,  0,  2,  2,  2,  0,  4,  5,  0,  0,
     &  2,  3,  0,  0,  0,  0,  0,  5,  5,  2,  2,  2,
     &  4,  0,  0,  0,  5,  2,  2, 10,  0,  0,  5,  5,
     &  1,  2,  0,  5,  0, 10,  0,  0,  0,  5,  1,  1,
     &  0,  2,  0,  2, 10,  0,  5,  1,  1,  5,  4,  0,
     &  0,  2,  0,  2,  0,  5,  0, 10,  5,  2,  3,  3,
     &  6,  0,  5, 10,  0,  1, 10,  0,  0,  0,  5,  0,
     &  2,  4,  5,  0,  0,  1,  5,  0,  0,  0, 10, 10,
     &  1,  5,  2,  0,  5,  5,  2,  0,  0,  0,  5,  0,
     &  1,  0,  2,  5,  1,  4,  3,  5, 10,  5,  0,  2,
     &  1,  0,  2,  5,  1,  0,  3,  0, 10,  0,  2,  0 /

      data perm_correct /
     & 1, 10,  5,  6,  2, 11,  4,  8,  3,  9, 12,  7 /

      unendl = 1000000000
      olwert_correct = 620
      rep = 15
      seed = 123456789

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST09'
      write ( *, '(a)' ) '  QAPH1 is a heuristic solver for the '
      write ( *, '(a)' ) '  the Quadratic Assignment Problem.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Instead of a "correct" solution, we have'
      write ( *, '(a)' ) '  an "example" solution computed with a'
      write ( *, '(a)' ) '  small number of repetitions.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix order N = ', n
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Distance matrix A:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(12(2x,i4))' ) ( a(i,j), j = 1, n )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Connection matrix B:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(12(2x,i4))' ) ( b(i,j), j = 1, n )
      end do

      call qaph1 ( n, a, b, unendl, rep, perm, olwert, menge, phiofm,
     &  iperm, bool, bool1, seed )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Optimal permutation:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I           P          P'
      write ( *, '(a)' ) '        (computed)  (example)'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i4,4x,i8,4x,i8)' )
     &  i, perm(i), perm_correct(I)
      end do
      write ( *, '(a)' ) ' '
      write ( *, '(a,I8)' )
     &  '  Objective functional value (computed) = ', olwert
      write ( *, '(a,I8)' )
     &  '  Objective functional value (example)  = ', olwert_correct

      return
      end
      subroutine test10 ( )

c*********************************************************************72
c
cc TEST10 tests QAPH2.
c
c  Modified:
c
c    02 December 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 12 )

      integer miter
      parameter ( miter = 3 * n )

      integer a(n,n)
      integer ap(n)
      real areal(n,n)
      integer b(n,n)
      logical bool(n)
      integer c(n,n)
      real dd(n)
      real eps
      integer h1(n,n)
      integer h2(n)
      integer h3(n)
      integer h4(n)
      real hr(n)
      integer i
      integer j
      real lambda(n,n)
      integer ol
      integer ol_correct
      integer op(n)
      integer ope(n)
      integer ope_correct(n)
      integer rep
      integer rest(miter,n,n)
      integer rest0(miter)
      integer seed
      logical startp
      real u(n)
      integer unendl
      real v(n)

      data a /
     &  0, 1, 2, 3, 1, 2, 3, 4, 2, 3, 4, 5,
     &  1, 0, 1, 2, 2, 1, 2, 3, 3, 2, 3, 4,
     &  2, 1, 0, 1, 3, 2, 1, 2, 4, 3, 2, 3,
     &  3, 2, 1, 0, 4, 3, 2, 1, 5, 4, 3, 2,
     &  1, 2, 3, 4, 0, 1, 2, 3, 1, 2, 3, 4,
     &  2, 1, 2, 3, 1, 0, 1, 2, 2, 1, 2, 3,
     &  3, 2, 1, 2, 2, 1, 0, 1, 3, 2, 1, 2,
     &  4, 3, 2, 1, 3, 2, 1, 0, 4, 3, 2, 1,
     &  2, 3, 4, 5, 1, 2, 3, 4, 0, 1, 2, 3,
     &  3, 2, 3, 4, 2, 1, 2, 3, 1, 0, 1, 2,
     &  4, 3, 2, 3, 3, 2, 1, 2, 2, 1, 0, 1,
     &  5, 4, 3, 2, 4, 3, 2, 1, 3, 2, 1, 0 /

      data b /
     &  0,  5,  2,  4,  1,  0,  0,  6,  2,  1,  1,  1,
     &  5,  0,  3,  0,  2,  2,  2,  0,  4,  5,  0,  0,
     &  2,  3,  0,  0,  0,  0,  0,  5,  5,  2,  2,  2,
     &  4,  0,  0,  0,  5,  2,  2, 10,  0,  0,  5,  5,
     &  1,  2,  0,  5,  0, 10,  0,  0,  0,  5,  1,  1,
     &  0,  2,  0,  2, 10,  0,  5,  1,  1,  5,  4,  0,
     &  0,  2,  0,  2,  0,  5,  0, 10,  5,  2,  3,  3,
     &  6,  0,  5, 10,  0,  1, 10,  0,  0,  0,  5,  0,
     &  2,  4,  5,  0,  0,  1,  5,  0,  0,  0, 10, 10,
     &  1,  5,  2,  0,  5,  5,  2,  0,  0,  0,  5,  0,
     &  1,  0,  2,  5,  1,  4,  3,  5, 10,  5,  0,  2,
     &  1,  0,  2,  5,  1,  0,  3,  0, 10,  0,  2,  0 /

      data ope_correct /
     & 5, 6, 10, 2, 4, 8, 11, 1, 12, 7, 9, 3 /

      eps = 0.00001E+00
      ol_correct = 578
      rep = 4
      seed = 123456789
      startp = .false.
      unendl = 1000000000

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST10'
      write ( *, '(a)' ) '  QAPH2 is a heuristic solver for the '
      write ( *, '(a)' ) '  the Quadratic Assignment Problem.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Instead of a "correct" solution, we have'
      write ( *, '(a)' ) '  an "example" solution computed with a'
      write ( *, '(a)' ) '  small number of repetitions.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) '  Matrix order N = ', n
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Distance matrix A:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(12(2x,i4))' ) ( a(i,j), j = 1, n )
      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Connection matrix B:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(12(2x,i4))' ) ( b(i,j), j = 1, n )
      end do

      call qaph2 ( n, a, b, miter, rep, unendl, eps, ope, startp,
     &  ol, c, rest0, rest, h1, h2, h3, h4, ap, op, areal, lambda,
     &  u, v, dd, hr, bool, seed )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Optimal permutation:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I           P          P'
      write ( *, '(a)' ) '        (computed)  (example)'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i4,4x,i8,4x,i8)' )
     &  i, ope(i), ope_correct(I)
      end do
      write ( *, '(a)' ) ' '
      write ( *, '(a,I8)' )
     &  '  Objective functional value (computed) = ', ol
      write ( *, '(a,I8)' )
     &  '  Objective functional value (example)  = ', ol_correct

      return
      end
      subroutine test11 ( )

c*********************************************************************72
c
cc TEST11 tests SMP.
c
c  Modified:
c
c    06 December 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 8 )

      integer basis(n)
      integer c(n,n)
      integer cc(((n-1)*n)/2)
      real dminus(n)
      real dplus(n)
      real eps
      integer i
      integer j
      integer k
      integer ka(n)
      integer kb(n)
      integer m1(n)
      integer max
      integer mem(n)
      integer nmatch(n)
      integer nmatch_correct(n)
      integer p(n)
      integer sm(n)
      real sup
      integer tma(n)
      integer tmb(n)
      real y1(n)
      real y2(n)
      integer zfw
      integer zfw_correct

      data c /
     &  0, 33, 55, 46, 29, 68, 38, 37,
     & 33,  0, 57, 95, 46, 30, 38, 28,
     & 55, 57,  0, 43, 71, 60, 51, 42,
     & 46, 95, 43,  0, 20, 37, 14, 57,
     & 29, 46, 71, 20,  0, 48, 46, 93,
     & 68, 30, 60, 37, 48,  0, 68, 77,
     & 38, 38, 51, 14, 46, 68,  0, 61,
     & 37, 28, 42, 57, 93, 77, 61,  0 /

      data nmatch_correct /
     &  5, 6, 8, 7, 1, 2, 4, 3 /

      eps = 0.000001E+00

      max = 0
      do i = 1, n
        do j = 1, n
          if ( max .lt. c(i,j) ) then
            max = c(i,j)
          end if
        end do
      end do

      sup = float ( max )
      zfw_correct = 115

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST11'
      write ( *, '(a)' ) '  SMP solves the Sum Matching Problem. '

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Cost matrix C:'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(12(2x,i4))' ) ( c(i,j), j = 1, n )
      end do

      k = 0
      do j = 2, n
        do i = 1, j - 1
          k = k + 1
          cc(k) = c(i,j)
        end do
      end do

      call smp ( n, zfw, nmatch, cc, p, basis, mem, ka, kb, sm, tma,
     &  tmb, m1, y1, y2, dplus, dminus, sup, eps )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Optimal matching:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '     I           P          P'
      write ( *, '(a)' ) '        (computed)  (correct)'
      write ( *, '(a)' ) ' '
      do i = 1, n
        write ( *, '(2x,i4,4x,i8,4x,i8)' )
     &  i, nmatch(i), nmatch_correct(I)
      end do
      write ( *, '(a)' ) ' '
      write ( *, '(a,I8)' )
     &  '  Cost of optimal matching (computed) = ', zfw
      write ( *, '(a,I8)' )
     &  '  Cost of optimal matching (correct)  = ', zfw_correct

      return
      end
      subroutine test12 ( )

c*********************************************************************72
c
cc TEST12 tests ZUFALL.
c
c  Modified:
c
c    27 November 2007
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 10 )

      logical bool(n)
      integer i
      integer j
      integer perm(n)
      integer seed
      integer seed_save

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST12'
      write ( *, '(a)' )
     &  '  ZUFALL generates pseudorandom permutations.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,a)' )
     &  '     I          SEED     1   2   3   4   5',
     &  '   6   7   8   9  10'
      write ( *, '(a)' ) ' '
      seed = 123456789

      do i = 1, 10

        seed_save = seed
        call zufall ( n, perm, bool, seed )

        write ( *, '(2x,i4,2x,i12,2x,10i4)' )
     &    i, seed_save, ( perm(j), j = 1, n )

      end do

      return
      end
