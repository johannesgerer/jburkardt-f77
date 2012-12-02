      program main

c*********************************************************************72
c
cc MAIN is the main program for ASA007_PRB.
c
c  Discussion:
c
c    ASA007_PRB calls the ASA007 routines.
c
c  Modified:
c
c    31 January 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA007_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the ASA007 library.'

      call test01 ( )
      call test02 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'ASA007_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 demonstrates the use of SYMINV.
c
c  Modified:
c
c    31 January 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 15 )

      double precision a((n_max*(n_max+1))/2)
      double precision afull(n_max,n_max)
      double precision c((n_max*(n_max+1))/2)
      double precision cfull(n_max,n_max)
      double precision cta
      double precision diff
      integer i
      integer ifault
      integer j
      integer k
      integer l
      integer n
      integer nullty
      double precision w(n_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) '  SYMINV computes the inverse of a positive'
      write ( *, '(a)' ) '  definite symmetric matrix.'
      write ( *, '(a)' ) '  A compressed storage format is used.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Here we look at the matrix A which is'
      write ( *, '(a)' ) '  N+1 on the diagonal and'
      write ( *, '(a)' ) '  N   on the off diagonals.'

      do n = 1, n_max
!
!  Set A to the lower triangle of the matrix which is N+1 on the diagonal
!  and N on the off diagonals.
!
        k = 0
        do i = 1, n
          do j = 1, i - 1
            k = k + 1
            a(k) = dble ( n )
          end do
          k = k + 1
          a(k) = dble ( n + 1 )
        end do

        call syminv ( a, n, c, w, nullty, ifault )

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Matrix order N = ', n
        write ( *, '(a,i8)' ) '  Maxtrix nullity NULLTY = ', nullty

        k = 0
        do j = 1, n
          do i = 1, j - 1
            k = k + 1
            cfull(i,j) = c(k)
            cfull(j,i) = c(k)
          end do
          k = k + 1
          cfull(j,j) = c(k)
        end do

        k = 0
        do j = 1, n
          do i = 1, j - 1
            k = k + 1
            afull(i,j) = a(k)
            afull(j,i) = a(k)
          end do
          k = k + 1
          afull(j,j) = a(k)
        end do
!
!  Compute C * A - I.
!
        diff = 0.0D+00
        do i = 1, n
          do j = 1, n
            cta = 0.0D+00
            do k = 1, n
              cta = cta + cfull(i,k) * afull(k,j)
            end do
            if ( i .eq. j ) then
              diff = diff + ( 1.0D+00 - cta )**2
            else
              diff = diff + cta**2
            end if
          end do
        end do

        diff = sqrt ( diff )

        write ( *, '(a,g14.6)' ) '  RMS ( C * A - I ) = ', diff

      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 demonstrates the use of SYMINV.
c
c  Modified:
c
c    31 January 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 15 )

      double precision a((n_max*(n_max+1))/2)
      double precision afull(n_max,n_max)
      double precision c((n_max*(n_max+1))/2)
      double precision cfull(n_max,n_max)
      double precision cta
      double precision diff
      integer i
      integer ifault
      integer j
      integer k
      integer l
      integer n
      integer nullty
      double precision w(n_max)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01:'
      write ( *, '(a)' ) '  SYMINV computes the inverse of a positive'
      write ( *, '(a)' ) '  definite symmetric matrix.'
      write ( *, '(a)' ) '  A compressed storage format is used.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Here we look at the Hilbert matrix'
      write ( *, '(a)' ) '  A(I,J) = 1/(I+J-1).'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  For this particular matrix, we expect the'
      write ( *, '(a)' ) '  errors to grow rapidly.'

      do n = 1, n_max
!
!  Set A to the lower triangle of the matrix which is N+1 on the diagonal
!  and N on the off diagonals.
!
        k = 0
        do i = 1, n
          do j = 1, i
            k = k + 1
            a(k) = 1.0D+00 / dble ( i + j - 1 )
          end do
        end do

        call syminv ( a, n, c, w, nullty, ifault )

        write ( *, '(a)' ) ' '
        write ( *, '(a,i8)' ) '  Matrix order N = ', n
        write ( *, '(a,i8)' ) '  Maxtrix nullity NULLTY = ', nullty

        k = 0
        do j = 1, n
          do i = 1, j - 1
            k = k + 1
            cfull(i,j) = c(k)
            cfull(j,i) = c(k)
          end do
          k = k + 1
          cfull(j,j) = c(k)
        end do

        k = 0
        do j = 1, n
          do i = 1, j - 1
            k = k + 1
            afull(i,j) = a(k)
            afull(j,i) = a(k)
          end do
          k = k + 1
          afull(j,j) = a(k)
        end do
!
!  Compute C * A - I.
!
        diff = 0.0D+00
        do i = 1, n
          do j = 1, n
            cta = 0.0D+00
            do k = 1, n
              cta = cta + cfull(i,k) * afull(k,j)
            end do
            if ( i .eq. j ) then
              diff = diff + ( 1.0D+00 - cta )**2
            else
              diff = diff + cta**2
            end if
          end do
        end do

        diff = sqrt ( diff )

        write ( *, '(a,g14.6)' ) '  RMS ( C * A - I ) = ', diff

      end do

      return
      end
