      program main

c*********************************************************************72
c
cc MAIN is the main program for LATIN_COVER_PRB.
c
c  Discussion:
c
c    LATIN_COVER_PRB tests the LATIN_COVER library.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 June 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LATIN_COVER_PRB:'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the LATIN_COVER library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'LATIN_COVER_PRB:'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests LATIN_COVER.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 June 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 9 )

      integer a(n_max,n_max)
      integer base
      integer n
      integer p(n_max)
      integer seed
      integer test

      base = 1

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  LATIN_COVER:'

      do n = 3, 9, 2

        seed = 123456789

        do test = 1, 3

          call perm_uniform ( n, base, seed, p )
     
          call perm_print ( n, p, '  Permutation' )

          call latin_cover ( n, p, a )

          call i4mat_print ( n, n, a, '  Latin cover' )

        end do

      end do

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests LATIN_COVER_2D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 June 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 9 )

      integer a(n_max,n_max)
      integer base
      integer n
      integer p(2,n_max)
      integer seed
      integer test

      base = 1

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  LATIN_COVER_2D:'

      do n = 3, 9, 2

        seed = 123456789

        do test = 1, 3

          call perm_uniform ( n, base, seed, p(1,1:n) )
     
          call perm_print ( n, p(1,1:n), '  Permutation 1' )

          call perm_uniform ( n, base, seed, p(2,1:n) ) 
     
          call perm_print ( n, p(2,1:n), '  Permutation 2' )
          call latin_cover_2d ( n, p, a )

          call i4mat_print ( n, n, a, '  Latin cover' )

        end do

      end do

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests LATIN_COVER_3D.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    26 June 2012
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n_max
      parameter ( n_max = 9 )

      integer a(n_max,n_max,n_max)
      integer base
      integer n
      integer p(3,n_max)
      integer seed
      integer test

      base = 1

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  LATIN_COVER_3D:'

      do n = 3, 9, 2

        seed = 123456789

        do test = 1, 3

          call perm_uniform ( n, base, seed, p(1,1:n) )
     
          call perm_print ( n, p(1,1:n), '  Permutation 1' )

          call perm_uniform ( n, base, seed, p(2,1:n) ) 
     
          call perm_print ( n, p(2,1:n), '  Permutation 2' )

          call perm_uniform ( n, base, seed, p(3,1:n) ) 
     
          call perm_print ( n, p(3,1:n), '  Permutation 3' )

          call latin_cover_3d ( n, p, a )

          call i4block_print ( n, n, n, a, '  Latin cover' )

        end do

      end do

      return
      end
