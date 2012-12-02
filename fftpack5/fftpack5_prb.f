      program main

c*********************************************************************72
c
cc FFTPACK5_PRB calls the FFTPACK5 test routines.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      call timestamp ( )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FFTPACK5_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the FFTPACK5 library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
      call test06 ( )
      call test07 ( )
      call test08 ( )
      call test09 ( )
      call test10 ( )
      call test11 ( )
      call test12 ( )
      call test13 ( )
      call test14 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FFTPACK5_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'

      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 tests CFFT1B, CFFT1F and CFFT1I.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4096 )

      integer lenwrk
      parameter ( lenwrk = 2 * n )

c     lensav = 2 * n + int ( log ( real ( n ) ) ) + 4
      integer lensav
      parameter ( lensav = 68144 )

      complex c(n)
      integer ier
      integer inc
      integer lenc
      integer seed
      real work(lenwrk)
      real wsave(lensav)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  For complex fast Fourier transforms, 1D,'
      write ( *, '(a)' ) '  CFFT1I initializes the transform,'
      write ( *, '(a)' ) '  CFFT1F does a forward transform;'
      write ( *, '(a)' ) '  CFFT1B does a backward transform.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  The number of data items is N = ', n
c
c  Set the data values.
c
      seed = 1973

      call c4vec_uniform_01 ( n, seed, c )

      call c4vec_print_some ( n, c, 10, '  The original data:' )
c
c  Initialize the WSAVE array.
c
      call cfft1i ( n, wsave, lensav, ier )
c
c  Compute the FFT coefficients.
c
      inc = 1
      lenc = n

      call cfft1f ( n, inc, c, lenc, wsave, lensav, work, lenwrk, ier )

      call c4vec_print_some ( n, c, 10, '  The FFT coefficients:' )
c
c  Compute inverse FFT of coefficients.  Should get back the
c  original data.
c
      call cfft1b ( n, inc, c, lenc, wsave, lensav, work, lenwrk, ier )

      call c4vec_print_some ( n, c, 10, '  The retrieved data:' )

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 tests CFFT2B, CFFT2F and CFFT2I.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer l
      integer m

      parameter ( l = 32 )
      parameter ( m = 64 )

      integer lenwrk
      parameter ( lenwrk = 2 * l * m )

c     parameter ( lensav = 2 * ( l + m ) + int ( log ( real ( l ) ) ) 
c   &  + int ( log ( real ( m ) ) ) + 8 )
      integer lensav
      parameter ( lensav = 208 )

      complex c(l,m)
      integer ier
      integer ldim
      integer seed
      real work(lenwrk)
      real wsave(lensav)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' ) '  For complex fast Fourier transforms, 2D,'
      write ( *, '(a)' ) '  CFFT2I initializes the transform,'
      write ( *, '(a)' ) '  CFFT2F does a forward transform;'
      write ( *, '(a)' ) '  CFFT2B does a backward transform.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The data is stored in an L by M array, with'
      write ( *, '(a,i4)' ) '  L = ', l
      write ( *, '(a,i4)' ) '  M = ', m
c
c  Set the data values.
c
      seed = 1973

      call c4mat_uniform_01 ( l, m, seed, c )

      call c4mat_print_some ( l, m, c, 1, 1, 5, 5, 
     &  '  Part of the original data:' )
c
c  Initialize the WSAVE array.
c
      call cfft2i ( l, m, wsave, lensav, ier )
c
c  Compute the FFT coefficients.
c
      ldim = l

      call cfft2f ( ldim, l, m, c, wsave, lensav, work, lenwrk, ier )

      call c4mat_print_some ( l, m, c, 1, 1, 5, 5, 
     &  '  Part of the FFT coefficients:' )
c
c  Compute inverse FFT of coefficients.  Should get back the
c  original data.
c
      call cfft2b ( ldim, l, m, c, wsave, lensav, work, lenwrk, ier )

      call c4mat_print_some ( l, m, c, 1, 1, 5, 5, 
     &  '  Part of the retrieved data:' )

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 tests CFFTMB, CFFTMF and CFFTMI.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lot

      parameter ( n = 32 )
      parameter ( lot = 6 )

      integer lenc
      integer lenwrk

      parameter ( lenc = n * lot )
      parameter ( lenwrk = 2 * lot * n )

c     lensav = 2 * n + int ( log ( real ( n ) ) ) + 4
      integer lensav
      parameter ( lensav = 72 )

      complex c(lenc)
      integer ier
      integer inc
      integer jump
      integer seed
      real work(lenwrk)
      real wsave(lensav)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  For complex fast Fourier transforms,'
      write ( *, '(a)' ) '  1D, multiple'
      write ( *, '(a)' ) '  CFFTMI initializes the transform,'
      write ( *, '(a)' ) '  CFFTMF does a forward transform;'
      write ( *, '(a)' ) '  CFFTMB does a backward transform.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  The number of sequences is LOT = ', lot
      write ( *, '(a,i4)' ) '  The length of each sequence is N = ', n
c
c  Set the data values.
c
      seed = 1973

      call c4mat_uniform_01 ( n, lot, seed, c )

      call c4mat_print_some ( n, lot, c, 1, 1, 5, 5, 
     &  '  Part of the original data:' )
c
c  Initialize the WSAVE array.
c
      call cfftmi ( n, wsave, lensav, ier )
c
c  Compute the FFT coefficients.
c
      jump = n
      inc = 1

      call cfftmf ( lot, jump, n, inc, c, lenc, wsave, lensav, 
     &  work, lenwrk, ier )

      call c4mat_print_some ( n, lot, c, 1, 1, 5, 5, 
     &  '  Part of the FFT coefficients:' )
c
c  Compute inverse FFT of coefficients.  Should get back the
c  original data.
c
      call cfftmb ( lot, jump, n, inc, c, lenc, wsave, lensav, 
     &  work, lenwrk, ier )

      call c4mat_print_some ( n, lot, c, 1, 1, 5, 5, 
     &  '  Part of the retrieved data:' )

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 tests COSQ1B, COSQ1F and COSQ1I.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4096 )

      integer lenwrk
      parameter ( lenwrk = n )

c     lensav = 2 * n + int ( log ( real ( n ) ) ) + 4
      integer lensav
      parameter ( lensav = 8205 )

      integer ier
      integer inc
      integer lenr
      real r(n)
      integer seed
      real work(lenwrk)
      real wsave(lensav)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  For real fast cosine transforms, 1D,'
      write ( *, '(a)' ) '  COSQ1I initializes the transform,'
      write ( *, '(a)' ) '  COSQ1F does a forward transform;'
      write ( *, '(a)' ) '  COSQ1B does a backward transform.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  The number of data items is N = ', n
c
c  Set the data values.
c
      seed = 1973

      call r4vec_uniform_01 ( n, seed, r )

      call r4vec_print_some ( n, r, 10, '  The original data:' )
c
c  Initialize the WSAVE array.
c
      call cosq1i ( n, wsave, lensav, ier )
c
c  Compute the FFT coefficients.
c
      inc = 1
      lenr = n

      call cosq1f ( n, inc, r, lenr, wsave, lensav, work, 
     &  lenwrk, ier )

      call r4vec_print_some ( n, r, 10, '  The FFT coefficients:' )
c
c  Compute inverse FFT of coefficients.  Should get back the
c  original data.
c
      call cosq1b ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

      call r4vec_print_some ( n, r, 10, '  The retrieved data:' )

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 tests COSQMB, COSQMF and COSQMI.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lot

      parameter ( n = 32 )
      parameter ( lot = 6 )

      integer lenr
      parameter ( lenr = n * lot )

      integer lenwrk
      parameter ( lenwrk = lot * n )

c      lensav = 2 * n + int ( log ( real ( n ) ) ) + 4
      integer lensav
      parameter ( lensav = 72 )

      integer ier
      integer inc
      integer jump
      real r(lenr)
      integer seed
      real work(lenwrk)
      real wsave(lensav)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' ) '  For real fast cosine transforms, '
      write ( *, '(a)' ) '  1D, multiple'
      write ( *, '(a)' ) '  COSQMI initializes the transform,'
      write ( *, '(a)' ) '  COSQMF does a forward transform;'
      write ( *, '(a)' ) '  COSQMB does a backward transform.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  The number of sequences is LOT = ', lot
      write ( *, '(a,i4)' ) '  The length of each sequence is N = ', n
c
c  Set the data values.
c
      seed = 1973

      call r4mat_uniform_01 ( n, lot, seed, r )

      call r4mat_print_some ( n, lot, r, 1, 1, 5, 5, 
     &  '  Part of the original data:' )
c
c  Initialize the WSAVE array.
c
      call cosqmi ( n, wsave, lensav, ier )
c
c  Compute the FFT coefficients.
c
      jump = n
      inc = 1

      call cosqmf ( lot, jump, n, inc, r, lenr, wsave, lensav, 
     &  work, lenwrk, ier )

      call r4mat_print_some ( n, lot, r, 1, 1, 5, 5, 
     &  '  Part of the FFT coefficients:' )
c
c  Compute inverse FFT of coefficients.  Should get back the
c  original data.
c
      call cosqmb ( lot, jump, n, inc, r, lenr, wsave, lensav, 
     &  work, lenwrk, ier )

      call r4mat_print_some ( n, lot, r, 1, 1, 5, 5, 
     &  '  Part of the retrieved data:' )

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 tests COST1B, COST1F and COST1I.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4096 )

      integer lenwrk
      parameter ( lenwrk = n - 1 )

c      lensav = 2 * n + int ( log ( real ( n ) ) ) + 4
      integer lensav
      parameter ( lensav = 8205 )

      integer ier
      integer inc
      integer lenr
      real r(n)
      integer seed
      real work(lenwrk)
      real wsave(lensav)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' ) '  For real fast cosine transforms, 1D,'
      write ( *, '(a)' ) '  COST1I initializes the transform,'
      write ( *, '(a)' ) '  COST1F does a forward transform;'
      write ( *, '(a)' ) '  COST1B does a backward transform.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  The number of data items is N = ', n
c
c  Set the data values.
c
      seed = 1973

      call r4vec_uniform_01 ( n, seed, r )

      call r4vec_print_some ( n, r, 10, '  The original data:' )
c
c  Initialize the WSAVE array.
c
      call cost1i ( n, wsave, lensav, ier )
c
c  Compute the FFT coefficients.
c
      inc = 1
      lenr = n

      call cost1f ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

      call r4vec_print_some ( n, r, 10, '  The FFT coefficients:' )
c
c  Compute inverse FFT of coefficients.  Should get back the
c  original data.
c
      call cost1b ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

      call r4vec_print_some ( n, r, 10, '  The retrieved data:' )

      return
      end
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07 tests COSTMB, COSTMF and COSTMI.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lot

      parameter ( n = 32 )
      parameter ( lot = 6 )

      integer lenr
      parameter ( lenr = n * lot )

      integer lenwrk
      parameter ( lenwrk = lot * ( n + 1 ) )

c      lensav = 2 * n + int ( log ( real ( n ) ) ) + 4
      integer lensav
      parameter ( lensav = 72 )

      integer ier
      integer inc
      integer jump
      real r(lenr)
      integer seed
      real work(lenwrk)
      real wsave(lensav)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07'
      write ( *, '(a)' ) '  For real fast cosine transforms, '
      write ( *, '(a)' ) '  1D, multiple'
      write ( *, '(a)' ) '  COSTMI initializes the transform,'
      write ( *, '(a)' ) '  COSTMF does a forward transform;'
      write ( *, '(a)' ) '  COSTMB does a backward transform.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  The number of sequences is LOT = ', lot
      write ( *, '(a,i4)' ) '  The length of each sequence is N = ', n
c
c  Set the data values.
c
      seed = 1973

      call r4mat_uniform_01 ( n, lot, seed, r )

      call r4mat_print_some ( n, lot, r, 1, 1, 5, 5, 
     &  '  Part of the original data:' )
c
c  Initialize the WSAVE array.
c
      call costmi ( n, wsave, lensav, ier )
c
c  Compute the FFT coefficients.
c
      jump = n
      inc = 1

      call costmf ( lot, jump, n, inc, r, lenr, wsave, lensav, 
     &  work, lenwrk, ier )

      call r4mat_print_some ( n, lot, r, 1, 1, 5, 5, 
     &  '  Part of the FFT coefficients:' )
c
c  Compute inverse FFT of coefficients.  Should get back the
c  original data.
c
      call costmb ( lot, jump, n, inc, r, lenr, wsave, lensav, 
     &  work, lenwrk, ier )

      call r4mat_print_some ( n, lot, r, 1, 1, 5, 5, 
     &  '  Part of the retrieved data:' )

      return
      end
      subroutine test08 ( )

c*********************************************************************72
c
cc TEST08 tests RFFT1B, RFFT1F and RFFT1I.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4096 )

      integer lenwrk
      parameter ( lenwrk = 2 * n )

c     lensav = n + int ( log ( real ( n ) ) ) + 4
      integer lensav
      parameter ( lensav = 34074 )

      integer ier
      integer inc
      integer lenr
      real r(n)
      integer seed
      real work(lenwrk)
      real wsave(lensav)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST08'
      write ( *, '(a)' ) '  For real fast Fourier transforms, 1D,'
      write ( *, '(a)' ) '  RFFT1I initializes the transform,'
      write ( *, '(a)' ) '  RFFT1F does a forward transform;'
      write ( *, '(a)' ) '  RFFT1B does a backward transform.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  The number of data items is N = ', n
c
c  Set the data values.
c
      seed = 1973

      call r4vec_uniform_01 ( n, seed, r )

      call r4vec_print_some ( n, r, 10, '  The original data:' )
c
c  Initialize the WSAVE array.
c
      call rfft1i ( n, wsave, lensav, ier )
c
c  Compute the FFT coefficients.
c
      inc = 1
      lenr = n

      call rfft1f ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

      call r4vec_print_some ( n, r, 10, '  The FFT coefficients:' )
c
c  Compute inverse FFT of coefficients.  Should get back the
c  original data.
c
      call rfft1b ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

      call r4vec_print_some ( n, r, 10, '  The retrieved data:' )

      return
      end
      subroutine test09 ( )

c*********************************************************************72
c
cc TEST09 tests RFFT2B, RFFT2F and RFFT2I.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer l
      integer m
      parameter ( l = 32 )
      parameter ( m = 64 )

      integer ldim
      parameter ( ldim = 2 * ( l / 2 + 1 ) )

      integer lenwrk
      parameter ( lenwrk = 2 * ldim * m )

c     lensav = 2 * ( l + m ) + int ( log ( real ( l ) ) ) 
c    &  + int ( log ( real ( m ) ) ) + 8
      integer lensav
      parameter ( lensav = 208 )

      integer ier
      real r(ldim,m)
      integer seed
      real work(lenwrk)
      real wsave(lensav)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST09'
      write ( *, '(a)' ) '  For real fast Fourier transforms, 2D,'
      write ( *, '(a)' ) '  RFFT2I initializes the transform,'
      write ( *, '(a)' ) '  RFFT2F does a forward transform;'
      write ( *, '(a)' ) '  RFFT2B does a backward transform.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  The L by M data is stored in an '
      write ( *, '(a)' ) '  LDIM by M array, with'
      write ( *, '(a,i4)' ) '  L = ', l
      write ( *, '(a,i4)' ) '  LDIM = ', ldim
      write ( *, '(a,i4)' ) '  M = ', m
c
c  Set the data values.
c
      seed = 1973

      call r4mat_uniform_01 ( ldim, m, seed, r )

      call r4mat_print_some ( ldim, m, r, 1, 1, 5, 5, 
     &  '  Part of the original data:' )
c
c  Initialize the WSAVE array.
c
      call rfft2i ( l, m, wsave, lensav, ier )
c
c  Compute the FFT coefficients.
c
      call rfft2f ( ldim, l, m, r, wsave, lensav, work, lenwrk, ier )

      call r4mat_print_some ( ldim, m, r, 1, 1, 5, 5, 
     &  '  Part of the FFT coefficients:' )
c
c  Compute inverse FFT of coefficients.  Should get back the
c  original data.
c
      call rfft2b ( ldim, l, m, r, wsave, lensav, work, lenwrk, ier )

      call r4mat_print_some ( ldim, m, r, 1, 1, 5, 5, 
     &  '  Part of the retrieved data:' )

      return
      end
      subroutine test10 ( )

c*********************************************************************72
c
cc TEST10 tests RFFTMB, RFFTMF and RFFTMI.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lot

      parameter ( n = 32 )
      parameter ( lot = 6 )

      integer lenr
      parameter ( lenr = n * lot )

      integer lenwrk
      parameter ( lenwrk = lot * n )

c     lensav = n + int ( log ( real ( n ) ) ) + 4
      integer lensav
      parameter ( lensav = 40 )

      integer ier
      integer inc
      integer jump
      real r(lenr)
      integer seed
      real work(lenwrk)
      real wsave(lensav)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST10'
      write ( *, '(a)' ) '  For real fast Fourier transforms, '
      write ( *, '(a)' ) '  1D, multiple'
      write ( *, '(a)' ) '  RFFTMI initializes the transform,'
      write ( *, '(a)' ) '  RFFTMF does a forward transform;'
      write ( *, '(a)' ) '  RFFTMB does a backward transform.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  The number of sequences is LOT = ', lot
      write ( *, '(a,i4)' ) '  The length of each sequence is N = ', n
c
c  Set the data values.
c
      seed = 1973

      call r4mat_uniform_01 ( n, lot, seed, r )

      call r4mat_print_some ( n, lot, r, 1, 1, 5, 5, 
     &  '  Part of the original data:' )
c
c  Initialize the WSAVE array.
c
      call rfftmi ( n, wsave, lensav, ier )
c
c  Compute the FFT coefficients.
c
      jump = n
      inc = 1

      call rfftmf ( lot, jump, n, inc, r, lenr, wsave, lensav, 
     &  work, lenwrk, ier )

      call r4mat_print_some ( n, lot, r, 1, 1, 5, 5, 
     &  '  Part of the FFT coefficients:' )
c
c  Compute inverse FFT of coefficients.  Should get back the
c  original data.
c
      call rfftmb ( lot, jump, n, inc, r, lenr, wsave, lensav, 
     &  work, lenwrk, ier )

      call r4mat_print_some ( n, lot, r, 1, 1, 5, 5, 
     &  '  Part of the retrieved data:' )

      return
      end
      subroutine test11 ( )

c*********************************************************************72
c
cc TEST11 tests SINQ1B, SINQ1F and SINQ1I.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4096 )

      integer lenwrk
      parameter ( lenwrk = n )

c     lensav = 2 * n + int ( log ( real ( n ) ) ) + 4
      integer lensav
      parameter ( lensav = 8205 )

      integer ier
      integer inc
      integer lenr
      real r(n)
      integer seed
      real work(lenwrk)
      real wsave(lensav)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST11'
      write ( *, '(a)' ) '  For real fast sine transforms, 1D,'
      write ( *, '(a)' ) '  SINQ1I initializes the transform,'
      write ( *, '(a)' ) '  SINQ1F does a forward transform;'
      write ( *, '(a)' ) '  SINQ1B does a backward transform.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  The number of data items is N = ', n
c
c  Set the data values.
c
      seed = 1973

      call r4vec_uniform_01 ( n, seed, r )

      call r4vec_print_some ( n, r, 10, '  The original data:' )
c
c  Initialize the WSAVE array.
c
      call sinq1i ( n, wsave, lensav, ier )
c
c  Compute the FFT coefficients.
c
      inc = 1
      lenr = n

      call sinq1f ( n, inc, r, lenr, wsave, lensav, work, 
     &  lenwrk, ier )

      call r4vec_print_some ( n, r, 10, '  The FFT coefficients:' )
c
c  Compute inverse FFT of coefficients.  Should get back the
c  original data.
c
      call sinq1b ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

      call r4vec_print_some ( n, r, 10, '  The retrieved data:' )

      return
      end
      subroutine test12 ( )

c*********************************************************************72
c
cc TEST12 tests SINQMB, SINQMF and SINQMI.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lot

      parameter ( n = 32 )
      parameter ( lot = 6 )

      integer lenr
      parameter ( lenr = n * lot )

      integer lenwrk
      parameter ( lenwrk = lot * n )

c     lensav = 2 * n + int ( log ( real ( n ) ) ) + 4
      integer lensav
      parameter ( lensav = 72 )

      integer ier
      integer inc
      integer jump
      real r(lenr)
      integer seed
      real work(lenwrk)
      real wsave(lensav)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST12'
      write ( *, '(a)' ) '  For real fast sine transforms, 1D, multiple'
      write ( *, '(a)' ) '  SINQMI initializes the transform,'
      write ( *, '(a)' ) '  SINQMF does a forward transform;'
      write ( *, '(a)' ) '  SINQMB does a backward transform.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  The number of sequences is LOT = ', lot
      write ( *, '(a,i4)' ) '  The length of each sequence is N = ', n
c
c  Set the data values.
c
      seed = 1973

      call r4mat_uniform_01 ( n, lot, seed, r )

      call r4mat_print_some ( n, lot, r, 1, 1, 5, 5, 
     &  '  Part of the original data:' )
c
c  Initialize the WSAVE array.
c
      call sinqmi ( n, wsave, lensav, ier )
c
c  Compute the FFT coefficients.
c
      jump = n
      inc = 1

      call sinqmf ( lot, jump, n, inc, r, lenr, wsave, lensav, 
     &  work, lenwrk, ier )

      call r4mat_print_some ( n, lot, r, 1, 1, 5, 5, 
     &  '  Part of the FFT coefficients:' )
c
c  Compute inverse FFT of coefficients.  Should get back the
c  original data.
c
      call sinqmb ( lot, jump, n, inc, r, lenr, wsave, lensav, 
     &  work, lenwrk, ier )

      call r4mat_print_some ( n, lot, r, 1, 1, 5, 5, 
     &  '  Part of the retrieved data:' )

      return
      end
      subroutine test13 ( )

c*********************************************************************72
c
cc TEST13 tests SINT1B, SINT1F and SINT1I.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      parameter ( n = 4096 )

      integer lenwrk
      parameter ( lenwrk = 2 * ( n + 1 ) )

c     lensav = n/2 + n + int ( log ( real ( n ) ) ) + 4
      integer lensav
      parameter ( lensav = 6157 )

      integer ier
      integer inc
      integer lenr
      real r(n)
      integer seed
      real work(lenwrk)
      real wsave(lensav)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST13'
      write ( *, '(a)' ) '  For real fast sine transforms, 1D,'
      write ( *, '(a)' ) '  SINT1I initializes the transform,'
      write ( *, '(a)' ) '  SINT1F does a forward transform;'
      write ( *, '(a)' ) '  SINT1B does a backward transform.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  The number of data items is N = ', n
c
c  Set the data values.
c
      seed = 1973

      call r4vec_uniform_01 ( n, seed, r )

      call r4vec_print_some ( n, r, 10, '  The original data:' )
c
c  Initialize the WSAVE array.
c
      call sint1i ( n, wsave, lensav, ier )
c
c  Compute the FFT coefficients.
c
      inc = 1
      lenr = n

      call sint1f ( n, inc, r, lenr, wsave, lensav, work, 
     &  lenwrk, ier )

      call r4vec_print_some ( n, r, 10, '  The FFT coefficients:' )
c
c  Compute inverse FFT of coefficients.  Should get back the
c  original data.
c
      call sint1b ( n, inc, r, lenr, wsave, lensav, work, lenwrk, ier )

      call r4vec_print_some ( n, r, 10, '  The retrieved data:' )

      return
      end
      subroutine test14 ( )

c*********************************************************************72
c
cc TEST14 tests SINTMB, SINTMF and SINTMI.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    27 March 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer n
      integer lot

      parameter ( n = 32 )
      parameter ( lot = 6 )

      integer lenr
      parameter ( lenr = n * lot )

      integer lenwrk
      parameter ( lenwrk = lot * 2 * ( n + 2 ) )

c     lensav = n / 2 + n + int ( log ( real ( n ) ) ) + 4
      integer lensav
      parameter ( lensav = 56 )

      integer ier
      integer inc
      integer jump
      real r(lenr)
      integer seed
      real work(lenwrk)
      real wsave(lensav)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST14'
      write ( *, '(a)' ) '  For real fast sine transforms, 1D, multiple'
      write ( *, '(a)' ) '  SINTMI initializes the transform,'
      write ( *, '(a)' ) '  SINTMF does a forward transform;'
      write ( *, '(a)' ) '  SINTMB does a backward transform.'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i4)' ) '  The number of sequences is LOT = ', lot
      write ( *, '(a,i4)' ) '  The length of each sequence is N = ', n
c
c  Set the data values.
c
      seed = 1973

      call r4mat_uniform_01 ( n, lot, seed, r )

      call r4mat_print_some ( n, lot, r, 1, 1, 5, 5, 
     &  '  Part of the original data:' )
c
c  Initialize the WSAVE array.
c
      call sintmi ( n, wsave, lensav, ier )
c
c  Compute the FFT coefficients.
c
      jump = n
      inc = 1

      call sintmf ( lot, jump, n, inc, r, lenr, wsave, lensav, 
     &  work, lenwrk, ier )

      call r4mat_print_some ( n, lot, r, 1, 1, 5, 5, 
     &  '  Part of the FFT coefficients:' )
c 
c  Compute inverse FFT of coefficients.  Should get back the
c  original data.
c
      call sintmb ( lot, jump, n, inc, r, lenr, wsave, lensav, 
     &  work, lenwrk, ier )

      call r4mat_print_some ( n, lot, r, 1, 1, 5, 5, 
     &  '  Part of the retrieved data:' )

      return
      end
      subroutine c4mat_print ( m, n, a, title )

c*********************************************************************72
c
cc C4MAT_PRINT prints a C4MAT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 March 2005
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns in the matrix.
c
c    Input, complex A(M,N), the matrix.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer m
      integer n

      complex a(m,n)
      character ( len = * ) title

      call c4mat_print_some ( m, n, a, 1, 1, m, n, title )

      return
      end
      subroutine c4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, title )

c*********************************************************************72
c
cc C4MAT_PRINT_SOME prints some of a C4MAT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 March 2005
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns in the matrix.
c
c    Input, complex A(M,N), the matrix.
c
c    Input, integer ILO, JLO, IHI, JHI, the first row and
c    column, and the last row and column to be printed.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 4 )
      integer m
      integer n

      complex a(m,n)
      character ( len = 20 ) ctemp(incx)
      integer i
      integer i2hi
      integer i2lo
      integer ihi
      integer ilo
      integer inc
      integer j
      integer j2
      integer j2hi
      integer j2lo
      integer jhi
      integer jlo
      character ( len = * ) title
      complex zero

      zero = cmplx ( 0.0E+00, 0.0E+00 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) title
c
c  Print the columns of the matrix, in strips of INCX.
c
      do j2lo = jlo, jhi, incx

        j2hi = j2lo + incx - 1
        j2hi = min ( j2hi, n )
        j2hi = min ( j2hi, jhi )

        inc = j2hi + 1 - j2lo

        write ( *, '(a)' ) ' '

        do j = j2lo, j2hi
          j2 = j + 1 - j2lo
          write ( ctemp(j2), '(i10,10x)' ) j
        end do

        write ( *, '(a,4a20)' ) '  Col: ', ( ctemp(j2), j2 = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) '  ---'
c
c  Determine the range of the rows in this strip.
c
        i2lo = max ( ilo, 1 )
        i2hi = min ( ihi, m )

        do i = i2lo, i2hi
c
c  Print out (up to) INCX entries in row I, that lie in the current strip.
c
          do j2 = 1, inc

            j = j2lo - 1 + j2

            if ( a(i,j) == zero ) then
              ctemp(j2) = '    0.0'
            else if ( imag ( a(i,j) ) == 0.0E+00 ) then
              write ( ctemp(j2), '(g10.3,10x)' ) real ( a(i,j) )
            else
              write ( ctemp(j2), '(2g10.3)' ) a(i,j)
            end if

          end do

          write ( *, '(i5,1x,4a20)' ) i, ( ctemp(j2), j2 = 1, inc )

        end do

      end do

      return
      end
      subroutine c4mat_uniform_01 ( m, n, seed, c )

c*********************************************************************72
c
cc C4MAT_UNIFORM_01 returns a unit pseudorandom C4MAT.
c
c  Discussion:
c
c    The angles should be uniformly distributed between 0 and 2 * PI,
c    the square roots of the radius uniformly distributed between 0 and 1.
c
c    This results in a uniform distribution of values in the unit circle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2005
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns in the matrix.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, complex C(M,N), the pseudorandom complex matrix.
c
      implicit none

      integer m
      integer n

      complex c(m,n)
      integer i
      integer j
      real r
      integer k
      real pi
      parameter ( pi = 3.1415926E+00 )
      integer seed
      real theta

      do j = 1, n
        do i = 1, m

          k = seed / 127773

          seed = 16807 * ( seed - k * 127773 ) - k * 2836

          if ( seed < 0 ) then
            seed = seed + 2147483647
          end if

          r = sqrt ( real ( dble ( seed ) * 4.656612875D-10 ) )

          k = seed / 127773

          seed = 16807 * ( seed - k * 127773 ) - k * 2836

          if ( seed < 0 ) then
            seed = seed + 2147483647
          end if

          theta = 2.0D+00 * pi * real ( dble ( seed ) 
     &      * 4.656612875D-10 )

          c(i,j) = r * cmplx ( cos ( theta ), sin ( theta ) )

        end do

      end do

      return
      end
      subroutine c4vec_print_some ( n, x, max_print, title )

c*********************************************************************72
c
cc C4VEC_PRINT_SOME prints some of a C4VEC.
c
c  Discussion:
c
c    The user specifies MAX_PRINT, the maximum number of lines to print.
c
c    If N, the size of the vector, is no more than MAX_PRINT, then
c    the entire vector is printed, one entry per line.
c
c    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
c    followed by a line of periods suggesting an omission,
c    and the last entry.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    17 December 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries of the vector.
c
c    Input, complex X(N), the vector to be printed.
c
c    Input, integer MAX_PRINT, the maximum number of lines to print.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer n

      integer i
      integer max_print
      character ( len = * ) title
      complex x(n)

      if ( max_print <= 0 ) then
        return
      end if

      if ( n <= 0 ) then
        return
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '

      if ( n <= max_print ) then

        do i = 1, n
          write ( *, '(i6,2x,2g14.6)' ) i, x(i)
        end do

      else if ( 3 <= max_print ) then

        do i = 1, max_print-2
          write ( *, '(i6,2x,2g14.6)' ) i, x(i)
        end do
        write ( *, '(a)' ) '......  ..............'
        i = n
        write ( *, '(i6,2x,2g14.6)' ) i, x(i)

      else

        do i = 1, max_print - 1
          write ( *, '(i6,2x,2g14.6)' ) i, x(i)
        end do
        i = max_print
        write ( *, '(i6,2x,2g14.6,2x,a)' ) i, x(i), '...more entries...'

      end if

      return
      end
      subroutine c4vec_uniform_01 ( n, seed, c )

c*********************************************************************72
c
cc C4VEC_UNIFORM_01 returns a unit pseudorandom C4VEC.
c
c  Discussion:
c
c    The angles should be uniformly distributed between 0 and 2 * PI,
c    the square roots of the radius uniformly distributed between 0 and 1.
c
c    This results in a uniform distribution of values in the unit circle.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    15 March 2005
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of values to compute.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, complex C(N), the pseudorandom complex vector.
c
      implicit none

      integer n

      complex c(n)
      integer i
      real r
      integer k
      real pi
      parameter ( pi = 3.1415926E+00 )
      integer seed
      real theta

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed < 0 ) then
          seed = seed + 2147483647
        end if

        r = sqrt ( real ( dble ( seed ) * 4.656612875D-10 ) )

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed < 0 ) then
          seed = seed + 2147483647
        end if

        theta = 2.0E+00 * pi * real ( dble ( seed ) * 4.656612875D-10 )

        c(i) = r * cmplx ( cos ( theta ), sin ( theta ) )

      end do

      return
      end
      subroutine r4mat_print_some ( m, n, a, ilo, jlo, ihi, jhi, 
     &  title )

c*********************************************************************72
c
cc R4MAT_PRINT_SOME prints some of an R4MAT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    12 September 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns.
c
c    Input, real A(M,N), an M by N matrix to be printed.
c
c    Input, integer ILO, JLO, the first row and column to print.
c
c    Input, integer IHI, JHI, the last row and column to print.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer incx
      parameter ( incx = 5 )

      integer m
      integer n

      real a(m,n)
      character ( len = 14 ) ctemp(incx)
      integer i
      integer i2hi
      integer i2lo
      integer ihi
      integer ilo
      integer inc
      integer j
      integer j2
      integer j2hi
      integer j2lo
      integer jhi
      integer jlo
      character ( len = * ) title

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )

      do j2lo = max ( jlo, 1 ), min ( jhi, n ), incx

        j2hi = j2lo + incx - 1
        j2hi = min ( j2hi, n )
        j2hi = min ( j2hi, jhi )

        inc = j2hi + 1 - j2lo

        write ( *, '(a)' ) ' '

        do j = j2lo, j2hi
          j2 = j + 1 - j2lo
          write ( ctemp(j2), '(i7,7x)') j
        end do

        write ( *, '(''  Col   '',5a14)' ) ( ctemp(j), j = 1, inc )
        write ( *, '(a)' ) '  Row'
        write ( *, '(a)' ) ' '

        i2lo = max ( ilo, 1 )
        i2hi = min ( ihi, m )

        do i = i2lo, i2hi

          do j2 = 1, inc

            j = j2lo - 1 + j2

            if ( a(i,j) == real ( int ( a(i,j) ) ) ) then
              write ( ctemp(j2), '(f8.0,6x)' ) a(i,j)
            else
              write ( ctemp(j2), '(g14.6)' ) a(i,j)
            end if

          end do

          write ( *, '(i5,1x,5a14)' ) i, ( ctemp(j), j = 1, inc )

        end do

      end do

      return
      end
      subroutine r4mat_uniform_01 ( m, n, seed, r )

c*********************************************************************72
c
cc R4MAT_UNIFORM_01 returns a unit pseudorandom R4MAT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 March 2005
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, L E Schrage,
c    A Guide to Simulation,
c    Springer Verlag, pages 201-202, 1983.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, pages 362-376, 1986.
c
c    P A Lewis, A S Goodman, J M Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, pages 136-143, 1969.
c
c  Parameters:
c
c    Input, integer M, N, the number of rows and columns in the array.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, real R(M,N), the array of pseudorandom values.
c
      implicit none

      integer m
      integer n

      integer i
      integer j
      integer k
      integer seed
      real r(m,n)

      do j = 1, n

        do i = 1, m

          k = seed / 127773

          seed = 16807 * ( seed - k * 127773 ) - k * 2836

          if ( seed < 0 ) then
            seed = seed + 2147483647
          end if

          r(i,j) = real ( seed ) * 4.656612875E-10

        end do
      end do

      return
      end
      subroutine r4vec_print_some ( n, a, max_print, title )

c*********************************************************************72
c
cc R4VEC_PRINT_SOME prints "some" of an R4VEC.
c
c  Discussion:
c
c    The user specifies MAX_PRINT, the maximum number of lines to print.
c
c    If N, the size of the vector, is no more than MAX_PRINT, then
c    the entire vector is printed, one entry per line.
c
c    Otherwise, if possible, the first MAX_PRINT-2 entries are printed,
c    followed by a line of periods suggesting an omission,
c    and the last entry.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 December 2001
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the number of entries of the vector.
c
c    Input, real A(N), the vector to be printed.
c
c    Input, integer MAX_PRINT, the maximum number of lines to print.
c
c    Input, character ( len = * ) TITLE, a title.
c
      implicit none

      integer n

      real a(n)
      integer i
      integer max_print
      character ( len = * ) title

      if ( max_print <= 0 ) then
        return
      end if

      if ( n <= 0 ) then
        return
      end if

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) trim ( title )
      write ( *, '(a)' ) ' '

      if ( n <= max_print ) then

        do i = 1, n
          write ( *, '(2x,i6,2x,g14.6)' ) i, a(i)
        end do

      else if ( 3 <= max_print ) then

        do i = 1, max_print-2
          write ( *, '(2x,i6,2x,g14.6)' ) i, a(i)
        end do
        write ( *, '(a)' ) '  ......  ..............'
        i = n
        write ( *, '(2x,i6,2x,g14.6)' ) i, a(i)

      else

        do i = 1, max_print - 1
          write ( *, '(2x,i6,2x,g14.6)' ) i, a(i)
        end do
        i = max_print
        write ( *, '(2x,i6,2x,g14.6,2x,a)' ) 
     &    i, a(i), '...more entries...'

      end if

      return
      end
      subroutine r4vec_uniform_01 ( n, seed, r )

c*********************************************************************72
c
cc R4VEC_UNIFORM_01 returns a unit pseudorandom R4VEC.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    22 March 2005
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Paul Bratley, Bennett Fox, L E Schrage,
c    A Guide to Simulation,
c    Springer Verlag, pages 201-202, 1983.
c
c    Bennett Fox,
c    Algorithm 647:
c    Implementation and Relative Efficiency of Quasirandom
c    Sequence Generators,
c    ACM Transactions on Mathematical Software,
c    Volume 12, Number 4, pages 362-376, 1986.
c
c    P A Lewis, A S Goodman, J M Miller,
c    A Pseudo-Random Number Generator for the System/360,
c    IBM Systems Journal,
c    Volume 8, pages 136-143, 1969.
c
c  Parameters:
c
c    Input, integer N, the number of entries in the vector.
c
c    Input/output, integer SEED, the "seed" value, which should NOT be 0.
c    On output, SEED has been updated.
c
c    Output, real R(N), the vector of pseudorandom values.
c
      implicit none

      integer n

      integer i
      integer k
      integer seed
      real r(n)

      do i = 1, n

        k = seed / 127773

        seed = 16807 * ( seed - k * 127773 ) - k * 2836

        if ( seed < 0 ) then
          seed = seed + 2147483647
        end if

        r(i) = real ( seed ) * 4.656612875E-10

      end do

      return
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
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 September 2005
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

      character ( len = 8 ) date
      character ( len = 10 ) time

      call date_and_time ( date, time )

      write ( *, '(a8,2x,a10)' ) date, time

      return
      end
