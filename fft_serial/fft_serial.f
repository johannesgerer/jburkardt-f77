      program main

c*********************************************************************72
c
cc MAIN is the main program for FFT_SERIAL.
c
c  Discussion:
c
c    The complex data in an N vector is stored as pairs of values in a
c    double precision real vector of length 2*N.
c
c  Modified:
c
c    23 March 2009
c
c  Author:
c
c    Original C version by Wesley Petersen.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Wesley Petersen, Peter Arbenz, 
c    Introduction to Parallel Computing - A practical guide with examples in C,
c    Oxford University Press,
c    ISBN: 0-19-851576-6,
c    LC: QA76.58.P47.
c
      implicit none

      integer n_max
      parameter ( n_max = 131072 )

      double precision ctime
      double precision ctime1
      double precision ctime2
      double precision error
      logical first
      double precision flops
      double precision fnm1
      double precision ggl
      integer i
      integer icase
      integer id
      integer it
      integer ln2
      double precision mflops
      integer n
      integer nits
      double precision seed
      double precision sgn
      double precision w(n_max)
      double precision x(2*n_max)
      double precision y(2*n_max)
      double precision z(2*n_max)
      double precision z0
      double precision z1

      nits = 10000

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FFT_SERIAL'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 
     &  '  Demonstrate an implementation of the Fast Fourier Transform'
      write ( *, '(a)' ) '  of a complex data vector,'
c
c  Prepare for tests.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Accuracy check:'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '    FFT ( FFT ( X(1:N) ) ) == N * X(1:N)'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '             N      NITS    Error         ' //
     &  'Time          Time/Call     MFLOPS'
      write ( *, '(a)' ) ' '

      seed  = 331.0D+00
      n = 1
c
c  LN2 is the log base 2 of N.  Each increase of LN2 doubles N.
c
      do ln2 = 1, 17

        n = 2 * n
c
c  Allocate storage for the complex arrays W, X, Y, Z.  
c
c  We handle the complex arithmetic,
c  and store a complex number as a pair of floats, a complex vector as a doubly
c  dimensioned array whose second dimension is 2. 
c
        first = .true.

        do icase = 0, 1

          if ( first ) then

            do i = 1, 2 * n - 1, 2

              z0 = ggl ( seed )
              z1 = ggl ( seed )
              x(i) = z0
              z(i) = z0
              x(i+1) = z1
              z(i+1) = z1
            end do

          else

            do i = 1, 2 * n - 1, 2
              z0 = 0.0D+00
              z1 = 0.0D+00
              x(i) = z0
              z(i) = z0
              x(i+1) = z1
              z(i+1) = z1
            end do

          end if
c
c  Initialize the sine and cosine tables.
c
          call cffti ( n, w )
c
c  Transform forward, back 
c
          if ( first ) then

            sgn = + 1.0D+00
            call cfft2 ( n, x, y, w, sgn )
            sgn = - 1.0D+00
            call cfft2 ( n, y, x, w, sgn )
c 
c  Results should be same as initial multiplied by n 
c
            fnm1 = 1.0D+00 / dble ( n )

            error = 0.0D+00
            do i = 1, 2 * n - 1, 2
              error = error 
     &        + ( z(i)   - fnm1 * x(i) )**2 
     &        + ( z(i+1) - fnm1 * x(i+1) )**2
            end do
            error = sqrt ( fnm1 * error )

            first = .false.

          else

            call cpu_time ( ctime1 )

            do it = 1, nits

              sgn = + 1.0D+00
              call cfft2 ( n, x, y, w, sgn )
              sgn = - 1.0D+00
              call cfft2 ( n, y, x, w, sgn )

            end do

            call cpu_time ( ctime2 )
            ctime = ctime2 - ctime1

            flops = 2.0D+00 * dble ( nits ) 
     &        * ( 5.0D+00 * dble ( n ) * dble ( ln2 ) )

            mflops = flops / 1.0D+06 / ctime

            write ( *, 
     &        '(2x,i12,2x,i8,2x,g12.4,2x,g14.6,2x,g12.4,2x,g12.4)' ) 
     &        n, nits, error, ctime, ctime / dble ( 2 * nits ), mflops

          end if

        end do

        if ( mod ( ln2, 4 ) .eq. 0 ) then
          nits = nits / 10
        end if

        if ( nits .lt. 1 ) then
          nits = 1
        end if

      end do

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'FFT_SERIAL:'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine cfft2 ( n, x, y, w, sgn )

c*********************************************************************72
c
cc CFFT2 performs a complex Fast Fourier Transform.
c
c  Modified:
c
c    23 March 2009
c
c  Author:
c
c    Original C version by Wesley Petersen.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Wesley Petersen, Peter Arbenz, 
c    Introduction to Parallel Computing - A practical guide with examples in C,
c    Oxford University Press,
c    ISBN: 0-19-851576-6,
c    LC: QA76.58.P47.
c
c  Parameters:
c
c    Input, integer N, the size of the array to be transformed.
c
c    Input/output, double precision X(2*N), the data to be transformed.  
c    On output, the contents of X have been overwritten by work information.
c
c    Output, double precision Y(2*N), the forward or backward FFT of X.
c
c    Input, double precision W(N), a table of sines and cosines.
c
c    Input, double precision SGN, is +1 for a "forward" FFT 
C    and -1 for a "backward" FFT.
c
      implicit none

      integer n

      integer i
      integer j
      integer m
      integer mj
      double precision sgn
      logical tgle
      double precision w(n)
      double precision x(2*n)
      double precision y(2*n)

       m = int ( ( log ( dble ( n ) ) / log ( 1.99D+00 ) ) )
       mj = 1
c
c  Toggling switch for work array.
c
      tgle = .true.
      call step ( n, mj, x(1), x((n/2)*2+1), y(1), y(mj*2+1), w, sgn )

      if ( n .eq. 2 ) then
        return
      end if

      do j = 1, m - 2  

        mj = mj * 2

        if ( tgle ) then
          call step ( n, mj, y(1), y((n/2)*2+1), x(1), x(mj*2+1), 
     &      w, sgn )
          tgle = .false.
        else
          call step ( n, mj, x(1), x((n/2)*2+1), y(1), y(mj*2+1), 
     &      w, sgn )
          tgle = .true.
        end if

      end do
c
c  Last pass through data: move Y to X if needed. 
c
      if ( tgle ) then
        do i = 1, 2 * n
          x(i) = y(i)
        end do
      end if

      mj = n / 2
      call step ( n, mj, x(1), x((n/2)*2+1), y(1), y(mj*2+1), w, sgn )

      return
      end
      subroutine cffti ( n, w )

c*********************************************************************72
c
cc CFFTI sets up sine and cosine tables needed for the FFT calculation.
c
c  Modified:
c
c    23 March 2009
c
c  Author:
c
c    Original C version by Wesley Petersen.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Wesley Petersen, Peter Arbenz, 
c    Introduction to Parallel Computing - A practical guide with examples in C,
c    Oxford University Press,
c    ISBN: 0-19-851576-6,
c    LC: QA76.58.P47.
c
c  Parameters:
c
c    Input, integer N, the size of the array to be transformed.
c
c    Output, double precision W(N), a table of sines and cosines.
c
      implicit none

      integer n

      double precision arg
      double precision aw
      integer i
      integer n2
      double precision pi
      parameter ( pi = 3.141592653589793D+00 )
      double precision w(n)

      n2 = n / 2
      aw = 2.0D+00 * pi / dble ( n )

      do i = 1, n2
        arg = aw * dble ( i - 1 )
        w(2*i-1) = cos ( arg )
        w(2*i)   = sin ( arg )
      end do

      return
      end
      function ggl ( seed )

c*********************************************************************72
c
cc GGL generates uniformly distributed pseudorandom numbers. 
c
c  Modified:
c
c    23 March 2009
c
c  Author:
c
c    Original C version by Wesley Petersen, M Troyer, I Vattulainen.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Wesley Petersen, Peter Arbenz, 
c    Introduction to Parallel Computing - A practical guide with examples in C,
c    Oxford University Press,
c    ISBN: 0-19-851576-6,
c    LC: QA76.58.P47.
c
c  Parameters:
c
c    Input/output, double precision SEED, used as a seed for the sequence.
c
c    Output, double precision GGL, the next pseudorandom value.
c
      implicit none

      double precision d2
      parameter ( d2 = 0.2147483647D+10 )
      double precision ggl
      double precision seed
      double precision t

      t = mod ( 16807.0D+00 * seed, d2 )
      seed = t
      ggl = dble ( ( t - 1.0D+00 ) / ( d2 - 1.0D+00 ) )

      return
      end
      subroutine step ( n, mj, a, b, c, d, w, sgn )

c*********************************************************************72
c
cc STEP carries out one step of the workspace version of CFFT2.
c
c  Modified:
c
c    23 March 2009
c
c  Author:
c
c    Original C version by Wesley Petersen.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    Wesley Petersen, Peter Arbenz, 
c    Introduction to Parallel Computing - A practical guide with examples in C,
c    Oxford University Press,
c    ISBN: 0-19-851576-6,
c    LC: QA76.58.P47.
c
c  Parameters:
c
      implicit none

      integer n

      double precision a(n)
      double precision ambr
      double precision ambu
      double precision b(n)
      double precision c(n)
      double precision d(n)
      integer j
      integer ja
      integer jb
      integer jc
      integer jd
      integer jw
      integer k
      integer lj
      integer mj
      integer mj2
      double precision sgn
      double precision w(n)
      double precision wjw(2)

      mj2 = 2 * mj
      lj = n / mj2

      do j = 0, lj - 1

        jw = j * mj
        ja  = jw
        jb  = ja
        jc  = j * mj2
        jd  = jc

        wjw(1) = w(jw*2+1) 
        wjw(2) = w(jw*2+2)

        if ( sgn .lt. 0.0D+00 ) then
          wjw(2) = - wjw(2)
        end if

        do k = 0, mj - 1

          c((jc+k)*2+1) = a((ja+k)*2+1) + b((jb+k)*2+1)
          c((jc+k)*2+2) = a((ja+k)*2+2) + b((jb+k)*2+2)

          ambr = a((ja+k)*2+1) - b((jb+k)*2+1)
          ambu = a((ja+k)*2+2) - b((jb+k)*2+2)

          d((jd+k)*2+1) = wjw(1) * ambr - wjw(2) * ambu
          d((jd+k)*2+2) = wjw(2) * ambr + wjw(1) * ambu

        end do
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

      write ( *, 
     &  '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) 
     &  d, month(m), y, h, ':', n, ':', s, '.', mm, ampm

      return
      end
