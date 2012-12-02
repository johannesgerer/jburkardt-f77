      function r4_exp ( jsr, ke, fe, we )

c*********************************************************************72
c
cc R4_EXP returns an exponentially distributed single precision real value.
c
c  Discussion:
c
c    The underlying algorithm is the ziggurat method.
c
c    Before the first call to this function, the user must call R4_EXP_SETUP
c    to determine the values of KE, FE and WE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 May 2008
c
c  Author:
c
c    Original C version by George Marsaglia, Wai Wan Tsang.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    George Marsaglia, Wai Wan Tsang,
c    The Ziggurat Method for Generating Random Variables,
c    Journal of Statistical Software,
c    Volume 5, Number 8, October 2000, seven pages.
c
c  Parameters:
c
c    Input/output, integer JSR, the seed.
c
c    Input, integer KE(256), data computed by R4_EXP_SETUP.
c
c    Input, real FE(256), WE(256), data computed by R4_EXP_SETUP.
c
c    Output, real R4_EXP, an exponentially distributed random value.
c
      implicit none

      real fe(256)
      integer iz
      integer jsr
      integer jz
      integer ke(256)
      real r4_exp
      real r4_uni
      integer shr3
      real value
      real we(256)
      real x

      jz = shr3 ( jsr )
      iz = iand ( jz, 255 )

      if ( abs ( jz  ) .lt. ke(iz+1) ) then

        value = abs ( jz ) * we(iz+1)

      else

10      continue

          if ( iz .eq. 0 ) then
            value = 7.69711E+00 - log ( r4_uni ( jsr ) )
            go to 20
          end if

          x = abs ( jz ) * we(iz+1)

          if ( fe(iz+1) + r4_uni ( jsr ) * ( fe(iz) - fe(iz+1) ) .lt. 
     &      exp ( - x ) )  then
            value = x
            go to 20
          end if

          jz = shr3 ( jsr )
          iz = iand ( jz, 255 )

          if ( abs ( jz ) .lt. ke(iz+1) ) then
            value = abs ( jz ) * we(iz+1)
            go to 20
          end if

        go to 10

20      continue
        
      end if

      r4_exp = value

      return
      end
      subroutine r4_exp_setup ( ke, fe, we )

c*********************************************************************72
c
cc R4_EXP_SETUP sets data needed by R4_EXP.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 May 2008
c
c  Author:
c
c    Original C version by George Marsaglia, Wai Wan Tsang.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    George Marsaglia, Wai Wan Tsang,
c    The Ziggurat Method for Generating Random Variables,
c    Journal of Statistical Software,
c    Volume 5, Number 8, October 2000, seven pages.
c
c  Parameters:
c
c    Output, integer KE(256), data needed by R4_EXP.
c
c    Output, real FE(256), WE(256), data needed by R4_EXP.
c
      implicit none

      double precision de
      real fe(256)
      integer i
      integer ke(256)
      double precision m2
      parameter ( m2 = 2147483648.0D+00 )
      double precision q
      double precision te
      double precision ve
      parameter ( ve = 3.949659822581572D-03 )
      real we(256)

      de = 7.697117470131487D+00
      te = 7.697117470131487D+00

      q = ve / exp ( - de )
      ke(1) = int ( ( de / q ) * m2 )
      ke(2) = 0

      we(1) = q / m2
      we(256) = de / m2

      fe(1) = 1.0
      fe(256) = exp ( - de )

      do i = 255, 2, -1
        de = - log ( ve / de + exp ( - de ) )
        ke(i+1) = int ( ( de / te ) * m2 )
        te = de
        fe(i) = exp ( - de )
        we(i) = de / m2
      end do

      return
      end
      function r4_nor ( jsr, kn, fn, wn )

c*********************************************************************72
c
cc R4_NOR returns a normally distributed single precision real value.
c
c  Discussion:
c
c    The value returned is generated from a distribution with mean 0 and 
c    variance 1.
c
c    The underlying algorithm is the ziggurat method.
c
c    Before the first call to this function, the user must call R4_NOR_SETUP
c    to determine the values of KN, FN and WN.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 May 2008
c
c  Author:
c
c    Original C version by George Marsaglia, Wai Wan Tsang.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    George Marsaglia, Wai Wan Tsang,
c    The Ziggurat Method for Generating Random Variables,
c    Journal of Statistical Software,
c    Volume 5, Number 8, October 2000, seven pages.
c
c  Parameters:
c
c    Input/output, integer JSR, the seed.
c
c    Input, integer KN(128), data computed by R4_NOR_SETUP.
c
c    Input, real FN(128), WN(128), data computed by R4_NOR_SETUP.
c
c    Output, real R4_NOR, a normally distributed random value.
c
      implicit none

      real fn(128)
      integer hz
      integer iz
      integer jsr
      integer kn(128)
      real r
      parameter ( r = 3.442620E+00 )
      real r4_nor
      real r4_uni
      integer shr3
      real value
      real wn(128)
      real x
      real y

      hz = shr3 ( jsr )
      iz = iand ( hz, 127 )

      if ( abs ( hz ) .lt. kn(iz+1) ) then

        value = hz * wn(iz+1)

      else

10      continue

          if ( iz .eq. 0 ) then

20          continue
              x = - 0.2904764E+00 * log ( r4_uni ( jsr ) )
              y = - log ( r4_uni ( jsr ) )
              if ( x * x .le. y + y ) then
                go to 30
              end if
            go to 20

30          continue

            if ( hz .le. 0 ) then
              value = - r - x
            else
              value = + r + x
            end if

            go to 40

          end if

          x = hz * wn(iz+1)

          if ( fn(iz+1) + r4_uni ( jsr ) * ( fn(iz) - fn(iz+1) ) .lt. 
     &      exp ( - 0.5E+00 * x * x ) ) then
            value = x
            go to 40
          end if

          hz = shr3 ( jsr )
          iz = iand ( hz, 127 )

          if ( abs ( hz ) .lt. kn(iz+1) ) then
            value = hz * wn(iz+1)
            go to 40
          end if

        go to 10

40      continue

      end if

      r4_nor = value

      return
      end
      subroutine r4_nor_setup ( kn, fn, wn )

c*********************************************************************72
c
cc R4_NOR_SETUP sets data needed by R4_NOR.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 May 2008
c
c  Author:
c
c    Original C version by George Marsaglia, Wai Wan Tsang.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    George Marsaglia, Wai Wan Tsang,
c    The Ziggurat Method for Generating Random Variables,
c    Journal of Statistical Software,
c    Volume 5, Number 8, October 2000, seven pages.
c
c  Parameters:
c
c    Output, integer KN(128), data needed by R4_NOR.
c
c    Output, real FN(128), WN(128), data needed by R4_NOR.
c
      implicit none

      double precision dn
      real fn(128)
      integer i
      integer kn(128)
      double precision m1
      parameter ( m1 = 2147483648.0D+00 )
      double precision q
      double precision tn
      double precision vn
      parameter ( vn = 9.91256303526217D-03 )
      real wn(128)

      dn = 3.442619855899D+00
      tn = 3.442619855899D+00

      q = vn / exp ( - 0.5 * dn * dn )
      kn(1) = ( dn / q ) * m1
      kn(2) = 0

      wn(1) = q / m1
      wn(128) = dn / m1

      fn(1) = 1.0
      fn(128) = exp ( - 0.5E+00 * dn * dn )

      do i = 127, 2, -1
        dn = sqrt ( - 2.0D+00 * log ( vn / dn 
     &    + exp ( - 0.5D+00 * dn * dn ) ) )
        kn(i+1) = int ( ( dn / tn ) * m1 )
        tn = dn
        fn(i) = exp ( - 0.5 * dn * dn )
        wn(i) = dn / m1
      end do

      return
      end
      function r4_uni ( jsr )

c*********************************************************************72
c
cc R4_UNI returns a uniformly distributed real value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    20 May 2008
c
c  Author:
c
c    Original C version by George Marsaglia, Wai Wan Tsang.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    George Marsaglia, Wai Wan Tsang,
c    The Ziggurat Method for Generating Random Variables,
c    Journal of Statistical Software,
c    Volume 5, Number 8, October 2000, seven pages.
c
c  Parameters:
c
c    Input/output, integer JSR, the seed.
c
c    Output, real R4_UNI, a uniformly distributed random value in
c    the range [0,1].
c
      implicit none

      integer jsr
      integer jsr_input
      real r4_uni

      jsr_input = jsr

      jsr = ieor ( jsr, ishft ( jsr,   13 ) )
      jsr = ieor ( jsr, ishft ( jsr, - 17 ) )
      jsr = ieor ( jsr, ishft ( jsr,    5 ) )

      r4_uni = 0.5E+00 + 0.2328306E-09 * real ( jsr_input + jsr )

      return
      end
      function shr3 ( jsr )

c*********************************************************************72
c
cc SHR3 evaluates the SHR3 generator for integers.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 May 2008
c
c  Author:
c
c    Original C version by George Marsaglia, Wai Wan Tsang.
c    FORTRAN77 version by John Burkardt.
c
c  Reference:
c
c    George Marsaglia, Wai Wan Tsang,
c    The Ziggurat Method for Generating Random Variables,
c    Journal of Statistical Software,
c    Volume 5, Number 8, October 2000, seven pages.
c
c  Parameters:
c
c    Input/output, integer JSR, the seed.
c
c    Output, integer SHR3, the value of the SHR3 generator.
c
      implicit none

      integer jsr
      integer jsr_input
      integer shr3

      jsr_input = jsr

      jsr = ieor ( jsr, ishft ( jsr,   13 ) )
      jsr = ieor ( jsr, ishft ( jsr, - 17 ) )
      jsr = ieor ( jsr, ishft ( jsr,    5 ) )

      shr3 = jsr_input + jsr

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
