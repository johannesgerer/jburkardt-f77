      program main

c*********************************************************************72
c
cc MAIN is the main program for MATMUL.
c
c  Discussion:
c
c    MATMUL is a matrix multiplication performance test.
c
c    MATMUL is an interactive FORTRAN77 program that sets up, carries
c    out, and times a matrix multiplication.  MATMUL can use several
c    different algorithms and matrix sizes, and can be run on many
c    different computers.
c
c
c  Making a version for a given machine:
c
c  1) Modify MATMUL_CPU_TIMER
c
c    MATMUL_CPU_TIMER is a routine to compute the current reading of
c    the CPU timer in seconds.  Several sets of appropiate code
c    are listed in the routine, depending on the machine and compiler
c    used.  Uncomment the one appropriate for your machine, or
c    add a new one if necessary.
c
c  2) Set the value of "MACHINE" to be the name of your machine.
c
c    This occurs in routine INIT.
c
c  3) Set the maximum problem size, LENA.
c
c     The value of LENA should be adjusted to the memory available
c     on your machine.  MATMUL needs three real matrices of size
c     LENA by LENA, and three integer, three double precision,
c     and three complex.
c
c     On the Cray, six real matrices are needed, rather than three.
c
c     My choices for LENA have been:
c
c       65 for Apple Macintosh
c       150 for Apple PowerMac
c       513 for Cray YMP or Cray C90.
c       513 for DEC Alpha.
c       100 for DECstation/ULTRix.
c       65 for IBM PC
c       300 for SGI/IRIS.
c       300 for SUN.
c       300 for VAX/VMS or VECTOR/VAX.
c
c    LENA is defined only once, in a parameter statement in the
c    main program.
c
c  4) Make some special routines available for special machines:
c
c     For the Cray,
c       uncomment the calls to
c         MXMA
c         SAXPY
c         SDOT
c         SGEMM
c         SGEMMS
c         TIMEF.
c       uncomment the declaration of the array WORK in routine USGEMMS.
c
c     For the SGI/IRIS,
c       uncomment the calls to
c         SAXPY
c         SDOT
c         SECNDS
c         SGEMM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 November 2008
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer lena
      parameter ( lena = 1100 )

      character*20 command
      logical cshow
      logical fshow
      integer ido
      integer ierror
      integer itemp
      integer lchar
      integer lda
      integer lenc
      integer s_length
      logical s_eqi
      character*7 lingo
      logical lnshow
      logical lshow
      character*10 machine
      logical mshow
      integer n
      integer nhi
      integer ninc
      integer nlo
      integer nmult
      logical noshow
      integer nrep
      logical nrshow
      logical nshow
      character*6 order
      character*82 output
      logical tshow

      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'MATMUL'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  An interactive demonstration of the speed'
      write ( *, '(a)' ) '  of matrix multiplication.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  This is version 1.19'
      write ( *, '(a)' ) '  Last modified on 28 August 1999.'
c
c  Initialize the data.
c
      call init ( command, cshow, fshow, lda, lena, lingo, lnshow, 
     &  lshow, machine, mshow, n, nhi, ninc, nlo, nmult, noshow,
     &  nrep, nrshow, nshow, order, tshow )
c
c  Say hello to the user.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  This is the version for ' // machine
      write ( *, '(a,i8)' ) 
     &  '  The maximum matrix order allowed is N = ', lena

10    continue
c
c  Print the prompt.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'Command?  (Type H for help)'
c
c  Read the command.
c
      read ( *, '(a)', end = 40 ) command
c
c  Capitalize the command.
c
      call s_cap ( command )
c
c  Remove all blanks to make interpretation easier.
c
      call chrdb1 ( command )

      ierror = 0
 
20    continue
c
c  Command "A" means abort the run.
c
      if ( s_eqi ( command(1:1), 'A' ) ) then
 
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'MATMUL:'
        write ( *, '(a)' ) '  Normal end of execution.'
        write ( *, '(a)' ) ' '
        call timestamp ( )
 
        stop
c
c  Command "H" means print out the help message.
c
      else if ( s_eqi ( command(1:1), 'H' ) ) then
 
        call help
c
c  Command "LDA = " means the user wants to set lda.
c
      else if ( s_eqi ( command(1:4), 'LDA=' ) ) then
 
        call chrcti ( command(5:), itemp, ierror, lchar )
 
        if ( ierror .ne. 0 ) then

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 
     &    'I did not understand your definition of LDA.'

        else if ( itemp .le. 0 ) then

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'The assignment of LDA was not acceptable!'
          write ( *, '(a)' ) 'LDA must be positive.'

        else if ( itemp .gt. lena ) then

          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'The assignment of LDA was not acceptable!'
          write ( *, '(a,i8)' ) 'LDA must be no greater than ', lena

        else

          lda = itemp
          write ( *, '(a)' ) ' '
          write ( *, '(a,i8)' ) 'LDA has been set to ', lda

        end if
 
        if ( lda .lt. n ) then
          n = lda
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 
     &      'Note: Since N must be no greater than LDA,'
          write ( *, '(a)' ) 'MATMUL has decreased the value of N.'
          write ( *, '(a,i8)' ) 'N has been set to ', n
        end if
c
c  Command "M" means the user wants the multiplication to be carried out.
c
      else if ( s_eqi ( command(1:1), 'M' ) ) then
c
c  Carry out multiplication for one, or many values of N.
c
        n = nlo
 
30      continue
 
        call mult ( cshow, fshow, lda, lena, lingo, lnshow, lshow, 
     &    machine, mshow, n, noshow, nrep, nrshow, nshow, order, 
     &    output, tshow )

        call nstep ( ido, n, nhi, ninc, nlo, nmult )

        if ( ido .eq. 1 ) then
          go to 30
        end if
 
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 
     &    'The matrix multiplication has been carried out.'
c
c  Command "N=" means the user wants to set the matrix size N.
c
      else if ( s_eqi ( command(1:2), 'N=' ) ) then
 
        call getn ( command(3:), ierror, lda, lena, n, nhi, ninc, nlo, 
     &    nmult )
c
c  Command "NOSHOW" means the user wants to turn off the display of all
c  quantities.
c
      else if ( s_eqi ( command, 'NOSHOW' ) ) then
 
        command = 'NOSHOW=ALL'
 
        call getsho ( command(8:), cshow, fshow, lnshow,
     &    lshow, .false., mshow, noshow, nrshow, nshow, tshow )
c
c  COMMAND "NOSHOW=" means the user wants to turn off the display of a
c  particular quantity.
c
      else if ( s_eqi ( command(1:7), 'NOSHOW=' ) ) then
 
        call getsho ( command(8:), cshow, fshow, lnshow, lshow,
     &    .false., mshow, noshow, nrshow, nshow, tshow )
c
c  Command "NREP=" sets the number of repetitions.
c
      else if ( s_eqi ( command(1:5), 'NREP=' ) ) then
 
        call chrcti( command(6:), nrep, ierror, lchar )
        write ( *, '(a,i8)' ) 
     &    'The repetition factor is now NREP = ', nrep
 
        if ( nrep .eq. 1 ) then
          nrshow = .false.
        else
          nrshow = .true.
        end if
c
c  Command "ORDER=" means the user wants to set the method.
c
      else if ( s_eqi ( command(1:6), 'ORDER=' ) ) then
 
        call get_order ( command(7:), ierror, machine, order )
c
c  Command "Q" means the user wants to quit.
c
      else if ( s_eqi ( command(1:1), 'Q' ) ) then
 
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Type "Y" to confirm that you want to quit.'

        read ( *, '(a)', end = 20 ) command
        call s_cap ( command )

        if ( s_eqi ( command(1:1), 'Y' ) ) then
          command = 'abort'
        end if
 
        go to 20
c
c  Command "P" means the user wants to print out the current settings.
c
      else if ( s_eqi ( command(1:1), 'P' ) ) then
 
        call printr ( lda, lena, lingo, n, nhi, ninc,
     &    nlo, nmult, nrep, order )
 
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'Valid choices for the order are:'
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'ALL, C4_IJK, R8_IJK, I4_IJK, UIJK, IUJK,'
        write ( *, '(a)' ) 'IJUK, IJK, IKJ, JIK, JKI, KIJ, KJI, L4_IJK'
        write ( *, '(a)' ) 'TAXPYC, TAXPYR, TDOT, TGEMM'
 
        if ( s_eqi ( machine(1:4), 'CRAY') ) then
          write ( *, '(a)' ) 'MIJK, MXMA, SAXPYC, SAXPYR, SDOT,'
          write ( *, '(a)' ) 'SGEMM, SGEMMS, SIJK.'
        else if ( s_eqi ( machine, 'SGI/IRIS' ) ) then
          write ( *, '(a)' ) 'MKJI, SAXPYC, SAXPYR, SDOT, SGEMM.'
        end if
c
c  Command "SHOW" means the user wants all items to be displayed.
c
      else if ( s_eqi ( command, 'SHOW' ) ) then
 
        command = 'SHOW=ALL'
 
        call getsho ( command(6:), cshow, fshow, lnshow, lshow,
     &    .true., mshow, noshow, nrshow, nshow, tshow )
c
c  Command "SHOW=" means the user wants a particular item displayed.
c
      else if ( s_eqi(command(1:5), 'SHOW=' ) ) then
 
        call getsho ( command(6:), cshow, fshow, lnshow, lshow,
     &    .true., mshow, noshow, nrshow, nshow, tshow )
c
c  The user's input did not match any acceptable command.
c
      else

        lenc = s_length ( command )

        if ( lenc .gt. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'Your command was not recognized.'
          write ( *, '(a)' ) 'You typed "' // trim ( command ) // '".'
          write ( *, '(a)' ) 'Type HELP for a list of commands.'
        end if

      end if
 
      go to 10
c
c  We jump here on certain input errors.
c
40    continue

      command = 'ABORT'
      go to 20
      end
      subroutine ch_cap ( c )

c*********************************************************************72
c
cc CH_CAP capitalizes a single character.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 July 1998
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character*1 C, the character to capitalize.
c
      implicit none

      character*1 c
      integer itemp

      itemp = ichar ( c )

      if ( 97 .le. itemp .and. itemp .le. 122 ) then
        c = char ( itemp - 32 )
      end if

      return
      end
      subroutine chrchp ( string, ilo, ihi )

c*********************************************************************72
c
cc CHRCHP "chops out" a portion of a string, and closes up the hole.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    06 July 1998
c
c  Example:
c
c    STRING = 'Fred is not a jerk!'
c
c    call chrchp ( STRING, 9, 12 )
c
c    STRING = 'Fred is a jerk!    '
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character*(*) STRING, the string to be transformed.
c
c    Input, integer ILO, IHI, the locations of the first and last
c    characters to be removed.  ILO must be at least 1, IHI must
c    be at least ILO, and no more than the length of the string.
c
      implicit none

      integer ihi
      integer ihi2
      integer ilo
      integer ilo2
      integer lens
      character*(*) string

      lens = len ( string )

      ilo2 = max ( ilo, 1 )
      ihi2 = min ( ihi, lens )

      if ( ilo2 .gt. ihi2 ) then
        return
      end if

      string(ilo2:lens+ilo2-ihi2-1) = string(ihi2+1:lens)
      string(lens+ilo2-ihi2:lens) = ' '

      return
      end
      subroutine chrcti ( string, intval, ierror, lchar )

c*********************************************************************72
c
cc CHRCTI reads an integer from a string.
c
c  Discussion:
c
c    CHRCTI will read as many characters as possible until it reaches
c    the end of the STRING, or encounters a character which cannot be
c    part of the number.
c
c    Legal input is
c
c      blanks,
c      initial sign,
c      blanks,
c      integer part,
c      blanks,
c      final comma or semicolon,
c
c    with most quantities optional.
c
c  Example:
c
c    STRING            INTVAL
c
c    '1'               1
c    '     1   '       1
c    '1A'              1
c    '12,34,56'        12
c    '  34 7'          34
c    '-1E2ABCD'        -100
c    '-1X2ABCD'        -1
c    ' 2E-1'           0
c    '23.45'           23
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 November 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) STRING, the string containing the
c    data to be read.  Reading will begin at position 1 and
c    terminate at the end of the string, or when no more
c    characters can be read to form a legal integer.  Blanks,
c    commas, or other nonnumeric data will, in particular,
c    cause the conversion to halt.
c
c    Output, integer INTVAL, the integer read from the string.
c
c    Output, integer IERROR, error flag.
c    0 if no errors,
c    Value of IHAVE when error occurred otherwise.
c
c    Output, integer LCHAR, number of characters read from
c    STRING to form the number.
c
      implicit none

      character*1 chrtmp
      integer ierror
      integer ihave
      integer intval
      integer isgn
      integer iterm
      integer itop
      integer lchar
      integer nchar
      integer ndig
      character*(*) string

      nchar = len ( string )

      ierror = 0
      intval = 0
      lchar = -1
      isgn = 1
      itop = 0
      ihave = 1
      iterm = 0

10    continue

      lchar = lchar + 1
      chrtmp = string(lchar+1:lchar+1)
c
c  Blank.
c
      if ( chrtmp .eq. ' ' ) then

        if ( ihave .eq. 2 ) then

        else if ( ihave .eq. 3 ) then
          ihave = 11
        end if
c
c  Comma.
c
      else if ( chrtmp .eq. ',' .or.
     &  chrtmp .eq. ';' ) then

        if ( ihave .ne. 1 ) then
          iterm = 1
          ihave = 12
          lchar = lchar + 1
        end if
c
c  Minus sign.
c
      else if ( chrtmp .eq. '-' ) then

        if ( ihave .eq. 1 ) then
          ihave = 2
          isgn = -1
        else
          iterm = 1
        end if
c
c  Plus sign.
c
      else if ( chrtmp .eq. '+' ) then

        if ( ihave .eq. 1 ) then
          ihave = 2
        else
          iterm = 1
        end if
c
c  Digit.
c
      else if ( ihave .lt. 11 .and.
     &  lge ( chrtmp, '0' ) .and.
     &  lle ( chrtmp, '9' ) ) then

        ihave = 3

        call digten ( chrtmp, ndig )

        itop = 10 * itop + ndig
c
c  Anything else is regarded as a terminator.
c
      else
        iterm = 1
      end if

      if ( iterm .ne. 1 .and. lchar+1 .lt. nchar ) then
        go to 10
      end if

      if ( iterm .ne. 1 .and. lchar+1 .eq. nchar ) then
        lchar = nchar
      end if
c
c  Number seems to have terminated.  Have we got a legal number?
c
      if ( ihave .eq. 1 .or. ihave .eq. 2 ) then
        ierror = ihave
        return
      end if
c
c  Number seems OK.  Form it.
c
      intval = isgn * itop

      return
      end
      subroutine chrdb1 ( string )

c*********************************************************************72
c
cc CHRDB1 removes blanks from a string, left justifying the remainder.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 November 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character*(*) STRING, the string to be transformed.
c
      implicit none

      character*1 chrtmp
      integer i
      integer j
      integer nchar
      character*(*) string

      nchar = len ( string )

      j = 0
      do i = 1, nchar
        chrtmp = string(i:i)
        string(i:i) = ' '

        if ( chrtmp .ne. ' ' ) then
          j = j + 1
          string(j:j) = chrtmp
        end if

      end do

      return
      end
      subroutine c4_ijk ( a, b, c, lda, n, nrep, ttime )

c*********************************************************************72
c
cc C4_IJK computes A = B*C using index order IJK and C4 arithmetic.
c
c  Discussion:
c
c    C4 arithmetic uses "single precision" complex values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Workspace, complex A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
c    used in the multiplication.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
c    Input, integer NREP, the number of times the multiplication should
c    be carried out.
c
c    Output, real TTIME, an estimate of the CPU time in seconds required
c    for the matrix multiplications.
c
      implicit none

      integer lda
      integer n

      complex a(lda,n)
      complex b(lda,n)
      complex c(lda,n)
      integer i
      integer irep
      integer j
      integer k
      integer nrep
      real time1
      real time2
      real ttime

      ttime = 0.0E+00
 
      do irep = 1, nrep
 
        call c4_set ( a, b, c, lda, n )
 
        call matmul_cpu_timer ( time1 )
 
        do i = 1, n
          do j = 1, n
            do k = 1, n
              a(i,k) = a(i,k) + b(i,j) * c(j,k)
            end do
          end do
        end do
 
        call matmul_cpu_timer ( time2 )
        ttime = ttime + time2 - time1

      end do
 
      return
      end
      subroutine c4_set ( a, b, c, lda, n )

c*********************************************************************72
c
cc C4_SET initializes the matrices for C4 arithmetic.
c
c  Discussion:
c
c    C4 arithmetic uses "single precision" complex values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, complex A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
c    used in the multiplication.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
      implicit none

      integer lda
      integer n

      complex a(lda,n)
      complex b(lda,n)
      complex c(lda,n)
      integer i
      integer j

      do i = 1, n
        do j = 1, n
          a(i,j) = cmplx ( 0.0E+00, 0.0E+00 )
          b(i,j) = cmplx ( 2.0E+00, 1.0E+00 )
          c(i,j) = cmplx ( 1.0E+00, 1.0E+00 )
        end do
      end do
 
      return
      end
      subroutine digten ( tenval, intval )

c*********************************************************************72
c
cc DIGTEN returns the integer value of a base 10 digit.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    19 November 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*1 TENVAL, the decimal digit, '0' through '9'.
c
c    Output, integer INTVAL, the corresponding integer value.
c
      implicit none

      integer intval
      character*1 tenval

      if ( lge ( tenval, '0' ) .and. lle ( tenval, '9' ) ) then

        intval = ichar ( tenval ) - 48

      else if ( tenval .eq. ' ' ) then

        intval = 0

      else

        write ( *, * ) ' '
        write ( *, * ) 'DIGTEN - Serious error!'
        write ( *, * ) '  Illegal decimal digit = ' // tenval
        write ( *, * ) '  ASCII number ', ichar ( tenval )
        intval = 0
        stop

      end if

      return
      end
      subroutine domethod ( lda, lena, n, nrep, order, ttime )

c*********************************************************************72
c
cc DOMETHOD calls a specific multiplication routine.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 March 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer LENA, the maximum dimension.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
c    Input, integer NREP, the number of times the multiplication should
c    be carried out.
c
c    Input, character*6 ORDER, specifies the method to be used.
c
c    Output, real TTIME, an estimate of the CPU time in seconds required
c    for the matrix multiplications.
c
      implicit none

      integer lda
      integer lena
      integer n

      real a(lena,lena)
      real b(lena,lena)
      real c(lena,lena)
      complex ca(lena,lena)
      complex cb(lena,lena)
      complex cc(lena,lena)
      double precision da(lena,lena)
      double precision db(lena,lena)
      double precision dc(lena,lena)
      integer ia(lena,lena)
      integer ib(lena,lena)
      integer ic(lena,lena)
      logical la(lena,lena)
      logical lb(lena,lena)
      logical lc(lena,lena)
      logical s_eqi
      integer nrep
      character*6 order
      real ttime

      if ( s_eqi ( order, 'C4_IJK' ) ) then

        call c4_ijk ( ca, cb, cc, lda, n, nrep, ttime )

      else if ( s_eqi ( order, 'I4_IJK' ) ) then

        call i4_ijk ( ia, ib, ic, lda, n, nrep, ttime )

      else if ( s_eqi ( order, 'IJK' ) ) then
 
        call ijk ( a, b, c, lda, n, nrep, ttime )

      else if ( s_eqi ( order, 'IJUK' ) ) then

        call ijuk ( a, b, c, lda, n, nrep, ttime )
 
      else if ( s_eqi ( order, 'IKJ' ) ) then
 
        call ikj ( a, b, c, lda, n, nrep, ttime )
 
      else if ( s_eqi ( order, 'IUJK' ) ) then
 
        call iujk ( a, b, c, lda, n, nrep, ttime )

      else if ( s_eqi ( order, 'JIK' ) ) then
 
        call jik ( a, b, c, lda, n, nrep, ttime )

      else if ( s_eqi ( order, 'JKI' ) ) then
 
        call jki ( a, b, c, lda, n, nrep, ttime )

      else if ( s_eqi ( order, 'KIJ' ) ) then
 
        call kij ( a, b, c, lda, n, nrep, ttime )

      else if ( s_eqi ( order, 'KJI' ) ) then
 
        call kji ( a, b, c, lda, n, nrep, ttime )

      else if ( s_eqi ( order, 'L4_IJK' ) ) then

        call l4_ijk ( la, lb, lc, lda, n, nrep, ttime )

      else if ( s_eqi ( order, 'MIJK' ) ) then

        call mijk ( a, b, c, lda, n, nrep, ttime )
 
      else if ( s_eqi ( order, 'MKJI' ) ) then

        call mkji ( a, b, c, lda, n, nrep, ttime )

      else if ( s_eqi ( order, 'MXMA' ) ) then

        call umxma ( a, b, c, lda, n, nrep, ttime )

      else if ( s_eqi ( order, 'R8_IJK' ) ) then

        call r8_ijk ( da, db, dc, lda, n, nrep, ttime )

      else if ( s_eqi ( order, 'SAXPYC' ) ) then

        call usaxpyc ( a, b, c, lda, n, nrep, ttime )

      else if ( s_eqi ( order, 'SAXPYR' ) ) then

        call usaxpyr ( a, b, c, lda, n, nrep, ttime )

      else if ( s_eqi ( order, 'SDOT' ) ) then

        call usdot ( a, b, c, lda, n, nrep, ttime )

      else if ( s_eqi ( order, 'SGEMM' ) ) then

        call usgemm ( a, b, c, lda, n, nrep, ttime )

      else if ( s_eqi ( order, 'SGEMMS' ) ) then

        call usgemms ( a, b, c, lda, n, nrep, 3*lena*lena, ttime )

      else if ( s_eqi ( order, 'SIJK' ) ) then

        call sijk ( a, b, c, lda, n, nrep, ttime )

      else if ( s_eqi ( order, 'TAXPYC' ) ) then

        call utaxpyc ( a, b, c, lda, n, nrep, ttime )

      else if ( s_eqi ( order, 'TAXPYR' ) ) then

        call utaxpyr ( a, b, c, lda, n, nrep, ttime )

      else if ( s_eqi ( order, 'TDOT' ) ) then

        call utdot ( a, b, c, lda, n, nrep, ttime )

      else if ( s_eqi ( order, 'TGEMM' ) ) then

        call utgemm ( a, b, c, lda, n, nrep, ttime )

      else if ( s_eqi ( order, 'UIJK' ) ) then

        call uijk ( a, b, c, lda, n, nrep, ttime )

      else
 
        ttime = 909.0

      end if

      return
      end
      subroutine getn ( string, ierror, lda, lena, n, nhi, ninc, nlo, 
     &  nmult )

c*********************************************************************72
c
cc GETN determines the problem sizes desired by the user.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) STRING, a string containing the user's
c    data for an "N" command.  This command might have one of the
c    following forms, where "n", "nlo", "nhi", "ninc" and "nmult"
c    are actually numeric values:
c
c      n                   Solve a single problem of size N.
c      nlo, nhi            Solve every problem from N = NLO to N = NHI.
c      nlo, nhi, ninc      Solve from N = NLO to N = NHI, incrementing by NINC.
c      nlo, nhi, *nmult    Solve from N = NLO to N = NHI, multiplying by NMULT.
c
c    Output, integer IERROR, an error flag.
c    0, no error occurred.
c    nonzero, an error occurred, and the operation could not be done.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer LENA, the maximum matrix order allowed.
c
c    Output, integer N, the new value for the number of rows and columns 
c    to use in the matrices.
c
c    Output, integer NHI, the maximum value of N to use.
c
c    Output, integer NINC, the additive increment to use, if additive
c    steps are being taken.
c
c    Output, integer NLO, the smallest value of N to use.
c
c    Output, integer NMULT, the multiplier to use, if multiplicative
c    steps are being taken.
c
      implicit none

      integer ierror
      integer imult
      integer itemp
      integer lchar
      integer lda
      integer lena
      integer n
      integer next
      integer nhi
      integer ninc
      integer nlo
      integer nmult
      character*(*) string

      nhi = 0
      ninc = 0
      nmult = 0
      nlo = 0
c
c  Read the first number, N or NLO.
c
      call chrcti ( string, itemp, ierror, lchar )
      next = lchar + 1

      if ( ierror .ne. 0 ) then
        write ( *, * ) ' '
        write ( *, * ) 'I could not understand your definition of N.'
        return
      else if ( itemp .le. 0 ) then
        write ( *, * ) ' '
        write ( *, * ) 'This value of N is not acceptable!'
        write ( *, * ) 'N must be positive.'
        return
      else if ( itemp .gt. lena ) then
        write ( *, * ) ' '
        write ( *, * ) 'This value of N is not acceptable!'
        write ( *, * ) 'N must be no greater than ', lena
        return
      else
        n = itemp
        nlo = itemp
        nhi = itemp
        write ( *, * ) 'N has been set to ', n
      end if
c
c  Read the second number, NHI.
c
      call chrcti ( string(next:), itemp, ierror, lchar )
      next = next + lchar

      if ( ierror .ne. 0 ) then
c       write ( *, * ) ' '
c       write ( *, * ) 'I could not understand your definition of NHI.'
        go to 10
      else if ( itemp .le. 0 ) then
        write ( *, * ) ' '
        write ( *, * ) 'This value of NHI is not acceptable!'
        write ( *, * ) 'NHI must be positive.'
        go to 10
      else if ( itemp .gt. lena ) then
        write ( *, * ) ' '
        write ( *, * ) 'This value of NHI is not acceptable!'
        write ( *, * ) 'NHI must be no greater than ', lena
        go to 10
      else
        nhi = itemp
        write ( *, * ) 'NHI has been set to ', nhi
        ninc = 1
      end if
 
      if ( string(next:next) .eq. '*' ) then
        next = next+1
        imult = 1
      else
        imult = 0
      end if
c
c  Read third number, ninc or nmult
c
      call chrcti ( string(next:), itemp, ierror, lchar )

      if ( ierror .ne. 0 ) then
c       write ( *, * ) ' '
c       write ( *, * ) 'I could not understand your definition of NINC.'
      else
        if ( imult .eq. 0 ) then
          ninc = itemp
          nmult = 0
          write ( *, * ) 'NINC has been set to ', ninc
        else
          ninc = 0
          nmult = itemp
          write ( *, * ) 'NMULT has been set to ', nmult
        end if
      end if
c
c  Check that LDA is no less than NLO and NHI.
c
10    continue

      if ( lda .lt. max ( nlo, nhi ) ) then
        lda = max ( nlo, nhi )
        write ( *, * ) ' '
        write ( *, * ) 'Note: Since LDA must be at least as large as'
        write ( *, * ) 'N, MATMUL has increased the value of LDA.'
        write ( *, * ) 'LDA has been set to ', lda
      end if
 
      return
      end
      subroutine get_order ( string, ierror, machine, order )

c*********************************************************************72
c
cc GET_ORDER reads a new value of order from the user.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) STRING, a string of characters, containing
c    the user's choice for the matrix multiplication method to employ.
c
c    Output, integer IERROR, an error flag.
c    0, no error occurred.
c    nonzero, an error occurred and the operation could not be completed.
c
c    Input, character*10 MACHINE, the name of the machine for which
c    the package is prepared.
c
c    Output, character*6 ORDER, specifies the method to be used.
c
      implicit none

      integer nchoice1
      parameter ( nchoice1 = 18 )

      character*6 choices1(nchoice1)
      integer i
      integer ierror
      logical s_eqi
      character*10 machine
      character*6 order
      character*(*) string

      data  ( choices1(i), i = 1, nchoice1 ) /
     &  'ALL'   , 'C4_IJK', 'R8_IJK'  , 'IJK'   , 'IJUK'  ,
     &  'IKJ'   , 'IUJK'  , 'JIK'   , 'JKI'   , 'KIJ'   ,
     &  'KJI'   , 'L4_IJK'  , 'I4_IJK'  , 'TAXPYC', 'TAXPYR',
     &  'TDOT'  , 'TGEMM' , 'UIJK' /

      do i = 1, nchoice1
        if ( s_eqi ( order, choices1(i) ) ) then
          go to 10
        end if
      end do

      if ( s_eqi(machine(1:4),'cray') ) then

        if ( s_eqi(string(1:4),'mijk'))go to 10
        if ( s_eqi(string(1:4),'mxma'))go to 10
        if ( s_eqi(string(1:6),'saxpyc'))go to 10
        if ( s_eqi(string(1:6),'saxpyc'))go to 10
        if ( s_eqi(string(1:4),'sdot'))go to 10
        if ( s_eqi(string(1:5),'sgemm'))go to 10
        if ( s_eqi(string(1:6),'sgemms'))go to 10
        if ( s_eqi(string(1:4),'sijk'))go to 10

      end if
 
      if (s_eqi(machine(1:8),'sgi/iris') ) then

        if (s_eqi(string(1:4),'mkji'))go to 10
        if (s_eqi(string(1:6),'saxpyc'))go to 10
        if (s_eqi(string(1:6),'saxpyc'))go to 10
        if (s_eqi(string(1:4),'sdot'))go to 10
        if (s_eqi(string(1:5),'sgemm'))go to 10

      end if
 
      write ( *, * ) ' '
      write ( *, * ) 'The order you chose was not a valid choice.'
      write ( *, * ) 'This order was "' // string // '".'
 
      write ( *, * ) ' '
      write ( *, * ) 'Valid choices for the order are:'
      write ( *, * ) ' '
      write ( *, * ) 'ALL, C4_IJK, R8_IJK, IJK, IKJ, IJUK, IUJK,'
      write ( *, * ) 'JIK, JKI, KIJ, KJI, L4_IJK, I4_IJK,'
      write ( *, * ) 'TAXPYC, TAXPYR, TDOT, TGEMM, UIJK'
 
      if (s_eqi(machine(1:4),'cray') ) then
        write ( *, * ) 'MIJK, MXMA, SAXPYC, SAXPYR, SDOT, '
        write ( *, * ) 'SGEMM, SGEMMS, SIJK.'
      else if (s_eqi(machine,'sgi/iris') ) then
        write ( *, * ) 'MKJI, SAXPYC, SAXPYR, SDOT, SGEMM.'
      end if
 
      write ( *, * ) ' '
      write ( *, * ) 'Your command was not carried out.'
 
      ierror = 1
      return
 
10    continue
 
      order = string(1:6)

      write ( *, * ) ' '
      write ( *, * ) 'The order has been set to ', order
 
      return
      end
      subroutine getsho ( string, cshow, fshow, lnshow, lshow,
     &  lval, mshow, noshow, nrshow, nshow, tshow )

c*********************************************************************72
c
cc GETSHO determines what items the user wishes to print out.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character STRING*(*), a string which specifies which
c    variable is to be assigned the input value of LVAL.
c
c    Output, logical CSHOW, is TRUE if the multiplication method
c    is to be shown.
c
c    Output, logical FSHOW, is TRUE if the MegaFLOP rate is to be shown.
c
c    Output, logical LNSHOW, is TRUE if the software's programming
c    language, stored in LINGO, is to be shown.
c
c    Output, logical LSHOW, is TRUE if the variable LDA is to be shown.
c
c    Input, logical LVAL, a logical value, which is to be assigned
c    to one of the variables.
c
c    Output, logical MSHOW, is TRUE if the machine name is to be shown.
c
c    Output, logical NOSHOW, is TRUE if the number of operations is
c    to be shown.
c
c    Output, logical NRSHOW, is TRUE if the number of repetitions is
c    to be shown.
c
c    Output, logical NSHOW, is TRUE if the matrix size N is to be shown.
c
c    Output, logical TSHOW, is TRUE if the execution time is to be shown.
c
      implicit none

      logical cshow
      logical fshow
      logical s_eqi
      logical lnshow
      logical lshow
      logical lval
      logical mshow
      logical noshow
      logical nrshow
      logical nshow
      character*(*) string
      logical tshow

      if (s_eqi(string(1:3),'all') ) then
        cshow = lval
        fshow = lval
        lnshow = lval
        lshow = lval
        mshow = lval
        noshow = lval
        nrshow = lval
        nshow = lval
        tshow = lval
      else if (s_eqi(string(1:3),'cpu') ) then
        tshow = lval
      else if (s_eqi(string(1:8),'language') ) then
        lnshow = lval
      else if (s_eqi(string(1:3),'lda') ) then
        lshow = lval
      else if (s_eqi(string(1:7),'machine') ) then
        mshow = lval
      else if (s_eqi(string(1:6),'mflops') ) then
        fshow = lval
      else if (s_eqi(string(1:1),'n') ) then
        nshow = lval
      else if (s_eqi(string(1:4),'nrep') ) then
        nrshow = lval
      else if (s_eqi(string(1:3),'ops') ) then
        noshow = lval
      else if (s_eqi(string(1:5),'order') ) then
        cshow = lval
      else if (s_eqi(string(1:4),'time') ) then
        tshow = lval
      else
        write ( *, * ) ' '
        write ( *, * ) 'That is not a legal name!'
        write ( *, * )
     &  'Legal names are CPU, LANGUAGE, LDA, MACHINE, MFLOPS,'//
     &  ' N, NREP, OPS, ORDER, and TIME.'
      end if
 
      return
      end
      subroutine header ( cshow, fshow, lnshow, lshow, mshow,
     &  noshow, nrshow, nshow, output, tshow )

c*********************************************************************72
c
cc HEADER prints out a header for the results.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, logical CSHOW, is TRUE if the multiplication method
c    is to be shown.
c
c    Input, logical FSHOW, is TRUE if the MegaFLOP rate is to be shown.
c
c    Input, logical LNSHOW, is TRUE if the software's programming
c    language, stored in LINGO, is to be shown.
c
c    Input, logical LSHOW, is TRUE if the variable LDA is to be shown.
c
c    Input, logical MSHOW, is TRUE if the machine name is to be shown.
c
c    Input, logical NOSHOW, is TRUE if the number of operations is
c    to be shown.
c
c    Input, logical NRSHOW, is TRUE if the number of repetitions is
c    to be shown.
c
c    Input, logical NSHOW, is TRUE if the matrix size N is to be shown.
c
c    Output, character*82 OUTPUT, a string containing a header for
c    all the variables which are to be printed.
c
c    Input, logical TSHOW, is TRUE if the execution time is to be shown.
c
      implicit none

      logical cshow
      logical fshow
      integer s_length
      logical lnshow
      logical lshow
      logical mshow
      integer next
      logical noshow
      logical nrshow
      logical nshow
      character*82 output
      logical tshow
c
c  Prepare the header string.  
c
      next = 1
  
      if ( cshow ) then
        output(next:) = ' Order'
        next = s_length(output) + 1
      end if

      if ( lshow ) then
        output(next:) = ' LDA'
        next = s_length(output) + 1
      end if

      if ( nshow ) then
        output(next:) = '   N'
        next = s_length(output) + 1
      end if

      if ( cshow ) then
        output(next:) = '      CPU'
        next = s_length(output) + 1
      end if

      if ( tshow ) then
        output(next:) = ' Secs'
        next = s_length(output) + 1
      end if

      if ( noshow ) then
        output(next:) = '       Ops'
        next = s_length(output) + 1
      end if

      if ( nrshow ) then
        output(next:) = ' NREP'
        next = s_length(output) + 1
      end if

      if ( fshow ) then
        output(next:) = '    MFLOPS'
        next = s_length(output) + 1
      end if

      if ( mshow ) then
        output(next:) = '  Machine'
        next = s_length(output) + 1
      end if
 
      if ( lnshow ) then
        output(next:) = '  Language'
        next = s_length(output) + 1
      end if

      write ( *, * ) ' '
      write ( *, '(a)' ) output(1:next-1)
      write ( *, * ) ' '
 
      return
      end
      subroutine help

c*********************************************************************72
c
cc HELP prints a list of the available commands.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
      implicit none

      write ( *, * ) ' '
      write ( *, * ) 'This is the list of legal commands.'
      write ( *, * ) ' '
      write ( *, * ) 'H                 Help. List the legal commands.'
      write ( *, * )
     &  'LDA=value         Set leading dimension of arrays.'
      write ( *, * ) 'M                 Multiply two matrices.'
      write ( *, * ) 'N=value           Assigns the size of the arrays.'
      write ( *, * ) 'N=nlo,nhi,ninc    Sets N=nlo, N=nlo+ninc, ....'
      write ( *, * ) 'N=nlo,nhi,*nmult  Sets N=nlo, N=nlo*nmult, ....'
      write ( *, * ) 'NREP=nrep         Sets the repetition factor.'
      write ( *, * ) 'ORDER=name        Chooses the algorithm.'
      write ( *, * ) 'P                 Prints out current results.'
      write ( *, * ) 'Q                 Quit.'
      write ( *, * ) 'SHOW=name         Include "name" in output.'
      write ( *, * )
     &  '                  "name" = ORDER, LDA, N, CPU, OPS,'
      write ( *, * ) '                  MFLOPS, MACHINE, or LANGUAGE.'
      write ( *, * )
     &  '                  If "name"=ALL, all items are shown.'
      write ( *, * ) 'NOSHOW=name       removes item from output list.'
      write ( *, * ) ' '
 
      return
      end
      subroutine ijk ( a, b, c, lda, n, nrep, ttime )

c*********************************************************************72
c
cc IJK multiplies A = B*C using index order IJK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Workspace, real A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
c    used in the multiplication.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
c    Input, integer NREP, the number of times the multiplication should
c    be carried out.
c
c    Output, real TTIME, an estimate of the CPU time in seconds required
c    for the matrix multiplications.
c
      implicit none

      integer lda
      integer n

      real a(lda,n)
      real b(lda,n)
      real c(lda,n)
      integer i
      integer irep
      integer j
      integer k
      integer nrep
      real time1
      real time2
      real ttime

      ttime = 0.0E+00
 
      do irep = 1, nrep
 
        call r4_set ( a, b, c, lda, n )
 
        call matmul_cpu_timer ( time1 )
 
        do i = 1, n
          do j = 1, n
            do k = 1, n
              a(i,k) = a(i,k) + b(i,j) * c(j,k)
            end do
          end do
        end do
 
        call matmul_cpu_timer ( time2 )
 
        ttime = ttime + time2 - time1
 
      end do
 
      return
      end
      subroutine ijuk ( a, b, c, lda, n, nrep, ttime )

c*********************************************************************72
c
cc IJUK multiplies A = B*C using index order IJK and unrolling.
c
c  Discussion:
c
c    The K loop is unrolled to a depth of NROLL = 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Workspace, real A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
c    used in the multiplication.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
c    Input, integer NREP, the number of times the multiplication should
c    be carried out.
c
c    Output, real TTIME, an estimate of the CPU time in seconds required
c    for the matrix multiplications.
c
      implicit none

      integer nroll
      parameter (nroll = 4)

      integer lda
      integer n

      real a(lda,n)
      real b(lda,n)
      real c(lda,n)
      integer i
      integer irep
      integer j
      integer k
      integer khi
      integer nrep
      real time1
      real time2
      real ttime

      ttime = 0.0E+00
 
      do irep = 1, nrep
 
        call r4_set ( a, b, c, lda, n )
 
        call matmul_cpu_timer ( time1 )
 
        khi = ( n / nroll ) * nroll

        do i = 1, n
          do j = 1, n
            do k = 1, khi, nroll
              a(i,k)   = a(i,k)   + b(i,j) * c(j,k)
              a(i,k+1) = a(i,k+1) + b(i,j) * c(j,k+1)
              a(i,k+2) = a(i,k+2) + b(i,j) * c(j,k+2)
              a(i,k+3) = a(i,k+3) + b(i,j) * c(j,k+3)
            end do
          end do
        end do
c
c  Take care of the few cases we missed if N is not a multiple of 4.
c
        do i = 1, n
          do j = 1, n
            do k = khi+1, n
              a(i,k) = a(i,k) + b(i,j) * c(j,k)
            end do
          end do
        end do
 
        call matmul_cpu_timer ( time2 )
        ttime = ttime + time2 - time1
 
      end do
 
      return
      end
      subroutine ikj ( a, b, c, lda, n, nrep, ttime )

c*********************************************************************72
c
cc IKJ multiplies A = B*C using index order IKJ.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Workspace, real A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
c    used in the multiplication.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
c    Input, integer NREP, the number of times the multiplication should
c    be carried out.
c
c    Output, real TTIME, an estimate of the CPU time in seconds required
c    for the matrix multiplications.
c
      implicit none

      integer lda
      integer n

      real a(lda,n)
      real b(lda,n)
      real c(lda,n)
      integer i
      integer irep
      integer j
      integer k
      integer nrep
      real time1
      real time2
      real ttime

      ttime = 0.0E+00
 
      do irep = 1, nrep
 
        call r4_set ( a, b, c, lda, n )
 
        call matmul_cpu_timer ( time1 )
 
        do i = 1, n
          do k = 1, n
            do j = 1, n
              a(i,k) = a(i,k) + b(i,j) * c(j,k)
            end do
          end do
        end do
 
        call matmul_cpu_timer ( time2 )
        ttime = ttime + time2 - time1
 
      end do
 
      return
      end
      subroutine init ( command, cshow, fshow, lda, lena, lingo, lnshow,
     &  lshow, machine, mshow, n, nhi, ninc, nlo, nmult, noshow,
     &  nrep, nrshow, nshow, order, tshow )

c*********************************************************************72
c
cc INIT initializes data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character*6 COMMAND, the most recent user command.
c
c    Output, logical CSHOW, is TRUE if the multiplication method
c    is to be shown.
c
c    Output, logical FSHOW, is TRUE if the MegaFLOP rate is to be shown.
c
c    Output, integer LDA, the leading dimension used for arrays.
c
c    Input, integer LENA, the maximum matrix order allowed.
c
c    Output, character*7 LINGO, the language in which MATMUL is written.
c
c    Output, logical LNSHOW, is TRUE if the software's programming
c    language, stored in LINGO, is to be shown.
c
c    Output, logical LSHOW, is TRUE if the variable LDA is to be shown.
c
c    Output, character*10 MACHINE, the name of the computer on which
c    MATMUL has been compiled.
c
c    Output, logical MSHOW, is TRUE if the machine name is to be shown.
c
c    Output, integer N, the number of rows and columns in the matrices.
c
c    Output, integer NHI, the maximum value of N to use.
c
c    Output, integer NINC, the additive increment to use, if additive
c    steps are being taken.
c
c    Output, integer NLO, the smallest value of N to use.
c
c    Output, integer NMULT, the multiplier to use, if multiplicative
c    steps are being taken.
c
c    Output, logical NOSHOW, is TRUE if the number of operations is
c    to be shown.
c
c    Output, integer NREP, the number of times the multiplication should
c    be carried out.
c
c    Output, logical NRSHOW, is TRUE if the number of repetitions is
c    to be shown.
c
c    Output, logical NSHOW, is TRUE if the matrix size N is to be shown.
c
c    Output, character*6 ORDER, specifies the method to be used.
c
c    Output, logical TSHOW, is TRUE if the execution time is to be shown.
c
      implicit none

      integer lena

      character*6 command
      logical cshow
      logical fshow
      integer lda
      character*7 lingo
      logical lnshow
      logical lshow
      character*10 machine
      logical mshow
      integer n
      integer nhi
      integer ninc
      integer nlo
      integer nmult
      logical noshow
      integer nrep
      logical nrshow
      logical nshow
      character*6 order
      logical tshow

      command = ' '
      cshow = .true.
      fshow = .true.
      lda = lena
      lingo = 'Fortran'
      lnshow = .true.
      lshow = .true.

c     machine = 'CM-2'
c     machine = 'Cray C90'
c     machine = 'Cray YMP'
c     machine = 'DEC Alpha'
c     machine = 'DECstation'
c     machine = 'IBM PC'
c     machine = 'Macintosh'
c     machine = 'Mac/6881'
      machine = 'Mac/G5'
c     machine = 'NEXT'
c     machine = 'SGI/IRIS'
c     machine = 'SUN'
c     machine = 'VAX/VMS'
c     machine = 'VECTOR/VAX'
c
      mshow = .true.
      n = 16
      nhi = n
      ninc = 0
      nlo = n
      nmult = 1
      noshow = .true.
      nrep = 1
      nrshow = .false.
      nshow = .true.
      order = 'IJK'
      tshow = .true.
 
      return
      end
      subroutine i4_set ( a, b, c, lda, n )

c*********************************************************************72
c
cc I4_SET initializes the I4 matrices.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Workspace, integer A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
c    used in the multiplication.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
      implicit none

      integer lda
      integer n

      integer a(lda,n)
      integer b(lda,n)
      integer c(lda,n)
      integer i
      integer j

      do i = 1, n
        do j = 1, n
          a(i,j) = 0
          b(i,j) = 1
          c(i,j) = 1
        end do
      end do
 
      return
      end
      subroutine iujk ( a, b, c, lda, n, nrep, ttime )

c*********************************************************************72
c
cc IUJK multiplies A = B*C using index order IJK, and unrolling on J.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Workspace, real A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
c    used in the multiplication.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
c    Input, integer NREP, the number of times the multiplication should
c    be carried out.
c
c    Output, real TTIME, an estimate of the CPU time in seconds required
c    for the matrix multiplications.
c
      implicit none

      integer nroll
      parameter ( nroll = 4 )

      integer lda
      integer n

      real a(lda,n)
      real b(lda,n)
      real c(lda,n)
      integer i
      integer irep
      integer j
      integer jhi
      integer k
      integer nrep
      real time1
      real time2
      real ttime

      ttime = 0.0E+00
 
      do irep = 1, nrep
 
        call r4_set ( a, b, c, lda, n )
 
        call matmul_cpu_timer ( time1 )
        jhi = ( n / nroll ) * nroll
 
        do i = 1, n
          do j = 1, jhi, nroll
            do k = 1, n
              a(i,k) = a(i,k)
     &          + b(i,j  ) * c(j,k)
     &          + b(i,j+1) * c(j+1,k)
     &          + b(i,j+2) * c(j+2,k)
     &          + b(i,j+3) * c(j+3,k)
            end do
          end do
        end do
c
c  Take care of the few cases we missed if N is not a multiple of 4.
c
        do i = 1, n
          do j = jhi+1, n
            do k = 1, n
              a(i,k) = a(i,k) + b(i,j) * c(j,k)
            end do
          end do
        end do
 
        call matmul_cpu_timer ( time2 )
        ttime = ttime + time2 - time1
 
      end do
 
      return
      end
      subroutine jik ( a, b, c, lda, n, nrep, ttime )

c*********************************************************************72
c
cc JIK multiplies A = B*C using index order JIK.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Workspace, real A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
c    used in the multiplication.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
c    Input, integer NREP, the number of times the multiplication should
c    be carried out.
c
c    Output, real TTIME, an estimate of the CPU time in seconds required
c    for the matrix multiplications.
c
      implicit none

      integer lda
      integer n

      real a(lda,n)
      real b(lda,n)
      real c(lda,n)
      integer i
      integer irep
      integer j
      integer k
      integer nrep
      real time1
      real time2
      real ttime

      ttime = 0.0E+00
 
      do irep = 1, nrep
 
        call r4_set ( a, b, c, lda, n )
 
        call matmul_cpu_timer ( time1 )
 
        do j = 1, n
          do i = 1, n
            do k = 1, n
              a(i,k) = a(i,k) + b(i,j) * c(j,k)
            end do
          end do
        end do
 
        call matmul_cpu_timer ( time2 )
        ttime = ttime + time2 - time1
 
      end do
 
      return
      end
      subroutine jki ( a, b, c, lda, n, nrep, ttime )

c*********************************************************************72
c
cc JKI multiplies A = B*C using index order JKI.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Workspace, real A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
c    used in the multiplication.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
c    Input, integer NREP, the number of times the multiplication should
c    be carried out.
c
c    Output, real TTIME, an estimate of the CPU time in seconds required
c    for the matrix multiplications.
c
      implicit none

      integer lda
      integer n

      real a(lda,n)
      real b(lda,n)
      real c(lda,n)
      integer i
      integer irep
      integer j
      integer k
      integer nrep
      real time1
      real time2
      real ttime

      ttime = 0.0E+00
 
      do irep = 1, nrep
 
        call r4_set ( a, b, c, lda, n )
 
        call matmul_cpu_timer ( time1 )
 
        do j = 1, n
          do k = 1, n
            do i = 1, n
              a(i,k) = a(i,k) + b(i,j) * c(j,k)
            end do
          end do
        end do
 
        call matmul_cpu_timer ( time2 )
        ttime = ttime + time2 - time1
 
      end do
 
      return
      end
      subroutine kij ( a, b, c, lda, n, nrep, ttime )

c*********************************************************************72
c
cc KIJ multiplies A = B*C using index order KIJ.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Workspace, real A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
c    used in the multiplication.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
c    Input, integer NREP, the number of times the multiplication should
c    be carried out.
c
c    Output, real TTIME, an estimate of the CPU time in seconds required
c    for the matrix multiplications.
c
      implicit none

      integer lda
      integer n

      real a(lda,n)
      real b(lda,n)
      real c(lda,n)
      integer i
      integer irep
      integer j
      integer k
      integer nrep
      real time1
      real time2
      real ttime

      ttime = 0.0E+00
 
      do irep = 1, nrep
 
        call r4_set ( a, b, c, lda, n )
 
        call matmul_cpu_timer ( time1 )
 
        do k = 1, n
          do i = 1, n
            do j = 1, n
              a(i,k) = a(i,k) + b(i,j) * c(j,k)
            end do
          end do
        end do
 
        call matmul_cpu_timer ( time2 )
        ttime = ttime + time2 - time1
 
      end do
 
      return
      end
      subroutine kji ( a, b, c, lda, n, nrep, ttime )

c*********************************************************************72
c
cc KJI multiplies A = B*C using index order KJI.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Workspace, real A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
c    used in the multiplication.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
c    Input, integer NREP, the number of times the multiplication should
c    be carried out.
c
c    Output, real TTIME, an estimate of the CPU time in seconds required
c    for the matrix multiplications.
c
      implicit none

      integer lda
      integer n

      real a(lda,n)
      real b(lda,n)
      real c(lda,n)
      integer i
      integer irep
      integer j
      integer k
      integer nrep
      real time1
      real time2
      real ttime

      ttime = 0.0E+00
 
      do irep = 1, nrep
 
        call r4_set ( a, b, c, lda, n )
 
        call matmul_cpu_timer ( time1 )
 
        do k = 1, n
          do j = 1, n
            do i = 1, n
              a(i,k) = a(i,k) + b(i,j) * c(j,k)
            end do
          end do
        end do
 
        call matmul_cpu_timer ( time2 )
        ttime = ttime + time2 - time1
 
      end do
 
      return
      end
      subroutine l4_ijk ( a, b, c, lda, n, nrep, ttime )

c*********************************************************************72
c
cc L4_IJK "multiplies" A = B*C using index order IJK, using logical data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Workspace, logical A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
c    used in the multiplication.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
c    Input, integer NREP, the number of times the multiplication should
c    be carried out.
c
c    Output, real TTIME, an estimate of the CPU time in seconds required
c    for the matrix multiplications.
c
      implicit none

      integer lda
      integer n

      logical a(lda,n)
      logical b(lda,n)
      logical c(lda,n)
      integer i
      integer irep
      integer j
      integer k
      integer nrep
      real time1
      real time2
      real ttime

      ttime = 0.0E+00
 
      do irep = 1, nrep
 
        call l4_set ( a, b, c, lda, n )
 
        call matmul_cpu_timer ( time1 )
 
        do i = 1, n
          do j = 1, n
            do k = 1, n
              a(i,k) = a(i,k) .or. ( b(i,j) .and. c(j,k) )
            end do
          end do
        end do
 
        call matmul_cpu_timer ( time2 )
        ttime = ttime + time2 - time1

      end do
 
      return
      end
      subroutine l4_set ( a, b, c, lda, n )

c*********************************************************************72
c
cc L4_SET initializes the L4 matrices.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, logical A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
c    used in the multiplication.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
      implicit none

      integer lda
      integer n

      logical a(lda,n)
      logical b(lda,n)
      logical c(lda,n)
      integer i
      integer j

      do i = 1, n
        do j = 1, n
          a(i,j) = .false.
          b(i,j) = .true.
          c(i,j) = .true.
        end do
      end do
 
      return
      end
      subroutine mijk ( a, b, c, lda, n, nrep, ttime )

c*********************************************************************72
c
cc MIJK multiplies A = B*C using index order IJK.
c
c  Discussion:
c
c    MIJK uses a Cray directive to run the triple loop using
c    multitasking.  The benefit of such a directive depends on the
c    algorithm and the load on the machine.
c
c    Except on the Cray, this routine should not be used, and in
c    particular, the call to TIMEF should be commented out.
c
c    Note that the Cray routine TIMEF must be called, rather than
c    SECOND.  TIMEF reports elapsed "real" time or "wallclock" time,
c    which should go down with multitasking, whereas CPU time should
c    remain roughly constant.
c
c    In order for parallel processing to occur, this routine must be
c    compiled on the Cray with the directive "-Zu"; moreover, the user must
c    set the environment variable NCPUS to the number of processors the
c    user would like.  For instance, a C shell user would type:
c
c      setenv NCPUS 8
c
c    while a Bourne shell user would type
c
c      NCPUS = 8
c      export NCPUS
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Workspace, real A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
c    used in the multiplication.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
c    Input, integer NREP, the number of times the multiplication should
c    be carried out.
c
c    Output, real TTIME, an estimate of the CPU time in seconds required
c    for the matrix multiplications.
c
      implicit none

      integer lda
      integer n

      real a(lda,n)
      real b(lda,n)
      real c(lda,n)
      integer i
      integer irep
      integer j
      integer k
      integer nrep
      real time1
      real time2
      real ttime

      time1 = 0.0E+00
      time2 = 0.0E+00
      ttime = 0.0E+00
 
      do irep = 1, nrep
 
        call r4_set ( a, b, c, lda, n )
c
c       call timef ( time1 )
c
cMIC$ do GLOBAL
c
        do i = 1, n
          do j = 1, n
            do k = 1, n
              a(i,k) = a(i,k) + b(i,j) * c(j,k)
            end do
          end do
        end do
 
c       call timef ( time2 )
        ttime = ttime + (time2-time1) / 1000.0E+00
 
      end do
 
      return
      end
      subroutine mkji ( a, b, c, lda, n, nrep, ttime )

c*********************************************************************72
c
cc MKJI multiplies A = B*C using index order KJI and multitasking.
c
c  Discussion:
c
c    MKJI uses an SGI/IRIS parallel processing directive to run the
c    triple loop using multitasking.
c
c    The benefit of such a directive depends on the algorithm and the
c    load on the machine.
c
c    Except on the SGI/IRIS, this routine should not be used, and in
c    particular, the call to SECNDS should be commented out.
c
c    Note that the SGI/IRIS routine SECNDS must be called, rather than
c    SECOND.  SECNDS reports elapsed "real" time or "wallclock" time,
c    which should go down with multitasking, whereas CPU time should
c    remain roughly constant.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Workspace, real A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
c    used in the multiplication.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
c    Input, integer NREP, the number of times the multiplication should
c    be carried out.
c
c    Output, real TTIME, an estimate of the CPU time in seconds required
c    for the matrix multiplications.
c
      implicit none

      integer lda
      integer n

      real a(lda,n)
      real b(lda,n)
      real c(lda,n)
      integer i
      integer irep
      integer j
      integer k
      integer nrep
      real time1
      real time2
      real ttime

      time1 = 0.0E+00
      time2 = 0.0E+00
      ttime = 0.0E+00
 
      do irep = 1, nrep
 
        call r4_set ( a, b, c, lda, n )
 
        call matmul_real_timer ( time1 )
c
c$DOACROSS LOCAL(I, J, K)
c
        do k = 1, n
          do j = 1, n
            do i = 1, n
              a(i,k) = a(i,k) + b(i,j) * c(j,k)
            end do
          end do
        end do
 
        call matmul_real_timer ( time2 )

        ttime = ttime + time2 - time1
 
      end do
 
      return
      end
      subroutine mult ( cshow, fshow, lda, lena, lingo, lnshow, lshow, 
     &  machine, mshow, n, noshow, nrep, nrshow, nshow, order, 
     &  output, tshow )

c*********************************************************************72
c
cc MULT carries out the matrix multiplication, using the requested method.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 March 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, logical CSHOW, is TRUE if the multiplication method
c    is to be shown.
c
c    Input, logical FSHOW, is TRUE if the MegaFLOP rate is to be shown.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer LENA, the maximum dimension.
c
c    Input, character*7 LINGO, the language in which MATMUL is written.
c
c    Input, logical LNSHOW, is TRUE if the software's programming
c    language, stored in LINGO, is to be shown.
c
c    Input, logical LSHOW, is TRUE if the variable LDA is to be shown.
c
c    Input, character*10 MACHINE, the computer on which MATMUL has
c    been compiled.
c
c    Input, logical MSHOW, is TRUE if the machine name is to be shown.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
c    Input, logical NOSHOW, is TRUE if the number of operations is
c    to be shown.
c
c    Input, integer NREP, the number of times the multiplication should
c    be carried out.
c
c    Input, logical NRSHOW, is TRUE if the number of repetitions is
c    to be shown.
c
c    Input, logical NSHOW, is TRUE if the matrix size N is to be shown.
c
c    Input, character*6 ORDER, specifies the method to be used.
c
c    Output, character*82 OUTPUT, a string containing the data for
c    this multiplication.
c
c    Input, logical TSHOW, is TRUE if the execution time is to be shown.
c
      implicit none

      integer order_num
      parameter ( order_num = 24 )

      integer lda
      integer n

      logical cshow
      logical fshow
      integer i
      logical s_eqi
      integer lena
      character*7 lingo
      logical lnshow
      logical lshow
      character*10 machine
      logical mshow
      logical noshow
      integer nrep
      logical nrshow
      logical nshow
      character*6 order
      character*6 order2(order_num)
      character*82 output
      logical tshow
      real ttime

      data order2 /
     &  'C4_IJK',
     &  'R8_IJK',
     &  'IJK',
     &  'IJUK',
     &  'IKJ',
     &  'IUJK',
     &  'JIK',
     &  'JKI',
     &  'L4_IJK',
     &  'MIJK',
     &  'MKJI',
     &  'MXMA',
     &  'I4_IJK',
     &  'SAXPYC',
     &  'SAXPYR',
     &  'SDOT',
     &  'SGEMM',
     &  'SGEMMS',
     &  'SIJK',
     &  'TAXPYC',
     &  'TAXPYR',
     &  'TDOT',
     &  'TGEMM',
     &  'UIJK' /

      call header ( cshow, fshow, lnshow, lshow, mshow, noshow,
     &  nrshow, nshow, output, tshow )

      if ( s_eqi ( order, 'ALL' ) ) then

        do i = 1, order_num

          order = order2(i)

          call domethod ( lda, lena, n, nrep, order, ttime )

          call report ( cshow, fshow, lda, lingo, lnshow, lshow,  
     &      machine, mshow, n, noshow, nrep, nrshow, nshow, order, 
     &      output, tshow, ttime )
 
        end do

        order = 'ALL'

      else

        call domethod ( lda, lena, n, nrep, order, ttime )

        call report ( cshow, fshow, lda, lingo, lnshow, lshow,  
     &    machine, mshow, n, noshow, nrep, nrshow, nshow, order, 
     &    output, tshow, ttime )

      end if

      return
      end
      subroutine i4_ijk ( a, b, c, lda, n, nrep, ttime )
cDIR$ integer = 64
c
c  The above line is a Cray compiler directive, which requests that
c  integers be stored as 64 bit quantities, and that 64 bit integer
c  arithmetic be used.
c
c*********************************************************************72
c
cc I4_IJK uses I4 arithmetic and IJK order.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Workspace, integer A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
c    used in the multiplication.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
c    Input, integer NREP, the number of times the multiplication should
c    be carried out.
c
c    Output, real TTIME, an estimate of the CPU time in seconds required
c    for the matrix multiplications.
c
      implicit none

      integer lda
      integer n

      integer a(lda,n)
      integer b(lda,n)
      integer c(lda,n)
      integer i
      integer irep
      integer j
      integer k
      integer nrep
      real time1
      real time2
      real ttime

      ttime = 0.0E+00
 
      do irep = 1, nrep
 
        call i4_set ( a, b, c, lda, n )

        call matmul_cpu_timer ( time1 )
 
        do i = 1, n
          do j = 1, n
            do k = 1, n
              a(i,k) = a(i,k) + b(i,j) * c(j,k)
            end do
          end do
        end do
 
        call matmul_cpu_timer ( time2 )
        ttime = ttime + time2 - time1

      end do
 
      return
      end
      subroutine nstep ( ido, n, nhi, ninc, nlo, nmult )

c*********************************************************************72
c
cc NSTEP is used when a set of values of N is being generated.
c
c  Discussion:
c
c    The routine checks whether addition or multiplication is being used,
c    and increases the value of N.  It also checks whether the set of
c    values is done, or whether the input values are inconsistent.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IDO.
c
c    0, N has reached or surpassed NHI, and no further
c    increase should be carried out.  N has been reset to NLO.
c
c    1, N had not reached or surpassed NHI, and so N has been
c    incremented by NINC, or multiplied by NMULT.
c
c    2, N had not reached or surpassed NHI, but we're never
c    going to get there!  Either:
c
c      NINC is negative but NHI is greater than NLO, or
c
c      NINC is positive but NHI is less than NLO, or
c
c      NMULT is positive but NLO is greater than NHI.
c
c    3, NINC is 0 and NMULT is less than or equal to 1.
c
c    Input/output, integer N, the number of rows and columns in the matrices.
c    This routine modifies the value of N as appropriate.
c
c    Input, integer NHI, the maximum value of N to use.
c
c    Input, integer NINC, the additive increment to use, if additive
c    steps are being taken.
c
c    Input, integer NLO, the smallest value of N to use.
c
c    Input, integer NMULT, the multiplier to use, if multiplicative
c    steps are being taken.
c
      implicit none

      integer ido
      integer n
      integer nhi
      integer ninc
      integer nlo
      integer nmult
c
c  If NINC is not 0, then
c    if it's pointing in the right direction, then
c      add NINC to N,
c      set a continuation flag
c      if N+NINC exceeds NHI, then
c        reset N, and
c        set a completion flag
c    else
c      set an error flag.
c
      if ( ninc .ne. 0 ) then
        if ( ( nlo .lt. nhi .and. ninc .gt. 0 )  .or. 
     &     ( nlo .gt. nhi .and. ninc .lt. 0 ) ) then
          n = n + ninc
          ido = 1
          if ( ( n .gt. nhi .and. nhi .ge. nlo )  .or. 
     &       ( n .lt. nhi .and. nhi .le. nlo ) ) then
            ido = 0
            n = nlo
          end if
        else
          ido = 2
        end if
 
        return
      end if
c
c  If NMULT is greater than 1, then
c    if it's pointing in the right direction, then
c      multiply N by NMULT,
c      set a continuation flag
c      if N*NMULT exceeds NHI, then
c        reset N, and
c        set a completion flag
c    else
c      set an error flag.
c
      if ( nmult .gt. 1 ) then
        if ( nlo .lt. nhi ) then
          n = n * nmult
          ido = 1
          if ( 
     &      ( n .gt. nhi .and. nhi .ge. nlo )  .or. 
     &      ( n .lt. nhi .and. nhi .le. nlo ) ) then
            ido = 0
            n = nlo
          end if
        else
          ido = 2
        end if
 
        return

      end if
c
c  NINC was 0, and NMULT wasn't greater than 1.
c
      ido = 3
 
      return
      end
      subroutine printr ( lda, lena, lingo, n, nhi, ninc, nlo, nmult, 
     &  nrep, order )

c*********************************************************************72
c
cc PRINTR prints out those parameters the user wants to see.
c
c  Discussion:
c
c    These parameters include:
c
c      the language MATMUL is written in,
c      the algorithm,
c      the leading dimension,
c      the maximum allowable dimension,
c      the actual size of arrays,
c      the number of multiplications carried out.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer LENA, the maximum matrix order allowed.
c
c    Input, character*7 LINGO, the language in which MATMUL was written.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
c    Input, integer NHI, the maximum value of N to use.
c
c    Input, integer NINC, the additive increment to use, if additive
c    steps are being taken.
c
c    Input, integer NLO, the smallest value of N to use.
c
c    Input, integer NMULT, the multiplier to use, if multiplicative
c    steps are being taken.
c
c    Input, integer NREP, the number of times the multiplication should
c    be carried out.
c
c    Input, character*6 ORDER, specifies the method to be used.
c
      implicit none

      integer lda
      integer lena
      character*7 lingo
      integer n
      integer nhi
      integer ninc
      integer nlo
      integer nmult
      integer nrep
      character*6 order

      write ( *, * ) ' '
      write ( *, * ) 'This version of MATMUL is written in ' // lingo
      write ( *, * ) 'The algorithm chosen is ' // order
      write ( *, * ) 'The leading dimension of arrays, lda, is ', lda
      write ( *, * ) 'The maximum legal choice for LDA is ', lena
      write ( *, * ) 'The actual size of the arrays, N, is ', n

      if ( nhi .ne. nlo ) then
        write ( *, * ) 'Several problem sizes will be solved in order.'
        write ( *, * ) 'The final size of arrays, NHI, will be ', nhi
      end if
 
      if ( ninc .ne. 0 ) then
        write ( *, * ) 'Array size will be incremented by NINC = ', ninc
      end if
 
      if ( nmult .ne. 1 ) then
        write ( *, * ) 'Array size will be multiplied by NMULT = ', 
     &    nmult
      end if
 
      if ( nrep .ne. 1 ) then
        write ( *, * ) 'Multiplications repeated NREP = ', 
     &    nrep, ' times.'
      end if
 
      return
      end
      subroutine r4_set ( a, b, c, lda, n )

c*********************************************************************72
c
cc R4_SET initializes the R4 matrices.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, real A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
c    used in the multiplication.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
      implicit none

      integer lda
      integer n

      real a(lda,n)
      real b(lda,n)
      real c(lda,n)
      integer i
      integer j

      do i = 1, n
        do j = 1, n
          a(i,j) = 0.0E+00
          b(i,j) = 1.0E+00
          c(i,j) = 1.0E+00
        end do
      end do
 
      return
      end
      subroutine r8_ijk ( a, b, c, lda, n, nrep, ttime )

c*********************************************************************72
c
cc R8_IJK using R8 arithmetic and IJK order.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Workspace, double precision A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
c    used in the multiplication.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
c    Input, integer NREP, the number of times the multiplication should
c    be carried out.
c
c    Output, real TTIME, an estimate of the CPU time in seconds required
c    for the matrix multiplications.
c
      implicit none

      integer lda
      integer n

      double precision a(lda,n)
      double precision b(lda,n)
      double precision c(lda,n)
      integer i
      integer irep
      integer j
      integer k
      integer nrep
      real time1
      real time2
      real ttime

      ttime = 0.0E+00
 
      do irep = 1, nrep
 
        call r8_set ( a, b, c, lda, n )
 
        call matmul_cpu_timer ( time1 )
 
        do i = 1, n
          do j = 1, n
            do k = 1, n
              a(i,k) = a(i,k) + b(i,j) * c(j,k)
            end do
          end do
        end do
 
        call matmul_cpu_timer ( time2 )
        ttime = ttime + time2 - time1

      end do
 
      return
      end
      subroutine r8_set ( a, b, c, lda, n )

c*********************************************************************72
c
cc R8_SET initializes the R8 matrices.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, double precision A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
c    used in the multiplication.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
      implicit none

      integer lda
      integer n

      double precision a(lda,n)
      double precision b(lda,n)
      double precision c(lda,n)
      integer i
      integer j

      do i = 1, n
        do j = 1, n
          a(i,j) = 0.0D+00
          b(i,j) = 1.0D+00
          c(i,j) = 1.0D+00
        end do
      end do
 
      return
      end
      subroutine report ( cshow, fshow, lda, lingo, lnshow, lshow,  
     &  machine, mshow, n, noshow, nrep, nrshow, nshow, order, output, 
     &  tshow, ttime )

c*********************************************************************72
c
cc REPORT reports the results for each multiplication experiment.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, logical CSHOW, is TRUE if the multiplication method
c    is to be shown.
c
c    Input, logical FSHOW, is TRUE if the MegaFLOP rate is to be shown.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, character*7 LINGO, the language in which MATMUL is written.
c
c    Input, logical LNSHOW, is TRUE if the software's programming
c    language, stored in LINGO, is to be shown.
c
c    Input, logical LSHOW, is TRUE if the variable LDA is to be shown.
c
c    Input, character*10 MACHINE, the computer on which MATMUL has
c    been compiled.
c
c    Input, logical MSHOW, is TRUE if the machine name is to be shown.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
c    Input, logical NOSHOW, is TRUE if the number of operations is
c    to be shown.
c
c    Input, integer NREP, the number of times the multiplication should
c    be carried out.
c
c    Input, logical NRSHOW, is TRUE if the number of repetitions is
c    to be shown.
c
c    Input, logical NSHOW, is TRUE if the matrix size N is to be shown.
c
c    Input, character*6 ORDER, specifies the method to be used.
c
c    Output, character*82 OUTPUT, a string containing the data for
c    this multiplication.
c
c    Input, logical TSHOW, is TRUE if the execution time is to be shown.
c
c    Output, real TTIME, an estimate of the CPU time in seconds required
c    for the matrix multiplications.
c
      implicit none

      integer lda
      integer n

      logical cshow
      logical fshow
      real ftemp
      integer ihi
      integer ilo
      character*7 lingo
      logical lnshow
      logical lshow
      character*10 machine
      logical mshow
      logical noshow
      integer nrep
      logical nrshow
      logical nshow
      integer ntemp
      character*6 order
      character*82 output
      logical tshow
      real ttime

      output = ' '
      ihi = 0

      if ( cshow ) then
        ilo = ihi + 1
        ihi = ilo + 5
        output(ilo:ihi) = order
      end if

      if ( lshow ) then
        ilo = ihi + 1
        ihi = ilo + 3
        write ( output(ilo:ihi), '(i4)' ) lda
      end if

      if ( nshow ) then
        ilo = ihi + 1
        ihi = ilo + 3
        write ( output(ilo:ihi), '(i4)' ) n
      end if

      if ( tshow ) then
        ilo = ihi + 1
        ihi = ilo + 13
        write ( output(ilo:ihi), '(f14.6)' ) ttime
      end if

      if ( noshow ) then
        ntemp = 2 * n**3
        ilo = ihi + 1
        ihi = ilo + 9
        write ( output(ilo:ihi), '(i10)' ) ntemp
      end if

      if ( nrshow ) then
        ilo = ihi + 1
        ihi = ilo + 4
        write ( output(ilo:ihi), '(i5)' ) nrep
      end if

      if ( fshow ) then

        if ( ttime .eq. 0.0E+00 ) then
          ftemp = 0.0E+00
        else
          ftemp = real ( ntemp * nrep ) / ( 1.0E+06 * ttime )
        end if

        ilo = ihi + 1
        ihi = ilo + 9
        write ( output(ilo:ihi), '(f10.4)' ) ftemp

      end if

      if ( mshow ) then
        ilo = ihi + 1
        ilo = ilo + 1
        ihi = ilo + 9
        output(ilo:ihi) = machine
      end if

      if ( lnshow ) then
        ilo = ihi + 1
        ilo = ilo + 1
        ihi = ilo + 6
        output(ilo:ihi) = lingo
      end if

      if ( ihi .gt. 0 ) then
        write ( *, '(a)' ) output(1:ihi)
      end if
 
      return
      end
      subroutine s_cap ( string )

c*********************************************************************72
c
cc S_CAP replaces any lowercase letters by uppercase ones in a string.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    16 May 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character*(*) STRING, the string to be transformed.
c
      implicit none

      character*1 c
      integer i
      integer nchar
      character*(*) string

      nchar = len ( string )

      do i = 1, nchar

        c = string(i:i)
        call ch_cap ( c )
        string(i:i) = c

      end do

      return
      end
      function s_eqi ( strng1, strng2 )

c*********************************************************************72
c
cc S_EQI is a case insensitive comparison of two strings for equality.
c
c  Example:
c
c    S_EQI ( 'Anjana', 'ANJANA' ) is .TRUE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    14 April 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) STRNG1, STRNG2, the strings to compare.
c
c    Output, logical S_EQI, the result of the comparison.
c
      implicit none

      integer i
      integer len1
      integer len2
      integer lenc
      logical s_eqi
      character*1 s1
      character*1 s2
      character*(*) strng1
      character*(*) strng2

      len1 = len ( strng1 )
      len2 = len ( strng2 )
      lenc = min ( len1, len2 )

      s_eqi = .false.

      do i = 1, lenc

        s1 = strng1(i:i)
        s2 = strng2(i:i)
        call ch_cap ( s1 )
        call ch_cap ( s2 )

        if ( s1 .ne. s2 ) then
          return
        end if

      end do

      do i = lenc + 1, len1
        if ( strng1(i:i) .ne. ' ' ) then
          return
        end if
      end do

      do i = lenc + 1, len2
        if ( strng2(i:i) .ne. ' ' ) then
          return
        end if
      end do

      s_eqi = .true.

      return
      end
      function s_length ( string )

c*********************************************************************72
c
cc S_LENGTH returns the length of a string up to the last nonblank.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    14 April 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) STRING, the string to be measured.
c
c    Output, integer S_LENGTH, the location of the last nonblank in STRING.
c
      integer i
      integer s_length
      character*(*) string

      do i = len ( string ), 1, -1

        if ( string(i:i) .ne. ' ' ) then
          s_length = i
          return
        end if

      end do

      s_length = 0

      return
      end
      subroutine matmul_cpu_timer ( cpu )

c*********************************************************************72
c
cc MATMUL_CPU_TIMER computes total CPU seconds.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, real CPU, the total CPU time, in seconds, since the
c    program began running.
c
      implicit none

      real cpu
c
c
c  Cray
c
c
c     cpu = second ( )
c
c
c  VAX/VMS
c
c
c     integer*4 code
c     parameter ( code = 2 )
c
c     integer*4 addr
c     integer*4 cputim
c     real oldtim
c     real temp
c
c     save addr
c     save oldtim
c
c     data addr / 0 /
c     data oldtim / 0.0E+00 /
c
c     if ( addr .eq. 0 ) then
c       call lib$init_timer ( addr )
c     end if
c
c     call lib$stat_timer ( code, cputim, addr )
c     temp = real ( cputim ) / 100.0E+00
c
c     oldtim = temp
c     cpu = temp
c
c
c  UNIX systems
c
c
      real tarray(2)
      call etime ( tarray, cpu )
c
c
c  PowerMac using Absoft FORTRAN.
c
c
c     include "OSUtils.inc"
c
c     RECORD /DateTimeRec/ DateTime
c
c     call GetTime ( DateTime )
c
c     cpu = 3600 * DateTime.hour + 60 * DateTime.minute 
c    &  + DateTime.second
c
c
c  Apple Macintosh using Absoft FORTRAN.
c
c
c     integer isecnd
c
c     call time ( isecnd )
c     cpu = real ( isecnd )
c
c
c  Apple Macintosh using LS FORTRAN.
c
c
c     cpu = secnds ( 0.0E+00 )
c
c
c  IBM PC using Microsoft FORTRAN.
c
c
c     integer*2 i100th
c     integer*2 ihr
c     integer*2 imin
c     integer*2 isec
c     real oldtim
c     real temp
c
c     save oldtim
c
c     call gettim ( ihr, imin, isec, i100th )
c
c     temp = 3600 * ihr + 60 * imin + isec + real ( i100th ) / 100.0E+00
c
c     oldtim = temp
c     cpu = temp
c
      return
      end
      subroutine matmul_real_timer ( seconds )

c*********************************************************************72
c
cc MATMUL_REAL_TIMER returns a reading of the real time clock.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, real SECONDS, a current real time in seconds.
c
      implicit none

      real seconds

      seconds = 0.0E+00
c
c  Cray code:
c
c     seconds = secnds ( 0.0E+00 )
c
      return
      end
      subroutine sijk ( a, b, c, lda, n, nrep, ttime )

c*********************************************************************72
c
cc SIJK multiplies A = B*C using index order IJK, and no Cray vectorization.
c
c  Discussion:
c
c    SIJK uses a Cray directive to run the inner do loop WITHOUT
c    vectorization.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Workspace, real A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
c    used in the multiplication.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
c    Input, integer NREP, the number of times the multiplication should
c    be carried out.
c
c    Output, real TTIME, an estimate of the CPU time in seconds required
c    for the matrix multiplications.
c
      implicit none

      integer lda
      integer n

      real a(lda,n)
      real b(lda,n)
      real c(lda,n)
      integer i
      integer irep
      integer j
      integer k
      integer nrep
      real time1
      real time2
      real ttime

      ttime = 0.0E+00
 
      do irep = 1, nrep
 
        call r4_set ( a, b, c, lda, n )
 
        call matmul_cpu_timer ( time1 )
 
        do i = 1, n
          do j = 1, n
cDIR$ NEXTSCALAR
            do k = 1, n
              a(i,k) = a(i,k) + b(i,j) * c(j,k)
            end do
          end do
        end do
 
        call matmul_cpu_timer ( time2 )
        ttime = ttime + time2 - time1
 
      end do
 
      return
      end
      subroutine taxpy ( n, sa, sx, incx, sy, incy )

c*********************************************************************72
c
cc TAXPY is unoptimized standard BLAS routine SAXPY.
c
c  Discussion:
c
c    Roughly, TAXPY adds SA * SX(I) to SY(I) for I = 1 to N.
c    However, the increments INCX and INCY allow this to be done
c    even when SX or SY is a row or column of a matrix.  
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Parameters:
c
c    Input, integer N, the "logical" number of items in the vectors.
c
c    Input, real SA, the multiplier.
c
c    Input, real SX(*), a vector, a multiple of which is to be added to SY.
c
c    Input, integer INCX, the increment in SX between the successive
c    elements that we will use.
c
c    Input/output, real SY(*), a vector, to which is to be added SA 
c    times an entry of SX.
c
c    Input, integer INCY, the increment in SY between the successive
c    elements that we will use.
c
      implicit none

      integer i
      integer incx
      integer incy
      integer ix
      integer iy
      integer m
      integer n
      real sa
      real sx(*)
      real sy(*)

      if (n.le.0) return
 
      if (sa.eq.0.0E+00)return
 
      if (incx.ne.1 .or. incy.ne.1 ) then
 
        if (incx.lt.0 ) then
          ix = (-n+1)*incx+1
        else
          ix = 1
        end if
 
        if (incy.lt.0 ) then
          iy = (-n+1)*incy+1
        else
          iy = 1
        end if
 
        do i = 1, n
          sy(iy) = sy(iy)+sa*sx(ix)
          ix = ix+incx
          iy = iy+incy
        end do
 
      else
 
        m = mod(n,4)
 
        do i = 1,m
          sy(i) = sy(i)+sa*sx(i)
        end do
 
        do i = m+1,n,4
          sy(i) = sy(i)+sa*sx(i)
          sy(i+1) = sy(i+1)+sa*sx(i+1)
          sy(i+2) = sy(i+2)+sa*sx(i+2)
          sy(i+3) = sy(i+3)+sa*sx(i+3)
        end do
 
      end if
 
      return
      end
      function tdot ( n, sx, incx, sy, incy )

c*********************************************************************72
c
cc TDOT computes the inner product of two vectors.
c
c  Discussion:
c
c    TDOT is a source code version of the BLAS level 1 routine SDOT,
c    which can be used to compare performance with optimized versions
c    of SDOT supplied by the compiler or operating system.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Parameters:
c
c    Input, integer N, the "logical" number of items in the vectors.
c
c    Input, real SX(*), the first vector.
c
c    Input, integer INCX, the increment in SX between the successive
c    elements that we will use.
c
c    Input/output, real SY(*), the second vector.
c
c    Input, integer INCY, the increment in SY between the successive
c    elements that we will use.
c
c    Output, real TDOT, the sum of the products of the appropriate
c    entries of SX and SY.
c
      implicit none

      integer i
      integer incx
      integer incy
      integer ix
      integer iy
      integer m
      integer n
      real stemp
      real sx(*)
      real sy(*)
      real tdot

      tdot = 0.0E+00
 
      if ( n .le. 0 ) then
        return
      end if
 
      if ( incx .ne. 1 .or. incy .ne. 1 ) then

        if ( incx .lt. 0 ) then
          ix = ( - n + 1 ) * incx + 1
        else
          ix = 1
        end if
 
        if ( incy .lt. 0 ) then
          iy = ( - n + 1 ) * incy + 1
        else
          iy = 1
        end if
 
        stemp = 0.0E+00
        do i = 1, n
          stemp = stemp + sx(ix) * sy(iy)
          ix = ix + incx
          iy = iy + incy
        end do
 
        tdot = stemp
 
      else
 
        m = mod ( n, 5)
 
        stemp = 0.0E+00
 
        do i = 1, m
          stemp = stemp 
     &      + sx(i) * sy(i)
        end do
 
        do i = m+1, n, 5
          stemp = stemp
     &      + sx(i)   * sy(i)
     &      + sx(i+1) * sy(i+1)
     &      + sx(i+2) * sy(i+2)
     &      + sx(i+3) * sy(i+3)
     &      + sx(i+4) * sy(i+4)
        end do
 
      end if
 
      tdot = stemp
 
      return
      end
      subroutine tgemm(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,
     & ldc)

c*********************************************************************72
c
cc TGEMM is a source code copy of SGEMM, a BLAS matrix * matrix routine.
c
c  TGEMM performs one of the matrix-matrix operations
c
c     C : =  alpha*op( A )*op( B ) + beta*C,
c
c  where  op( X ) is one of
c
c     op( X ) = X   or   op( X ) = X',
c
c  alpha and beta are scalars, anda, B and C are matrices, with op( A )
c  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
c
c  Parameters:
c
c  transa - character*1.
c           On entry, transa specifies the form of op( A ) to be used in
c           the matrix multiplication as follows:
c
c              transa = 'N' or 'n',  op( A ) = A.
c
c              transa = 'T' or 't',  op( A ) = A'.
c
c              transa = 'C' or 'c',  op( A ) = A'.
c
c           Unchanged on exit.
c
c  transb - character*1.
c           On entry, transb specifies the form of op( B ) to be used in
c           the matrix multiplication as follows:
c
c              transa = 'N' or 'n',  op( B ) = B.
c
c              transa = 'T' or 't',  op( B ) = B'.
c
c              transa = 'C' or 'c',  op( B ) = B'.
c
c           Unchanged on exit.
c
c  m      - integer.
c           On entry,  m  specifies  the number  of rows  of the  matrix
c           op( A )  and of the  matrix  C.  m  must  be at least  zero.
c           Unchanged on exit.
c
c  N      - integer.
c           On entry,  N  specifies the number  of columns of the matrix
c           op( B ) and the number of columns of the matrix C. N must be
c           at least zero.
c           Unchanged on exit.
c
c  K      - integer.
c           On entry,  K  specifies  the number of columns of the matrix
c           op( A ) and the number of rows of the matrix op( B ). K must
c           be at least  zero.
c           Unchanged on exit.
c
c  alpha  - real .
c           On entry, alpha specifies the scalar alpha.
c           Unchanged on exit.
c
c  A      - real array of DIMENSION ( lda, ka ), where ka is
c           k  when  transa = 'N' or 'n',  and is  m  otherwise.
c           Before entry with  transa = 'N' or 'n',  the leading  m by k
c           part of the array  A  must contain the matrix a,  otherwise
c           the leading  k by m  part of the array  A  must contain  the
c           matrix A.
c           Unchanged on exit.
c
c  lda    - integer.
c           On entry, lda specifies the first dimension of A as declared
c           in the calling (sub) program. When  transa = 'N' or 'n' then
c           lda must be at least  max( 1, m ), otherwise  lda must be at
c           least  max( 1, k ).
c           Unchanged on exit.
c
c  B      - real array of DIMENSION ( ldb, kb ), where kb is
c           n  when  transb = 'N' or 'n',  and is  k  otherwise.
c           Before entry with  transb = 'N' or 'n',  the leading  k by n
c           part of the array  B  must contain the matrix  B,  otherwise
c           the leading  n by k  part of the array  B  must contain  the
c           matrix B.
c           Unchanged on exit.
c
c  ldb    - integer.
c           On entry, ldb specifies the first dimension of B as declared
c           in the calling (sub) program. When  transb = 'N' or 'n' then
c           ldb must be at least  max( 1, k ), otherwise  LDB must be at
c           least  max( 1, n ).
c           Unchanged on exit.
c
c  beta   - real .
c           On entry,  beta  specifies the scalar  beta.  When  BETA  is
c           supplied as zero then C need not be set on input.
c           Unchanged on exit.
c
c  c      - real array of DIMENSION ( LDC, n ).
c           Before entry, the leading  m by n  part of the array  C must
c           contain the matrix  C,  except when  beta  is zero, in which
c           case C need not be set on entry.
c           On exit, the array  C  is overwritten by the  m by n  matrix
c           ( alpha*op( A )*op( B ) + beta*C ).
c
c  LDC    - integer.
c           On entry, LDC specifies the first dimension of C as declared
c           in  the  calling  (sub)  program.   LDC  must  be  at  least
c           max( 1, m ).
c           Unchanged on exit.
c
      implicit none

      real alpha, beta
      integer m, n, k, lda, ldb, ldc
      character*1 transa, transb
      real a(lda,*),b(ldb,*),c(ldc,*)

      real one, zero
      parameter ( one  = 1.0E+0, zero = 0.0E+0 )

      integer i, info, j, ncola, nrowa, nrowb
      logical nota, notb

      logical tlsame
c
c     Set  nota  and  notb  as  true if  A  and  B  respectively are not
c     transposed, and set  nrowa, ncola and  nrowb as the number of rows
c     and  columns  of  A  and the  number of rows  of  B  respectively.
c
      nota = tlsame( transa, 'N' )
      notb = tlsame( transb, 'N' )
      if ( nota  ) then
         nrowa = m
         ncola = k
      else
         nrowa = k
         ncola = m
      end if
      if ( notb  ) then
         nrowb = k
      else
         nrowb = n
      end if
c
c     Test the input parameters.
c
      info = 0
      if (      ( .not.nota                 ).and.
     &         ( .not.tlsame( transa, 'T' ) ).and.
     &         ( .not.tlsame( transa, 'C' ) )       ) then
         info = 1
      else if ( ( .not.notb                 ).and.
     &         ( .not.tlsame( transb, 'T' ) ).and.
     &         ( .not.tlsame( transb, 'C' ) )       ) then
         info = 2
      else if ( m.lt.0  ) then
         info = 3
      else if ( n.lt.0  ) then
         info = 4
      else if ( k.lt.0  ) then
         info = 5
      else if ( lda.lt.max( 1,nrowa )  ) then
         info = 8
      else if ( ldb.lt.max( 1,nrowb )  ) then
         info = 10
      else if ( ldc.lt.max( 1, m )  ) then
         info = 13
      end if
 
      if ( info .ne. 0  ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TGEMM - Fatal error!'
        write ( *, '(a,i8)' ) 
     &    '  Illegal value for input argument #', info
        return
      end if
c
c  Quick return if possible.
c
      if ( 
     &  m .eq. 0 .or. 
     &  n.eq.0 .or.
     &  ( ( alpha.eq.zero .or. k.eq.0 ) .and. beta.eq.one ) ) then
        return
      end if
c
c  Start the operations.
c
      if ( k.eq.0  ) then
c
c  Form  C : =  beta*C.
c
         if ( beta.eq.zero  ) then
           do j = 1, n
             do i = 1,m
               c(i,j) = zero
             end do
           end do
         else
           do j = 1, n
             do i = 1,m
               c(i,j) = beta*c(i,j)
             end do
           end do
         end if
      else if ( notb  ) then
c
c  Form  C : =  alpha*op( A )*B + beta*C.
c
        do j = 1, n
          call tgemvf(transa, nrowa, ncola,
     &                   alpha,a, lda, b( 1,j), 1,
     &                   beta, c( 1,j), 1 )
        end do
 
      else
c
c  Form  C : =  alpha*op( A )*B' + beta*C.
c
        do j = 1, n
          call tgemvf(transa,nrowa,ncola,alpha,a,lda,b(j,1),ldb,
     &      beta,c(1,j),1)
        end do
 
      end if
 
      return
      end
      subroutine tgemvf(trans,m,n,alpha,a,lda,x,incx,beta,y,incy)

c*********************************************************************72
c
cc TGEMVF is a source code copy of BLAS SGEMVF, a matrix * vector routine.
c
c  TGEMVF performs one of the matrix-vector operations
c
c     y : =  alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
c
c  where alpha and beta are scalars, x and y are vectors and A is an
c  m by n matrix.
c
c  Parameters:
c
c  trans  - character*1.
c           On entry, trans specifies the operation to be performed as
c           follows:
c
c              trans = 'N' or 'n'   y := alpha*A*x + beta*y.
c
c              trans = 'T' or 't'   y := alpha*A'*x + beta*y.
c
c              trans = 'C' or 'c'   y := alpha*A'*x + beta*y.
c
c           Unchanged on exit.
c
c  m      - integer.
c           On entry, m specifies the number of rows of the matrix A.
c           m must be at least zero.
c           Unchanged on exit.
c
c  N      - integer.
c           On entry, N specifies the number of columns of the matrix A.
c           N must be at least zero.
c           Unchanged on exit.
c
c  alpha  - real .
c           On entry, alpha specifies the scalar alpha.
c           Unchanged on exit.
c
c  A      - real array of DIMENSION ( lda, n ).
c           Before entry, the leading m by n part of the array A must
c           contain the matrix of coefficients.
c           Unchanged on exit.
c
c  lda    - integer.
c           On entry, lda specifies the first dimension of A as declared
c           in the calling (sub) program. lda must be at least
c           max( 1, m ).
c           Unchanged on exit.
c
c  X      - real array of DIMENSION at least
c           ( 1 + ( n - 1 )*abs( incx ) ) when trans = 'N' or 'n'
c           and at least
c           ( 1 + ( m - 1 )*abs( incx ) ) otherwise.
c           Before entry, the incremented array X must contain the
c           vector x.
c           Unchanged on exit.
c
c  incx   - integer.
c           On entry, incx specifies the increment for the elements of
c           X. incx must not be zero.
c           Unchanged on exit.
c
c  beta   - real .
c           On entry, beta specifies the scalar beta. When BETA is
c           supplied as zero then Y need not be set on input.
c           Unchanged on exit.
c
c  Y      - real array of DIMENSION at least
c           ( 1 + ( m - 1 )*abs( incy ) ) when trans = 'N' or 'n'
c           and at least
c           ( 1 + ( n - 1 )*abs( incy ) ) otherwise.
c           Before entry with beta non-zero, the incremented array Y
c           must contain the vector y. On exit, Y is overwritten by the
c           updated vector y.
c
c  incy   - integer.
c           On entry, incy specifies the increment for the elements of
c           Y. incy must not be zero.
c           Unchanged on exit.
c
      implicit none

      real alpha, beta
      integer incx, incy, lda, m, n
      character*1 trans

      real a( lda, * ), x( * ), y( * )

      real one, zero
      parameter( one = 1.0E+0, zero = 0.0E+0 )

      real temp
      integer i, info, ix, iy, j, jx, jy, kx, ky, lenx, leny

      logical tlsame
c
c  Test the input parameters.
c
      info = 0
      if ( .not.tlsame( trans, 'N' ).and.
     &         .not.tlsame( trans, 'T' ).and.
     &         .not.tlsame( trans, 'C' )       ) then
         info = 1
      else if ( m.lt.0  ) then
         info = 2
      else if ( n.lt.0  ) then
         info = 3
      else if ( lda.lt.max( 1, m )  ) then
         info = 6
      else if ( incx.eq.0  ) then
         info = 8
      else if ( incy.eq.0  ) then
         info = 11
      end if

      if ( info .ne. 0  ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'TGEMVF - Fatal error!'
        write ( *, '(a,i8)' ) 
     &    '  Illegal value for input argument #', info
         return
      end if
c
c  Quick return if possible.
c
      if ( ( m.eq.0 ) .or. ( n.eq.0 ).or.
     &    ( ( alpha.eq.zero ).and.( beta.eq.one ) ) )
     &   return
c
c  Set  lenx  and  leny, the lengths of the vectors x and y, and set
c  up the start points in  X  and  Y.
c
      if ( tlsame( trans, 'N' )  ) then
         lenx = n
         leny = m
      else
         lenx = m
         leny = n
      end if
      if ( incx.gt.0  ) then
         kx = 1
      else
         kx = 1 - ( lenx - 1 )*incx
      end if
      if ( incy.gt.0  ) then
         ky = 1
      else
         ky = 1 - ( leny - 1 )*incy
      end if
c
c  Start the operations. In this version the elements of A are
c  accessed sequentially with one pass through A.
c
c  First form  y : =  beta*y.
c
      if ( beta.ne.one  ) then
         if ( incy.eq.1  ) then
            if ( beta.eq.zero  ) then
              do i = 1,leny
                y(i) = zero
              end do
            else
              do i = 1, leny
                y(i) = beta*y(i)
              end do
            end if
         else
            iy = ky
            if ( beta.eq.zero  ) then
              do i = 1, leny
                y(iy) = zero
                iy      = iy   + incy
              end do
            else
              do i = 1,leny
                y(iy) = beta*y(iy)
                iy      = iy           + incy
               end do
            end if
         end if
      end if
      if ( alpha.eq.zero )
     &   return
      if ( tlsame( trans, 'N' )  ) then
c
c  Form  y : =  alpha*A*x + y.
c
         jx = kx
         if ( incy.eq.1  ) then
            do j = 1, n
               if ( x( jx ).ne.zero  ) then
                  temp = alpha*x( jx )
                  do i = 1, m
                     y(i) = y( i ) + temp*a(i,j)
                  end do
               end if
               jx = jx + incx
            end do
         else
            do j = 1, n
               if ( x( jx ).ne.zero  ) then
                  temp = alpha*x( jx )
                  iy   = ky
                  do i = 1, m
                     y(iy) = y(iy) + temp*a(i,j)
                     iy      = iy      + incy
                 end do
               end if
               jx = jx + incx
            end do
         end if
      else
c
c  Form  y : =  alpha*A'*x + y.
c
         jy = ky
         if ( incx.eq.1  ) then
           do j = 1, n
               temp = zero
               do i = 1, m
                  temp = temp + a(i,j)*x(i)
               end do
               y( jy ) = y( jy ) + alpha*temp
               jy      = jy      + incy
           end do
         else
            do j = 1, n
               temp = zero
               ix   = kx
               do i = 1,m
                  temp = temp + a(i,j)*x(ix)
                  ix   = ix   + incx
               end do
               y( jy ) = y( jy ) + alpha*temp
               jy      = jy      + incy
          end do
        end if
      end if
 
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
      function tlsame(ca,cb)

c*********************************************************************72
c
cc TLSAME is a source code copy of BLAS LSAME, testing character equality.
c
c  Parameters:
c
c  cb is assumed to be an upper case letter. tlsame returns .true. if
c  CA is either the same as cb or the equivalent lower case letter.
c
c  CA     - character*1
c  cb     - character*1
c           On entry, CA and cb specify characters to be compared.
c           Unchanged on exit.
c
      implicit none

      integer ioff
      parameter (ioff = 32)

      character*1 ca
      character*1 cb
      logical tlsame
c
c  Test if the characters are equal
c
      tlsame = ca.eq.cb
c
c  Now test for equivalence
c
      if ( .not. tlsame ) then
        tlsame = ichar(ca)-ioff.eq.ichar(cb)
      end if
 
      return
      end
      subroutine uijk ( a, b, c, lda, n, nrep, ttime )

c*********************************************************************72
c
cc UIJK multiplies A = B*C using index order IJK and I unrolling to depth 4.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Workspace, real A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
c    used in the multiplication.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
c    Input, integer NREP, the number of times the multiplication should
c    be carried out.
c
c    Output, real TTIME, an estimate of the CPU time in seconds required
c    for the matrix multiplications.
c
      implicit none

      integer nroll
      parameter (nroll = 4)

      integer lda
      integer n

      real a(lda,n)
      real b(lda,n)
      real c(lda,n)
      integer i
      integer ihi
      integer irep
      integer j
      integer k
      integer nrep
      real time1
      real time2
      real ttime

      ttime = 0.0E+00
 
      do irep = 1, nrep
 
        call r4_set ( a, b, c, lda, n )
 
        call matmul_cpu_timer ( time1 )
 
        ihi = ( n / nroll ) * nroll
        do i = 1, ihi, nroll
          do j = 1, n
            do k = 1, n
              a(i,k)   = a(i,k)   + b(i,j)   * c(j,k)
              a(i+1,k) = a(i+1,k) + b(i+1,j) * c(j,k)
              a(i+2,k) = a(i+2,k) + b(i+2,j) * c(j,k)
              a(i+3,k) = a(i+3,k) + b(i+3,j) * c(j,k)
            end do
          end do
        end do
c
c  Take care of the few cases we missed if N is not a multiple of 4.
c
        do i = ihi+1, n
          do j = 1, n
            do k = 1, n
              a(i,k) = a(i,k) + b(i,j) * c(j,k)
            end do
          end do
        end do
 
        call matmul_cpu_timer ( time2 )
        ttime = ttime + time2 - time1

      end do
 
      return
      end
      subroutine umxma ( a, b, c, lda, n, nrep, ttime )

c*********************************************************************72
c
cc UMXMA multiplies A = B*C using optimized MXMA.
c
c  Discussion:
c
c    Since the routine MXMA is only available on the Cray, in the
c    SCILIB library, the statement
c
c      call mxma(...)
c
c    should be commented out in versions of MATMUL that are to run
c    on other machines.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Workspace, real A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
c    used in the multiplication.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
c    Input, integer NREP, the number of times the multiplication should
c    be carried out.
c
c    Output, real TTIME, an estimate of the CPU time in seconds required
c    for the matrix multiplications.
c
      implicit none

      integer lda
      integer n

      real a(lda,n)
      real b(lda,n)
      real c(lda,n)
      integer inca
      integer incb
      integer incc
      integer irep
      integer nrep
      real time1
      real time2
      real ttime

      inca = 1
      incb = 1
      incc = 1

      ttime = 0.0E+00
 
      do irep = 1, nrep
 
        call r4_set ( a, b, c, lda, n )
 
        call matmul_cpu_timer ( time1 )

c       call mxma ( b, incb, lda, c, incc, lda, a, inca, lda, n, n, n )

        call matmul_cpu_timer ( time2 )
        ttime = ttime + time2 - time1

      end do
 
      return
      end
      subroutine usaxpyc ( a, b, c, lda, n, nrep, ttime )

c*********************************************************************72
c
cc USAXPYC multiplies A = B*C columnwise, using optimized SAXPY.
c
c  Discussion:
c
c    SAXPY is used to carry out the multiplication "columnwise".
c
c    This is equivalent to the following "JKI" code:
c
c       do j = 1, n
c         do k = 1, n
c           do i = 1, n
c              a(i,k) = a(i,k)+b(i,j)*c(j,k)
c            end do
c          end do
c        end do
c
c    Except on the Cray and SGI/IRIS, the statement "call saxpy" below
c    should be commented out.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Workspace, real A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
c    used in the multiplication.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
c    Input, integer NREP, the number of times the multiplication should
c    be carried out.
c
c    Output, real TTIME, an estimate of the CPU time in seconds required
c    for the matrix multiplications.
c
      implicit none

      integer lda
      integer n

      real a(lda,n)
      real b(lda,n)
      real c(lda,n)
      integer irep
      integer j
      integer k
      integer nrep
      real time1
      real time2
      real ttime

      ttime = 0.0E+00
 
      do irep = 1, nrep
 
        call r4_set ( a, b, c, lda, n )
 
        call matmul_cpu_timer ( time1 )
 
        do j = 1, n
          do k = 1, n
            call saxpy ( n, c(j,k), b(1,j), 1, a(1,k), 1 )
          end do
        end do
 
        call matmul_cpu_timer ( time2 )
        ttime = ttime + time2 - time1

      end do
 
      return
      end
      subroutine usaxpyr ( a, b, c, lda, n, nrep, ttime )

c*********************************************************************72
c
cc USAXPYR multiplies A = B*C "rowwise", using optimized SAXPY.
c
c  Discussion:
c
c    This is equivalent to the following "IJK" code:
c
c     do i = 1, n
c       do j = 1, n
c          do k = 1, n
c            a(i,k) = a(i,k)+b(i,j)*c(j,k)
c          end do
c        end do
c      end do
c
c    Except on the Cray and SGI/IRIS, the statement "call saxpy" below
c    should be commented out.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Workspace, real A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
c    used in the multiplication.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
c    Input, integer NREP, the number of times the multiplication should
c    be carried out.
c
c    Output, real TTIME, an estimate of the CPU time in seconds required
c    for the matrix multiplications.
c
      implicit none

      integer lda
      integer n

      real a(lda,n)
      real b(lda,n)
      real c(lda,n)
      integer i
      integer irep
      integer j
      integer nrep
      real time1
      real time2
      real ttime

      ttime = 0.0E+00
 
      do irep = 1, nrep
 
        call r4_set ( a, b, c, lda, n )
 
        call matmul_cpu_timer ( time1 )
 
        do i = 1, n
          do j = 1, n
            call saxpy ( n, b(i,j), c(j,1), lda, a(i,1), lda )
          end do
        end do
 
        call matmul_cpu_timer ( time2 )
        ttime = ttime + time2 - time1

      end do
 
      return
      end
      subroutine usdot ( a, b, c, lda, n, nrep, ttime )

c*********************************************************************72
c
cc USDOT multiplies A = B*C using optimized SDOT.
c
c  Discussion:
c
c    This is equivalent to the following "IKJ" code:
c
c       do i = 1, n
c         do k = 1, n
c           do j = 1, n
c             a(i,k) = a(i,k) + b(i,j) * c(j,k)
c          end do
c        end do
c      end do
c
c    Except on the Cray or SGI/IRIS, the call to SDOT should be commented out.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Workspace, real A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
c    used in the multiplication.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
c    Input, integer NREP, the number of times the multiplication should
c    be carried out.
c
c    Output, real TTIME, an estimate of the CPU time in seconds required
c    for the matrix multiplications.
c
      implicit none

      integer lda
      integer n

      real a(lda,n)
      real b(lda,n)
      real c(lda,n)
      integer i
      integer irep
      integer k
      integer nrep
      real sdot
      real time1
      real time2
      real ttime

      ttime = 0.0E+00
 
      do irep = 1, nrep
 
        call r4_set ( a, b, c, lda, n )
 
        call matmul_cpu_timer ( time1 )
 
        do i = 1, n
          do k = 1, n
            a(i,k) = sdot ( n, b(i,1), lda, c(1,k), 1 )
          end do
        end do
 
        call matmul_cpu_timer ( time2 )
        ttime = ttime + time2 - time1

      end do
 
      return
      end
      subroutine usgemm ( a, b, c, lda, n, nrep, ttime )

c*********************************************************************72
c
cc USGEMM multiplies A = B*C using optimized SGEMM.
c
c  Discussion:
c
c    Except on the Cray or SGI/IRIS, the call to SGEMM should be commented out.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Workspace, real A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
c    used in the multiplication.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
c    Input, integer NREP, the number of times the multiplication should
c    be carried out.
c
c    Output, real TTIME, an estimate of the CPU time in seconds required
c    for the matrix multiplications.
c
      implicit none

      integer lda
      integer n

      real a(lda,n)
      real alpha
      real b(lda,n)
      real beta
      real c(lda,n)
      integer irep
      integer nrep
      real time1
      real time2
      real ttime

      ttime = 0.0E+00
 
      do irep = 1, nrep
 
        call r4_set ( a, b, c, lda, n )
 
        alpha = 1.0E+00
        beta = 0.0E+00
 
        call matmul_cpu_timer ( time1 )

        call sgemm ( 'n', 'n', n, n, n, alpha, b, lda, c, lda, 
     &    beta, a, lda )

        call matmul_cpu_timer ( time2 )
        ttime = ttime + time2 - time1

      end do
 
      return
      end
      subroutine usgemms ( a, b, c, lda, n, nrep, nwork, ttime )

c*********************************************************************72
c
cc USGEMMS multiplies A = B*C using optimized SGEMMS.
c
c  Discussion:
c
c    SGEMMS is the Cray SCILIB variant of the BLAS3 routine SGEMM.  
c    The difference is that SGEMMS uses Strassen's algorithm.
c
c    Except on the Cray, the call to SGEMMS should be commented out,
c    as well as the statement "REAL WORK(NWORK)".
c
c    Notice that the dimensioning of WORK is illegal in standard FORTRAN.
c    A Cray extension allows this creation of scratch arrays.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Workspace, real A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
c    used in the multiplication.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
c    Input, integer NREP, the number of times the multiplication should
c    be carried out.
c
c    Input, integer NWORK, the size to be used for the locally
c    declared workspace array WORK.  Apparently, a recommended value
c    is 3*N*N, rather a lot, really.
c
c    Output, real TTIME, an estimate of the CPU time in seconds required
c    for the matrix multiplications.
c
      implicit none

      integer lda
      integer n

      real a(lda,n)
      real alpha
      real b(lda,n)
      real beta
      real c(lda,n)
      integer irep
      integer nrep
      integer nwork
      real time1
      real time2
      real ttime
c     real work(nwork)

      ttime = 0.0E+00
 
      do irep = 1, nrep
 
        call r4_set ( a, b, c, lda, n )
 
        alpha = 1.0E+00
        beta = 0.0E+00
 
        call matmul_cpu_timer ( time1 )

c       call sgemms ( 'n', 'n', n, n, n, alpha, b, lda, c, lda, beta, a, 
c    &    lda, work )

        call matmul_cpu_timer ( time2 )
        ttime = ttime + time2 - time1

      end do
 
      return
      end
      subroutine utaxpyc ( a, b, c, lda, n, nrep, ttime )

c*********************************************************************72
c
cc UTAXPYC multiplies A = B*C columnwise, using unoptimized SAXPY.
c
c  Discussion:
c
c    This is equivalent to the following "JKI" code:
c
c       do j = 1, n
c         do k = 1, n
c           do i = 1, n
c             a(i,k) = a(i,k)+b(i,j)*c(j,k)
c           end do
c         end do
c       end do
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Workspace, real A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
c    used in the multiplication.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
c    Input, integer NREP, the number of times the multiplication should
c    be carried out.
c
c    Output, real TTIME, an estimate of the CPU time in seconds required
c    for the matrix multiplications.
c
      implicit none

      integer lda
      integer n

      real a(lda,n)
      real b(lda,n)
      real c(lda,n)
      integer irep
      integer j
      integer k
      integer nrep
      real time1
      real time2
      real ttime

      ttime = 0.0E+00
 
      do irep = 1, nrep
 
        call r4_set ( a, b, c, lda, n )
 
        call matmul_cpu_timer ( time1 )
 
        do j = 1, n
          do k = 1, n
            call taxpy ( n, c(j,k), b(1,j), 1, a(1,k), 1 )
          end do
        end do
 
        call matmul_cpu_timer ( time2 )
        ttime = ttime + time2 - time1

      end do
 
      return
      end
      subroutine utaxpyr ( a, b, c, lda, n, nrep, ttime )

c*********************************************************************72
c
cc UTAXPYR multiplies A = B*C rowwise using source code SAXPY.
c
c  Discussion:
c
c    This is equivalent to the following "IJK" code:
c
c        do i = 1, n
c          do j = 1, n
c            do k = 1, n
c              a(i,k) = a(i,k)+b(i,j)*c(j,k)
c            end do
c          end do
c        end do
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Workspace, real A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
c    used in the multiplication.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
c    Input, integer NREP, the number of times the multiplication should
c    be carried out.
c
c    Output, real TTIME, an estimate of the CPU time in seconds required
c    for the matrix multiplications.
c
      implicit none

      integer lda
      integer n

      real a(lda,n)
      real b(lda,n)
      real c(lda,n)
      integer i
      integer irep
      integer j
      integer nrep
      real time1
      real time2
      real ttime

      ttime = 0.0E+00
 
      do irep = 1, nrep
 
        call r4_set ( a, b, c, lda, n )
 
        call matmul_cpu_timer ( time1 )
 
        do i = 1, n
          do j = 1, n
            call taxpy ( n, b(i,j), c(j,1), lda, a(i,1), lda )
          end do
        end do
 
        call matmul_cpu_timer ( time2 )
        ttime = ttime + time2 - time1

      end do
 
      return
      end
      subroutine utdot ( a, b, c, lda, n, nrep, ttime )

c*********************************************************************72
c
cc UTDOT multiplies A = B * C using source code SDOT.
c
c  Discussion:
c
c    This is equivalent to the following "IKJ" code:
c
c      do i = 1, n
c        do k = 1, n
c          do j = 1, n
c            a(i,k) = a(i,k) + b(i,j) * c(j,k)
c          end do
c        end do
c      end do
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Workspace, real A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
c    used in the multiplication.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
c    Input, integer NREP, the number of times the multiplication should
c    be carried out.
c
c    Output, real TTIME, an estimate of the CPU time in seconds required
c    for the matrix multiplications.
c
      implicit none

      integer lda
      integer n

      real a(lda,n)
      real b(lda,n)
      real c(lda,n)
      integer i
      integer irep
      integer k
      integer nrep
      real tdot
      real time1
      real time2
      real ttime

      ttime = 0.0E+00
 
      do irep = 1, nrep
 
        call r4_set ( a, b, c, lda, n )
 
        call matmul_cpu_timer ( time1 )
 
        do i = 1, n
          do k = 1, n
            a(i,k) = tdot ( n, b(i,1), lda, c(1,k), 1 )
          end do
        end do
 
        call matmul_cpu_timer ( time2 )
        ttime = ttime + time2 - time1
 
      end do
 
      return
      end
      subroutine utgemm ( a, b, c, lda, n, nrep, ttime )

c*********************************************************************72
c
cc UTGEMM multiplies A = B*C using SGEMM.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    28 August 1999
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Workspace, real A(LDA,N), B(LDA,N), C(LDA,N), the three matrices
c    used in the multiplication.
c
c    Input, integer LDA, the leading dimension used for arrays.
c
c    Input, integer N, the number of rows and columns in the matrices.
c
c    Input, integer NREP, the number of times the multiplication should
c    be carried out.
c
c    Output, real TTIME, an estimate of the CPU time in seconds required
c    for the matrix multiplications.
c
      implicit none

      integer lda
      integer n

      real a(lda,n)
      real alpha
      real b(lda,n)
      real beta
      real c(lda,n)
      integer irep
      integer nrep
      real time1
      real time2
      real ttime

      ttime = 0.0E+00
 
      do irep = 1, nrep
 
        call r4_set ( a, b, c, lda, n )
 
        alpha = 1.0E+00
        beta = 0.0E+00
 
        call matmul_cpu_timer ( time1 )

        call tgemm ( 'n', 'n', n, n, n, alpha, b, lda, c, lda, beta,
     &    a, lda )

        call matmul_cpu_timer ( time2 )
        ttime = ttime + time2 - time1
 
      end do
 
      return
      end
