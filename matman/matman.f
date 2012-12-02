c  matman.f  Version 1.59  06 July 1998
c
      program matman

c*********************************************************************72
c
cc MATMAN is a program for interactive linear algebra demonstrations.
c
c
c  Version 1.59   06 July 1998
c
c    * Made the linear algebra sample problems random.
c    * Forced A(I,J) = A(J,I) = 0 exactly in Jacobi method.
c    * Replaced "NEXPER" by newer version from SUBSET named "PERNEX".
c    * Used new format for argument documentation.
c    * Inserted latest versions of CHRPAK/SUBPAK/SUBSET routines.
c    * Passed MAXINT to DEC and RAT routines, allowed user to set it.
c    * Used a better RANDOM routine, and updated MATKEY.DAT.
c
c  Version 1.58  10 April 1996
c
c    CGC requested the following changes:
c
c    1) The confirmation of the A command should be
c
c       Row 5 <= 3 Row 2 + Row 5, with "Row 5" coming last.  (FIXED)
c
c    2) In LP mode, the feasibility ratios are printed out in
c       both real and rational values, but they disagree.  (FIXED)
c
c    I made the following changes:
c
c    3) I also tried to add some more comments at the beginnings of
c       routines, to define the variables.
c
c    4) I made CHRREL print out in G14.7 format, to try to get seven
c       digit accuracy where possible, and modified RELPRN to print out
c       7 decimals as well.
c
c    5) I also altered SETDIG to allow the user to exceed the recommended
c       maximum of MAXDIG digits, with a warning.
c
c    6) I modified the linear algebra optimization checker to print out
c       the column, as well as the row, where rule 1 is violated.
c
c    7) I changed RELPRN so that if the printed quantities are all integers,
c       the printout is more compact.
c
c    8) I modified AUTERO to eliminate unnecessary row operations,
c       where an entry is already zero.
c
c    9) Added an extra check in ROWADD to skip out immediately if the
c       multiplier is zero.
c
c    10) Added documentation for the DECIMAL, RATIONAL, and real
c       commands, which will eventually replace the "F" command.
c
c    11) Updated my address and EMAIL address.
c
c
c  Version 1.57  15 December 1995
c
c    CGC requested the following changes:
c
c    1) In LP mode, when doing a two phase problem, if you convert
c       from one arithmetic form to another, the objective function
c       data in row NROW+1 is not converted.  So I modified FORM
c       to convert all the rows and columns, not just NROW by NCOL.
c       FIXED.
c
c    2) Move the L command to the short menu, and add to the short
c       menu a description of how to get the long menu.  DONE.
c
c    3) Replace the prompt which follows the "Enter command" prompt
c       with ("H" for short menu, "HELP" for full menu, ? for full help).
c       DONE.
c
c    4) Stumbled across a slight problem.  When in LP mode, using
c       fractional arithmetic, and you enter a problem with artificial
c       variables, the last column of the auxilliary objective function
c       was set to 0/0 instead of 0/1.  FIXED.
c
c    5) Instead of printing the ERO determinant after every operation,
c       I added an EDET command to print it out only on demand.
c       This is a request CGC made earlier, and which we had both
c       forgotten.
c
c    6) Cosmetic change: rational matrices were printed out with two
c       trailing blank lines, which I cut back to one.
c
c    7) Fixed an obscure error in RELREA, which only caused problems
c       on the ALPHA.  CHRCTR was reading numbers all the way to the
c       last blank, setting LCHAR=NCHAR=80, and then RELREA was asking
c       if character LCHAR+1 was a '/'.
c
c    8) I replaced all the "WRITE(*,*)" statements by the more robust
c       "OUTPUT=...", "CALL CHRWRT()" pair.  Some error messages were
c       only going to the screen, and not the permanent output file.
c
c    9) Discovered that when a sample problem is chosen, only NROW
c       by NCOL of the matrix area was set, leaving garbage possibly
c       in other areas.  I rectified this, zeroing out all the
c       rest of the matrix area, via a routine INIMAT.
c
c
c  Version 1.56  12 October 1995
c
c    I fixed the program, so that you can type
c
c      Row 2 <=> Row 3
c
c    (spelling out the word "Row") if you want to.
c
c    I found a logic mistake in DECMUL which occurs if one of the
c    input quantities is the same as the output, and I fixed it.
c
c    I modified the program so that, if you are working with 4
c    digit decimals, any decimal input is automatically truncated
c    to 4 decimals as well.
c
c    Replaced CHLDEC by a routine which is exact.
c
c    Corrected DECRAT, and many other decimal discrepancies.
c
c    Corrected DECADD, so that decimal addition is exact.
c
c    Added the DECIMAL, RATIONAL and real commands, although I
c    did not mention them.
c
c    Corrected the phrase "row reduced echelon form" by replacing
c    it with "reduced row echelon form".
c
c    Added the BASIC command to allow the user to assign a row
c    to a basic variable without using the Change command.
c
c    An unneeded change to CHRINP disabled the "<" command, but
c    I fixed that.
c
c    I modified routine PASS so that its default key corresponds
c    to the value currently stored in MATKEY.DAT.
c
c    Corrected TRANSC so that VMS output files would have correct
c    carriagecontrol.
c
c
c  Version 1.55  06 October 1995
c
c    Experimenting with allowing longer command names, so that I
c    can have more reasonable names.  Changed COMNEW to 4 characters
c    in length.  Now I can request a determinant with DET, and
c    a transpose with TR.  Note that TR is potentially in conflict
c    with the "Type a Row" command, unless the user puts a space there.
c
c    Added the determinant command.
c
c    Dropped the T C, T R and T E commands.
c
c    Added a square matrix example for the determinant problem.
c
c    Print out the determinant of the ERO's.
c
c
c  Notes:
c
c    The "<" has a logical flaw.  If you use the "X" command in
c    the input file, you will return to the user input, not the
c    file input, once you're done.  You need to save each input
c    unit number and file name in a stack in order to properly
c    recover.
c
c    If you add NCON, you're going to have to save it in RESTORE
c    and in read/write examples as well...
c
c    CGC requested the ability to add a row or column for linear
c    programming!  In this case, it would correspond to a slack
c    variable.
c
c
c  Version 1.54  05 September 1994
c
c    EVJACO allows the user to type an integer or a character.
c    But INTREA was allowing an error message to appear if the
c    user typed a "Q", because the IHUSH parameter was reset
c    to 0.  I took out the resetting.
c
c    Because of capability of entering several commands, separated
c    by a semicolon, comments weren't being ignored properly.
c    So now, once a comment is recognized, it's blanked out.
c
c
c  Version 1.53  25 July 1994
c
c    Replaced constant "0" by variable "ITERM" in all calls to
c    CHRREA.
c
c    Added "Y" command, so user can turn autoprinting off or on.
c
c    Autoprint after "V" command.
c
c    Added capability to enter several commands, separated by a semicolon.
c
c
c  Version 1.52  12 May 1994
c
c
c    12 May 1994:
c
c      Corrected SETDIG so that the value the user typed in does
c      not immediately overwrite the current value, until it has
c      been checked.
c
c    11 May 1994:
c
c      Struggling with a "final" problem, in which DECPRN does
c      not work when trying to print out the linear programming
c      solution for the advanced sample problem.  Think it's
c      fixed now.
c
c    10 May 1994:
c
c      Program freezes when going into decimal fraction mode,
c      with simple linear programming problem.
c
c      It would be preferable to delete trailing zeroes from the
c      printouts of real numbers by RELPRN.
c
c      Tracking down a bug in DECREA, I think.
c
c      The linear algebra stuff seems to be working properly with
c      decimal arithmetic.  I still need to check linear
c      programming.
c
c      Added an option to use a default key.
c
c      Inserted new version of CHLDEC which does not try to convert
c      IVAL*10**JVAL into a real first, and so should be able to print
c      out the exact representation, with no trailing blanks or
c      roundoff problems.  Looks a lot better.
c
c      I corrected the sample problem routines, to take account of
c      the 3 different arithmetic modes.
c
c      How about an option to generate a random test problem,
c      similar to the sample?
c
c    09 May 1994:
c
c      Fixed DECPRN.
c      Fixed RELDEC and DECREL.
c      Now, CHLDEC is printing 0.4 as 0.39999999999, and
c      DECREA is reading too many digits, and overflowing.
c
c    08 May 1994:
c
c      Wrote DBLDEC.  Updated DECMUL.
c      Rewrote FORM.  Created new RELDEC, DECREL, RATDEC,
c      DECRAT.
c
c      I checked all the lines where IFORM.EQ.2 occurred,
c      and tried to correct them, but gave up in the LP routines.
c      Check them later!
c
c      I still need to write CHLDEC and fix up DECPRN.
c
c    07 May 1994:
c
c      Rewrote DECDIV, adding argument NDIG, and requiring
c      new routine DBLDEC, converts double precision quantity
c      to decimal.
c
c      Added "N" command allowing user to set NDIG.
c
c    06 May 1994:
c
c      Proposal: the decimals should be stored
c      as SMANT (an integer representing the signed mantissa) and
c      IEXP (an integer representing the power of 10).
c
c      All calculations should be carried out by converting to a
c      real value first.  So I need to write just two routines,
c      RELDEC and DECREL.  Oh, and I forgot to mention that the
c      restriction on NDIG simply restricts the size of SMANT.
c
c    05 May 1994:
c
c      Modified AUTERO to use SCADIV rather than
c      SCAMUL.  However, I still have something to complain about.
c      The interim values in the matrix (using decimals)
c      are not themselves decimals.  Once that happens,
c      the whole point is lost.  Is this DECDIV's fault?
c
c    03 May 1994:
c
c      What's with this DECPRN routine?
c
c    03 May 1994:
c
c      Fixed error found on 23 April.  It was a tiny mistake in
c      RELDEC.
c
c    23 April 1994:
c
c      Error:  I specified "B" for sample problem,
c      specified "F" and converted to "Decimal", chose "1" digit.
c      4 by 4 matrix got divided by 10, while RHS was correctly
c      rounded.
c
c      Fixed a mistake in FORM which meant that conversions were
c      not being done when going from real to rational.
c
c      Replaced all labeled DO loops by DO/ENDDO pairs.
c
c      Added a "D" command that does division, renaming old
c      D (disk file) command to "K".
c
c
c  Version 1.51  19 September 1993
c
c    Changed from WRITE(* to WRITE(6 so that output redirection
c    works on DOS machines.
c
c    Also changed READ(* to READ(5.
c
c    I removed all occurrences of real(...), to make it easier
c    to convert the code to double precision if desired.
c    To convert this code to double precision, the only change
c    needed is a global substitute of "double precision" for
c    "real".
c
c  Version 1.50  06 June 1993
c
c    Added a sample linear system solve problem.
c
c    Dropped the "N" command.
c
c    Reversed order of arguments in call to CHRWRT.
c
c  Version 1.49  04 May 1993
c
c    Well guess what, the new version of Language Systems FORTRAN
c    seems to have fixed that bug!
c
c    Renamed SCALE to SCALER, because of a possible conflict
c    with a Macintosh routine (like this is going to help!).
c
c    A continuing bug that occurs on the Macintosh has led me
c    to add all sorts of checks for "NULL" characters in strings.
c    Right now it's just a hunch.
c
c    The "Z" command would claim an error had occurred if you
c    had not yet set up a matrix.
c
c    Changed CHRWRT to print an explicit blank as carriage control
c    on the Mac, since Language Systems "FORTRAN" will otherwise
c    print a null.
c
c    To make life easier for the "<" command, I dropped the
c    initial demand for arithmetic specification.
c
c    Added "<" command to allow user to specify an input file.
c    This is because it's hard to do on VMS via system commands,
c    and impossible on a Macintosh.
c
c    Added autoprint after pivoting.
c
c    Program should no longer fail if using rational arithmetic
c    and an overflow or underflow occurs.  Right now, MATMAN will
c    catch this problem, and halt the computation.  A better
c    solution would allow the user to request that overflows and
c    underflows be "rounded" to decimal and recomputed as ratios.
c
c    Noticed that Macintosh requires FORTRAN carriage control (1X)
c    for output to console, so had to modify CHRWRT.
c
c    Modified advanced LP problem, and corrected label.
c
c    Added FLUSHL, and forced CHRREA to flush the string left
c    once it has been read.  This is so that, if I like, I can
c    type "B S" and have it mean the same as "BS".
c
c    Renamed SETSOL to LPSOL.
c
c    Modified INTREA to have IHUSH parameter, so that I can
c    type a "Q" to quit in the Jacobi iteration.
c
c    Added the # feature, which allows a one line comment
c    beginning with the sharp symbol.
c
c    Added the $ and % commands, which allows me to turn paging off
c    and back on.
c
c    A row or column can be added to the matrix, or deleted
c    from it, in linear algebra mode.  An added row or column
c    can be inserted anywhere in the matrix.
c
c    CGC requested automatic printout of a matrix or tableau
c    when it is entered, and I have added that.
c
c    In linear algebra mode, the "o" command will now check
c    whether the matrix is in row echelon or reduced row echelon
c    form.
c
c    I added LEQI to this program, to avoid having to capitalize
c    everything.
c
c    Set up new paging routine, which, if it works, will be
c    added to MATALG also.  This also allows me to drop that
c    stupid ICOMP parameter from CHRWRT.
c
c    Cleaned up the main program a bit, so that the commands
c    are all part of one big IF/ELSIF block.
c
c    Ran the CLEANER program, to standardize indentation and
c    statement numbering.
c
c  Version 1.48  16 April 1993
c
c    Ran the STRIPPER program, to lower case all statements.
c
c    Removed last argument of CHRREA.
c
c    Replaced old version of CHRCTR by a newer, better one which
c    does not suffer from integer "wrap around" when a large number
c    of decimal places are entered.
c
c  Version 1.47  25 April 1992
c
c    Cleaning up program format:
c      Continuation character is now always "&".
c      Marked the beginning of each routine with a "C****..." line.
c      Routines placed in alphabetical order.
c      Declared all variables.
c    Set SOL, ISLTOP and ISLBOT to be vectors, rather than 2D
c      arrays with a first dimension of 1.
c    Made JACOBI automatically iterative.
c    JACOBI will now accept nonsymmetric matrices.
c    CHRCTI and CHRCTR reset LCHAR = 0 on error.
c
c  Minor modification, 25 September 1991
c
c    Dropped "MAXCOL" from the argument list of JACOBI.
c
c  Version 1.46  01 February 1991
c
c    Tried a modification in LOOP 20 in HLPVMS.
c
c  Version 1.45  29 January 1991
c
c    Corrected error in example file reading.
c
c  Version 1.44  24 December 1990
c
c    Moved initializations to INIT.
c
c  Version 1.43  07 December 1990
c
c    ROWADD will now refuse to add a multiple of a row to itself.
c
c  Version 1.42  03 December 1990
c
c    Added Jacobi example.
c
c    CGC complained that "1" in the P column was missing, for
c      problems with artificial variables.  I tried to restore
c      that.
c
c    CGC complained that in problems with artificial variables, the
c      artificial objective function picks up the constant term of
c      the original objective.
c
c    Also added sample problem with artificial variables.
c
c    Added forced typeout of matrix after each Jacobi step.
c
c  Version 1.41  06 November 1990
c
c    FR and FF commands force real and rational arithmetic.
c
c    Added JACOBI routine to do Jacobi rotations.
c
c    Moved "hello" stuff to routine HELLO.
c
c    Moved transcript stuff to a routine TRANSC.
c
c  Version 1.40  15 October 1990
c
c    Modification made to allow matrix to be entered in entirety,
c    rather than a row at a time.
c
c  Version 1.39  12 October 1990
c
c    CHRREA modified so that output need not be capitalized.  This
c    keeps filenames from being forcibly capitalized.
c
c  Version 1.38
c
c    * command allows you to transpose the matrix, LA mode
c      only.
c
c  Version 1.37
c
c    In non linear programming mode, you may enter the entire
c      matrix, or several rows at a time, on one line, if you like.
c
c  Version 1.36
c
c    Updated interface to CHRPAK.
c
c    Inserted obsolete CHRPAK routines into MATMAN source.
c
c  Version 1.35
c
c    Correctly writes out hidden objective row for problems using
c    artificial variables.
c
c  Version 1.34
c
c    Corrected a transposition of variables that meant that, for
c    linear programming problems, constraints and variables were
c    interchanged.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer maxcol
      parameter (maxcol=30)

      integer maxrow
      parameter (maxrow = 16)

      real a(maxrow,maxcol)
      logical autop
      real b(maxrow,maxcol)
      real c(maxrow,maxcol)
      character*1 chineq(maxrow)
      character*22 chldec
      character*22 chlrat
      character*14 chrrel
      character*22 chrtmp
c
c  Experimentally setting COMNEW to more than just 1 letter.
c
      character*20 comnew
      character*20 comold
      real det
      real dete
      character*60 filex
      character*60 filhlp
      character*60 filinp
      character*60 filkey
      character*60 filtrn
      integer iabot(maxrow,maxcol)
      integer iarray(maxrow)
      integer iatop(maxrow,maxcol)
      integer iauthr
      integer iauto
      integer ibase(maxrow)
      integer ibaseb(maxrow)
      integer ibasec(maxrow)
      integer ibbot(maxrow,maxcol)
      integer ibtop(maxrow,maxcol)
      integer icbot(maxrow,maxcol)
      integer icol
      integer ictop(maxrow,maxcol)
      integer idbot
      integer idebot
      integer idtop
      integer idetop
      integer ierror
      integer iform
      integer imat
      integer iopti
      integer iounit(4)
      integer iprint
      integer irow
      integer irow1
      integer irow2
      character*1 isay
      integer isbot
      integer iseed
      integer islbot(maxcol)
      integer isltop(maxcol)
      integer istop
      integer iterm
      integer jform
      logical ldigit
      integer lenchr
      logical leqi
      character*80 line
      character*80 line2
      integer lpage
      integer lpmoda
      integer lpmodb
      integer lpmodc
      integer maxdig
      integer maxint
      integer nart
      integer nartb
      integer nartc
      integer ncol
      integer ncolb
      integer ncolc
      integer ncon
      integer ndig
      integer nline
      integer nrow
      integer nrowb
      integer nrowc
      integer nslak
      integer nslakb
      integer nslakc
      integer nvar
      integer nvarb
      integer nvarc
      character*100 output
      character*80 prompt
      real sol(maxcol)
      real sval
c
c  Initializations.
c
      call init ( a, autop, chineq, comnew, comold, dete, filex,
     &  filhlp, filinp, filkey, filtrn, iabot, iatop, iauthr, ibase,
     &  idebot, idetop, ierror, iform, imat, iounit, iprint, iseed,
     &  islbot,
     &  isltop, line, lpage, lpmoda, maxcol, maxdig, maxint, maxrow,
     &  nart, ncol, ncon, ndig, nline, nrow, nslak, nvar, sol )
 
      call copmat ( a, b, iatop, iabot, ibtop, ibbot, ibase, ibaseb, 
     &  lpmoda, lpmodb, maxcol, maxrow, nart, nartb, ncol, ncolb, 
     &  nrow, nrowb, nslak, nslakb, nvar, nvarb )
 
      call copmat ( a, c, iatop, iabot, ictop, icbot, ibase, ibasec, 
     &  lpmoda, lpmodc, maxcol, maxrow, nart, nartc, ncol, ncolc, 
     &  nrow, nrowc, nslak, nslakc, nvar, nvarc )
c
c  Say hello.
c
      call hello ( iounit, output )
c
c  Print out arithmetic warning.
c
      if ( iform .eq. 0 ) then
        call ratwrn ( iounit, maxint, output )
      else if ( iform .eq. 1 ) then
        call relwrn ( iounit, output )
      else if ( iform .eq. 2 ) then
        call decwrn ( iounit, output )
      end if
c
c  Get the next command from the user.
c
10    continue
 
      line = ' '
      nline = 0
 
      if ( ierror .ne. 0 ) then

        output = ' '
        call chrwrt ( iounit, output )
        output = 'Because of an error, your command was not completed.'
        call chrwrt ( iounit, output )
        output = 'We return to the main menu.'
        call chrwrt ( iounit, output )
        ierror = 0
c
c  Wipe out the offending command line.
c
        nline = 0
 
        if ( iounit(1) .ne. 0 ) then
          close ( unit = iounit(1) )
          iounit(1) = 0
          output = ' '
          call chrwrt ( iounit, output )
          output = 'Because an error occurred, we are closing'
          call chrwrt ( iounit, output )
          output = 'the input file, and requiring you to respond'
          call chrwrt ( iounit, output )
          output = 'directly!'
          call chrwrt ( iounit, output )
        end if
 
      end if
c
c  Insert a blank line.
c
      if ( comnew .ne. '#' ) then
        output = ' '
        call chrwrt ( iounit, output )
      end if
c
c  Save the name of the previous command as COMOLD, in case we need
c  to undo it.  But only save "interesting" commands.
c
      if ( .not.
     &  (leqi ( comnew, 'DET' ) .or. 
     &   leqi ( comnew, 'H' ) .or. 
     &   leqi ( comnew, 'HELP' ) .or. 
     &   leqi ( comnew, 'N' ) .or. 
     &   leqi ( comnew, 'O' ) .or. 
     &   leqi ( comnew, 'S' ) .or. 
     &   leqi ( comnew, 'T' ) .or. 
     &   leqi ( comnew, '$' ) .or. 
     &   leqi ( comnew, '?' ) .or. 
     &   leqi ( comnew, '%' ) .or. 
     &   leqi ( comnew, '<' ) .or. 
     &   leqi ( comnew, '#' ) ) ) then
        comold = comnew
      end if
 
      if ( comnew .ne. '#' ) then
        prompt = 'command? ("H" for short menu, ' //
     &    '"HELP" for full menu, ? for full help)'
      else
        prompt = ' '
      end if
c
c  No check for terminators.
c
      iterm = 0
      call chrrea ( line2, line, nline, prompt, iounit, ierror,
     &  iterm )

      nline = lenchr ( line2 )
      line = line2
 
      if ( line2 .eq. ' ' .or. nline .eq. 0 ) go to 10
 
      if ( ierror .ne. 0 ) then
        ierror = 0
        comnew = 'Q'
      end if
c
c  Check to see if the command is an ERO, in a special format.
c
      call chrdb1 ( line2 )
 
      if ( leqi ( line2(1:1), 'R' ) .and. 
     &  ldigit ( line2(2:2) ) ) then

        call chkero ( comnew, ierror, iounit, line2, output )
        if ( ierror .ne. 0 ) go to 10
        line = line2
        nline = lenchr(line)

      else if ( leqi ( line2(1:3), 'ROW' ) .and. 
     &  ldigit ( line2(4:4) ) ) then

        call chkero ( comnew, ierror, iounit, line2, output )
        if ( ierror .ne. 0 ) go to 10
        line = line2
        nline = lenchr(line)

      else

        comnew = ' '

      end if
c
c  If command was not an ERO that had to be translated, read it
c  the regular way.
c
c  Blank, slash, comma, semicolon, equals terminate COMNEW input.
c
      if ( comnew .eq. ' ' ) then

        nline = lenchr(line)
        iterm = 1

        call chrrea ( comnew, line, nline, prompt, iounit, ierror, 
     &    iterm )

        if ( ierror .ne. 0 ) then
          ierror = 0
          comnew = 'Q'
        end if

      end if
c
c  If the "Z" command was issued, the user must give authorization
c  the first time.
c
      if ( leqi ( comnew, 'Z' ) .and. iauthr .eq. 0 ) then

        call pass ( filkey, iauthr, ierror, iounit, line, nline, output,
     &    prompt )

        if ( iauthr .eq. 0 .or. imat .eq. 0 ) then
          go to 10
        end if

      end if
c
c  Jump here when one command needs to switch to another.
c
20    continue
c
c  Save a copy of the matrix A in B before the operation, but only
c  for certain commands.
c
      if ( 
     &  leqi ( comnew, 'A' ) .or. 
     &  leqi ( comnew, 'B' ) .or. 
     &  leqi ( comnew, 'BASIC' ) .or. 
     &  leqi ( comnew, 'C' ) .or. 
     &  leqi ( comnew, 'D' ) .or. 
     &  leqi ( comnew, 'E' ) .or. 
     &  leqi ( comnew, 'F' ) .or. 
     &  leqi ( comnew, 'G' ) .or. 
     &  leqi ( comnew, 'I' ) .or. 
     &  leqi ( comnew, 'J' ) .or. 
     &  leqi ( comnew, 'L' ) .or. 
     &  leqi ( comnew, 'M' ) .or. 
     &  leqi ( comnew, 'P' ) .or. 
     &  leqi ( comnew, 'R' ) .or. 
     &  leqi ( comnew, 'TR' ) .or. 
     &  leqi ( comnew, 'V' ) .or. 
     &  leqi ( comnew, 'X' ) .or. 
     &  leqi ( comnew, 'Z' ) ) then
  
        call copmat ( a, b, iatop, iabot, ibtop, ibbot, ibase, ibaseb,
     &    lpmoda, lpmodb, maxcol, maxrow, nart, nartb, ncol, ncolb, 
     &    nrow, nrowb, nslak, nslakb, nvar, nvarb )
 
      end if
c
c  A=Add a multiple of one row to another.
c
      if ( leqi ( comnew, 'A' ) ) then
 
        call chkadd ( ierror, iform, imat, iounit, irow1, irow2, istop,
     &    isbot, line, maxdig, ndig, nline, nrow, output, prompt, sval )
 
        if ( ierror .eq. 0 ) then
 
          call rowadd ( a, iatop, iabot, ierror, iform, iounit, irow1,
     &      irow2, maxcol, maxint, maxrow, ncol, ndig, output, sval, 
     &      istop, isbot )
 
          iprint = 1
 
        end if
c
c  B=Set up sample problem.
c
      else if ( leqi ( comnew, 'B' ) ) then
 
        call sample ( a, chineq, iatop, iabot, ibase, ierror, iform,
     &    imat, iounit, iseed, line, lpmoda, maxcol, maxrow, nart, 
     &    ncol, nline, nrow, nslak, nvar, output, prompt )
 
        if ( ierror .eq. 0 ) then
 
          if ( iform .eq. 0 ) then
            idetop = 1
            idebot = 1
          else if ( iform .eq. 1 ) then
            dete = 1.0
          else if ( iform .eq. 2 ) then
            idetop = 1
            idebot = 0
          end if
 
          call copmat ( a, c, iatop, iabot, ictop, icbot, ibase,
     &      ibasec, lpmoda, lpmodc, maxcol, maxrow, nart, nartc, ncol,
     &      ncolc, nrow, nrowc, nslak, nslakc, nvar, nvarc )
 
          iprint = 1
 
      end if
c
c  BASIC = Assign row I to basic variable J.
c
      else if ( leqi ( comnew, 'BASIC' ) ) then
 
        call basic ( ibase, ierror, imat, iounit, line, lpmoda,
     &    maxrow, nart, nline, nrow, nslak, nvar, output, prompt )
c
c  C=Change entry.
c
      else if ( leqi ( comnew, 'C' ) ) then
 
        call change ( a, iatop, iabot, ierror, iform, imat, iounit,
     &    line, maxcol, maxdig, maxrow, ncol, ndig, nline, nrow,
     &    output, prompt )
 
        iprint = 1
 
        if ( lpmoda .eq. 0 ) then

          if ( 
     &      ( iform .eq. 0 .and. idetop .ne. idebot ) .or. 
     &      ( iform .eq. 1 .and. dete .ne. 1.0 ) .or. 
     &      ( iform .eq. 2 .and. idetop .ne. idebot ) ) then

            output = 'Warning!  Changing the matrix has probably made'
            call chrwrt ( iounit, output )
            output = 'the ERO determinant incorrect.'
            call chrwrt ( iounit, output )

          end if

        end if
c
c  D = Divide row by scalar.
c
      else if ( leqi ( comnew, 'D' ) ) then
 
        call divide ( a, dete, iatop, iabot, idetop, idebot, ierror,
     &    iform, imat, iounit, line, maxcol, maxdig, maxint, maxrow,
     &    ncol, ndig, nline, nrow, output, prompt )
 
        if ( ierror .eq. 0 ) then
          iprint = 1
        end if
c
c  DECimal = use decimal arithmetic
c
      else if ( leqi ( comnew(1:3), 'DEC' ) ) then
 
        jform = 2
 
        call form ( a, b, c, dete, iatop, iabot, ibtop, ibbot, ictop,
     &    icbot, idetop, idebot, iform, imat, iounit, jform, maxcol,
     &    maxint, maxrow, ndig, output )
 
        iprint = 1
c
c  DET = Determinant of the matrix.
c
      else if ( leqi ( comnew(1:3), 'DET' ) ) then
 
        call chkdet ( ierror, imat, iounit, lpmoda, ncol, nrow, 
     &    output )
 
        if ( ierror .eq. 0 ) then
 
          if ( iform .eq. 0 ) then

            call ratdet ( iatop, iabot, idtop, idbot, iarray, ierror,
     &        maxrow, maxint, nrow )

            chrtmp = chlrat ( idtop, idbot )

          else if ( iform .eq. 1 ) then

            call reldet ( a, det, iarray, maxrow, nrow )

            chrtmp = chrrel(det)

          else if ( iform .eq. 2 ) then

            call decdet ( iarray, iatop, iabot, idtop, idbot, ierror,
     &        maxint, maxrow, nrow, ndig )

              chrtmp = chldec ( idtop, idbot )

          end if

          if ( ierror .eq. 0 ) then
            output = ' '
            call chrwrt ( iounit, output )
            output = 'The determinant is ' // chrtmp
            call chrwrt ( iounit, output )
          end if
 
        end if
c
c  E=Enter problem definition
c
      else if ( leqi ( comnew, 'E' ) ) then
 
        if ( lpmoda .eq. 0 ) then
 
          call lainp0 ( a, iatop, iabot, ierror, iform, iounit, line,
     &      maxcol, maxrow, ncol, nline, nrow, nvar, output, prompt )

          irow = 1
          icol = 1
 
          call lainp1 ( a, iabot, iatop, icol, ierror, iform, iounit,
     &      irow, line, maxcol, maxdig, maxrow, ncol, ndig, nline, nrow,
     &      output, prompt )
 
          if ( iform .eq. 0 ) then
            idetop = 1
            idebot = 1
          else if ( iform .eq. 1 ) then
            dete = 1.0
          else if ( iform .eq. 2 ) then
            idetop = 1
            idebot = 0
          end if
 
        else
 
          call lpinp ( a, chineq, iatop, iabot, ibase, ierror, iform,
     &      iounit, line, maxcol, maxdig, maxrow, nart, ncol, ncon,
     &      ndig, nline, nrow, nslak, nvar, output, prompt )
 
        end if
 
        if ( iounit(1) .eq. 41 ) then
          close(unit = iounit(1))
          iounit(1) = 0
          output = 'The example has been read.'
          call chrwrt ( iounit, output )
        end if
 
        if ( ierror .ne. 0 ) go to 10
 
        imat = 1
 
        call copmat ( a, c, iatop, iabot, ictop, icbot, ibase, ibasec,
     &    lpmoda, lpmodc, maxcol, maxrow, nart, nartc, ncol, ncolc, 
     &    nrow, nrowc, nslak, nslakc, nvar, nvarc )
 
        output = 'A copy of this matrix is being saved.'
        call chrwrt ( iounit, output )
 
        output = 'The "R" command can bring it back.'
        call chrwrt ( iounit, output )
 
        iprint = 1
c
c  EDET = ERO matrix determinant.
c
      else if ( leqi ( comnew, 'EDET' ) ) then
 
        if ( lpmoda .eq. 0 ) then
 
          output = ' '
          call chrwrt ( iounit, output )

          if ( iform .eq. 0 ) then
            chrtmp = chlrat ( idetop, idebot )
          else if ( iform .eq. 1 ) then
            chrtmp = chrrel ( dete )
          else if ( iform .eq. 2 ) then
            chrtmp = chldec ( idetop, idebot )
          end if
 
          output = 'The ERO determinant is ' // chrtmp
          call chrwrt ( iounit, output )

        end if
c
c  F=Form of arithmetic.
c
      else if ( leqi ( comnew, 'F' ) ) then
 
        prompt = '"F" for fractions, "R" for real, "D" for decimal.'
        iterm = 0
        call chrrea ( isay, line, nline, prompt, iounit, ierror, 
     &    iterm )
        if ( ierror .ne. 0 ) go to 10
 
        if ( leqi ( isay, 'F' ) ) then
          jform = 0
        else if ( leqi ( isay, 'R' ) ) then
          jform = 1
        else if ( leqi ( isay, 'D' ) ) then
          jform = 2
        else
          ierror = 1
          output = 'Your choice was not acceptable!'
          call chrwrt ( iounit, output )
          go to 10
        end if
 
        call form ( a, b, c, dete, iatop, iabot, ibtop, ibbot, ictop,
     &    icbot, idetop, idebot, iform, imat, iounit, jform, maxcol,
     &    maxint, maxrow, ndig, output )
 
        iprint = 1
c
c  G=Add/delete a row or column of the matrix,
c    Add a constraint to the tableau.
c
      else if ( leqi ( comnew, 'G' ) ) then
 
        call deladd ( a, iabot, iatop, ibase, ierror, iform, imat,
     &    iounit, line, lpmoda, maxcol, maxdig, maxrow, ncol, ncon,
     &    ndig, nline, nrow, nslak, nvar, output, prompt )
 
        iprint = 1
c
c  H=Help.
c
      else if ( leqi ( comnew, 'H' ) ) then
 
        if ( lpmoda .eq. 0 ) then
          call lahlp1 ( iounit, output )
        else
          call lphlp1 ( iounit, output )
        end if
 
      else if ( leqi ( comnew(1:4), 'HELP' ) ) then
 
        call help ( iounit, output )
c
c  I=Interchange rows I and J.
c
      else if ( leqi ( comnew, 'I' ) ) then
 
        if ( imat .eq. 0 ) then
          ierror = 1
          output = 'You must set up a matrix first!'
          call chrwrt ( iounit, output )
          go to 10
        end if
 
        prompt = 'row I, row J.'
        call intrea ( irow1, line, nline, prompt, iounit, ierror )
        if ( ierror .ne. 0 ) go to 10
 
        call intrea ( irow2, line, nline, prompt, iounit, ierror )
        if ( ierror .ne. 0 ) go to 10
 
        call swprow ( a, iatop, iabot, ibase, ierror, iform, iounit,
     &    irow1, irow2, lpmoda, maxcol, maxrow, ncol, nrow, output )
 
          if ( iform .eq. 0 ) then
            idetop = - idetop
          else if ( iform .eq. 1 ) then
            dete = - dete
          else if ( iform .eq. 2 ) then
            idetop = - idetop
          end if
 
        iprint = 1
c
c  J=Jacobi pre and post multiplication by (I,J) plane rotation.
c
      else if ( leqi ( comnew, 'J' ) ) then
 
        call evjaco ( a, ibase, ierror, iform, imat, iounit, line,
     &    lpmoda, maxcol, maxrow, ncol, nline, nrow, output, prompt )
 
        iprint = 0
c
c  K=Disk file is to be opened or closed.
c
      else if ( leqi ( comnew, 'K' ) ) then
 
        call transc ( filtrn, ierror, iounit, line, nline, output,
     &    prompt )
c
c  L=Change between linear algebra and linear programming modes.
c
      else if ( leqi ( comnew, 'L' ) ) then
 
        call lpset ( ierror, imat, iounit, line, lpmoda, nart, ncol,
     &    ncon, nline, nrow, nslak, nvar, output, prompt )
c
c  M=Multiply row by scalar.
c
      else if ( leqi ( comnew, 'M' ) ) then
 
        call chkmul ( ierror, iform, imat, iounit, irow, istop, isbot,
     &    line, maxdig, ndig, nline, output, prompt, sval )
 
        if ( ierror .eq. 0 ) then
 
          call mulply ( a, dete, iatop, iabot, idetop, idebot, ierror,
     &      iform, iounit, irow, maxcol, maxint, maxrow, ncol, ndig,
     &      nrow, output, sval, istop, isbot )
 
          iprint = 1
 
        end if
c
c  MAXINT = Set maximum integer.
c
      else if ( leqi ( comnew, 'MAXINT' ) ) then
 
        prompt = 'maximum integer for rational representations.'
        call chrdb2 ( prompt )
 
        call intrea ( maxint, line, nline, prompt, iounit, ierror )
c
c  N=Set number of digits.
c
      else if ( leqi ( comnew, 'N' ) ) then
 
        call setdig ( ierror, iounit, line, maxdig, ndig, nline,
     &    output)
c
c  O=Optimality check.
c
      else if ( leqi ( comnew, 'O' ) ) then
 
        if ( lpmoda .eq. 1 ) then
 
          call lpopt ( a, iatop, iabot, ibase, ierror, iform, imat,
     &      iopti, iounit, isltop, islbot, lpmoda, maxcol, maxrow,
     &      nart, ncol, nrow, nslak, nvar, output, sol )
 
        else
 
          call laopt ( a, iabot, iatop, ierror, iform, imat, iounit,
     &      maxcol, maxrow, ncol, nrow, output )
 
        end if
c
c  P=Pivot.
c
      else if ( leqi ( comnew, 'P' ) ) then
 
        iauto = 0
 
        call lppiv ( a, iatop, iabot, iauto, ibase, ierror, iform,
     &    imat, iounit, isltop, islbot, line, lpmoda, maxcol,
     &    maxint, maxrow, nart, ncol, ndig, nline, nrow, nslak, nvar,
     &    output, prompt, sol )
 
        iprint = 1
c
c  Q=Quit.
c  QY= QUIT NOW!
c
      else if ( leqi ( comnew(1:1), 'Q' ) ) then
 
        if ( leqi ( comnew(2:2), 'Y' ) ) then
          isay='Y'
        else
          nline = 0
          prompt = '"Y" to confirm you want to quit.'
          iterm = 0
          call chrrea ( isay, line, nline, prompt, iounit, ierror,
     &      iterm )
          if ( ierror .ne. 0) then
            isay = 'Y'
          end if
        end if
 
        if ( leqi ( isay, 'Y' ) ) then
          output = 'MATMAN is stopping now.'
          call chrwrt ( iounit, output )
 
          if ( iounit(3) .ne. -1 ) then
            call transc ( filtrn, ierror, iounit, line, nline, output,
     &        prompt )
          end if
 
          stop
        end if
c
c  R=Restore matrix.
c
      else if ( leqi ( comnew, 'R' ) ) then
 
        call restor ( a, c, iabot, iatop, ibase, ibasec, icbot, ictop,
     &    ierror, imat, iounit, lpmoda, lpmodc, maxcol, maxrow, nart,
     &    nartc, ncol, ncolc, nrow, nrowc, nslak, nslakc, nvar, nvarc,
     &    output )
 
        if ( ierror .eq. 0 ) then
          iprint = 1
          if ( iform .eq. 0 ) then
            idetop = 1
            idebot = 1
          else if ( iform .eq. 1 ) then
            dete = 1.0
          else if ( iform .eq. 2 ) then
            idetop = 1
            idebot = 0
          end if
        end if
c
c  RATional = use rational arithmetic
c
      else if ( leqi ( comnew(1:3), 'RAT' ) ) then
 
        jform = 0
 
        call form ( a, b, c, dete, iatop, iabot, ibtop, ibbot, ictop,
     &    icbot, idetop, idebot, iform, imat, iounit, jform, maxcol,
     &    maxint, maxrow, ndig, output )
 
        iprint = 1
c
c  REAl = use real arithmetic
c
      else if ( leqi ( comnew(1:3), 'REA' ) ) then
 
        jform = 1
 
        call form ( a, b, c, dete, iatop, iabot, ibtop, ibbot, ictop,
     &    icbot, idetop, idebot, iform, imat, iounit, jform, maxcol,
     &    maxint, maxrow, ndig, output )
 
        iprint = 1
c
c  S=Store a matrix.
c
      else if ( leqi ( comnew, 'S' ) ) then
 
        if ( imat .ne. 1 ) then
          output = 'No matrix has been defined yet!'
          call chrwrt ( iounit, output )
          go to 10
        end if
 
        call copmat ( a, c, iatop, iabot, ictop, icbot, ibase, ibasec,
     &    lpmoda, lpmodc, maxcol, maxrow, nart, nartc, ncol, ncolc,
     &    nrow, nrowc, nslak, nslakc, nvar, nvarc )
 
        output = 'A copy of the matrix has been stored.'
        call chrwrt ( iounit, output )
c
c  T=Type matrix or tableau.
c
      else if ( leqi ( comnew, 'T' ) ) then
 
        call type ( a, iabot, iatop, ibase, ierror, iform, imat,
     &    iounit, lpmoda, maxcol, maxrow, nart, ncol, nrow, output )
c
c  TR = Transpose matrix.
c
      else if ( comnew(1:2) .eq. 'TR' ) then
 
        call chktrn ( ierror, imat, iounit, lpmoda, maxcol, maxrow,
     &    ncol, nrow, output )
 
        if ( ierror .eq. 0 ) then
 
          if ( iform .eq. 0 .or. iform .eq. 2 ) then
            call rattrn ( iatop, iabot, maxcol, maxrow, ncol, nrow )
          else
            call reltrn ( a, maxcol, maxrow, ncol, nrow )
          end if
 
          iprint = 1
 
        end if
c
c  TS = Type linear programming solution.
c
      else if ( comnew(1:2) .eq. 'TS' ) then
 
        call types ( a, iabot, iatop, ibase, ierror, iform, imat,
     &    iounit, islbot, isltop, lpmoda, maxcol, maxrow, nart, ncol,
     &    nrow, nslak, nvar, output, sol )
c
c  U=Undo last command.
c
      else if ( leqi ( comnew, 'U' ) ) then
 
        if ( leqi ( comold, 'K' ) ) then
          comold = 'U'
          comnew = 'K'
          go to 20
        end if
 
        if ( leqi ( comold, 'H' ) ) go to 10
        if ( leqi ( comold, 'HELP' ) ) go to 10
        if ( leqi ( comold, 'N' ) ) go to 10
 
        if ( leqi ( comold, 'L' ) ) then
          comold = 'U'
          comnew = 'L'
          go to 20
        end if
 
        if ( leqi ( comold, 'O' ) ) go to 10
        if ( leqi ( comold, 'T' ) ) go to 10
        if ( leqi ( comold, 'W' ) ) go to 10
 
        call copmat ( b, a, ibtop, ibbot, iatop, iabot, ibaseb,
     &    ibase, lpmodb, lpmoda, maxcol, maxrow, nartb, nart, ncolb,
     &    ncol, nrowb, nrow, nslakb, nslak, nvarb, nvar )
 
        iprint = 1
c
c  V=Remove artificial variables.
c
      else if ( leqi ( comnew, 'V' ) ) then
 
        call lprem ( a, iabot, iatop, ibase, ierror, iform, imat,
     &    iounit, lpmoda, maxcol, maxrow, nart, ncol, nrow, nslak,
     &    nvar, output )
 
        iprint = 1
c
c  W=Write example to file.
c
      else if ( leqi ( comnew, 'W' ) ) then
 
        call filwrt ( a, chineq, filex, iatop, iabot, ierror, iform,
     &    imat, iounit, line, lpmoda, maxcol, maxrow, nart, ncol, 
     &    nline, nrow, nvar, output, prompt )
c
c  X=Read example from file.
c
      else if ( leqi ( comnew, 'X' ) ) then
 
        call filred ( filex, ierror, iform, iounit, line, lpmoda,
     &    nline, output, prompt )
 
        if ( ierror .ne. 0 ) go to 10
 
        comnew = 'e'
        go to 20
c
c  Y=Turn autoprint off or on.
c
      else if ( leqi ( comnew, 'Y' ) ) then
 
        autop = .not. autop
 
        if ( autop ) then
          output = 'Autoprinting turned ON.'
        else
          output = 'Autoprinting turned OFF.'
        end if

        call chrwrt ( iounit, output )
c
c  Z=Automatic reduction.
c
      else if ( leqi ( comnew, 'Z' ) ) then
 
        if ( lpmoda .eq. 0 ) then
 
          call autero ( a, dete, iatop, iabot, ibase, idetop,
     &      idebot, ierror, iform, imat, iounit, maxcol, maxint, 
     &      maxrow, ncol, ndig, nrow, output )
 
        else if ( lpmoda .eq. 1 ) then
 
          iauto = 1
 
          call lppiv ( a, iatop, iabot, iauto, ibase, ierror, iform,
     &      imat, iounit, isltop, islbot, line, lpmoda, maxcol,
     &      maxint, maxrow, nart, ncol, ndig, nline, nrow, nslak, nvar,
     &      output, prompt, sol )

        end if
 
        iprint = 1
c
c  # = Comment.
c  Blank out the input line so MATMAN doesn't reparse it, looking for commands.
c
      else if ( comnew .eq. '#' ) then

        line = ' '
        nline = 0
c
c  $ sign means no paging.
c
      else if ( comnew .eq. '$' ) then

        lpage = 0
        call setpag ( lpage )
        output = 'Paging turned OFF.'
        call chrwrt ( iounit, output )
c
c  % means restore paging.
c
      else if ( comnew .eq. '%' ) then
 
        prompt = 'number of lines to print before pausing.'
        call intrea ( lpage, line, nline, prompt, iounit, ierror )
        if ( ierror .ne. 0 ) go to 10
 
        call setpag ( lpage )
        output = 'Paging turned ON.'
        call chrwrt ( iounit, output )
c
c  < means input from a file.
c
      else if ( comnew .eq. '<' ) then
 
        call infile ( filinp, ierror, iounit, line, lpage, nline,
     &    output, prompt )
c
c ? Extensive help from file.
c
      else if ( comnew .eq. '?' ) then
 
        call hlpvms ( filhlp, iounit, line, nline, output, prompt )
c
c  No match!
c
      else if ( comnew .ne. ' ' ) then

        output = 'You issued the command "' // comnew // '",'
        call chrwrt ( iounit, output )
        output = 'which is not a legal command to MATMAN!'
        call chrwrt ( iounit, output )
        ierror = 1

      end if
c
c  After certain operations, print out the matrix.
c
      if ( 
     &  autop .and. 
     &  ierror .eq. 0 .and.
     &  imat .eq. 1 .and.
     &  iprint .eq. 1 ) then
 
        call type ( a, iabot, iatop, ibase, ierror, iform, imat,
     &    iounit, lpmoda, maxcol, maxrow, nart, ncol, nrow, output )
 
        iprint = 0
 
      end if
 
      go to 10
      end
      subroutine addlin ( )

c*********************************************************************72
c
cc ADDLIN is called whenever a new line is printed.  
c
c  Discussion:
c
c    It simply updates an internal count of the number of lines printed 
c    since the last pause.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
      implicit none

      integer nline

      nline = 0
      call indata ( 'GET', 'NLINE', nline )
      nline = nline + 1
      call indata ( 'SET', 'NLINE', nline )
 
      return
      end
      subroutine autero ( a, dete, iatop, iabot, ibase, idetop, idebot,
     &  ierror, iform, imat, iounit, maxcol, maxint, maxrow, ncol, 
     &  ndig, nrow, output )

c*********************************************************************72
c
cc AUTERO automatically row reduces the current matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, real A(MAXROW,MAXCOL).  A is the matrix to
c    which elementary row operations will be applied.
c
c    Input/output, real DETE, the determinant of the product of the
c    elementary row operations applied to the current matrix.
c
c    Input/output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
c    IATOP and IABOT represent the rational or decimal matrix
c    to which elementary row operations will be applied.
c
c    Input/output, integer IBASE(MAXROW).  IBASE is information
c    really only used by the linear programming routines.
c    AUTERO only needs it because some lower level routines are
c    shared with the linear programming routines.
c
c    Input/output, integer IDETOP, IDEBOT, the rational or
c    decimal representation of the determinant of the product of
c    the elementary row operations applied to the current matrix.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IFORM, specifies the arithmetic being used.
c    0=rational, 1=real, 2=decimal.
c
c    Input, integer IMAT.
c    0, no matrix has been defined by the user.
c    1, a matrix has been defined by the user.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXINT, the maximum integer representable, which
c    should probably be 2147483647.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input, integer NCOL, the number of columns in the matrix.
c
c    Input, integer NDIG, the number of decimal digits used.
c
c    Input, integer NROW, the number of rows in the matrix.
c
c    Workspace, character*100 OUTPUT.
c
      implicit none

      integer maxcol
      integer maxrow

      real a(maxrow,maxcol)
      real amax
      real atemp
      real dete
      integer i
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer ibase(maxrow)
      integer idebot
      integer idetop
      integer ierror
      integer iform
      integer imat
      integer imax
      integer iounit(4)
      integer irow
      integer isbot
      integer istop
      integer j
      integer jcol
      integer krow
      integer l
      integer lpmoda
      integer lrow
      integer maxint
      integer ncol
      integer ndig
      integer nrow
      character*100 output
      real sval

      ierror = 0
 
      if ( imat .ne. 1 ) then
        ierror = 1
        output = 'You must define a matrix first!'
        call chrwrt ( iounit, output )
        return
      end if
 
      do i = 1, nrow
 
        irow = i
 
        do j = 1, ncol
 
          jcol = j
c
c  In column JCOL, seek the row between IROW and NROW with
c  maximum nonzero entry AMAX.
c
          imax = 0
          amax = 0.0
 
          do krow = irow, nrow
 
            if ( iform .eq. 0 ) then
              call ratrel ( atemp, iatop(krow,jcol), iabot(krow,jcol) )
            else if ( iform .eq. 1 ) then
              atemp = a(krow,jcol)
            else if ( iform .eq. 2 ) then
              call decrel ( atemp, iatop(krow,jcol), iabot(krow,jcol) )
            end if
 
            atemp = abs ( atemp )
 
            if ( atemp .gt. amax ) then
              amax = atemp
              imax = krow
            end if
 
          end do
 
          if ( imax .ne. 0 ) then
            krow = imax
            go to 10
          end if
 
        end do
 
        return
 
10      continue
 
        output = ' '
        call chrwrt ( iounit, output )
c
c  Interchange the IROW-th and the pivot rows.
c
        if ( krow .ne. irow ) then
          lpmoda = 0
          call swprow ( a, iatop, iabot, ibase, ierror, iform, iounit,
     &      krow, irow, lpmoda, maxcol, maxrow, ncol, nrow, output )
          dete = - dete
          idetop = - idetop
        end if
c
c  Divide the pivot row by A(IROW,JCOL) so that A(IROW,JCOL) = 1.
c
        if ( iform .eq. 0 ) then
 
          istop = iatop(irow,jcol)
          isbot = iabot(irow,jcol)
 
          call ratdiv ( idebot, idebot, isbot, ierror, idetop,
     &      idetop, istop, maxint )
 
        else if ( iform .eq. 1 ) then
 
          sval = a(irow,jcol)
 
          dete = dete / sval
 
        else if ( iform .eq. 2 ) then
 
          istop = iatop(irow,jcol)
          isbot = iabot(irow,jcol)
 
          call decdiv ( idebot, idebot, isbot, ierror, idetop, idetop,
     &      istop, ndig )
 
        end if
 
        call scadiv ( a, iatop, iabot, ierror, iform, iounit, irow,
     &    maxcol, maxint, maxrow, ncol, ndig, nrow, output, sval, 
     &    istop, isbot )
c
c  Annihilate A(L,JCOL) for L not equal to IROW.
c
        do l = 1, nrow
 
          lrow = l
 
          if ( lrow .ne. irow ) then
 
            if ( iform .eq. 0 ) then

              if ( iatop(lrow,jcol) .ne. 0 ) then

                istop = - iatop(lrow,jcol)
                isbot = iabot(lrow,jcol)

                call rowadd ( a, iatop, iabot, ierror, iform, iounit,
     &            lrow, irow, maxcol, maxint, maxrow, ncol, ndig, 
     &            output, sval, istop, isbot )

                iatop(lrow,jcol) = 0
                iabot(lrow,jcol) = 1

              end if

            else if ( iform .eq. 1 ) then

              if ( a(lrow,jcol) .ne. 0.0 ) then

                sval = -a(lrow,jcol)

                call rowadd ( a, iatop, iabot, ierror, iform, iounit,
     &            lrow, irow, maxcol, maxint, maxrow, ncol, ndig, 
     &            output, sval, istop, isbot )

                a(lrow,jcol) = 0.0

              end if

            else if ( iform .eq. 2 ) then

              if ( iatop(lrow,jcol) .ne. 0 ) then

                istop = -iatop(lrow,jcol)
                isbot = iabot(lrow,jcol)

                call rowadd ( a, iatop, iabot, ierror, iform, iounit,
     &            lrow, irow, maxcol, maxint, maxrow, ncol, ndig, 
     &            output, sval, istop, isbot )

                iatop(lrow,jcol) = 0
                iabot(lrow,jcol) = 0

              end if

            end if
 
          end if
 
        end do
 
      end do
 
      return
      end
      subroutine basic ( ibase, ierror, imat, iounit, line, lpmoda, 
     &  maxrow, nart, nline, nrow, nslak, nvar, output, prompt )

c*********************************************************************72
c
cc BASIC assigns a row of the tableau to one of the basic variables.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer IBASE(MAXROW), keeps track of the basic
c    variables.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IMAT.
c    0, no matrix has been defined by the user.
c    1, a matrix has been defined by the user.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Workspace, character*80 LINE.
c    Used to hold the user's input.
c
c    Input, integer LPMODA.
c    0, the program is in linear algebra mode.
c    1, the program is in linear programming mode.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input, integer NART, the number of artificial variables.
c
c    Input/output, integer NLINE.
c    Keeps track of the number of useful characters in LINE.
c
c    Input, integer NROW, the number of rows in the matrix.
c
c    Input, integer NSLAK, the number of slack variables.
c
c    Input, integer NVAR, the number of basic variables.
c
c    Workspace, character*100 OUTPUT.
c
c    Workspace, character*80 PROMPT.
c
      implicit none

      integer maxrow

      character*6 chrint
      character*22 chrtmp
      integer ibase(maxrow)
      integer ierror
      integer imat
      integer iounit(4)
      integer irow
      integer ivar
      character*80 line
      integer lpmoda
      integer nart
      integer nline
      integer nrow
      integer nslak
      integer nvar
      character*100 output
      character*80 prompt

      if ( imat .eq. 0 ) then
        ierror = 1
        output = 'You must set up a matrix first!'
        call chrwrt ( iounit, output )
        return
      end if
 
      if ( lpmoda .ne. 1 ) then
        ierror = 1
        output = 'Error!  You must be in linear programming mode!'
        call chrwrt ( iounit, output )
        return
      end if
 
      prompt = 'row I, basic variable J.'
c
c  Get the row number I.
c
      call intrea ( irow, line, nline, prompt, iounit, ierror )
      if ( ierror .ne. 0 ) return
 
      if ( irow .lt. 1 .or. irow .gt. nrow ) then
        output = 'Error!  Illegal row number!'
        call chrwrt ( iounit, output )
        ierror = 1
        return
      end if
c
c  Get the basic variable index J.
c
      call intrea ( ivar, line, nline, prompt, iounit, ierror )
      if ( ierror .ne. 0 ) return
 
      if ( ivar .lt. 1 .or. ivar .gt. nvar+nslak+nart ) then
        output = 'Error!  Illegal basic variable number!'
        call chrwrt ( iounit, output )
        ierror = 1
        return
      end if
 
      ibase(irow) = ivar
 
      chrtmp = chrint(ivar)
      output = 'Assigning row ' // chrint(irow) // ' to basic variable'
     &  // chrtmp
      call chrdb2 ( output )
      call chrwrt ( iounit, output )
 
      return
      end
      subroutine capchr ( string )
c
c*********************************************************************72
c
cc CAPCHR replaces any lowercase letters by uppercase ones.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character*(*) STRING, is the string of characters to
c    be transformed.
c
      implicit none

      integer i
      integer itemp
      integer nchar
      character*(*) string

      nchar = len ( string )
 
      do i = 1, nchar
 
        itemp = ichar ( string(i:i) )
 
        if ( 97 .le. itemp .and. itemp .le. 122 ) then
          string(i:i) = char ( itemp-32 )
        end if
 
      end do
 
      return
      end
      subroutine change ( a, iatop, iabot, ierror, iform, imat, iounit,
     &  line, maxcol, maxdig, maxrow, ncol, ndig, nline, nrow, output,
     &  prompt )

c*********************************************************************72
c
cc CHANGE allows the user to change an entry in the array.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, real A(MAXROW,MAXCOL).  A is the matrix
c    whose entry is to be changed.
c
c    Input/output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
c    IATOP and IABOT represent the rational or decimal matrix
c    whose entry is to be changed.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IFORM, specifies the arithmetic being used.
c    0=rational, 1=real, 2=decimal.
c
c    Input, integer IMAT.
c    0, no matrix has been defined by the user.
c    1, a matrix has been defined by the user.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Workspace, character*80 LINE.
c    Used to hold the user's input.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXDIG, the maximum number of decimal digits
c    allowed.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input, integer NCOL, the number of columns in the matrix.
c
c    Input, integer NDIG, the number of decimal digits used.
c
c    Input/output, integer NLINE.
c    Keeps track of the number of useful characters in LINE.
c
c    Input, integer NROW, the number of rows in the matrix.
c
c    Workspace, character*100 OUTPUT.
c
c    Workspace, character*80 PROMPT.
c
      implicit none

      integer maxcol
      integer maxrow

      real a(maxrow,maxcol)
      character*22 chldec
      character*22 chlrat
      character*6 chrint
      character*14 chrrel
      character*22 chrtmp
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer icol
      integer ierror
      integer iform
      integer igcf
      integer imat
      integer iounit(4)
      integer irow
      integer isbot
      integer istop
      integer itemp
      character*80 line
      integer maxdig
      integer ncol
      integer ndig
      integer nline
      integer nrow
      character*100 output
      character*80 prompt
      real rval

      if ( imat .eq. 0 ) then
        ierror = 1
        output = 'You must set up a matrix first!'
        call chrwrt ( iounit, output )
        return
      end if
 
      prompt = 'row I, column J, new value S.'
c
c  Get the row number.
c
      call intrea ( irow, line, nline, prompt, iounit, ierror )
      if ( ierror .ne. 0 ) return
 
      if ( irow .lt. 1 .or. irow .gt. nrow ) then
        output = 'Error!  Illegal row value!'
        call chrwrt ( iounit, output )
        ierror = 1
        return
      end if
c
c  Get the column number.
c
      call intrea ( icol, line, nline, prompt, iounit, ierror )
      if ( ierror .ne. 0 ) return
 
      if ( icol .lt. 1 .or. icol .gt. ncol ) then
        output = 'Error!  Illegal column value!'
        call chrwrt ( iounit, output )
        ierror = 1
        return
      end if
c
c  Read the value.
c
      if ( iform .eq. 0 ) then
 
        call ratrea ( istop, isbot, rval, line, nline, prompt, iounit,
     &    ierror )
        if ( ierror .ne. 0 ) return
 
        chrtmp = chlrat ( istop, isbot )
        output = 'Change entry ' // chrint(irow) // ',' // chrint(icol)
     &    // ' to ' // chrtmp
        call chrdb2 ( output )
        call chrwrt ( iounit, output )
        itemp = igcf ( istop, isbot )
        iatop(irow,icol) = istop / itemp
        iabot(irow,icol) = isbot / itemp
 
      else if ( iform .eq. 1 ) then
 
        call relrea ( rval, line, nline, prompt, iounit, ierror )
        if ( ierror .ne. 0 ) return
        a(irow,icol) = rval
        output = 'Change entry '//chrint(irow)//','//chrint(icol)//
     &    ' to '//chrrel(rval)
        call chrdb2 ( output )
        call chrwrt ( iounit, output )
 
      else if ( iform .eq. 2 ) then
 
        call decrea ( istop, isbot, rval, line, maxdig, nline, prompt,
     &    iounit, ierror )
 
        if ( ierror .ne. 0 ) return
 
        call deccut ( istop, isbot, ndig )
 
        chrtmp = chldec ( istop, isbot )
        output = 'Change entry ' // chrint(irow) // ',' // 
     &    chrint(icol) // ' to ' // chrtmp
        call chrdb2 ( output )
        call chrwrt ( iounit, output )
        iatop(irow,icol) = istop
        iabot(irow,icol) = isbot
 
      end if
 
      return
      end
      subroutine chkadd ( ierror, iform, imat, iounit, irow1, irow2, 
     &  istop, isbot, line, maxdig, ndig, nline, nrow, output, prompt, 
     &  sval )

c*********************************************************************72
c
cc CHKADD checks the command to add a multiple of one row to another.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IFORM, specifies the arithmetic being used.
c    0=rational, 1=real, 2=decimal.
c
c    Input, integer IMAT.
c    0, no matrix has been defined by the user.
c    1, a matrix has been defined by the user.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Input, integer IROW1, the row to which the multiple is to be added.
c
c    Input, integer IROW2, the row which is to be multiplied and
c    added to another row.
c
c    Output, integer ISTOP, ISBOT, the parts of the rational
c    or decimal fraction of the multiplier, if that is the
c    arithmetic being used.
c
c    Workspace, character*80 LINE.
c    Used to hold the user's input.
c
c    Input, integer MAXDIG, the maximum number of decimal digits allowed.
c
c    Input, integer NDIG, the number of decimal digits used.
c
c    Input/output, integer NLINE.
c    Keeps track of the number of useful characters in LINE.
c
c    Input, integer NROW, the number of rows in the matrix.
c
c    Workspace, character*100 OUTPUT.
c
c    Workspace, character*80 PROMPT.
c
c    Output, real SVAL, the multiplier, if real arithmetic is used.
c
      implicit none

      integer ierror
      integer iform
      integer imat
      integer iounit(4)
      integer irow1
      integer irow2
      integer isbot
      integer istop
      character*80 line
      integer maxdig
      integer ndig
      integer nline
      integer nrow
      character*100 output
      character*80 prompt
      real sval

      if ( imat .eq. 0 ) then
        ierror = 1
        output = 'You must set up a matrix first!'
        call chrwrt ( iounit, output )
        return
      end if
 
      prompt = 'multiplier S, row I to add, target row J.'
c
c  Get the multiplier, SVAL or ISTOP/ISBOT.
c
      if ( iform .eq. 0 ) then
 
        call ratrea ( istop, isbot, sval, line, nline, prompt, iounit,
     &    ierror )
        if ( ierror .ne. 0 ) return
 
      else if ( iform .eq. 1 ) then
 
        call relrea ( sval, line, nline, prompt, iounit, ierror )
        if ( ierror .ne. 0 ) return
 
      else if ( iform .eq. 2 ) then
 
        call decrea ( istop, isbot, sval, line, maxdig, nline, prompt,
     &    iounit, ierror )
 
        if ( ierror .ne. 0 ) return
 
        call deccut ( istop, isbot, ndig )
 
      end if
c
c  Get the row to add, IROW2.
c
      call intrea ( irow2, line, nline, prompt, iounit, ierror )
      if ( ierror .ne. 0 ) return
 
      if ( irow2 .lt. 1 .or. irow2 .gt. nrow ) then
        ierror = 1
        output = 'Error!  Row index was not acceptable!'
        call chrwrt ( iounit, output )
        return
      end if
c
c  Get the row to which we are adding, IROW1.
c
      call intrea ( irow1, line, nline, prompt, iounit, ierror )
      if ( ierror .ne. 0 ) return
 
      if ( irow1 .lt. 1 .or. irow1 .gt. nrow ) then
        ierror = 1
        output = 'Error!  Row index was not acceptable!'
        call chrwrt ( iounit, output )
        return
      end if
c
c  Make sure the rows are different.
c
      if ( irow1 .eq. irow2 ) then
        output = 'Error!  The rows should not be the same!'
        call chrwrt ( iounit, output )
        ierror = 1
        return
      end if
 
      ierror = 0
 
      return
      end
      subroutine chkdet ( ierror, imat, iounit, lpmoda, ncol, nrow, 
     &  output )

c*********************************************************************72
c
cc CHKDET checks that a request for a determinant can be carried out.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IMAT.
c    0, no matrix has been defined by the user.
c    1, a matrix has been defined by the user.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Input, integer LPMODA.
c    0, the program is in linear algebra mode.
c    1, the program is in linear programming mode.
c
c    Input, integer NCOL, the number of columns in the matrix.
c
c    Input, integer NROW, the number of rows in the matrix.
c
c    Workspace, character*100 OUTPUT.
c
      implicit none

      integer ierror
      integer imat
      integer iounit(4)
      integer lpmoda
      integer ncol
      integer nrow
      character*100 output

      if ( lpmoda .ne. 0 ) then
        ierror = 1
        output = 'Error!  You must get into linear algebra mode with'
        call chrwrt ( iounit, output )
        output = 'the "L" command before asking for a determinant.'
        call chrwrt ( iounit, output )
        return
      end if
 
      if ( imat .eq. 0 ) then
        ierror = 1
        output = 'Error!  You must define a matrix with the "E" command'
        call chrwrt ( iounit, output )
        output = 'before asking for its determinant.'
        call chrwrt ( iounit, output )
        return
      end if
 
      if ( nrow .ne. ncol ) then
        ierror = 1
        output = 'Error!  A matrix must be square in order for you to'
        call chrwrt ( iounit, output )
        output = 'ask for its determinant.'
        call chrwrt ( iounit, output )
        return
      end if
 
      ierror = 0
 
      return
      end
      subroutine chkero ( comnew, ierror, iounit, line2, output )

c*********************************************************************72
c
cc CHKERO checks for commands given in the form of ERO's.
c
c
c  1) The row interchange command  RI1 <=> RI2
c     Note that this will fail if user types "R I1 <=> R I2"
c
c  2a) The scalar multiply command
c     RI1 <= S * RI1
c  with or without the "*".
c
c  2b) The scalar divide command
c     RI1 <= RI1 / S
c
c  3) The add row command:
c     RI1 <= RI1 + S * RI2
c  or
c     RI1 <= S * RI2 + RI1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character*4 COMNEW.
c    If CHKERO decides that the user has input an ERO in the
c    natural format, then COMNEW contains the necessary
c    one letter MATMAN command to carry out the ERO.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Workspace, character*80 LINE2.
c    Used to hold a copy of the user input normally kept in LINE.
c
c    Workspace, character*100 OUTPUT.
c
      character*22 chlrat
      character*6 chrint
      character*20 comnew
      integer idbot2
      integer idbot3
      integer idtop2
      integer idtop3
      integer ierror
      integer iounit(4)
      integer irow1
      integer irow2
      integer irow3
      integer isbot2
      integer isbot3
      integer istop2
      integer istop3
      integer itemp
      integer lchar
      logical ldiv
      integer lenchr
      character*80 line2
      integer nline
      character*100 output
      character*80 string

      comnew = ' '
c
c  1. Remove all blanks from the line, and capitalize it.
c
      call chrdb1 ( line2 )
      call capchr ( line2 )
      nline = lenchr ( line2 )
c
c  2. Is the first character an "R" or "ROW"?
c
      if ( line2(1:1) .ne. 'R' ) return
 
      if ( line2(1:3) .eq. 'ROW' ) then
        call chrchp ( line2, 1, 3 )
      else
        call chrchp ( line2, 1, 1 )
      end if

      nline = lenchr ( line2 )
c
c  3. The next item should be a row number, IROW1.
c
      call chrcti ( line2, irow1, ierror, lchar )
 
      if ( ierror .ne. 0 ) then
        output = 'Your ERO command could not be understood.'
        call chrwrt ( iounit, output )
        output = 'The first row number "R1" did not make sense.'
        call chrwrt ( iounit, output )
        return
      end if
 
      call chrchp ( line2, 1, lchar )
      nline = lenchr ( line2 )
c
c  4. Check for the row interchange string "=", "<>", "<=>" or "<->".
c
      if ( line2(1:2) .eq. '<>' ) then
        string = '<>'
      else if ( line2(1:3) .eq. '<=>' ) then
        string = '<=>'
      else if ( line2(1:3) .eq. '<->' ) then
        string = '<->'
      else if ( line2(1:2) .eq. '<=' ) then
        string = '<='
      else if ( line2(1:2) .eq. '<-' ) then
        string = '<-'
      else if ( line2(1:2) .eq. '=>' ) then
        string = '=>'
      else if ( line2(1:2) .eq. '->' ) then
        string = '->'
      else if ( line2(1:1) .eq. '=' ) then
        string = '='
      else if ( line2(1:2) .eq. ':=' ) then
        string = ':='
      else
        ierror = 1
        output = 'Your ERO command could not be understood.'
        call chrwrt ( iounit, output )
        output = 'The assignment symbol <=> was missing.'
        call chrwrt ( iounit, output )
        return
      end if
 
      lchar = lenchr ( string )
 
      call chrchp ( line2, 1, lchar )
      nline = lenchr ( line2 )
c
c  5. The next quantity could be a possible signed scalar, S2,
c     or an implicit +-1.
c
      if ( line2(1:1) .eq. 'R' ) then

        istop2 = 1.0
        isbot2 = 1.0

      else

        if ( line2(1:2) .eq. '+R' ) then
          istop2 = 1.0
          isbot2 = 1.0
          call chrchp ( line2, 1, 1 )
          nline = lenchr ( line2 )
        else if ( line2(1:2) .eq. '-R' ) then
          istop2 = - 1.0
          isbot2 = 1.0
          call chrchp ( line2, 1, 1 )
          nline = lenchr ( line2 )
        else
          call chrctg ( line2, istop2, isbot2, ierror, lchar )
          call chrchp ( line2, 1, lchar )
          nline = lenchr ( line2 )
 
          if ( ierror .ne. 0 ) then
            output = 'Your ERO command could not be understood.'
            call chrwrt ( iounit, output )
            output = 'The multiplier S2 did not make sense.'
            call chrwrt ( iounit, output )
            ierror = 1
            return
          end if
 
        end if
      end if
c
c  6. Is the next character an optional "*"?
c
      if ( line2(1:1) .eq. '*' ) then
        call chrchp ( line2, 1, 1 )
        nline = lenchr ( line2 )
      end if
c
c  7. Is the next character an "R"?
c
      if ( line2(1:3) .eq. 'ROW' ) then
        call chrchp ( line2, 1, 3 )
        nline = lenchr ( line2 )
      else if ( line2(1:1) .eq. 'R' ) then
        call chrchp ( line2, 1, 1 )
        nline = lenchr ( line2 )
      else
        ierror = 1
        output = 'Your ERO command could not be understood.'
        call chrwrt ( iounit, output )
        output = 'Could not find the second row index.'
        call chrwrt ( iounit, output )
        return
      end if
c
c  8. The next item should be a row number, IROW2.
c
      call chrcti ( line2, irow2, ierror, lchar )
 
      if ( ierror .ne. 0 ) then
        output = 'Your ERO command could not be understood.'
        call chrwrt ( iounit, output )
        output = 'The second row number "R2" did not make sense.'
        call chrwrt ( iounit, output )
        ierror = 1
        return
      else
        call chrchp ( line2, 1, lchar )
        nline = lenchr ( line2 )
      end if
c
c  9. If there's nothing more, this must be an interchange
c     or a scaling.  Form the equivalent MATMAN M or I command.
c
      if ( nline .eq. 0 ) then
 
        if ( irow1 .eq. irow2 ) then
          comnew = 'm'
          if ( isbot2 .lt. 0 ) then
            isbot2 = - isbot2
            istop2 = - istop2
          end if
          line2 = chrint(irow1) // ' ' // chlrat(istop2,isbot2)
          nline = lenchr ( line2 )
          return
        end if
 
        if ( istop2 .eq. 1 .and. isbot2 .eq. 1 ) then
          comnew = 'i'
          line2 = chrint(irow1) // ' ' // chrint(irow2)
          call chrdb2 ( line2 )
          nline = lenchr ( line2 )
          return
        end if
 
        ierror = 1
        output = 'Your ERO command could not be understood.'
        call chrwrt ( iounit, output )
        output = 'A MULTIPLY command must have R1 and R2 the same.'
        call chrwrt ( iounit, output )
        output = 'An INTERCHANGE command cannot have a multiplier.'
        call chrwrt ( iounit, output )
        return
      end if
c
c  10. Is the next quantity a '/', or perhaps a '*'?
c
      ldiv = .false.
 
      if ( line2(1:1) .eq. '/' ) then
 
        ldiv = .true.
        call chrchp ( line2, 1, 1 )
        nline = lenchr ( line2 )
        call chrctg ( line2, idtop2, idbot2, ierror, lchar )
 
        if ( ierror .ne. 0 ) then
          output = 'Your ERO command could not be understood.'
          call chrwrt ( iounit, output )
          output = 'The divisor of row 2 did not make sense.'
          call chrwrt ( iounit, output )
          return
        end if
 
        istop2 = istop2 * idbot2
        isbot2 = isbot2 * idtop2
 
        if ( irow1 .eq. irow2 ) then

          if ( ldiv ) then
            comnew = 'd'
            itemp = istop2
            istop2 = isbot2
            isbot2 = itemp
          else
            comnew = 'm'
          end if

          if ( isbot2 .lt. 0 ) then
            isbot2 = - isbot2
            istop2 = - istop2
          end if

          line2 = chrint(irow1) // ' ' // chlrat(istop2,isbot2)
          nline = lenchr(line2)
          return

        end if
 
      else if ( line2(1:1) .eq. '*' ) then
 
        call chrchp ( line2, 1, 1 )
        nline = lenchr ( line2 )
        call chrctg ( line2, idtop2, idbot2, ierror, lchar )
 
        if ( ierror .ne. 0 ) then
          output = 'Your ERO command could not be understood.'
          call chrwrt ( iounit, output )
          output = 'The multiplier of row 2 did not make sense.'
          call chrwrt ( iounit, output )
          return
        end if
 
        istop2 = istop2 * idtop2
        isbot2 = isbot2 * idbot2
 
        if ( irow1 .eq. irow2 ) then
          comnew = 'm'
          if ( isbot2 .lt. 0 ) then
            isbot2 = - isbot2
            istop2 = - istop2
          end if
          line2 = chrint(irow1)//' '//chlrat(istop2,isbot2)
          nline = lenchr(line2)
          return
        end if
 
      end if
c
c  11. Is the next quantity a scalar, S3?
c
      if ( line2(1:2) .eq. '+R' ) then
 
        istop3 = 1.0
        isbot3 = 1.0
        call chrchp ( line2, 1, 1) 
        nline = lenchr(line2)
 
      else if ( line2(1:2) .eq. '-R' ) then
 
        istop3 = - 1.0
        isbot3 = 1.0
        call chrchp ( line2, 1, 1 )
        nline = lenchr ( line2 )
 
      else
 
        call chrctg ( line2, istop3, isbot3, ierror, lchar )
 
        if ( ierror .ne. 0 ) then
          output = 'Your ERO command could not be understood.'
          call chrwrt ( iounit, output )
          output = 'The multiplier S2 did not make sense.'
          call chrwrt ( iounit, output )
          ierror = 1
          return
        end if
 
        call chrchp ( line2, 1, lchar )
        nline = lenchr ( line2 )
 
      end if
c
c  12. Is the next quantity an optional "*"?
c
      if ( line2(1:1) .eq. '*' ) then
        call chrchp ( line2, 1, 1 )
        nline = lenchr(line2)
      end if
c
c  13. Is the next quantity an "R" or ROW?
c
      if ( line2(1:3) .eq. 'ROW' ) then
        call chrchp ( line2, 1, 3 )
        nline = lenchr(line2)
      else if ( line2(1:1) .eq. 'R' ) then
        call chrchp ( line2, 1, 1) 
        nline = lenchr(line2)
      else
        ierror = 1
        output = 'Your ERO command could not be understood.'
        call chrwrt ( iounit, output )
        output = 'The "R" marking the third row was misplaced.'
        call chrwrt ( iounit, output )
        return
      end if
c
c  14. The next item should be a row number, IROW3.
c
      call chrcti ( line2, irow3, ierror, lchar )
 
      if ( ierror .ne. 0 ) then
        output = 'Your ERO command could not be understood.'
        call chrwrt ( iounit, output )
        output = 'The third row number "R3" did not make sense.'
        call chrwrt ( iounit, output )
        ierror = 1
        return
      end if
 
      call chrchp ( line2, 1, lchar )
c
c  15. Is the next quantity a '/', or perhaps a '*'?
c
      if ( line2(1:1) .eq. '/' ) then
 
        call chrchp ( line2, 1, 1 )
        nline = lenchr ( line2 )
        call chrctg ( line2, idtop3, idbot3, ierror, lchar )
 
        if ( ierror .ne. 0 ) then
          output = 'Your ERO command could not be understood.'
          call chrwrt ( iounit, output )
          output = 'The divisor of row 3 did not make sense.'
          call chrwrt ( iounit, output )
          return
        end if
 
        istop3 = istop3 * idbot3
        isbot3 = isbot3 * idtop3
 
      else if ( line2(1:1) .eq. '*' ) then
 
        call chrchp ( line2, 1, 1 )
        nline = lenchr ( line2 )
        call chrctg ( line2, idtop3, idbot3, ierror, lchar )
 
        if ( ierror .ne. 0 ) then
          output = 'Your ERO command could not be understood.'
          call chrwrt ( iounit, output )
          output = 'The multiplier of row 3 did not make sense.'
          call chrwrt ( iounit, output )
          return
        end if
 
        istop3 = istop3 * idtop3
        isbot3 = isbot3 * idbot3
 
      end if
c
c  16. Form the equivalent MATMAN ADD command.
c
      if ( irow1 .eq. irow2 ) then

        comnew = 'a'

        if ( isbot3 .lt. 0 ) then
          isbot3 = - isbot3
          istop3 = - istop3
        end if

        line2 = chlrat(istop3,isbot3)//' '//chrint(irow3)//' '
     &    //chrint(irow1)
        call chrdb2 ( line2 )
        nline = lenchr ( line2 )

      else if ( irow1 .eq. irow3 ) then

        comnew = 'a'

        if ( isbot2 .lt. 0 ) then
          isbot2 = - isbot2
          istop2 = - istop2
        end if

        line2 = chlrat(istop2,isbot2)//' '//chrint(irow2)//' '
     &    //chrint(irow1)

        call chrdb2 ( line2 )
        nline = lenchr ( line2 )

      else

        ierror = 1
        output = 'Your ERO command could not be understood.'
        call chrwrt ( iounit, output )
        output = 'R2 or R3 must equal R1 in an ERO command.'
        call chrwrt ( iounit, output )
      end if
 
      return
      end
      subroutine chkmul ( ierror, iform, imat, iounit, irow, istop, 
     &  isbot, line, maxdig, ndig, nline, output, prompt, rval )

c*********************************************************************72
c
cc CHKMUL checks a command to multiply a row by a scalar.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IFORM, specifies the arithmetic being used.
c    0=rational, 1=real, 2=decimal.
c
c    Input, integer IMAT.
c    0, no matrix has been defined by the user.
c    1, a matrix has been defined by the user.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Output, integer IROW, the row to be multiplied.
c
c    Output, integer ISTOP, ISBOT, the multiplier to use for
c    fractional or decimal arithmetic.
c
c    Workspace, character*80 LINE.
c    Used to hold the user's input.
c
c    Input, integer MAXDIG, the maximum number of decimal digits
c    allowed.
c
c    Input, integer NDIG, the number of decimal digits used.
c
c    Input/output, integer NLINE.
c    Keeps track of the number of useful characters in LINE.
c
c    Workspace, character*100 OUTPUT.
c
c    Workspace, character*80 PROMPT.
c
c    Output, real RVAL, the multiplier to use for real arithmetic.
c
      implicit none

      integer ierror
      integer iform
      integer imat
      integer iounit(4)
      integer irow
      integer isbot
      integer istop
      character*80 line
      integer maxdig
      integer ndig
      integer nline
      character*100 output
      character*80 prompt
      real rval

      if ( imat .eq. 0 ) then
        ierror = 1
        output = 'You must set up a matrix first!'
        call chrwrt ( iounit, output )
        return
      end if
 
      prompt = 'row I, multiplier S.'
c
c  Read the row number to be multiplied.
c
      call intrea ( irow, line, nline, prompt, iounit, ierror )
      if ( ierror .ne. 0 ) return
c
c  Read the multiplier, either RVAL or ISTOP/ISBOT.
c
      if ( iform .eq. 0 ) then
 
        call ratrea ( istop, isbot, rval, line, nline, prompt, iounit,
     &    ierror )
 
      else if ( iform .eq. 1 ) then
 
        call relrea ( rval, line, nline, prompt, iounit, ierror )
 
      else if ( iform .eq. 2 ) then
 
        call decrea ( istop, isbot, rval, line, maxdig, nline, prompt,
     &    iounit, ierror )
 
        call deccut ( istop, isbot, ndig )
 
      end if
 
      return
      end
      subroutine chktrn ( ierror, imat, iounit, lpmoda, maxcol, maxrow,
     &  ncol, nrow, output )

c*********************************************************************72
c
cc CHKTRN checks that a request to transpose a matrix can be carried out.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IMAT.
c    0, no matrix has been defined by the user.
c    1, a matrix has been defined by the user.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Input, integer LPMODA.
c    0, the program is in linear algebra mode.
c    1, the program is in linear programming mode.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input, integer NCOL, the number of columns in the matrix.
c
c    Input, integer NROW, the number of rows in the matrix.
c
c    Workspace, character*100 OUTPUT.
c
      implicit none

      integer maxcol
      integer maxrow

      integer ierror
      integer imat
      integer iounit(4)
      integer lpmoda
      integer ncol
      integer nrow
      character*100 output

      ierror = 0
c
c  The user must have entered a matrix.
c
      if ( imat .ne. 1 ) then
 
        output = ' '
        call chrwrt ( iounit, output )
        output = 'Error!  You must set up a matrix with the'
        call chrwrt ( iounit, output )
        output = '"E" command before you can transpose it!'
        call chrwrt ( iounit, output )
 
        ierror = 1
 
        return
      end if
c
c  The user must be in linear algebra mode.
c
      if ( lpmoda .ne. 0 ) then
        output = 'Error!  You must be in linear algebra mode'
        call chrwrt ( iounit, output )
        output = 'in order to transpose a matrix!'
        call chrwrt ( iounit, output )
        output = 'Use the "L" command to switch over!'
        call chrwrt ( iounit, output )
        ierror = 1
        return
      end if
c
c  The matrix can't have too many columns or rows.
c
      if ( nrow .gt. maxcol ) then
        output = 'Error!'
        call chrwrt ( iounit, output )
        output = 'The matrix has too many rows to transpose.'
        call chrwrt ( iounit, output )
        ierror = 1
        return
      end if
 
      if ( ncol .gt. maxrow ) then
        output = 'Error!'
        call chrwrt ( iounit, output )
        output = 'The matrix has too many columns to transpose!'
        call chrwrt ( iounit, output )
        ierror = 1
        return
      end if
 
      return
      end
      function chldec ( ival, jval )

c*********************************************************************72
c
cc CHLDEC returns a left-justified representation of IVAL * 10**JVAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer IVAL, JVAL, the two integers which represent the 
c    decimal.
c
c    Output, character*22 CHLDEC, a left-justified string
c    containing the representation of the decimal.
c
      character*22 chldec
      character*22 chrtmp
      integer i
      integer icopy
      integer ival
      integer ival1
      integer ival2
      integer jval
      integer l0
      integer l1
      integer l2
      integer l3
      integer l4
      integer nd
      integer nd1
      integer nd2

      chrtmp = ' '
c
c  Take care of case where IVAL is 0.
c
      if ( ival .eq. 0 ) then
        chldec = '0'
        return
      end if
c
c  Take care of case where JVAL is 0.
c
      if ( jval .eq. 0 ) then
        write ( chrtmp, '(i22)' ) ival
        call chrdb1 ( chrtmp )
        call chrdb1 ( chrtmp )
        chldec = chrtmp
        return
      end if
c
c  Count the digits in IVAL.
c
      nd = 0
      icopy = abs ( ival )
 
10    continue
 
      if ( icopy .gt. 0 ) then
        icopy = icopy / 10
        nd = nd + 1
        go to 10
      end if
c
c  If JVAL is greater than 0:
c
      if ( jval .gt. 0 ) then

        if ( ival .gt. 0 ) then
          l1 = nd
        else
          l1 = nd + 1
        end if

        write ( chrtmp, '(i22)' ) ival
        call chrdb1 ( chrtmp )

        do i = l1 + 1, l1 + jval
          chrtmp(i:i) = '0'
        end do

        chldec = chrtmp
        return

      end if
c
c  JVAL is negative.
c
      ival1 = abs ( ival ) / 10**(-jval)
      ival2 = abs ( ival ) - 10**(-jval) * ival1
 
      nd1 = 0
      icopy = abs ( ival1 )
 
20    continue
 
      if ( icopy .gt. 0 ) then
        icopy = icopy/10
        nd1 = nd1 + 1
        go to 20
      end if
 
      if ( nd1 .eq. 0 ) then
        nd1 = 1
      end if
 
      nd2 = 0
      icopy = abs ( ival2 )
 
30    continue
 
      if ( icopy .gt. 0 ) then
        icopy = icopy/10
        nd2 = nd2 + 1
        go to 30
      end if
 
      chrtmp = ' '
 
      if ( ival .lt. 0 ) then
        l0 = 2
      else
        l0 = 1
      end if
 
      l1 = l0 + nd1 - 1
      l2 = l1 + 1
      l3 = l2 + 1
      l4 = l3 + abs ( jval ) - 1
  
      if ( ival .lt. 0 ) then
        chrtmp(1:1) = '-'
      end if
 
      call chritc0 ( chrtmp(l0:l1), ival1 )

      chrtmp(l2:l2) = '.'

      call chritc0 ( chrtmp(l3:l4), ival2 )
 
      chldec = chrtmp
 
      return
      end
      function chlint ( intval )

c*********************************************************************72
c
cc CHLINT returns a 6-character representation of an integer.
c
c
c  The representation is left justified.  The representation is 
c  '******' if the integer is too large or negative to fit in six 
c  positions.
c
c  Compare CHRINT and CHR0NT.
c
c  Examples:
c
c    INTVAL  STRING
c
c         1  1
c        -1  -1
c         0  0
c      1952  1952
c    123456  123456
c   1234567  ******  <-- Not enough room!
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer INTVAL, an integer variable to be converted.
c
c    Output (through function value), character*6 CHLINT,
c    a 6 character representation of the integer, left
c    justified.  Thus, if INTVAL = 1, CHLINT = '1     '.
c    CHLINT must be declared "character*6 CHLINT" in the
c    calling program.
c
      character*6 chlint
      character*6 chrtmp
      integer i
      integer idig
      integer intval
      integer ipos
      integer ival

      chlint = ' '
 
      do i = 1, 6
        chrtmp(i:i) = ' '
      end do
 
      ival = abs ( intval )
 
      ipos = 6
 
10    continue
c
c  If we ever run out of room in STRING, then overwrite it
c  with stars and return.
c
      if ( ipos .le. 0 ) then
        do i = 1, 6
          chrtmp(i:i) = '*'
        end do
        return
      end if
c
c  Read the digits of INTVAL, from right to left.
c
      idig = mod ( ival, 10 )
      chrtmp(ipos:ipos) = char ( 48+idig )
      ipos = ipos - 1
      ival = ival / 10

      if ( ival .ne. 0 ) then
        go to 10
      end if
 
      if ( intval .lt. 0 ) then
 
        if ( ipos .le. 0 ) then
          do i = 1, 6
            chrtmp(i:i) = '*'
          end do
          return
        else
          chrtmp(ipos:ipos) = '-'
        end if
 
      end if
 
      call chrdb1 ( chrtmp )
      chlint = chrtmp
 
      return
      end
      function chlrat ( ival, jval )

c*********************************************************************72
c
cc CHLRAT returns a left-justified representation of IVAL/JVAL.
c
c
c  If the ratio is negative, a minus sign precedes IVAL.
c  A slash separates IVAL and JVAL.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer IVAL, JVAL, the two integers whose
c    ratio IVAL/JVAL is to be represented.
c
c    Note that if IVAL is nonzero and JVAL is 0, CHLRAT will
c    be returned as "Inf" or "-Inf" (Infinity), and if both
c    IVAL and JVAL are zero, CHLRAT will be returned as "NaN"
c    (Not-a-Number).
c
c    Output, character*22 CHLRAT, a left-justified string
c    containing the representation of IVAL/JVAL.
c
      character*22 chlrat
      character*22 chrtmp
      integer ival
      integer ival2
      integer jval
      integer jval2
c
c  Take care of simple cases right away.
c
      if ( ival .eq. 0 ) then
 
        if ( jval .ne. 0 ) then
          chrtmp = '0'
        else
          chrtmp = 'NaN'
        end if
 
      else if ( jval .eq. 0 ) then
 
        if ( ival .gt. 0 ) then
          chrtmp = 'Inf'
        else
          chrtmp = '-Inf'
        end if
c
c  Make copies of IVAL and JVAL.
c
      else
 
        ival2 = ival
        jval2 = jval
 
        if ( jval2 .eq. 1 ) then
          write ( chrtmp, '(i11)' ) ival2
        else
          write ( chrtmp, '(i11, ''/'', i10)' ) ival2, jval2
        end if
 
        call chrdb1 ( chrtmp )
 
      end if
 
      chlrat = chrtmp
 
      return
      end
      subroutine chrchp ( string, ilo, ihi )

c*********************************************************************72
c
cc CHRCHP "chops out" a portion of a string, and closes up the hole.
c
c
c  Using quotes to denote the beginning and end of the string, then
c  calling CHRCHP with STRING = 'Fred is not a jerk!' and ILO = 9 and
c  IHI = 12 will result in the output STRING = 'Fred is a jerk!    '
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character*(*) STRING, the character string
c    to be transformed.
c
c    Input, integer ILO, the location of the first character
c    to be removed.
c
c    Input, integer IHI, the location of the last character
c    to be removed.
c
      character*1 chrtmp
      integer i
      integer ihi
      integer ilo
      integer inew
      integer nchar
      character*(*) string

      nchar = len ( string )
 
      if ( ilo .gt. ihi ) then
        return
      end if
 
      do i = ilo, nchar
 
        inew = i-ilo+ihi+1
 
        if ( inew .le. nchar ) then
          chrtmp = string(inew:inew)
          string(i:i) = chrtmp
        else
          string(i:i) = ' '
        end if
 
      end do
 
      return
      end
      subroutine chrctf ( string, itop, ibot, ierror, lchar )

c*********************************************************************72
c
cc CHRCTF reads an integer or rational fraction from a string.
c
c
c  The integer may be in real format, for example '2.25'.  It
c  returns ITOP and IBOT.  If the input number is an integer, ITOP
c  equals that integer, and IBOT is 1.  But in the case of 2.25,
c  the program would return ITOP = 225, IBOT = 100.
c
c  Legal input is
c
c    blanks,
c    initial sign,
c    blanks,
c    integer part,
c    decimal point,
c    fraction part,
c    'E' or 'e' or 'D' or 'd', exponent marker,
c    exponent sign,
c    exponent integer part,
c    blanks,
c    final comma or semicolon,
c
c  with most quantities optional.
c
c  Example:
c
c  STRING            ITOP      IBOT
c
c  '1'               1         1
c  '     1   '       1         1
c  '1A'              1         1
c  '12,34,56'        12        1
c  '  34 7'          34        1
c  '-1E2ABCD'        -100      1
c  '-1X2ABCD'        -1        1
c  ' 2E-1'           2         10
c  '23.45'           2345      100
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) STRING, the string containing the
c    data to be read.  Reading will begin at position 1 and
c    terminate when no more characters
c    can be read to form a legal integer.  Blanks, commas,
c    or other nonnumeric data will, in particular, cause
c    the conversion to halt.
c
c    Output, integer ITOP, the integer read from the string,
c    assuming that no negative exponents or fractional parts
c    were used.  Otherwise, the 'integer' is ITOP/IBOT.
c
c    Output, integer IBOT, the integer divisor required to
c    represent numbers which are in real format or have a
c    negative exponent.
c
c    Output, integer IERROR, error flag.
c    0 if no errors,
c    Value of IHAVE when error occurred otherwise.
c
c    Output, integer LCHAR, number of characters read from
c    STRING to form the number.
c
      character*1 chrtmp
      integer ibot
      integer ierror
      integer ihave
      integer isgn
      integer iterm
      integer itop
      integer jsgn
      integer jtop
      integer lchar
      logical leqi
      integer nchar
      integer ndig
      character*(*) string

      nchar = len ( string )
 
      ierror = 0
      lchar = -1
      isgn = 1
      itop = 0
      ibot = 1
      jsgn = 1
      jtop = 0
      ihave = 1
      iterm = 0

10    continue

      lchar = lchar+1
      chrtmp = string(lchar+1:lchar+1)
c
c  Blank.
c
      if ( chrtmp .eq. ' ' ) then
 
        if ( ihave .eq. 2 ) then
 
        else if ( ihave .eq. 6 .or. ihave .eq. 7 ) then
          iterm = 1
        else if ( ihave .gt. 1 ) then
          ihave = 11
        end if
c
c  Comma.
c
      else if ( chrtmp .eq. ',' .or. chrtmp .eq. ';' ) then
 
        if ( ihave .ne. 1 ) then
          iterm = 1
          ihave = 12
          lchar = lchar+1
        end if
c
c  Minus sign.
c
      else if ( chrtmp .eq. '-' ) then
 
        if ( ihave .eq. 1 ) then
          ihave = 2
          isgn = -1
        else if ( ihave .eq. 6 ) then
          ihave = 7
          jsgn = -1
        else
          iterm = 1
        end if
c
c  Plus sign.
c
      else if ( chrtmp .eq. '+' ) then
 
        if ( ihave .eq. 1 ) then
          ihave = 2
        else if ( ihave .eq. 6 ) then
          ihave = 7
        else
          iterm = 1
        end if
c
c  Decimal point.
c
      else if ( chrtmp .eq. '.' ) then
 
        if ( ihave .lt. 4 ) then
          ihave = 4
        else
          iterm = 1
        end if
c
c  Exponent marker.
c
      else if ( leqi ( chrtmp, 'E' ) .or.
     &       leqi ( chrtmp, 'D' ) ) then
 
        if ( ihave .lt. 6 ) then
          ihave = 6
        else
          iterm = 1
        end if
c
c  Digit.
c
      else if ( 
     &  lge ( chrtmp, '0' ) .and.
     &  lle ( chrtmp, '9' ) .and.
     &  ihave .lt. 11 ) then
 
        if ( ihave .le. 2 ) then
          ihave = 3
        else if ( ihave .eq. 4 ) then
          ihave = 5
        else if ( ihave .eq. 6 .or. ihave .eq. 7 ) then
          ihave = 8
        end if
 
        call digten ( chrtmp, ndig )
 
        if ( ihave .eq. 3 ) then
          itop = 10 * itop + ndig
        else if ( ihave .eq. 5 ) then
          itop = 10 * itop + ndig
          ibot = 10*ibot
        else if ( ihave .eq. 8 ) then
          jtop = 10 * jtop + ndig
        end if
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
      if ( 
     &  (ihave .eq. 1) .or.
     &  (ihave .eq. 2) .or.
     &  (ihave .eq. 6) .or.
     &  (ihave .eq. 7) ) then
        ierror = ihave
        write ( *, * ) ' '
        write ( *, * ) 'CHRCTF - Serious error!'
        write ( *, * ) '  Illegal or nonnumeric input:'
        write ( *, '(a)' ) string
        return
      end if
c
c  Number seems OK.  Form it.
c
      if ( jsgn .eq. 1 ) then
        itop = itop * 10**jtop
      else
        ibot = ibot * 10**jtop
      end if
 
      itop = isgn * itop
 
      return
      end
      subroutine chrctg ( string, itop, ibot, ierror, lchar )

c*********************************************************************72
c
cc CHRCTG reads an integer, decimal fraction or a ratio from a string.
c
c
c  CHRCTG returns an equivalent ratio (ITOP/IBOT).
c
c  If the input number is an integer, ITOP equals that integer, and
c  IBOT is 1.   But in the case of 2.25, the program would return
c  ITOP = 225, IBOT = 100.
c
c  A ratio is either
c    a number
c  or
c    a number, "/", a number.
c
c  A "number" is defined as:
c
c    blanks,
c    initial sign,
c    integer part,
c    decimal point,
c    fraction part,
c    E,
c    exponent sign,
c    exponent integer part,
c    blanks,
c    final comma or semicolon,
c
c  with most quantities optional.
c
c  Examples of a number:
c
c    15, 15.0, -14E-7, E2, -12.73E-98, etc.
c
c  Examples of a ratio:
c
c    15, 1/7, -3/4.9, E2/-12.73
c
c  Sample results:
c
c  STRING            ITOP      IBOT
c
c  '1'               1         1
c  '     1   '       1         1
c  '1A'              1         1
c  '12,34,56'        12        1
c  '  34 7'          34        1
c  '-1E2ABCD'        -100      1
c  '-1X2ABCD'        -1        1
c  ' 2E-1'           2         10
c  '23.45'           2345      100
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) STRING, the string containing the
c    data to be read.  Reading will begin at position 1 and
c    terminate when no more characters
c    can be read to form a legal integer.  Blanks, commas,
c    or other nonnumeric data will, in particular, cause
c    the conversion to halt.
c
c    Output, integer ITOP, the integer read from the string,
c    assuming that no negative exponents or fractional parts
c    were used.  Otherwise, the 'integer' is ITOP/IBOT.
c
c    Output, integer IBOT, the integer divisor required to
c    represent numbers which are in decimal format or have a
c    negative exponent.
c
c    Output, integer IERROR, error flag.
c    0 if no errors,
c    Value of IHAVE in CHRCTF when error occurred otherwise.
c
c    Output, integer LCHAR, the number of characters read from
c    STRING to form the number.
c
      integer i
      integer ibot
      integer ibotb
      integer ierror
      integer igcf
      integer itemp
      integer itop
      integer itopb
      integer lchar
      integer lchar2
      integer lenchr
      integer nchar
      character*(*) string

      itop = 0
      ibot = 1
      lchar = 0
 
      call chrctf ( string, itop, ibot, ierror, lchar )

      if ( ierror .ne. 0) then
        return
      end if
c
c  The number is represented as a fraction.
c  If the next nonblank character is "/", then read another number.
c
      nchar = lenchr ( string )
 
      do i = lchar+1, nchar-1
 
        if ( string(i:i) .eq. '/' ) then
 
          call chrctf ( string(i+1:), itopb, ibotb, ierror, lchar2 )

          if ( ierror .ne. 0 ) then
            return
          end if
 
          itop = itop * ibotb
          ibot = ibot * itopb
 
          itemp = igcf ( itop, ibot )
 
          itop = itop / itemp
          ibot = ibot / itemp
 
          lchar = i + lchar2
 
          return
 
        else if ( string(i:i) .ne. ' ' ) then
 
          return
 
        end if
 
      end do
 
      return
      end
      subroutine chrcti ( string, intval, ierror, lchar )

c*********************************************************************72
c
cc CHRCTI reads an integer from a string.
c
c
c  CHRCTI will read as many characters as possible until it reaches
c  the end of the STRING, or encounters a character which cannot be
c  part of the number.
c
c  Legal input is
c
c    blanks,
c    initial sign,
c    blanks,
c    integer part,
c    blanks,
c    final comma or semicolon,
c
c  with most quantities optional.
c
c  Sample results:
c
c  STRING            INTVAL
c
c  '1'               1
c  '     1   '       1
c  '1A'              1
c  '12,34,56'        12
c  '  34 7'          34
c  '-1E2ABCD'        -100
c  '-1X2ABCD'        -1
c  ' 2E-1'           0
c  '23.45'           23
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
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
      subroutine chrctr ( string, rval, ierror, lchar )

c*********************************************************************72
c
cc CHRCTR reads a real number from a string.
c
c
c  CHRCTR will read as many characters as possible until it reaches
c  the end of the string, or encounters a character which cannot be
c  part of the real number.
c
c  Legal input is:
c
c     1 blanks,
c     2 '+' or '-' sign,
c     2.5 spaces
c     3 integer part,
c     4 decimal point,
c     5 fraction part,
c     6 'E' or 'e' or 'D' or 'd', exponent marker,
c     7 exponent sign,
c     8 exponent integer part,
c     9 exponent decimal point,
c    10 exponent fraction part,
c    11 blanks,
c    12 final comma or semicolon.
c
c  with most quantities optional.
c
c  Examples:
c
c    STRING            RVAL
c
c    '1'               1.0
c    '     1   '       1.0
c    '1A'              1.0
c    '12,34,56'        12.0
c    '  34 7'          34.0
c    '-1E2ABCD'        -100.0
c    '-1X2ABCD'        -1.0
c    ' 2E-1'           0.2
c    '23.45'           23.45
c    '-4.2E+2'         -420.0
c    '17d2'            1700.0
c    '-14e-2'         -0.14
c    'e2'              100.0
c    '-12.73e-9.23'   -12.73 * 10.0**(-9.23)
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
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
c    characters can be read to form a legal real.  Blanks,
c    commas, or other nonnumeric data will, in particular,
c    cause the conversion to halt.
c
c    Output, real RVAL, the real value that was read from the string.
c
c    Output, integer IERROR, error flag.
c
c    0, no errors occurred.
c
c    1, 2, 6 or 7, the input number was garbled.  The
c    value of IERROR is the last type of input successfully
c    read.  For instance, 1 means initial blanks, 2 means
c    a plus or minus sign, and so on.
c
c    Output, integer LCHAR, the number of characters read from
c    STRING to form the number, including any terminating
c    characters such as a trailing comma or blanks.
c
      character*1 chrtmp
      integer ierror
      integer ihave
      integer isgn
      integer iterm
      integer jbot
      integer jsgn
      integer jtop
      integer lchar
      logical leqi
      integer nchar
      integer ndig
      real rbot
      real rexp
      real rtop
      real rval
      character*(*) string

      nchar = len ( string )
 
      ierror = 0
      rval = 0.0
      lchar = -1
      isgn = 1
      rtop = 0.0
      rbot = 1.0
      jsgn = 1
      jtop = 0
      jbot = 1
      ihave = 1
      iterm = 0
 
10    continue

      lchar = lchar+1
      chrtmp = string(lchar+1:lchar+1)
c
c  Blank character.
c
      if ( chrtmp .eq. ' ' ) then
c
c  20 November 1993
c
c  I would like to allow input like "+ 2", where there is a space
c  between the plus and the number.  So I am going to comment out
c  this line, because I think that's all that's keeping me from
c  doing this.
c
c       if ( ihave .eq. 2 .or.
c    &     ihave .eq. 6 .or.
c    &     ihave .eq. 7 ) then
 
        if ( ihave .eq. 2 ) then
 
        else if ( ihave .eq. 6 .or. ihave .eq. 7 ) then
          iterm = 1
        else if ( ihave .gt. 1 ) then
          ihave = 11
        end if
c
c  Comma.
c
      else if ( chrtmp .eq. ',' .or. chrtmp .eq. ';' ) then
 
        if ( ihave .ne. 1 ) then
          iterm = 1
          ihave = 12
          lchar = lchar+1
        end if
c
c  Minus sign.
c
      else if ( chrtmp .eq. '-' ) then
 
        if ( ihave .eq. 1 ) then
          ihave = 2
          isgn = -1
        else if ( ihave .eq. 6 ) then
          ihave = 7
          jsgn = -1
        else
          iterm = 1
        end if
c
c  Plus sign.
c
      else if ( chrtmp .eq. '+' ) then
 
        if ( ihave .eq. 1 ) then
          ihave = 2
        else if ( ihave .eq. 6 ) then
          ihave = 7
        else
          iterm = 1
        end if
c
c  Decimal point.
c
      else if ( chrtmp .eq. '.' ) then
 
        if ( ihave .lt. 4 ) then
          ihave = 4
        else if ( ihave .ge. 6 .and. ihave .le. 8 ) then
          ihave = 9
        else
          iterm = 1
        end if
c
c  Exponent marker.
c
      else if ( leqi ( chrtmp, 'E' ) .or. 
     &  leqi ( chrtmp, 'D' ) ) then
 
        if ( ihave .lt. 6 ) then
          ihave = 6
        else
          iterm = 1
        end if
c
c  Digit.
c
      else if ( ihave .lt. 11 .and.
     &  lge ( chrtmp, '0' ) .and.
     &  lle ( chrtmp, '9' ) ) then
 
        if ( ihave .le. 2 ) then
          ihave = 3
        else if ( ihave .eq. 4 ) then
          ihave = 5
        else if ( ihave .eq. 6 .or. ihave .eq. 7 ) then
          ihave = 8
        else if ( ihave .eq. 9 ) then
          ihave = 10
        end if
 
        call digten ( chrtmp, ndig )
 
        if ( ihave .eq. 3 ) then
          rtop = 10.0 * rtop + real(ndig)
        else if ( ihave .eq. 5 ) then
          rtop = 10.0 * rtop + real(ndig)
          rbot = 10.0 * rbot
        else if ( ihave .eq. 8 ) then
          jtop = 10 * jtop + ndig
        else if ( ihave .eq. 10 ) then
          jtop = 10 * jtop + ndig
          jbot = 10 * jbot
        end if
c
c  Anything else is regarded as a terminator.
c
      else
        iterm = 1
      end if
c
c  If we haven't seen a terminator, and we haven't examined the
c  entire string, go get the next character.
c
      if ( iterm .ne. 1 .and. lchar+1 .lt. nchar ) then
        go to 10
      end if
c
c  If we haven't seen a terminator, and we have examined the
c  entire string, then we're done, and LCHAR is equal to NCHAR.
c
      if ( iterm .ne. 1 .and. lchar+1 .eq. nchar ) then
        lchar = nchar
      end if
c
c  Number seems to have terminated.  Have we got a legal number?
c  Not if we terminated in states 1, 2, 6 or 7!
c
      if ( 
     &  ihave .eq. 1 .or.
     &  ihave .eq. 2 .or.
     &  ihave .eq. 6 .or.
     &  ihave .eq. 7 ) then
 
        ierror = ihave
 
        return
      end if
c
c  Number seems OK.  Form it.
c
      if ( jtop .eq. 0 ) then
        rexp = 1.0
      else
 
        if ( jbot .eq. 1 ) then
          rexp = 10.0**( jsgn * jtop )
        else
          rexp = jsgn * jtop
          rexp = rexp / jbot
          rexp = 10.0**rexp
        end if
 
      end if
 
      rval = isgn * rexp * rtop / rbot
 
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
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character*(*) STRING, the string to be transformed.
c
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
      subroutine chrdb2 ( string )

c*********************************************************************72
c
cc CHRDB2 replaces consecutive blanks by one.
c
c
c  CHRDB2 left justifies the remainder and right pads with blanks.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character*(*) STRING, the string to be transformed.
c
      integer i
      integer j
      integer nchar
      character*1 newchr
      character*1 oldchr
      character*(*) string

      nchar = len ( string )
 
      j = 0
      newchr = ' '
 
      do i = 1, nchar
 
        oldchr = newchr
        newchr = string(i:i)
        string(i:i) = ' '
 
        if ( oldchr .ne. ' ' .or. newchr .ne. ' ' ) then
          j = j + 1
          string(j:j) = newchr
        end if
 
      end do
 
      return
      end
      subroutine chrinp ( ierror, iounit, line, nline, output, prompt )

c*********************************************************************72
c
cc CHRINP requests new input if the LINE buffer is empty.
c
c
c  CHRINP checks to see whether there is any more information in
c  the buffer array LINE.  If so, it simply updates the prompt
c  and returns.  Otherwise, it prints the prompt string out,
c  reads the input from the user, and reprints the prompt and
c  the user input on those I/O units where it is appropriate.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer IERROR.
c
c    If IERROR is nonzero on input, CHRINP stops.  It is the
c    calling routine's responsibility to make sure IERROR is
c    zero on input.  This is because CHRINP signals problems
c    to the calling routine using IERROR.  If the routine
c    does not take the trouble to reset IERROR, then it
c    is likely not to have addressed the problem itself.
c    These problems can include things like end of input,
c    so a failure to act can be catastrophic.
c
c    On output, IERROR = 
c      0 if no errors were detected,
c      1 if there was an error in the read,
c      2 if there was an end-of-file in the read.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Input/output, character*80 LINE.
c
c    On input, LINE may contain information that the calling
c    program can use, or LINE may be empty.
c
c    On output, LINE is unchanged if it contained information
c    on input.  But if the input LINE was empty, then the
c    output LINE contains whatever information the user typed.
c
c    Input/output, integer NLINE.
c
c    On input, if NLINE is zero, CHRINP assumes LINE is
c    empty, and asks the user for more information.
c
c    If NLINE is greater than zero, than NLINE and LINE
c    are left unchanged.
c
c    On output, NLINE is reset to the length of LINE.
c
c    Workspace, character*100 OUTPUT.
c
c    Input/output, character*80 PROMPT.
c    On input, the prompt string to be printed.
c    On output, PROMPT has been blanked out, up to the first comma.
c
      character*6 chrint
      integer i
      integer icomma
      integer ierror
      integer iosave
      integer iounit(4)
      integer lchar
      integer lenchr
      character*80 line
      integer nline
      character*100 output
      character*80 prompt
c
c  Catch nasty errors in calling routines.
c
      if ( ierror .ne. 0 ) then
        output = 'Error!'
        call chrwrt ( iounit, output )
        output = 'Nonzero input value of IERROR='//chrint(ierror)
        call chrdb2 ( output )
        call chrwrt ( iounit, output )
        stop
      end if
 
10    continue
c
c  If there is nothing in the LINE buffer, then:
c    "turn off" the automatic echo for units between 30 and 39,
c    print the prompt line,
c    "turn on" the automatic echo for units between 30 and 39,
c    read the input line,
c    remove double blanks,
c    set NLINE to the length of the LINE,
c    don't print a copy of the input on units between 40 and 49.
c
      if ( nline .le. 0 ) then
 
        do i = 2, 4
          if ( iounit(i) .ge. 30 .and. 
     &       iounit(i) .le. 39 ) then
            iounit(i)=-iounit(i)
          end if
        end do
 
        lchar = lenchr ( prompt )
        if ( lchar .gt. 0 ) then
          output = 'Enter ' // prompt(1:lchar)
          call chrwrt ( iounit, output )
        end if
 
        do i = 2, 4
          if ( iounit(i) .le. -30 .and. 
     &       iounit(i) .ge. -39) then
            iounit(i)=-iounit(i)
          end if
        end do
 
        if ( iounit(1) .le. 0 ) then
          read ( *, '(a80)', end = 50, err = 40 ) line
        else
          read ( iounit(1), '(a80)', end = 50, err = 40 ) line
        end if

        call chrdb2 ( line )
c
c  Don't echo input to IOUNIT(2).
c
        if ( iounit(1) .lt. 40 .or. iounit(1) .gt. 49 ) then
          iosave = iounit(2)
          if ( iounit(1) .le. 0 ) then
            iounit(2) = - 1
          end if
          output = line
          call chrwrt ( iounit, output )
          iounit(2) = iosave
        end if
 
      end if
c
c  Reset NLINE.
c
      nline = lenchr(line)
c
c  If the user typed something in, reset the line position to 0.
c
      if ( iounit(1) .eq. 0 ) then
        call setlin(0)
      end if
c
c  If item was read, remove item from PROMPT list.
c
      if ( nline .gt. 0 ) then

        icomma = index(prompt,',')

        if ( icomma .gt. 0 .and. 
     &    icomma .lt. 80 .and. 
     &    prompt(icomma+1:icomma+1) .eq. ' ' ) then

          icomma = icomma+1

        end if

        call chrchp ( prompt, 1, icomma )

      end if
 
      return
c
c  Error in input.
c
40    continue

      ierror = 1
      output = 'Error in input format.'
      call chrwrt ( iounit, output )
      output = 'Input line follows:'
      call chrwrt ( iounit, output )
      output = line
      call chrwrt ( iounit, output )
 
      if ( iounit(1) .le. 0 ) then
        nline = 0
        go to 10
      end if
 
      return
c
c  End of input.
c
c  If we are reading from a file, then set IERROR=2 and return.
c  But if we are reading from the user, something is seriously
c  wrong, and we must stop.
c
50    continue

      ierror = 2
      nline = 0
      output = 'End of input!'
      call chrwrt ( iounit, output )
 
      if ( iounit(1) .eq. 0 ) then
        output = 'The program is being stopped now!'
        call chrwrt ( iounit, output )
      else
        close ( unit = iounit(1) )
        iounit(1) = 0
        output = 'Closing current input file!'
        call chrwrt ( iounit, output )
      end if
 
      return
      end
      function chrint ( intval )

c*********************************************************************72
c
cc CHRINT returns a 6-character representation of an integer.
c
c
c  The representation is right justified, or '******' if
c  the integer is too large or negative to fit in six positions.
c
c  Compare CHLINT and CHR0NT.
c
c  Example:
c
c    INTVAL  STRING
c
c         1       1
c        -1      -1
c         0       0
c      1952    1952
c    123456  123456
c   1234567  ******  <-- Not enough room!
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer INTVAL, an integer variable to be converted.
c
c    Output (through function value), character*6 CHRINT, a 6
c    character representation of the integer, right justified.
c    Thus, if INTVAL = 1, CHRINT = '     1'.  CHRINT must be
c    declared "character*6 CHRINT" in the calling program.
c
      character*6 chrint
      character*6 chrtmp
      integer intval

      if ( intval .gt. 999999 ) then
        chrtmp = '******'
      else if ( intval .lt. -99999 ) then
        chrtmp = '-*****'
      else
        write ( chrtmp, '(i6)' ) intval
      end if
 
      chrint = chrtmp
 
      return
      end
      subroutine chritc0 ( string, intval )

c*********************************************************************72
c
cc CHRITC0 stores an integer into a string, with left zero padding.
c
c
c  The integer is left-justified in the string.  If the integer is
c  negative, the first character of the string is a minus sign.
c  All other entries of the string are filled with zero characters.
c  However, if the integer is too large to fit into the string, 
c  the entire string is filled with asterisks.
c
c
c  Example:
c
c  INTVAL   STRING
c  ------   ------
c
c       1   00001
c     -10   -0010
c   12345   12345
c  123456   *****
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character*(*) STRING, the string into which the
c    integer is to be stored.
c
c    Input, integer INTVAL, the integer to be stored in STRING.
c
      character*1 digit(0:9)
      integer i
      integer ihi
      integer ilo
      integer intval
      integer ival
      integer j
      character*(*) string

      save digit

      data digit / '0', '1', '2', '3', '4', '5,', '6', '7', '8', '9' /

      ihi = len ( string )
 
      if ( intval .lt. 0 ) then
        string(1:1) = '-'
        ilo = 2
      else
        ilo = 1
      end if
 
      ival = abs ( intval )
 
      do i = ihi, ilo, -1
        j = mod ( ival, 10 )
        string(i:i) = digit(j)
        ival = ival / 10 
      end do
 
      if ( ival .ne. 0 ) then
        do i = 1, ihi
          string(i:i) = '*'
        end do
      end if

      return
      end
      subroutine chrrea ( string, line, nline, prompt, iounit, ierror,
     &  iterm )

c*********************************************************************72
c
cc CHRREA extracts a character string from the input buffer.
c
c
c  CHRREA accepts LINE, which is assumed to contain NLINE user
c  input characters, where NLINE may be less than 1, and a PROMPT
c  line.
c
c  If NLINE is less than 1, the PROMPT is printed and user input
c  read from IOUNIT(1) into LINE, and NLINE updated.
c
c  In either case, enough characters are read from LINE to fill
c  STRING and the positions read are removed, and NLINE updated.
c
c  PROMPT is also updated.  On satisfactory input of STRING,
c  everything in PROMPT up to and including the first comma is
c  removed.
c
c  NOTE:
c
c  IOUNIT is assumed to have the following properties, which
c  also apply to routines CHRWRT, CHRINP, RELREA, RELWRT, INTREA
c  and RATREA:
c
c  IOUNIT(1) represents the input unit.  0 is taken to be the user
c  and we READ(*,format) the input.
c
c  IOUNIT(2) is taken to be a standard output unit.  Input is never
c  echoed to IOUNIT(2), but may be to other units.
c
c  Later units:  If their values is between 30 and 39, user input
c  is copied to them, but no output.
c  If between 40 and 49, output is copied to them, but no input.
c  If the unit number is negative, no input is read, nor output
c  written.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, character*(*) STRING.
c    The user's response to the PROMPT, as read from LINE.
c
c    Input/output, character*80 LINE.
c    A buffer containing the user's input.
c
c    Input/output, integer NLINE.
c    The number of characters of information in LINE.
c
c    Input/output, character*80 PROMPT.
c    On input, a prompt string that will be printed if NLINE is not positive.
c
c    On output, if STRING has been read, then PROMPT is cleared out
c    up to, and including, the first comma.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Output, integer IERROR, 
c    0, No error occurred.
c    1, Format error during read.
c    2, End of file during read.
c
c    Input, integer ITERM,
c    0, No check for terminators.
c    1, Blank, slash, comma, semicolon, equals, greater or
c       lesser signs terminate input.
c    2, Nonalphabetic terminates input
c    3, Nonalphanumeric terminates input
c    4, Blank, slash, comma, semicolon, equals, greater or
c       lesser signs or nonalphabetic characters terminate input.
c
      character*1 chrtmp
      integer i
      integer ierror
      integer iounit(4)
      integer iterm
      integer lchar
      integer lenchr
      logical let
      character*80 line
      integer nchar
      integer nline
      character*1 null
      logical num
      character*100 output
      character*80 prompt
      character*(*) string

      ierror = 0
      null = char(0)
      string = ' '
      call chrinp ( ierror, iounit, line, nline, output, prompt )
      if ( ierror .ne. 0 ) return
c
c  Remove double blanks.
c
      if ( iterm .eq. 2 .or. iterm .eq. 3 ) then
        call chrdb2 ( line )
      end if
c
c  Null input acceptable for character input only.
c
      if ( nline .le. 0 ) return
      lchar = 0
      nchar = len(string)
 
      do i = 1, nchar
 
        if ( lchar .ne. 0 ) go to 10
 
        chrtmp = line(i:i)
 
        if ( iterm .eq. 1 ) then
 
          if ( chrtmp .eq. ' ' .or. 
     &       chrtmp .eq. null .or. 
     &       chrtmp .eq. '/' .or. 
     &       chrtmp .eq. ',' .or. 
     &       chrtmp .eq. ';' .or. 
     &       chrtmp .eq. '=') then
            lchar = i
          end if
 
        else if ( iterm .eq. 2 ) then
 
          let = (lge(chrtmp,'a') .and. lle(chrtmp,'z')) .or. 
     &        (lge(chrtmp,'A') .and. lle(chrtmp,'Z'))

          if ( .not. let ) then
            lchar = i
          end if
 
        else if ( iterm .eq. 3 ) then
 
          let = (lge(chrtmp,'a') .and. lle(chrtmp,'z')) .or. 
     &        (lge(chrtmp,'A') .and. lle(chrtmp,'Z'))

          num=lge(chrtmp,'0') .and. lle(chrtmp,'9')

          if ( (.not.let) .and. (.not.num) ) then
            lchar = i
          end if
 
        else if ( iterm .eq. 4 ) then
 
          let = (lge(chrtmp,'a') .and. lle(chrtmp,'z')) .or. 
     &        (lge(chrtmp,'A') .and. lle(chrtmp,'Z'))

          if ( .not. let ) then
            lchar = i
          end if
 
          if ( chrtmp .eq. ' ' .or. 
     &       chrtmp .eq. null .or. 
     &       chrtmp .eq. '/' .or. 
     &       chrtmp .eq. ',' .or. 
     &       chrtmp .eq. ';' .or. 
     &       chrtmp .eq. '<' .or. 
     &       chrtmp .eq. '>' .or. 
     &       chrtmp .eq. '=' ) then
            lchar = i
          end if
 
        end if
 
        if ( lchar .eq. 0 ) then
          string(i:i) = chrtmp
        end if
 
      end do
 
10    continue
c
c  Chop out the character positions that have been used.
c
      if ( lchar .eq. 0 ) then
        lchar = nchar
      end if

      call chrchp ( line, 1, lchar )
c
c  Force the string to be flush left by removing leading blanks.
c
      call flushl ( line )
c
c  Update the line length.
c
      nline = lenchr ( line )
 
      return
      end
      function chrrel ( rval )

c*********************************************************************72
c
cc CHRREL represents a real using 14 right justified characters.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real RVAL, a real number.
c
c    Output (through function value), character*14 CHRREL,
c    a right-justified character variable containing the
c    representation of RVAL, using a G14.6 format.
c
      character*14 chrrel
      character*14 chrtmp
      real rval
c
c  We can't seem to write directly into CHRREL because of compiler
c  quibbles.
c
      if ( real ( int ( rval ) ) .eq. rval .and.
     &  abs ( rval ) .lt. 1.0e+13 ) then
 
        write ( chrtmp, '(i14)' ) int ( rval )
 
      else
 
        write ( chrtmp, '(g14.6)' ) rval
 
      end if
 
      chrrel = chrtmp
      return
      end
      subroutine chrup3 ( string, first, middle, last )

c*********************************************************************72
c
cc CHRUP3 divides a string into three parts, given the middle.
c
c
c  *  The FIRST part, up to the first occurrence of MIDDLE,
c  *  the MIDDLE part,
c  *  the LAST part.
c
c  For instance, if, on input,
c
c    STRING = 'aBCdEfgh'
c    MIDDLE = 'eF'
c
c  then the output will be:
c
c    FIRST = 'aBCd'
c    LAST =  'gh'
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) STRING, the string to be analyzed.
c
c    Output, character*(*) FIRST, the entries in STRING, up
c    to, but not including, the first occurrence, if any,
c    of MIDDLE.  If MIDDLE occurs immediately, then FIRST = ' '.
c    If FIRST is not long enough, trailing entries will be lost.
c    Unless FIRST is all blank, FIRST will not have leading blanks.
c
c    Input, character*(*) MIDDLE, the string to be searched
c    for in STRING.  Trailing blanks in MIDDLE are ignored,
c    unless MIDDLE is entirely blank, in which case it will
c    be treated as a single blank.  The check for the occurrence
c    of MIDDLE in STRING ignores case.
c
c    Output, character*(*) LAST, any characters of STRING which
c    occur after the first occurrence of MIDDLE.  If MIDDLE
c    occurs at the end of STRING, then LAST = ' '.  If LAST is not
c    long enough, then trailing entries will be lost.
c    Unless LAST is all blank, LAST will not have leading blanks.
c
      character*(*) first
      integer i
      integer indexi
      character*(*) last
      integer lenm
      integer lens
      character*(*) middle
      character*(*) string
c
c  If the divider is the blank, then remove multiple blanks.
c
      if ( middle .eq. ' ' ) then
        call chrdb2 ( string )
      end if

      lens = len ( string )
      lenm = len ( middle )
      i = indexi ( string, middle )
 
      if ( i .eq. 0 ) then
        first = string
        last = ' '
      else if ( i .eq. 1 ) then
        first = ' '
        last = string(lenm+1:)
      else if ( i + lenm .gt. lens ) then
        first = string
        last = ' '
      else
        first = string(1:i-1)
        last = string(i+lenm: )
      end if
c
c  Drop leading blanks from FIRST and LAST.
c
      call flushl ( first )
      call flushl ( last )
 
      return
      end
      subroutine chrwrt ( iounit, string )

c*********************************************************************72
c
cc CHRWRT writes a character STRING to one or more output units.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Input, character*(*) STRING, the string to be printed.
c
      character*1 cc
      integer i
      integer iounit(4)
      integer k
      integer lchar
      integer lenchr
      integer npage
      character*(*) string
c
c  If output is to the user, rather than to a file, then
c  see if we need to pause for a new page.
c
      if ( iounit(2) .eq. 0 .and. npage() .gt. 0 ) then
        write ( 6, * ) '(more)'
        read ( 5, '(a1)', end = 10, err = 10 )
10      continue
      end if
 
      lchar = lenchr(string)
      if ( lchar .le. 0 ) then
        lchar = 1
      end if
 
      do i = 2, 4

        if ( iounit(i) .eq. 0 ) then
c
c  Use the following line for UNIX machines, and the IBM PC:
c
c         write ( 6, '(80a1)' ) ( string(k:k), k = 1, lchar )
c
c  Use the following lines for Macintosh and VAX/VMS systems:
c
          cc=' '
          write ( 6, '(a1,80a1)' ) cc, ( string(k:k), k = 1, lchar )
 
        else if ( iounit(i) .gt. 0 ) then
          write( iounit(i), '(80a1)' ) ( string(k:k), k = 1, lchar )
        end if
 
      end do
c
c  Update the line count.
c
      call addlin ( )
 
      return
      end
      subroutine copmat ( a, b, iatop, iabot, ibtop, ibbot, ibase,
     &  ibaseb, lpmoda, lpmodb, maxcol, maxrow, nart, nartb, ncol,
     &  ncolb, nrow, nrowb, nslak, nslakb, nvar, nvarb )

c*********************************************************************72
c
cc COPMAT makes a copy of the current problem information.
c
c
c  A      --> B
c  IATOP  --> IBTOP
c  IABOT  --> IBBOT
c  IBASE  --> IBASEB
c  LPMODA --> LPMODB
c  NART   --> NARTB
c  NCOL   --> NCOLB
c  NROW   --> NROWB
c  NSLAK  --> NSLAKB
c  NVAR   --> NVARB
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real A(MAXROW,MAXCOL).  A is the current matrix.
c
c    Output, real B(MAXROW,MAXCOL), a copy of the input value of A.
c
c    Input, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
c    IATOP and IABOT represent the rational or decimal matrix.
c
c    Output, integer IBBOT(MAXROW,MAXCOL), IBBOT(MAXROW,MAXCOL),
c    a copy of the input values of IATOP and IABOT.
c
c    Input, integer IBASE(MAXROW).
c    Keeps track of the basic variables.
c
c    Output, integer IBASEB(MAXROW), a copy of the input value of IBASE.
c
c    Input, integer LPMODA.
c    0, the program is in linear algebra mode.
c    1, the program is in linear programming mode.
c
c    Output, integer LPMODB.
c    A copy of the input value of LPMODA.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input, integer NART, the number of artificial variables.
c
c    Output, integer NARTB, a copy of the input value of NART.
c
c    Input, integer NCOL, the number of columns in the matrix.
c
c    Output, integer NCOLB, a copy of the input value of NCOL.
c
c    Input, integer NROW, the number of rows in the matrix.
c
c    Output, integer NROWB, a copy of the input value of NROW.
c
c    Input, integer NSLAK, the number of slack variables.
c
c    Output, integer NSLAKB, a copy of the input value of NSLAK.
c
c    Input, integer NVAR, the number of basic variables.
c
c    Output, integer NVARB, a copy of the input value of NVAR.
c
      integer maxcol
      integer maxrow

      real a(maxrow,maxcol)
      real b(maxrow,maxcol)
      integer i
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer ibbot(maxrow,maxcol)
      integer ibtop(maxrow,maxcol)
      integer ibase(maxrow)
      integer ibaseb(maxrow)
      integer j
      integer lpmoda
      integer lpmodb
      integer nart
      integer nartb
      integer ncol
      integer ncolb
      integer nrow
      integer nrowb
      integer nslak
      integer nslakb
      integer nvar
      integer nvarb

      lpmodb = lpmoda
      nartb = nart
      ncolb = ncol
      nrowb = nrow
      nslakb = nslak
      nvarb = nvar
 
      do i = 1, maxrow
 
        ibaseb(i) = ibase(i)
 
        do j = 1, maxcol
 
          ibtop(i,j) = iatop(i,j)
          ibbot(i,j) = iabot(i,j)
          b(i,j) = a(i,j)
 
        end do
      end do
 
      return
      end
      subroutine dbldec ( dval, itop, ibot, ndig )

c*********************************************************************72
c
cc DBLDEC converts a double precision quantity to a decimal representation.
c
c
c  Given the double precision value DVAL, DBLDEC computes integers ITOP and 
c  IBOT so that it is approximately true that:
c
c    DVAL = ITOP * 10 ** IBOT
c
c  In particular, only NDIG digits of DVAL are used in constructing the 
c  representation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, double precision DVAL, the value whose decimal representation 
c    is desired.
c
c    Output, integer ITOP, IBOT, the approximate decimal representation of DVAL.
c    ITOP is an integer, strictly between -10**NDIG and 10**NDIG.
c    IBOT is an integer exponent of 10.
c
c    Input, integer NDIG, the number of decimal digits of DVAL
c    to be used in constructing the decimal representation.
c    Rounding is used.  NDIG should normally be 1 or greater.
c    Because of limited computer accuracy, NDIG should normally
c    be no more than 7.
c
      integer ibot
      integer itop
      integer ndig
      double precision dtop
      double precision dval
      real ten1
      real ten2
c
c  Special cases.
c
      if ( dval .eq. 0.0 ) then
        itop = 0
        ibot = 0
        return
      end if
c
c  Factor DVAL = DTOP * 10**IBOT
c
      dtop = dval
      ibot = 0
c
c  Now normalize so that 10**(NDIG-1) <= ABS(DTOP) < 10**(NDIG)
c
      ten1 = 10**(ndig-1)
      ten2 = 10**ndig
 
10    continue
 
      if ( abs(dtop) .lt. ten1 ) then
        dtop = dtop * 10.0
        ibot = ibot - 1
        go to 10
      else if ( abs(dtop) .ge. ten2 ) then
        dtop = dtop / 10.0
        ibot = ibot + 1
        go to 10
      end if
c
c  ITOP is the integer part of DTOP, rounded.
c
      itop = nint ( dtop )
c
c  Now divide out any factors of ten from ITOP.
c
20    continue
 
      if ( itop .ne. 0 ) then
        if ( 10 * ( itop / 10 ) .eq. itop ) then
          itop = itop / 10
          ibot = ibot + 1
          go to 20
        end if
      end if
 
      return
      end
      subroutine decadd ( ibot, ibot1, ibot2, itop, itop1, itop2, ndig )

c*********************************************************************72
c
cc DECADD adds two decimal quantities.
c
c
c  DECADD computes
c
c    ITOP * 10**IBOT = ITOP1 * 10**IBOT1 + ITOP2 * 10**IBOT2
c
c  while trying to avoid integer overflow.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IBOT, the exponent of the result.
c
c    Input, integer IBOT1, IBOT2, the exponents of the two
c    numbers to be added.
c
c    Output, integer ITOP, the coefficient of the result.
c
c    Input, integer ITOP1, ITOP2, the coefficients of the two
c    numbers to be added.
c
c    Input, integer NDIG, the number of decimal digits used.
c
      integer ibot
      integer ibot1
      integer ibot2
      integer itop
      integer itop1
      integer itop2
      integer jtop1
      integer jtop2
      integer ndig

      if ( itop1 .eq. 0 ) then
        itop = itop2
        ibot = ibot2
        return
      else if ( itop2 .eq. 0 ) then
        itop = itop1
        ibot = ibot1
        return
      else if ( ibot1 .eq. ibot2 ) then
        itop = itop1 + itop2
        ibot = ibot1
        call deccut ( itop, ibot, ndig )
        return
      end if
c
c  Line up the exponents.
c
      jtop1 = itop1
      jtop2 = itop2
 
      if ( ibot1 .lt. ibot2 ) then
        jtop2 = jtop2 * 10**(ibot2-ibot1)
      else
        jtop1 = jtop1 * 10**(ibot1-ibot2)
      end if
c
c  Add the coefficients.
c
      itop = jtop1 + jtop2
      ibot = min ( ibot1, ibot2 )
c
c  Clean up the result.
c
      call deccut ( itop, ibot, ndig )
 
      return
      end
      subroutine deccut ( itop, ibot, ndig )

c*********************************************************************72
c
cc DECCUT rounds a decimal fraction to a given number of digits.
c
c
c  DECCUT takes an arbitrary decimal fraction represented by
c
c    ITOP * 10**IBOT
c
c  and makes sure that ITOP has no more than NDIG digits.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer ITOP, IBOT, the coefficient and exponent
c    of a decimal fraction.  On return, ITOP has no more than 
c    NDIG decimal digits.
c
c    Input, integer NDIG, the number of decimal digits used.
c
      integer ibot
      integer itop
      integer ndig

      if ( itop .eq. 0 ) then
        ibot = 0
        return
      end if
 
10    continue
 
      if ( abs ( itop ) .ge. 10**ndig ) then
        itop = itop / 10
        ibot = ibot + 1
        go to 10
      end if
 
20    continue
 
      if ( ( itop / 10 ) * 10 .eq. itop ) then
        itop = itop / 10
        ibot = ibot + 1
        go to 20
      end if
 
      return
      end
      subroutine decdet ( iarray, iatop, iabot, idtop, idbot, ierror,
     &  lda, maxint, n, ndig )

c*********************************************************************72
c
cc DECDET finds the determinant of an N by N matrix of decimal entries.
c
c
c  The brute force method is used.
c
c  DECDET should only be used for small matrices, since this calculation
c  requires the summation of N! products of N numbers.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Workspace, integer IARRAY(N).
c
c    Input, integer IATOP(LDA,N), IABOT(LDA,N), the decimal
c    representation of the matrix.
c
c    Output, integer IDTOP, IDBOT, the decimal determinant of the matrix.
c
c    Output, integer IERROR.
c    0, no error occurred.
c    1, an error occurred (probably overflow).
c
c    Input, integer LDA, the leading dimension of A.
c
c    Input, integer MAXINT, the maximum integer representable, which
c    should probably be 2147483647.
c
c    Input, integer N, the number of rows and columns of A.
c
c    Input, integer NDIG, the number of decimal digits used.
c
      integer lda
      integer n

      logical even
      integer i
      integer iabot(lda,n)
      integer iatop(lda,n)
      integer iarray(n)
      integer ibot
      integer ibot1
      integer ibot2
      integer idbot
      integer idtop
      integer ierror
      integer itop
      integer itop1
      integer itop2
      integer maxint
      logical more
      integer ndig

      ierror = 0
      more = .false.
      idtop = 0
      idbot = 1
c
c  Compute the next permutation.
c
10    continue
 
      call pernex ( n, iarray, more, even )
c
c  The sign of this term depends on the sign of the permutation.
c
      if ( even ) then
        itop = 1
      else
        itop = -1
      end if
c
c  Choose one item from each row, as specified by the permuation,
c  and multiply them
c
      ibot = 0
 
      do i = 1, n
 
        itop1 = itop
        ibot1 = ibot
        itop2 = iatop(i,iarray(i))
        ibot2 = iabot(i,iarray(i))

        call decmul ( ibot, ibot1, ibot2, itop, itop1, itop2, 
     &    maxint, ndig )
 
      end do
c
c  Add this term to the total.
c
      itop1 = itop
      ibot1 = ibot
 
      itop2 = idtop
      ibot2 = idbot
 
      call decadd ( ibot, ibot1, ibot2, itop, itop1, itop2, ndig )
 
      idtop = itop
      idbot = ibot
 
      if ( more ) then
        go to 10
      end if
 
      return
      end
      subroutine decdiv ( ibot, ibot1, ibot2, ierror, itop, itop1,
     &  itop2, ndig )

c*********************************************************************72
c
cc DECDIV divides two NDIG digit decimals.
c
c
c  DECDIV computes 
c
c    ITOP * 10**IBOT = (ITOP1 * 10**IBOT1) / (ITOP2 * 10**IBOT2)
c
c                    = (ITOP1/ITOP2) * 10**(IBOT1-IBOT2)
c
c  while avoiding integer overflow.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IBOT, the exponent of the result.
c
c    Input, integer IBOT1, the exponent of the dividend.
c
c    Input, integer IBOT2, the exponent of the divisor.
c
c    Output, integer IERROR.
c    0, no error occurred.
c    1, an error occurred.
c
c    Output, integer ITOP, the coefficient of the result.
c
c    Input, integer ITOP1, the coefficient of the dividend.
c
c    Input, integer ITOP2, the coefficient of the divisor.
c
c    Input, integer NDIG, the number of decimal digits used.
c
      double precision dval
      integer ibot
      integer ibot1
      integer ibot2
      integer ibot3
      integer ierror
      integer itop
      integer itop1
      integer itop2
      integer itop3
      integer ndig
c
c  First special case, top fraction is 0.
c
      if ( itop1 .eq. 0 ) then
        itop = 0
        ibot = 0
        return
      end if
c
c  First error, bottom of fraction is 0.
c
      if ( itop2 .eq. 0 ) then
        ierror = 1
        itop = 0
        ibot = 0
        return
      end if
c
c  Second special case, result is 1.
c
      if ( itop1 .eq. itop2 .and. ibot1 .eq. ibot2 ) then
        itop = 1
        ibot = 0
        return
      end if
c
c  Third special case, result is power of 10.
c
      if ( itop1 .eq. itop2 ) then
        itop = 1
        ibot = ibot1 - ibot2
        return
      end if
c
c  Fourth special case: ITOP1/ITOP2 is exact.
c
      if ( (itop1/itop2) * itop2 .eq. itop1 ) then
        itop = itop1 / itop2
        ibot = ibot1 - ibot2
        return
      end if
c
c  General case.
c
      dval = dble ( itop1 ) / dble ( itop2 )
 
      call dbldec ( dval, itop3, ibot3, ndig )
 
      itop = itop3
      ibot = ibot3 + ibot1 - ibot2
 
      return
      end
      subroutine decmul ( ibot, ibot1, ibot2, itop, itop1, itop2, 
     &  maxint, ndig )

c*********************************************************************72
c
cc DECMUL multiplies two decimals.
c
c
c  DECMUL computes
c
c    ITOP * 10**IBOT = (ITOP1 * 10**IBOT1) * (ITOP2 * 10**IBOT2)
c                    = (ITOP1*ITOP2) * 10**(IBOT1+IBOT2)
c
c  while avoiding integer overflow.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IBOT, the exponent of the result.
c
c    Input, integer IBOT1, the exponent of the first factor.
c
c    Input, integer IBOT2, the exponent of the second factor.
c
c    Output, integer ITOP, the coefficient of the result.
c
c    Input, integer ITOP1, the coefficient of the first factor.
c
c    Input, integer ITOP2, the coefficient of the second factor.
c
c    Input, integer MAXINT, the maximum integer representable, which
c    should probably be 2147483647.
c
c    Input, integer NDIG, the number of decimal digits used.
c
      double precision dval
      integer ibot
      integer ibot1
      integer ibot2
      integer ibot3
      integer itop
      integer itop1
      integer itop2
      integer itop3
      integer maxint
      integer ndig
      real temp
c
c  The result is zero if either ITOP1 or ITOP2 is zero.
c
      if ( itop1 .eq. 0 .or. itop2 .eq. 0 ) then
        itop = 0
        ibot = 0
        return
      end if
c
c  The result is simple if either ITOP1 or ITOP2 is one.
c
      if ( itop1 .eq. 1 .or. itop2 .eq. 1 ) then
        itop = itop1 * itop2
        ibot = ibot1 + ibot2
        return
      end if
 
      temp =
     &  log ( real ( abs ( itop1 ) ) ) +
     &  log ( real ( abs ( itop2 ) ) )
 
      if ( temp .lt. log ( real ( maxint ) ) ) then
 
        itop = itop1 * itop2
        ibot = ibot + ibot2
 
      else
 
        dval = dble(itop1) * dble(itop2)
 
        call dbldec ( dval, itop3, ibot3, ndig )
 
        itop = itop3
        ibot = ibot3 + ( ibot1 + ibot2 )
 
      end if
 
      return
      end
      subroutine decprn ( iatop, iabot, ibase, iounit, ihi, ilo, jhi,
     &  jlo, lpmoda, maxcol, maxrow, ncol, nrow, output, title )

c*********************************************************************72
c
cc DECPRN prints out decimal vectors and matrices.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
c    IATOP and IABOT represent the decimal matrix.
c
c    Input, integer IBASE(MAXROW), keeps track of basic variables.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Input, integer IHI, ILO, the last and first rows to print.
c
c    Input, integer JHI, JLO, the last and first columns to print.
c
c    Input, integer LPMODA.
c    0, the program is in linear algebra mode.
c    1, the program is in linear programming mode.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input, integer NCOL, the number of columns in the matrix.
c
c    Input, integer NROW, the number of rows in the matrix.
c
c    Workspace, character*100 OUTPUT.
c
c    Input, character*(*) TITLE, a label for the object being printed.
c
      integer ncolum
      parameter ( ncolum = 80 )

      integer maxcol
      integer maxrow

      character*22 chldec
      character*6 chrint
      character*22 chrtmp
      character*40 fortwo
      integer i
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer ibase(maxrow)
      integer ichi
      integer iclo
      integer ihi
      integer ilo
      integer imax
      integer imin
      integer iounit(4)
      integer izhi
      integer izlo
      integer j
      integer jhi
      integer jlo
      integer jmax
      integer jmin
      integer khi
      integer klo
      integer kmax
      character*4 lab
      integer lenc
      integer lenchr
      integer llab
      integer lpmoda
      integer ncol
      integer npline
      integer nrow
      character*100 output
      character*(*) title

      if ( lpmoda .eq. 1 ) then
        llab = 4
      else
        llab = 0
      end if
c
c  Figure out how wide we must make each column.
c
      imax = 0
      jmax = 0
 
      do i = ilo, ihi
        do j = jlo, jhi
 
          chrtmp = chldec ( iatop(i,j), iabot(i,j) )
          lenc = lenchr ( chrtmp )
 
          jmax = max ( jmax, lenc )
 
        end do
      end do
 
      kmax = 2 + imax + 1 + jmax
      npline = ( ncolum - llab ) / kmax
c
c  Set up the format for the heading.
c
      if ( lpmoda .eq. 1 ) then
        fortwo = '(' // chrint(llab) // 'x,' // chrint(npline) // 'i'
     &    // chrint(kmax) // ')'
      else
        fortwo = '(' // chrint(npline) // 'i' // chrint(kmax) // ')'
      end if
 
      call chrdb1 ( fortwo )
 
      do jmin = jlo, jhi, npline
 
        jmax = min ( jmin+npline-1, jhi )
 
        lab = '    '
c
c  Handle a column vector.
c
        if ( jlo .eq. jhi .and. ilo .ne. ihi ) then
 
          output = ' '
          call chrwrt ( iounit, output )
 
          if ( ilo .eq. 1 ) then
            output = title
            call chrwrt ( iounit, output )
            output = 'Column ' // chrint(jlo) // ' transposed.'
            call chrdb2 ( output )
            call chrwrt ( iounit, output )
          end if
 
          do imin = ilo, ihi, npline
 
            imax = min ( imin+npline-1, ihi )
 
            output = ' '
            call chrwrt ( iounit, output )
 
            do i = imin, imax
              ilo = 4 + ( i - imin ) * kmax + 1
              ihi = 4 + ( i - imin ) * kmax + kmax
              chrtmp = chldec ( iatop(i,jlo), iabot(i,jlo) )
              call flushr ( chrtmp(1:kmax) )
              output(ilo:ihi) = chrtmp(1:kmax)
            end do
 
            call chrwrt ( iounit, output )
 
          end do
 
          go to 90
        end if
 
        output = ' '
        call chrwrt ( iounit, output )
 
        if ( jmin .eq. 1 ) then
          output = title
          call chrwrt ( iounit, output )
          output = ' '
          call chrwrt ( iounit, output )
        end if
c
c  Print heading for linear programming tableau.
c
        if ( lpmoda .eq. 1 ) then
 
          write ( output, fortwo ) ( j, j = jmin, jmax )
 
          if ( jmin .le. ncol-1 .and. ncol-1 .le. jmax ) then
            izlo = llab + ((ncol-1)-jmin) * kmax + kmax - 2
            izhi = izlo + 2
            output(izlo:izhi) = '  P'
          end if
 
          if ( jmin .le. ncol .and. ncol .le. jmax ) then
            iclo = llab + (ncol-jmin) * kmax + kmax - 2
            ichi = iclo+2
            output(iclo:ichi) = '  C'
          end if
 
          call chrwrt ( iounit, output )
 
          output = ' '
          call chrwrt ( iounit, output )
c
c  Print heading for linear algebra matrix.
c
        else
 
          if ( jmin .gt. 1 .or. jmax .lt. ncol .or.
     &       ilo .gt. 1 .or. ihi .lt. nrow ) then

            output = 'Columns ' // chrint(jmin) // ' to ' // 
     &        chrint(jmax)
            call chrdb2 ( output )
            call chrwrt ( iounit, output )
            output = ' '
            call chrwrt ( iounit, output )
          end if
 
        end if
 
        do i = ilo, ihi
 
          if ( lpmoda .eq. 1 ) then

            if ( i .lt. nrow ) then
              if ( ibase(i) .lt. 10 ) then
                write ( lab, '(a1,i1)' ) 'X', ibase(i)
              else
                write ( lab, '(a1,i2)' ) 'X', ibase(i)
              end if
            else if ( i .lt. ihi ) then
              lab = 'Obj2'
            else
              lab = 'Obj '
            end if

            if ( maxrow .eq. 1 ) then
              lab = '    '
            end if

          end if
 
          if ( lpmoda .eq. 1 ) then
            output(1:4) = lab
          else
            output(1:4) = '    '
          end if
 
          do j = jmin, jmax
            klo = 4 + (j-jmin) * kmax + 1
            khi = 4 + (j-jmin) * kmax + kmax
            chrtmp = chldec ( iatop(i,j), iabot(i,j) )
            call flushr ( chrtmp(1:kmax) )
            output(klo:khi) = chrtmp(1:kmax)
          end do
 
          call chrwrt ( iounit, output )
 
        end do
 
90      continue
 
      end do
 
      return
      end
      subroutine decrat ( iatop, iabot )

c*********************************************************************72
c
cc DECRAT converts a decimal to a rational representation.
c
c
c  On input, a value is represented as IATOP * 10**IABOT.
c
c  On output, approximately the same value is represented as IATOP / IABOT.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer IATOP, IABOT.
c    On input, these quantities represent the value IATOP * 10 ** IABOT.
c    On output, these quantities represent the value IATOP / IABOT.
c
      integer iabot
      integer iatop
      integer igcf
      integer itmp

      external igcf

      if ( iabot .ge. 0 ) then
        iatop = iatop * 10**iabot
        iabot = 1
      else
        iabot = 10**(-iabot)
        itmp = igcf ( iatop, iabot )
        iatop = iatop / itmp
        iabot = iabot / itmp
      end if
 
      return
      end
      subroutine decrea ( itop, ibot, rval, line, maxdig, nline,
     &  prompt, iounit, ierror )

c*********************************************************************72
c
cc DECREA reads a decimal, rational or integer, and returns a decimal fraction.
c
c
c  Right now, this routine uses a lousy method, converting an arbitrary
c  fraction to a real value, then taking a MAXDIG digit decimal
c  expansion of that.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer ITOP, IBOT, represents the decimal fraction.
c
c    Output, real RVAL, the real value equivalent to the decimal fraction.
c
c    Workspace, character*80 LINE.
c    Used to hold the user's input.
c
c    Input, integer MAXDIG, the maximum number of decimals to use.
c
c    Input/output, integer NLINE.
c    Keeps track of the number of useful characters in LINE.
c
c    Input/output, character*80 PROMPT, the prompt string.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Output, integer IERROR.
c    0, no error occurred.
c    1, an error occurred.
c
      integer ibot
      integer ibot1
      integer ibot2
      integer ierror
      integer iounit(4)
      integer itop
      integer itop1
      integer itop2
      integer lchar
      integer lenchr
      character*80 line
      integer llchar
      integer maxdig
      integer nline
      character*100 output
      character*80 prompt
      real rval
c
      itop = 0
      ibot = 1
      rval = 0
      llchar = len(line)
 
10    continue
 
      call chrinp ( ierror, iounit, line, nline, output, prompt )
      if ( ierror .ne. 0 ) return
 
      if ( nline .le. 0 ) go to 10
 
      call chrctf ( line, itop1, ibot1, ierror, lchar )
 
      if ( lchar .ge. llchar ) then
        itop = itop1
        ibot = ibot1
      else if ( line(lchar+1:lchar+1) .ne. '/' ) then
        itop = itop1
        ibot = ibot1
      else
        lchar = lchar+1
        call chrchp ( line, 1, lchar )
        call chrctf ( line, itop2, ibot2, ierror, lchar )
        itop = itop1 * ibot2
        ibot = ibot1 * itop2
      end if
 
      rval = itop
      rval = rval / ibot
 
      call reldec ( rval, itop, ibot, maxdig )
 
      call chrchp ( line, 1, lchar )
      nline = lenchr ( line )
 
      return
      end
      subroutine decrel ( a, itop, ibot )
c
c*********************************************************************72
c
cc DECREL converts a decimal ITOP * 10**IBOT to a real value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, real A, the equivalent real value.
c
c    Input, integer ITOP, IBOT, the coefficient and exponent
c    of the decimal value.
c
      real a
      integer ibot
      integer itop

      a = itop * 10.0**ibot
 
      return
      end
      subroutine decwrn ( iounit, output )

c*********************************************************************72
c
cc DECWRN prints out a warning about using decimal arithmetic.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Workspace, character*100 OUTPUT.
c
      integer iounit(4)
      character*100 output
      logical said

      save said

      data said /.false./

      if ( .not. said ) then
 
        output = ' '
        call chrwrt ( iounit, output )
        output = 'Note:  The representation of decimals is exact.'
        call chrwrt ( iounit, output )
        output = ' '
        call chrwrt ( iounit, output )
        output = 'However, this representation will break down'
        call chrwrt ( iounit, output )
        output = 'if any exponent becomes too large or small.'
        call chrwrt ( iounit, output )
        output = ' '
        call chrwrt ( iounit, output )
        output = 'Calculations with decimals are NOT exact.'
        call chrwrt ( iounit, output )
 
        said = .true.
 
      end if
 
      return
      end
      subroutine deladd ( a, iabot, iatop, ibase, ierror, iform, 
     &  imat, iounit, line, lpmoda, maxcol, maxdig, maxrow, ncol,
     &  ncon, ndig, nline, nrow, nslak, nvar, output, prompt )

c*********************************************************************72
c
cc DELADD deletes or adds a row or column to the matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, real A(MAXROW,MAXCOL).  A is the matrix
c    to be changed.
c
c    Input/output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
c    IATOP and IABOT represent the rational or decimal matrix
c    to be changed.
c
c    Input/output, integer IBASE(MAXROW), keeps track of the basic variables.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IFORM, specifies the arithmetic being used.
c    0=rational, 1=real, 2=decimal.
c
c    Input, integer IMAT.
c    0, no matrix has been defined by the user.
c    1, a matrix has been defined by the user.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Workspace, character*80 LINE.
c    Used to hold the user's input.
c
c    Input, integer LPMODA.
c    0, the program is in linear algebra mode.
c    1, the program is in linear programming mode.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXDIG, the maximum number of decimals to use.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input/output, integer NCOL, the number of columns in the matrix.
c
c    Input/output, integer NCON, the number of constraints.
c
c    Input, integer NDIG, the number of decimal digits used.
c
c    Input/output, integer NLINE.
c    Keeps track of the number of useful characters in LINE.
c
c    Input/output, integer NROW, the number of rows in the matrix.
c
c    Input/output, integer NSLAK, the number of slack variables.
c
c    Input, integer NVAR, the number of basic variables.
c
c    Workspace, character*100 OUTPUT.
c
c    Workspace, character*80 PROMPT.
c
      integer maxcol
      integer maxrow
c
      real a(maxrow,maxcol)
      character*6 chrint
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer ibase(maxrow)
      integer icol
      integer ierror
      integer iform
      integer imat
      integer iounit(4)
      integer irow
      character*1 isay
      integer iterm
      logical leqi
      character*80 line
      integer lpmoda
      integer maxdig
      integer ncol
      integer ncon
      integer ndig
      integer nline
      integer nrow
      integer nslak
      integer nvar
      character*100 output
      character*80 prompt
c
      ierror = 0

      if ( imat .eq. 0 ) then
        ierror = 1
        output = 'You must set up a matrix first!'
        call chrwrt ( iounit, output )
        return
      end if
 
      if ( lpmoda .eq. 0 ) then
 
        prompt = '"+" to add a row or column, "-" to delete one.'
        iterm = 0
        call chrrea ( isay, line, nline, prompt, iounit, ierror,
     &    iterm )
 
        if ( isay .eq. '-' ) then
 
          prompt = 'R or C to delete a row or column.'

          iterm = 0
          call chrrea ( isay, line, nline, prompt, iounit, ierror,
     &      iterm )
c
c  -R: Delete a row in linear algebra mode.
c  --  -----------------------------------
c
          if ( leqi ( isay, 'R' ) ) then
c
c  Get row index.
c
            prompt = 'row to delete, between 1 and ' // chrint(nrow)
            call chrdb2 ( prompt )
            call intrea ( irow, line, nline, prompt, iounit, ierror )
            if ( ierror .ne. 0 ) return
 
            if ( irow .lt. 1 .or. irow .gt. nrow ) then
              ierror = 1
              output = 'Your row index was not acceptable!'
              return
            end if
c
c  Shift matrix rows.
c
            call delrow ( a, iabot, iatop, irow, maxcol, maxrow, ncol,
     &        nrow )
 
            nrow = nrow - 1
 
            output = 'The row has been deleted!'
            call chrwrt ( iounit, output )
c
c  -C:  Delete a column in linear algebra mode.
c  --   --------------------------------------
c
          else if ( leqi ( isay, 'C' ) ) then
c
c  Get column index.
c
            prompt = 'column to delete, between 1 and ' //
     &        chrint(ncol)
            call chrdb2 ( prompt )
            call intrea ( icol, line, nline, prompt, iounit, ierror )
            if ( ierror .ne. 0 ) return
 
            if ( icol .lt. 1 .or. icol .gt. ncol ) then
              ierror = 1
              output = 'Your column index was not acceptable!'
              return
            end if
c
c  Shift matrix columns.
c
            call delcol ( a, iabot, iatop, icol, maxcol, maxrow, ncol,
     &        nrow )
 
            ncol = ncol - 1
 
            output = 'The column has been deleted!'
            call chrwrt ( iounit, output )
 
          else
 
            ierror = 1
 
          end if
 
        else if ( isay .eq. '+' ) then
 
          prompt = 'R or C to add a row or column.'
          iterm = 0
          call chrrea ( isay, line, nline, prompt, iounit, ierror,
     &      iterm )
c
c  +R:  Add a row in linear programming mode.
c  --   ------------------------------------
c
          if ( leqi ( isay, 'R' ) ) then
 
            if ( nrow .lt. maxrow ) then
 
              nrow=nrow+1
c
c  Get row index.
c
              prompt = 'index for new row between 1 and '
     &          // chrint(nrow)
              call chrdb2 ( prompt )
              call intrea ( irow, line, nline, prompt, iounit, ierror )
              if ( ierror .ne. 0 ) return
 
              if ( irow .lt. 1 .or. irow .gt. nrow ) then
                ierror = 1
                output = 'Your row index was not acceptable!'
                return
              end if
c
c  Shift matrix rows.
c
              call shfrow ( a, iabot, iatop, irow, maxcol, maxrow, ncol,
     &          nrow )
c
c  Read in values for new row.
c
              icol = 0
 
              call lainp1 ( a, iabot, iatop, icol, ierror, iform, 
     &          iounit, irow, line, maxcol, maxdig, maxrow, ncol, ndig,  
     &          nline, nrow, output, prompt )
 
            else
              ierror = 1
              output = 'There is no space for more rows!'
              call chrwrt ( iounit, output )
            end if
          end if
c
c  +C: Add a column in linear programming mode.
c  --  ---------------------------------------
c 
          if ( leqi ( isay, 'C' ) ) then
 
            if ( ncol .lt. maxcol ) then
 
              ncol=ncol+1
c
c  Get column index.
c
              prompt = 'index for new column between 1 and '//
     &          chrint(ncol)
              call chrdb2 ( prompt )
              call intrea ( icol, line, nline, prompt, iounit, ierror )
              if ( ierror .ne. 0 ) return
 
              if ( icol .lt. 1 .or. icol .gt. ncol ) then
                ierror = 1
                output = 'Your column index was not acceptable!'
                return
              end if
c
c  Shift matrix columns.
c
              call shfcol ( a, iabot, iatop, icol, maxcol, maxrow,
     &          ncol, nrow )
c
c  Read in values for new column.
c
              irow = 0
 
              call lainp1 ( a, iabot, iatop, icol, ierror, iform,
     &          iounit, irow, line, maxcol, maxdig, maxrow, ncol, ndig,
     &          nline, nrow, output, prompt )
 
            else
              output = 'Error!  There is no space for more rows!'
              call chrwrt ( iounit, output )
            end if
 
          end if
        end if
c
c  Add new constraint and slack variable for linear programming.
c
      else
 
        if ( nrow .ge. maxrow-2 ) then
          ierror = 1
          output = 'Error!'
          call chrwrt ( iounit, output )
          output = 'The tableau cannot be increased in size to'
          call chrwrt ( iounit, output )
          output = 'make room for the new constraint!'
          call chrwrt ( iounit, output )
          return
        end if
 
        if ( nvar .ge. maxcol ) then
          ierror = 1
          output = 'The tableau cannot be increased in size to'
          call chrwrt ( iounit, output )
          output = 'make room for the new slack variable!'
          call chrwrt ( iounit, output )
          return
        end if
 
        output = 'Add a new constraint and slack variable!'
        call chrwrt ( iounit, output )
c
c  Shift last row down, shift last column to right.
c
c
c  DOES THIS DEPEND ON WHETHER ARTIFICIAL VARIABLES ARE INVOLVED?
c
        ncon = ncon + 1
        irow = ncon
        nrow = nrow + 1
        call shfrow ( a, iabot, iatop, irow, maxcol, maxrow, ncol,
     &    nrow )
        nslak = nslak + 1
        icol = nvar + nslak
        ncol = ncol + 1
        call shfcol ( a, iabot, iatop, icol, maxcol, maxrow, ncol,
     &    nrow )
        ibase(irow) = nvar + nslak
c
c  Read in values of constraint.
c
      end if
 
      return
      end
      subroutine delcol ( a, iabot, iatop, icol, maxcol, maxrow, ncol,
     &  nrow )

c*********************************************************************72
c
cc DELCOL deletes a column by shifting other columns to the left.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, real A(MAXROW,MAXCOL).  A is the matrix to be changed.
c
c    Output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
c    IATOP and IABOT represent the rational or decimal matrix
c    to be changed.
c
c    Input, integer ICOL, the column to be deleted.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input, integer NCOL, the number of columns in the matrix.
c
c    Input, integer NROW, the number of rows in the matrix.
c
      integer maxcol
      integer maxrow

      real a(maxrow,maxcol)
      integer i
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer icol
      integer j
      integer ncol
      integer nrow

      do j = icol, ncol-1
        do i = 1, nrow
 
          a(i,j) = a(i,j+1)
          iatop(i,j) = iatop(i,j+1)
          iabot(i,j) = iabot(i,j+1)
 
        end do
      end do
 
      return
      end
      subroutine delrow ( a, iabot, iatop, irow, maxcol, maxrow, ncol,
     &  nrow )

c*********************************************************************72
c
cc DELROW deletes a row by shifting other rows up.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, real A(MAXROW,MAXCOL).  A is the matrix to be changed.
c
c    Output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
c    IATOP and IABOT represent the rational or decimal matrix
c    to be changed.
c
c    Input, integer IROW, the row to be deleted.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input, integer NCOL, the number of columns in the matrix.
c
c    Input, integer NROW, the number of rows in the matrix.
c
      integer maxcol
      integer maxrow
c
      real a(maxrow,maxcol)
      integer i
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer irow
      integer j
      integer ncol
      integer nrow
c
      do i = irow, nrow-1
        do j = 1, ncol
          a(i,j) = a(i+1,j)
          iatop(i,j) = iatop(i+1,j)
          iabot(i,j) = iabot(i+1,j)
        end do
      end do
 
      return
      end
      subroutine digten ( tenval, intval )
c
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
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*1 TENVAL, the decimal digit, '0'
c    through '9'.
c
c    Output, integer INTVAL, the corresponding integer value.
c
      integer intval
      character*1 tenval
c
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
      subroutine divide ( a, dete, iatop, iabot, idetop, idebot, 
     &  ierror, iform, imat, iounit, line, maxcol, maxdig, maxint, 
     &  maxrow, ncol, ndig, nline, nrow, output, prompt)
c
c*********************************************************************72
c
cc DIVIDE divides one row of the matrix by a nonzero value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, real A(MAXROW,MAXCOL).  A is the matrix
c    whose row is to be divided.
c
c    Input/output, real DETE, the determinant of the product of the
c    elementary row operations applied to the current matrix.
c
c    Input/output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
c    IATOP and IABOT represent the rational or decimal matrix
c    whose row is to be divided.
c
c    Input/output, integer IDETOP, IDEBOT, the rational or
c    decimal representation of the determinant of the product of
c    the elementary row operations applied to the current matrix.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IFORM, specifies the arithmetic being used.
c    0=rational, 1=real, 2=decimal.
c
c    Input, integer IMAT.
c    0, no matrix has been defined by the user.
c    1, a matrix has been defined by the user.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Workspace, character*80 LINE.
c    Used to hold the user's input.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXDIG, the maximum number of decimals to use.
c
c    Input, integer MAXINT, the maximum integer representable, which
c    should probably be 2147483647.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input, integer NCOL, the number of columns in the matrix.
c
c    Input, integer NDIG, the number of decimal digits used.
c
c    Input/output, integer NLINE.
c    Keeps track of the number of useful characters in LINE.
c
c    Input, integer NROW, the number of rows in the matrix.
c
c    Workspace, character*100 OUTPUT.
c
c    Workspace, character*80 PROMPT.
c
      integer maxcol
      integer maxrow
c
      real a(maxrow,maxcol)
      real dete
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer idebot
      integer idetop
      integer ierror
      integer iform
      integer imat
      integer iounit(4)
      integer irow
      integer isbot
      integer istop
      character*80 line
      integer maxdig
      integer maxint
      integer ncol
      integer ndig
      integer nline
      integer nrow
      character*100 output
      character*80 prompt
      real sval
c
      if ( imat .eq. 0 ) then
        ierror = 1
        output = 'You must set up a matrix first!'
        call chrwrt ( iounit, output )
        return
      end if

      prompt = 'row I, divisor S.'
c
c  Read the row number to be divided.
c
      call intrea ( irow, line, nline, prompt, iounit, ierror )
      if ( ierror .ne. 0 ) return
c
c  Read the divisor, either RVAL or ISTOP/ISBOT.
c
      if ( iform .eq. 0 ) then
 
        call ratrea ( istop, isbot, sval, line, nline, prompt, iounit,
     &    ierror )
 
      else if ( iform .eq. 1 ) then
 
        call relrea ( sval, line, nline, prompt, iounit, ierror )
 
      else if ( iform .eq. 2 ) then
 
        call decrea ( istop, isbot, sval, line, maxdig, nline, prompt,
     &    iounit, ierror )
 
        call deccut ( istop, isbot, ndig )
 
      end if
 
      if ( ierror .ne. 0 ) return
c
c  Divide the row by the divisor.
c
      call scadiv ( a, iatop, iabot, ierror, iform, iounit, irow,
     &  maxcol, maxint, maxrow, ncol, ndig, nrow, output, sval,
     &  istop, isbot )

      if ( ierror .ne. 0 ) return
c
c  Update the ERO determinant.
c
      if ( iform .eq. 0 ) then

        call ratdiv ( idebot, idebot, isbot, ierror, idetop,
     &    idetop, istop, maxint )

      else if ( iform .eq. 1 ) then

        dete = dete / sval

      else if ( iform .eq. 2 ) then

        call decdiv ( idebot, idebot, isbot, ierror, idetop, idetop,
     &    istop, ndig )

      end if
 
      return
      end
      subroutine evjaco ( a, ibase, ierror, iform, imat, iounit, line,
     &  lpmoda, maxcol, maxrow, ncol, nline, nrow, output, prompt )
c
c*********************************************************************72
c
cc EVJACO carries out a Jacobi rotation on a square matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, real A(MAXROW,MAXCOL).  A is the matrix
c    on which Jacobi rotation is carried out.
c
c    Input, integer IBASE(MAXROW), keeps track of basic variables.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IFORM, specifies the arithmetic being used.
c    0=rational, 1=real, 2=decimal.
c
c    Input, integer IMAT.
c    0, no matrix has been defined by the user.
c    1, a matrix has been defined by the user.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Workspace, character*80 LINE.
c    Used to hold the user's input.
c
c    Input, integer LPMODA.
c    0, the program is in linear algebra mode.
c    1, the program is in linear programming mode.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input, integer NCOL, the number of columns in the matrix.
c
c    Input/output, integer NLINE.
c    Keeps track of the number of useful characters in LINE.
c
c    Input, integer NROW, the number of rows in the matrix.
c
c    Workspace, character*100 OUTPUT.
c
c    Workspace, character*80 PROMPT.
c
      integer maxcol
      integer maxrow
c
      real a(maxrow,maxcol)
      real cj
      integer i
      integer ibase(maxrow)
      integer ierror
      integer iform
      integer ihi
      integer ilo
      integer imat
      integer iounit(4)
      integer j
      integer jhi
      integer jlo
      integer k
      logical leqi
      character*80 line
      integer lpmoda
      integer ncol
      integer nline
      integer nrow
      character*100 output
      character*80 prompt
      real sj
      logical sym
      real t1
      real t2
      real temp
      character*80 title
      real tj
      real u
c
c  Return if no matrix has been entered yet.
c
      if ( imat .ne. 1 ) then
        output = 'You must first enter a matrix with the "E" command!'
        call chrwrt ( iounit, output )
        ierror = 1
        return
      end if
c
c  Return if in linear programming mode.
c
      if ( lpmoda .ne. 0 ) then
        output = 'Please exit linear programming mode with the "L" '
     &    // 'command!'
        call chrwrt ( iounit, output )
        ierror = 1
        return
      end if
c
c  Return if using rational arithmetic.
c
      if ( iform .ne. 1 ) then
        ierror = 1
        output = 'Jacobi iteration requires real arithmeticc'
        call chrwrt ( iounit, output )
        output = 'Please issue the command "F R" first!'
        call chrwrt ( iounit, output )
        return
      end if
c
c  Return if matrix is not square.
c
      if ( nrow .ne. ncol ) then
        output = 'Jacobi iteration requires a square matrix!'
        ierror = 1
        call chrwrt ( iounit, output )
        return
      end if
c
c  Test for symmetry.
c
      sym = .true.
      do i = 1, min ( nrow, ncol )
        do j = 1, i-1
          if ( a(i,j) .ne. a(j,i) ) then
            sym = .false.
          end if
        end do
      end do
 
      if ( .not. sym ) then
        output = 'Warning!  Because the matrix is not symmetric,'
        call chrwrt ( iounit, output )
        output = 'Jacobi''s method may not converge!'
        call chrwrt ( iounit, output )
      end if
c
c  Here is where repetition will begin.
c
30    continue
c
c  Print the current matrix.
c
      ilo = 1
      ihi = nrow
      jlo = 1
      jhi = ncol
      title = 'The current matrix'

      call relprn ( a, ibase, iounit, ihi, ilo, jhi, jlo, lpmoda,
     &  maxcol, maxrow, ncol, nrow, output, title )
c
c  Get the row of the entry.
c
40    continue

      output = ' '
      call chrwrt ( iounit, output )
      prompt = 'row I, column J, or "Q" to quit.'
      call intrea ( i, line, nline, prompt, iounit, ierror )
 
      if ( ierror .ne. 0 ) then
        if ( leqi ( line(1:1), 'Q' ) ) then
          nline = 0
          ierror = 0
        end if
        return
      end if
 
      if ( i .le. 0 .or. i .gt. nrow ) then
        output = 'The value of I, the row index, is illegal.'
        call chrwrt ( iounit, output )
        go to 40
      end if
c
c  Get the column of the entry.
c
50    continue

      prompt = 'column J or "Q" to quit.'
      call intrea ( j, line, nline, prompt, iounit, ierror )
 
      if ( ierror .ne. 0 ) then
        if ( leqi ( line(1:1), 'Q' ) ) then
          nline = 0
          ierror = 0
        end if
        return
      end if
 
      if ( j .le. 0 .or. j .gt. ncol ) then
        output = 'The value of J, the column index, is illegal.'
        call chrwrt ( iounit, output )
        go to 50
      end if
c
c  I and J must not be equal.
c
      if ( i .eq. j ) then
        output = 'Jacobi rotations require I and J to be distinct!'
        call chrwrt ( iounit, output )
        go to 40
      end if
c
c  A(I,J) should not already be zero.
c
      if ( a(i,j) .eq. 0.0 ) then
        output = 'A(I,J) is already zero!'
        call chrwrt ( iounit, output )
        go to 40
      end if
c
c  If matrix is nonsymmetric, we require that A(I,J)+A(J,I)
c  not be zero.
c
      if ( a(i,j) + a(j,i) .eq. 0.0 ) then
        output = 'A(I,J)+A(J,I) is zero!'
        call chrwrt ( iounit, output )
        go to 40
      end if
c
c  Compute CJ and SJ.
c
      u = ( a(j,j) - a(i,i) ) / ( a(i,j) + a(j,i) )
 
      if ( u .ge. 0 ) then
        temp = 1.0
      else
        temp = -1.0
      end if
 
      tj = temp / ( abs(u) + sqrt(u**2+1.0) )
      cj = 1.0 / sqrt ( tj**2 + 1.0 )
      sj = tj * cj
c
c  Premultiply by Q transpose.
c
      do k = 1, ncol
        t1 = a(i,k)
        t2 = a(j,k)
        a(i,k) = cj * t1 - sj * t2
        a(j,k) = sj * t1 + cj * t2
      end do
c
c  Postmultiply by Q.
c
      do k = 1, nrow
        t1 = a(k,i)
        t2 = a(k,j)
        a(k,i) = cj * t1 - sj * t2
        a(k,j) = sj * t1 + cj * t2
      end do
c
c  Force A(I,J) and A(J,I) to be zero.
c
      a(i,j) = 0.0
      a(j,i) = 0.0
 
      go to 30
 
      end
      subroutine filadd ( filnam, ierror, inew, iold, iounit, nrec,
     &  output )
c
c*********************************************************************72
c
cc FILADD allows us to append information to a pre-existing file.
c
c
c  FILADD was created to address the fact that ANSI FORTRAN
c  does not let one easily append information to a sequential
c  access file once it has been closed.  In order to allow a user
c  to append new information, we create a new, writeable copy
c  of the file by means of a temporary copy.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*60 FILNAM, the name of the old file.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer INEW, the unit number on which the new copy
c    should be opened.
c
c    Input, integer IOLD, the unit number on which the old file
c    should be opened.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Output, integer NREC, the number of records in the old file.
c
c    Workspace, character*100 OUTPUT.
c
      character*6 chrint
      character*60 filnam
      character*60 filtmp
      integer ierror
      integer inew
      integer iold
      integer iounit(4)
      character*80 line
      integer nrec
      character*100 output
c
      filtmp = 'tmpfil.dat'
      ierror = 0
      nrec = 0
c
c  Open old file as readable.  If it doesn't exist, we can
c  skip ahead.  Otherwise, also open new file as writeable.
c
      open ( unit = iold, file = filnam, status = 'old', err = 50 )
      rewind iold
      open ( unit = inew, file = filtmp, status = 'new', err = 60 )
c
c  Copy old into temporary, then delete old.
c
10    continue

      read ( iold, '(a80)', end = 20 ) line
      nrec = nrec + 1
      write ( inew, '(a80)' ) line
      go to 10
 
20    continue

      output = 'The file contains ' // chrint(nrec) // ' lines.'
      call chrdb2 ( output )
      call chrwrt ( iounit, output )
      close ( unit = iold, status = 'delete' )
      close ( unit = inew )
c
c  Reopen old as writeable, write copy of temporary into it.
c
      open ( unit = iold, file = filnam, status = 'new', err = 60 )
      open ( unit = inew, file = filtmp, status = 'old', err = 60 )
      rewind inew
 
30    continue

      read ( inew, '(a80)', end = 40 ) line
      write ( iold, '(a80)' ) line
      go to 30
 
40    continue
c
c  Delete temporary file and return.
c
      close ( unit = inew, status = 'delete' )
      return
c
c  The file does not exist.  We may write into it immediately.
c
50    continue

      output = 'Creating a new file.'
      call chrwrt ( iounit, output )
      open ( unit = iold, file = filnam, status = 'new', err = 60 )
      return
c
c  Delete old copy of FILTMP.
c
60    continue

      ierror = 1
      output = 'The file could not be opened!'
      call chrwrt ( iounit, output )
 
      return
      end
      subroutine filred ( filex, ierror, iform, iounit, line, lpmoda,
     &  nline, output, prompt )
c
c*********************************************************************72
c
cc FILRED reads an example from a file.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character*60 FILEX.
c    On input, the default name of the example file.
c    On output, the chosen name of the example file.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IFORM, specifies the arithmetic being used.
c    0=rational, 1=real, 2=decimal.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Workspace, character*80 LINE.
c    Used to hold the user's input.
c
c    Output, integer LPMODA.
c    0, the program is in linear algebra mode.
c    1, the program is in linear programming mode.
c
c    Input/output, integer NLINE.
c    Keeps track of the number of useful characters in LINE.
c
c    Workspace, character*100 OUTPUT.
c
c    Input, character*80 PROMPT, the prompt string.
c
      character*6 chrint
      character*60 filex
      character*60 filnam
      integer ierror
      integer iform
      integer ilabel
      integer iounit(4)
      integer iterm
      integer jform
      integer jlabel
      integer jpmoda
      integer lchar
      integer lenchr
      logical leqi
      character*80 line
      integer lpmoda
      integer nline
      character*100 output
      character*80 prompt
      character*60 ylabel
c
      ierror = 0
      lchar = lenchr(filex)
 
      prompt = 'filename to read, default= "' // filex(1:lchar) 
     &  // '".'
      call chrdb2 ( prompt )
      iterm = 0

      call chrrea ( filnam, line, nline, prompt, iounit, ierror, iterm )

      if ( ierror .ne. 0 ) return
 
      if ( filnam(1:1) .ne. ' ' ) then
         filex = filnam
      end if
 
      open ( unit = 41, file = filex, status = 'old', err = 70 )
      iounit(1) = 41
c
c  Read and print labels.
c
      ilabel = 0
      ylabel = 'to cancel.'
      output = chrint(ilabel) // ' ' // ylabel
      call chrdb2 ( output )
      call chrwrt ( iounit, output )
      prompt = ' '
 
10    continue

      nline = 0
      iterm = 0

      call chrrea ( ylabel, line, nline, prompt, iounit, ierror, 
     &  iterm )

      if ( ierror .ne. 0 ) go to 20
      call chrdb2 ( ylabel )
      call capchr ( ylabel(1:6) )
 
      if ( leqi ( ylabel(1:6), 'LABEL:') ) then
        ilabel = ilabel + 1
        ylabel(1:6) = ' '
        output = chrint (ilabel) // ' ' // ylabel
        call chrdb2 ( output )
        call chrwrt ( iounit, output )
      end if
 
      go to 10
c
c  Close file.
c
20    continue

      ierror = 0
      close ( unit = iounit(1) )
      iounit(1) = 0
c
c  Get example number from user.
c
30    continue

      nline = 0
      prompt = 'example number.'
      call intrea ( jlabel, line, nline, prompt, iounit, ierror )
      if ( ierror .ne. 0 ) return
 
      if ( jlabel .le. 0 .or. 
     &   jlabel .gt. ilabel ) then
        output = 'Your choice was not acceptable.'
        call chrwrt ( iounit, output )
        go to 30
      end if
c
c  Reopen file, seek that example number.
c
      ilabel = 0
      open ( unit = 41, file = filex, status = 'old', err = 70 )
      iounit(1) = 41
      prompt = ' '
 
40    continue

      nline = 0
      iterm = 0

      call chrrea ( ylabel, line, nline, prompt, iounit, ierror, 
     &  iterm )

      if ( ierror .ne. 0 ) go to 50
      call chrdb2 ( ylabel )
      call capchr ( ylabel(1:6) )
 
      if ( leqi ( ylabel(1:6), 'LABEL:' ) ) then
        ilabel = ilabel + 1
        if ( ilabel .eq. jlabel ) go to 60
      end if
 
      go to 40
 
50    continue

      ierror = 1
      close ( unit = iounit(1) )
      iounit(1) = 0
      output = 'Could not retrieve example.'
      call chrwrt ( iounit, output )
      return
 
60    continue
c
c  See if the arithmetic mode should be changed.
c
      nline = 0
      call intrea ( jform, line, nline, prompt, iounit, ierror )

      if ( ierror .ne. 0 ) then
        iounit(1) = 0
        return
      end if
 
      if ( iform .ne. jform ) then
 
        if ( jform .eq. 0 ) then
          output = 'Arithmetic switched to rational form.'
          iform = jform
        else if ( jform .eq. 1 ) then
          output = 'Arithmetic switched to real form.'
          iform = jform
        else if ( jform .eq. 2 ) then
          output = 'Arithmetic switched to decimal form.'
          iform = jform
        else
          output = 'Illegal value for arithmetic switch:' //
     &      chrint(jform)
          call chrwrt ( iounit, output )
          ierror = 1
          iounit(1) = 0
          return
        end if
 
        call chrwrt ( iounit, output )
 
      end if
c
c  See if linear programming mode should be changed.
c
      nline = 0
      call intrea ( jpmoda, line, nline, prompt, iounit, ierror )
      if ( ierror .ne. 0 ) then
        iounit(1) = 0
        return
      end if
 
      if ( lpmoda .ne. jpmoda ) then
        if ( jpmoda .eq. 0 ) then
          output = 'Switching to linear algebra mode.'
          lpmoda = jpmoda
        else if ( jpmoda .eq. 1 ) then
          output = 'Switching to linear programming mode.'
          lpmoda = jpmoda
        else
          output = 'Illegal value for linear programming mode:'
     &      // chrint(jpmoda)
          call chrwrt ( iounit, output )
          iounit(1) = 0
          ierror = 1
          return
        end if
        call chrwrt ( iounit, output )
      end if
 
      nline = 0
      return
 
70    continue

      ierror = 1
      output = 'Error!  The example file could not be opened!'
      call chrwrt ( iounit, output )
 
      return
      end
      subroutine filwrt ( a, chineq, filex, iatop, iabot, ierror,
     &  iform, imat, iounit, line, lpmoda, maxcol, maxrow, nart,
     &  ncol, nline, nrow, nvar, output, prompt )
c
c*********************************************************************72
c
cc FILWRT writes an example to a file.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real A(MAXROW,MAXCOL).  A is the current matrix.
c
c    Input, character*1 CHINEQ(MAXROW), the '<', '=', or '>'
c    sign for each linear programming constraint.
c
c    Input/output, character*60 FILEX.
c    On input, the default name of the example file.
c    On output, the chosen name of the example file.
c
c    Input, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
c    IATOP and IABOT represent the current rational or decimal matrix.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IFORM, specifies the arithmetic being used.
c    0=rational, 1=real, 2=decimal.
c
c    Input, integer IMAT.
c    0, no matrix has been defined by the user.
c    1, a matrix has been defined by the user.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Workspace, character*80 LINE.
c    Used to hold the user's input.
c
c    Input, integer LPMODA.
c    0, the program is in linear algebra mode.
c    1, the program is in linear programming mode.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input, integer NART, the number of artificial variables.
c
c    Input, integer NCOL, the number of columns in the matrix.
c
c    Input/output, integer NLINE.
c    Keeps track of the number of useful characters in LINE.
c
c    Input, integer NROW, the number of rows in the matrix.
c
c    Input, integer NVAR, the number of basic variables.
c
c    Workspace, character*100 OUTPUT.
c
c    Workspace, character*80 PROMPT.
c
      integer maxcol
      integer maxrow
c
      real a(maxrow,maxcol)
      character*1 chineq(maxrow)
      character*6 chrint
      character*60 filex
      character*60 filnam
      integer i
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer ierror
      integer iform
      integer ii
      integer imat
      integer inew
      integer iold
      integer iosave
      integer iounit(4)
      character*1 isay
      integer iterm
      integer j
      integer jinc
      integer k
      integer khi
      integer lchar
      integer lenchr
      character*80 line
      integer lpmoda
      integer nart
      integer ncol
      integer nline
      integer nrec
      integer nrow
      integer nvar
      character*100 output
      character*80 prompt
      character*60 xlabel
c
c  Quit immediately if there are no matrices to write.
c
      if ( imat .ne. 1 ) then
        ierror = 1
        output = 'Operation canceled!'
        call chrwrt ( iounit, output )
        output = 'There is no matrix to write!'
        call chrwrt ( iounit, output )
        return
      end if
c
c  Get the filename to use.
c
      lchar = lenchr(filex)
      prompt = 'file to use, default= "' // filex(1:lchar) // '".'
      call chrdb2 ( prompt )
      iterm = 0
      call chrrea ( filnam, line, nline, prompt, iounit, ierror,
     &  iterm )
      if ( ierror .ne. 0 ) return
      if ( filnam(1:1) .ne. ' ' ) then
        filex = filnam
      end if
c
c  Get the label to use.
c
      nline = 0
      prompt = 'label.'
      iterm = 0
      call chrrea ( xlabel, line, nline, prompt, iounit, ierror, 
     &  iterm )
      if ( ierror .ne. 0 ) return
c
c  Open the file, whether it is new or old, and prepare to write
c  the new information at the end of the file.
c
      iold = 31
      inew = 32
      call filadd ( filex, ierror, inew, iold, iounit, nrec, output )
      if ( ierror .ne. 0 ) return

      iounit(4) = 31
      output = 'label:    '//xlabel
      call chrwrt ( iounit, output )
      iounit(2) = -1
      iosave = iounit(3)
      iounit(3) = -1

      output = chrint(iform) // 
     &  ', iform (0 fraction, 1 real, 2 decimal)'
      call chrwrt ( iounit, output )

      output = chrint(lpmoda) // ', lpmode (1 for linear programming)'
      call chrwrt ( iounit, output )

      if ( lpmoda .eq. 0 ) then
        output = chrint(nrow) // ',' // chrint(ncol)
      else
        output = chrint(nrow-1) // ',' // chrint(nvar)
      end if

      call chrdb1 ( output )
      call chrwrt ( iounit, output )
 
      if ( lpmoda .eq. 0 ) then
 
        do i = 1, nrow
 
          isay = ' '
 
          if ( iform .eq. 0 ) then
            jinc = 3
          else if ( iform .eq. 1 ) then
            jinc = 5
          else if ( iform .eq. 2 ) then
            jinc = 3
          end if
 
          do j = 1, ncol, jinc
 
            khi = min ( j+jinc-1, ncol )
 
            if ( iform .eq. 0 ) then
              write ( output, 70 ) isay,
     &          ( iatop(i,k), iabot(i,k), k = j, khi )
            else if ( iform .eq. 1 ) then
              write ( output, 80 ) isay, ( a(i,k), k = j, khi )
            else if ( iform .eq. 2 ) then
              write ( output, 70 ) isay, 
     &          ( iatop(i,k), iabot(i,k), k = j, khi )
            end if
 
            call chrdb2 ( output )
            call chrwrt ( iounit, output )
 
          end do
 
        end do
 
      else
 
        do i = 1, nrow
 
          ii = i

          if ( i .eq. nrow .and. nart .gt. 0) then
            ii = nrow + 1
          end if
 
          if ( i .eq. nrow ) then
 
            do j = 1, nvar
              iatop(ii,j) = - iatop(ii,j)
              a(ii,j) = - a(ii,j)
            end do
 
          end if
 
          isay = chineq(i)
 
          if ( iform .eq. 0 ) then
            jinc = 3
          else if ( iform .eq. 1 ) then
            jinc = 5
          else if ( iform .eq. 2 ) then
            jinc = 3
          end if
 
          do j = 1, nvar, jinc
 
            khi = min ( j+jinc-1, nvar )
 
            if ( iform .eq. 0 ) then
              write ( output, 70 )isay,
     &          ( iatop(ii,k), iabot(ii,k), k = j, khi )
            else if ( iform .eq. 1 ) then
              write ( output, 80 ) isay, ( a(ii,k), k = j, khi )
            else if ( iform .eq. 2 ) then
              write(output,70)isay, 
     &          ( iatop(ii,k), iabot(ii,k), k= j, khi )
            end if
 
            isay=' '
            call chrdb2 ( output )
            call chrwrt ( iounit, output )
 
          end do
 
          if ( iform .eq. 0 ) then
            write ( output, 70 ) isay, iatop(ii,ncol), iabot(ii,ncol)
          else if ( iform .eq. 1 ) then
            write ( output, 80 ) isay, a(ii,ncol)
          else if ( iform .eq. 2 ) then
            write ( output, 70 ) isay, iatop(ii,ncol), iabot(ii,ncol)
          end if
 
          call chrdb2 ( output )
          call chrwrt ( iounit, output )
 
          if ( i .eq. nrow ) then
 
            do j = 1, nvar
              iatop(ii,j) = - iatop(ii,j)
              a(ii,j) = - a(ii,j)
            end do
 
          end if
 
        end do
 
      end if
c
c  Write one blank line to avoid end-of-file problems on Macintosh.
c
      output = ' '
      call chrwrt ( iounit, output )
 
      iounit(2) = 0
      iounit(3) = iosave
      close(unit = iounit(4))
      iounit(4) = - 1
      output = 'The problem has been stored.'
      call chrwrt ( iounit, output )
 
      return
70    format(a1,3(i12,'/',i12,','))
80    format(a1,5(g14.6,','))
      end
      subroutine flushl ( string )
c
c*********************************************************************72
c
cc FLUSHL flushes a string left.
c
c
c  Both blanks and tabs are treated as "white space".
c
c  Example:
c
c    Input             Output
c
c    '     Hello'      'Hello     '
c    ' Hi there!  '    'Hi there!   '
c    'Fred  '          'Fred  '
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character*(*) STRING.
c
c    On input, STRING is a string of characters.
c
c    On output, any initial blank or tab characters in STRING
c    have been cut.
c
      character*1 TAB
      parameter ( TAB = char(9) )
c
      integer i
      integer lchar
      integer lenchr
      integer nonb
      character*(*) string
c
c  Check the length of the string to the last nonblank.
c  If nonpositive, return.
c
      lchar = lenchr ( string )

      if ( lchar .le. 0 ) then
        return
      end if
c
c  Find NONB, the location of the first nonblank, nontab.
c
      do i = 1, lchar

        nonb = i

        if ( string(i:i) .ne. ' ' .and. 
     &    string(i:i) .ne. TAB ) then
          go to 10
        end if

      end do

      string = ' '
 
      return
 
10    continue
c
c  Shift the string left.
c
      if ( nonb .gt. 1 ) then
        do i = 1, lchar + 1 - nonb
          string(i:i) = string(i+nonb-1:i+nonb-1)
        end do
      end if
c
c  Blank out the end of the string.
c
      string(lchar+2-nonb:lchar) = ' '
 
      return
      end
      subroutine flushr ( string )
c
c*********************************************************************72
c
cc FLUSHR flushes a string right.
c
c
c  Example:
c
c    Input             Output
c
c    'Hello     '      '     Hello'
c    ' Hi there!  '    '   Hi there!'
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character*(*) STRING.
c
c    On input, STRING is a string of characters.
c
c    On output, any trailing blank characters in STRING
c    have been cut, and pasted back onto the front.
c
      integer i
      integer lchar
      integer lenchr
      integer nonb
      character*(*) string
c
c  Check the full length of the string.
c
      lchar = len ( string )
c
c  Find the occurrence of the last nonblank.
c
      nonb = lenchr ( string )
c
c  Shift the string right.
c
      do i = lchar, lchar + 1 - nonb, -1
        string(i:i) = string(i-lchar+nonb:i-lchar+nonb)
      end do
c
c  Blank out the beginning of the string.
c
      string(1:lchar-nonb) = ' '
 
      return
      end
      subroutine form ( a, b, c, dete, iatop, iabot, ibtop, ibbot,
     &  ictop, icbot, idetop, idebot, iform, imat, iounit, jform,
     &  maxcol, maxint, maxrow, ndig, output )
c
c*********************************************************************72
c
cc FORM converts from one arithmetic form to another.
c
c
c  On input, IFORM contains a code for the current arithmetic
c  form, and JFORM contains the code for the new form.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, real A(MAXROW,MAXCOL), B(MAXROW,MAXCOL),
c    C(MAXROW,MAXCOL), the current real matrix, and its two
c    backup copies.
c
c    Input/output, real DETE, the determinant of the product of the
c    elementary row operations applied to the current matrix.
c
c    Input/output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL),
c    IBTOP(MAXROW,MAXCOL), IBBOT(MAXROW,MAXCOL), ICTOP(MAXROW,MAXCOL),
c    ICBOT(MAXROW,MAXCOL), the current fractional or decimal matrix
c    and its two backup copies.
c
c    Input/output, integer IDETOP, IDEBOT, the rational or
c    decimal representation of the determinant of the product of
c    the elementary row operations applied to the current matrix.
c
c    Input/output, integer IFORM, specifies the arithmetic being
c    used.  0=rational, 1=real, 2=decimal.
c
c    Input, integer IMAT.
c    0, no matrix has been defined by the user.
c    1, a matrix has been defined by the user.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Input, integer JFORM, the arithmetic to be converted to.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXINT, the maximum integer representable, which
c    should probably be 2147483647.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input, integer NDIG, the number of decimal digits used.
c
c    Workspace, character*100 OUTPUT.
c
      integer maxcol
      integer maxrow
c
      real a(maxrow,maxcol)
      real b(maxrow,maxcol)
      real c(maxrow,maxcol)
      real dete
      integer i
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer ibbot(maxrow,maxcol)
      integer ibtop(maxrow,maxcol)
      integer icbot(maxrow,maxcol)
      integer ictop(maxrow,maxcol)
      integer idebot
      integer idetop
      integer iform
      integer imat
      integer iounit(4)
      integer j
      integer jform
      integer maxint
      integer ndig
      character*100 output
c
      if ( iform .eq. jform ) then
        output = 'You are already using the arithmetic type that'
        call chrwrt ( iounit, output )
        output = 'you have requested.'
        call chrwrt ( iounit, output )
        return
      end if
c
c  Tell the user what we think we're doing.
c
      if ( jform .eq. 0 ) then
 
        output = 'Converting to fractional arithmetic.'
        call chrwrt ( iounit, output )
 
        call ratwrn ( iounit, maxint, output )
 
      else if ( jform .eq. 1 ) then
 
        output = 'Converting to real arithmetic.'
        call chrwrt ( iounit, output )
 
        call relwrn ( iounit, output )
 
      else if ( jform .eq. 2 ) then
 
        output = 'Converting to decimal arithmetic.'
        call chrwrt ( iounit, output )
 
        call decwrn ( iounit, output )
 
      end if
c
c  If there's no matrix, then just set the arithmetic mode
c  and return.
c
      if ( imat .eq. 0 ) then
        iform = jform
        return
      end if
c
c  Convert the matrix data.
c  In a special case, there is data in the matrix beyond row NROW.
c  (Linear programming, with artificial variables).
c  So just convert MAXROW by MAXCOL, rather than the more modest
c  NROW by NCOL.
c
      if ( jform .eq. 0 .and. iform .eq. 1 ) then
 
        do i = 1, maxrow
          do j = 1, maxcol
 
            call relrat ( a(i,j), iatop(i,j), iabot(i,j), ndig )
            call relrat ( b(i,j), ibtop(i,j), ibbot(i,j), ndig )
            call relrat ( c(i,j), ictop(i,j), icbot(i,j), ndig )
 
          end do
        end do
 
        call relrat ( dete, idetop, idebot, ndig )
 
      else if ( jform .eq. 0 .and. iform .eq. 2 ) then
 
        do i = 1, maxrow
          do j = 1, maxcol
 
            call decrat ( iatop(i,j), iabot(i,j) )
            call decrat ( ibtop(i,j), ibbot(i,j) )
            call decrat ( ictop(i,j), icbot(i,j) )
  
          end do
        end do
 
        call decrat ( idetop, idebot )
 
      else if ( jform .eq. 1 .and. iform .eq. 0 ) then
 
        do i = 1, maxrow
          do j = 1, maxcol
 
            call ratrel ( a(i,j), iatop(i,j), iabot(i,j) )
            call ratrel ( b(i,j), ibtop(i,j), ibbot(i,j) )
            call ratrel ( c(i,j), ictop(i,j), icbot(i,j) )
 
          end do
        end do
 
        call ratrel ( dete, idetop, idebot )
 
      else if ( jform .eq. 1 .and. iform .eq. 2 ) then
 
        do i = 1, maxrow
          do j = 1, maxcol
 
            call decrel ( a(i,j), iatop(i,j), iabot(i,j) )
            call decrel ( b(i,j), ibtop(i,j), ibbot(i,j) )
            call decrel ( c(i,j), ictop(i,j), icbot(i,j) )
 
          end do
        end do
 
        call decrel ( dete, idetop, idebot )
 
      else if ( jform .eq. 2 .and. iform .eq. 0 ) then
 
        do i = 1, maxrow
          do j = 1, maxcol
 
            call ratdec ( iatop(i,j), iabot(i,j ), ndig )
            call ratdec ( ibtop(i,j), ibbot(i,j ), ndig )
            call ratdec ( ictop(i,j), icbot(i,j ), ndig )
 
          end do
        end do
 
        call ratdec ( idetop, idebot, ndig )
 
      else if ( jform .eq. 2 .and. iform .eq. 1 ) then
 
        do i = 1, maxrow
          do j = 1, maxcol
 
            call reldec ( a(i,j), iatop(i,j), iabot(i,j), ndig )
            call reldec ( b(i,j), ibtop(i,j), ibbot(i,j), ndig )
            call reldec ( c(i,j), ictop(i,j), icbot(i,j), ndig )
 
          end do
        end do
 
        call reldec ( dete, idetop, idebot, ndig )
 
      end if
c
c  Update the arithmetic form.
c
      iform = jform
 
      return
      end
      subroutine hello ( iounit, output )
c
c*********************************************************************72
c
cc HELLO greets the user on startup.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Workspace, character*100 OUTPUT.
c
      integer iounit(4)
      character*100 output
c
      output = ' '
      call chrwrt ( iounit, output )
      output = 'MATMAN, version 1.59'
      call chrwrt ( iounit, output )
      output = 'Last modified on 06 July 1998.'
      call chrwrt ( iounit, output )
      output = ' '
      call chrwrt ( iounit, output )
      output = 'An interactive program which carries out'
      call chrwrt ( iounit, output )
      output = 'elementary row operations on a matrix, or'
      call chrwrt ( iounit, output )
      output = 'the simplex method of linear programming.'
      call chrwrt ( iounit, output )
      output = ' '
      call chrwrt ( iounit, output )
      output = 'Developed by Charles Cullen and John Burkardt.'
      call chrwrt ( iounit, output )
      output = 'All rights reserved by the authors.  This program may'
      call chrwrt ( iounit, output )
      output =
     &  'not be reproduced in any form without written permission.'
      call chrwrt ( iounit, output )
      output = ' '
      call chrwrt ( iounit, output )
      output = 'Send comments to burkardt@psc.edu.'
      call chrwrt ( iounit, output )
 
      return
      end
      subroutine help ( iounit, output )
c
c*********************************************************************72
c
cc HELP prints out a brief list of the available commands.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Workspace, character*100 OUTPUT.
c
      integer iounit(4)
      character*100 output

      output = ' '
      call chrwrt ( iounit, output )
      output = 'Here is a list of all MATMAN commands:'
      call chrwrt ( iounit, output )
      output = ' '
      call chrwrt ( iounit, output )
 
      output = 'A      Add S times row I to row J.'
      call chrwrt ( iounit, output )
      output = 'B      Set up sample problem.'
      call chrwrt ( iounit, output )
      output = 'BASIC I, J changes basic variable I to J.'
      call chrwrt ( iounit, output )
      output = 'C      Change entry I, J to S.'
      call chrwrt ( iounit, output )
      output = 'D/M    Divide/Multiply row I by S.'
      call chrwrt ( iounit, output )
      output = 'DEC    Use decimal arithmetic.'
      call chrwrt ( iounit, output )
      output = 'DET    Print the determinant of the matrix.'
      call chrwrt ( iounit, output )
      output = 'E      Enter matrix with I rows and J columns.'
      call chrwrt ( iounit, output )
      output = 'E      Enter a linear programming problem, ' //
     &  'I constraints, J variables.'
      call chrwrt ( iounit, output )
      output = 'EDET   Print ERO determinant.'
      call chrwrt ( iounit, output )
      output = 'F      Choose arithmetic (Real, Fraction, or Decimal).'
      call chrwrt ( iounit, output )
      output = 'G      Add/delete a row or column of the matrix.'
      call chrwrt ( iounit, output )
      output = 'H      for quick help.'
      call chrwrt ( iounit, output )
      output = 'HELP   for full help (this list).'
      call chrwrt ( iounit, output )
      output = 'I      Interchange rows I and J.'
      call chrwrt ( iounit, output )
      output = 'J      Jacobi rotation in (I,J) plane.'
      call chrwrt ( iounit, output )
      output = 'K      Open/close the transcript file.'
      call chrwrt ( iounit, output )
      output = 'L      To switch between linear algebra and linear '
     &  //'programming.'
      call chrwrt ( iounit, output )
      output = 'MAXINT Set the maximum integer.'
      call chrwrt ( iounit, output )
      output = 'N      Set the number of decimal digits.'
      call chrwrt ( iounit, output )
      output = 'O      Check matrix for reduced row echelon form.'
      call chrwrt ( iounit, output )
      output = 'O      Check linear program tableau for optimality.'
      call chrwrt ( iounit, output )
      output = 'P      Pivot linear program, entering I, departing J.'
      call chrwrt ( iounit, output )
      output = 'Q      Quit.'
      call chrwrt ( iounit, output )
      output = 'R      Restore a saved matrix or tableau'
      call chrwrt ( iounit, output )
      output = 'RAT    Use rational arithmetic.'
      call chrwrt ( iounit, output )
      output = 'REAL   Use real arithmetic.'
      call chrwrt ( iounit, output )
      output = 'S      Store the current matrix or tableau.'
      call chrwrt ( iounit, output )
      output = 'T      Type out the matrix'
      call chrwrt ( iounit, output )
      output = 'TR     Transpose the matrix.'
      call chrwrt ( iounit, output )
      output = 'TS     Type linear programming solution.'
      call chrwrt ( iounit, output )
      output = 'U      Undo last operation.'
      call chrwrt ( iounit, output )
      output = 'V      Remove LP artificial variables.'
      call chrwrt ( iounit, output )
      output = 'W/X    Write/read example to/from file.'
      call chrwrt ( iounit, output )
      output = 'Y      Turn automatic printing ON or OFF.'
      call chrwrt ( iounit, output )
      output = 'Z      Automatic operation (requires password).'
      call chrwrt ( iounit, output )
      output = '#      Begins a comment line.'
      call chrwrt ( iounit, output )
      output = '<      Get input from a file.'
      call chrwrt ( iounit, output )
      output = '%/$    Turn paging on/off.'
      call chrwrt ( iounit, output )
      output = ' '
      call chrwrt ( iounit, output )
      output = 'R1 <=> R2  interchanges two rows'
      call chrwrt ( iounit, output )
      output = 'R1 <= S R1 multiplies a row by S.'
      call chrwrt ( iounit, output )
      output = 'R1 <= R1 + S R2 adds a multiple of another row.'
      call chrwrt ( iounit, output )
 
      return
      end
      subroutine hlpvms ( filhlp, iounit, line, nline, output, prompt )
c
c*********************************************************************72
c
cc HLPVMS provides extensive help from the MATMAN help file.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*60 FILHLP, the name of the help file.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Workspace, character*80 LINE, used to hold the user's input.
c
c    Input/output, integer NLINE.
c    Keeps track of the number of useful characters in LINE.
c
c    Workspace, character*100 OUTPUT.
c
c    Workspace, character*80 PROMPT.
c
      integer maxtop
      parameter (maxtop=40)
c
      character*75 choice
      character*75 ctemp
      character*75 ctemp2
      character*60 filhlp
      integer i
      integer ierror
      integer iline
      character*75 inline
      integer iounit(4)
      integer iterm
      integer itop
      integer jerror
      character*1 lab
      integer lchar
      integer lenc
      integer lenchr
      logical leqi
      integer level
      character*75 levelc(maxtop)
      integer levelm(10)
      integer levelo
      integer levelt(maxtop)
      integer lhunit
      character*80 line
      integer move
      integer nline
      integer ntop
      integer num
      character*100 output
      character*80 prompt
c
      ierror = 0
      lhunit = 55
      call setlin ( 0 )
c
c  Open help file
c
c  These lines work for a "private" copy of MATMAN on VAX/VMS,
c  UNIX, IBM PC or Macintosh
c
      open ( unit = lhunit, file = filhlp, status = 'old', err = 100 )
c
c  These lines work for a "shared" copy of MATMAN on a VAX/VMS
c  system:
c
c     open ( unit = lhunit, file = filhlp, status = 'old',err = 100,
c    &  shared, readonly )
c
      levelo = 0
      level = 1
      iline = 1
c
c  Move to beginning of current topic by reading MOVE lines from
c  the top of the file.  Record this position, corresponding to
c  the current LEVEL, in LEVELM, in case we later want to back up.
c
c  Print out the heading line of this topic.
c
10    continue

      jerror = 0
      move = iline
      levelm(level) = iline
 
      do i = 1, move-1
        read ( lhunit, '(1x)', end = 110, err = 110 )
      end do
 
      output = ' '
      call chrwrt ( iounit, output )
      read ( lhunit, '(a1,a75)', end = 110, err = 110 ) lab, inline
      output = inline
      call chrwrt ( iounit, output )
c
c  If 'going down' or redisplaying, (as opposed to backing up),
c  display information available under the current topic.
c
c  We stop printing when we hit a numeric label.
c
c  If this label is less than or equal to current level, there are
c  no subtopics.
c
c  Otherwise, we now move ahead to print out the list of subtopics
c  available for this topic.
c
      if ( level .ge. levelo ) then

        ntop = -1
 
30      continue
 
        read ( lhunit, '(a1,a75)', end = 50 ) lab, inline
        move = move + 1
 
        if ( lge(lab,'0') .and. lle(lab,'9') ) then
          read ( lab, '(i1)' ) num
          if ( num .le. level ) go to 50
          ntop = 0
          go to 40
        end if
 
        output = inline
        call chrwrt ( iounit, output )
        go to 30
      else
        ntop = 0
        inline = ' '
        lab = ' '
      end if
c
c  Locate each subtopic by examining column 1, searching for
c  integer label.
c
c  Assuming we are at level LEVEL, we are searching for labels
c  equal to LEVEL+1.  As we encounter each such label, we want to
c  store the rest of the line as a subtopic.  We ignore labels
c  greater than LEVEL+1 because these are sub-subtopics, and we
c  cease our search when we reach a label less than or equal to
c  LEVEL.
c
40    continue
 
      if ( lge(lab,'0') .and. lle(lab,'9') ) then

        read ( lab, '(i1)' ) num
        if ( num .le. level ) go to 50

        if ( num .eq. level+1 ) then

          ntop = ntop + 1
 
          if ( ntop .eq. 1 ) then
            output = ' '
            call chrwrt ( iounit, output )
            output = 'Help is available on:'
            call chrwrt ( iounit, output )
            output = ' '
            call chrwrt ( iounit, output )
          end if
 
          output = inline
          call chrwrt ( iounit, output )
          levelt(ntop) = move
          levelc(ntop) = inline

        end if

      end if

      read ( lhunit, '(a1,a75)', end = 50, err = 50 ) lab, inline
      move = move + 1
      go to 40
 
50    continue
c
c  Display subtopics.
c
      output = ' '
      call chrwrt ( iounit, output )
      output = 'Return to back up, ? to redisplay.'
      call chrwrt ( iounit, output )
c
c  Prompt for user choice of new topic, exit, or back up.
c
60    continue
 
      ierror = 0
      nline = 0
 
      if ( ntop .gt. 0 ) then
        prompt = 'topic you want help on, or RETURN or ?.'
      else
        prompt = 'RETURN or ?.'
      end if
 
      iterm = 0
      call chrrea ( choice, line, nline, prompt, iounit, ierror,
     &  iterm )

      if ( ierror .ne. 0 ) then
        ierror = 0
        close ( unit = lhunit )
        return
      end if
 
      call setlin ( 0 )
      call chrdb2 ( choice )
      lenc = lenchr ( choice )
      if ( lenc .le. 0 ) then
        choice = '!'
      end if
      ctemp = choice
c
c  Two errors in a row, OK, but three suggests that something is wrong.
c
      if ( ierror .ne. 0 ) then
        jerror = jerror + 1
        if ( jerror .le. 4) go to 60
        output = 'Too many input errors in a row!'
        call chrwrt ( iounit, output )
      end if
c
c  Consider ending this help session.
c
      if ( 
     &  (ctemp .eq. '!' .and. level .eq. 1 ) .or. 
     &  jerror .gt. 4 ) then
        close ( unit = lhunit )
        return
      end if
c
c  User wants to back up to a supertopic.  We must rewind.
c
      rewind lhunit
      levelo = level

      if ( ctemp .eq. '!' ) then

        level = level-1
        iline = levelm(level)
c
c  Redisplay current topic.
c
      else if ( ctemp .eq. '?' ) then

        go to 10
c
c  User wants to go down to a subtopic.
c
      else
 
        do i = 1, ntop

          ctemp2 = levelc(i)
          call chrdb2 ( ctemp2 )
          itop = i

          if ( leqi ( ctemp(1:lenc), ctemp2(1:lenc) ) ) then
            go to 90
          end if

        end do
 
        lchar = lenchr(choice)
        output = 'Sorry, no help available on "' //
     &    choice(1:lchar) // '".'
        call chrdb2 ( output )
        call chrwrt ( iounit, output )
        jerror = jerror + 1
        go to 60
 
90      continue

        level = level + 1
        iline = levelt(itop)

      end if
 
      go to 10
c
c  Error opening help file.
c
100   continue

      ierror = 1
      lchar = lenchr ( filhlp ) 
      output = 'Could not open the help file "' //
     &  filhlp(1:lchar) // '".'
      call chrdb2 ( output )
      call chrwrt ( iounit, output )
      return
c
c  Error reading help file.
c
110   continue

      ierror = 1
      lchar = lenchr(filhlp)
      output = 'Unexpected error while reading "' //
     &  filhlp(1:lchar) // '".'
      call chrdb2 ( output )
      call chrwrt ( iounit, output )
      close ( unit = lhunit )
 
      return
      end
      function igcf ( i, j )
c
c*********************************************************************72
c
cc IGCF finds the greatest common factor of I and J.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer I, J, two numbers whose greatest common factor 
c    is desired.
c
c    Output, integer IGCF, the greatest common factor of I and J.
c
c    Note that only the absolute values of I and J are
c    considered, so that IGCF is always nonnegative.
c
c    If I or J is 0, IGCF is returned as max ( 1, abs(I), abs(J) ).
c
c    If I and J have no common factor, IGCF is returned as 1.
c
c    Otherwise, using the Euclidean algorithm, IGCF is the
c    largest common factor of I and J.
c
      integer i
      integer igcf
      integer ip
      integer iq
      integer ir
      integer j
c
      igcf = 1
c
c  Return immediately if either I or J is zero.
c
      if ( i .eq. 0 ) then
        igcf = max ( 1, abs(j) )
        return
      else if ( j .eq. 0 ) then
        igcf = max ( 1, abs(i) )
        return
      end if
c
c  Set IP to the larger of I and J, IQ to the smaller.
c  This way, we can alter IP and IQ as we go.
c
      ip = max ( abs ( i ), abs ( j ) )
      iq = min ( abs ( i ), abs ( j ) )
c
c  Carry out the Euclidean algorithm.
c
10    continue
 
      ir = mod ( ip, iq )
 
      if ( ir .ne. 0 ) then
        ip = iq
        iq = ir
        go to 10
      end if
 
      igcf = iq
 
      return
      end
      subroutine indata ( op, var, ival )
c
c*********************************************************************72
c
cc INDATA stores and retrieves common data items.
c
c
c  INDATA works like a sort of COMMON block.  It stores or returns
c  the values of certain variables.  Thus, it allows routines
c  to "communicate" without having to have data passed up and
c  down the calling tree in argument lists.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) OP, describes the operation to be done.
c    'SET' means set a value.
c    'GET' means get a value.
c
c    Input, character*(*) VAR, the name of the variable to be set
c    or gotten.
c    VAR may have the value 'NLINE' or 'LPAGE'.
c
c    Input/output, integer IVAL.
c    If OP is 'SET', then the variable named in VAR is set to the
c    value IVAL.
c    If OP is 'GET', then the value of IVAL is set to the value of
c    the variable named in VAR.
c
      integer ival
      logical leqi
      integer lpage
      integer nline
      character*(*) op
      character*(*) var
c
      save lpage
      save nline
c
      data lpage /24/
      data nline /0/
c
      if ( leqi ( op, 'SET' ) .and. leqi ( var, 'NLINE' ) ) then
        nline = ival
      else if ( leqi ( op, 'SET' ) .and. leqi ( var, 'LPAGE' ) ) then
        lpage = ival
      else if ( leqi ( op, 'GET' ) .and. leqi ( var, 'NLINE' ) ) then
        ival = nline
      else if ( leqi ( op, 'GET' ) .and. leqi ( var, 'LPAGE' ) ) then
        ival = lpage
      end if
 
      return
      end
      function indexi ( string, sub )
c
c*********************************************************************72
c
cc INDEXI is a case-insensitive INDEX function.
c
c
c  It returns the location in STRING at which the substring SUB is
c  first found, or 0 if the substring does not occur at all.
c
c  INDEXI is also trailing blank insensitive.  This is very
c  important for those cases where you have stored information in
c  larger variables.  If STRING is of length 80, and SUB is of
c  length 80, then if STRING = 'FRED' and SUB = 'RED', a match would
c  not be reported by the standard FORTRAN INDEX, because it treats
c  both variables as being 80 characters long!  INDEXI assumes that
c  trailing blanks represent garbage!
c
c  Because of the suppression of trailing blanks, INDEXI cannot be
c  used to find, say, the first occurrence of the two-character
c  string 'A '.  However, INDEXI treats as a special case the
c  occurrence where STRING or SUB is entirely blank.  Thus you can
c  use INDEXI to search for occurrences of double or triple blanks
c  in a string, for example, although INDEX itself would be just as
c  suitable for that problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) STRING, the string to be searched.
c
c    Input, character*(*) SUB, the substring to search for.
c
c    Output, integer INDEXI.  0 if SUB does not occur in
c    STRING.  Otherwise STRING(INDEXI:INDEXI+LENS-1) = SUB,
c    where LENS is the length of SUB, and is the first place
c    this happens.  However, note that INDEXI ignores case,
c    unlike the standard FORTRAN INDEX function.
c
      integer i
      integer indexi
      integer lenchr
      logical leqi
      integer llen1
      integer llen2
      character*(*) string
      character*(*) sub
c
      indexi = 0

      llen1 = lenchr ( string )
      llen2 = lenchr ( sub )
c
c  In case STRING or SUB is blanks, use LEN.
c
      if ( llen1 .eq. 0 ) then
        llen1 = len ( string )
      end if

      if ( llen2 .eq. 0 ) then
        llen2 = len ( sub )
      end if

      if ( llen2 .gt. llen1 ) then
        return
      end if

      do i = 1, llen1 + 1 - llen2

        if ( leqi ( string(i:i+llen2-1), sub ) ) then
          indexi = i
          return
        end if

      end do

      return
      end
      subroutine infile ( filinp, ierror, iounit, line, lpage, nline,
     &  output, prompt )
c
c*********************************************************************72
c
cc INFILE handles a new input file.
c
c
c  INFILE reads the name of an input file, and changes the internal
c  values of IOUNIT(1), and opens the file.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character*(*) FILINP, the input file name.
c    On input, this is a default value, or the name of a previously
c    used input file.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Workspace, character*80 LINE.
c    Used to hold the user's input.
c
c    Output, integer LPAGE.
c    The number of lines per page.  Reset to 24 if the user is
c    now to be typing the input.
c
c    Input/output, integer NLINE.
c    Keeps track of the number of useful characters in LINE.
c
c    Workspace, character*100 OUTPUT.
c
c    Workspace, character*80 PROMPT.
c
      character*60 filnam
      character*60 filinp
      integer ierror
      integer iounit(4)
      integer iterm
      integer lchar
      integer lenchr
      character*80 line
      integer lpage
      integer nline
      character*100 output
      character*80 prompt
c
c  If we were already reading an input file, close it!
c
      if ( iounit(1) .ne. 0 ) then
        lchar = lenchr(filinp)
        output = 'Closing previous input file "'//filinp(1:lchar)//'".'
        call chrdb2 ( output )
        call chrwrt ( iounit, output )
        close ( unit = iounit(1) )
        iounit(1) = 0
      end if
c
c  Get the name of the input file.
c
      lchar = lenchr(filinp)
      prompt = 'file name, default= "' // filinp(1:lchar) // '".'
      call chrdb2 ( prompt )
      iterm = 0
      call chrrea ( filnam, line, nline, prompt, iounit, ierror,
     &  iterm )
      if ( ierror .ne. 0 ) return
c
c  If the input file is "*", then the user is typing input.
c
      if ( filnam .eq. '*' ) then
        output = 'MATMAN now expects input directly from the user.'
        call chrwrt ( iounit, output )
        lpage = 24
        call setpag ( lpage )
        output = 'Paging turned ON.'
        call chrwrt ( iounit, output )
        return
      end if
 
      if ( filnam .ne. ' ' ) then
        filinp = filnam
      end if
 
      iounit(1) = 11
      open ( unit = iounit(1), file = filinp, status = 'old', 
     &  err = 10 )
c
c  Turn paging off.
c
      lpage = 0
      call setpag ( lpage )
      output = 'Paging turned OFF.'
      call chrwrt ( iounit, output )
 
      lchar = lenchr(filinp)
      output = 'MATMAN now expects input from "'//filinp(1:lchar)//'".'
      call chrdb2 ( output )
      call chrwrt ( iounit, output )
      return
c
c  Opening failed.
c
10    continue

      ierror = 1
      iounit(1) = 0
      output = 'MATMAN could not open the input file!'
      call chrwrt ( iounit, output )
 
      return
      end
      subroutine inimat ( a, iabot, iatop, iform, maxcol, maxrow )
c
c*********************************************************************72
c
cc INIMAT initializes the matrix by zeroing it out.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, real A(MAXROW,MAXCOL).  A is the current matrix.
c
c    Output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
c    IATOP and IABOT represent the current rational or decimal
c    matrix.
c
c    Output, integer IFORM, specifies the arithmetic being used.
c    0=rational, 1=real, 2=decimal.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
      integer maxcol
      integer maxrow

      real a(maxrow,maxcol)
      integer i
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer iform
      integer j

      do i = 1, maxrow
        do j = 1, maxcol

          if ( iform .eq. 0 ) then
            iatop(i,j) = 0
            iabot(i,j) = 1
          else if ( iform .eq. 1 ) then
            a(i,j) = 0.0
          else if ( iform .eq. 2 ) then
            iatop(i,j) = 0
            iabot(i,j) = 1
          end if

        end do
      end do
 
      return
      end
      subroutine init ( a, autop, chineq, comnew, comold, dete,
     &  filex, filhlp, filinp, filkey, filtrn, iabot, iatop, iauthr,
     &  ibase, idebot, idetop, ierror, iform, imat, iounit, iprint, 
     &  iseed, islbot, isltop, line, lpage, lpmoda, maxcol, maxdig, 
     &  maxint, 
     &  maxrow, nart, ncol, ncon, ndig, nline, nrow, nslak, nvar, 
     &  sol )
c
c*********************************************************************72
c
cc INIT initializes the data.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, real A(MAXROW,MAXCOL).  A is the current matrix.
c
c    Output, LOGICAL AUTOP, .TRUE. if the matrix should be
c    automatically printed after most operations, .FALSE. otherwise.
c
c    Output, character*1 CHINEQ(MAXROW), the '<', '=', or '>'
c    sign for each linear programming constraint.
c
c    Output, character*20 COMNEW, the newest command from the user.
c
c    Output, character*20 COMOLD, the previous command from the user.
c
c    Output, real DETE, the determinant of the product of the
c    elementary row operations applied to the current matrix.
c
c    Output, character*60 FILEX, the default examples file.
c
c    Output, character*60 FILHLP, the default help file.
c
c    Output, character*60 FILINP, the default input file.
c
c    Output, character*60 FILKEY, the default password file.
c
c    Output, character*60 FILTRN, the default transcript file.
c
c    Output, integer IABOT(MAXROW,MAXCOL), IATOP(MAXROW,MAXCOL).
c    IATOP and IABOT represent the current rational or decimal
c    matrix.
c
c    Output, integer IAUTHR,
c    0 if the user has typed the correct password.
c    1 if the user has not typed the correct password.
c
c    Output, integer IBASE(MAXROW), keeps track of basic variables.
c
c    Output, integer IDETOP, IDEBOT, the rational or
c    decimal representation of the determinant of the product of
c    the elementary row operations applied to the current matrix.
c
c    Output, integer IERROR.
c    The error flag, which is initialized to zero by this routine.
c
c    Output, integer IFORM, specifies the arithmetic being used.
c    0=rational, 1=real, 2=decimal.
c
c    Output, integer IMAT.
c    0, no matrix has been defined by the user.
c    1, a matrix has been defined by the user.
c
c    Output, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Output, integer IPRINT.
c    0, if the most recent command does not require printout.
c    1, if the most recent command requires printout.
c
c    Output, integer ISEED, a random number generator seed.
c
c    Output, integer ISLBOT(MAXROW), ISLTOP(MAXROW), the decimal
c    or fractional representation of the linear programming solution.
c
c    Output, character*80 LINE,
c    a buffer used to hold the user's input.
c
c    Output, INTEGE LPAGE,
c    the number of lines per page.
c
c    Output, integer LPMODA.
c    0, the program is in linear algebra mode.
c    1, the program is in linear programming mode.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c         in the matrices used by MATMAN.
c
c    Output, integer MAXDIG, the maximum number of decimal digits to use.
c
c    Output, integer MAXINT, the maximum legal integer, which
c    should probably be 2147483647.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Output, integer NART, the number of artificial variables.
c
c    Output, integer NCOL, the number of columns in the matrix.
c
c    Output, integer NCON, the number of constraints.
c
c    Output, integer NDIG, the number of decimal digits used.
c
c    Output, integer NLINE, the number of characters of input in LINE.
c
c    Output, integer NROW, the number of rows in the matrix.
c
c    Output, integer NSLAK, the number of slack variables.
c
c    Output, integer NVAR, the number of basic variables.
c
c    Output, real SOL(MAXROW), the current linear programming solution.
c
      integer maxcol
      integer maxrow
c
      real a(maxrow,maxcol)
      logical autop
      character*1 chineq(maxrow)
      character*20 comnew
      character*20 comold
      real dete
      character*60 filex
      character*60 filhlp
      character*60 filinp
      character*60 filkey
      character*60 filtrn
      integer i
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer iauthr
      integer ibase(maxrow)
      integer idebot
      integer idetop
      integer ierror
      integer iform
      integer imat
      integer iounit(4)
      integer iprint
      integer iseed
      integer islbot(maxcol)
      integer isltop(maxcol)
      character*(*) line
      integer lpage
      integer lpmoda
      integer maxdig
      integer maxint
      integer nart
      integer ncol
      integer ncon
      integer ndig
      integer nline
      integer nrow
      integer nslak
      integer nvar
      real sol(maxcol)
c
      call inimat ( a, iabot, iatop, iform, maxcol, maxrow )
 
      autop = .true.
 
      do i = 1, maxrow
        chineq(i) = ' '
      end do
 
      comnew = ' '
      comold = ' '
 
      dete = 1.0
c
c  The file names, as given here, assume that the MATMAN files
c  are in the directory where the user is working.
c
c  If MATMAN is installed on a computer in a special directory,
c  but the user wishes to run it while working in another
c  directory, then the names of FILHLP and FILKEY must be
c  changed to include the directory information.
c
c  Similarly, if a single copy of MATMAN is installed on a
c  multi-user computer, then the file names for FILHLP and FILKEY
c  would need to be changed to include the directory information.
c
      filex = 'matman.dat'
      filhlp = 'matman.hlp'
      filinp = 'matman.inp'
      filkey = 'matkey.dat'
      filtrn = 'matman.lpt'
 
      iauthr = 0
 
      do i = 1, maxrow
        ibase(i) = i
      end do
 
      idebot = 1
      idetop = 1
 
      ierror = 0
      iform = 0
      imat = 0
 
      iounit(1) = 0
      iounit(2) = 0
      iounit(3) = -1
      iounit(4) = -1
 
      iprint = 0
 
      iseed = 10031952

      do i = 1, maxcol
        isltop(i) = 0
        islbot(i) = 1
      end do
 
      line = ' '
      lpage = 24
      call setpag ( lpage )
 
      lpmoda = 0
c
c  Use MAXDIG=7 for 32 bit IEEE machines (whether or not double
c  precision is used!)
c
      maxdig = 7
c
c  Use MAXDIG = 14 (approximately) for a Cray, or other machines
c  with 64 bit integers.
c
c     maxdig = 14
c
c  Use MAXINT = 2147483647 for 32 bit IEEE machines.
c
      maxint = 2147483647
c
c  Use MAXINT = 9223372036854775807 for 64 bit integer machines.
c
c     maxint = 9223372036854775807
c
      nart = 0
      ncol = 0
      ncon = 0
      ndig = 6
      nline = 0
      nrow = 0
      nslak = 0
      nvar = 0
 
      do i = 1, maxcol
        sol(i) = 0.0
      end do
 
      return
      end
      subroutine intran ( i, ilo, ihi, iseed )
c
c***********************************************************************
c
cc INTRAN returns a random integer in a given range.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer I, the randomly chosen integer.
c
c    Input, integer ILO, IHI, the minimum and maximum values acceptable
c    for I.
c
c    Input/output, integer ISEED, a seed for the random number generator.
c
      integer i
      integer ihi
      integer ilo
      integer iseed
      real r
      real random
      real rhi
      real rlo
c
c  Pick a random number in (0,1).
c
      r = random ( iseed )
c
c  Set a real interval [RLO,RHI] which contains the integers [ILO,IHI],
c  each with a "neighborhood" of width 1.
c
      rlo = real(ilo) - 0.5
      rhi = real(ihi) + 0.5
c
c  Set I to the integer that is nearest the scaled value of R.
c
      i = nint ( (1.0-r) * rlo + r * rhi ) 
c
c  In case of oddball events at the boundary, enforce the limits.
c
      i = max ( i, ilo )
      i = min ( i, ihi )

      return
      end
      subroutine intrea ( intval, line, nline, prompt, iounit, ierror )
c
c***********************************************************************
c
cc INTREA reads an integer from the input buffer.
c
c
c  INTREA accepts LINE which contains NLINE characters (NLINE may be
c  less than 1) and a PROMPT line.  If NLINE is less than 1, the
c  PROMPT will be printed and LINE read from the input unit,
c  IOUNIT(1), and NLINE will be updated.
c
c  In either case, the integer INTVAL will be read from LINE,
c  beginning at character 1 and ending at the first comma, slash,
c  blank, or the end of LINE.
c
c  The PROMPT should consist of a string of names of data items,
c  separated by commas, with the current one first.
c
c  The program will print 'ENTER' PROMPT and after reading LINE
c  will strip the characters corresponding to INTVAL from LINE,
c  and the part of PROMPT up to the first comma, leaving LINE and
c  PROMPT ready for another call to INTREA, CHRREA, RATREA or
c  RELREA.
c
c  If NLINE is greater than 0, but no characters can be read into
c  INTVAL, IERROR = 1 and we return
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer INTVAL, the integer that was read.
c
c    Workspace, character*80 LINE.
c    Used to hold the user's input.
c
c    Input/output, integer NLINE.
c    Keeps track of the number of useful characters in LINE.
c
c    Input, character*80 PROMPT, the prompt string.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
      integer ierror
      integer intval
      integer iounit(4)
      integer lchar
      integer lenchr
      character*80 line
      integer nline
      character*100 output
      character*80 prompt
c
      intval = 0
c
c  Retrieve a likely character string from input.
c
10    continue

      call chrinp ( ierror, iounit, line, nline, output, prompt )
      if ( ierror .ne. 0 ) return

      if ( nline .le. 0 ) then
        go to 10
      end if
c
c  Convert the character string to an integer.
c
      call chrcti ( line, intval, ierror, lchar )
      if ( ierror .ne. 0 ) return
c
c  Remove the character string from the input line.
c
      call chrchp ( line, 1, lchar )
      nline = lenchr ( line )
 
      return
      end
      subroutine lahlp1 ( iounit, output )
c
c***********************************************************************
c
cc LAHLP1 prints out a brief list of useful linear algebra commands.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Workspace, character*100 OUTPUT.
c
      integer iounit(4)
      character*100 output
c
      output = 'C    I,J,S changes matrix entry I, J to S.'
      call chrwrt ( iounit, output )
      output = 'E    enters a matrix to work on.'
      call chrwrt ( iounit, output )
      output = 'HELP for full help.'
      call chrwrt ( iounit, output )
      output = 'L    switches to linear programming.'
      call chrwrt ( iounit, output )
      output = 'O    checks if the matrix is row reduced.'
      call chrwrt ( iounit, output )
      output = 'Q    quits.'
      call chrwrt ( iounit, output )
      output = 'Z    automatic row reduction (requires password).'
      call chrwrt ( iounit, output )
      output = '?    for interactive help.'
      call chrwrt ( iounit, output )
      output = ' '
      call chrwrt ( iounit, output )
      output = 'R1 <=> R2  interchanges two rows'
      call chrwrt ( iounit, output )
      output = 'R1 <= S R1 multiplies a row by S.'
      call chrwrt ( iounit, output )
      output = 'R1 <= R1 + S R2 adds a multiple of another row.'
      call chrwrt ( iounit, output )
      output = ' '
      call chrwrt ( iounit, output )
 
      return
      end
      subroutine lainp0 ( a, iatop, iabot, ierror, iform, iounit, line,
     &  maxcol, maxrow, ncol, nline, nrow, nvar, output, prompt )
c
c***********************************************************************
c
cc LAINP0 begins the process of receiving a matrix from the user.
c
c
c  It finds out the dimensions of the array, and zeroes it out.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, real A(MAXROW,MAXCOL).  A is the current matrix.
c
c    Output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
c    IATOP and IABOT represent the current rational or decimal
c    matrix.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IFORM, specifies the arithmetic being used.
c    0=rational, 1=real, 2=decimal.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Workspace, character*80 LINE.
c    Used to hold the user's input.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Output, integer NCOL, the number of columns in the matrix.
c
c    Input/output, integer NLINE.
c    Keeps track of the number of useful characters in LINE.
c
c    Output, integer NROW, the number of rows in the matrix.
c
c    Output, integer NVAR, the number of basic variables.
c
c    Workspace, character*100 OUTPUT.
c
c    Workspace, character*80 PROMPT.
c
      integer maxcol
      integer maxrow
c
      real a(maxrow,maxcol)
      character*6 chrint
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer ierror
      integer iform
      integer iounit(4)
      character*80 line
      integer ncol
      integer nline
      integer nrow
      integer nvar
      character*100 output
      character*80 prompt
c
      nrow = 0
      ncol = 0
      nvar = 0
 
      prompt = 'number of rows, number of columns.'
c
c  Get NROW, the number of rows.
c
      call intrea ( nrow, line, nline, prompt, iounit, ierror )
      if ( ierror .ne. 0 ) return
 
      if ( nrow .lt. 1 ) then
        output = 'Error!  Negative number of rows not allowed!'
        call chrwrt ( iounit, output )
        ierror = 1
        return
      else if ( nrow .gt. maxrow ) then
        output = 'Number of rows must be less than '//chrint(maxrow)
        call chrdb2 ( output )
        ierror = 1
        return
      end if
c
c  Get NCOL, the number of columns.
c
      call intrea ( ncol, line, nline, prompt, iounit, ierror )
      if ( ierror .ne. 0 ) return
 
      if ( ncol .lt. 1 ) then
        output = 'Error!  Negative number of columns not allowed!'
        call chrwrt ( iounit, output )
        ierror = 1
        return
      else if ( ncol .gt. maxcol ) then
        output = 'Number of columns must be less than '
     &    // chrint(maxcol)
        call chrdb2 ( output )
        ierror = 1
        return
      end if
c
c  Zero out the matrix.
c
      call inimat ( a, iabot, iatop, iform, maxcol, maxrow )
 
      return
      end
      subroutine lainp1 ( a, iabot, iatop, icol, ierror, iform, iounit,
     &  irow, line, maxcol, maxdig, maxrow, ncol, ndig, nline, nrow,
     &  output, prompt )
c
c***********************************************************************
c
cc LAINP1 accepts the values of the entries of a matrix from the user.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, real A(MAXROW,MAXCOL).  A is the current matrix.
c
c    Output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
c    IATOP and IABOT represent the current rational or decimal
c    matrix.
c
c    Input, integer ICOL.
c    0, enter a single row of the matrix.
c    1, enter all rows of the matrix.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IFORM, specifies the arithmetic being used.
c    0=rational, 1=real, 2=decimal.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Input, integer IROW.
c    0, enter a single column of the matrix.
c    1, enter all columns of the matrix.
c
c    Workspace, character*80 LINE.
c    Used to hold the user's input.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXDIG, the maximum number of decimal digits
c    to use.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input, integer NCOL, the number of columns in the matrix.
c
c    Input, integer NDIG, the number of decimal digits to use.
c
c    Input/output, integer NLINE.
c         Keeps track of the number of useful characters in LINE.
c
c    Input, integer NROW, the number of rows in the matrix.
c
c    Workspace, character*100 OUTPUT.
c
c    Workspace, character*80 PROMPT.
c
      integer maxcol
      integer maxrow
c
      real a(maxrow,maxcol)
      character*6 chrint
      integer i
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer ibot
      integer icol
      integer ierror
      integer iform
      integer iounit(4)
      integer irow
      integer itop
      integer j
      character*80 line
      integer maxdig
      integer ncol
      integer ndig
      integer nline
      integer nrow
      character*100 output
      character*80 prompt
      real rval
c
      if ( irow .eq. 0 .and. icol.eq.0 ) then
        output = 'LAINP1 - Programming error!'
        call chrwrt ( iounit, output )
        output = 'IROW=ICOL = 0'
        call chrwrt ( iounit, output )
        ierror = 1
        return
      end if
c
c  Enter a single row.
c
      if ( icol .eq. 0 ) then
 
c       nline = 0
 
        do j = 1, ncol
 
          prompt = 'entries '//chrint(j)//' to '//chrint(ncol)//
     &      ' of row '//chrint(irow)
          call chrdb2 ( prompt )
 
          if ( iform .eq. 0 ) then
 
            call ratrea ( itop, ibot, rval, line, nline, prompt, iounit,
     &        ierror )
            if ( ierror .ne. 0 ) return
            iatop(irow,j) = itop
            iabot(irow,j) = ibot
 
          else if ( iform .eq. 1 ) then
 
            call relrea ( rval, line, nline, prompt, iounit, ierror )
 
            if ( ierror .ne. 0 ) return
            a(irow,j) = rval
 
          else if ( iform .eq. 2 ) then
 
            call decrea ( itop, ibot, rval, line, maxdig, nline, prompt,
     &        iounit, ierror )
 
            if ( ierror .ne. 0 ) return
 
            call deccut ( itop, ibot, ndig )
 
            iatop(irow,j) = itop
            iabot(irow,j) = ibot
 
          end if
 
        end do
c
c  Enter a single column.
c
      else if ( irow .eq. 0 ) then
 
c       nline = 0
 
        do i = 1, nrow
 
          prompt = 'entries '//chrint(i)//' to '//chrint(nrow)//
     &      ' of column ' // chrint(icol)
          call chrdb2 ( prompt )
 
          if ( iform .eq. 0 ) then
 
            call ratrea ( itop, ibot, rval, line, nline, prompt, iounit,
     &        ierror )
            if ( ierror .ne. 0 ) return
            iatop(i,icol) = itop
            iabot(i,icol) = ibot
 
          else if ( iform .eq. 1 ) then
 
            call relrea ( rval, line, nline, prompt, iounit, ierror )
            if ( ierror .ne. 0 ) return
            a(i,icol) = rval
 
          else if ( iform .eq. 2 ) then
 
            call decrea ( itop, ibot, rval, line, maxdig, nline, prompt,
     &        iounit, ierror )
 
            if ( ierror .ne. 0 ) return
 
            call deccut ( itop, ibot, ndig )
 
            iatop(i,icol) = itop
            iabot(i,icol) = ibot
 
          end if
 
        end do
c
c  Enter an entire matrix.
c
      else
 
c       nline = 0
 
        do i = 1, nrow
          do j = 1, ncol
 
            prompt = 'entries '//chrint(j)//' to '//chrint(ncol)//
     &        ' of row '//chrint(i)
            call chrdb2 ( prompt )
 
            if ( iform .eq. 0 ) then
 
              call ratrea ( itop, ibot, rval, line, nline, prompt,
     &          iounit, ierror )
              if ( ierror .ne. 0 ) return
              iatop(i,j) = itop
              iabot(i,j) = ibot
 
            else if ( iform .eq. 1 ) then
 
              call relrea ( rval, line, nline, prompt, iounit, ierror )
              if ( ierror .ne. 0 ) return
              a(i,j) = rval
 
            else if ( iform .eq. 2 ) then
 
              call decrea ( itop, ibot, rval, line, maxdig, nline,
     &          prompt, iounit, ierror )
 
              if ( ierror .ne. 0 ) return
 
              call deccut ( itop, ibot, ndig )
 
              iatop(i,j) = itop
              iabot(i,j) = ibot
 
            end if
 
          end do
        end do
 
      end if
 
      return
      end
      subroutine laopt ( a, iabot, iatop, ierror, iform, imat, iounit,
     &  maxcol, maxrow, ncol, nrow, output )
c
c***********************************************************************
c
cc LAOPT checks for row echelon or reduced row echelon form.
c
c
c  A matrix is in row echelon form if:
c
c    1.  The first nonzero entry in each row is 1.
c
c    2.  The leading 1 in a given row occurs in a column to
c        the right of the leading 1 in the previous row.
c
c    3.  Rows which are entirely zero must occur last.
c
c  The matrix is in reduced row echelon form if, in addition to
c  the first three conditions, it also satisfies:
c
c    4.  Each column containing a leading 1 has no other nonzero
c        entries.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real A(MAXROW,MAXCOL).  A is the current matrix.
c
c    Input, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
c    IATOP and IABOT represent the current rational or decimal
c    matrix.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IFORM, specifies the arithmetic being used.
c    0=rational, 1=real, 2=decimal.
c
c    Input, integer IMAT.
c    0, no matrix has been defined by the user.
c    1, a matrix has been defined by the user.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input, integer NCOL, the number of columns in the matrix.
c
c    Input, integer NROW, the number of rows in the matrix.
c
c    Workspace, character*100 OUTPUT.
c
      integer maxcol
      integer maxrow
c
      real a(maxrow,maxcol)
      character*22 chldec
      character*22 chlrat
      character*6 chrint
      character*14 chrrel
      character*22 chrtmp
      integer i
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer ierror
      integer iform
      integer ii
      integer imat
      integer iounit(4)
      integer izer
      integer j
      integer lead
      integer leadp
      integer ncol
      integer nrow
      character*100 output
c
      if ( imat .eq. 0 ) then
        ierror = 1
        output = 'You must set up a matrix first!'
        call chrwrt ( iounit, output )
        return
      end if
c
c  Check rule 1.
c
      do i = 1, nrow
        do j = 1, ncol
 
          if ( iform .eq. 0 ) then
 
            if ( iatop(i,j) .eq. 0 ) then
              go to 10
            else if ( iatop(i,j) .eq. iabot(i,j) ) then
              go to 20
            end if
 
          else if ( iform .eq. 1 ) then
 
            if ( a(i,j) .eq. 0 ) then
              go to 10
            else if ( a(i,j) .eq. 1 ) then
              go to 20
            end if
 
          else if ( iform .eq. 2 ) then
 
            if ( iatop(i,j) .eq. 0 ) then
              go to 10
            else if ( iatop(i,j) .eq. 1 .and. iabot(i,j).eq.0 ) then
              go to 20
            end if
          end if
 
          output = 'This matrix is NOT in row echelon form.'
          call chrwrt ( iounit, output )
 
          output = 'The first nonzero entry in row ' // chrint(i)
          call chrdb2 ( output )
          call chrwrt ( iounit, output )
 
          output = 'which occurs in column '//chrint(j)
          call chrdb2 ( output )
          call chrwrt ( iounit, output )
 
          if ( iform .eq. 0 ) then
            chrtmp = chlrat ( iatop(i,j), iabot(i,j) )
          else if ( iform .eq. 1 ) then
            chrtmp = chrrel ( a(i,j) )
          else if ( iform .eq. 2 ) then
            chrtmp = chldec ( iatop(i,j), iabot(i,j) )
          end if
 
          output = 'is ' // chrtmp // ' rather than 1.'
          call chrdb2 ( output )
          call chrwrt ( iounit, output )
          return
 
10        continue
 
          end do
 
20      continue
 
        end do
c
c  Check rule 2.
c
      lead = 0
 
      do i = 1, nrow
        do j = 1, ncol
 
          if ( iform .eq. 0 ) then
 
            if ( iatop(i,j) .eq. 0 ) then
              go to 30
            else if ( iatop(i,j) .eq. iabot(i,j) ) then
              leadp = lead
              lead = j
              if ( lead .gt. leadp ) go to 40
            end if
 
          else if ( iform .eq. 1 ) then
 
            if ( a(i,j) .eq. 0 ) then
              go to 30
            else if ( a(i,j) .eq. 1 ) then
              leadp = lead
              lead = j
              if ( lead .gt. leadp ) go to 40
            end if
 
          else if ( iform .eq. 2 ) then
 
            if ( iatop(i,j) .eq. 0 ) then
              go to 30
            else if ( iatop(i,j) .eq. 1 .and. iabot(i,j).eq.0 ) then
              leadp = lead
              lead = j
              if ( lead .gt. leadp ) go to 40
            end if
 
          end if
 
          output = 'This matrix is NOT in row echelon form.'
          call chrwrt ( iounit, output )
          output = 'The first 1 in row ' // chrint(i) // ' does NOT'
          call chrdb2 ( output )
          call chrwrt ( iounit, output )
          output = 'occur to the right of the first 1 in row'
          call chrwrt ( iounit, output )
          output = chrint(i-1)
          call chrdb2 ( output )
          call chrwrt ( iounit, output )
          return
 
30      continue
 
        end do
 
40    continue
 
      end do
c
c  Check rule 3.
c
      izer = 0
 
      do i = 1, nrow
 
        if ( izer .eq. 0 ) then
 
          do j = 1, ncol
 
            if ( iform .eq. 0 ) then
              if ( iatop(i,j) .ne. 0 ) go to 70
            else if ( iform .eq. 1 ) then
              if ( a(i,j) .ne. 0 ) go to 70
            else if ( iform .eq. 2 ) then
              if ( iatop(i,j) .ne. 0 ) go to 70
            end if
 
          end do
 
          izer = i
 
        else
 
          do j = 1, ncol
 
            if ( 
     &        ( iform .eq. 0 .and. iatop(i,j) .ne. 0 ) .or.
     &        ( iform .eq. 1 .and. a(i,j) .ne. 0 ) .or.
     &        ( iform .eq. 2 .and. iatop(i,j) .ne. 0 ) ) then
 
              output = 'This matrix is NOT in row echelon form.'
              call chrwrt ( iounit, output )
              output = 'Row '//chrint(izer)//' is entirely zero.'
              call chrdb2 ( output )
              call chrwrt ( iounit, output )
              output = 'Row '//chrint(i)//' occurs afterwards, and has'
              call chrdb2 ( output )
              call chrwrt ( iounit, output )
              output = 'nonzero entries in it!'
              call chrwrt ( iounit, output )
              return

            end if
 
          end do
 
        end if
 
70      continue
 
      end do
 
      output = 'This matrix is in row echelon form.'
      call chrwrt ( iounit, output )
c
c  Check rule 4.
c
      do i = 1, nrow
 
        do j = 1, ncol
c
c  We know first nonzero in this row will be 1.
c
          if ( iform .eq. 0 ) then
            if ( iatop(i,j) .eq. 0 ) go to 90
          else if ( iform .eq. 1 ) then
            if ( a(i,j) .eq. 0 ) go to 90
          else if ( iform .eq. 2 ) then
            if ( iatop(i,j) .eq. 0 ) go to 90
          end if
c
c  The leading 1 of this row is entry (i,j).
c
          do ii = 1, nrow
 
            if ( ii .ne. i ) then
 
              if ( iform .eq. 0 ) then
                if ( iatop(ii,j) .eq. 0 ) go to 80
              else if ( iform .eq. 1 ) then
                if ( a(ii,j) .eq. 0 ) go to 80
              else if ( iform .eq. 2 ) then
                if ( iatop(ii,j) .eq. 0 ) go to 80
              end if
 
              output = ' '
              call chrwrt ( iounit, output )
              output = 'This matrix is NOT in reduced row echelon form.'
              call chrwrt ( iounit, output )
              output = 'Row '//chrint(i)//' has its leading 1 in '//
     &          'column '//chrint(j)
              call chrdb2 ( output )
              call chrwrt ( iounit, output )
              output = 'This means that all other entries of that '//
     &          'column should be zero.'
              call chrwrt ( iounit, output )
 
              if ( iform .eq. 0 ) then
                chrtmp = chlrat ( iatop(ii,j), iabot(ii,j) )
                output = 'But the entry in row '//chrint(ii)//' is '//
     &            chrtmp
              else if ( iform .eq. 1 ) then
                output = 'But the entry in row '//chrint(ii)//' is '//
     &            chrrel(a(ii,j))
              else if ( iform .eq. 2 ) then
                chrtmp = chldec ( iatop(ii,j), iabot(ii,j) )
                output = 'But the entry in row '//chrint(ii)//' is '//
     &            chrtmp
              end if
 
              call chrdb2 ( output )
              call chrwrt ( iounit, output )
              return
            end if
 
80          continue
 
          end do
 
          go to 100
 
90        continue
 
        end do
 
100     continue

      end do
 
      output = ' '
      call chrwrt ( iounit, output )
      output = 'In fact, this matrix is in reduced row echelon form.'
      call chrwrt ( iounit, output )
 
      return
      end
      function ldigit ( string )
c
c***********************************************************************
c
cc LDIGIT returns .TRUE. if STRING contains only digits or blanks.
c
c
c  Note that LDIGIT MUST be declared logical in your calling program.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) STRING, the string to be checked.
c
c    Output, logical LDIGIT, .TRUE. if STRING contains only digits and
c    blanks, .FALSE. otherwise.
c
      character*1 chrtmp
      integer i
      logical ldigit
      integer lenc
      character*(*) string
c
      lenc = len ( string )
 
      ldigit = .false.
 
      do i = 1, lenc
 
        chrtmp = string(i:i)
 
        if ( chrtmp .ne. ' ' ) then
          if ( llt ( chrtmp, '0' ) .or. lgt ( chrtmp, '9' ) ) then
            return
          end if
        end if
 
      end do
 
      ldigit = .true.

      return
      end
      function lenchr ( string )
c
c***********************************************************************
c
cc LENCHR returns the length of STRING up to the last nonblank.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) STRING, the string to be measured.
c
c    Output, integer LENCHR, the location of the last nonblank
c    character in STRING.
c
      integer i
      integer lenchr
      character*(*) string
c
      do i = len ( string ), 1, -1
 
        if ( string(i:i) .ne. ' ' ) then
          lenchr = i
          return
        end if
 
      end do
 
      lenchr = 0
 
      return
      end
      function leqi ( strng1, strng2 )
c
c***********************************************************************
c
cc LEQI is a case insensitive comparison of two strings for equality.  
c
c
c  Thus, LEQI('Anjana','ANJANA') is .TRUE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) STRNG1, STRNG2, the strings to compare.
c
c    Output, logical LEQI, the result of the comparison.
c
      integer i
      integer len1
      integer len2
      integer lenc
      logical leqi
      character*1 s1
      character*1 s2
      character*(*) strng1
      character*(*) strng2
c
      len1 = len ( strng1 )
      len2 = len ( strng2 )
      lenc = min ( len1, len2 )
 
      leqi = .false.

      do i = 1, lenc
        s1 = strng1(i:i)
        s2 = strng2(i:i)
        call capchr ( s1 )
        call capchr ( s2 )
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
 
      leqi = .true.
 
      return
      end
      subroutine lphlp1 ( iounit, output )
c
c***********************************************************************
c
cc LPHLP1 prints out a brief list of useful linear programming commands.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Workspace, character*100 OUTPUT.
c
      integer iounit(4)
      character*100 output

      output = 'C     I, J, S changes tableau entry I, J to S.'
      call chrwrt ( iounit, output )
      output = 'E     Enters a tableau to work on.'
      call chrwrt ( iounit, output )
      output = 'HELP  for full help.'
      call chrwrt ( iounit, output )
      output = 'L     switches to linear algebra.'
      call chrwrt ( iounit, output )
      output = 'O     checks if the solution is optimal.'
      call chrwrt ( iounit, output )
      output = 'P     I, J performs a pivot operation.'
      call chrwrt ( iounit, output )
      output = 'Q     quits.'
      call chrwrt ( iounit, output )
      output = 'TS    types the linear programming solution.'
      call chrwrt ( iounit, output )
      output = 'V     removes artificial variables.'
      call chrwrt ( iounit, output )
      output = 'Z     automatic solution (requires password).'
      call chrwrt ( iounit, output )
      output = '?     interactive help.'
      call chrwrt ( iounit, output )
      output = ' '
      call chrwrt ( iounit, output )
      output = 'R1 <=> R2  interchanges two rows'
      call chrwrt ( iounit, output )
      output = 'R1 <= S R1 multiplies a row by S.'
      call chrwrt ( iounit, output )
      output = 'R1 <= R1 + S R2 adds a multiple of another row.'
      call chrwrt ( iounit, output )
      output = ' '
      call chrwrt ( iounit, output )
 
      return
      end
      subroutine lpinp ( a, chineq, iatop, iabot, ibase, ierror,
     &  iform, iounit, line, maxcol, maxdig, maxrow, nart, ncol, ncon,
     &  ndig, nline, nrow, nslak, nvar, output, prompt )
c
c***********************************************************************
c
cc LPINP allows the user to enter a linear programming problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, real A(MAXROW,MAXCOL), the current tableau.
c
c    Output, character*1 CHINEQ(MAXROW), the '<', '=', or '>'
c    sign for each linear programming constraint.
c
c    Input, integer IABOT(MAXROW,MAXCOL), IATOP(MAXROW,MAXCOL).
c    IATOP and IABOT represent the current rational or decimal
c    matrix.
c
c    Output, integer IBASE(MAXROW), keeps track of basic variables.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IFORM, specifies the arithmetic being used.
c    0=rational, 1=real, 2=decimal.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Workspace, character*80 LINE.
c    Used to hold the user's input.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXDIG, the maximum number of decimal digits
c    to use.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Output, integer NART, the number of artificial variables.
c
c    Output, integer NCOL, the number of columns in the matrix.
c
c    Output, integer NCON, the number of constraints.
c
c    Input, integer NDIG, the number of decimal digits in use.
c
c    Input/output, integer NLINE.
c     Keeps track of the number of useful characters in LINE.
c
c    Output, integer NROW, the number of rows in the matrix.
c
c    Output, integer NSLAK, the number of slack variables.
c
c    Output, integer NVAR, the number of basic variables.
c
c    Workspace, character*100 OUTPUT.
c
c    Workspace, character*80 PROMPT.
c
      integer maxcol
      integer maxrow

      real a(maxrow,maxcol)
      character*1 chineq(maxrow)
      character*6 chrint
      integer i
      integer iabot(maxrow,maxcol)
      integer iart
      integer iatop(maxrow,maxcol)
      integer ibase(maxrow)
      integer ibot
      integer ierror
      integer iform
      integer iounit(4)
      character*1 isay
      integer islak
      integer iterm
      integer itop
      integer j
      integer j1
      integer j2
      integer jcol
      integer jhi
      character*80 line
      integer maxdig
      integer nart
      integer ncol
      integer ncon
      integer ndig
      integer nline
      integer nrow
      integer nslak
      integer nvar
      character*100 output
      character*80 prompt
      real rval
c
c  Zero out the matrix.
c
      call inimat ( a, iabot, iatop, iform, maxcol, maxrow )

      nrow = 0
      ncol = 0
      ncon = 0
      nvar = 0
      nslak = 0
      nart = 0
c
c  Read two integers defining problem.
c
      prompt = 'number of constraints, number of variables.'
c
c  Get number of constraints.
c
      call intrea ( ncon, line, nline, prompt, iounit, ierror )
      if ( ierror .ne. 0 ) return
 
      if ( ncon .lt. 0 ) then
        output = 'Number of constraints must be positive!'
        call chrwrt ( iounit, output )
        ierror = 1
        return
      else if ( ncon .gt. maxrow-2 ) then
        output = 'Number of constraints must be no greater than '//
     &    chrint(maxrow-2)
        call chrdb2 ( output )
        call chrwrt ( iounit, output )
        ierror = 1
        return
      end if
 
      nrow = ncon+1
c
c  Get number of variables.
c
      call intrea ( nvar, line, nline, prompt, iounit, ierror )
      if ( ierror .ne. 0 ) return
 
      if ( nvar .lt. 1 ) then
        output = 'A negative number of variables is not allowed!'
        call chrwrt ( iounit, output )
        ierror = 1
        return
      else if ( nvar .gt. maxcol ) then
        output = 'Number of variables must be no greater than '//
     &    chrint(maxcol)
        call chrwrt ( iounit, output )
        ierror = 1
        return
      end if
c
c  Zero out the matrix.
c
      call inimat ( a, iabot, iatop, iform, maxcol, maxrow )
 
      nline = 0
 
      do i = 1, nrow
 
        chineq(i) = ' '
 
        if ( i .le. ncon ) then
          nline = 0
30        continue
          prompt = 'sign < > or = and coefficients and RHS of '//
     &      'constraint'//chrint(i)
          call chrdb2 ( prompt )
          iterm = 0
          call chrrea ( isay, line, nline, prompt, iounit, ierror,
     &      iterm ) 
          if ( ierror .ne. 0 ) return
 
          if ( isay .eq. '<' ) then
            ibase(i) = - 1
            nslak = nslak+1
          else if ( isay .eq. '=' ) then
            ibase(i) = 1
            nart=nart+1
          else if ( isay .eq. '>' ) then
            ibase(i) = 0
            nslak = nslak+1
            nart = nart+1
          else
            output = 'Huh?  Try again.'
            go to 30
          end if
 
          chineq(i) = isay
        else
          nline = 0
          prompt = 'coefficients and constant of objective function.'
        end if
 
        jhi = nvar+1
 
        do j = 1, jhi
 
          jcol = j
          if ( j .eq. jhi ) then
            jcol = maxcol
          end if

          if ( i .le. ncon ) then
 
            if ( j .lt. jhi ) then
              prompt = 'entries '//chrint(j)//' to '//chrint(nvar)//
     &          ' and RHS of constraint '//chrint(i)
              call chrdb2 ( prompt )
            else
              prompt = 'right hand side of constraint '//chrint(i)
            end if
 
          end if
 
          call chrdb2 ( prompt )
 
          if ( iform .eq. 0 ) then
 
            call ratrea ( itop, ibot, rval, line, nline, prompt, iounit,
     &        ierror )
            if ( ierror .ne. 0 ) return

            iatop(i,jcol) = itop
            iabot(i,jcol) = ibot
            if ( (i .eq. nrow) .and. (jcol .le. nvar) ) then
              iatop(i,jcol) = - iatop(i,jcol)
            end if
 
          else if ( iform .eq. 1 ) then
 
            call relrea ( rval, line, nline, prompt, iounit, ierror )
            if ( ierror .ne. 0 ) return

            a(i,jcol) = rval

            if ( i .eq. nrow .and. jcol .le. nvar ) then
              a(i,jcol) = -a(i,jcol)
            end if
 
          else if ( iform .eq. 2 ) then
 
            call decrea ( itop, ibot, rval, line, maxdig, nline, prompt,
     &        iounit, ierror )
 
            if ( ierror .ne. 0 ) return
 
            call deccut ( itop, ibot, ndig )
 
            iatop(i,jcol) = itop
            iabot(i,jcol) = ibot
            if ( i .eq. nrow .and. jcol .le. nvar ) then
              iatop(i,jcol) = - iatop(i,jcol)
            end if
 
          end if
 
        end do
 
      end do
c
c  Move the right hand sides to the proper column.
c
      j2 = nvar + nslak + nart + 2

      do i = 1, nrow
 
        if ( iform .eq. 0 ) then
          iatop(i,j2) = iatop(i,maxcol)
          iabot(i,j2) = iabot(i,maxcol)
        else if ( iform .eq. 1 ) then
          a(i,j2) = a(i,maxcol)
        else if ( iform .eq. 2 ) then
          iatop(i,j2) = iatop(i,maxcol)
          iabot(i,j2) = iabot(i,maxcol)
        end if
 
      end do
c
c  Place the 1 in the bottom of the P column.
c
      j1 = nvar + nslak + nart + 1

      if ( iform .eq. 0 ) then
        iatop(nrow,j1) = 1
        iabot(nrow,j1) = 1
      else if ( iform .eq. 1 ) then
        a(nrow,j1) = 1.0
      else if ( iform .eq. 2 ) then
        iatop(nrow,j1) = 1
        iabot(nrow,j1) = 0
      end if
c
c  For artificial variable problems, move the objective row down
c  one row to a "hidden" row, and set up a dummy objective row.
c
      if ( nart .gt. 0 ) then
        output = ' '
        call chrwrt ( iounit, output )
        output = 'Because we have artificial variables, the objective'
        call chrwrt ( iounit, output )
        output = 'function is also "artificial".'
        call chrwrt ( iounit, output )
        output = 'The true objective will be stored away until the'
        call chrwrt ( iounit, output )
        output = 'artificial variables are gone.'
        call chrwrt ( iounit, output )
        output = ' '
        call chrwrt ( iounit, output )
 
        do j = 1, nvar+nslak
 
          if ( iform .eq. 0 ) then
            iatop(nrow+1,j) = iatop(nrow,j)
            iabot(nrow+1,j) = iabot(nrow,j)
            iatop(nrow,j) = 0
            iabot(nrow,j) = 1
          else if ( iform .eq. 1 ) then
            a(nrow+1,j) = a(nrow,j)
            a(nrow,j) = 0.0
          else if ( iform .eq. 2 ) then
            iatop(nrow+1,j) = iatop(nrow,j)
            iabot(nrow+1,j) = iabot(nrow,j)
            iatop(nrow,j) = 0
            iabot(nrow,j) = 0
          end if
 
        end do
c
c  Move the last two entries of the original row to where they
c  would properly line up in the full problem.
c
        j1 = nvar+nslak+nart+1
        j2 = nvar+nslak+nart+2

        if ( iform .eq. 0 ) then
 
          iatop(nrow+1,j1) = 1
          iabot(nrow+1,j1) = 1
 
          iatop(nrow+1,j2) = iatop(nrow,j2)
          iabot(nrow+1,j2) = iabot(nrow,j2)
 
          iatop(nrow,j2) = 0
          iabot(nrow,j2) = 1
 
        else if ( iform .eq. 1 ) then
 
          a(nrow+1,j1) = 1.0
          a(nrow+1,j2) = a(nrow,j2)
          a(nrow,j2) = 0.0
 
        else if ( iform .eq. 2 ) then
 
          iatop(nrow+1,j1) = 1
          iabot(nrow+1,j1) = 0
 
          iatop(nrow+1,j2) = iatop(nrow,j2)
          iabot(nrow+1,j2) = iabot(nrow,j2)
 
          iatop(nrow,j2) = 0
          iabot(nrow,j2) = 0
 
        end if
 
      end if
c
c  Set entries corresponding to slack and artificial variables.
c
      islak = 0
      iart = 0
      ncol = nvar+nslak+nart+2
 
      do i = 1, ncon
 
        if ( ibase(i) .eq. -1 ) then
          islak = islak+1
          ibase(i) = nvar+islak
 
          if ( iform .eq. 0 ) then
            iatop(i,nvar+islak) = 1
            iabot(i,nvar+islak) = 1
          else if ( iform .eq. 1 ) then
            a(i,nvar+islak) = 1.0
          else if ( iform .eq. 2 ) then
            iatop(i,nvar+islak) = 1
            iabot(i,nvar+islak) = 0
          end if
 
        else if ( ibase(i) .eq. 0 ) then
 
          islak = islak+1
          iart = iart+1
          j = nvar+nslak+iart
          ibase(i) = j
 
          if ( iform .eq. 0 ) then
            iatop(i,nvar+islak) = - 1
            iabot(i,nvar+islak) = 1
            iatop(i,j) = 1
            iabot(i,j) = 1
            iatop(nrow,j) = 1
            iabot(nrow,j) = 1
          else if ( iform .eq. 1 ) then
            a(i,nvar+islak) = - 1.0
            a(i,j) = 1.0
            a(nrow,j) = 1.0
          else if ( iform .eq. 2 ) then
            iatop(i,nvar+islak) = - 1
            iabot(i,nvar+islak) = 0
            iatop(i,j) = 1
            iabot(i,j) = 0
            iatop(nrow,j) = 1
            iabot(nrow,j) = 0
          end if
 
        else if ( ibase(i) .eq. 1 ) then
 
          iart = iart+1
          j = nvar+nslak+iart
          ibase(i) = j
 
          if ( iform .eq. 0 ) then
            iatop(i,j) = 1
            iabot(i,j) = 1
            iatop(nrow,j) = 1
            iabot(nrow,j) = 1
          else if ( iform .eq. 1 ) then
            a(i,j) = 1.0
            a(nrow,j) = 1.0
          else if ( iform .eq. 2 ) then
            iatop(i,j) = 1
            iabot(i,j) = 0
            iatop(nrow,j) = 1
            iabot(nrow,j) = 0
          end if
 
        end if
 
      end do
 
      return
      end
      subroutine lpopt ( a, iatop, iabot, ibase, ierror, iform,
     &  imat, iopti, iounit, isltop, islbot, lpmoda, maxcol,
     &  maxrow, nart, ncol, nrow, nslak, nvar, output, sol )
c
c*******************************************************************
c
cc LPOPT checks the current linear programming tableau for optimality.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real A(MAXROW,MAXCOL).  A is the current matrix.
c
c    Input, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
c    IATOP and IABOT represent the current rational or decimal matrix.
c
c    Input, integer IBASE(MAXROW), keeps track of the basic variables.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IFORM, specifies the arithmetic being used.
c    0=rational, 1=real, 2=decimal.
c
c    Input, integer IMAT.
c    0, no matrix has been defined by the user.
c    1, a matrix has been defined by the user.
c
c    Output, integer IOPTI.
c    0, the current solution is not optimal.
c    1, the current solution is optimal.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Output, integer ISLTOP(MAXROW), ISLBOT(MAXROW), the fractional
c    or decimal representation of the linear programming solution.
c
c    Input, integer LPMODA.
c    0, the program is in linear algebra mode.
c    1, the program is in linear programming mode.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input, integer NART, the number of artificial variables.
c
c    Input, integer NCOL, the number of columns in the matrix.
c
c    Input, integer NROW, the number of rows in the matrix.
c
c    Input, integer NSLAK, the number of slack variables.
c
c    Input, integer NVAR, the number of basic variables.
c
c    Workspace, character*100 OUTPUT.
c
c    Output, real SOL(MAXROW), the real representation of the
c    linear programming solution.
c
      integer maxcol
      integer maxrow
c
      real a(maxrow,maxcol)
      character*22 chldec
      character*22 chlrat
      character*6 chrint
      character*14 chrrel
      character*22 chrtmp
      integer i
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer ibase(maxrow)
      integer ierror
      integer iform
      integer ihi
      integer ilo
      integer imat
      integer iopti
      integer iounit(4)
      integer islbot(maxcol)
      integer isltop(maxcol)
      integer jhi
      integer jlo
      integer lpmoda
      integer nart
      integer ncol
      integer nrow
      integer nslak
      integer nvar
      character*100 output
      real sol(maxcol)
      real temp
      character*80 title

      if ( imat .eq. 0 ) then
        ierror = 1
        output = 'You must set up a tableau first!'
        call chrwrt ( iounit, output )
        return
      end if
 
      output = ' '
      call chrwrt ( iounit, output )
      output = 'Optimality test'
      call chrwrt ( iounit, output )
      output = ' '
      call chrwrt ( iounit, output )
      iopti = 1
      output = ' '
      call chrwrt ( iounit, output )
      output = 'Are all objective entries nonnegative?'
      call chrwrt ( iounit, output )
 
      do i = 1, nslak+nvar+nart
 
        if ( iform .eq. 0 ) then
          call ratrel ( temp, iatop(nrow,i), iabot(nrow,i) )
        else if ( iform .eq. 1 ) then
          temp = a(nrow,i)
        else if ( iform .eq. 2 ) then
          call decrel(temp,iatop(nrow,i),iabot(nrow,i))
        end if
 
        if ( temp .lt. 0 ) then
          output = 'Negative objective coefficient, entry '//chrint(i)
          call chrdb2 ( output )
          call chrwrt ( iounit, output )
          iopti = 0
        end if
 
      end do
 
      output = ' '
      call chrwrt ( iounit, output )
 
      if ( iopti .eq. 0 ) then
        output = 'The current solution is NOT optimal.'
        call chrwrt ( iounit, output )
      else
        output = 'Yes.  The current solution is optimal.'
        call chrwrt ( iounit, output )
      end if
c
c  Print the current linear programming solution.
c
      call lpsol ( a, iatop, iabot, ibase, iform, isltop, islbot,
     &  maxcol, maxrow, ncol, nrow, sol )
 
      title = 'The linear programming solution:'
 
      jhi = nvar + nslak + nart
      jlo = 1
      ilo = 1
      ihi = 1
 
      if ( iform .eq. 0 ) then

        call ratprn ( isltop, islbot, ibase, iounit, ihi, ilo, jhi,
     &    jlo, lpmoda, jhi, 1, ncol, nrow, output, title )

      else if ( iform .eq. 1 ) then

        call relprn ( sol, ibase, iounit, ihi, ilo, jhi, jlo, lpmoda,
     &    jhi, 1, ncol, nrow, output, title )

      else if ( iform .eq. 2 ) then

        call decprn ( isltop, islbot, ibase, iounit, ihi, ilo, jhi,
     &    jlo, lpmoda, jhi, 1, ncol, nrow, output, title )

      end if
 
      output = ' '
      call chrwrt ( iounit, output )
 
      if ( iform .eq. 0 ) then
        chrtmp = chlrat ( iatop(nrow,ncol), iabot(nrow,ncol) )
      else if ( iform .eq. 1 ) then
        chrtmp = chrrel ( a(nrow,ncol) )
      else if ( iform .eq. 2 ) then
        chrtmp = chldec ( iatop(nrow,ncol), iabot(nrow,ncol) )
      end if

      output = 'Objective = ' // chrtmp
      call chrdb2 ( output )
      call chrwrt ( iounit, output )
c
c  Warn user if artificial variables must be deleted.
c
      if ( nart .gt. 0 ) then
        output = ' '
        call chrwrt ( iounit, output )
        output = 'This problem has artificial variables.'
        call chrwrt ( iounit, output )
        output = 'Use the "V" command to remove them.'
        call chrwrt ( iounit, output )
      end if
 
      return
      end
      subroutine lppiv ( a, iatop, iabot, iauto, ibase, ierror, 
     &  iform, imat, iounit, isltop, islbot, line, lpmoda, maxcol,
     &  maxint, maxrow, nart, ncol, ndig, nline, nrow, nslak, nvar,
     &  output, prompt, sol )
c
c*******************************************************************
c
cc LPPIV carries out pivoting for a linear programming problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, real A(MAXROW,MAXCOL).  A is the current matrix.
c
c    Input/output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
c    IATOP and IABOT represent the current rational or decimal matrix.
c
c    Input, integer IAUTO.
c    0, automatic processing is not being carried out.
c    1, automatic processing is being carried out.
c
c    Input/output, integer IBASE(MAXROW), keeps track of basic variables.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IFORM, specifies the arithmetic being used.
c    0=rational, 1=real, 2=decimal.
c
c    Input, integer IMAT.
c    0, no matrix has been defined by the user.
c    1, a matrix has been defined by the user.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Output, integer ISLTOP(MAXROW), ISLBOT(MAXROW), the fractional
c    or decimal representation of the linear programming solution.
c
c    Workspace, character*80 LINE.
c    Used to hold the user's input.
c
c    Input, integer LPMODA.
c    0, the program is in linear algebra mode.
c    1, the program is in linear programming mode.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXINT, the maximum integer representable, which
c    should probably be 2147483647.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input, integer NART, the number of artificial variables.
c
c    Input, integer NCOL, the number of columns in the matrix.
c
c    Input, integer NDIG, the number of decimal digits used.
c
c    Input/output, integer NLINE.
c    Keeps track of the number of useful characters in LINE.
c
c    Input, integer NROW, the number of rows in the matrix.
c
c    Input, integer NSLAK, the number of slack variables.
c
c    Input, integer NVAR, the number of basic variables.
c
c    Workspace, character*100 OUTPUT.
c
c    Workspace, character*80 PROMPT.
c
c    Output, real SOL(MAXROW), the real representation of the
c    linear programming solution.
c
      integer maxcol
      integer maxrow

      real a(maxrow,maxcol)
      character*22 chldec
      character*22 chlrat
      character*14 chrrel
      character*22 chrtmp
      character*22 chrtmp2
      integer i
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer iauto
      integer ibase(maxrow)
      integer ierror
      integer iform
      integer ihi
      integer ilo
      integer imat
      integer iobbot
      integer iobtop
      integer iopti
      integer iounit(4)
      integer ipiv
      integer irow
      integer irow1
      integer irow2
      integer isbot
      integer islbot(maxcol)
      integer isltop(maxcol)
      integer istop
      integer j
      integer jhi
      integer jlo
      integer jpiv
      character*80 line
      integer lpmoda
      integer maxint
      integer nart
      integer ncol
      integer ndig
      integer nline
      integer nrow
      integer nslak
      integer nvar
      real objnew
      real objold
      character*100 output
      character*80 prompt
      real sol(maxcol)
      real sval
      real temp
      character*80 title

      if ( lpmoda .ne. 1 ) then
        ierror = 1
        output = 'This command should only be given during'
        call chrwrt ( iounit, output )
        output = 'linear programming!'
        call chrwrt ( iounit, output )
        return
      end if
 
      if ( imat .eq. 0 ) then
        ierror = 1
        output = 'You must set up a tableau first!'
        call chrwrt ( iounit, output )
        return
      end if
c
c  For each basic variable, check that the objective row entry is zero.
c  If not, then if notautomatic, complain, else fix it.
c
10    continue
 
      call lppiv1 ( a, iabot, iatop, iauto, ibase, ierror, iform,
     &  iounit, maxcol, maxint, maxrow, ncol, ndig, nrow, output )

      if ( ierror .ne. 0 ) return
c
c  Save current value of objective function.
c
      if ( iform .eq. 0 ) then

        iobtop = iatop(nrow,ncol)
        iobbot = iabot(nrow,ncol)
        call ratrel ( objold, iatop(nrow,ncol), iabot(nrow,ncol) )

      else if ( iform .eq. 1 ) then

        objold = a(nrow,ncol)

      else if ( iform .eq. 2 ) then

        iobtop = iatop(nrow,ncol)
        iobbot = iabot(nrow,ncol)
        call decrel ( objold, iatop(nrow,ncol), iabot(nrow,ncol) )

      end if
c
c  Print out objective row.
c
      if ( iauto .eq. 0 ) then
        title = 'Objective row'
        ilo = nrow
        ihi = nrow
        jlo = 1
        jhi = ncol
 
        if ( iform .eq. 0 ) then

          call ratprn ( iatop, iabot, ibase, iounit, ihi, ilo, jhi,
     &      jlo, lpmoda, maxcol, maxrow, ncol, nrow, output, title )

        else if ( iform .eq. 1 ) then

          call relprn ( a, ibase, iounit, ihi, ilo, jhi, jlo, lpmoda,
     &      maxcol, maxrow, ncol, nrow, output, title )

        else if ( iform .eq. 2 ) then

          call decprn ( iatop, iabot, ibase, iounit, ihi, ilo, jhi,
     &      jlo, lpmoda, maxcol, maxrow, ncol, nrow, output, title )

        end if
 
      end if
c
c  Check for optimality.
c
      do j = 1, nvar+nslak+nart
 
        if ( iform .eq. 0 ) then
          call ratrel ( temp, iatop(nrow,j), iabot(nrow,j) )
        else if ( iform .eq. 1 ) then
          temp = a(nrow,j)
        else if ( iform .eq. 2 ) then
          call decrel ( temp, iatop(nrow,j), iabot(nrow,j) )
        end if
 
        if ( temp .lt. 0.0 ) go to 30
 
      end do
 
      call lpopt ( a, iatop, iabot, ibase, ierror, iform, imat,
     &  iopti, iounit, isltop, islbot, lpmoda, maxcol, maxrow, nart,
     &  ncol, nrow, nslak, nvar, output, sol )
 
      return
c
c  Choose the entering variable.
c
30    continue
 
      call lppiv2 ( a, iabot, iatop, iauto, ierror, iform, iounit,
     &  jpiv, line, maxcol, maxrow, nart, nline, nrow, nslak, nvar,
     &  output, prompt )
      if ( ierror .ne. 0 ) return
c
c  Choose the departing variable.
c
      call lppiv3 ( a, iabot, iatop, iauto, ibase, ierror, iform,
     &  iounit, ipiv, jpiv, line, maxcol, maxint, maxrow, ncol, nline,
     &  nrow, output, prompt )

      if ( ierror .ne. 0 ) return
c
c  Pivot on entry (IPIV,JPIV).
c
      irow=ipiv
 
      if ( iform .eq. 0 ) then
        istop = iatop(ipiv,jpiv)
        isbot = iabot(ipiv,jpiv)
      else if ( iform .eq. 1 ) then
        sval = a(ipiv,jpiv)
      else if ( iform .eq. 2 ) then
        istop = iatop(ipiv,jpiv)
        isbot = iabot(ipiv,jpiv)
      end if
 
      call scadiv ( a, iatop, iabot, ierror, iform, iounit, irow,
     &  maxcol, maxint, maxrow, ncol, ndig, nrow, output, sval,
     &  istop, isbot )
 
      irow2 = ipiv
 
      do i = 1, nrow
 
        irow1 = i
 
        if ( irow1 .ne. ipiv ) then
 
          if ( iform .eq. 0 ) then
            istop = - iatop(irow1,jpiv)
            isbot = iabot(irow1,jpiv)
          else if ( iform .eq. 1 ) then
            sval = -a(irow1,jpiv)
          else if ( iform .eq. 2 ) then
            istop = -iatop(irow1,jpiv)
            isbot = iabot(irow1,jpiv)
          end if
 
          call rowadd ( a, iatop, iabot, ierror, iform, iounit,
     &      irow1, irow2, maxcol, maxint, maxrow, ncol, ndig, 
     &      output, sval, istop, isbot )

        end if
 
      end do
c
c  Print out change in objective.
c
      output = ' '
      call chrwrt ( iounit, output )
 
      output = 'No change in objective.'
 
      if ( iform .eq. 0 ) then
 
        call ratrel ( objnew, iatop(nrow,ncol), iabot(nrow,ncol) )
 
        if ( objold .ne. objnew ) then
 
          chrtmp = chlrat ( iobtop, iobbot )
 
          if ( iobbot .ne. 1 ) then
            chrtmp2 = chrrel(objold)
            output = 'The objective changed from ' // chrtmp // ' = ' //
     &        chrtmp2
          else
            output = 'The objective changed from ' // chrtmp
          end if
 
          call chrdb2 ( output )
          call chrwrt ( iounit, output )
 
          chrtmp = chlrat ( iatop(nrow,ncol), iabot(nrow,ncol) )
 
          if ( iabot(nrow,ncol) .ne. 1 ) then
            chrtmp2 = chrrel ( objnew )
            output = 'to ' // chrtmp // ' = ' // chrtmp2
          else
            output = 'to ' // chrtmp
          end if
 
          call chrdb2 ( output )
          call chrwrt ( iounit, output )
 
        end if
 
      else if ( iform .eq. 1 ) then
 
        objnew = a(nrow,ncol)
 
        if ( objold .ne. objnew ) then
          chrtmp = chrrel(objold)
          chrtmp2 = chrrel(objnew)
          output = 'The objective changed from '//chrtmp//' to '
     &      //chrtmp2
          call chrdb2 ( output )
          call chrwrt ( iounit, output )
        end if
 
      else if ( iform .eq. 2 ) then
 
        call decrel ( objnew, iatop(nrow,ncol), iabot(nrow,ncol) )
 
        if ( objold .ne. objnew ) then
 
          chrtmp = chldec ( iobtop, iobbot )
          chrtmp2 = chldec ( iatop(nrow,ncol), iabot(nrow,ncol) )
 
          output = 'The objective changed from ' // chrtmp // ' to '
     &      // chrtmp2
          call chrdb2 ( output )
          call chrwrt ( iounit, output )
 
        end if
 
      end if
 
      if ( iauto .eq. 1 ) go to 10
 
      return
      end
      subroutine lppiv1 ( a, iabot, iatop, iauto, ibase, ierror, iform,
     &  iounit, maxcol, maxint, maxrow, ncol, ndig, nrow, output )
c
c***********************************************************************
c
cc LPPIV1 zeroes out objective row entries for basic variables.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, real A(MAXROW,MAXCOL).  A is the current matrix.
c
c    Input/output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
c    IATOP and IABOT represent the current rational or decimal matrix.
c
c    Input, integer IAUTO.
c    0, automatic processing is not being carried out.
c    1, automatic processing is being carried out.
c
c    Input, integer IBASE(MAXROW), keeps track of basic variables.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IFORM, specifies the arithmetic being used.
c    0=rational, 1=real, 2=decimal.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXINT, the maximum integer representable, which
c    should probably be 2147483647.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input, integer NCOL, the number of columns in the matrix.
c
c    Input, integer NDIG, the number of decimal digits used.
c
c    Input, integer NROW, the number of rows in the matrix.
c
c    Workspace, character*100 OUTPUT.
c
      integer maxcol
      integer maxrow

      real a(maxrow,maxcol)
      real amax
      character*6 chrint
      integer i
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer iauto
      integer ibase(maxrow)
      integer ierror
      integer iform
      integer imax
      integer iounit(4)
      integer irow
      integer irow1
      integer irow2
      integer isbot
      integer istop
      integer jcol
      integer maxint
      integer ncol
      integer ndig
      integer nrow
      character*100 output
      real sval
      real temp
c
      do i = 1, nrow-1
 
        jcol = ibase(i)
 
        if ( jcol .lt. 1 .or. jcol .gt. ncol ) then
          output = 'Error in the IBASE vector!'
          call chrwrt ( iounit, output )
          output = 'Entry ' // chrint(i) // ' of IBASE = ' 
     &      // chrint(jcol)
          call chrwrt ( iounit, output )
          ierror = 1
          return
        end if
c
c  Check the objective entry in column JCOL.
c
        if ( iform .eq. 0 ) then
          call ratrel ( temp, iatop(nrow,jcol), iabot(nrow,jcol) )
        else if ( iform .eq. 1 ) then
          temp = a(nrow,jcol)
        else if ( iform .eq. 2 ) then
          call decrel ( temp, iatop(nrow,jcol), iabot(nrow,jcol) )
        end if
 
        if ( temp .ne. 0 ) then
 
          irow1 = nrow
          output = ' '
          call chrwrt ( iounit, output )
          output = 'The objective entry in column '//chrint(jcol)//
     &      ' is not zero,'
          call chrdb2 ( output )
          call chrwrt ( iounit, output )
          output = 'but this corresponds to a basic variable.'
          call chrwrt ( iounit, output )
          output = ' '
          call chrwrt ( iounit, output )
 
          if ( iauto .eq. 0 ) then
            output = 'Use the "A" command to zero out this entry.'
            call chrwrt ( iounit, output )
            output = 'THEN you may use the "P" command!'
            call chrwrt ( iounit, output )
            ierror = 1
            return
          end if
c
c  Search for the IMAX, the row of the maximum entry in column JCOL.
c
          amax = 0.0
          imax = 0
 
          do irow= 1, nrow-1
 
            if ( iform .eq. 0 ) then
              call ratrel ( temp, iatop(irow,jcol), iabot(irow,jcol) )
            else if ( iform .eq. 1 ) then
              temp = a(irow,jcol)
            else if ( iform .eq. 2 ) then
              call decrel ( temp, iatop(irow,jcol), iabot(irow,jcol) )
            end if
 
            temp = abs(temp)
 
            if ( temp .gt. amax ) then
              amax = temp
              imax = irow
            end if
 
          end do
 
          if ( amax .eq. 0.0 ) then
            output = 'The artificial variable cannot be eliminated!'
            call chrwrt ( iounit, output )
            go to 20
          end if
c
c  Divide row IMAX by entry (IMAX,JCOL) to normalize it.
c
          irow = imax
 
          if ( iform .eq. 0 ) then
            istop = iatop(irow,jcol)
            isbot = iabot(irow,jcol)
          else if ( iform .eq. 1 ) then
            sval = a(irow,jcol)
          else if ( iform .eq. 2 ) then
            istop = iatop(irow,jcol)
            isbot = iabot(irow,jcol)
          end if
 
          call scadiv ( a, iatop, iabot, ierror, iform, iounit, irow,
     &      maxcol, maxint, maxrow, ncol, ndig, nrow, output, sval,
     &      istop, isbot )
c
c  Add a multiple of row IMAX to row NROW, to eliminate entry (NROW,JCOL).
c
          irow2 = imax
          irow1 = nrow
 
          if ( iform .eq. 0 ) then
            istop = -iatop(irow1,jcol)
            isbot = iabot(irow1,jcol)
          else if ( iform .eq. 1 ) then
            sval = -a(irow1,jcol)
          else if ( iform .eq. 2 ) then
            istop = -iatop(irow1,jcol)
            isbot = iabot(irow1,jcol)
          end if
 
          call rowadd ( a, iatop, iabot, ierror, iform, iounit,
     &      irow1, irow2, maxcol, maxint, maxrow, ncol, ndig, 
     &      output, sval, istop, isbot )

          if ( iform .eq. 0 ) then
            iatop(nrow,jcol) = 0
            iabot(nrow,jcol) = 1
          else if ( iform .eq. 1 ) then
            a(nrow,jcol) = 0.0
          else if ( iform .eq. 2 ) then
            iatop(nrow,jcol) = 0
            iabot(nrow,jcol) = 0
          end if
 
        end if
 
20      continue
 
      end do
 
      return
      end
      subroutine lppiv2 ( a, iabot, iatop, iauto, ierror, iform,
     &  iounit, jpiv, line, maxcol, maxrow, nart, nline, nrow, nslak,
     &  nvar, output, prompt )
c
c***********************************************************************
c
cc LPPIV2 chooses the entering variable for pivoting.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real A(MAXROW,MAXCOL).  A is the current matrix.
c
c    Input, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
c    IATOP and IABOT represent the current rational or decimal matrix.
c
c    Input, integer IAUTO.
c    0, automatic processing is not being carried out.
c    1, automatic processing is being carried out.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IFORM, specifies the arithmetic being used.
c    0=rational, 1=real, 2=decimal.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Output, integer JPIV, the entering variable.
c
c    Workspace, character*80 LINE.
c    Used to hold the user's input.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input, integer NART, the number of artificial variables.
c
c    Input/output, integer NLINE.
c    Keeps track of the number of useful characters in LINE.
c
c    Input, integer NROW, the number of rows in the matrix.
c
c    Input, integer NSLAK, the number of slack variables.
c
c    Input, integer NVAR, the number of basic variables.
c
c    Workspace, character*100 OUTPUT.
c
c    Workspace, character*80 PROMPT.
c
      integer maxcol
      integer maxrow

      real a(maxrow,maxcol)
      character*6 chrint
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer iauto
      integer ierror
      integer iform
      integer iounit(4)
      integer j
      integer jmin
      integer jpiv
      character*80 line
      integer nart
      integer nline
      integer nrow
      integer nslak
      integer nvar
      character*100 output
      character*80 prompt
      real temp
      real tmin
c
10    continue
c
c  Set TMIN to the (NROW,1) entry, and JMIN to 1.
c
      if ( iform .eq. 0 ) then
        call ratrel ( tmin, iatop(nrow,1), iabot(nrow,1) )
      else if ( iform .eq. 1 ) then
        tmin = a(nrow,1)
      else if ( iform .eq. 2 ) then
        call decrel ( tmin, iatop(nrow,1), iabot(nrow,1) )
      end if
 
      jmin = 1
c
c  Set TMIN to the (NROW,J) entry, if it is smaller.
c
      do j = 2, nvar+nslak+nart
 
        if ( iform .eq. 0 ) then
          call ratrel ( temp, iatop(nrow,j), iabot(nrow,j) )
        else if ( iform .eq. 1 ) then
          temp = a(nrow,j)
        else if ( iform .eq. 2 ) then
          call decrel ( temp, iatop(nrow,j), iabot(nrow,j) )
        end if
 
        if ( temp .le. tmin ) then
          tmin = temp
          jmin = j
        end if
 
      end do
c
c  Now get pivot index JPIV.
c
      if ( iauto .eq. 1 ) then
 
        jpiv = jmin
 
      else
 
        output = ' '
        call chrwrt ( iounit, output )
        output = 'Variable with most negative objective coefficient?'
        call chrwrt ( iounit, output )
 
        prompt = 'column (=variable number)'
        call intrea ( jpiv, line, nline, prompt, iounit, ierror )
        if ( ierror .ne. 0 ) return
 
        if ( jpiv .lt. 1 .or. jpiv .gt. nvar+nslak+nart ) then
          output = 'Your input was out of bounds.'
          call chrwrt ( iounit, output )
          go to 10
        end if
 
        if ( iform .eq. 0 ) then
          call ratrel ( temp, iatop(nrow,jpiv), iabot(nrow,jpiv) )
        else if ( iform .eq. 1 ) then
          temp = a(nrow,jpiv)
        else if ( iform .eq. 2 ) then
          call decrel ( temp, iatop(nrow,jpiv), iabot(nrow,jpiv) )
        end if
 
        if ( temp .gt. tmin+0.0001 ) then
          output = 'Not acceptable.'
          call chrwrt ( iounit, output )
          go to 10
        end if
 
      end if
 
      output = ' '
      call chrwrt ( iounit, output )
      output = 'The entering variable is ' // chrint(jpiv)
      call chrdb2 ( output )
      call chrwrt ( iounit, output )
 
      return
      end
      subroutine lppiv3 ( a, iabot, iatop, iauto, ibase, ierror, iform,
     &  iounit, ipiv, jpiv, line, maxcol, maxint, maxrow, ncol, nline,
     &  nrow, output, prompt )
c
c***********************************************************************
c
cc LPPIV3 chooses the departing variable for pivoting.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real A(MAXROW,MAXCOL).  A is the current matrix.
c
c    Input, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
c    IATOP and IABOT represent the current rational or decimal matrix.
c
c    Input, integer IAUTO.
c    0, automatic processing is not being carried out.
c    1, automatic processing is being carried out.
c
c    Input, integer IBASE(MAXROW), keeps track of basic variables.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IFORM, specifies the arithmetic being used.
c    0=rational, 1=real, 2=decimal.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Output, integer IPIV, the row of the departing variable.
c
c    Input, integer JPIV, the entering variable.
c
c    Workspace, character*80 LINE.
c    Used to hold the user's input.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXINT, the maximum integer representable, which
c    should probably be 2147483647.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input, integer NCOL, the number of columns in the matrix.
c
c    Input/output, integer NLINE.
c    Keeps track of the number of useful characters in LINE.
c
c    Input, integer NROW, the number of rows in the matrix.
c
c    Workspace, character*100 OUTPUT.
c
c    Workspace, character*80 PROMPT.
c
      integer maxcol
      integer maxrow

      real a(maxrow,maxcol)
      real bot
      character*22 chlrat
      character*6 chrint
      character*14 chrrel
      character*22 chrtmp
      integer i
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer iauto
      integer ibase(maxrow)
      integer ibot
      integer ierror
      integer iform
      integer imin
      integer iounit(4)
      integer ipiv
      integer itop
      integer jpiv
      character*80 line
      integer maxint
      integer ncol
      integer nline
      integer nrow
      character*100 output
      character*80 prompt
      real ratio
      real ratj
      real ratmin
      real temp1
      real temp2
      real top
c
      if ( iauto .eq. 0 ) then
        output = ' '
        call chrwrt ( iounit, output )
        output = 'Variable with smallest nonnegative feasibility ratio?'
        call chrwrt ( iounit, output )
      end if
 
10    continue
 
      imin = 0

      if ( iauto .eq. 0 ) then
        output = ' '
        call chrwrt ( iounit, output )
        output = 'Nonnegative feasibility ratios:'
        call chrwrt ( iounit, output )
        output = ' '
        call chrwrt ( iounit, output )
      end if
 
      ratmin = - 1.0
 
      do i = 1, nrow-1
 
        if ( iform .eq. 0 ) then
 
          if ( iabot(i,jpiv) .lt. 0 ) then
            iatop(i,jpiv) = - iatop(i,jpiv)
            iabot(i,jpiv) = - iabot(i,jpiv)
          end if
 
          if ( iatop(i,jpiv) .le. 0 ) go to 20
 
          call ratrel ( top, iatop(i,ncol), iabot(i,ncol) )
          call ratrel ( bot, iatop(i,jpiv), iabot(i,jpiv) )
 
        else if ( iform .eq. 1 ) then
 
          if ( a(i,jpiv) .le. 0.0 ) go to 20
 
          top = a(i,ncol)
          bot = a(i,jpiv)
 
        else if ( iform .eq. 2 ) then
 
          if ( iatop(i,jpiv) .le. 0 ) go to 20
 
          call decrel ( top, iatop(i,ncol), iabot(i,ncol) )
          call decrel ( bot, iatop(i,jpiv), iabot(i,jpiv) )
 
        end if
 
        if ( bot .eq. 0.0 ) go to 20
 
        ratio = top / bot
 
        if ( iauto .eq. 0 ) then
 
          if ( iform .eq. 0 ) then
 
            call ratdiv ( ibot, iabot(i,ncol), iabot(i,jpiv), ierror,
     &        itop, iatop(i,ncol), iatop(i,jpiv), maxint )
 
            if ( ibot .ne. 1 ) then
              chrtmp = chlrat ( itop, ibot )
              write ( output, 120 ) i, ibase(i), ratio, chrtmp
            else
              write ( output, 111 ) i, ibase(i), ratio
            end if
 
          else if ( iform .eq. 1 ) then
 
            write ( output, 111 ) i, ibase(i), ratio
 
          else if ( iform .eq. 2 ) then
 
            write ( output, 111 ) i, ibase(i), ratio
 
          end if
 
          call chrdb2 ( output )
          call chrwrt ( iounit, output )
 
        end if
 
        if ( imin .eq. 0 .and. 0.0 .le. ratio ) then
          imin = i
          ratmin = ratio
        end if
 
        if ( 0.0 .le. ratio .and. ratio .lt. ratmin ) then
          ratmin = ratio
          imin = i
        end if
 
20      continue
 
      end do
 
      if ( imin .eq. 0 ) then
        output = ' '
        call chrwrt ( iounit, output )
        output = 'Cannot find a departing variable.'
        call chrwrt ( iounit, output )
        output = 'Presumably, the feasible set is unbounded.'
        call chrwrt ( iounit, output )
        ierror = 1
        return
      end if
 
30    continue
 
      if ( iauto .eq. 1 ) then

        ipiv = imin

      else

        output = ' '
        call chrwrt ( iounit, output )
        prompt = 'the row of the departing variable.'
        call intrea ( ipiv, line, nline, prompt, iounit, ierror )
        if ( ierror .ne. 0 ) return
 
        if ( ipiv .le. 0 .or. ipiv .gt. nrow-1 ) then
          output = 'Illegal row.'
          call chrwrt ( iounit, output )
          go to 30
        end if
 
        if ( iform .eq. 0 ) then
 
          if ( iatop(ipiv,jpiv) .eq. 0 ) then
            output = 'Illegal zero divisor.'
            call chrwrt ( iounit, output )
            go to 10
          end if
 
        else if ( iform .eq. 1 ) then
 
          if ( a(ipiv,jpiv) .eq. 0 ) then
            output = 'Illegal zero divisor.'
            call chrwrt ( iounit, output )
            go to 10
          end if
 
        else if ( iform .eq. 2 ) then
 
          if ( iatop(ipiv,jpiv) .eq. 0 ) then
            output = 'Illegal zero divisor.'
            call chrwrt ( iounit, output )
            go to 10
          end if
 
        end if
 
        if ( iform .eq. 0 ) then
 
          call ratrel ( temp1, iatop(ipiv,ncol), iabot(ipiv,ncol) )
          call ratrel ( temp2, iatop(ipiv,jpiv), iabot(ipiv,jpiv) )
 
        else if ( iform .eq. 1 ) then
 
          temp1 = a(ipiv,ncol)
          temp2 = a(ipiv,jpiv)
 
        else if ( iform .eq. 2 ) then
 
          call decrel ( temp1, iatop(ipiv,ncol), iabot(ipiv,ncol) )
          call decrel ( temp2, iatop(ipiv,jpiv), iabot(ipiv,jpiv) )
 
        end if
 
        ratj = temp1 / temp2
 
        if ( ratj .lt. 0.0 ) then
          output = 'The pivot ratio is not acceptable because ' //
     &      'it is negative.'
          call chrwrt ( iounit, output )
          go to 10
        else if ( ratmin .lt. ratj ) then
          output = 'The pivot ratio is not acceptable because ' //
     &      'it is not the smallest nonnegative ratio.'
          call chrwrt ( iounit, output )
          go to 10
        end if
 
      end if
 
      output = ' '
      call chrwrt ( iounit, output )
 
      output = 'The departing variable is ' // chrint(ibase(ipiv))
     &  // ' with feasibility ratio ' // chrrel(ratmin)
      call chrdb2 ( output )
      call chrwrt ( iounit, output )
 
      output = ' '
      call chrwrt ( iounit, output )
 
      ibase(ipiv) = jpiv
 
  111 format('Row ',i2,', variable ',i2,', ratio = ',g14.6)
  120 format('Row ',i2,', variable ',i2,', ratio = ',g14.6,' = ',
     &  a22)
 
      return
      end
      subroutine lprem ( a, iabot, iatop, ibase, ierror, iform, imat,
     &  iounit, lpmoda, maxcol, maxrow, nart, ncol, nrow, nslak, nvar,
     &  output )
c
c***********************************************************************
c
cc LPREM removes the artificial variables.
c
c
c  LPREM should be called once the artificial objective function has 
c  reached zero.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, real A(MAXROW,MAXCOL).  A is the current matrix.
c
c    Input/output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
c    IATOP and IABOT represent the current rational or decimal matrix.
c
c    Input, integer IBASE(MAXROW), keeps track of basic variables.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IFORM, specifies the arithmetic being used.
c    0=rational, 1=real, 2=decimal.
c
c    Input, integer IMAT.
c    0, no matrix has been defined by the user.
c    1, a matrix has been defined by the user.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Input, integer LPMODA.
c    0, the program is in linear algebra mode.
c    1, the program is in linear programming mode.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input, integer NART, the number of artificial variables.
c
c    Input/output, integer NCOL, the number of columns in the matrix.
c    On output, this number may have changed because of the elimination
c    of artificial variables.
c
c    Input, integer NROW, the number of rows in the matrix.
c
c    Input, integer NSLAK, the number of slack variables.
c
c    Input, integer NVAR, the number of basic variables.
c
c    Workspace, character*100 OUTPUT.
c
      integer maxcol
      integer maxrow
c
      real a(maxrow,maxcol)
      character*6 chrint
      integer i
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer ibase(maxrow)
      integer ierror
      integer iform
      integer imat
      integer inext
      integer iounit(4)
      integer j
      integer jhi
      integer jvar
      integer lpmoda
      integer mart
      integer nart
      integer ncol
      integer nrow
      integer nslak
      integer nvar
      character*100 output

      if ( lpmoda .ne. 1 ) then
        ierror = 1
        output = 'This command should only be given during'
        call chrwrt ( iounit, output )
        output = 'linear programming!'
        call chrwrt ( iounit, output )
        return
      end if
 
      if ( imat .eq. 0 ) then
        ierror = 1
        output = 'You must set up a tableau first!'
        call chrwrt ( iounit, output )
        return
      end if
 
      if ( nart .eq. 0 ) then
        output = 'There aren''t any artificial variables to delete!'
        call chrwrt ( iounit, output )
        return
      end if
 
      if ( 
     &  (iform .eq. 0 .and. iatop(nrow,ncol) .ne. 0) .or. 
     &  (iform .eq. 1 .and. a(nrow,ncol) .ne. 0.0) .or.
     &  (iform .eq. 2 .and. iatop(nrow,ncol) .ne. 0) ) then

        output = 'The phase 1 objective function is nonzero.'
        call chrwrt ( iounit, output )
        output = 'Hence, this problem may have no solution.'
        call chrwrt ( iounit, output )

      end if
 
      jhi = nvar+nslak+nart+2
      mart = nart
      nart = 0
      inext = nvar+nslak
 
      do jvar = nvar+nslak+1, jhi
 
        if ( jvar .gt. jhi-2 ) go to 30
 
        do i = 1, nrow-1
          if ( ibase(i) .eq. jvar ) then
            nart = nart+1
            go to 30
          end if
        end do
 
        do i = 1, nrow-1
          if ( ibase(i) .gt. jvar) then
            ibase(i) = ibase(i)-1
          end if
        end do
 
        go to 50
 
30      continue

        inext = inext+1
 
        if ( iform .eq. 0 ) then
 
          do i = 1, nrow
            iatop(i,inext) = iatop(i,jvar)
            iabot(i,inext) = iabot(i,jvar)
          end do
 
        else if ( iform .eq. 1 ) then
 
          do i = 1, nrow
            a(i,inext) = a(i,jvar)
          end do
 
        else if ( iform .eq. 2 ) then
 
          do i = 1, nrow
            iatop(i,inext) = iatop(i,jvar)
            iabot(i,inext) = iabot(i,jvar)
          end do
 
        end if
 
50      continue
 
      end do
c
c  If possible, restore the original objective function.
c
      ncol = nvar+nslak+2
 
      output = ' '
      call chrwrt ( iounit, output )
 
      if ( nart .ne. 0 ) then
        output = chrint(nart)//' artificial variables were not deleted.'
        call chrdb2 ( output )
        call chrwrt ( iounit, output )
        output = 'You must revise the objective row by hand!'
        call chrwrt ( iounit, output )
      else
        output = 'All the artificial variables were deleted.'
        call chrwrt ( iounit, output )
        output = 'The original objective function is restored.'
        call chrwrt ( iounit, output )
 
        do j = 1, nvar+nslak
 
          if ( iform .eq. 0 ) then
            iatop(nrow,j) = iatop(nrow+1,j)
            iabot(nrow,j) = iabot(nrow+1,j)
          else if ( iform .eq. 1 ) then
            a(nrow,j) = a(nrow+1,j)
          else if ( iform .eq. 2 ) then
            iatop(nrow,j) = iatop(nrow+1,j)
            iabot(nrow,j) = iabot(nrow+1,j)
          end if
 
        end do
 
        do j = nvar+nslak+1, nvar+nslak+2
 
          if ( iform .eq. 0 ) then
            iatop(nrow,j) = iatop(nrow+1,j+mart)
            iabot(nrow,j) = iabot(nrow+1,j+mart)
          else if ( iform .eq. 1 ) then
            a(nrow,j) = a(nrow+1,j+mart)
          else if ( iform .eq. 2 ) then
            iatop(nrow,j) = iatop(nrow+1,j+mart)
            iabot(nrow,j) = iabot(nrow+1,j+mart)
          end if
 
        end do
 
      end if
 
      output = ' '
      call chrwrt ( iounit, output )
      output = 'You must now use the "A" command to zero out'
      call chrwrt ( iounit, output )
      output = 'objective row entries for all basic variables.'
      call chrwrt ( iounit, output )
 
      return
      end
      subroutine lpsama ( a, chineq, iatop, iabot, ibase, iform, imat,
     &  iounit, maxcol, maxrow, nart, ncol, nrow, nslak, nvar, 
     &  output )
c
c*********************************************************************
c
cc LPSAMA sets up an advanced linear programming problem.
c
c
c  The problem includes artificial variables.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, real A(MAXROW,MAXCOL).  A is the current matrix.
c
c    Output, character*1 CHINEQ(MAXROW), the '<', '=', or '>'
c    sign for each linear programming constraint.
c
c    Output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
c    IATOP and IABOT represent the current rational or decimal matrix.
c
c    Output, integer IBASE(MAXROW), keeps track of basic variables.
c
c    Input, integer IFORM, specifies the arithmetic being used.
c    0=rational, 1=real, 2=decimal.
c
c    Output, integer IMAT.
c    0, no matrix has been defined by the user.
c    1, a matrix has been defined by the user.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Output, integer NART, the number of artificial variables.
c
c    Output, integer NCOL, the number of columns in the matrix.
c
c    Output, integer NROW, the number of rows in the matrix.
c
c    Output, integer NSLAK, the number of slack variables.
c
c    Output, integer NVAR, the number of basic variables.
c
c    Workspace, character*100 OUTPUT.
c
      integer maxcol
      integer maxrow
c
      real a(maxrow,maxcol)
      character*1 chineq(maxrow)
      integer i
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer ibase(maxrow)
      integer iform
      integer imat
      integer iounit(4)
      integer j
      integer nart
      integer ncol
      integer nrow
      integer nslak
      integer nvar
      character*100 output
c
c  Zero out the matrix.
c
      call inimat ( a, iabot, iatop, iform, maxcol, maxrow )

      nvar = 2
      nslak = 4
      nart = 2
      nrow = nslak+1
      ncol = nvar+nslak+nart+2
 
      iatop(1,1) = 1
      iatop(1,2) = 2
      iatop(1,3) = -1
      iatop(1,4) = 0
      iatop(1,5) = 0
      iatop(1,6) = 0
      iatop(1,7) = 1
      iatop(1,8) = 0
      iatop(1,9) = 0
      iatop(1,10) = 6
 
      iatop(2,1) = 2
      iatop(2,2) = 1
      iatop(2,3) = 0
      iatop(2,4) = -1
      iatop(2,5) = 0
      iatop(2,6) = 0
      iatop(2,7) = 0
      iatop(2,8) = 1
      iatop(2,9) = 0
      iatop(2,10) = 4
 
      iatop(3,1) = 1
      iatop(3,2) = 1
      iatop(3,3) = 0
      iatop(3,4) = 0
      iatop(3,5) = 1
      iatop(3,6) = 0
      iatop(3,7) = 0
      iatop(3,8) = 0
      iatop(3,9) = 0
      iatop(3,10) = 5
 
      iatop(4,1) = 2
      iatop(4,2) = 1
      iatop(4,3) = 0
      iatop(4,4) = 0
      iatop(4,5) = 0
      iatop(4,6) = 1
      iatop(4,7) = 0
      iatop(4,8) = 0
      iatop(4,9) = 0
      iatop(4,10) = 8
 
      iatop(5,1) = 0
      iatop(5,2) = 0
      iatop(5,3) = 0
      iatop(5,4) = 0
      iatop(5,5) = 0
      iatop(5,6) = 0
      iatop(5,7) = 1
      iatop(5,8) = 1
      iatop(5,9) = 1
      iatop(5,10) = 0
 
      iatop(6,1) = -40
      iatop(6,2) = -30
      iatop(6,3) = 0
      iatop(6,4) = 0
      iatop(6,5) = 0
      iatop(6,6) = 0
      iatop(6,7) = 0
      iatop(6,8) = 0
      iatop(6,9) = 1
      iatop(6,10) = 0
 
      do i = 1, nrow+1
        do j = 1, ncol
          if ( iform .eq. 0 ) then
            iabot(i,j) = 1
          else if ( iform .eq. 2 ) then
            iabot(i,j) = 0
          end if
        end do
      end do
 
      do i = 1, nrow+1
        do j = 1, ncol
          a(i,j) = real ( iatop(i,j) )
        end do
      end do
 
      ibase(1) = 7
      ibase(2) = 8
      ibase(3) = 5
      ibase(4) = 6
 
      chineq(1) = '>'
      chineq(2) = '>'
      chineq(3) = '<'
      chineq(4) = '<'
      chineq(5) = ' '
 
      imat = 1
 
      output = ' '
      call chrwrt ( iounit, output )
      output = 'Advanced linear programming problem:'
      call chrwrt ( iounit, output )
      output = ' '
      call chrwrt ( iounit, output )
      output = 'Maximize'
      call chrwrt ( iounit, output )
      output = '  Z=40 X + 30 Y'
      call chrwrt ( iounit, output )
      output = 'subject to'
      call chrwrt ( iounit, output )
      output = '  X + 2 Y > 6'
      call chrwrt ( iounit, output )
      output = '2 X +   Y > 4'
      call chrwrt ( iounit, output )
      output = '  X +   Y < 5'
      call chrwrt ( iounit, output )
      output = '2 X +   Y < 8'
      call chrwrt ( iounit, output )
      output = ' '
      call chrwrt ( iounit, output )
 
      return
      end
      subroutine lpsams ( a, chineq, iatop, iabot, ibase, iform,
     &  imat, iounit, maxcol, maxrow, nart, ncol, nrow, nslak,
     &  nvar, output )
c
c*********************************************************************
c
cc LPSAMS sets up a simple linear programming problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, real A(MAXROW,MAXCOL).  A is the current matrix.
c
c    Output, character*1 CHINEQ(MAXROW), the '<', '=', or '>'
c    sign for each linear programming constraint.
c
c    Output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
c    IATOP and IABOT represent the current rational or decimal matrix.
c
c    Output, integer IBASE(MAXROW), keeps track of basic variables.
c
c    Input, integer IFORM, specifies the arithmetic being used.
c    0=rational, 1=real, 2=decimal.
c
c    Output, integer IMAT.
c    0, no matrix has been defined by the user.
c    1, a matrix has been defined by the user.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Output, integer NART, the number of artificial variables.
c
c    Output, integer NCOL, the number of columns in the matrix.
c
c    Output, integer NROW, the number of rows in the matrix.
c
c    Output, integer NSLAK, the number of slack variables.
c
c    Output, integer NVAR, the number of basic variables.
c
c    Workspace, character*100 OUTPUT.
c
      integer maxcol
      integer maxrow
c
      real a(maxrow,maxcol)
      character*1 chineq(maxrow)
      integer i
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer ibase(maxrow)
      integer iform
      integer imat
      integer iounit(4)
      integer j
      integer nart
      integer ncol
      integer nrow
      integer nslak
      integer nvar
      character*100 output
c
c  Zero out the matrix.
c
      call inimat ( a, iabot, iatop, iform, maxcol, maxrow )
c
      nvar = 2
      nslak = 2
      nart = 0
      nrow = nslak+1
      ncol = nvar+nslak+nart+2
 
      iatop(1,1) = 2
      iatop(1,2) = 2
      iatop(1,3) = 1
      iatop(1,4) = 0
      iatop(1,5) = 0
      iatop(1,6) = 8
 
      iatop(2,1) = 5
      iatop(2,2) = 3
      iatop(2,3) = 0
      iatop(2,4) = 1
      iatop(2,5) = 0
      iatop(2,6) = 15
 
      iatop(3,1) = -120
      iatop(3,2) = -100
      iatop(3,3) = 0
      iatop(3,4) = 0
      iatop(3,5) = 1
      iatop(3,6) = 70
 
      do i = 1, nrow
        do j = 1, ncol
          if ( iform .eq. 0 ) then
            iabot(i,j) = 1
          else if ( iform .eq. 2 ) then
            iabot(i,j) = 0
          end if
        end do
      end do
 
      do i = 1, nrow
        do j = 1, ncol
          a(i,j) = real(iatop(i,j))
        end do
      end do
 
      ibase(1) = 3
      ibase(2) = 4
 
      chineq(1) = '<'
      chineq(2) = '<'
      chineq(3) = ' '
 
      imat = 1
 
      output = ' '
      call chrwrt ( iounit, output )
      output = 'Simple linear programming problem:'
      call chrwrt ( iounit, output )
      output = ' '
      call chrwrt ( iounit, output )
      output = 'Maximize:'
      call chrwrt ( iounit, output )
      output = '  Z = 120 X + 100 Y + 70'
      call chrwrt ( iounit, output )
      output = 'subject to'
      call chrwrt ( iounit, output )
      output = '  2 X + 2 Y < 8'
      call chrwrt ( iounit, output )
      output = '  5 X + 3 Y < 15'
      call chrwrt ( iounit, output )
      output = ' '
      call chrwrt ( iounit, output )
 
      return
      end
      subroutine lpset ( ierror, imat, iounit, line, lpmoda, nart,
     &  ncol, ncon, nline, nrow, nslak, nvar, output, prompt )

c***********************************************************************
c
cc LPSET switches the linear programming mode.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IERROR.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IMAT.
c    0, no matrix has been defined by the user.
c    1, a matrix has been defined by the user.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Workspace, character*80 LINE.
c    Used to hold the user's input.
c
c    Input/output, integer LPMODA.
c    0, the program is in linear algebra mode.
c    1, the program is in linear programming mode.
c
c    Input, integer NART, the number of artificial variables.
c
c    Input, integer NCOL, the number of columns in the matrix.
c
c    Output, integer NCON, the number of constraints.
c
c    Input/output, integer NLINE.
c    Keeps track of the number of useful characters in LINE.
c
c    Input, integer NROW, the number of rows in the matrix.
c
c    Output, integer NSLAK, the number of slack variables.
c
c    Output, integer NVAR, the number of basic variables.
c
c    Workspace, character*100 OUTPUT.
c
c    Workspace, character*80 PROMPT.
c
      integer ierror
      integer imat
      integer iounit(4)
      character*1 isay
      integer iterm
      logical leqi
      character*80 line
      integer lpmoda
      integer nart
      integer ncol
      integer ncon
      integer nline
      integer nrow
      integer nslak
      integer nvar
      character*100 output
      character*80 prompt

      lpmoda = 1 - lpmoda
 
      if ( lpmoda .eq. 0 ) then
        output = 'Switching to linear algebra mode.'
        call chrwrt ( iounit, output )
        return
      end if
 
      output = 'Switching to linear programming mode.'
      call chrwrt ( iounit, output )
 
      if ( imat .eq. 0 ) return
 
      prompt = '"Y" to use current matrix in linear programming.'
      nline = 0
      iterm = 0
      call chrrea ( isay, line, nline, prompt, iounit, ierror, iterm )
      if ( ierror .ne. 0 ) return
 
      if ( .not. leqi ( isay, 'Y' ) ) then
        imat = 0
        nrow = 0
        ncol = 0
        return
      end if
 
      output = ' '
      call chrwrt ( iounit, output )

10    continue

      nline = 0
      prompt = '# of slack variables, # of artificial variables.'
      call intrea ( nslak, line, nline, prompt, iounit, ierror )
      if ( ierror .ne. 0 ) return
 
      call intrea ( nart, line, nline, prompt, iounit, ierror )
      if ( ierror .ne. 0 ) return
 
      nvar = ncol - 2 - nart - nslak
      ncon = nrow - 1
 
      if ( nvar .le. 0 ) then
        output = 'Values too large or too small!'
        call chrwrt ( iounit, output )
        return
      end if
 
      output = ' '
      call chrwrt ( iounit, output )
      output = 'Now please set the row labels (=basic variables)'
      call chrwrt ( iounit, output )
      output = 'using the "C" command, with I2 = 0.'
      call chrwrt ( iounit, output )
 
      return
      end
      subroutine lpsol ( a, iatop, iabot, ibase, iform, isltop, islbot,
     &  maxcol, maxrow, ncol, nrow, sol )
c
c***********************************************************************
c
cc LPSOL determines the current linear programming solution.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real A(MAXROW,MAXCOL).  A is the current matrix.
c
c    Input, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
c    IATOP and IABOT represent the current rational or decimal matrix.
c
c    Input, integer IBASE(MAXROW), keeps track of basic variables.
c
c    Input, integer IFORM, specifies the arithmetic being used.
c    0=rational, 1=real, 2=decimal.
c
c    Output, integer ISLTOP(MAXROW), ISLBOT(MAXROW), the fractional
c    or decimal representation of the linear programming solution.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input, integer NCOL, the number of columns in the matrix.
c
c    Input, integer NROW, the number of rows in the matrix.
c
c    Output, real SOL(MAXROW), the real representation of the
c    linear programming solution.
c
      integer maxcol
      integer maxrow

      real a(maxrow,maxcol)
      integer i
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer ibase(maxrow)
      integer iform
      integer islbot(maxcol)
      integer isltop(maxcol)
      integer j
      integer jbase
      integer ncol
      integer nrow
      real sol(maxcol)

      do i = 1, ncol-2
 
        jbase = 0
 
        do j = 1, nrow-1
          if ( ibase(j) .eq. i ) then
            jbase = j
          end if
        end do
 
        if ( jbase .ne. 0 ) then
 
          if ( iform .eq. 0 ) then
            isltop(i) = iatop(jbase,ncol)
            islbot(i) = iabot(jbase,ncol)
          else if ( iform .eq. 1 ) then
            sol(i) = a(jbase,ncol)
          else if ( iform .eq. 2 ) then
            isltop(i) = iatop(jbase,ncol)
            islbot(i) = iabot(jbase,ncol)
          end if
 
        else
 
          if ( iform .eq. 0 ) then
            isltop(i) = 0
            islbot(i) = 1
          else if ( iform .eq. 1 ) then
            sol(i) = 0.0
          else if ( iform .eq. 2 ) then
            isltop(i) = 0
            islbot(i) = 0
          end if
 
        end if
 
      end do
 
      return
      end
      subroutine mulply ( a, dete, iatop, iabot, idetop, idebot,
     &  ierror, iform, iounit, irow, maxcol, maxint, maxrow, ncol,
     &  ndig, nrow, output, sval, istop, isbot )
c
c***********************************************************************
c
cc MULPLY multiplies a row of the A matrix by a scale factor.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, real A(MAXROW,MAXCOL).  A is the current matrix.
c
c    Input, real DETE, the determinant of the product of the
c    elementary row operations applied to the current matrix.
c
c    Input/output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
c    IATOP and IABOT represent the current rational or decimal matrix.
c
c    Input/output, integer IDETOP, IDEBOT, the rational or
c    decimal representation of the determinant of the product of
c    the elementary row operations applied to the current matrix.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IFORM, specifies the arithmetic being used.
c    0=rational, 1=real, 2=decimal.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Input, integer IROW, the row that is to be multiplied.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXINT, the maximum integer representable, which
c    should probably be 2147483647.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input, integer NCOL, the number of columns in the matrix.
c
c    Input, integer NDIG, the number of decimal digits used.
c
c    Input, integer NROW, the number of rows in the matrix.
c
c    Workspace, character*100 OUTPUT.
c
c    Input, real SVAL, the real row multiplier.
c
c    Input, integer ISTOP, ISBOT, the decimal or fractional row multiplier.
c
      integer maxcol
      integer maxrow
c
      real a(maxrow,maxcol)
      character*22 chldec
      character*22 chlrat
      character*6 chrint
      character*14 chrrel
      character*22 chrtmp
      real dete
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer ibot
      integer idebot
      integer idetop
      integer ierror
      integer iform
      integer iounit(4)
      integer irow
      integer isbot
      integer istop
      integer itop
      integer j
      integer maxint
      integer ncol
      integer ndig
      integer nrow
      character*100 output
      real sval
c
c  Make sure row number is OK.
c
      if ( irow .lt. 1 .or. irow .gt. nrow ) then
        output = 'Error!  The row number is out of range!'
        call chrwrt ( iounit, output )
        ierror = 1
        return
      end if
c
c  For rational arithmetic, make sure bottom of scale factor
c  is not 0.
c
      if ( iform .eq. 0 ) then
        if ( isbot .eq. 0 ) then
          output = 'Error!  Illegal 0 divisor in multiplier!'
          call chrwrt ( iounit, output )
          ierror = 1
          return
        end if
      end if
c
c  Check for multiplication by 0 or 1.
c
      if ( iform .eq. 0 ) then
        if ( istop .eq. 0 ) then
          output = 'Warning - Multiplication by zero is not an ERO.'
          call chrwrt ( iounit, output )
          ierror = 1
          return
        else if ( istop .eq. isbot ) then
          return
        end if
      else if ( iform .eq. 1 ) then
        if ( sval .eq. 0.0 ) then
          output = 'Warning - Multiplication by zero is not an ERO.'
          call chrwrt ( iounit, output )
          ierror = 1
          return
        else if ( sval .eq. 1.0 ) then
          return
        end if
      else if ( iform .eq. 2 ) then
        if ( istop .eq. 0 ) then
          output = 'Warning - Multiplication by zero is not an ERO.'
          call chrwrt ( iounit, output )
          ierror = 1
          return
        else if ( istop .eq. 1 .and. isbot.eq.0 ) then
          return
        end if
      end if
c
c  Carry out multiplication.
c
      if ( iform .eq. 0 ) then
 
        do j = 1, ncol
 
          call ratmul ( ibot, iabot(irow,j), isbot, ierror, 
     &      itop, iatop(irow,j), istop, maxint )

          if ( ierror .ne. 0 ) return
 
          iatop(irow,j) = itop
          iabot(irow,j) = ibot
 
        end do
 
        chrtmp = chlrat ( istop, isbot )
        output = 'ERO: Row '//chrint(irow)//' <= '//chrtmp//
     &    ' Row '//chrint(irow)
 
      else if ( iform .eq. 1 ) then
 
        do j = 1, ncol
          a(irow,j) = sval * a(irow,j)
        end do
 
        output = 'ERO: Row '//chrint(irow)//' <= '//chrrel(sval)//
     &    ' Row ' // chrint(irow)
 
      else if ( iform .eq. 2 ) then
 
        do j = 1, ncol
 
          call decmul ( ibot, iabot(irow,j), isbot,
     &      itop, iatop(irow,j), istop, maxint, ndig )

          if ( ierror .ne. 0 ) return
 
          iatop(irow,j) = itop
          iabot(irow,j) = ibot
 
        end do
 
        chrtmp = chldec ( istop, isbot )
        output = 'ERO: Row '//chrint(irow)//' <= '//chrtmp//
     &    ' Row ' // chrint(irow)
 
      end if
 
      call chrdb2 ( output )
      call chrwrt ( iounit, output )
 
      if ( iform .eq. 0 ) then

        call ratmul ( idebot, idebot, isbot, ierror, idetop, idetop,
     &    istop, maxint )

      else if ( iform .eq. 1 ) then

        dete = dete * sval

      else if ( iform .eq. 2 ) then

        call decmul ( idebot, idebot, isbot, idetop, idetop, istop,
     &    maxint, ndig )

      end if
 
      return
      end
      function npage ( )
c
c***********************************************************************
c
cc NPAGE determines whether it's time to pause before more printing.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer NPAGE.
c    The current number of pages completed, defined as the
c    number of lines printed, divided by the number of pages per line.
c
      integer lpage
      integer nline
      integer npage
c
      lpage = 0
      nline = 0
c
c  Get the page length.
c
      call indata ( 'GET', 'LPAGE', lpage )
 
      if ( lpage .le. 0 ) then
        npage = 0
        return
      end if
c
c  Get the current line number.
c
      call indata ( 'GET', 'NLINE', nline )
 
      npage = nline / lpage
      nline = nline - npage * lpage
      call setlin ( nline )
 
      return
      end
      subroutine pass ( filkey, iauthr, ierror, iounit, line, nline,
     &  output, prompt )
c
c***********************************************************************
c
cc PASS authenticates the password used for the automatic option.
c
c
c  The automatic option allows the user to do row elimination, or
c  solve a linear programming problem, using a single command.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*60 FILKEY, the name of the file which
c    contains the MATMAN authorization key.
c
c    Input/output, integer IAUTHR,
c    0 if the user has not yet typed the correct key,
c    1 if the user has typed the key.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Workspace, character*80 LINE.
c    Used to hold the user's input.
c
c    Input/output, integer NLINE.
c    Keeps track of the number of useful characters in LINE.
c
c    Workspace, character*100 OUTPUT.
c
c    Workspace, character*80 PROMPT.
c
      character*60 filkey
      integer iauthr
      integer ierror
      integer inew
      integer inkey
      integer iounit(4)
      integer ips
      integer key
      integer lchar
      integer lenchr
      character*80 line
      integer nline
      character*100 output
      character*80 prompt
      real random
      real x
c
c  If authorization has already been given, refuse to demand it
c  a second time.
c
      if ( iauthr .eq. 1 ) then
        return
      end if
c
c  Get the key from the user.
c
      prompt = 'authorization key for "Z" command.'
      call intrea ( inkey, line, nline, prompt, iounit, ierror )
      if ( ierror .ne. 0 ) return
c
c  If the user key is negative, this may be an attempt to change
c  the password.
c
      if ( inkey .gt. 0 ) then
        inew = 0
      else
        inew = 1
      end if
 
      inkey = abs ( inkey )
      x = random ( inkey )
c
c  Here's how to find out what your input key gets turned into.
c
c     write ( *, * ) 'RANDOM(INKEY) = ', inkey
c
c  This line works for a "private" copy of MATMAN on VAX/VMS,
c  IBM PC or Macintosh.
c
      open ( unit = 32, file = filkey, status = 'old', err = 50 )
c
c  This line works for a "shared" copy of MATMAN on VAX/VMS:
c
c     open ( unit = 32, file = filkey, status = 'old', err = 50,
c    &  shared, readonly )
c
      read ( 32, *, end = 30, err = 30 ) key
      close ( unit = 32 )
 
10    continue
 
      if ( inkey .ne. key ) then
        output = 'Authorization denied.  See your instructor for help.'
        call chrwrt ( iounit, output )
        return
      end if
 
      output = 'Authorization confirmed.'
      call chrwrt ( iounit, output )
      iauthr = 1
      if ( inew .eq. 0 ) return
      nline = 0
      prompt = '5 digit positive integer password.'
      call intrea ( ips, line, nline, prompt, iounit, ierror )
 
      if ( ierror .ne. 0 ) then
        ierror = 0
        output = 'Sorry, could not accept your new key.'
        call chrwrt ( iounit, output )
        return
      end if
 
      ips = abs ( ips )
      ips = mod ( ips, 100000 )
      inkey = ips
      x = random ( inkey )
      open ( unit = 32, file = filkey, status = 'old', err = 20 )
      close ( unit = 32, status = 'delete' )
 
20    continue
 
      open ( unit = 32, file = filkey, status = 'new', err = 60 )
      write ( 32, * ) inkey
      close ( unit = 32 )
      output = 'Password file updated.'
      call chrwrt ( iounit, output )
      return
 
30    continue
 
      close ( unit = 32 )
 
      output = 'Authorization denied.  See your instructor for help.'
      call chrwrt ( iounit, output )
      return
c
c  If the key file can't be found, use a default value.
c
50    continue
 
      output = 'The usual key file cannot be found.'
      call chrwrt ( iounit, output )
      output = 'MATMAN will use the default key.'
      call chrwrt ( iounit, output )
      key = 14897
      go to 10
c
c  The password file could not be opened as a "NEW" file.
c
60    continue
 
      lchar = lenchr(filkey)
      output = 'Problems opening the file "'//filkey(1:lchar)//'".'
      call chrdb2 ( output )
      call chrwrt ( iounit, output )
 
      return
      end
      subroutine pernex ( n, iarray, more, even )
c
c***********************************************************************
c
cc PERNEX computes all of the permutations of N objects, one at a time.
c
c
c  The routine is initialized by calling with MORE = TRUE, in which case
c  it returns the identity permutation.
c
c  If the routine is called with MORE = FALSE, then the successor of the
c  permutation in IARRAY is computed.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Reference:
c
c    Combinatorial Algorithms,
c    A Nijenhuis and H Wilf,
c    Academic Press, 1978, second edition,
c    ISBN 0-12-519260-6
c
c  Parameters:
c
c    Input, integer N, the number of objects being permuted.
c
c    Input/output, integer IARRAY(N).
c
c    If MORE is .TRUE., then IARRAY is assumed to contain the
c    "previous" permutation, and on IARRAY(I) is the value
c    of the I-th object under the next permutation.
c
c    Otherwise, IARRAY(I) will be set to the "first" permutation.
c
c    Input/output, logical MORE.
c
c    Set MORE = .FALSE. before the first call.
c
c    MORE will be reset to .TRUE. and a permutation will be returned.
c
c    Each new call produces a new permutation until
c    MORE is returned .FALSE.
c
c    Output, logical EVEN.
c
c    EVEN is .TRUE. if the output permutation is even, that is,
c    involves an even number of transpositions.
c
c    EVEN is .FALSE. otherwise.
c
      integer n
c
      integer i
      integer i1
      integer ia
      integer iarray(n)
      integer id
      integer is
      integer j
      integer l
      integer m
      logical more
      logical even
c
      if ( .not. more ) then

        do i = 1, n
          iarray(i) = i
        end do

        more = .true.
        even = .true.

        if ( n .eq. 1 ) then
          more = .false.
          return
        end if

        if ( iarray(n) .ne. 1 .or. 
     &    iarray(1) .ne. 2 + mod ( n, 2 ) ) then
          return
        end if

        do i = 1, n-3
          if ( iarray(i+1) .ne. iarray(i)+1 ) then
            return
          end if
        end do

        more = .false.

      else

        if ( n .eq. 1 ) then
          iarray(1) = 0
          more = .false.
          return
        end if

        if ( even ) then

          ia = iarray(1)
          iarray(1) = iarray(2)
          iarray(2) = ia
          even = .false.

          if ( iarray(n) .ne. 1 .or. 
     &      iarray(1) .ne. 2 + mod ( n, 2 ) ) then
            return
          end if

          do i = 1, n-3
            if ( iarray(i+1) .ne. iarray(i)+1 ) then
              return
            end if
          end do

          more = .false.
          return

        else

          is = 0

          do i1 = 2, n

            ia = iarray(i1)
            i = i1 - 1
            id = 0

            do j = 1, i
              if ( iarray(j) .gt. ia ) then
                id = id + 1
              end if
            end do

            is = id + is
            if ( id .ne. i * mod ( is, 2 ) ) then
              go to 10
            end if

          end do

          iarray(1) = 0
          more = .false.
          return

        end if

10      continue

        m = mod ( is+1, 2 ) * (n+1)

        do j = 1, i

          if ( isign(1,iarray(j)-ia) .ne. isign(1,iarray(j)-m) ) then
            m = iarray(j)
            l = j
          end if

        end do

        iarray(l) = ia
        iarray(i1) = m
        even = .true.

      end if

      return
      end
      function random ( iseed )
c
c***********************************************************************
c
cc RANDOM is a portable random number generator.
c
c
c  The recursion used has the form:
c
c    ISEED = ISEED*IA mod IP
c
c  RANDOM was extracted from ACM algorithm 570, "LOPSI".
c
c  IA = 7**5, IB=2**15, IB16=2**16, IP=2**31-1
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer ISEED, the integer "seed" used to generate
c    the output random number, and updated in preparation for the
c    next one.
c
c    Output, real RANDOM, a random value between 0 and 1.
c
      integer ia
      parameter ( ia = 16807 )
c
      integer ib15
      parameter ( ib15 = 32768 )
c
      integer ib16
      parameter ( ib16 = 65536 )
c
      integer ip
      parameter ( ip = 2147483647 )
c
      integer iprhi
      integer iseed
      integer ixhi
      integer k
      integer leftlo
      integer loxa
      real random
c
c  Get the 15 high order bits of ISEED.
c
      ixhi = iseed / ib16
c
c  Get the 16 low bits of ISEED and form the low product.
c
      loxa = ( iseed - ixhi * ib16 ) * ia
c
c  Get the 15 high order bits of the low product.
c
      leftlo = loxa / ib16
c
c  Form the 31 highest bits of the full product.
c
      iprhi = ixhi * ia + leftlo
c
c  Get overflow past the 31st bit of full product.
c
      k = iprhi / ib15
c
c  Assemble all the parts and presubtract IP.  The parentheses are
c  essential.
c
      iseed = ( ( ( loxa - leftlo * ib16 ) - ip ) +
     &  ( iprhi - k * ib15 ) * ib16 ) + k
c
c  Add IP back in if necessary.
c
      if ( iseed .lt. 0 ) then
        iseed = iseed + ip
      end if
c
c  Multiply by 1 / (2**31-1).
c
      random = real ( iseed ) * 4.656612875e-10
 
      return
      end
      subroutine ratadd ( ibot, ibot1, ibot2, ierror, itop, itop1, 
     &  itop2, maxint )
c
c***********************************************************************
c
cc RATADD adds two rational values.
c
c
c  RATADD computes
c
c    ITOP/IBOT = ITOP1/IBOT1 + ITOP2/IBOT2
c
c  while trying to avoid integer overflow.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IBOT, the denominator of the result.
c
c    Input, integer IBOT1, IBOT2, the denominators of the
c    two rational values to be added.
c
c    Output, integer IERROR.
c    0, no error occurred.
c    1, an error occurred.  The addition of the two values
c    requires a numerator or denominator larger than the
c    maximum legal integer.
c
c    Output, integer ITOP, the numerator of the result.
c
c    Input, integer ITOP1, ITOP2, the numerators of the
c    two rational values to be added.
c
c    Input, integer MAXINT, the maximum integer representable, which
c    should probably be 2147483647.
c
      integer ibot
      integer ibot1
      integer ibot2
      integer ierror
      integer igcf
      integer itemp
      integer itop
      integer itop1
      integer itop2
      integer jbot1
      integer jbot2
      integer jbot3
      integer jtop1
      integer jtop2
      integer maxint
      real temp1
      real temp2
c
      ierror = 0
 
      if ( itop1 .eq. 0 ) then
        itop = itop2
        ibot = ibot2
        return
      else if ( itop2 .eq. 0 ) then
        itop = itop1
        ibot = ibot1
        return
      end if
c
c  Make copies of the input arguments, since we will change them.
c
      jbot1 = ibot1
      jbot2 = ibot2
      jtop1 = itop1
      jtop2 = itop2
c
c  Compute the greatest common factor of the two denominators,
c  and factor it out.
c
      jbot3 = igcf ( jbot1, jbot2 )
      jbot1 = jbot1 / jbot3
      jbot2 = jbot2 / jbot3
c
c  The fraction may now be formally written as:
c
c    (jtop1*jbot2 + jtop2*jbot1) / (jbot1*jbot2*jbot3)
c
c  Check the tops for overflow.
c
      temp1 = jtop1
      temp1 = abs ( temp1 * jbot2 )
 
      temp2 = maxint
 
      if ( temp1 .gt. temp2 ) then
        ierror = 1
        write ( *, * ) ' '
        write ( *, * ) 'RATADD - Fatal error!'
        write ( *, * ) '  Overflow of top of rational sum.'
        itop = 0
        stop
      else
        jtop1 = jtop1 * jbot2
      end if
 
      temp1 = jtop2
      temp1 = abs ( temp1 * jbot1 )
 
      temp2 = maxint
 
      if ( temp1 .gt. temp2 ) then
        ierror = 1
        write ( *, * ) ' '
        write ( *, * ) 'RATADD - Fatal error!'
        write ( *, * ) '  Overflow of top of rational sum.'
        itop = 0
        stop
      else
        jtop2 = jtop2 * jbot1
      end if
 
      temp1 = jtop1
      temp1 = abs ( temp1 + jtop2 )
 
      temp2 = maxint
 
      if ( temp1 .gt. temp2 ) then
        ierror = 1
        write ( *, * ) ' '
        write ( *, * ) 'RATADD - Fatal error!'
        write ( *, * ) '  Overflow of top of rational sum.'
        itop = 0
        stop
      else
        itop = jtop1 + jtop2
      end if
c
c  Check the bottom for overflow.
c
      temp1 = jbot1
      temp1 = temp1 * jbot2
      temp1 = abs ( temp1 * jbot3 )
 
      temp2 = maxint
 
      if ( temp1 .gt. temp2 ) then
        ierror = 1
        write ( *, * ) ' '
        write ( *, * ) 'RATADD - Fatal error!'
        write ( *, * ) '  Overflow of bottom of rational sum.'
        ibot = 1
        stop
      else
        ibot = jbot1 * jbot2 * jbot3
      end if
c
c  Put the fraction in lowest terms.
c
      itemp = igcf ( itop, ibot )
      itop = itop / itemp
      ibot = ibot / itemp
c
c  Sign of bottom should be positive.
c
      if ( ibot .lt. 0 ) then
        ibot = - ibot
        itop = - itop
      end if
 
      return
      end
      subroutine ratdec ( iatop, iabot, ndig )
c
c***********************************************************************
c
cc RATDEC converts a rational value to a decimal value.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer IATOP, IABOT.
c
c    On input, the rational value (IATOP/IABOT) which is to be
c    converted.
c
c    On output, the rational decimal value IATOP * 10**IABOT.
c
c    Input, integer NDIG, the number of decimal digits used.
c
      double precision dval
      integer iabot
      integer iatop
      integer ndig
c
      dval = dble ( iatop ) / dble ( iabot )
 
      call dbldec ( dval, iatop, iabot, ndig )
 
      return
      end
      subroutine ratdet ( iatop, iabot, idtop, idbot, iarray, ierror,
     &  lda, maxint, n )
c
c***********************************************************************
c
cc RATDET finds the determinant of an N by N matrix of rational entries.
c
c
c  The brute force method is used.
c
c  RATDET should only be used for small matrices, since this calculation
c  requires the summation of N! products of N numbers.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer IATOP(LDA,N), IABOT(LDA,N), the numerators
c    and denominators of the entries of the matrix.
c
c    Output, integer IDTOP, IDBOT, the determinant of the matrix,
c    expressed as IDTOP/IDBOT.
c
c    Workspace, integer IARRAY(N).
c
c    Output, integer IERROR.
c    0, the determinant was computed.
c    1, an overflow error occurred, and the determinant was not
c    computed.
c
c    Input, integer LDA, the leading dimension of A.
c
c    Input, integer MAXINT, the maximum integer representable, which
c    should probably be 2147483647.
c
c    Input, integer N, the number of rows and columns of A.
c
      integer lda
      integer n
c
      logical even
      integer i
      integer iabot(lda,n)
      integer iatop(lda,n)
      integer iarray(n)
      integer ibot
      integer ibot1
      integer ibot2
      integer idbot
      integer idtop
      integer ierror
      integer itop
      integer itop1
      integer itop2
      integer maxint
      logical more
c
      ierror = 0
 
      more = .false.
      idtop = 0
      idbot = 1
 
10    continue
 
      call pernex ( n, iarray, more, even )
 
      if ( even ) then
        itop = 1
      else
        itop = -1
      end if
 
      ibot = 1
 
      do i = 1, n
 
        itop1 = itop
        ibot1 = ibot
        itop2 = iatop(i,iarray(i))
        ibot2 = iabot(i,iarray(i))
 
        call ratmul ( ibot, ibot1, ibot2, ierror, itop, itop1, itop2,
     &    maxint )
 
        if ( ierror .ne. 0 ) then
          write ( *, * ) ' '
          write ( *, * ) 'RATDET - Fatal error!'
          write ( *, * ) '  An overflow occurred.'
          write ( *, * ) '  The determinant calculation cannot be done'
          write ( *, * ) '  for this matrix.'
          idtop = 0
          idbot = 1
          return
        end if
 
      end do
 
      itop1 = itop
      ibot1 = ibot
 
      itop2 = idtop
      ibot2 = idbot
 
      call ratadd ( ibot, ibot1, ibot2, ierror, itop, itop1, itop2,
     &  maxint )
 
      if ( ierror .eq. 0 ) then
        idtop = itop
        idbot = ibot
      else
        write ( *, * ) ' '
        write ( *, * ) 'RATDET - Fatal error!'
        write ( *, * ) '  An overflow occurred.'
        write ( *, * ) '  The determinant calculation cannot be done'
        write ( *, * ) '  for this matrix.'
        idtop = 0
        idbot = 1
        return
      end if
 
      if ( more ) then
        go to 10
      end if
c
c  Sign of bottom should be positive.
c
      if ( idbot .lt. 0 ) then
        idbot = - idbot
        idtop = - idtop
      end if

      return
      end
      subroutine ratdiv ( ibot, ibot1, ibot2, ierror, itop, itop1,
     &  itop2, maxint )
c
c***********************************************************************
c
cc RATDIV divides one rational value by another.
c
c
c  RATDIV computes
c
c    ITOP / IBOT = ( ITOP1 / IBOT1 ) / ( ITOP2 / IBOT2 ).
c
c  while avoiding integer overflow.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IBOT, the denominator of the result.
c
c    Input, integer IBOT1, IBOT2, the denominators of the
c    two rational values.
c
c    Output, integer IERROR.
c    0, no error occurred.
c    1, an error occurred.  One of the quantities IBOT1, IBOT2,
c    or ITOP2 is zero, or the result of the division
c    requires a numerator or denominator larger than the
c    maximum legal integer.
c
c    Output, integer ITOP, the numerator of the result.
c
c    Input, integer ITOP1, ITOP2, the numerators of the
c    two rational values.
c
c    Input, integer MAXINT, the maximum integer representable, which
c    should probably be 2147483647.
c
      integer ibot
      integer ibot1
      integer ibot2
      integer ierror
      integer igcf
      integer itemp
      integer itop
      integer itop1
      integer itop2
      integer jbot1
      integer jbot2
      integer jtop1
      integer jtop2
      integer maxint
      real temp1
      real temp2
c
      ierror = 0
 
      if ( ibot1 .eq. 0 .or. itop2 .eq. 0 .or. ibot2 .eq. 0 ) then
        ierror = 1
        return
      end if
 
      if ( itop1 .eq. 0 ) then
        itop = 0
        ibot = 1
        return
      end if
c
c  Make copies of the input arguments, since we will change them.
c  Implicitly invert the divisor fraction here.  The rest of
c  the code will be a multiply operation.
c
      jbot1 = ibot1
      jbot2 = itop2
      jtop1 = itop1
      jtop2 = ibot2
c
c  Get rid of all common factors in top and bottom.
c
      itemp = igcf ( jtop1, jbot1 )
      jtop1 = jtop1 / itemp
      jbot1 = jbot1 / itemp
      itemp = igcf ( jtop1, jbot2 )
      jtop1 = jtop1 / itemp
      jbot2 = jbot2 / itemp
      itemp = igcf ( jtop2, jbot1 )
      jtop2 = jtop2 / itemp
      jbot1 = jbot1 / itemp
      itemp = igcf ( jtop2, jbot2 )
      jtop2 = jtop2 / itemp
      jbot2 = jbot2 / itemp
c
c  The fraction (ITOP1*IBOT2)/(IBOT1*ITOP2) is in lowest terms.
c
c  Check the top for overflow.
c
      temp1 = jtop1
      temp1 = abs ( temp1 * jtop2 )
 
      temp2 = maxint
 
      if ( temp1 .gt. temp2 ) then
        ierror = 1
        write ( *, * ) ' '
        write ( *, * ) 'RATDIV - Fatal error!'
        write ( *, * ) '  Overflow of top of rational fraction.'
        itop = 0
        return
      else
        itop = jtop1 * jtop2
      end if
c
c  Check the bottom IBOT1*ITOP2 for overflow.
c
      temp1 = jbot1
      temp1 = abs ( temp1 * jbot2 )
 
      temp2 = maxint
 
      if ( temp1 .gt. temp2 ) then
        ierror = 1
        write ( *, * ) ' '
        write ( *, * ) 'RATDIV - Fatal error!'
        write ( *, * ) '  Overflow of bottom of rational fraction.'
        ibot = 1
        return
      else
        ibot = jbot1 * jbot2
      end if
c
c  Sign of bottom should be positive.
c
      if ( ibot .lt. 0 ) then
        ibot = - ibot
        itop = - itop
      end if
c
c  The fraction is ITOP/IBOT with no loss of accuracy.
c
      return
      end
      subroutine ratmul ( ibot, ibot1, ibot2, ierror, itop, itop1,
     &  itop2, maxint )
c
c***********************************************************************
c
cc RATMUL multiplies two fractions.
c
c
c  RATMUL computes
c
c    ITOP/IBOT = ITOP1/IBOT1 * ITOP2/IBOT2.
c
c  while avoiding integer overflow.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IBOT, the denominator of the result.
c
c    Input, integer IBOT1, IBOT2, the denominators of the
c    two rational values to be multiplied.
c
c    Output, integer IERROR.
c    0, no error occurred.
c    1, an error occurred.  The multiplication of the two values
c    requires a numerator or denominator larger than the
c    maximum legal integer.
c
c    Output, integer ITOP, the numerator of the result.
c
c    Input, integer ITOP1, ITOP2, the numerators of the
c    two rational values to be multiplied.
c
c    Input, integer MAXINT, the maximum integer representable, which
c    should probably be 2147483647.
c
      integer ibot
      integer ibot1
      integer ibot2
      integer ierror
      integer igcf
      integer itemp
      integer itop
      integer itop1
      integer itop2
      integer jbot1
      integer jbot2
      integer jtop1
      integer jtop2
      integer maxint
      real temp1
      real temp2
c
      ierror = 0
 
      if ( itop1 .eq. 0 .or. itop2 .eq. 0 ) then
        itop = 0
        ibot = 1
        return
      end if
c
c  Make copies of the input arguments, since we will change them.
c
      jbot1 = ibot1
      jbot2 = ibot2
      jtop1 = itop1
      jtop2 = itop2
c
c  Get rid of all common factors in top and bottom.
c
      itemp = igcf ( jtop1, jbot1 )
      jtop1 = jtop1 / itemp
      jbot1 = jbot1 / itemp
      itemp = igcf ( jtop1, jbot2 )
      jtop1 = jtop1 / itemp
      jbot2 = jbot2 / itemp
      itemp = igcf ( jtop2, jbot1 )
      jtop2 = jtop2 / itemp
      jbot1 = jbot1 / itemp
      itemp = igcf ( jtop2, jbot2 )
      jtop2 = jtop2 / itemp
      jbot2 = jbot2 / itemp
c
c  The fraction (ITOP1*ITOP2)/(IBOT1*IBOT2) is in lowest terms.
c
c  Check the top ITOP1*ITOP2 for overflow.
c
      temp1 = jtop1
      temp1 = abs ( temp1 * jtop2 )
 
      temp2 = maxint
 
      if ( temp1 .gt. temp2 ) then
        ierror = 1
        write ( *, * ) ' '
        write ( *, * ) 'RATMUL - Fatal error!'
        write ( *, * ) '  Overflow of top of rational product.'
        itop = 0
        return
      else
        itop = jtop1 * jtop2
      end if
c
c  Check the bottom IBOT1*IBOT2 for overflow.
c
      temp1 = jbot1
      temp1 = abs ( temp1 * jbot2 )
 
      temp2 = maxint
 
      if ( temp1 .gt. temp2 ) then
        ierror = 1
        write ( *, * ) ' '
        write ( *, * ) 'RATMUL - Fatal error!'
        write ( *, * ) '  Overflow of bottom of rational product.'
        ibot = 1
        return
      else
        ibot = jbot1 * jbot2
      end if
c
c  Sign of bottom should be positive.
c
      if ( ibot .lt. 0 ) then
        ibot = - ibot
        itop = - itop
      end if
c
c  The fraction is ITOP/IBOT with no loss of accuracy.
c
      return
      end
      subroutine ratprn ( iatop, iabot, ibase, iounit, ihi, ilo, jhi,
     &  jlo, lpmoda, maxcol, maxrow, ncol, nrow, output, title )
c
c***********************************************************************
c
cc RATPRN prints out rational vectors or matrices.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
c    IATOP and IABOT represent the current rational or decimal
c    matrix.
c
c    Input, integer IBASE(MAXROW), keeps track of basic variables.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Input, integer IHI, ILO, the last and first rows to print.
c
c    Input, integer JHI, JLO, the last and first columns to print.
c
c    Input, integer LPMODA.
c    0, the program is in linear algebra mode.
c    1, the program is in linear programming mode.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input, integer NCOL, the number of columns in the matrix.
c
c    Input, integer NROW, the number of rows in the matrix.
c
c    Workspace, character*100 OUTPUT.
c
c    Input, character*(*) TITLE, a label for the object being printed.
c
      integer ncolum
      parameter ( ncolum = 80 )
c
      integer maxcol
      integer maxrow
c
      character*6 chrint
      character*40 fornam
      character*40 fortwo
      integer i
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer ibase(maxrow)
      integer ichi
      integer iclo
      integer ihi
      integer ilo
      integer imax
      integer imin
      integer ione
      integer iounit(4)
      integer itemp
      integer izhi
      integer izlo
      integer j
      integer jhi
      integer jlo
      integer jmax
      integer jmin
      integer kmax
      character*4 lab
      integer llab
      integer lpmoda
      integer ncol
      integer none
      integer npline
      integer nrow
      character*100 output
      character*(*) title
c
      if ( lpmoda .eq. 1 ) then
        llab = 4
      else
        llab = 0
      end if
c
c  Figure out how many rationals we can get in (NCOLUM-LLAB) columns.
c
      lab = '    '
      kmax = 3
 
      do i = ilo, ihi
        do j = jlo, jhi
 
          itemp = abs ( iatop(i,j) )
 
10        continue

          if ( itemp .ge. 10**(kmax-2) ) then
            kmax = kmax + 1
            go to 10
          end if
 
          itemp = abs ( iabot(i,j) )
 
20        continue

          if ( itemp .gt. 10**(kmax-2) ) then
            kmax = kmax + 1
            go to 20
          end if
 
        end do
      end do
 
      kmax = kmax + 1
      npline = (ncolum-llab) / kmax
c
c  Create the formats.
c
      if ( lpmoda .eq. 1 ) then
        fornam = '(a'//chrint(llab)//','//chrint(npline)//'i'
     &    //chrint(kmax)//')'
      else
        fornam = '('//chrint(npline)//'i'//chrint(kmax)//')'
      end if
 
      call chrdb1 ( fornam )
 
      if ( lpmoda .eq. 1 ) then
        fortwo='('//chrint(llab)//'x,'//chrint(npline)//'i'
     &    //chrint(kmax)//')'
      else
        fortwo='('//chrint(npline)//'i'//chrint(kmax)//')'
      end if
 
      call chrdb1 ( fortwo )
 
      do jmin = jlo, jhi, npline
 
        jmax = min ( jmin+npline-1, jhi )
        lab = '    '
c
c  Handle a column vector.
c
        if ( jlo .eq. jhi .and. ilo .ne. ihi ) then
 
          output = ' '
          call chrwrt ( iounit, output )
 
          if ( ilo .eq. 1 ) then
            output = title
            call chrwrt ( iounit, output )
            output = ' '
            call chrwrt ( iounit, output )
            output = 'Column '//chrint(jlo)//' (transposed).'
            call chrdb2 ( output )
            call chrwrt ( iounit, output )
          end if
 
          do imin = ilo, ihi, npline
 
            imax = min ( imin+npline-1, ihi )
 
            output = ' '
            call chrwrt ( iounit, output )
 
            none = 0
 
            do i = imin, imax
              if ( iabot(i,jlo) .eq. 1 ) then
                ione = 3+(i-imin+1)*kmax
                output(ione:ione) = ' '
              else
                none = 1
              end if
            end do
 
            if ( lpmoda .eq. 1 ) then
              write(output,fornam) lab, (iatop(i,jlo),i=imin,imax)
              call chrwrt ( iounit, output )
              if ( none .eq. 1 ) then
                write(output,fornam)lab,(iabot(i,jlo),i=imin,imax)
                call chrwrt ( iounit, output )
              end if
            else
              write(output,fornam)(iatop(i,jlo),i=imin,imax)
              call chrwrt ( iounit, output )
              write(output,fornam)(iabot(i,jlo),i=imin,imax)
              if ( none .eq. 1 ) then
                write(output,fornam)lab,(iabot(i,jlo),i=imin,imax)
                call chrwrt ( iounit, output )
              end if
            end if
 
          end do
 
          go to 90
 
        end if
c
c  Handle a 2D array or tableau.
c
        output = ' '
        call chrwrt ( iounit, output )
 
        if ( jmin .eq. 1 ) then
          output = title
          call chrwrt ( iounit, output )
          output = ' '
          call chrwrt ( iounit, output )
        end if
 
        if ( lpmoda .eq. 1 ) then
 
          write ( output, fortwo ) ( j, j = jmin, jmax )
 
          if ( jmin .le. ncol-1 .and. ncol-1.le.jmax ) then
            izlo = llab+((ncol-1)-jmin)*kmax+kmax-2
            izhi = izlo+2
            output(izlo:izhi)='  P'
          end if
 
          if ( jmin .le. ncol .and. ncol.le.jmax ) then
            iclo = llab+(ncol-jmin)*kmax+kmax-2
            ichi = iclo+2
            output(iclo:ichi)='  C'
          end if
 
          call chrwrt ( iounit, output )
 
          output = ' '
          call chrwrt ( iounit, output )
 
        else
 
          if ( jmin .gt. 1 .or. jmax .lt. ncol.or.
     &       ilo .gt. 1 .or. ihi .lt. nrow ) then
            output = 'Columns '//chrint(jmin)//' to '//chrint(jmax)
            call chrdb2 ( output )
            call chrwrt ( iounit, output )
            output = ' '
            call chrwrt ( iounit, output )
          end if
 
        end if
 
        do i = ilo, ihi
 
          if ( lpmoda .eq. 1 ) then
 
            if ( i .lt. nrow ) then
 
              if ( ibase(i) .lt. 10 ) then
                write ( lab, '(''X'',i1)' ) ibase(i)
              else
                write ( lab, '(''X'',i2)' ) ibase(i)
              end if
 
            else if ( i .lt. ihi ) then
              lab = 'Obj2'
            else
              lab = 'Obj '
            end if
 
            if ( maxrow .eq. 1 ) then
              lab='    '
            end if
 
          end if
 
          if ( lpmoda .eq. 1 ) then
            write(output,fornam)lab,(iatop(i,j),j=jmin,jmax)
            call chrwrt ( iounit, output )
            lab = '    '
            write(output,fornam)lab,(iabot(i,j),j=jmin,jmax)
          else
            write(output,fornam)(iatop(i,j),j=jmin,jmax)
            call chrwrt ( iounit, output )
            write(output,fornam)(iabot(i,j),j=jmin,jmax)
          end if
c
c  Delete each denominator that is 1.  If all are 1, don't
c  even print out the line.
c
          none = 0
 
          do j = jmin, jmax
 
            if ( iabot(i,j) .eq. 1 ) then
              ione = llab + (j-jmin+1) * kmax
              output(ione:ione) = ' '
            else
              none = 1
            end if
 
          end do
 
          if ( none .eq. 1 ) then
            call chrwrt ( iounit, output )
          end if
 
          if ( jmax .eq. jhi .and. i.eq.ihi ) then
          else
            output = ' '
            call chrwrt ( iounit, output )
          end if
 
        end do
 
90      continue
 
      end do
 
      return
      end
      subroutine ratrea ( itop, ibot, rval, line, nline, prompt,
     &  iounit, ierror )
c
c***********************************************************************
c
cc RATREA reads a rational value, expressed as integer, decimal or fraction.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer ITOP, IBOT, the top and bottom of the
c    fraction that was read.
c
c    Output, real RVAL, the real value that approximates ITOP/IBOT.
c
c    Workspace, character*80 LINE.
c    Used to hold the user's input.
c
c    Input/output, integer NLINE.
c    Keeps track of the number of useful characters in LINE.
c
c    Workspace, character*80 PROMPT.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Output, integer IERROR.
c    0, no error occurred.
c    1, an error occurred.
c
      integer ibot
      integer ibot1
      integer ibot2
      integer igcf
      integer ierror
      integer iounit(4)
      integer itemp
      integer itop
      integer itop1
      integer itop2
      integer lchar
      integer lenchr
      character*80 line
      integer llchar
      integer nline
      character*100 output
      character*80 prompt
      real rval
c
      itop = 0
      ibot = 1
      rval = 0
      llchar = len(line)
 
10    continue
 
      call chrinp ( ierror, iounit, line, nline, output, prompt )
      if ( ierror .ne. 0 ) return
 
      if ( nline .le. 0 ) go to 10
 
      call chrctf ( line, itop1, ibot1, ierror, lchar )
 
      if ( lchar .ge. llchar ) then
        itop = itop1
        ibot = ibot1
      else if ( line(lchar+1:lchar+1) .ne. '/' ) then
        itop = itop1
        ibot = ibot1
      else
        lchar = lchar+1
        call chrchp ( line, 1, lchar )
        call chrctf ( line, itop2, ibot2, ierror, lchar )
        itop = itop1 * ibot2
        ibot = ibot1 * itop2
      end if
 
      call chrchp ( line, 1, lchar )
c
c  Make sure fraction is in lowest terms.
c
      itemp = igcf(itop,ibot)
      itop = itop / itemp
      ibot = ibot / itemp
 
      rval = itop
      rval = rval / ibot
 
      nline = lenchr(line)
 
      return
      end
      subroutine ratrel ( a, iatop, iabot )
c
c***********************************************************************
c
cc RATREL converts rational values to real values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, real A, the value of the rational quantity.
c
c    Input, integer IATOP, IABOT, the rational quantity
c    (IATOP/IABOT) that is to be converted.
c
      real a
      integer iabot
      integer iatop
c
      if ( iabot .eq. 0 ) then
        write ( *, * ) ' '
        write ( *, * ) 'RATREL - Warning!'
        write ( *, * ) '  The input fraction to be converted had a'
        write ( *, * ) '  zero denominator.'
        a = 0.0
      else
        a = real ( iatop) / real ( iabot )
      end if
 
      return
      end
      subroutine rattrn ( iatop, iabot, maxcol, maxrow, ncol, nrow )
c
c***********************************************************************
c
cc RATTRN transposes a rational matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
c    IATOP and IABOT represent the current rational or decimal
c    matrix.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input, integer NCOL, the number of columns in the matrix.
c
c    Input, integer NROW, the number of rows in the matrix.
c
      integer maxcol
      integer maxrow
c
      integer i
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer itemp
      integer j
      integer ncol
      integer nhigh
      integer nrow
c
      nhigh = min ( maxrow, maxcol )
 
      do i = 1, nhigh
        do j = i+1, nhigh
 
          itemp = iatop(i,j)
          iatop(i,j) = iatop(j,i)
          iatop(j,i) = itemp

          itemp = iabot(i,j)
          iabot(i,j) = iabot(j,i)
          iabot(j,i) = itemp
 
        end do
      end do
c
c  Swap the dimensions of the matrix.
c
      itemp = nrow
      nrow = ncol
      ncol = itemp
 
      return
      end
      subroutine ratwrn ( iounit, maxint, output )
c
c***********************************************************************
c
cc RATWRN prints out a warning about using rational arithmetic.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Input, integer MAXINT, the maximum integer representable, which
c    should probably be 2147483647.
c
c    Workspace, character*100 OUTPUT.
c
      integer iounit(4)
      integer maxint
      character*100 output
      logical said
c
      save said
c
      data said /.false./
c
      if ( .not. said ) then
 
        output = ' '
        call chrwrt ( iounit, output )
        output = 'Note:  The representation of fractions is exact.'
        call chrwrt ( iounit, output )
        output = 'Calculations with fractions are exact.'
        call chrwrt ( iounit, output )
        output = ' '
        call chrwrt ( iounit, output )
        output = 'However, this representation will break down'
        call chrwrt ( iounit, output )
        output = 'if any numerator or denominator becomes larger'
        call chrwrt ( iounit, output )
        output = 'than the maximum legal integer: '  
        call chrwrt ( iounit, output )
        write ( output, '(i20)' ) maxint
        call chrwrt ( iounit, output )
 
        said = .true.
 
      end if
 
      return
      end
      subroutine reldec ( rval, itop, ibot, ndig )
c
c***********************************************************************
c
cc RELDEC converts a real value to a decimal fraction form.
c
c
c  RELDEC is given RVAL, and computes ITOP and IBOT, so that approximately:
c
c    RVAL = ITOP * 10 ** IBOT
c
c  In particular, only NDIG digits of RVAL are used
c  in constructing the representation.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real RVAL, the real number whose decimal
c    representation is desired.
c
c    Output, integer ITOP, IBOT, form the decimal
c    representation of RVAL, approximately.
c
c    ITOP is an integer, strictly between -10**NDIG and 10**NDIG.
c    IBOT is an integer exponent of 10.
c
c    Input, integer NDIG, the number of decimal digits used.
c
      integer ibot
      integer itop
      integer ndig
      real rtop
      real rval
      real ten1
      real ten2
c
c  Special cases.
c
      if ( rval .eq. 0.0 ) then
        itop = 0
        ibot = 0
        return
      end if
c
c  Factor RVAL = RTOP * 10**IBOT
c
      rtop = rval
      ibot = 0
c
c  Now normalize so that 10**(NDIG-1) <= ABS(RTOP) < 10**(NDIG)
c
      ten1 = 10.0**(ndig-1)
      ten2 = 10.0**ndig
 
10    continue
 
      if ( abs(rtop) .lt. ten1 ) then
        rtop = rtop * 10.0
        ibot = ibot - 1
        go to 10
      else if ( abs(rtop) .ge. ten2 ) then
        rtop = rtop / 10.0
        ibot = ibot + 1
        go to 10
      end if
c
c  ITOP is the integer part of RTOP, rounded.
c
      itop = nint ( rtop )
c
c  Now divide out any factors of ten from ITOP.
c
20    continue
 
      if ( itop .ne. 0 ) then
        if ( 10*(itop/10) .eq. itop ) then
          itop = itop / 10
          ibot = ibot + 1
          go to 20
        end if
      end if
 
      return
      end
      subroutine reldet ( a, det, iarray, lda, n )
c
c***********************************************************************
c
cc RELDET finds the determinant of a real N by N matrix.
c
c
c  The brute force method is used.
c
c  RELDET should only be used for small matrices, since this calculation
c  requires the summation of N! products of N numbers.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real A(LDA,N), the matrix whose determinant is desired.
c
c    Output, real DET, the determinant of the matrix.
c
c    Workspace, integer IARRAY(N).
c
c    Input, integer LDA, the leading dimension of A.
c
c    Input, integer N, the number of rows and columns of A.
c
      integer lda
      integer n
c
      real a(lda,n)
      real det
      logical even
      integer i
      integer iarray(n)
      logical more
      real term
c
      more = .false.
      det = 0.0
 
10    continue
 
      call pernex ( n, iarray, more, even )
 
      if ( even ) then
        term = 1.0
      else
        term = -1.0
      end if
 
      do i = 1, n
        term = term * a(i,iarray(i))
      end do
 
      det = det + term
 
      if ( more ) then
        go to 10
      end if
 
      return
      end
      subroutine relprn ( a, ibase, iounit, ihi, ilo, jhi, jlo, 
     &  lpmoda, maxcol, maxrow, ncol, nrow, output, title )
c
c***********************************************************************
c
cc RELPRN prints out real vectors and matrices.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real A(MAXROW,MAXCOL).  A is the current matrix.
c
c    Input, integer BASE(MAXROW), keeps track of basic variables.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Input, integer IHI, ILO, the last and first rows to print.
c
c    Input, integer JHI, JLO, the last and first columns to print.
c
c    Input, integer LPMODA.
c    0, the program is in linear algebra mode.
c    1, the program is in linear programming mode.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input, integer NCOL, the number of columns in the matrix.
c
c    Input, integer NROW, the number of rows in the matrix.
c
c    Workspace, character*100 OUTPUT.
c
c    Input, character*(*) TITLE, a label for the object being printed.
c
      integer ncolum
c
      parameter ( ncolum = 80 )
c
      integer maxcol
      integer maxrow
c
      real a(maxrow,maxcol)
      logical allint
      real atemp
      real btemp
      character*6 chrint
      character*40 fornam
      character*40 fortwo
      integer i
      integer ibase(maxrow)
      integer ichi
      integer iclo
      integer ihi
      integer ilo
      integer imax
      integer imin
      integer iounit(4)
      integer itemp
      integer izhi
      integer izlo
      integer j
      integer jhi
      integer jlo
      integer jmax
      integer jmin
      integer kmax
      integer kmin
      character*4 lab
      integer llab
      integer lpmoda
      integer ncol
      integer npline
      integer nrow
      character*100 output
      character*(*) title
c
      if ( lpmoda .eq. 1 ) then
        llab = 4
      else
        llab = 0
      end if
c
c  Figure out how many numbers we can fit in (NCOLUM-LLAB) columns.
c
      kmin = 10
      kmax = kmin
 
      do i = ilo, ihi
        do j = jlo, jhi
 
          atemp = abs(a(i,j))
 
10        continue
 
          if ( atemp .ge. 10.0**(kmax-kmin) ) then
            kmax = kmax + 1
            go to 10
          end if
 
        end do
      end do
 
      npline = (ncolum-llab) / kmax
c
c  Create the formats used to print out the data.
c
      allint = .true.
 
      do i = ilo, ihi
        do j = jlo, jhi
 
          atemp = a(i,j)
          itemp = 10 * int(atemp)
          btemp = itemp
          btemp = btemp / 10.0
          if ( atemp .ne. btemp ) then
            allint=.false.
          end if

        end do
      end do
 
      if ( allint ) then
c
c  If all integers, cut down KMAX, the width of each number,
c  and update NPLINE, the number of numbers we can print on one line.
c
        kmax = kmax - 7
        npline = (ncolum-llab) / kmax
 
        if ( lpmoda .eq. 1 ) then
          fornam='(a'//chrint(llab)//','//chrint(npline)//'f'
     &      //chrint(kmax)//'.0)'
        else
          fornam='('//chrint(npline)//'f'//chrint(kmax)//'.0)'
        end if
 
      else
 
        if ( lpmoda .eq. 1 ) then
          fornam='(a'//chrint(llab)//','//chrint(npline)//'f'
     &      //chrint(kmax)//'.7)'
        else
          fornam='('//chrint(npline)//'f'//chrint(kmax)//'.7)'
        end if
 
      end if
 
      call chrdb1(fornam)
 
      if ( lpmoda .eq. 1 ) then
        fortwo='('//chrint(llab)//'x,'//chrint(npline)//'i'
     &    //chrint(kmax)//')'
      else
        fortwo='('//chrint(npline)//'i'//chrint(kmax)//')'
      end if

      call chrdb1 ( fortwo )
 
      do jmin = jlo, jhi, npline
 
        jmax = min ( jmin+npline-1, jhi )
        lab = '    '
c
c  Handle a column vector.
c
        if ( jlo .eq. jhi .and. ilo .ne. ihi ) then
 
          output = ' '
          call chrwrt ( iounit, output )
 
          if ( ilo .eq. 1 ) then
            output = title
            call chrwrt ( iounit, output )
            output = 'Column '//chrint(jlo)//' transposed.'
            call chrdb2 ( output )
            call chrwrt ( iounit, output )
          end if
 
          do imin = ilo, ihi, npline
 
            imax = min ( imin+npline-1, ihi )
 
            output = ' '
            call chrwrt ( iounit, output )
 
            if ( lpmoda .eq. 0 ) then
              write(output,fornam)(a(i,jlo),i=imin,imax)
              call chrwrt ( iounit, output )
            else
              write(output,fornam)lab,(a(i,jlo),i=imin,imax)
              call chrwrt ( iounit, output )
            end if
 
          end do
 
          go to 90
        end if
 
        output = ' '
        call chrwrt ( iounit, output )
 
        if ( jmin .eq. 1 ) then
          output = title
          call chrwrt ( iounit, output )
          output = ' '
          call chrwrt ( iounit, output )
        end if
c
c  Print heading for linear programming tableau.
c
        if ( lpmoda .eq. 1 ) then
 
          write ( output, fortwo ) ( j, j = jmin, jmax )
 
          if ( jmin .le. ncol-1 .and. ncol-1.le.jmax ) then
            izlo = llab+((ncol-1)-jmin)*kmax+kmax-2
            izhi = izlo+2
            output(izlo:izhi) = '  P'
          end if
 
          if ( jmin .le. ncol .and. ncol .le. jmax ) then
            iclo = llab+(ncol-jmin)*kmax+kmax-2
            ichi = iclo+2
            output(iclo:ichi) = '  C'
          end if
 
          call chrwrt ( iounit, output )
 
          output = ' '
          call chrwrt ( iounit, output )
c
c  Print heading for linear algebra matrix.
c
        else
          if ( jmin .gt. 1 .or. jmax .lt. ncol.or.
     &       ilo .gt. 1 .or. ihi .lt. nrow ) then
 
            output = 'Columns '//chrint(jmin)//' to '//chrint(jmax)
            call chrdb2 ( output )
            call chrwrt ( iounit, output )
            output = ' '
            call chrwrt ( iounit, output )
          end if
        end if
 
        do i = ilo, ihi
 
          if ( lpmoda .eq. 1 ) then

            if ( i .lt. nrow ) then
              if ( ibase(i) .lt. 10 ) then
                write(lab,'(a1,i1)')'X',ibase(i)
              else
                write(lab,'(a1,i2)')'X',ibase(i)
              end if
            else if ( i .lt. ihi ) then
              lab='Obj2'
            else
              lab='Obj '
            end if

            if ( maxrow .eq. 1 )then
              lab='    '
            end if

          end if
 
          if ( lpmoda .eq. 1 ) then
            write ( output, fornam ) lab, ( a(i,j), j = jmin, jmax )
          else
            write ( output, fornam ) ( a(i,j), j = jmin, jmax )
          end if
          call chrwrt ( iounit, output )
 
        end do
 
90      continue
 
      end do
 
      return
      end
      subroutine relrat ( a, iatop, iabot, ndig )
c
c***********************************************************************
c
cc RELRAT converts a real value to a rational value.  
c
c
c  The rational value (IATOP/IABOT) is essentially computed by truncating 
c  the decimal representation of the real value after a given number of
c  decimal digits.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real A, the real value to be converted.
c
c    Output, integer IATOP, IABOT, the numerator and denominator
c    of the rational value that approximates A.
c
c    Input, integer NDIG, the number of decimal digits used.
c
      real a
      real factor
      integer iabot
      integer iatop
      integer ibot
      integer ifac
      integer igcf
      integer itemp
      integer itop
      integer jfac
      integer ndig
c
      factor = 10.0**ndig
 
      if ( ndig .gt. 0 ) then
        ifac = 10**ndig
        jfac = 1
      else
        ifac = 1
        jfac = 10**(-ndig)
      end if
 
      itop = nint ( a * factor ) * jfac
      ibot = ifac
c
c  Factor out the greatest common factor.
c
      itemp = igcf ( itop, ibot )
 
      iatop = itop / itemp
      iabot = ibot / itemp
 
      return
      end
      subroutine relrea ( rval, line, nline, prompt, iounit, ierror )
c
c***********************************************************************
c
cc RELREA "reads" a real value from a line of text.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c  RELREA accepts a LINE of NLINE characters which may contain some
c  user input.  If not, it prints out the PROMPT and reads new
c  information into LINE, seeking to find a real number RVAL to
c  return.
c
c  RELREA will accept integers, decimals, and ratios of the
c  form R1/R2.  Real numbers may be in scientific notation, as
c  +12.34E-56.78
c
c
c    Output, real RVAL, the real value found in LINE.
c
c    Workspace, character*80 LINE.
c    Used to hold the user's input.
c
c    Input/output, integer NLINE.
c    Keeps track of the number of useful characters in LINE.
c
c    Workspace, character*80 PROMPT.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Output, integer IERROR.
c    0, no error occurred.
c    nonzero, an error occurred while trying to read the value.
c
      real bot
      integer ierror
      integer iounit(4)
      integer lchar
      integer lenchr
      character*80 line
      integer nline
      character*100 output
      character*80 prompt
      real rval
      real top
c
      rval = 0.0
      top = 0.0
      bot = 1.0
c
c  Read a character string.
c
10    continue
 
      call chrinp ( ierror, iounit, line, nline, output, prompt )
      if ( ierror .ne. 0 ) return
      if ( nline .le. 0 ) go to 10
c
c  Convert the character string to a decimal value, TOP.
c
      call chrctr ( line, top, ierror, lchar )
c
c  If we haven't used up all our characters,
c  and if the next character is '/',
c  then the user means to input the value as a ratio,
c  so prepare to read BOT as well.
c
      if ( lchar+1 .lt. nline ) then
        if ( line(lchar+1:lchar+1) .eq. '/' ) then
          lchar = lchar + 1
          call chrchp ( line, 1, lchar )
          call chrctr ( line, bot, ierror, lchar )
          if ( bot .eq. 0.0 ) then
            bot = 1.0
          end if
        end if
      end if
c
c  Set the value of RVAL.
c
      rval = top / bot
c
c  Chop out the characters that were used.
c
      call chrchp ( line, 1, lchar )
      nline = lenchr ( line )
 
      return
      end
      subroutine reltrn ( a, maxcol, maxrow, ncol, nrow )
c
c***********************************************************************
c
cc RELTRN transposes the real matrix A.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, real A(MAXROW,MAXCOL), the matrix which is
c    transposed.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input, integer NCOL, the number of columns in the matrix.
c
c    Input, integer NROW, the number of rows in the matrix.
c
      integer maxcol
      integer maxrow
c
      real a(maxrow,maxcol)
      integer i
      integer itemp
      integer j
      integer ncol
      integer nhigh
      integer nrow
      real temp
c
      nhigh = min ( maxrow, maxcol )
 
      do i = 1, nhigh
        do j = i+1, nhigh
 
          temp = a(i,j)
          a(i,j) = a(j,i)
          a(j,i) = temp
 
        end do
      end do
c
c  Swap the dimensions of the matrix.
c
      itemp = nrow
      nrow = ncol
      ncol = itemp
 
      return
      end
      subroutine relwrn ( iounit, output )
c
c***********************************************************************
c
cc RELWRN prints out, just once, a warning about using real arithmetic.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Workspace, character*100 OUTPUT.
c
      integer iounit(4)
      character*100 output
      logical said
c
      save said
c
      data said /.false./
c
c  If real arithmetic is being used, print a warning message,
c  but only once.
c
      if ( .not. said ) then
 
        output = ' '
        call chrwrt ( iounit, output )
        output = 'Note:  Real arithmetic can be inaccurate.'
        call chrwrt ( iounit, output )
        output = ' '
        call chrwrt ( iounit, output )
        output = 'In particular, a singular matrix may be'
        call chrwrt ( iounit, output )
        output = 'incorrectly found to be nonsingular.'
        call chrwrt ( iounit, output )
 
        said = .true.
 
      end if
 
      return
      end
      subroutine restor ( a, c, iabot, iatop, ibase, ibasec, icbot,
     &  ictop, ierror, imat, iounit, lpmoda, lpmodc, maxcol, maxrow,
     &  nart, nartc, ncol, ncolc, nrow, nrowc, nslak, nslakc, nvar,
     &  nvarc, output )
c
c***********************************************************************
c
cc RESTOR restores a matrix that was saved earlier.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, real A(MAXROW,MAXCOL), the input value of C.
c
c    Input, real C(MAXROW,MAXCOL), a saved matrix.
c
c    Output, integer IABOT(MAXROW,MAXCOL), IATOP(MAXROW,MAXCOL).
c    The input value of ICTOP, ICBOT.
c
c    Output, integer IBASE(MAXROW), the input value of IBASEC.
c
c    Input, integer IBASEC(MAXROW), a saved vector to keep track
c    of basic variables.
c
c    Input, integer ICBOT(MAXROW,MAXCOL), ICTOP(MAXROW,MAXCOL).
c    A saved matrix.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IMAT.
c    0, no matrix has been defined by the user.
c    1, a matrix has been defined by the user.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Output, integer LPMODA.
c    The input value of LPMODC.
c
c    Input, integer LPMODC, a saved linear programming switch.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c         in the matrices used by MATMAN.
c
c    Output, integer NART, the input value of NARTC.
c
c    Input, integer NARTC, a saved number of artificial variables.
c
c    Output, integer NCOL, the input value of NCOLC.
c
c    Input, integer NCOLC, a saved number of columns.
c
c    Output, integer NROW, the input value of NROWC.
c
c    Input, integer NROWC, a saved number of rows.
c
c    Output, integer NSLAK, the input value of NSLAKC.
c
c    Input, integer NSLAKC, a saved number of slack variables.
c
c    Output, integer NVAR, the input value of NVARC.
c
c    Input, integer NVARC, a saved number of basic variables.
c
c    Workspace, character*100 OUTPUT.
c
      integer maxcol
      integer maxrow

      real a(maxrow,maxcol)
      real c(maxrow,maxcol)
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer ibase(maxrow)
      integer ibasec(maxrow)
      integer icbot(maxrow,maxcol)
      integer ictop(maxrow,maxcol)
      integer ierror
      integer imat
      integer iounit(4)
      integer lpmoda
      integer lpmodc
      integer lpmods
      integer nart
      integer nartc
      integer ncol
      integer ncolc
      integer nrow
      integer nrowc
      integer nslak
      integer nslakc
      integer nvar
      integer nvarc
      character*100 output

      if ( imat .eq. 0 ) then
        ierror = 1
        output = 'You must set up a matrix with the "E" command'
        call chrwrt ( iounit, output )
        output = 'before using the "R" command to restore it!'
        call chrwrt ( iounit, output )
        return
      end if
c
c  Is there a saved matrix to restore?
c
      if ( ncolc .eq. 0 ) then
        output = 'There is no saved matrix to restore!'
        call chrwrt ( iounit, output )
        ierror = 1
        return
      end if
c
c  Save a copy of the current linear programming mode.
c
      lpmods = lpmoda
c
c  Overwrite the current information by the old information.
c
      call copmat ( c, a, ictop, icbot, iatop, iabot, ibasec, ibase,
     &  lpmodc, lpmoda, maxcol, maxrow, nartc, nart, ncolc, ncol, 
     &  nrowc, nrow, nslakc, nslak, nvarc, nvar )
 
      output = 'The saved matrix has been restored.'
      call chrwrt ( iounit, output )
c
c  Print a warning if linear programming mode has been switched.
c
      if ( lpmods .ne. lpmoda ) then
        output = 'Note: The linear programming mode has been switched.'
        call chrwrt ( iounit, output )
      end if
 
      return
      end
      subroutine rowadd ( a, iatop, iabot, ierror, iform, iounit,
     &  irow1, irow2, maxcol, maxint, maxrow, ncol, ndig, output, 
     &  sval, istop, isbot )
c
c***********************************************************************
c
cc ROWADD adds a multiple of one matrix row to another.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, real A(MAXROW,MAXCOL).  A is the current matrix.
c
c    Input/output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
c    IATOP and IABOT represent the current rational or decimal
c    matrix.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IFORM, specifies the arithmetic being used.
c    0=rational, 1=real, 2=decimal.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Input, integer IROW1, the row which is to be modified.
c
c    Input, INPUT IROW2, the row which is to be multiplied by
c    a given value and added to row IROW1.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXINT, the maximum integer representable, which
c    should probably be 2147483647.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input, integer NCOL, the number of columns in the matrix.
c
c    Input, integer NDIG, the number of decimal digits used.
c
c    Workspace, character*100 OUTPUT.
c
c    Input, real SVAL, the real multiplier to use.
c
c    Input, integer ISTOP, ISBOT, the fractional or decimal
c    multiplier to use.
c
      integer maxcol
      integer maxrow
c
      real a(maxrow,maxcol)
      character*22 chldec
      character*22 chlrat
      character*6 chrint
      character*14 chrrel
      character*22 chrtmp
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer ibot
      integer ierror
      integer iform
      integer iounit(4)
      integer irow1
      integer irow2
      integer isbot
      integer isbot2
      integer istop
      integer istop2
      integer itop
      integer j
      integer maxint
      integer ncol
      integer ndig
      character*100 output
      real sval
c
      if ( iform .eq. 0 ) then
 
        if ( istop .eq. 0 ) return
 
        do j = 1, ncol
 
          call ratmul ( isbot2, isbot, iabot(irow2,j), ierror,
     &      istop2, istop, iatop(irow2,j), maxint )
 
          call ratadd ( ibot, iabot(irow1,j), isbot2, ierror, itop,
     &      iatop(irow1,j), istop2, maxint )
 
          iatop(irow1,j) = itop
          iabot(irow1,j) = ibot
 
        end do
 
      else if ( iform .eq. 1 ) then
 
        if ( sval .eq. 0.0 ) return
 
        do j = 1, ncol
          a(irow1,j) = a(irow1,j) + sval * a(irow2,j)
        end do
 
      else if ( iform .eq. 2 ) then
 
        if ( istop .eq. 0 ) return
 
        do j = 1, ncol
 
          call decmul ( isbot2, isbot, iabot(irow2,j), istop2, istop,
     &      iatop(irow2,j), maxint, ndig )
 
          call decadd ( ibot, iabot(irow1,j), isbot2, itop,
     &      iatop(irow1,j), istop2, ndig )
 
          iatop(irow1,j) = itop
          iabot(irow1,j) = ibot
 
        end do
 
      end if
 
      if ( iform .eq. 0 ) then
 
        chrtmp = chlrat ( istop, isbot )
 
      else if ( iform .eq. 1 ) then
 
        chrtmp = chrrel ( sval )
 
      else if ( iform .eq. 2 ) then
 
        chrtmp = chldec ( istop, isbot )
 
      end if
 
      output = 'ERO: Row ' // chrint(irow1) // ' <= '
     &  //chrtmp//' Row '//chrint(irow2)//' + Row '//chrint(irow1)
      call chrdb2 ( output )
      call chrwrt ( iounit, output )
 
      return
      end
      subroutine sample ( a, chineq, iatop, iabot, ibase, ierror,
     &  iform, imat, iounit, iseed, line, lpmoda, maxcol, maxrow, nart, 
     &  ncol, nline, nrow, nslak, nvar, output, prompt )
c
c***********************************************************************
c
cc SAMPLE allows the user to choose a particular sample problem.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, real A(MAXROW,MAXCOL).  A is the current matrix.
c
c    Output, character*1 CHINEQ(MAXROW), the '<', '=', or '>'
c    sign for each linear programming constraint.
c
c    Output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
c    IATOP and IABOT represent the current rational or decimal
c    matrix.
c
c    Output, integer IBASE(MAXROW), keeps track of basic variables.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Output, integer IFORM, specifies the arithmetic being used.
c    0=rational, 1=real, 2=decimal.
c
c    Input, integer IMAT.
c    0, no matrix has been defined by the user.
c    1, a matrix has been defined by the user.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Input/output, integer ISEED, a random number generator seed.
c
c    Workspace, character*80 LINE.
c    Used to hold the user's input.
c
c    Input, integer LPMODA.
c    0, the program is in linear algebra mode.
c    1, the program is in linear programming mode.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Output, integer NART, the number of artificial variables.
c
c    Output, integer NCOL, the number of columns in the matrix.
c
c    Input/output, integer NLINE.
c    Keeps track of the number of useful characters in LINE.
c
c    Output, integer NROW, the number of rows in the matrix.
c
c    Output, integer NSLAK, the number of slack variables.
c
c    Output, integer NVAR, the number of basic variables.
c
c    Workspace, character*100 OUTPUT.
c
c    Workspace, character*80 PROMPT.
c
      integer maxcol
      integer maxrow
c
      real a(maxrow,maxcol)
      character*1 chineq(maxrow)
      character*6 chrint
      integer i
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer ibase(maxrow)
      integer ierror
      integer iform
      integer ihi
      integer ilo
      integer imat
      integer iounit(4)
      character*1 isay
      integer iseed
      integer iterm
      integer ival
      integer ival2
      integer j
      logical leqi
      character*80 line
      integer lpmoda
      integer nart
      integer ncol
      integer nline
      integer nrow
      integer nslak
      integer nvar
      character*100 output
      character*80 prompt
c
      if ( lpmoda .eq. 0 ) then
 
        output = ' '
        call chrwrt ( iounit, output )
        output = 'The following examples are available:'
        call chrwrt ( iounit, output )
        output = '  "D" for determinant;'
        call chrwrt ( iounit, output )
        output = '  "E" for eigenvalues;'
        call chrwrt ( iounit, output )
        output = '  "I" for inverse;'
        call chrwrt ( iounit, output )
        output = '  "S" for linear solve.'
        call chrwrt ( iounit, output )
        output = ' '
        call chrwrt ( iounit, output )
        output = '  "C" to cancel.'
        call chrwrt ( iounit, output )
        output = ' '
        call chrwrt ( iounit, output )
 
        prompt = '"D", "E", "I", "S" or "C" to cancel.'
        iterm = 0
        call chrrea ( isay, line, nline, prompt, iounit, ierror,
     &    iterm )
 
        if ( ierror .ne. 0 ) return
c
c  Linear System problem, random square matrix, with RHS vector appended.
c
        if ( leqi ( isay, 'S' ) ) then
 
10        continue
 
          prompt = 'number of rows desired.'
 
          call intrea ( nrow, line, nline, prompt, iounit, ierror )
          if ( ierror .ne. 0 ) return
 
          ncol = nrow + 1

          if ( nrow .lt. 1 ) then
            output = 'Error!  Negative number of rows not allowed!'
            call chrwrt ( iounit, output )
            nline = 0
            go to 10
          else if ( nrow .gt. maxrow ) then
            output = 'Number of rows must be less than ' 
     &        // chrint(maxrow)
            call chrdb2 ( output )
            call chrwrt ( iounit, output )
            nline = 0
            go to 10
          else if ( ncol .gt. maxcol ) then
            output = 'Please ask for fewer rows NROW, so that '
            call chrwrt ( iounit, output )
            output = 'NROW+1 is no more than ' // chrint(maxcol)
            call chrdb2 ( output )
            nline = 0
            go to 10
          end if
 
          ilo = -10
          ihi = 10

          do i = 1, nrow
            do j = 1, nrow

              call intran ( ival, ilo, ihi, iseed )

              if ( iform .eq. 0 ) then
                iatop(i,j) = ival
                iabot(i,j) = 1
              else if ( iform .eq. 1 ) then
                a(i,j) = real ( ival )
              else if ( iform .eq. 2 ) then
                iatop(i,j) = ival
                iabot(i,j) = 0
              end if

            end do
          end do

          do i = 1, nrow

            call intran ( ival, ilo, ihi, iseed )
            ibase(i) = ival

          end do

          do i = 1, nrow

            ival = 0
            do j = 1, nrow
              ival = ival + iatop(i,j) * ibase(j)
            end do

            if ( iform .eq. 0 ) then
              iatop(i,nrow+1) = ival
              iabot(i,nrow+1) = 1
            else if ( iform .eq. 1 ) then
              a(i,nrow+1) = real ( ival )
            else if ( iform .eq. 2 ) then
              iatop(i,nrow+1) = ival
              iabot(i,nrow+1) = 0
            end if
 
          end do

          do i = 1, nrow
            ibase(i) = 0
          end do
 
          imat = 1
c
c  Inverse problem, random square matrix, with identity appended.
c
        else if ( leqi ( isay, 'I' ) ) then
 
20        continue
 
          nrow = 0
          ncol = 0
 
          prompt = 'number of rows desired.'
 
          call intrea ( nrow, line, nline, prompt, iounit, ierror )
          if ( ierror .ne. 0 ) return
 
          if ( nrow .lt. 1 ) then
            output = 'Error!  Negative number of rows not allowed!'
            call chrwrt ( iounit, output )
            nline = 0
            go to 20
          else if ( nrow .gt. maxrow ) then
            output = 'Number of rows must be less than ' 
     &        // chrint(maxrow)
            call chrdb2 ( output )
            call chrwrt ( iounit, output )
            nline = 0
            go to 20
          else if ( 2*nrow .gt. maxcol ) then
            output = 'Please ask for fewer rows NROW, so that '
            call chrwrt ( iounit, output )
            output = '2 * NROW is no more than ' // chrint(maxcol)
            call chrdb2 ( output )
            nline = 0
            go to 20
          end if
 
          ncol = 2 * nrow
 
          ilo = -10
          ihi = 10

          do i = 1, nrow
            do j = 1, nrow

              call intran ( ival, ilo, ihi, iseed )

              if ( i .eq. j ) then
                ival2 = 1
              else
                ival2 = 0
              end if

              if ( iform .eq. 0 ) then
                iatop(i,j) = ival
                iabot(i,j) = 1
                iatop(i,j+nrow) = ival2
                iabot(i,j+nrow) = 1
              else if ( iform .eq. 1 ) then
                a(i,j) = real ( ival )
                a(i,j+nrow) = real ( ival2 )
              else if ( iform .eq. 2 ) then
                iatop(i,j) = ival
                iabot(i,j) = 0
                iatop(i,j+nrow) = ival2
                iabot(i,j+nrow) = 0
              end if

            end do
          end do
 
          imat = 1
c
c  Determinant problem, a random square matrix.
c
        else if ( leqi ( isay, 'D' ) ) then
 
30        continue
 
          nrow = 0
          ncol = 0
 
          prompt = 'number of rows desired.'
 
          call intrea ( nrow, line, nline, prompt, iounit, ierror )
          if ( ierror .ne. 0 ) return
 
          if ( nrow .lt. 1 ) then
            output = 'Error!  Negative number of rows not allowed!'
            call chrwrt ( iounit, output )
            nline = 0
            go to 30
          else if ( nrow .gt. maxrow ) then
            output = 'Number of rows must be less than ' 
     &        // chrint(maxrow)
            call chrdb2 ( output )
            call chrwrt ( iounit, output )
            nline = 0
            go to 30
          end if
 
          ncol = nrow
 
          ilo = -10
          ihi = 10

          do i = 1, nrow
            do j = 1, ncol

              call intran ( ival, ilo, ihi, iseed )

              if ( iform .eq. 0 ) then
                iatop(i,j) = ival
                iabot(i,j) = 1
              else if ( iform .eq. 1 ) then
                a(i,j) = real ( ival )
              else if ( iform .eq. 2 ) then
                iatop(i,j) = ival
                iabot(i,j) = 0
              end if

            end do
          end do
 
          imat = 1
c
c  Eigenvalue sample problem, a random square symmetric matrix.
c
        else if ( leqi ( isay, 'E' ) ) then
 
31        continue
 
          nrow = 0
          ncol = 0
 
          prompt = 'number of rows desired.'
 
          call intrea ( nrow, line, nline, prompt, iounit, ierror )
          if ( ierror .ne. 0 ) return
 
          if ( nrow .lt. 1 ) then
            output = 'Error!  Negative number of rows not allowed!'
            call chrwrt ( iounit, output )
            nline = 0
            go to 31
          else if ( nrow .gt. maxrow ) then
            output = 'Number of rows must be less than ' //
     &        chrint(maxrow)
            call chrdb2 ( output )
            call chrwrt ( iounit, output )
            nline = 0
            go to 31
          end if
 
          ncol = nrow
 
          ilo = -10
          ihi = 10

          do i = 1, nrow
            do j = i, ncol

              call intran ( ival, ilo, ihi, iseed )

              if ( iform .eq. 0 ) then
                iatop(i,j) = ival
                iabot(i,j) = 1
                iatop(j,i) = ival
                iabot(j,i) = 1
              else if ( iform .eq. 1 ) then
                a(i,j) = real ( ival )
                a(j,i) = real ( ival )
              else if ( iform .eq. 2 ) then
                iatop(i,j) = ival
                iabot(i,j) = 0
                iatop(j,i) = ival
                iabot(j,i) = 0
              end if

            end do
          end do
 
          imat = 1
 
          output =
     &      'We have a symmetric matrix, whose eigenvalues'
          call chrwrt ( iounit, output )
          output =
     &      'can be found with the "J" command.'
          call chrwrt ( iounit, output )
          output = 'Be sure to use real arithmeticc'
          call chrwrt ( iounit, output )

        else
 
          output = 'No problem was selected.'
          call chrwrt ( iounit, output )
 
        end if
c
c  Linear programming.
c
      else
 
        output = ' '
        call chrwrt ( iounit, output )
        output = 'The following examples are available:'
        call chrwrt ( iounit, output )
        output = '  "S" a simple linear programming problem;'
        call chrwrt ( iounit, output )
        output = '  "A" an advanced linear programming problem.'
        call chrwrt ( iounit, output )
        output = ' '
        call chrwrt ( iounit, output )
        output = '  "C" to cancel.'
        call chrwrt ( iounit, output )
        output = ' '
        call chrwrt ( iounit, output )
 
        prompt = '"S", "A", or "C" to cancel.'
        iterm = 0
        call chrrea ( isay, line, nline, prompt, iounit, ierror,
     &    iterm )
 
        if ( ierror .ne. 0 ) then
          return
        end if
 
        if ( leqi ( isay, 'S' ) ) then

          call lpsams ( a, chineq, iatop, iabot, ibase, iform, imat,
     &      iounit, maxcol, maxrow, nart, ncol, nrow, nslak, nvar,
     &      output)

        else if ( leqi ( isay, 'A' ) ) then

          call lpsama ( a, chineq, iatop, iabot, ibase, iform, imat,
     &      iounit, maxcol, maxrow, nart, ncol, nrow, nslak, nvar,
     &      output ) 

        else
 
          output = 'No problem was selected.'
          call chrwrt ( iounit, output )
 
        end if
 
      end if
 
      return
      end
      subroutine scadiv ( a, iatop, iabot, ierror, iform, iounit, irow,
     &  maxcol, maxint, maxrow, ncol, ndig, nrow, output, sval, istop, 
     &  isbot )

c***********************************************************************
c
cc SCADIV divides row IROW of the A matrix by a scale factor.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, real A(MAXROW,MAXCOL).  A is the current matrix.
c
c    Input/output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
c    IATOP and IABOT represent the current rational or decimal
c    matrix.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IFORM, specifies the arithmetic being used.
c    0=rational, 1=real, 2=decimal.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Input, integer IROW, the row to be divided.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXINT, the maximum integer representable, which
c    should probably be 2147483647.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input, integer NCOL, the number of columns in the matrix.
c
c    Input, integer NDIG, the number of decimal digits used.
c
c    Input, integer NROW, the number of rows in the matrix.
c
c    Workspace, character*100 OUTPUT.
c
c    Input, real SVAL, the real divisor.
c
c    Input, integer ISTOP, ISBOT, the fractional or decimal
c    divisor.
c
      integer maxcol
      integer maxrow

      real a(maxrow,maxcol)
      character*22 chldec
      character*6 chlint
      character*22 chlrat
      character*6 chrint
      character*14 chrrel
      character*24 chrtmp
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer ibot
      integer ierror
      integer iform
      integer iounit(4)
      integer irow
      integer isbot
      integer istop
      integer itop
      integer j
      integer maxint
      integer ncol
      integer ndig
      integer nrow
      character*100 output
      real sval
c
c  Make sure that the row number is legal.
c
      if ( irow .lt. 1 .or. irow .gt. nrow ) then
        output = 'Error!  The row number is out of range!'
        call chrwrt ( iounit, output )
        ierror = 1
        return
      end if
c
c  Check for an illegal divisor of 0, or a pointless divisor of 1.
c
      if ( iform .eq. 0 ) then
        if ( istop .eq. 0 ) then
          output = 'Error!  It is illegal to divide by 0!'
          call chrwrt ( iounit, output )
          ierror = 1
          return
        else if ( istop .eq. isbot ) then
          return
        end if
      else if ( iform .eq. 1 ) then
        if ( sval .eq. 0.0 ) then
          output = 'Error!  It is illegal to divide by 0!'
          call chrwrt ( iounit, output )
          ierror = 1
          return
        else if ( sval .eq. 1.0 ) then
          return
        end if
      else if ( iform .eq. 2 ) then
        if ( istop .eq. 0 ) then
          output = 'Error!  It is illegal to divide by 0!'
          call chrwrt ( iounit, output )
          ierror = 1
          return
        else if ( istop .eq. 1 .and. isbot.eq.0 ) then
          return
        end if
      end if
c
c  Carry out the division.
c
      if ( iform .eq. 0 ) then
 
        do j = 1, ncol
 
          call ratdiv ( ibot, iabot(irow,j), isbot, ierror,
     &      itop, iatop(irow,j), istop, maxint )
 
          if ( ierror .ne. 0 ) return
 
          iatop(irow,j) = itop
          iabot(irow,j) = ibot
 
        end do
 
      else if ( iform .eq. 1 ) then
 
        do j = 1, ncol
          a(irow,j) = a(irow,j) / sval
        end do
 
      else if ( iform .eq. 2 ) then
 
        do j = 1, ncol
 
          call decdiv ( ibot, iabot(irow,j), isbot, ierror, itop,
     &      iatop(irow,j), istop, ndig )
 
          iatop(irow,j) = itop
          iabot(irow,j) = ibot
 
        end do
 
      end if
c
c  Print out a statement about what has been done.
c
      if ( iform .eq. 0 ) then
 
        if ( isbot .eq. 1 ) then
          chrtmp = chlint(istop)
        else
          chrtmp = '(' // chlrat(istop,isbot) // ')'
          call chrdb1(chrtmp)
        end if
 
      else if ( iform .eq. 1 ) then
 
        chrtmp = chrrel(sval)
 
      else if ( iform .eq. 2 ) then
 
        chrtmp = chldec(istop,isbot)
        call chrdb1(chrtmp)
 
      end if
 
      output = 'ERO: Row ' // chrint(irow) // ' <=  Row '
     &  // chrint(irow) // ' / ' // chrtmp
      call chrdb2 ( output )
      call chrwrt ( iounit, output )
 
      return
      end
      subroutine setdig ( ierror, iounit, line, maxdig, ndig, nline,
     &  output )

c***********************************************************************
c
cc SETDIG allows the user to specify the maximum digits in a decimal.
c
c  Discussion:
c
c    NDIG is
c
c    * the number of digits used when converting a real number
c      to a fraction using the "FI" or "FD" command;
c
c    * the maximum number of digits in a decimal.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Workspace, character*80 LINE.
c    Used to hold the user's input.
c
c    Input, integer MAXDIG, the maximum number of decimal
c    digits allowed.
c
c    Output, integer NDIG, the number of decimal digits to use
c    in constructing decimal representations.  NDIG should normally
c    be between 1 and 7.
c
c    Input/output, integer NLINE.
c    Keeps track of the number of useful characters in LINE.
c
c    Workspace, character*100 OUTPUT.
c
      character*6 chrint
      integer ierror
      integer iounit(4)
      integer itemp
      character*80 line
      integer maxdig
      integer ndig
      integer nline
      character*100 output
      character*80 prompt

      output = 'How many decimal places should be used in '
      call chrwrt ( iounit, output )
      output = 'converting real results to a decimal?'
      call chrwrt ( iounit, output )
      output = ' '
      call chrwrt ( iounit, output )
      output = ' 1 means 123.45 becomes 1 * 10**2'
      call chrwrt ( iounit, output )
      output = ' 2 means 123.45 becomes 12 * 10**1'
      call chrwrt ( iounit, output )
      output = ' 3 means 123.45 becomes 123'
      call chrwrt ( iounit, output )
      output = 'and so on.'
      call chrwrt ( iounit, output )
 
      prompt = 'number of decimals (1 to '//chrint(maxdig)//').'
      call chrdb2 ( prompt )
 
      call intrea ( itemp, line, nline, prompt, iounit, ierror )
 
      if ( ierror .ne. 0 ) then
        output = 'Your choice was not acceptable!'
        call chrwrt ( iounit, output )
        return
      end if
c
c  Absolutely do not let NDIG be less than 1.
c
      if ( itemp .lt. 1 ) then
        output = 'The number of decimals must be positive!'
        call chrwrt ( iounit, output )
        ierror = 1
        return
      end if
c
c  Allow user to exceed MAXDIG, with a warning.
c
      if ( itemp .gt. maxdig ) then
        output = ' '
        call chrwrt ( iounit, output )
        output = 'Warning!'
        call chrwrt ( iounit, output )
        output = 'Your choice is larger than the recommended maximum!'
        call chrwrt ( iounit, output )
        output = 'which is '//chrint(maxdig)
        call chrdb2 ( output )
        call chrwrt ( iounit, output )
        output = 'It is possible that calculations will break down'
        call chrwrt ( iounit, output )
        output = 'at any time!  Be careful!'
        call chrwrt ( iounit, output )
      end if
 
      output = ' '
      call chrwrt ( iounit, output )
      ndig=itemp
      output = 'The number of decimal digits will now be '//chrint(ndig)
      call chrdb2 ( output )
      call chrwrt ( iounit, output )
 
      return
      end
      subroutine setlin ( nline )

c***********************************************************************
c
cc SETLIN sets the current output page line number.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer NLINE, the current line number.
c
      integer nline

      call indata ( 'SET', 'NLINE', nline )
 
      return
      end
      subroutine setpag ( lpage )

c***********************************************************************
c
cc SETPAG sets the number of lines per output page.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer LPAGE, the desired number of lines per page.
c
      integer lpage

      call indata ( 'SET', 'LPAGE', lpage )
 
      return
      end
      subroutine shfcol ( a, iabot, iatop, icol, maxcol, maxrow, ncol,
     &  nrow )

c***********************************************************************
c
cc SHFCOL allows a new column to be inserted by shifting others right.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, real A(MAXROW,MAXCOL).  A is the current matrix.
c
c    Input/output, integer IABOT(MAXROW,MAXCOL), IATOP(MAXROW,MAXCOL).
c    IATOP and IABOT represent the current rational or decimal
c    matrix.
c
c    Input, integer ICOL, the position of the new column.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input, integer NCOL, the number of columns in the matrix.
c
c    Input, integer NROW, the number of rows in the matrix.
c
      integer maxcol
      integer maxrow

      real a(maxrow,maxcol)
      integer i
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer icol
      integer j
      integer ncol
      integer nrow

      do j = ncol, icol+1, -1
        do i = 1, nrow
 
          a(i,j) = a(i,j-1)
          iatop(i,j) = iatop(i,j-1)
          iabot(i,j) = iabot(i,j-1)
 
        end do
      end do
 
      return
      end
      subroutine shfrow ( a, iabot, iatop, irow, maxcol, maxrow, ncol,
     &  nrow )
c
c*********************************************************************72
c
cc SHFROW allows a new row to be inserted by shifting other rows down.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, real A(MAXROW,MAXCOL).  A is the current matrix.
c
c    Input/output, integer IABOT(MAXROW,MAXCOL), IATOP(MAXROW,MAXCOL).
c    IATOP and IABOT represent the current rational or decimal matrix.
c
c    Input, integer IROW, the position of the new row.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input, integer NCOL, the number of columns in the matrix.
c
c    Input, integer NROW, the number of rows in the matrix.
c
      implicit none

      integer maxcol
      integer maxrow

      real a(maxrow,maxcol)
      integer i
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer irow
      integer j
      integer ncol
      integer nrow

      do i = nrow, irow+1, -1
        do j = 1, ncol
 
          a(i,j) = a(i-1,j)
          iatop(i,j) = iatop(i-1,j)
          iabot(i,j) = iabot(i-1,j)
 
        end do
      end do
 
      return
      end
      subroutine swprow ( a, iatop, iabot, ibase, ierror, iform, 
     &  iounit, irow1, irow2, lpmoda, maxcol, maxrow, ncol, nrow,
     &  output )

c*********************************************************************72
c
cc SWPROW swaps two rows of a matrix.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, real A(MAXROW,MAXCOL).  A is the current matrix.
c
c    Input/output, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
c    IATOP and IABOT represent the current rational or decimal matrix.
c
c    Input/output, integer IBASE(MAXROW), keeps track of basic
c    variables.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IFORM, specifies the arithmetic being used.
c    0=rational, 1=real, 2=decimal.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Input, integer IROW1, IROW2, the numbers of the two rows
c    to be swapped.
c
c    Input, integer LPMODA.
c    0, the program is in linear algebra mode.
c    1, the program is in linear programming mode.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input, integer NCOL, the number of columns in the matrix.
c
c    Input, integer NROW, the number of rows in the matrix.
c
c    Workspace, character*100 OUTPUT.
c
      implicit none

      integer maxcol
      integer maxrow

      real a(maxrow,maxcol)
      character*6 chrint
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer ibase(maxrow)
      integer ierror
      integer iform
      integer iounit(4)
      integer irow1
      integer irow2
      integer itemp
      integer j
      integer lpmoda
      integer ncol
      integer nrow
      character*100 output
      real temp
c
c  Skip out if the two rows are the same.
c
      if ( irow1 .eq. irow2 ) then
        output = 'You have asked to swap a row with itself!'
        call chrwrt ( iounit, output )
        return
      end if
c
c  Refuse to continue if a row is out of bounds.
c
      if ( (irow1 .lt. 1 .or. irow1 .gt. nrow).or.
     &   (irow2 .lt. 1 .or. irow2 .gt. nrow) ) then
        ierror = 1
        output = 'One of the rows is illegal!'
        call chrwrt ( iounit, output )
        return
      end if
c
c  Refuse to swap the last row in linear programming mode.
c
      if ( lpmoda .eq. 1 ) then
        if ( irow1 .eq. nrow .or. irow2.eq.nrow ) then
          ierror = 1
          output = 'You are in linear programming mode.'
          call chrwrt ( iounit, output )
          output = 'You may not swap the last row!'
          call chrwrt ( iounit, output )
          return
        end if
      end if
c
c  Swap the rows.
c
      do j = 1, ncol
 
        if ( iform .eq. 0 ) then
 
          itemp = iatop(irow1,j)
          iatop(irow1,j) = iatop(irow2,j)
          iatop(irow2,j) = itemp

          itemp = iabot(irow1,j)
          iabot(irow1,j) = iabot(irow2,j)
          iabot(irow2,j)=itemp
 
        else if ( iform .eq. 1 ) then
 
          temp = a(irow1,j)
          a(irow1,j) = a(irow2,j)
          a(irow2,j) = temp
 
        else if ( iform .eq. 2 ) then
 
          itemp = iatop(irow1,j)
          iatop(irow1,j) = iatop(irow2,j)
          iatop(irow2,j) = itemp

          itemp = iabot(irow1,j)
          iabot(irow1,j) = iabot(irow2,j)
          iabot(irow2,j) = itemp
 
        end if
 
      end do
 
      itemp = ibase(irow1)
      ibase(irow1) = ibase(irow2)
      ibase(irow2) = itemp
 
      output = 'ERO: Row '//chrint(irow1)//' <=> Row '//chrint(irow2)
      call chrdb2 ( output )
      call chrwrt ( iounit, output )
 
      return
      end
      subroutine transc ( filtrn, ierror, iounit, line, nline,
     &  output, prompt )

c*********************************************************************72
c
cc TRANSC opens or closes a transcript file.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character*60 FILTRN.
c    On input, FILTRN is the current or default transcript file.
c    On output, FILTRN is the file name chosen by the user.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Workspace, character*80 LINE.
c    Used to hold the user's input.
c
c    Input/output, integer NLINE.
c    Keeps track of the number of useful characters in LINE.
c
c    Workspace, character*100 OUTPUT.
c
c    Workspace, character*80 PROMPT.
c
      implicit none

      character*60 filnam
      character*60 filtrn
      integer ierror
      integer iounit(4)
      integer iosave
      integer iterm
      integer lchar
      integer lenchr
      character*80 line
      integer nline
      character*100 output
      character*80 prompt
c
c  Get the name of the file.
c
      if ( iounit(3) .eq. -1 ) then
        lchar = lenchr(filtrn)
        prompt = 'file name, default= "'//filtrn(1:lchar)//'".'
        call chrdb2 ( prompt )
        iterm = 0
        call chrrea ( filnam, line, nline, prompt, iounit, ierror,
     &    iterm )
        if ( ierror .ne. 0 ) return

        if ( filnam .ne. ' ' ) then
          filtrn = filnam
        end if

        iounit(3) = 21
c
c  This command works for non-VMS systems.
c
        open ( unit = iounit(3), file = filtrn, status = 'new',
     &    form = 'formatted', err = 10 )
c
c  This command is preferable for VMS systems.
c
c       open ( unit = iounit(3), file = filtrn, status = 'new',
c    &    form = 'formatted', carriagecontrol = 'list', err = 10 )
 
        go to 20
c
c  Opening with STATUS='NEW' failed.  Try opening with STATUS='OLD'.
c
10      continue

        ierror = 1
        iounit(3) = - 1
        open ( unit = 21, file = filtrn, status = 'old', err = 30 )
        close ( unit = 21, status = 'delete' )
        ierror = 0
        iounit(3) = 21
c
c  This command works for non-VMS systems.
c
        open ( unit = iounit(3), file = filtrn, status = 'new',
     &    form = 'formatted' )
c
c  This command is preferrable for VMS systems.
c
c       open ( unit = iounit(3), file = filtrn, status = 'new',
c    &    form = 'formatted', carriagecontrol = 'list' )
 
20      continue

        lchar = lenchr(filtrn)
        output = 'Opening the transcript file "'//filtrn(1:lchar)//'".'
        call chrwrt ( iounit, output )
 
        iosave = iounit(2)
        iounit(2) = - 1
        call hello ( iounit, output )
        iounit(2) = iosave
 
      else
 
        lchar = lenchr(filtrn)
        output = 'Closing the transcript file "'//filtrn(1:lchar)//'".'
        call chrwrt ( iounit, output )
        close ( unit = iounit(3) )
        iounit(3) = - 1
 
      end if
 
      return
 
30    continue

      ierror = 1
      output = 'Unable to open transcript file!'
      call chrwrt ( iounit, output )
 
      return
      end
      subroutine type ( a, iabot, iatop, ibase, ierror, iform, imat,
     &  iounit, lpmoda, maxcol, maxrow, nart, ncol, nrow, output )

c*********************************************************************72
c
cc TYPE prints out the matrix or tableau.
c
c  Discussion:
c
c    It also can print out the linear programming solution.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real A(MAXROW,MAXCOL).  A is the current matrix.
c
c    Input, integer IABOT(MAXROW,MAXCOL), IATOP(MAXROW,MAXCOL).
c    IATOP and IABOT represent the current rational or decimal matrix.
c
c    Input, integer IBASE(MAXROW), keeps track of basic variables.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IFORM, specifies the arithmetic being used.
c    0=rational, 1=real, 2=decimal.
c
c    Input, integer IMAT.
c    0, no matrix has been defined by the user.
c    1, a matrix has been defined by the user.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Input, integer LPMODA.
c    0, the program is in linear algebra mode.
c    1, the program is in linear programming mode.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input, integer NART, the number of artificial variables.
c
c    Input, integer NCOL, the number of columns in the matrix.
c
c    Input, integer NROW, the number of rows in the matrix.
c
c    Workspace, character*100 OUTPUT.
c
      implicit none

      integer maxcol
      integer maxrow

      real a(maxrow,maxcol)
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer ibase(maxrow)
      integer ierror
      integer iform
      integer ihi
      integer ilo
      integer imat
      integer iounit(4)
      integer jhi
      integer jlo
      integer lpmoda
      integer nart
      integer ncol
      integer nrow
      character*100 output
      character*80 title
c
c  Make sure there is something to print.
c
      if ( imat .ne. 1 ) then
        ierror = 1
        output = ' '
        call chrwrt ( iounit, output )
        output = 'Error!  There is no data to print out.'
        call chrwrt ( iounit, output )
        return
      end if
 
      ilo = 1
      ihi = nrow
      jlo = 1
      jhi = ncol
 
      if ( lpmoda .eq. 1 .and. nart .gt. 0 .and. ihi .eq. nrow ) then
        ihi = nrow + 1
      end if
 
      if ( lpmoda .eq. 0 ) then
        title = 'The current matrix:'
      else if ( lpmoda .eq. 1 ) then
        title = 'The linear programming tableau:'
      else
        title = ' '
      end if
 
      if ( iform .eq. 0 ) then
 
        call ratprn ( iatop, iabot, ibase, iounit, ihi, ilo, jhi, jlo,
     &    lpmoda, maxcol, maxrow, ncol, nrow, output, title )
 
      else if ( iform .eq. 1 ) then
 
        call relprn ( a, ibase, iounit, ihi, ilo, jhi, jlo, lpmoda,
     &    maxcol, maxrow, ncol, nrow, output, title )
 
      else if ( iform .eq. 2 ) then
 
        call decprn ( iatop, iabot, ibase, iounit, ihi, ilo, jhi,
     &    jlo, lpmoda, maxcol, maxrow, ncol, nrow, output, title )
 
      end if
 
      return
      end
      subroutine types ( a, iabot, iatop, ibase, ierror, iform, imat,
     &  iounit, islbot, isltop, lpmoda, maxcol, maxrow, nart, ncol,
     &  nrow, nslak, nvar, output, sol )

c*********************************************************************72
c
cc TYPES prints out the linear programming solution.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    23 January 2009
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, real A(MAXROW,MAXCOL).  A is the current matrix.
c
c    Input, integer IATOP(MAXROW,MAXCOL), IABOT(MAXROW,MAXCOL).
c    IATOP and IABOT represent the current rational or decimal matrix.
c
c    Input, integer IBASE(MAXROW), records the basic variables.
c
c    Output, integer IERROR, error flag.
c    0, no error occurred.
c    1, an error occurred.
c
c    Input, integer IFORM, specifies the arithmetic being used.
c    0=rational, 1=real, 2=decimal.
c
c    Input, integer IMAT.
c    0, no matrix has been defined by the user.
c    1, a matrix has been defined by the user.
c
c    Input, integer IOUNIT(4).
c    IOUNIT(1) is the FORTRAN input unit.
c    IOUNIT(2) is the standard output unit, while IOUNIT(3) and
c    IOUNIT(4), if nonzero, are auxilliary output units.
c
c    Output, integer ISLBOT, ISLTOP.
c    Represents the linear programming solution, if fractional
c    or decimal arithmetic is used.
c
c    Input, integer LPMODA.
c    0, the program is in linear algebra mode.
c    1, the program is in linear programming mode.
c
c    Input, integer MAXCOL, the maximum number of columns allowed
c    in the matrices used by MATMAN.
c
c    Input, integer MAXROW, the maximum number of rows allowed
c    in the matrices used by MATMAN.
c
c    Input, integer NART, the number of artificial variables.
c
c    Input, integer NCOL, the number of columns in the matrix.
c
c    Input, integer NROW, the number of rows in the matrix.
c
c    Input, integer NSLAK, the number of slack variables.
c
c    Input, integer NVAR, the number of basic variables.
c
c    Workspace, character*100 OUTPUT.
c
c    Output, real SOL, represents the linear programming solution,
c    if real arithmetic is used.
c
      implicit none

      integer maxcol
      integer maxrow

      real a(maxrow,maxcol)
      integer iabot(maxrow,maxcol)
      integer iatop(maxrow,maxcol)
      integer ibase(maxrow)
      integer ierror
      integer iform
      integer ihi
      integer ilo
      integer imat
      integer iounit(4)
      integer islbot(maxcol)
      integer isltop(maxcol)
      integer jhi
      integer jlo
      integer lpmoda
      integer nart
      integer ncol
      integer nrow
      integer nslak
      integer nvar
      character*100 output
      real sol(maxcol)
      character*80 title
c
c  Make sure there is something to print.
c
      if ( imat .ne. 1 ) then
        ierror = 1
        output = 'Error!  You haven''t set up a tableau yet!'
        call chrwrt ( iounit, output )
        return
      end if
 
      if ( lpmoda .ne. 1 ) then
        ierror = 1
        output = 'Error!  There is no solution to print!'
        call chrwrt ( iounit, output )
        output = 'because we are not doing linear programming!'
        call chrwrt ( iounit, output )
        return
      end if
 
      call lpsol ( a, iatop, iabot, ibase, iform, isltop, islbot, 
     &  maxcol, maxrow, ncol, nrow, sol )
 
      title = 'The linear programming solution'
      ilo = 1
      ihi = 1
      jhi = nvar + nslak + nart
      jlo = 1
 
      if ( iform .eq. 0 ) then

        call ratprn ( isltop, islbot, ibase, iounit, ihi, ilo, jhi, jlo,
     &    lpmoda, jhi, 1, ncol, nrow, output, title )

      else if ( iform .eq. 1 ) then

        call relprn ( sol, ibase, iounit, ihi, ilo, jhi, jlo, lpmoda,
     &    jhi, 1, ncol, nrow, output, title )

      else if ( iform .eq. 2 ) then

        call decprn ( isltop, islbot, ibase, iounit, ihi, ilo, jhi, jlo,
     &    lpmoda, jhi, 1, ncol, nrow, output, title )

      end if
 
      return
      end
