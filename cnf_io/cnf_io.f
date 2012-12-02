      subroutine ch_cap ( ch )

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
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character CH, the character to capitalize.
c
      implicit none

      character ch
      integer itemp

      itemp = ichar ( ch )
 
      if ( 97 .le. itemp .and. itemp .le. 122 ) then
        ch = char ( itemp - 32 )
      end if
 
      return
      end
      function ch_eqi ( c1, c2 )

c*********************************************************************72
c
cc CH_EQI is a case insensitive comparison of two characters for equality.  
c
c  Example:
c
c    CH_EQI ( 'A', 'a' ) is TRUE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character C1, C2, the characters to compare.
c
c    Output, logical CH_EQI, the result of the comparison.
c
      implicit none
   
      character c1
      character c1_cap
      character c2
      character c2_cap
      logical ch_eqi

      c1_cap = c1
      c2_cap = c2

      call ch_cap ( c1_cap )
      call ch_cap ( c2_cap )

      if ( c1_cap == c2_cap ) then
        ch_eqi = .true.
      else
        ch_eqi = .false.
      end if

      return
      end
      function ch_is_space ( ch )

c*********************************************************************72
c
cc CH_IS_SPACE is TRUE if a character is a whitespace character.
c
c  Discussion:
c
c    A whitespace character is a space, a form feed, a newline,
c    a carriage return, a tab, or a vertical tab.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character CH, a character to check.
c
c    Output, logical CH_IS_SPACE is TRUE if the character is a whitespace
c    character.
c
      implicit none

      character ch
      logical ch_is_space

      if ( ch .eq. ' ' ) then
        ch_is_space = .true.
      else if ( ch .eq. char ( 12 ) ) then
        ch_is_space = .true.
      else if ( ch .eq. char ( 10 ) ) then
        ch_is_space = .true.
      else if ( ch .eq. char ( 13 ) ) then
        ch_is_space = .true.
      else if ( ch .eq. char ( 9 ) ) then
        ch_is_space = .true.
      else if ( ch .eq. char ( 11 ) ) then
        ch_is_space = .true.
      else
        ch_is_space = .false.
      end if

      return
      end
      subroutine cnf_data_read ( cnf_file_name, v_num, c_num, l_num, 
     &  l_c_num, l_val )

c*********************************************************************72
c
cc CNF_DATA_READ reads the data of a CNF file.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) CNF_FILE_NAME, the name of the CNF file.
c
c    Input, integer V_NUM, the number of variables.
c
c    Input, integer C_NUM, the number of clauses.
c
c    Input, integer L_NUM, the number of signed literals.
c
c    Output, integer L_C_NUM(C_NUM), the number of signed
c    literals occuring in each clause.
c
c    Output, integer L_VAL(L_NUM), a list of all the signed 
c    literals in all the clauses, ordered by clause.
c
      implicit none

      integer  c_num
      integer l_num
      integer v_num

      integer c_num2
      logical ch_eqi
      logical ch_is_space
      character * ( * ) cnf_file_name
      integer cnf_file_status
      integer cnf_file_unit
      integer ierror
      integer l_c_num(c_num)
      integer l_c_num2
      integer l_num2
      integer l_val(l_num)
      integer l_val2
      integer length
      character * ( 255 ) line
      logical s_eqi
      character * ( 20 ) word

      call get_unit ( cnf_file_unit )

      open ( unit = cnf_file_unit, file = cnf_file_name, 
     &  status = 'old', iostat = cnf_file_status )

      if ( cnf_file_status .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CNF_DATA_READ - Fatal error!'
        write ( *, '(a)' ) '  Could not open filec'
        stop
      end if
c
c  Read lines until you find one that is not blank and does not begin
c  with a "c".  This should be the header line.
c
      line = ' '

10    continue

        read ( cnf_file_unit, '(a)', iostat = cnf_file_status ) line

        if ( cnf_file_status .ne. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'CNF_DATA_READ - Fatal error!'
          write ( *, '(a)' ) '  Error while reading the file,'
          write ( *, '(a)' ) '  searching for TITLE line.'
          stop
        end if

        if ( line(1:1) .eq. 'c' .or. line(1:1) .eq. 'C' ) then
          go to 10
        end if

        if ( 0 .lt. len_trim ( line ) ) then
          go to 20
        end if

      go to 10

20    continue
c
c  We expect to be reading the line "p cnf V_NUM C_NUM"
c
      if ( .not. ch_eqi ( line(1:1), 'p' ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CNF_DATA_READ - Fatal error!'
        write ( *, '(a)' ) '  First non-comment non-blank line does not'
        write ( *, '(a)' ) '  start with "p " marker.'
        stop
      end if

      if ( .not. ch_is_space ( line(2:2) ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CNF_DATA_READ - Fatal error!'
        write ( *, '(a)' ) '  Character after "p" must be whitespace.'
        stop
      end if
c
c  Remove the first two characters and shift left.
c
      line(1:2) = ' '
      line = adjustl ( line )
c
c  Expect the string 'CNF'
c
      if ( .not. s_eqi ( line(1:3), 'cnf' ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CNF_DATA_READ - Fatal error!'
        write ( *, '(a)' ) '  First non-comment non-blank line does not'
        write ( *, '(a)' ) '  start with "p cnf" marker.'
        stop
      end if

      if ( .not. ch_is_space ( line(4:4) ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CNF_DATA_READ - Fatal error!'
        write ( *, '(a)' ) 
     &    '  Character after "p cnf" must be whitespace.'
        stop
      end if
c
c  Remove the first four characters and shift left.
c
      line(1:4) = ' '
      line = adjustl ( line )
c
c  Extract the next word, which is the number of variables.
c
      call s_word_extract_first ( line, word )

      if ( len_trim ( word ) .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CNF_DATA_READ - Fatal error!'
        write ( *, '(a)' ) '  Unexpected End of input.'
        stop
      end if

      call s_to_i4 ( word, v_num, ierror, length )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CNF_DATA_READ - Fatal error!'
        write ( *, '(a)' ) '  Unexpected End of input.'
        stop
      end if
c
c  Extract the next word, which is the number of clauses.
c
      call s_word_extract_first ( line, word )

      if ( len_trim ( word ) .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CNF_DATA_READ - Fatal error!'
        write ( *, '(a)' ) '  Unexpected End of input.'
        stop
      end if

      call s_to_i4 ( word, c_num, ierror, length )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CNF_DATA_READ - Fatal error!'
        write ( *, '(a)' ) '  Unexpected End of input.'
        stop
      end if
c
c  Read remaining lines, counting the literals while ignoring occurrences of '0'.
c
      l_num2 = 0
      c_num2 = 0
      l_c_num2 = 0
      line = ' '

30    continue

        read ( cnf_file_unit, '(a)', iostat = cnf_file_status ) line

        if ( cnf_file_status .ne. 0 ) then
          go to 60
        end if

        if ( line(1:1) .eq. 'c' ) then
          go to 30
        end if

        if ( len_trim ( line ) .lt. 0 ) then
          go to 60
        end if

40      continue

          call s_word_extract_first ( line, word )

          if ( len_trim ( word ) .le. 0 ) then
            go to 50
          end if

          call s_to_i4 ( word, l_val2, ierror, length )

          if ( ierror .ne. 0 ) then
            go to 50
          end if

          if ( l_val2 .ne. 0 ) then
            l_num2 = l_num2 + 1
            l_val(l_num2) = l_val2
            l_c_num2 = l_c_num2 + 1
          else
            c_num2 = c_num2 + 1
            l_c_num(c_num2) = l_c_num2
            l_c_num2 = 0
          end if

        go to 40

50      continue

      go to 30

60    continue
c
c  At the end:
c
c    C_NUM2 should equal C_NUM, 
c    L_NUM2 should equal L_NUM.
c
c  Close file and return.
c
      close ( unit = cnf_file_unit )

      return
      end
      subroutine cnf_data_write ( c_num, l_num, l_c_num, l_val, 
     &  output_unit )

c*********************************************************************72
c
cc CNF_DATA_WRITE writes data to a CNF file.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer C_NUM, the number of clauses.
c
c    Input, integer L_NUM, the total number of signed literals.
c
c    Input, integer L_C_NUM(C_NUM), the number of signed
c    literals occuring in each clause.
c
c    Input, integer L_VAL(L_NUM), a list of all the signed 
c    literals in all the clauses, ordered by clause.
c
c    Input, integer OUTPUT_UNIT, the output unit.
c
      implicit none

      integer c_num
      integer l_num

      integer c
      integer i1
      integer i2
      integer l
      integer l_c
      integer l_c_num(c_num)
      integer l_val(l_num)
      integer output_unit
      character * ( 80 ) string

      l = 0
      string = ' '

      do c = 1, c_num

        i1 = 1
        i2 = 10

        do l_c = 1, l_c_num(c)

          l = l + 1

          write ( string(i1:i2), '(1x,i7)' ) l_val(l)

          i1 = i1 + 10
          i2 = i2 + 10

          if ( mod ( l_c, 10 ) .eq. 0 ) then
            call s_blanks_delete ( string )
            write ( output_unit, '(a)' ) string(1:len_trim(string))
            string = ' '
          end if

        end do

        string(i2+1:i2+2) = ' 0'
        call s_blanks_delete ( string )
        write ( output_unit, '(a)' ) string(1:len_trim(string))
        string = ' '

      end do

      return
      end
      function cnf_evaluate ( v_num, c_num, l_num, l_c_num, l_val, 
     &  v_val )

c*********************************************************************72
c
cc CNF_EVALUATE evaluates a formula in CNF form.
c
c  Discussion:
c
c    The formula is in conjunctive normal form.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer V_NUM, the number of variables.
c
c    Input, integer C_NUM, the number of clauses.
c
c    Input, integer L_NUM, the total number of signed literals.
c
c    Input, integer L_C_NUM(C_NUM), the number of signed
c    literals occuring in each clause.
c
c    Input, integer L_VAL(L_NUM), a list of all the signed 
c    literals in all the clauses, ordered by clause.
c
c    Input, logical V_VAL(V_NUM), the values assigned to the variables.
c
c    Output, logical CNF_EVALUATE, the value of the CNF formula for the
c    given variable values.
c
      implicit none

      integer c_num
      logical cnf_evaluate
      integer l_num
      integer v_num

      integer c
      logical c_val
      logical f_val
      integer l
      integer l_c
      integer l_c_num(c_num)
      integer l_val(l_num)
      logical s_val
      integer v_index
      logical v_val(v_num)
      logical value

      f_val = .true.

      l = 0

      do c = 1, c_num
c
c  The clause is false unless some signed literal is true.
c
        c_val = .false.
        do l_c = 1, l_c_num(c)
          l = l + 1
          s_val = ( 0 .lt. l_val(l) )
          v_index = abs ( l_val(l) )
c
c  The signed literal is true if the sign "equals" the value.
c  Note that we CANNOT exit the loop because we need to run out the 
c  L indexc
c
          if ( v_val(v_index) .eqv. s_val ) then
            c_val = .true.
          end if
        end do
c
c  The formula is false if any clause is false.
c
        if ( .not. c_val ) then
          f_val = .false.
          go to 10
        end if

      end do

10    continue

      cnf_evaluate = f_val

      return
      end
      subroutine cnf_header_read ( cnf_file_name, v_num, c_num, l_num )

c*********************************************************************72
c
cc CNF_HEADER_READ reads the header of a CNF file.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) CNF_FILE_NAME, the name of the CNF file.
c
c    Output, integer V_NUM, the number of variables.
c
c    Output, integer C_NUM, the number of clauses.
c
c    Output, integer L_NUM, the number of signed literals.
c
      implicit none

      integer c_num
      logical ch_eqi
      logical ch_is_space
      character * ( * ) cnf_file_name
      integer cnf_file_status
      integer cnf_file_unit
      integer ierror
      integer l_num
      integer l_val
      integer length
      character * ( 255 ) line
      logical s_eqi
      integer v_num
      character * ( 20 ) word

      call get_unit ( cnf_file_unit )

      open ( unit = cnf_file_unit, file = cnf_file_name, 
     &  status = 'old', iostat = cnf_file_status )

      if ( cnf_file_status .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CNF_HEADER_READ - Fatal error!'
        write ( *, '(a)' ) '  Could not open filec'
        stop
      end if
c
c  Read lines until you find one that is not blank and does not begin
c  with a "c".  This should be the header line.
c
      line = ' '

10    continue

        read ( cnf_file_unit, '(a)', iostat = cnf_file_status ) line

        if ( cnf_file_status .ne. 0 ) then
          write ( *, '(a)' ) ' '
          write ( *, '(a)' ) 'CNF_HEADER_READ - Fatal error!'
          write ( *, '(a)' ) '  Error while reading the file,'
          write ( *, '(a)' ) '  searching for TITLE line.'
          stop
        end if

        if ( line(1:1) .eq. 'c' .or. line(1:1) .eq. 'C' ) then
          go to 10
        end if

        if ( 0 .lt. len_trim ( line ) ) then
          go to 20
        end if

      go to 10

20    continue
c
c  We expect to be reading the line "p cnf V_NUM C_NUM"
c
      if ( .not. ch_eqi ( line(1:1), 'p' ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CNF_HEADER_READ - Fatal error!'
        write ( *, '(a)' ) '  First non-comment non-blank line does '
        write ( *, '(a)' ) '  not start with "p " marker.'
        stop
      end if

      if ( .not. ch_is_space ( line(2:2) ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CNF_HEADER_READ - Fatal error!'
        write ( *, '(a)' ) '  Character after "p" must be whitespace.'
        stop
      end if
c
c  Remove the first two characters and shift left.
c
      line(1:2) = ' '
      line = adjustl ( line )
c
c  Expect the string 'CNF'
c
      if ( .not. s_eqi ( line(1:3), 'cnf' ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CNF_HEADER_READ - Fatal error!'
        write ( *, '(a)' ) '  First non-comment non-blank line does'
        write ( *, '(a)' ) '  not start with "p cnf" marker.'
        stop
      end if

      if ( .not. ch_is_space ( line(4:4) ) ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CNF_HEADER_READ - Fatal error!'
        write ( *, '(a)' ) 
     &    '  Character after "p cnf" must be whitespace.'
        stop
      end if
c
c  Remove the first four characters and shift left.
c
      line(1:4) = ' '
      line = adjustl ( line )
c
c  Extract the next word, which is the number of variables.
c
      call s_word_extract_first ( line, word )

      if ( len_trim ( word ) .le. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CNF_HEADER_READ - Fatal error!'
        write ( *, '(a)' ) '  Unexpected End of input.'
        stop
      end if

      call s_to_i4 ( word, v_num, ierror, length )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CNF_HEADER_READ - Fatal error!'
        write ( *, '(a)' ) '  Unexpected End of input.'
        stop
      end if
c
c  Extract the next word, which is the number of clauses.
c
      call s_word_extract_first ( line, word )

      if ( len_trim ( word ) .eq. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CNF_HEADER_READ - Fatal error!'
        write ( *, '(a)' ) '  Unexpected End of input.'
        stop
      end if

      call s_to_i4 ( word, c_num, ierror, length )

      if ( ierror .ne. 0 ) then
        write ( *, '(a)' ) ' '
        write ( *, '(a)' ) 'CNF_HEADER_READ - Fatal error!'
        write ( *, '(a)' ) '  Unexpected End of input.'
        stop
      end if
c
c  Read remaining lines, counting the literals while ignoring occurrences of '0'.
c
      l_num = 0
      line = ' '

30    continue

        read ( cnf_file_unit, '(a)', iostat = cnf_file_status ) line

        if ( cnf_file_status .ne. 0 ) then
          go to 60
        end if

        if ( line(1:1) .eq. 'c' ) then
          go to 30
        end if

        if ( len_trim ( line ) .lt. 0 ) then
          go to 60
        end if

40      continue

          call s_word_extract_first ( line, word )

          if ( len_trim ( word ) .le. 0 ) then
            go to 50
          end if

          call s_to_i4 ( word, l_val, ierror, length )

          if ( ierror .ne. 0 ) then
            go to 50
          end if

          if ( l_val .ne. 0 ) then
            l_num = l_num + 1
          end if

        go to 40

50      continue

      go to 30

60    continue
c
c  Close file and return.
c
      close ( unit = cnf_file_unit )

      return
      end
      subroutine cnf_header_write ( v_num, c_num, output_name, 
     &  output_unit )

c*********************************************************************72
c
cc CNF_HEADER_WRITE writes the header for a CNF file.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer V_NUM, the number of variables.
c
c    Input, integer C_NUM, the number of clauses.
c
c    Input, character * ( * ) OUTPUT_NAME, the name of the output file.
c
c    Input, integer OUTPUT_UNIT, the output unit.
c
      implicit none

      integer c_num
      character * ( * ) output_name
      integer output_unit
      character * ( 80 ) string
      integer v_num

      write ( output_unit, '(a)' ) 'c ' // trim ( output_name )
      write ( output_unit, '(a)' ) 'c'
      write ( string, '(a,1x,i7,1x,i7)' ) 'p cnf', v_num, c_num
      call s_blanks_delete ( string )
      write ( output_unit, '(a)' ) string(1:len_trim(string))

      return
      end
      subroutine cnf_print ( v_num, c_num, l_num, l_c_num, l_val )

c*********************************************************************72
c
cc CNF_PRINT prints CNF information.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer V_NUM, the number of variables.
c
c    Input, integer C_NUM, the number of clauses.
c
c    Input, integer L_NUM, the total number of signed literals.
c
c    Input, integer L_C_NUM(C_NUM), the number of signed
c    literals occuring in each clause.
c
c    Input, integer L_VAL(L_NUM), a list of all the signed 
c    literals in all the clauses, ordered by clause.
c
      implicit none

      integer c_num
      integer l_num

      integer c
      integer l
      integer l_c
      integer l_c_num(c_num)
      integer l_val(l_num)
      integer v_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CNF data printout:'
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' ) 
     &  '  The number of variables       V_NUM  = ', v_num
      write ( *, '(a,i8)' ) 
     &  '  The number of clauses         C_NUM  = ', c_num
      write ( *, '(a,i8)' ) 
     &  '  The number of signed literals L_NUM  = ', l_num
      l = 0
      do c = 1, c_num
        write ( *, '(a)' ) ' '
        write ( *, '(a,i8,a,i8,a)' ) 
     &    '  Clause ', c, ' includes ', l_c_num(c), ' signed literals'
        do l_c = 1, l_c_num(c)
          l = l + 1
          write ( *, '(i4)' ) l_val(l)
        end do
      end do

      return
      end
      subroutine cnf_write ( v_num, c_num, l_num, l_c_num, l_val, 
     &  output_name )

c*********************************************************************72
c
cc CNF_WRITE writes the header and data of a CNF file.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer V_NUM, the number of variables.
c
c    Input, integer C_NUM, the number of clauses.
c
c    Input, integer L_NUM, the total number of signed literals.
c
c    Input, integer L_C_NUM(C_NUM), the number of signed
c    literals occuring in each clause.
c
c    Input, integer L_VAL(L_NUM), a list of all the signed 
c    literals in all the clauses, ordered by clause.
c
c    Input, character * ( * ) OUTPUT_NAME, the name of the output file.
c
      implicit none

      integer c_num
      integer l_num

      integer l_c_num(c_num)
      integer l_val(l_num)
      character * ( * ) output_name
      integer output_unit
      integer v_num

      call get_unit ( output_unit )

      open ( unit = output_unit, file = output_name, 
     &  status = 'replace' )

      call cnf_header_write ( v_num, c_num, output_name, output_unit )

      call cnf_data_write ( c_num, l_num, l_c_num, l_val, output_unit )

      close ( unit = output_unit )

      return
      end
      subroutine get_unit ( unit )

c*********************************************************************72
c
cc GET_UNIT returns a free FORTRAN unit number.
c
c  Discussion:
c
c    A "free" FORTRAN unit number is an integer between 1 and 99 which
c    is not currently associated with an I/O device.  A free FORTRAN unit
c    number is needed in order to open a file with the OPEN command.
c
c    If UNIT = 0, then no free FORTRAN unit could be found, although
c    all 99 units were checked (except for units 5, 6 and 9, which
c    are commonly reserved for console I/O).
c
c    Otherwise, UNIT is an integer between 1 and 99, representing a
c    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
c    are special, and will never return those values.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    26 Ocober 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Output, integer UNIT, the free unit number.
c
      implicit none

      integer i
      integer unit

      unit = 0

      do i = 1, 99

        if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

          open ( unit = i, err = 10, status = 'scratch' )
          close ( unit = i )

          unit = i

          return
        end if

10      continue

      end do

      return
      end
      subroutine lvec_next ( n, lvec )

c*********************************************************************72
c
cc  Purpose:
c
c    LVEC_NEXT generates the next logical vector.
c
c  Discussion:
c
c    In the following discussion, we will let '0' stand for FALSE and
c    '1' for TRUE.
c
c    The vectors have the order
c
c      (0,0,...,0),
c      (0,0,...,1),
c      ...
c      (1,1,...,1)
c
c    and the "next" vector after (1,1,...,1) is (0,0,...,0).  That is,
c    we allow wrap around.
c
c  Example:
c
c    N = 3
c
c    Input      Output
c    -----      ------
c    0 0 0  =>  0 0 1
c    0 0 1  =>  0 1 0
c    0 1 0  =>  0 1 1
c    0 1 1  =>  1 0 0
c    1 0 0  =>  1 0 1
c    1 0 1  =>  1 1 0
c    1 1 0  =>  1 1 1
c    1 1 1  =>  0 0 0
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    31 May 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, integer N, the dimension of the vectors.
c
c    Input/output, logical LVEC(N), on output, the successor to the
c    input vector.
c
      implicit none

      integer n

      integer i
      logical lvec(n)

      do i = n, 1, -1

        if ( .not. lvec(i) ) then
          lvec(i) = .true.
          return
        end if

        lvec(i) = .false.

      end do

      return
      end
      subroutine s_adjustl ( s )

c*********************************************************************72
c
cc S_ADJUSTL flushes a string left.
c
c  Discussion:
c
c    Both blanks and tabs are treated as "white space".
c
c    This routine is similar to the FORTRAN90 ADJUSTL routine.
c
c  Example:
c
c    Input             Output
c
c    '     Hello'      'Hello     '
c    ' Hi therec  '    'Hi therec   '
c    'Fred  '          'Fred  '
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 Jun3 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character * ( * ) S.
c    On input, S is a string of characters.
c    On output, any initial blank or tab characters have been cut.
c
      implicit none

      integer i
      integer nonb
      character * ( * ) s
      integer s_length
      character tab

      tab = char ( 9 )
c
c  Check the length of the string to the last nonblank.
c  If nonpositive, return.
c
      s_length = len_trim ( s )

      if ( s_length .le. 0 ) then
        return
      end if
c
c  Find NONB, the location of the first nonblank, nontab.
c
      nonb = 0

      do i = 1, s_length

        if ( s(i:i) .ne. ' ' .and. s(i:i) .ne. tab ) then
          nonb = i
          go to 10
        end if

      end do

10    continue

      if ( nonb .eq. 0 ) then
        s = ' '
        return
      end if
c
c  Shift the string left.
c
      if ( 1 .lt. nonb ) then
        do i = 1, s_length + 1 - nonb
          s(i:i) = s(i+nonb-1:i+nonb-1)
        end do
      end if
c
c  Blank out the end of the string.
c
      s(s_length+2-nonb:s_length) = ' '
     
      return
      end
      subroutine s_blanks_delete ( s )

c*********************************************************************72
c
cc S_BLANKS_DELETE replaces consecutive blanks by one blank.
c
c  Discussion:
c
c    Thanks to Bill Richmond for pointing out a programming flaw which
c    meant that, as characters were slid to the left through multiple
c    blanks, their original images were not blanked out.  This problem
c    is easiest resolved by using a copy of the string.
c
c    The remaining characters are left justified and right padded with blanks.
c    TAB characters are converted to spaces.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character*(*) S, the string to be transformed.
c
      implicit none

      integer i
      integer j
      character newchr
      character oldchr
      character*(*) s
      character*255 s_copy
      integer s_length
      character tab

      tab = char ( 9 )

      s_length = len ( s )

      j = 0
      s_copy(1:s_length) = s(1:s_length)
      s(1:s_length) = ' '

      newchr = ' '

      do i = 1, s_length

        oldchr = newchr
        newchr = s_copy(i:i)

        if ( newchr == tab ) then
          newchr = ' '
        end if
  
        if ( oldchr .ne. ' ' .or. newchr .ne. ' ' ) then
          j = j + 1
          s(j:j) = newchr
        end if

      end do

      return
      end
      function s_eqi ( s1, s2 )

c*********************************************************************72
c
cc S_EQI is a case insensitive comparison of two strings for equality.
c
c  Example:
c
c    S_EQI ( 'Anjana', 'ANJANA' ) is TRUE.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    03 January 2007
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) S1, S2, the strings to compare.
c
c    Output, logical S_EQI, the result of the comparison.
c
      implicit none

      character c1
      character c2
      integer i
      integer lenc
      logical s_eqi
      character*(*) s1
      integer s1_length
      character*(*) s2
      integer s2_length

      s1_length = len ( s1 )
      s2_length = len ( s2 )
      lenc = min ( s1_length, s2_length )

      s_eqi = .false.

      do i = 1, lenc

        c1 = s1(i:i)
        c2 = s2(i:i)
        call ch_cap ( c1 )
        call ch_cap ( c2 )

        if ( c1 .ne. c2 ) then
          return
        end if

      end do

      do i = lenc + 1, s1_length
        if ( s1(i:i) .ne. ' ' ) then
          return
        end if
      end do

      do i = lenc + 1, s2_length
        if ( s2(i:i) .ne. ' ' ) then
          return
        end if
      end do

      s_eqi = .true.

      return
      end
      function s_len_trim ( s )

c*********************************************************************72
c
cc S_LEN_TRIM returns the length of a string to the last nonblank.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    05 March 2004
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character*(*) S, a string.
c
c    Output, integer S_LEN_TRIM, the length of the string to the last nonblank.
c
      implicit none

      integer i
      character*(*) s
      integer s_len_trim

      do i = len ( s ), 1, -1

        if ( s(i:i) .ne. ' ' ) then
          s_len_trim = i
          return
        end if

      end do

      s_len_trim = 0

      return
      end
      subroutine s_to_i4 ( s, ival, ierror, length )

c*********************************************************************72
c
cc S_TO_I4 reads an I4 from a string.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license.
c
c  Modified:
c
c    28 April 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input, character * ( * ) S, a string to be examined.
c
c    Output, integer IVAL, the integer value read from the string.
c    If the string is blank, then IVAL will be returned 0.
c
c    Output, integer IERROR, an error flag.
c    0, no error.
c    1, an error occurred.
c
c    Output, integer LENGTH, the number of characters of S
c    used to make IVAL.
c
      implicit none

      character c
      integer i
      integer ierror
      integer isgn
      integer istate
      integer ival
      integer length
      character * ( * ) s
      integer s_len_trim

      ierror = 0
      istate = 0
      isgn = 1
      ival = 0

      do i = 1, s_len_trim ( s )

        c = s(i:i)
c
c  Haven't read anything.
c
        if ( istate .eq. 0 ) then

          if ( c .eq. ' ' ) then

          else if ( c .eq. '-' ) then
            istate = 1
            isgn = -1
          else if ( c .eq. '+' ) then
            istate = 1
            isgn = + 1
          else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
            istate = 2
            ival = ichar ( c ) - ichar ( '0' )
          else
            ierror = 1
            return
          end if
c
c  Have read the sign, expecting digits.
c
        else if ( istate .eq. 1 ) then

          if ( c .eq. ' ' ) then

          else if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
            istate = 2
            ival = ichar ( c ) - ichar ( '0' )
          else
            ierror = 1
            return
          end if
c
c  Have read at least one digit, expecting more.
c
        else if ( istate .eq. 2 ) then

          if ( lle ( '0', c ) .and. lle ( c, '9' ) ) then
            ival = 10 * ival + ichar ( c ) - ichar ( '0' )
          else
            ival = isgn * ival
            length = i - 1
            return
          end if

        end if

      end do
c
c  If we read all the characters in the string, see if we're OK.
c
      if ( istate .eq. 2 ) then
        ival = isgn * ival
        length = s_len_trim ( s )
      else
        ierror = 1
        length = 0
      end if

      return
      end
      subroutine s_word_extract_first ( s, w )

c*********************************************************************72
c
cc S_WORD_EXTRACT_FIRST extracts the first word from a string.
c
c  Discussion:
c
c    A "word" is a string of characters terminated by a blank or
c    the end of the string.
c
c  Licensing:
c
c    This code is distributed under the GNU LGPL license. 
c
c  Modified:
c
c    03 June 2008
c
c  Author:
c
c    John Burkardt
c
c  Parameters:
c
c    Input/output, character * ( * ) S, the string.  On output, the first
c    word has been removed, and the remaining string has been shifted left.
c
c    Output, character * ( * ) W, the leading word of the string.
c
      implicit none

      integer get1
      integer get2
      character * ( * ) s
      integer s_len_trim
      integer s_length
      character * ( * ) w

      w = ' '

      s_length = s_len_trim ( s )

      if ( s_length .lt. 1 ) then
        return
      end if
c
c  Find the first nonblank.
c
      get1 = 0

10    continue

        get1 = get1 + 1

        if ( s_length .lt. get1 ) then
          return
        end if

        if ( s(get1:get1) .ne. ' ' ) then
          go to 20
        end if

      go to 10

20    continue
c
c  Look for the last contiguous nonblank.
c
      get2 = get1

30    continue

        if ( s_length .le. get2 ) then
          go to 40
        end if

        if ( s(get2+1:get2+1) .eq. ' ' ) then
          go to 40
        end if

        get2 = get2 + 1

      go to 30

40    continue
c
c  Copy the word.
c
      w = s(get1:get2)
c
c  Shift the string.
c
      s(1:get2) = ' '
      s = adjustl ( s )

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
