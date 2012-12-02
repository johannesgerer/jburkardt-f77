      program main

c*********************************************************************72
c
cc MAIN is the main program for CNF_IO_PRB.
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
      call timestamp ( )
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CNF_IO_PRB'
      write ( *, '(a)' ) '  FORTRAN77 version'
      write ( *, '(a)' ) '  Test the CNF_IO library.'

      call test01 ( )
      call test02 ( )
      call test03 ( )
      call test04 ( )
      call test05 ( )
      call test06 ( )
      call test07 ( )
c
c  Terminate.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'CNF_IO_PRB'
      write ( *, '(a)' ) '  Normal end of execution.'
      write ( *, '(a)' ) ' '
      call timestamp ( )

      stop
      end
      subroutine test01 ( )

c*********************************************************************72
c
cc TEST01 calls CNF_WRITE to write a small CNF example to a CNF file.
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
      integer c_num
      parameter ( c_num = 2 )
      integer l_num
      parameter ( l_num = 5 )

      character * ( 80 ) cnf_file_name
      integer l_c_num(c_num)
      data l_c_num / 2, 3 /
      integer l_val(l_num)
      data l_val / 1, -3, 2, 3, -1 /
      integer v_num
      parameter ( v_num = 3 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST01'
      write ( *, '(a)' ) '  CNF_WRITE can write CNF data to a CNF file.'

      cnf_file_name = 'cnf_io_v3_c2.cnf'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Here is the data:'

      call cnf_print ( v_num, c_num, l_num, l_c_num, l_val )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  Now we call CNF_WRITE to store this information'
      write ( *, '(a)' )
     &  '  in the file "' // trim ( cnf_file_name ) // '".'

      call cnf_write ( v_num, c_num, l_num, l_c_num, l_val,
     &  cnf_file_name )

      return
      end
      subroutine test02 ( )

c*********************************************************************72
c
cc TEST02 calls CNF_HEADER_READ to read the header of a small CNF example file.
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
      integer c_num
      character * ( 80 ) cnf_file_name
      integer l_num
      integer v_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST02'
      write ( *, '(a)' )
     &  '  CNF_HEADER_READ reads the header of a CNF file.'

      cnf_file_name = 'cnf_io_v3_c2.cnf'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  Read the header of "' // trim ( cnf_file_name ) // '".'

      call cnf_header_read ( cnf_file_name, v_num, c_num, l_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' )
     &  '  The number of variables       V_NUM  = ', v_num
      write ( *, '(a,i8)' )
     &  '  The number of clauses         C_NUM  = ', c_num
      write ( *, '(a,i8)' )
     &  '  The number of signed literals L_NUM  = ', l_num

      return
      end
      subroutine test03 ( )

c*********************************************************************72
c
cc TEST03 calls CNF_DATA_READ to read the header of a small CNF example file.
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
      integer c_num_max
      parameter ( c_num_max = 100 )
      integer l_num_max
      parameter ( l_num_max = 100 )

      integer c_num
      character * ( 80 ) cnf_file_name
      integer l_c_num(c_num_max)
      integer l_num
      integer l_val(l_num_max)
      integer v_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST03'
      write ( *, '(a)' ) '  CNF_DATA_READ reads the data of a CNF file.'

      cnf_file_name = 'cnf_io_v3_c2.cnf'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  Read the header of "' // trim ( cnf_file_name ) // '".'

      call cnf_header_read ( cnf_file_name, v_num, c_num, l_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' )
     &  '  The number of variables       V_NUM  = ', v_num
      write ( *, '(a,i8)' )
     &  '  The number of clauses         C_NUM  = ', c_num
      write ( *, '(a,i8)' )
     &  '  The number of signed literals L_NUM  = ', l_num

      call cnf_data_read ( cnf_file_name, v_num, c_num, l_num,
     &  l_c_num, l_val )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Here is the data as read from the file:'

      call cnf_print ( v_num, c_num, l_num, l_c_num, l_val )

      return
      end
      subroutine test04 ( )

c*********************************************************************72
c
cc TEST04 calls CNF_WRITE to write a CNF example to a CNF file.
c
c  Discussion:
c
c    This formula is used as an example in the Quinn reference.
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
      integer c_num
      parameter ( c_num = 18 )
      integer l_num
      parameter ( l_num = 36 )

      character * ( 80 ) cnf_file_name
      integer l_c_num(c_num)
      data l_c_num /
     &  2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
     &  2, 2, 2, 2, 2, 2, 2, 2 /
      integer l_val(l_num)
      data l_val /
     &  1,    2,
     & -2,   -4,
     &  3,    4,
     & -4,   -5,
     &  5,   -6,
     &  6,   -7,
     &  6,    7,
     &  7,  -16,
     &  8,   -9,
     & -8,  -14,
     &  9,   10,
     &  9,  -10,
     &-10,  -11,
     & 10,   12,
     & 11,   12,
     & 13,   14,
     & 14,  -15,
     & 15,   16 /
      integer v_num
      parameter ( v_num = 16 )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST04'
      write ( *, '(a)' ) '  CNF_WRITE can write CNF data to a CNF file.'

      cnf_file_name = 'cnf_io_v16_c18.cnf'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  Here is the data to be written to the file:'

      call cnf_print ( v_num, c_num, l_num, l_c_num, l_val )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  Now we call CNF_WRITE to store this information'
      write ( *, '(a)' )
     &  '  in the file "' // trim ( cnf_file_name ) // '".'

      call cnf_write ( v_num, c_num, l_num, l_c_num, l_val,
     &  cnf_file_name )

      return
      end
      subroutine test05 ( )

c*********************************************************************72
c
cc TEST05 calls CNF_HEADER_READ to read the header of a CNF example file.
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
      integer c_num
      character * ( 80 ) cnf_file_name
      integer l_num
      integer v_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST05'
      write ( *, '(a)' )
     &  '  CNF_HEADER_READ reads the header of a CNF file.'

      cnf_file_name = 'cnf_io_v16_c18.cnf'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  Read the header of "' // trim ( cnf_file_name ) // '".'

      call cnf_header_read ( cnf_file_name, v_num, c_num, l_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' )
     &  '  The number of variables       V_NUM  = ', v_num
      write ( *, '(a,i8)' )
     &  '  The number of clauses         C_NUM  = ', c_num
      write ( *, '(a,i8)' )
     &  '  The number of signed literals L_NUM  = ', l_num

      return
      end
      subroutine test06 ( )

c*********************************************************************72
c
cc TEST06 calls CNF_DATA_READ to read the header of a small CNF example file.
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
      integer c_num_max
      parameter ( c_num_max = 100 )
      integer l_num_max
      parameter ( l_num_max = 100 )

      integer c_num
      character * ( 80 ) cnf_file_name
      integer l_c_num(c_num_max)
      integer l_num
      integer l_val(l_num_max)
      integer v_num

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST06'
      write ( *, '(a)' )
     &  '  CNF_DATA_READ reads the data of a CNF file.'

      cnf_file_name = 'cnf_io_v16_c18.cnf'

      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '  Read the header of "' // trim ( cnf_file_name ) // '".'

      call cnf_header_read ( cnf_file_name, v_num, c_num, l_num )

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' )
     &  '  The number of variables       V_NUM  = ', v_num
      write ( *, '(a,i8)' )
     &  '  The number of clauses         C_NUM  = ', c_num
      write ( *, '(a,i8)' )
     &  '  The number of signed literals L_NUM  = ', l_num

      call cnf_data_read ( cnf_file_name, v_num, c_num, l_num,
     &  l_c_num, l_val )

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Here is the data as read from the file:'

      call cnf_print ( v_num, c_num, l_num, l_c_num, l_val )

      return
      end
      subroutine test07 ( )

c*********************************************************************72
c
cc TEST07 calls CNF_EVALUATE to evaluate a formula.
c
c  Discussion:
c
c    This formula is used as an example in the Quinn reference.
c    Here, we seek the logical inputs that make the formula true.
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
      integer c_num
      parameter ( c_num = 18 )
      integer l_num
      parameter ( l_num = 36 )
      integer v_num
      parameter ( v_num = 16 )

      logical cnf_evaluate
      logical f_val
      integer i
      integer ihi
      integer l_c_num(c_num)
      data l_c_num /
     &  2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
     &  2, 2, 2, 2, 2, 2, 2, 2 /
      integer l_val(l_num)
      data l_val /
     &  1,    2,
     & -2,   -4,
     &  3,    4,
     & -4,   -5,
     &  5,   -6,
     &  6,   -7,
     &  6,    7,
     &  7,  -16,
     &  8,   -9,
     & -8,  -14,
     &  9,   10,
     &  9,  -10,
     &-10,  -11,
     & 10,   12,
     & 11,   12,
     & 13,   14,
     & 14,  -15,
     & 15,   16 /
      integer solution_num
      logical v_val(v_num)

      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'TEST07'
      write ( *, '(a)' )
     &  '  Seek inputs to a circuit that produce a 1 (TRUE) output.'
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) '  Here is the CNF data defining the formula:'

      call cnf_print ( v_num, c_num, l_num, l_c_num, l_val )
c
c  Initialize the logical vector.
c
      v_val(1:v_num) = .false.
c
c  Compute the number of binary vectors to check.
c
      ihi = 2**v_num

      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' )
     &  '  Number of input vectors to check is  ', ihi
      write ( *, '(a)' ) ' '
      write ( *, '(a)' )
     &  '   #       Index    ---------Input Values----------'
      write ( *, '(a)' ) ' '
c
c  Check every possible input vector.
c
      solution_num = 0

      do i = 1, ihi

        f_val = cnf_evaluate ( v_num, c_num, l_num, l_c_num, l_val,
     &    v_val )

        if ( f_val ) then
          solution_num = solution_num + 1
          write ( *, '(2x,i2,2x,i10,3x,16l1)' )
     &      solution_num, i - 1, v_val(1:v_num)
        end if

        call lvec_next ( v_num, v_val )

      end do
c
c  Report.
c
      write ( *, '(a)' ) ' '
      write ( *, '(a,i8)' )
     &  '  Number of solutions found was ', solution_num

      return
      end
