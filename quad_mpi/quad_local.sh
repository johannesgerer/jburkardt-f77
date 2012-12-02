#!/bin/bash
#
mpif77 quad_mpi.f
#
if [ $? -ne 0 ]; then
  echo "Errors compiling quad_mpi.f"
  exit
fi
#
#  Rename the executable.
#
mv a.out quad
#
#  Ask MPI to use 8 processes to run your program.
#
mpirun -np 8 ./quad > quad_local_output.txt
#
#  Clean up.
#
rm quad
#
echo "Program output written to quad_local_output.txt"

