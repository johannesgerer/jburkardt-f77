#!/bin/bash
#
#  Compile the program with GFORTRAN.
#
gfortran -fopenmp helmholtz.f
#
#  Compile the program with IFORT.
#
#ifort -openmp -parallel -fpp helmholtz.f
#
mv a.out helmholtz
#
#  Run with 1, 2, and 4 threads.
#
echo "Run with 1 thread."
export OMP_NUM_THREADS=1
./helmholtz > helmholtz_local_output.txt
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./helmholtz >> helmholtz_local_output.txt
#
echo "Run with 4 threads."
export OMP_NUM_THREADS=4
./helmholtz >> helmholtz_local_output.txt
#
#  Discard the executable file.
#
rm helmholtz
#
echo "Program output written to helmholtz_local_output.txt"
