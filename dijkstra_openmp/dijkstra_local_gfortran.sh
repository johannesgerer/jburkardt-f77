#!/bin/bash
#
#  Compile the program with GFORTRAN.
#
gfortran -fopenmp dijkstra_openmp.f
#
mv a.out dijkstra
#
#  Run with 1, 2, and 4 threads.
#
echo "Run with 1 thread."
export OMP_NUM_THREADS=1
./dijkstra > dijkstra_local_gfortran_output.txt
#
echo "Run with 2 threads."
export OMP_NUM_THREADS=2
./dijkstra >> dijkstra_local_gfortran_output.txt
#
echo "Run with 4 threads."
export OMP_NUM_THREADS=4
./dijkstra >> dijkstra_local_gfortran_output.txt
#
#  Discard the executable file.
#
rm dijkstra
#
echo "Program output written to dijkstra_local_gfortran_output.txt"
