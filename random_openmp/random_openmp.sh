#!/bin/bash
#
gfortran -fopenmp random_openmp.f
mv a.out random_openmp
#
./random_openmp > random_openmp_output.txt
rm random_openmp
#
echo "Program output written to random_openmp_output.txt."
