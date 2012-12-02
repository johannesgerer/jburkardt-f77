#!/bin/bash
#
gfortran ../hello_openmp/hello_openmp.f -L$HOME/libf77/$ARCH -lopenmp_stubs
mv a.out hello
#
#  Run the program.
#
./hello > hello_output.txt
rm hello
#
echo "Normal end of execution."
