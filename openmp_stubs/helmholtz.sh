#!/bin/bash
#
gfortran ../openmp/helmholtz.f -L$HOME/libf77/$ARCH -lopenmp_stubs
mv a.out helmholtz
#
#  Run the program.
#
./helmholtz > helmholtz_output.txt
rm helmholtz
#
echo "Normal end of execution."
