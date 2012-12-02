#!/bin/bash
#
cp mpi_stubs_f77.h ~/include/mpi_stubs_f77.h
#
mkdir temp
cd temp
~/binc/$ARCH/f77split ../mpi_stubs.f
cp ../mpi_stubs_f77.h .
#
for FILE in `ls -1 *.f`;
do
  gfortran -c -g $FILE >& compiler.txt
  if [ $? -ne 0 ]; then
    echo "Errors compiling " $FILE
    exit
  fi
  rm compiler.txt
done
rm *.f
rm mpi_stubs_f77.h
#
ar qc libmpi_stubs.a *.o
rm *.o
#
mv libmpi_stubs.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/libmpi_stubs.a."
