#!/bin/bash
#
mkdir temp
cd temp
~/binc/$ARCH/f77split ../../blas1_c/blas1_c.f
~/binc/$ARCH/f77split ../../blas1_d/blas1_d.f
~/binc/$ARCH/f77split ../../blas1_s/blas1_s.f
~/binc/$ARCH/f77split ../../blas1_z/blas1_z.f
~/binc/$ARCH/f77split ../../blas2/blas2.f
~/binc/$ARCH/f77split ../../blas3/blas3.f
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
#
#  Create the library.
#
ar qc libblas.a *.o
rm *.o
#
mv libblas.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/libblas.a."
