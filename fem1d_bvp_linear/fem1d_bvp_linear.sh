#!/bin/bash
#
mkdir temp
cd temp
rm *
f77split ../fem1d_bvp_linear.f
#
for FILE in `ls -1 *.f`
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
ar qc libfem1d_bvp_linear.a *.o
rm *.o
#
mv libfem1d_bvp_linear.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/libfem1d_bvp_linear.a"
