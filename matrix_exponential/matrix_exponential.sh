#!/bin/bash
#
mkdir temp
cd temp
rm *
f77split ../matrix_exponential.f
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
ar qc libmatrix_exponential.a *.o
rm *.o
#
mv libmatrix_exponential.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/libmatrix_exponential.a"
