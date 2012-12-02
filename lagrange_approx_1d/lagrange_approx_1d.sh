#!/bin/bash
#
mkdir temp
cd temp
rm *
f77split ../lagrange_approx_1d.f
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
ar qc liblagrange_approx_1d.a *.o
rm *.o
#
mv liblagrange_approx_1d.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/liblagrange_approx_1d.a"
