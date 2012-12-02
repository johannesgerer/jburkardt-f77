#!/bin/bash
#
mkdir temp
cd temp
rm *
f77split ../lagrange_interp_2d.f
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
ar qc liblagrange_interp_2d.a *.o
rm *.o
#
mv liblagrange_interp_2d.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/liblagrange_interp_2d.a"
