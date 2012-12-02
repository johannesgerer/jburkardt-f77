#!/bin/bash
#
mkdir temp
cd temp
rm *
~/binc/$ARCH/f77split ../nl2sol.f
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
ar qc libnl2sol.a *.o
rm *.o
#
mv libnl2sol.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/libnl2sol.a."
