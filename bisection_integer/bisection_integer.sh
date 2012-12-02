#!/bin/bash
#
mkdir temp
cd temp
rm *
f77split ../bisection_integer.f
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
ar qc libbisection_integer.a *.o
rm *.o
#
mv libbisection_integer.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/libbisection_integer.a"
