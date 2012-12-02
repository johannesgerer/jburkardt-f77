#!/bin/bash
#
mkdir temp
cd temp
rm *
f77split ../test_int_hermite.f
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
ar qc libtest_int_hermite.a *.o
rm *.o
#
mv libtest_int_hermite.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/libtest_int_hermite.a"
