#!/bin/bash
#
mkdir temp
cd temp
rm *
f77split ../subset_sum.f
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
ar qc libsubset_sum.a *.o
rm *.o
#
mv libsubset_sum.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/libsubset_sum.a"
