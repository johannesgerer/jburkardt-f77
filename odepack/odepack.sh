#!/bin/bash
#
mkdir temp
cd temp
rm *
~/binc/$ARCH/f77split ../odepack.f
~/binc/$ARCH/f77split ../odepack_sub1.f
~/binc/$ARCH/f77split ../odepack_sub2.f
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
ar qc libodepack.a *.o
rm *.o
#
mv libodepack.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/libodepack.a."
