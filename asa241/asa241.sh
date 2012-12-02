#!/bin/bash
#
mkdir temp
cd temp
rm *
~/binc/$ARCH/f77split ../asa241.f
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
ar qc libasa241.a *.o
rm *.o
#
mv libasa241.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/libasa241.a."
