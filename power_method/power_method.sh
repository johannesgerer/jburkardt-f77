#!/bin/bash
#
mkdir temp
cd temp
~/binc/$ARCH/f77split ../power_method.f
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
ar qc libpower_method.a *.o
rm *.o
#
mv libpower_method.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/libpower_method.a."
