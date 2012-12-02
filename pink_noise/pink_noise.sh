#!/bin/bash
#
mkdir temp
cd temp
rm *
~/binc/$ARCH/f77split ../pink_noise.f
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
ar qc libpink_noise.a *.o
rm *.o
#
mv libpink_noise.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/libpink_noise.a."
