#!/bin/bash
#
mkdir temp
cd temp
~/binc/$ARCH/f77split ../fftpack5.1.f
#
for FILE in `ls -1 *.f`;
do
  gfortran -c -O3 $FILE >& compiler.txt
  if [ $? -ne 0 ]; then
    echo "Errors compiling " $FILE
    exit
  fi
  rm compiler.txt
done
rm *.f
#
ar cr libfftpack5.1.a *.o
rm *.o
#
mv libfftpack5.1.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/libfftpack5.1.a."
