#!/bin/bash
#
mkdir temp
cd temp
~/binc/$ARCH/f77split ../toms352.f
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
ar qc libtoms352.a *.o
rm *.o
#
mv libtoms352.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/libtoms352.a."
