#!/bin/bash
#
mkdir temp
cd temp
~/binc/$ARCH/f77split ../sandia_cvt.f
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
ar qc libsandia_cvt.a *.o
rm *.o
#
mv libsandia_cvt.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/libsandia_cvt.a."
