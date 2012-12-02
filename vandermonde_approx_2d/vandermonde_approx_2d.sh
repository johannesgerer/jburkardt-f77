#!/bin/bash
#
mkdir temp
cd temp
rm *
f77split ../vandermonde_approx_2d.f
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
ar qc libvandermonde_approx_2d.a *.o
rm *.o
#
mv libvandermonde_approx_2d.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/libvandermonde_approx_2d.a"
