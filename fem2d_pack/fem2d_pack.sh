#!/bin/bash
#
mkdir temp
cd temp
~/binc/$ARCH/f77split ../fem2d_pack.f
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
ar qc libfem2d_pack.a *.o
rm *.o
#
mv libfem2d_pack.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/libfem2d_pack.a."
