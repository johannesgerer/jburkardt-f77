#!/bin/bash
#
mkdir temp
cd temp
rm *
~/binc/$ARCH/f77split ../linpack_z.f
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
ar qc liblinpack_z.a *.o
rm *.o
#
mv liblinpack_z.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/liblinpack_z.a."
