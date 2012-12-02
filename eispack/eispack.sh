#!/bin/bash
#
mkdir temp
cd temp
rm *
~/binc/$ARCH/f77split ../eispack.f
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
ar qc libeispack.a *.o
rm *.o
#
mv libeispack.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/libeispack.a."
