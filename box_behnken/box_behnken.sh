#!/bin/bash
#
mkdir temp
cd temp
~/binc/$ARCH/f77split ../box_behnken.f
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
ar qc libbox_behnken.a *.o
rm *.o
#
mv libbox_behnken.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/libbox_behnken.a."
