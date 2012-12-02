#!/bin/bash
#
mkdir temp
cd temp
~/binc/$ARCH/f77split ../steam_nbs.f
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
ar qc libsteam_nbs.a *.o
rm *.o
#
mv libsteam_nbs.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/libsteam_nbs.a."
