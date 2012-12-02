#!/bin/bash
#
mkdir temp
cd temp
~/binc/$ARCH/f77split ../beta_nc.f
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
ar qc libbeta_nc.a *.o
rm *.o
#
mv libbeta_nc.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/libbeta_nc.a."
