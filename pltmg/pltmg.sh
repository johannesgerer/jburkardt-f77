#!/bin/bash
#
mkdir temp
cd temp
rm *
~/binc/$ARCH/f77split ../pltmg.f
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
gcc -c ../mgxdr.c >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling mgxdr.c"
  exit
fi
rm compiler.txt
#
ar qc libpltmg.a *.o
rm *.o
#
mv libpltmg.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/libpltmg.a"
