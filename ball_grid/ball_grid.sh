#!/bin/bash
#
mkdir temp
cd temp
rm *
~/binc/$ARCH/f77split ../ball_grid.f
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
ar qc libball_grid.a *.o
rm *.o
#
mv libball_grid.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/libball_grid.a"
