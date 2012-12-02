#!/bin/bash
#
mkdir temp
cd temp
rm *
f77split ../doomsday.f
#
for FILE in `ls -1 *.f`
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
ar qc libdoomsday.a *.o
rm *.o
#
mv libdoomsday.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/libdoomsday.a"
