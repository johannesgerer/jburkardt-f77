#!/bin/bash
#
mkdir temp
cd temp
rm *
f77split ../sine_transform.f
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
ar qc libsine_transform.a *.o
rm *.o
#
mv libsine_transform.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/libsine_transform.a"
