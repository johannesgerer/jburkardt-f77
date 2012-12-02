#!/bin/bash
#
mkdir temp
cd temp
rm *
~/binc/$ARCH/f77split ../colored_noise.f
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
ar qc libcolored_noise.a *.o
rm *.o
#
mv libcolored_noise.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/libcolored_noise.a."
