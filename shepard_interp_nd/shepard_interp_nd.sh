#!/bin/bash
#
mkdir temp
cd temp
rm *
f77split ../shepard_interp_nd.f
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
ar qc libshepard_interp_nd.a *.o
rm *.o
#
mv libshepard_interp_nd.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/libshepard_interp_nd.a"