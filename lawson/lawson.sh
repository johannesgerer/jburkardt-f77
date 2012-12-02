#!/bin/bash
#
mkdir temp
cd temp
rm *
~/binc/$ARCH/f77split ../lawson.f
#
foreach FILE (`ls -1 *.f`)
  gfortran -c -g $FILE >& compiler.txt
  if [ $? -ne 0 ]; then
    echo "Errors compiling " $FILE
    exit
  fi
  rm compiler.txt
end
rm *.f
#
ar qc liblawson.a *.o
rm *.o
#
mv liblawson.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/liblawson.a."
