#!/bin/bash
#
mkdir temp
cd temp
rm *
f77split ../ice_io.f
#
for FILE in `ls -1 *.f`;
do
  g95 -c -g -I /usr/local/include $FILE >& compiler.txt
  if [ $? -ne 0 ]; then
    echo "Errors compiling " $FILE
    exit
  fi
  rm compiler.txt
done
rm *.f
#
ar qc libice_io.a *.o
rm *.o
#
mv libice_io.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/libice_io.a"
