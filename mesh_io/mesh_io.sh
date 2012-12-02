#!/bin/bash
#
mkdir temp
cd temp
rm *
~/binc/$ARCH/f77split ../mesh_io.f
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
ar qc libmesh_io.a *.o
rm *.o
#
mv libmesh_io.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/libmesh_io.a"
