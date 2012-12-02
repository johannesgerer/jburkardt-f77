#!/bin/bash
#
mkdir temp
cd temp
rm *
f77split ../brownian_motion_simulation.f
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
ar qc libbrownian_motion_simulation.a *.o
rm *.o
#
mv libbrownian_motion_simulation.a ~/libf77/$ARCH
cd ..
rmdir temp
#
echo "Library installed as ~/libf77/$ARCH/libbrownian_motion_simulation.a"
