#!/bin/bash
#
gfortran -c -g sphere_lebedev_rule_prb.f >& compiler.txt
if [ $? -ne 0 ]; then
  echo "Errors compiling sphere_lebedev_rule_prb.f"
  exit
fi
rm compiler.txt
#
gfortran sphere_lebedev_rule_prb.o -L$HOME/libf77/$ARCH -lsphere_lebedev_rule
if [ $? -ne 0 ]; then
  echo "Errors linking and loading sphere_lebedev_rule_prb.o"
  exit
fi
rm sphere_lebedev_rule_prb.o
#
mv a.out sphere_lebedev_rule_prb
./sphere_lebedev_rule_prb > sphere_lebedev_rule_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running sphere_lebedev_rule_prb"
  exit
fi
rm sphere_lebedev_rule_prb
#
echo "Test results written to sphere_lebedev_rule_prb_output.txt."
