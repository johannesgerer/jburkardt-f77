#!/bin/bash
#
~/binf77/$ARCH/toms431 < toms431_prb_input.txt > toms431_prb_output.txt
if [ $? -ne 0 ]; then
  echo "Errors running toms431"
  exit
fi
#
echo "Test results written to toms431_prb_output.txt."
