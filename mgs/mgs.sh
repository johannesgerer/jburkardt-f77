#!/bin/bash
#
gfortran -c mgs.f
if [ $? -ne 0 ]; then
  echo "Errors compiling mgs.f"
  exit
fi
#
echo "The mgs.f file was compiled."
