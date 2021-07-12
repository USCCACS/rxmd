#!/bin/sh

for f in dimer-*.xyz;do
  file=${f%xyz}bin
  ./geninit -i $f -n
  ./geninit -i norm.xyz 
  mv rxff.bin -v ${f%xyz}bin
done
