#!/bin/bash

echo CPU
echo ===========================
for i in {500..20000..500}
do
  export j=`echo $i/20 | bc`
  ./cpu.out $i $j
  ./cpu.out $i $j
  ./cpu.out $i $j
  echo
  
  export j=`echo $i/10 | bc`
  ./cpu.out $i $j
  ./cpu.out $i $j
  ./cpu.out $i $j
  echo

  export j=`echo $i/5 | bc`
  ./cpu.out $i $j
  ./cpu.out $i $j
  ./cpu.out $i $j
  echo
done
echo ===========================

echo Phi
echo ===========================
for i in {500..20000..500}
do
  export j=`echo $i/20 | bc`
  ./phi.out $i $j
  ./phi.out $i $j
  ./phi.out $i $j
  echo
  
  export j=`echo $i/10 | bc`
  ./phi.out $i $j
  ./phi.out $i $j
  ./phi.out $i $j
  echo

  export j=`echo $i/5 | bc`
  ./phi.out $i $j
  ./phi.out $i $j
  ./phi.out $i $j
  echo
done
echo ===========================

