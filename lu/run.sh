#!/bin/bash

echo CPU Banded
echo ===========================
for i in {500..20000..500}
do
  export j=`echo $i/20 | bc`
  ./cpu_banded.out $i $j
  ./cpu_banded.out $i $j
  ./cpu_banded.out $i $j
  echo
  
  export j=`echo $i/10 | bc`
  ./cpu_banded.out $i $j
  ./cpu_banded.out $i $j
  ./cpu_banded.out $i $j
  echo

  export j=`echo $i/5 | bc`
  ./cpu_banded.out $i $j
  ./cpu_banded.out $i $j
  ./cpu_banded.out $i $j
  echo
done
echo ===========================

echo Phi Banded
echo ===========================
for i in {500..20000..500}
do
  export j=`echo $i/20 | bc`
  ./phi_banded.out $i $j
  ./phi_banded.out $i $j
  ./phi_banded.out $i $j
  echo
  
  export j=`echo $i/10 | bc`
  ./phi_banded.out $i $j
  ./phi_banded.out $i $j
  ./phi_banded.out $i $j
  echo

  export j=`echo $i/5 | bc`
  ./phi_banded.out $i $j
  ./phi_banded.out $i $j
  ./phi_banded.out $i $j
  echo
done
echo ===========================

echo CPU DENSE
echo ===========================
for i in {500..20000..500}
do
  ./cpu_dense.out $i
  ./cpu_dense.out $i
  ./cpu_dense.out $i
  echo
done
echo ===========================

echo Phi DENSE
echo ===========================
for i in {500..20000..500}
do
  ./phi_dense.out $i
  ./phi_dense.out $i
  ./phi_dense.out $i
  echo
done
echo ===========================
