#!/bin/bash

export OMP_NUM_THREADS=40 
export KMP_AFFINITY=granularity=fine,compact,1,0 
export MIC_OMP_NUM_THREAD=236 
export MIC_KMP_AFFINITY=explicit,granularity=fine,proclist=[1-236:1] 
export MIC_ENV_PREFIX=MIC_ 

echo ===========================
echo CPU DENSE
echo ===========================
for i in {10500..20000..500}
do
  ./cpu_dense.out $i
  ./cpu_dense.out $i
  ./cpu_dense.out $i
  echo
done
echo ===========================

echo ===========================
echo Phi DENSE
echo ===========================
for i in {10500..20000..500}
do
  ./phi_dense.out $i
  ./phi_dense.out $i
  ./phi_dense.out $i
  echo
done
echo ===========================
