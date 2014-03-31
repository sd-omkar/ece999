#!/bin/bash

export KMP_AFFINITY=compact

for i in {1..40}
do
export OMP_NUM_THREADS=$i
echo OMP_NUM_THREADS = $i
./speedometer ./dgetrf.out 1000
done

unset OMP_NUM_THREADS
unset KMP_AFFINITY
