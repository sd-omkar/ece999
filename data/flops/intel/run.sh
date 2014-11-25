#!/bin/bash

for i in {500..5000..500}
do
  vtune -collect -allow-multiple-runs snb-general-exploration -- ./a.out $i
done
