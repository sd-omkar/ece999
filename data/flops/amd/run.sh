#!/bin/bash

for i in {500..5000..500}
do
  sudo opcontrol --reset
  sudo opcontrol -e=DISPATCHED_FPU_OPS:500 --no-vmlinux
  sudo opcontrol --start
  ./a.out $i
  sudo opcontrol --stop
  sudo opcontrol --dump
  opreport -l ./a.out > oprofile_$i
done
