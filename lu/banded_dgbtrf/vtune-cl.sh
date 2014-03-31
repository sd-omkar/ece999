#!/bin/bash

export OMP_NUM_THREADS=40 
export KMP_AFFINITY=granularity=fine,compact,1,0 
export MIC_OMP_NUM_THREAD=236 
export MIC_KMP_AFFINITY=explicit,granularity=fine,proclist=[1-236:1] 
export MIC_ENV_PREFIX=MIC_ 

amplxe-cl -allow-multiple-runs -search-dir=src=./ -collect concurrency -r results/concurrency ./cpu.out 5000 500
amplxe-cl -report hotspots -report-output=results/concurrency.rep -r results/concurrency -search-dir=src=./

amplxe-cl -allow-multiple-runs -search-dir=src=./ -collect hotspots -r results/hostspots ./cpu.out 5000 500
amplxe-cl -report hotspots -report-output=results/hotspots.rep -r results/hostspots -search-dir=src=./

