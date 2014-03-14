export OFFLOAD_REPORT=2

export MKL_MIC_DISABLE_HOST_FALLBACK=1 

export OMP_NUM_THREADS=1 

export KMP_AFFINITY=granularity=fine,compact,1,0 

export MIC_OMP_NUM_THREAD=236 

export MIC_KMP_AFFINITY=explicit,granularity=fine,proclist=[1-236:1] 

export MIC_ENV_PREFIX=MIC_ 

