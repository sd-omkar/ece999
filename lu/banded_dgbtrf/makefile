all: cpu phi

cpu:
	icc dgbtrf_cpu.c -o cpu.out -xHOST -O3 -no-prec-div -mkl -parallel -opt-prefetch=2 -opt-streaming-stores auto -fp-model fast

phi:
	icc dgbtrf_phi.c -o phi.out -xHOST -O3 -no-prec-div -mkl -parallel -opt-prefetch=2 -opt-streaming-stores auto -fp-model fast

offload:
	icc dgbtrf.c -o offload.out -openmp -mkl -O3
