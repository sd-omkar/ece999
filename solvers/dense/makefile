all: cpu phi

cpu:
	icc dgesv_cpu.c -o cpu.out -xHOST -O3 -no-prec-div -mkl -parallel -opt-prefetch=2 -opt-streaming-stores auto -fp-model fast

phi:
	icc dgesv_phi.c -o phi.out -xHOST -O3 -no-prec-div -mkl -parallel -opt-prefetch=2 -opt-streaming-stores auto -fp-model fast
