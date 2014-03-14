all: cpu phi

cpu:
	icc dgbsv_cpu.c -o cpu.out -xHOST -O3 -no-prec-div -mkl -parallel -opt-prefetch=2 -opt-streaming-stores auto -fp-model fast

phi:
	icc dgbsv_phi.c -o phi.out -xHOST -O3 -no-prec-div -mkl -parallel -opt-prefetch=2 -opt-streaming-stores auto -fp-model fast
