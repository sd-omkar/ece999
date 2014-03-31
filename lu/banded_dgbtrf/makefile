all: cpu phi native

cpu:
	icc dgbtrf_cpu.c -o cpu.out -xHOST -O3 -no-prec-div -mkl -parallel -opt-prefetch=2 -opt-streaming-stores auto -fp-model fast

phi:
	icc dgbtrf_phi.c -o phi.out -xHOST -O3 -no-prec-div -mkl -parallel -opt-prefetch=2 -opt-streaming-stores auto -fp-model fast

native:
	icc dgbtrf_native.c -o native.out -mmic -O3 -no-prec-div -mkl -parallel -opt-prefetch=2 -opt-streaming-stores auto -fp-model fast

profile:
	icc dgbtrf_cpu.c -o profile_cpu.out -xHOST -O3 -no-prec-div -mkl -opt-prefetch=2 -opt-streaming-stores auto -fp-model fast -profile-functions -profile-loops=all
