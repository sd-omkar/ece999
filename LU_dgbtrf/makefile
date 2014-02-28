all: host native

trf:
	icc dgbtrf.c -o trf.out -xHOST -O3 -no-prec-div -mkl -parallel -opt-prefetch=2 -opt-streaming-stores auto -fp-model fast
