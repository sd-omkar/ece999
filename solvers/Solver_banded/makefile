all: host native

trf:
	icc -03 -xhost -w -mkl -opt-prefetch=2 dgbtrf.c -o trf.out

trf2:
	icc dgbtrf.c -o trf.out -xHOST -O3 -no-prec-div -mkl -parallel -opt-prefetch=2 -opt-streaming-stores auto -fp-model fast

host:
	icc -03 -xhost -w -mkl -vec-report1 dgbsv.c -o host.out

native:
	icc -O3 -mmic -w -mkl -vec-report1 dgbsv.c -o native.out
