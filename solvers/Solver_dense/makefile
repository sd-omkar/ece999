all: host native

host:
	icc -03 -xhost -w -mkl -vec-report1 dgesv.c -o host.out

native:
	icc -O3 -mmic -w -mkl -vec-report1 dgesv.c -o native.out
