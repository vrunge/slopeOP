all: coeff.update.so prune2R.so
.PHONY: all clean

coeff.update.so: coeff.updateR.c
	R CMD SHLIB coeff.updateR.c

prune2R.so: prune2R.c
	R CMD SHLIB prune2R.c      

clean:
	\rm -f *.so *.o
