CC = gcc

CFLAGS = -Wall -O3 -std=c99 -pedantic

OBJ = cephes_eigens.o GenData.o GenPars.o Matrix.o Overlap.o Print.o QF.o Random.o symmEigValDec.o

CARP:	CARP.c $(OBJ)
	$(CC) $(CFLAGS) CARP.c -o CARP $(OBJ) -lm
	$(CC) $(CFLAGS) C-MixSim.c -o C-MixSim $(OBJ) -lm
	$(CC) $(CFLAGS) AdjRand.c -o AdjRand -lm
	$(CC) $(CFLAGS) hierclust.c -o hierclust -lm
	$(CC) $(CFLAGS) unittest.c -o unittest $(OBJ) -lm

clean:
	rm -rf *.o CARP C-MixSim AdjRand hierclust unittest

check:
	./unittest
