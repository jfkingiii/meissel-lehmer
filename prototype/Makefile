CC = gcc
CFLAGS = -Wall -lm -lgmp -O3 -funroll-all-loops

.PHONY: all
all: ~/bin/ML ML_test validate_ML

~/bin/ML: ML.c
	$(CC) ML.c $(CFLAGS) -o ~/bin/ML

ML_test:
	$(CC) ML_test.c $(CFLAGS) -o ./ML_test

validate_ML:
	$(CC) validate_ML.c $(CFLAGS) -o ./validate_ML

.PHONY: clean
clean:
	rm ML_test validate_ML ~/bin/ML
