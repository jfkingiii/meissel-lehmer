CC = gcc
CFLAGS = -Wall -lm -lgmp -O3 -funroll-all-loops -mpopcnt

.PHONY: all
all: S2 num_blocks S2_parallel_blocks S2_total

S2:
	$(CC) S2.c $(CFLAGS) -o ~/bin/S2

num_blocks:
	$(CC) num_blocks.c $(CFLAGS) -o ~/bin/num_blocks

S2_parallel_blocks:
	$(CC) S2_parallel_blocks.c $(CFLAGS) -o ~/bin/S2_parallel_blocks

S2_total:
	$(CC) S2_total.c $(CFLAGS) -o ~/bin/S2_total

.PHONY: clean
clean:
	rm  ~/bin/S2 ~/bin/num_blocks ~/bin/S2_parallel_blocks ~/bin/S2_total
