Meissel Lehmer Sieve

This project, long stopped, was an attempt to improve the Meissel-Lehmer sieve to that large values of pi(x)
could be computed quickly. As I recall 10^20 was the highest I got. The main thing I think is of interest here
is the program ML.c in the prototype folder. It keeps all the data in memory so can't go much beyond 10^14,
but I think reading this code would give a good understanding of how the algorithm works.

-------------------------------------------------------------------------------------------------------------

Implementation of calculation of pi(x) by the Meissel-method as described in 
Lagarias, Miller, and Odlyzko, Computing π(x): the Meissel-Lehmer Method.

www.dtc.umn.edu/~odlyzko/doc/arch/meissel.lehmer.pdf

All programs assume a 64 bit machine and are expected to fail on a 32 bit
target. They have been successfully compiled and run under Ubuntu and Redhat.
