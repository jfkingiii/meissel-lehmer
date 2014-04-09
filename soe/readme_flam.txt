Calling prime_sieve with no parameter prints a default usage message and exits.

It needs at least one parameter, called 'stop', the largest number to check for
to-be-prime. Then it will generate all prime numbers starting from 1 up to
this 'stop'. If a second parameter is given, called 'start', the interval
searched for primes will be only [start,stop].
A further, third parameter, called 'size', can be given indicating the size of
the sieve to be used. This will strongly influence the run time of the program.
The optimal value is mostly a value slightly smaller than the 1st level data
chache size of the CPU. For very large 'stop' values this has to be extended to
the minimum sieve size which is about Sqrt(stop)/24 but the program will check
on start whether the given 'size' is sufficient and if not aborts and prints
the minimal needed size.


If you compile prime_sieve.c there are a few important defines to pay attention:
LONG
_ANSI_EXTENSION
true_64bit_words

If you define nothing of these three then all internal variables will be
'unsigned' which will probably mean they are 32-bit variables.
Then you can use 'stop' at most 2^32-1.

If you want to handle primes >= 2^32 you need to define  LONG=uns64
This pseudo datatype uns64 is either  unsigned long (the default)
or if you define  _ANSI_EXTENSION  then it will become  unsigned long long.
What you need will depend on your compiler.
If your CPU and compiler handles true 64-bit register variables you may
also want to define  true_64bit_words

Lastly, there is a macro  use(x)
Its argument x is the next found prime number.
Its default function is to print this x to stdout and increase the number of
total found primes by 1. You can replace the macro body by code of your likeness
to process the found prime as you want.


Runtime Storage requirements:

max{sieve_size,Sqrt(stop)/24} + (sizeof(LONG)+1)*Pi(Sqrt(stop))  bytes
