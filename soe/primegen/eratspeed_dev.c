#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#define MOD 30					//Mod 30 wheel
#define GROUPS 8				//8 residue classes coprime to 30
#define B32 10001
#define B (B32 * 32)				//Segment size

uint32_t a[GROUPS][B];				//The sieve
uint64_t dtab[GROUPS] = { 1, 7, 11, 13, 17, 19, 23, 29};
uint64_t *next[GROUPS];
uint64_t total = 10;				// 2, 3, 5, 7, 11, 13, 17, 19, 23, 29 

uint64_t small_sieve(uint64_t N, uint64_t* p)	//Generates primes to sieve with
{
						//Initialize
    assert(N > 1);
    uint64_t M = (N-1)/2 + 1;
    uint64_t pi;
    uint64_t i, f;
    uint8_t* S = malloc(sizeof(uint8_t) * M);
    assert(S != NULL);
    memset(S, 1, sizeof(uint8_t)*M);
    S[0]=0;                 

						//Sift
    for(i=1; i*i<=N; i+=2)
	if (S[(i-1)/2] == 1)			//Only sieve by primes
	    for(f = i; f*i <= N; f += 2) S[(f*i - 1)/2] = 0;

						//Count and record primes
    pi = 0;
    for(i=3; 2*i+1<=N; i++)			//Treat 7 as the first prime since wheel is mod 30
	if(S[i]==1)
	{
//	    printf("%lu\t%lu\n",pi,N);
	    p[pi]= 2*i+1;
	    pi++; 
	}

    free(S);
    return(pi);
}


void init(uint64_t L, uint64_t num_primes, uint64_t* primes)
{
/*
constructs next[8][3509]
The 8 rows are the 8 numbers d less than and coprime to 30
the 3509 columns are the first 3509 primes >= 7

Calculate qz = smallest multiple of q (>= q*q) so that qz%30 == d
next[i][j] = number of blocks of size 30 that qz - d is above L

Example: residue = 13, prime = 17

q*q = 289
(289 + 12*17) = 493  = 13 (mod 30)
(493 - 13)/30 - 1 = 15

493 = (15 + 1) * 30 + 13

So to sieve multiples of 17 which are congruent to 13 mod 30, sieve

493, 493 + 30*17, 493 + 60 * 17, etc.

*/
  uint64_t i, j, d, q, qz;
  
  for (i = 0; i < GROUPS; ++i) 
  {
    d = dtab[i];
    for (j = 0;j < num_primes;++j) 
    {
      q = primes[j];
      qz = q * q;
      while (qz % MOD != d) qz += q + q;
      next[i][j] = (qz - d)/MOD - L;
    }
  }
#if 0
  for(i=0; i<8; i++)
      for(j=0; j<10; j++)
	  printf("%lu			\t%lu\tnext: %lu\n", dtab[i], primes[j], next[i][j]);
#endif
}

void sift_segment(uint64_t num_primes, uint64_t* primes)
{
  int i,j;
  uint64_t *nexti;
  register uint64_t k;
  register uint32_t *buf;
  register uint64_t q;
						//array a is 8 rows (the residues mod 30) by 1001 columns
  for (i = 0;i < GROUPS;++i) 
  {
      buf = a[i];				//buf is the address of the row of a with data on the ith residue
      nexti = next[i]; 				//nexti is the row of next giving the next multiple congruent to dtab[i] mod 30
      memset(buf, 0, B32*sizeof(uint32_t));
      for (j = 0; j < num_primes; ++j)		//loop over the 3509 primes in qtab
      {	
	  k = nexti[j];				//??(k+1)*q is first multiple of qtab[j] that must be sieved
	  if (k < B) 
	  {					//B=32,032
	      q = primes[j];
	      do 
	      {
		  buf[k/32] |= 1 << k%32;       //Or  buf[k/32] |= 1 << k%32; 
		  k += q;			//next multiple of q to be sieved
	      } 
	      while (k < B);
	  }
	  nexti[j] = k - B;
      }
  }
}

void count_segment()
{
  int i;
  register uint32_t *ai;
  register uint64_t pos;
  register uint32_t bits;
  register uint64_t result;

/* To print the primes (slowly), given L:
  for (k = 0;k < B;++k)
    for (i = 0;i < 8;++i)
      if (!(a[i][k / 32] & two[k & 31]))
	printf("%d			\n",30 * (L + k) + dtab[i]);
*/

  result = 0;
  uint64_t s_size = 0;
  for (i = 0; i < GROUPS; ++i) 
  {
    ai = a[i];
    for (pos = 0; pos < B32; ++pos) 
    {
      s_size += 32;
      bits = ~ai[pos];
      result += __builtin_popcountl(bits);
    }
  }
  total += result;
  //printf("Sieve size calcluated: %lu\n", s_size);
}

int main(int argc, char *argv[])
{
    uint64_t U;
    if(argc == 2) U = atol(argv[1]);
    else if(argc == 3) U = roundl(pow(atol(argv[1]), atol(argv[2])));
    else if(argc == 4) U = atol(argv[1]) * roundl(pow(atol(argv[2]), atol(argv[3])));
    else
  {
      printf("Usage: %s x returns pi(x)	\n", argv[0]);
      exit(1);
  }

  uint64_t L = 1;
  uint64_t sqrt_U = ceil(sqrt(U));
  size_t psize = 1.26 * sqrt_U/log(sqrt_U) * sizeof(uint64_t);
  uint64_t* sift_primes = malloc(psize);
  assert(sift_primes != NULL);
  uint64_t I = small_sieve(sqrt_U, sift_primes);
  int i;

/*
  if(U < 9600990) 
  {
      printf("Upper limit is smaller than sieve size\n");
      exit(1);
  }
*/

  for(i=0; i<GROUPS; i++) 
  {
     next[i] = malloc(I * sizeof(uint64_t));
     assert(next[i] != NULL);
  } 
  init(L, I, sift_primes);

  do 
  {
    sift_segment(I, sift_primes);
    count_segment();
    L += B;
  } 
  while (L < ((U/30)/B)*B);

  printf("%lu primes up to %lu\n", total, 30 * L);

  exit(0);
}
