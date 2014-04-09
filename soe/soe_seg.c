#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>

#define NEXT_ODD(u) u + (1 - u%2)
#define PREV_ODD(u) u - (1 - u%2)

uint64_t sieve(uint64_t N, uint64_t* p)
{
  assert(N > 1);
  uint64_t M = (N-1)/2 + 1;
  uint64_t i, f;

  /*contents of S[i] will be 1 if 2*i + 1 is prime, zero otherwise*/
  uint8_t* S = malloc(sizeof(uint8_t) * M);
  assert(S!=NULL);
  memset(S, 1, sizeof(uint8_t)*M);
  S[0]=0;                                         //1 is not prime
  for(i=1; i*i<=N; i+=2)
    if (S[(i-1)/2]==1)                            //sift only by primes
      for(f=i; f*i<=N; f+=2) S[(f*i-1)/2]=0;

  p[1]=2;
  uint64_t pi = 1;
  for(i=1; 2*i+1<=N; i++)
    if(S[i]==1)
      {pi++; p[pi]= 2*i+1;}
  free(S);
  return(pi);
}

uint64_t sieve_segment(uint64_t L, uint64_t U, const uint64_t* sift_primes, uint64_t I, uint8_t* s)
{
    /*                         
    L = lower bound of interval
    U = upper bound of interval
    (L and U must be odd)		
    sift_primes = array of primes to sieve with, with p[1] = 2, p[2] = 3, etc.
    I = length of 
    s = sieving array - reused for each segment  
    */

    assert(L%2==1 && U%2==1);				//for simplicity assume upper and lower limits are always odd
    if(L>U) return 0;

    uint64_t p;						//sieving prime
    uint64_t i;						//sieving prime index in sift_primes
    uint64_t f;						//multiple of p to be crossed out
    uint64_t start;					//starting value for f
    uint64_t pi=0;					//accumulator for value of pi(x)

    memset(s, 1, sizeof(uint8_t) * ((U - L)/2+1));

    /*L and U are assumed to be odd 
      L+2*0, L+2*1, L+2*2, ..., L+ 2*(U-L)/2 -> s[0], s[1], ...,s[(U - L)/2]                                                          s[i] is an indicator for L + 2*i*/

    for(i=2; i<=I; i++)					//loop through the sieving primes
    {
	p = sift_primes[i];
	if (p*p > U) break;

	/*Determine first factor to sift by:
	  start = smallest value of f so that f*p >= L*/

	start = L%p==0 ? L/p : L/p +1;
	assert(start*p>=L); assert((start-1)*p<L);
	if(start%2==0) start++;
	if (start<p) start=p;				//can start crossing out at p^2                                                                                                 
	for(f=start*p; f<=U; f+= 2*p) s[(f - L)/2]=0;	//sift the prime p
    }

    for (i = 0; i<= (U - L)/2; i++) pi += s[i];		//count the unsieved values
    return pi;
}


#define SIEVE_INCR 256000

int main(int argc, char *argv[])
{
  uint64_t L, U;
  if(argc==3)
  {
    L = atol(argv[1]);
    U = atol(argv[2]);
  }
  else if(argc==2)
  {
    L=1;
    U = atol(argv[1]);
  }
  else
  {
      printf("usage\n");
      exit(1);
  }

  assert(L>=1);
  assert(U>=L);
  uint64_t sqrt_U = ceil(sqrt(U));
  size_t psize = (1.26*sqrt_U/log(sqrt_U)+1)*sizeof(uint64_t);
  uint64_t* sift_primes = malloc(psize);
  assert(sift_primes != NULL);
  uint64_t I = sieve(sqrt_U, sift_primes);

  uint64_t sieve_to, sieve_from, pi = 1;
  uint8_t* s = malloc((SIEVE_INCR/2 + 1) * sizeof(uint8_t));
  assert(s != NULL);

  sieve_from = L;
  sieve_to = L + SIEVE_INCR -1;
  if (sieve_to > U) sieve_to = U;
  uint64_t cnt = 0;
  for(;;)
  {
    #if DEBUG 
      printf("%lu\t%lu\n",sieve_from, sieve_to);
    #endif
    pi += sieve_segment(NEXT_ODD(sieve_from), PREV_ODD(sieve_to), sift_primes, I, s);
    cnt++;
    if(sieve_to == U) break;
    sieve_from += SIEVE_INCR;
    sieve_to += SIEVE_INCR;
    if(sieve_to > U) sieve_to = U;
  }
  printf("%lu segments sieved\n", cnt);
  free(s);
  free(sift_primes);
  printf("%lu\n",pi - 1);
  return(0);
}

