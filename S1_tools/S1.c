#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include "gmp.h"

uint64_t sqr(uint64_t n) {return n*n;}

uint64_t lpow(mpz_t n, uint64_t k, uint64_t l)
{
  //Returns floor(n^(k/l))
  mpz_t x;
  mpz_init(x);
  mpz_set(x,n);
  mpz_pow_ui(x, x, k);
  mpz_root(x, x, l);
  uint64_t  m = mpz_get_ui(x);
  //Check
  /*  mpz_t M, M1, N, C1, C2;
 mpz_inits(M, M1, N, C1, C2, NULL);
 mpz_set_ui(M, m); mpz_set_ui(M1, m+1); mpz_set_ui(N, n);
 mpz_pow_ui(C1, M, l); mpz_pow_ui(C2, N, k);
 assert(mpz_cmp(C1, C2) <= 0);
 mpz_pow_ui(C1, M1, l);
 assert(mpz_cmp(C1, C2) > 0);

 mpz_clears(x, M, M1, N, C1, C2, NULL );*/
  return m;
}

uint64_t sieve(uint64_t N, uint64_t* p)
{
  uint64_t M = N/2 + N%2;
  uint64_t i, j;
  uint64_t sieve_to = ceil((sqrt(N)-1)/2);
  assert(sqr(2*sieve_to + 1 )>= N);
  uint64_t pi = 1;

  /*contents of S[i] will be 1 if 2*i + 1 is prime, zero otherwise*/
  uint8_t* S = malloc(sizeof(uint8_t) * M);
  assert(S!=NULL);
  memset(S, 1, sizeof(uint8_t)*M);
  S[0]=0;     /*1 is not prime*/

  for(i=1; i<=sieve_to; i++)
    if (S[i]==1)
      for(j=2*i*(i+1); j<M; j+=2*i+1) S[j]=0;

  p[1]=2;
  for(i=1; 2*i+1<=N; i++)
    if(S[i]==1)
      {pi++; p[pi]= 2*i+1;}
  free(S);
  return(pi);
}

void S1(mpz_t n, uint64_t crn, int8_t* mu, mpz_t result)
{
  mpz_t tmp;
  mpz_inits(result, tmp, NULL);
  uint64_t i;
  for(i=1; i<=crn; i++)
    {
      mpz_tdiv_q_ui(tmp, n, i);
      mpz_mul_si(tmp, tmp, mu[i]);
      mpz_add(result, result, tmp);
    }
}


main(long argc, char *argv[])
{
  mpz_t n, buffer1, buffer2;
  mpz_init(n); mpz_init(buffer1); mpz_init(buffer2);
  if(argc == 2)
    mpz_set_str(n, argv[1], 10);
  if(argc == 3)
    mpz_ui_pow_ui(n, atol(argv[1]), atol(argv[2]));
  if (argc<2 || argc>3)
    {
      printf("Incorrect input format\n");
      exit(1);
    }
  printf("Using GMP...\n");
  uint64_t crn = lpow(n, 1, 3);
  size_t psize = (1.26*crn)/log(crn)*sizeof(uint64_t) + 1;                  //add 1 because primes array will run from 1 to pi(crN) with primes[0] not used
  uint64_t* primes = malloc(psize);
  assert(primes!=NULL);
  uint64_t i,j;

  uint64_t J = sieve(crn ,primes);

  int8_t* mu = malloc(sizeof(int8_t) * (crn+1));                        //mu will run from 1 to crN with s[0] not used
  assert(mu != NULL);
  memset(mu, 1, sizeof(int8_t) * (crn+1));

  uint64_t* lpf = malloc(sizeof(int64_t) * (crn+1));                        //lpf will run from 1 to crN with s[0] not used
  assert(lpf != NULL);
  memset(lpf, 0, sizeof(int64_t) * (crn+1));

  for(i=1; i<=J; i++)
    for(j=primes[i]; j<=crn; j+=primes[i])
      {
	mu[j] = -mu[j];
	if(lpf[j] == 0) lpf[j] = primes[i];
      }
  for(i=1; i<=J; i++)
    for(j=sqr(primes[i]); j<=crn; j+=sqr(primes[i])) mu[j]=0;           //remove numbers that are not squarefree

  mpz_t S1_result;
  mpz_init(S1_result);
  printf("Calculating S1..."); fflush(stdout);
  S1(n, crn, mu, S1_result);
  printf("\nS1 = ");
  mpz_out_str (stdout, 10, S1_result);
  printf("\n\n");

  free(primes); free(mu); free(lpf);
}
