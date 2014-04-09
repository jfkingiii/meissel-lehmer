//This is the parallel implementation. A range of block numbers is supplied in the arguments and multiple
//parallel instances are run to calculate the contributions of all the blocks
//and the results are printed to a file. S2_short is run to calculate the short block contribution
//and S2_total then reads all the files and calculates S2.
//TO DO: Change the sieveing method to a bit sieve mod 6.

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include "gmp.h"

uint64_t sqr(uint64_t n) {return n*n;}

uint64_t max(uint64_t a, uint64_t b) {return a>b ? a:b;}
uint64_t min(uint64_t a, uint64_t b) {return a<b ? a:b;}


uint64_t lpow(mpz_t n, uint64_t k, uint64_t l)
{
  mpz_t x;
  mpz_init(x);
  mpz_set(x,n);
  mpz_pow_ui(x, x, k);
  mpz_root(x, x, l);
  uint64_t  m = mpz_get_ui(x);
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

int main(int argc, char *argv[])
{
  if(argc != 5)
  {
    printf("Four arguments required\n");
    exit(1);
  }
  mpz_t n, buffer1, buffer2;
  mpz_init(n); mpz_init(buffer1); mpz_init(buffer2);
  mpz_ui_pow_ui(n, atol(argv[1]), atol(argv[2]));
  uint64_t start_block = atol(argv[3]);
  uint64_t finish_block = atol(argv[4]);
  //uint64_t frn = lpow(n, 1, 4);
  uint64_t crn = lpow(n, 1, 3);
  //uint64_t srn = lpow(n, 1, 2);
  uint64_t ttrn = lpow(n, 2, 3);

  size_t psize = (size_t)(((1.26*crn)/log(crn) + 1) * sizeof(uint64_t));
  uint64_t* primes = (uint64_t*)malloc(psize);
  assert(primes!=NULL);
  uint64_t pi_crn = sieve(crn, primes);

  int8_t* mu = malloc(sizeof(int8_t) * (crn+1));                            //mu will run from 1 to crn with s[0] not used
  assert(mu != NULL);
  memset(mu, 1, sizeof(int8_t) * (crn+1));

  uint64_t* lpf = malloc(sizeof(int64_t) * (crn+1));                        //lpf will run from 1 to crn with s[0] not used
  assert(lpf != NULL);
  memset(lpf, 0, sizeof(int64_t) * (crn+1));

  uint64_t i, j;

  for(i=1; i<=pi_crn; i++)
    for(j=primes[i]; j<=crn; j+=primes[i])
    {
      mu[j] = -mu[j];
      if(lpf[j] == 0) lpf[j] = primes[i];
    }
  for(i=1; i<=pi_crn; i++)
    for(j=sqr(primes[i]); j<=crn; j+=sqr(primes[i])) mu[j]=0;               //remove numbers that are not squarefree

  uint64_t blocksize = crn;
  uint64_t num_blocks = ttrn/blocksize;
  assert(start_block <= finish_block);
  assert(finish_block <= num_blocks);
  //  printf("Calculating blocks %lu to %lu of %lu...\n", start_block, finish_block, num_blocks);

  mpz_t S2_result;
  mpz_init(S2_result);

  uint8_t* s = malloc((crn+1)*sizeof(uint8_t));
  assert(s != NULL);

  uint64_t* phi = malloc((pi_crn+1)*sizeof(uint64_t));
  memset(phi, 0, (pi_crn+1)*sizeof(uint64_t));

  int64_t* cum_mu = malloc((pi_crn+1)*sizeof(int64_t));
  memset(cum_mu,0, (pi_crn+1)*sizeof(int64_t));

  uint64_t p, f, m, k, start, q, block_id, first_m, last_m, L, U, first_f;

  for(block_id = start_block; block_id <= finish_block; block_id++)
  {
    memset(s, 1, sizeof(uint8_t) * (crn+1));
    L = (block_id -1)*blocksize + 1;
    U = L + blocksize;
    for(i=1; i<=pi_crn; i++)
    {
      p = primes[i];
      mpz_fdiv_q_ui(buffer1, n, L*p);
      mpz_fdiv_q_ui(buffer2, n, U*p);
      first_m = mpz_get_ui(buffer1);
      last_m = mpz_get_ui(buffer2);
      if(first_m > crn) first_m = crn;
      if(last_m > first_m) last_m = first_m;
      if(p > first_m) break;
      start = L;
      for(m = first_m; m > last_m && p*m > crn ; m--)
      {
	if(mu[m] != 0 && p < lpf[m] && p*m > crn && m <= crn)
	{
	  cum_mu[i] += -mu[m];
	  mpz_fdiv_q_ui(buffer1, n, p*m);
	  q = mpz_get_ui(buffer1);
	  assert(q>=L);
	  assert(q<=U);
	  for(k=start; k<=q; k++) phi[i] += s[k-L];
	  start = q+1;
	  if(mu[m] > 0) mpz_sub_ui(S2_result, S2_result, mu[m] * phi[i]);
	  else mpz_add_ui(S2_result, S2_result, -mu[m] * phi[i]);
	}
      }
      for(k = start; k < U; k++) phi[i] += s[k - L];
      first_f = L%p==0 ? L/p : L/p +1;
      for(f = first_f; p*f < U; f++)  s[p*f - L]=0;
    }
  }
  //  printf("\n");

  // printf("\nS2 = ");
  mpz_out_str (stdout, 10, S2_result);
  printf("\n");

  for(i=1; i<=pi_crn; i++) printf("%lu\t%lu\t%ld\n", primes[i], phi[i], cum_mu[i]);

  free(s); free(phi);
  free(primes); free(mu);
  free(lpf); free(cum_mu);
  return(0);
}
