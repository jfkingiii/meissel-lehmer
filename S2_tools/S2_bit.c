#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include "gmp.h"

#define DEBUG 0

uint64_t sqr(uint64_t n) {return n*n;}

uint64_t lpow(mpz_t n, uint64_t k, uint64_t l)
{
  //Returns floor(n^(k/l))
  assert(n>0); assert(l>0);
  mpz_t x;
  mpz_init(x);
  mpz_set(x,n);
  mpz_pow_ui(x, x, k);
  mpz_root(x, x, l);
  uint64_t y = mpz_get_ui(x);
  mpz_clear(x);
  return y;
}

uint64_t sieve(uint64_t N, uint64_t* p)
{
  assert(N>0); assert(p!=NULL);
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


void clear_bit(uint64_t* s, uint64_t j)
{
  uint64_t ind = j/64;
  uint64_t bit = j%64;
  s[ind] &= ~(1ul<<bit);
}

uint64_t count_bits(uint64_t *s, uint64_t first_bit, uint64_t last_bit)
{
  if(first_bit > last_bit) return 0;
  uint64_t first_block = first_bit/64;
  uint64_t last_block = last_bit/64;
  uint64_t i, j=0;
  uint64_t mask1 = (1ul<<(last_bit%64+1))-1;
  if (last_bit%64 + 1 == 64) mask1 = ~0;
  uint64_t mask2 = ~((1ul<<first_bit%64)-1);
  if(first_block==last_block) return __builtin_popcountl(s[first_block] & mask1 & mask2);
  if(last_block - first_block >=2) for(i=first_block+1; i<=last_block-1; i++) j += __builtin_popcountl(s[i]);
  j +=  __builtin_popcountl(mask1 & s[last_block]);
  j +=  __builtin_popcountl(mask2 & s[first_block]);
  return j;
}


void sieve_step(uint64_t start, uint64_t U, uint64_t* phi, uint64_t* s, uint64_t L, uint64_t p)
{
  uint64_t f;
  uint64_t first_f = L%p==0 ? L/p : L/p +1;
  for(f = first_f; p*f < U; f++)
  {
      assert(p*f-L>=0);
      s[(p*f -L)/64] &= ~(1ul<<((p*f-L)%64));
  }
}


int main(long argc, char *argv[])
{
  if(argc != 4)
  {
      printf("three arguments required\n");
      exit(1);
  }
  mpz_t n, buffer1, buffer2;
  mpz_init(n); mpz_init(buffer1); mpz_init(buffer2);
  mpz_ui_pow_ui(n, atol(argv[2]), atol(argv[3]));
  mpz_mul_ui(n, n, atol(argv[1]));

  uint64_t crn = lpow(n, 1, 3);
  uint64_t srn = lpow(n, 1, 2);
  uint64_t ttrn = lpow(n, 2, 3);
  uint64_t start_block = 1;
  uint64_t finish_block = ttrn/crn;

  size_t psize = (size_t)(((1.26*crn)/log(crn) + 1) * sizeof(uint64_t));
  uint64_t* primes = (uint64_t*)malloc(psize);
  assert(primes!=NULL);
  uint64_t pi_crn = sieve(crn, primes);

  int8_t* mu = malloc(sizeof(int8_t) * (crn+1));                            //mu will run from 1 to crn with s[0] not used
  assert(mu != NULL);
  memset(mu, 1, sizeof(int8_t) * (crn+1));

  uint64_t* lpf = malloc(sizeof(uint64_t) * (crn+1));                        //lpf will run from 1 to crn with s[0] not used
  assert(lpf != NULL);
  memset(lpf, 0, sizeof(uint64_t) * (crn+1));

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

  mpz_t S2_result;
  mpz_init(S2_result);

  uint64_t* s = malloc((crn/64+2)*sizeof(uint64_t));
  assert(s != NULL);

  uint64_t* phi = malloc((pi_crn+1)*sizeof(uint64_t));
  memset(phi, 0, (pi_crn+1)*sizeof(uint64_t));

  uint64_t p, f, m, k, start, q, block_id, first_m, last_m, L, U, first_f;

  for(block_id = start_block; block_id <= finish_block; block_id++)
    {
      memset(s, ~0, sizeof(uint64_t) * (crn/64+2));
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
		  mpz_fdiv_q_ui(buffer1, n, p*m);
		  q = mpz_get_ui(buffer1);
		  assert(q>=L);
		  assert(q<=U);

		  phi[i] += count_bits(s, start - L, q - L);

		  start = q+1;
		  if(mu[m] > 0) mpz_sub_ui(S2_result, S2_result, mu[m] * phi[i]);
		  else mpz_add_ui(S2_result, S2_result, -mu[m] * phi[i]);
		}
	    }
	  phi[i] += count_bits(s, start - L, U -1 - L);
	  sieve_step(start, U, phi,s, L, p);
	}
    }

  uint64_t short_block_length = ttrn % blocksize;

  if(short_block_length>0)
    {
      //short block
      L = L + blocksize;
      U = L + short_block_length; 

      memset(s, ~0, sizeof(uint64_t) * (crn/64+2));
      for(i=1; i<=pi_crn; i++)
	{
	  start = L;
	  p = primes[i];
	  mpz_fdiv_q_ui(buffer1, n, L*p);
	  mpz_fdiv_q_ui(buffer2, n, U*p);
	  first_m = mpz_get_ui(buffer1);
	  last_m = mpz_get_ui(buffer2);
	  if(first_m > crn) first_m = crn;
	  if(last_m > first_m) last_m = first_m;
	  if(p > first_m) break;
	  for(m=first_m; m>last_m && p*m>crn ; m--)
	  {
	      if(mu[m]!=0 && p<lpf[m] && p*m>crn && m<=crn)
	      {
		  mpz_fdiv_q_ui(buffer1, n, p*m);
		  q = mpz_get_ui(buffer1);
		  assert(q>=L);
		  assert(q<=U);
		  phi[i] += count_bits(s, start - L, q - L);
		  start = q+1;
		  if(mu[m] > 0) mpz_sub_ui(S2_result, S2_result, mu[m]*phi[i]);
		  else mpz_add_ui(S2_result, S2_result, -mu[m]*phi[i]);
	      }
	  }
	  sieve_step(start, U, phi,s, L, p);
	}
    }

  mpz_out_str (stdout, 10, S2_result);
  printf("\n");
  mpz_clears(buffer1, buffer2, n, S2_result, NULL);
  free(s); free(phi);
  free(primes); free(mu);
  free(lpf); 
  return(0);
}
