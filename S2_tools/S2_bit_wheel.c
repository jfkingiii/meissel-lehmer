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

//experimental mod 30 wheel
const int64_t PREV[30] = {-1, 0, -1, -2, -3, -4, -5, 0, -1, -2,
			  -3, 0, -1,  0, -1, -2, -3, 0, -1,  0,
			  -1, -2, -3, 0, -1, -2, -3, -4, 0, -1};

const int64_t NEXT[30] = {1,0,5,4,3,2,1,0,3,2,1,0,1,0,3,2,1,0,1,0,3,2,1,0,5,4,3,2,1,0};

const int64_t POS[30] = {8, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3,
		     	 4, 4, 4, 4, 5, 5, 6, 6, 6, 6, 7, 7, 7,
			 7, 7, 7, 8};
const int64_t DIST[30] = {1,6,5,4,3,2,1,4,3,2,1,2,1,4,3,2,1,2,1,4,3,2,1,6,5,4,3,2,1,2};


const int64_t MOD = 30;
const int64_t GROUPS = 8;


uint64_t sieve_count(uint64_t* s, uint64_t prime_num, uint64_t first, uint64_t last, uint64_t L, uint64_t U, uint64_t* s2)
{
  /*use Legendre sieve for 2, 3, and 5 so sieve s does not have to store numbers
   divisible by 2, 3, or 5*/
  if      (prime_num == 1) return last - (first -1);
  else if (prime_num == 2) return last - last/2 - (first - 1 - (first -1)/2);
  else if (prime_num == 3) return last - last/2 - last/3 + last/6 -((first - 1) - (first - 1)/2 - (first - 1)/3 + (first - 1)/6);
  
  uint64_t a, b, c, ans;
  if(first > 0) a = (first  + NEXT[first % MOD])/MOD * GROUPS + POS[first % MOD];
  else a = 0;
  if (last > 0) b = (last   + PREV[last  % MOD])/MOD * GROUPS + POS[last % MOD];
  else b = 0;

  if (L > 1) c = (L-1 + PREV[(L-1) % MOD])/MOD * GROUPS + POS[(L-1) % MOD];
  else c = 0;


  if(a - c > 0 && b - c > 0) ans =  count_bits(s2, a - c, b - c);
  else ans = 0;  

  //  printf("L = %lu is bit %lu U is %lu\n", L, c, U);
  //  printf("sieve count args: prime_num = %lu first = %lu last = %lu L = %lu, U = %lu\n",prime_num,first,last,L,U);
  //printf("a = %lu b = %lu c = %lu\n",a,b,c);
  //  printf("start bit = %lu finish_bit = %lu\n", a-c, b-c);
  printf("L = %lu\ti = %lu\twheel = %lu\ttraditional = %lu\n",L, prime_num, ans, count_bits(s, first - L, last - L));

  return count_bits(s, first - L, last - L);
}

void sieve_step(uint64_t* s, uint64_t U, uint64_t L, uint64_t p, uint64_t* s2)
{
  uint64_t f, c;
  uint64_t first_f = L%p==0 ? L/p : L/p +1;
  for(f = first_f; p*f < U; f++)
    {
      assert(p*f-L>=0);
      clear_bit(s, p*f - L);
    }

  //version 2
  if(p > 5)
  {
    if (L > 1) c = (L-1 + PREV[(L-1) % MOD])/MOD * GROUPS + POS[(L-1) % MOD];
    else c = 0;
    assert(c >= 0);

    if (first_f > 1) first_f += DIST[first_f % MOD];
    assert(first_f > 0); assert(first_f == 1 || first_f >= 7);
    for(f = first_f; p*f < U; f += DIST[f % MOD])
    {
      //  printf("wheel sifting p = %lu f = %lu  p*f = %lu\n", p, f, p*f);
      clear_bit(s2, GROUPS*((p*f)/MOD) + POS[(p*f)%MOD] - c);      
    }
   }

  /*  if (L > 1) c = (L-1 + PREV[(L-1) % MOD])/MOD * GROUPS + POS[(L-1) % MOD];
  else c = 0;
  printf("L = %lu is bit %lu U is %lu p is %lu\n", L, c, U, p);
  if(p > 5)
  {  
     printf("first_f = %lu\n",first_f);
    if(first_f == 0) first_f = 1;
    if (first_f>1 && first_f < 7) first_f = 7;
    first_f += NEXT[first_f%MOD];
    printf("adjusted first_f = %lu\n",first_f);
    for(f=first_f; p*f<U; f+=DIST[f%MOD])
    {
      printf("wheel sifting %lu bit %lu\n", p*f, GROUPS*((p*f)/MOD) + POS[(p*f)%MOD] - c);
      clear_bit(s2, GROUPS*((p*f)/MOD) + POS[(p*f)%MOD] - c);
    }
    }*/
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

  int64_t* cum_mu = malloc((pi_crn+1)*sizeof(int64_t));
  memset(cum_mu, 0, (pi_crn+1)*sizeof(int64_t));

  uint64_t p, f, m, k, start, q, block_id, first_m, last_m, L, U, first_f;

  //Experimental mod 30 sieve
  uint64_t wheel_size = (crn  + PREV[crn % MOD])/MOD * GROUPS + POS[crn % MOD];
  uint64_t* s2 = malloc((wheel_size/64 + 2) * sizeof(uint64_t));

  for(block_id = start_block; block_id <= finish_block; block_id++)
    {
      memset(s, ~0, sizeof(uint64_t) * (crn/64+2));
      memset(s2, ~0, sizeof(uint64_t) * (wheel_size/64+2));
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
		  phi[i] += sieve_count(s, i, start, q, L, U, s2);
		  start = q+1;
		  if(mu[m] > 0) mpz_sub_ui(S2_result, S2_result, mu[m] * phi[i]);
		  else mpz_add_ui(S2_result, S2_result, -mu[m] * phi[i]);
		}
	    }
	  phi[i] += sieve_count(s, i, start, U - 1, L, U, s2);
	  sieve_step(s, U,  L, p, s2);
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
	  sieve_step(s, U,  L, p, s2);
	}
    }

  mpz_out_str (stdout, 10, S2_result);
  printf("\n");
  mpz_clears(buffer1, buffer2, n, S2_result, NULL);
  free(s); free(phi);
  free(primes); free(mu);
  free(lpf); free(cum_mu);
  free(s2);
  return(0);
}
