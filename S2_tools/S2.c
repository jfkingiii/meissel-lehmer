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


#if 0
void S2(uint64_t n, uint64_t crn, uint64_t* primes, uint64_t J, int8_t* mu, uint64_t* lpf, mpz_t S2_result)
{
  /*BLOCKING IDEA:

    Maintain an array phi(i), i=1, 2,... pi(N^(1/3)) to store phi(U[j],i) for each prime i. After each block B[j]
    is processed, increment the array phi(i) += number of elements of B[j] remaining after sieveing by the
    first i primes. Use the results to calculate phi for the special leaves in block B[j+1].

    For a given prime p, and block B[j] of numbers between L[j] and U[j]:
    Find the special leaves x associated with p[i] in B[j] in order from lowest to highest.
    Calculate phi(x1, i) buy summing the contents of the sieveing array s up to x1. s has already been
    sieved by p[1]...p[i].
    Add contribution of x1 to S2. phi(L[j]-1, i) has been stored in the previous step.
    Proceed to the next special leaf x2, and sum the elements of s from x1+1 to x2. Add phi(x1, i) + phi(L[j]-1, i)
    and iterate until all special leaves in B[j] have been accounted for.
    Now continue summing elements of s until we reach U[j] = L[j+1]-1. Store phi(U[j], i) in the array.
    Repeat with the next block.
  */
}
#endif


int main(int argc, char *argv[])
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

  //uint64_t frn = lpow(n, 1, 4);
  uint64_t crn = lpow(n, 1, 3);
  //uint64_t srn = lpow(n, 1, 2);
  uint64_t ttrn = lpow(n, 2, 3);

  size_t psize = (size_t)((1.26*(crn)/log(crn) + 1)*sizeof(uint64_t));
  assert(psize < 1000000000);
  uint64_t* primes = (uint64_t*)malloc(psize);
  assert(primes!=NULL);
  uint64_t pi_crn = sieve(crn, primes);

  int8_t* mu = malloc(sizeof(int8_t) * (crn+1));                        //mu will run from 1 to crn with s[0] not used
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
    for(j=sqr(primes[i]); j<=crn; j+=sqr(primes[i])) mu[j]=0;           //remove numbers that are not squarefree


  //  for(i=1; i<=crn; i++) printf("%lu\t%d\t%lu\n",i, mu[i], lpf[i]);


  mpz_t S2_result;
  mpz_init(S2_result);
  uint8_t* s = malloc((ttrn+1)*sizeof(uint8_t));
  assert(s != NULL);
  memset(s, 1, sizeof(uint8_t) * (ttrn+1));
  uint64_t p, f, m, k, phi, start, q;

  for(i=0; i<pi_crn; i++)
  {
    p = primes[i+1];
    start=1;
    phi = 0;
    for(m=crn; m>=crn/p; m--)
    {
      if(mu[m]!=0 && p<lpf[m] && m>crn/p)
      {
	/*We've found a special leaf at p*m, so compute its contribution -mu(m)*phi(n/(p*m)) to S2*/
	printf("%lu\t%lu\n", p, m);
	mpz_tdiv_q_ui(buffer1, n, p*m);
	q = mpz_get_ui(buffer1);
	switch(i)
	{
	   case 0: phi = q; break;
	   case 1: phi = q - q/2; break;
	   case 2: phi = q - q/2 - q/3 + q/6; break;
	   case 3: phi = q - q/2 - q/3 -q/5 + q/6 + q/10 + q/15 - q/30;  break;
	   case 4: phi = q - q/2 - q/3 -q/5 - q/7 + q/6 + q/10 + q/15 + q/14 + q/21 +q/35 - q/30 - q/105 - q/70 -q/42 + q/210;  break;
	   default:
	   {
	     for(k=start; k<=q; k++) phi += s[k];
	     start = q+1;
	   }
	}
	if(mu[m] > 0)  mpz_sub_ui(S2_result, S2_result, mu[m]*phi);
	else mpz_add_ui(S2_result, S2_result, -mu[m]*phi);
      }
    }
    /*Now sift the i+1th prime*/
    if (p>2) for(f=1; p*f <= ttrn; f+=2) s[p*f]=0;
    else for(f=1; p*f<=ttrn; f++)  s[p*f]=0;
  }
  printf("S2 = ");
  mpz_out_str (stdout, 10, S2_result);
  printf("\n");
  free(s);
  free(primes); free(mu); free(lpf);
  return 0;
}
