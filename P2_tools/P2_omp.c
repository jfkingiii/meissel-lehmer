#define NUM_THREADS 16

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <omp.h>
#include "gmp.h"

uint64_t blocksize;

uint64_t sqr(uint64_t n) {return n*n;}

uint64_t cube(uint64_t n) {return n*n*n;}

uint64_t lpow(mpz_t n, uint64_t k, uint64_t l)
{
  //Returns floor(n^(k/l))
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
  assert(N > 1);
  uint64_t M = (N-1)/2 + 1;
  uint64_t i, f;

  /*contents of S[i] will be 1 if 2*i + 1 is prime, zero otherwise*/
  uint8_t* S = (uint8_t*)malloc(sizeof(uint8_t) * M);
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

uint64_t sieve_segment(uint64_t L, uint64_t U, const uint64_t* sift_primes, uint64_t I, uint8_t* s, uint64_t* found_primes)
{
  assert(L>=1);
  assert(L<=U);
  if(L==1) L++;

  /*
    L = lower bound of interval
    U = upper bound of interval
    sift_primes = array of primes to sieve with, with p[1] = 2, p[2] = 3, etc.
    I = length of p
    s = sieving array - reused for each segment
  */

  uint64_t i,p, start, pi=0;
  int64_t f;
  memset(s, 1, sizeof(uint8_t) * (U - L +1));

  /*s[i] is an indicator for L + i
    L,L+1,L+2,..., L + (U - L) -> s[0],s[1],s[2],s[3]...s[U]*/

  for(i=L%2; i<=U-L; i+=2) s[i]=0;
  if(L==2) s[0]=1;

  for(i=2; i<=I; i++)
  {
      p = sift_primes[i];
      if (p*p > U) break;

      /*Determine first factor to sift by
	start = smallest value of f so that f*p >= L*/

      start = L%p==0 ? L/p : L/p +1;
      assert(start*p>=L); assert((start-1)*p<L);
      if(start%2==0) start++;
      if (start<p) start=p;       //can start crossing out at p^2
      for(f=start*p; f<=U; f+= 2*p) s[f - L]=0;                                    
  }

  if (found_primes==NULL)
    for (i=0; i<= U - L; i++) pi += s[i];
  else
    for (i=0; i<= U - L; i++) if(s[i]==1) {pi++; found_primes[pi]=i+L;}
  return pi;
}


uint64_t sieve_segment2(uint64_t L, uint64_t U, const uint64_t* sift_primes, uint64_t I, uint8_t* s)
{
  assert(L>2);
  assert(L<=U);

  assert(L%2==1 && U%2==1);
  /*
    L = lower bound of interval
    U = upper bound of interval
    sift_primes = array of primes to sieve with, with p[1] = 2, p[2] = 3, etc.
    I = length of p
    s = sieving array - reused for each segment
  */

  uint64_t i,p, start, pi=0;
  int64_t f;
  memset(s, 1, sizeof(uint8_t) * ((U - L)/2+1));

  /*L and U are assumed to be odd
    L+2*0, L+2*1, L+2*2, ..., L+ 2*(U-L)/2 -> s[0], s[1], ...,s[(U - L)/2]
    s[i] is an indicator for L + 2*i*/

  for(i=2; i<=I; i++)
  {
      p = sift_primes[i];
      if (p*p > U) break;

      /*Determine first factor to sift by
	start = smallest value of f so that f*p >= L*/

      start = L%p==0 ? L/p : L/p +1;
      assert(start*p>=L); assert((start-1)*p<L);
      if(start%2==0) start++;
      if (start<p) start=p;       //can start crossing out at p^2
      for(f=start*p; f<=U; f+= 2*p) s[(f - L)/2]=0;
      }

  for (i=0; i<= (U - L)/2; i++) pi += s[i];
  return pi;
  }

#define SEGMENT_SIZE_KB 64
#define SEGMENT_SIZE (SEGMENT_SIZE_KB*1024)

uint64_t sieve_incrementally(uint64_t L, uint64_t U,  const uint64_t* sift_primes, uint64_t I)
{
  assert(L>=1);
  assert(L<=U);
  uint64_t sieve_to, sieve_from = L, pi = 0;
  uint8_t* s =(uint8_t*) malloc(SEGMENT_SIZE * sizeof(uint8_t));
  assert(s!=NULL);
  if(L==1) L++;
  sieve_from = L;
  sieve_to = L + SEGMENT_SIZE -1;
  if (sieve_to > U) sieve_to = U;

  for(;;)
  {
    pi += sieve_segment(sieve_from, sieve_to, sift_primes, I, s, NULL);
    if(sieve_to == U) break;
    sieve_from += SEGMENT_SIZE;
    sieve_to += SEGMENT_SIZE;
    if(sieve_to > U) sieve_to = U;
  }
  free(s);
  return pi;
}

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

  uint64_t frn = lpow(n, 1, 4);
  uint64_t crn = lpow(n, 1, 3);
  uint64_t srn = lpow(n, 1, 2);
  uint64_t ttrn = lpow(n, 2, 3);

  size_t psize = (size_t)((1.26*(crn+1)/log(crn+1) + 1)*sizeof(uint64_t));
  assert(psize < 1000000000);
  uint64_t* primes = (uint64_t*)malloc(psize);
  assert(primes!=NULL);
  uint64_t pi_crn = sieve(crn, primes);

  uint64_t  pi_frn;
  int64_t i;
  for(i=1; i<=pi_crn && primes[i]<=frn; i++); 
  pi_frn = i-1;

  uint64_t pi_srn =  sieve_incrementally(1, srn, primes, pi_frn);

  printf("pi(x^(1/4)) = %lu\n",pi_frn);
  printf("pi(x^(1/3)) = %lu\n", pi_crn);
  printf("pi(x^(1/2)) = %lu\n", pi_srn); fflush(stdout);

  //L= smallest integer >= sqrt(n)
  uint64_t L = srn;
  mpz_ui_pow_ui(buffer1, srn, 2);
  if(mpz_cmp(buffer1, n) < 0) L++;
  mpz_ui_pow_ui(buffer1, L, 2);
  mpz_ui_pow_ui(buffer2, L-1, 2);
  assert(mpz_cmp(buffer1, n)>=0 && mpz_cmp(buffer2, n)<0);

  uint64_t pi_L_1 = sieve_incrementally(1, L-1 , primes, pi_frn);

  //U = largest integer < n^(2/3)
  uint64_t U = ttrn;
  mpz_ui_pow_ui(buffer1, crn, 3);
  if(mpz_cmp(buffer1, n) == 0) U--;

  mpz_ui_pow_ui(buffer1, U, 3); mpz_ui_pow_ui(buffer2, U+1, 3);
  mpz_root(buffer1, buffer1, 2); mpz_root(buffer2, buffer2, 2);

  assert(mpz_cmp(buffer1, n)<0 && mpz_cmp(buffer2, n)>=0);

  blocksize = crn;
  uint64_t num_blocks = (U - L + 1)/blocksize;
  uint64_t short_block_length = (U - L + 1) % blocksize;
  assert(num_blocks * blocksize + short_block_length == U - L + 1);
  printf("L = %lu\nU = %lu\nnum_blocks = %lu\nblocksize = %lu\nshort_block_length = %lu\nlen1 = %lu\nlen2 = %lu\n",
	 L, U, num_blocks, blocksize, short_block_length, num_blocks * blocksize + short_block_length, U - L + 1);

  uint64_t sieve_from, sieve_to, start, finish, pi_seg2, pi_seg1, j, q;
  mpz_t P2, P2_total; 
  mpz_init(P2_total); mpz_set_ui(P2_total,0);
  uint8_t* s;
  uint64_t* pi;
  uint64_t* found_primes;
  uint64_t* pi_array = malloc((num_blocks+1)*sizeof(uint64_t));
  uint64_t* count_array = malloc((num_blocks+1)*sizeof(uint64_t));
  uint64_t blocks_complete=0;
  uint64_t first_odd, last_odd;

  //Begin parallel block
#pragma omp parallel num_threads(NUM_THREADS)  private(sieve_from, sieve_to, start, finish, pi_seg1, pi_seg2, j, q, P2, s, pi, found_primes, buffer1, buffer2, first_odd, last_odd)
  {
    if (omp_get_thread_num() == 0) printf("Running %d threads\n", omp_get_num_threads());

    mpz_init(buffer1); 
    mpz_init(buffer2);
    mpz_init(P2); 
    mpz_set_ui(P2,0);
    s = malloc(blocksize * sizeof(uint8_t));
    s = malloc(blocksize * sizeof(uint8_t));
    pi  = malloc(blocksize * sizeof(uint64_t));
    found_primes = malloc(blocksize * sizeof(uint64_t));

    #pragma omp for schedule(dynamic)
    for(i=1; i<=num_blocks; i++)
    {
      sieve_from = L + (i -1)*blocksize;
      sieve_to = sieve_from + blocksize - 1;
      mpz_tdiv_q_ui(buffer1, n, sieve_to);
      start = mpz_get_ui(buffer1);                                                    //start = n/sieve_to;
      mpz_tdiv_q_ui(buffer1, n, start);
      if(mpz_cmp_ui(buffer1, sieve_to) > 0 ) start++;                                 //if(n/start > sieve_to) start++;
      if (start<=crn) start++; 
      mpz_tdiv_q_ui(buffer1, n, start); mpz_tdiv_q_ui(buffer2, n, start - 1);
      assert(mpz_get_ui(buffer1)<=sieve_to); assert(mpz_get_ui(buffer2)>sieve_to);    //assert(n/start<=sieve_to); assert(n/(start-1)>sieve_to);
      mpz_tdiv_q_ui(buffer1, n, sieve_from);
      finish = mpz_get_ui(buffer1);                                                   //finish = n/sieve_from;
      mpz_tdiv_q_ui(buffer1, n, finish);                                              //assert(n/finish>=sieve_from); assert(n/(finish+1)<sieve_from);
      mpz_tdiv_q_ui(buffer2, n, finish + 1);   
      assert(mpz_get_ui(buffer1)>=sieve_from); assert(mpz_get_ui(buffer2)<sieve_from);
      assert(start <= finish); assert(sieve_from <= sieve_to);

      pi_seg1 = sieve_segment(start, finish, primes, pi_frn, s, found_primes);
      first_odd = sieve_from + 1 - sieve_from%2;
      last_odd = sieve_to - 1 + sieve_to%2;
      assert(last_odd>first_odd);
      pi_seg2 = sieve_segment2(first_odd ,last_odd, primes, pi_crn,  s);
      pi_array[i] = pi_seg2;
      count_array[i] = pi_seg1;

      if (sieve_from%2==0){pi[0]=0; pi[1]=s[0];}
      else {pi[0] = s[0]; pi[1] = s[0];}
      for(j=1; j<= (last_odd - first_odd)/2; j++)
      {
	pi[first_odd - sieve_from +2*j - 1] = pi[first_odd - sieve_from +2*j - 2] ;
	pi[first_odd - sieve_from +2*j] = pi[first_odd - sieve_from +2*j - 2] + s[j];
      }	
      if(sieve_to%2==0) pi[sieve_to - sieve_from]= pi[sieve_to - sieve_from -1];
      assert(pi[sieve_to - sieve_from]==pi_seg2);

      for(j=1;j<=pi_seg1;j++)
      {
	mpz_tdiv_q_ui(buffer1, n, found_primes[j]); 
	q = mpz_get_ui(buffer1);
	assert(q >= sieve_from && q/found_primes[j] <= sieve_to);
	mpz_add_ui(P2, P2, pi[q - sieve_from]);
      }
      blocks_complete++;
      if(blocks_complete % (num_blocks/1000) == 0) {printf("\rPercent complete %2.1f\n", ((blocks_complete*1000)/(num_blocks))/10.0); fflush(stdout);}

    }      //end parallel for
    #pragma omp critical
    {
     mpz_add(P2_total, P2_total, P2);
     free(s);  free(pi); free(found_primes);
     mpz_clear(buffer1); mpz_clear(buffer2);
    }
  }                //end parallel block


  pi_array[0] = pi_L_1;
  for(i=1; i<=num_blocks; i++) pi_array[i] += pi_array[i-1];
  for(i=1; i<=num_blocks; i++) mpz_add_ui(P2_total, P2_total, pi_array[i-1]*count_array[i]);

  s = (uint8_t*) malloc(blocksize * sizeof(uint8_t));
  pi  = (uint64_t*) malloc(blocksize * sizeof(uint64_t));
  found_primes = (uint64_t*) malloc(blocksize * sizeof(uint64_t));

  if(short_block_length>0)
  {
    printf("\nstarting short block\n"); fflush(stdout);
    //short block
    sieve_from = L + num_blocks * blocksize;
    sieve_to = sieve_from + short_block_length -1;
    mpz_tdiv_q_ui(buffer1, n, sieve_to);
    start = mpz_get_ui(buffer1);      
    if (start<=crn) start++;
    mpz_tdiv_q_ui(buffer1, n, start); mpz_tdiv_q_ui(buffer2, n, start - 1);
    assert(mpz_get_ui(buffer1)<=sieve_to); assert(mpz_get_ui(buffer2)>sieve_to);      //assert(n/start<=sieve_to); assert(n/(start-1)>sieve_to);
    mpz_tdiv_q_ui(buffer1, n, sieve_from);
    finish = mpz_get_ui(buffer1);                                                     //finish = n/sieve_from;
    if(start<=finish)
    {
      pi_seg1 = sieve_segment(start, finish, primes, pi_frn,  s, found_primes);
      pi_seg2 = sieve_segment(sieve_from, sieve_to, primes, pi_crn,  s, NULL);
      pi[0]=s[0];
      for(j=1; j<=sieve_to - sieve_from ; j++) pi[j]= pi[j-1] + s[j];
      for(j=1;j<=pi_seg1;j++)
      {
	mpz_tdiv_q_ui(buffer1, n, found_primes[j]);
	q = mpz_get_ui(buffer1);
	assert(q>=sieve_from && q<=sieve_to);
	mpz_add_ui(P2_total, P2_total,pi[q - sieve_from] + pi_L_1);
      }
    }
  }

  mpz_set_ui(buffer1, pi_crn); mpz_set_ui(buffer2, pi_srn);
  mpz_mul_ui(buffer1, buffer1, pi_crn-1); mpz_mul_ui(buffer2, buffer2, pi_srn-1);
  mpz_fdiv_q_2exp (buffer1, buffer1, 1); mpz_fdiv_q_2exp (buffer2, buffer2, 1);
  mpz_add(P2_total, P2_total, buffer1);
  mpz_sub(P2_total, P2_total, buffer2);
  printf("P2 = ");
  mpz_out_str(stdout, 10, P2_total);
  printf("\n");
  free(s); free(pi); free(found_primes);
  return 0;
}
