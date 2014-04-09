#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include "gmp.h"
#include <inttypes.h>

/*
  A simple sieve of Erastothenes using bit arrays and a mod 6 wheel.

  For machines with hardware popcnt instructions bit arrays are superior to arrays because much of the run time can be spent
  looping through the sieve array to count the primes.

  The mod 6 wheel saves time by avoiding sieving even numbers and numbers which are divisible by 3.

  s is the sieve array. The first bit of s[0] represents 0, the 64th bit represents 63, the 20th bit of s[1] represents 83 and so on.
  The bits in s are set to 1 if the associated number is prime, zero if it is composite
*/

#define CLEAR(f,x)  f[x/64] &= ~(1ul<<(x%64))			//Clear the bit associated with x from array f
#define TEST(f,x)    (f[x/64] & (1ul<<(x%64))) 			//Returns zero if the bit for x is not set, otherwise returns a positive number

#if 0
void printbitssimple(uint64_t n) 				//Prints the 64 bits of its argument in binary, useful for debugging
{
  uint64_t i;
  i = 1;

  while (i!=0) 
  {
    if (n & i)
      printf("1");
    else
      printf("0");
    i <<= 1;
    printf(" ");
  }
  printf("\n");
}
#endif

uint64_t count_bits(uint64_t *s, uint64_t first_bit, uint64_t last_bit)
{
  if(first_bit > last_bit) return 0;

  uint64_t first_block = first_bit/64;
  uint64_t last_block = last_bit/64;
  uint64_t i,j = 0;
  uint64_t mask1 = (1ul<<(last_bit%64+1))-1;
  if (last_bit%64 + 1 == 64) mask1 = ~0;
  uint64_t mask2 = ~((1ul<<first_bit%64)-1);

  if(first_block==last_block) return __builtin_popcountl(s[first_block] & mask1 & mask2);
  
  if(last_block - first_block > 1) 
    for(i=first_block+1; i<=last_block-1; i++) j += __builtin_popcountl(s[i]);

  j +=  __builtin_popcountl(mask1 & s[last_block]);
  j +=  __builtin_popcountl(mask2 & s[first_block]);
  return j;
}

int main(long argc, char *argv[])
{
  if(argc != 4)
  {
      printf("Three arguments required\n");
      printf("Usage: test6 a b c calculates pi(a*b^c)\n");
      exit(1);
  }

  //Legacy GMP code below; the input must fit in a 64 bit unsigned integer
  mpz_t n;
  mpz_init(n);

  mpz_ui_pow_ui(n,atol(argv[2]), atol(argv[3]));
  mpz_mul_ui(n, n, atol(argv[1]));
  uint64_t n2 = mpz_get_ui(n);

  uint64_t* s = malloc((n2/64+3)*sizeof(uint64_t));
  assert(s != NULL);


  //for wheel mod 6
  uint64_t word1 = 0x28A28A28A28A28A2;
  uint64_t word2 = 0xA28A28A28A28A28A;
  uint64_t word3 = 0x8A28A28A28A28A28;

  uint64_t k;

  for(k=0; k <= n2/64; k+=3)
  {
    s[k]=word1;
    s[k+1]=word2;
    s[k+2]=word3;
  }

  uint64_t sieve_to=floor(sqrt(n2));

  uint64_t i,j;

  i=5;								//5 is the first number which is relatively prime to 6
  int next = 4;             					//If i%6 = 5 then add 2 to i; if i%6 = 1 then add 4 to i;
  while (i<=sieve_to)
  {
    next = 6 - next;						//Alternate jumps of 2 and 4 to move through all the numbers congruent to +/- 1 mod 6
    if(TEST(s,i)!=0)						//Only sieve by primes
    {
      for(j=i*i; j<=n2; j+= 6*i)  CLEAR(s,j); 			//since i*i is always 1 (mod 6), this sieves out all the multiples of i congruent to 1 (mod 6)
      for(j=i*(i + next); j<=n2; j+= 6*i)  CLEAR(s,j); 		//i*(i+next) is the first multiple of i which exceeds i*i and is congruent to -1 (mod 6)
    }
    i += next;							//Advance pointer to the next number which is congruent to +/- 1
  }

  s[0] &= (~0 - 0xF);     					//Remove sieve entries for numbers less than 5;
  s[0] += 0xC;            					//Put in the correct entries; 2 and 3 are prime, 1100b = 4d + 8d = 12d = 0xC

 #if 0 
  for(i=0; i<=n2; i++)
    if(TEST(s,i)>0) printf("%" PRIu64 "\n",i);
#endif 

  //Note ridiculous macro format required to print unit64_t without generating a compiler warning                                              
  printf("%" PRIu64 "\n", count_bits(s,0,n2));

  return(0);
}
