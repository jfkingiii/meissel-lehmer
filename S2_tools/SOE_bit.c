#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include "gmp.h"
#include <inttypes.h>

#define CLEAR(f,x)  f[x/64] &= ~(1ul<<(x%64))
#define TEST(f,x)    (f[x/64] & (1ul<<(x%64)))

void printbitssimple(uint64_t n) 
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
      printf("three arguments required\n");
      exit(1);
  }
  mpz_t n;
  mpz_init(n);

  mpz_ui_pow_ui(n,atol(argv[2]), atol(argv[3]));
  mpz_mul_ui(n, n, atol(argv[1]));
  uint64_t n2 = mpz_get_ui(n);

  uint64_t* s = malloc((n2/64+2)*sizeof(uint64_t));
  assert(s != NULL);
  uint64_t init = 0xAAAAAAAAAAAAAAAA;

  //for wheel mod 6
  //uint64_t word1 = 0x28A28A28A28A28A2;
  //uint64_t word2 = 0xA28A28A28A28A28A;
  //uint64_t word3 = 0x8A28A28A28A28A28;

  //printbitssimple(word1); printbitssimple(word2); printbitssimple(word3);
  memset(s, init, sizeof(uint64_t) * (n2/64+2));
  uint64_t sieve_to=floor(sqrt(n2));


  CLEAR(s,0);
  CLEAR(s,1);

  uint64_t i,j;
  for(i=3; i<=sieve_to; i+=2)
    if(TEST(s,i)!=0)
      for(j=i*i; j<=n2; j+= 2*i) CLEAR(s,j);

  //Need to add one to count since all even numbers, including 2, were zeroed in the bit array
  //Note ridiculous macro format required to print unit64_t without generating a compiler warning                                              
  printf("%" PRIu64 "\n", count_bits(s,0,n2) + 1);

  return(0);
}
