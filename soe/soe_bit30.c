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

#define CLEAR(f, x)  f[x/64] &= ~(1ul<<(x%64))	  //Clear the bit associated with x from array f
#define SET(f, x)    f[x/64] |= (1ul<<(x%64))     //Set the bit associated with x from array f
#define TEST(f, x)    (f[x/64] & (1ul<<(x%64)))   //Returns zero if the bit for x is not set, otherwise returns a positive number


#if 0
void printbitssimple(uint64_t n) 		  //Prints the 64 bits of its argument in binary, useful for debugging
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

int get_next(uint64_t x) 
{
    uint64_t  y;
    switch (x%30) 
    {
    case 1:
	y=6; 
	break; 
    case 7:
	y=4; 
	break; 
    case 11:
	y=2; 
	break; 
    case 13:
	y=4; 
	break; 
    case 17:
	y=2; 
	break; 
    case 19:
	y=4; 
	break; 
    case 23:
	y=6; 
	break; 
    case 29:
	y=2; 
	break;
    default:
	printf("Error: could not find residue %lu\n", x%30);
	exit(1);
	break; 
    }
    return y;
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

void init(uint64_t* pattern)
{
    int rel_prime[8] = {1, 7, 11, 13, 17, 19, 23, 29};
    int bit_list[30];
    int i;
    for(i=0; i<30; i++) bit_list[i] = 0;
    for(i=0; i<8; i++) bit_list[rel_prime[i]] = 1;
    for(i=0; i<15; i++) pattern[i] = 0;
    for(i=0; i<64 * 15; i++) if(bit_list[i%30] == 1) SET(pattern, i);
}

int main(int argc, char *argv[])
{
  if(argc != 4)
  {
      printf("Three arguments required\n");
      printf("Usage: %s a b c calculates pi(a*b^c)\n", argv[0]);
      exit(1);
  }

  //Legacy GMP code below; the input must fit in a 64 bit unsigned integer
  mpz_t n;
  mpz_init(n);

  mpz_ui_pow_ui(n, atol(argv[2]), atol(argv[3]));
  mpz_mul_ui(n, n, atol(argv[1]));
  uint64_t n2 = mpz_get_ui(n);
  mpz_clear(n);

  uint64_t* s = calloc(n2/64+3, sizeof(uint64_t));
  assert(s != NULL);

  //for wheel mod 30
  uint64_t wheel_pattern[15];
  init(wheel_pattern);
  uint64_t k;
  for(k=0; k <= n2/64; k++) s[k] = wheel_pattern[k%15];

  uint64_t sieve_to = floor(sqrt(n2));

  uint64_t i,j;


  i=7;							//7 is the first number > 1 which is relatively prime to 30
  while (i <= sieve_to)
  {
    if(TEST(s, i) != 0)					//Only sieve by primes
    {
	switch(i%30)
	{
	case 1:
	    for(j=i*i; j<=n2; j+= 30*i)  CLEAR(s, j);        
	    for(j=i*(i+6); j<=n2; j+= 30*i)  CLEAR(s, j);    
	    for(j=i*(i+10); j<=n2; j+= 30*i)  CLEAR(s, j);   
	    for(j=i*(i+12); j<=n2; j+= 30*i)  CLEAR(s, j);   
	    for(j=i*(i+16); j<=n2; j+= 30*i)  CLEAR(s, j);   
	    for(j=i*(i+18); j<=n2; j+= 30*i)  CLEAR(s, j);  
	    for(j=i*(i+22); j<=n2; j+= 30*i)  CLEAR(s, j);  
	    for(j=i*(i+28); j<=n2; j+= 30*i)  CLEAR(s, j);  
        break;

	case 7:
	    for(j=i*i; j<=n2; j+= 30*i)  CLEAR(s, j);    
	    for(j=i*(i+4); j<=n2; j+= 30*i)  CLEAR(s, j);   
	    for(j=i*(i+6); j<=n2; j+= 30*i)  CLEAR(s, j);   
	    for(j=i*(i+10); j<=n2; j+= 30*i)  CLEAR(s, j);    
	    for(j=i*(i+12); j<=n2; j+= 30*i)  CLEAR(s, j);  
	    for(j=i*(i+16); j<=n2; j+= 30*i)  CLEAR(s, j);  
	    for(j=i*(i+22); j<=n2; j+= 30*i)  CLEAR(s, j);  
	    for(j=i*(i+24); j<=n2; j+= 30*i)  CLEAR(s, j);  
        break;

	case 11:
	    for(j=i*(i+0); j<=n2; j+= 30*i)  CLEAR(s, j);    
	    for(j=i*(i+2); j<=n2; j+= 30*i)  CLEAR(s, j);   
	    for(j=i*(i+6); j<=n2; j+= 30*i)  CLEAR(s, j);   
	    for(j=i*(i+8); j<=n2; j+= 30*i)  CLEAR(s, j);    
	    for(j=i*(i+12); j<=n2; j+= 30*i)  CLEAR(s, j);  
	    for(j=i*(i+18); j<=n2; j+= 30*i)  CLEAR(s, j);  
	    for(j=i*(i+20); j<=n2; j+= 30*i)  CLEAR(s, j);  
	    for(j=i*(i+26); j<=n2; j+= 30*i)  CLEAR(s, j);  
        break;

	case 13:
	    for(j=i*(i+0); j<=n2; j+= 30*i)  CLEAR(s, j);    
	    for(j=i*(i+4); j<=n2; j+= 30*i)  CLEAR(s, j);   
	    for(j=i*(i+6); j<=n2; j+= 30*i)  CLEAR(s, j);   
	    for(j=i*(i+10); j<=n2; j+= 30*i)  CLEAR(s, j);    
	    for(j=i*(i+16); j<=n2; j+= 30*i)  CLEAR(s, j);  
	    for(j=i*(i+18); j<=n2; j+= 30*i)  CLEAR(s, j);  
	    for(j=i*(i+24); j<=n2; j+= 30*i)  CLEAR(s, j);  
	    for(j=i*(i+28); j<=n2; j+= 30*i)  CLEAR(s, j);  
        break;

	case 17:
	    for(j=i*(i+0); j<=n2; j+= 30*i)  CLEAR(s, j);    
	    for(j=i*(i+2); j<=n2; j+= 30*i)  CLEAR(s, j);   
	    for(j=i*(i+6); j<=n2; j+= 30*i)  CLEAR(s, j);   
	    for(j=i*(i+12); j<=n2; j+= 30*i)  CLEAR(s, j);    
	    for(j=i*(i+14); j<=n2; j+= 30*i)  CLEAR(s, j);  
	    for(j=i*(i+20); j<=n2; j+= 30*i)  CLEAR(s, j);  
	    for(j=i*(i+24); j<=n2; j+= 30*i)  CLEAR(s, j);  
	    for(j=i*(i+26); j<=n2; j+= 30*i)  CLEAR(s, j);  
        break;

	case 19:
	    for(j=i*(i+0); j<=n2; j+= 30*i)  CLEAR(s, j);    
	    for(j=i*(i+4); j<=n2; j+= 30*i)  CLEAR(s, j);   
	    for(j=i*(i+10); j<=n2; j+= 30*i)  CLEAR(s, j);   
	    for(j=i*(i+12); j<=n2; j+= 30*i)  CLEAR(s, j);    
	    for(j=i*(i+18); j<=n2; j+= 30*i)  CLEAR(s, j);  
	    for(j=i*(i+22); j<=n2; j+= 30*i)  CLEAR(s, j);  
	    for(j=i*(i+24); j<=n2; j+= 30*i)  CLEAR(s, j);  
	    for(j=i*(i+28); j<=n2; j+= 30*i)  CLEAR(s, j);  
        break;

	case 23:
	    for(j=i*(i+0); j<=n2; j+= 30*i)  CLEAR(s, j);    
	    for(j=i*(i+6); j<=n2; j+= 30*i)  CLEAR(s, j);   
	    for(j=i*(i+8); j<=n2; j+= 30*i)  CLEAR(s, j);  
	    for(j=i*(i+14); j<=n2; j+= 30*i)  CLEAR(s, j);   
	    for(j=i*(i+18); j<=n2; j+= 30*i)  CLEAR(s, j);    
	    for(j=i*(i+20); j<=n2; j+= 30*i)  CLEAR(s, j);  
	    for(j=i*(i+24); j<=n2; j+= 30*i)  CLEAR(s, j);  
	    for(j=i*(i+26); j<=n2; j+= 30*i)  CLEAR(s, j);  
        break;

	case 29:
            for(j=i*(i+0); j<=n2; j+= 30*i)  CLEAR(s, j);
            for(j=i*(i+2); j<=n2; j+= 30*i)  CLEAR(s, j);
            for(j=i*(i+8); j<=n2; j+= 30*i)  CLEAR(s, j);
            for(j=i*(i+12); j<=n2; j+= 30*i)  CLEAR(s, j);
            for(j=i*(i+14); j<=n2; j+= 30*i)  CLEAR(s, j);
            for(j=i*(i+18); j<=n2; j+= 30*i)  CLEAR(s, j);
            for(j=i*(i+20); j<=n2; j+= 30*i)  CLEAR(s, j);
            for(j=i*(i+24); j<=n2; j+= 30*i)  CLEAR(s, j);
	break;

	default:
	    printf("Error\n");
	    exit(1);

	}
    }
    i += get_next(i);						//Advance pointer to the next number
  }



  #if 0
  for(i=0; i<=n2; i++)
    if(TEST(s,i)>0) printf("%" PRIu64 "\n",i);
  #endif 

  //Note ridiculous macro format required to print unit64_t without generating a compiler warning                                              
  printf("%" PRIu64 "\n", count_bits(s, 0, n2)+2);
  free(s);

  return(0);
}


