#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>


void printbitssimple( uint64_t n )
{
  uint64_t i;
  i = 1ul << ( 64 - 1 );

  while( i > 0 )
    {
      if( n & i )
	printf( "1" );
      else
	printf( "0" );
      i >>= 1;
      printf( " " );
    }
  printf( "\n" );
}



void clear_bit( uint64_t* s, uint64_t j )
{
  /*                                                                                                                                                      Clear the jth bit of array s, first bit of the 0th                                                                                                    element of s is zero, bits are numbered 0, 1, ...63                                                                                                   Examples:                                                                                                                                          
    clear_bit(s, 0) clears the first bit of s[0]                                                                                                       
    clear_bit(s, 64) clears the first bit of s[1]                                                                                                       */

  uint64_t ind = j / 64;
  uint64_t bit = j % 64;
  s[ind] &= ~( 1ul << bit );
}


int main(long argc, char *argv[])
{
  long s_length = atol(argv[1]);
  long bit = atol(argv[2]);
  uint64_t* s = malloc(s_length * sizeof(uint64_t));
  memset(s, ~0, s_length* sizeof(uint64_t));
  clear_bit(s, bit);
  long i;
  for(i=63; i>=0; i--) printf("%lu ", i/10);
  printf("\n");
  for(i=63; i>=0; i--) printf("%lu ", i%10);
  printf("\n");
  for(i=0; i<s_length; i++)
    printbitssimple(s[i]);
  return 0;
}
