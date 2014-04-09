#define WHEEL 4

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

uint64_t sqr(uint64_t y) {return(y*y);}

uint64_t sieve(uint64_t N, uint64_t* p)
{
  uint64_t M = N/2+N%2;
  uint8_t* S = malloc(sizeof(uint8_t) * M);
  if(S==NULL) {printf("malloc failed\n"); exit(1);}
  uint64_t i, j;
  uint64_t result = 0;
  uint64_t sieve_to = ceil((sqrt(N)-1)/2);

  /*contents of S[i] will be 1 if 2*i + 1 is prime, zero otherwise*/
  memset(S, 1, M*sizeof(uint8_t));
  S[0]=0;     /*1 is not prime*/

  for(i=1; i<=sieve_to; i++)
    if (S[i]==1)
     	  for(j=2*i*(i+1); j<M; j+=2*i+1) S[j]=0;
  p[1]=2;
  j=2;
  for(i=1; 2*i+1<=N; i++) 
  {
    if(S[i]==1)
    {
      result++;
      p[j]= 2*i+1; 
      j++;
    }
  }
  free(S);
  return(result+1); /*+1 to account for 2, the only even prime*/
}

#if WHEEL == 4
#define k 4
#define MOD 210
#define GROUPS 48
#define DIST {1,10,9,8,7,6,5,4,3,2,1,2,1,4,3,2,1,2,1,4,3,2,1,6,5,4,3,2,1,2,1,6,5,4,3,2,1,4,3,2,1,2,1,4,3,2,1,6,5,4,3,2,1,6,5,4,3,2,1,2,1,6,5,4,3,2,1,4,3,2,1,2,1,6,5,4,3,2,1,4,3,2,1,6,5,4,3,2,1,8,7,6,5,4,3,2,1,4,3,2,1,2,1,4,3,2,1,2,1,4,3,2,1,8,7,6,5,4,3,2,1,6,5,4,3,2,1,4,3,2,1,6,5,4,3,2,1,2,1,4,3,2,1,6,5,4,3,2,1,2,1,6,5,4,3,2,1,6,5,4,3,2,1,4,3,2,1,2,1,4,3,2,1,6,5,4,3,2,1,2,1,6,5,4,3,2,1,4,3,2,1,2,1,4,3,2,1,2,1,10,9,8,7,6,5,4,3,2,1,2}
#define POS {0,1,0,0,0,0,0,0,0,0,0,2,0,3,0,0,0,4,0,5,0,0,0,6,0,0,0,0,0,7,0,8,0,0,0,0,0,9,0,0,0,10,0,11,0,0,0,12,0,0,0,0,0,13,0,0,0,0,0,14,0,15,0,0,0,0,0,16,0,0,0,17,0,18,0,0,0,0,0,19,0,0,0,20,0,0,0,0,0,21,0,0,0,0,0,0,0,22,0,0,0,23,0,24,0,0,0,25,0,26,0,0,0,27,0,0,0,0,0,0,0,28,0,0,0,0,0,29,0,0,0,30,0,0,0,0,0,31,0,32,0,0,0,33,0,0,0,0,0,34,0,35,0,0,0,0,0,36,0,0,0,0,0,37,0,0,0,38,0,39,0,0,0,40,0,0,0,0,0,41,0,42,0,0,0,0,0,43,0,0,0,44,0,45,0,0,0,46,0,47,0,0,0,0,0,0,0,0,0,48}
#endif


#if WHEEL==3
#define k 3
#define MOD 30
#define GROUPS 8
#define DIST {1,6,5,4,3,2,1,4,3,2,1,2,1,4,3,2,1,2,1,4,3,2,1,6,5,4,3,2,1,2}
#define POS  {0,1,0,0,0,0,0,2,0,0,0,3,0,4,0,0,0,5,0,6,0,0,0,7,0,0,0,0,0,8}
#endif


int main(int argc, char *argv[])
{
  //Input
  uint64_t N;
  if(argc == 2)
    N = atol(argv[1]);
  else if(argc == 3)
    N = pow(atol(argv[1]), atol(argv[2]));
  else if (argc<2 || argc>3)
  {
      printf("Incorrect input format\n");
      exit(1);
  }
  else 
  {
      printf("Usage\n");
      exit(1);
  }

  //Small sieve
  uint64_t sqrt_N = sqrt(N);
  size_t psize = ((1.26*sqrt_N)/log(sqrt_N)+1)*sizeof(long);
  uint64_t* primes = malloc(psize);
  uint64_t J = sieve(sqrt_N ,primes);

  //Large sieve
  uint64_t i, j, f, p;
  uint8_t dist[MOD] = DIST;
  uint8_t pos[MOD] = POS;
  uint64_t array_size = GROUPS*(N/MOD) + pos[(N+dist[N%MOD])%MOD];
  uint8_t* s =  malloc(array_size*sizeof(uint8_t));
  memset(s, 1, array_size*sizeof(uint8_t));

  for(j=k+1; j<=J; j++)
  {
      p=primes[j];
      for(f=p; f<=N/p; f+=dist[f%MOD])
	    s[GROUPS*((p*f)/MOD) + pos[(p*f)%MOD]]=0;
  }

  uint64_t pi_N = k;
  for(i=2; i<array_size; i++) pi_N += s[i];
  printf("%lu\n", pi_N );
  free(s);
  free(primes);
  exit(0);
}
