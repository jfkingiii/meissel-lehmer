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

  //uint64_t frn = lpow(n, 1, 4);
  uint64_t crn = lpow(n, 1, 3);
  //uint64_t srn = lpow(n, 1, 2);
  uint64_t ttrn = lpow(n, 2, 3);


  uint64_t blocksize = crn;
  uint64_t num_blocks = ttrn/blocksize;

  printf("\nnum_blocks = %lu",num_blocks);
  printf("\n\n");
  return(0);
}
