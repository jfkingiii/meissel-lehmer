#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "gmp.h"

#define FILE_NAME "S2_test_cases.txt"
#define NUM_LINES 67147

int main ()
{
    FILE* myfile = fopen(FILE_NAME,"r");
    assert(myfile != NULL);
    int a, b, c;
    uint64_t N, S2;
    int i;
    char str[30];
    for(i=0; i<NUM_LINES; i++)
    {
      fscanf(myfile,"%d %d %d %lu\n", &a, &b, &c, &S2);
      if( a*pow(b,c) > pow(10,14)) continue;
      printf("%d\t%d\t%d\t%lu\t",a, b, c, S2);
      fflush(stdout);
      sprintf(str, "./S2 %d %d %d", a, b, c);
      system(str);
      fflush(stdout);
    }
    fclose(myfile);
    return 0;
}
