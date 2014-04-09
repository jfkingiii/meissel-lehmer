#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include "gmp.h"

#define ROW_CNT 159555
#define NUM_FILES 16
#define STRING_WIDTH 20
#define FILE_NAME "./s2_19/time.o%d"

int main ()
{
  uint64_t phi[ROW_CNT+1];
  int64_t mu[ROW_CNT+1];
  uint64_t cum_phi[ROW_CNT+1];
  memset(cum_phi, 0, (ROW_CNT+1)*sizeof(uint64_t));
  uint64_t tmp1;
  int64_t S2_chunk;
  int64_t S2_total = 0;
  int i,j;
  char infile[NUM_FILES+1][STRING_WIDTH];
  for(i=1; i<=NUM_FILES; i++)
    sprintf(infile[i],FILE_NAME,i);
  for(i=1; i<=NUM_FILES; i++)
  {
    printf("Reading file %s\n",infile[i]);
    FILE* myfile = fopen(infile[i],"r");
    assert(myfile!=NULL);

    assert(fscanf(myfile,"%ld\n", &S2_chunk)==1);
    S2_total += S2_chunk;
    for(j=1; j<=ROW_CNT;j++)
	assert(fscanf(myfile,"%lu %lu %ld", &tmp1, &phi[j], &mu[j])==3);
    if(i>1)
      for(j=1; j<=ROW_CNT; j++) S2_total += cum_phi[j]*mu[j];
    for(j=1; j<=ROW_CNT; j++) cum_phi[j] += phi[j];
    printf("File number %d read\n\n",i);fflush(stdout);
    fclose(myfile);
  }
  printf("\nS2 total = %ld\n", S2_total);
  return 0;
}
