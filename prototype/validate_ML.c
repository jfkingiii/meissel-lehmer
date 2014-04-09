/*
  validate_ML.c: Reads and prints records from the test case file, runs the prototype ML calculation,
  and prints the results for comparison to make sure ML is working correctly.

  Once ML is validated it is used to produce values of P2, S1, and S2 to validate the high performance programs
  validate_ML must be launched from the same directory that contains the test case files and the program ML_test
*/
 
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "gmp.h"

#define FILE_NAME "test_cases.txt"
#define NUM_LINES 161473

int main ()
{
  FILE* myfile = fopen(FILE_NAME,"r");
  assert(myfile != NULL);
  int a, b, c;
  uint64_t pi;
  uint64_t i;
  char str[30];
  for(i=0; i<NUM_LINES; i++)
    {
      /*Read one record from test case file*/
      if(fscanf(myfile,"%d %d %d %lu\n", &a, &b, &c, &pi) !=4 )
	printf("Failed to read file!\n");
      if( a*pow(b,c) > pow(10,11)) continue;
      
      /*Print record*/
      printf("%d\t%d\t%d\t%lu\t",a, b, c, pi);
      fflush(stdout);
      
      /*Execute ML program and print calculated value of pi(x)*/
      sprintf(str, "./ML_test %d %d %d", a, b, c);
      if(system(str) !=0 )
	printf("Program execution failed!\n");
      fflush(stdout);
    }
  fclose(myfile);
  return 0;
}

