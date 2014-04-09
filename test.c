#include <stdio.h>
#include <malloc.h>
#include <time.h>
#include <stdint.h>



#define TEST(f,x)(*(f+(x)/16)&(1<<(((x)%16L)/2)))
#define SET(f,x)*(f+(x)/16)|=1<<(((x)%16L)/2)

int
main(int argc, char *argv[])
{
  unsigned char *feld=NULL, *zzz;
  unsigned long teste=1, max, mom, hits=1, count, alloc, s=0, e=1;
  time_t begin;

  if (argc > 1)
    max = atol (argv[1]) + 10000;
  else
    max = 14010000L;

  printf("max=%lu\n", max);

  while (feld==NULL)
    zzz = feld = malloc (alloc=(((max-=10000L)>>4)+1L));

  for (count=0; count<alloc; count++) *zzz++ = 0x00;

  printf ("Searching prime numbers to : %ld\n", max);
  
  uint64_t loop_cnt=0;

  begin = time (NULL);
  while ((teste+=2) < max)
    if (!TEST(feld, teste)) {
      if  (++hits%2000L==0) {printf (" %ld. prime number\x0d", hits); fflush(stdout);}
      for (mom=3L*teste; mom<max; mom+=teste<<1) {SET (feld, mom); loop_cnt++;}
    }

  printf("loop count = %lu\n", loop_cnt);
  printf (" %ld prime numbers foundn %ld secs.\n\nShow prime numbers", 
          hits, time(NULL)-begin);
  /*
  while (s<e) {
    printf ("\n\nStart of Area : "); fflush (stdout); scanf ("%ld", &s);
    printf ("End   of Area : ");     fflush (stdout); scanf ("%ld", &e);

    count=s-2; if (s%2==0) count++;
    while ((count+=2)<e) if (!TEST(feld,count)) printf ("%ld\t", count);
    }*/
  free (feld);
  return 0;
}
