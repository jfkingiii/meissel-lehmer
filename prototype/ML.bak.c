/*

    DON'T TOUCH THIS CODE - IT MATCHES THE TEST SUITE AND WILL BE
    USED FOR ALL FUTURE VALIDATION

    Prototype implementation of calculation of pi(x) by the Meissel-Lehmer
    method as described in Lagarias, Miller, and Odlyzko, Computing π(x): the
    Meissel-Lehmer Method. The prototype is designed to clarify the algorithm,
    not for high performance. This implementation keeps all numbers up to x^(2/3)
    in memory, so it will fail if x is too large, say > 10^15

    Usage: ML c b k to compute pi(c*b^k) or ML b k to compute pi(b^k).
    This program is stable and has been validated against the values of pi(x)
    from the web page of Tomás Oliveira e Silva.

    http://www.ieeta.pt/~tos/primes.html

    For best performance, compile with gcc ML.c -lm -lgmp -O3 -funroll-all-loops.

    Only 64 bit architectures are supported. The executable has been verified to run correctly 
    on Ubuntu and Red Hat. No Windows version is available but the code might compile and run
    under Cygwin (not verified). The GNU multiple precision arithmetic library (GMP) is required,
    see http://gmplib.org/.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#define OUTPUT_TYPE 0          //0 for full output, 1 to print pi(x) only, 2 to print pi(x) components only

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include "gmp.h"

uint64_t sqr (uint64_t n)
{
    return n * n;
}

uint64_t n_choose2 (uint64_t n)
{
    return n * (n - 1) / 2;
}

/*GMP is needed to compute floor(n^(2/3)) accurately*/
uint64_t lpow (uint64_t n, uint64_t k, uint64_t l)
{
    //Returns floor(n^(k/l))
    mpz_t x;
    mpz_init (x);
    mpz_set_ui (x, n);
    mpz_pow_ui (x, x, k);
    mpz_root (x, x, l);
    uint64_t  m = mpz_get_ui (x);
    mpz_clear(x);
    return m;
}


//Basic Sieve of Erastothenes sieving only odd numbers
//N is the number to be sieved to, the array p will be populated with primes up to N
uint64_t sieve (uint64_t N, uint64_t* p)    
{
    uint64_t M = N / 2 + N % 2;
    uint64_t i, j;
    uint64_t sieve_to = ceil ( (sqrt (N) - 1) / 2);
    assert (sqr (2 * sieve_to + 1) >= N);
    uint64_t pi = 1;

    /*contents of S[i] will be 1 if 2*i + 1 is prime, zero otherwise*/
    uint8_t* S = malloc (sizeof (uint8_t) * M);
    assert (S != NULL);
    memset (S, 1, sizeof (uint8_t) *M);
    S[0] = 0;   /*1 is not prime*/

    for (i = 1; i <= sieve_to; i++)
        if (S[i] == 1)
            for (j = 2 * i * (i + 1); j < M; j += 2 * i + 1) S[j] = 0;

    p[1] = 2;
    for (i = 1; 2 * i + 1 <= N; i++)
        if (S[i] == 1)
        {
            pi++;
            p[pi] = 2 * i + 1;
        }
    free (S);
    return (pi);
}

//For definitions of S1, S2 and P1 see Lagarias, Miller, and Odlyzko, Computing π(x): the Meissel-Lehmer Method
void S1 (uint64_t n, uint64_t crn, int8_t* mu, mpz_t result)
{
    mpz_t tmp;
    mpz_init (tmp);
    uint64_t i;
    for (i = 1; i <= crn; i++)
    {
        mpz_set_si (tmp, mu[i] * (n / i));
        mpz_add (result, result, tmp);
    }
    mpz_clear(tmp);
}

void P2 (uint64_t n, mpz_t result, uint64_t pi_crn, uint64_t crn, uint64_t srn, uint64_t ttn)
{
    size_t psize = (1.26 * ttn) / log (ttn) * sizeof (uint64_t) + 1;
    uint64_t* primes = malloc (psize);
    assert (primes != NULL);

    uint64_t J = sieve (ttn, primes);

    uint64_t i, pi_sqrtn;
    for(i=1; i <= J; i++) if(primes[i] > srn) break;
    pi_sqrtn = i - 1;

    uint64_t  j, k;
    i = 1;
    while (primes[i] <= crn) i++;                                      //primes[i] is the first prime > floor(n^(1/3))
    j = i;
    while (primes[j + 1] <= srn) j++;                                  //primes[j+1] is the first prime > floor(sqrt(n))
    assert (primes[i] > crn && primes[i - 1] <= crn);                  //so primes[j] is the largest prime <= sqrt(n)
    assert (primes[j] <= srn && primes[j + 1] > srn);

    uint64_t P2 = 0, pi = pi_sqrtn;
    for (k = j; k >= i; k--)
    {
      while (primes[pi] <= n / primes[k] && pi<=J ) pi++;              //Ugly, should be a cleaner way
      pi--;
      P2 += pi;
    }

    P2 += n_choose2 (pi_crn) - n_choose2 (pi_sqrtn);

    mpz_set_ui (result, P2);

    free (primes);
}

/*S2 is the most complex function, and uses most of the execution time*/
void S2 (uint64_t n, uint64_t crn, uint64_t ttn, uint64_t* primes, uint64_t J, int8_t* mu, uint64_t* lpf, mpz_t S2_result)
{
    uint8_t* s = malloc ( (ttn + 1) * sizeof (uint8_t));
    assert (s != NULL);
    memset (s, 1, sizeof (uint8_t) * (ttn + 1));
    uint64_t i, p, f, m, k, phi, start;
    int64_t answer;
    answer = 0;
    for (i = 0; i < J; i++)
    {
        p = primes[i + 1];
        start = 1;
        phi = 0;
        for (m = crn; m >= crn / p; m--)
        {
            if (mu[m] != 0 && p < lpf[m] && p * m > crn)
            {
                /*We've found a special leaf at p*m, so compute its contribution -mu(m)*phi(n/(p*m)) to S2*/
                if (i == 0) phi = n / (m * p);
                else
                {
                    for (k = start; k <= n / (m * p); k++) phi += s[k];
                    start = n / (m * p) + 1;
                }
                answer += -mu[m] * phi;
            }
        }
        /*Now sift the i+1th prime*/
        for (f = 1; p*f <= ttn; f++) s[p * f] = 0;
    }
    mpz_set_si (S2_result, answer);
    free (s);
}

int main (int argc, char *argv[])
{
    uint64_t multiplier, base, exponent, N;
    mpz_t tmp;
    mpz_init (tmp);

    /*Get the input*/
    if (argc > 4 || argc < 3)
    {
        printf ("Usage: ML c b k to compute pi(c*b^k) or ML b k to compute pi(b^k).\n");
        exit (1);
    }
    else if (argc == 3)
    {
        base = atol (argv[1]);
        exponent = atol (argv[2]);
        multiplier = 1;
    }
    else if (argc == 4)
    {
        multiplier = atol (argv[1]);
        base = atol (argv[2]);
        exponent = atol (argv[3]);
    }
    else
    {
      printf("Error reading input.\n");
      exit(1);
    }

    mpz_ui_pow_ui (tmp, base, exponent);
    N = multiplier * mpz_get_ui (tmp);
    if(N > pow(10, 15))
    {
      printf("Enter an argument less than or equal to 10^15.\nFor larger arguments use the high performance implementation.\n");
      exit(1);
    }
    mpz_clear(tmp);

    /*Deal with primes up to 11 for which the algorithm may not give the right answer*/
    uint8_t small[12] = {0, 0, 1, 2, 2, 3, 3, 4, 4, 4, 4, 5};
    if (N <= 11)
    {
        printf ("%d\n", small[N]);
        exit (0);
    }

    /*Processing starts here*/
    uint64_t crN = lpow (N, 1, 3);                                              //cube root of N
    uint64_t srN = lpow (N, 1, 2);                                              //square root of N
    uint64_t ttN = lpow (N, 2, 3);                                              //N to the two-thirds power

    size_t psize = (1.26 * crN) / log (crN) * sizeof (uint64_t) + 1;            //An upper bound on the number of primes <= crN
    uint64_t* primes = malloc (psize);
    assert (primes != NULL);

    uint64_t i, j;
    uint64_t J = sieve (crN , primes);                                          //We will need the primes up to crN. J = pi(crN).

    /*mu holds the value of the Möbius function*/
    int8_t* mu = malloc (sizeof (int8_t) * (crN + 1));
    assert (mu != NULL);
    memset (mu, 1, sizeof (int8_t) * (crN + 1));

    /*lpf holds the least prime factor of the associated integer*/
    uint64_t* lpf = malloc (sizeof (int64_t) * (crN + 1)); 
    assert (lpf != NULL);
    memset (lpf, 0, sizeof (int64_t) * (crN + 1));
    
    /*Both mu and lpf are needed to compute S2*/

    /*Compute mu and lpf up to crN*/
    for (i = 1; i <= J; i++)
        for (j = primes[i]; j <= crN; j += primes[i])
        {
            mu[j] = -mu[j];                                                     //Get the right sign for mu by flipping for each divisor
            if (lpf[j] == 0) lpf[j] = primes[i];
        }
    for (i = 1; i <= J; i++)
        for (j = sqr (primes[i]); j <= crN; j += sqr (primes[i])) mu[j] = 0;    //Set mu to zero for numbers that are not squarefree

    mpz_t S1_result, P2_result, S2_result, total;
    mpz_inits(S1_result, P2_result, S2_result, total, NULL);

    /*Calculate the three components of pi(x)*/
    S1 (N, crN, mu, S1_result);
    P2 (N, P2_result, J, crN, srN, ttN);
    S2 (N, crN, ttN, primes, J, mu, lpf, S2_result);

    /*pi(x) = pi(x^(1/3)) - P2 + S1 + S2 - 1 */
    mpz_set (total, S2_result);
    mpz_add (total, total, S1_result);
    mpz_sub (total, total, P2_result);
    mpz_sub_ui (total, total, 1);
    mpz_add_ui (total, total, J);

    /*OUTPUT_TYPE==0 is normal interactive output*/
    #if OUTPUT_TYPE==0
    printf ("\nS1 = ");
    mpz_out_str (stdout, 10, S1_result);
    printf ("\n\n");

    printf ("\nP2 = ");
    mpz_out_str (stdout, 10, P2_result);
    printf ("\n\n");

    printf ("\nS2 = ");
    mpz_out_str (stdout, 10, S2_result);
    printf ("\n\n");

    printf ("pi(x) = ");
    mpz_out_str (stdout, 10, total);
    printf("\n");

    /*OUTPUT_TYPE==1 is for automated validation of the calculated pi(x) against test cases*/
    #elif OUTPUT_TYPE==1
    mpz_out_str (stdout, 10, total);
    printf("\n");

    /*OUTPUT_TYPE==2 produces output to validate the high performance routines*/
    #elif OUTPUT_TYPE==2
    mpz_out_str (stdout, 10, S1_result); printf("\t");
    mpz_out_str (stdout, 10, P2_result); printf("\t");
    mpz_out_str (stdout, 10, S2_result); printf("\t");
    printf ("%lu\t", J);
    mpz_out_str (stdout, 10, total);
    printf ("\n");

    #else
    printf("Invalid output type\n");
    #endif

    mpz_clears(S1_result,P2_result, S2_result, total, NULL);
    free (primes);
    free (mu);
    free (lpf);
    return 0;
}
