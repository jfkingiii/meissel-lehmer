/* VERSION 1.2
This program requires the GNU GMP, MPFR, and MPC libraries.
Compile with gcc zeta.c -lmpc -lmpfr -lgmp
(the libraries must be referenced in this order or linking will fail)

11/10/2011. Added MPC library to calculate the complex functions
11/12/2011. Added a more precise estimate of the number of terms needed.
11/12/2011. Added progress indicator.
11/12/2011. Added a function to calculate Z(t)

Calculation of the Riemann Zeta function using the alternating series
approximation, very early draft, by James F. King

See:

H. Cohen, F. Rodriguez Villegas, D. Zagier, Convergence acceleration of
alternating series, Bonn, (1991)

Xavier Gourdon and Pascal Sebah, Numerical evaluation of the Riemann
Zeta-function, July 23, 2003, http://numbers.computation.free.fr/Constants/constants.html

Copyright 1995, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2005, 2009
Free Software Foundation, Inc.

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program.  If not, see http://www.gnu.org/licenses/.  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "gmp.h"
#include "mpfr.h"
#include "mpc.h"


#define PRECISION 512
#define ZETA_VARS d, d_inverse, b, c, x                            /*local mpfr variables for zeta calculation*/
#define PROGRESS 1                                                 /*Report the progress of the calculation?*/
#define PROGRESS_INTERVAL 1000                                     /*How often to report progress*/

double sqr(double y) {return y*y;}

void zeta(double sigma, double t, int decimals, mpc_t z_result)
{
  mpfr_t ZETA_VARS;
  mpc_t minus_s, z_accum, z, two, kplus1, z1;                      /*local mpc variables for zeta calculation*/

  mpfr_inits2(PRECISION, ZETA_VARS, NULL);                         /*Initialize mpfr variables*/
  
  mpc_init2(minus_s, PRECISION);  mpc_init2(z_accum, PRECISION);   /*Initialize mpc variables*/
  mpc_init2(z, PRECISION);  mpc_init2(two, PRECISION);  
  mpc_init2(kplus1, PRECISION); mpc_init2(z1, PRECISION); 

  mpc_set_ui(two, 2, MPC_RNDNN);                                   /*Set inital values*/
  mpc_set_d_d (minus_s, -sigma, -t, MPC_RNDNN);
  mpc_add_ui(z1, minus_s, 1, MPC_RNDNN);
  mpc_pow(z1, two, z1, MPC_RNDNN);    
  mpc_ui_sub(z1, 1, z1, MPC_RNDNN);                                /*z1 = 1 - 2^(1 - s)*/
                                                                   /*z1 is used in the calculation of both n and zeta*/
  mpc_abs(x,z1, MPFR_RNDN);
  double pole_term =  log(mpfr_get_d (x, MPFR_RNDN));              /*This is the term that gets large near s = 1*/
  unsigned long n = ceil(                                          /*Number of terms to use in expansion*/
			 (log(10)* decimals +
			 log(3) +
			 (fabs(t)*M_PI/2) +
			 log(1 + 2*fabs(t)) -
			 pole_term)/log(3 + sqrt(8))
						     );          
  
  #if PROGRESS > 0
  fprintf(stderr, "Computing %d terms...\n", n);
  #endif

  mpc_set_ui(z_accum, 0, MPC_RNDNN); 
  mpfr_sqrt_ui(d, 8,  MPFR_RNDN);
  mpfr_add_ui(d, d, 3,  MPFR_RNDN);
  mpfr_pow_ui(d, d, n, MPFR_RNDN);                                 /*d = (3 + sqrt(8))^n;*/
  mpfr_ui_div(d_inverse, 1, d, MPFR_RNDN);
  mpfr_add(d, d, d_inverse, MPFR_RNDN);
  mpfr_div_2exp(d, d, 1, MPFR_RNDN);                               /*d = (d + 1.0/d)/2.0;*/
  mpfr_set_si(b, -1, MPFR_RNDN);                                   /*b = -1*/
  mpfr_set(c, d, MPFR_RNDN);
  mpfr_neg(c, c, MPFR_RNDN);                                       /*c = -d*/
  long k;                                                          /*Loop variable for summing series*/

  for(k=0; k<n; k++){
    mpfr_sub(c, b, c, MPFR_RNDN);                                  /*c = b - c*/
    mpc_set_ui(kplus1, k+1, MPC_RNDNN);
    mpc_pow(z,  kplus1, minus_s, MPC_RNDNN);                       /*z = (k + 1)^(-s)*/
    mpc_mul_fr(z, z, c, MPC_RNDNN); 
    mpc_add(z_accum, z_accum, z, MPC_RNDNN);                       /*z_accum = z_accum + c*(k + 1)^(-s)*/
    mpfr_mul_si(b, b, (k + n)*(k - n),  MPFR_RNDN);
    mpfr_div_ui(b, b, (2*k + 1)*(k + 1),  MPFR_RNDN);
    mpfr_mul_2ui(b, b, 1, MPFR_RNDN);                              /*b = (k + n)*(k - n)*b/((k + 0.5)*(k+1))*/

    #if PROGRESS > 0
    if (k>0 && k % PROGRESS_INTERVAL == 0) {
      fprintf(stderr, "\r%d",k);
      fflush(stderr);
    }
    #endif

  }

  #if PROGRESS > 0
    fprintf(stderr, "\r%d\n\n",n);
  #endif

  mpc_div_fr(z_accum, z_accum, d, MPC_RNDNN);                     /*Divide result by d to complete alternating series calc*/
  mpc_div(z_result, z_accum, z1, MPC_RNDNN);                      /*z_result = value of zeta*/

  mpfr_clears(ZETA_VARS, NULL);                                   /*Clear all locals before returning*/
  mpc_clear(minus_s);  mpc_clear(z_accum); mpc_clear(z);  
  mpc_clear(two);  mpc_clear(kplus1); mpc_clear(z1);
}

#define DEFAULT_DECIMALS 10

mpfr_t pi, two_pi, pi_over_8;
unsigned int decimals;


void init_Z() {
  mpfr_init2(pi, PRECISION);
  mpfr_init2(two_pi, PRECISION);
  mpfr_init2(pi_over_8, PRECISION);
  mpfr_const_pi(pi, MPFR_RNDN);

  mpfr_const_pi(pi, MPFR_RNDN);
  mpfr_mul_2ui(two_pi, pi, 1, MPFR_RNDN);
  mpfr_div_2ui(pi_over_8, pi, 3, MPFR_RNDN);
}

void Z(mpfr_t Z_value, mpfr_t t,  mpc_t zeta_value)
{
    mpfr_t term, theta, zero;  
    mpc_t Z_complex;
    mpfr_init2(theta, PRECISION);
    mpfr_init2(term, PRECISION);
    mpfr_init2(zero, PRECISION);
    mpc_init2(Z_complex, PRECISION);


    mpfr_div(term, t, two_pi, MPFR_RNDN);
    mpfr_log(term, term, MPFR_RNDN);                /*log(t/(2*pi))*/
    mpfr_mul(term, t, term, MPFR_RNDN);
    mpfr_div_2ui(theta, term, 1,MPFR_RNDN);         /*(t/2) * log(t/(2*pi))*/
    
    mpfr_div_2ui(term,t,1,MPFR_RNDN);
    mpfr_sub(theta, theta, term, MPFR_RNDN);        /*(t/2) * log(t/(2*pi)) - t/2*/
    
    mpfr_sub(theta, theta, pi_over_8, MPFR_RNDN);
    
    mpfr_mul_ui(term, t, 48, MPFR_RNDN);
    mpfr_ui_div(term, 1, term, MPFR_RNDN);  

    mpfr_add(theta, theta, term, MPFR_RNDN);       /*theta = (t/2) * log(t/(2*pi)) - t/2 - pi/8* + 1/(48t)*/

    mpfr_set_ui(zero, 0, MPFR_RNDN);
    mpc_set_fr_fr (Z_complex, zero, theta, MPC_RNDNN);
    mpc_exp(Z_complex, Z_complex, MPC_RNDNN);
    mpc_mul(Z_complex, Z_complex, zeta_value, MPC_RNDNN);
    mpc_real(Z_value, Z_complex, MPC_RNDNN);
}

     
main(int argc, char *argv[])
{
  double sigma;
  double t;

  if(argc < 3){
    printf("Must specify two arguments, the real and imaginary part of s.\n");
    exit(1);
  }

  sigma = atof(argv[1]);
  t = atof(argv[2]);

  if (argc >= 4) decimals = atoi(argv[3]);
  else decimals = DEFAULT_DECIMALS;
  if (sigma==1 && t==0){
    fprintf(stderr, "WARNING: Attempt to evaluate zeta(1) in %s\n", argv[0]);
    printf("Inf");
    exit(1);
  }

  if(sqrt(sqr(sigma-1) + sqr(t)) < 0.1) 
    fprintf(stderr, "WARNING: The requested accuracy may not be acheived this close to the pole.\n");
  
  mpc_t answer;
  mpc_init2(answer, PRECISION);

  zeta(sigma, t, decimals, answer);
  mpc_out_str(stdout, 10, decimals, answer, MPC_RNDNN);


  /*  init_Z();

  mpfr_t val, T;
  mpfr_inits2(PRECISION,val, T, NULL);


  mpfr_set_d(T, t, MPFR_RNDN);
  Z(val, T,  answer);

  printf("\n\n");
  mpfr_out_str(stdout, 10, decimals, val, MPFR_RNDN);

  mpfr_clears(pi,two_pi, pi_over_8, NULL);          */
}
