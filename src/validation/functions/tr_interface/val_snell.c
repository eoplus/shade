/*
 gcc val_snell.c ../../../tr_interface.c -lm -lgsl -lgslcblas -Winline -O3 -o val_snell.o

 ./val_snell.o
*/


 #define COMPLEX_N

 #include <stdio.h>
 #include <math.h>
 #include <complex.h>
 #include "../../../constants.h"
 #include "../../../tr_interface.h"

 void main (void)
 {
  int    i;
  double CMPLX_T na = 1.000278 + CMPLX_0(0.000E-00 * I);
  double CMPLX_T nw = 1.342033 + CMPLX_0(2.442E-09 * I);
  double CMPLX_T nrat_a[3];
  double CMPLX_T nrat_w[3];
  nrat_a[0] = na / nw;
  nrat_a[1] = nrat_a[0] * nrat_a[0];
  nrat_a[2] = nrat_a[0] * nrat_a[0] * nrat_a[0];
  nrat_w[0] = nw / na;
  nrat_w[1] = nrat_w[0] * nrat_w[0];
  nrat_w[2] = nrat_w[0] * nrat_w[0] * nrat_w[0];

  double thetai, mui;
  double thetar;
  double step = M_PI_2 / 90.0;

  printf("Set 01 ==========================================================\n");
  printf("  Refractive indexes given by: \n");
  printf("    (na - 1): %.3E + %.3Ei (dry air at 550 nm)\n", CMPLX_F(creal)(na)-1, CMPLX_0(cimag(na)));
  printf("    (nw - 1): %.3E + %.3Ei (average seawater at 550 nm)\n\n", CMPLX_F(creal)(nw)-1, CMPLX_0(cimag(nw)));

  printf("  Test 1:------------------------------------------------------\n");
  printf("     Critical angle for total internal reflection: \n");
  printf("     Incidence in air   Incidence in water \n");
  printf("         %06.2f              %06.2f\n\n", asin(nrat_w[0]) * DEG, 
    asin(nrat_a[0]) * DEG);


  printf("  Test 2:------------------------------------------------------\n");
  printf("     Incidence in air   Incidence in water \n");
  printf("    __________________  __________________ \n");
  printf("    Incident Refracted  Incident Refracted \n");
  thetai = 0.0;
  for (i = 0; i < 91; i++)
  {
    mui = cos( thetai );
    thetar = acos( snell( mui, nrat_a ) );
    printf("      %06.2f    %06.2f", thetai * DEG, thetar*DEG);
    thetar = acos( snell( mui, nrat_w ) );
    printf("    %06.2f    %06.2f\n", thetai * DEG, thetar*DEG);
    thetai += step;
  }
 }


