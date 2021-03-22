
/*
gcc val_fresnel.c ../../../tr_interface.c -lm -lgsl -lgslcblas -O3 -o val_fresnel.o
./val_fresnel.o
*/

 //#define COMPLEX_N

 #include <stdio.h>
 #include <stdlib.h>
 #include <math.h>
 #include <complex.h>
 #include "../../../constants.h"
 #include "../../../tr_interface.h"

 void main (void)
 {

  int    i;
  double CMPLX_T na = 1.000278 - CMPLX_0(0.000E-00 * I);
  double CMPLX_T nw = 1.342033 - CMPLX_0(2.442E-09 * I);
  double CMPLX_T nrat_a[3];
  double CMPLX_T nrat_w[3];
  nrat_a[0] = na / nw;
  nrat_a[1] = nrat_a[0] * nrat_a[0];
  nrat_a[2] = nrat_a[1] * nrat_a[0];
  nrat_w[0] = nw / na;
  nrat_w[1] = nrat_w[0] * nrat_w[0];
  nrat_w[2] = nrat_w[1] * nrat_w[0];

  double thetai, thetar;
  double mui;
  double CMPLX_T mur;
  double step = M_PI_2 / 360.0;
  struct Fm Fmat = {0.0};
  FILE   *fpo;

  printf("Set 01 ==========================================================\n");
  printf("  Refractive indexes given by: \n");
  printf("    (na - 1): %.3E + %.3Ei (dry air at 550 nm)\n", CMPLX_F(creal)(na)-1, CMPLX_0(cimag(na)));
  printf("    (nw - 1): %.3E + %.3Ei (average seawater at 550 nm)\n\n", CMPLX_F(creal)(nw)-1, CMPLX_0(cimag(nw)));

  printf("  Test 1:------------------------------------------------------\n");
  printf("    Incidence at the Brewster's angle. \n");
  printf("    theta_b = %06.2f\n", atan(nrat_w[0]) * DEG);
  printf("\n");
  mui = CMPLX_F(creal)(CMPLX_P(cos)(CMPLX_P(atan)(nrat_w[0])));
  mur = snell( mui, nrat_a );
  thetar = acos( CMPLX_F(creal)(mur) );
  fresnel (mui, mur, na, nw, nrat_w, &Fmat);
  //print_Fmat (&Fmat);
  printf("\n");

  fpo = fopen ("res/val_fresnel/air_water_550nm.txt", "w");
  printf("  Test 2:------------------------------------------------------\n");
  printf("                                                     Incidence in air \n");
  printf("    ____________________________________________________________________________________________________________________ \n");
  printf("    Incident Refracted        R(1,1)      R(1,2)      R(3,3)      R(3,4)      T(1,1)      T(1,2)      T(3,3)      T(3,4) \n");
  fprintf(fpo, "    Incident Refracted        R(1,1)      R(1,2)      R(3,3)      R(3,4)      T(1,1)      T(1,2)      T(3,3)      T(3,4) \n");
  thetai = 0.0;
  for (i = 0; i < 360; i++)
  {
    mui = cos( thetai );
    mur = snell( mui, nrat_a );
    thetar = acos( CMPLX_F(creal)(mur) );
    fresnel (mui, mur, na, nw, nrat_w, &Fmat);
    if( (i % 4) == 0.0)
    {
      printf("      %06.2f    %06.2f    %10.3E  %10.3E  %10.3E  %10.3E  %10.3E  %10.3E  %10.3E  %10.3E\n", 
        thetai * DEG, thetar * DEG, Fmat.R[0], Fmat.R[1], Fmat.R[2], Fmat.R[3], Fmat.T[0], Fmat.T[1], Fmat.T[2], Fmat.T[3]);
    }
    fprintf(fpo, "      %06.2f    %06.2f    %10.3E  %10.3E  %10.3E  %10.3E  %10.3E  %10.3E  %10.3E  %10.3E\n", 
      thetai * DEG, thetar * DEG, Fmat.R[0], Fmat.R[1], Fmat.R[2], Fmat.R[3], Fmat.T[0], Fmat.T[1], Fmat.T[2], Fmat.T[3]);
    thetai += step;
  }
  fclose(fpo);

  fpo = fopen ("res/val_fresnel/water_air_550nm.txt", "w");
  printf("\n");
  printf("  Test 3:------------------------------------------------------\n");
  printf("                                                    Incidence in water \n");
  printf("    ____________________________________________________________________________________________________________________ \n");
  printf("    Incident Refracted        R(1,1)      R(1,2)      R(3,3)      R(3,4)      T(1,1)      T(1,2)      T(3,3)      T(3,4) \n");
  fprintf(fpo, "    Incident Refracted        R(1,1)      R(1,2)      R(3,3)      R(3,4)      T(1,1)      T(1,2)      T(3,3)      T(3,4) \n");
  thetai = 0.0;
  for (i = 0; i < 360; i++)
  {
    mui = cos( thetai );
    mur = snell( mui, nrat_w );
    thetar = acos( CMPLX_F(creal)(mur) );
    fresnel (mui, mur, nw, na, nrat_a, &Fmat);
    if( (i % 4) == 0.0)
    {
      printf("      %06.2f    %06.2f    %10.3E  %10.3E  %10.3E  %10.3E  %10.3E  %10.3E  %10.3E  %10.3E\n", 
        thetai * DEG, thetar * DEG, Fmat.R[0], Fmat.R[1], Fmat.R[2], Fmat.R[3], Fmat.T[0], Fmat.T[1], Fmat.T[2], Fmat.T[3]);
    }
    fprintf(fpo, "      %06.2f    %06.2f    %10.3E  %10.3E  %10.3E  %10.3E  %10.3E  %10.3E  %10.3E  %10.3E\n", 
      thetai * DEG, thetar * DEG, Fmat.R[0], Fmat.R[1], Fmat.R[2], Fmat.R[3], Fmat.T[0], Fmat.T[1], Fmat.T[2], Fmat.T[3]);
    thetai += step;
  }
  fclose(fpo);

  printf("\nRunning R script val_fresnel.R:\n");
  system("Rscript val_fresnel.R");
  printf("Done!\n");
 }

