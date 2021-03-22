
#include <stdio.h>
#include <stdlib.h>	// FUnctions: calloc
#include <math.h>
#include "aux.h"


/* Function to find the interval index where a value lies. It is open at the 
 * left and closed in the right. brks extremes can be be -INFINITE and INFINITE.
 *
 * However, is not very efficient... The function finBin below adress this 
 * problem, requiring that the minimum of the breaks and "k" 
 * (= n_bins / (max - min)) be passed to the function. Min is fixed but not 
 * directly available as would always be -INFINITY for spatial dimension, so it
 *
 *
 * val  = scalar value
 * brks = vector of breaks
 * n    = length of brks
 */

 int findInterv (double val, 
                 double *p_brks, 
                 size_t n)
 {
   int i;

   if(val <= p_brks[1])
   {
     return 0;
   } else if (val > p_brks[n - 2]) {
     return n - 2;
   }

   for(i = 1; i < (n - 2); i++)
   {
     if(((val - p_brks[i+1]) * (val - p_brks[i])) <= 0)
     {
       return i;
     }
   }
 }

/*
 * Find the bin number (interval number)
 *
 * This is an attempt towards a more efficient function to find the interval in 
 * which a number lies (used to find the cell for a given direction of sky 
 * radiance or the cell position in the spatially resolved integration), with 
 * the downside that the intervals must be of equal width.
 *
 * The logic is that an integer interval can be found with the following 
 * formulation:
 *
 * interval = floor(val - min / k),
 * k = (max - min) / n_intervals,
 *
 * if val is in the range [min, max]. It is easy to accomodate the situation 
 * when val can be outside that range.
 *
 *  0 1 2 3 4 5 6 7 8  -> Bins (intervals)
 * | | | | | | | | | |
 * 0 1 2 3 4 5 6 7 8 9 -> Break points
 * 
 * Floor is used instead of ceiling since in C the index starts at 0. So in the 
 * example above, for a value of 3.2, the bin index is 3 (= floor(3.2 - 0) / 1).
 * This also results that all intervals are closed on the left.
 *
 * Since k is fixed for any vector of break points, it only needs to be 
 * calculated once and is stored in the data structure for reuse.
 *
 * Note however that spatially resolved accumulators have as min and max 
 * -INFINITY and INFINITY. In this case the values should logically be:
 *
 * min = vec[1] - res
 * max = vec[n - 2] + res
 *
 * To avoid those repetitive calculations of  a pratical min and max and to make
 * the function more generic, the needed min can be stored in the data structure
 * while max is only necessary for the calculation of k.
 *
 * This is also a potential use for the phase function or cumulative phase 
 * function sampling. Instead of using an interpolation function, tabulated data 
 * can be used directly if at sufficient resolution. But note again that this
 * function requires that the intervals have the same width. If a tabulated 
 * phase function has unequal spacing, findInterv should be used instead or 
 * interpolation.
 *
 */
/* 
 int findBin (double val,	// Value to find interval for
              double min,	// (Practical) Minimum of the break points
              double k_inv,		// Scaling
              size_t n)		// Number of bins
 {
   int i;

   i = floor((val - min) * k_inv);

   if (i < 0)
   {
     i = 0;
   } else if (i >= n) {
     i = n - 1;
   }

   return i;
 }
*/
 int update_cdir (double psi, 
                  double phi, 
                  double *p_cdir) 
 {
   double mus = cos(psi);
   double sint = sin(acos(p_cdir[2]));
   double sins = sin(psi);
   double scatmb[3] = {0.0};
   double scatma[3][3] = {0.0};

   scatmb[0] = sins * cos(phi);
   scatmb[1] = sins * sin(phi);
   scatmb[2] = mus;

   if(ABS(sint) > 1.0E-012)
   {
     scatma[0][0] =  p_cdir[0] * p_cdir[2] / sint; 
     scatma[0][1] = -p_cdir[1] / sint; 
     scatma[0][2] =  p_cdir[0];
     scatma[1][0] =  p_cdir[1] * p_cdir[2] / sint;
     scatma[1][1] =  p_cdir[0] / sint;
     scatma[1][2] =  p_cdir[1];
     scatma[2][0] = -sint;
     scatma[2][1] =  0.0;
     scatma[2][2] =  p_cdir[2]; 
     p_cdir[0] = (scatma[0][0] * scatmb[0]) + 
                 (scatma[0][1] * scatmb[1]) + 
                 (scatma[0][2] * scatmb[2]); 
     p_cdir[1] = (scatma[1][0] * scatmb[0]) + 
                 (scatma[1][1] * scatmb[1]) + 
                 (scatma[1][2] * scatmb[2]);
     p_cdir[2] = (scatma[2][0] * scatmb[0]) + 
                 (scatma[2][1] * scatmb[1]) + 
                 (scatma[2][2] * scatmb[2]); 
   } else {
     p_cdir[0] = SIGN(p_cdir[2]) * scatmb[0];
     p_cdir[1] = SIGN(p_cdir[2]) * scatmb[1];
     p_cdir[2] = SIGN(p_cdir[2]) * scatmb[2];
   }
   return 0;
 }


 int chk_intsct (double *ppos,
                 double *cpos_i, 
                 double *cdir, 
                 double s, 
                 double *ptfm, 
                 double *sspos, 
                 double radius, 
                 double sshgt, 
                 int    virt)
 {

   double delta = 0.0; 
   double z1 = 0.0, z2 = 0.0;
   double x1 = 0.0, x2 = 0.0;
   double dx1 = 0.0, dx2 = 0.0;
   double dx = 0.0, dy = 0.0, dz = 0.0;
   double dr2 = 0.0, D = 0.0;
   double root = 0.0;
   double cpos[3] = {0.0};

   int res = 0;


   // Propagate back to surface if above. This is a way to garantee that function 
   // will capture photon exiting from inside the cilinder. 
   if(cpos_i[2] < 0) 
   {
     s = cpos_i[2] / cdir[2];
     cpos[0] = cpos_i[0] + (-cdir[0] * s);
     cpos[1] = cpos_i[1] + (-cdir[1] * s);
     cpos[2] = 0.0; // to avoid rounding errors
   } else {
     cpos[0] = cpos_i[0];
     cpos[1] = cpos_i[1];
     cpos[2] = cpos_i[2];
   }
 
   // Check platform intersection:
   if(cpos[0] >= ptfm[0] && cpos[0] <= ptfm[1] &&
      cpos[1] >= ptfm[2] && cpos[1] <= ptfm[3] &&
      cpos[2] < 1E-12) res = 1;

   if ( sqrt(cpos[0] * cpos[0] + cpos[1] * cpos[1]) < radius && 
       (cpos[2] < sspos[2] && cpos[2] >= (sspos[2] - sshgt)) && 
       radius > 0.0) res = 2;

   dx  = cpos[0] - ppos[0];
   dy  = cpos[1] - ppos[1];
   dz  = cpos[2] - ppos[2];
   dr2 = (dx * dx) + (dy * dy);
   D   = (ppos[0] * cpos[1]) - (cpos[0] * ppos[1]);
   delta = ((radius * radius) * dr2) - (D * D);
   if (delta >= 0)
   {
     root = SIGN(dy) * dx * sqrt(delta);
     x1 = ((D * dy) - root) / dr2;
     x2 = ((D * dy) + root) / dr2;
     dx1 = (x1 - ppos[0]) / dx;
     dx2 = (x2 - ppos[0]) / dx;
     z1 = ppos[2] + dx1 * dz;
     z2 = ppos[2] + dx2 * dz;
     if ( (dx1 >= 0.0 && dx1 <= 1.0 & z1 <= sspos[2] & z1 >= sshgt) || 
          (dx2 >= 0.0 & dx2 <= 1.0 & z2 <= sspos[2] & z2 >= sshgt)) res = 3;
   }

   return res;
 }


