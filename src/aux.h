
 #ifndef AUXILIARY
 #define AUXILIARY

 #include <stddef.h>

 #define ABS(x) \
         (((x) >= 0.0) ? (x) : -(x))

 #define NUM_EQU(x, y, pre) \
         (ABS((x) - (y)) <= (pre))

 #define SIGN(x) \
         (NUM_EQU((x), 0.0, 1.0E-012) ? 0.0 : (((x) > 0.0) ? 1.0 : -1.0))

 #define SIGN_N0(x) \
         (((x) >= 0) ? 1.0 : -1.0)

/* Find bin of constant spacing that contains a value: *************************

 findBin

 Find the bin number (interval number)

 This is an attempt towards a more efficient function to find the interval in 
 which a number lies (used to find the cell for a given direction of sky 
 radiance or the cell position in the spatially resolved integration), with 
 the downside that the intervals must be of equal width.

 The logic is that an integer interval can be found with the following 
 formulation:

 interval = floor(val - min / k),
 k = (max - min) / n_intervals,

 if val is in the range [min, max]. It is easy to accomodate the situation 
 when val can be outside that range.

  0 1 2 3 4 5 6 7 8  -> Bins (intervals)
 | | | | | | | | | |
 0 1 2 3 4 5 6 7 8 9 -> Break points
 
 Floor is used instead of ceiling since in C the index starts at 0. So in the 
 example above, for a value of 3.2, the bin index is 3 (= floor(3.2 - 0) / 1).
 This also results that all intervals are closed on the left.

 Since k is fixed for any vector of break points, it only needs to be 
 calculated once and is stored in the data structure for reuse.

 Note however that spatially resolved accumulators have as min and max 
 -INFINITY and INFINITY. In this case the values should logically be:

 min = vec[1] - res
 max = vec[n - 2] + res

 To avoid those repetitive calculations of  a pratical min and max and to make
 the function more generic, the needed min can be stored in the data structure
 while max is only necessary for the calculation of k.

 This is also a potential use for the phase function or cumulative phase 
 function sampling. Instead of using an interpolation function, tabulated data 
 can be used directly if at sufficient resolution. But note again that this
 function requires that the intervals have the same width. If a tabulated 
 phase function has unequal spacing, findInterv should be used instead or 
 interpolation.

 INPUT:
 val - Value to find interval for;
 min - (Practical) Minimum of the break points;
 k   - Scaling;
 n   - Number of bins;

 OUTPUT:
 The index of the bin (interval) that contains the value.

*******************************************************************************/
 
 static inline __attribute__((always_inline)) int
 findBin
 (
   double const val,
   double const min,
   double const k_inv,
   size_t const n
 )
 {
   int i = floor((val - min) * k_inv);

   if ( i < 0 )
   {
     i = 0;
   } else if ( i >= n ) {
     i = n - 1;
   }

   return i;
 }

 #define FINDBIN(I_val, I_min, I_k_inv, I_n, I_i) \
         I_i = floor((I_val - I_min) * I_k_inv); \
         if ( I_i < 0 ) \
         { \
           I_i = 0; \
         } else if ( I_i >= I_n ) { \
           I_i = I_n - 1; \
         }

 /* Function prototypes: ******************************************************/

 int findInterv (double val, 
                 double *p_brks, 
                 size_t n);
/*
 int findBin (double val, 
              double min,
              double k,
              size_t n);
*/

 int update_cdir (double psi, 
                  double phi, 
                  double *p_cdir);

 int chk_intsct (double *ppos, 
                 double *cpos, 
                 double *cdir, 
                 double s, 
                 double *ptfm,
                 double *sspos, 
                 double radius, 
                 double sshgt, 
                 int    virt);


 #endif // AUXILIARY

