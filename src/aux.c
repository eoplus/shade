
/*******************************************************************************
 aux.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 Version: 2021-01-30

 Auxiliary functions. findInterv is defined here, finBin is declared in the
 aux.h header file.

*******************************************************************************/

 #include <stddef.h>	// size_t
 #include <math.h>	// floor

 #include "aux.h"

/* Find interval that contain a value: *****************************************

 findInterv

 Function to find the interval index where a value lies. It is open at the left
 and closed in the right. Breaks extremes can be be -INFINITE and INFINITE.

 However, is not very efficient... The function finBin (aux.h header file) below
 adress this problem, requiring that the minimum of the breaks and "k"
 (= n_bins / (max - min)) be passed to the function. Min is fixed but not 
 directly available as would always be -INFINITY for spatial dimension, so it
 
 INPUT: 
 val  - Scalar value
        Constant double;
 brks - Vector of breaks
        Pointer to constant double;
 n    - Number of break points
        Constant size_t.

 OUTPUT:
 The index of the interval in which the value lies.

*******************************************************************************/

 int
 findInterv
 (
   double const val,
   double const * brks,
   size_t const n
 )
 {
   if ( val <= brks[1] )
   {
     return 0;
   }
   else if ( val > brks[n - 2] )
   {
     return n - 2;
   }

   for (size_t i = 1; i < (n - 2); i++)
   {
     if( ( (val - brks[i+1]) * (val - brks[i]) ) <= 0 )
     {
       return i;
     }
   }
 }


