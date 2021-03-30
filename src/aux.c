
/*******************************************************************************
 aux.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 Version: 1.6
 Date: 2021-03-25
 License: GPL-3.0

 This source file and the homonymous header file provide functions auxiliary
 functions.

 FUNCTIONS AND MACROS LIST:
 findInterv - Retrieve the index of the interval in which a value lies
 findBin    - Retrieve the index of the regular interval in which a value lies
 FINDBIN    - Retrieve the index of the regular interval in which a value lies
              (macro)
 ABS        - Retrieve the absolute value of a number
 SIGN       - Retrieve the sign of a number (zeros == 0)
 SIGN_N0    - Retrieve the sign of a number (zeros == +)
 NUM_EQU    - Near equivalence of two real numbers

*******************************************************************************/

 #include <math.h>	// floor

 #include "aux.h"

/* Find interval that contain a value: *****************************************

 findInterv

 Function to find the interval index where a value lies. It is open at the left
 and closed in the right. Breaks extremes can be be -INFINITE and INFINITE.

 However, is not very efficient... if the break points form regular intervals, 
 the function finBin (aux.h header file) is many times faster.
 
 INPUT: 
 val  - Scalar value
        Range: (-Inf,Inf), undefined units
        Constant double;
 brks - Vector of breaks
        Range: [-Inf,Inf], undefined units
        Pointer to constant double;
 n    - Number of break points
        Integer range: (0,Inf), unitless
        Constant int.

 OUTPUT:
 The index of the interval in which the value lies.

*******************************************************************************/

 int
 findInterv
 (
   double const val,
   double const * brks,
   int const n
 )
 {
   /*
   In this code, all sequence of breaks start at -INFINITY and end at INFINITY.
   That means that a value is always in some interval. If lower that the 2nd
   value it is between -INFINITY and the 2nd value. If larger than the nth-1
   value it is between the nth-1 and INFINITY.
   */
   if ( val <= brks[1] )
   {
     return 0;
   }
   else if ( val > brks[n - 2] )
   {
     return n - 1;
   }

   /*
   If passing first condition, loop over the break points to find the interval
   that contains the value.
   */
   for (int i = 1; i < (n - 2); i++)
   {
     if( ( (val - brks[i + 1]) * (val - brks[i]) ) <= 0 )
     {
       return i;
     }
   }

   /*
   All possibilities are covered, but gcc flag '-Wreturn-type' will complain
   that the control ended without return for a non void function. The following 
   will never happen.
   */
   return 0;
 }


