
/*******************************************************************************
 bottom.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 2021-03-03

 Provides functions describing the bottom.

 In the current version, horizonthal homogeneity is assumed and bottom is
 described by a single reflectance model (though, several reflectance values
 are possible in a single run) and single depth.

 This is expected to be expanded in the future to incorporate spatially resolved
 simulations, with the bottom represented by a triangular mesh, and so with
 variable depth, orientation and possibility to associate different reflectance
 "classes" to each cell.

*******************************************************************************/

 #include <stdlib.h>
 #include <stdio.h>
 #include <string.h>
 #include <math.h>

 #include "config.h"
 #include "reflectance.h"
 #include "bottom.h"

/* Allocate and deallocate bottom struct: **************************************

 btt_alloc, btt_free

 Allocate and deallocate memory for a bottom struct. When allocating, values are
 initialized to zero and when deallocating, the pointer is set to NULL.

 INPUT:
 btt - Bottom parameters
       Pointer to pointer to bottom struct;

 OUTPUT:
 The adress of a bottom struct memory block or none, when deallocating.
 
*******************************************************************************/

 struct bottom *
 btt_alloc
 ( void )
 { 
   struct bottom * btt = (struct bottom*) calloc(1, sizeof(struct bottom));
   #ifdef CHCK_MEM_ALLOC
   if ( !btt )
   {
     printf("\nERROR: Could not allocate memory for btt in btt_alloc"
       " (bottom.c)\n");
     exit(-1);
   }
   #endif // CHCK_MEM_ALLOC

   return btt;
 }

 void
 btt_free
 ( struct bottom ** btt )
 {
   free( *btt );
   *btt = NULL;
 }

/* Setup bottom properties: ****************************************************

 btt_setup

 Setup a bottom struct.

 INPUT:
 btt  - Bottom parameters
         Pointer to bottom struct;
 depth - Single bottom depth
         Range: (0,Inf), meters
         Constant double;
 type  - Reflectance model of the bottom
         Options: See the reflectance.c source file
         Constant char;
 nbhr  - Number of bi-hemispherical reflectances
         Range: (0,Inf), unitless
         Constant int;
 bhr   - Bi-hemispherical reflectances
         Range: [0,1], unitless
         Pointer to constant double;
 k     - Expoent for the Minnaert BRDF model
         Range: unspecified
         Constant double;

 OUTPUT:
 None. Updates the values pointed by btt.
 
*******************************************************************************/

 void
 btt_setup
 (
   struct bottom * btt,
   double const depth,	// Depth, meters
   char const * type,	// Reflectance model
   int const nbhr,	// Number of bi-hemispherical reflectances
   double const * bhr,	// bi-hemispherical reflectances
   double const k	// Minnaert exponent
 )
 {
   btt->depth = depth;
   refl_setup(&btt->refl, type, nbhr, bhr, k);
 }

/* Print bottom properties: ****************************************************

 btt_fprintf

 Prints the bottom struct.

 INPUT:
 btt   - Bottom parameters
          Pointer to bottom struct;
 indent - Number of indentations at the header level (each indent = 2 spaces)
          Constant int;

 OUTPUT:
 None. Prints values pointed by btt.
 
*******************************************************************************/

 void
 btt_fprintf
 (
   FILE * odv,
   struct bottom const * btt,
   int const indent 
 )
 {
   char pre_0[STRMXLEN] = "";
   char pre_1[STRMXLEN] = "";
   for (size_t i = 0; i < indent; i++)
     strcat(pre_0, "  ");
   strncpy(pre_1, pre_0, STRMXLEN);
   strcat(pre_1, "  ");

   fprintf(odv, "%sBottom parameters:\n", pre_0);
   fprintf(odv, "%sDepth: %.2lf (m)\n", pre_1, btt->depth);
   if ( isfinite( btt->depth ) )
   {
     refl_fprintf(odv, &btt->refl, indent + 1);
   }
   else
   {
     fprintf(odv, "%sReflectance parameters: Undefined\n", pre_1);
   }
 }


