
/******************************************************************************
 memory.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 2021-03-03

 Provides functions to allocate and free n-dimensional arrays of doubles, where
 n is from 1 to 5. The functions test for memory allocation and updated the 
 pointer to NULL when releasing memory.

*******************************************************************************/

 #include <stdio.h>
 #include <stdlib.h>

 #include "memory.h"

/* One-dimensional array alloc: ************************************************

 calloc_1d, free_1d

 Functions to allocate and free one-dimensional arrays of doubles. Performs a 
 test for memory allocation.

 INPUT:
 nc    - Number of columns.
 varnm - Name of the variable or NULL;

 OUTPUT:
 Pointer to an array of doubles (alloc) or none (free), updating the pointer to 
 NULL.

*******************************************************************************/

 double *
 calloc_1d
 (
   int const nc, 
   char const * varnm
 ) 
 {
   double * array = (double*) calloc(nc, sizeof (double));

   #ifdef CHCK_MEM_ALLOC
   if(array == NULL)
   {
     printf("\nERROR: Could not allocate memory for %s.\n", varnm);
     exit(-2);
   }
   #endif // CHCK_MEM_ALLOC

   return array;
 }

 void
 free_1d
 ( double ** array )
 {
   if( *array )
   {
     free( *array );
     *array = NULL;
   }
 }

/* Two-dimensional array alloc: ************************************************

 calloc_2d, free_2d

 Functions to allocate and free bi-dimensional arrays of doubles. Performs a 
 test for memory allocation.

 INPUT:
 nr    - Number of rows.
 nc    - Number of columns.
 varnm - Name of the variable or NULL;

 OUTPUT:
 Pointer to an array of pointers of doubles (alloc) or none (free), updating 
 the pointer to NULL.

*******************************************************************************/

 double **
 calloc_2d
 (
   int const nr,
   int const nc,
   char const * varnm
 ) 
 {
   double ** array = (double **) calloc(nr, sizeof (double *));

   #ifdef CHCK_MEM_ALLOC
   if( !array )
   {
     printf("\nERROR: Could not allocate memory for %s.\n", varnm);
     exit(-2);
   }
   #endif // CHCK_MEM_ALLOC

   for (size_t cr = 0; cr < nr; cr++)
   {
     array[cr] = (double*) calloc(nc, sizeof(double));
     #ifdef CHCK_MEM_ALLOC
     if( !array[cr] )
     {
       printf("\nERROR: Could not allocate memory for index [%d] of %s.\n", 
         cr, varnm);
       exit(-2);
     }
     #endif // CHCK_MEM_ALLOC
   }

   return array;
 }

 void
 free_2d
 (
   int const nr,
   double *** array
 )
 {
   if( *array )
   {
     for (size_t cr = 0; cr < nr; cr++)
     {
       free( (*array)[cr] );
     }
     free( *array );
     *array = NULL;
   }
 }

/* Three-dimensional array alloc: **********************************************

 calloc_3d, free_3d

 Functions to allocate and free three-dimensional arrays of doubles. Performs a 
 test for memory allocation.

 INPUT:
 nl    - Number of layers;
 nr    - Number of rows;
 nc    - Number of columns;
 varnm - Name of the variable or NULL.

 OUTPUT:
 Pointer to an array of pointers of an array of pointer of doubles (alloc) or 
 none (free), updating the pointer to NULL.

*******************************************************************************/

 double ***
 calloc_3d
 (
   int const nl, 
   int const nr, 
   int const nc,
   char const * varnm
 )
 {
   double *** array = (double ***) calloc(nl, sizeof (double **));

   #ifdef CHCK_MEM_ALLOC
   if( !array )
   {
     printf("\nERROR: Could not allocate memory for %s.\n", varnm);
     exit(-2);
   }
   #endif // CHCK_MEM_ALLOC

   for (size_t cl = 0; cl < nl; cl++)
   {
     array[cl] = (double**) calloc(nr, sizeof (double*));

     #ifdef CHCK_MEM_ALLOC
     if ( !array[cl] )
     {
       printf("\nERROR: Could not allocate memory for index [%d] of %s.\n", 
         cl, varnm);
       exit(-2);
     }
     #endif // CHCK_MEM_ALLOC

     for (size_t cr = 0; cr < nr; cr++)
     {
       array[cl][cr] = (double*) calloc(nc, sizeof (double));

       #ifdef CHCK_MEM_ALLOC
       if ( !array[cl][cr] )
       {
         printf("\nERROR: Could not allocate memory for index [%d][%d] of "
           "%s.\n", cl, cr, varnm);
         exit(-2);
       }
       #endif // CHCK_MEM_ALLOC
     }
   }

   return array;
 }

 void
 free_3d
 (
   int const nl,
   int const nr,
   double **** array
 )
 {
   if( *array )
   {
     for (size_t cl = 0; cl < nl; cl++)
     {
       for (size_t cr = 0; cr < nr; cr++)
       {
         free( (*array)[cl][cr] );
       }
       free( (*array)[cl] );
     }
     free( *array );
     *array = NULL;
   }
 }

/* Four-dimensional array alloc: ***********************************************

 calloc_4d, free_4d

 Functions to allocate and free four-dimensional arrays of doubles. Performs a 
 test for memory allocation.

 INPUT:
 nm    - Number of members;
 nl    - Number of layers;
 nr    - Number of rows;
 nc    - Number of columns;
 varnm - Name of the variable or NULL.

 OUTPUT:
 Pointer to an array of pointers of an array of pointer of an array of pointers 
 to doubles (alloc) or none (free), updating the pointer to NULL.

*******************************************************************************/

 double ****
 calloc_4d
 (
   int const nm, 
   int const nl, 
   int const nr, 
   int const nc,
   char const * varnm
 )
 {
   double **** array = (double ****) calloc(nm, sizeof (double ***));

   #ifdef CHCK_MEM_ALLOC
   if ( !array )
   {
     printf("\nERROR: Could not allocate memory for %s.\n", varnm);
     exit(-2);
   }
   #endif // CHCK_MEM_ALLOC

   for (size_t cm = 0; cm < nm; cm++)
   {
     array[cm] = (double***) calloc(nl, sizeof (double**));

     #ifdef CHCK_MEM_ALLOC
     if ( !array[cm] )
     {
       printf("\nERROR: Could not allocate memory for index [%d] of %s.\n", 
         cm, varnm);
       exit(-2);
     }
     #endif // CHCK_MEM_ALLOC

     for (size_t cl = 0; cl < nl; cl++)
     {
       array[cm][cl] = (double**) calloc(nr, sizeof (double*));

       #ifdef CHCK_MEM_ALLOC
       if ( !array[cm][cl] )
       {
         printf("\nERROR: Could not allocate memory for index [%d][%d] of "
           "%s.\n", cm, cl, varnm);
         exit(-2);
       }
       #endif // CHCK_MEM_ALLOC

       for (size_t cr = 0; cr < nr; cr++)
       {
         array[cm][cl][cr] = (double*) calloc(nc, sizeof (double));

         #ifdef CHCK_MEM_ALLOC
         if ( !array[cm][cl][cr] )
         {
           printf("\nERROR: Could not allocate memory for index [%d][%d][%d] of"
             " %s.\n", cm, cl, cr, varnm);
           exit(-2);
         }
         #endif // CHCK_MEM_ALLOC
       }
     }
   }

   return array;
 }

 void
 free_4d
 (
   int const nm,
   int const nl,
   int const nr,
   double ***** array
 )
 {
   if( *array )
   {
     for (size_t cm = 0; cm < nm; cm++)
     {
       for (size_t cl = 0; cl < nl; cl++)
       {
         for (size_t cr = 0; cr < nr; cr++)
         {
           free( (*array)[cm][cl][cr] );
         }
         free( (*array)[cm][cl] );
       }
       free( (*array)[cm] );
     }
     free( *array );
     *array = NULL;
   }
 }

/* Five-dimensional array alloc: ***********************************************

 calloc_5d, free_5d

 Functions to allocate and free five-dimensional arrays of doubles. Performs a 
 test for memory allocation.

 INPUT:
 ng    - Number of groups;
 nm    - Number of members;
 nl    - Number of layers;
 nr    - Number of rows;
 nc    - Number of columns;
 varnm - Name of the variable or NULL.

 OUTPUT:
 Pointer to an array of pointers of an array of pointer of an array of pointers 
 of an array of pointers to doubles (alloc) or none (free), updating the pointer
 to NULL.

*******************************************************************************/

 double *****
 calloc_5d
 (
   int const ng, 
   int const nm, 
   int const nl, 
   int const nr, 
   int const nc,
   char const * varnm
 )
 {
   double ***** array = (double *****) calloc(ng, sizeof (double ****));

   #ifdef CHCK_MEM_ALLOC
   if ( !array )
   {
     printf("\nERROR: Could not allocate memory for %s.\n", varnm);
     exit(-2);
   }
   #endif // CHCK_MEM_ALLOC

   for (size_t cg = 0; cg < ng; cg++)
   {
     array[cg] = (double ****) calloc(nm, sizeof (double ***));

     #ifdef CHCK_MEM_ALLOC
     if ( !array[cg] )
     {
       printf("\nERROR: Could not allocate memory for index [%d] of %s.\n", 
         cg, varnm);
       exit(-2);
     }
     #endif // CHCK_MEM_ALLOC

     for (size_t cm = 0; cm < nm; cm++)
     {
       array[cg][cm] = (double ***) calloc(nl, sizeof (double **));

       #ifdef CHCK_MEM_ALLOC
       if ( !array[cg][cm] )
       {
         printf("\nERROR: Could not allocate memory for index [%d][%d] of "
           "%s.\n", cg, cm, varnm);
         exit(-2);
       }
       #endif // CHCK_MEM_ALLOC

       for (size_t cl = 0; cl < nl; cl++)
       {
         array[cg][cm][cl] = (double **) calloc(nr, sizeof (double *));

         #ifdef CHCK_MEM_ALLOC
         if ( !array[cg][cm][cl] )
         {
           printf("\nERROR: Could not allocate memory for index [%d][%d][%d] of"
             " %s.\n", cg, cm, cl, varnm);
           exit(-2);
         }
         #endif // CHCK_MEM_ALLOC

         for (size_t cr = 0; cr < nr; cr++)
         {
           array[cg][cm][cl][cr] = (double *) calloc(nc, sizeof (double));

           #ifdef CHCK_MEM_ALLOC
           if ( !array[cg][cm][cl][cr] )
           {
             printf("\nERROR: Could not allocate memory for index [%d][%d][%d]"
               "[%d] of %s.\n", cg, cm, cl, cr, varnm);
             exit(-2);
           }
           #endif // CHCK_MEM_ALLOC
         }
       }
     }
   }

   return array;
 }

 void
 free_5d
 (
   int const ng,
   int const nm,
   int const nl,
   int const nr,
   double ****** array
 )
 {
   if( *array )
   {
     for (size_t cg = 0; cg < ng; cg++)
     {
       for (size_t cm = 0; cm < nm; cm++)
       {
         for (size_t cl = 0; cl < nl; cl++)
         {
           for (size_t cr = 0; cr < nr; cr++)
           {
             free( (*array)[cg][cm][cl][cr] );
           }
           free( (*array)[cg][cm][cl] );
         }
         free( (*array)[cg][cm] );
       }
       free( (*array)[cg] );
     }
     free( *array );
     *array = NULL;
   }
 }


