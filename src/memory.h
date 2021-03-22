
 #ifndef MEMORY
 #define MEMORY

 /* Function protopytes: ******************************************************/

 double *
 calloc_1d
 (
   int const nc, 
   char const * varnm
 ); 

 void
 free_1d
 ( double ** array );

 double **
 calloc_2d
 (
   int const nr,
   int const nc,
   char const * varnm
 ); 

 void
 free_2d
 (
   int const nr,
   double *** array
 );

 double ***
 calloc_3d
 (
   int const nl, 
   int const nr, 
   int const nc,
   char const * varnm
 );

 void
 free_3d
 (
   int const nl,
   int const nr,
   double **** array
 );

 double ****
 calloc_4d
 (
   int const nm, 
   int const nl, 
   int const nr, 
   int const nc,
   char const * varnm
 );

 void
 free_4d
 (
   int const nm,
   int const nl,
   int const nr,
   double ***** array
 );

 double *****
 calloc_5d
 (
   int const ng, 
   int const nm, 
   int const nl, 
   int const nr, 
   int const nc,
   char const * varnm
 );

 void
 free_5d
 (
   int const ng,
   int const nm,
   int const nl,
   int const nr,
   double ****** array
 );

 #endif // MEMORY

