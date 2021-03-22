
 #ifndef SKYRAD
 #define SKYRAD

/*******************************************************************************
 skyrad

 COMPONENTS:
 ns     - Number of sky radiance distributions;
 ny     - Number of rows;
 nx     - Number of columns;
 ybrks  - Pointer to row coordinate break points;
 xbrks  - Pointer to column coordinate break points;
 grid;  - Pointer to three dimensional array [ns][ny][nx];
 vector - Pointer to unidimensional array [ns];
 kx     - Scaling for use in the findBin function;
 ky     - Scaling for use in the findBin function;
 minx   - Practical minimum for use in the findBin function;
 miny   - Practical minimum for use in the findBin function;

*******************************************************************************/

 struct skyradiance
 {
   int    ns;
   int    ny;
   int    nx;
   double *ybrks;
   double *xbrks;
   double ***grid;
   double *vector;
   double kx_inv;
   double ky_inv;
   double minx;
   double miny;
 };

 // Function prototypes: *******************************************************

 struct skyradiance*
 skr_alloc
 ( void );

 void
 skr_free
 ( struct skyradiance **skr );

 void
 skr_setup
 (
   struct skyradiance *skr,
   int const sim_ns,
   int const skr_ny,
   int const skr_nx,
   double const skr_resy,
   double const skr_resx,
   char const **skr_fls,
   double const *sim_sza
 );

 void
 skr_read
 (
   struct skyradiance *skr,
   char const **skr_fls
 );

 #endif // SKYRAD

