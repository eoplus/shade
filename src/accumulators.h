
 #ifndef ACUMULATORS
 #define ACUMULATORS

 #include "skyrad.h"
 #include "config.h"

/*******************************************************************************
 accm_bmc

 This struct implements the accumulators for Backward Monte Carlo integration.
 It contains a three-dimensional array that accumulates from -INFINITY to
 INFINITY spatial range. If using spatially resolved integral (pre-processor 
 flag: SPATIALLY_RESOLVED), a five-dimensional array contains the spatially
 resolved accumulators.

 Since parameters that do not affect directionality can be included in a single
 simulation, the accumulator resolves all combinations of single scattering 
 albedos and bi-hemispherical reflectance of the bottom. Additionally, because
 this version of the code only resolves scattering in water (that is, the sky 
 radiance is prescribed), the directionality of the incoming light field does
 not afffect the simulation per se, only its scaling during accumulation. 
 Therefore, the accumulators also resolve all combinations of the previous
 properties with illumination condition.

 COMPONENTS (accumulator_bmc):
 ns     - Number of sun zenith angles (or sky radiance distributions);
 nw0    - Number of single scattering albedos;
 nbr    - Number of bottom reflectances;
 ny     - Number of rows;
 nx     - Number of columns;
 ybrks  - Pointer to row coordinate break points;
 xbrks  - Pointer to column coordinate break points;
 grid   - Pointer to five dimensional array [ns][nw0][nbr][ny][nx];
 vector - Pointer to three array [ns][nw0][nbr];
 kx_inv - Scaling for use in the findBin function;
 ky_inv - Scaling for use in the findBin function;
 minx   - Practical minimum for use in the findBin function;
 miny   - Practical minimum for use in the findBin function;
 fun    - accumulator_fun struct;
 type   - type of the accumulator;

*******************************************************************************/

 struct accumulator_bmc;
 struct accumulator_fun;

 typedef void
 accm_add_fun
 (
   struct accumulator_bmc *,
   struct skyradiance const *, 
   double const *, 
   double const *, 
   double const ***,
   double const, 
   int const,
   int const
 );

 struct accumulator_fun
 {
   accm_add_fun * add;
 };

 struct accumulator_bmc
 {
   int    ns;
   int    nw0;
   int    nbr;
   int    ny;
   int    nx;
   double * ybrks;
   double * xbrks;
   double ***** grid;
   double *** vector;
   double kx_inv;
   double ky_inv;
   double minx;
   double miny;
   struct accumulator_fun fun;
   char   type[STRMXLEN];
 };

/* Add to backward accumulator: ************************************************

 accm_b_add

 Add contribution of rays to a backward accumulator. This function is a generic 
 function to call the specific function for the accumulator type. It is provided 
 as an inline function as as a macro.

 For details of the functions, see the functions 'accm_b_add_sect' and 
 'accm_b_add_grid' in the accumulators.c source file.

 INPUT:
 acm   - Accumulator parameters
         Pointer to accumulator_bmc struct;
 skr   - Sky radiance distribution parameters
         Pointer to constant skyradiance struct;
 p     - Terminal point of the ray
         Range: (-Inf,Inf), meters
         [0]X, [1]Y, [2]Z
         Pointer to constant double;
 s     - Refracted spherical directions of the ray
         Range: [0,PI], [0,2PI], radians
         [0]Theta, [1]Phi, [2]Radius (always 1.0 and never used)
         Pointer to constant double;
 stks  - Diffuse Stokes vector
         [N w0][N bhr][4]: [][][0]I, [][][1]Q, [][][2]U, [][][3]V
         Pointer to pointer to pointer to constant double;
 scale - Arbitrary multiplier (e.g., Fresnel transmittance or probability to 
         Sun)
         Range: (0,Inf), unitless
         Constant double;
 dirf  - Flag to indicate if the contribution is for the direct component
         Values: (0) No, (1) Yes
         Constant int
 cs    - Index of the accumulator layer for the SZA in case of direct component
         Values: [0, N SZA]
         Constant int.

 OUTPUT:
 None. Updates the accumulator structure.

*******************************************************************************/

 static inline __attribute__((always_inline)) void
 accm_b_add
 (
   struct accumulator_bmc * accm,
   struct skyradiance const * skr,
   double const * pos,
   double const * s,
   double const *** stks,
   double const scale,
   int const dirf,
   int const cs
 )
 { accm->fun.add(accm, skr, pos, s, stks, scale, dirf, cs); }

 #define ACCM_B_ADD(I_accm, I_skr, I_pos, I_s, I_stks, I_scale, I_dirf, I_cs) \
         (I_accm)->fun.add( (I_accm), (I_skr), (I_pos), (I_s), \
         (double const ***) (I_stks), (I_scale), (I_dirf), (I_cs) );

/* Function prototypes: *******************************************************/

 // Setup functions:

 struct accumulator_bmc *
 accm_b_alloc
 ( void );

 void
 accm_b_free
 ( struct accumulator_bmc ** accm );

 void
 accm_setup
 (
   struct accumulator_bmc * accm,
   char const * accm_tp,
   int const sim_ns,
   int const iop_nw0,
   int const btt_nbr,
   double const acc_ext,
   double const acc_resy,
   double const acc_resx
 );


 // Sectorial functions:

 void
 accm_b_setup_sect
 (
   struct accumulator_bmc * accm,
   int const sim_ns,
   int const iop_nw0,
   int const btt_nbr,
   double const acc_ext,
   double const acc_resy,
   double const acc_resx
 );

 void
 accm_b_add_sect
 (
   struct accumulator_bmc * accm,
   struct skyradiance const * skr,
   double const * p,
   double const * s,
   double const *** stks,
   double const scale,
   int const dirf,
   int const cs
 );


 // Regular grid functions:

 void
 accm_b_setup_grid
 (
   struct accumulator_bmc * accm,
   int const sim_ns,
   int const iop_nw0,
   int const btt_nbr,
   double const acc_ext,
   double const acc_resy,
   double const acc_resx
 );

 void
 accm_b_add_grid
 (
   struct accumulator_bmc * accm,
   struct skyradiance const * skr,
   double const * p,
   double const * s,
   double const *** stks,
   double const scale,
   int const dirf,
   int const cs
 );


 // Statistical functions:

 void
 accm_b_norm
 (
   struct accumulator_bmc * accm,
   double const f0,
   double const * normf
 );

 void 
 accm_b_add_k
 (
   struct accumulator_bmc * accm,
   double const val,
   int const cs,
   int const cr,
   int const cc
 );

 void
 accm_b_sum
 (
   struct accumulator_bmc * accm_sum,
   struct accumulator_bmc const * accm_in,
   double const scale
 );

 void
 accm_b_reset
 ( struct accumulator_bmc * accm );


 // Output functions:

 void
 accm_b_write_grid
 (
   struct accumulator_bmc * accm,
   int const sim_ns,
   char const * ofbn,
   char const * sufx,
   double const * iop_w0,
   double const * btt_bh,
   double const * sim_sza
 );

 void
 accm_b_write_vec
 (
   struct accumulator_bmc * accm,
   int const sim_ns,
   char const * ofbn,
   char const * sufx,
   double const * iop_w0,
   double const * btt_bhr,
   double const * sim_sza
 );

 #endif // ACUMULATORS

