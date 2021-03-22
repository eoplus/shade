
 #ifndef SOURCES
 #define SOURCES

 #include <stdio.h>
 #include "config.h"
 #include "constants.h"
 
/*******************************************************************************
 source

 This structure defines the light source properties.

 COMPONENTS:
 lnf      - Normalization factor for conical sources (e.g., scalar irradiance and
            radiance sensors as sources)
            Double;
 o        - Origin (center) of the source
            Range: (-Inf,Inf), meters
            Array of double: [0]X, [1]Y, [2]Z;
 rotate_f - Flag to indicate if emitted ray needs to be rotated to system's
            reference frame.
            0 - No rotation, 1 - Rotate;
 M        - Rotation matrix from the source reference frame to the system's
            reference frame
            [3][3]
            Pointer to pointer to double;
 fun      - Set of functions for the source
            Pointer to source_fun struct;
 u        - Cartesian coordinates of source axis
            [0]X, [1]Y, [2]Z
            Range: [-1, 1], unitless
            Pointer to double;
 s        - Spherical coordinates of the axis
            [0]Theta, [1]Phi, [2]Radius (= 1.0);
            Range: [0,PI], [0,2PI], radians;
 fov      - Filed-of-view of the sensor
            Range: [0,2 PI], radians
            Double;
 hfov     - Half-angle of the field-of-view of the sensor
            Range: [0,PI], radians
            Double;
 stks     - Stokes parameters of the source
            [0]I, [1]Q, [2]U, [3]V
            Array of length 1 if scalar and 4 if vector RTE;
 tp       - Source type name;


 To allow flexibility and easy of extension, generic functions are defined, wich
 will call the correct function for the source type, as the function pointers
 are included in the source structure. The function pointers are defined with
 typedef and contained in an auxiliary structure source_fun. The correct
 function adresses are copied by the setup function.

 COMPONENTS:
 qtl - Quantile function pointer;
 cdf - CDF function pointer;
 pdf - PDF function pointer;
 ef  - Emission function pointer;

*******************************************************************************/

 struct source;
 struct source_fun;

 typedef void src_qtl_fun (struct source const *, double const *, double *);
 typedef void src_cdf_fun (struct source const *, double *, double const *);
 typedef void src_pdf_fun (struct source const *, double *, double const *);
 typedef void src_ef_fun  (struct source const *, double *, double const *);

 struct source_fun
 {
   src_qtl_fun * qtl;
   src_cdf_fun * cdf;
   src_pdf_fun * pdf;
   src_ef_fun  * ef;
 };

 struct source
 {
   double o[3];
   int    rotate_f;
   double ** M;
   double lnf;
   double lnf_inv;
   struct source_fun fun;
   #ifdef VECTOR_RT
   double stks[4];
   #else
   double stks[1];
   #endif // VECTOR_RT
   double u[3];
   double s[3];
   double fov;
   double hfov;
   char   tp[STRMXLEN];
 };

/* Generic source functions: ***************************************************

 src_qtl, src_cdf, src_pdf, src_ef

 Generic functions to retrieve the quantile, CDF, PDF and emission function for
 a given source.

 INPUT:
 src    - Source parameters
          Pointer to constant source struct; 
 prob   - Cumulative probability associated with a polar and a azimuth angle
          [0]Theta, [1]Phi
          Range: [0,1] unitless 
          Pointer to constant double or Pointer to double;
 stks   - Stokes parameters
          If VECTOR_RT: [0]I, [1]Q, [2]U, [3]V, else: [0]I
          Pointer to double;
 s_src  - Spherical direction
          [0]Theta, [1]Phi, [2]Radius
          Range: [0,PI] radians, [0,2PI] radians, 1.0 unitless
          Pointer to constant double or Pointer to double;

 OUTPUT:
 None. Updates values pointed by s_src or prob.

*******************************************************************************/

 static inline __attribute__((always_inline)) void
 src_qtl
 (
   struct source const *src,
   double *prob,
   double *s_src
 )
 { src->fun.qtl(src, prob, s_src); }

 #define SRC_QTL(I_src, I_prob, I_s_src) \
         (I_src)->fun.qtl( (I_src), (I_prob), (I_s_src) );

 static inline __attribute__((always_inline)) void
 src_cdf
 (
   struct source const *src,
   double *prob,
   double const *s_src
 )
 { src->fun.cdf(src, prob, s_src); }

 #define SRC_CDF(I_src, I_prob, I_s_src) \
         (I_src)->fun.cdf( (I_src), (I_prob), (I_s_src) );

 static inline __attribute__((always_inline)) void
 src_pdf
 (
   struct source const *src,
   double *prob,
   double const *s_src
 )
 { src->fun.pdf(src, prob, s_src); }

 #define SRC_PDF(I_src, I_dens, I_s_src) \
         (I_src)->fun.pdf( (I_src), (I_dens), (I_s_src) );

 static inline __attribute__((always_inline)) void
 src_ef
 (
   struct source const *src,
   double *stks,
   double const *s_src
 )
 { src->fun.ef(src, stks, s_src); }

 #define SRC_EF(I_src, I_stks, I_s_src) \
         (I_src)->fun.ef( (I_src), (I_stks), (I_s_src) );

/* Function prototypes: *******************************************************/

 // Setup functions:
 struct source*
 src_alloc
 ( void );

 void 
 src_free
 ( struct source **src );

 void
 src_setup
 (
   struct source *src,
   double const src_fov,
   double const *src_s,
   double const *src_ref_o,
   double const *src_rel_o,
   double const *src_stks,
   char const   *src_tp
 );

 void
 src_fprintf
 ( 
   FILE *odv,
   struct source const *src,
   int const indent
 );

 // Point sources, scalar irradiance sensors as sources:
 void
 src_qtl_pnt_sir
 (
   struct source const *src,
   double const *prob,
   double *s_src
 );

 void
 src_cdf_pnt_sir
 (
   struct source const *src,
   double *prob,
   double const *s_src
 );

 void
 src_pdf_pnt_sir
 (
   struct source const *src,
   double *prob,
   double const *s_src
 );

 void
 src_ef_pnt_sir
 (
   struct source const *src,
   double *stks,
   double const *s_src
 );


 // Radiance sensors as sources:
 void
 src_qtl_pnt_rad
 (
   struct source const *src,
   double const *prob,
   double *s_src
 );

 void
 src_cdf_pnt_rad
 (
   struct source const *src,
   double *prob,
   double const *s_src
 );

 void
 src_pdf_pnt_rad
 (
   struct source const *src,
   double *prob,
   double const *s_src
 );

 void
 src_ef_pnt_rad
 (
   struct source const *src,
   double *stks,
   double const *s_src
 );


 // Plane sources (e.g., irradiance sensors as sources):
 void
 src_qtl_pnt_pir
 (
   struct source const *src,
   double const *prob,
   double *s_src
 );

 void
 src_cdf_pnt_pir
 (
   struct source const *src,
   double *prob,
   double const *s_src
 );

 void
 src_pdf_pnt_pir
 (
   struct source const *src,
   double *prob,
   double const *s_src
 );

 void
 src_ef_pnt_pir
 (
   struct source const *src,
   double *stks,
   double const *s_src
 );

 #endif // SOURCES

