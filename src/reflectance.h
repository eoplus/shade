
 #ifndef REFLECTANCE
 #define REFLECTANCE

 #include "config.h"

/*******************************************************************************


*******************************************************************************/

 struct reflectance;
 struct reflectance_fun;
 struct reflectance_minn;
 struct reflectance_lamb;

 typedef void
 refl_qtl_fun
 (
   struct reflectance const *,
   double const *,
   double *, 
   double const *
 );

 typedef void
 refl_cdf_fun
 (
   struct reflectance const *,
   double const *, 
   double const *,
   double *
 );

 typedef void
 refl_pdf_fun
 (
   struct reflectance const *,
   double const *, 
   double const *,
   double *
 );

 typedef void 
 refl_dhr_fun
 (
   struct reflectance const *,
   double const *,
   double *
 );

 typedef void
 refl_brdf_fun
 (
   struct reflectance const *,
   double const *, 
   double const *,
   double *
 );

 struct reflectance_fun
 {
   refl_qtl_fun * qtl;
   refl_cdf_fun * cdf;
   refl_pdf_fun * pdf;
   refl_dhr_fun * dhr;
   refl_brdf_fun * brdf;
 };

 struct reflectance_minn
 {
   int    nbhr;
   double * bdr;		// Fixed part of the bi-directional reflectance
   double * bhr;		// Fixed part of the bi-hemispherical reflectance
   double k;		// Exponential parameter
   double kp1;		// k + 1
   double kp2;		// k + 2
   double invkp2;	// 1 / (k + 2)
 };

 struct reflectance_lamb
 {
   int    nbhr;
   double * bdr;		// Bi-directional reflectance
   double * bhr;		// Bi-hemispherical reflectance
 };

 struct reflectance
 {
   struct reflectance_fun  fun;
   struct reflectance_minn minn;
   struct reflectance_lamb lamb;
   char tp[STRMXLEN];
 };

/* Generic source functions: ***************************************************



*******************************************************************************/

 static inline __attribute__((always_inline)) void
 refl_qtl
 (
   struct reflectance const *refl,
   double const *s_scat_i,
   double       *s_scat_r,
   double const *prob
 )
 {
   refl->fun.qtl(refl, s_scat_i, s_scat_r, prob);
 }

 #define REFL_QTL(I_refl, I_s_scat_i, I_s_scat_r, I_prob) \
         (I_refl)->fun.qtl( (I_refl), (I_s_scat_i), (I_s_scat_r), (I_prob) );

 static inline __attribute__((always_inline)) void
 refl_cdf
 (
   struct reflectance const *refl,
   double const *s_scat_i,
   double const *s_scat_r,
   double       *prob
 )
 {
   refl->fun.cdf(refl, s_scat_i, s_scat_r, prob);
 }

 #define REFL_CDF(I_refl, I_s_scat_i, I_s_scat_r, I_prob) \
         (I_refl)->fun.cdf( (I_refl), (I_s_scat_i), (I_s_scat_r), (I_prob) );

 static inline __attribute__((always_inline)) void
 refl_pdf
 (
   struct reflectance const *refl,
   double const *s_scat_i,
   double const *s_scat_r,
   double       *dens
 )
 {
   refl->fun.pdf(refl, s_scat_i, s_scat_r, dens);
 }

 #define REFL_PDF(I_refl, I_s_scat_i, I_s_scat_r, I_dens) \
         (I_refl)->fun.pdf( (I_refl), (I_s_scat_i), (I_s_scat_r), (I_dens) );

 static inline __attribute__((always_inline)) void
 refl_dhr
 (
   struct reflectance const *refl,
   double const *s_scat_i,
   double       *dhr
 )
 {
   refl->fun.dhr(refl, s_scat_i, dhr);
 }

 #define REFL_DHR(I_refl, I_s_scat_i, I_dhr) \
         (I_refl)->fun.dhr( (I_refl), (I_s_scat_i), (I_dhr) );

 static inline __attribute__((always_inline)) void
 refl_brdf
 (
   struct reflectance const *refl,
   double const *s_scat_i,
   double const *s_scat_r,
   double       *bdr
 )
 {
   refl->fun.brdf(refl, s_scat_i, s_scat_r, bdr);
 }

 #define REFL_BRDF(I_refl, I_s_scat_i, I_s_scat_r, I_bdr) \
         (I_refl)->fun.brdf( (I_refl), (I_s_scat_i), (I_s_scat_r), (I_bdr) );


/* Function prototypes: *******************************************************/

 // Setup functions:

 struct reflectance *
 refl_alloc
 ( void );

 void
 refl_free
 ( struct reflectance ** refl );

 void
 refl_setup
 (
   struct reflectance * refl,
   char const * tp,
   int const nbhr,
   double const * bhr,
   double const k
 );

 void
 refl_fprintf
 (
   FILE * odv,
   struct reflectance const * refl,
   int const indent
 );


 // Minnaert:

 struct reflectance_minn *
 refl_alloc_minn
 ( void );

 void
 refl_free_minn
 ( struct reflectance_minn ** minn );

 void
 refl_setup_minn
 ( 
   struct reflectance * refl,
   int const nbhr, 
   double const * bhr,
   double const k
 );

 void
 refl_qtl_minn
 (
   struct reflectance const * refl,
   double const * s_scat_i,
   double * s_scat_r,
   double const * prob
 );

 void
 refl_cdf_minn
 (
   struct reflectance const * refl,
   double const * s_scat_i,
   double const * s_scat_r,
   double * prob
 );

 void
 refl_pdf_minn
 (
   struct reflectance const * refl,
   double const * s_scat_i,
   double const * s_scat_r,
   double * dens
 );

 void
 refl_dhr_minn
 (
   struct reflectance const * refl,
   double const * s_scat_i,
   double * dhr
 );

 void
 refl_brdf_minn
 (
   struct reflectance const * refl,
   double const * s_scat_i,
   double const * s_scat_r,
   double * bdr
 );


 // Lambert:

 struct reflectance_lamb *
 refl_alloc_lamb
 ( void );

 void 
 refl_free_lamb
 ( struct reflectance_lamb ** lamb );

 void
 refl_setup_lamb
 (
   struct reflectance * refl,
   int const nbhr, 
   double const * bhr
 );

 void
 refl_qtl_lamb
 (
   struct reflectance const * refl, 
   double const * s_scat_i,
   double * s_scat_r,
   double const * prob
 );

 void
 refl_cdf_lamb
 (
   struct reflectance const * refl,
   double const * s_scat_i,
   double const * s_scat_r,
   double * prob
 );

 void
 refl_pdf_lamb
 (
   struct reflectance const * refl,
   double const * s_scat_i,
   double const * s_scat_r, 
   double * dens
 );

 void 
 refl_dhr_lamb
 (
   struct reflectance const * refl,
   double const * s_scat_i,
   double * dhr
 );

 void
 refl_brdf_lamb
 (
   struct reflectance const * refl, 
   double const * s_scat_i,
   double const * s_scat_r,
   double * bdr
 );

/* Unidirectional model */
/*
 struct brdf_par_unid
 {
   int    nbhr;		// Number of bottom reflectances
   double *p_bhr;	// Bi-hemispherical reflectance
   double *p_dhr;	// Directional-hemispherical reflectance
   double *p_brdf; 	// Delta function BRDF
   double *p_brdf0; 	// Delta eval to 0
   double mutf;		// Reference fixed cosine of polar angle
   double phif;		// Reference fixed cosine of azimuth angle
 };

 void brdf_build_unid (void  **pp_brdf_par,
                       const int    nbhr,
                       const double *p_bhr,
                       const double thetaf,
                       const double phif);

 double* brdf_unid (const double mut_i,
                    const double phi_i,
                    const double mut_r,
                    const double phi_r,
                    void   *p_brdf_par);

 double* brdf_dhr_unid (const double mut_i,
                        const double phi_i,
                        void   *p_brdf_par);

 double brdf_cdf_phi_unid (const double mut_i,
                           const double phi_i,
                           const double mut_r,
                           const double phi_r,
                           const void   *p_brdf_par);

 double brdf_cdf_psi_unid (const double mut_i,
                           const double phi_i,
                           const double mut_r,
                           const double phi_r,
                           const void   *p_brdf_par);

 void brdf_qtl_unid (const double mut_i,
                     const double phi_i,
                     double       *p_mut_r,
                     double       *p_phi_r,
                     const void   *p_brdf_par,
                     const double prob_psi,
                     const double prob_phi);
*/
 #endif // REFLECTANCE

