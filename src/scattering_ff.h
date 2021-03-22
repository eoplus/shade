
 #ifndef SCATTERING_FF
 #define SCATTERING_FF

 #include "scattering.h"

/* Function prototypes: *******************************************************/

 void scat_setup_ff (struct scattering *scat, double const scat_fbb,
   int const truncate, char const *method);

 void scat_findpar_ff (double const scat_fbb, double *j, double *iop_nh);


 void scat_cdf_ff (struct scattering const *scat, double const *s_scat,
   double *prob);

 void scat_pdf_ff (struct scattering const *scat, double const *s_scat,
   double *dens);

 void scat_pf_ff (struct scattering const *scat, double const *s_scat,
   double *pf);


 void scat_qtl_ff_itp (struct scattering const *scat, double *s_scat,
   double const *prob);

 void scat_cdf_ff_itp (struct scattering const *scat, double const *s_scat,
   double *prob);

 void scat_pdf_ff_itp (struct scattering const *scat, double const *s_scat,
   double *dens);

 void scat_pf_ff_itp (struct scattering const *scat, double const *s_scat,
   double *pf);


 void scat_qtl_ff_lut (struct scattering const *scat, double *s_scat,
   double const *prob);

 void scat_cdf_ff_lut (struct scattering const *scat, double const *s_scat,
   double *prob);

 void scat_pdf_ff_lut (struct scattering const *scat, double const *s_scat,
   double *dens);

 void scat_pf_ff_lut (struct scattering const *scat, double const *s_scat,
   double *pf);


 void scat_qtl_ff_min (struct scattering const *scat, double *s_scat,
   double const *prob);

 void scat_cdf_ff_min (struct scattering const *scat, double const *s_scat,
   double *prob);

 void scat_pdf_ff_min (struct scattering const *scat, double const *s_scat,
   double *dens);

 void scat_pf_ff_min (struct scattering const *scat, double const *s_scat,
   double *pf);

 #endif // SCATTERING_FF

