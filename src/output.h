
 #ifndef OUTPUT
 #define OUTPUT


 void write_vec (struct str_acm *acm,
                 const int    sim_ns,
                 const char   *ofbn,
                 const char   *sufx,
                 const int    ST, //STRMXLEN,
                 const double *iop_w0,
                 const double *btt_bh,
                 const double *sim_sza);

 void write_grid (struct str_acm *acm,
                  int    NT,
                  char   *ofbn,
                  char   *sufx,
                  int    ST, //STRMXLEN,
                  double *iop_w0,
                  double *btthdr,
                  double *sun);

 void write_acm (struct str_acm *acm,
                  int    NT,
                  char   *ofbn,
                  char   *sufx,
                  int    ST, //STRMXLEN,
                  double *iop_w0,
                  double *btthdr,
                  double *sun);

 void print_sim (const char   *ofbn,
                 const char   *version,
                 const int    sim_nr,
                 const int    sim_ns,
                 const int    iop_nw0,
                 const int    btt_nr,
                 const int    str_ncl,
                 const int    str_ncn,
                 const int    str_nbx,
                 const int    ST, //STRMXLEN,
                 const char   *sns_tp,
                 const double sns_fov,
                 const double *sns_pos,
                 const double *sim_sza,
                 const double sim_f0,
                 const double iop_na,
                 const double iop_nw,
                 const double iop_c,
                 const double *iop_w0,
                 const double btt_d,
                 const char   *btt_tp,
                 const double *btt_bh,
                 const char   *scat_tp,
                 const double scat_g,
                 const double scat_fbb,
                 const int    acc_geom,
                 const int    acc_fgeom,
                 const double acc_ext,
                 const double acc_resx,
                 const double acc_resy,
                 const struct str_cylnd *str_cls,
                 const struct str_cone *str_cns,
                 double **str_bxs,
                 const struct str_acm *acm);

 void print_cls (struct str_cylnd *cls, 
                 int    NCL);

 void print_cns (struct str_cone *cns, 
                 int    NCN);

 #endif

