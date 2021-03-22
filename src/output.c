
/*
 * Two functions are available to write the integrals to disk, one for each 
 * component of the accumulators.
 *
 * To save in the number of files, the function writing the total spatial 
 * integral writes down a matrix where each line represents an illumination and
 * each column represents a single scattering albedo. Different files are 
 * created per bottom relfectance. The first row and the first column are row 
 * and column names.
 * 
 * For the spatially resolved integrals, one file is written for each 
 * illumination, single scattering albedo and bottom reflectance. Rows and 
 * column names are not provided. The simulation printout file will provide the
 * actual extent and number of cells for the grid.
 *
 */

 #include <stdio.h>	// Functions:   fprintf, printf,
 #include <string.h>	// Functions:   strcat, strcpy, sprintf
 #include <stdlib.h>	// Functions:   exit
 #include <math.h>      // Functions:   ceil

 #include "config.h"    // Definitions: SPATIALLY_RESOLVED, SHADOWING
 #include "constants.h"     // Macros:      DEG
 #include "strc.h"	// Structures:  str_acm


 void write_vec (struct str_acm *acm,	// Pointer to accumulator 
                 const int    sim_ns,	// Number of Sun zenith angles
                 const char   *ofbn,	// Base file name
                 const char   *sufx,	// Suffix
                 const int    ST, //STRMXLEN,	// Maximum string size
                 const double *iop_w0,	// Pointer to single scattering albedos
                 const double *btt_bh,	// Pointer to bottom reflectance
                 const double *sim_sza)	// Pointer to Sun zenith angles

 {
   int  cw, cb, cz;
   char wf[10], bf[10];
   char ofn[200];
   FILE *fpo;

   for (cb = 0; cb < acm->nbr; cb++)
   {
     sprintf(bf, "_B%05.3f", btt_bh[cb]);
     strncpy(ofn, ofbn, STRMXLEN);
     strcat(ofn, sufx);
     strcat(ofn, bf);
     strcat(ofn, ".txt");

     fpo = fopen (ofn, "w");
     if (fpo == NULL)
     {
       printf("ERROR: Failed to create output file %s\n", ofn);
       exit(-3);
     }

     fprintf(fpo, "\t");
     for (cw = 0; cw < (acm->nw0 - 1); cw++)
     {
       fprintf(fpo, "W_%6.4f ", iop_w0[cw]);
     }
     fprintf(fpo, "W_%6.4f\n", iop_w0[cw]);
     for (cz = 0; cz < acm->ns; cz++)
     {
       if (cz < sim_ns)
       {
         fprintf(fpo, "I_%04.1f\t", sim_sza[cz] * DEG);
       } else if (cz == sim_ns) {
         fprintf(fpo, "UNIF\t");
       } else {
         fprintf(fpo, "CARD\t");
       }
       for (cw = 0; cw < (acm->nw0 - 1); cw++)
       {
         fprintf(fpo, "%.6E ", acm->ppp_vector[cz][cw][cb]);
       }
       fprintf(fpo, "%.6E\n", acm->ppp_vector[cz][cw][cb]);
     }
     fclose (fpo);
   }
 }

 void write_grid (struct str_acm *acm,	// Pointer to accumulator
                  int    sim_ns,		// Number of Sun zenith angles
                  char   *ofbn,		// Base file name
                  char   *sufx,		// Suffix
                  int    ST, //STRMXLEN,	// Maximum string size
                  double *iop_w0,	// Pointer to single scattering albedos
                  double *btt_bh,	// Pointer to bottom reflectance
                  double *sim_sza)		// Pointer to Sun zenith angles
 {
   int cw, cb, cz, cr, cc;
   char wf[12], bf[12], sz[12];
   char ofn[200];
   FILE *fpo;

   for (cz = 0; cz < acm->ns; cz++)
   {
     if (cz < sim_ns)
     {
       sprintf(sz, "_S%03d", (int) ceil(sim_sza[cz] * DEG));
     } else if (cz == sim_ns) {
       sprintf(sz, "_Isotropic");
     } else {
       sprintf(sz, "_Cardioidal");
     }
     for (cw = 0; cw < acm->nw0; cw++)
     {
       sprintf(wf, "_W%05.2f", iop_w0[cw]);
       for (cb = 0; cb < acm->nbr; cb++)
       {
         sprintf(bf, "_B%05.3f", btt_bh[cb]);
         strncpy(ofn, ofbn, STRMXLEN);
         strcat(ofn, sufx);
         strcat(ofn, sz);
         strcat(ofn, wf);
         strcat(ofn, bf);
         strcat(ofn, ".txt");

         fpo = fopen (ofn, "w");
         if (fpo == NULL)
         {
           printf("ERROR: Failed to create output file %s\n", ofn);
           exit(-3);
         }
         for (cr = 0; cr < acm->ny; cr++)
         {
           for (cc = 0; cc < (acm->nx - 1); cc++)
           {
             fprintf(fpo, "%.6E ", acm->ppppp_grid[cz][cw][cb][cr][cc]);
           }
           fprintf(fpo, "%.6E\n", acm->ppppp_grid[cz][cw][cb][cr][cc]);
         }
         fclose (fpo);
       }
     }
   }
 }

 void write_acm (struct str_acm *acm,	// Pointer to accumulator
                  int    sim_ns,		// Number of Sun zenith angles
                  char   *ofbn,		// Base file name
                  char   *sufx,		// Suffix
                  int    ST, //STRMXLEN,	// Maximum string size
                  double *iop_w0,	// Pointer to single scattering albedos
                  double *btt_bh,	// Pointer to bottom reflectance
                  double *sim_sza)		// Pointer to Sun zenith angles
 {
   char *sufx_1, *sufx_2;

   sufx_1 = (char*) calloc(STRMXLEN, sizeof(char));
   sufx_2 = (char*) calloc(STRMXLEN, sizeof(char));

   strncpy(sufx_1, sufx, STRMXLEN);
   strcat(sufx_1, "_spi");
   write_vec(acm, sim_ns, ofbn, sufx_1, STRMXLEN, iop_w0, btt_bh, sim_sza);

   strncpy(sufx_2, sufx, STRMXLEN);
   strcat(sufx_2, "_spr");
   write_grid(acm, sim_ns, ofbn, sufx_2, STRMXLEN, iop_w0, btt_bh, sim_sza);

   free(sufx_1);
   free(sufx_2);
 }

 void print_cls (struct str_cylnd *str_cls, int str_ncl)
 {
   int ci;

   printf("  Number of cylinders: %d\n", str_ncl);
   printf("                               base                            axis "
     "     height   closed     radius\n");
   for(ci = 0; ci < str_ncl; ci++)
   {
     printf("  %02d  %9.2e,%9.2e,%9.2e   %9.2e,%9.2e,%9.2e   %9.2e      "
       "%1d,%1d  %9.2e\n",
       ci+1,
       str_cls[ci].base[0], 
       str_cls[ci].base[1], 
       str_cls[ci].base[2], 
       str_cls[ci].caxis[0], 
       str_cls[ci].caxis[1], 
       str_cls[ci].caxis[2], 
       str_cls[ci].height[1],
       str_cls[ci].closed[0], 
       str_cls[ci].closed[1], 
       str_cls[ci].radius);
   }
   printf("\n");
 }

 void print_cns (struct str_cone *str_cns, int str_ncn)
 {
   int ci;

   printf("  Number of cones: %d\n", str_ncn);
   printf("                             vertex                            axis "
     "              height  angle closed               radius\n");
   for(ci = 0; ci < str_ncn; ci++)
   {
     printf("  %02d  %9.2e,%9.2e,%9.2e   %9.2e,%9.2e,%9.2e  %9.2e,%9.2e   "
       "%3.1f    %1d,%1d  %9.2e,%9.2e \n", 
       ci+1,
       str_cns[ci].vertex[0], 
       str_cns[ci].vertex[1], 
       str_cns[ci].vertex[2], 
       str_cns[ci].caxis[0], 
       str_cns[ci].caxis[1], 
       str_cns[ci].caxis[2], 
       str_cns[ci].height[0], 
       str_cns[ci].height[1], 
       acos(str_cns[ci].mu) * DEG, 
       str_cns[ci].closed[0], 
       str_cns[ci].closed[1],
       str_cns[ci].radius[0], 
       str_cns[ci].radius[1]);
   }
   printf("\n");
 }

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
                 const struct str_acm *acm)
 {
   int  ci, cj;
   char ofn[200];
   FILE *fpo;

   strncpy(ofn, ofbn, STRMXLEN);
   strcat(ofn, "_out_summary.txt");

   fpo = fopen (ofn, "w");
   if (fpo == NULL)
   {
     printf("ERROR: Failed to create output file %s\n", ofn);
     exit(-3);
   }

   fprintf (fpo, "\nBackward Monte Carlo radiative transfer code\n"); 
   fprintf (fpo, "Version: %s\n", version);
   fprintf (fpo, "Author: Alexandre Castagna\n\n");

   fprintf (fpo, "System properties:\n");
   fprintf (fpo, "  Rays:    %.2e\n", (double) sim_nr);
   if (strcmp (sns_tp, "lu") == 0)
   {
     fprintf (fpo, "  Source:  %s with %.2fº FOV at z = %.2e m\n", sns_tp, 
       sns_fov * DEG, sns_pos[2]);
   } else {
     fprintf (fpo, "  Source:  %s at z = %.2f m\n", sns_tp, sns_pos[2]);
   }
   fprintf (fpo, "  SZA:     (%d): ", sim_ns);
   for (ci = 0; ci < (sim_ns - 1); ci++)
   {
     fprintf(fpo, "%.1f, ", DEG * sim_sza[ci]);
   }
   fprintf(fpo, "%.1f (º)\n", DEG * sim_sza[ci]);
   fprintf (fpo, "  f0:      %.2f (W/m2)\n", sim_f0);
   fprintf (fpo, "  n_air:   %.2f (unitless)\n", iop_na);
   fprintf (fpo, "  n_water: %.2f (unitless)\n", iop_nw);
   fprintf (fpo, "  c coef.: %.2f (1/m)\n", iop_c);
   fprintf (fpo, "  w0:      (%d): ", iop_nw0);
   for (ci = 0; ci < (iop_nw0 - 1); ci++)
   {
     fprintf(fpo, "%.3f, ", iop_w0[ci]);
   }
   fprintf(fpo, "%.3f (unitless)\n", iop_w0[ci]);
   if (strcmp (scat_tp, "isotropic") == 0) 
   {
     fprintf(fpo, "  Scat.:   Isotropic (Henyey-Greenstein with g = 0)\n");
   }
   if (strcmp (scat_tp, "hg") == 0) 
   {
     fprintf(fpo, "  Scat.:   Henyey-Greenstein, with asymmetry parameter of "
       "%.3e (unitless)\n", scat_g);
   }
   if (strcmp (scat_tp, "ff") == 0) 
   {
     fprintf(fpo, "  Scat.:   Fournier-Forand, with fbb of %.3e (unitless)\n", 
       scat_fbb);
   }
   if (strcmp (scat_tp, "petzold") == 0) 
   {
     fprintf(fpo, "  Scat.:   Average Petzold (fbb of 0.0183, unitless)\n");
   }

   if ( isfinite(btt_d) )
   {
     if ( strcmp (btt_tp, "unidir") == 0 ) 
     {
       fprintf(fpo, "  Bottom:  Unidirectional, at %.3e m, and with %d BHR: ", 
         btt_d, btt_nr);
     }
     if ( strcmp (btt_tp, "minnaert") == 0 ) 
     {
       fprintf(fpo, "  Bottom:  Minnaert, at %.3e m, and with %d BHR: ", btt_d, 
         btt_nr);
     }
     if ( strcmp (btt_tp, "lambert") == 0 ) 
     {
       fprintf(fpo, "  Bottom:  Lambert (Minnaert with k = 0), at %.3e m, and "
        "with %d BHR: ", btt_d, btt_nr);
     }
     for (ci = 0; ci < (btt_nr - 1); ci++)
     {
       fprintf(fpo, "%.3f, ", btt_bh[ci]);
     }
     fprintf(fpo, "%.3f (unitless)\n", btt_bh[ci]);
   } else {
     fprintf(fpo,"  Bottom:  Undefined, at INFINITE depth\n");
   } // Logic: btt_d

   #ifdef SPATIALLY_RESOLVED
   if (!acc_fgeom)
   {
     fprintf (fpo, "\nSpatially resolved integral:\n");
     if (acc_geom == 1)
     {
       fprintf (fpo, "  Geometry: Sectorial\n");
     } else {
       fprintf (fpo, "  Geometry: Grid\n");
     }
     fprintf (fpo, "  Extent: %.3e (m)\n", acc_ext);
     if (acc_geom == 1)
     {
       fprintf (fpo, "  Resolution (x): %.3e (º)\n", acc_resx);
     } else {
       fprintf (fpo, "  Resolution (x): %.3e (m)\n", acc_resx);
     }
     fprintf (fpo, "  Resolution (y): %.3e (m)\n", acc_resy);
     fprintf (fpo, "  (grid breakpoints in x and y in file 'out_grid_breakpoints.txt')\n", acc_resy);
   } else {
     fprintf(fpo, "\nSpatially resolved integral:\n");
     fprintf(fpo, "  Code compiled with spatially resolved integral, but acc_geom was\n"
            "  set to 0. For a simulation without spatially resolved integral,\n"
            "  it will give a gain in performance to comment out the the \n"
            "  definition SPATIALLY_RESOLVED in config.h\n");
   }
   #endif // SPATIALLY_RESOLVED

   #ifdef SHADOWING
   fprintf(fpo, "\nStructures:\n");
   if (str_ncl > 0)
   {
     fprintf(fpo, "  Number of cylinders: %d\n", str_ncl);
     fprintf(fpo, "                               base                            axis "
       "     height   closed     radius\n");
     for(ci = 0; ci < str_ncl; ci++)
     {
       fprintf(fpo, "  %02d  %9.2e,%9.2e,%9.2e   %9.2e,%9.2e,%9.2e   %9.2e      "
         "%1d,%1d  %9.2e\n",
         ci+1,
         str_cls[ci].base[0], 
         str_cls[ci].base[1], 
         str_cls[ci].base[2], 
         str_cls[ci].caxis[0], 
         str_cls[ci].caxis[1], 
         str_cls[ci].caxis[2], 
         str_cls[ci].height[1],
         str_cls[ci].closed[0], 
         str_cls[ci].closed[1], 
         str_cls[ci].radius);
     }
     printf("\n");
   }
   if (str_ncn > 0)
   {
     fprintf(fpo, "  Number of cones: %d\n", str_ncn);
     fprintf(fpo, "                             vertex                            axis "
       "              height  angle closed               radius\n");
     for(ci = 0; ci < str_ncn; ci++)
     {
       fprintf(fpo, "  %02d  %9.2e,%9.2e,%9.2e   %9.2e,%9.2e,%9.2e  %9.2e,%9.2e   "
         "%3.1f    %1d,%1d  %9.2e,%9.2e \n", 
         ci+1,
         str_cns[ci].vertex[0], 
         str_cns[ci].vertex[1], 
         str_cns[ci].vertex[2], 
         str_cns[ci].caxis[0], 
         str_cns[ci].caxis[1], 
         str_cns[ci].caxis[2], 
         str_cns[ci].height[0], 
         str_cns[ci].height[1], 
         acos(str_cns[ci].mu) * DEG, 
         str_cns[ci].closed[0], 
         str_cns[ci].closed[1],
         str_cns[ci].radius[0], 
         str_cns[ci].radius[1]);
     }
     printf("\n");
   }
   if (str_nbx > 0)
   {
     if (str_nbx > 2) {
       fprintf(fpo, "  Number of axis-aligned boxes: %d\n", str_nbx - 1);
       cj = 1;
     } else {
       fprintf(fpo, "  Number of axis-aligned boxes: %d\n", str_nbx);
       cj = 0;
     }
     fprintf(fpo, "     xmn       xmx       ymn       ymx       zmn       zmx\n");
     for(ci = cj; ci < str_nbx; ci++)
     {
       fprintf(fpo, "  %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e\n", str_bxs[ci][0], 
         str_bxs[ci][1], str_bxs[ci][2], str_bxs[ci][3], str_bxs[ci][4], str_bxs[ci][5]);
     } 
   }
   if ( str_ncl == 0 && str_ncn == 0 && str_nbx == 0 ) 
   {
     fprintf(fpo, "  Code compiled with shadow tracking, but no structure was\n"
            "  defined. For a simulation without shadowing, it will give\n"
            "  a gain in performance to comment out the the definition\n"
            "  SDHW in config.h\n");
   }
   #endif // SHADOWING

   #ifdef AIRTIGHT_ALEX
   fprintf(fpo,"\nCAUTION! Code compiled with AIRTIGHT_ALEX defined in config.h - this results"
          "in a specific circunstances for Alex's system and is not general, "
          "even without shadowing! Comment out AIRTIGHT_ALEX in config.h for generic "
          "simulation.\n");
   #endif // AIRTIGHT_ALEX

   fclose (fpo);

   #ifdef SPATIALLY_RESOLVED
   if (acc_geom > 0 && acc_fgeom == 0)
   {
     strncpy(ofn, ofbn, STRMXLEN);
     strcat(ofn, "_out_grid_breakpoints.txt");

     fpo = fopen (ofn, "w");
     if (fpo == NULL)
     {
       printf("ERROR: Failed to create output file %s\n", ofn);
       exit(-3);
     }

     fprintf(fpo, "y ");
     for (ci = 0; ci < acm->ny; ci++)
     {
       fprintf(fpo, "%.3f, ", acm->p_ybrks[ci]);
     }
     fprintf(fpo, "%.3f\n", acm->p_ybrks[ci]);

     fprintf(fpo, "x ");
     for (ci = 0; ci < acm->nx; ci++)
     {
       fprintf(fpo, "%.3f, ", acm->p_xbrks[ci]);
     }
     fprintf(fpo, "%.3f\n", acm->p_xbrks[ci]);
     fclose (fpo);
   }
   #endif // SPATIALLY_RESOLVED
 }

