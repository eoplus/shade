
 #include <stdio.h>
 #include <stddef.h>
 #include <string.h>

 #include "config.h"
 #include "memory.h"
 #include "constants.h"
 #include "structures.h"
 #include "geometry.h"
 #include "rotation.h"
 #include "sources.h"
 #include "scattering.h"
 #include "bottom.h"

 void
 input_read
 (
   const char *flnm,
   int      *sim_nr,
   int      *sim_ns,
   double   **sim_sza,
   double   **sim_saa,
   double   *sim_f0,
   double   *skr_resx,
   double   *skr_resy,
   int      *skr_nx,
   int      *skr_ny,
   char     ***skr_fls,
   double   *iop_na,
   double   *iop_nw,
   double   *iop_c,
   int      *iop_nw0,
   double   **iop_w0,
   char     *scat_tp,
   int      *scat_trc,
   char     *scat_mtd,
   double   *scat_d,
   double   *scat_g,
   double   *scat_fbb,
   double   *btt_d,
   char     *btt_tp,
   double   *btt_par_thetaf,
   double   *btt_par_phif, 
   double   *btt_par_rho,
   double   *btt_par_k,
   int      *btt_nbr,
   double   **btt_bhr,
   char     *src_tp,
   double   *src_fov,
   double   *src_saxs,
   double   *src_ref_o,
   double   *src_rel_o,
   int      *acc_fgeom,
   char     *acc_geom,
   double   *acc_ext,
   double   *acc_resx,
   double   *acc_resy,
   int      *str_def,
   int      *str_ncl,
   str_cyln ***cylns,
   int      *str_ncn,
   str_cone ***cones,
   int      *str_ncb,
   str_cubd ***cubds,
   int      *expand_str_ncb,
   int const cnfg_spres,
   int const cnfg_shdw
 )
 {
   FILE* fip;			// Pointer for input file
   long int fppos;		// Store the current position in the file stream
   char varnm[STRMXLEN];	// Testing variables names
   char hcheck[STRMXLEN];	// String to test presence of header line

   // Temporary variables to read each structure:
   double origin[3];
   double axis[3];
   double alpha;
   double lengths[3];
   double radius;
   double height;
   double s_base[3] = {1.0};
   double s_top[3] = {1.0};
   int    closed[2];

   int ci;

   fip = fopen(flnm, "r");
   if(fip == NULL)
   {
     printf("Failed to open input file %s\n", flnm);
     exit(-1);
   }
   fppos = ftell(fip);

   // Skip input file header lines:
   while (fgets (hcheck, 200, fip) != NULL && hcheck[0] == ';')
   {
     fppos = ftell(fip);
   }
   fseek(fip, fppos, SEEK_SET);
 
   // Number of rays to trace:
   fscanf(fip, "%s %d", varnm, sim_nr);
   if ( strcmp (varnm, "sim_nr:") )
   {
     printf ("\nERROR: Expected 'sim_nr:' (number of rays to trace) from input,"
       " but got '%s'\n", varnm);
     exit (-1);    
   }
   if((*sim_nr) < SUBDV)
   {
     printf ("\nERROR: sim_nr (number of rays to trace) must be at least %d\n", 
       SUBDV);
     exit (-1);
   }

   // Number of sun zenith angles:
   fscanf(fip, "%s %d", varnm, sim_ns);
   if ( strcmp (varnm, "sim_ns:") )
   {
     printf ("\nERROR: Expected 'sim_ns:' (number of Sun positions) from input,"
       " but got '%s'\n", varnm);
     exit (-1);    
   }
   if((*sim_ns) <= 0)
   {
     printf ("\nERROR: sim_ns (number of Sun positions) must be > 0\n");
     exit (-1);
   }
   *sim_sza = calloc_1d ((*sim_ns), "sim_sza in input_read (read_input.c)");
   *sim_saa = calloc_1d ((*sim_ns), "sim_saa in input_read (read_input.c)");
   fscanf(fip, "%s", varnm);
   if ( strcmp (varnm, "sim_sza(s):") )
   {
     printf ("\nERROR: Expected 'sim_sza(s):' (Sun zenith angles) from input, "
       "but got '%s'\n", varnm);
     exit (-1);    
   }
   for(size_t cs = 0; cs < (*sim_ns); cs++)
   {
     fscanf(fip, "%lf", &(*sim_sza)[cs]);
     (*sim_sza)[cs] *= RAD;
   }
   if ( (*sim_sza)[ci] < 0.0 || (*sim_sza)[ci] > M_PI_2 )
   {
     printf ("\nERROR: sim_sza (Sun zenith angles) must be in the range "
       "[0º,90º]\n");
     exit (-1);
   }
   fscanf(fip, "%s", varnm);
   if ( strcmp (varnm, "sim_saa(s):") )
   {
     printf ("\nERROR: Expected 'sim_saa(s):' (Sun azimuth angles) from input, "
       "but got '%s'\n", varnm);
     exit (-1);    
   }
   for(size_t cs = 0; cs < (*sim_ns); cs++)
   {
     fscanf(fip, "%lf", &(*sim_saa)[cs]);
     (*sim_saa)[ci] *= RAD;
   }
   if ( (*sim_saa)[ci] < 0.0 || (*sim_saa)[ci] > K_2PI )
   {
     printf ("\nERROR: sim_saa(s) (Sun azimuth angles) must be in the range "
       "[0º,360º]\n");
     exit (-1);
   }

   // Incident irradiance:
   fscanf(fip, "%s %lf", varnm, sim_f0);
   if ( strcmp (varnm, "sim_f0:") )
   {
     printf ("\nERROR: Expected 'sim_f0:' (incident irradiance) from input, but"
       " got '%s'\n", varnm);
     exit (-1);    
   }
   if ( (*sim_f0) <= 0.0 )
   {
     printf("ERROR: sim_f0 (incident irradiance) must be in the > 0.0\n");
     exit(-1);
   }

   // Sky radiance distribution parameters:
   fscanf(fip, "%s %lf", varnm, skr_resx);
   if ( strcmp (varnm, "skr_resx:") )
   {
     printf ("\nERROR: Expected 'skr_resx:' (azimuthal angle resolution of the "
       "sky radiance data) from input, but got '%s'\n", varnm);
     exit (-1);    
   }
   if ( (*skr_resx) < 0.0 || (*skr_resx) > 360.0 || 
        fmod(360.0, (*skr_resx)) != 0.0 )
   {
     printf("ERROR: skr_resx (azimuthal angle resolution of the sky radiance "
       "data) must be in the > 0.0 and < 360 degrees. It also needs to divide "
       "360 in equal and complete parts.\n");
     exit(-1);
   }
   *skr_nx = (int) (360.0 / (*skr_resx));
   *skr_resx *= RAD;

   fscanf(fip, "%s %lf", varnm, skr_resy);
   if ( strcmp (varnm, "skr_resy:") )
   {
     printf ("\nERROR: Expected 'skr_resy:' (polar angle resolution of the sky "
       "radiance data) from input, but got '%s'\n", varnm);
     exit (-1);    
   }
   if ( (*skr_resy) < 0.0 || (*skr_resy) > 90.0 || 
        fmod(90.0, (*skr_resy)) != 0.0 )
   {
     printf("ERROR: skr_resx (polar angle resolution of the sky radiance data) "
       "must be in the > 0.0 and < 90 degrees. It also needs to divide 90 in "
       "equal and complete parts.\n");
     exit(-1);
   }
   *skr_ny = (int) (90.0 / (*skr_resy));
   *skr_resy *= RAD;

   fscanf(fip, "%s", varnm);
   if ( strcmp (varnm, "skr_fls:") )
   {
     printf ("\nERROR: Expected 'skr_fls:' (Sky radiance file path) from input,"
       "but got '%s'\n", varnm);
     exit (-1);    
   }
   *skr_fls = (char**) calloc((*sim_ns) + 2, sizeof(char*));
   for (size_t cs = 0; cs < (*sim_ns); cs++)
   {
     (*skr_fls)[cs] = (char*) calloc(STRMXLEN, sizeof(char));
     fscanf(fip, "%s", (*skr_fls)[cs]);
   }
   (*skr_fls)[(*sim_ns)] = (char*) calloc(STRMXLEN, sizeof(char));
   strncpy((*skr_fls)[(*sim_ns)], "anc/skyrad/isotropic_sky.txt", STRMXLEN);
   (*skr_fls)[(*sim_ns) + 1] = (char*) calloc(STRMXLEN, sizeof(char));
   strncpy((*skr_fls)[(*sim_ns) + 1], "anc/skyrad/cardioidal_sky.txt", STRMXLEN);

   // Inherent optical properties: *********************************************
   fscanf(fip, "%s %lf", varnm, iop_na);
   if ( strcmp (varnm, "iop_na:") )
   {
     printf ("\nERROR: Expected 'iop_na:' (refractive index of air) from input,"
       "but got '%s'\n", varnm);
     exit (-1);    
   }
   if ( (*iop_na) < 1.0 )
   {
     printf("ERROR: iop_na (refractive index of air) must be in the >= 1.0\n");
     exit(-1);
   }

   fscanf(fip, "%s %lf", varnm, iop_nw); 
   if ( strcmp (varnm, "iop_nw:") )
   {
     printf ("\nERROR: Expected 'iop_nw:' (refractive index of water) from "
       "input, but got '%s'\n", varnm);
     exit (-1);    
   }
   if ( (*iop_nw) < 1.0 )
   {
     printf("ERROR: iop_nw (refractive index of water) must be in the >= 1.0\n");
     exit(-1);
   }
  
   fscanf(fip, "%s %lf", varnm, iop_c);
   if ( strcmp (varnm, "iop_c:") )
   {
     printf ("\nERROR: Expected 'iop_c:' (beam attenuation coefficient) from "
       "input, but got '%s'\n", varnm);
     exit (-1);    
   }
   if ( (*iop_c) <= 0.0 )
   {
     printf("ERROR: iop_c (beam attenuation coefficient) must be > 0.0\n");
     exit(-1);
   }

   fscanf(fip, "%s %d", varnm, iop_nw0);
   if ( strcmp (varnm, "iop_nw0:") )
   {
     printf ("\nERROR: Expected 'iop_nw0:' from input, but got '%s'\n", varnm);
     exit (-1);    
   }
   if ( (*iop_nw0) <= 0 )
   {
     printf("ERROR: iop_nw0 (number of single scattering albedos) must be > "
       "0\n");
     exit(-1);
   }
   *iop_w0 = calloc_1d ((*iop_nw0), "iop_w0 in input_read (read_input.c)");
   fscanf(fip, "%s", varnm);
   if ( strcmp (varnm, "iop_w0:") )
   {
     printf ("\nERROR: Expected 'iop_w0:' from input, but got '%s'\n", varnm);
     exit (-1);    
   }
   for(ci = 0; ci < (*iop_nw0); ci++)
   {
     fscanf(fip, "%lf", &(*iop_w0)[ci]);
     if ( ((*iop_w0)[ci] < 0.0) || ((*iop_w0)[ci] > 1.0) )
     {
       printf("ERROR: iop_w0 (single scattering albedo) must be in the range "
         "[0,1]\n");
       exit(-1);
     }
   }

   fscanf(fip, "%s %s", varnm, scat_tp);
   if ( strcmp (varnm, "scat_tp:") )
   {
     printf ("\nERROR: Expected 'scat_tp:' (scattering phase function type) "
       "from input, but got '%s'\n", varnm);
     exit (-1);    
   }
   if ( strcmp (scat_tp, "hg") != 0 && 
        strcmp (scat_tp, "ff") != 0 &&
        strcmp (scat_tp, "rayleigh") != 0 &&
        strcmp (scat_tp, "petzold") != 0 && 
        strcmp (scat_tp, "isotropic") != 0)
   {
     printf ("\nERROR: scat_tp (scattering phase function type) must be one of:"
       " 'isotropic', 'rayleigh', 'hg', 'ff', 'petzold'\n");
     exit (-1);
   }

   if ( strcmp (scat_tp, "rayleigh") == 0 )
   {
     fscanf(fip, "%s %lf", varnm, scat_d);
     if ( strcmp (varnm, "scat_d:") )
     {
       printf ("\nERROR: Expected 'scat_d:' (linear depolarization ratio) from "
         "input, but got '%s'\n", varnm);
       exit (-1);    
     }
     if ( ((*scat_d) < 0.0) || ((*scat_d) > 1.0) )
     {
       printf("ERROR: scat_d (depolarization ratio) must be in the range "
         "[0,1]\n");
       exit(-1);
     }
   }

   if ( strcmp (scat_tp, "hg") == 0 )
   {
     fscanf(fip, "%s %lf", varnm, scat_g);
     if ( strcmp (varnm, "scat_g:") )
     {
       printf ("\nERROR: Expected 'scat_g:' (average scattering angle cosine) "
         "from input, but got '%s'\n", varnm);
       exit (-1);    
     }
     if ( ((*scat_g) < 0.0) || ((*scat_g) > 1.0) )
     {
       printf("ERROR: scat_g (assymetry parameter for the Henyey-Greestein "
         "phase function) must be in the range [0,1]\n");
       exit(-1);
     }
     // Current implementation of hg will not work for scat_g = 0, but this is 
     // the same as isotropic scattering:
     if (strcmp (scat_tp, "hg") == 0 && (*scat_g) == 0)
     {
       strcpy (scat_tp, "isotropic");
     }
   }

   if ( strcmp (scat_tp, "ff") == 0 )
   {
     fscanf(fip, "%s %d", varnm, scat_trc);
     if ( strcmp (varnm, "scat_trc:") )
     {
       printf ("\nERROR: Expected 'scat_trc:' (logical flag: should the "
         "diffraction peak of the phase function be truncated?) from input, but"
         " got '%s'\n", varnm);
       exit (-1);    
     }
     fscanf(fip, "%s %s", varnm, scat_mtd);
     if ( strcmp (varnm, "scat_mtd:") )
     {
       printf ("\nERROR: Expected 'scat_mtd:' (method used to retrieve "
         "scattering values) from input, but got '%s'\n", varnm);
       exit (-1);
     }
     if ( strcmp (scat_mtd, "lut") != 0 && 
          strcmp (scat_mtd, "itp") != 0 &&
          strcmp (scat_mtd, "min") != 0 )
     {
       printf ("\nERROR: scat_mtd (method used to retrieve scattering values) "
         "must be one of: 'lut' (Look-up table direct search), 'itp' (Look-up "
         "table interpolation), 'min' (function minimization)\n");
       exit (-1);
     }
     fscanf(fip, "%s %lf", varnm, scat_fbb);
     if ( strcmp (varnm, "scat_fbb:") )
     {
       printf ("\nERROR: Expected 'scat_fbb:' (backscattering fraction for the "
         "Fournier-Forand phase function) from input, but got '%s'\n", varnm);
       exit (-1);    
     }
     if ( ((*scat_fbb) < 0.0) || ((*scat_fbb) > 1.0) )
     {
       printf("ERROR: scat_fbb (backscattering fraction for the Fournier-Forand"
         " phase function) must be in the range [0,1]\n");
       exit(-1);
     }
   }

   fscanf(fip, "%s %lf", varnm, btt_d);
   if ( strcmp (varnm, "btt_d:") )
   {
     printf ("\nERROR: Expected 'btt_d:' (bottom depth) from input, but got "
       "'%s'\n", varnm);
     exit (-1);    
   }
   if ( (*btt_d) <= 0.0 )
   {
     printf("ERROR: btt_d (bottom depth) must be in the range (0,Inf)\n");
     exit(-1);
   }
   if ( isfinite(*btt_d) )
   {
     fscanf(fip, "%s %s", varnm, btt_tp);
     if ( strcmp (varnm, "btt_tp:") )
     {
       printf ("\nERROR: Expected 'btt_tp:' (bottom BRDF type) from input, but "
         "got '%s'\n", varnm);
       exit (-1);
     }

     if ( strcmp (btt_tp, "unidir") == 0 )
     {
       fscanf(fip, "%s %lf %lf %lf", varnm, btt_par_thetaf, btt_par_phif, 
         btt_par_rho);
       if ( strcmp (varnm, "btt_par:") )
       {
         printf ("\nERROR: Expected 'btt_par:' (parameters for the bottom BRDF)"
           " from input, but got '%s'\n", varnm);
         exit (-1);    
       }
     }
     
     if ( strcmp (btt_tp, "minnaert") == 0 )
     {
       fscanf(fip, "%s %lf", varnm, btt_par_k);
       if ( (*btt_par_k) == 0.0 )
       {
         strcpy (btt_tp, "lambert");
       }
       if ( strcmp (varnm, "btt_par:") )
       {
         printf ("\nERROR: Expected 'btt_par:' (parameters for the bottom BRDF)"
           " from input, but got '%s'\n", varnm);
         exit (-1);    
       }
     }

     // if btt_tp == "lambert", no parameters are necessary.

     fscanf(fip, "%s %d", varnm, btt_nbr);
     if ( strcmp (varnm, "btt_nbr:") )
     {
       printf ("\nERROR: Expected 'btt_nbr:' (number of bottom reflectances) "
         "from input, but got '%s'\n", varnm);
       exit (-1);    
     }
     if ( (*btt_nbr) <= 0 )
     {
       printf("ERROR: btt_nbr (number of bottom reflectances) must be > 1\n");
       exit(-1);
     }
     *btt_bhr = calloc_1d ((*btt_nbr), "btt_bhr in input_read (read_input.c)");
     fscanf(fip, "%s", varnm);
     if ( strcmp (varnm, "btt_bhr:") )
     {
       printf ("\nERROR: Expected 'btt_bhr:' (bottom bi-hemispherical "
         "reflectance) from input, but got '%s'\n", varnm);
       exit (-1);    
     }
     for(ci = 0; ci < (*btt_nbr); ci++)
     {
       fscanf(fip, "%lf", &(*btt_bhr)[ci]);
       if ( ((*btt_bhr)[ci] < 0.0) || ((*btt_bhr)[ci] > 1.0) )
       {
         printf("ERROR: btt_bhr (bottom bi-hemispherical reflectance) must be "
           "in the range [0,1]\n");
         exit(-1);
       }
     }
   } else {
     *btt_nbr = 1;
     *btt_bhr = calloc_1d ((*btt_nbr), "btt_bhr in input_read (read_input.c)");
   }


  /* Sensor Variables: ********************************************************/
   fscanf(fip, "%s %s", varnm, src_tp);
   if ( strcmp (varnm, "src_tp:") )
   {
     printf ("\nERROR: Expected 'src_tp:' (sensor type) from input, but got "
       "'%s'\n", varnm);
     exit (-1);    
   }
   if ( strcmp (src_tp, "rad")  != 0 && 
        strcmp (src_tp, "pir")  != 0 &&
        strcmp (src_tp, "sir") != 0 )
   {
     printf("ERROR: src_tp (sensor type) must be rad (radiance), pir (plane "
       "irradiance) or sir (scalar irradiance)\n");
       exit(-1);
   }
   if ( (strcmp (src_tp, "rad") == 0) || (strcmp (src_tp, "sir") == 0) )
   {
     fscanf(fip, "%s %lf", varnm, src_fov);
     if ( strcmp (varnm, "src_fov:") )
     {
       printf ("\nERROR: Expected 'src_fov:' (sensor Field-of-View) from input,"
         " but got '%s'\n", varnm);
       exit (-1);
     }
     if ( (strcmp (src_tp, "rad") == 0) && 
          ( ((*src_fov) < 0.0) || ((*src_fov) > 180.0) ) )
     {
       printf("ERROR: src_fov (sensor Field-of-View) must be in the range "
         "[0º,180º] for a radiance sensor\n");
       exit(-1);
     }
     if ( (strcmp (src_tp, "sir") == 0) && 
          ( ((*src_fov) < 0.0) || ((*src_fov) > 360.0) ) )
     {
       printf("ERROR: src_fov (sensor Field-of-View) must be in the range "
         "[0º,360º] for a scalar radiance sensor\n");
       exit(-1);
     }
     *src_fov *= RAD;
   } else {
     *src_fov = M_PI;
   } 

   fscanf(fip, "%s %lf %lf", varnm, &src_saxs[0], &src_saxs[1]);
   if ( strcmp (varnm, "src_axs:") )
   {
     printf ("\nERROR: Expected 'src_axs:' (sensor orientation axis) from "
       "input, but got '%s'\n", varnm);
     exit (-1);    
   }
   #ifdef AIRTIGHT_ALEX
   if ( src_saxs[0] < 0.0 || src_saxs[0] > 45.0 ||
        src_saxs[1] < 0.0 || src_saxs[1] > 360.0 )
   {
     printf ("\nERROR: Value of 'src_axs:' (sensor orientation axis) must be in" 
       " the range [0º,45º] for the polar angle under airtight conditions and "
       "[0º,360º] for the azimuthal angle.\n");
     exit (-1);
   }
   #else
   if ( src_saxs[0] < 0.0 || src_saxs[0] > 180.0 ||
        src_saxs[1] < 0.0 || src_saxs[1] > 360.0 )
   {
     printf ("\nERROR: Value of 'src_axs:' (sensor orientation axis) must be in"
       " the range [0º,180º] for the polar angle and [0º,360º] for the azimuthal"
       " angle.\n");
     exit (-1);    
   }
   #endif // AIRTIGHT_ALEX
   src_saxs[0] *= RAD;
   src_saxs[1] *= RAD;

   fscanf(fip, "%s %lf %lf %lf", varnm, &src_ref_o[0], &src_ref_o[1], 
     &src_ref_o[2]);
   if ( strcmp (varnm, "src_ref_o:") )
   {
     printf ("\nERROR: Expected 'src_ref_o:' (sensor reference position) from "
       "input, but got '%s'\n", varnm);
     exit (-1);    
   }
   fscanf(fip, "%s %lf %lf %lf", varnm, &src_rel_o[0], &src_rel_o[1], 
     &src_rel_o[2]);
   if ( strcmp (varnm, "src_rel_o:") )
   {
     printf ("\nERROR: Expected 'src_rel_o:' (sensor relative position) from "
       "input, but got '%s'\n", varnm);
     exit (-1);    
   }

  /* Accumulator parameters: **************************************************/
   fscanf(fip, "%s %s", varnm, acc_geom);
   if ( strcmp (varnm, "acc_geom:") )
   {
     printf ("\nERROR: Expected 'acc_geom:' (MC accumulator geometry) from "
       "input, but got '%s'\n", varnm);
     exit (-1);    
   }
   if ( strcmp (acc_geom, "none") != 0 &&
        strcmp (acc_geom, "sectorial") != 0 &&
        strcmp (acc_geom, "grid") != 0 )
   {
     printf("\nERROR: acc_geom (accumulator geometry) must be 'none' (none "
       "reported), 'sectorial' or 'grid' when SPRES (spatially resolved "
       "integral) definition is active in 'config.h'. If spatially resolved "
       "integrals are not necessary, comment out SPRES in 'config.h' and "
       "recompile for a gain in speed.\n");
     exit(-1);
   }
   if ( ( strcmp (acc_geom, "sectorial") == 0 || 
          strcmp (acc_geom, "grid") == 0 ) && !cnfg_spres )
   {
     printf("\nERROR: acc_geom (accumulator geometry) is 'grid' or 'sectorial'"
       " but code was compiled without spatially resolved integral activated. "
       "Set acc_geom to 0 or uncomment the definition SPRES in 'config.h' and "
       "recompile.\n");
     exit(-1);
   } else {
     fscanf(fip, "%s %lf", varnm, acc_ext);
     if ( strcmp (varnm, "acc_ext:") )
     {
       printf ("\nERROR: Expected 'acc_ext:' (MC accumulator extent) from "
         "input, but got '%s'\n", varnm);
       exit (-1);    
     }
     if ( (*acc_ext) <= 0.0 )
     {
       printf("ERROR: acc_ext (accumulator extent) must be > 0.\n");
       exit(-1);
     }

     fscanf(fip, "%s %lf", varnm, acc_resx);
     if ( strcmp (varnm, "acc_resx:") )
     {
       printf ("\nERROR: Expected 'acc_resx:' (MC accumulator x resolution) "
         "from input, but got '%s'\n", varnm);
       exit (-1);    
     }
     if ( (*acc_resx) <= 0.0 )
     {
       printf("ERROR: acc_resx (accumulator X-axis resolution) must be > 0.\n");
       exit(-1);
     }

     fscanf(fip, "%s %lf", varnm, acc_resy);
     if ( strcmp (varnm, "acc_resy:") )
     {
       printf ("\nERROR: Expected 'acc_resy:' (MC accumulator y resolution) "
         "from input, but got '%s'\n", varnm);
       exit (-1);    
     }
     if ( (*acc_resy) <= 0.0 )
     {
       printf("ERROR: acc_resy (accumulator Y-axis resolution) must be > 0.\n");
       exit(-1);
     }
   }
   if ( strcmp (acc_geom, "none") == 0 )
   {
     strncpy(acc_geom, "grid", STRMXLEN);
     *acc_ext   = 1.0;
     *acc_resx  = 0.2;
     *acc_resy  = 0.2;
     *acc_fgeom = 1;
   }

  /* Structures parameters: ***************************************************/
   fscanf(fip, "%s %d", varnm, str_def);
   if ( strcmp (varnm, "str_def:") )
   {
     printf ("\nERROR: Expected 'str_def:' from input, but got '%s'\n", 
       varnm);
     exit (-1);    
   }

   #ifdef AIRTIGHT_ALEX
   if ( !str_def )
   {
     printf("\nERROR: str_def (define structures flag) == 0 but code was "
       "compiled with special airtight definition that requires a structure. "
       "Set str_def to 1 and add a structure or uncomment the definition ALEX "
       "in 'config.h' and recompile.\n");
     exit(-1);
   }
   #endif // AIRTIGHT_ALEX

   if ( (*str_def) )
   {
     if ( !cnfg_shdw )
     {
       printf("\nERROR: str_def (define structures flag) != 0 but code was "
       "compiled without tracking shadows activated. Set str_def to 0 or "
       "uncomment the definition SHDW in 'config.h' and recompile.\n");
       exit(-1);
     }

     fscanf(fip, "%s", varnm);
     if ( strcmp (varnm, "str_ncl:") == 0 )
     {

       fscanf(fip, "%d", str_ncl);
       if ( (*str_ncl) > 0 )
       {

         *cylns = (str_cyln**) calloc ((*str_ncl), sizeof(str_cyln*));
         if ( (*cylns) == NULL )
         {
           printf("\nERROR: Failed to allocate memory for cylns in input_read"
             " (read_input.c).\n");
           exit(-2);
         }
         fppos = ftell(fip);
         for(ci = 0; ci < (*str_ncl); ci++)
         {

           (*cylns)[ci] = str_cyln_alloc();
           str_cyln_read(fip, &fppos, origin, axis, &alpha, &radius, &height, 
             s_base, s_top, closed);

           str_cyln_setup((*cylns)[ci], origin, axis, alpha, radius, height, 
             s_base, s_top, closed);  
         }

       }
       fscanf(fip, "%s", varnm);
     }

     if ( strcmp (varnm, "str_ncn:") == 0 )
     {
       fscanf(fip, "%d", str_ncn);
/*      if (str_ncn > 0) // NEED TO RE-IMPLEMENT CONES!
       {
         p_str_cns = (struct str_cone*) calloc (str_ncn, sizeof(struct str_cone));
         if ( p_str_cns == NULL )
         {
           printf("ERROR: Failed to allocate memory for p_str_cns (main.c).\n");
           exit(-2);
         }
         for(ci = 0; ci < str_ncn; ci++)
         {
           fscanf(fip, "%lf %lf %lf %lf %lf %lf %lf %lf %d %d", &str_cn_vertex[0], 
             &str_cn_vertex[1], &str_cn_vertex[2], &str_cn_axis[0], 
             &str_cn_axis[1], &str_cn_height[0], &str_cn_height[1], 
             &str_cn_theta, &str_cn_closed[0], &str_cn_closed[1]);
           str_cn_axis[0] *= RAD;
           str_cn_axis[1] *= RAD;
           str_cn_theta *= RAD;
           if ( (str_cn_axis[0] < 0.0) || (str_cn_axis[0] > M_PI) || 
                (str_cn_axis[1] < 0.0) || (str_cn_axis[1] > K_2PI) )
           {
             printf("ERROR: Cone axis must be between 0 and 180 for the polar "
               "angle and 0 and 360 for the azimuthal angle.\n");
             exit(-1);
           }
           if ( (str_cn_height[0] < 0) || (str_cn_height[0] > str_cn_height[1]) ) 
           {
             printf("ERROR: Cone height must be greater or equal to zero and "
               "height max must be higher than height min.\n");
             exit(-1);
           }
           if ( (str_cn_theta < 0.0) || (str_cn_theta > M_PI_2) )
           {
             printf("ERROR: Cone half angle must be between 0 and 90 degrees.\n");
             exit(-1);
           }
           setup_cone (&p_str_cns[ci], str_cn_vertex, str_cn_axis, str_cn_height, 
             str_cn_theta, str_cn_closed);
         } // Loop: ci
       } // Logic: str_ncn
*/       fscanf(fip, "%s", varnm);
     } // Logic: varnm "str_ncn"

     if ( strcmp (varnm, "str_ncb:") == 0 )
     {
       // One extra bounding box is created, to be the bounding box of all  
       // boxes, if the number of boxes is > 1:
       fscanf(fip, "%d", str_ncb);
       if ( (*str_ncb) > 0 )
       {
         (*expand_str_ncb) = 0;
         if ( (*str_ncb) > 1 ) *expand_str_ncb = 1;
         (*str_ncb) += (*expand_str_ncb);
         *cubds = (str_cubd**) calloc ((*str_ncb), sizeof(str_cubd*));
         if ( (*cubds ) == NULL )
         {
           printf("\nERROR: Failed to allocate memory for cubds in input_read"
             " (read_input.c).\n");
           exit(-2);
         }
         fppos = ftell(fip);
         for(ci = (*expand_str_ncb); ci < (*str_ncb); ci++)
         {
           (*cubds)[ci] = str_cubd_alloc();
           str_cubd_read(fip, &fppos, origin, axis, &alpha, lengths, s_base, 
             s_top, closed);
           str_cubd_setup((*cubds)[ci], origin, axis, alpha, lengths, s_base,
             s_top, closed); 
         } // Loop: ci
         if ( (*expand_str_ncb) )
         {
           // If a bounding box is to be created, it is necessary to find the 
           // limits of all boxes in the systems reference frame. This requires
           // projecting from box reference frame.
           double bbx_mn[3];
           double bbx_mx[3];
           bbx_mn[0] =  INFINITY;
           bbx_mn[1] =  INFINITY;
           bbx_mn[2] =  INFINITY;
           bbx_mx[0] = -INFINITY;
           bbx_mx[1] = -INFINITY;
           bbx_mx[2] = -INFINITY;
           double **Mt = calloc_2d(3, 3, "Mt in input_read (read_input.c)");
           double v[3];
           double v_rot[3];
           for (ci = 1; ci < (*str_ncb); ci++)
           {
             mat_transpose(Mt, (double const**) (*cubds)[ci]->M, 3, 3);
             // Test all vertices:
             for (size_t cx = 0; cx < 2; cx++)             
             {
               for (size_t cy = 0; cy < 2; cy++)             
               {
                 for (size_t cz = 0; cz < 2; cz++)             
                 {
                   v[0] = (*cubds)[ci]->x[cx];
                   v[1] = (*cubds)[ci]->y[cy];
                   v[2] = (*cubds)[ci]->z[cz];

                   rot_vec_unit(v_rot, v, (double const**) Mt);

                   for (size_t cg = 0; cg < 3; cg++)
                   {
                     v_rot[cg] += (*cubds)[ci]->o[cg];
                     if (v_rot[cg] < bbx_mn[cg])
                     {
                       bbx_mn[cg] = v_rot[cg];
                     }
                     if (v_rot[cg] > bbx_mx[cg])
                     {
                       bbx_mx[cg] = v_rot[cg];
                     }
                   } // cg
                 } // cz
               } // cy
             } // cx
           } // ci
           free_2d(3, &Mt);
           (*cubds)[0] = str_cubd_alloc();
           for (size_t cg = 0; cg < 3; cg++)
           {
             lengths[cg] = bbx_mx[cg] - bbx_mn[cg];
             origin[cg]  = bbx_mn[cg] + lengths[cg] / 2.0;
             axis[cg]    = 0.0;
             s_base[cg]  = 0.0;
             s_top[cg]   = 0.0;
           }
           alpha = 0.0;
           str_cubd_setup((*cubds)[0], origin, axis, alpha, lengths, s_base,
             s_top, closed);
         } // Bounding box calculation
       } // Setup cuboids
     } // Cuboids

     if ( (*str_ncl) == 0 && (*str_ncn) == 0 && (*str_ncb) == 0 )
     {
       printf("ERROR: str_def indicates that structures are defined, but no "
         "structure was provided.\n");
       exit(-1);
     }
   }

   fclose(fip);
 }

 void
 input_printf
 (
   int const    sim_nr,
   char const   *src_tp,
   double const src_fov,
   double const *src_pos,
   int const    sim_ns,
   double const *sim_sza,
   double const sim_f0,
   double const iop_na,
   double const iop_nw,
   double const iop_c,
   int const    iop_nw0,
   double const *iop_w0,
   char const   *scat_tp,
   double const scat_g,
   double const scat_fbb,
   double const btt_d,
   char const   *btt_tp,
   int const    btt_nbr,
   double const *btt_bhr,
   int const    acc_fgeom,
   char const   *acc_geom,
   double const acc_ext,
   double const acc_resx,
   double const acc_resy,
   int const    str_ncl,
   str_cyln     **cylns,
   int const    str_ncn,
   str_cone     **cones,
   int const    str_ncb,
   str_cubd     **cubds,
   struct source const * src,
   struct scattering const * scat,
   struct bottom const * btt,
   char const   *version
 )
 {
   int ci;

   // Terminal printout:
   printf ("\n  Backward Monte Carlo radiative transfer code\n"); 
   printf ("  Version: %s\n", version);
   printf ("  Author: Alexandre Castagna\n\n");

   printf("\n  Summary: ***************************************************\n");
   printf("\n");
   printf ("  System properties:\n");
   printf ("    Rays:    %.2e\n", (double) sim_nr);
   if (strcmp (src_tp, "rad") == 0 || strcmp (src_tp, "sir") == 0)
   {
     printf ("    Source:  %s with %.2fº FOV at z = %.2e m\n", src_tp, 
       src_fov * DEG, src_pos[2]);
   } else {
     printf ("    Source:  %s at z = %.2f m\n", src_tp, src_pos[2]);
   }
   printf ("    SZA:     (%d): ", sim_ns);
   for (ci = 0; ci < (sim_ns - 1); ci++)
   {
     if ( ci && !(ci % 6) )
       printf("\n                   ");
     printf("%.1f, ", DEG * sim_sza[ci]);
   }
   printf("%.1f (º)\n", DEG * sim_sza[ci]);
   printf ("    f0:      %.2f (W/m2)\n", sim_f0);
   printf ("    n_air:   %.2f (unitless)\n", iop_na);
   printf ("    n_water: %.2f (unitless)\n", iop_nw);
   printf ("    c coef.: %.2f (1/m)\n", iop_c);
   printf ("    w0:      (%d): ", iop_nw0);
   for (ci = 0; ci < (iop_nw0 - 1); ci++)
   {
     if ( ci && !(ci % 6) )
       printf("\n                   ");
     printf("%.3f, ", iop_w0[ci]);
   }
   printf("%.3f (unitless)\n", iop_w0[ci]);
   if (strcmp (scat_tp, "isotropic") == 0) 
   {
     printf("    Scat.:   Isotropic (Henyey-Greenstein with g = 0)\n");
   }
   if (strcmp (scat_tp, "hg") == 0) 
   {
     printf("    Scat.:   Henyey-Greenstein, with asymmetry parameter of %.3e "
       "(unitless)\n", scat_g);
   }
   if (strcmp (scat_tp, "ff") == 0) 
   {
     printf("    Scat.:   Fournier-Forand, with fbb of %.3e (unitless)\n", 
       scat_fbb);
   }
   if (strcmp (scat_tp, "petzold") == 0) 
   {
     printf("    Scat.:   Average Petzold (fbb of 0.0183, unitless)\n");
   }

   if ( isfinite(btt_d) )
   {
     if ( strcmp (btt_tp, "fixed") == 0 ) 
     {
       printf("    Bottom:  Unidirectional, at %.3e m, and with %d BHR: ", btt_d, 
         btt_nbr);
     }
     if ( strcmp (btt_tp, "minnaert") == 0 ) 
     {
       printf("    Bottom:  Minnaert, at %.3e m, and with %d BHR: ", btt_d, 
         btt_nbr);
     }
     if ( strcmp (btt_tp, "lambert") == 0 ) 
     {
       printf("    Bottom:  Lambert (Minnaert with k = 0), at %.3e m, and with %d"
        " BHR: ", btt_d, btt_nbr);
     }
     for (ci = 0; ci < (btt_nbr - 1); ci++)
     {
       printf("%.3f, ", btt_bhr[ci]);
     }
     printf("%.3f (unitless)\n", btt_bhr[ci]);
   } else {
     printf("    Bottom:  Undefined, at INFINITE depth\n");
   } // Logic: btt_d

   #ifdef SPATIALLY_RESOLVED
   if ( !acc_fgeom )
   {
     printf ("\n  Spatially resolved integral:\n");
     if ( strcmp (acc_geom, "sectorial") == 0 )
     {
       printf ("    Geometry: Sectorial\n");
     } else {
       printf ("    Geometry: Grid\n");
     }
     printf ("    Extent: %.3e (m)\n", acc_ext);
     if ( strcmp (acc_geom, "grid") == 0 )
     {
       printf ("    Resolution (x): %.3e (º)\n", acc_resx);
     } else {
       printf ("    Resolution (x): %.3e (m)\n", acc_resx);
     }
     printf ("    Resolution (y): %.3e (m)\n", acc_resy);
   } else {
     printf ("\n  Spatially resolved integral:\n");
     printf("    Code compiled with spatially resolved integral, but acc_geom was\n"
            "    set to 0. For a simulation without spatially resolved integral,\n"
            "    it will give a gain in performance to comment out the the \n"
            "    definition SPRES in config.h\n");
   } // Logic: acc_fgeom
   #endif // SPATIALLY_RESOLVED

   #ifdef SHADOWING
   printf("\n  Structures:\n");
   if ( str_ncl > 0 )
   {
     str_cylns_printf((str_cyln const **) cylns, str_ncl);
     printf("\n");
   }
   if ( str_ncn > 0 )
   {
//     str_cones_printf(cones, str_ncn);
//     printf("\n");
   }
   if ( str_ncb > 0 )
   {
     str_cubds_printf((str_cubd const **) cubds, str_ncb);
     printf("\n");
   }
   if (str_ncl == 0 && str_ncn == 0 && str_ncb == 0) 
   {
     printf("    Code compiled with shadow tracking, but no structure was\n"
            "    defined. For a simulation without shadowing, it will give\n"
            "    a gain in performance to comment out the the definition\n"
            "    SDHW in config.h\n");
   }
   #endif // SHADOWING

   #ifdef AIRTIGHT_ALEX
   printf("\nCAUTION! Code compiled with ALEX defined in config.h - this results"
          "in a specific circunstances for Alex's system and is not general, "
          "even without shadowing! Comment out ALEX in config.h for generic "
          "simulation.\n");
   #endif // AIRTIGHT_ALEX

   printf("\n  Detailed parameters: ***************************************\n");
   printf("\n");
   src_fprintf(stdout, src, 2);
   printf("\n");
   scat_fprintf(stdout, scat, 2);
   printf("\n");
   btt_fprintf(stdout, btt, 2);
   printf("\n");

   #ifdef SHADOWING
   if ( str_ncl > 0 )
   {
     for (size_t i = 0; i < str_ncl; i++)
     {
       str_cyln_fprintf(stdout, (str_cyln const *) cylns[i], 2);
       printf("\n");
     }
   }
   if ( str_ncn > 0 )
   {
   }
   if ( str_ncb > 0 )
   {
     for (size_t i = 0; i < str_ncb; i++)
     {
       str_cubd_fprintf(stdout, (str_cubd const *) cubds[i], 2);
       printf("\n");
     }
   }
   #endif // SHADOWING

   printf("\n");
 }

