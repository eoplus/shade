
 #ifndef TR_INTERFACE
 #define TR_INTERFACE

 #include "config.h"	// CMPLX_T, CMPLX_P, CMPLX_F

 struct Fm
 {
   double CMPLX_T rp; 
   double CMPLX_T tp; 
   double CMPLX_T rs; 
   double CMPLX_T ts; 
   double R[STKS_N];
   double T[STKS_N];
 };

/* Snell refraction: ***********************************************************

 snell_s, snell_u

 The angle of refraction between two media with different complex refractive
 indexes is given by the Snell formulation:

 theta_r = asin(sin(theta_i) * n_i / n_r),                                   (1)

 where theta is the angle to the surface normal (upper and lower) and n is the 
 complex refractive index of the medium, with the subscript i and r used to 
 indicate "incident" and "refracted" variables. It is possible to express Eq. 1 
 as a function of the cosines of the angles:

 mu_r = sqrt(1 - (sqrt(1 - mu_i^2) * n_i / n_r)^2)

 Since the refracted angle will be complex if the incidence angle is larger than
 the angle of total internal reflection, the function will deal with general 
 case of complex numbers and in this program refractive indexes will always take
 the full complex expression. Treatment of complex variables is important for
 correct reflection including polarization and when the imaginary part of the
 refractive index is large (e.g., water at long wavelengths), as the refracted
 angle will be impacted (Liu et al., 2003). 

 The refraction angle in the case that the refractive indexes are real is equal
 to the refracted angle of constant phase (the real part of the complex
 refracted angle).
 
 INPUT:
 mui  - (snell_u) Pointer to the cosine of the incident angle, untiless, [0,1].
 si   - (snell_s) Incident angle in radians
 n_rat - Pointer to the ratio of the complex refractive indexes as: n_i / n_r, 
         where i is the medium of the incident ray and r is the one of the 
         refracted ray. Note that in this program, n_rat is vector of length 3,
         with [0] n_rat, [1] n_rat^2, [2] n_rat^3.

 OUTPUT:
 The (complex) cosine of the (complex) refracted angle for snell_u or the
 (complex) refracted angle for snell_s.
 
 References:
 Bohren, C. F.; Huffman, D. R. 1983. Absorption and Scattering of Light by 
   Small Particles. Wiley, New York.
 Liu, Y.; Qian, J.; Tian, Y. 2003. Succinct formulas for decomposition of 
   complex refraction angle. IEEE Antennas and Propagation Society 
   International Symposium. Digest. Held in conjunction with: USNC/CNC/URSI 
   North American Radio Sci. Meeting (Cat. No.03CH37450), Columbus, OH, pp. 
   487-490, vol. 3. DOI: 10.1109/APS.2003.1219892

*******************************************************************************/

 static inline __attribute__((always_inline)) double CMPLX_T
 snell_s
 ( 
   double const si,
   double CMPLX_T const *n_rat
 )
 { return CMPLX_P( asin )( sin(si) * n_rat[0] ); }

 #define SNELL_S(I_si, I_sr, I_n_rat) \
         (I_sr) = CMPLX_P( asin )( sin( (I_si) ) * (I_n_rat)[0] );

 static inline __attribute__((always_inline)) double CMPLX_T
 snell_u
 ( 
   double const mui,
   double CMPLX_T const *n_rat
 )
 {
   double CMPLX_T const ssin = sqrt( 1.0 - mui * mui ) * n_rat[0];
   return CMPLX_P(sqrt) ( 1.0 - ssin * ssin );
 }

 #define SNELL_U(I_mui, I_mur, I_n_rat, I_ssin) \
         (I_ssin) = sqrt( 1.0 - (I_mui) * (I_mui) ) * (I_n_rat)[0]; \
         (I_mur)  = CMPLX_P(sqrt)( 1.0 - (I_ssin) * (I_ssin) );

/* Fresnel Reflection and Transmission matrices: *******************************

 fresnel

 As the Stokes column vector is represented in this program in the common form:
 S = | I |
     | Q |
     | U |
     | V |,
 
 The reflection and transmission matrices for the diffuse Stokes column vector
 are given by:

 R = | R(1,1)  R(1,2)    0       0    |
     | R(2,1)  R(2,2)    0       0    |
     |   0       0     R(3,3)  R(3,4) |
     |   0       0     R(4,3)  R(4,4) |,

 T = | T(1,1)  T(1,2)    0       0    | * (n_r/n_i)^3 * mu_r / mu_i,
     | T(2,1)  T(2,2)    0       0    |
     |   0       0     T(3,3)  T(3,4) |
     |   0       0     T(4,3)  T(4,4) |

 where:
 X(1,1) =           0.5 * (r_p * Conj(r_p) + r_s * Conj(r_s))
 X(1,2) =  X(2,1) = 0.5 * (r_p * Conj(r_p) - r_s * Conj(r_s))
 X(3,3) =  X(4,4) = Real(t_p * Conj(t_s))
 X(3,4) = -X(4,3) = Imag(t_p * Conj(t_s)),

 with X being R or T. The subscripts p and s identify parallel and perpendicular
 polarization and are given by:
 r_p = ( n_i * mu_r - n_r * mu_i ) / ( n_i * mu_r + n_r * mu_i )
 t_p = ( 2.0 * n_i * mu_i ) / ( n_i * mu_r + n_r * mu_i )
 r_s = ( n_i * mu_i - n_r * mu_r ) / ( n_i * mu_i + n_r * mu_r )
 t_s = ( 2.0 * n_i * mu_i ) / ( n_i * mu_i + n_r * mu_r ),

 where mu is the cosine and n is the complex index of refraction, with the 
 subscript i and r used to indicate "incident" and "refracted" variables.

 The body of the function below will depend on how the Fresnel matrices are
 represented. In Version 1.6, the code is optimized for directionally 
 isotropic media with mirror-symnmetric particles, such that, only four 
 elements of the Mueller and Fresnel matrices need to be tracked: (1,1), 
 (1,2), (3,3), (3,4), the other being 0 or equal (or the negative) to these,
 as in the definitions above.

 If solving in the scalar approximation to the RTE, the reflection and 
 transmittance coefficients are given by the first element of the R and T 
 matrices.

 The Fresnel function is provided for the diffuse Stokes transmittance (taking
 into account the refractive index change) and as a colimated Stokes
 transmittance, which does not take into account the refractive index change.
 This last function is used on the MC code since the MC code tracks L/n^2 and
 the refractive index factor is corrected for in the end.

 During initial tests, inlining was performed but performance was slower than
 directly using macros, though no difference was expected. Therefore the
 MC function is also provided as a macro. Since a series of intermediate
 variables are necessary, another macro is provided that expands into the
 definition of those variables.

 References:
 Bohren, C. F.; Huffman, D. R. 1983. Absorption and Scattering of Light by 
   Small Particles. Wiley, New York.
 Garcia, R. D. M. 2012. Fresnel boundary and interface conditions for 
   polarized radiative transfer in a multilayer medium. Journal of 
   Quantitative Spectroscopy & Radiative Transfer, 113, 306–317. 
   DOI: 10.1016/j.jqsrt.2011.11.015
 Zhai, P.-W.; Kattawar, G. W.; Hu, Y. 2012. Comment on the transmission matrix  
   for a dielectric interface. Journal of Quantitative Spectroscopy &  
   Radiative Transfer 113, 1981-1984. DOI: 10.1016/j.jqsrt.2012.07.001

 INPUT:
 mui  - Pointer to the cosine of the incident angle, untiless, [0,1].
 mur  - Pointer to the complex cosine of the refracted angle, unitless, [].
 ni   - Pointer to the complex index of refraction of the medium of incidence.
 nr   - Pointer to the complex index of refraction of the medium of 
        refraction.
 n_rat - Pointer to the ratio of the complex refractive indexes as: n_r / n_i, 
          where i is the medium of the incident ray and r is the one of the 
          refracted ray. Note that this is the inverse of the ratio input for 
          the snell function. Note also that in this program, n_rat is vector of 
          length 3, with [0] n_rat, [1] n_rat^2, [2] n_rat^3.
 Fmat - Pointer to the Fresnel matrix structure object.
 
 OUTPUT:
 Updates the Fresnel matrix structure object (Fmat).
 
*******************************************************************************/

 static inline __attribute__((always_inline)) void
 fresnel
 ( 
   double const mui,
   double CMPLX_T const mur,
   double CMPLX_T const ni,
   double CMPLX_T const nr,
   double CMPLX_T const *n_rat,
   struct Fm *Fmat
 )
 {
   double CMPLX_T const nimi = ni * mui;
   double CMPLX_T const nimr = ni * mur;
   double CMPLX_T const nrmi = nr * mui;
   double CMPLX_T const nrmr = nr * mur;

   Fmat->rp = ( nrmi - nimr ) / ( nrmi + nimr );
   Fmat->tp = ( 2.0 * nimi  ) / ( nrmi + nimr );
   Fmat->rs = ( nimi - nrmr ) / ( nimi + nrmr );
   Fmat->ts = ( 2.0 * nimi  ) / ( nimi + nrmr );

   double CMPLX_T xpcxp = Fmat->rp * CMPLX_F( conj )( Fmat->rp );
   double CMPLX_T xscxs = Fmat->rs * CMPLX_F( conj )( Fmat->rs );

   #ifdef VECTOR_RT
   double CMPLX_T xpcxs = Fmat->rp * CMPLX_F( conj )( Fmat->rs );
   Fmat->R[0] =  0.5 * ( xpcxp + xscxs );
   Fmat->R[1] =  0.5 * ( xpcxp - xscxs );
   Fmat->R[2] =  CMPLX_F( creal )( xpcxs );
   Fmat->R[3] =  CMPLX_0( cimag( xpcxs ) );

   if ( CMPLX_F(creal)(mur) > TOLERANCE ) // SEEMS IT SHOULD BE MU CRITICAL AND NOT TOLERANCE... CHECK!
   {
     xpcxp = Fmat->tp * CMPLX_F( conj )( Fmat->tp );
     xscxs = Fmat->ts * CMPLX_F( conj )( Fmat->ts );
     xpcxs = Fmat->tp * CMPLX_F( conj )( Fmat->ts );

     const double kt2 = CMPLX_F( creal )(n_rat[2] * mur / mui);

     Fmat->T[0] =  kt2 * 0.5 * ( xpcxp + xscxs );
     Fmat->T[1] =  kt2 * 0.5 * ( xpcxp - xscxs );
     Fmat->T[2] =  kt2 * CMPLX_F( creal )( xpcxs );
     Fmat->T[3] =  kt2 * CMPLX_0( cimag )( xpcxs );
   } else {
     Fmat->T[0] =  0.0;
     Fmat->T[1] =  0.0;
     Fmat->T[2] =  0.0;
     Fmat->T[3] =  0.0;
   }
   #else
   Fmat->R[0] = 0.5 * ( xpcxp + xscxs );
   Fmat->T[0] = 1.0 - Fmat->R[0];
   #endif // VETOR_RT
 }

 static inline __attribute__((always_inline)) void
 fresnel_mc
 ( 
   double const mui,
   double CMPLX_T const mur,
   double CMPLX_T const ni,
   double CMPLX_T const nr,
   struct Fm *Fmat
 )
 {
   double CMPLX_T const nimi = ni * mui;
   double CMPLX_T const nimr = ni * mur;
   double CMPLX_T const nrmi = nr * mui;
   double CMPLX_T const nrmr = nr * mur;

   Fmat->rp = ( nrmi - nimr ) / ( nrmi + nimr );
   Fmat->tp = ( 2.0 * nimi  ) / ( nrmi + nimr );
   Fmat->rs = ( nimi - nrmr ) / ( nimi + nrmr );
   Fmat->ts = ( 2.0 * nimi  ) / ( nimi + nrmr );

   double CMPLX_T xpcxp = Fmat->rp * CMPLX_F( conj )( Fmat->rp );
   double CMPLX_T xscxs = Fmat->rs * CMPLX_F( conj )( Fmat->rs );

   #ifdef VECTOR_RT
   double CMPLX_T xpcxs = Fmat->rp * CMPLX_F( conj )( Fmat->rs );

   Fmat->R[0] =  0.5 * ( xpcxp + xscxs );
   Fmat->R[1] =  0.5 * ( xpcxp - xscxs );
   Fmat->R[2] =  CMPLX_F( creal )( xpcxs );
   Fmat->R[3] =  CMPLX_0( cimag( xpcxs ) );

   if ( CMPLX_F(creal)(mur) > TOLERANCE )
   {
     /*
     In the current version L/n^2 is tracked and the refractive index factor is
     removed at the end. Therefore, for now the correct form for radiance
     transmittance below is not included
     */
     Fmat->T[0] = 1.0 - Fmat->R[0];
     Fmat->T[1] = 1.0 - Fmat->R[1];
     Fmat->T[2] = 1.0 - Fmat->R[2];
     Fmat->T[3] = 1.0 - Fmat->R[3];
   } else {
     Fmat->T[0] =  0.0;
     Fmat->T[1] =  0.0;
     Fmat->T[2] =  0.0;
     Fmat->T[3] =  0.0;
   }
   #else
   Fmat->R[0] = 0.5 * ( xpcxp + xscxs );
   Fmat->T[0] = 1.0 - Fmat->R[0];
   #endif // VETOR_RT
 }

 #define FRESNEL_MC_VAR \
         double CMPLX_T V_nimi;  \
         double CMPLX_T V_nimr;  \
         double CMPLX_T V_nrmi;  \
         double CMPLX_T V_nrmr;  \
         double CMPLX_T V_xpcxp; \
         double CMPLX_T V_xscxs; \
         double CMPLX_T V_xpcxs;

 #ifdef VECTOR_RT
 #define FRESNEL_MC_U(I_mui, I_mur, I_ni, I_nr, I_Fmat) \
         V_nimi = (I_ni) * (I_mui); \
         V_nimr = (I_ni) * (I_mur); \
         V_nrmi = (I_nr) * (I_mui); \
         V_nrmr = (I_nr) * (I_mur); \
         (I_Fmat).rp = ( V_nrmi - V_nimr ) / ( V_nrmi + V_nimr ); \
         (I_Fmat).tp = ( 2.0 * V_nimi  ) / ( V_nrmi + V_nimr ); \
         (I_Fmat).rs = ( V_nimi - V_nrmr ) / ( V_nimi + V_nrmr ); \
         (I_Fmat).ts = ( 2.0 * V_nimi  ) / ( V_nimi + V_nrmr ); \
         V_xpcxp = (I_Fmat).rp * CMPLX_F( conj )( (I_Fmat).rp ); \
         V_xscxs = (I_Fmat).rs * CMPLX_F( conj )( (I_Fmat).rs ); \
         (I_Fmat).R[0] =  0.5 * ( V_xpcxp + V_xscxs ); \
         V_xpcxs = (I_Fmat)->rp * CMPLX_F( conj )( (I_Fmat)->rs ); \
         (I_Fmat).R[1] =  0.5 * ( V_xpcxp - V_xscxs ); \
         (I_Fmat).R[2] =  CMPLX_F( creal )( V_xpcxs ); \
         (I_Fmat).R[3] =  CMPLX_0( cimag( V_xpcxs ) ); \
         if ( CMPLX_F(creal)( I_mur ) > TOLERANCE ) \
         { \
           for(size_t i = 0; i < STKS_N; i++) \
           { \
             (I_Fmat).T[i] = 1.0 - (I_Fmat).R[i]; \
           } \
        } \
        else \
        { \
           for(size_t i = 0; i < STKS_N; i++) \
           { \
             (I_Fmat).T[i] = 0.0; \
           } \
        }
 #else
 #define FRESNEL_MC_U(I_mui, I_mur, I_ni, I_nr, I_Fmat) \
         V_nimi = (I_ni) * (I_mui); \
         V_nimr = (I_ni) * (I_mur); \
         V_nrmi = (I_nr) * (I_mui); \
         V_nrmr = (I_nr) * (I_mur); \
         (I_Fmat).rp = ( V_nrmi - V_nimr ) / ( V_nrmi + V_nimr ); \
         (I_Fmat).tp = ( 2.0 * V_nimi  ) / ( V_nrmi + V_nimr ); \
         (I_Fmat).rs = ( V_nimi - V_nrmr ) / ( V_nimi + V_nrmr ); \
         (I_Fmat).ts = ( 2.0 * V_nimi  ) / ( V_nimi + V_nrmr ); \
         V_xpcxp = (I_Fmat).rp * CMPLX_F( conj )( (I_Fmat).rp ); \
         V_xscxs = (I_Fmat).rs * CMPLX_F( conj )( (I_Fmat).rs ); \
         (I_Fmat).R[0] =  0.5 * ( V_xpcxp + V_xscxs ); \
         if ( CMPLX_F(creal)( I_mur ) >= TOLERANCE ) \
         { \
           (I_Fmat).T[0] = 1.0 - (I_Fmat).R[0]; \
         } \
        else \
        { \
           (I_Fmat).T[0] = 0.0; \
        }
 #endif // VECTOR_RT

/* void
 diel_refl
 (
   double *s;
   double *u;
   double ***stks;
   struct Fm const *Fmat
 );
*/

 #endif // TR_INTERFACE

