
/*******************************************************************************
 tr_interface.c

 Alexandre Castagna Mourão e Lima (alexandre.castagna@ugent.be)
 Version: 1.6
 Date: 2021-03-25
 License: GPL-3.0

 Transmittance and reflection at dielectric interfaces

 Snell refraction function is defined and declared in the tr_interface.h header
 file.

*******************************************************************************/

 #include <stdio.h>
 #include <math.h>
 #include <complex.h>

 #include "config.h"
 #include "tr_interface.h"


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
 nrat - Pointer to the ratio of the complex refractive indexes as: n_r / n_i, 
          where i is the medium of the incident ray and r is the one of the 
          refracted ray. Note that this is the inverse of the ratio input for 
          the snell function. Note also that in this program, nrat is vector of 
          length 3, with [0] nrat, [1] nrat^2, [2] nrat^3.
 Fmat - Pointer to the Fresnel matrix structure object.
 
 OUTPUT:
 Updates the Fresnel matrix structure object (Fmat).
 
*******************************************************************************/

 void
 fresnel
 ( 
   double const mui,
   double CMPLX_T const mur,
   double CMPLX_T const ni,
   double CMPLX_T const nr,
   double CMPLX_T const *nrat,
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

   #ifdef VECTOR_RT
   double CMPLX_T xpcxp = Fmat->rp * CMPLX_F( conj )( Fmat->rp );
   double CMPLX_T xscxs = Fmat->rs * CMPLX_F( conj )( Fmat->rs );
   double CMPLX_T xpcxs = Fmat->rp * CMPLX_F( conj )( Fmat->rs );

   Fmat->R[0] =  0.5 * ( xpcxp + xscxs );
   Fmat->R[1] =  0.5 * ( xpcxp - xscxs );
   Fmat->R[2] =  CMPLX_F( creal )( xpcxs );
   Fmat->R[3] =  CMPLX_0( cimag( xpcxs ) );

   if ( CMPLX_F(creal)(mur) > TOLERANCE )
   {

/* In the current version L/n^2 is tracked and the refractive index factor is
   removed at the end. Therefore, for now the correct form for radiance
   transmittance below is not included;

     const double kt2 = CMPLX_F( creal )(nrat[2] * mur / mui);

     xpcxp = Fmat->tp * CMPLX_F( conj )( Fmat->tp );
     xscxs = Fmat->ts * CMPLX_F( conj )( Fmat->ts );
     xpcxs = Fmat->tp * CMPLX_F( conj )( Fmat->ts );

     Fmat->T[0] =  kt2 * 0.5 * ( xpcxp + xscxs );
     Fmat->T[1] =  kt2 * 0.5 * ( xpcxp - xscxs );
     Fmat->T[2] =  kt2 * CMPLX_F( creal )( xpcxs );
     Fmat->T[3] =  kt2 * CMPLX_0( cimag )( xpcxs );
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
   Fmat->R[0] = ( Fmat->rp * CMPLX_F( conj )( Fmat->rp ) + 
                  Fmat->rs * CMPLX_F( conj )( Fmat->rs ) ) * 0.5;
   Fmat->T[0] = 1.0 - Fmat->R[0];
   #endif // VETOR_RT
 }

/* reflectance out of dieletric interface */

/* void
 diel_refl
 (
   double *s;
   double *u;
   double ***stks;
   struct Fm const *Fmat,
 )
 {

 }
*/
