
 #ifndef CONFIG
 #define CONFIG

// #define SHADOWING
// #define SPATIALLY_RESOLVED
// #define AIRTIGHT_ALEX

 #ifdef AIRTIGHT_ALEX
   #ifndef SHADOWING
     #error AIRTIGHT_ALEX defined will only calculate shadowing, so SHADOWING must also be defined. 
   #endif
 #endif 

 #define SCALAR_RT
 #ifndef SCALAR_RT
   #define VECTOR_RT
 #endif

 #ifdef SCALAR_RT
   #define STKS_N 1
 #else
   #define STKS_N 4
 #endif // SCALAR_RT

 #define CHCK_MEM_ALLOC		// Check each memory allocation

 #define TOLERANCE 1E-12
 #define CRTCW 1E-6		// Critical weight below which a ray is terminated
 #define STRMXLEN 201		// Maximum length of strings (including reading)
 #define SUBDV 10		// Subdivisions of the total number ray simulations
 #define SCAT_NPSI 200001	// Number of scattering angles in look-up table
 #define SCAT_DIR "anc/pf/"	// Directory for scattering phase functions

 // If complex refractive indexes are used the interface functions for Snell and
 // Fresnell will use complex numbers and complex versions of functions. 
 // CMPLX_T will add the qualifier 'complex' where appropriate and CMPLX_P will
 // add the profix 'c' as a prefix to functions, as the complex.h library use a
 // c prefix for the complex versions of the functions in the math.h library.

 #ifdef COMPLEX_N
   #define CMPLX_T complex		// Adds the 'complex' qualifier to variables (e.g., double CMPLX_T const x -> double complex const x)
   #define CMPLX_P(EXP) c ## EXP	// concatenates a 'c' before a function name (e.g., CMPLX_P(sqrt)(x) -> csqrt(x))
   #define CMPLX_F(EXP) EXP		// No effect (e.g., CMPLX_F(creal)(x) -> creal(x) )
   #define CMPLX_0(EXP) EXP		// No effect (e.g., CMPLX_0(cimag(x)) -> cimag(x) )
 #else
   #define CMPLX_T			// No effect (e.g., double CMPLX_T const x -> double const x)
   #define CMPLX_P(EXP) EXP		// No effect (e.g., CMPLX_P(sqrt)(x) -> sqrt(x))
   #define CMPLX_F(EXP)			// Remove a complex only function (e.g., CMPLX_F(creal)(x) -> (x) )
   #define CMPLX_0(EXP) 0.0		// Substitute a complex expression by 0 (e.g., CMPLX_0(cimag(x)) -> 0.0 )
 #endif // COMPLEX_N

 #endif

