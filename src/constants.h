
 #ifndef CONSTANTS
 #define CONSTANTS

 #include <math.h>

 #define RAD (M_PI / 180.0)			// Degree to radian conversion factor

 #define DEG (1 / RAD)				// Radian to degree conversion factor

 #define K_2PI 2.0 * M_PI			// 2 PI

 #define K_1_2PI 1.0 / (2.0 * M_PI)		// 1 / (2 PI)

 #define WATER_DEPOLR 0.039			// Water linear depolarization ratio

 #define SUN_OMG_A K_2PI * (1.0 - cos(0.00464))	// Subtended angle of the Sun as seen from the surface of the Earth at the mean Sun-Earth distance

 #define SUN_OMG_W SUN_OMG_A / (1.34 * 1.34)	// Subtended angle of the Sun as seen from underwater

 // Colors:
 #define ANSI_BOLD   "\033[1m"
 #define ANSI_RED    "\033[31m"
 #define ANSI_YELLOW "\033[33m"
 #define ANSI_BLUE   "\033[34m"
 #define ANSI_RESET  "\033[0m"

 #endif // CONSTANTS
