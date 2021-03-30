
 #ifndef CONSTANTS
 #define CONSTANTS

 // Geometry:
 #define K_PI 3.14159265358979323846						// PI
 #define K_1_PI 1.0 / K_PI							// 1 / PI
 #define K_2PI 2.0 * K_PI							// 2 PI
 #define K_PI_2 0.5 * K_PI							// PI / 2
 #define K_1_2PI 1.0 / (2.0 * K_PI)						// 1 / (2 PI)
 #define RAD (K_PI / 180.0)							// Degree to radian conversion factor
 #define DEG (1 / RAD)								// Radian to degree conversion factor

 // Rayleigh scattering:
 #define WATER_DEPOLR 0.039							// Water linear depolarization ratio

 #define SUN_OMG_A K_2PI * (1.0 - cos(0.00464))					// Subtended angle of the Sun as seen from the surface of the Earth at the mean Sun-Earth distance
 #define SUN_OMG_W SUN_OMG_A / (1.34 * 1.34)					// Subtended angle of the Sun as seen from underwater

 // Wave-slope wind-speed law:
 #define WIND_K_U 0.00316				 			// Upwind wave-slope coefficient s / m
 #define WIND_K_C 0.00192 							// Crosswind wave-slope coefficient s / m
 #define WIND_DELTA 1.0
 #define WIND_EPSILON WIND_DELTA * sqrt(3.0 * WIND_K_U / WIND_K_C / 4.0)

 // Colors:
 #define ANSI_BOLD   "\033[1m"
 #define ANSI_RED    "\033[31m"
 #define ANSI_YELLOW "\033[33m"
 #define ANSI_BLUE   "\033[34m"
 #define ANSI_RESET  "\033[0m"

 #endif // CONSTANTS
