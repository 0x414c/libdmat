#pragma once

//#define WITH_SPINNER

#ifdef WITH_SPINNER
#define MIN_TICKS ( 64 )
void spinActivityIndicator (void);
void clearActivityIndicator (void);
#else //WITH_SPINNER
#define spinActivityIndicator() 
#define clearActivityIndicator() 
#endif //WITH_SPINNER
