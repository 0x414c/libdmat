#pragma once

//#define SPINNER

#ifdef SPINNER
#define MIN_TICKS 64
void spinActivityIndicator (void);
void clearActivityIndicator (void);
#else
#define spinActivityIndicator() 
#define clearActivityIndicator() 
#endif
