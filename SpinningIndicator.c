#include "SpinningIndicator.h"

#ifdef SPINNER
#include <stdio.h>
#include <stdint.h>


#ifdef __GNUC__
#pragma message "Got GCC..."
#include <sys/time.h>
uint64_t GetTickCounts() {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (tv.tv_sec * 1000 + tv.tv_usec / 1000);
}
#endif // __GNUC__

#ifdef _MSC_VER
#pragma message("Got MSVC...")
#include <Windows.h>
#if WINVER >= 0x0600
#define GetTickCounts() GetTickCount64()
#else
#define GetTickCounts() GetTickCount()
#endif // WINVER
#endif // _MSC_VER


static uint8_t currentPos = 0;
static char chars[] = { '-', '\\', '|', '/', '-', '\\', '|', '/' };

void spinActivityIndicator (void) {
	static uint64_t lastTickTime = 0, currentTickTime = 0;

	currentTickTime = GetTickCounts();
	if (currentTickTime < lastTickTime + MIN_TICKS) {
		return;
	} else {
		lastTickTime = currentTickTime;
	}

	if (currentPos >= sizeof(chars)) {
		currentPos = 0;
	}

	putchar(chars[currentPos++]);
	putchar('\b');

	return;
}

void clearActivityIndicator (void) {
	putchar(' ');
	putchar('\b');

	return;
}
#endif // SPINNER
