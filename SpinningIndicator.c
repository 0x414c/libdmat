#include "SpinningIndicator.h"

#ifdef SPINNER
#include <stdio.h>
#include <stdint.h>
#ifdef _MSC_VER
#include <Windows.h>
#if WINVER >= 0x0600
#define GetTickCount() GetTickCount64();
#endif // WINVER
#endif // _MSC_VER

#ifdef __GNUC__
#include <sys/time.h>
uint64_t GetTickCount() {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (tv.tv_sec * 1000 + tv.tv_usec / 1000);
}
#endif // __GNUC__

static uint8_t currentPos = 0;
static char chars[] = { '-', '\\', '|', '/', '-', '\\', '|', '/' };

void spinActivityIndicator (void) {
	static uint64_t lastTickTime = 0, currentTickTime = 0;

	currentTickTime = GetTickCount();
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
