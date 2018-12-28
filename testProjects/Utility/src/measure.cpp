#include "measure.h"

size_t Memory(){
#if defined(__linux__)
	pid_t pid = getpid();
	FILE *fpm;
	char command_str[1024];
	char buf[1024];
	double mem;

	sprintf(command_str, "ps up %d|awk '{print $5 }' | tail -1", pid);
	if ((fpm = popen(command_str, "r")) == nullptr) {
		fprintf(stderr, "ps failed\n");
		exit(1);
	}

	fgets(buf, 1024, fpm);
	mem = atoi(buf);
	mem /= 1024;
	pclose(fpm);
	return(mem);
#else
	typedef struct {
		DWORD  cb;
		DWORD  PageFaultCount;
		SIZE_T PeakWorkingSetSize;
		SIZE_T WorkingSetSize;
		SIZE_T QuotaPeakPagedPoolUsage;
		SIZE_T QuotaPagedPoolUsage;
		SIZE_T QuotaPeakNonPagedPoolUsage;
		SIZE_T QuotaNonPagedPoolUsage;
		SIZE_T PagefileUsage;
		SIZE_T PeakPagefileUsage;
	} info_t;
	typedef BOOL(WINAPI*func_t)(HANDLE, info_t*, DWORD);
	static func_t func; static enum { FIRST, OK, ERR } flag;
	static HANDLE proc; info_t info; HINSTANCE dll;
	if (flag == FIRST) {
		if ((dll = LoadLibraryA("psapi")) == 0) { flag = ERR; return 0; }
		func = (func_t)GetProcAddress(dll, "GetProcessMemoryInfo");
		if (func == 0) { flag = ERR; return 0; }
		proc = GetCurrentProcess(); flag = OK;
	}
	if (flag == ERR) return 0;
	func(proc, &info, sizeof info);
	return (int)info.WorkingSetSize;
#endif
}

//return time[ms]
double GetCPUTime(){
#if defined(__linux__)
	struct rusage RU;
	getrusage(RUSAGE_SELF, &RU);
	return (double)RU.ru_utime.tv_sec*1e3 + (double)RU.ru_utime.tv_usec*1e-3;
#else
	FILETIME dummyC, dummyE;
	__int64 kernelTime1, userTime1;
	GetThreadTimes(GetCurrentThread(),
		&dummyC,
		&dummyE,
		(FILETIME*)&kernelTime1,
		(FILETIME*)&userTime1);
	return userTime1*1e-3;
#endif
}
