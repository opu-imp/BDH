/**
* @file    Measure.h
* @author  T.Sato
* @date    2015.05.03
* @version 1.0
*/

#if !defined(__MEASURE__)
#define __MEASURE__
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#if defined(__linux__)
#include <sys/types.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <unistd.h>
#else
#include <time.h>
#include <Windows.h>
#endif

using namespace std;

/******************************** system ********************************************/
/**
* @brief Check Memory useage[byte]
* @return memory usage[byte]
*/
size_t Memory();

/**
* @brief get cpu time[ms]
* @return cpu time[ms]
*/
double GetCPUTime();

#endif /*__MEASURE__*/
