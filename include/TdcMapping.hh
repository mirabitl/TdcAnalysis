#ifndef TDC_MAPPING_HH
#define TDC_MAPPING_HH
#include<map>
#include<stdint.h>

//const int LEMO2PR[32]={14,13,12,11,10,9,8,7,6,4,2,0,1,3,5,24,
//		 26,28,30,31,29,27,25,23,22,21,20,19,18,17,16,15};
//const int PR2LEMO[32]={11,12,10,13,9,14,8,7,6,5,4,3,2,1,0,31,
//		 30,29,28,27,26,25,24,23,15,22,16,21,17,20,18,19};
//const int PR2LEMO[32]={11,12,10,13,9,14,8,7,6,5,4,3,2,1,0,31,
//		 30,29,28,27,26,25,24,23,15,22,16,21,17,20,18,19};


// From injection tests
const int TDC2LEMO[32]={19,11,18,12,20,10,17,13,21,9,16,14,22,8,15,7,23,6,24,5,25,3,26,4,-1,-1,-1,-1,-1,-1,-1,-1};

const int LEMO2TDC[32]={-1,-1,-1,21,23,19,17,15,13,9,5,1,3,7,11,
			14,10,6,2,0,4,8,12,16,18,20,22,-1,-1,-1,-1,-1};
#define NEWFIRMWARE
#ifndef NEWFIRMWARE
// From Xushian firmware

const int PR2TDC[32]={1,3,5,7,9,11,13,15,
		17,19,21,23,-1,-1,-1,-1,
		      -1,-1,-1,-1,22,20,18,16,
		14,12,10,8,6,4,2,0};

const int TDC2PR[32]={31,0,30,1,29,2,28,3,27,4,26,5,25,6,24,7,
		 23,8,22,9,21,10,20,11,-1,-1,-1,-1,-1,-1,-1,-1,};
#else
// From Sirley firmware
const int PR2TDC[32]={1,3,5,7,9,11,13,15,
		17,19,21,-1,-1,-1,-1,-1,
		      -1,-1,-1,23,22,20,18,16,
		14,12,10,8,6,4,2,0};

const int TDC2PR[32]={31,0,30,1,29,2,28,3,27,4,26,5,25,6,24,7,
		      23,8,22,9,21,10,20,19,-1,-1,-1,-1,-1,-1,-1,-1,};
#endif
// So we can build

const int PR2LEMO[32]={11,12,10,13,9,14,8,7,6,5,3,4,-1,-1,-1,-1,-1,-1,-1,-1,26,25,24,23,15,22,16,21,17,20,18,19
};

const int LEMO2PR[32]={-1,-1,-1,10,11,9,8,7,6,4,2,0,1,3,5,24,26,28,30,31,29,27,25,23,22,21,20,-1,-1,-1,-1,-1};

// Example connection to be put in JSON geometry 
const int LEMO2STRIP[32]={86,85,84,83,82,81,80,79,78,77,76,75,74,73,72,71,
		    71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86};

const int FEB2STRIP[24]={0,0,0,0,0,0,24,36,12,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};


// #define DEBUG_PRINT_ENABLED 1  // uncomment to enable DEBUG statements
#define INFO_PRINT_ENABLED 1
#if DEBUG_PRINT_ENABLED
#define INFO_PRINT_ENABLED 1
#define DEBUG_PRINTF printf
#else
#define DEBUG_PRINTF(format, args...) ((void)0)
#endif
#if INFO_PRINT_ENABLED
#define INFO_PRINTF printf
#else
#define INFO_PRINTF(format, args...) ((void)0)
#endif


#endif
