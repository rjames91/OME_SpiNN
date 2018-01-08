/*
 * IHC_AN.h
 *
 *  Created on: 1 Nov 2016
 *      Author: rjames
 */

#ifndef OME_SpiNN_H_
#define OME_SpiNN_H_

#define REAL float
#define REAL_CONST(x) x##f
#define REAL_MAX 1e38
#define PI 3.14159265358979323846f

#define MAX_CHIPX 1//255
#define MAX_CHIPY 1//255
#define MAX_COREID 16
//#define SEED_SEL_SIZE (((MAX_CHIPX << 8) | MAX_CHIPY) << 8) | (MAX_COREID << 3)   
#define SEED_SEL_SIZE 1024

#define SEGSIZE 96//100

//#define SAMPLING_FREQUENCY 44100.f

#define MAX_SIGNAL_S 1

#define TIMER2_CONF        0x82
#define TIMER2_LOAD        0

#define max(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

#define ABS(x) (((x)<0) ? -(x) : (x))
#define SIGN(x) (((x)<0) ? -(1.0f) : (1.0f))

typedef union
{
	uint32_t u;
	REAL f;
} uint_float_union;


#endif /* IHC_AN_H_ */
