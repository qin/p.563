/********************************************************************
ITU-T Draft Recommendation P.563
Version 1.0 - 23 March 2004
  
NOTICE 
  
The Single Ended Assessment Model P.563 algorithm and the copyright therein 
is the joint property of Psytechnics Limited, OPTICOM GmbH and SwissQual AG 
and is protected by UK, US and other patents, either applied for or registered. 
Permission is granted to use this source code solely for the purpose of 
evaluation of ITU-T recommendation P.563. 
Any other use of this software requires a licence, which may be obtained from: 
  
OPTICOM GmbH 
Am Weichselgarten 7, D- 91058 Erlangen, Germany 
Phone: +49 9131 691 160			Fax: +49 9131 691 325  
E-mail: info@opticom.de         www.3sqm.com  
  
Psytechnics Limited 
Fraser House, 23 Museum Street, Ipswich, IP1 1HN, UK 
Phone: +44 1 473 261 800		Fax: +44 1 473 261 880 
E-mail: info@psytechnics.com    www.psytechnics.com
  
SwissQual AG 
Gewerbestrasse 2 CH-4528 Zuchwil, Switzerland 
Phone: +41 32 685 08 30			Fax: +41 32 685 08 31   
E-mail: sales@swissqual.com     www.swissqual.com
  
Psytechnics, SwissQual or Opticom can provide licences and further information. 
  
Authors: 
      Ludovic Malfait ludovic.malfait@psytechnics.com 
      Roland Bitto rb@opticom.de 
      Pero Juric pero.juric@swissqual.com

********************************************************************/


#ifndef __GENERIC_TYPEDEFS_H
#define __GENERIC_TYPEDEFS_H

#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
 #define CPP "C"
#else
 #define CPP
#endif


/* Type definitions */

typedef struct 
{
	FLOAT	SpecLevelDev;
	FLOAT	SpecLevelRange;
	FLOAT	EstSegSNR;
	FLOAT	RelNoiseFloor;
} typSegSNR_results;


typedef struct 
{
	FILE   	*Fptr;
    INT32	filetype;
	INT32	frame_len;
	INT32	frame_shft;
	INT16	useWholeFile;
    INT32	start;
	INT32	stop;
	INT32	FileSize;
	FLOAT	sumsq;
	FLOAT	AverageRMS;
	FLOAT	DCOffset;
} fs;


typedef struct 
{
	INT16	 packetmove;
	INT32	 start;
	INT32	 length;
	INT32	 SummOfAllInterrupts;
	FLOAT  rms;
	FLOAT  EstimatedRms;
	FLOAT  startSlope;
} typInterruption;



typedef struct 
{
	INT32    Frequency;
	INT32    bitResolution;
	FLOAT   MOS;
	FLOAT  *ReferData;
	FLOAT  *CodedData;
} typInputParameter;


typedef struct 
{
	fs			file[2];
	INT32		fileSize;
	INT32		FFTSize;
	INT32		FrameCnt;
	INT32		MaxNrFrames;
	FLOAT		NoiseLevelLog;
	FLOAT		sampleLength;
	FLOAT		RefLevelCorrection;
	FLOAT	   *dataBuffer;
	FLOAT	   *dataBufferPrev;
	FLOAT	   *dataBufferNext;
} typChannel;


typedef struct 
{
	FLOAT	 meanHist;
	FLOAT  maxHist;
	FLOAT  minHist;
	FLOAT  minInputVal;
	FLOAT  maxInputVal;
	FLOAT	*Histogram;
	FLOAT	*NormIndex;
	INT32   *index;
	INT32    histLength;
	INT32    HistMode;
} histogram;


typedef struct
{
	FLOAT* pdDenomCoefA;
	FLOAT* pdNumCoefB;
	FLOAT* pdState;

	INT32 lDenomCoefALength;
	INT32 lNumCoefBLength;
	INT32 lStateLength;
} sctFilterInfo;


#endif
