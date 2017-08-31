/********************************************************************
ITU-T Draft Recommendation P.563
Version 1.0 - 23 March 2004

NOTICE

The Single Ended Assessment Model P.563 algorithm and the copyright therein
is the joint property of Psytechnics Limited, OPTICOM GmbH and SwissQual AG
and is protected by UK, US and other patents, either applied for or
registered.
Permission is granted to use this source code solely for the purpose of
evaluation of ITU-T recommendation P.563.
Any other use of this software requires a licence, which may be obtained
from:

OPTICOM GmbH
Am Weichselgarten 7, D- 91058 Erlangen, Germany
Phone: +49 9131 691 160   Fax: +49 9131 691 325
E-mail: info@opticom.de         www.3sqm.com

Psytechnics Limited
Fraser House, 23 Museum Street, Ipswich, IP1 1HN, UK
Phone: +44 1 473 261 800  Fax: +44 1 473 261 880
E-mail: info@psytechnics.com    www.psytechnics.com

SwissQual AG
Gewerbestrasse 2 CH-4528 Zuchwil, Switzerland
Phone: +41 32 685 08 30   Fax: +41 32 685 08 31
E-mail: sales@swissqual.com     www.swissqual.com

Psytechnics, SwissQual or Opticom can provide licences and further
information.

Authors:
      Ludovic Malfait ludovic.malfait@psytechnics.com
      Roland Bitto rb@opticom.de
      Pero Juric pero.juric@swissqual.com

********************************************************************/


#include <stdio.h>
#include <math.h>
// #include <malloc.h>
#include <stdlib.h>

#include "defines.h"
#include "dsp.h"
#include "Statistics.h"

/*********************
* FUNCTION: MultiplyWith
*
* DESCRIPTION: Vector Multiplication with constant
*
***********************/

void MultiplyWith (FLOAT           *pTimeSeries,
                   INT32             pNumberOfSamples,
                   FLOAT			  pScale)
{

    INT32 i;
    for (i = 0; i < pNumberOfSamples; i++)
	 {
        pTimeSeries [i] *= pScale;
    }
}

/*********************
* FUNCTION: LocalPower
*
* DESCRIPTION: Evaluate the power of an arry in
*					the range started by pStartIndex to pStopIdex
*
***********************/

FLOAT LocalPower (const FLOAT  *pTimeSeries,
                  INT32         pStartIndex,
                  INT32         pStopIndex)
{

    INT32     i;
    FLOAT  power;


	 power = 0;


    for (i = pStartIndex; i < pStopIndex; i++)
	 {
        FLOAT h = pTimeSeries [i];
        power += h * h;
    }

    power /= (pStopIndex - pStartIndex);
    return (FLOAT) power;
}

/*********************
* FUNCTION: GlobalPower
*
* DESCRIPTION: Evaluate the power of an array
*
***********************/

FLOAT GlobalPower (const FLOAT *pTimeSeries,
                   INT32       pNumberOfSamples)
{
    return LocalPower (pTimeSeries,
                       0,
                       pNumberOfSamples);
}


/*********************
*
* FUNCTION: QuickSortIncreasing
*
* DESCRIPTION:
*	sort elements in increasing order
*
***********************/


void QuickSortIncreasing (FLOAT *pfThis,
                          INT32 iLeftIndex,
                          INT32 iRightIndex,
                          INT32 *piThat)
{
    INT32    i, j;
    FLOAT    x, w;
    INT32    ww;

    i = iLeftIndex;
    j = iRightIndex;
    x = pfThis [(iLeftIndex + iRightIndex)/2];

    do {
        while (pfThis [i] < x) {i++;}
        while (pfThis [j] > x) {j--;}
        if (i > j) break;
        w = pfThis [i]; pfThis [i] = pfThis [j]; pfThis [j] = w;
        ww = piThat [i]; piThat [i] = piThat [j]; piThat [j] = ww;
    } while (++i <= --j);

    if (iLeftIndex < j) {QuickSortIncreasing (pfThis, iLeftIndex, j, piThat);}
    if (i < iRightIndex) {QuickSortIncreasing (pfThis, i, iRightIndex, piThat);}
}

void SortIncreasing (FLOAT * pfThis,
							INT32    pN,
                     INT32   *piThat)
{

    QuickSortIncreasing (pfThis, 0, pN - 1, piThat);
}


/*********************
*
* FUNCTION: interpolate
*
* DESCRIPTION:
*	Interpolation filter
*
*
***********************/
static FLOAT interpolate (FLOAT    freq,
									const FLOAT   fppFilterCurve [][2],
									INT32      iNumberOfPoints) {

FLOAT  result;
INT32     i;
FLOAT  freqLow, freqHigh;
FLOAT  curveLow, curveHigh;

    if (freq <= fppFilterCurve [0][0]) {
        freqLow = fppFilterCurve [0][0];
        curveLow = fppFilterCurve [0][1];
        freqHigh = fppFilterCurve [1][0];
        curveHigh = fppFilterCurve [1][1];

        result = ((freq - freqLow) * curveHigh + (freqHigh - freq) * curveLow)/ (freqHigh - freqLow);

        return (float) result;
    }

    if (freq >= fppFilterCurve [iNumberOfPoints-1][0]) {
        freqLow = fppFilterCurve [iNumberOfPoints-2][0];
        curveLow = fppFilterCurve [iNumberOfPoints-2][1];
        freqHigh = fppFilterCurve [iNumberOfPoints-1][0];
        curveHigh = fppFilterCurve [iNumberOfPoints-1][1];

        result = ((freq - freqLow) * curveHigh + (freqHigh - freq) * curveLow)/ (freqHigh - freqLow);

        return (float) result;
    }

    i = 1;
    freqHigh = fppFilterCurve [i][0];
    while (freqHigh < freq) {
        i++;
        freqHigh = fppFilterCurve [i][0];
    }
    curveHigh = fppFilterCurve [i][1];

    freqLow = fppFilterCurve [i-1][0];
    curveLow = fppFilterCurve [i-1][1];

    result = ((freq - freqLow) * curveHigh + (freqHigh - freq) * curveLow)/ (freqHigh - freqLow);

    return (float) result;
}



/*********************
*
* FUNCTION: FftFilter
*
* DESCRIPTION:
*	Filter the signal using a fft Filter
*
***********************/

void FftFilter( FLOAT * pfData, INT32 iNrSamples, INT32 iNumberOfPoints,const FLOAT Fs, const FLOAT fppFilterCurve [][2] ,FLOAT *fpFilterdData)
{
    INT32    pow_of_2    = nextpow2 (iNrSamples);
    FLOAT    *x = NULL;

    FLOAT    factorDb, factor;

    FLOAT   overallGainFilter = interpolate ((FLOAT) 1000, fppFilterCurve, iNumberOfPoints);
    FLOAT   freq_resolution;
    INT32    i;


	 x = (FLOAT *) calloc ((pow_of_2 + 2) , sizeof (FLOAT));


    for (i = 0; i < iNrSamples; i++) {
        x [i] = pfData[i];
    }

    RealFFT (x, pow_of_2);

    freq_resolution = (FLOAT) Fs / (FLOAT) pow_of_2;

    for (i = 0; i <= pow_of_2/2; i++)
	 {
        factorDb = interpolate (i * freq_resolution, fppFilterCurve, iNumberOfPoints) - overallGainFilter;
        factor = (FLOAT) pow ((FLOAT) 10, factorDb / (FLOAT) 20);

        x [2 * i] *= factor;
        x [2 * i + 1] *= factor;
    }


    RealIFFT (x, pow_of_2);

    for (i = 0; i < iNrSamples; i++) {
        fpFilterdData[i] = x[i];
    }

    free (x);
}



/*********************
* FUNCTION: StandardIRSFilter
*
* DESCRIPTION: Filter signal using IRS like Curve
*
***********************/


void StandardIRSFilter(FLOAT *pfSignal,
						INT32 lSignalLen,
						FLOAT *pfFilterdSignal)
{


const FLOAT IRSLikeCurvedB [26][2] ={{0.0F,-100.0F},
												{ 50.0F, -27.5F},
												{100.0F, -27.5F},
												{125.0F, -18.8F},
												{160.0F, -10.8F},
												{200.0F, -2.7F},
												{250.0F, 2.7F},
												{300.0F, 6.4F},
												{315.0F, 7.2F},
												{400.0F, 9.9F},
												{500.0F, 11.3F},
												{600.0F, 11.8F},
												{630.0F, 11.9F},
												{800.0F, 12.3F},
												{1000.0F, 12.6F},
												{1250.0F, 12.5F},
												{1600.0F, 13.0F},
												{2000.0F, 13.1F},
												{2200.0F, 13.1F},
												{2700.0F, 12.5F},
												{2850.0F, 12.6F},
												{3200.0F, -8.1F},
												{3700.0F, -43.6F},
												{4700.0F, -66.9F},
												{6000.0F, -79.5F},
												{8000.0F, -102.0F}};


	FftFilter( pfSignal, lSignalLen, 26 ,SAMPLE_FREQUENCY, IRSLikeCurvedB ,pfFilterdSignal);

}





/*********************
* FUNCTION: GetFilteredLevel
*
* DESCRIPTION: Evaluate singal level in the range of 250-3000 Hz
*
***********************/


FLOAT  GetFilteredLevel(FLOAT *pfSignal,
								INT32 lSignalLen)
{

	FLOAT *pFilterdSignal=NULL;
	FLOAT fLevel=0;

	const FLOAT gAlignFilterdB [26][2] = {{0.0F,-500.0F},
													{50.0F, -500},
													{100.0F, -500},
													{125.0F, -500},
													{160.0F, -500},
													{200.0F, -500},
													{250.0F,  0},
													{300.0F,  0},
													{315.0F,  0},
													{400.0F,  0},
													{500.0F,  0},
													{600.0F,  0},
													{630.0F,  0},
													{800.0F,  0},
													{1000.0F, 0},
													{1250.0F, 0},
													{1600.0F, 0},
													{2000.0F, 0},
													{2500.0F, -5},
													{3000.0F, -10},
													{3150.0F, -20},
													{3500.0F, -50},
													{4000.0F, -500},
													{5000.0F, -500},
													{6300.0F, -500},
													{8000.0F, -500}};

	pFilterdSignal = (FLOAT *) calloc (lSignalLen , sizeof (FLOAT));

	FftFilter( pfSignal, lSignalLen, 26 ,SAMPLE_FREQUENCY, gAlignFilterdB ,pFilterdSignal);

   fLevel = GlobalPower (pFilterdSignal, lSignalLen);

	free(pFilterdSignal);

	return(fLevel);

}
