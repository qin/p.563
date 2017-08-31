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

#include <stdlib.h>
#include <math.h>
#include <assert.h>

#include "defines.h"
#include "beeprob.h"
#include "Statistics.h"
#include "dsp.h"
#include "tools1.h"


#include <stdio.h>



/*********************
*
* FUNCTION: L1
*
* DESCRIPTION: Mean absolute value computation over an interval 
*	
*
***********************/

FLOAT L1 (const FLOAT  *pTimeSeries,
          INT32           pStartIndex, 
          INT32           pStopIndex) 
{


    INT32    i;
    const FLOAT  *p = &(pTimeSeries [pStartIndex]);
    FLOAT  result = 0;

    for (i = pStopIndex - pStartIndex - 1; i >= 0; i--) {
        result += (FLOAT) fabs (*(p++));
    }
        
    result /= (pStopIndex - pStartIndex);
    return (FLOAT) result;
}


/*********************
*
* FUNCTION: Envelope
*
* DESCRIPTION:  sine windowed wheighted energy
*	
*
***********************/

FLOAT Envelope (const FLOAT *pTimeSeries, INT32 pNumberOfSamples, INT32 pStartIndex, INT32 pEnvelopeFilterLength) 
{
    INT32            i;
    FLOAT          envelope;
    static  INT32     storedEnvelopeFilterLength = -1;
    static  FLOAT  *filterCoeffs = NULL;

  
    if (pEnvelopeFilterLength != storedEnvelopeFilterLength) 
	 {
        
        if (filterCoeffs != NULL) 
		  {
            free (filterCoeffs);
        }
        filterCoeffs = (FLOAT *) malloc (sizeof (FLOAT) * pEnvelopeFilterLength);

        for (i = 0; i < pEnvelopeFilterLength; i++) 
		  {
            filterCoeffs [i] = (FLOAT) sin (PI * i/(FLOAT) pEnvelopeFilterLength);
        }

        storedEnvelopeFilterLength = pEnvelopeFilterLength;
    }

    envelope = 0;

    for (i = 0; i < pEnvelopeFilterLength; i++) 
	 {
        if (pStartIndex + i < pNumberOfSamples) 
		  {
            FLOAT h;
            h = filterCoeffs [i] * pTimeSeries [pStartIndex + i];
            envelope += h * h;
        }
    }
    
    envelope /= (pEnvelopeFilterLength / 2);
    envelope =(FLOAT) sqrt (envelope);
    return (FLOAT) envelope;
}


/*********************
*
* FUNCTION: WindowedSample
*
* DESCRIPTION: scale time series using a hann window 
*	
*
***********************/


static FLOAT WindowedSample (const FLOAT    *pTimeSeries,
                      INT32             pStartSample, 
                      INT32             pSampleIndex, 
                      INT32             pWindowSize) 
{


    static  INT32      storedWindowSize = -1;
    static  FLOAT   *storedWeights = NULL;    
    FLOAT            result;
    INT32              i;
    
    assert ((0 <= pSampleIndex) && (pSampleIndex < pWindowSize));
    
    if (pWindowSize != storedWindowSize) {
        INT32 sampleIndex;


        if (storedWeights != NULL) {
            free(storedWeights);
        }
        storedWeights =(FLOAT*)calloc(pWindowSize + 1,sizeof(FLOAT));
 
    
        for (sampleIndex = 0; sampleIndex <= pWindowSize; sampleIndex++) { 
            storedWeights [sampleIndex] = (FLOAT) (0.5 - 0.5 * cos (2 * PI * sampleIndex / pWindowSize));
        }
    
        storedWindowSize = pWindowSize;
    }

    i = pStartSample + pSampleIndex;

    result = storedWeights [pSampleIndex] * pTimeSeries [i];

    return result;
}

/*********************
*
* FUNCTION: SpectralUnFlatness
*
* DESCRIPTION: spectral Flatness evaluation 
* using computation in the frequency domain
*	
*
***********************/

FLOAT SpectralUnFlatness (const FLOAT *  pTimeSeries,
                          INT32           pStartSample, 
                          FLOAT         pLowerFrequencyHz, 
                          FLOAT         pUpperFrequencyHz,
                          FLOAT         pConstant,
                          FLOAT *        x1,
                          FLOAT *        y) 
{

    INT32              i;
    FLOAT            frequencyResolution; 
    FLOAT           sumPower = 0;
    FLOAT           sumLog = 0;
    INT32              count = 0;
    FLOAT           avgPower;
    FLOAT           avgLog;
    FLOAT            spectralUnFlatness;

    if (pStartSample < 0) {
        pStartSample = 0;
    }

    for (i = 0; i < TRANSFORM_LENGTH + 2; i++) {
        x1 [i] = 0.;
    }

    for (i = 0; i < TRANSFORM_LENGTH; i++) {
        x1 [i] = WindowedSample (pTimeSeries, pStartSample, i, TRANSFORM_LENGTH);
    }

	 {

		INT32 Len=0x01;
		Len=Len << (LOG2_TRANSFORM_LENGTH);
		RealFFT(x1, Len);

	 }
    
    frequencyResolution = (FLOAT) SAMPLE_FREQUENCY_HZ / (FLOAT) TRANSFORM_LENGTH;
   
    for (i = 0; i <= TRANSFORM_LENGTH / 2; i++) 
	 { 
        FLOAT freq = i * frequencyResolution;
        y [i] = x1 [2*i] * x1 [2*i] + x1 [2*i + 1] * x1 [2*i + 1]; 
        if ((freq >= pLowerFrequencyHz) && (freq <= pUpperFrequencyHz)) 
		  {
            sumPower += y[i] + pConstant;
            sumLog +=(FLOAT) log (y[i] + pConstant);
            count++;
        }
    }  
    
    assert (count > 0);
    assert (sumPower > 0);
    
    avgPower = sumPower / count;
    avgLog = sumLog / count;

    spectralUnFlatness = (FLOAT) (log (avgPower) - avgLog) / (FLOAT) log (10);
    return spectralUnFlatness;
}

/*********************
*
* FUNCTION: Periodicity
*
* DESCRIPTION: Periodicity evaluation. 
* using computation in the frequency domain
*	
*
***********************/

FLOAT Periodicity(const FLOAT *fpTimeSeriesA,
						const INT32 iSizeTimeSerA,
						const	FLOAT *fpTimeSeriesB,
						const INT32 iSizeTimeSerB,
						const	INT32	 iTransformLength,
						const	FLOAT  fLowerFrequencyHz, 
						const	FLOAT  fUpperFrequencyHz,
						const	FLOAT  fConstant,
						const INT32	 iMode,
								FLOAT *x1,
								FLOAT *x2,
								FLOAT *y,
								INT32 *plBestDelayIndex) 
{

    FLOAT  energy1 = 0, energy2 = 0, normalization;
    INT32  i;
    INT32  bestDelay;
    FLOAT h;
    FLOAT frequencyResolution; 
    FLOAT maximumAbsoluteCorrelation = (FLOAT) 0.;
	 INT32 iFFTLen=iTransformLength;
 
    for (i = 0; i < iTransformLength + 2; i++) 
	 {
        x1 [i] = 0.0F;
        x2 [i] = 0.0F;
		   y [i] = 0.0F;
    }

	 for (i = 0; i < iSizeTimeSerA; i++) 
	 {
        x1 [i] = fpTimeSeriesA[i];
    }
	 for (i = 0; i < iSizeTimeSerB; i++) 
	 {
        x2 [i] = fpTimeSeriesB[i];
    }


	RealFFT(x1, iFFTLen);
	RealFFT(x2, iFFTLen);


    frequencyResolution = (FLOAT) SAMPLE_FREQUENCY_HZ / (FLOAT) iTransformLength;

    for (i = 0; i <= iTransformLength / 2; i++) 
	 {  
        FLOAT freq = i * frequencyResolution;
    
        if ((freq >= fLowerFrequencyHz) && (freq <= fUpperFrequencyHz)) 
		  {
             
            y [2*i] = x1 [2*i] * x2 [2*i] + x1 [2*i + 1] * x2 [2*i + 1];
    
            y [2*i + 1] = -x1 [2*i + 1] * x2 [2*i] + x1 [2*i] * x2 [2*i + 1];
        } 
		  else
		  {
            y [2 * i] = 0;   
            y [2 * i + 1] = 0;   

            x1 [2 * i] = 0;   
            x1 [2 * i + 1] = 0;   
            
            x2 [2 * i] = 0;   
            x2 [2 * i + 1] = 0;   
        }

    }
    
		RealIFFT(y, iFFTLen);
		RealIFFT(x1, iFFTLen);
		RealIFFT(x2, iFFTLen);

 
    for (i = 0; i < iTransformLength; i++) 
	 {
        energy1 += x1 [i] * x1 [i];
        energy2 += x2 [i] * x2 [i];
    }

	 normalization=1;
	 if(iMode==0)
	 {
		energy1 += fConstant * iTransformLength;
		energy2 += fConstant * iTransformLength;
		normalization =(FLOAT) sqrt (energy1 * energy2);
	 }
	 if(iMode==1)
	 {
		if(energy1==0)	energy1 += fConstant * iTransformLength;
		if(energy2==0)energy2 += fConstant * iTransformLength;

		normalization =(FLOAT) sqrt (energy1 * energy2);

	 }


    bestDelay = 0;

    for (i = -iTransformLength / 2; i <= -1; i++) 
	 { 
		  h = (FLOAT)(y [(i + iTransformLength)]) / normalization;
		  if(iMode==0)h=(FLOAT)fabs(h);
       
        if (h > (FLOAT) maximumAbsoluteCorrelation) 
		  {
            maximumAbsoluteCorrelation = h;
            bestDelay = i;
        }
    }

    for (i = 0; i < iTransformLength / 2; i++) 
	 {
        h = (FLOAT)(y [i]) / normalization;
		  if(iMode==0)h=(FLOAT)fabs(h);

        if (h >= (FLOAT) maximumAbsoluteCorrelation) 
		  {
            maximumAbsoluteCorrelation = h;
            bestDelay = i;
        }
    }
    
	 *plBestDelayIndex=bestDelay;

    return maximumAbsoluteCorrelation;
}

/*********************
*
* FUNCTION: UnnaturalBeeps
*
* DESCRIPTION: The dedector for unnatural beeps considers the peroidicities in time singnal frames of 32 ms
*					over a total lenth of 160 ms whith shifts of 16ms. The periodicity itself is drived from the 
*					spectral unflatnes measure n a frequency range between 250 and 1200 Hz.
*	
*
***********************/


INT16 UnnaturalBeeps(FLOAT *pProcessVector,
                     INT32  pBlockLengthInSamples,
							FLOAT	*fpUnnaturalBeeps,
							FLOAT	*fpMeanUnnaturalBeeps,
							FLOAT	*fpUnBeepsMeanDistSamp
							)
{

    FLOAT   windowShift = TRANSFORM_LENGTH / 10.0f; 
    INT32     lWindowNum = (INT32) (1 + (pBlockLengthInSamples - TRANSFORM_LENGTH)/ windowShift);
    FLOAT *  periodicity = calloc(lWindowNum,sizeof(FLOAT));
    INT32     k;
    INT32     numberOfBeeps = 0;

	FLOAT *pHelperBufferX1;
	FLOAT *pHelperBufferY;

	SISTATISTIC_HANDLE hStatistics =NULL;

	SiStatisticsCreate(&hStatistics,"UnnaturalBeeps");

	pHelperBufferX1 = (FLOAT*)calloc(TRANSFORM_LENGTH + 2,sizeof(FLOAT));
	pHelperBufferY = (FLOAT*)calloc(TRANSFORM_LENGTH + 2,sizeof(FLOAT));


    for (k = 0; k < lWindowNum; k++) 
	 {
        INT32 startSample = (INT32) (k * windowShift);
    
        periodicity [k] = SpectralUnFlatness (pProcessVector,
                                              startSample, 
                                              250.0f, 
                                              1200.0f, 
                                              1,
                                              pHelperBufferX1,
                                              pHelperBufferY);  
    }  

    for (k = 10; k < lWindowNum - 10; k++) {
        INT32 startSample = (INT32) (k * windowShift);
        INT32 stopSample = startSample + TRANSFORM_LENGTH;
        FLOAT l11 = (FLOAT) L1 (pProcessVector, startSample, startSample + TRANSFORM_LENGTH);

        FLOAT periodicity1 = periodicity [k ];
        FLOAT periodicity2 = periodicity [k - 5];
        FLOAT periodicity3 = periodicity [k + 5];


        if ((l11 > 1E3) && (periodicity1 > 1.54) && (periodicity2 < 0.89) && (periodicity3 < 0.89)) 
		  
		  {
            
            numberOfBeeps++;
	
				if(hStatistics != NULL)
				{

					SiStatisticsSetDefaultVector(hStatistics,pProcessVector,startSample,stopSample,1);
				}
        }
    }



	
	*fpUnnaturalBeeps = (FLOAT) numberOfBeeps * 1000 / (FLOAT) pBlockLengthInSamples; 

	{
		 FLOAT dMean=0;
		 FLOAT dStd=0;
		 INT32 lNrEntries=0;

		SiStatisticsGetMoments(hStatistics,&dMean,&dStd,&lNrEntries);
		dMean=(FLOAT)(fabs(dMean));
		*fpMeanUnnaturalBeeps=dMean;

		*fpUnBeepsMeanDistSamp=(FLOAT)hStatistics->lNrEntries/(FLOAT)pBlockLengthInSamples;
	 }

	SiStatisticsDelete(&hStatistics);

	free(periodicity);
	free(pHelperBufferX1);
	free(pHelperBufferY);

    return ((INT16)(numberOfBeeps > 0));
} 


/*********************
*
* FUNCTION: Robotization
*
* DESCRIPTION: Using the frequencies in the range between 2200Hz and 3300 the 
*					the periodicitie is calculated over adjacent time singal frames of 32ms
*					if the periodic frame among the nonsilenc frame exceeds 3.4% at least 3.4% frame 
*					are declared as robotic
*	
***********************/

#define KK          2.0 /* Best 2 */
#define KK_SPLIT    20

BOOL Robotization(const FLOAT *pProcessVector,
								INT32  pBlockLengthInSamples,
                        FLOAT *pFractionActive,
								FLOAT *pfRobotization)
{

    INT32     periodWindowLength = (INT32) (KK * TRANSFORM_LENGTH); 
    INT32     periodWindowShift = (INT32) (TRANSFORM_LENGTH * KK / (FLOAT) KK_SPLIT); 
    INT32     numberOfPeriodWindows = 1 + (pBlockLengthInSamples - periodWindowLength)/ periodWindowShift;
    FLOAT *  periodicity =(FLOAT *)calloc(numberOfPeriodWindows,sizeof(FLOAT));

    INT16   *isPeriodic = (INT16 *)malloc(numberOfPeriodWindows * sizeof (INT16));
    INT32    *sortedPeriodk = (INT32 *)calloc(numberOfPeriodWindows , sizeof (INT32));
    INT32     numberOfPerodicWindows = 0;
    INT32     numberOfNonSilentWindows = 0;
    INT32     periodk;
	 INT32   iDelayIndex=0;

    FLOAT   fractionPeriodic;
    BOOL    result = FALSE;
	 FLOAT		fMOVValue=0;


	FLOAT *pHelperBufferX1;
	FLOAT *pHelperBufferX2;
	FLOAT *pHelperBufferY;

	SISTATISTIC_HANDLE hStatistics =NULL;
	SiStatisticsCreate(&hStatistics,"Robotization");

	pHelperBufferX1 = (FLOAT*)calloc(TRANSFORM_LENGTH + 2,sizeof(FLOAT));
	pHelperBufferX2 = (FLOAT*)calloc(TRANSFORM_LENGTH + 2,sizeof(FLOAT));
	pHelperBufferY = (FLOAT*)calloc(TRANSFORM_LENGTH + 2,sizeof(FLOAT));



    for (periodk = 0; periodk < numberOfPeriodWindows; periodk++) {
        INT32 startSample = periodk * periodWindowShift;
        INT32 stopSample = startSample + periodWindowLength;
        FLOAT power = LocalPower (pProcessVector, startSample, stopSample);
        FLOAT l1 = L1 (pProcessVector,  startSample, stopSample);

        if (power > 1E6) {
            numberOfNonSilentWindows++;
        }

        sortedPeriodk [periodk] = (INT32) periodk;
        
        periodicity [periodk] = Periodicity(&pProcessVector[startSample],
                                            TRANSFORM_LENGTH,
														  &pProcessVector[startSample+TRANSFORM_LENGTH],
														  TRANSFORM_LENGTH,
                                            TRANSFORM_LENGTH, 
                                            2200.0f, 
                                            3300.0f, 
                                            1E3,
														  0,
                                            pHelperBufferX1,
                                            pHelperBufferX2,
                                            pHelperBufferY,
														  &iDelayIndex);       


#define PERIODICITY_LIMIT  0.84

        if ((periodicity [periodk] > PERIODICITY_LIMIT) && (l1 > 2.5E2)) 
		  { 
            isPeriodic [periodk] = TRUE;
            numberOfPerodicWindows++;
        } 
		  else 
		  {
            isPeriodic [periodk] = FALSE;
        }
    } 

    fractionPeriodic = (FLOAT) numberOfPerodicWindows / (FLOAT) numberOfNonSilentWindows;
    *pFractionActive = (FLOAT) numberOfNonSilentWindows / (FLOAT) numberOfPeriodWindows;
    
#define MAX_FRAC_PERIODIC   0.034 

    if (fractionPeriodic > MAX_FRAC_PERIODIC) 
	 {
        INT32 i;
        INT32 j = numberOfPeriodWindows - 1;

        result = TRUE;

        SortIncreasing (periodicity, numberOfPeriodWindows, sortedPeriodk);

        while (periodicity [j] > PERIODICITY_LIMIT) {
            j--;
        }

        j += (INT32) (MAX_FRAC_PERIODIC * numberOfNonSilentWindows);
        
        j = (INT32) (0.45 * j + 0.55 * numberOfPeriodWindows);       

        
			for (i = j; i < numberOfPeriodWindows; i++) 
			{
          
            INT32 startSamplePeriodWindow = sortedPeriodk [i] * periodWindowShift + TRANSFORM_LENGTH / 2;
            INT32 stopSamplePeriodWindow = startSamplePeriodWindow + periodWindowLength;
			
            if (startSamplePeriodWindow < 0) { startSamplePeriodWindow = 0;}
            if (stopSamplePeriodWindow >= pBlockLengthInSamples) {stopSamplePeriodWindow = pBlockLengthInSamples - 1;}

				if(hStatistics != NULL)
				{
						
							SiStatisticsSetDefaultVector(hStatistics,pProcessVector,startSamplePeriodWindow,stopSamplePeriodWindow,1);
				}


			} 

	 	  fMOVValue = (FLOAT) (((j - i) - 1) * periodWindowShift + periodWindowLength) / (FLOAT) pBlockLengthInSamples;

    }

	if (fMOVValue < 0){fMOVValue = 0;}
	*pfRobotization=fMOVValue * 1000;


	free(periodicity);
	free(isPeriodic);
	free(sortedPeriodk);

	free(pHelperBufferX1);
	free(pHelperBufferX2);
	free(pHelperBufferY);


	SiStatisticsDelete(&hStatistics);


    return result;

} 

/*********************
*
* FUNCTION: FrameRepeats
*
* DESCRIPTION: By evaluation of the usualy high cross correlation of repeated frames 
*             this function determines the nr of repeted frames of 256 samples and the energy of the 
*					repeated parts in a signa.l
*	
*
***********************/


void FrameRepeats (const FLOAT *  pProcessVector,
                         INT32    pBlockLengthInSamples,
								 FLOAT *pfFrRepNum,
								 FLOAT *pfFrRepRelEnergy)
{

   INT32	lwindowShift =32; 
	INT32	lWindowLength = 256;
   INT32  lWindowNum = (INT32) (1 + (pBlockLengthInSamples - lWindowLength)/ lwindowShift);
   INT32  k = 0;
	INT32	iBestWindowIdx = 0;
	FLOAT	fMaximumKorr = 0;
	FLOAT	sumEnvelopeEnergy = 0;
	INT32	numFrameRepeats = 0;
   INT32   iDelayIndex=0;

	FLOAT *  periodicity = (FLOAT*)calloc(lWindowNum,sizeof(FLOAT));
	FLOAT*  pEnvelope = (FLOAT*)calloc(lWindowNum,sizeof(FLOAT));
	INT32*	pBestDelay = (INT32*)calloc(lWindowNum,sizeof(INT32));

	FLOAT *pHelperBufferX1;
	FLOAT *pHelperBufferX2;
	FLOAT *pHelperBufferY;
	SISTATISTIC_HANDLE hStatistics =NULL;


	INT32 iTrLength1=TRANSFORM_LENGTH/2;
	INT32 iTrLength2=1024;

	SiStatisticsCreate(&hStatistics,"FrameRepeats");



	pHelperBufferX1 = (FLOAT*)calloc(max(iTrLength1,iTrLength2) + 2,sizeof(FLOAT));
	pHelperBufferX2 = (FLOAT*)calloc(max(iTrLength1,iTrLength2) + 2,sizeof(FLOAT));
	pHelperBufferY = (FLOAT*)calloc(max(iTrLength1,iTrLength2) + 2,sizeof(FLOAT));



    for (k = 0; k < lWindowNum; k++)
	{
			INT32		startSample = (INT32) (k * lwindowShift);
			pEnvelope[k] = Envelope(pProcessVector, pBlockLengthInSamples, startSample, lWindowLength);   		
		
			if ((pBlockLengthInSamples > (startSample + 5*TRANSFORM_LENGTH )) & (pEnvelope[k] > 100.0f))
			{


			if(Periodicity(&pProcessVector[startSample],
                        iTrLength1,
								&pProcessVector[startSample+iTrLength1],
								iTrLength1,
                        iTrLength1, 
                        0.0f, 
                        4000.0f,
								1,
								0,
                        pHelperBufferX1,
                        pHelperBufferX2,
                        pHelperBufferY, 
								&iDelayIndex) <0.6)
														 
				{

					periodicity [k] =	Periodicity(&pProcessVector[startSample],
									      256,
											&pProcessVector[startSample+256],
											512,
											1024, 
											0.0f, 
											4000.0f,
											1,
											1,
											pHelperBufferX1,
											pHelperBufferX2,
											pHelperBufferY, 
											&(pBestDelay[k]));



				if ((fmod((FLOAT)(k),8.f) != 0) || (k == 0))
				{
					if (periodicity [k] > fMaximumKorr)
					{
						fMaximumKorr = periodicity [k];
						iBestWindowIdx = k;
					}
				}
				else
				{

					if ((fMaximumKorr > 0.72f))/*0.72f */
					{
						INT32 idx = iBestWindowIdx + (INT32)((256 + pBestDelay[iBestWindowIdx])/ (FLOAT)(lwindowShift));
						
						if ((pBestDelay[iBestWindowIdx] > -17.0f) & (pBestDelay[iBestWindowIdx] < 100.0f))
						{
							sumEnvelopeEnergy += pEnvelope[idx];
							numFrameRepeats++;
	  
							if(hStatistics != NULL)
							{
						
									SiStatisticsSetDefaultVector(hStatistics,pProcessVector,startSample,startSample+lWindowLength,1);
							}
							
						}
						else
						{
							if ((fMaximumKorr > 0.75f) & (pEnvelope[idx] > 100.0f))
							{
								sumEnvelopeEnergy += pEnvelope[idx];
								numFrameRepeats++;

								if(hStatistics != NULL)
								{
				
									SiStatisticsSetDefaultVector(hStatistics,pProcessVector,startSample,startSample+lWindowLength,1);
								}
							}

						}
					}
					fMaximumKorr = periodicity [k];
					iBestWindowIdx = k;
				}
			}
		}

	}

	*pfFrRepNum=(FLOAT)(numFrameRepeats);
	*pfFrRepRelEnergy=(FLOAT)hStatistics->dSummValues/(FLOAT)pBlockLengthInSamples;
	free(pEnvelope);

	free(periodicity);
	free(pBestDelay);

	free(pHelperBufferX1);
	free(pHelperBufferX2);
	free(pHelperBufferY);

	SiStatisticsDelete(&hStatistics);

} 



/*********************
*
* FUNCTION: SharpDeclines
*
* DESCRIPTION:this dedector works in several steps:
*					The signal power of two frames each with a duration of 10 ms are compared to each other.
*					if the ratio is larger than 10.5 two further indicators are calculated.
*					These indicators are derived from a periodicity measure. They determine if the periodicity in the 
*					past blocks are below a certain threshold, if the power ration is larger than 24 
*					and if the signal containes enough enrgy.
*	
*
***********************/

#define LL          0.33333333333
#define LL_SPLIT    10

INT16 SharpDeclines (const FLOAT *   pProcessVector,
									INT32    pBlockLengthInSamples,
									FLOAT *fpSharpDeclines)
                          
{

    INT32         shortWindowLength = (INT32) (TRANSFORM_LENGTH * LL);
    INT32         shortWindowShift = (INT32) (TRANSFORM_LENGTH * LL / LL_SPLIT); 
    INT32         numberOfShortWindows = 1 + (pBlockLengthInSamples - shortWindowLength)/ shortWindowShift;
    FLOAT       currentShortPower;
    FLOAT       pastShortPower;
    INT32         shortk; 

    BOOL        result = FALSE;
    INT32         numberOfSharpDeclines = 0;
	 INT32   iDelayIndex=0;
	
	FLOAT *pHelperBufferX1;
	FLOAT *pHelperBufferX2;
	FLOAT *pHelperBufferY;

	pHelperBufferX1 = (FLOAT*)calloc(TRANSFORM_LENGTH + 2,sizeof(FLOAT));
	pHelperBufferX2 = (FLOAT*)calloc(TRANSFORM_LENGTH + 2,sizeof(FLOAT));
	pHelperBufferY = (FLOAT*)calloc(TRANSFORM_LENGTH + 2,sizeof(FLOAT));

	*fpSharpDeclines=0;


    for (shortk = LL_SPLIT / 2; shortk < numberOfShortWindows; shortk++) 
	 {
        INT32 currentStartSample = shortk * shortWindowShift * 1;
        INT32 currentStopSample = currentStartSample + shortWindowLength;
        INT32 pastStartSample = (shortk - LL_SPLIT / 2) * shortWindowShift * 1;
        INT32 pastStopSample = pastStartSample + shortWindowLength;
        FLOAT factor;
            
	 	  currentShortPower = (FLOAT) sqrt (LocalPower (pProcessVector, currentStartSample, currentStopSample));
        currentShortPower += 10; 

        pastShortPower = (FLOAT) sqrt (LocalPower (pProcessVector, pastStartSample, pastStopSample));
        pastShortPower += 10; 

        factor = pastShortPower/ currentShortPower;

        if (factor > 10.5) 
		  {
            FLOAT longEnvelope;
            FLOAT periodicity1 =0;
				FLOAT periodicity2 =0;
				INT32 farStartSample = currentStartSample - TRANSFORM_LENGTH; 
				INT32 iPerStart=0;


				
				iPerStart=currentStartSample - TRANSFORM_LENGTH;
				if(iPerStart < 0 )iPerStart=0;

				periodicity1=Periodicity(&pProcessVector[iPerStart],
									            TRANSFORM_LENGTH/2,
													&pProcessVector[iPerStart+TRANSFORM_LENGTH/2],
													TRANSFORM_LENGTH/2,
													TRANSFORM_LENGTH/2, 
												   250.0f, 
                                       1000.0f, 
                                       1,
													0,
													pHelperBufferX1,
													pHelperBufferX2,
													pHelperBufferY,
													&iDelayIndex);
				
				
				iPerStart=currentStartSample - 2*TRANSFORM_LENGTH;
				if(iPerStart < 0 )iPerStart=0;

				periodicity2=Periodicity(&pProcessVector[iPerStart],
									            TRANSFORM_LENGTH/2,
													&pProcessVector[iPerStart+TRANSFORM_LENGTH/2],
													TRANSFORM_LENGTH/2,
													TRANSFORM_LENGTH/2, 
												   250.0f, 
                                       1000.0f, 
                                       1,
													0,
													pHelperBufferX1,
													pHelperBufferX2,
													pHelperBufferY,
													&iDelayIndex);

                  
            if (farStartSample < 0) 
				{
                farStartSample = 0;
            }
            longEnvelope = (FLOAT) L1 (pProcessVector,  farStartSample, currentStopSample); 
            

            if ((longEnvelope > 150 ) && 
                ((factor > 24 ) || ((periodicity1 < 0.70) && (periodicity2 < 0.65)))) 
				{
                result = TRUE;
                numberOfSharpDeclines++;
            } 
        }
    }

	
	*fpSharpDeclines = (FLOAT) numberOfSharpDeclines * 1000 / (FLOAT) pBlockLengthInSamples;
	
	free(pHelperBufferX1);
	free(pHelperBufferX2);
	free(pHelperBufferY);


    return result;
}   


#define SUBSAMPLE_FACTOR_ENVELOPES_1c       2560                                                               
#define UUUU                                8        
#define MAX_NUMBER_OF_DIPS                  50

#define IN_DIP_ENVELOPE                     100.f  
#define OUT_OF_DIP_ENVELOPE                 300.f  
    
#define DURATION_INTERVAL_1c                2.6f


/*********************
*
* FUNCTION: UnnaturalSilences
*
* DESCRIPTION: Detection of unnatural silenced speech signals
*				using an hysteresis approach the algorithm iterates through the envelopes of 320ms frames 
*				with a shift of 40ms, and determines the number of transitions within low to high 
*				and vice versa within an interval of 2.6 seconds
*	 
***********************/


INT32 UnnaturalSilences (const FLOAT *   pProcessVector,
									INT32      pBlockLengthInSamples,
									FLOAT	*fpMeanUnnatSilence
								 )
												 
{   
                                            
    INT32     numberOfEnvelopes = 1 + (pBlockLengthInSamples - SUBSAMPLE_FACTOR_ENVELOPES_1c)/ (SUBSAMPLE_FACTOR_ENVELOPES_1c / UUUU);
    FLOAT *  envelopes = (FLOAT *) calloc (numberOfEnvelopes , sizeof (FLOAT)); 
    INT32     envelopeIndex;
    INT32     numberOfDips = 0;
    INT32    *dipStartEnvelope = (INT32 *) calloc (MAX_NUMBER_OF_DIPS , sizeof (INT32));
    INT32    *dipStopEnvelope = (INT32 *) calloc (MAX_NUMBER_OF_DIPS , sizeof (INT32));
    INT32    *dipTransitionEnvelope = (INT32 *) calloc (MAX_NUMBER_OF_DIPS , sizeof (INT32));
    BOOL    statusInDip = TRUE;
    INT32     numberOfEnvelopesInAnInterval = 1 + (INT32) ((DURATION_INTERVAL_1c  * (FLOAT) SAMPLE_FREQUENCY - SUBSAMPLE_FACTOR_ENVELOPES_1c) / (FLOAT) (SUBSAMPLE_FACTOR_ENVELOPES_1c /UUUU));
    INT32     i;
    BOOL    result = FALSE;

	 INT32 iNumChangedSamples = 0;

	SISTATISTIC_HANDLE hStatistics =NULL;

	SiStatisticsCreate(&hStatistics,"UnnaturalSilences");


    /* Fill envelope signal */
    for (envelopeIndex = 0; envelopeIndex < numberOfEnvelopes; envelopeIndex++) 
	 {
        FLOAT sum = 0;
        FLOAT  h;
        INT32    i;
        INT32   sampleIndex = envelopeIndex * SUBSAMPLE_FACTOR_ENVELOPES_1c / UUUU;
            
        for (i = 0; i < SUBSAMPLE_FACTOR_ENVELOPES_1c; i++) {
            h = pProcessVector [sampleIndex];    
            sum += h * h;
            sampleIndex++;
        }

        sum /= SUBSAMPLE_FACTOR_ENVELOPES_1c;

        envelopes [envelopeIndex] = (FLOAT) sqrt (sum);
    }

    /* Find dips of energy */
    for (envelopeIndex = 0; envelopeIndex < numberOfEnvelopes; envelopeIndex++) {
        if (statusInDip) {
            if (envelopes [envelopeIndex] > OUT_OF_DIP_ENVELOPE) {
                statusInDip = FALSE;
                dipStopEnvelope [numberOfDips] = envelopeIndex+UUUU;
                dipTransitionEnvelope [numberOfDips] = envelopeIndex+UUUU;
                numberOfDips++;
                if (numberOfDips >= MAX_NUMBER_OF_DIPS) {
                    exit (1);
                }

            }
        } else {
            if (envelopes [envelopeIndex] < IN_DIP_ENVELOPE) {             
                statusInDip = TRUE;
                dipStartEnvelope [numberOfDips] = envelopeIndex;
                dipTransitionEnvelope [numberOfDips] = envelopeIndex;
                numberOfDips++;
                if (numberOfDips >= MAX_NUMBER_OF_DIPS) {
                    exit (1);
                }
            }
        }
    }        
    /* */
    for (i = 0; i < numberOfDips; i++) {

		if (dipTransitionEnvelope [i] > 0) {
            INT32 numberOfTransitionsInIntervalStartingAti = 1;

            INT32 j = i + 1;
            while ((j < numberOfDips - 1) && (dipTransitionEnvelope [j] < dipTransitionEnvelope[i] + numberOfEnvelopesInAnInterval)) {
                j++;
                numberOfTransitionsInIntervalStartingAti++;
            }

            if (numberOfTransitionsInIntervalStartingAti >= 4) 
				{
                /* Add noise */
                INT32 startEnvelope = dipTransitionEnvelope [(i/2) * 2 + 1]; /* odd indices are starting env's of dips */
                INT32 startSample = startEnvelope * SUBSAMPLE_FACTOR_ENVELOPES_1c / UUUU;
                INT32 stopSample = startSample + SUBSAMPLE_FACTOR_ENVELOPES_1c;
                
                if (startSample < 0) {startSample = 0;}
                if (stopSample > pBlockLengthInSamples) {stopSample = pBlockLengthInSamples;}

                result = TRUE;

					if(hStatistics != NULL)
					{
						/* add data to statistics */
						SiStatisticsSetDefaultVector(hStatistics,pProcessVector,startSample,stopSample,1);
					}
				/* I assume that there is no overlapping of the frames! */
				iNumChangedSamples += stopSample - startSample;
            }
        }
    }	 
	 
	 
	 {
		 FLOAT dMean=0;
		 FLOAT dStd=0;
		 INT32 lNrEntries=0;

		SiStatisticsGetMoments(hStatistics,&dMean,&dStd,&lNrEntries);
		dMean=(FLOAT)(fabs(dMean)+1.0);
		dMean=(FLOAT)(10*log10(dMean));
		*fpMeanUnnatSilence=dMean;
		
	 }


	 free(dipStartEnvelope);
    free(dipStopEnvelope);
    free(dipTransitionEnvelope);
	
    free(envelopes);
	 SiStatisticsDelete(&hStatistics);

    return result;
} /* end UnnaturalSilences */
