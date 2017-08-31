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


#include "defines.h"


#include "back_noise.h"
#include "Statistics.h"
#include "dsp.h"
#include "tools1.h"


#define MAX_NUMBER_OF_VOWELS            1000
#define SUBSAMPLE_FACTOR_ENVELOPES      160 
  

/*********************
*
* FUNCTION: LocalBackgroundNoise
*
* DESCRIPTION:
*	Local BackgroundNoise evaluation
*
***********************/

void LocalBackgroundNoise (const FLOAT *pProcessVector,
									INT32  iNrSamples,
									INT32  iSampSkipStart,
									INT32  iSampSkipEnd,
									FLOAT *pFractionZeroedNoVowels,

									FLOAT *pfLocalBackNoise,
									FLOAT *pfLocalMeanDistSamp,
									FLOAT *pfLocalBackMean,
									FLOAT *pfLocalBackLog
								)
{

    INT32     numberOfEnvelopes = iNrSamples / SUBSAMPLE_FACTOR_ENVELOPES;
    FLOAT *  envelopes = calloc(numberOfEnvelopes,sizeof(FLOAT));
    INT16 *madeZero = (INT16 *) calloc(iNrSamples, sizeof (INT16)); 
    FLOAT *  vowelStartTime = calloc(MAX_NUMBER_OF_VOWELS,sizeof(FLOAT));
    FLOAT *  vowelStopTime = calloc(MAX_NUMBER_OF_VOWELS,sizeof(FLOAT));
    INT32     envelopeIndex;
    INT32     numberOfSamplesZeroed = 0;
    FLOAT   endTime;
    FLOAT   offset;
    INT32    sample;
	 INT32     numberOfVowels = 0;


	SISTATISTIC_HANDLE hStatistics =NULL;
	SiStatisticsCreate(&hStatistics,"LocalBackgroundNoise ");



    for (sample = 0; sample < iNrSamples; sample++) 
	 {
        madeZero [sample] = FALSE;
    }

    for (envelopeIndex = 0; envelopeIndex < numberOfEnvelopes; envelopeIndex++)
	 {
        FLOAT sum = 0;
        INT32    i;
        INT32    sampleIndex = envelopeIndex * SUBSAMPLE_FACTOR_ENVELOPES;
            
        for (i = 0; i < SUBSAMPLE_FACTOR_ENVELOPES; i++) {
            FLOAT h = pProcessVector [sampleIndex++];    
            sum += h * h;
        }

        sum /= SUBSAMPLE_FACTOR_ENVELOPES;

        envelopes [envelopeIndex] = (FLOAT) sqrt (sum);

    }

#define     MIN_ENV_DURING_SPEECH   700.f  


    for (envelopeIndex = 0; envelopeIndex < numberOfEnvelopes; envelopeIndex++) 
	 {
        envelopes [envelopeIndex] += MIN_ENV_DURING_SPEECH;
    }


    for (envelopeIndex = 0; envelopeIndex < numberOfEnvelopes; envelopeIndex++)
	 {


#define MIN_PEAK_GAIN               2.2f
#define MAX_TIME_TO_PEAK_SECS       0.1f

        FLOAT currentEnvelope = envelopes [envelopeIndex];
        INT32 peakIndex = envelopeIndex + 1;
        INT32 maxEnvelopesToPeak = (INT32) (SAMPLE_FREQUENCY * MAX_TIME_TO_PEAK_SECS / SUBSAMPLE_FACTOR_ENVELOPES);
        INT32 envelopesToPeak = 0;
        FLOAT maxEnvelope = envelopes [0];

        while ((peakIndex < numberOfEnvelopes) 
                && (envelopes [peakIndex] < MIN_PEAK_GAIN * currentEnvelope) 
                && (envelopesToPeak < maxEnvelopesToPeak)) {
            peakIndex++;
            envelopesToPeak++;
            if ((peakIndex < numberOfEnvelopes) && (envelopes [peakIndex] > maxEnvelope)) 
				{
                maxEnvelope = envelopes [peakIndex];
            }
        }

#define MAX_VOWEL_DURATION_SECS       0.4f

        if ((peakIndex < numberOfEnvelopes) && (envelopes [peakIndex] >= MIN_PEAK_GAIN * currentEnvelope)) 
		  {

            INT32 maxEnvelopesInVowel = (INT32) (SAMPLE_FREQUENCY * MAX_VOWEL_DURATION_SECS / SUBSAMPLE_FACTOR_ENVELOPES);
            INT32 envelopesInVowel = envelopesToPeak + 1;
            
            FLOAT maxRecedeGain = (FLOAT)sqrt (maxEnvelope / currentEnvelope);

            INT32 recedeIndex = peakIndex + 1;
            if ((recedeIndex < numberOfEnvelopes) && (envelopes [recedeIndex] > maxEnvelope)) 
				{
                maxEnvelope = envelopes [recedeIndex];
            }
            
            while ((recedeIndex < numberOfEnvelopes) 
                    && (envelopes [recedeIndex] > maxRecedeGain * currentEnvelope) 
                    && (envelopesInVowel < maxEnvelopesInVowel)) {
                recedeIndex++;
                envelopesInVowel++;
                if ((recedeIndex < numberOfEnvelopes) && (envelopes [recedeIndex] > maxEnvelope)) 
					 {
                    maxEnvelope = envelopes [recedeIndex];
                }                
                maxRecedeGain =(FLOAT)sqrt (maxEnvelope / currentEnvelope);
            }
            
            if ((recedeIndex < numberOfEnvelopes) && (envelopes [recedeIndex] <= maxRecedeGain * currentEnvelope)) 
				{
          
					
                FLOAT startTime = envelopeIndex * SUBSAMPLE_FACTOR_ENVELOPES / (FLOAT) SAMPLE_FREQUENCY; 
                FLOAT stopTime = recedeIndex * SUBSAMPLE_FACTOR_ENVELOPES / (FLOAT) SAMPLE_FREQUENCY; 
                
              
                vowelStartTime [numberOfVowels] = startTime;
                vowelStopTime [numberOfVowels] = stopTime;                         
                numberOfVowels++;

                if (numberOfVowels >= MAX_NUMBER_OF_VOWELS) 
					 {

                    exit (1);
                }

                envelopeIndex = recedeIndex;
            }
        }
    }        

 
    endTime = (numberOfEnvelopes - 0.1f) * SUBSAMPLE_FACTOR_ENVELOPES / (FLOAT) SAMPLE_FREQUENCY;

#define DURATION_INTERVAL           1.f

    
    for (offset = 0; offset < 0.99 * DURATION_INTERVAL; offset += DURATION_INTERVAL / 4) {
        FLOAT startTimeInterval;
        for (startTimeInterval = offset; 
             startTimeInterval + DURATION_INTERVAL < endTime; 
             startTimeInterval += DURATION_INTERVAL) 
				 {
            INT32 count = 0;
            INT32 i;

            for (i = 0; i < numberOfVowels; i++) {
                FLOAT t1 = vowelStartTime [i];
                FLOAT t2 = vowelStopTime [i];
                if ((t1 >= startTimeInterval) && (t1 < startTimeInterval + DURATION_INTERVAL)) {
                    count++;
                }
                if ((t2 >= startTimeInterval) && (t2 < startTimeInterval + DURATION_INTERVAL)) {
                    count++;
                }
            }
        
           
#define MIN_COUNT                   2
#define DURATION_SILENCE            0.2f 
             
            if (count < 2 * MIN_COUNT) 
				{

          
                INT32 startEnvelopeInterval = (INT32) (startTimeInterval * SAMPLE_FREQUENCY / SUBSAMPLE_FACTOR_ENVELOPES);
                INT32 stopEnvelopeInterval = (INT32) ((startTimeInterval + DURATION_INTERVAL)  * SAMPLE_FREQUENCY / SUBSAMPLE_FACTOR_ENVELOPES);
                INT32 numberOfEnvelopesInSilence = (INT32) (DURATION_SILENCE * SAMPLE_FREQUENCY / SUBSAMPLE_FACTOR_ENVELOPES);
            
                FLOAT minRMS = (FLOAT)1E38;
                INT32 bestStartEnvelopeSilence = -1;
                INT32 startEnvelopeSilence;

                for (startEnvelopeSilence = startEnvelopeInterval;
                     startEnvelopeSilence + numberOfEnvelopesInSilence < stopEnvelopeInterval;
                     startEnvelopeSilence++) {

                    FLOAT sum = 0;
                    INT32    i;
                    FLOAT  rms;
                    
                    for (i = startEnvelopeSilence; i < startEnvelopeSilence + numberOfEnvelopesInSilence; i++) {
                        FLOAT h = envelopes [i];
                        sum += h * h;
                    }
                    sum /= (numberOfEnvelopesInSilence + 1);

                    rms = (FLOAT) sqrt (sum);

                    if (rms < minRMS) {
                        minRMS = rms;
                        bestStartEnvelopeSilence = startEnvelopeSilence;
                    }
                }

                if (envelopes [bestStartEnvelopeSilence] > 0) 
					 {
                    INT32 startSampleSilence = bestStartEnvelopeSilence * SUBSAMPLE_FACTOR_ENVELOPES;
                    INT32 stopSampleSilence = startSampleSilence + numberOfEnvelopesInSilence * SUBSAMPLE_FACTOR_ENVELOPES;
                    INT32 sample;

						  for (sample = startSampleSilence; sample < stopSampleSilence; sample++) 
						  {
							  madeZero [sample] = TRUE;
                    }

						  if(hStatistics != NULL)
						  {
						
								SiStatisticsSetDefaultVector(hStatistics,pProcessVector,startSampleSilence,stopSampleSilence,1);
						  }

                }

           
                for (i = bestStartEnvelopeSilence;
                     i < bestStartEnvelopeSilence + numberOfEnvelopesInSilence;
                     i++) {
                    envelopes [i] = 0;
                }                    
            }
        }
    }

    for (sample = iSampSkipStart; sample < iNrSamples - iSampSkipEnd; sample++) 
	 {
        if (madeZero [sample]) {
            numberOfSamplesZeroed++;
        }
    }
    (*pFractionZeroedNoVowels) = (FLOAT) numberOfSamplesZeroed / (FLOAT) (iNrSamples - iSampSkipStart - iSampSkipEnd);


	{
		 FLOAT dMean=0;
		 FLOAT dMov=0;
		 FLOAT dStd=0;
		 INT32 lNrEntries=0;

		 *pfLocalBackNoise=0;
		if (*pFractionZeroedNoVowels > 0){*pfLocalBackNoise = *pFractionZeroedNoVowels;}

		SiStatisticsGetMoments(hStatistics,&dMean,&dStd,&lNrEntries);

		*pfLocalMeanDistSamp=(FLOAT)lNrEntries/(FLOAT)(iNrSamples-iSampSkipStart-iSampSkipEnd);
		*pfLocalBackMean=dMean;
			
		dMov=dMean;
		dMov=(FLOAT)(fabs(dMov)+1.0);
		dMov=(FLOAT)(10*log10(dMov));
		*pfLocalBackLog=dMov;

	 }


    free(envelopes);
    free(vowelStartTime);
    free(vowelStopTime);    
    free(madeZero);

	 SiStatisticsDelete(&hStatistics);
}

/*********************
*
* FUNCTION: GlobalBackgroundNoise
*
* DESCRIPTION:
*	Global BackgroundNoise evaluation
*
***********************/


void GlobalBackgroundNoise(const FLOAT *pProcessVector,
									const FLOAT *pEvalVector,
									INT32     iNrSamples,
									FLOAT     pFractionZeroedNoVowels,
							
								  FLOAT *pfGlobalBackNoise
					  
							   )

{   

    INT32     windowLength = TRANSFORM_LENGTH;
    INT32     lWindowNum = (INT32) (2 * iNrSamples / windowLength - 1);
    FLOAT  *envelope = calloc(lWindowNum,sizeof(FLOAT));

    INT32     maxDistanceToTheRightInWindows = (INT32) ceil ((14 * TRANSFORM_LENGTH) / windowLength);
    INT32     maxDistanceToTheLeftInWindows = (INT32) ceil ((14 * TRANSFORM_LENGTH) / windowLength);
    FLOAT  *smearedEnvelope = calloc(lWindowNum,sizeof(FLOAT));
    INT32   *sortedk = (INT32 *) calloc(lWindowNum , sizeof (INT32));
    INT32     k;
    INT32     firstActiveWindow = 0;
    INT32     lastActiveWindow = lWindowNum - 1;
    FLOAT   fractionToBeSilenced; 
    INT32     numberOfToBeSilencedWindows;
    INT32     windowsNoBreath;

	 FLOAT	fTotalEnergy;

	if(pEvalVector==0)pEvalVector=pProcessVector;

    fTotalEnergy = GlobalPower (pEvalVector, iNrSamples);


    for (k = 0; k < lWindowNum; k++) 
	 {
        INT32 startSample = k * windowLength / 2;
        INT32 stopSample = startSample + windowLength;
        envelope [k] = (FLOAT) sqrt (LocalPower (pProcessVector, startSample, stopSample));
        smearedEnvelope [k] = envelope [k];
    }

    while ((firstActiveWindow < lastActiveWindow) && (envelope [firstActiveWindow] < 200)) 
	 {
        firstActiveWindow++;
    }
    while ((firstActiveWindow < lastActiveWindow) && (envelope [lastActiveWindow] < 200)) 
	 {
        lastActiveWindow--;
    }
    if (lastActiveWindow < lWindowNum - 1) 
	 {
        lastActiveWindow++;
    }
    

    for (k = 0; k < lWindowNum; k++)
	 {
        FLOAT h = envelope [k];
        INT32 lowk = k - maxDistanceToTheLeftInWindows;
        INT32 highk = k + maxDistanceToTheRightInWindows;
        INT32 i;
        if (lowk < 0) {
            lowk = 0;
        }
        if (highk >= lWindowNum) {
            highk = lWindowNum - 1;
        }
        sortedk [k] = (INT32) k;
        for (i = lowk; i <= highk; i++) 
		  {
        
            if ((smearedEnvelope [i] > 0) && (smearedEnvelope [i] < h)) {
                smearedEnvelope [i] = h;
            }
        }
    }

   
    QuickSortIncreasing (smearedEnvelope, firstActiveWindow, lastActiveWindow, sortedk);

#define SILENT_TAIL_FRACTION                    0.07f
#define DO_NOT_NEED_TO_TAKE_A_BREATH_YET_SECS   1.f    

    pFractionZeroedNoVowels -= 0.05f;
    if (pFractionZeroedNoVowels < 0) {
        pFractionZeroedNoVowels = 0;
    }

    fractionToBeSilenced = SILENT_TAIL_FRACTION - 0.5f * pFractionZeroedNoVowels;
    windowsNoBreath = (INT32) (DO_NOT_NEED_TO_TAKE_A_BREATH_YET_SECS * SAMPLE_FREQUENCY) / (TRANSFORM_LENGTH / 2);
    numberOfToBeSilencedWindows = (INT32) (fractionToBeSilenced *  (lastActiveWindow - firstActiveWindow + 1 - windowsNoBreath));
    if (numberOfToBeSilencedWindows < 0) {
        numberOfToBeSilencedWindows = 0;
    }

	{

		*pfGlobalBackNoise=0;
		if (numberOfToBeSilencedWindows > 1)
		{
			*pfGlobalBackNoise = (FLOAT) ((numberOfToBeSilencedWindows - 1) * windowLength * 1.5) / (FLOAT) iNrSamples;
		}
				
	 }


    free(envelope);
    free(smearedEnvelope);
    free(sortedk);

} 



