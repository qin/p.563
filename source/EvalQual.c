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



#include <assert.h>
#include <math.h>

#include <stdlib.h>
#include <stdio.h>

#include "defines.h"
#include "tools1.h"
#include "dsp.h"
#include "EvalQual.h"
#include "SignalsPercept.h"

/*********************
*
* FUNCTION: WindowedSample
*
* DESCRIPTION: scale time series using a hann window 
*	
***********************/

static FLOAT WindowedSample (const FLOAT    *pTimeSeries,
										INT32          pStartSample, 
										INT32          pSampleIndex, 
										INT32          pWindowSize) 
{


    static  INT32   storedWindowSize = -1;
    static  FLOAT   *storedWeights = NULL;    
    FLOAT            result;
    INT32              i;
    
    assert ((0 <= pSampleIndex) && (pSampleIndex < pWindowSize));
    
    if (pWindowSize != storedWindowSize) 
	 {
        INT32 sampleIndex;


        if (storedWeights != NULL) 
		  {
            free(storedWeights);
        }
        storedWeights = calloc(pWindowSize + 1,sizeof(FLOAT));

        for (sampleIndex = 0; sampleIndex <= pWindowSize; sampleIndex++) 
		  { 
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
* FUNCTION: PowerSpectrumOf
*
* DESCRIPTION: evaluation of power spectrum of (pFrameIndex * TRANSFORM_LENGTH / 2) 
*					samples taken from beginning of block of samples 
*	
*
***********************/

void PowerSpectrumOf (const FLOAT * pBlockOfSamples, 
                      INT32          pFrameIndex,
                      FLOAT *       pHzSpectrum,
                      FLOAT *       x1) 
{
    
    INT32                 i;
    FLOAT               a, b;
    INT32                 bandIndex;

    for (i = 0; i < TRANSFORM_LENGTH; i++) 
	 {
        x1 [i] = WindowedSample (pBlockOfSamples, pFrameIndex * TRANSFORM_LENGTH / 2, i, TRANSFORM_LENGTH); 

    }

	 {
		INT32 Len=0x01;
		Len=Len << (LOG2_TRANSFORM_LENGTH);
		RealFFT(x1, Len);
	 }
   
    pHzSpectrum[0] = 0;

    for (bandIndex = 1; bandIndex < NUMBER_OF_HZ_BANDS; bandIndex++) 
	 {
        a = x1 [2 * bandIndex];        
        b = x1 [2 * bandIndex + 1];    
        pHzSpectrum [bandIndex] = a * a + b * b; 
    }
}

/*********************
*
* FUNCTION: ComputeBasicVoiceQual
*
* DESCRIPTION: Basic Voice Quality evaluation 
*	
*
***********************/


FLOAT ComputeBasicVoiceQual(FLOAT * fpDistortedSignal,
								    FLOAT * fpEnhancedSignal,
									 INT32   iNrSamples,
									 INT32	iSampSkipBegin,
									 INT32	iSampSkipEnd,
									 FLOAT	fractionActive,
									 FLOAT	*pfABD,
									 FLOAT	*pfACRMOS)
{

    FLOAT   blockPredictionACRMOS;
    FLOAT   fractionLow;
  
    FLOAT   symmetricBlockDisturbance;
    FLOAT   asymmetricBlockDisturbance;

    FLOAT   fDistortedPower;
    FLOAT   fEnhancedPower;
	 INT32     startFrameIndex = iSampSkipBegin / (TRANSFORM_LENGTH / 2);
 	 INT32     stopFrameIndex = (((iNrSamples - iSampSkipEnd) / (TRANSFORM_LENGTH /2)) - 2); 
	 INT32     frameIndex;

	FLOAT   calibrationFactorPower = 1.0f;
	FLOAT   calibrationFactorSone = 1.0f;

	FLOAT *  symmetricFrameDisturbance;
	FLOAT *  asymmetricFrameDisturbance;
	
	FLOAT   s, a;

        {   
			  	FLOAT *  pHelperBufferX1 = calloc(TRANSFORM_LENGTH + 2,sizeof(FLOAT));
            FLOAT *  calibrationTimeSignal = calloc(TRANSFORM_LENGTH,sizeof(FLOAT));
            FLOAT *  calibrationHzSpectrum = calloc(NUMBER_OF_HZ_BANDS,sizeof(FLOAT));
            FLOAT *  calibrationPitchPowerDensity = calloc(NUMBER_OF_BARK_BANDS,sizeof(FLOAT));
            FLOAT *  calibrationLoudnessDensity = calloc(NUMBER_OF_BARK_BANDS,sizeof(FLOAT));
            FLOAT   peak;
            FLOAT   totalSone;

            FLOAT   approxPeriodInSamples = (FLOAT) SAMPLE_FREQUENCY_HZ / (FLOAT) ONE_KILO_HZ;
            INT32     numberOfPeriodsPerFrame = (INT32) floor (0.5f + TRANSFORM_LENGTH
															/ approxPeriodInSamples);
            FLOAT   periodInSamples = (FLOAT) TRANSFORM_LENGTH / (FLOAT) numberOfPeriodsPerFrame;
            FLOAT   omega = 2.0f * (FLOAT) PI / (FLOAT) periodInSamples;
            INT32     i;

            for (i = 0; i < TRANSFORM_LENGTH; i++) 
				{
                calibrationTimeSignal [i] = 29.54f * (FLOAT) sin (i * omega);    
            }

            PowerSpectrumOf (calibrationTimeSignal,
                             0,
                             calibrationHzSpectrum,
                             pHelperBufferX1);
            
            FrequencyWarpingOf (calibrationHzSpectrum, calibrationPitchPowerDensity);

            peak = Maximum (calibrationPitchPowerDensity, NUMBER_OF_BARK_BANDS);

            calibrationFactorPower = 10000.f / peak;

            MultiplyWith (calibrationPitchPowerDensity, NUMBER_OF_BARK_BANDS, calibrationFactorPower);

            IntensityWarpingOf (calibrationPitchPowerDensity, calibrationLoudnessDensity);
            
            totalSone = BarkIntegral (calibrationLoudnessDensity);

            calibrationFactorSone = 1.f/ totalSone;  

            free (calibrationTimeSignal);
            free (calibrationHzSpectrum);
            free (calibrationPitchPowerDensity);
            free (calibrationLoudnessDensity);
				free (pHelperBufferX1);
			
        }

 
            symmetricFrameDisturbance= calloc(stopFrameIndex + 1,sizeof(FLOAT));
            asymmetricFrameDisturbance= calloc(stopFrameIndex + 1,sizeof(FLOAT));



        if (fractionActive > 0.35) 
		  {
            fractionActive = 0.35f;
        }

        { 

			  	FLOAT *  pHelperBufferX1 = calloc(TRANSFORM_LENGTH + 2,sizeof(FLOAT));

            FLOAT *  fpEnhancedHzSpectrum = calloc(NUMBER_OF_HZ_BANDS,sizeof(FLOAT));
            FLOAT *  fpDistortedHzSpectrum = calloc(NUMBER_OF_HZ_BANDS,sizeof(FLOAT));
            FLOAT *  fpEnhancedPitchPowerDensity = calloc(NUMBER_OF_BARK_BANDS,sizeof(FLOAT));
            FLOAT *  fpDistortedPitchPowerDensity = calloc(NUMBER_OF_BARK_BANDS,sizeof(FLOAT));
            FLOAT *  fpEnhancedLoudnessDensity = calloc(NUMBER_OF_BARK_BANDS,sizeof(FLOAT));
            FLOAT *  fpDistortedLoudnessDensity = calloc(NUMBER_OF_BARK_BANDS,sizeof(FLOAT));
            FLOAT *  fpEnhancedPitchPowerDensityAvg = calloc(NUMBER_OF_BARK_BANDS,sizeof(FLOAT));
            FLOAT *  fpDistortedPitchPowerDensityAvg = calloc(NUMBER_OF_BARK_BANDS,sizeof(FLOAT));
            FLOAT *  fpSymDistDensity = calloc(NUMBER_OF_BARK_BANDS,sizeof(FLOAT));
            FLOAT *  fpASymDistDensity = calloc(NUMBER_OF_BARK_BANDS,sizeof(FLOAT));
            FLOAT *  mask = calloc(NUMBER_OF_BARK_BANDS,sizeof(FLOAT));
            
            FLOAT *  avgDistortedHzSpectrum = calloc(NUMBER_OF_HZ_BANDS,sizeof(FLOAT));
            INT32    numberOfDistortedHzSpectraAveraged = 0;
            
            FLOAT  fCurrDistortedTotalPower;
            FLOAT  fCurrEnhancedTotalPower;
            FLOAT  currentScale;
            FLOAT  oldScale = 1.0f;

            FLOAT  multiplier1;
            INT32    i;


            for (i = 0; i < NUMBER_OF_HZ_BANDS; i++) avgDistortedHzSpectrum [i] = 0;
         

            fDistortedPower = GlobalPower (fpDistortedSignal, iNrSamples);
            fEnhancedPower = GlobalPower (fpEnhancedSignal, iNrSamples);        
                       

  
            for (frameIndex = startFrameIndex; frameIndex <= stopFrameIndex; frameIndex++) 
				{                
                FLOAT reduceAvg;
                
 
                PowerSpectrumOf (fpDistortedSignal,
                                 frameIndex,
                                 fpDistortedHzSpectrum,
                                 pHelperBufferX1);
                PowerSpectrumOf (fpEnhancedSignal,
                                 frameIndex,
                                 fpEnhancedHzSpectrum,
                                 pHelperBufferX1);
            
  
                numberOfDistortedHzSpectraAveraged++; 
                reduceAvg = (FLOAT) (numberOfDistortedHzSpectraAveraged - 1)
									/ (FLOAT) numberOfDistortedHzSpectraAveraged;
                
					 for (i = 0; i < NUMBER_OF_HZ_BANDS; i++) 
					 {
                    avgDistortedHzSpectrum [i] *= reduceAvg;
                    avgDistortedHzSpectrum [i] += fpDistortedHzSpectrum [i]
												/ numberOfDistortedHzSpectraAveraged;
                }
                
      
                FrequencyWarpingOf (fpDistortedHzSpectrum, fpDistortedPitchPowerDensity);
                FrequencyWarpingOf (fpEnhancedHzSpectrum, fpEnhancedPitchPowerDensity);

       
                MultiplyWith (fpDistortedPitchPowerDensity, NUMBER_OF_BARK_BANDS, calibrationFactorPower);
                MultiplyWith (fpEnhancedPitchPowerDensity, NUMBER_OF_BARK_BANDS, calibrationFactorPower);

          
                fCurrDistortedTotalPower = TotalAudible (fpDistortedPitchPowerDensity, 100.f); 
                fCurrEnhancedTotalPower = TotalAudible (fpEnhancedPitchPowerDensity, 100.f);  

           
                currentScale = (fCurrDistortedTotalPower + (FLOAT) 1E4) / (fCurrEnhancedTotalPower + (FLOAT) 1E4); 

#define MAX_SCALE   20.0f            

                if (currentScale > (FLOAT) MAX_SCALE) currentScale = (FLOAT) MAX_SCALE;

#define MIN_SCALE   0.10f

                if (currentScale < (FLOAT) MIN_SCALE) currentScale = (FLOAT) MIN_SCALE;            
                
#define ZETA        0.16f
        
                if (frameIndex > startFrameIndex) 
					 {
                    currentScale = (FLOAT) ZETA * oldScale + (FLOAT) (1.0 - ZETA) * currentScale;
                }
                oldScale = currentScale;

     
                MultiplyWith (fpEnhancedPitchPowerDensity, NUMBER_OF_BARK_BANDS, currentScale);                

  
                IntensityWarpingOf (fpDistortedPitchPowerDensity, fpDistortedLoudnessDensity);
                IntensityWarpingOf (fpEnhancedPitchPowerDensity, fpEnhancedLoudnessDensity);

                MultiplyWith (fpDistortedLoudnessDensity, NUMBER_OF_BARK_BANDS, calibrationFactorSone);
                MultiplyWith (fpEnhancedLoudnessDensity, NUMBER_OF_BARK_BANDS, calibrationFactorSone);                
    
 
                DifferenceOf (fpEnhancedLoudnessDensity, fpDistortedLoudnessDensity, fpSymDistDensity);
    

                MinimumWeightedOf (0.05f, fpEnhancedLoudnessDensity, 0.075f, fpDistortedLoudnessDensity, mask);
                MaxWith (-1, 0.040f, fpDistortedLoudnessDensity, mask);   
                MaxWith (-2, 0.025f, fpDistortedLoudnessDensity, mask);   
    
                MaskWith (mask, fpSymDistDensity);
                
                CopyOf (fpSymDistDensity, NUMBER_OF_BARK_BANDS, fpASymDistDensity);
    
 
                MultiplyWithAsymmetryFactorAddOf (fpDistortedPitchPowerDensity, 
                                                  fpEnhancedPitchPowerDensity, 
                                                  fpASymDistDensity);    

                symmetricFrameDisturbance[frameIndex] = BarkLp (fpSymDistDensity, POWER_DISTURBANCE_BARK_SYM );
                asymmetricFrameDisturbance[frameIndex] = BarkLp (fpASymDistDensity, POWER_DISTURBANCE_BARK_ASYM );
              
                
                multiplier1 = (FLOAT) pow ((fCurrDistortedTotalPower + 1E2) / 1E7, 0.20);  
                
                symmetricFrameDisturbance[frameIndex] *= multiplier1;
                

#define MAX_DIST        88
                if (symmetricFrameDisturbance[frameIndex] > MAX_DIST) symmetricFrameDisturbance[frameIndex] = MAX_DIST;
                
					 asymmetricFrameDisturbance[frameIndex] *= multiplier1;            

#define MAX_ADIST       15.6       

                 if (asymmetricFrameDisturbance[frameIndex] > MAX_ADIST) asymmetricFrameDisturbance[frameIndex] = (FLOAT) MAX_ADIST;           
                  
					  asymmetricFrameDisturbance[frameIndex] /= (FLOAT) pow (fractionActive / 0.3, 0.06f); 
               

            }

  
            fractionLow = FractionInBetween (avgDistortedHzSpectrum, 20, 170, SAMPLE_FREQUENCY_HZ);           
            
			
            free(fpEnhancedHzSpectrum);
            free(fpDistortedHzSpectrum);
            free(fpEnhancedPitchPowerDensity);
				free(fpDistortedPitchPowerDensity);
            free(fpEnhancedLoudnessDensity);
            free(fpDistortedLoudnessDensity);
            free(fpEnhancedPitchPowerDensityAvg);
            free(fpDistortedPitchPowerDensityAvg);
            free(fpSymDistDensity);
				free(fpASymDistDensity);
            free(mask);
				free(avgDistortedHzSpectrum);		
				free(pHelperBufferX1);
          
        }  
    



		symmetricBlockDisturbance = LpqWeight (symmetricFrameDisturbance, 
                                   POWER_SPLIT_SECOND_SYM,
                                   POWER_BLOCK_SYM, 
                                   startFrameIndex,
                                   stopFrameIndex);
            
      
       asymmetricBlockDisturbance = LpqWeight (asymmetricFrameDisturbance, 
                                    POWER_SPLIT_SECOND_ASYM,
                                    POWER_BLOCK_ASYM,
                                    startFrameIndex,
                                    stopFrameIndex);
                
	
        s = symmetricBlockDisturbance;
        a = asymmetricBlockDisturbance;


			blockPredictionACRMOS=SBDFACTOR*s+ADBFACTOR*a+FLFACTOR*fractionLow;
			blockPredictionACRMOS=ACRMOSOFFSET + ACRMOSSCALE*blockPredictionACRMOS;

			*pfABD=a;
			*pfACRMOS=blockPredictionACRMOS;

         free(symmetricFrameDisturbance);
         free(asymmetricFrameDisturbance);


	return blockPredictionACRMOS;
}
