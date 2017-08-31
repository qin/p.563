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
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "defines.h"
#include "EvalQual.h"
#include "SignalsPercept.h"

#ifndef MAX
#define  MAX(a,b)  (((a) > (b)) ? (a) : (b))
#endif

#ifndef MIN
#define  MIN(a,b)  (((a) < (b)) ? (a) : (b))
#endif


FLOAT gCentreOfBandBark [NUMBER_OF_BARK_BANDS] 
                                        = { 0.078672f,
			                                0.316341f,
			                                0.636559f,
			                                0.961246f,
			                                1.290450f,
			                                1.624217f,
			                                1.962597f,
			                                2.305636f,
			                                2.653383f,
			                                3.005889f,
			                                3.363201f,
			                                3.725371f,
			                                4.092449f,
			                                4.464486f,
			                                4.841533f,
			                                5.223642f,
			                                5.610866f,
			                                6.003256f,
			                                6.400869f,
			                                6.803755f,
			                                7.211971f,
			                                7.625571f,
			                                8.044611f,
			                                8.469146f,
			                                8.899232f,
			                                9.334927f,
			                                9.776288f,
			                                10.223374f,
			                                10.676242f,
			                                11.134952f,
			                                11.599563f,
			                                12.070135f,
			                                12.546731f,
			                                13.029408f,
			                                13.518232f,
			                                14.013264f,
			                                14.514566f,
			                                15.022202f,
			                                15.536238f,
			                                16.056736f,
			                                16.583761f,
                                            17.117382f};

FLOAT gWidthOfBandBark [NUMBER_OF_BARK_BANDS]
                                       = {  0.157344f,
			                                0.317994f,
			                                0.322441f,
			                                0.326934f,
			                                0.331474f,
			                                0.336061f,
			                                0.340697f,
			                                0.345381f,
			                                0.350114f,
			                                0.354897f,
			                                0.359729f,
			                                0.364611f,
			                                0.369544f,
			                                0.374529f,
			                                0.379565f,
			                                0.384653f,
			                                0.389794f,
			                                0.394989f,
			                                0.400236f,
			                                0.405538f,
			                                0.410894f,
			                                0.416306f,
			                                0.421773f,
			                                0.427297f,
			                                0.432877f,
			                                0.438514f,
			                                0.444209f,
			                                0.449962f,
			                                0.455774f,
			                                0.461645f,
			                                0.467577f,
			                                0.473569f,
			                                0.479621f,
			                                0.485736f,
			                                0.491912f,
			                                0.498151f,
			                                0.504454f,
			                                0.510819f,
			                                0.517250f,
			                                0.523745f,
			                                0.530308f,
                                            0.536934f};

FLOAT gPowerDensityCorrectionFactor [NUMBER_OF_BARK_BANDS]
                                         = {100.000000f,
			                                99.999992f,
			                                100.000000f,
			                                100.000008f,
			                                100.000008f,
			                                100.000015f,
			                                99.999992f,
			                                99.999969f,
			                                50.000027f,
			                                100.000000f,
			                                99.999969f,
			                                100.000015f,
			                                99.999947f,
			                                100.000061f,
			                                53.047077f,
			                                110.000046f,
			                                117.991989f,
			                                65.000000f,
			                                68.760147f,
			                                69.999931f,
			                                71.428818f,
			                                75.000038f,
			                                76.843384f,
			                                80.968781f,
			                                88.646126f,
			                                63.864388f,
			                                68.155350f,
			                                72.547775f,
			                                75.584831f,
			                                58.379192f,
			                                80.950836f,
			                                64.135651f,
			                                54.384785f,
			                                73.821884f,
			                                64.437073f,
			                                59.176456f,
			                                65.521278f,
			                                61.399822f,
			                                58.144047f,
			                                57.004543f,
			                                64.126297f,
                                            59.248363f};

INT32 gNumberOfHzBandsInBarkBand [NUMBER_OF_BARK_BANDS] 
                                        = { 1,
			                                1,
			                                1,
			                                1,
			                                1,
			                                1,
			                                1,
			                                1,
			                                2,
			                                1,
			                                1,
			                                1,
			                                1,
			                                1,
			                                2,
			                                1,
			                                1,
			                                2,
			                                2,
			                                2,
			                                2,
			                                2,
			                                2,
			                                2,
			                                2,
			                                3,
			                                3,
			                                3,
			                                3,
			                                4,
			                                3,
			                                4,
			                                5,
			                                4,
			                                5,
			                                6,
			                                6,
			                                7,
			                                8,
			                                9,
			                                9,
                                            11};

FLOAT gAbsoluteThresholdPower [NUMBER_OF_BARK_BANDS] 
                                        = { 407380384.000000f,
						                    19498450.000000f,
						                    562341.437500f,
						                    38904.519531f,
						                    9332.543945f,
						                    3090.295898f,
						                    831.763855f,
						                    363.078094f,
						                    141.253769f,
						                    77.624718f,
						                    38.904518f,
						                    24.547091f,
						                    15.135613f,
						                    10.000000f,
						                    7.762471f,
						                    5.754399f,
						                    4.466836f,
						                    3.630781f,
						                    3.090296f,
						                    2.630268f,
						                    2.344229f,
						                    2.137962f,
						                    2.041738f,
						                    1.995262f,
						                    1.995262f,
						                    1.995262f,
						                    1.995262f,
						                    2.089296f,
						                    2.290868f,
						                    2.454709f,
						                    2.691535f,
						                    2.951209f,
						                    3.162278f,
						                    3.467369f,
						                    3.715352f,
						                    3.890451f,
						                    3.981072f,
						                    3.981072f,
						                    4.073803f,
						                    4.168694f,
						                    4.168694f,
                                            4.168694f};


/*********************
*
* FUNCTION: FrequencyWarpingOf
*
* DESCRIPTION: Mapping of spectrum to Pitch power densities
*                 ( mapp to Bark frequency axis)
*					
*
***********************/

void FrequencyWarpingOf (const FLOAT * pHzSpectrum, FLOAT * pPitchPowerDensity) 
{

    INT32   i, n;
    INT32   hzBandIndex = 0;
    INT32   barkBandIndex;
    FLOAT   sum;
   
    for (barkBandIndex = 0; barkBandIndex < NUMBER_OF_BARK_BANDS; barkBandIndex++) 
	 {
        n = gNumberOfHzBandsInBarkBand [barkBandIndex];
        sum = 0;
        for (i = 0; i < n; i++) {
            sum += pHzSpectrum [hzBandIndex++];
            assert (hzBandIndex <= TRANSFORM_LENGTH / 2);
        }
        pPitchPowerDensity [barkBandIndex] = sum * gPowerDensityCorrectionFactor [barkBandIndex];;
    }
}


/*********************
*
* FUNCTION: TotalAudible
*
* DESCRIPTION: 
*					The following function is used to measure the instaneous total audible power. 
*					A momentary difference between input and output power is reduced by scaling
*					
*
***********************/

FLOAT TotalAudible (const FLOAT * pPitchPowerDensity, FLOAT pFactor) 
{
    
    INT32         bandIndex;
    FLOAT      threshold, result;

    result = 0.;
    for (bandIndex = 1; bandIndex < NUMBER_OF_BARK_BANDS; bandIndex++) 
	 {
        threshold = pFactor * gAbsoluteThresholdPower [bandIndex];
        if (pPitchPowerDensity [bandIndex] > threshold) 
		  {
            result += pPitchPowerDensity [bandIndex];
        }
    }
    return (FLOAT) result;
}

/*********************
*
* FUNCTION: IntensityWarpingOf
*
* DESCRIPTION: 
*					The following function is used to implement the mapping of power to Sone
*				
*	
*
***********************/


#define ZWICKER_POWER       0.23f

void IntensityWarpingOf (const FLOAT * pPitchPowerDensity, FLOAT * pLoudnessDensity) 
{


    INT32         bandIndex;
    FLOAT       h;
    FLOAT       modifiedZwickerPower;

    for (bandIndex = 0; bandIndex < NUMBER_OF_BARK_BANDS; bandIndex++) 
	 {

        if (gCentreOfBandBark [bandIndex] < 7.0f) 
		  {
            h =  9.0f / (gCentreOfBandBark [bandIndex] + 2.0f);
        }
		  else 
		  {
            h = 1.f;
        }
        h = (FLOAT) pow (h, (FLOAT) 0.15f);
        modifiedZwickerPower = ZWICKER_POWER * h;

        if (pPitchPowerDensity [bandIndex] > gAbsoluteThresholdPower [bandIndex])
		  {
            pLoudnessDensity [bandIndex] = (FLOAT) ((pow (gAbsoluteThresholdPower [bandIndex] / 0.5, modifiedZwickerPower)
                    * (pow (0.5f + 0.5f * pPitchPowerDensity [bandIndex] / gAbsoluteThresholdPower [bandIndex], modifiedZwickerPower) - 1)));
        } 
		  else 
		  {
            pLoudnessDensity [bandIndex] = 0;
        }
    }
}

/*********************
*
* FUNCTION: DifferenceOf
*
* DESCRIPTION: 
*					The following function is used to deduct the enhanced bark spectrum and the 
*					distorted bark spectrum
*				
*	
***********************/


void DifferenceOf(const FLOAT * pSignal1, 
                  const FLOAT * pSignal2,  
                  FLOAT *       pResult) 
{ 

    INT32        bandIndex;

    for (bandIndex = 0; bandIndex < NUMBER_OF_BARK_BANDS; bandIndex++) 
	 {
        pResult [bandIndex] = pSignal1 [bandIndex]- pSignal2 [bandIndex];
    }
}

/*********************
*
* FUNCTION: MinimumWeightedOf
*
* DESCRIPTION: 
*					The following function is used to compute the deadzone which determines an 
*					insensitivity of the disturbance
*				
*	
***********************/


void MinimumWeightedOf (FLOAT        pFactor1, 
                        const FLOAT * pSignal1, 
                        FLOAT        pFactor2, 
                        const FLOAT * pSignal2, 
                        FLOAT *       pResult) {

    INT32       bandIndex;

    for (bandIndex = 0; bandIndex < NUMBER_OF_BARK_BANDS; bandIndex++) 
	 {
        FLOAT h1 = pFactor1 * pSignal1 [bandIndex];
        FLOAT h2 = pFactor2 * pSignal2 [bandIndex];

        if (h1 < h2) 
		  {
            pResult [bandIndex] = h1;
        } 
		  else 
		  {
            pResult [bandIndex] = h2;
        }
    }
}


/*********************
*
* FUNCTION: MaxWith
*
* DESCRIPTION: 
*					The following function is used to implement smearing of the masking array
*				
*	
***********************/

void MaxWith (INT32 pShift, FLOAT pFactor1, const FLOAT * pSignal1, FLOAT * pResult) 
{


    INT32       bandIndex;

    for (bandIndex = MAX (-pShift, 0); bandIndex < NUMBER_OF_BARK_BANDS - MAX (pShift, 0); bandIndex++) {
        FLOAT h1 = pFactor1 * pSignal1 [bandIndex + pShift];
        
        pResult [bandIndex] = MAX (pResult[bandIndex], h1);
    }
}


/*********************
*
* FUNCTION: MaxWith
*
* DESCRIPTION: 
*		The following function pulls the signed disturbance values towards zero, thus implementing a deadzone.
*		 As a result, small disturbance values are mapped onto zero.
*				
*	
***********************/


void MaskWith (const FLOAT * pMaskSpectrum, FLOAT * pDisturbance) 
{


    FLOAT  maskLevel;
    INT32    bandIndex;

    for (bandIndex = 0; bandIndex < NUMBER_OF_BARK_BANDS; bandIndex++) 
	 {
        FLOAT h = pDisturbance [bandIndex];
        maskLevel = pMaskSpectrum [bandIndex];

        if (h > maskLevel) 
		  {
            pDisturbance [bandIndex] -= maskLevel;
        } 
		  else 
		  {
            if (h < -maskLevel) 
				{
                pDisturbance [bandIndex] += maskLevel;
            } 
				else 
				{
                pDisturbance [bandIndex] = 0;
            }
        }
    }
}

/*********************
*
* FUNCTION: Maximum
*
* DESCRIPTION: 
*   find maximum value of an array
*			
*	
***********************/


FLOAT Maximum (const FLOAT * pSignal, INT32 pNumberOfSamples) 
{

    FLOAT result = pSignal [0];
    INT32   i;
    FLOAT h;

    for (i = 0; i < pNumberOfSamples; i++) 
	 {
        h = pSignal [i];
        if (h > result) 
		  {
            result = h;
        }
    }

    return result;
}

/*********************
*
* FUNCTION: BarkIntegral
*
* DESCRIPTION: 
*   The following function is used in the calibration procedure at the beginning of the disturbance section
*			
*	
***********************/


FLOAT BarkIntegral (const FLOAT * pSignal) {


    FLOAT result = 0;
    INT32   i;
    
    for (i = 0; i < NUMBER_OF_BARK_BANDS; i++) {
        result += pSignal [i] * gWidthOfBandBark [i];
    }

    return (FLOAT) result;
}

/*********************
*
* FUNCTION: CopyOf
*
* DESCRIPTION: 
*   The following function is used to make of copy of the disturbance array
*			
*	
***********************/

void CopyOf (const FLOAT * pInputSignal, 
             INT32          pNumberOfSamples, 
             FLOAT *       pOutputSignal) 
{

    INT32 i;

    for (i = 0; i < pNumberOfSamples; i++) 
	 {        
        pOutputSignal [i] = pInputSignal [i];
    }
}
    
/*********************
*
* FUNCTION: MultiplyWithAsymmetryFactorAddOf
*
* DESCRIPTION: 
*   The following function is used to make the disturbance values asymmetric
*			
*	
***********************/


void MultiplyWithAsymmetryFactorAddOf (const FLOAT * fpDisortedPitchPowerDensity, 
                                       const FLOAT * fpEnhancedPitchPowerDensity,
                                       FLOAT *       fpDisturbanceDensity) 
{
    INT32   i;
    FLOAT threshold;
    FLOAT ratio, h;

    for (i = 0; i < NUMBER_OF_BARK_BANDS; i++) 
	 {
        threshold = 1 * gAbsoluteThresholdPower [i];
        ratio = (fpEnhancedPitchPowerDensity [i] + 1.f) / (fpDisortedPitchPowerDensity [i] + 1.f); 
        h = (FLOAT) pow (ratio, (FLOAT) 0.23f);  
        
#define MAX_H       5

        if (h > (FLOAT) MAX_H) {h = (FLOAT) MAX_H;} 
        if (h < (FLOAT) 0.5f) {h = (FLOAT) 0.0f;}
        fpDisturbanceDensity [i] *= h;
    }
}

/*********************
*
* FUNCTION: BarkLp
*
* DESCRIPTION: 
*   The following function is used to aggregate the disturbance over the frequency axis
*			
*	
***********************/


FLOAT BarkLp (const FLOAT * pDisturbance, FLOAT pPower) 
{
    INT32    bandIndex;
    FLOAT totalWeight = 0;
    FLOAT result = 0;
    FLOAT  h, w;
    
    for (bandIndex = 1; bandIndex < NUMBER_OF_BARK_BANDS; bandIndex++) 
	 {
        h = (FLOAT) fabs (pDisturbance [bandIndex]);        
        w = gWidthOfBandBark [bandIndex];
        result +=(FLOAT) pow (h * w, pPower);
        totalWeight += w;
    }

    result /= totalWeight;
    result =(FLOAT) pow (result, 1.f/ pPower);
    result *= totalWeight;
    return (FLOAT) result;
}

/*********************
*
* FUNCTION: BarkLp
*
* DESCRIPTION: 
*   The following function is used for the  two-step time aggregation
*			
*	
***********************/

#define NUMBER_OF_PSQM_FRAMES_PER_SplitSecond       20                                                


FLOAT LpqWeight (const FLOAT *            pFrameDisturbance,
                 FLOAT                   pPowerSplitSecond,
                 FLOAT                   pPowerTime,
                 INT32                     pStartFrameIndex,
                 INT32                     pStopFrameIndex) 
{

    FLOAT      resultTime = 0;
    FLOAT      totalTimeWeightTime = 0;
    INT32         startFrameOfSplitSecond;
    
    for (startFrameOfSplitSecond = pStartFrameIndex; 
         startFrameOfSplitSecond <= pStopFrameIndex; 
         startFrameOfSplitSecond += NUMBER_OF_PSQM_FRAMES_PER_SplitSecond/2) 
	{

        FLOAT  resultSplitSecond = 0;
        INT32     countSplitSecond = 0;
        INT32     frameIndex;

        for (frameIndex = startFrameOfSplitSecond;
             frameIndex < startFrameOfSplitSecond + NUMBER_OF_PSQM_FRAMES_PER_SplitSecond;
             frameIndex++) 
			{
            if (frameIndex <= pStopFrameIndex) 
				{
                resultSplitSecond += (FLOAT) pow (pFrameDisturbance [frameIndex], pPowerSplitSecond); 
                countSplitSecond++;                
            }
        }

        resultSplitSecond /= countSplitSecond;
        resultSplitSecond =(FLOAT) pow (resultSplitSecond, (FLOAT) 1/pPowerSplitSecond);        
     
        resultTime += (FLOAT) pow (resultSplitSecond, pPowerTime); 
        totalTimeWeightTime ++;
    }

    resultTime /= totalTimeWeightTime;
    resultTime =(FLOAT) pow (resultTime, (FLOAT) 1 / pPowerTime);

    return (FLOAT) resultTime;
}

/*********************
*
* FUNCTION: FractionInBetween
*
* DESCRIPTION: 
*   The following function for evaluation of the relative energy in a 
*		given frequency range
*			
*	
***********************/

FLOAT FractionInBetween (const FLOAT *  pHzSpectrum,
                         INT32           pLowerFrequencyHz,
                         INT32           pUpperFrequencyHz,
                         INT32           pSampleFrequencyHz) 
{

    FLOAT  frequencyResolution = (FLOAT) pSampleFrequencyHz / (FLOAT) NUMBER_OF_HZ_BANDS;
    FLOAT totalPower = 0;
    FLOAT selectedPower = 0;
    INT32    i;
    FLOAT  fraction;
    

    for (i = 0; i < NUMBER_OF_HZ_BANDS; i++) 
	 {
        FLOAT h = pHzSpectrum [i];
        FLOAT freq = i * frequencyResolution;
        
        totalPower += h;

        if ((freq >= pLowerFrequencyHz) && (freq <= pUpperFrequencyHz)) 
		  {
            selectedPower += h;
        }
    }

    fraction = (FLOAT) (10 * log10 (selectedPower / totalPower));
    return fraction;
}
