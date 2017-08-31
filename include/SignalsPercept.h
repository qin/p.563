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

#ifndef SIGNALS_PERCEPT_H
#define SIGNALS_PERCEPT_H

#define NUMBER_OF_BARK_BANDS        42
#define NUMBER_OF_HZ_BANDS          ((TRANSFORM_LENGTH / 2) + 1)

void    FrequencyWarpingOf (const FLOAT * pHzSpectrum, 
									 FLOAT * pPitchPowerDensity);                

FLOAT   TotalAudible (const FLOAT * pPitchPowerDensity, 
							 FLOAT pFactor);

void    IntensityWarpingOf (const FLOAT * pPitchPowerDensity, 
									 FLOAT * pLoudnessDensity);

void    DifferenceOf (const FLOAT * pSignal1, 
                      const FLOAT * pSignal2, 
                      FLOAT *       pResult);
 
void    MinimumWeightedOf (FLOAT        pFactor1, 
                           const FLOAT * pSignal1, 
                           FLOAT        pFactor2, 
                           const FLOAT * pSignal2,
                           FLOAT * pResult);


void    MaxWith (INT32 pShift, 
					  FLOAT pFactor1, 
					  const FLOAT * pSignal1, 
					  FLOAT * pResult);

void    MaskWith (const FLOAT * pMaskSpectrum, 
						FLOAT * pDisturbance);

FLOAT   Maximum (const FLOAT * pSignal,
					  INT32 pNumberOfSamples);

FLOAT   BarkIntegral (const FLOAT * pSignal);
    
void    CopyOf (const FLOAT * pInputSignal, 
					 INT32 pNumberOfSamples, 
					 FLOAT * pOutputSignal);

void    MultiplyWithAsymmetryFactorAddOf (const FLOAT *fpDistortedPitchPowerDensity, 
                                          const FLOAT *fpEnhancedPitchPowerDensity,
                                          FLOAT *fpDisturbanceDensity);

FLOAT   BarkLp (const FLOAT * pDisturbance, FLOAT pPower);

FLOAT   LpqWeight (const FLOAT *            pFrameDisturbance,
                   FLOAT                   pPowerSyllabe,
                   FLOAT                   pPowerTime,
                   INT32                     pStartFrameIndex,
                   INT32                     pStopFrameIndex);

FLOAT   FractionInBetween (const FLOAT *  pHzSpectrum,
                          INT32           pLowerFrequencyHz,
                          INT32           pUpperFrequencyHz,
                          INT32           pSampleFrequencyHz);


#endif
