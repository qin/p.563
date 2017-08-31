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


#ifndef BEEBROB_H
#define BEEBROB_H

#ifndef FLT_MAX
#define FLT_MAX         (FLOAT)1.e38   /* largest floating point number */ 
#endif
   
/* Define 8 kHz sampling frequency and a 32 ms local FFT window for the disturbance processing */
#define     SAMPLE_FREQUENCY_HZ             8000
#define     TRANSFORM_LENGTH                256
#define     LOG2_TRANSFORM_LENGTH           8

INT16 UnnaturalBeeps(FLOAT *pProcessVector,
                     INT32  pBlockLengthInSamples,
							FLOAT	*fpUnnaturalBeeps,
							FLOAT	*fpMeanUnnaturalBeeps,
							FLOAT	*fpUnBeepsMeanDistSamp
							);


BOOL Robotization(const FLOAT *pProcessVector,
								INT32  pBlockLengthInSamples,
                        FLOAT *pFractionActive,
								FLOAT *pfRobotization);

void FrameRepeats (const FLOAT *  pProcessVector,
                         INT32    pBlockLengthInSamples,
								 FLOAT *pfFrRepNum,
								 FLOAT *pfFrRepRelEnergy);

INT16 SharpDeclines (const FLOAT * pProcessVector,
									INT32  pBlockLengthInSamples,
									FLOAT *fpSharpDeclines);
                         


INT32 UnnaturalSilences (const FLOAT *   pOrigDistBlockOfSamples,
									INT32      pBlockLengthInSamples,
									FLOAT	*MeanUnatSilence
								 );                      

#endif
