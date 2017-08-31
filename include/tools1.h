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


#ifndef MYTOOLS
#define MYTOOLS

void MultiplyWith (FLOAT           *pTimeSeries, 
                   INT32             pNumberOfSamples, 
                   FLOAT			  pScale);

FLOAT LocalPower (const FLOAT  *pTimeSeries,
                  INT32         pStartIndex, 
                  INT32         pStopIndex);

FLOAT GlobalPower (const FLOAT *pTimeSeries,
                   INT32         pNumberOfSamples);

void QuickSortIncreasing (FLOAT *pfThis,
                          INT32 iLeftIndex, 
                          INT32 iRightIndex,
                          INT32 *piThat);

void SortIncreasing (FLOAT * pfThis,
							INT32    pN,
                     INT32   *piThat);



void FftFilter( FLOAT * pfData, 
					INT32 iNrSamples, 
					INT32 iNumberOfPoints,
					const FLOAT Fs, 
					const FLOAT fppFilterCurve [][2] ,
					FLOAT *fpFilterdData);

void StandardIRSFilter(FLOAT *pfSignal, 
						INT32 lSignalLen,
						FLOAT *pfFilterdSignal);

FLOAT  GetFilteredLevel(FLOAT *pfSignal, 
								INT32 lSignalLen);

#endif
