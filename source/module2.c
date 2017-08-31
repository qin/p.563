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
#include <memory.h>
#include <limits.h>
#include <math.h>

#include "defines.h"
#include "dsp.h"
#include "module2.h"
#include "tools1.h"
#include "Enhance.h"

#include "back_noise.h"
#include "beeprob.h"
#include "EvalQual.h"



/*********************
* FUNCTION: get_StartStop
*
* DESCRIPTION:
*   Compute a row estimation of the sart and stop point of speech in signal
*    
***********************/

static void get_StartStop(FLOAT *pfSignal,INT32 lSignalLen,INT32 *plStartSignal,INT32 *plStopSignal)
{
	const INT32 lWinSize=5;
	const FLOAT fSilenceLevel=100.0F;


	FLOAT fSum=0;
	INT32 lStart=0;
	INT32 lStop=0;
	  
	lStart=0;
    do {INT32 i;
        fSum = 0.0f;
        for (i = 0; i < lWinSize; i++) {
            fSum += (FLOAT) fabs (pfSignal [lStart + i]);
        }
        if (fSum < lWinSize * fSilenceLevel) {
			lStart++;         
        }        
    } while ((fSum < lWinSize * fSilenceLevel) && (lStart < lSignalLen));
    


	lStop=0;
    do {INT32 i;
        fSum = 0.0f;
        for (i = 0; i < lWinSize; i++) {
            fSum += (FLOAT) fabs (pfSignal[lSignalLen - 1 - lStop - i]);
        }
        if (fSum < lWinSize * fSilenceLevel) {
            lStop++;         
        }        
    } while ((fSum < lWinSize * fSilenceLevel) && (lStop < lSignalLen));


	*plStartSignal=lStart;
	*plStopSignal=lSignalLen-lStop;
  
}

/*********************
* FUNCTION: AlignAndFilter
*
* DESCRIPTION:
*   LevelAlignement and IRS filter
*    
*
***********************/

#define     TARGET_AVG_POWER                1E7 
static void AlignAndFilter(FLOAT *pfSignal,INT32 lSignalLen)
{

		FLOAT fLevel;
		INT32 i;
		
		fLevel=GetFilteredLevel(pfSignal,lSignalLen);
		  
		fLevel += TARGET_AVG_POWER / 1E3;
		fLevel = (FLOAT)sqrt(TARGET_AVG_POWER/fLevel);

		for(i=0;i<lSignalLen;i++)
		{
			pfSignal[i]=pfSignal[i]*fLevel;
		}

		StandardIRSFilter(pfSignal, lSignalLen,pfSignal);

}

/*********************
* FUNCTION: module2float
*
* DESCRIPTION: Main distribution function for evalutation of following MOVS:
*
*	 fLocalBGNoise,
*	LocalMeanDistSamp
*	LocalBGNoiseMean
*	LocalBGNoiseLog
*	GlobalBGNoise
*	UBeeps
*	UBeepsMean
*	UnBeepsMeanDistSamp
*	Robotisation
*	FrameRepeats
*	FrameRepeatsMean
*	UnnaturalSilenceMean
*	SharpDeclines
*	BasicVoiceQualityAsym
*	BasicVoiceQuality
*
***********************/
static void module2float(FLOAT *fpBufferSignal, UINT32 iSignalLen, p563Results_struct *ptResults)
{

	FLOAT *fpBufferEnhance=NULL;
	INT32 lDelay=0;
	INT32 iBufferLen=0;
	INT32 i;

	INT32 lStartSignal=0;
	INT32 lStopSignal=0;

	FLOAT fZeroFraction=0;
   FLOAT fFractionActive=0;	



	fpBufferEnhance =(FLOAT *)calloc(iSignalLen,sizeof(FLOAT));


	 SpeechEnhaceFilter(fpBufferSignal,fpBufferEnhance,iSignalLen,&iBufferLen,&lDelay);

    for (i = 0; i < iBufferLen - lDelay; i++) {
            fpBufferEnhance[i] = fpBufferEnhance[i + lDelay];
    }
    iBufferLen -= lDelay;


	get_StartStop(fpBufferEnhance,iBufferLen,&lStartSignal,&lStopSignal);

	AlignAndFilter(fpBufferSignal,iBufferLen);
	AlignAndFilter(fpBufferEnhance,iBufferLen);

	{


   LocalBackgroundNoise (fpBufferSignal,
								 iBufferLen,                    
								 lStartSignal, 
								 max(0,iBufferLen-lStopSignal),

                         &fZeroFraction,
								 &ptResults->tNoise.fLocalBGNoise,
								 &ptResults->tNoise.fLocalMeanDistSamp,
								 &ptResults->tNoise.fLocalBGNoiseMean,
								 &ptResults->tNoise.fLocalBGNoiseLog
								);



	  GlobalBackgroundNoise(	fpBufferEnhance,
							  			fpBufferSignal,
										iBufferLen,
										fZeroFraction,
										&ptResults->tNoise.fGlobalBGNoise			  
									);


	  UnnaturalBeeps(fpBufferSignal,
							iBufferLen,
							&(ptResults->tUnnatural.fUBeeps),	
							&(ptResults->tUnnatural.fUBeepsMean),
							&(ptResults->tUnnatural.fUnBeepsMeanDistSamp)
						);


	  Robotization( fpBufferSignal,
						 iBufferLen,
					    &fFractionActive,
						 &(ptResults->tUnnatural.fRobotisation)
					 );

	  FrameRepeats( fpBufferSignal,
						 iBufferLen,
							&(ptResults->tUnnatural.fFrameRepeats),
							&(ptResults->tUnnatural.fFrameRepeatsMean)
						);


	   UnnaturalSilences (fpBufferSignal,
							   iBufferLen,
								&ptResults->tMutes.fUnnaturalSilenceMean);

	   SharpDeclines(fpBufferSignal,
					  	  iBufferLen,
						  &ptResults->tMutes.fSharpDeclines);



	ComputeBasicVoiceQual(fpBufferSignal,
								 fpBufferEnhance,
								 iBufferLen,
								 lStartSignal, 
								 max(0,iBufferLen-lStopSignal),
								fFractionActive,
								&ptResults->tSpeechExtract.fBasicVoiceQualityAsym,
								&ptResults->tSpeechExtract.fBasicVoiceQuality);


	}






	free(fpBufferEnhance);



}


/*********************
* FUNCTION: module2
*
* DESCRIPTION: Main interface module.
*					Transfer INT16RawSpeech signals to FLOAT format
*
*					calculation of theese MOVS
*							 fLocalBGNoise,
*							LocalMeanDistSamp
*							LocalBGNoiseMean
*							LocalBGNoiseLog
*							GlobalBGNoise
*							UBeeps
*							UBeepsMean
*							UnBeepsMeanDistSamp
*							Robotisation
*							FrameRepeats
*							FrameRepeatsMean
*							UnnaturalSilenceMean
*							SharpDeclines
*							BasicVoiceQualityAsym
*							BasicVoiceQuality
*
***********************/


void module2(const INT16 *psRawSpeech, const UINT32 ulRawSpeechSz, p563Results_struct *tResults)
{

	UINT32 i;
   FLOAT *pfDistSpeech=NULL;

	pfDistSpeech = calloc(ulRawSpeechSz,sizeof(FLOAT));

		
	for (i = 0; i < ulRawSpeechSz; i++) 
	{
           pfDistSpeech[i] = (FLOAT)psRawSpeech[i];
	}

	module2float(pfDistSpeech,ulRawSpeechSz,tResults);
		
	free(pfDistSpeech);

}

