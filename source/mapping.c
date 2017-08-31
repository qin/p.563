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

#include "defines.h"
#include "mapping.h"

static FLOAT Max(FLOAT x, FLOAT y)
{
	if (x>y) return x;
	return y;
}

void GetPartitionNumber( int *lPartition, p563Results_struct *ptResults)
{
    *lPartition = 0;

    if (ptResults->tNoise.fSnr <= 15) *lPartition+=2;
	if ((ptResults->tNoise.fEstSegSNR <= 25) && (ptResults->tUnnatural.fLPCCurt > 2) 
		&& (ptResults->tUnnatural.fLPCCurt < 9)) *lPartition+=4;

	if ( (ptResults->tMutes.fMuteLength>0) || (ptResults->tMutes.fSpeechInterruptions>0) || 
			(ptResults->tMutes.fSharpDeclines>0) ) *lPartition+=8;

	if (ptResults->tUnnatural.fRobotisation>0) *lPartition+=16;

	if ( (ptResults->tBasicDesc.fPitchAverage<=50) && (*lPartition==0) ) *lPartition+=1;

	if ((*lPartition) & 2) *lPartition=((*lPartition) & 3);
    if ((*lPartition) & 8) *lPartition=((*lPartition) & 9);
	if ((*lPartition) & 4) *lPartition=((*lPartition) & 5);
    if ((*lPartition) & 16) *lPartition=((*lPartition) & 17);
}

void CalculateOverallMapping(int *PartitionNumber, p563Results_struct *ptResults, FLOAT *PredictedMos)
{
	FLOAT MOV1, MOV2, MOV3, MOV4, MOV5, MOV6, MOV7, MOV8, MOV9, MOV10, MOV11, MOV12, PARTMAP=0;

	switch (*PartitionNumber)
	{		case 0:

			 MOV1 =  Max(0, ptResults->tNoise.fSnr - (FLOAT)15.04100);
			 MOV2 =  Max(0, ptResults->tNoise.fSpectralClarity - (FLOAT)2.83680);
			 MOV3 =  Max(0, ptResults->tBasicDesc.fSpeechSectionLevelVar - (FLOAT)0.00018);
			 MOV4 =  Max(0, ptResults->tUnnatural.fVtpMaxTubeSection - (FLOAT)1.53670);
			 MOV5 =  Max(0, ptResults->tMutes.fUnnaturalSilenceMean - (FLOAT).395720E-08);
			 MOV6 =  Max(0, ptResults->tUnnatural.fLPCSkewAbs - (FLOAT)0.00038);
			 MOV7 =  Max(0, ptResults->tNoise.fHiFreqVar - (FLOAT)0.30711);
			 MOV8 =  Max(0, ptResults->tSpeechExtract.fBasicVoiceQuality - (FLOAT)1.75450);
			 MOV9 =  Max(0, ptResults->tBasicDesc.fPitchAverage - (FLOAT)51.00000);
			 MOV10 = Max(0, ptResults->tUnnatural.fArtAverage - (FLOAT)1.95340);
			 MOV11 = Max(0, ptResults->tUnnatural.fCepADev - (FLOAT)0.03368);
			 MOV12 = Max(0, ptResults->tNoise.fEstBGNoise + (FLOAT)0.24619);
			 
			 PARTMAP = (FLOAT)(2.06246 + 0.01643 * MOV1 + 0.05904 * MOV2 - 0.11457 * MOV3
			             + 0.04918 * MOV4 - 0.02067 * MOV5 + 0.19825 * MOV6
			             - 2.25281 * MOV7 + 0.29361 * MOV8 + 0.01670 * MOV9
			             + 0.08901 * MOV10 - 69.98210 * MOV11 - 1.83748 * MOV12);
			
		break;

		case 1:

			 MOV1 = Max(0, ptResults->tNoise.fSpectralClarity - (FLOAT)3.31100);
			 MOV2 = Max(0, ptResults->tNoise.fSnr - (FLOAT)15.19900);
			 MOV3 = Max(0, ptResults->tUnnatural.fFinalVtpAverage - (FLOAT)0.02293);
			 MOV4 = Max(0, ptResults->tSpeechExtract.fBasicVoiceQuality - (FLOAT)2.60990);
			 MOV5 = Max(0, ptResults->tNoise.fEstSegSNR - (FLOAT)5.40810);
			 MOV6 = Max(0, ptResults->tBasicDesc.fPitchAverage - (FLOAT)27.00000);
			 MOV7 = Max(0, ptResults->tNoise.fEstBGNoise + (FLOAT)0.23730);
			 MOV8 = Max(0, ptResults->tMutes.fUnnaturalSilenceMean + (FLOAT).460037E-09);
			 MOV9 = Max(0, ptResults->tNoise.fLocalBGNoiseMean + (FLOAT)0.00233);
			 MOV10 = Max(0, ptResults->tUnnatural.fLPCSkewAbs - (FLOAT)0.00001);
			 MOV11 = Max(0, ptResults->tUnnatural.fVtpVadOverlap - (FLOAT)0.16647);
			 MOV12 = Max(0, ptResults->tBasicDesc.fSpeechSectionLevelVar + (FLOAT).660578E-08);
			 
			 PARTMAP = (FLOAT)(0.32070 + 0.09140 * MOV1 + 0.02205 * MOV2 + 1.79735 * MOV3
			             + 0.24701 * MOV4 + 0.03224 * MOV5 + 0.02063 * MOV6
			             - 3.06822 * MOV7 - 0.02178 * MOV8 - .514851E-06 * MOV9
			             + 0.30707 * MOV10 - 1.08290 * MOV11 - 0.10127 * MOV12);
			
		break;

		case 16:

			 MOV1 = Max(0, ptResults->tNoise.fSpectralClarity - (FLOAT)4.71610);
			 MOV2 = Max(0, ptResults->tUnnatural.fCepCurt - (FLOAT)24.11000);
			 MOV3 = Max(0, ptResults->tUnnatural.fVtpPeakTracker - (FLOAT)0.19831);
			 MOV4 = Max(0, ptResults->tMutes.fUnnaturalSilenceMean + (FLOAT).301314E-07);
			 MOV5 = Max(0, ptResults->tUnnatural.fLPCCurt - (FLOAT)0.48799);
			 MOV6 = Max(0, ptResults->tUnnatural.fUBeeps + (FLOAT).961321E-12);
			 MOV7 = Max(0, ptResults->tUnnatural.fPitchCrossPower - (FLOAT)2.50520);
			 MOV8 = Max(0, ptResults->tBasicDesc.fSpeechSectionLevelVar - (FLOAT)0.00096);
			 MOV9 = Max(0, ptResults->tUnnatural.fLPCSkew + (FLOAT)0.74012);
			 MOV10 = Max(0, ptResults->tNoise.fHiFreqVar - (FLOAT)0.43271);
			 MOV11 = Max(0, ptResults->tUnnatural.fArtAverage - (FLOAT)1.97970);
			 MOV12 = Max(0, ptResults->tUnnatural.fVtpMaxTubeSection - (FLOAT)1.65670);
			 
			 PARTMAP = (FLOAT)(2.39146 + 0.04377 * MOV1 + 0.03865 * MOV2 - 1.32388 * MOV3
			             - 0.01387 * MOV4 + 0.09580 * MOV5 - 20.17612 * MOV6
			             + 0.03013 * MOV7 - 0.17167 * MOV8 - 0.39843 * MOV9
			             - 1.81732 * MOV10 + 0.15903 * MOV11 - 0.12303 * MOV12);
			
		break;

		case 2:

			 MOV1 = Max(0, ptResults->tUnnatural.fLPCSkewAbs - (FLOAT)0.00005);
			 MOV2 = Max(0, ptResults->tUnnatural.fFinalVtpAverage - (FLOAT)0.06182);
			 MOV3 = Max(0, ptResults->tUnnatural.fPitchCrossCorrelOffset - (FLOAT)1.12850);
			 MOV4 = Max(0, ptResults->tBasicDesc.fPitchAverage - (FLOAT)32.00000);
			 MOV5 = Max(0, ptResults->tUnnatural.fCepADev - (FLOAT)0.04023);
			 MOV6 = Max(0, ptResults->tNoise.fSpecLevelRange - (FLOAT)14.97200);
			 MOV7 = Max(0, ptResults->tNoise.fSpecLevelDev - (FLOAT)6.74700);
			 MOV8 = Max(0, ptResults->tUnnatural.fLPCCurt - (FLOAT)0.07322);
			 MOV9 = Max(0, ptResults->tSpeechExtract.fBasicVoiceQuality - (FLOAT)3.18260);
			 MOV10 = Max(0, ptResults->tUnnatural.fFrameRepeatsMean - (FLOAT)0.00067);
			 MOV11 = Max(0, ptResults->tNoise.fLocalBGNoiseLog - (FLOAT)53.60100);
			 MOV12 = Max(0, ptResults->tBasicDesc.fSpeechLevel + (FLOAT)36.13300);
			 
			 PARTMAP = (FLOAT)(0.83407 + 1.90631 * MOV1 + 2.08569 * MOV2 - 0.06350 * MOV3
			             + 0.02800 * MOV4 + 303.48782 * MOV5 - 0.15274 * MOV6
			             + 0.27536 * MOV7 - 0.14481 * MOV8 + 0.29626 * MOV9
			             + .396050E-05 * MOV10 - 0.04222 * MOV11 + 0.03798 * MOV12);
			
		break;

		case 4:

			 MOV1 = Max(0, ptResults->tNoise.fSpecLevelDev - (FLOAT)5.36230);
			 MOV2 = Max(0, ptResults->tUnnatural.fLPCSkewAbs - (FLOAT)0.00010);
			 MOV3 = Max(0, ptResults->tUnnatural.fPitchCrossPower - (FLOAT)3.52670);
			 MOV4 = Max(0, ptResults->tUnnatural.fLPCSkew + (FLOAT)1.48220);
			 MOV5 = Max(0, ptResults->tNoise.fHiFreqVar - (FLOAT)0.28086);
			 MOV6 = Max(0, ptResults->tUnnatural.fLPCCurt - (FLOAT)2.00530);
			 MOV7 = Max(0, ptResults->tNoise.fRelNoiseFloor - (FLOAT)19.21700);
			 MOV8 = Max(0, ptResults->tUnnatural.fUBeeps + (FLOAT).414795E-11);
			 MOV9 = Max(0, ptResults->tUnnatural.fUnBeepsMeanDistSamp + (FLOAT).406471E-11);
			 MOV10 = Max(0, ptResults->tNoise.fSpectralClarity - (FLOAT)4.80410);
			 MOV11 = Max(0, ptResults->tUnnatural.fVtpMaxTubeSection - (FLOAT)1.33150);
			 MOV12 = Max(0, ptResults->tUnnatural.fCepCurt - (FLOAT)2.19830);
			 
			 PARTMAP = (FLOAT)(0.79317 + 0.67259 * MOV1 + 0.37375 * MOV2 - 0.00782 * MOV3
			             - 0.16638 * MOV4 + 4.67619 * MOV5 + 0.08335 * MOV6
			             + 0.01999 * MOV7 - 87.92210 * MOV8 + 362.28912 * MOV9
			             + 0.01873 * MOV10 + 0.06747 * MOV11 - 0.01138 * MOV12);
			
		break;

		case 8:

			 MOV1 = Max(0, ptResults->tNoise.fSpectralClarity - (FLOAT)3.69380);
			 MOV2 = Max(0, ptResults->tMutes.fMuteLength - (FLOAT).284187E-05);
			 MOV3 = Max(0, ptResults->tSpeechExtract.fBasicVoiceQualityAsym - (FLOAT)0.71223);
			 MOV4 = Max(0, ptResults->tUnnatural.fLPCSkewAbs - (FLOAT)0.00251);
			 MOV5 = Max(0, ptResults->tNoise.fEstBGNoise + (FLOAT)0.19972);
			 MOV6 = Max(0, ptResults->tUnnatural.fArtAverage - (FLOAT)2.25600);
			 MOV7 = Max(0, ptResults->tUnnatural.fUBeepsMean + (FLOAT)0.03415);
			 MOV8 = Max(0, ptResults->tBasicDesc.fSpeechSectionLevelVar - (FLOAT)0.00105);
			 MOV9 = Max(0, ptResults->tNoise.fLocalMeanDistSamp - (FLOAT)0.03353);
			 MOV10 = Max(0, ptResults->tSpeechExtract.fBasicVoiceQuality - (FLOAT)1.86190);
			 MOV11 = Max(0, ptResults->tBasicDesc.fSpeechLevel + (FLOAT)34.96900);
			 MOV12 = Max(0, ptResults->tNoise.fGlobalBGNoise - (FLOAT).301834E-08);
			 
			 PARTMAP = (FLOAT)(2.45541 + 0.08748 * MOV1 - 0.00342 * MOV2 + 0.20826 * MOV3
			             + 0.57269 * MOV4 - 5.80198 * MOV5 + 0.15964 * MOV6
			             - .647241E-08 * MOV7 - 0.19399 * MOV8 - 3.93634 * MOV9
			             + 0.14741 * MOV10 + 0.04647 * MOV11 - 3.02823 * MOV12);
		break;

	}

	 MOV1 = Max(0, PARTMAP - (FLOAT)0.59975);
	 MOV2 = Max(0, ptResults->tBasicDesc.fSpeechLevel + (FLOAT)38.68900);
	 MOV3 = Max(0, ptResults->tNoise.fLocalMeanDistSamp - (FLOAT)0.01735);
	 MOV4 = Max(0, ptResults->tNoise.fGlobalBGNoise + (FLOAT).165375E-08);
	 MOV5 = Max(0, ptResults->tNoise.fLocalBGNoise + (FLOAT).667780E-08);
	 MOV6 = Max(0, ptResults->tSpeechExtract.fBasicVoiceQualityAsym - (FLOAT)0.59026);
	 MOV7 = Max(0, ptResults->tUnnatural.fPitchCrossPower - (FLOAT)2.50520);
	 MOV8 = Max(0, ptResults->tUnnatural.fCepSkew + (FLOAT)0.96454);
	 MOV9 = Max(0, ptResults->tUnnatural.fLPCCurt + (FLOAT)0.02441);
	 MOV10 = Max(0, ptResults->tNoise.fLocalBGNoiseMean - (FLOAT)0.00152);
	 MOV11 = Max(0, ptResults->tUnnatural.fFrameRepeats + (FLOAT).157036E-07);
	 MOV12 = Max(0, ptResults->tUnnatural.fUBeeps + (FLOAT).343165E-10);
	 
	 *PredictedMos = (FLOAT)(1.44818 + 0.94472 * MOV1 + 0.01934 * MOV2 - 2.65333 * MOV3
	             - 5.35764 * MOV4 - 1.84121 * MOV5 + 0.13724 * MOV6
	             - 0.00670 * MOV7 - 0.02235 * MOV8 + 0.01935 * MOV9
	             - .168446E-06 * MOV10 - 0.03036 * MOV11 - 2.23638 * MOV12);

	if (*PredictedMos<1) *PredictedMos =1;
	if (*PredictedMos>5) *PredictedMos =5;
	
}
