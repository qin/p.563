/********************************************************************
ITU-T Draft Recommendation P.563
Version 1.0 - 23 March 2004
  
NOTICE 
  
The Single Ended Assessment Model P.563 algorithm and the copyright therein 
is the joint property of Psytechnics Limited, OPTICOM GmbH and SwissQual AG 
and is protected by UK, US and other patents, either applied for or registered. 
Permission is granted to use this source code solely for the purpose of 
evaluation of ITU-T recommendation P.563. 
Any other use of this software requires a licence, which may be obtained from: 
  
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
  
Psytechnics, SwissQual or Opticom can provide licences and further information. 
  
Authors: 
      Ludovic Malfait ludovic.malfait@psytechnics.com 
      Roland Bitto rb@opticom.de 
      Pero Juric pero.juric@swissqual.com
 
********************************************************************/


#ifndef MODULE1_H
	#define MODULE1_H

	void module1(const INT16 *psRawSpeech, const UINT32 ulRawSpeechSz, p563Results_struct *tResults);

	
	void vtp_extract(INT16 average_pitch, INT32 * pitch, INT32 pitch_length, 
					 INT16 * speech, INT32 speech_length, FLOAT * vtp[13], INT32 * Nframes);

	FLOAT FFT_gaps(const INT16 *speech, const INT32 speech_Nelements, const INT32 *final_pitch, const INT32 final_pitch_Nelements);

	FLOAT Tract_Averages(FLOAT *VTP[8], INT32 VTP_row, INT32 VTP_tubes, INT32 index);
	
	void averages(const FLOAT *in, const UINT32 in_Nelements, FLOAT *out);

	void VTPtoART(FLOAT *VTP[8], const UINT32 VTP_Nelements, FLOAT *ART[3]);

	FLOAT ConsistentArtTracker(FLOAT *ART[3], INT32 ART_Nelements);

	void Correlations_Zero_Cross( INT16 *speech, INT32 *pitch, INT32 pitch_Nelements, 
								  FLOAT *PtoMean_SA_out, FLOAT *corr_offset_out);

	FLOAT VTPPeakTracker(FLOAT **VTP, INT32 VTP_Nelements, FLOAT *pfPeakDiffOut, FLOAT *pfFinalVtpOut);

	void DC_remove_float(const FLOAT *spa_Speech, FLOAT *spa_DCRemoved, UINT32 l_SpeechNelements);
	void DC_remove_int16(const INT16 *spa_Speech, FLOAT *spa_DCRemoved, UINT32 l_SpeechNelements);

	INT32 SNR(const FLOAT *pfSignal, const INT32 lSignalLen, const FLOAT *pfVad, 
				 const FLOAT fSpeechActivity, FLOAT *pfSpeechLevel, FLOAT *pfNoiseLevel, FLOAT *pfSnr);

	void vad_process( INT16 * X, INT32 Nsamples, FLOAT * VAD,
		FLOAT *LevelThresh, FLOAT *LevelSpeech, FLOAT *LevelNoise, FLOAT *Activity, INT32 Downsample );

	void normalise(INT16 *spa_NormSpeech, UINT32 l_NormSpeechNelements, const FLOAT f_SpeechLevel);

	FLOAT SpeechSectionLevelVar(const FLOAT *VAD, const UINT32 l_VADNelements);

	FLOAT hi_freq_variations(const INT16 *in, const INT32 in_Nelements, const FLOAT *pfVad, const INT32 lVadSz);

	UINT32 mute_spotter(const INT16 *psRawSpeech, const INT32 psRawSpeechSz);

	FLOAT VtpVadOverlap(FLOAT *LeftRight[2], INT32 VTP_Nelements, FLOAT *VAD, INT32 VAD_Nelements, 
			INT32 speech_Nelements);
#endif
