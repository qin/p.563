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

#include <stdlib.h>
#include <memory.h>
#include <limits.h>
#include <math.h>

#include "defines.h"
#include "dsp.h"
#include "pitch.h"
#include "tools.h"
#include "lpc.h"
#include "module1.h"


/*********************
*
*	FUNCTION: module1
*
*	DESCRIPTION:
*		This module calculates the following parameters:
*			fPitchAverage
*			fSpeechLevel
*			fSpeechSectionLevelVar
*			fMuteLength
*			fHiFreqVar
*			fNoiseLevel
*			fSnr
*			fSpectralClarity
*			fArtAverage
*			fConsistentArtTracker
*			fFinalVtpAverage
*			fPitchCrossCorrelOffset
*			fPitchCrossPower
*			fVtpMaxTubeSection
*			fVtpPeakTracker
*			fVtpVadOverlap
*
***********************/
void module1(const INT16 *psRawSpeech, const UINT32 ulRawSpeechSz, p563Results_struct *tResults)
{
    INT16   *psNormSpeech = NULL;
	FLOAT   *pfNormSpeech = NULL;
    UINT32 lNormSpeechSz;
    FLOAT   *pfVad = NULL;
    INT32	lVadSz;
    FLOAT	fSpeechActivity = 0;
	FLOAT	fNoiseLevel = 0;
    FLOAT	fSpeechLevel = 0;
	FLOAT	fSpeechLevel1 = 0;
    FLOAT	fLevelThresh = 0;
	FLOAT	fSnr = 0;
    FLOAT	*pfNormSpeechNoDc = NULL;

	FLOAT	*ppdVtp[VTP_ORDER+2] ;
	FLOAT	*ppdArt[3];

	INT32 *plPitch=NULL;
	INT32 lPitchSz;
	INT16 sPitchAv=0;

	UINT32 i;

	lNormSpeechSz = ulRawSpeechSz + 192;
	lVadSz = lNormSpeechSz / VAD_FRAME_LEN;
    psNormSpeech = (INT16 *) calloc(lNormSpeechSz,2);
	pfNormSpeech = (FLOAT *) malloc(lNormSpeechSz * sizeof(FLOAT));
    memcpy(psNormSpeech+96,psRawSpeech, 2*ulRawSpeechSz);
    pfVad = (FLOAT *) calloc(lVadSz , sizeof(FLOAT));
    pfNormSpeechNoDc = (FLOAT *) malloc(lNormSpeechSz * sizeof(FLOAT));
	plPitch = (INT32 *) malloc( (lNormSpeechSz/2) * sizeof(INT32));

	vad_process(psNormSpeech,lNormSpeechSz, pfVad, &fLevelThresh,
			&fSpeechLevel, &fNoiseLevel, &fSpeechActivity,VAD_FRAME_LEN);

    DC_remove_int16(psNormSpeech,pfNormSpeechNoDc,lNormSpeechSz);
    normalise(psNormSpeech,lNormSpeechSz,fSpeechLevel);
	for (i=0;i<lNormSpeechSz;i++) pfNormSpeech[i] = psNormSpeech[i];

	SNR(pfNormSpeechNoDc,lNormSpeechSz, pfVad, fSpeechActivity,&fSpeechLevel,&fNoiseLevel,&fSnr);

	vad_process(psNormSpeech,lNormSpeechSz, pfVad, &fLevelThresh, &fSpeechLevel1, &fNoiseLevel,
		&fSpeechActivity,VAD_FRAME_LEN);
    SNR(pfNormSpeech, lNormSpeechSz, pfVad, fSpeechActivity, &fSpeechLevel1, &fNoiseLevel, &fSnr);

	pitch_extract(psNormSpeech,lNormSpeechSz,pfVad, lVadSz, plPitch, &lPitchSz,lNormSpeechSz/2,
		&sPitchAv);

	for (i=0; i<(VTP_ORDER+2); i++)
		ppdVtp[i] = (FLOAT  *) calloc((INT32 )(lPitchSz * 1.5), sizeof(FLOAT));

	for (i=0; i<3; i++)
		ppdArt[i] = (FLOAT *) malloc (lPitchSz * sizeof(FLOAT));

	if (lPitchSz>0)
	{
		vtp_extract(sPitchAv, plPitch, lPitchSz, psNormSpeech, lNormSpeechSz, ppdVtp, &lPitchSz );

		VTPtoART(ppdVtp, lPitchSz, ppdArt);

		tResults->tNoise.fSpectralClarity = FFT_gaps(psNormSpeech, lNormSpeechSz, plPitch, lPitchSz);

		tResults->tUnnatural.fConsistentArtTracker = ConsistentArtTracker(ppdArt, lPitchSz);

		Correlations_Zero_Cross( psNormSpeech, plPitch, lPitchSz,
							  &(tResults->tUnnatural.fPitchCrossPower),
							  &(tResults->tUnnatural.fPitchCrossCorrelOffset));

		VTPPeakTracker(ppdVtp, lPitchSz, &(tResults->tUnnatural.fVtpPeakTracker),
					&(tResults->tUnnatural.fFinalVtpAverage));

		tResults->tUnnatural.fVtpVadOverlap =
			VtpVadOverlap(&ppdVtp[VTP_ORDER],lPitchSz,pfVad,lVadSz,lNormSpeechSz);

		tResults->tUnnatural.fVtpMaxTubeSection = Tract_Averages(ppdVtp, lPitchSz, VTP_NB_TUBES, 3);

		tResults->tUnnatural.fArtAverage = Tract_Averages(ppdArt, lPitchSz, 3, 5);

		tResults->tBasicDesc.fPitchAverage = sPitchAv;
	}
	else
	{
		tResults->tNoise.fSpectralClarity=0;
		tResults->tUnnatural.fConsistentArtTracker=0;
		tResults->tUnnatural.fPitchCrossPower=0;
		tResults->tUnnatural.fPitchCrossCorrelOffset=0;
		tResults->tUnnatural.fFinalVtpAverage=0;
		tResults->tBasicDesc.fPitchAverage=0;
		tResults->tUnnatural.fVtpVadOverlap=0;
		tResults->tUnnatural.fVtpMaxTubeSection=0;
		tResults->tUnnatural.fArtAverage=0;
	}

	tResults->tBasicDesc.fSpeechLevel = fSpeechLevel;
	tResults->tNoise.fNoiseLevel= fSpeechLevel - fSnr;
	tResults->tNoise.fSnr = fSnr;

	tResults->tBasicDesc.fSpeechSectionLevelVar = SpeechSectionLevelVar(pfVad,lVadSz);

	tResults->tNoise.fHiFreqVar = hi_freq_variations(psNormSpeech, lNormSpeechSz,pfVad, lVadSz);

	tResults->tMutes.fMuteLength = (FLOAT) mute_spotter(psRawSpeech,ulRawSpeechSz);




	for (i=0; i<(VTP_ORDER+2); i++) free(ppdVtp[i]) ;
	for (i=0; i<3; i++) free(ppdArt[i]) ;
	free(plPitch);
	free(pfVad);
	free(psNormSpeech);
	free(pfNormSpeech);
	free(pfNormSpeechNoDc);
    FFTFree();
}


/*********************
*
*	FUNCTION: vtp_extract
*
*	DESCRIPTION:
*		This function calculates the vocal tract parameters for all pitch sections
*		in the speech file being processed. The pitch marks are first broken up
*		into voiced sections and then within these sections they are paired before
*		presentation to the coefficient generation routine. Parameters are extracted
*		from the reflection coefficients.
*
***********************/
void vtp_extract(INT16 average_pitch, INT32 * pitch, INT32 pitch_length,
				 INT16 * speech, INT32 speech_length, FLOAT * vtp[13], INT32 * Nframes)
{

	INT32 * first_array ;
	INT32 * last_array ;
	INT32 * frames ;
	FLOAT * double_speech ;
	FLOAT * coefs ;
	INT32 * left_array ;
	INT32 * right_array ;
	FLOAT frame[VTP_ORDER] ;
	INT32 frames_length = 0;
	INT32 first_last_length = 0;
	INT32 pitch_index = 0 ;
	INT32 build_index = 0 ;

	INT32 section_index ;
	INT32 cycle_length ;
	INT32 frame_index = 0 ;
	INT32 loop ;
	INT32 frame_loop ;

		first_array = (INT32 *) calloc(pitch_length, sizeof(INT32)) ;
		last_array = (INT32 *) calloc(pitch_length, sizeof(INT32)) ;
		frames = (INT32 *) calloc((pitch_length * 2), sizeof(INT32)) ;
		double_speech = (FLOAT *) calloc(speech_length, sizeof(FLOAT)) ;
		coefs = (FLOAT *) calloc(((*Nframes) * (VTP_ORDER + 0)), sizeof(FLOAT)) ;
		left_array = (INT32 *) calloc(*Nframes, sizeof(INT32)) ;
		right_array = (INT32 *) calloc(*Nframes, sizeof(INT32)) ;

		find_pos_neg_first_last(average_pitch, pitch, pitch_length, first_array, last_array, &first_last_length) ;

		for(section_index = 0; section_index < first_last_length; section_index++)
		{

			while(first_array[section_index] != pitch[pitch_index])
			{
				pitch_index++ ;
			}
			while(pitch[pitch_index] < last_array[section_index])
			{
				cycle_length = pitch[pitch_index + 1] - pitch[pitch_index] ;
				left_array[frame_index] = pitch[pitch_index] ;
				right_array[frame_index] = pitch[pitch_index + 1] ;
				frame_index++ ;
				frames[frames_length] = pitch[pitch_index] - (INT32)(cycle_length/2) ;
				frames_length++ ;
				frames[frames_length] = pitch[pitch_index + 1] + (INT32)(cycle_length/2) ;
				frames_length++ ;
				pitch_index++ ;
			}

			pitch_index++ ;
		}

		*Nframes = frame_index;

		for(loop = 0; loop < speech_length; loop++)
		{
			double_speech[loop] = (FLOAT)speech[loop] ;
		}

		calc_gen_coef(double_speech, frames, *Nframes, coefs, VTP_ORDER) ;

		section_index = 0 ;

		for(frame_index =0; frame_index < *Nframes; frame_index++)
		{
			for(frame_loop = 0; frame_loop < VTP_ORDER; frame_loop++)
			{
				frame[frame_loop] = coefs[(frame_index * VTP_ORDER) + frame_loop] ;
			}

			reverse_array(frame, VTP_ORDER, FLOAT_TYPE) ;
			for(loop = 0; loop < VTP_ORDER; loop++)
			{
				vtp[loop][build_index] = frame[loop] ;
			}
			vtp[VTP_ORDER][build_index] = (FLOAT) left_array[frame_index] ;
			vtp[VTP_ORDER+1][build_index] = (FLOAT) right_array[frame_index] ;
			build_index++ ;

			if(right_array[frame_index] == last_array[section_index])
			{
				for(loop = 0; loop < (VTP_ORDER+2); loop++)
				{
					vtp[loop][build_index] = -9999;
				}
				build_index++ ;
				section_index++ ;
			}

		}

		*Nframes = build_index ;

		free(first_array) ;
		free(last_array) ;
		free(double_speech) ;
		free(left_array) ;
		free(right_array) ;
		free(frames) ;
		free(coefs) ;

}

/*********************
*
*	FUNCTION: FFT_gaps
*
*	DESCRIPTION:
*		Energy analysis of harmonics in pitched sections, centered on each pitch mark.
*		The difference between the first harmonics and gaps energies
*		is calculated and averaged.
*
***********************/
FLOAT FFT_gaps(const INT16 *speech, const INT32 speech_Nelements, const INT32 *final_pitch, const INT32 final_pitch_Nelements)
{
	const FFT_FRAME_SIZE = 512;
	const HALF_FFT_FRAME_SIZE = 256;
	const NFORMANTS = 4;
	FLOAT *FFT_result;
	UINT32 last_FFT_pos = 0;
	UINT32 pos_temp = 0, pos_pitch=0;
	UINT32 loop = 0, loop1 = 0;
	INT32 pitch_length = 0;
	FLOAT formant_nrg = 0, gap_nrg= 0;
	INT32 spec_unit = 0,half_spec_unit = 0, quarter_spec_unit = 0;
	FLOAT mean_diff_nrg=0;

	FLOAT *blackman_window = NULL;

	FFT_result = (FLOAT *) malloc ( (FFT_FRAME_SIZE + 2) * sizeof(FLOAT) );
	blackman_window = (FLOAT *) malloc (FFT_FRAME_SIZE * sizeof(FLOAT) );

	blackmanHarris_window(blackman_window, FFT_FRAME_SIZE);

	last_FFT_pos = speech_Nelements;

	for (pos_pitch=1; pos_pitch<(UINT32)(final_pitch_Nelements-1); pos_pitch++)
	{
		pos_temp = (UINT32) round2(( ((FLOAT) final_pitch[pos_pitch]) / HALF_FFT_FRAME_SIZE) );
		if ((speech_Nelements - (INT32)(pos_temp*HALF_FFT_FRAME_SIZE))<FFT_FRAME_SIZE)
			pos_temp = (speech_Nelements - FFT_FRAME_SIZE-1)/HALF_FFT_FRAME_SIZE;

		if ((last_FFT_pos != pos_temp))
		{
			loop1 = pos_temp*HALF_FFT_FRAME_SIZE;
			for (loop = 0; loop<(UINT32) FFT_FRAME_SIZE; loop ++,loop1++)
				FFT_result[loop] = ((FLOAT) speech[loop1]) * blackman_window[loop];

			RealFFT(FFT_result,FFT_FRAME_SIZE);

			for (loop = 0,loop1=0; loop<(UINT32) FFT_FRAME_SIZE; loop+=2,loop1++)
			{
				FFT_result[loop1] = (FLOAT) sqrt(FFT_result[loop] * FFT_result[loop] + FFT_result[loop+1] * FFT_result[loop+1]);
				if (FFT_result[loop1] <1) FFT_result[loop1] = 1;
				FFT_result[loop1] = (FLOAT) (20*log10(FFT_result[loop1]));
			}
			last_FFT_pos = pos_temp;
		}

		pitch_length = min(final_pitch[pos_pitch] - final_pitch[pos_pitch-1],
			final_pitch[pos_pitch+1] - final_pitch[pos_pitch]);

		spec_unit = (INT32) round2((FFT_FRAME_SIZE / (FLOAT) pitch_length));
		half_spec_unit = round2( (FLOAT) spec_unit /2);
		quarter_spec_unit = round2( (FLOAT) spec_unit /4);

		formant_nrg=0;
		gap_nrg=0;

		if (((NFORMANTS+1)*spec_unit+quarter_spec_unit) > FFT_FRAME_SIZE)
		{
			for (loop=1; loop<(UINT32) (NFORMANTS); loop++)
				for (loop1=(loop*spec_unit - quarter_spec_unit);loop1<(loop*spec_unit + quarter_spec_unit);loop1++ )
				{
					formant_nrg+=FFT_result[loop1];
					gap_nrg+=FFT_result[loop1 - half_spec_unit];
				}
			mean_diff_nrg += (formant_nrg - gap_nrg) / ((NFORMANTS-1) * half_spec_unit);
		}
		else
		{
			for (loop=1; loop<(UINT32) (NFORMANTS+1); loop++)
				for (loop1=(loop*spec_unit - quarter_spec_unit);loop1<(loop*spec_unit + quarter_spec_unit);loop1++ )
				{
					formant_nrg+=FFT_result[loop1];
					gap_nrg+=FFT_result[loop1 - half_spec_unit];
				}
			mean_diff_nrg += (formant_nrg - gap_nrg) / (NFORMANTS * half_spec_unit);
		}
	}
	free( FFT_result);
	free(blackman_window);
	return (FLOAT) (mean_diff_nrg/(final_pitch_Nelements-2));
}

/*********************
*
*	FUNCTION: Tract_Averages
*
*	DESCRIPTION:
*		Generate statistic measures from the tract parameters.
*
***********************/
FLOAT Tract_Averages(FLOAT *VTP[8], INT32 VTP_row, INT32 VTP_tubes, INT32 index)
{

	INT32 i,j,cpt = 0;
	INT32 rmv=0;

	FLOAT *temp;

	FLOAT VTP_av[105];

	for (i=0; i<VTP_row; i++)
	{
		if (VTP[0][i]==-9999) rmv++;
		else for (j=0;j<VTP_tubes;j++) VTP[j][i-rmv] = VTP[j][i];
	}
	VTP_row-=rmv;
	for (i=0; i<VTP_tubes; i++)
	{
		cpt = 10*i;
		VTP_av[cpt+3]=VTP[i][0];

		VTP_av[cpt+1] = (FLOAT) sqrt(variance(VTP[i],VTP_row, FLOAT_TYPE, &VTP_av[cpt]));
		VTP_av[cpt+2] = ((1000000)*VTP_av[cpt])/VTP_av[cpt+1];
		for (j=0; j<VTP_row;j++) if (VTP_av[cpt+3]<VTP[i][j]) VTP_av[cpt+3]=VTP[i][j];
		VTP_av[cpt+4] = VTP_av[cpt+3]/VTP_av[cpt];
		averages(VTP[i], VTP_row, &VTP_av[cpt+5]);
	}

	cpt+=10;
	temp = (FLOAT *) malloc(VTP_tubes * sizeof(FLOAT));

	for (i=0;i<5;i++,cpt+=5)
	{
		for (j=0;j<VTP_tubes;j++) temp[j] = VTP_av[5+i+10*j];

		averages(temp,VTP_tubes,&VTP_av[cpt]);
	}

	free(temp);
	return VTP_av[index];
}

/*********************
*
*	FUNCTION: averages
*
*	DESCRIPTION:
*		Calculate power averages of the input array.
*
***********************/
void averages(const FLOAT *in, const UINT32 in_Nelements, FLOAT *out)
{
	UINT32 cpt;
	FLOAT calc;

	for (cpt = 0; cpt<5; cpt++) out[cpt]=0;

	if (in_Nelements)
	{
		for (cpt = 0; cpt<in_Nelements; cpt++)
		{
			out[0]+=(FLOAT) sqrt(in[cpt]);
			out[1]+=in[cpt];
			calc=in[cpt];
			calc*=calc;
			out[2]+=calc;
			calc*=calc;
			out[3]+=calc;
			calc*=calc;
			out[4]+=calc;
		}

		for (cpt = 0; cpt<5; cpt++) out[cpt]/=in_Nelements;

		out[0] *= out[0];
		out[2] = (FLOAT) sqrt(out[2]);
		out[3] = (FLOAT) pow(out[3],0.25);
		out[4] = (FLOAT) pow(out[4],0.125);
	}
}

/*********************
*
*	FUNCTION: VTPtoART
*
*	DESCRIPTION:
*		Combine vocal tract parameters to generate articulator parameters.
*
***********************/
void VTPtoART(FLOAT *VTP[8], const UINT32 VTP_Nelements, FLOAT *ART[3])
{
	UINT32 i;

	for(i=0;i<VTP_Nelements;i++)
	{
		if (VTP[0][i]!=-9999)
		{
			ART[0][i] = VTP[0][i] + VTP[1][i] + VTP[2][i];
			ART[1][i] = VTP[3][i] + VTP[4][i] + VTP[5][i];
			ART[2][i] = VTP[6][i] + VTP[7][i];
		}
		else
		{
			ART[0][i] = -9999;
			ART[1][i] = -9999;
			ART[2][i] = -9999;
		}
	}
}

/*********************
*
*	FUNCTION: ConsistentArtTracker
*
*	DESCRIPTION:
*		For each pitch frame, the difference in section areas between
*		the back and middle cavity is calculated.
*		The length of smooth sections in the generated array is then calculated.
*		Sections are considered to be consistent when consecutive elements do
*		not vary by more than �0.25. The length of each section is then averaged
*		over the number of sections found and over the number of pitch cycles.
*
***********************/
FLOAT ConsistentArtTracker(FLOAT *ART[3], INT32 ART_Nelements)
{

	FLOAT *diff, *deriv,count1=0;
	INT32 diff_Nelements=0;
	INT32 i,cpt=0,count=0;
	char test=0;
	FLOAT output;

	diff = (FLOAT *) malloc(ART_Nelements * sizeof(FLOAT));

	for (i=0;i<ART_Nelements;i++)
	{
		diff[diff_Nelements]=(ART[0][i]-ART[1][i]);

		if ((diff[diff_Nelements])!=0)
		{
			if (diff[diff_Nelements]>0) diff[diff_Nelements] = (FLOAT)((INT32)(100*diff[diff_Nelements]+0.5))/100;
			else diff[diff_Nelements] = (FLOAT)((INT32)(100*diff[diff_Nelements]-0.5))/100;
			diff_Nelements++;
		}
	}

	deriv = (FLOAT *) malloc(diff_Nelements * sizeof(FLOAT));

	derivative(diff, diff_Nelements, deriv);

	for (i=0;i<diff_Nelements;i++)
	{
		if ((deriv[i]<=0.25) && (deriv[i]>-0.25))
		{
			if (test) count++;
			else
			{
				count=1;
				test=1;
			}
		}
		else if (test)
		{
			deriv[cpt] = (FLOAT) count;
			test = 0;
			cpt++;
		}
	}

	if (cpt)
	{
		for (i=0;i<cpt;i++) count1+= deriv[i];
		output = count1 / cpt;
		output*= (count1 / diff_Nelements);
	}
	else
	{
		output=0;
	}
	free(diff);
	free(deriv);

	return output;
}

/*********************
*
*	FUNCTION: Correlations_Zero_Cross
*
*	DESCRIPTION:
*		This function generates two parameters:
*		PitchCrossPower:
*			Cross power between 2 frames for each consecutive frame.
*			The peak to mean ratio is calculated for each cross power computation.
*			The peak to mean array is then averaged.
*
*		PitchCorrelationOffset:
*			Cross-correlation of consecutive frames.
*			For each pair of frames, the distance between the position of
*			the maximum value and the middle of the frame is calculated.
*			The averaged distance is then calculated.
*
***********************/
void Correlations_Zero_Cross( INT16 *speech, INT32 *pitch, INT32 pitch_Nelements,
							  FLOAT *PtoMean_SA_out, FLOAT *corr_offset_out)

{
	INT32 i=0,j=0,k=0,l=0;
	INT32 *first_last[2], fl_Nelements=0;
	INT32 section[2], section_Nelements[2];
	INT32 global_count=-1;
	FLOAT *cpwr=NULL;
	INT32 cpwr_Nelements=0;
	FLOAT max=0,sum=0;
	FLOAT *PtoMean=NULL, *maxi=NULL;
	INT32 PtoMean_Nelements=0, maxi_Nelements=0;
	FLOAT *correl=NULL;
	INT32 correl_Nelements = 0;
	FLOAT *correl_offset=NULL;
	INT32 correl_offset_Nelements;

	First_Last(pitch, &pitch_Nelements, first_last, &fl_Nelements);

	maxi			= (FLOAT *) calloc( (pitch_Nelements-2),sizeof(FLOAT) );
	PtoMean			= (FLOAT *) calloc( (pitch_Nelements-2),sizeof(FLOAT) );
	correl_offset	= (FLOAT *) calloc( (pitch_Nelements-2),sizeof(FLOAT) );

	for (i=0; i<fl_Nelements; i++)
	{
		for (j=0; j<first_last[1][i]-first_last[0][i]-2; j++)
		{
			global_count++;

			for (k = first_last[0][i],l=0; l<2; k++,l++)
			{
				section_Nelements[l]	= pitch[ k+j+1 ] - pitch[ k+j ];
				section[l]				= pitch[ k+j ] - (section_Nelements[l]>>1);
				section_Nelements[l]  <<= 1;
			}

			cross_power(&speech[section[0]], section_Nelements[0], &speech[section[1]], section_Nelements[1],
						&cpwr, &cpwr_Nelements);

			for ( k=1, sum=max=cpwr[0]; k<cpwr_Nelements; k++)
			{
				if (max<cpwr[k]) max = cpwr[k];
				sum+=cpwr[k];
			}
			sum-=max;

			if (max && sum) PtoMean[global_count]= max / (sum/(cpwr_Nelements-1));
			maxi[global_count]= max;

			if (j && (j!=(first_last[1][i]-first_last[0][i]-3)))
			{
				correl = cross_correl(&speech[section[0]], section_Nelements[0], &speech[section[1]], section_Nelements[1],
					&correl_Nelements);
				for (k=0, max=correl[0];k<correl_Nelements;k++)
						if (max<correl[k])
						{
							max=correl[k];
							correl_offset[global_count] = (FLOAT) k;
						}
				correl_offset[global_count] = (FLOAT)fabs(correl_offset[global_count] - ( (FLOAT)correl_Nelements)/2);
			}
			else
			{
				correl_offset[global_count] =0;
			}
			if (correl) free(correl);
			if (cpwr) free(cpwr);
			cpwr = NULL;
			correl= NULL;
		}
	}

	maxi_Nelements = correl_offset_Nelements = PtoMean_Nelements = pitch_Nelements-2;
	zero_removerf(PtoMean,&PtoMean_Nelements);
	zero_removerf(maxi,&maxi_Nelements);
	zero_removerf(correl_offset,&correl_offset_Nelements);

	*PtoMean_SA_out = smooth_averagesf(PtoMean, PtoMean_Nelements, 10);

	*corr_offset_out =0;
	for (i=0; i< correl_offset_Nelements; i++)
		*corr_offset_out +=(FLOAT) sqrt(correl_offset[i]);
	*corr_offset_out /=correl_offset_Nelements;
	*corr_offset_out *=(*corr_offset_out);

	free(first_last[0]);
	free(first_last[1]);
	free(correl_offset);
	free(PtoMean);
	free(maxi);
	first_last[0] = first_last[1] = NULL;
	PtoMean = maxi = NULL;

}

/*********************
*
*	FUNCTION: Peak_tracker
*
*	DESCRIPTION:
*		This function generates  the VTPPeakTracker parameter.
*		It looks for amplitude variations within the vocal tract array.
*		For each frame, the largest of the eight tubes is searched.
*		The derivative of the position array in calculated to determine the variations.
*		The variation array is then averaged.
*
***********************/
FLOAT VTPPeakTracker(FLOAT **VTP, INT32 VTP_Nelements, FLOAT *fpPeakDiffOut, FLOAT *fpFinalVtpOut)
{
	INT32 i,j,k;
	FLOAT *temp, *temp1;
	FLOAT calc;

	temp  = (FLOAT*) malloc (VTP_Nelements*sizeof(FLOAT));
	temp1 = (FLOAT*) malloc (VTP_Nelements*sizeof(FLOAT));

	for (i=0,k=0;i<VTP_Nelements;i++)
	{
		if (VTP[0][i]!=-9999)
		{
			for (j=1, temp[i-k]=0; j<8; j++)
				if ( VTP[(INT32) temp[i-k] ][i] < VTP[j][i] ) temp[i-k]=(FLOAT) j;
		}
		else
		{
			k++;
		}
	}

	derivative(temp, VTP_Nelements-k,temp1);

	for (j=0 ; j<VTP_Nelements-k; j++) temp1[j] = (FLOAT) fabs(temp1[j]);

	*fpPeakDiffOut = 0;
	for (i=0; i<(VTP_Nelements-k); i++) *fpPeakDiffOut += temp1[i];
	*fpPeakDiffOut/=(VTP_Nelements-k);

	for (i=0,k=0;i<VTP_Nelements;i++)
	{
		if (VTP[0][i]!=-9999)
		{
			for (j=1, temp[i-k]=(VTP[1][i]-VTP[0][i]); j<8; j++)
			{
				calc = (FLOAT) fabs(VTP[j][i]-VTP[j-1][i]);

				if ( temp[i-k] < calc) temp[i-k]=calc;
			}
			if (temp[i-k]==0) k++;
		}
		else k++;
	}
	for (i=0,k=0;i<VTP_Nelements;i++)
	{
		if ((VTP[0][i]!=-9999) && (VTP[0][i])) temp[i-k]=VTP[8][i]/10;
		else k++;
	}

	*fpFinalVtpOut=0;
	for (i=0; i<(VTP_Nelements-k); i++) *fpFinalVtpOut+=temp[i];
	*fpFinalVtpOut/=VTP_Nelements-k;

	free(temp);
	free(temp1);

	return 0;
}


/*********************
*
*	FUNCTION: VtpVadOverlap
*
*	DESCRIPTION:
*		This function calculates the ratio of
*		the total length of voiced sections within speech sections
*		over the total length of speech sections.
*
***********************/
FLOAT VtpVadOverlap(FLOAT *LeftRight[2], INT32 VTP_Nelements, FLOAT *VAD, INT32 VAD_Nelements,
		INT32 speech_Nelements)
{
	INT32	VIO_Nelements=0;
	FLOAT	*VIO[2];
	INT32	i=0, k;

	INT16	*vio_bool ;
	INT16	*vad_bool ;
	INT32	vio_index ;
	INT32	place_loop ;
	INT32	vad_index ;
	INT32	bool_index ;
	INT32	in_total = 0 ;
	INT32	vad_total = 0 ;

	for (i=0; i<2; i++) VIO[i] = (FLOAT*) malloc(VTP_Nelements * sizeof(FLOAT));
	vio_bool = (INT16 *) calloc(speech_Nelements, sizeof(INT16)) ;
	vad_bool = (INT16 *) calloc(speech_Nelements, sizeof(INT16)) ;

	for (i=0,k=0; i<VTP_Nelements; i++)
	{
		if (LeftRight[0][i]!=-9999)
		{
			VIO[0][i-k]=LeftRight[0][i];
			VIO[1][i-k]=LeftRight[1][i];
		}
		else k++;
	};

	VIO_Nelements = i-k;

	for(vio_index = 0; vio_index < VIO_Nelements; vio_index++)
	{
		for(place_loop = (INT32)VIO[0][vio_index]; place_loop < VIO[1][vio_index]; place_loop++)
		{
			vio_bool[place_loop] = 1 ;
		}
	}

	for(vad_index = 0; vad_index < VAD_Nelements; vad_index++)
	{
		for(place_loop = 0; place_loop < VAD_FRAME_LEN; place_loop++)
		{
			if(VAD[vad_index] != 0)
			{
				vad_bool[(vad_index*VAD_FRAME_LEN) + place_loop] = 1 ;
			}
			else
			{
				vad_bool[(vad_index*VAD_FRAME_LEN) + place_loop] = 0 ;
			}
		}
	}

	for(bool_index = 0; bool_index < speech_Nelements; bool_index++)
	{
		if(vad_bool[bool_index]) vad_total++ ;
		if(vad_bool[bool_index] && vio_bool[bool_index]) in_total++ ;
	}

	free(vio_bool);
	free(vad_bool);

	for (i=0; i<2; i++) free(VIO[i]) ;

	return (FLOAT)in_total/(FLOAT)vad_total;
}


/*********************
*
*	FUNCTION: SNR
*
*	DESCRIPTION:
*		Noise level, speech level and SNR.
*		Two different method is used depending on the speech activity.
*		If the speech activity is lower or equal to 80%, then the level of
*		noise sections and speech sections are calculated using a modified
*		version of the VAD profile.
*		If the speech activity is greate than 80%, then a statistical
*		approach is used. The level of each VAD frame is calculated,
*		and the average of the lowest RMS values is retained to be the
*		noise level. The average of all other frames gives the speech level.
*
*		In both case, the SNR is the difference betzeen the calculated
*		speech level and noise level.
*
***********************/
INT32 SNR(const FLOAT *pfSignal, const INT32 lSignalLen, const FLOAT *pfVad,
             const FLOAT fSpeechActivity, FLOAT *pfSpeechLevel, FLOAT *pfNoiseLevel, FLOAT *pfSnr)
{
	UINT32 i=0, l_StartNelements =0, l_StopNelements = 0;
	UINT32 l_NoiseLength = 0, l_SpeechLength = 0;
	UINT32 *lpa_Start = NULL, *lpa_Stop = NULL, *lpa_Basket=NULL;
    const UINT32 l_VADNelements = lSignalLen / VAD_FRAME_LEN;
	const UINT32 l_MaxNbSections = l_VADNelements / 2;
	FLOAT *fpa_RMSValues = NULL;
	FLOAT f_MaxRMS = 0, f_Ratio = 0, f_Stepm;

	FLOAT f_RMSNoise=0,f_RMSSpeech=0, f_RMSSum=0, f_RMSSum1 = 0;
	if (fSpeechActivity>0.80)
	{
		lpa_Basket = (UINT32 *) calloc(5001, sizeof(UINT32));
		fpa_RMSValues = (FLOAT *) malloc(l_VADNelements * sizeof(FLOAT));

		l_NoiseLength  = (UINT32) (l_VADNelements * 0.05) ;
		l_SpeechLength = l_VADNelements - l_NoiseLength;
		for (i=0 ; i<l_VADNelements ; i++)
		{
			fpa_RMSValues[i] = (FLOAT) rms( &pfSignal[i*VAD_FRAME_LEN], VAD_FRAME_LEN, FLOAT_TYPE);
			if (fpa_RMSValues[i]>f_MaxRMS) f_MaxRMS = fpa_RMSValues[i];
			f_RMSSum+=fpa_RMSValues[i];
		}
		f_Ratio = f_MaxRMS / 5000;

		for (i=0 ; i<l_VADNelements ; i++)
            lpa_Basket[(short) (fpa_RMSValues[i]/f_Ratio) ]++;
		i=0;
		f_Stepm=f_Ratio/2;

		while (l_NoiseLength>0)
		{
			if (l_NoiseLength>lpa_Basket[i])
			{
				f_RMSSum1+=lpa_Basket[i] * f_Stepm;
				l_NoiseLength-=lpa_Basket[i];
			}
			else
			{
				f_RMSSum1+=l_NoiseLength * f_Stepm;
				l_NoiseLength=0;
			}
			f_Stepm+=f_Ratio;
			i++;
		}
		*pfNoiseLevel = (FLOAT) (20*log10( 1 + f_RMSSum1 / (l_VADNelements-l_SpeechLength) )) - (FLOAT) REF_DBOV;
		*pfSpeechLevel = (FLOAT) (20*log10( 1 + (f_RMSSum-f_RMSSum1) / l_SpeechLength )) - (FLOAT) REF_DBOV;
		*pfSnr = *pfSpeechLevel - *pfNoiseLevel;
		free(lpa_Basket);
		free(fpa_RMSValues);
	}
	else
	{
		lpa_Start = (UINT32 *) malloc (l_MaxNbSections * sizeof(UINT32));
		lpa_Stop  = (UINT32 *) malloc (l_MaxNbSections * sizeof(UINT32));

		lpa_Stop[0] = 0;
		l_StopNelements++;

		if (pfVad[0])
		{
			lpa_Start[l_StartNelements] = 1;
			l_StartNelements++;
		}
		for (i=1; i<l_VADNelements ; i++)
		{
			if ( (pfVad[i-1] == 0) && (pfVad[i] != 0) )
			{
				lpa_Start[l_StartNelements] = i * VAD_FRAME_LEN;
				l_StartNelements++;
			}
			else if ( (pfVad[i-1] != 0) && (pfVad[i] == 0) )
			{
				lpa_Stop[l_StopNelements] = i * VAD_FRAME_LEN;
				l_StopNelements++;
			}
		}

		lpa_Start[l_StartNelements] = l_VADNelements * VAD_FRAME_LEN;
		l_StartNelements++;

		for (i=0; i<l_StartNelements;i++)
		{
			if ((lpa_Start[i]-lpa_Stop[i])>16000)
			{
				if ( (i!=0) && ((lSignalLen-lpa_Stop[i])>4000) )
						lpa_Stop[i]+= 4000;
				if ((i!=l_StartNelements-1) && (lpa_Start[i]>4000) )
						lpa_Start[i]-=4000;
			}
			else
			{
				if ( (i!=0) && ( (lSignalLen-lpa_Stop[i])>((lpa_Start[i]-lpa_Stop[i])/8) ))
						lpa_Stop[i]+= (lpa_Start[i]-lpa_Stop[i])/8;
				if ( (i!=l_StartNelements-1) && ( lpa_Start[i]>((lpa_Start[i]-lpa_Stop[i])/8) ) )
						lpa_Start[i]-= (lpa_Start[i]-lpa_Stop[i])/8;
			}
		}

		for (i=0;i<l_StartNelements;i++)
		{
			f_RMSNoise += rms_nosqrt( &pfSignal[lpa_Stop[i]], lpa_Start[i] - lpa_Stop[i],FLOAT_TYPE);
			l_NoiseLength+=lpa_Start[i] - lpa_Stop[i];
		}
		for (i=0;i<l_StopNelements-1;i++)
		{
			f_RMSSpeech+= rms_nosqrt( &pfSignal[lpa_Start[i]], lpa_Stop[i+1] - lpa_Start[i],FLOAT_TYPE);
			l_SpeechLength+=lpa_Stop[i+1] - lpa_Start[i];
		}

		f_RMSNoise  = (FLOAT)sqrt(f_RMSNoise/l_NoiseLength);
		f_RMSSpeech = (FLOAT)sqrt(f_RMSSpeech/l_SpeechLength);

		if (f_RMSNoise)
		{
			*pfNoiseLevel = (FLOAT) (20*log10(1 + f_RMSNoise)) - REF_DBOV ;
		}
		else *pfNoiseLevel = - REF_DBOV;

		if (f_RMSSpeech)
		{
			*pfSpeechLevel = (FLOAT) (20*log10(1 + f_RMSSpeech)) - REF_DBOV;
		}
		else *pfSpeechLevel = - REF_DBOV;

		free(lpa_Start);
		free(lpa_Stop);
		*pfSnr = *pfSpeechLevel - *pfNoiseLevel;
	}

	return 0;
}

/*********************
*
*	FUNCTION: DC_remove
*
*	DESCRIPTION:
*		Time variing DC offset removal by a Notch filter.
*
***********************/
void DC_remove_float(const FLOAT *spa_Speech, FLOAT *spa_DCRemoved, UINT32 l_SpeechNelements)
{
	UINT32 i;
	const FLOAT f_coef = (FLOAT) 32735 / 32768;

	spa_DCRemoved[0] = spa_Speech[0];
    for (i = 1; i < l_SpeechNelements; i++)
    {
        spa_DCRemoved[i] = (spa_Speech[i] - spa_Speech[i-1] +
			((FLOAT) f_coef * spa_DCRemoved[i-1] ));
    }
}

/*********************
*
*	FUNCTION: DC_remove1
*
*	DESCRIPTION:
*		Time variing DC offset removal by a Notch filter.
*
***********************/
void DC_remove_int16(const INT16 *spa_Speech, FLOAT *spa_DCRemoved, UINT32 l_SpeechNelements)
{
	UINT32 i;
	const FLOAT f_coef = (FLOAT) 32735 / 32768;

	spa_DCRemoved[0] = spa_Speech[0];
    for (i = 1; i < l_SpeechNelements; i++)
    {
        spa_DCRemoved[i] = (spa_Speech[i] - spa_Speech[i-1] +
			((FLOAT) f_coef * spa_DCRemoved[i-1] ));
    }
}

/*********************
*
*	FUNCTION: normalise
*
*	DESCRIPTION:
*		This function normalises the input signal.
*		It first adjusts the signal level to -26dBov,
*		then a 4th oder Butterworth high-pass filter
*		at 100Hz cut-off frequency is applied.
*
***********************/
void normalise(INT16 *spa_NormSpeech, UINT32 l_NormSpeechNelements, const FLOAT f_SpeechLevel)
{
    UINT32 i;
    FLOAT f_coef;
    FLOAT *fpa_Speech;

    int Butter_4_100_8k_Nsos = 2;
    FLOAT Butter_4_100_8k_Hsos[10] = {
        0.9309753f, -1.8619506f, 0.9309753f, -1.8590763f, 0.8648249f, 0.96935381f, -1.9387076f, 0.96935381f,
        -1.9357148f, 0.94170045f};

    fpa_Speech = (FLOAT *) malloc(l_NormSpeechNelements * sizeof(FLOAT));

	f_coef = 1640 / (FLOAT) sqrt(f_SpeechLevel);

	for (i =0; i<l_NormSpeechNelements; i++)
        fpa_Speech[i] = (FLOAT) (spa_NormSpeech[i] * f_coef);

	IIRFilt( Butter_4_100_8k_Hsos, Butter_4_100_8k_Nsos, NULL, (FLOAT *) fpa_Speech, l_NormSpeechNelements, NULL );
	for (i=0; i<l_NormSpeechNelements; i++)
        spa_NormSpeech[i] = (INT16) fpa_Speech[i] ;
	free(fpa_Speech);
}

/*********************
*
*	FUNCTION: vad_process
*
*	DESCRIPTION:
*		This function generates the VAD profile based on an adaptative
*		power threshold, using an iterative approach.
*
***********************/
void vad_process( INT16 * X, INT32 Nsamples, FLOAT * VAD,
    FLOAT *LevelThresh, FLOAT *LevelSpeech, FLOAT *LevelNoise, FLOAT *Activity, INT32 Downsample )
{
    FLOAT g;
    FLOAT O_LevelThresh;
    FLOAT O_LevelNoise;
    FLOAT StDNoise;
    FLOAT O_LevelSig;
    FLOAT O_LevelMin;
    INT32  count;
    INT32  iteration;
    INT32  length;
    INT32  length_noise;
    INT32  start;
    INT32  finish;
    INT32  Nwindows;
    INT32 k;
    int VAD_alloc = 0;
    INT32 is_speech;

    Nwindows = Nsamples / Downsample;
    if( Nwindows < 10 ) return;

    for( count = 0; count < Nwindows; count++ )
    {
        VAD[count] = 0.0f;

        for( iteration = 0; iteration < Downsample; iteration++ )
        {
            g = X[count * Downsample + iteration];
            VAD[count] += (g * g);
        }
        VAD[count] /= Downsample;
    }

    O_LevelThresh = 0.0f;
    for( count = 0; count < Nwindows; count++ )
        O_LevelThresh += VAD[count];
    O_LevelThresh /= Nwindows;

    O_LevelMin = 0.0f;
    for( count = 0; count < Nwindows; count++ )
        if( VAD[count] > O_LevelMin )
            O_LevelMin = VAD[count];
    if( O_LevelMin > 0.0f )
        O_LevelMin *= 1.0e-4f;
    else
        O_LevelMin = 1.0f;
    for( count = 0; count < Nwindows; count++ )
        if( VAD[count] < O_LevelMin )
            VAD[count] = O_LevelMin;

    for( iteration = 0; iteration < 12; iteration++ )
    {
        O_LevelNoise = 0.0f;
        StDNoise = 0.0f;
        length = 0L;
        for( count = 0; count < Nwindows; count++ )
            if( VAD[count] <= O_LevelThresh )
            {
                O_LevelNoise += VAD[count];
                length++;
            }
        if( length > 0L )
        {
            O_LevelNoise /= length;
            for( count = 0; count < Nwindows; count++ )
                if( VAD[count] <= O_LevelThresh )
                {
                    g = VAD[count] - O_LevelNoise;
                    StDNoise += g * g;
                }
            StDNoise = (FLOAT)sqrt(StDNoise / length);
        }

        O_LevelThresh = 1.001f * (O_LevelNoise + 2.0f * StDNoise);
    }

    O_LevelNoise = 0.0f;
    O_LevelSig = 0.0f;
    length = 0L;
    for( count = 0; count < Nwindows; count++ )
    {
        if( VAD[count] > O_LevelThresh )
        {
            O_LevelSig += VAD[count];
            length++;
        }
        else
            O_LevelNoise += VAD[count];
    }
    if( length > 0L )
        O_LevelSig /= length;
    else
        O_LevelThresh = -1.0f;
    if( length < Nwindows )
        O_LevelNoise /= (Nwindows - length);
    else
        O_LevelNoise = 1.0f;

    for( count = 0; count < Nwindows; count++ )
        if( VAD[count] <= O_LevelThresh )
            VAD[count] = -VAD[count];

    VAD[0] = -O_LevelMin;
    VAD[Nwindows-1] = -O_LevelMin;

    start = 0L;
    finish = 0L;
    for( count = 1; count < Nwindows; count++ )
    {
        if( (VAD[count] > 0.0f) && (VAD[count-1] <= 0.0f) )
            start = count;
        if( (VAD[count] <= 0.0f) && (VAD[count-1] > 0.0f) )
        {
            finish = count;
            if( (finish - start) <= VAD_MIN_SPEECH_SECTION_LEN	 )
                for( iteration = start; iteration < finish; iteration++ )
                    VAD[iteration] = -VAD[iteration];
        }
    }

    if( O_LevelSig >= (O_LevelNoise * 1000.0f) )
    {
        for( count = 1; count < Nwindows; count++ )
        {
            if( (VAD[count] > 0.0f) && (VAD[count-1] <= 0.0f) )
                start = count;
            if( (VAD[count] <= 0.0f) && (VAD[count-1] > 0.0f) )
            {
                finish = count;
                g = 0.0f;
                for( iteration = start; iteration < finish; iteration++ )
                    g += VAD[iteration];
                if( g < 3.0f * O_LevelThresh * (finish - start) )
                    for( iteration = start; iteration < finish; iteration++ )
                        VAD[iteration] = -VAD[iteration];
            }
        }
    }

    start = 0L;
    finish = 0L;
    for( count = 1; count < Nwindows; count++ )
    {
        if( (VAD[count] > 0.0f) && (VAD[count-1] <= 0.0f) )
        {
            start = count;
            if( (finish > 0L) && ((start - finish) <= VAD_JOIN_SPEECH_SECTION_LEN	) )
                for( iteration = finish; iteration < start; iteration++ )
                    VAD[iteration] = O_LevelMin;
        }
        if( (VAD[count] <= 0.0f) && (VAD[count-1] > 0.0f) )
            finish = count;
    }

    if( start == 0L )
    {
        for( count = 0; count < Nwindows; count++ )
            VAD[count] = (FLOAT) fabs(VAD[count]);
        VAD[0] = -O_LevelMin;
        VAD[Nwindows-1] = -O_LevelMin;
    }

    count = 2;
    while( count < (Nwindows-2) )
    {
        if( (VAD[count] > 0.0f) && (VAD[count-2] <= 0.0f) )
        {
            VAD[count-2] = (FLOAT) (VAD[count] * 0.1);
            VAD[count-1] = (FLOAT) (VAD[count] * 0.3);
            count++;
        }
        if( (VAD[count] <= 0.0f) && (VAD[count-1] > 0.0f) )
        {
            VAD[count] = (FLOAT) (VAD[count-1] * 0.3);
            VAD[count+1] = (FLOAT) (VAD[count-1] * 0.1);
            count += 3;
        }
        count++;
    }

    for( count = 0; count < Nwindows; count++ )
        if( VAD[count] < 0.0f ) VAD[count] = 0.0f;

    length = 0;
    O_LevelSig = 0.0f;
    for( count = 0; count < Nwindows; count++ ) {
        if( VAD[count] > 0.0f ) {
            for( iteration = 0; iteration < Downsample; iteration++ )
            {
                g = X[count * Downsample + iteration];
                O_LevelSig += (g * g);
            }
            length++;
        }
    }

    length_noise = 0;
    O_LevelNoise = 0.0f;
    start = 0;
    is_speech = 0;
    for( count = 0; count < Nwindows; count++ ) {
        if( (!is_speech) && ((VAD[count] > 0.0f) || (count == (Nwindows-1))) )
        {
            finish = count;

            if( (finish - start) > (3 * VAD_MIN_SPEECH_SECTION_LEN	) ) {
                start += VAD_MIN_SPEECH_SECTION_LEN	;
                finish -= VAD_MIN_SPEECH_SECTION_LEN	;
                k = (finish - start) / 10;
                if( k > (VAD_JOIN_SPEECH_SECTION_LEN	/2) ) k = (VAD_JOIN_SPEECH_SECTION_LEN	/2);
                start += k * 2;
                finish -= k;

                length_noise += (finish - start);
                for( iteration = (start * Downsample);
                    iteration < (finish * Downsample); iteration++ )
                {
                    g = X[iteration];
                    O_LevelNoise += (g * g);
                }
            }
        }
        if( is_speech && (VAD[count] <= 0.0f) )
            start = count;
        is_speech = (VAD[count] > 0.0f);
    }
    if( length_noise == 0 ) {
        length_noise = 0;
        O_LevelNoise = 0.0f;
        for( count = 0; count < Nwindows; count++ ) {
            if( VAD[count] <= 0.0f ) {
                for( iteration = 0; iteration < Downsample; iteration++ )
                {
                    g = X[count * Downsample + iteration];
                    O_LevelNoise += (g * g);
                }
                length_noise++;
            }
        }
    }

    if( length > 0 ) O_LevelSig /= (length * Downsample);
    if( length_noise > 0 ) O_LevelNoise /= (length_noise * Downsample);

    if( Activity != NULL ) *Activity = (FLOAT)length / (FLOAT)Nwindows;

    if( LevelThresh != NULL ) *LevelThresh = O_LevelThresh;
    if( LevelSpeech != NULL ) *LevelSpeech = O_LevelSig;
    if( LevelNoise != NULL ) *LevelNoise = O_LevelNoise;

    if( VAD_alloc ) free(VAD);
}

/*********************
*
*	FUNCTION: SpeechSectionLevelVar
*
*	DESCRIPTION:
*		Computes the averaged energy for each speech section,
*		then calculate the difference between the max and min section.
*
***********************/
FLOAT SpeechSectionLevelVar(const FLOAT *VAD, const UINT32 l_VADNelements)
{
	UINT32 cpt,cpt1;
	FLOAT max=0,min=0,count;

	for (cpt=0; cpt<l_VADNelements; cpt++)
	{
		if (VAD[cpt]!=0)
		{
			for (cpt1=cpt,count=0; cpt1<l_VADNelements; cpt1++)
			{
				if (VAD[cpt1]!=0) count+=VAD[cpt1];
				else break;
			}
			if ((max==0)&&(min==0)) max=min=(count/(cpt1-cpt));
			else
			{
				count/=(cpt1-cpt);
				if (count<min) min=count;
				if (count>max) max=count;
			}
			cpt = cpt1;
		}
	}
	return (FLOAT) (sqrt(max)-sqrt(min))/1000;
}

/*********************
*
*	FUNCTION: hi_freq_variations
*
*	DESCRIPTION:
*		This function calculates the averaged energy
*		in the frequency band of 2.5-3.5kHz in speech sections.
*
***********************/
FLOAT hi_freq_variations(const INT16 *psNormSpeech, const INT32 lNormSpeechSz,
							const FLOAT *pfVad, const INT32 lVadSz)
{
	INT32 cnt,cnt1,cnt2;
	FLOAT *temp;
	FLOAT *hamming;
	FLOAT w;
	FLOAT value;
	INT32 frame_size = 1024;
	INT16 *psSpeechExtract=NULL;
	INT32 lSpeechExtractSz=0;
	INT16 min=0,max=0;
	FLOAT coef;

	psSpeechExtract = (INT16 *) calloc(lNormSpeechSz,sizeof(INT16));

	for (cnt=0; cnt<lVadSz; cnt++)
	{
		if (pfVad[cnt]>0)
		{
			for (cnt1=0; cnt1<VAD_FRAME_LEN; cnt1++)
				psSpeechExtract[lSpeechExtractSz + cnt1] = psNormSpeech[cnt*VAD_FRAME_LEN + cnt1];
			lSpeechExtractSz+=VAD_FRAME_LEN;
		}
	}

	for (cnt=0; cnt<lNormSpeechSz; cnt++)
	{
		if (psSpeechExtract[cnt]>max) max=psSpeechExtract[cnt];
		if (psSpeechExtract[cnt]<min) min=psSpeechExtract[cnt];
	}

	if (max<abs(min)) max = (INT16) abs(min);

	coef = (FLOAT) SHRT_MAX/max;

	for (cnt=0; cnt<lSpeechExtractSz; cnt++)
		psSpeechExtract[cnt] = (INT16) (psSpeechExtract[cnt]*coef);

	value = 0;
	temp = (FLOAT*) malloc( (frame_size +2) * sizeof(FLOAT) );

	hamming = (FLOAT*) malloc( (frame_size) * sizeof(FLOAT) );

	w = PI*2/frame_size;
	for (cnt=0;cnt<frame_size;cnt++) hamming[cnt] = (FLOAT) (0.54-0.46*cos(w*cnt));

	for (cnt=0; cnt<(lSpeechExtractSz-frame_size) ; cnt+=frame_size)
	{
		if ( standard_deviation(&psSpeechExtract[cnt],frame_size, SHORT_TYPE)>1 )
		{
			for (cnt1=0; cnt1<frame_size; cnt1++)
				temp[cnt1] = (FLOAT)(psSpeechExtract[cnt + cnt1]*hamming[cnt1]);

			RealFFT(temp,frame_size);

			for (cnt1=frame_size/8,cnt2=0; cnt1<7*frame_size/8; cnt1+=2,cnt2++)
				temp[cnt2] = (FLOAT) log10( sqrt(temp[cnt1]*temp[cnt1] + temp[cnt1+1]*temp[cnt1+1]) );

			value+=standard_deviation(temp,cnt2,FLOAT_TYPE);
		}
	}
	if (lSpeechExtractSz>0) value/=lSpeechExtractSz/frame_size;
	else value=0;
	FFTFree();

	free(psSpeechExtract);
	free(temp);
	free(hamming);

	return value;
}

/*********************
*
*	FUNCTION: mute_spotter
*
*	DESCRIPTION:
*		This function searches in the input signal for potential
*		mutes. It returns sum of the length of each signal drop.
*
***********************/
UINT32 mute_spotter(const INT16 *psRawSpeech, const INT32 psRawSpeechSz)
{
	FLOAT *Fspeech=NULL;
	FLOAT *Sspeech1=NULL,*Sspeech2=NULL;
	FLOAT *speech_RMS=NULL;
	INT32  RMS_Nelements=0;
	INT32 cnt=0,cnt1=0;
	INT32 *start=NULL,*stop=NULL;
	INT32 start_Nelements=0, stop_Nelements=0;
	char mute_detected=0;
	INT32 mute_amount = 0;
	FLOAT *start_slope=NULL, *stop_slope=NULL;
	INT32 start_RMS=0,stop_RMS=0;
	int IIR_high_1_500_8000_Nsos = 1;
	FLOAT IIR_high_1_500_8000_Hsos[5] = { 0.83408931895965f,  -0.83408931895965f, 0.0f, -0.66817863791930f, 0.0f };

	Sspeech1 = (FLOAT*) malloc ( psRawSpeechSz * sizeof(FLOAT) );
	Sspeech2= (FLOAT*) malloc(psRawSpeechSz * sizeof(FLOAT) );
	Fspeech = (FLOAT*) malloc ( psRawSpeechSz * sizeof(FLOAT) );

	for (cnt=0;cnt<psRawSpeechSz;cnt++)	Fspeech[cnt] = psRawSpeech[cnt];
	IIRFilt( IIR_high_1_500_8000_Hsos, IIR_high_1_500_8000_Nsos, NULL, Fspeech, psRawSpeechSz, NULL );

	for (cnt=0;cnt<psRawSpeechSz;cnt++) Sspeech1[cnt] = (FLOAT)round2(Fspeech[cnt]);

	free(Fspeech);

	DC_remove_float(Sspeech1,Sspeech2,psRawSpeechSz);
	DC_remove_float(Sspeech2,Sspeech1,psRawSpeechSz);
	DC_remove_float(Sspeech1,Sspeech2,psRawSpeechSz);
	DC_remove_float(Sspeech2,Sspeech1,psRawSpeechSz);
	DC_remove_float(Sspeech1,Sspeech2,psRawSpeechSz);

	RMS_Nelements = (INT32) (psRawSpeechSz / 32 );
	speech_RMS = (FLOAT *) calloc( RMS_Nelements + 4, sizeof(FLOAT) );
	for (cnt=0,cnt1=0;cnt<RMS_Nelements;cnt1+=32)
		speech_RMS[2+cnt++] = (FLOAT) (20*log10( rms( &(Sspeech2[cnt1]) ,32,FLOAT_TYPE) +1 ));

	start_RMS = 2;
	while( (speech_RMS[start_RMS]<50) && (start_RMS<RMS_Nelements) ) start_RMS++;
	stop_RMS = RMS_Nelements-1;
	while( (speech_RMS[stop_RMS]<50) && (stop_RMS>start_RMS) ) stop_RMS--;

	free(Sspeech1);
	free(Sspeech2);

	if (start_RMS<stop_RMS)
	{
		start = (INT32 *) malloc( (stop_RMS - start_RMS) * sizeof(INT32));
		stop =  (INT32 *) malloc( (stop_RMS - start_RMS) * sizeof(INT32));
	}

	for(cnt = (start_RMS+2); cnt<stop_RMS; cnt++)
		if((speech_RMS[cnt] < 30) && (speech_RMS[cnt-1] >30) && (speech_RMS[cnt-2] <30)) speech_RMS[cnt -1] = 0 ;
	for(cnt = (start_RMS+1); cnt<stop_RMS+1; cnt++)
	{
		if((speech_RMS[cnt] < 30) && (speech_RMS[cnt-1] > 30)) start[start_Nelements++] = cnt;
		if((speech_RMS[cnt] > 30) && (speech_RMS[cnt-1] < 30)) stop[stop_Nelements++] = cnt;
	}

	if (start_Nelements != stop_Nelements)
	{
		if (start_Nelements>stop_Nelements) start_Nelements = stop_Nelements;
		else stop_Nelements = start_Nelements;
	}

	start_slope = (FLOAT *) malloc ((start_Nelements+1) * sizeof(FLOAT));
	stop_slope = (FLOAT *) malloc  ((stop_Nelements+1) * sizeof(FLOAT) );

	for (cnt=0;cnt<start_Nelements;cnt++)
	{
		if (stop[cnt]-start[cnt] > 1)
		{
			start_slope[cnt] = speech_RMS[start[cnt]-2] - speech_RMS[start[cnt]];
			stop_slope[cnt] = speech_RMS[stop[cnt]+2] - speech_RMS[stop[cnt]-2];
			if ((start_slope[cnt]>30) && (stop_slope[cnt]>30)) mute_detected=1;
		}
		else stop[cnt]=start[cnt]=-1;
	}

	if (mute_detected)
	{
		for (cnt=0,cnt1=0;cnt<start_Nelements;cnt++)
			if ((start[cnt]!=-1) && (start_slope[cnt]>20) && (stop_slope[cnt]>20))
			{
				start[cnt1] = start[cnt];
				stop[cnt1] = stop[cnt];
				cnt1++;
			}
		start_Nelements = cnt1;

		for (cnt=0,cnt1=0;cnt<start_Nelements;cnt++)
		{
			start[cnt1]=start[cnt];
			while ( ((cnt+1)<start_Nelements) && ((start[cnt+1]-stop[cnt])<90) ) cnt++;
			stop[cnt1]=stop[cnt];
			cnt1++;
		}
		for (cnt=0;cnt<cnt1;cnt++) mute_amount+=stop[cnt]-start[cnt];
		if (mute_amount>250) mute_amount=250;
		if (mute_amount<0) mute_amount=0;
	}

	free(start_slope);
	free(stop_slope);
	free(start);
	free(stop);
	free(speech_RMS);
	return mute_amount;
}
