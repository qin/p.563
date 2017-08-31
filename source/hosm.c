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
Phone: +49 9131 691 160			Fax: +49 9131 691 325  
E-mail: info@opticom.de         www.3sqm.com  
  
Psytechnics Limited 
Fraser House, 23 Museum Street, Ipswich, IP1 1HN, UK 
Phone: +44 1 473 261 800		Fax: +44 1 473 261 880 
E-mail: info@psytechnics.com    www.psytechnics.com
  
SwissQual AG 
Gewerbestrasse 2 CH-4528 Zuchwil, Switzerland 
Phone: +41 32 685 08 30			Fax: +41 32 685 08 31   
E-mail: sales@swissqual.com     www.swissqual.com
  
Psytechnics, SwissQual or Opticom can provide licences and further information. 
  
Authors: 
      Ludovic Malfait ludovic.malfait@psytechnics.com 
      Roland Bitto rb@opticom.de 
      Pero Juric pero.juric@swissqual.com

********************************************************************/


/* Standard prototypes */ 
#include <time.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

/* specific prototypes */ 
#include "defines.h"
#include "hosm.h"
#include "dsp.h"
#include "vector_lib.h"
#include "generic_typedefs.h"


/***** local functions *******/
INT32 calc_seg_snr( FLOAT *pdA_Speech, INT32 lSpeechLen, typSegSNR_results *SegSNR_results);
INT32 find_noise_floor(UINT32 *peak, UINT32 *bound, UINT32 *min, UINT32 *Hist, FLOAT *HistValues, 
                             UINT32 hist_size, FLOAT max_above_peak);

INT32 energy_histogram(typInputParameter *par,typChannel *channel, INT32 file_id, UINT32 frame_size, 
					   UINT32 resolution, UINT32 *Hist, FLOAT *HistValues, UINT32 hist_size, FLOAT *hist_mean);

INT32 dCalcSpecPowDensFFT( FLOAT *pdA_InX, FLOAT *pdA_XSpecDens, INT32 iFrameLen);

INT32 calc_seg_snr( FLOAT *pdA_Speech,					/* pointer to speech samples */
                   INT32 lSpeechLen,					/* length of speech in samples */
                   typSegSNR_results *SegSNR_results	/* result structure */
                   );


/********************* 
* 
* FUNCTION: hosm 
* 
* DESCRIPTION: 
*	A main function for speech statistics analysis. 
*   The parameters used for speech statistics are mainly based on high order 
*	statistical evaluation of cepstral and LPC analyses. Kurtosis measures 
*   the degree of fat tails of a distribution whereas skewness measures the                   
* 	coefficient of asymmetry of a distribution.
*
*   Following features are calculated in this function:
*	- calculation of non-static noise parameters
*	- speech statistics
*	- signal mutes
*	- background noise
*
*    Following parameters are evaluated in this function:
*	- EstSegSNR;   						
*	- RelNoiseFloor;   					
*	- SpecLevelDev;   					
*	- SpecLevelRange;   				
*	- LPCskewAbs  							
*	- CepCurt
*	- CepSkew
*	- CepADev
*	- LPCskew
*	- EstBGNoise
*	- SpeechInterruptions	
* 
***********************/ 

INT32 hosm ( typInputParameter *par,typChannel *channel, FLOAT *HOS_MOVs)
{  
	FLOAT *dDegradedCopy=NULL;
	FLOAT *CepsOut=NULL;
	FLOAT dCurrentProcBuffer[BUFFSIZE];
	FLOAT dInputDataBuffers[4][BUFFSIZE/2];
	FLOAT dRMS_CurrentBuffer=0;
	FLOAT dCurrent_Mean=0;
	FLOAT lpcResidual[520];
	FLOAT LpcFilteredData[520];
	FLOAT msqd=0;
	FLOAT A_coeff=1.0;
	FLOAT lpcRMS=0;
	FLOAT lpcCoeffs[LPC_ORDER_HOSM+2];
	FLOAT lpcAver=0;
	FLOAT lpcADev=0;
	FLOAT lpcSDev=0;
	FLOAT lpcVar=0;
	FLOAT lpcSkew=0;
	FLOAT lpcCurt=0;
	FLOAT cepAver=0;
	FLOAT cepADev=0;
	FLOAT cepSDev=0;
	FLOAT cepVar=0;
	FLOAT cepSkew=0;
	FLOAT cepCurt=0;
	FLOAT StatNoiseDegradation;
	FLOAT PercOfInterruption;
	FLOAT *DataBuffX=NULL;
	FLOAT SumOfInterruptions=0;
	
	INT32 lpc_order = LPC_ORDER_HOSM;
	INT32 IntCnt=0;
	INT32 CodecDegCnt=0;
	INT32 iProcBufferCounter=0;
	INT32 MaxInterrLength=0;
	INT32 retVal=0;
			 
	UINT32 InterrCnt=0;

 	sctFilterInfo* psctFilter=NULL;
	typInterruption *Interruptions=NULL;

	typInterr_detect int_det;
    typSegSNR_results SegSNR_results;

	retVal = CalcScalingParams(par, channel);

	if(par->Frequency <= 0 || channel->FFTSize <= 0 || channel->sampleLength <= 0)
		return SQERR_WRONG_LENGTH;

	channel->MaxNrFrames = ((channel->file[REFER_FILE].stop - channel->file[REFER_FILE].start) / channel->FFTSize)+1;
    Interruptions = (typInterruption *)calloc(channel->MaxNrFrames/2, sizeof(typInterruption));

    dDegradedCopy = (FLOAT *)calloc(channel->FFTSize+2, sizeof(FLOAT));
    CepsOut = (FLOAT *)calloc(2*channel->FFTSize, sizeof(FLOAT));

	if(CepsOut==NULL || dDegradedCopy==NULL || Interruptions==NULL)	
	{
		fprintf(stderr,"can't allocate memory in \"hosm\"\r\n");
		return SQERR_ALLOC_ERROR;	
	}

	retVal = (INT32) FindNoiseFloors (par, channel);
	
	if(retVal < 0)	
		return retVal;

	retVal=(INT16)calc_seg_snr(par->ReferData, (INT32) channel->file[REFER_FILE].FileSize/2, &SegSNR_results);

	if(retVal < 0)  
		return retVal;
   
	CodecDegCnt	=0;
	InterrCnt	=0;
    channel->FrameCnt	=0;
	iProcBufferCounter	=0;
	SumOfInterruptions	=0;

	Interruptions[InterrCnt].SummOfAllInterrupts=0;
	MaxInterrLength = (INT32)((FLOAT)(MAX_LEN_FOR_INTERPOLATION)/channel->sampleLength);

	vclr(lpcCoeffs, 1, lpc_order+2);

	StatNoiseDegradation = MAX_MOS_VALUE;

	vclr(HOS_MOVs, 1, NR_OF_RESULTS);

	vclr(dInputDataBuffers[0], 1, BUFFSIZE/2);
	vclr(dInputDataBuffers[1], 1, BUFFSIZE/2);
	vclr(dInputDataBuffers[2], 1, BUFFSIZE/2);
	vclr(dInputDataBuffers[3], 1, BUFFSIZE/2);

	dInputDataBuffers[0][0]	= 0;
	dInputDataBuffers[1][0]	= 0;
	dInputDataBuffers[2][0]	= 0;
	dInputDataBuffers[3][0]	= 0;

	retVal=initDetectInterruptions(&int_det, par);
	
	if (retVal<0)
		return retVal;

	DataBuffX = par->ReferData;
	channel->MaxNrFrames = ((channel->file[REFER_FILE].stop - channel->file[REFER_FILE].start) / channel->FFTSize)-1;

	while (channel->FrameCnt < channel->MaxNrFrames)
	{
		iProcBufferCounter      = (channel->FrameCnt) % 4;
		channel->dataBuffer     = dInputDataBuffers[iProcBufferCounter];
		channel->dataBufferNext = dInputDataBuffers[(iProcBufferCounter+1) % 4];
		channel->dataBufferPrev = dInputDataBuffers[(iProcBufferCounter+3) % 4];

		vmov(DataBuffX + channel->FrameCnt*channel->FFTSize/2, 1, dCurrentProcBuffer, 1, channel->FFTSize);

		vmov(dCurrentProcBuffer, 1, dInputDataBuffers[(iProcBufferCounter+1) % 4], 1, channel->FFTSize/2);                                   	
		RemoveDCOffset(channel, dCurrentProcBuffer, dCurrentProcBuffer, (UINT32)(REFER_FILE), (UINT32)(channel->FFTSize));                   	
		mve(dCurrentProcBuffer, 1, &dCurrent_Mean, channel->FFTSize);                                                                        	
		                                                                                                                                     
        vsmul (dCurrentProcBuffer, 1, channel->RefLevelCorrection, dCurrentProcBuffer, 1, channel->FFTSize);                                     
		rmvesq (dCurrentProcBuffer, 1, &dRMS_CurrentBuffer, channel->FFTSize);                                                               	
		                                                                                                                                     
		retVal = detectInterruptions(par, channel, Interruptions+InterrCnt, &int_det);

		IntCnt=0;

		if (retVal==INTERR_FOUND && Interruptions[InterrCnt+IntCnt].start > 0)	
			for (IntCnt=0; IntCnt<int_det.num_int; IntCnt++) 
			{
				SumOfInterruptions += (FLOAT)(Interruptions[InterrCnt+IntCnt].length);

				InterrCnt++;
				if(Interruptions[InterrCnt+IntCnt].length < MaxInterrLength)
					Interruptions[InterrCnt+IntCnt].SummOfAllInterrupts=InterrCnt;
			}

		rmvesq (dCurrentProcBuffer, 1, &dRMS_CurrentBuffer, channel->FFTSize);
			
		if(dRMS_CurrentBuffer > (FLOAT)(MIN_SPEECH_ACTIV_RMS)) 
		{
			vmov(dCurrentProcBuffer,1,dDegradedCopy,1,channel->FFTSize);
			dDegradedCopy[channel->FFTSize+1] = dDegradedCopy[channel->FFTSize-1];
			dDegradedCopy[channel->FFTSize]   = dDegradedCopy[channel->FFTSize-1];

			EvalHigherMoments (dDegradedCopy, &lpcAver, &lpcADev, &lpcSDev, &lpcVar, &lpcSkew, &lpcCurt, channel->FFTSize);

			vsadd (dDegradedCopy,1,-lpcAver,dDegradedCopy,1,channel->FFTSize);

			if(lpcSDev>0)
			{
				vsmul (dDegradedCopy,1, (FLOAT)(1.0/lpcSDev), dDegradedCopy,1,channel->FFTSize);
				lpcCoeffs[0]=0;                                                                                             
				LPC_Burg(dDegradedCopy, lpc_order, lpcCoeffs, &msqd, channel->FFTSize-1);      

				psctFilter = sqft_NewFilter(lpcCoeffs+1, &A_coeff, lpc_order, 1);                                           
					                                                                                                            
				sqft_Filter(psctFilter, dDegradedCopy, LpcFilteredData, channel->FFTSize-1);                                      
				sqft_DeleteFilter(&psctFilter);                                                                              
					                                                                                                            
				vsub(dDegradedCopy+lpc_order+1, 1, LpcFilteredData+lpc_order, 1, lpcResidual, 1, channel->FFTSize-2*lpc_order);   
				rmvesq(lpcResidual, 1, &lpcRMS, channel->FFTSize-2*lpc_order);                                              
				EvalHigherMoments (lpcResidual, &lpcAver, &lpcADev, &lpcSDev, &lpcVar, &lpcSkew, &lpcCurt, channel->FFTSize-2*lpc_order);
					                                                                                                            
				cepstrum(dDegradedCopy, CepsOut, channel->FFTSize);                                                               
				EvalHigherMoments(CepsOut+1, &cepAver, &cepADev, &cepSDev, &cepVar, &cepSkew, &cepCurt, channel->FFTSize/2-1);      

				HOS_MOVs[0] += cepADev;	  
				HOS_MOVs[1] += cepSkew;	            
				HOS_MOVs[2] += cepCurt;	            
				HOS_MOVs[3] += lpcCurt;	            
				HOS_MOVs[4] += lpcSkew;	            

				CodecDegCnt++;                        
			}
		}
		channel->FrameCnt++;
	}

	if(channel->fileSize)
		SumOfInterruptions /= (FLOAT)(channel->fileSize/2);
	PercOfInterruption = SumOfInterruptions*(FLOAT)(100.0);
	
	/* scaling of the noise level */
	channel->NoiseLevelLog += (FLOAT)(20.0*log10(channel->RefLevelCorrection));

	/* stationary noise weighting */
	StatNoiseDegradation = (FLOAT)(-0.0625*channel->NoiseLevelLog + 1.3);

	if(StatNoiseDegradation > MAX_MOS_VALUE) 
		StatNoiseDegradation = MAX_MOS_VALUE;
	else if(StatNoiseDegradation < 0.5) 
		StatNoiseDegradation = 0.5;

	if(CodecDegCnt)  
	{
		vsdiv(HOS_MOVs, 1, (FLOAT)(CodecDegCnt), HOS_MOVs, 1, 6);

		HOS_MOVs[5] += (FLOAT)(fabs(HOS_MOVs[4]));  							
		HOS_MOVs[6]  = PercOfInterruption;							
		HOS_MOVs[7]  = StatNoiseDegradation/(FLOAT)(MAX_MOS_VALUE) - (FLOAT)(1.0);	
		HOS_MOVs[8]  = SegSNR_results.EstSegSNR;   						
		HOS_MOVs[9]	 = SegSNR_results.RelNoiseFloor;   					
		HOS_MOVs[10] = SegSNR_results.SpecLevelDev;   					
		HOS_MOVs[11] = SegSNR_results.SpecLevelRange;   				
	}
	
	DeinitDetectInterruptions(&int_det);
	Free((void**)&dDegradedCopy);

	free(Interruptions); 
	free(CepsOut);
	  
	return retVal;
}


/********************* 
* 
* FUNCTION: RemoveDCOffset
* 
* DESCRIPTION: 
*	Remove a mean value of the signal. Mean is calculated over the whole signal length 
* 
***********************/ 

void RemoveDCOffset(typChannel *channel, FLOAT *in, FLOAT *out, UINT32 file, UINT32 length)
{
	FLOAT dcOff=0;
	
	dcOff=(FLOAT)(channel->file[file].DCOffset*MAX_16BIT_VALUE/100);
	vsadd(in, 1, -dcOff, out, 1, length);
}




/********************* 
* 
* FUNCTION: initStructures 
* 
* DESCRIPTION: 
*	A main initialisation function for hosm(). All relevant parameters are initialized here.
* 
***********************/ 

INT32 initStructures(typInputParameter *par, typChannel *channel)
{
	channel->FFTSize = (INT32)(BUFFSIZE/2);
	channel->sampleLength = (FLOAT)(1.0)/(par->Frequency*(FLOAT)(1000.0));
	channel->file[REFER_FILE].frame_len  = channel->FFTSize;
	channel->file[REFER_FILE].frame_shft = (INT32)((FLOAT)(channel->file[REFER_FILE].frame_len)*FRAME_OVERLAP);
	channel->file[REFER_FILE].sumsq = 0.0;
	channel->file[REFER_FILE].AverageRMS = 0.0;
	channel->file[REFER_FILE].DCOffset = 0.0;
	channel->file[REFER_FILE].FileSize = 0;
	channel->file[REFER_FILE].start = 0;
	channel->file[REFER_FILE].stop  = 0;
	channel->file[REFER_FILE].useWholeFile = 1;	
	channel->fileSize=0;
	channel->FrameCnt=0;
	channel->RefLevelCorrection=0.0;
    par->bitResolution = 16;

	return 0;
}



/********************* 
* 
* FUNCTION: FindNoiseFloors
* 
* DESCRIPTION: 
*	Evaluates a background noise floor expressed in dBov. 
*	A calculation is based on a histogram of the RMS amplitude values of the signal. 
*	A main parameter calculated here is a "EstBGnoise".
* 
***********************/ 

INT32 FindNoiseFloors ( typInputParameter *par, typChannel *channel )
{
	UINT32	 EnDetFrameSize=0;
	UINT32   pos, bound_pos=0;
	UINT32   min_pos=0;
  	UINT32  *Histogram=NULL;    
  	UINT32   BitResolution=0;
	FLOAT	 hist_mean=0;
	FLOAT	 NoiseLevelReference=0;
  	FLOAT  *HistValues=0;  
	INT32    posMean=0;
  	INT32    err_code=0;     

	EnDetFrameSize = (UINT32)(channel->FFTSize);
	BitResolution = par->bitResolution;

    HistValues = (FLOAT *) calloc(ENERGY_WINDOWS, sizeof(FLOAT));
    if(HistValues == NULL)
		fprintf(stderr,"can't allocate memory for \"HistValues\"\r\n");

    Histogram = (UINT32 *)calloc(ENERGY_WINDOWS, sizeof(UINT32));
    if(Histogram == NULL) 
		fprintf(stderr,"can't allocate memory for \"Histogram\"\r\n");

	/*find peak in histogram, which should be the noise floor */
	posMean = energy_histogram(par, channel, REFER_FILE, EnDetFrameSize/16, BitResolution, 
		Histogram, HistValues, ENERGY_WINDOWS, &hist_mean);

	if(posMean < 0)
	{
		Free((void**)&HistValues);
		Free((void**)&Histogram);
		return posMean;
	}

	/*find peak in a histogram, which is a mean value of a noise */
	err_code=find_noise_floor(&pos, &bound_pos, &min_pos, Histogram, HistValues, ENERGY_WINDOWS, MAX_BOUND_ABOVE_PEAK);
	if(err_code < 0)
	{
		Free((void**)&HistValues);
		Free((void**)&Histogram);
		return err_code;
	}

	if(HistValues[bound_pos] > hist_mean)
		NoiseLevelReference = (hist_mean + HistValues[pos]) / 2;
	else
		NoiseLevelReference = HistValues[bound_pos];
			
	/* if NoiseBound less then Noise Peak */
	if(HistValues[bound_pos] < HistValues[pos])
		NoiseLevelReference = HistValues[pos];

	channel->NoiseLevelLog = NoiseLevelReference;

	Free((void**)&HistValues);
	Free((void**)&Histogram);

	return SQ_NO_ERRORS;
}



/********************* 
* 
* FUNCTION: energy_histogram
* 
* DESCRIPTION: 
*	A Noise floor, also known as the "bound position", is calculated 
*	as a minimum value between a noise peak and a mean. 
*   This function is called from FindNoiseFloors() which results in a parameter "EstBGnoise"
* 
***********************/ 

INT32 energy_histogram(typInputParameter *par, 
					   typChannel *channel, 
					   INT32 file_id, 
					   UINT32 frame_size, 
					   UINT32 resolution,
					   UINT32 *Hist, 
					   FLOAT *HistValues, 
					   UINT32 hist_size, 
					   FLOAT *hist_mean)
{
    INT32	i=0;
	INT32   maxFrames=0;
	INT32   eofcode=0;
	INT32   tempStop=0;
	INT32   tempStart=0;
	UINT32	index=0;
	UINT32  entry=0;
	UINT32  pos=0;
	
	FLOAT  min_en=0;
	FLOAT  max_en=0;
	FLOAT  total_energy=0;
	FLOAT  frame_en=0;
	FLOAT  EnDiff=0;
	FLOAT  Offset=0;
	FLOAT  frame_en_at0dB=0;
	FLOAT *Energy=NULL;
	FLOAT *wrkspace=NULL;
	FLOAT *DataBuffX=NULL;

	frame_en_at0dB = (FLOAT)pow(2.0, 2 * (resolution - 1)) * frame_size;
	tempStop  = channel->file[file_id].stop;
	tempStart = channel->file[file_id].start;

	channel->file[file_id].stop  = min(tempStop+(tempStop-tempStart)/3, channel->file[file_id].FileSize);

	if(frame_size > 0)
		channel->MaxNrFrames = (channel->file[file_id].stop - channel->file[file_id].start) / (2 * frame_size) + frame_size;
	else
	{
		channel->file[file_id].stop =tempStop ;
		channel->file[file_id].start=tempStart;
		return SQERR_WRONG_LENGTH;
	}
    
    wrkspace = (FLOAT *) calloc(max(channel->FFTSize, channel->MaxNrFrames), sizeof(FLOAT));
    if(wrkspace == NULL)
	{
		channel->file[file_id].stop =tempStop ;
		channel->file[file_id].start=tempStart;
		fprintf(stderr,"can't allocate memory for \"wrkspace\"\r\n");
		return SQERR_ALLOC_ERROR;
	}

	
	/* quit for unreasonable values */
	if(hist_size < 2) {
		fprintf(stderr,"Not enough binary trays in EnergyLevels Histogram!\n");
		channel->file[file_id].stop =tempStop ;
		channel->file[file_id].start=tempStart;
		return SQERR_HISTOGRAM_NOT_EN_TRAYS;
	}

	/* allocate moemory for frame energies */
	Energy = (FLOAT *)calloc(channel->MaxNrFrames + 10, sizeof(FLOAT));
    if(Energy == NULL) fprintf(stderr,"can't allocate memory for \"Energy\"\r\n");

    /* clear arrays before use */
    vfill ( 0.0, HistValues, 1, hist_size);
	for(index = 0; index < hist_size; index++)
		Hist[index] = 0;
	
    channel->FrameCnt = 0;
	total_energy	  = 0.0;
	eofcode = 0;
	Offset  = 0;
	EnDiff  = 0;

	DataBuffX = par->ReferData+channel->file[REFER_FILE].start;
	maxFrames = ((channel->file[REFER_FILE].stop - channel->file[REFER_FILE].start) / channel->FFTSize)-1;

	i=0;
    while (i < maxFrames-(INT32)frame_size)
	{
	 	if((UINT32)channel->FrameCnt >= channel->file[REFER_FILE].stop/(2*frame_size))
	 		break;

		/*	Read one frame of the data */
		vmov(DataBuffX + i*channel->FFTSize/2, 1, wrkspace, 1, channel->FFTSize);
		i++;

	    RemoveDCOffset(channel, wrkspace, wrkspace, (INT32)(REFER_FILE), channel->FFTSize);

   		frame_en = 0.0;
		for(index = 0; index < (UINT32)channel->FFTSize; index+=frame_size)
		{
			svesq(wrkspace+index, 1, &frame_en, frame_size);
			frame_en = max(frame_en, (FLOAT)frame_size);	/*prevent log10(0) errors */

			if(frame_en_at0dB > 0)
				Energy[channel->FrameCnt] = (FLOAT)10.0 * (FLOAT)log10(frame_en/frame_en_at0dB);
			channel->FrameCnt++;
		}
	}
    
	/*find min and max energies for setting the histogram boundries */
	maxv(Energy, 1, &max_en, &pos, channel->FrameCnt);
	minv(Energy, 1, &min_en, &pos, channel->FrameCnt);
 
	if(max_en == min_en)
	{
		fprintf(stderr,"Couldn\'t estimate the histogram boundries!\n");
		Free((void**)&Energy);
		channel->file[file_id].stop =tempStop ;
		channel->file[file_id].start=tempStart;
		return SQERR_HISTOGRAM_NOT_CREATED;
	}

	/*create histogram */
	for(index = 0; index < (UINT32)(channel->FrameCnt); index++) 
	{
		entry = max(0, min(hist_size - 1, (UINT32)((hist_size - 1) * (Energy[index] - min_en)/(max_en - min_en))));
		total_energy += Energy[index];		/* for creating a mean over the whole sample */
		Hist[entry] += 1;
	}

	/*frame energies in each bin tray */
	for(index = 0; index < hist_size; index++)
		HistValues[index] = min_en + (FLOAT)index * (max_en - min_en)/(hist_size - 1);

	Free((void**)&Energy);
	
	if(channel->FrameCnt<=0)
	{
		fprintf(stderr,"Number of frame wrong (%d)!\n", channel->FrameCnt);
		channel->file[file_id].stop =tempStop ;
		channel->file[file_id].start=tempStart;

		return SQERR_HISTOGRAM_NOT_CREATED;
	}

	/* return histogram mean */
	*hist_mean = total_energy / channel->FrameCnt;

	channel->file[file_id].stop =tempStop ;
	channel->file[file_id].start=tempStart;

	Free((void**)&wrkspace);

	return max(0, min(hist_size - 1, (UINT32)((hist_size - 1) * (*hist_mean - min_en)/(max_en - min_en))));
}




/********************* 
* 
* FUNCTION: find_noise_floor
* 
* DESCRIPTION: 
*	This function is called from FindNoiseFloors() which results in a parameter "EstBGnoise".
*   The resulting histogram of rms values calculated in energy_histogram() is analysed here.
* 
***********************/ 

INT32 find_noise_floor(UINT32 *peak, 
					   UINT32 *bound, 
					   UINT32 *min, 
					   UINT32 *Hist, 
					   FLOAT *HistValues,
					   UINT32 hist_size, 
					   FLOAT max_above_peak)
{
	UINT32 max_pos=0;
	UINT32 bound_pos=0;
	UINT32 min_pos=0;
	UINT32 mean_pos=0;
	UINT32 index=0;
	UINT32 frames=0;
		
	FLOAT mean_energy=0;
	FLOAT histMax=0;
	FLOAT histMin=0;
	INT32 Quotient=0;

	/* determine mean energy */
	mean_energy = 0.0;
	frames = 0;

	for(index = 0; index < hist_size; index++) 
	{
		mean_energy += Hist[index] * HistValues[index];
		frames += Hist[index];
	}

	if(!frames)
	{
		fprintf(stderr,"Couldn\'t estimate noise floor!\n");
		return SQERR_HISTOGRAM_NOT_CREATED;
	}

	Quotient = div(hist_size, 3).quot;
	mean_energy /=(FLOAT)frames;

	if(HistValues[hist_size - 1] == HistValues[0])
	{
		fprintf(stderr,"Couldn\'t estimate noise floor!\n");
		return SQERR_HISTOGRAM_NOT_CREATED;
	}

	/* use mean_pos as an upper bound for noise floor seraching */
	mean_pos = max(0, min(hist_size - 1, (UINT32)((hist_size - 1) * (mean_energy - HistValues[0])/
		(HistValues[hist_size - 1] - HistValues[0]))));

	/* find peak in histogram, which is the noise floor */
	histMax = 0;
	index = Quotient + 1;

	while(*((INT32 *)&histMax) < (INT32)(0.01*frames)) /* peak must be above 1% value */
	{
		maxi((INT32 *)Hist, 1, (INT32 *)&histMax, &max_pos, index);
		index++;
		if(index>= mean_pos) break;
	}
	
	if(mean_pos < max_pos)
		mean_pos=max(0, max_pos-1);
		
	/*find minimum beteween peak and mean */
	mini((INT32 *)(Hist + max_pos), 1, (INT32 *)&histMin, &min_pos, mean_pos - max_pos + 1);
	min_pos += max_pos;	/*started from max_pos */

	/*find the "right" point above noise */
	for(bound_pos = max_pos + 1; bound_pos < hist_size; bound_pos++) 
	{
		if(Hist[max_pos]>=0)
		{
			/*find MAX_BOUND_TO_PEAK bound_to_peak drop off point */
			if((FLOAT)Hist[bound_pos]/(FLOAT)Hist[max_pos] < MAX_BOUND_TO_PEAK) 
				break;
		}
		/*if the above fails, stop search when we are 'max_above_peak' in dB above noise peak */
		if(HistValues[bound_pos] > HistValues[max_pos]  + max_above_peak)
			break;
	}					   	
	
	bound_pos = min(bound_pos, hist_size - 1);
	
	if(bound_pos == hist_size)
		fprintf(stderr,"Warning, required speech to noise ratio not detected!");

	*bound = bound_pos;
	*peak = max_pos;
	*min = min_pos;
	return SQ_NO_ERRORS;
}




/********************* 
* 
* FUNCTION: sqft_NewFilter
* 
* DESCRIPTION: 
* 
*	Initialises One-dimensional digital filter.
    Y = FILTER(B,A,X) filters the data in vector X with the
    filter described by vectors A and B to create the filtered
    data Y.  The filter is a "Direct Form II Transposed"
    implementation of the standard difference equation:
 
    a(1)*y(n) = b(1)*x(n) + b(2)*x(n-1) + ... + b(nb+1)*x(n-nb)
                          - a(2)*y(n-1) - ... - a(na+1)*y(n-na)* 
***********************/ 

sctFilterInfo* sqft_NewFilter(FLOAT* pdNumCoefB,
							  FLOAT* pdDenomCoefA,
							  INT32 lNumCoefBLength,
							  INT32 lDenomCoefALength)
{
	sctFilterInfo* psctNewFilter = NULL;

	/* create new sturcture and fill it with the intitial parameters */
	if (psctNewFilter == NULL)
	{
		psctNewFilter = (sctFilterInfo*)calloc(1, sizeof(sctFilterInfo));
		psctNewFilter->pdDenomCoefA = NULL;
		psctNewFilter->pdNumCoefB = NULL;
		psctNewFilter->pdState = NULL;
	}

	if (psctNewFilter != NULL)
	{
		psctNewFilter->lNumCoefBLength = lNumCoefBLength;
		psctNewFilter->lDenomCoefALength = lDenomCoefALength;
		psctNewFilter->lStateLength = MAX(lNumCoefBLength, lDenomCoefALength) - 1;

		psctNewFilter->pdNumCoefB = (FLOAT*)calloc(lNumCoefBLength, sizeof(FLOAT));
		psctNewFilter->pdDenomCoefA = (FLOAT*)calloc(lDenomCoefALength, sizeof(FLOAT));
		psctNewFilter->pdState = (FLOAT*)calloc(psctNewFilter->lStateLength, sizeof(FLOAT));

		/* copy coefficients and scale them if a[0] != 0 */
		if (sqft_CoefCopyScale(psctNewFilter, pdNumCoefB, pdDenomCoefA) == -1)
		{
			sqft_DeleteFilter(&psctNewFilter);
			return NULL;
		}
		
		return psctNewFilter;
	}
	else
		return NULL;
}



/********************* 
* 
* FUNCTION: sqft_DeleteFilter
* 
* DESCRIPTION: 
*	Delete a filter structure created with sqft_NewFilter()
* 
***********************/ 

void sqft_DeleteFilter(sctFilterInfo** ppsctFilter)
{
	if (*ppsctFilter != NULL)
	{
		if ((*ppsctFilter)->pdDenomCoefA != NULL)
		{
			free((*ppsctFilter)->pdDenomCoefA);
			(*ppsctFilter)->pdDenomCoefA = NULL;
		}
		if ((*ppsctFilter)->pdNumCoefB != NULL)
		{
			free((*ppsctFilter)->pdNumCoefB);
			(*ppsctFilter)->pdNumCoefB = NULL;
		}
		if ((*ppsctFilter)->pdState != NULL)
		{
			free((*ppsctFilter)->pdState);
			(*ppsctFilter)->pdState = NULL;
		}

		free(*ppsctFilter);
		*ppsctFilter = NULL;
	}
}



/********************* 
* 
* FUNCTION: sqft_CoefCopyScale
* 
* DESCRIPTION: 
*	Copies values in coefficient vectors in filter structure,
*	if pdDenomCoefA[0] is not equal to 1, sqft_CoefCopyScale normalizes the filter
*	coefficients by pdDenomCoefA[0].
* 
***********************/ 

INT16 sqft_CoefCopyScale(sctFilterInfo* psctFilter,
						 FLOAT* pdNumCoefB,
						 FLOAT* pdDenomCoefA)
{
	INT32 li = 0;

	if (pdDenomCoefA != 0)
	{
		/* copy and normalize if a[0] != 1 */
		for (li = 0; li < psctFilter->lNumCoefBLength; li++)
			psctFilter->pdNumCoefB[li] = pdNumCoefB[li] / pdDenomCoefA[0];

		for (li = 0; li < psctFilter->lDenomCoefALength; li++)
			psctFilter->pdDenomCoefA[li] = pdDenomCoefA[li] / pdDenomCoefA[0];

		return 0;
	}
	else 
		return -1; /* error a[0] = 0 */
}




/********************* 
* 
* FUNCTION: sqft_Filter
* 
* DESCRIPTION: 
* 
*	One-dimensional digital filter of type pdOutData = sqft_Filter(B,A,pdInData).
*   filters the data in vector pdInData with the filter described by vectors A and B
*   to create the filtered data pdOutData.
*	The filter implements the following difference equation:
*
*                M                       N
*               ___                     ___
*	a[0] * y[n] = \   (b[k] * x[n - k]) - \   (a[k] * y[n - k])
*               /__                     /__
*               k=0                     k=1
*
*	where a's and b's are the filter feedback and feedforward coefficients and
*	x[n] and y[n] are the input and the output, respectively.  N+1 is the number of
*	feedback coefficients and M+1 is the number of feedforward
*	coefficients (N >= M).
* 
***********************/ 

INT16 sqft_Filter(sctFilterInfo* psctFilter,
				  FLOAT* pdInData,
				  FLOAT* pdOutData,
				  INT32 lDataLength)
{
	FLOAT* pdInPtrX  = NULL; /* input array */
	FLOAT* pdOutPtrY = NULL; /* output array */

	INT32 lSample = 0;
	INT32 li = 0;

	/* init input and output pointers */
	pdInPtrX = pdInData;
	pdOutPtrY = pdOutData;

	/* special case: 0th-order filter */
	if (psctFilter->lStateLength == 0)
	{
		for (lSample = 0; lSample < lDataLength; lSample++)
			pdOutPtrY[lSample] = psctFilter->pdNumCoefB[0] * pdInPtrX[lSample];
	}
	else
	{
		/* loop through blocks, processing samples */
		for (lSample = 0; lSample < lDataLength; lSample++)
		{
			/* output */
			pdOutPtrY[lSample] = psctFilter->pdNumCoefB[0] * pdInPtrX[lSample] + psctFilter->pdState[0];

			/* move filter timer */
			for (li = 1; li < psctFilter->lStateLength; li++)
				psctFilter->pdState[li - 1] = psctFilter->pdState[li];

			psctFilter->pdState[psctFilter->lStateLength - 1] = (FLOAT)0.0;

			/* calc b coefs */
			for (li = 1; li < psctFilter->lNumCoefBLength; li++)
				psctFilter->pdState[li - 1] += psctFilter->pdNumCoefB[li] * pdInPtrX[lSample];

			/* calc a coefs */
			for (li = 1; li < psctFilter->lDenomCoefALength; li++)
				psctFilter->pdState[li - 1] -= psctFilter->pdDenomCoefA[li] * pdOutPtrY[lSample];
		}
	}

	return 0;
}




/********************* 
* 
* FUNCTION: rmvesqth
* 
* DESCRIPTION: 
*	Normalized root mean squared of vector is calculated if the absolute value 
*	of an input data is above a threshold.
* 
***********************/ 

void rmvesqth(FLOAT *inp, 
			  INT16 astr, 
			  FLOAT *mean, 
			  FLOAT threshold, 
			  UINT32 len)
{
	FLOAT sum=0;
	UINT32 cnt=0;
	UINT32 i=0;                                         
	                                                    
  	sum = 0.0;                                          
  	*mean = *inp;                                       
  	if(len==0)                                          
  		  return;                                         
	                                                    
  	i=0;                                                
  	for(cnt=0; cnt<len; cnt++)                          
  	{                                                   
  		  if(fabs(*(inp+cnt*astr))> threshold)            
  		  {                                               
  			  sum += (*(inp+cnt*astr) * *(inp+cnt*astr)); 
  			  i++;                                        
  		  }                                               
  	}                                                   
	                                                    
  	if(i)                                               
  	{		                                              
  		  sum /= i;                                       
  		  sum = (FLOAT)sqrt(sum);                         
  	}                                                   
	                                                    
  	*mean = sum;                                        
}




/********************* 
* 
* FUNCTION: CalcScalingParams
* 
* DESCRIPTION: 
*	Opens input speech data and analyses the following properties:
*	- start position
*	- stop position 
*	- length of the data for analysis
*	- rms level
*	- DC Offset
*	- Nr of processing frames 
* 
***********************/ 

INT32 CalcScalingParams( typInputParameter *par, typChannel *channel)
{
  FLOAT	 EnDiff=0;
  FLOAT	 OffsetTmp=0;
  FLOAT	 DCOffset=0;
  FLOAT	 rmsBuff=0;
  FLOAT	 meanBuff=0;
  FLOAT	 *DataBuff=NULL;
  INT32	 BitResolution=0;
  INT32	 actflag=0;
  INT32	 startflag=0;
  INT32 cnt1=0;
  INT32 sampfact=0;
  INT32 cnt2=0;
  INT32 Delay=0;
  INT32 fileType=0;
  INT32 cnt3=0;
  INT32  id=REFER_FILE;

  INT32 frame_len  =(UINT32)(channel->FFTSize); 
  INT32 frame_shift=(UINT32)(channel->FFTSize)/2;

  if(id == REFER_FILE)
		DataBuff = par->ReferData;
  else
		DataBuff = par->CodedData;

  fileType = channel->file[id].filetype;

  if (fileType==0)
    sampfact=1;
  else
    sampfact=2;

  switch(channel->file[REFER_FILE].filetype)	
  {
  	case 0:	
		BitResolution = 8; 
		break;
  	case 1: case 2:
  		BitResolution = 12; 
		break;
  	default:
  		BitResolution = 16;
  }

  actflag = 0;
  startflag = 0;
  Delay = 0;
						 
  channel->file[id].sumsq   = 0.0;
  channel->file[id].start   = 0;
  channel->file[id].stop    = channel->file[id].FileSize;
  channel->file[id].filetype   = fileType;
  channel->file[id].frame_len  = frame_len;
  channel->file[id].frame_shft = frame_shift*(INT32)sampfact;

  /* Find start of speech activity in a reference file */
  cnt1 = cnt2 = cnt3 = 0;
  while(cnt1 < channel->file[id].FileSize/2 - channel->FFTSize/4)    
  {                                                                  
  	rmvesq(DataBuff+cnt1, 1, &rmsBuff, channel->FFTSize/4);          
  	mve(DataBuff+cnt1, 1, &meanBuff, channel->FFTSize/4);	/* calculate DC-offset in a signal     */
  	rmsBuff -= (FLOAT)fabs(meanBuff);						/* substract DC-Offset from rms value  */

  	if(rmsBuff > MIN_SPEECH_ACTIV_RMS)                               
  		cnt2++;                                                      
  	else                                                             
  		cnt2=0;                                                      
  	if(cnt2 >= 4)	
	{                                                
  		channel->file[id].start = cnt1 - cnt2*(channel->FFTSize/4);  
  		break;                                                       
  	}                                                                
  	cnt1+=(channel->FFTSize/4);                                      
  }                                                                  
  
  channel->file[id].start *= 2;                                   
  if(channel->file[id].start < 0)
	  channel->file[id].start = 0;

  /* Find end of speech activity in a reference file */
  cnt2 = 0;
  cnt1 = channel->file[id].FileSize/2 - channel->FFTSize/4;
  while(cnt1 > channel->FFTSize/4)    
  {                                                                  
  	rmvesq(DataBuff+cnt1, 1, &rmsBuff, channel->FFTSize/4);          
  	mve(DataBuff+cnt1, 1, &meanBuff, channel->FFTSize/4);	/* calculate DC-offset in a signal */
  	rmsBuff -= (FLOAT)fabs(meanBuff);						/* remove DC-Offset from rms value */

  	if(rmsBuff > MIN_SPEECH_ACTIV_RMS)                               
  		cnt2++;                                                      
  	else                                                             
  		cnt2=0;                                                      
  	if(cnt2 >= 4)	
	{                                                
  		channel->file[id].stop = cnt1 + cnt2*(channel->FFTSize/4);  
  		break;                                                       
  	}                                                                
  	cnt1-=(channel->FFTSize/4);                                      
  }                                                                  
  
  channel->file[id].stop *= 2;
  if(channel->file[id].stop < 0)
	  channel->file[id].stop = 0;

  if(channel->file[id].stop > channel->file[id].FileSize)
	  channel->file[id].stop = channel->file[id].FileSize;

  OffsetTmp = 0;
  DCOffset = 0;
  EnDiff = 0;
  cnt2 = 0;

  for (cnt1 = channel->file[id].start/2; ; cnt1+=channel->FFTSize)
  {
    if (cnt1 >= channel->file[id].stop/2-channel->FFTSize)
		break;

   	mve(DataBuff+cnt1, 1, &EnDiff, channel->FFTSize);
	OffsetTmp = EnDiff = (FLOAT)((0.9 * OffsetTmp + EnDiff)/1.9);
	DCOffset+=OffsetTmp;

	rmvesqth(DataBuff+cnt1, 1, &rmsBuff, (FLOAT)(MIN_SPEECH_ACTIV_RMS), channel->FFTSize);

	if(rmsBuff > MIN_SPEECH_ACTIV_RMS)	
	{							
		channel->file[id].sumsq += rmsBuff;
		cnt3++;
	}
	cnt2++;
  } 

  if(cnt2>0) 
	  channel->file[id].DCOffset=(FLOAT)(DCOffset/cnt2/MAX_16BIT_VALUE*100);

  if(cnt2)
  {
	if(channel->file[id].sumsq>0 && cnt3 > 0)
	{
		channel->file[id].AverageRMS = channel->file[id].sumsq / cnt3;
		channel->file[id].AverageRMS = (FLOAT)(20.0*log10(channel->file[id].AverageRMS/pow(2.0,BitResolution-1.0))+3);

		/* calculate Level difference of Nominal rms Level to Overload point (-26dBov) */
		if(id==REFER_FILE)
			channel->RefLevelCorrection = (FLOAT)pow(10.0, (MEAN_RMS_LEVEL - channel->file[id].AverageRMS)/20.0);
 	}
	else
		channel->file[id].AverageRMS = 0;
  }
  else
		channel->file[id].AverageRMS = 0;

	channel->fileSize = channel->file[REFER_FILE].stop - channel->file[REFER_FILE].start;

	if(channel->file[REFER_FILE].useWholeFile) 
	{
		channel->file[REFER_FILE].start	= 0;
		channel->file[REFER_FILE].stop 	= channel->file[REFER_FILE].FileSize;
	}

 	channel->MaxNrFrames = ((channel->file[REFER_FILE].stop - channel->file[REFER_FILE].start) / channel->FFTSize)+1;

  return SQ_NO_ERRORS;
}




/********************* 
* 
* FUNCTION: calc_seg_snr 
* 
* DESCRIPTION: 
*	This function estimates the segmental signal-to-noise ratio
*	by analysing the range of the levels in the spectral domain for active speech frames
*  
*	This function includes an own pre-processing part for calculating the
*	speech activity threshold as well as a filtering with a weak telephony bandpass
*  
*	Besides of the main result
*  - Est SegSNR 
*
*  additional intermeadiate results will be provided to the frame function
* 
*	- RelNoiseFloor;   					
*	- SpecLevelDev;   					
*	- SpecLevelRange;   				
* 
***********************/ 


INT32 calc_seg_snr( FLOAT *pdA_Speech, 
				  INT32 lSpeechLen, 
				  typSegSNR_results *SegSNR_results)
{
   
   INT32      i, j, k;
   INT32      iRetVal;
   FLOAT      dEstSegSNR;       
   FLOAT      dSNRCheck;

   FLOAT      dSpeech20PercThresh;
   FLOAT      dSpeech80PercThresh;
   FLOAT      dStatSNR;

   FLOAT      dActiveSpeechTresh;
   FLOAT      dSpeechLevel;
   
   FLOAT      dFilterOut, dMean, dSqareSum, dSpecDev;

   INT32      iBlockLen;
   FLOAT      dNoisePowBlock;
   FLOAT      dNoiseLevel;  
   FLOAT      dNoisePowSum;  

   INT32      iFrameLen;
   INT32      iStepPerFrame;
   FLOAT      dFrameOverLap;

   INT32     lNbrOfFrames;
   INT32     lSpeechFrameCtr;

   FLOAT      dSpec15PercThresh;
   FLOAT      dSpec80PercThresh;

   FLOAT      dSpecDevSum;
   FLOAT      dSpecRangeSum;
   FLOAT      dSpecPowFrameSum;

   FLOAT      *pdA_SpeechBPfilt;

   FLOAT      *pdA_SpeechEnv;
   FLOAT      *pdA_Distribution;

   FLOAT      *pdA_SpecPowFrame;
   FLOAT      *pdA_SpecPowTotal;
   FLOAT      *pdA_SpecPowMin;

   INT32        iAllocErrCtr;



   INT32         iOrderBandPass8kHz = 64;
   FLOAT      pdA_CoeffBandPass8kHz[65]  =  
   {(FLOAT) 0.0000531, (FLOAT) 0.0003059, (FLOAT) 0.00028,   (FLOAT) 0.0005082, (FLOAT) 0.000634,  (FLOAT) 0.0008517, (FLOAT) 0.0010684,
    (FLOAT) 0.0014012, (FLOAT) 0.0013436, (FLOAT) 0.002116,  (FLOAT) 0.0011138, (FLOAT) 0.0026455, (FLOAT) 0.0001507, (FLOAT) 0.0021984,
	(FLOAT)-0.001554,  (FLOAT)-0.0002023, (FLOAT)-0.0040758, (FLOAT)-0.00498,   (FLOAT)-0.0081556, (FLOAT)-0.0111685, (FLOAT)-0.0155649,
	(FLOAT)-0.0162199, (FLOAT)-0.0285923, (FLOAT)-0.0169833, (FLOAT)-0.0485369, (FLOAT)-0.0115695, (FLOAT)-0.0739747, (FLOAT)-0.0009768,
	(FLOAT)-0.100091,  (FLOAT) 0.0107859, (FLOAT)-0.1199559, (FLOAT) 0.0184921, (FLOAT) 0.8725996, (FLOAT) 0.0184921, (FLOAT)-0.1199559,
	(FLOAT) 0.0107859, (FLOAT)-0.100091,  (FLOAT)-0.0009768, (FLOAT)-0.0739747, (FLOAT)-0.0115695, (FLOAT)-0.0485369, (FLOAT)-0.0169833,
	(FLOAT)-0.0285923, (FLOAT)-0.0162199, (FLOAT)-0.0155649, (FLOAT)-0.0111685, (FLOAT)-0.0081556, (FLOAT)-0.00498,   (FLOAT)-0.0040758,
	(FLOAT)-0.0002023, (FLOAT)-0.001554,  (FLOAT) 0.0021984, (FLOAT) 0.0001507, (FLOAT) 0.0026455, (FLOAT) 0.0011138, (FLOAT) 0.002116,
	(FLOAT) 0.0013436, (FLOAT) 0.0014012, (FLOAT) 0.0010684, (FLOAT) 0.0008517, (FLOAT) 0.000634,  (FLOAT) 0.0005082, (FLOAT) 0.00028, 
    (FLOAT) 0.0003059, (FLOAT)0.0000531 };


	/* **************** Initializing and pre-processing ************************************ */
   iFrameLen         = 512;   
   dFrameOverLap     = 0.75; 

   iStepPerFrame     = iFrameLen - (INT32) (dFrameOverLap*iFrameLen);
   lNbrOfFrames      = (INT32) ((FLOAT) (lSpeechLen-iFrameLen) / (iFrameLen - (INT32) (dFrameOverLap*iFrameLen)));

   iAllocErrCtr   = 0; 
   
   if( (pdA_SpeechBPfilt   = (FLOAT *) calloc( lSpeechLen+1,  sizeof( FLOAT))) == NULL) iAllocErrCtr++;
   if( (pdA_SpecPowFrame   = (FLOAT *) calloc( iFrameLen+1,   sizeof( FLOAT))) == NULL) iAllocErrCtr++;
   if( (pdA_SpecPowTotal   = (FLOAT *) calloc( iFrameLen+1,   sizeof( FLOAT))) == NULL) iAllocErrCtr++;
   if( (pdA_SpecPowMin     = (FLOAT *) calloc( iFrameLen+1,   sizeof( FLOAT))) == NULL) iAllocErrCtr++;
   if( (pdA_SpeechEnv      = (FLOAT *) calloc( lNbrOfFrames+1,sizeof( FLOAT))) == NULL) iAllocErrCtr++;
   if( (pdA_Distribution   = (FLOAT *) calloc( 101,           sizeof( FLOAT))) == NULL) iAllocErrCtr++;


   if( iAllocErrCtr) 
   {
      if( pdA_SpeechBPfilt)   free( pdA_SpeechBPfilt);
      if( pdA_SpecPowFrame)   free( pdA_SpecPowFrame);
      if( pdA_SpecPowTotal)   free( pdA_SpecPowTotal);
      if( pdA_SpecPowMin)     free( pdA_SpecPowMin);
      if( pdA_SpeechEnv)      free( pdA_SpeechEnv);
      if( pdA_Distribution)   free( pdA_Distribution);
      return( SQERR_ALLOC_ERROR);
   }

   for( i=0; i<iOrderBandPass8kHz; i++)
   {
		dFilterOut  =  0.0;
		for( j=0; j<iOrderBandPass8kHz && (i-j)>=0; j++)
		{  
			dFilterOut  += pdA_Speech[i-j] * pdA_CoeffBandPass8kHz[j];
		}
		pdA_SpeechBPfilt[i] = (FLOAT) dFilterOut;
   }


   for( i=iOrderBandPass8kHz; i<lSpeechLen; i++)
   {
		dFilterOut  =  0.0;
		for( j=0; j<iOrderBandPass8kHz; j++)
		{  
			dFilterOut  += pdA_Speech[i-j] * pdA_CoeffBandPass8kHz[j];
		}
		pdA_SpeechBPfilt[i] = dFilterOut;
   }
  
   
   for( i=0; i<lNbrOfFrames; i++)
   {
      for( j=0, dMean=0, dSqareSum=0; j<iFrameLen; j++)	
      {
         dSqareSum	+=	(pdA_SpeechBPfilt[i*iStepPerFrame+j]) * (pdA_SpeechBPfilt[i*iStepPerFrame+j]); 
         dMean	    +=  pdA_SpeechBPfilt[i*iStepPerFrame+j];
      }
	  pdA_SpeechEnv[i]	= (dSqareSum - dMean*dMean / iFrameLen) / iFrameLen / 32768 / 32768;
      pdA_SpeechEnv[i]  =  10 * (FLOAT)(log10( pdA_SpeechEnv[i] + 1e-16));
   }   
   
   Calc_DistrOfVector( pdA_SpeechEnv, lNbrOfFrames, -100, 0, 100, pdA_Distribution);


   dSpeech20PercThresh  = Calc_PercentilOfDistrVector( 20.0, -100, 0, 100, pdA_Distribution);
   dSpeech80PercThresh  = Calc_PercentilOfDistrVector( 80.0, -100, 0, 100, pdA_Distribution);

   dStatSNR             = dSpeech20PercThresh-dSpeech80PercThresh;
   dActiveSpeechTresh   = dSpeech20PercThresh - (FLOAT)(4.0); 

   lSpeechFrameCtr   = 0;
   dSpeechLevel      = 0;
   for( i=0; i<lNbrOfFrames; i++)
   {
      if( pdA_SpeechEnv[i] > dActiveSpeechTresh) 
      {
         dSpeechLevel  +=  (FLOAT)(pow( 10.0, pdA_SpeechEnv[i] / 10.0));
         lSpeechFrameCtr++;
      }
   }

   dSpeechLevel   = 10 * (FLOAT)(log10( dSpeechLevel / lSpeechFrameCtr + 1e-16));

   
/* ==> Start spectral evaluation #1: Calculates only total average spectrum (PSD) -------- */
   for( i=0; i<lNbrOfFrames; i++)       
   {
      if( pdA_SpeechEnv[i]>dActiveSpeechTresh)   
      {  

         iRetVal  =  dCalcSpecPowDensFFT( &pdA_SpeechBPfilt[i*iStepPerFrame], pdA_SpecPowFrame, iFrameLen);
         if( iRetVal == SQERR_ALLOC_ERROR)   return( SQERR_ALLOC_ERROR);

         for( k=0; k<iFrameLen/2; k++) 
         {
            pdA_SpecPowTotal[k] += pdA_SpecPowFrame[k]  / (FLOAT)(32768.0) / (FLOAT)(32768.0);
         }
      }
   }

   for( k=0; k<iFrameLen/2; k++) 
   {
      pdA_SpecPowTotal[k] = pdA_SpecPowTotal[k] / lNbrOfFrames;
   }



/* *************** SpcLevelDev and SpcLevelRange ****************** */
   
   lSpeechFrameCtr      =  0;
   dSpecDevSum          =  0.0;
   dSpecRangeSum        =  0.0;

   for( i=0; i<lNbrOfFrames; i++)      
   {
      if( pdA_SpeechEnv[i]>dActiveSpeechTresh)   
      {  

         iRetVal  =  dCalcSpecPowDensFFT( &pdA_SpeechBPfilt[i*iStepPerFrame], pdA_SpecPowFrame, iFrameLen);
         if( iRetVal == SQERR_ALLOC_ERROR)   return( SQERR_ALLOC_ERROR);
         
         for( k=0; k<(iFrameLen/2); k++) 
         {
            pdA_SpecPowFrame[k]  = pdA_SpecPowFrame[k]  / (FLOAT)(32768.0) / (FLOAT)(32768.0);
            pdA_SpecPowFrame[k]  = pdA_SpecPowFrame[k]  / (pdA_SpecPowTotal[k] + (FLOAT)(1e-06));
            pdA_SpecPowFrame[k]  = 10 * (FLOAT)(log10( pdA_SpecPowFrame[k] + 1e-16));

         }

         for( j=(INT32) (iFrameLen/8)+1, dMean=0, dSqareSum=0; j<(INT32)(iFrameLen/8)+(iFrameLen/4); j++)	
         {
            dSqareSum	+=	(pdA_SpecPowFrame[j]) * (pdA_SpecPowFrame[j]); 
            dMean	      +=	 pdA_SpecPowFrame[j];
         }
	      dSpecDev	= (dSqareSum - dMean*dMean / (FLOAT) (iFrameLen/4) ) / (FLOAT)  (iFrameLen/4);
         if( dSpecDev>0)  dSpecDevSum    +=  (FLOAT)sqrt( dSpecDev);


         Calc_DistrOfVector( pdA_SpecPowFrame, (iFrameLen/2),         -80, 0, 40, pdA_Distribution);
         dSpec15PercThresh       = Calc_PercentilOfDistrVector( 15.0, -80, 0, 40, pdA_Distribution);
         dSpec80PercThresh       = Calc_PercentilOfDistrVector( 80.0, -80, 0, 40, pdA_Distribution);
         dSpecRangeSum    +=  dSpec15PercThresh - dSpec80PercThresh;

         lSpeechFrameCtr ++;
      }
   }

   dSpecDevSum    = dSpecDevSum   / lSpeechFrameCtr;      /*  SpcLevelDev */
   dSpecRangeSum  = dSpecRangeSum / lSpeechFrameCtr;      /*  SpcLevelRange */


   
/* ******************* RelNoiseFloor ********************************** */

   iBlockLen       =  10;     
   dNoisePowBlock  =  0.0;
   dNoisePowSum    =  0.0;

   for( i=0; i<(lNbrOfFrames-1)/iBlockLen; i++)       
   {
      for( k=0; k<iFrameLen/2; k++) 
      {
         pdA_SpecPowMin[k]    = 32768.0 * 32768.0;
      }

      for( j=0; j<iBlockLen; j++)
      {  

         iRetVal =   dCalcSpecPowDensFFT( &pdA_SpeechBPfilt[(i*iBlockLen+j)*iStepPerFrame], pdA_SpecPowFrame, iFrameLen);
         if( iRetVal == SQERR_ALLOC_ERROR)   return( SQERR_ALLOC_ERROR);
         
         for( k=0, dSpecPowFrameSum=0.0; k<iFrameLen/2; k++) 
         {
            pdA_SpecPowFrame[k]  = pdA_SpecPowFrame[k]  / (FLOAT)(32768.0) / (FLOAT)(32768.0);
            if( pdA_SpecPowFrame[k]<1e-8)  pdA_SpecPowFrame[k]  =  (FLOAT)(1e-8);
            if( pdA_SpecPowFrame[k]<pdA_SpecPowMin[k])   
            {  
               pdA_SpecPowMin[k]   =   pdA_SpecPowFrame[k];
            }

         }
  
      }

      dNoisePowBlock=0.0;
      for( k=iFrameLen/8; k<(iFrameLen/4); k++) 
      {
         dNoisePowBlock    += pdA_SpecPowMin[k];
      }

      for( j=0; j<iBlockLen; j++)
      {
         if( pdA_SpeechEnv[i*iBlockLen+j] > dActiveSpeechTresh)   dNoisePowSum   +=  dNoisePowBlock;  /* valid for all 10 Blocks */
      }
   }

   dNoiseLevel   = 10 * (FLOAT)(log10( dNoisePowSum / lSpeechFrameCtr  + 1e-16));

   dEstSegSNR    = (dSpecRangeSum - (FLOAT)(10.5)) * (FLOAT)(2.7);
   dSNRCheck    = (dSpeechLevel - dNoiseLevel - 22) * 3 - dEstSegSNR;

   if( dEstSegSNR>25 || dStatSNR<38 || dSpecDevSum>8.1 || fabs( dSNRCheck)>20)  dEstSegSNR    =  40;

   SegSNR_results->EstSegSNR      =  dEstSegSNR;
   SegSNR_results->SpecLevelDev   =  dSpecDevSum;
   SegSNR_results->SpecLevelRange =  dSpecRangeSum;
   SegSNR_results->RelNoiseFloor  =  dSpeechLevel - dNoiseLevel;

   free( pdA_SpeechBPfilt);
   free( pdA_SpecPowFrame);
   free( pdA_SpecPowTotal);
   free( pdA_SpecPowMin);
   free( pdA_SpeechEnv);
   free( pdA_Distribution);
 
 
   return( 1);

}




/********************* 
* 
* FUNCTION: Calc_DistrOfVector 
* 
* DESCRIPTION: 
*	This function calculates the distribution of a given vector.
*	It is used within calc_seg_snr to provide thresholds for speech activity 
*	or upper and lower bounds in the spectral analysis.
* 
***********************/ 

INT32 Calc_DistrOfVector( FLOAT *pdA_Vector, 
						 INT32 lVecLen, 
						 FLOAT dMinValue, 
						 FLOAT dMaxValue, 
						 INT32 iLenDistrVec, 
						 FLOAT *pdA_DistrVec)
{
   INT32   i, j, iCtrStart;
   FLOAT   dStepSize;

   dStepSize   =  (dMaxValue - dMinValue) / (FLOAT) iLenDistrVec;

   for( i=0; i<lVecLen; i++)
   {
      if( pdA_Vector[i]>dMaxValue)	pdA_Vector[i] = dMaxValue;
      if( pdA_Vector[i]<dMinValue)	pdA_Vector[i] = dMinValue;

      iCtrStart   =  (INT32) ( (pdA_Vector[i] - dMinValue) / dStepSize );

      for( j=iCtrStart; j<iLenDistrVec; j++)  
		  pdA_DistrVec[j] += 1.0;
   }

   for( j=0; j<iLenDistrVec; j++)  
	   pdA_DistrVec[j] = pdA_DistrVec[j] / (FLOAT) lVecLen;

   return( 0);
   
}




/********************* 
* 
* FUNCTION: Calc_PercentilOfDistrVector 
* 
* DESCRIPTION: 
*	This function is directly related to Calc_DistrOfVector(). It is calculating 
*	percentile values from a given distribtion vector.
* 
***********************/ 

FLOAT Calc_PercentilOfDistrVector( FLOAT dPercentil, 
								  FLOAT dMinValue, 
								  FLOAT dMaxValue, 
								  INT32 iLenDistrVec, 
								  FLOAT *pdA_DistrVec)
{
   INT32  i;
   INT32  iPercIndex;
   FLOAT  dDistrSum;
   FLOAT  dTargetPercentil;
   FLOAT  dPercentilValue;
   FLOAT  dStepSize;

   dTargetPercentil  =  (FLOAT)(1.0) - (dPercentil / 100);
   dStepSize         =  (dMaxValue - dMinValue) / (FLOAT) iLenDistrVec;
   dDistrSum         =  0.0;

   for( i=0; i<iLenDistrVec && pdA_DistrVec[i]<dTargetPercentil; i++);  

   iPercIndex  =  i;
   
   dPercentilValue   =  ((FLOAT) iPercIndex * dStepSize) + dMinValue;

   return( dPercentilValue);

}



/********************* 
* 
* FUNCTION: dCalcSpecPowDensFFT 
* 
* DESCRIPTION: 
*	This function calculates the (estimated) spectral power density (periodogram) of a given (FLOAT) vector
*	with the length iFrameLen (limited to INT32)
*	The input vector will not be affected, the results are stored in pdA_XSpecDens.
*  
*	Please note: 
*	- A Hann-Window is applied
*	- The length of the power density vector is only the half of iFrameLen
* 
***********************/ 

INT32	dCalcSpecPowDensFFT( FLOAT *pdA_InX, 
							FLOAT *pdA_XSpecDens, 
							INT32 iFrameLen)
{	
	INT32 	iLenHalf, i;
	FLOAT   dHannCorrection;
	FLOAT   dHannValue;

	FLOAT  *pfA_DummyZ;

	dHannCorrection	=  (FLOAT)(2.66666666);  					/* Correction for HANN-Window  => 1 / 0.375 */

	iLenHalf = (INT32) ((FLOAT) iFrameLen / 2);

	pfA_DummyZ  =  (FLOAT *) malloc( sizeof( FLOAT)*(2*iFrameLen+2));
	if( pfA_DummyZ==NULL) 
		return( SQERR_ALLOC_ERROR);
   
	for( i=0; i<iFrameLen; i++)
	{  
		pfA_DummyZ[2*i]   =  (FLOAT) pdA_InX[i];
		pfA_DummyZ[2*i+1] =  0.0;
	}

	for (i=1; i<=iLenHalf; i++)
	{	
		dHannValue	=	(FLOAT)(0.5 * ( 1.0 + cos( 3.1415926 * ((2 * i) - 1) / (iFrameLen-1))));

		pfA_DummyZ[2*(iLenHalf + i - 1)]		*=		(FLOAT)dHannValue;
		pfA_DummyZ[2*(iLenHalf -i)]			*=		(FLOAT)dHannValue;
	}

	FFT( pfA_DummyZ, iFrameLen);

	pdA_XSpecDens[0]  =  0;

	for(i=1; i<iLenHalf; i++)
	{	
		pdA_XSpecDens[i]  =  (FLOAT) 2.0 * (pfA_DummyZ[2*i]*pfA_DummyZ[2*i] + pfA_DummyZ[2*i+1]*pfA_DummyZ[2*i+1]);
		pdA_XSpecDens[i]  =  pdA_XSpecDens[i] * dHannCorrection / iFrameLen / iFrameLen;
	}

	free(pfA_DummyZ);

	return( 1);
}


