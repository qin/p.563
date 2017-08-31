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
#include <math.h>
#include <string.h>
#include <float.h>

/* specific prototypes */ 
#include "defines.h"
#include "generic_typedefs.h"
#include "interr_detect.h"
#include "vector_lib.h"


/* Local prototypes */ 

/*---------------------------------------------------------------------------------------------- */
/* Build fundamental frames and its averages */
/*---------------------------------------------------------------------------------------------- */
INT32 fund_frames(typInputParameter *par, typChannel *channel, typInterruption *interruption,
								typInterr_detect *int_det, INT32 word_begin);

/*---------------------------------------------------------------------------------------------- */
/* Test fundamental frame for interruption */
/*----------------------------------------------------------------------------------------------	 */
INT32 testframe(typInputParameter *par, typChannel *channel, FLOAT average_frame[3][5], 
								typInterruption *interruption, typInterr_detect *int_det, 
								INT32 num_of_frames[3], INT32 *found_in_frame, INT32 which_frame,
								INT32 *word_end, INT32 *tested_all_frame);

/*---------------------------------------------------------------------------------------------- */
/* Test interruptions */
/*----------------------------------------------------------------------------------------------	 */
INT32 TestInterruptions(typInputParameter *par, typChannel *channel, FLOAT average_frame[3][5], 
									typInterruption *interruption, typInterr_detect *int_det, 
									INT32 num_of_frames[3]);

/*---------------------------------------------------------------------------------------------- */
/* Get start and end position of the detected interruption (startend==1 start, startend==2 end) */
/*----------------------------------------------------------------------------------------------	 */
INT32 getstartend(typInputParameter *par, typChannel *channel, typInterruption *interruption, 
									typInterr_detect *int_det, INT32 frame_cnt, INT32 startend,
									INT32 which_frame);

/*---------------------------------------------------------------------------------------------- */
/* Test detected interruption for minimum and maximum length and  */
/* write detected interruption into interruption structure. */
/*----------------------------------------------------------------------------------------------	 */
INT32 testwriteInterruption(typInputParameter *par, typChannel *channel, 
								typInterruption *interruption, typInterr_detect *int_det);

/*---------------------------------------------------------------------------------------------- */
/* Removes DC-offset */
/*---------------------------------------------------------------------------------------------- */
void RemoveDC_Offset(FLOAT *in, FLOAT *out, typInterr_detect *int_det, typChannel *channel); 

/*---------------------------------------------------------------------------------------------- */
/* Linear fit at the start- and end of interruption */
/*---------------------------------------------------------------------------------------------- */
FLOAT linearfit(typInterr_detect *int_det, INT32 end_pos, INT32 small_frame);

/*---------------------------------------------------------------------------------------------- */
/* Test for packet-move */
/*---------------------------------------------------------------------------------------------- */
INT32 packetmove(typInterr_detect *int_det, typInterruption *interruption, INT32 str_num,
											INT32 start_end_pos, INT32 start_end, INT32 samp_freq);

/*---------------------------------------------------------------------------------------------- */
/* Set all data of packet-switching test back */
/*---------------------------------------------------------------------------------------------- */
void setpacketback(typInterr_detect *int_det, INT32 beg_end);

/*---------------------------------------------------------------------------------------------- */
/* Get threshold for fundamental frame test */
/*---------------------------------------------------------------------------------------------- */
void getthreshold(FLOAT average_frame, FLOAT *threshold);

/*---------------------------------------------------------------------------------------------- */
/* Get threshold for frame test */
/*---------------------------------------------------------------------------------------------- */
void getthressamll(typInputParameter *par, FLOAT sum_value_prev, FLOAT *sum_value_thres);



/********************* 
* 
* FUNCTION: detectInterruptions 
* 
* DESCRIPTION: 
* 
*	Detects interruptions in an input speech signal
*	It can detect the interruptions from a couple of samples up to 80 ms. 
*	An input signal has to be normalized with the regards to level and its DC offset.
*
*	Main result: SpeechInterruptions
*
*	Input:		*par			Pointer to input parameter (example: sampling frequency)
*				*channel		Pointer to channel parameter (example: datas)
*				*interruption	Pointer to interruption parameter (example: interruption start)
*				*int_det		Pointer to interruption detection structure (see above)
*				*ImpulseNoise	Pointer to impulse noise positions
*				
*	Output:		*interruption	Interruption start, length and type of interruption (Frame move, 
*								no frame move)
*
***********************/ 

INT32 detectInterruptions(typInputParameter *par, typChannel *channel, 
						  typInterruption *interruption, typInterr_detect *int_det) 
{
	INT32	retVal=SQ_NO_ERRORS;
	INT32	word_begin=0;
	INT32	open_cnt=0;
	FLOAT	value_power=0;

	
	/* Reset counter */
	int_det->num_int=0;
	
	/* End of file reached */
	if (channel->FrameCnt>=channel->MaxNrFrames-1) {
		memcpy(int_det->buffer, channel->dataBufferPrev, int_det->num_samples*sizeof(FLOAT));
		memcpy(int_det->buffer+int_det->num_samples, channel->dataBuffer, int_det->num_samples*sizeof(FLOAT));
		vclr (int_det->buffer+int_det->num_samples*2, 1, int_det->num_samples);
	}
	else {
		memcpy(int_det->buffer, channel->dataBufferPrev, int_det->num_samples*sizeof(FLOAT));
		memcpy(int_det->buffer+int_det->num_samples, channel->dataBuffer, int_det->num_samples*sizeof(FLOAT));
		memcpy(int_det->buffer+int_det->num_samples*2, channel->dataBufferNext, int_det->num_samples*sizeof(FLOAT));

		memcpy(int_det->packetmem, int_det->packetmem+int_det->num_samples,	int_det->num_samples*sizeof(FLOAT));
		memcpy(int_det->packetmem+int_det->num_samples, int_det->packetmem+int_det->num_samples*2, int_det->num_samples*sizeof(FLOAT));
		memcpy(int_det->packetmem+int_det->num_samples*2, channel->dataBufferPrev, int_det->num_samples*sizeof(FLOAT));
		memcpy(int_det->packetmem+int_det->num_samples*3, channel->dataBuffer, int_det->num_samples*sizeof(FLOAT));

		/* Test interruption type */
		if ((int_det->open_packet_cnt>0) && (int_det->strpacket[0].end_frame>0)) 
		{
			for (open_cnt=0; open_cnt<int_det->open_packet_cnt; open_cnt++) 
			{
				int_det->strpacket[open_cnt].end_frame--;
				if (int_det->strpacket[open_cnt].end_frame==0) 
				{
					retVal=packetmove(int_det, interruption, 0, int_det->strpacket[open_cnt].end_pos, 2, par->Frequency);
					if (retVal<0)
						return retVal;
					setpacketback(int_det, 1);
					open_cnt--;
				}
			} /* end for */
		}
	}
	
	/* Remove DC-offset */
	RemoveDC_Offset(int_det->buffer, int_det->buffer, int_det, channel);

	/* Get start point of the next speech active frame */
	if (int_det->speech==0) 
	{	
		mvesq(int_det->buffer+int_det->num_samples, 1, &value_power, int_det->num_samples);

		/* Signal greater then 50 dBov */
		if (sqrt(value_power+1e-15)>MAX_16BIT_VALUE*pow(10, (FLOAT)-50/20)) 
		{
			int_det->speech=1;
			word_begin=1;
			retVal=fund_frames(par, channel, interruption, int_det, word_begin);
		}
	}
	else /* Speech active frame has allready been detected */
		retVal=fund_frames(par, channel, interruption, int_det, word_begin);

	return retVal;
}




/********************* 
* 
* FUNCTION: initDetectInterruptions
* 
* DESCRIPTION: 
*	Initialize structure of typ typInterr_detect
*				
* Input:		*int_det	Pointer to structure
*				*par		Pointer to structure to get sampling frequency
*				
* Output:		*int_det	Pointer to initialized structure
*
***********************/ 


INT32 initDetectInterruptions(typInterr_detect *int_det, typInputParameter *par) 
{

	int_det->speech=0;						/* Flag used for speech activity detection */
											/* (0 no speech, 1 speech) */
	int_det->num_samples=(INT32)			/* Number of samples within 16 ms */
			(FRAME_TIME*FRAME_OVERLAP*
			par->Frequency*1000);
	int_det->num_samples_frame=(INT32)		/* Number of samples in fundamental frame */
			(0.005*par->Frequency*1000);
	int_det->no_interruption[0]=1;			/* Interruption status */
	int_det->no_interruption[1]=1;			/* (0 interruption, 1 no interruption) */
	int_det->no_interruption[2]=1;
	int_det->last_position[0]=0;			/* Last checked sample */
	int_det->last_position[1]=0;
	int_det->last_position[2]=0;
	int_det->num_int=0;						/* Number of interruptions in a processing buffer */
	int_det->num_int_all=0;					/* Total number of interruptions */
	int_det->open_packet_cnt=0;				/* Counter for not checked packet-switching */
	int_det->max_open=0;					/* Maximum actual number of not tested interruptions */
	int_det->start_position=0;				/* Start position of interruption before writing into interruption buffer */
	int_det->end_position=0;				/* End position of interruption before writing  into interruption buffer */
	int_det->last_pos=0;					/* End position for calculation of the signal power, */
											/* 1ms after start of interruption */

	int_det->local_end=0;					/* End position of 5ms frames during interruption */
	int_det->frame_offset[0]=0;				/* Data offset in dataBuffer	 */
	int_det->frame_offset[1]=0;
	int_det->frame_offset[2]=0;
	int_det->energy_rate=0;					/* Short-time-energy value during interruption ocurence */
	int_det->est_rms=0;						/* Estimated RMS value to the left- and right side of an interruption */
	int_det->value_level=(FLOAT)			/* Maximum noise value. Used for threshold */
			(MAX_16BIT_VALUE*(FLOAT)pow(10,	
			(FLOAT)MAX_NOISE_DB/20));
	int_det->pow_af_st=0;					/* Signal power 1ms after interruption start */
	int_det->DCint=0;						/* Sum of signal values during interruption for DC removal */
	int_det->DCOffset[0]=0;					/* DC-Offset value of 16ms buffers */
	int_det->DCOffset[1]=0;
	int_det->DCOffset[2]=0;
	int_det->DCOffset[3]=0;
	int_det->DCOffset[4]=0;
	int_det->DCOffset[5]=0;
	int_det->DCOffset[6]=0;
	int_det->DCOffset[7]=0;
	int_det->DCOffset[8]=0;
	int_det->DCOffset[9]=0;
	int_det->DCOffset[10]=0;
	int_det->DCOffset[11]=0;
	int_det->local_pow[0]=0;				/* Power of 5ms long frames during interruption */
	int_det->local_pow[1]=0;
	int_det->packetmem=NULL;				/* Data buffer for packet-switching test */
	
	/* Allocate memory for packet-switching test */
	int_det->packetmem=(FLOAT *) calloc(int_det->num_samples*4,sizeof(FLOAT));

	/* Buffer with 3 times 16ms datas (50% overlapping) */
	int_det->buffer=(FLOAT *) calloc(int_det->num_samples*3, sizeof(FLOAT));	

	if ((int_det->packetmem==NULL) || (int_det->buffer==NULL))
		return SQERR_ALLOC_ERROR;

	return SQ_NO_ERRORS;
}



/********************* 
* 
* FUNCTION: DeinitDetectInterruptions
* 
* DESCRIPTION: 
*	Free memory for initialized structures
*	Input:		*int_det Pointer to structure
*				
*	Output:		*int_det Pointer to deinitialized structure
* 
***********************/ 

void DeinitDetectInterruptions(typInterr_detect *int_det) 
{
	INT32 cnt=0;

	Free((void**)&int_det->buffer);
	Free((void**)&int_det->packetmem);
	if (int_det->max_open>=0) {
		for (cnt=0; cnt<int_det->max_open; cnt++) {
			Free((void**)&int_det->strpacket[cnt].samples_middle);
			Free((void**)&int_det->strpacket[cnt].ceps_left);
		}
	}
}



/********************* 
* 
* FUNCTION: fund_frames
* 
* DESCRIPTION: 
*	Build fundamental frames and averages of fundamental frames
* 
***********************/ 

INT32 fund_frames(typInputParameter *par, typChannel *channel, typInterruption *interruption,
				  typInterr_detect *int_det, INT32 word_begin) 
{
	/* Initialize */
	INT32	num_of_frames[3]={0};	
	INT32	frame_cnt=0;
	INT32	num_of_values[3]={0};
	INT32	num_sam_data[3]={0};
	INT32	retVal=SQ_NO_ERRORS;	
	INT32	begin_cnt=0;
	INT32	max_value=0;
	FLOAT	average_frame[3][5]={0};

	/* Search for word begin */
	if (word_begin==1) {
		for (begin_cnt=int_det->num_samples; begin_cnt<int_det->num_samples*2; begin_cnt++) {
			if (fabs(*(int_det->buffer+begin_cnt))>int_det->value_level*2.095) { 
				int_det->last_position[0]=int_det->num_samples+begin_cnt;
				int_det->last_position[1]=int_det->num_samples+begin_cnt+ 
													(INT32)	floor(int_det->num_samples_frame/3);
				int_det->last_position[2]=int_det->num_samples+begin_cnt+
													(INT32)	floor(int_det->num_samples_frame/3*2);
				break;
			}
		}
		/* No begin found */
		if (begin_cnt==int_det->num_samples*2) {
			int_det->speech=0;
			return SQ_NO_ERRORS;
		}
	}

	/* Number of samples that have not been checked in a previous frame */
	num_sam_data[0]=(INT32)(int_det->num_samples*2-int_det->last_position[0]);
	num_sam_data[1]=(INT32)(int_det->num_samples*2-int_det->last_position[1]);
	num_sam_data[2]=(INT32)(int_det->num_samples*2-int_det->last_position[2]);

	/* Number of fundamental frames in 16 ms plus not checked data in previous frame */
	if (int_det->num_samples_frame!=0) {
		num_of_frames[0]=(INT32)(floor(((int_det->num_samples+num_sam_data[0])/
																int_det->num_samples_frame)))+1;
		num_of_frames[1]=(INT32)(floor(((int_det->num_samples+num_sam_data[1])/
																int_det->num_samples_frame)))+1;
		num_of_frames[2]=(INT32)(floor(((int_det->num_samples+num_sam_data[2])/
																int_det->num_samples_frame)))+1;
	}
	else
		return SQERR_UNKNOWN_ERROR;

	num_of_values[0]=num_of_frames[0]*int_det->num_samples_frame;
	num_of_values[1]=num_of_frames[1]*int_det->num_samples_frame;
	num_of_values[2]=num_of_frames[2]*int_det->num_samples_frame;
	
	/* Frame A */
	int_det->frame_offset[0]=int_det->last_position[0]-
												(int_det->num_samples+int_det->num_samples_frame);
	/* Frame B */
	int_det->frame_offset[1]=int_det->last_position[1]-
												(int_det->num_samples+int_det->num_samples_frame);
	/* Frame C */
	int_det->frame_offset[2]=int_det->last_position[2]-
												(int_det->num_samples+int_det->num_samples_frame);

	/* Save last positions */
	int_det->last_position[0]=(int_det->last_position[0]-int_det->num_samples)+
												(num_of_frames[0]-1)*int_det->num_samples_frame;
	int_det->last_position[1]=(int_det->last_position[1]-int_det->num_samples)+
												(num_of_frames[1]-1)*int_det->num_samples_frame;
	int_det->last_position[2]=(int_det->last_position[2]-int_det->num_samples)+
												(num_of_frames[2]-1)*int_det->num_samples_frame;

	/* Find max number of frames */
	for (frame_cnt=0; frame_cnt<3; frame_cnt++) {
		if (max_value<num_of_frames[frame_cnt])
			max_value=num_of_frames[frame_cnt];
	}

	/* Build of each fundamental frames first the absolute value of each sample, sum and then  */
	/* calculate a mean value */
	for (frame_cnt=0; frame_cnt<max_value; frame_cnt++) {

		if (frame_cnt<num_of_frames[0])
			mvemg(int_det->buffer+int_det->frame_offset[0]+frame_cnt*int_det->num_samples_frame, 
									1, &average_frame[0][frame_cnt], int_det->num_samples_frame);

		if (frame_cnt<num_of_frames[1])
			mvemg(int_det->buffer+int_det->frame_offset[1]+frame_cnt*int_det->num_samples_frame, 
									1, &average_frame[1][frame_cnt], int_det->num_samples_frame);
		if (frame_cnt<num_of_frames[2])
			mvemg(int_det->buffer+int_det->frame_offset[2]+frame_cnt*int_det->num_samples_frame, 
									1, &average_frame[2][frame_cnt], int_det->num_samples_frame);
	}
	
	retVal=TestInterruptions(par, channel, average_frame, interruption, int_det, num_of_frames);
	
	return retVal;
}



							

/********************* 
* 
* FUNCTION: TestInterruptions
* 
* DESCRIPTION: 
*	Test fundamental frames for interruptions
* 
***********************/ 

INT32 TestInterruptions(typInputParameter *par, typChannel *channel, FLOAT average_frame[3][5], 
						typInterruption *interruption, typInterr_detect *int_det, INT32 num_of_frames[3]) 
{
	INT32	retVal=0;
	INT32	retVal2=0;
	INT32	word_end[3]={0, 0, 0};
	INT32	found_in_frame[3]={-1, -1, -1};
	INT32	tested_all_frame[3]={0, 0, 0};

	/* Frame A */
	if ((int_det->no_interruption[1]!=0) && (int_det->no_interruption[2]!=0)) {
		retVal=testframe(par, channel, average_frame, interruption, int_det, num_of_frames, 
						&found_in_frame[0], 0, &word_end[0], &tested_all_frame[0]);
		if (retVal<0)
			return retVal;
	}

	/* Frame B */
	if ((int_det->no_interruption[0]!=0) && (int_det->no_interruption[2]!=0)) {
		if (retVal==INTERR_FOUND)
		{
			retVal2=testframe(par, channel, average_frame, interruption, int_det, num_of_frames, 
						&found_in_frame[0], 1, &word_end[0], &tested_all_frame[1]);
		}
		else
		{
			retVal=testframe(par, channel, average_frame, interruption, int_det, num_of_frames, 
						&found_in_frame[0], 1, &word_end[0], &tested_all_frame[1]);
		}

		if (retVal<0)
			return retVal;
		else if  (retVal2<0)
			return retVal2;
	
		if ((retVal==INTERR_FOUND) && (tested_all_frame[0]==0) && 
															(int_det->no_interruption[1]!=0)) {
			retVal2=testframe(par, channel, average_frame, interruption, int_det, num_of_frames,
						&found_in_frame[0], 0, &word_end[0], &tested_all_frame[0]);
			if (retVal2<0)
				return retVal;
		}
	}

	/* Frame C */
	if ((int_det->no_interruption[0]!=0) && (int_det->no_interruption[1]!=0)) {
		if (retVal==INTERR_FOUND)
			retVal2=testframe(par, channel, average_frame, interruption, int_det, num_of_frames,
						&found_in_frame[0], 2, &word_end[0], &tested_all_frame[2]);
		else
			retVal=testframe(par, channel, average_frame, interruption, int_det, num_of_frames, 
						&found_in_frame[0], 2, &word_end[0], &tested_all_frame[2]);
		if (retVal<0)
			return retVal;
		else if (retVal2<0)
			return retVal2;

		if ((retVal==INTERR_FOUND) && (int_det->no_interruption[2]!=0)) {
			/* Test not tested frames A */
			if (tested_all_frame[0]==0) {
				retVal2=testframe(par, channel, average_frame, interruption, int_det, 
											num_of_frames, &found_in_frame[0], 0, &word_end[0],
											&tested_all_frame[0]);
				if (retVal2<0)
					return retVal;
			}
		   
			if ((tested_all_frame[1]==0) && (int_det->no_interruption[0]!=0)) {
				retVal2=testframe(par, channel, average_frame, interruption, int_det, 
											num_of_frames, &found_in_frame[0], 1, &word_end[0],
											&tested_all_frame[1]);		
				if (retVal2<0)
					return retVal2;
			}
		}
	}
	return retVal;
}



/********************* 
* 
* FUNCTION: testframe
* 
* DESCRIPTION: 
*	Test fundamental frame for interruption
* 
***********************/ 

INT32 testframe(typInputParameter *par, typChannel *channel, FLOAT average_frame[3][5], 
								typInterruption *interruption, typInterr_detect *int_det, 
								INT32 num_of_frames[3], INT32 *found_in_frame, INT32 which_frame,
								INT32 *word_end, INT32 *tested_all_frame) 
{
	INT32	frame_cnt=0;
	INT32	cnt=0;
	INT32	found=-1;
	INT32	where=-1;
	INT32	start_pos=0;
	INT32	end_pos=0;
	INT32	retVal=0;
	INT32	retVal2=0;
	FLOAT	threshold=0;
	FLOAT	energy=0;
	FLOAT	DCsum=0;
	FLOAT	mean_en=0;


	/* More than one interruption found within 16ms. Get start frame to avoid FLOAT testing. */
	if (int_det->num_int!=0) {
		/* Search for last detected interruption */
		for (cnt=0; cnt<3; cnt++) {
			if (found<found_in_frame[cnt]) {
				found=found_in_frame[cnt];
				where=cnt;
			}
		}

		switch (which_frame) {
		case 0:
			if (int_det->frame_offset[where]-int_det->frame_offset[0]>0)
				frame_cnt=found+1;
			else
				frame_cnt=found;
			break;
		case 1:
			if (int_det->frame_offset[where]-int_det->frame_offset[1]>0)
				frame_cnt=found+1;
			else
				frame_cnt=found;
			break;
		case 2:
			if (int_det->frame_offset[where]-int_det->frame_offset[2]>0)
				frame_cnt=found+1;
			else
				frame_cnt=found;
			break;
		}
	}

	for (frame_cnt; frame_cnt<(num_of_frames[which_frame]-1); frame_cnt++) {
		if (int_det->no_interruption[which_frame]!=0) 
		{
			getthreshold(average_frame[which_frame][frame_cnt], &threshold);

			/* No interruption detected yet */
			if ((average_frame[which_frame][frame_cnt]*threshold>
													average_frame[which_frame][frame_cnt+1])
													&& (average_frame[which_frame][frame_cnt]>
													int_det->value_level/1.5)) 
			{
				int_det->no_interruption[which_frame]=0;
				if (retVal==INTERR_FOUND) 
					retVal2=getstartend(par, channel, interruption, int_det, frame_cnt+1, 1, which_frame);
				else
					retVal=getstartend(par, channel, interruption, int_det, frame_cnt+1, 1, which_frame);
				if (retVal<0)
					return retVal;
				else if  (retVal2<0)
					return retVal2;
			}
	
			/* Condition for end of speech activity recognition */
			else if (average_frame[which_frame][frame_cnt+1]<int_det->value_level/1.3) 
			{
				frame_cnt=num_of_frames[which_frame];
				word_end[which_frame]=1;
				if ((word_end[0]==1) && (word_end[1]==1) && (word_end[2]==1)) {
					int_det->speech=0;
					break;
				}
			}			
		}
		else if (int_det->value_level*1.52<average_frame[which_frame][frame_cnt+1]) 
		{
			/* More than one interruption detected */
			if (retVal==INTERR_FOUND) 
			{
				retVal2=getstartend(par, channel, interruption, int_det, frame_cnt, 2, which_frame);
				if (retVal2<0)
					return retVal2;
				if (retVal2==INTERR_FOUND)
					found_in_frame[which_frame]=frame_cnt;
			}
			else 
			{
				retVal=getstartend(par, channel, interruption, int_det, frame_cnt, 2, which_frame);
				if (retVal<0)
					return retVal;
				if (retVal==INTERR_FOUND)
					found_in_frame[which_frame]=frame_cnt;
			}
		}
		else {
			start_pos=int_det->frame_offset[which_frame]+(frame_cnt+1)*int_det->num_samples_frame;
			end_pos=start_pos+int_det->num_samples_frame;

			svesq(int_det->buffer+start_pos, 1, &energy, end_pos-start_pos);
			int_det->energy_rate+=energy;
			sve(int_det->buffer+start_pos, 1, &DCsum, end_pos-start_pos);
			int_det->DCint+=DCsum;
			/* Test if power during interruption is too high */
			mvesq(int_det->buffer+start_pos, 1, &mean_en, int_det->num_samples_frame);
			if (int_det->local_pow[0]<mean_en) {
				int_det->local_pow[1]=int_det->local_pow[0];
				int_det->local_pow[0]=mean_en;
				int_det->local_end=(channel->FrameCnt-2)*(int_det->num_samples)+end_pos;
			}
		}
	}

	*tested_all_frame=1;
	return retVal;
}




/********************* 
* 
* FUNCTION: getstartend
* 
* DESCRIPTION: 
*	Get start and end position of the detected interruption 
*	(startend==1 start, startend==2 end)
* 
***********************/ 

INT32 getstartend(typInputParameter *par, typChannel *channel, typInterruption *interruption, 
				  typInterr_detect *int_det, INT32 frame_cnt, INT32 startend,	INT32 which_frame) 
{
	INT32	position_left=0;
	INT32	position_middle=0;
	INT32	position_cnt=0;
	INT32	big_cnt=0;
	INT32	little_cnt=0;
	INT32	min_cnt=0;
	INT32	retVal=0;
	INT32	start_pos=0;
	INT32	end_pos=0;
	INT32	samp=0;
	INT32	num_values=3*int_det->num_samples;
	INT32	small_frame=(INT32) ((FLOAT) int_det->num_samples_frame)/4;
	INT32	small_interr=(INT32) (0.00625*par->Frequency*1000);
	INT32	true_start_end=0;
	FLOAT	sum_value_thres=0;
	FLOAT	calc_power_rate=0;
	FLOAT	sum_value_prev=0;
	FLOAT	sum_value_next=0;
	FLOAT	act_sum=0;
	FLOAT	minimum=FLT_MAX;
	FLOAT	rms_right=0;
	FLOAT	steep_before=0;
	FLOAT	steep_after=0;
	FLOAT	energy=0;
	FLOAT	DCsum=0;
	FLOAT	mean_start_end;
	FLOAT	start_end_power=(FLOAT)(MAX_16BIT_VALUE * (FLOAT)pow(10, (FLOAT) (MAX_NOISE_DB+5)/20));
	FLOAT	thres_power=0;

	if (startend==1) {
		position_left=int_det->frame_offset[which_frame]+frame_cnt*int_det->num_samples_frame;

		/* Search for lowest sum(abs(values)) in splitted fundamental frame to find start position */
		if (small_frame!=0)
			samp=(INT32)floor((FLOAT)int_det->num_samples_frame/small_frame);
		else
			return SQERR_UNKNOWN_ERROR;
		for (big_cnt=0; big_cnt<samp; big_cnt++) {
			act_sum=0;

			for (little_cnt=0; little_cnt<small_frame; little_cnt++) {
				act_sum+=(FLOAT)fabs(*(int_det->buffer+position_left+big_cnt*small_frame+little_cnt));
			}

			if (act_sum<minimum) {
				minimum=act_sum;
				min_cnt=big_cnt;
			}
		}
		position_middle=int_det->frame_offset[which_frame]+frame_cnt*
												int_det->num_samples_frame+min_cnt*small_frame+
												(INT32) ((FLOAT)small_frame*3.0/4);

		/* Move towards negativ time axis */
		for (position_cnt=position_middle; position_cnt>=(FLOAT) small_frame/2; position_cnt--) {
			mvesq(int_det->buffer+position_cnt-(INT32) ((FLOAT) small_frame/2)-1, 1, 
												&mean_start_end, (INT32) (FLOAT) small_frame/2);
			if (sqrt(mean_start_end)>start_end_power) {

				for (little_cnt=position_cnt; little_cnt>=0; little_cnt--) {
					if ((fabs(*(int_det->buffer+little_cnt))+
														fabs(*(int_det->buffer+little_cnt+1)))/2>
														int_det->value_level*2.06) {

						if (fabs(*(int_det->buffer+little_cnt+1))>int_det->value_level*2.06) {
							int_det->start_position=(channel->FrameCnt-2)
														*(int_det->num_samples)+little_cnt+2;
							little_cnt+=1;
							break;
						}
						else if (fabs(*(int_det->buffer+little_cnt))>int_det->value_level*2.06) {
							int_det->start_position=(channel->FrameCnt-2)
														*(int_det->num_samples)+little_cnt+1;

							break;
						}
					}
				}

				if ((little_cnt>=small_frame) && (little_cnt<num_values-small_frame)) {
					mvemg(int_det->buffer+little_cnt-small_frame, 1, &sum_value_prev, 
																					small_frame);
					mvemg(int_det->buffer+little_cnt+1, 1, &sum_value_next, small_frame);

					getthressamll(par, sum_value_prev, &sum_value_thres);

					/* Test if it is an interruption or end of word */
					if ((sum_value_prev<sum_value_next*sum_value_thres) || 
														(sum_value_prev<int_det->value_level)) {
						int_det->no_interruption[which_frame]=1;
						return INTERR_NOT_FOUND;
					}

					start_pos=little_cnt+1;
					end_pos=int_det->frame_offset[which_frame]+
								frame_cnt*int_det->num_samples_frame+int_det->num_samples_frame;

					if (little_cnt>=int_det->num_samples_frame) {
						/* Calculate estimated RMS left of an interruption */
						rmvesq((int_det->buffer+little_cnt-int_det->num_samples_frame), 1, 
												&(int_det->est_rms), int_det->num_samples_frame);
					}
					else {
						int_det->est_rms=0;
					}
					
					int_det->energy_rate=0;
					/* Calculate energy */
					svesq(int_det->buffer+start_pos, 1, &energy, end_pos-start_pos);

					int_det->energy_rate+=energy;
					/* Calculate sum of vector values during interruption for DC removal */
					sve(int_det->buffer+start_pos, 1, &DCsum, end_pos-start_pos);
					int_det->DCint+=DCsum;

					/* Calculate power 1 ms after start */
					int_det->pow_af_st=0;
					if (num_values-small_interr>start_pos) {
						mvesq(int_det->buffer+start_pos+small_frame, 1, &int_det->pow_af_st, 
																	int_det->num_samples_frame);

						int_det->last_pos=(channel->FrameCnt-2)
														*(int_det->num_samples)+little_cnt+1+
														small_frame+int_det->num_samples_frame;
					}

					/* Steepness at the start position */
					steep_before=linearfit(int_det, start_pos-small_frame, small_frame);
					steep_after=linearfit(int_det, start_pos, small_frame);
					if ((fabs(steep_before)>58) && (fabs(steep_after)>10) && 
																(steep_before*steep_after>0)) {
						int_det->start_position=0;
						int_det->no_interruption[which_frame]=1;
						return INTERR_NOT_FOUND;
					}

					/* Initialization for local power test during interruption */
					int_det->local_pow[1]=0;
					int_det->local_end=(channel->FrameCnt-2)*(int_det->num_samples)+end_pos;
					mvesq(int_det->buffer+start_pos, 1, &int_det->local_pow[0], 
																			end_pos-start_pos);
					
					/* Test for packet-move */
					if (int_det->packetmem!=NULL) 
					{
						/* Allocate memory if necessary */
						if (int_det->max_open<=int_det->open_packet_cnt) {
							int_det->max_open++;
							int_det->strpacket[int_det->open_packet_cnt].end_pos=0;
							int_det->strpacket[int_det->open_packet_cnt].end_frame=0;
							int_det->strpacket[int_det->open_packet_cnt].samples_middle=
										(FLOAT *) calloc(2*int_det->num_samples, sizeof(FLOAT));
							int_det->strpacket[int_det->open_packet_cnt].ceps_left=
										(FLOAT *) calloc(2*int_det->num_samples, sizeof(FLOAT));

							if ((int_det->strpacket[int_det->open_packet_cnt].samples_middle
														==NULL) || (int_det->strpacket[int_det->
														open_packet_cnt].ceps_left==NULL))
								return SQERR_ALLOC_ERROR;
						}

						retVal=packetmove(int_det, interruption, int_det->open_packet_cnt,
																	start_pos, 1, par->Frequency);
						if (retVal<0)
							return retVal;
						int_det->strpacket[int_det->open_packet_cnt].inter_number=
																	int_det->num_int_all;							
						int_det->strpacket[int_det->open_packet_cnt].end_frame=-1;
						int_det->open_packet_cnt++;
					}
					return SQ_NO_ERRORS;
				}
				/* Start is out of range */
				else {
					int_det->start_position=0;
					int_det->no_interruption[which_frame]=1;
					return INTERR_NOT_FOUND;
				}				
			}
		}
		/* No start detected */
		if (position_cnt==(INT32) ((FLOAT) small_frame/2)-1) {
			int_det->no_interruption[which_frame]=1;
			return INTERR_NOT_FOUND;
		}
	}
	else if (startend==2) 
	{
		/* Get possible start position */
		position_left=int_det->frame_offset[which_frame]+frame_cnt*int_det->num_samples_frame;

		/* Search for lowest sum(abs(values)) in splitted fundamental frame to find end position */
		if (small_frame!=0)
			samp=(INT32)floor((FLOAT)int_det->num_samples_frame/small_frame);
		else
			return SQERR_UNKNOWN_ERROR;
		for (big_cnt=0; big_cnt<samp; big_cnt++) {
			act_sum=0;

			for (little_cnt=0; little_cnt<small_frame; little_cnt++) {
				act_sum+=(FLOAT)fabs(*(int_det->buffer+position_left+big_cnt*small_frame+little_cnt));
			}
			if (act_sum<minimum) {
				minimum=act_sum;
				min_cnt=big_cnt;
			}
		}
		position_middle=int_det->frame_offset[which_frame]+
										frame_cnt*int_det->num_samples_frame+min_cnt*small_frame+
										(INT32) ((FLOAT)small_frame/4);
	
		for (position_cnt=position_middle; position_cnt<num_values-(FLOAT) small_frame/2; position_cnt++) 
		{
			/* Move towards positiv time axis */
			mvesq(int_det->buffer+position_cnt, 1, &mean_start_end, 
																(INT32) ((FLOAT) small_frame/2));
			if (sqrt(mean_start_end)>start_end_power) {
				for (little_cnt=position_cnt; little_cnt<num_values; little_cnt++) {
					if ((fabs(*(int_det->buffer+little_cnt))+
														fabs(*(int_det->buffer+little_cnt-1)))/2>
														int_det->value_level*2.06) {
						if (fabs(*(int_det->buffer+little_cnt-1))>int_det->value_level*2.06) {
							int_det->end_position=
										(channel->FrameCnt-2)*(int_det->num_samples)+little_cnt-1;
							little_cnt-=1;
							true_start_end=1;
							break;
						}
						else if (fabs(*(int_det->buffer+little_cnt))>int_det->value_level*2.06) {
							int_det->end_position=
										(channel->FrameCnt-2)*(int_det->num_samples)+little_cnt;
							true_start_end=1;
							break;
						}
					}
				}

				if ((little_cnt>=small_frame) && (little_cnt<num_values-small_frame)) {
					start_pos=int_det->frame_offset[which_frame]+
														(frame_cnt+1)*int_det->num_samples_frame;
					int_det->no_interruption[which_frame]=1;

					end_pos=little_cnt;
					
					/* Start position is less than end position */
					if (start_pos<end_pos) {
						svesq(int_det->buffer+start_pos, 1, &energy, end_pos-start_pos);
						int_det->energy_rate+=energy;
						/* Calculate sum of vector values during interruption for DC removal */
						sve(int_det->buffer+start_pos, 1, &DCsum, end_pos-start_pos);
						int_det->DCint+=DCsum;
					}

					else if (start_pos>end_pos) {
						svesq(int_det->buffer+end_pos, 1, &energy, start_pos-end_pos);
						int_det->energy_rate-=energy;
						/* Calculate sum of vector values during interruption for DC removal */
						sve(int_det->buffer+end_pos, 1, &DCsum, start_pos-end_pos);
						int_det->DCint-=DCsum;

					}

					/* Test if interruption result is OK */
					if ((int_det->end_position-int_det->start_position)>0) {

						/* Calculate power during interruption */
						calc_power_rate=(FLOAT) int_det->energy_rate/
												(int_det->end_position-int_det->start_position);

						/* Hightest local power during interruption interval */
						if (int_det->local_end>int_det->end_position)
							int_det->local_pow[0]=int_det->local_pow[1];

						if (int_det->last_pos>int_det->end_position)
							int_det->pow_af_st=0;

						/* Test if power during interruption is too high (greater than maximum  */
						/* noise +3 [dB]) and power 1 ms after interruption start */
						thres_power=(FLOAT)pow(MAX_16BIT_VALUE*pow(10, ((FLOAT)MAX_NOISE_DB+3)/20), 2);
						if ((calc_power_rate<thres_power)	&& 
									(int_det->pow_af_st<thres_power) &&
									(int_det->local_pow[0]<thres_power)) {

							if (little_cnt<num_values-int_det->num_samples_frame) {
								/* Calculates estimated RMS right of an interruption position */
								rmvesq(int_det->buffer+little_cnt, 1, &rms_right, 
																	int_det->num_samples_frame);
				
								/* Calculates an estimated rms as an average of left and right values */
								int_det->est_rms=(rms_right+int_det->est_rms)/2;
							}
							else {
								int_det->est_rms=0;
							}

							/* Special test for critical interruptions using linear fitting */
							if ((int_det->end_position-int_det->start_position<small_interr) &&
																	(end_pos-small_frame>=0)) {
								steep_before=linearfit(int_det, end_pos-small_frame, small_frame);
								steep_after=linearfit(int_det, end_pos, small_frame);

								switch (par->Frequency) {
								case 8:
									if ((fabs(steep_before)>18) && (fabs(steep_after)>18) && 
																(steep_before*steep_after>0)) {
										int_det->start_position=0;
										int_det->energy_rate=0;
										setpacketback(int_det, 2);
										return INTERR_NOT_FOUND;
									}
									break;
								case 16:
									if ((fabs(steep_before)>11) && (fabs(steep_after)>59) && 
																(steep_before*steep_after>0)) {
										int_det->start_position=0;
										int_det->energy_rate=0;
										setpacketback(int_det, 2);
										return INTERR_NOT_FOUND;
									}
									break;
								}
							}
							retVal=testwriteInterruption(par, channel, interruption, int_det);
							/* DC removal of interruption values */
							if ((retVal==INTERR_FOUND) && (int_det->end_position-
										int_det->start_position>=int_det->num_samples_frame-
										small_frame)) {
								
								/* Get DC-Offset old  */
								mve(&(int_det->DCOffset[0]), 1, &DCsum, 12);

								if (int_det->end_position-int_det->start_position!=0)
									int_det->DCint/=(int_det->end_position-int_det->start_position);
								else
									return SQERR_UNKNOWN_ERROR;
								int_det->DCint+=DCsum;
								int_det->DCOffset[0]=int_det->DCint;
								int_det->DCOffset[1]=int_det->DCint;
								int_det->DCOffset[2]=int_det->DCint;
								int_det->DCOffset[3]=int_det->DCint;
								int_det->DCOffset[4]=int_det->DCint;
								int_det->DCOffset[5]=int_det->DCint;
								int_det->DCOffset[6]=int_det->DCint;
								int_det->DCOffset[7]=int_det->DCint;
								int_det->DCOffset[8]=int_det->DCint;
								int_det->DCOffset[9]=int_det->DCint;
								int_det->DCOffset[10]=int_det->DCint;
								int_det->DCOffset[11]=int_det->DCint;
							}
							else
								int_det->DCint=0;

							/* Test packet-move */
							if ((retVal==INTERR_FOUND) && (int_det->packetmem!=NULL)) 
							{
								if(int_det->open_packet_cnt > 0)
								{
									if (little_cnt<int_det->num_samples*2)
										int_det->strpacket[int_det->open_packet_cnt-1].end_frame=2;
									else
										int_det->strpacket[int_det->open_packet_cnt-1].end_frame=3;

									int_det->strpacket[int_det->open_packet_cnt-1].end_pos=end_pos;
								}
							}
							else 
								setpacketback(int_det, 2);
						}
						else {
							int_det->start_position=0;
							int_det->energy_rate=0;
							int_det->DCint=0;
							setpacketback(int_det, 2);
						}

					}
					/* Interruption starts and ends at the same position */
					else {
						int_det->start_position=0;
						int_det->energy_rate=0;
						int_det->DCint=0;
						setpacketback(int_det, 2);
					}
				}
				else {
					int_det->start_position=0;
					int_det->no_interruption[which_frame]=1;
					setpacketback(int_det, 2);
					return INTERR_NOT_FOUND;
				}
				break;
			}
		}
	} 

	return retVal;
}




/********************* 
* 
* FUNCTION: testwriteInterruption
* 
* DESCRIPTION: 
*	Test detected interruption for minimum and maximum length and if an impulse is near and 
*	writes detected interruption into interruption structure.
* 
***********************/ 

INT32 testwriteInterruption(typInputParameter *par, typChannel *channel, 
								typInterruption *interruption, typInterr_detect *int_det) 
{
	/* Initialize */
	INT32 min_length=(INT32) (0.0035*par->Frequency*1000);
	INT32 interr_length=0;

	/* Get length of detected interruption */
	interr_length=int_det->end_position-int_det->start_position;

	if ((interr_length > min_length) && ((interruption[int_det->num_int].start+                  
									interruption[int_det->num_int].length-int_det->start_position<0) 
									|| int_det->num_int==0) &&                                         
									(interr_length<channel->FFTSize*10))                               
		{                                                                                              
			interruption[int_det->num_int].start=int_det->start_position;                              
			interruption[int_det->num_int].length=interr_length;                                       
			interruption[int_det->num_int].EstimatedRms=int_det->est_rms;                              
			interruption[int_det->num_int].packetmove=-1;                                              
			int_det->num_int++;                                                                        
			int_det->num_int_all++;                                                                    
		                                                                                               
		                                                                                               
			return INTERR_FOUND;                                                                       
		}                                                                                              
		else                                                                                           
			return INTERR_NOT_FOUND;                                                                   
}



/********************* 
* 
* FUNCTION: RemoveDC_Offset
* 
* DESCRIPTION: 
*	Remove DC-offset
* 
***********************/ 

void RemoveDC_Offset(FLOAT *in, FLOAT *out, typInterr_detect *int_det, typChannel *channel) {
	
	/* Initialize */
	FLOAT DC_new=0;
	FLOAT DC_old=0;

	/* Get old mean value */
	mve(&(int_det->DCOffset[0]), 1, &DC_old, 12);

	/* Get DC-Offset of data buffer next */
	mve(in+int_det->num_samples*2, 1, &DC_new, int_det->num_samples);

	/* only then if above 50 (50=1.5259e-3*32767) */
	if ((fabs(DC_old-DC_new)>0.0015259*MAX_16BIT_VALUE) && (channel->FrameCnt>13)) 
		DC_new+=((DC_old-DC_new)*(FLOAT)0.65);

	int_det->DCOffset[(channel->FrameCnt % 12)]=DC_new;

	/* Get mean value */
	mve(&(int_det->DCOffset[0]), 1, &DC_new, 12);

	if (fabs(DC_new)>5)
		vsadd(in, 1, -DC_new, out, 1, int_det->num_samples*3);
}



/********************* 
* 
* FUNCTION: linearfit
* 
* DESCRIPTION: 
*	Linear fit at the start and at the end of an interruption
* 
***********************/ 

FLOAT linearfit(typInterr_detect *int_det, INT32 end_pos, INT32 small_frame) {

	/* Initialize */
	INT32  cnt=0;
	FLOAT step=0;
	FLOAT sumxy=0;
	FLOAT sumx=0;
	FLOAT sumy=0;
	FLOAT squaresumx=0;
	FLOAT divider=0;

	for (cnt=0; cnt<small_frame; cnt++) {
		sumx+=cnt;
		sumxy+=(cnt**(int_det->buffer+end_pos+cnt));
		squaresumx+=cnt*cnt;
	}
	/* Sum of vector elements */
	sve(int_det->buffer+end_pos, 1, &sumy, small_frame);

	divider=(small_frame*squaresumx-(FLOAT)pow(sumx, 2));

	if (divider!=0)
		step=(small_frame*sumxy-sumx*sumy)/divider;
	else
		return SQERR_UNKNOWN_ERROR;
	return step;
}


/********************* 
* 
* FUNCTION: packetmove
* 
* DESCRIPTION: 
*	Test for packet-move
* 
***********************/ 

INT32 packetmove(typInterr_detect *int_det, typInterruption *interruption, INT32 str_num,
										INT32 start_end_pos, INT32 start_end, INT32 samp_freq) {
	
	/* Initialize */
	INT32	start_cep=(INT32) floor(0.5625*samp_freq);
	INT32	half_samp=(INT32) ((FLOAT) int_det->num_samples/2);
	FLOAT	correlCoeff_left_middle=0;
	FLOAT	correlCoeff_right_middle=0;
	FLOAT	*help_mem=NULL;
	FLOAT	*ceps_right=NULL;
	FLOAT	*ceps_middle=NULL;

	if (start_end==1) {

		help_mem=(FLOAT *) calloc(2*int_det->num_samples, sizeof(FLOAT));

		if (help_mem==NULL)
			return SQERR_ALLOC_ERROR;

		if(start_end_pos+int_det->num_samples+half_samp > int_det->num_samples*2)
			return SQERR_ALLOC_ERROR;


		
		/* Memory copy for overlapping test (to left and to right of interruption) */
			memcpy(int_det->strpacket[str_num].samples_middle, 
								int_det->packetmem+start_end_pos+int_det->num_samples+half_samp, 
								int_det->num_samples*sizeof(FLOAT));

		/* Multiply data with hamming window */
		hamm (int_det->packetmem+start_end_pos+int_det->num_samples, 1, help_mem, 1, 
																			int_det->num_samples);
		/* Calculate cepstrum of data left of interruption */
		cepstrum(help_mem, int_det->strpacket[str_num].ceps_left, int_det->num_samples);
	}
	else if (start_end==2) {

		help_mem=(FLOAT *) calloc(2*int_det->num_samples, sizeof(FLOAT));
		ceps_right=(FLOAT *) calloc(2*int_det->num_samples, sizeof(FLOAT));
		ceps_middle=(FLOAT *) calloc(2*int_det->num_samples, sizeof(FLOAT));

		if ((help_mem==NULL) || (ceps_right==NULL) || (ceps_middle==NULL))
			return SQERR_ALLOC_ERROR;

		/* Copy to data middle */
		memcpy(&(int_det->strpacket[str_num].samples_middle[half_samp]), 
									int_det->packetmem+start_end_pos, (half_samp)*sizeof(FLOAT));

		/* Multiply data right of interruption with hamming window */
		hamm (int_det->packetmem+start_end_pos, 1, help_mem, 1, int_det->num_samples);

		/* Calculate cepstrum of data right of interruption */
		cepstrum(help_mem, ceps_right, int_det->num_samples);
		
		/* Multiply data middle of interruption with hamming window */
		hamm (int_det->strpacket[str_num].samples_middle, 1, 
						int_det->strpacket[str_num].samples_middle, 1, int_det->num_samples);
		
		cepstrum(int_det->strpacket[str_num].samples_middle, ceps_middle, int_det->num_samples);
		
		/* Data left middle */
		mcorrel (&int_det->strpacket[0].ceps_left[start_cep], &ceps_middle[start_cep], 
												&correlCoeff_left_middle, half_samp-start_cep);
		/* Data right middle */
		mcorrel (&ceps_right[start_cep], &ceps_middle[start_cep], &correlCoeff_right_middle,
																			half_samp-start_cep);

		if(int_det->strpacket[str_num].inter_number-int_det->num_int_all >= 0)
		{
			/* Test Discontinuation                                                                      */
			switch (samp_freq) {                                                                        
			case 8:                                                                                     
				if ((correlCoeff_left_middle+correlCoeff_right_middle)/2>0.6617)                        
					interruption[int_det->strpacket[str_num].inter_number-int_det->num_int_all].        
																			packetmove=FRAME_MOVE;      
				else                                                                                    
					interruption[int_det->strpacket[str_num].inter_number-int_det->num_int_all].        
																			packetmove=NO_FRAME_MOVE;   
				break;                                                                                  
			case 16:                                                                                    
				if ((correlCoeff_left_middle+correlCoeff_right_middle)/2>0.58832)                       
					interruption[int_det->strpacket[str_num].inter_number-int_det->num_int_all].        
																			packetmove=FRAME_MOVE;      
				else                                                                                    
					interruption[int_det->strpacket[str_num].inter_number-int_det->num_int_all].        
																			packetmove=NO_FRAME_MOVE;   
				break;                                                                                  
			}
		}

		Free((void**)&ceps_right);
		Free((void**)&ceps_middle);
	}

	Free((void**)&help_mem);
	return SQ_NO_ERRORS;

}


										
/********************* 
* 
* FUNCTION: setpacketback
* 
* DESCRIPTION: 
*	Set all data of packet-switching test back
* 
***********************/ 

void setpacketback(typInterr_detect *int_det, INT32 beg_end) {

	/* Initialize */
	INT32	open_move_cnt=0;

	if (int_det->open_packet_cnt>0) {
		if (beg_end==1) {
			for (open_move_cnt=0; open_move_cnt<int_det->open_packet_cnt-1; open_move_cnt++) {
				int_det->strpacket[open_move_cnt].end_frame=
												int_det->strpacket[open_move_cnt+1].end_frame;
				int_det->strpacket[open_move_cnt].end_pos=
												int_det->strpacket[open_move_cnt+1].end_pos;
				int_det->strpacket[open_move_cnt].inter_number=
												int_det->strpacket[open_move_cnt+1].inter_number;

				/* Copy fft_left */
				memcpy(&(int_det->strpacket[open_move_cnt].ceps_left[0]), 
										&(int_det->strpacket[open_move_cnt+1].ceps_left[0]),
										int_det->num_samples*sizeof(FLOAT));
				/* Copy samples middle */
				memcpy(&(int_det->strpacket[open_move_cnt].samples_middle[0]), 
										&(int_det->strpacket[open_move_cnt+1].samples_middle[0]),
										int_det->num_samples*sizeof(FLOAT));
			}
		}
		int_det->open_packet_cnt--;							
	}
}



/********************* 
* 
* FUNCTION: getthreshold
* 
* DESCRIPTION: 
*	Get threshold for fundamental frame test
* 
***********************/ 

/*---------------------------------------------------------------------------------------------- */
void getthreshold(FLOAT average_frame, FLOAT *threshold) {

	if (average_frame>0.045778*MAX_16BIT_VALUE)
		*threshold=0.625;
	else if ((average_frame<=0.045778*MAX_16BIT_VALUE) && 
														(average_frame>0.025483*MAX_16BIT_VALUE))
		*threshold=0.25;
	else if ((average_frame<=0.025483*MAX_16BIT_VALUE) && 
														(average_frame>0.010681*MAX_16BIT_VALUE))
		*threshold=(FLOAT)0.117;
	else if ((average_frame<=0.010681*MAX_16BIT_VALUE) && 
														(average_frame>0.006104*MAX_16BIT_VALUE))
		*threshold=(FLOAT)0.093;
	else
		*threshold=(FLOAT)0.088;
}


/********************* 
* 
* FUNCTION: getthressamll
* 
* DESCRIPTION: 
*	Get threshold for small frame test
* 
***********************/ 

void getthressamll(typInputParameter *par, FLOAT sum_value_prev, FLOAT *sum_value_thres) {

	switch (par->Frequency) {
	case 8:
		if (sum_value_prev>0.024415*MAX_16BIT_VALUE)
			*sum_value_thres=(FLOAT)12.55;
		else if ((sum_value_prev<=0.024415*MAX_16BIT_VALUE) && 
														(sum_value_prev>0.021668*MAX_16BIT_VALUE))
			*sum_value_thres=(FLOAT)13.96;						 
		else if ((sum_value_prev<=0.021668*MAX_16BIT_VALUE) && 
														(sum_value_prev>0.019532*MAX_16BIT_VALUE))
			*sum_value_thres=(FLOAT)11.41;
		else if ((sum_value_prev<=0.019532*MAX_16BIT_VALUE) && 
														(sum_value_prev>0.009491*MAX_16BIT_VALUE))
			*sum_value_thres=(FLOAT)12.87;
		else
			*sum_value_thres=(FLOAT)10.9;
		break;
	case 16:
		if (sum_value_prev>0.024415*MAX_16BIT_VALUE)
			*sum_value_thres=(FLOAT)12.55;
		else if ((sum_value_prev<=0.024415*MAX_16BIT_VALUE) && 
														(sum_value_prev>0.021668*MAX_16BIT_VALUE))
			*sum_value_thres=(FLOAT)13.96;						 
		else if ((sum_value_prev<=0.021668*MAX_16BIT_VALUE) && 
														(sum_value_prev>0.018769*MAX_16BIT_VALUE))
			*sum_value_thres=(FLOAT)12.62;
		else if ((sum_value_prev<=0.018769*MAX_16BIT_VALUE) && 
														(sum_value_prev>0.009156*MAX_16BIT_VALUE))
			*sum_value_thres=(FLOAT)12.59;
		else if ((sum_value_prev<=0.009156*MAX_16BIT_VALUE) && 
														(sum_value_prev>0.006104*MAX_16BIT_VALUE))
			*sum_value_thres=(FLOAT)12.38;
		else
			*sum_value_thres=(FLOAT)13.63;
		break;
	}
}

