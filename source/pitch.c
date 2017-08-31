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

#include <math.h>
#include <stdlib.h>
#include <memory.h>

#include "defines.h"
#include "dsp.h"

#include "pitch.h"
#include "tools.h"



/*********************
*
*	FUNCTION: pitch_extract
*
*	DESCRIPTION:
*		This function processes the synchronous pitch extraction algorithm.
*		It calculates the instant pitch length, makes the voiced/unvoiced decision 
*		and places pitch marks.
*
***********************/
void pitch_extract(const INT16 *psNormSpeech, const UINT32 ulNormSpeechSz,
					   const FLOAT *pfVad, const UINT32 ulVadSz,
					   INT32 *pitch, INT32 *pitch_length, INT32 pitch_allocated, INT16 *average_pitch)
{

	INT32 neg_loop, index_loop ;	
	INT32 num_frames ;			
	INT32 length_threshold ;		
	INT32 av_loop ;				
	FLOAT av_tmp = 0.0;		
	INT32 num_voiced = 0	;		

	INT32 *pos_pitch = NULL;	
	INT32 *neg_pitch = NULL;	
	INT32 pos_pitch_length;	
	INT32 neg_pitch_length;	
	INT32 *period;			
	INT32 *voiced;			
	INT32 *labels;
	INT16 *psNormSpeechCopy=NULL;
	UINT32 i;

	num_frames = 2 * (ulNormSpeechSz / APITCH_FRAME_SIZE) - 1 ;

	pos_pitch = (INT32 *) calloc(pitch_allocated, sizeof(INT32)) ;
	neg_pitch = (INT32 *) calloc(pitch_allocated, sizeof(INT32)) ;
	period = (INT32 *) calloc(num_frames, sizeof(INT32)) ;
	voiced = (INT32 *) calloc(num_frames, sizeof(INT32)) ;
	labels = (INT32 *) calloc(ulNormSpeechSz, sizeof(INT32)) ;
	psNormSpeechCopy = (INT16*) malloc(ulNormSpeechSz * sizeof(INT16));

	for (i=0; i<ulNormSpeechSz; i++) psNormSpeechCopy[i] = psNormSpeech[i];

	pos_pitch_length = 0 ;	
	neg_pitch_length = 0 ;

	apitch_run(psNormSpeechCopy, ulNormSpeechSz, period, voiced, num_frames ) ;

	split_voiced(voiced, num_frames, pfVad, ulVadSz) ;

	for(av_loop = 0; av_loop < num_frames; av_loop++)
	{
		if(voiced[av_loop] != 0)
		{
			av_tmp = av_tmp + period[av_loop] ;
			num_voiced++ ;
		}
	}

	*average_pitch = (INT16)(av_tmp / num_voiced);

	if (*average_pitch)
	{
		apitch_label(psNormSpeechCopy, labels, ulNormSpeechSz, period, voiced, num_frames ) ;

		for(index_loop = 0 ; index_loop < (INT32) ulNormSpeechSz; index_loop ++)
		{
			if(labels[index_loop] == 1)
			{	
				pos_pitch[pos_pitch_length] = index_loop ;
				pos_pitch_length = pos_pitch_length + 1 ;
			}
		}

		inconsistency_handler(pos_pitch, &pos_pitch_length, pitch_allocated, psNormSpeechCopy, 
							  ulNormSpeechSz, *average_pitch) ;

		length_threshold =  ((*average_pitch) / PITCH_MED ) * PITCH_NOM  ;

		extend_pitch_controller(psNormSpeechCopy, ulNormSpeechSz, pos_pitch, &pos_pitch_length,
								length_threshold, pitch_allocated) ;

		for(neg_loop = 0; neg_loop <(INT32) ulNormSpeechSz; neg_loop++)
			psNormSpeechCopy[neg_loop] = (INT16) -psNormSpeechCopy[neg_loop] ;

		apitch_run(psNormSpeechCopy, ulNormSpeechSz, period, voiced, num_frames ) ;

		split_voiced(voiced, num_frames, pfVad, ulVadSz) ;

		apitch_label(psNormSpeechCopy, labels, ulNormSpeechSz, period, voiced, num_frames ) ;

		for(index_loop = 0 ; index_loop < (INT32) ulNormSpeechSz; index_loop ++)
			if(labels[index_loop] == 1)
			{
				neg_pitch[neg_pitch_length] = index_loop ;
				neg_pitch_length = neg_pitch_length + 1 ;
			}

		inconsistency_handler(neg_pitch, &neg_pitch_length, pitch_allocated, psNormSpeechCopy, 
							  ulNormSpeechSz, *average_pitch) ;

		extend_pitch_controller(psNormSpeechCopy, ulNormSpeechSz, neg_pitch, &neg_pitch_length, 
								length_threshold, pitch_allocated) ;

		for(neg_loop = 0; neg_loop < (INT32) ulNormSpeechSz; neg_loop++)
			psNormSpeechCopy[neg_loop] = (INT16) -psNormSpeechCopy[neg_loop] ;

		pitch_allocated = neg_pitch_length + pos_pitch_length ;

		if (pitch_allocated) pos_or_neg(*average_pitch, psNormSpeechCopy, pos_pitch, pos_pitch_length, 
				   neg_pitch, neg_pitch_length, pitch, pitch_length) ;
		else *pitch_length=0;
	}

	free(psNormSpeechCopy);
	free(pos_pitch);
	free(neg_pitch);
	free(period);
	free(voiced);
	free(labels);
}

/*********************
*
*	FUNCTION: apitch_run
*
*	DESCRIPTION:
*		Compute pitch period and voiced/unvoiced decision.
*		Nframes should be 2 * (Nsamples / APITCH_FRAME_SIZE) - 1.
*
*
***********************/
void apitch_run( const INT16 *psSpeech, const UINT32 ulSpeechSz, INT32 *period, INT32 *voiced, INT32 Nframes )
{
    FLOAT *Window;

    INT32 frame, band;

    FLOAT R[APITCH_FRAME_SIZE_PLUS2];
    FLOAT N[APITCH_FRAME_SIZE_2];
    FLOAT power;

    FLOAT Nmax;
    INT32 Tmax, Told;

	Window = (FLOAT *) calloc( APITCH_FRAME_SIZE, sizeof(FLOAT) );

    Nframes = min( Nframes, (INT32) (2 * (ulSpeechSz / APITCH_FRAME_SIZE) - 1 ));

	for( band = 0; band < APITCH_FRAME_SIZE; band++ )
		Window[band] = (FLOAT)(0.5 * (1.0 - cos( (TWOPI * band) / (FLOAT)APITCH_FRAME_SIZE)) );

    Told = 0L;
    for( frame = 0; frame < Nframes; frame++ ) 
	{
        for( band = 0; band < APITCH_FRAME_SIZE; band++ ) 
			R[band] = Window[band] * (FLOAT)psSpeech[frame*APITCH_FRAME_SIZE_2 + band];

        RealFFT( R, APITCH_FRAME_SIZE );
        for( band = 0; band <= APITCH_FRAME_SIZE_2; band++ ) 
		{
            R[band*2] = R[band*2] * R[band*2] + R[band*2+1] * R[band*2+1];
            R[band*2+1] = 0.0;
        }
        RealIFFT( R, APITCH_FRAME_SIZE );
        power = (FLOAT) max( R[0], 1.0e-10 );

        for( band = 0; band < APITCH_ROLL; band++ ) 
			N[band] = R[band] * Window[band * (APITCH_FRAME_SIZE_2/APITCH_ROLL)] / power;
        for( band = APITCH_ROLL; band < APITCH_FRAME_SIZE_2; band++ ) N[band] = R[band] / power;

        Nmax = 0.0;
        Tmax = 0;
        for( band = 20; band < 100; band++ ) if( N[band] > Nmax ) 
		{
            Nmax = N[band];
            Tmax = band;
        }

        if( Nmax <= 0.5 ) voiced[frame] = 0; 
        else
		{
			voiced[frame] = 1;

			if( Told > 0L ) 
			{
				while( (Told > 24L) && (N[Told-1] > N[Told]) ) Told--;
				while( (Told < 94L) && (N[Told+1] > N[Told]) ) Told++;
				if( N[Told] > (1.0 + Nmax) / 3.0 ) Tmax = Told;
			}

			period[frame] = Tmax;
			Told = Tmax;
		}
    }

    free( Window );
}

/*********************
*
*	FUNCTION: split_voiced
*
*	DESCRIPTION:
*		Function to split the very long voiced sections so 
*		that the find long	consistency routines do not have huge arrays to process.                    
*		The aim is to reduce the overall pitch extraction times and 
*		make them more  linear with relation to the file length.
*
***********************/
void split_voiced(INT32 *voiced, INT32 voiced_length, const FLOAT *pfVad, const UINT32 ulVadSz)
{
	FLOAT * min_vad ;				
	INT32 loop, long_loop ;		
	INT32 right_result = 0, left_result = 0 ; 
	INT32 tmp_search = 0 ;		
	INT32 num_long_voiced = 0;   
	INT32 * new_voiced ;			
	INT32 some_left = 0 ;		
	FLOAT tmp_min = 0 ;				
	INT32 passes ;				
	char hold=0;

	typedef struct					
	{
		INT32	left;				
		INT32	right;
	} long_voiced_struc;

	long_voiced_struc *long_voiced ; 

	min_vad = (FLOAT *)calloc(ulVadSz/8,sizeof(FLOAT)) ;

	for(loop = 0;loop < ((INT32)(ulVadSz/8) -1); loop++)
		min_vad[loop] = find_min(&pfVad[(loop * 8)], 16L) ;

	long_voiced = (long_voiced_struc *)calloc((voiced_length/2),sizeof(long_voiced_struc)) ;

	new_voiced = (INT32 *)calloc(voiced_length, sizeof(INT32)) ;

	for(loop = 2,hold=0; loop < (voiced_length-2); loop ++)
	{
		if( min_vad[loop] == 0 ) new_voiced[loop]=0;
		else if( (voiced[loop-2]) && (voiced[loop-1]) && (!voiced[loop]) && (voiced[loop+1]) && (voiced[loop+2]) )
		{
			new_voiced[loop] = 0;
			new_voiced[loop-1] = 0;	
			hold=1;
		}
		else
		{
			if (hold)
			{
				new_voiced[loop] = 0;
				hold=0;
			}
			else new_voiced[loop] = voiced[loop];
			if(voiced[loop]) some_left = 1 ;	
		}
	}

	new_voiced[0] = new_voiced[1] = new_voiced[(voiced_length-1)] = new_voiced[(voiced_length-2)] = 0 ;

	for(passes = 0; passes < 2; passes++)
	{
		num_long_voiced = 0 ;
		right_result = 0 ;
		left_result = 0 ;
		tmp_search = 0 ;

		do
		{
			tmp_search = find_value(1,&new_voiced[right_result],voiced_length - right_result) ;

			if (tmp_search!=-1)	
			{
				left_result = right_result + tmp_search ;

				right_result = left_result + find_value(0,&new_voiced[left_result],voiced_length - left_result) ;

				if((right_result - left_result) >= PITCH_LONG_VOICED_THRESH)
				{
					long_voiced[num_long_voiced].left = left_result ;
					long_voiced[num_long_voiced].right = right_result ;
					num_long_voiced++ ;
				}
			}
		} while(tmp_search != -1) ;

		for(loop=0;loop<num_long_voiced;loop++)
		{
			tmp_min = find_min(&min_vad[(long_voiced[loop].left) + 1], ((long_voiced[loop].right - long_voiced[loop].left) - 2));

			for(long_loop = long_voiced[loop].left; long_loop < long_voiced[loop].right; long_loop ++)
				if(min_vad[long_loop] < (1.1 * tmp_min)) new_voiced[long_loop] = 0 ;
		}
	}

	if(some_left == 1)	memcpy(voiced,new_voiced,(voiced_length * sizeof(INT32))) ;

	free(new_voiced) ;
	free(long_voiced) ;
	free(min_vad) ;
}

/*********************
*
*	FUNCTION: apitch_label
*
*	DESCRIPTION:
*		Process frames and label up pitch points with 1s (0s elsewhere).
*		Method used is cross-correlation of windowed speech segments with
*		an impulse train.  This is performed frame by frame to find the
*		optimum pitch mark locations.
*
***********************/
void apitch_label(const INT16 *psSpeech, INT32 * labels, INT32 Nsamples, INT32 * period, INT32 * voiced, INT32 Nframes ) 
{
    INT32 frame, offset, sample;
    FLOAT * Window = (FLOAT *) calloc( APITCH_FRAME_SIZE, sizeof(FLOAT));

    FLOAT X, Xmax;
    INT32 Tmax;

    Nframes = min( Nframes, 2 * (Nsamples / APITCH_FRAME_SIZE) - 1 );

    for( sample = 0; sample < Nsamples; sample++ ) labels[sample] = 0;

	for( sample = 0; sample < APITCH_FRAME_SIZE; sample++ )
		Window[sample] = (FLOAT)(0.5 * (1.0 - cos( (TWOPI * sample) / (FLOAT) APITCH_FRAME_SIZE)) );

    for( frame = 0; frame < Nframes; frame++ ) 
	{

        if( voiced[frame] > 0 ) 
		{
            Xmax = 0.0f;
            Tmax = 0;
            for( offset = 0; offset < period[frame]; offset++ ) {
                X = 0.0f;
                for( sample = offset; sample < APITCH_FRAME_SIZE; sample += period[frame] )
                    X += psSpeech[frame * APITCH_FRAME_SIZE_2 + sample] * Window[sample];
                if( X > Xmax ) {
                    Xmax = X;
                    Tmax = offset;
                }
            }

            for( sample = Tmax; sample < APITCH_FRAME_SIZE; sample += period[frame] )
                labels[frame * APITCH_FRAME_SIZE_2 + sample] = 1;
        }
    }

    free( Window);
}

/*********************
*
*	FUNCTION: inconsistency_handler
*
*	DESCRIPTION:
*		This function improves the initial pitch mack placement
*		by cleaning and refining their position until convergence.
*
***********************/
void inconsistency_handler(INT32 *pitch, INT32 *pitch_length, INT32 pitch_allocated, INT16 *speech, INT32 speech_length, INT16 average_pitch)
{

	INT32 loop, direction_loop;
	INT32 do_counter ;
	INT32 *inconsist ;
	INT32 inconsist_length = 0 ;
	INT32 first_last_length = 0 ;
	INT32 long_consist_length = 0 ;
	INT32 finished ;
	first_last_row_type *first_last ;
	long_consist_row_type *long_consist ;
	inconsist_row_type *inconsist_array ;

	local_search(pitch, *pitch_length, speech ) ;

	for(direction_loop = 0; direction_loop < 2; direction_loop++)
	{
		do_counter = 0 ;

		first_last = (first_last_row_type  *) calloc(*pitch_length, sizeof(first_last_row_type));
		long_consist = (long_consist_row_type  *) calloc(*pitch_length, sizeof(long_consist_row_type ));
		inconsist_array = (inconsist_row_type*) calloc(*pitch_length, sizeof(inconsist_row_type));

		inconsist = (INT32 *) calloc(*pitch_length, sizeof(INT32)) ;

		do
		{
			inconsist_length = 0 ;
			first_last_length = 0 ;
			long_consist_length = 0 ;
			remove_doubles(pitch, pitch_length) ;

			if (*pitch_length==0) break;

			generate_inconsist(pitch, *pitch_length, average_pitch, inconsist, &inconsist_length);

			find_long_consist(pitch, *pitch_length, inconsist, inconsist_array, inconsist_length, average_pitch,	
							  first_last, &first_last_length, 
							  long_consist, &long_consist_length) ;

			finished = move_pitch(first_last, &first_last_length, long_consist, 
								  long_consist_length, pitch, pitch_length, speech, pitch_allocated ) ;

			do_counter++ ;

		} while((do_counter != 500) && (!finished)) ;

		free(first_last) ;
		free(long_consist) ;
		free(inconsist_array) ; 
		free(inconsist) ;	
		reverse_array(speech, speech_length, SHORT_TYPE) ;

		for(loop = 0; loop<*pitch_length; loop++) pitch[loop] = speech_length - pitch[loop];

		reverse_array(pitch, *pitch_length, LONG_TYPE)  ;

	}
}

/*********************
*
*	FUNCTION: extend_pitch_controller
*
*	DESCRIPTION:
*		This function is a bundle of the boundaries reestimation
*		and the pitch marks cleaning.
*
***********************/
void extend_pitch_controller(INT16 *speech, INT32 speech_length, 
	INT32 *pitch, INT32 *pitch_length, INT32 length_threshold, INT32 pitch_allocated)
{
	INT32 loop, direction_loop ;

	for(direction_loop = 0; direction_loop < 2; direction_loop++)
	{
		extend_pitch(speech, speech_length, pitch, pitch_length, length_threshold, pitch_allocated)  ;

		reverse_array(speech, speech_length, SHORT_TYPE) ;

		for(loop = 0; loop<*pitch_length; loop++)
		{
			pitch[loop] = speech_length - pitch[loop] ;
		}
		reverse_array(pitch, *pitch_length, LONG_TYPE)  ;
	}

	remove_singles(pitch, pitch_length) ;
	remove_doubles(pitch, pitch_length) ; 
}

/*********************
*
*	FUNCTION: pos_or_neg
*
*	DESCRIPTION:
*		Combine positive and negative pitch arrays
*
*		Combines two pitch index arays to form a single array with best 
*		pitch mark estimates for each voicing section
*
***********************/
void pos_or_neg(INT32 average_pitch, INT16 *speech, 
				INT32 *pos_pitch, INT32 pos_pitch_length, INT32 *neg_pitch, 
				INT32 neg_pitch_length, INT32 *combined_pitch, 
				INT32 *combined_pitch_length)
{

INT32 loop=0 ;				  
INT32 *left_array ;			  
INT32 *right_array ;			  
INT32 first_last_length = 0 ;  
INT32 *tmp_pos ;				  
INT32 *tmp_neg ;				  
INT32 tmp_pos_length ;		  
INT32 tmp_neg_length ;		  
INT32 pos_start = 0 ;		  
INT32 neg_start = 0 ;		  
INT32 pos_end = -1 ;			  
INT32 neg_end = -1;			  
INT32 pos_loop ;				  
INT32 neg_loop ;				  
FLOAT pos_total ;				  
FLOAT neg_total ;				  

INT32 neg_cnt=0, pos_cnt=0;        
char pos_com=0, neg_com=0;

	if( (pos_pitch_length != 0) && (neg_pitch_length != 0) )
	{
		do
		{
			if ( !(pos_com || neg_com) )
			{
				if (pos_pitch[pos_cnt] < neg_pitch[neg_cnt])
				{
					combined_pitch[loop] = pos_pitch[pos_cnt];
					pos_cnt++;
					if (pos_cnt==pos_pitch_length) pos_com=1;
				}
				else if (pos_pitch[pos_cnt] > neg_pitch[neg_cnt])
				{
					combined_pitch[loop] = neg_pitch[neg_cnt];
					neg_cnt++;
					if (neg_cnt==neg_pitch_length) neg_com=1;
				}
				else 
				{
					combined_pitch[loop] = neg_pitch[neg_cnt];
					neg_cnt++;
					if (neg_cnt==neg_pitch_length) neg_com=1;
					pos_cnt++;
					if (pos_cnt==pos_pitch_length) pos_com=1;
				}
			}
			else
			{
				if (!pos_com) 
				{
					combined_pitch[loop] = pos_pitch[pos_cnt];
					pos_cnt++;
					if (pos_cnt==pos_pitch_length) pos_com=1;
				}
				else
				{
					combined_pitch[loop] = neg_pitch[neg_cnt];
					neg_cnt++;
					if (neg_cnt==neg_pitch_length) neg_com=1;
				}
			}
			loop++;

		} while (!(neg_com && pos_com) );
		*combined_pitch_length = loop;

		left_array = (INT32 *) calloc(*combined_pitch_length, sizeof(INT32)) ;
		right_array = (INT32 *) calloc(*combined_pitch_length, sizeof(INT32)) ;

		find_pos_neg_first_last(average_pitch, combined_pitch, *combined_pitch_length, 
								left_array, right_array, &first_last_length) ;

		tmp_pos = (INT32 *) calloc(pos_pitch_length, sizeof(INT32)) ;
		tmp_neg = (INT32 *) calloc(neg_pitch_length, sizeof(INT32)) ;

		*combined_pitch_length = 0 ;

		for(loop=0; loop<first_last_length; loop++)
		{

			tmp_pos_length = 0 ;
			tmp_neg_length = 0 ;

			neg_start = start_query(left_array[loop], neg_pitch, neg_end, neg_pitch_length) ;
			pos_start = start_query(left_array[loop], pos_pitch, pos_end, pos_pitch_length) ;
			neg_end = query(right_array[loop], neg_pitch, neg_start, neg_pitch_length) ;
			pos_end = query(right_array[loop], pos_pitch, pos_start, pos_pitch_length) ;

			pos_total = 0 ;
			neg_total = 0 ;

			if((pos_end - pos_start) != 0)
			{
				for(pos_loop = pos_start; pos_loop <= pos_end; pos_loop++)
				{
					pos_total += (FLOAT)speech[pos_pitch[pos_loop]] ;
				}

				if(pos_total != 0) pos_total = pos_total / (FLOAT)((pos_end + 1) - pos_start) ;
			}

			if((neg_end - neg_start) != 0)
			{
				for(neg_loop = neg_start; neg_loop <= neg_end; neg_loop++)
				{
					neg_total += -(FLOAT)speech[neg_pitch[neg_loop]] ;
				}
				if(neg_total != 0) neg_total = neg_total / (FLOAT)((neg_end + 1) - neg_start) ;
			}

			if(pos_total > neg_total)
			{

				if((pos_end - pos_start) != 0)
				{
					for(pos_loop = pos_start; pos_loop <= pos_end; pos_loop++)
					{
						combined_pitch[*combined_pitch_length]= pos_pitch[pos_loop] ;
						*combined_pitch_length += 1 ;
					}
				}
			}
			else
			{
				if((neg_end - neg_start) != 0)
				{
					for(neg_loop = neg_start; neg_loop <= neg_end; neg_loop++)
					{
						combined_pitch[*combined_pitch_length] = neg_pitch[neg_loop] ;
						*combined_pitch_length += 1 ;
					}
				}
			}
		}
		free(tmp_pos) ;
		free(tmp_neg) ;
		free(left_array);
		free(right_array) ;

		remove_singles(combined_pitch, combined_pitch_length) ;
	}
	else if (pos_pitch_length != 0)
	{
		memcpy(combined_pitch, pos_pitch, pos_pitch_length * sizeof(INT32));
	}
	else if (neg_pitch_length != 0)
	{
		memcpy(combined_pitch, neg_pitch, neg_pitch_length * sizeof(INT32));
	}
}

/*********************
*
*	FUNCTION: local_search
*
*	DESCRIPTION:
*		After initial pitch mark placements, which happen on a frame by frame basis, there can
*		be some improvement through a waveform based search to align the pitch candidates with 
*		the local waveform peaks. This function performs a search for the local waveform maximum
*		and then updates the position of the candidates.
*
***********************/
void local_search(INT32 *search_array, INT32 search_array_length, INT16 *speech_array)
{
	INT32 search_loop, fill_loop ;
	INT16 tmp_array[PITCH_SEARCH_DISTANCEx2] ;

	for(search_loop=0;search_loop<search_array_length;search_loop++)
	{

		if((search_array[search_loop] - PITCH_SEARCH_DISTANCE) >= 0)
		{
			for(fill_loop=0;fill_loop<(PITCH_SEARCH_DISTANCEx2);fill_loop++)
				tmp_array[fill_loop] = speech_array[((search_array[search_loop]) + fill_loop - PITCH_SEARCH_DISTANCE)];
			search_array[search_loop] =  find_max( tmp_array,0 ,(PITCH_SEARCH_DISTANCEx2), SHORT_TYPE) 
				+ search_array[search_loop] - PITCH_SEARCH_DISTANCE ;
		}
	}
}

/*********************
*
*	FUNCTION: remove_doubles
*
*	DESCRIPTION:
*		Find errored pitch marks corresponding to same pitch cycle.
*		It is based on a hard threshold which detects pitch marks which are too close together.
*
***********************/
void remove_doubles(INT32 *double_array, INT32 *double_array_length)
{
	INT32 loop ;					 
	INT32 *store_index = NULL;		
	INT32 store_index_Nelements=0;	

	store_index = (INT32*) malloc( (*double_array_length) * sizeof(INT32) );

	for (loop = 1;loop< (*double_array_length) ; loop++)
		if ( (double_array[loop]-double_array[loop-1]) < PITCH_DOUBLE_LIMIT )
		{
			store_index[store_index_Nelements] = loop-1;
			store_index_Nelements++;
		}

	if (store_index_Nelements) 
		multi_remove(store_index,store_index_Nelements, double_array,double_array_length);

	free(store_index);
}

/*********************
*
*	FUNCTION: generate_inconsist
*
*	DESCRIPTION:
*		The combination of overlapping frames leads to multiple pitch marks 
*		within close proximity due to the inaccuracies of the cross-correlation time alignment.
*		If two candidate pitch marks are separated by less than certain amount of samples 
*		they are considered to be derived from the same pitch peak and are replaced 
*		by a single candidate at the most likely position of the two. 
*
***********************/
void generate_inconsist(INT32 *current_pitch,INT32 current_pitch_length, 
							INT32 av_pitch, INT32 *inconsist, INT32 *inconsist_length)
{
	INT32 tmp_store = 1 ;
	FLOAT difference1, difference2 ;
	INT32 loop ;
	INT32 tmp_av ;
	INT32 current_length ;
	INT32 scaled_length_diff ;

	tmp_av = 3 * av_pitch ;

	for(loop=1; loop<current_pitch_length; loop++) 
	{
		current_length = current_pitch[loop] - current_pitch[loop - 1] ;
		scaled_length_diff = (100 * abs(current_length - tmp_store)) ;
		difference1 = (FLOAT) (scaled_length_diff / current_length);
		difference2 = (FLOAT) (scaled_length_diff / tmp_store) ;
		tmp_store = current_length ;

		if(((difference1 > PITCH_INCONSIST_LIMIT) || (difference2 > PITCH_INCONSIST_LIMIT)) || (tmp_store > tmp_av))
		{
			inconsist[*inconsist_length] = current_pitch[loop - 1] ;
			*inconsist_length += 1 ;		
		}	
	}
	inconsist[*inconsist_length] = current_pitch[current_pitch_length - 1] ;
	*inconsist_length += 1 ;
}

/*********************
*
*	FUNCTION: move_pitch
*
*	DESCRIPTION:
*		This function is used to move pitch marks which have caused inconsistencies in a
*		particular voiced section. Also checks to see if the run of consistent pitch marks has
*		been extended to the beginning of the voiced section and returns the correct value (1)
*		if this happy situation has been reached.
*
***********************/
INT32 move_pitch(first_last_row_type *first_last, INT32 *first_last_length, long_consist_row_type *long_consist,
				INT32 long_consist_length, INT32 *current_pitch, INT32 *current_pitch_length,
				INT16 *speech, INT32 pitch_allocated)
{
	INT32 loop, first_last_search_result, current_pitch_search_result, local_search_result, previous_mark, insertion_point ;
	INT32 finished = 1 ;
	INT32 tmp_index ;
	INT32 element_to_remove ;
	INT32 allocated ;

	allocated = *current_pitch_length ;

	for(loop=0;loop<long_consist_length;loop++)
	{
		first_last_search_result = bin_search_first_last(long_consist[loop].consist_left, first_last, *first_last_length);

		current_pitch_search_result = bin_search(long_consist[loop].consist_left, current_pitch, *current_pitch_length);

		if((first_last_search_result == -1) && (long_consist[loop].consist_left > PITCH_CLOSE_TO_START))
		{
			find_local_search_boundaries(current_pitch, current_pitch_search_result, speech, &previous_mark, &local_search_result);
			tmp_index = bin_search(long_consist[loop].consist_left, current_pitch, *current_pitch_length) ;
			element_to_remove = current_pitch[tmp_index - 1] ;

			if(bin_search(local_search_result, current_pitch, *current_pitch_length) != -1)
				insertion_point = bin_insert(previous_mark, current_pitch, current_pitch_length, &pitch_allocated) ;
			else
				insertion_point = bin_insert(local_search_result, current_pitch, current_pitch_length, &pitch_allocated) ;
			if (tmp_index != insertion_point)
				bin_remove(element_to_remove, current_pitch, current_pitch_length ) ;

			finished = 0 ;
		}
	}
	return(finished) ;
}

/*********************
*
*	FUNCTION: find_long_consist
*
*	DESCRIPTION:
*		This function searches for the boundaries of consistent pitched sections.
*		A consistent section is defined by a section where the instant pitch length
*		within the section is less than 3 time the average pitch length.
*
***********************/
void find_long_consist(INT32 *current_pitch, INT32 current_pitch_length, 
					   INT32 *inconsist, inconsist_row_type *inconsist_array, INT32 inconsist_length, 
					   INT16 average_pitch, first_last_row_type *first_last, INT32 *first_last_length, 
					   long_consist_row_type *long_consist, INT32 *long_consist_length)
{
	INT32 tmp_index=0, previous_pitch_mark = 0;
	INT32 *pitch_mark_gaps ;
	INT32 gaps_length = 0 ;
	INT32 last_time = -1 ;
	INT32 while_ctr = -1 ;
	INT32 max_pitch_array_index, max_inconsist_left, max_inconsist_right;
	INT16 loop ;
	for(loop=0;loop<inconsist_length;loop++)
	{
		tmp_index = tmp_index + (search(inconsist[loop], current_pitch + tmp_index, current_pitch_length - tmp_index)) ;

		inconsist_array[loop].pitch_array_index = tmp_index - (last_time + 1);
		last_time = tmp_index ;	

		if(loop==0)	inconsist_array[loop].inconsist_left = 0;
		else inconsist_array[loop].inconsist_left = inconsist[loop-1];

		inconsist_array[loop].inconsist_right = inconsist[loop];
	}

	pitch_mark_gaps = (INT32 *) calloc(2 * current_pitch_length, sizeof(INT32));

	if( pitch_mark_gaps != NULL ) gaps_length = 0;
	else free(pitch_mark_gaps);

	for (loop=0;loop<current_pitch_length;loop++)
	{
		if((current_pitch[loop]-previous_pitch_mark) > (average_pitch * PITCH_MULTIPLIER))
		{	
			pitch_mark_gaps[gaps_length] = previous_pitch_mark ;
			pitch_mark_gaps[gaps_length + 1] = current_pitch[loop] ;
			gaps_length += 2 ;
		}
		previous_pitch_mark = current_pitch[loop] ;
	}
	pitch_mark_gaps[gaps_length] = current_pitch[loop-1] ;
	gaps_length += 1 ;

	for (loop=0;loop<((gaps_length-1)/2);loop++)
	{
		if(pitch_mark_gaps[(((loop+1)*2)-1)] != pitch_mark_gaps[((loop+1)*2)])
			{
			first_last[*first_last_length].segment_left = pitch_mark_gaps[(((loop+1)*2)-1)] ;
			first_last[*first_last_length].segment_right = pitch_mark_gaps[((loop+1)*2)] ;
			*first_last_length += 1 ;
		}
	}

	for(loop=0;loop<*first_last_length;loop++)
	{
		max_pitch_array_index = 0 ;
		max_inconsist_left = 0 ;
		max_inconsist_right = 0 ;

		do 
		{
			while_ctr += 1 ;
			if (inconsist_array[while_ctr].pitch_array_index >= max_pitch_array_index)
			{
				max_pitch_array_index = inconsist_array[while_ctr].pitch_array_index ;
				max_inconsist_left = inconsist_array[while_ctr].inconsist_left ;
				max_inconsist_right = inconsist_array[while_ctr].inconsist_right ;
			}
		}while(first_last[loop].segment_right != inconsist_array[while_ctr].inconsist_right) ;
		long_consist[*long_consist_length].consist_left = max_inconsist_left ;
		long_consist[*long_consist_length].consist_right = max_inconsist_right ;
		*long_consist_length +=1 ;
	}

	free(pitch_mark_gaps) ;
}

/*********************
*
*	FUNCTION: extend_pitch
*
*	DESCRIPTION:
*		This function refine the boundaries of the pitched sections
*		according to a pitch length and energy test.
*
***********************/
void extend_pitch(INT16 *speech, INT32 speech_length, INT32 *tmp_pitch, INT32 *tmp_pitch_length, INT32 length_threshold, INT32 pitch_allocated)
{
	INT32 first_mark ;		
	INT32 candidate ;		
	INT32 previous ;			
	INT32 length ;			
	INT32 left_loop ;		
	INT32 *first ;			
	INT32 first_length;		
	INT32 *stable_test ;		
	INT32 stable_length = 0;	
	INT32 *new_insertions ;	
	INT32 new_insertions_length = 0; 
	INT32 *pitch ;			
	INT32 *first_pitch ;		
	INT32 pitch_length ;		
	INT32 first_pitch_length ; 
	INT16 *current_cycle ;	
	INT16 *previous_cycle ;	
	FLOAT cross_power_result ; 
	FLOAT rms_previous ;		 
	FLOAT rms_ratio ;			
	INT32 beginning_result ;	
	INT32 comparison_result = 0 ; 
	INT32 pitch_same = 1 ;	
	INT32 test_count ;		
	INT32 tmp_loop = 0 ;
	INT32 stable_allocated = 0 ;
	INT32 insertions_allocated = 0 ;

	insertions_allocated = (speech_length / 4) ;
	new_insertions = (INT32  *) calloc(insertions_allocated, sizeof(INT32));

	stable_allocated = *tmp_pitch_length ;
	stable_test = (INT32  *) calloc(stable_allocated, sizeof(INT32) );
	pitch = (INT32  *) calloc(speech_length / 4, sizeof(INT32));

	first = (INT32  *) calloc(speech_length/10, sizeof(INT32));
	first_pitch = (INT32  *) calloc(speech_length / 4, sizeof(INT32)) ;
	tmp_loop = 0 ;

	do
	{

		pitch_length = *tmp_pitch_length ;
		memcpy(pitch, tmp_pitch, (sizeof(INT32) * (*tmp_pitch_length))) ;
		remove_singles(tmp_pitch, tmp_pitch_length) ;
		first_length = 0 ;
		find_first_only(tmp_pitch, *tmp_pitch_length, first, &first_length) ;
		first_pitch_length = *tmp_pitch_length ;

		memcpy(first_pitch, tmp_pitch, (sizeof(INT32) * (*tmp_pitch_length))) ;

		for(left_loop = 0;left_loop<first_length; left_loop++) 
		{

			first_mark = first_pitch[first[left_loop]] ;
			length = first_pitch[(first[left_loop] + 1)] - first_mark ;

			if(first[left_loop] == 0) previous = -1 ;
			else previous = first_pitch[(first[left_loop] - 1)] ;

			candidate = first_mark - length ;

			if( (candidate < 96) || (candidate > speech_length-96))
			{
				if(bin_search(first_mark, new_insertions, new_insertions_length) != -1)
					bin_insert(first_mark, stable_test, &stable_length, &stable_allocated) ;
				else 
					bin_remove(first_mark, tmp_pitch, tmp_pitch_length) ;
			}
			else if((bin_search(first_mark, stable_test, stable_length)) == -1)
			{
				current_cycle = (INT16  *) calloc(length, sizeof(INT16));
				previous_cycle = (INT16  *) calloc(length, sizeof(INT16));
				split_arrays(&candidate, length, speech, previous_cycle, current_cycle) ;
				rms_previous = rms(previous_cycle, length, SHORT_TYPE) ;
				rms_ratio = (FLOAT) ((rms(current_cycle, length, SHORT_TYPE) + 0.01) / (rms_previous + 0.01));
				cross_power_result = cross_power_ratio(previous_cycle, length, current_cycle, length) ;
				free(current_cycle);
				free(previous_cycle);

				if( (rms_previous > PITCH_ENERGY_THRESHOLD) && (rms_ratio < PITCH_RMS_THRESHOLD) && 
					(cross_power_result > PITCH_CROSS_POWER_THRESHOLD) && (length > length_threshold) )
				{
					comparison_result = 1 ;
				}
				else comparison_result = 0 ;

				beginning_result = beginning_test(candidate, previous) ;

				if(beginning_result)
				{
					if(!comparison_result)
					{
						if(bin_search(first_mark, new_insertions, new_insertions_length) != -1)
							bin_insert(first_mark, stable_test, &stable_length, &stable_allocated) ;
						else 
							bin_remove(first_mark, tmp_pitch, tmp_pitch_length) ;
					}
				}
				else
				{
					if(comparison_result)
					{
						bin_insert(candidate, new_insertions, &new_insertions_length, &insertions_allocated) ;
						bin_insert(candidate, tmp_pitch, tmp_pitch_length, &pitch_allocated) ;
					}
					else
					{
						if(bin_search(first_mark, new_insertions, new_insertions_length) != -1)
							bin_insert(first_mark, stable_test, &stable_length, &stable_allocated) ;
						else 
							bin_remove(first_mark, tmp_pitch, tmp_pitch_length) ;
					}
				}
			}
		}

		pitch_same = 1 ;

		if(*tmp_pitch_length != pitch_length) pitch_same = 0 ;
		else
		{
			test_count = 0 ;
			while((pitch_same == 1) && (test_count < pitch_length))
			{
				if (tmp_pitch[test_count] != pitch[test_count]) pitch_same = 0 ;
				test_count += 1 ;
			}
		}

		tmp_loop += 1 ;

	}
	while(pitch_same == 0) ;

	free(pitch) ;
	free(first) ;
	free(first_pitch) ;
	free(new_insertions);
	free(stable_test);
}

/*********************
*
*	FUNCTION: remove_singles
*
*	DESCRIPTION:
*		This function remove single pitch marks.
*		A pitch mark is defined as single when the ratio 
*		of the distance from the considered ptich mark
*		to the previous and next pitch mark exceed a 
*		hard threshold.
*
***********************/
void remove_singles(INT32 *pitch, INT32 *pitch_length)
{
	INT32  last;                
	INT32  pos;                 
	INT32  *diff;               
	INT32  diff_Nelements = 0;  
	INT32  *epitch;             
	INT32  epitch_Nelements = 0;
	FLOAT ratio = 0;           
	INT32  count;               

	epitch_Nelements = *pitch_length + 2;
	epitch = (INT32 *) malloc(epitch_Nelements * sizeof(INT32));

	for (count=0;count<*pitch_length;count++) epitch[count+1] = pitch[count];

	epitch[0] = -500;	
	epitch[epitch_Nelements-1] = epitch[epitch_Nelements-2]+500;

	diff_Nelements = epitch_Nelements-1;
	diff = (INT32 *) malloc(diff_Nelements * sizeof(INT32));
	for (count=0;count<diff_Nelements;count++) diff[count] = epitch[count+1] - epitch[count];

	for (count=0, last=0, pos=0; count<diff_Nelements; count++)
	{
		if ( (last < 150) || (diff[count]< 150)) 
		{
			pos++;
			epitch[pos-1] = epitch[count];
		}
		last = diff[count];
	}
	epitch[pos] = epitch[epitch_Nelements-1];
	epitch_Nelements = pos+1;

	diff_Nelements = epitch_Nelements-1;
	for (count=0;count<diff_Nelements;count++) diff[count] = epitch[count+1] - epitch[count];

	diff_Nelements--;
	for (count=0;count<diff_Nelements;count++) 
	{
		if (diff[count]>diff[count+1]) ratio = ((FLOAT) diff[count] / (FLOAT) diff[count+1] );
		else ratio = ((FLOAT) diff[count+1] / (FLOAT) diff[count] );
		if (ratio < PITCH_MAX_CONSECUTIVE_PITCH_LENGTH_RATIO) diff[count] = 1;
		else diff[count] = 0;
	}

	for (count=1,last=0; count<diff_Nelements;count++)
	{
		if (diff[count]!=last)
		{
			if (diff[count]==0) 
			{
				diff[count]=1;
				last=0;
			}
			else if (diff[count]==1)
			{
				diff[count-1]=1;
				last=1;
			}
		}
		else last=diff[count];
	}

	for (count=0,last=0;count<diff_Nelements;count++)
	{
		pitch[last] = epitch[count+1];
		if (diff[count]) last++;
	}
	*pitch_length = last;

	free(epitch);
	free(diff);
}

/*********************
*
*	FUNCTION: find_pos_neg_first_last
*
*	DESCRIPTION:
*		Finds the first pitch mark and the last pitch mark in each voicing section.
*
***********************/
void find_pos_neg_first_last(INT32 average_pitch, INT32 *current_pitch, 
							 INT32 current_pitch_length, INT32 *first_last_left_array, 
							 INT32 *first_last_right_array, INT32 *first_last_length)
{
	typedef struct  
	{ 
		INT32 segment_left;  
		INT32 segment_right; 
	}first_last_row_type ;      

	INT32 loop ;                    
	INT32 previous_pitch_mark = 0 ; 
	INT32 *pitch_mark_gaps ;        
	INT32 gaps_length = 0 ;

	first_last_row_type *first_last ; 

	first_last = (first_last_row_type  *) calloc(current_pitch_length,sizeof(first_last_row_type));
	if( first_last != NULL )
	{

	}
	else free(first_last); 

	pitch_mark_gaps = (INT32 *) calloc(2 * current_pitch_length, sizeof(INT32));

	if( pitch_mark_gaps != NULL ) gaps_length = 0;
	else free(pitch_mark_gaps); 

	for (loop=0;loop<current_pitch_length;loop++)
	{
		if(((current_pitch[loop]-previous_pitch_mark) > (average_pitch * PITCH_TEST)) || (loop == 0))
		{	
			pitch_mark_gaps[gaps_length] = previous_pitch_mark ;
			pitch_mark_gaps[gaps_length + 1] = current_pitch[loop] ;
			gaps_length += 2 ;

		}
		previous_pitch_mark = current_pitch[loop] ;
	}
	pitch_mark_gaps[gaps_length] = current_pitch[loop-1] ;
	gaps_length += 1 ;

	for (loop=0;loop<((gaps_length-1)/2);loop++)
	{
		if(pitch_mark_gaps[(((loop+1)*2)-1)] != pitch_mark_gaps[((loop+1)*2)])
		{
			first_last[*first_last_length].segment_left = pitch_mark_gaps[(((loop+1)*2)-1)] ;
			first_last[*first_last_length].segment_right = pitch_mark_gaps[((loop+1)*2)] ;
			*first_last_length += 1 ;
		}
	}

	for(loop=0;loop<*first_last_length;loop++)
	{
		first_last_left_array[loop] = first_last[loop].segment_left ;
		first_last_right_array[loop] = first_last[loop].segment_right ;
	}
	free(pitch_mark_gaps) ;
	free(first_last) ;
}

/*********************
*
*	FUNCTION: start_query
*
*	DESCRIPTION:
*		Find the position of start index within the input array
*
***********************/
INT32 start_query(INT32 value, INT32 *array, INT32 start, INT32 array_length)
{
	INT32 loop ;				
	INT32 result = start;	

	for(loop = array_length; loop > start; loop --)
		if(array[loop] >= value) result = loop ;

	return result ;
}

/*********************
*
*	FUNCTION: query
*
*	DESCRIPTION:
*		Find the position of a particular pitch index within an existing array
*
***********************/
INT32 query(INT32 value, INT32 *array, INT32 start, INT32 array_length)
{
	INT32 loop ;				
	INT32 result = start;	

	for(loop = start; loop < array_length; loop ++)
		if(array[loop] <= value) result = loop ;

	return result ;
}

/*********************
*
*	FUNCTION: find_local_search_boundaries
*
*	DESCRIPTION:
*		Determines the boundaries for this search, and then 
*		updates the local_search_result to contain the position 
*		of the most likely candidate within this area.
*
***********************/
void find_local_search_boundaries(INT32 *current_pitch_array, INT32 current_index, INT16 *speech_array, 
								  INT32 *previous_mark_position, INT32 *local_search_result)
{
	INT32 cycle_length;
	INT32 local_search_start;
	INT32 search_length;

	cycle_length = current_pitch_array[current_index +1] - current_pitch_array[current_index] ;
	*previous_mark_position = current_pitch_array[current_index] - cycle_length;
	if(floor(PITCH_LENGTH_CHANGE*cycle_length) < PITCH_LOCAL_SEARCH_MINIMUM) 
		local_search_start = *previous_mark_position - (INT32)(floor(PITCH_LENGTH_CHANGE * ((FLOAT)cycle_length))) ;
	else local_search_start = *previous_mark_position - PITCH_LOCAL_SEARCH_MINIMUM ;	

	search_length = 2 * (*previous_mark_position - local_search_start) ;

	*local_search_result = find_max(speech_array, local_search_start, search_length,SHORT_TYPE) ;
}

/*********************
*
*	FUNCTION: split_arrays
*
*	DESCRIPTION:
*
*
***********************/
void split_arrays(INT32 *new_mark, INT32 length, INT16 *speech, INT16 *previous_cycle, INT16 *current_cycle)
{

	INT32 loop ;
	for(loop = 0;loop<length; loop++)
	{
		previous_cycle[loop] = speech[*new_mark + loop] ;
		current_cycle[loop] = speech[*new_mark + length + loop] ;
	}

	*new_mark = find_max(speech, (*new_mark - 5), 10, SHORT_TYPE ) ;
}

/*********************
*
*	FUNCTION: beginning_test
*
*	DESCRIPTION:
*
*
***********************/
INT32 beginning_test(INT32 candidate_mark, INT32 previous_mark) 
{
	if( ( (previous_mark - (candidate_mark - 25)) > 50) || (previous_mark > (candidate_mark - 25)) 
		|| (candidate_mark <= 0) )
	return(1) ;
	else return(0) ;
}

/*********************
*
*	FUNCTION: First_Last
*
*	DESCRIPTION:
*		Extract the pitch sections boundaries.
*
*
***********************/
void First_Last(INT32 *pitch, INT32 *pitch_Nelements, INT32 *first_last[2],INT32 *fl_Nelements)

{
	INT32 i,j;
	INT32 *temp;

	*fl_Nelements = 0;
	temp = (INT32 *) malloc ((*pitch_Nelements) * sizeof(INT32));
	temp[0] = 0;

	for (i=1, j=0; i<*pitch_Nelements; i++)
	{
		if ( ((FLOAT)((pitch[i+1] - pitch[i])/(pitch[i] - pitch[i-1]))) > 1.5 ) 
		{
			temp[++j] = i;
			temp[++j] = i+1;
		}
	}
	temp[++j] = *pitch_Nelements-1;

	*fl_Nelements = j+1; 

	for (i=0, j=0; i<*fl_Nelements; i+=2)
	{
		if ( (temp[i+1]-temp[i])<2 ) j+=2;
		else 
		{
			temp[i-j] = temp[i]; 
			temp[i-j+1] = temp[i+1];
		}
	}
	*fl_Nelements -=j;
	*fl_Nelements /=2;

	first_last[0] = (INT32 *) malloc (*fl_Nelements * sizeof(INT32));
	first_last[1] = (INT32 *) malloc (*fl_Nelements * sizeof(INT32));

	for (i=0,j=0; i<*fl_Nelements;i++,j+=2)
	{
		first_last[0][i] = temp[j];
		first_last[1][i] = temp[j+1];
	}

	free(temp);
}

/*********************
*
*	FUNCTION: find_first_only
*
*	DESCRIPTION:
*
*
***********************/
void find_first_only(INT32 *pitch, INT32 pitch_length, INT32 *first, INT32 *first_length)
{
	INT32 loop ;
	INT32 last_difference = 1;
	INT32 difference ;

	for(loop = 0; loop < pitch_length; loop++)
	{
		if(loop == 0) difference = pitch[0] ;
		else difference = pitch[loop] - pitch[loop-1];
		if((difference/last_difference) > PITCH_DIFFERENCE_THRESHOLD)
		{
			first[*first_length] = loop ;
			*first_length += 1 ;			
		}
		last_difference = difference ;
	}
}

/*********************
*
*	FUNCTION: bin_search_first_last
*
*	DESCRIPTION:
*		Embedded of bin_query_first_last() function.
*		This function searches in segment_left of an array of first_last_row_type structure
*		for a specific value.
*
***********************/
INT32 bin_search_first_last( INT32 value, first_last_row_type * list_array, INT32 Nelements ) 
{
    INT32 ptr = bin_query_first_last( value, list_array, Nelements );
    if( ptr == -1L ) return -1L;
    if( list_array[ptr].segment_left == value ) return ptr;
    return -1L;
}

/*********************
*
*	FUNCTION: bin_query_first_last
*
*	DESCRIPTION:
*		Search for a specific value of segment_left in an array of first_last_row_type structure.
*
***********************/
INT32 bin_query_first_last( INT32 value, first_last_row_type * list_array, INT32 Nelements ) 
{
    INT32 step_size;
    INT32 ptr;

    if( Nelements <= 0L ) return -1L;
    if( list_array == NULL ) return -1L;
    if( Nelements == 1L ) return 0L;

    step_size = 1L;
    while( step_size < Nelements ) step_size *= 2L;
    step_size /= 2L;

    ptr = step_size;
    while( step_size > 0L ) {
        if( list_array[ptr].segment_left > value ) {
            if( (ptr - step_size) >= 0L ) ptr -= step_size;

        } else if( list_array[ptr].segment_left < value ) {
            if( (ptr + step_size) < Nelements ) ptr += step_size;
        } else
		{
            return ptr;
		}
        step_size /= 2L;
    }
    return ptr;
}
