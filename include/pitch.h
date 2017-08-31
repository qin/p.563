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



#ifndef PITCH_H
    #define PITCH_H

	#define APITCH_FRAME_SIZE 512
	#define APITCH_FRAME_SIZE_PLUS2 514
	#define APITCH_FRAME_SIZE_2 256
	#define APITCH_ROLL 32
	#define PITCH_LONG_VOICED_THRESH 15 	
	#define PITCH_MED 63				
	#define PITCH_NOM 32				
	#define PITCH_SEARCH_DISTANCE 5		
	#define PITCH_SEARCH_DISTANCEx2	10	
	#define PITCH_DOUBLE_LIMIT 8		
	#define PITCH_INCONSIST_LIMIT 15
	#define PITCH_CLOSE_TO_START 200
	#define PITCH_MULTIPLIER 3
	#define PITCH_RMS_THRESHOLD 3
	#define PITCH_ENERGY_THRESHOLD 500 
	#define PITCH_CROSS_POWER_THRESHOLD 0.1
	#define PITCH_MAX_CONSECUTIVE_PITCH_LENGTH_RATIO 1.5
	#define PITCH_DIFFERENCE_THRESHOLD 1.5
	#define PITCH_TEST 2				
	#define PITCH_LENGTH_CHANGE 0.15
	#define PITCH_LOCAL_SEARCH_MINIMUM 5

	#define VTP_ORDER 9
	#define VTP_NB_TUBES 8

	typedef struct  
		{ 
			INT32 segment_left;
			INT32 segment_right;
		} first_last_row_type;

	typedef struct  
		{ 
			INT32 consist_left;
			INT32 consist_right;
		} long_consist_row_type;

	typedef struct  
		{ 
			INT32 pitch_array_index;
			INT32 inconsist_left;
			INT32 inconsist_right;
		} inconsist_row_type ;


	void pitch_extract(const INT16 *psNormSpeech, const UINT32 ulNormSpeechSz,
						   const FLOAT *pfVad, const UINT32 ulVadSz,
						   INT32 *pitch, INT32 *pitch_length, INT32 pitch_allocated, INT16 *average_pitch);

	void apitch_run( const INT16 *psSpeech, const UINT32 ulSpeechSz, INT32 *period, INT32 *voiced, INT32 Nframes ) ;

	void split_voiced(INT32 *voiced, INT32 voiced_length, const FLOAT *pfVad, const UINT32 ulVadSz);

	void apitch_label(const INT16 *psSpeech, INT32 * labels, INT32 Nsamples, INT32 * period, INT32 * voiced, INT32 Nframes ) ;

	void inconsistency_handler(INT32 *pitch, INT32 *pitch_length, INT32 pitch_allocated, INT16 *speech, INT32 speech_length, INT16 average_pitch);

	void extend_pitch_controller(INT16 *speech, INT32 speech_length, 
		INT32 *pitch, INT32 *pitch_length, INT32 length_threshold, INT32 pitch_allocated);

	void pos_or_neg(INT32 average_pitch, INT16 *speech, 
					INT32 *pos_pitch, INT32 pos_pitch_length, INT32 *neg_pitch, 
					INT32 neg_pitch_length, INT32 *combined_pitch, 
					INT32 *combined_pitch_length);

	void local_search(INT32 *search_array, INT32 search_array_length, INT16 *speech_array);

	void remove_doubles(INT32 *double_array, INT32 *double_array_length);

	void generate_inconsist(INT32 *current_pitch,INT32 current_pitch_length, 
								INT32 av_pitch, INT32 *inconsist, INT32 *inconsist_length);

	INT32 move_pitch(first_last_row_type *first_last, INT32 *first_last_length, long_consist_row_type *long_consist,
					INT32 long_consist_length, INT32 *current_pitch, INT32 *current_pitch_length,
					INT16 *speech, INT32 pitch_allocated);

	void find_long_consist(INT32 *current_pitch, INT32 current_pitch_length, 
						   INT32 *inconsist, inconsist_row_type *inconsist_array, INT32 inconsist_length, 
						   INT16 average_pitch, first_last_row_type *first_last, INT32 *first_last_length, 
						   long_consist_row_type *long_consist, INT32 *long_consist_length);
	void extend_pitch(INT16 *speech, INT32 speech_length, 
						  INT32 *tmp_pitch, INT32 *tmp_pitch_length, INT32 length_threshold, 
						  INT32 pitch_allocated);
	void remove_singles(INT32 *pitch, INT32 *pitch_length);

	void find_pos_neg_first_last(INT32 average_pitch, INT32 *current_pitch, 
								 INT32 current_pitch_length, INT32 *first_last_left_array, 
								 INT32 *first_last_right_array, INT32 *first_last_length);
	INT32 start_query(INT32 value, INT32 *array, INT32 start, INT32 array_length);

	INT32 query(INT32 value, INT32 *array, INT32 start, INT32 array_length);
	
	void find_local_search_boundaries(INT32 *current_pitch_array, INT32 current_index, INT16 *speech_array, 
									 INT32 *previous_mark_position, INT32 *local_search_result);
	
	void find_first_only(INT32 *pitch, INT32 pitch_length, INT32 *first, INT32 *first_length);
	
	void split_arrays(INT32 *new_mark, INT32 length, INT16 *speech, INT16 *previous_cycle, INT16 *current_cycle);
	
	INT32 beginning_test(INT32 candidate_mark, INT32 previous_mark);
	
	void First_Last(INT32 *pitch, INT32 *pitch_Nelements, INT32 *first_last[2],INT32 *fl_Nelements);
	
	INT32 bin_search_first_last( INT32 value, first_last_row_type * list_array, INT32 Nelements ) ;

	INT32 bin_query_first_last( INT32 value, first_last_row_type * list_array, INT32 Nelements );

#endif
