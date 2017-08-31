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


#ifndef __INTERR_DETECT_H__
#define __INTERR_DETECT_H__

/* specific prototypes */ 
#include <math.h>
#include "generic_typedefs.h"
#include "hosm.h"

/* Local definitions */ 
#define INTERR_FOUND		301
#define INTERR_NOT_FOUND	302

#define FRAME_TIME			32E-3
#define MAX_NOISE_DB		-55			/* maximum noise level in dBov */

#define FRAME_MOVE		1
#define NO_FRAME_MOVE	0



/* Local Structures */
typedef	struct {
	INT32	inter_number;				/* Number of actual interruption */
	INT32	end_frame;					/* Show when enough values are available to test  */
										/* packet-move */
	INT32	end_pos;					/* End position of interruption */
	FLOAT	*samples_middle;			/* Samples middle of packet-move test */
	FLOAT	*ceps_left;					/* Calculated cepstrum to left of interruption */
}
typpacketswitch;



typedef struct {

	INT32	speech;						/* Variable for speech recognition */
	INT32	num_samples;				/* Number of samples in 32 ms */
	INT32	no_interruption[3];			/* Interruption state (0 interruption, 1 no interruption) */
	INT32	last_position[3];			/* Last sample that has been used in average frame A	 */
	INT32	num_samples_frame;			/* Number of samples in fundamental frame */
	INT32	num_int;					/* Number of interruptions in actual buffer */
	INT32	num_int_all;				/* Number of all interruptions */
	INT32	open_packet_cnt;			/* Counter for not checked packet-switching */
	INT32	max_open;					/* Maximum actual number of not tested interruptions */
	INT32	start_position;				/* Start position of interruption before it is  */
										/* written into interruption buffer */
	INT32	end_position;				/* End position of interruption before it is written */
										/* into interruption buffer */
	INT32	last_pos;					/* End position to calculate power 1ms after start of  */
										/* interruption */
	INT32	local_end;					/* End position of 5ms frames during interruption */
	INT32	frame_offset[3];			/* Data offset in dataBuffer frame  */
	FLOAT	energy_rate;				/* Short-time-energy value during interruption */
	FLOAT	est_rms;					/* Estimated RMS value to the left and to the right side  */
										/* of interruption */
	FLOAT	value_level;				/* Maximum of noise calculated into real value. Used  */
										/* for threshold. */
	FLOAT	pow_af_st;					/* Power 1ms after interruption start */
	FLOAT	DCint;						/* Sum of values during interruption for DC removal */
	FLOAT	DCOffset[12];				/* DC-Offset value of 16ms buffers */
	FLOAT	local_pow[2];				/* Power of 5ms frames during interruption */
	FLOAT	*buffer;					/* Data buffer with 3*16ms datas	 */
	FLOAT	*packetmem;					/* Data buffer with datas for packet-switching test */
	typpacketswitch strpacket[19];		/* Data for packet switching test */
}
typInterr_detect;



/*/////////////////////////////////////////////////////////////////////////////////////////////// */
/* Function:	detectInterruptions	 */
/*----------------------------------------------------------------------------------------------- */
/* Description:	Detects interruptions in given audio signal. Signal can be 8 kHz or 16 kHz. */
/*				 */
/* Authors:		 */
/* */
/* Log:			27/07/01 */
/*				27/07/01 */
/*  */
/* Input:		*par			Pointer to input parameter (example: sampling frequency) */
/*				*channel		Pointer to channel parameter (example: datas) */
/*				*interruption	Pointer to interruption parameter (example: interruption start) */
/*				*int_det		Pointer to interruption detection structure (see above) */
/*				*ImpulseNoise	Pointer to impulse noise positions */
/*				 */
/* Output:		*interruption	Interruption start, length and type of interruption (Frame move,  */
/*								no frame move) */
/* */
/* Return:		INTERR_FOUND		 */
/*				INTERR_NOT_FOUND */
/*				If an error appears return is < 0. */
/* */
/*/////////////////////////////////////////////////////////////////////////////////////////////// */
INT32 detectInterruptions(typInputParameter *par, typChannel *channel, 
						  typInterruption *interruption, typInterr_detect *int_det);



/*/////////////////////////////////////////////////////////////////////////////////////////////// */
/* Function:	initDetectInterruptions	 */
/*----------------------------------------------------------------------------------------------- */
/* Description:	Initialize structure of typ typInterr_detect */
/*				 */
/* Authors:		 */
/* */
/* Log:			02/10/01 */
/*				02/10/01 */
/*  */
/* Input:		*int_det	Pointer to structure */
/*				*par		Pointer to structure to get sampling frequency */
/*				 */
/* Output:		*int_det	Pointer to initialized structure */
/* */
/* Return:		If an error appears return is < 0. */
/* */
/*/////////////////////////////////////////////////////////////////////////////////////////////// */
INT32 initDetectInterruptions(typInterr_detect *int_det, typInputParameter *par);




/*/////////////////////////////////////////////////////////////////////////////////////////////// */
/* Function:	DeinitDetectInterruptions	 */
/*----------------------------------------------------------------------------------------------- */
/* Description:	Free memory of initialized structure. */
/*				 */
/* Authors:		 */
/* */
/* Log:			06/11/01 */
/*				06/11/01 */
/*  */
/* Input:		*int_det Pointer to structure */
/*				 */
/* Output:		*int_det Pointer to deinitialized structure */
/* */
/* Return:		none */
/* */
/*/////////////////////////////////////////////////////////////////////////////////////////////// */
void DeinitDetectInterruptions(typInterr_detect *int_det);

#endif /*__INTERR_DETECT_H__ */
