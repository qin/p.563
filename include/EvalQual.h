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



#ifndef __EVALQUAL_H__
#define __EVALQUAL_H__

#define  _PCM_SAMPLES
#define  NUMBER_OF_BYTES_PER_SAMPLE     2
  

/* String length */
#define     STRLEN                          300


#define     FALSE                           0
#define     TRUE                            1

/* Define 8 kHz sampling frequency and a 32 ms local FFT window for the disturbance processing */
#define     SAMPLE_FREQUENCY_HZ             8000
#define     TRANSFORM_LENGTH                256
#define     LOG2_TRANSFORM_LENGTH           8

#define     ONE_KILO_HZ                     1000

/* Set the subsampling factor for the envelope series */
#define    SUBSAMPLING_FACTOR_ENVELOPES     4


/* Parameters to delineate active speech part in file */
#define     CRITERIUM_FOR_SILENCE_OF_SAMPLES        100.
#define     CRIT_NUM                                5

/* For scaling */
#define     TARGET_AVG_POWER                1E7


		
FLOAT ComputeBasicVoiceQual(FLOAT * fpDistortedSignal,
								    FLOAT * fpEnhancedSignal,
									 INT32   iNrSamples,
									 INT32	iSampSkipBegin,
									 INT32	iSampSkipEnd,
									 FLOAT	fractionActive,
									 FLOAT	*pfABD,
									 FLOAT	*pfACRMOS);




#define POWER_DISTURBANCE_BARK_SYM      0.8f  
#define POWER_SPLIT_SECOND_SYM          0.7f  
#define POWER_BLOCK_SYM                 0.5f  

#define POWER_DISTURBANCE_BARK_ASYM     1.0f  
#define POWER_SPLIT_SECOND_ASYM         0.6f  
#define POWER_BLOCK_ASYM                0.8f	 

	
#define SBDFACTOR		-0.99f
#define ADBFACTOR		1.0F
#define FLFACTOR		-0.021f
#define ACRMOSOFFSET -3.0f
#define ACRMOSSCALE  4.0f




#endif /* __EVALQUAL_H__ */
