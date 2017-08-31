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


#ifndef __LPCANALYSIS_H_
#define __LPCANALYSIS_H_

/* LPC Analysis buffer */
struct _LPCBUFFER_DATA;
typedef struct _LPCBUFFER_DATA *LPCBUFFER_HANDLE;

typedef enum
{
	LPCBUFFER_OK,
	LPCBUFFER_MEMORY_ALLOCATION_FAILED,
	LPCBUFFER_SIZEERROR,
	LPCBUFFER_NOTINITIALIZED

}LPCBUFFER_ERRORCODE;


typedef struct _LPCBUFFER_DATA
{

	INT32 lStepsize;		/* */
	INT32 lLPCWindowSize;   /* for LPC analysis */
	INT32 lBufferDelay;	/* delay */
	INT32 lLPCOrder;		/* have to be 10 in this implementation */
	/* local buffers */
	FLOAT *pdLPCTimeBuffer;
	FLOAT *pdLPCTimeWindow;		/* const: window for lpc analysis */
	FLOAT *pdLPCWindTimeSignal; /* windowd time singnal - for LPC analysis*/


	FLOAT *pdAnalysisSignal; /* pointer to actual time signal */


}LPCBUFFER_DATA;


void LPCBufferDelete(LPCBUFFER_HANDLE *phBuff);
void	UpdateLPCBuffer(LPCBUFFER_HANDLE phBuff,FLOAT* InputTimeSignals);
LPCBUFFER_ERRORCODE LPCBufferCreate(LPCBUFFER_HANDLE *phBuff,INT32 lStepsize);


/* Main handle to LPC Modules */
struct _LPCANALYSIS_DATA;
typedef struct _LPCANALYSIS_DATA *LPCANALYSIS_HANDLE;



typedef enum
{
	LPCANALYSIS_OK,
	LPCANALYSIS_MEMORY_ALLOCATION_FAILED,
	LPCANALYSIS_SIZEERROR,
	LPCANALYSIS_NOTINITIALIZED

}LPCANALYSIS_ERRORCODE;


typedef struct _LPCANALYSIS_DATA
{

	INT32 lLPCOrder;
	INT32 lWindowSize;
	FLOAT fSmoothFactor;
	FLOAT *pdAutoCorrVector;
	FLOAT *pdLagWindow;
	FLOAT *pdLpcCoeff;
	FLOAT *pdReflectCoeff;
	FLOAT *pdLspCoeff;
	
	FLOAT *pdLpcCoeffQuant;
	FLOAT *pdLspCoeffQuant;
	FLOAT *pdLspCoeffOld;
	FLOAT *fWorkBuffer1;
	FLOAT *fWorkBuffer2;

}LPCANALYSIS_DATA;


LPCANALYSIS_ERRORCODE LPCAnalysisCreate(LPCANALYSIS_HANDLE *phLpc,INT32 lLPCOrder,INT32 lWindowSize);
void LPCAnalysisDelete(LPCANALYSIS_HANDLE *phLpc);

void LPCAnalysis
(
 	LPCANALYSIS_HANDLE hLpc,
	FLOAT *InputTimeSignals,
	INT32 *lLpcOrder,											
	FLOAT *pfParcor,
	FLOAT *pfReflect,
	FLOAT *pfQuantParcor
);

void ConvertLpcToLsp( 
			FLOAT fLpc[], 
			FLOAT fLsp[], 
			FLOAT fOldLsp[],
			INT32 iOrder,
			FLOAT fWorkBuff1[],
			FLOAT fWorkBuff2[]
);

void ConvertLspToLpc( 
		FLOAT fLsp[], 
		FLOAT fLpc[],
		INT32 iOrder,
		FLOAT fWorkBuff1[],
		FLOAT fWorkBuff2[]
);

#endif
