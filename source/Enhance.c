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


#include <assert.h>
#include <math.h>

#include <stdlib.h>
#include <stdio.h>

#include "defines.h"


#include "SpeechLib.h"
#include "LpcAnalysis.h"
#include "Enhance.h"
#include "Quant.h"


struct _FILTERBUFFER_DATA;
typedef struct _FILTERBUFFER_DATA *FILTERBUFFER_HANDLE;

typedef enum
{
	FILTERBUFFER_OK,
	FILTERBUFFER_MEMORY_ALLOCATION_FAILED,
	FILTERBUFFER_SIZEERROR,
	FILTERBUFFER_NOTINITIALIZED

}FILTERBUFFER_ERRORCODE;
 
typedef struct _FRAME_DATA
{

	INT32 FrameNr;
	INT32 lLpcOrder;

	FLOAT fParcor[MAXLPCORDER+1];		
	FLOAT fReflect[MAXLPCORDER+1];		
	FLOAT fQuantParcor[MAXLPCORDER+1];

	INT32 lStartSample;  

	FLOAT *pdAnalysisSignal; 
	FLOAT *pdAnalysisResiduum;
	FLOAT *pdSynthSignal; 


}FRAME_DATA;

typedef struct _FILTERBUFFER_DATA
{

	INT32 lStepsize;		/* */
	INT32 lMaxNrFrames; /* */
	INT32 lBufferDelay;	/* delay */

	FLOAT *pdTimeBufferEnd;
	INT32 NrSamplesInBuffer;
	INT32 NrFramesInBuffer;

	FLOAT *pdTimeBuffer;		/* unfiltert time signals  - full data */
	INT32 lTimeBufferSize;
	FLOAT *pdResiduumBuffer;			/* buffer to store residuum */
	INT32 lResiduumBufferSize;
	FLOAT *pdSynthBuffer;			/* buffer to store residuum */
	INT32 lSynthBufferSize;

	FRAME_DATA *pFrameData; /* info for 1 frame */

}FILTERBUFFER_DATA;





void FilterBufferDelete(FILTERBUFFER_HANDLE *phBuff);
FILTERBUFFER_ERRORCODE FilterBufferCreate(FILTERBUFFER_HANDLE *phBuff,INT32 lSignalSize);


/*********************
*
* FUNCTION: FilterBufferCreate
*
* DESCRIPTION: Initialisation and Memory allocation
*					for LPC fiter 
*	
*
***********************/


FILTERBUFFER_ERRORCODE FilterBufferCreate(FILTERBUFFER_HANDLE *phBuff,INT32 lSignalSize )
{

INT32 lStepsize=40;

INT32 BufferDelay=0;

FILTERBUFFER_ERRORCODE ErrorCode=FILTERBUFFER_OK;
FILTERBUFFER_HANDLE PHand=*phBuff;


	if(PHand != NULL)
	{
		FilterBufferDelete(&PHand);
	}

	if(PHand == NULL)
	{
		PHand = (FILTERBUFFER_HANDLE) calloc( 1,sizeof(FILTERBUFFER_DATA) );

		if( PHand == NULL )
		{
			 ErrorCode=FILTERBUFFER_MEMORY_ALLOCATION_FAILED;
		}


	}

	if(ErrorCode==FILTERBUFFER_OK)
	{

		PHand->lStepsize=lStepsize;

		PHand->NrSamplesInBuffer=0;
		PHand->NrFramesInBuffer=0;
	
		PHand->lMaxNrFrames=lSignalSize/lStepsize;

		PHand->lTimeBufferSize=0;
		PHand->lResiduumBufferSize=0;


		PHand->pdTimeBuffer=NULL;
		PHand->pdTimeBuffer=calloc(lSignalSize,sizeof(FLOAT));
		if(PHand->pdTimeBuffer==NULL)ErrorCode=LPCANALYSIS_MEMORY_ALLOCATION_FAILED;
		PHand->lTimeBufferSize=lSignalSize;


		PHand->pFrameData=calloc(PHand->lMaxNrFrames,sizeof(FRAME_DATA));
		if(PHand->pFrameData==NULL)ErrorCode=LPCANALYSIS_MEMORY_ALLOCATION_FAILED;

		PHand->pdResiduumBuffer=calloc(lSignalSize,sizeof(FLOAT));
		if(PHand->pdResiduumBuffer==NULL)ErrorCode=LPCANALYSIS_MEMORY_ALLOCATION_FAILED;
		PHand->lResiduumBufferSize=lSignalSize;

		PHand->pdSynthBuffer=calloc(lSignalSize,sizeof(FLOAT));
		if(PHand->pdSynthBuffer==NULL)ErrorCode=LPCANALYSIS_MEMORY_ALLOCATION_FAILED;
		PHand->lSynthBufferSize=lSignalSize;
		

	}

	if(ErrorCode==FILTERBUFFER_OK)
	{		
		PHand->lBufferDelay=BufferDelay;
	}


	if(ErrorCode!=FILTERBUFFER_OK)
	{
		FilterBufferDelete(&PHand);
	}
	

	*phBuff=PHand;

	return ErrorCode;
}

/*********************
*
* FUNCTION: FilterBufferDelete
*
* DESCRIPTION:  Memory deallocation
*					for LPC fiter 
*	
*
***********************/

void FilterBufferDelete(FILTERBUFFER_HANDLE *phBuff)
{

	FILTERBUFFER_HANDLE PHand=*phBuff;

	if(PHand != NULL)
	{
		if(PHand->pdTimeBuffer!=NULL)free(PHand->pdTimeBuffer);
		if(PHand->pdResiduumBuffer!=NULL)free(PHand->pdResiduumBuffer);
		if(PHand->pdSynthBuffer!=NULL)free(PHand->pdSynthBuffer);
		if(PHand->pFrameData!=NULL)free(PHand->pFrameData);

		free(PHand);
	}

	*phBuff=NULL;

}

/*********************
*
* FUNCTION: UpdateFilterBuffer
*
* DESCRIPTION:  update filter buffer with next frame data
*	
*
***********************/

static void	UpdateFilterBuffer(FILTERBUFFER_HANDLE phBuff,FRAME_DATA **pFrameData)
{

INT32	lStepSize=phBuff->lStepsize;

FRAME_DATA *pData=&phBuff->pFrameData[phBuff->NrFramesInBuffer];

	assert(phBuff->NrFramesInBuffer < phBuff->lMaxNrFrames);

	pData->lStartSample=phBuff->NrSamplesInBuffer;

	pData->pdAnalysisResiduum=&phBuff->pdResiduumBuffer[phBuff->NrSamplesInBuffer];
	pData->pdAnalysisSignal=&phBuff->pdTimeBuffer[phBuff->NrSamplesInBuffer];
	pData->pdSynthSignal=&phBuff->pdSynthBuffer[phBuff->NrSamplesInBuffer];

	pData->FrameNr=phBuff->NrFramesInBuffer;
	pData->lLpcOrder=0; /* not set yet*/


	phBuff->NrSamplesInBuffer+=lStepSize;
	phBuff->NrFramesInBuffer++;

*pFrameData=pData;
}





/*********************
*
* FUNCTION: ExtractNewData
*
* DESCRIPTION:  extract next frame smples
*	
*
***********************/



INT32 ExtractNewData (FLOAT fSpeechBuffer[], 
							INT32     pNumberOfSamplesToRead,
							FLOAT **pPointerInInputSpeechArray, 
							INT32   *pNumberOfSamplesLeft) 
{

    INT32 numberOfSamplesRead = 0;
    INT32 i;

    for (i = 0; i < pNumberOfSamplesToRead; i++) {
        if (*pNumberOfSamplesLeft > 0) {
            fSpeechBuffer[i] = (FLOAT) *((*pPointerInInputSpeechArray)++);
            
            numberOfSamplesRead++;
            (*pNumberOfSamplesLeft)--;
        }
    }
    return numberOfSamplesRead;
}


/*********************
*
* FUNCTION: SpeechEnhaceFilter
*
* DESCRIPTION:  main speech enhancer module.
*				
*	
*
***********************/


#define MAXSTEPSIZE 128

INT32 SpeechEnhaceFilter ( FLOAT *fpDistortedSignal, 
									  FLOAT *fpEnhancedSignal, 
									  INT32	iNrDistortedSamples, 
									  INT32	*piNrEnhancedSamples,
									  INT32	*lDelay
           ) 
{

    FLOAT *pointerInInputSpeechArray = fpDistortedSignal;
    INT32   numberOfSamplesLeft = iNrDistortedSamples;
    FLOAT *pointerInOutputSpeechArray = fpEnhancedSignal;

	 INT32 lStepSize=0;
    INT32 i;
	 INT32 lFrame=0;

	 FLOAT fpInputBuffer[MAXSTEPSIZE];


	 IIRFILTER_HANDLE hPreprocFilter=NULL; /* preporcess higpas filter */

	 LPCBUFFER_HANDLE hLPCBuffer=NULL;
    LPCANALYSIS_HANDLE hLpcAnalysis=NULL;

	 FILTERBUFFER_HANDLE hFilterBuffer=NULL;
	
		
	 ALLZEROFILTER_HANDLE hResFilter=NULL;  /* residuum calculation */
	 ALLPOLEFILTER_HANDLE hSynthFilter=NULL;  /* residuum calculation */

	 FRAME_DATA *pFrameData;

    *piNrEnhancedSamples = 0;
 
    for (i = 0; i < iNrDistortedSamples; i++) { fpEnhancedSignal[i] = 0;}


	if(FilterBufferCreate(&hFilterBuffer,iNrDistortedSamples) != FILTERBUFFER_OK) return(1);
	lStepSize=hFilterBuffer->lStepsize; /* get Stepsize for easy handling */


	if(LPCBufferCreate(&hLPCBuffer,lStepSize) != LPCBUFFER_OK) return(1);	
	if(LPCAnalysisCreate(&hLpcAnalysis,hLPCBuffer->lLPCOrder,hLPCBuffer->lLPCWindowSize) != LPCANALYSIS_OK) return(1);	
 
	if(AllZeroFilterCreate(&hResFilter,hLPCBuffer->lLPCOrder,lStepSize) != FILTER_OK) return(FILTER_MEMORY_ALLOCATION_FAILED);
	if(AllPoleFilterCreate(&hSynthFilter,hLPCBuffer->lLPCOrder,lStepSize) != FILTER_OK) return(FILTER_MEMORY_ALLOCATION_FAILED);


	if(IIRFilterCreate(&hPreprocFilter,2,2,lStepSize) != FILTER_OK) return(FILTER_MEMORY_ALLOCATION_FAILED);

	{

		FLOAT cpPrepCoefEnum[]={(FLOAT)0.92727435E+00, (FLOAT)-0.18544941E+01, (FLOAT)0.92727435E+00};
 		FLOAT cpPrepCoefDenom[]={(FLOAT)1.00000000E+00, (FLOAT)-0.19059465E+01, (FLOAT)+0.91140240E+00};

		IIRFilterSet(hPreprocFilter,cpPrepCoefEnum,2,cpPrepCoefDenom,2,1);	
	}


	lFrame=0;
		 	
	while (ExtractNewData(fpInputBuffer, lStepSize, &pointerInInputSpeechArray, &numberOfSamplesLeft) == lStepSize) 
   {
         lFrame++;

		
			IIRFilter(hPreprocFilter,fpInputBuffer,fpInputBuffer,lStepSize);

		
			UpdateFilterBuffer(hFilterBuffer,&pFrameData);

			UpdateLPCBuffer(hLPCBuffer,fpInputBuffer);
			LPCAnalysis(hLpcAnalysis,hLPCBuffer->pdLPCWindTimeSignal,						
												&pFrameData->lLpcOrder,
												pFrameData->fParcor,
												pFrameData->fReflect,
												pFrameData->fQuantParcor);

			for(i=0;i<lStepSize;i++)pFrameData->pdAnalysisSignal[i]=hLPCBuffer->pdAnalysisSignal[i];

			AllZeroFilterSet(hResFilter,pFrameData->fParcor,pFrameData->lLpcOrder,0);
			AllZeroFilter(hResFilter, pFrameData->pdAnalysisSignal,pFrameData->pdAnalysisResiduum, lStepSize);

			AllPoleFilterSet(hSynthFilter,pFrameData->fQuantParcor,pFrameData->lLpcOrder,0);
			AllPoleFilter(hSynthFilter, pFrameData->pdAnalysisResiduum,pFrameData->pdSynthSignal, lStepSize);

	 
	}


    for (i = 0; i < hFilterBuffer->NrSamplesInBuffer; i++) 
	  {
           *(pointerInOutputSpeechArray++) = (FLOAT) hFilterBuffer->pdSynthBuffer[i];
           (*piNrEnhancedSamples)++;
	  }  

	 *lDelay=hFilterBuffer->lBufferDelay+hLPCBuffer->lBufferDelay;

	FilterBufferDelete(&hFilterBuffer) ;
	LPCBufferDelete(&hLPCBuffer);
   LPCAnalysisDelete(&hLpcAnalysis);
	AllZeroFilterDelete(&hResFilter);
	AllPoleFilterDelete(&hSynthFilter);

	IIRFilterDelete(&hPreprocFilter);

    return 0;
}


