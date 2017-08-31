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
#include <memory.h>
#include <stdlib.h>
#include "defines.h"
#include "SpeechLib.h"

/*********************
*
* FUNCTION: AllZeroFilterCreate
*
* DESCRIPTION: initalilzation and memory allocation function for a all zero filter
*	
*
***********************/


FILTER_ERRORCODE AllZeroFilterCreate(ALLZEROFILTER_HANDLE *phAllZeroFilter, INT32 lFilterOrder,INT32 MaxStepsize)
{
FILTER_ERRORCODE ErrorCode=FILTER_OK;
ALLZEROFILTER_HANDLE PHand=*phAllZeroFilter;


	if(PHand != NULL){AllZeroFilterDelete(&PHand);}

	if(PHand == NULL)
	{
		PHand = (ALLZEROFILTER_HANDLE) calloc( 1,sizeof(ALLZEROFILTER_DATA) );
		if( PHand == NULL )
		{
			 ErrorCode=FILTER_MEMORY_ALLOCATION_FAILED;
		}
	}

	if(ErrorCode==FILTER_OK)
	{

		PHand->lOrder=lFilterOrder;
		PHand->lMaxStepsize =MaxStepsize;

		PHand->pdFilterBuffer=NULL;
		PHand->pdFilterBuffer=calloc((lFilterOrder+MaxStepsize+1),sizeof(FLOAT));
		if(PHand->pdFilterBuffer==NULL)ErrorCode=FILTER_MEMORY_ALLOCATION_FAILED;

		PHand->pdFilterMemory=NULL;
		PHand->pdFilterMemory=calloc((lFilterOrder+1),sizeof(FLOAT));
		if(PHand->pdFilterMemory==NULL)ErrorCode=FILTER_MEMORY_ALLOCATION_FAILED;
	
		PHand->pdFilterCoeff=NULL;
		PHand->pdFilterCoeff=calloc((lFilterOrder+1),sizeof(FLOAT));
		if(PHand->pdFilterCoeff==NULL)ErrorCode=FILTER_MEMORY_ALLOCATION_FAILED;
	}

	if(ErrorCode!=FILTER_OK)
	{
		AllZeroFilterDelete(&PHand);
	}	
	*phAllZeroFilter=PHand;

	return ErrorCode;
}

/*********************
*
* FUNCTION: AllZeroFilterCreate
*
* DESCRIPTION:  memory deallocation function for a all zero filter
*	
*
***********************/


void AllZeroFilterDelete(ALLZEROFILTER_HANDLE *phAllZeroFilter)
{

	ALLZEROFILTER_HANDLE PHand=*phAllZeroFilter;

	if(PHand != NULL)
	{
		if(PHand->pdFilterBuffer!=NULL)free(PHand->pdFilterBuffer);
		if(PHand->pdFilterMemory!=NULL)free(PHand->pdFilterMemory);
		if(PHand->pdFilterCoeff!=NULL)free(PHand->pdFilterCoeff);

		free(PHand);
	}

	*phAllZeroFilter=NULL;
}

/*********************
*
* FUNCTION: AllZeroFilterSet
*
* DESCRIPTION: preset all zero filter coefficients
*	
*
***********************/


void AllZeroFilterSet(ALLZEROFILTER_HANDLE phAllZeroFilter,
							FLOAT a[],    
							INT32 iOrder,
							INT32 Clear  
)
{
INT32 i;
	assert(iOrder == phAllZeroFilter->lOrder);
	for(i=0;i<=iOrder;i++)
	{
		phAllZeroFilter->pdFilterCoeff[i]=a[i];
	}
	if(Clear!=0)
	{
		for(i=0;i<=iOrder;i++)
		{
			phAllZeroFilter->pdFilterMemory[i]=0;
		}
	}
}

/*********************
*
* FUNCTION: AllZeroFilter
*
* DESCRIPTION: All Zero Filtering of a singnal
*	
*
***********************/


void AllZeroFilter(ALLZEROFILTER_HANDLE phAllZeroFilter,    
						 FLOAT *pInputSignal,			
						 FLOAT *pOutputSignal,	   
						INT32  iNrSmples			
)
{
  FLOAT s;
  INT32  i, j;
  FLOAT *pCoeff=phAllZeroFilter->pdFilterCoeff;				
  INT32 lOrder1=phAllZeroFilter->lOrder+1;
  FLOAT *WorkBuffer=phAllZeroFilter->pdFilterBuffer;

  assert(pInputSignal!=0);
  assert(pOutputSignal!=0);
  assert(phAllZeroFilter->lMaxStepsize >= iNrSmples);

  memcpy(WorkBuffer,phAllZeroFilter->pdFilterMemory,lOrder1*sizeof(FLOAT));
  memcpy(WorkBuffer+lOrder1,pInputSignal,iNrSmples*sizeof(FLOAT));
  
  for (i = lOrder1; i < iNrSmples+lOrder1; i++)
  {
    s = WorkBuffer[i]*pCoeff[0];
    for (j = 1; j < lOrder1; j++)
	 {
		 s += pCoeff[j]*WorkBuffer[i-j];
	 }
    *pOutputSignal++ = s;
  }

  memcpy(phAllZeroFilter->pdFilterMemory,WorkBuffer+iNrSmples,lOrder1*sizeof(FLOAT));

  return;
}


/*********************
*
* FUNCTION: AllPoleFilterCreate
*
* DESCRIPTION: initalilzation and memory allocation function for a all Pole filter
*	
*
***********************/


FILTER_ERRORCODE AllPoleFilterCreate(ALLPOLEFILTER_HANDLE *phAllPoleFilter, INT32 lFilterOrder,INT32 MaxStepsize)
{
FILTER_ERRORCODE ErrorCode=FILTER_OK;
ALLPOLEFILTER_HANDLE PHand=*phAllPoleFilter;


	if(PHand != NULL){AllPoleFilterDelete(&PHand);}

	if(PHand == NULL)
	{
		PHand = (ALLPOLEFILTER_HANDLE) calloc( 1,sizeof(ALLPOLEFILTER_DATA) );
		if( PHand == NULL )
		{
			 ErrorCode=FILTER_MEMORY_ALLOCATION_FAILED;
		}
	}

	if(ErrorCode==FILTER_OK)
	{

		PHand->lOrder=lFilterOrder;
		PHand->lMaxStepsize =MaxStepsize;

		PHand->pdFilterBuffer=NULL;
		PHand->pdFilterBuffer=calloc((lFilterOrder+MaxStepsize+1),sizeof(FLOAT));
		if(PHand->pdFilterBuffer==NULL)ErrorCode=FILTER_MEMORY_ALLOCATION_FAILED;

		PHand->pdFilterMemory=NULL;
		PHand->pdFilterMemory=calloc((lFilterOrder+1),sizeof(FLOAT));
		if(PHand->pdFilterMemory==NULL)ErrorCode=FILTER_MEMORY_ALLOCATION_FAILED;
	
		PHand->pdFilterCoeff=NULL;
		PHand->pdFilterCoeff=calloc((lFilterOrder+1),sizeof(FLOAT));
		if(PHand->pdFilterCoeff==NULL)ErrorCode=FILTER_MEMORY_ALLOCATION_FAILED;

	}

	if(ErrorCode!=FILTER_OK)
	{
		AllPoleFilterDelete(&PHand);
	}	
	*phAllPoleFilter=PHand;

	return ErrorCode;
}


/*********************
*
* FUNCTION: AllPoleFilterDelete
*
* DESCRIPTION:  memory deallocation function for a all Pole filter
*	
*
***********************/

void AllPoleFilterDelete(ALLPOLEFILTER_HANDLE *phAllPoleFilter)
{

	ALLPOLEFILTER_HANDLE PHand=*phAllPoleFilter;

	if(PHand != NULL)
	{
		if(PHand->pdFilterBuffer!=NULL)free(PHand->pdFilterBuffer);
		if(PHand->pdFilterMemory!=NULL)free(PHand->pdFilterMemory);
		if(PHand->pdFilterCoeff!=NULL)free(PHand->pdFilterCoeff);

		free(PHand);
	}

	*phAllPoleFilter=NULL;
}

/*********************
*
* FUNCTION: AllPoleFilterSet
*
* DESCRIPTION:  preset all pole filter coefficients
*	
*
***********************/
void AllPoleFilterSet(ALLPOLEFILTER_HANDLE phAllPoleFilter,
							FLOAT faFilterCoeff[],   
							INT32 iOrder,
							INT32 Clear  
)
{
INT32 i;
	assert(iOrder == phAllPoleFilter->lOrder);
	for(i=0;i<=iOrder;i++)
	{
		phAllPoleFilter->pdFilterCoeff[i]=faFilterCoeff[i];
	}
	if(Clear!=0)
	{
		for(i=0;i<=iOrder;i++)
		{
			phAllPoleFilter->pdFilterMemory[i]=0;
		}
	}
}

/*********************
*
* FUNCTION: AllPoleFilter
*
* DESCRIPTION: All Pole Filtering of a singnal
*	
*
***********************/

void AllPoleFilter(ALLPOLEFILTER_HANDLE phAllPoleFilter,
							FLOAT x[],   
							FLOAT y[],    
							INT32  iNrSmples     
)
{
   INT32  i,j;

   FLOAT s, *yy, *py, *pa;

	FLOAT *WorkBuffer=phAllPoleFilter->pdFilterBuffer;
	FLOAT *Memory=phAllPoleFilter->pdFilterMemory;
	FLOAT *Coeff=phAllPoleFilter->pdFilterCoeff;
   INT32 lOrder=phAllPoleFilter->lOrder;

	assert(x!=0);
   assert(y!=0);
   assert(phAllPoleFilter->lMaxStepsize >= iNrSmples);

   yy = WorkBuffer;
   for (i = 0; i <lOrder; i++)  *yy++ =  *Memory++;

   for (i = 0; i < iNrSmples; i++)
     {
        py=yy;
        pa=Coeff;
        s = *x++;
        for (j = 0; j <lOrder; j++)  s -= (*++pa) * (*--py);
        *yy++ = s;
        *y++ = s;
     }

   for (i = 0; i <lOrder; i++)  *--Memory =*--yy;

   return;
}

/*********************
*
* FUNCTION: IIRFilterCreate
*
* DESCRIPTION: memory allocation an initializstion of an IIR filter
*	
*
***********************/

FILTER_ERRORCODE IIRFilterCreate(IIRFILTER_HANDLE *phIIRFilter, INT32 lOrderEnum,INT32 lOrderDenum,INT32 MaxStepsize)
{
FILTER_ERRORCODE ErrorCode=FILTER_OK;
IIRFILTER_HANDLE PHand=*phIIRFilter;


	if(PHand != NULL){IIRFilterDelete(&PHand);}

	if(PHand == NULL)
	{
		PHand = (IIRFILTER_HANDLE) calloc( 1,sizeof(IIRFILTER_DATA) );
		if( PHand == NULL )
		{
			 ErrorCode=FILTER_MEMORY_ALLOCATION_FAILED;
		}
	}

	if(ErrorCode==FILTER_OK)
	{

		PHand->lOrderEnum=lOrderEnum;
		PHand->lOrderDenum=lOrderDenum;
		PHand->lMaxStepsize =MaxStepsize;

		PHand->hAlZeroFilter=NULL;
		if(AllZeroFilterCreate( &PHand->hAlZeroFilter,lOrderEnum,MaxStepsize) != FILTER_OK)
		{
			ErrorCode=FILTER_MEMORY_ALLOCATION_FAILED;
		}
			PHand->hAlPoleFilter=NULL;
		if(AllPoleFilterCreate( &PHand->hAlPoleFilter,lOrderDenum,MaxStepsize) != FILTER_OK)
		{
			ErrorCode=FILTER_MEMORY_ALLOCATION_FAILED;
		}

	}

	if(ErrorCode!=FILTER_OK)
	{
		IIRFilterDelete(&PHand);
	}	
	*phIIRFilter=PHand;

	return ErrorCode;
}

/*********************
*
* FUNCTION: IIRFilterDelete
*
* DESCRIPTION: memory deallocation of IIR filter
*	
*
***********************/

void IIRFilterDelete(IIRFILTER_HANDLE *phIIRFilter)
{

	IIRFILTER_HANDLE PHand=*phIIRFilter;

	if(PHand != NULL)
	{

		AllZeroFilterDelete(&PHand->hAlZeroFilter);
		AllPoleFilterDelete(&PHand->hAlPoleFilter);

		free(PHand);
	}

	*phIIRFilter=NULL;
}

/*********************
*
* FUNCTION: IIRFilterSet
*
* DESCRIPTION: preset filter coefficient
*	
*
***********************/
void IIRFilterSet(IIRFILTER_HANDLE hIIRFilter,  
						 FLOAT *pCoeffEnum, 
						INT32 iOrderEnum, 
						  FLOAT *pCoeffDenom, 
						INT32 iOrderDenom,		 
						INT32  Clear				 
)
{
	AllZeroFilterSet(hIIRFilter->hAlZeroFilter,pCoeffEnum,iOrderEnum,Clear);
	AllPoleFilterSet(hIIRFilter->hAlPoleFilter,pCoeffDenom,iOrderDenom,Clear);
}

/*********************
*
* FUNCTION: IIRFilter
*
* DESCRIPTION: IIR filtering of a singal
*	
*
***********************/
void IIRFilter(IIRFILTER_HANDLE hIIRFilter,
					FLOAT x[],     
					FLOAT y[],     
					INT32  iNrSamples      
)
{
	AllZeroFilter(hIIRFilter->hAlZeroFilter,x,y,iNrSamples);
	AllPoleFilter(hIIRFilter->hAlPoleFilter,y,y,iNrSamples);
}


