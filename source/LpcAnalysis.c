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


#include <string.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>

#include "defines.h"
 
#include "LpcAnalysis.h"
#include "Quant.h"
 

/*********************
*
* FUNCTION: LPCBufferCreate
*
* DESCRIPTION: Initialisation and Memory allocation
*					for LPC fiter 
*	
*
***********************/

LPCBUFFER_ERRORCODE LPCBufferCreate(LPCBUFFER_HANDLE *phBuff,INT32 lStepsize)
{


INT32 lLPCWindowSize=240; 
INT32 lAnalysisStart=160;
INT32 lLPCOrder=LPCORDER;


INT32 BufferDelay=lLPCWindowSize-lStepsize-lAnalysisStart;

LPCBUFFER_ERRORCODE ErrorCode=LPCBUFFER_OK;
LPCBUFFER_HANDLE PHand=*phBuff;

INT32 i;

	if(PHand != NULL)
	{
		LPCBufferDelete(&PHand);
	}

	if(PHand == NULL)
	{
		PHand = (LPCBUFFER_HANDLE) calloc( 1,sizeof(LPCBUFFER_DATA) );

		if( PHand == NULL )
		{
			 ErrorCode=LPCBUFFER_MEMORY_ALLOCATION_FAILED;
		}


	}

	if(ErrorCode==LPCBUFFER_OK)
	{

		PHand->lStepsize=lStepsize;
		PHand->lLPCWindowSize=lLPCWindowSize;
		PHand->lLPCOrder=lLPCOrder;

		PHand->pdLPCTimeBuffer=NULL;
		PHand->pdLPCTimeBuffer=calloc(lLPCWindowSize,sizeof(FLOAT));
		if(PHand->pdLPCTimeBuffer==NULL)ErrorCode=LPCANALYSIS_MEMORY_ALLOCATION_FAILED;

		PHand->pdLPCTimeWindow=NULL;
		PHand->pdLPCTimeWindow=calloc(lLPCWindowSize,sizeof(FLOAT));
		if(PHand->pdLPCTimeWindow==NULL)ErrorCode=LPCANALYSIS_MEMORY_ALLOCATION_FAILED;

		PHand->pdLPCWindTimeSignal=NULL;
		PHand->pdLPCWindTimeSignal=calloc(lLPCWindowSize,sizeof(FLOAT));
		if(PHand->pdLPCWindTimeSignal==NULL)ErrorCode=LPCANALYSIS_MEMORY_ALLOCATION_FAILED;

	}

	if(ErrorCode==LPCBUFFER_OK)
	{		
		PHand->pdAnalysisSignal=PHand->pdLPCTimeBuffer+lAnalysisStart;

		for(i=0;i<200;i++)
		{
				PHand->pdLPCTimeWindow[i]=(FLOAT)(0.54-0.46*cos(2*PI*i/399));
		}
		for(i=200;i<lLPCWindowSize;i++)
		{
				PHand->pdLPCTimeWindow[i]=(FLOAT)(cos(2*PI*(i-200)/159));
		}		
		PHand->lBufferDelay=BufferDelay;
	}


	if(ErrorCode!=LPCBUFFER_OK)
	{
		LPCBufferDelete(&PHand);
	}
	

	*phBuff=PHand;

	return ErrorCode;
}

/*********************
*
* FUNCTION: LPCBufferDelete
*
* DESCRIPTION: Memory deallocation
*					for LPC fiter 
*	
*
***********************/

void LPCBufferDelete(LPCBUFFER_HANDLE *phBuff)
{

	LPCBUFFER_HANDLE PHand=*phBuff;

	if(PHand != NULL)
	{
		if(PHand->pdLPCTimeBuffer!=NULL)free(PHand->pdLPCTimeBuffer);
		if(PHand->pdLPCTimeWindow!=NULL)free(PHand->pdLPCTimeWindow);
		if(PHand->pdLPCWindTimeSignal!=NULL)free(PHand->pdLPCWindTimeSignal);

		free(PHand);
	}

	*phBuff=NULL;

}

/*********************
*
* FUNCTION: UpdateLPCBuffer
*
* DESCRIPTION: 
*	
*
***********************/
void	UpdateLPCBuffer(LPCBUFFER_HANDLE phBuff,FLOAT* InputTimeSignals)
{
INT32 i;
INT32 lLPCWindowSize=phBuff->lLPCWindowSize;
INT32	lStepSize=phBuff->lStepsize;


/* lpc handling */
	for(i=0;i<lLPCWindowSize-lStepSize;i++)
	{
		phBuff->pdLPCTimeBuffer[i]=phBuff->pdLPCTimeBuffer[i+lStepSize];
	}
	for(i=0;i<lStepSize;i++)
	{
		phBuff->pdLPCTimeBuffer[i+lLPCWindowSize-lStepSize]=InputTimeSignals[i];
	}

	for(i=0;i<lLPCWindowSize;i++)
	{
		phBuff->pdLPCWindTimeSignal[i]=phBuff->pdLPCTimeBuffer[i]*phBuff->pdLPCTimeWindow[i];
	}

}


/*********************
*
* FUNCTION: Autocorr
*
* DESCRIPTION: computes Autocorrelation of a vector
*	
*
***********************/

static void Autocorr(
     FLOAT TimeSigVector[],              
	  INT32 WindowSize,
                     
     FLOAT AutoCorrVector[],              
	  INT32 LpcOrder                
)
{  
   FLOAT sum;
   INT32 i, j;

   for (i = 0; i <= LpcOrder; i++)
   {
     sum = (FLOAT)0.0;
     for (j = 0; j < WindowSize-i; j++)
          sum += TimeSigVector[j]*TimeSigVector[j+i];
     AutoCorrVector[i] = (FLOAT)sum;
   }
   if (AutoCorrVector[0]<(FLOAT)1.0) AutoCorrVector[0]=(FLOAT)1.0;

   return;
}

/*********************
*
* FUNCTION: LagWindow
*
* DESCRIPTION: Autocorrelation vector windowing
*	
*
***********************/


static void LagWindow(                    
     FLOAT AutoCorrVector[],              
	  FLOAT LagWindow[],
	  INT32 LpcOrder               
)
{  
 
   INT32 i;

   for (i = 0; i <= LpcOrder; i++)
   {
		AutoCorrVector[i] *= LagWindow[i];
   }

   return;
}

/*********************
*
* FUNCTION: Levinson
*
* DESCRIPTION: Alevinson dourbin allogrithm
*	
*
***********************/
static FLOAT Levinson(       
 FLOAT AutoCorrCoeff[],  
 INT32 LpcOrder,
 FLOAT LpcCoeff[],       
 FLOAT ReflectCoeff[]             
)
{
   FLOAT s, at, err;
   INT32 i, j, l;

   ReflectCoeff[0] = (-AutoCorrCoeff[1])/AutoCorrCoeff[0];
   LpcCoeff[0] = (FLOAT)1.0;
   LpcCoeff[1] = ReflectCoeff[0];
   err = AutoCorrCoeff[0] + AutoCorrCoeff[1]*ReflectCoeff[0];
 
	for (i = 2; i <= LpcOrder; i++)
   {
     s = (FLOAT)0.0;
     for (j = 0; j < i; j++)
       s += AutoCorrCoeff[i-j]*LpcCoeff[j];
     ReflectCoeff[i-1]= (-s)/(err) ;
     for (j = 1; j <= (i/2); j++)
     {
       l = i-j;
       at = LpcCoeff[j] + ReflectCoeff[i-1]*LpcCoeff[l];
       LpcCoeff[l] += ReflectCoeff[i-1]*LpcCoeff[j];
       LpcCoeff[j] = at;
     }
     LpcCoeff[i] = ReflectCoeff[i-1];
     err += ReflectCoeff[i-1]*s;
     if (err <= (FLOAT)0.0)
        err = (FLOAT)0.001;
   }
   return (err);
}
/*********************
*
* FUNCTION: EvalChebyshev
*
* DESCRIPTION:Evaluation of  chebyshev polinom 
*	
***********************/


static FLOAT EvalChebyshev
(
  FLOAT x,        
  FLOAT *f,       
  INT32 n  
)
{
  FLOAT b1, b2, b0, x2;
  INT32 i;                           
                                 
  x2 = (FLOAT)2.0*x;        
  b2 = (FLOAT)1.0;       
  b1 = x2 + f[1];          
  for (i=2; i<n; i++) 
  {               
    b0 = x2*b1 - b2 + f[i];    
    b2 = b1;                   
    b1 = b0;                        
  }                                
  return (x*b1 - b2 + (FLOAT)0.5*f[n]); 
}


/*********************
*
* FUNCTION: ConvertLpcToLsp
*
* DESCRIPTION: Conversion of LPC coefficients to Line Spectrum Paires
*	
***********************/

#define GRID_POINTS     60 

void ConvertLpcToLsp
(
  FLOAT fLpc[],        
  FLOAT fLsp[],       
  FLOAT fOldLsp[],
  INT32 iOrder,
  FLOAT fWorkBuff1[],
  FLOAT fWorkBuff2[]
)
{
 INT32 i, j;
 INT32 iIndex;
 FLOAT fXlow,fXhigh,fXmid;
 FLOAT fYlow,fYhigh,fYmid;
 FLOAT *fBuffer;

 INT16 iFlag;

 fWorkBuff1[0] = (FLOAT)1.0;
 fWorkBuff2[0] = (FLOAT)1.0;
 for (i=1, j=iOrder; i<=iOrder/2; i++, j--)
 {
    fWorkBuff1[i] = fLpc[i]+fLpc[j]-fWorkBuff1[i-1];
    fWorkBuff2[i] = fLpc[i]-fLpc[j]+fWorkBuff2[i-1];
 }


 iIndex=0;      
 iFlag=0;  
 fBuffer = fWorkBuff1;  
 fXlow=(FLOAT)cos((PI/GRID_POINTS)*0);
 fYlow = EvalChebyshev(fXlow,fBuffer,iOrder/2);

 j = 0;
 while ( (iIndex < iOrder) && (j < GRID_POINTS) )
 {
   j++;
   fXhigh = fXlow;
   fYhigh = fYlow;

	
	fXlow=(FLOAT)cos((PI/GRID_POINTS)*j);

   fYlow = EvalChebyshev(fXlow,fBuffer,iOrder/2);

   if (fYlow*fYhigh <= (FLOAT)0.0) 
   {
     j--;
     for (i = 0; i < 4; i++)
     {
       fXmid = (FLOAT)0.5*(fXlow + fXhigh);
       fYmid = EvalChebyshev(fXmid,fBuffer,iOrder/2);
       if (fYlow*fYmid <= (FLOAT)0.0)
       {
         fYhigh = fYmid;
         fXhigh = fXmid;
       }
       else
       {
         fYlow = fYmid;
         fXlow = fXmid;
       }
     }

     fLsp[iIndex] = fXlow - fYlow*(fXhigh-fXlow)/(fYhigh-fYlow);  

	  if(iFlag==1)
	  {
			fBuffer=fWorkBuff1;
			iFlag=0;
	  }
	  else
	  {
		   fBuffer=fWorkBuff2;
		   iFlag=1;
	  }

     fXlow = fLsp[iIndex];
     fYlow = EvalChebyshev(fXlow,fBuffer,iOrder/2);
	  iIndex++;

   }
 }


 if ( iIndex < iOrder)
    for(i=0; i<iOrder; i++)  fLsp[i] = fOldLsp[i];

 return;
}


/*********************
*
* FUNCTION: FindLspPolinom
*
* DESCRIPTION: Evaluate Lsp Polinom
*	
*
***********************/

static void FindLspPolinom
(
   FLOAT fLsp[],           
   FLOAT fPoly[],
	INT32 iOrder
)
{
  FLOAT fTemp;
  INT32   i,j;
  

  fPoly[0] = (FLOAT)1.0;
  fTemp = (FLOAT)-2.0*fLsp[0];
  fPoly[1] = fTemp;
  for (i = 2; i <= iOrder/2; i++)
  {
    fTemp = (FLOAT)-2.0*fLsp[2*i-2];
    fPoly[i] = fTemp*fPoly[i-1] + (FLOAT)2.0*fPoly[i-2];
    for (j = i-1; j > 1; j--)
	 {
      fPoly[j] += fTemp*fPoly[j-1] + fPoly[j-2];
	 }
    fPoly[1] += fTemp;
  }
  return;
}

/*********************
*
* FUNCTION: ConvertLspToLpc
*
* DESCRIPTION:Conversion of Line spectrum paires to Lpc coefficients
*	
*
***********************/

void ConvertLspToLpc
(
 FLOAT fLsp[],            
 FLOAT fLpc[],
 INT32 iOrder,
 FLOAT fWorkBuff1[],
 FLOAT fWorkBuff2[]
)
{

  INT32 i,j;


  FindLspPolinom(&fLsp[0],fWorkBuff1,iOrder);
  FindLspPolinom(&fLsp[1],fWorkBuff2,iOrder);

  for (i = iOrder/2; i > 0; i--)
  {
    fWorkBuff1[i] += fWorkBuff1[i-1];
    fWorkBuff2[i] -= fWorkBuff2[i-1];
  }
  fLpc[0] = (FLOAT)1.0;
  for (i = 1, j = iOrder; i <= iOrder/2; i++, j--)
  {
    fLpc[i] = (FLOAT)0.5*(fWorkBuff1[i] + fWorkBuff2[i]);
    fLpc[j] = (FLOAT)0.5*(fWorkBuff1[i] - fWorkBuff2[i]);
  }

  return;
}




/*********************
*
* FUNCTION: LPCAnalysisDelete
*
* DESCRIPTION: buffer deallocation for LPC analysis
*	
*
***********************/
void LPCAnalysisDelete(LPCANALYSIS_HANDLE *phLpc)
{

	LPCANALYSIS_HANDLE PHand=*phLpc;

	if(PHand != NULL)
	{
		if(PHand->pdAutoCorrVector!=NULL)free(PHand->pdAutoCorrVector);
		if(PHand->pdLagWindow!=NULL)free(PHand->pdLagWindow);

		if(PHand->pdLpcCoeff!=NULL)free(PHand->pdLpcCoeff);
		if(PHand->pdLpcCoeffQuant!=NULL)free(PHand->pdLpcCoeffQuant);
		if(PHand->pdReflectCoeff!=NULL)free(PHand->pdReflectCoeff);	
		if(PHand->pdLspCoeff!=NULL)free(PHand->pdLspCoeff);	
		if(PHand->pdLspCoeffQuant!=NULL)free(PHand->pdLspCoeffQuant);	
		if(PHand->pdLspCoeffOld!=NULL)free(PHand->pdLspCoeffOld);
		if(PHand->fWorkBuffer1!=NULL)free(PHand->fWorkBuffer1);
		if(PHand->fWorkBuffer2!=NULL)free(PHand->fWorkBuffer2);

		free(PHand);
	}

	*phLpc=NULL;

}

/*********************
*
* FUNCTION: LPCAnalysisDelete
*
* DESCRIPTION: buffer allocation for LPC analysis
*	
*
***********************/

LPCANALYSIS_ERRORCODE LPCAnalysisCreate(LPCANALYSIS_HANDLE *phLpc,INT32 lLPCOrder,INT32 lWindowSize)
{

LPCANALYSIS_ERRORCODE ErrorCode=LPCANALYSIS_OK;
LPCANALYSIS_HANDLE PHand=*phLpc;
INT32 i;
FLOAT expFactor=1;


	assert(lLPCOrder == 10);
	if(lLPCOrder != 10 ) return LPCANALYSIS_SIZEERROR;

	if(PHand != NULL)
	{
		LPCAnalysisDelete(&PHand);
	}

	if(PHand == NULL)
	{
		PHand = (LPCANALYSIS_HANDLE) calloc( 1,sizeof(LPCANALYSIS_DATA) );

		if( PHand == NULL )
		{
			 ErrorCode=LPCANALYSIS_MEMORY_ALLOCATION_FAILED;
		}


	}

	if(ErrorCode==LPCANALYSIS_OK)
	{

		PHand->lLPCOrder=lLPCOrder;
		PHand->lWindowSize=lWindowSize;
		PHand->fSmoothFactor=0;  /* 0 - swith off */

		PHand->pdAutoCorrVector=NULL;
		PHand->pdAutoCorrVector=calloc((lLPCOrder+1),sizeof(FLOAT));
		if(PHand->pdAutoCorrVector==NULL)ErrorCode=LPCANALYSIS_MEMORY_ALLOCATION_FAILED;
		
		PHand->pdLagWindow=NULL;
		PHand->pdLagWindow=calloc((lLPCOrder+1),sizeof(FLOAT));
		if(PHand->pdLagWindow==NULL)
		{
			ErrorCode=LPCANALYSIS_MEMORY_ALLOCATION_FAILED;
		}
		else
		{
			PHand->pdLagWindow[0]=(FLOAT)1.0001;
			for(i=1;i<lLPCOrder+1;i++)
			{
				expFactor=(2*PI*60/8000)*i;
				expFactor=(FLOAT)(-0.5*expFactor*expFactor);
				PHand->pdLagWindow[i]=(FLOAT)exp(expFactor);
			}
		}

		PHand->pdLpcCoeff=NULL;
		PHand->pdLpcCoeff=calloc((lLPCOrder+1),sizeof(FLOAT));
		if(PHand->pdLpcCoeff==NULL)ErrorCode=LPCANALYSIS_MEMORY_ALLOCATION_FAILED;
	
		PHand->pdLpcCoeffQuant=NULL;
		PHand->pdLpcCoeffQuant=calloc((lLPCOrder+1),sizeof(FLOAT));
		if(PHand->pdLpcCoeffQuant==NULL)ErrorCode=LPCANALYSIS_MEMORY_ALLOCATION_FAILED;

		PHand->pdReflectCoeff=NULL;
		PHand->pdReflectCoeff=calloc((lLPCOrder+1),sizeof(FLOAT));
		if(PHand->pdReflectCoeff==NULL)ErrorCode=LPCANALYSIS_MEMORY_ALLOCATION_FAILED;

		PHand->pdLspCoeff=NULL;
		PHand->pdLspCoeff=calloc((lLPCOrder+1),sizeof(FLOAT));
		if(PHand->pdLspCoeff==NULL)ErrorCode=LPCANALYSIS_MEMORY_ALLOCATION_FAILED;
		if(PHand->pdLspCoeff==NULL)
		{
			ErrorCode=LPCANALYSIS_MEMORY_ALLOCATION_FAILED;
		}
		else
		{ 
			PHand->pdLspCoeff[0]=(FLOAT)0.9595;
	      PHand->pdLspCoeff[1]=(FLOAT)0.8413;
			PHand->pdLspCoeff[2]=(FLOAT)0.6549;
			PHand->pdLspCoeff[3]=(FLOAT)0.4154;
			PHand->pdLspCoeff[4]=(FLOAT)0.1423,
			PHand->pdLspCoeff[5]=(FLOAT)-0.1423;
			PHand->pdLspCoeff[6]=(FLOAT)-0.4154;
			PHand->pdLspCoeff[7]=(FLOAT)-0.6549;
			PHand->pdLspCoeff[8]=(FLOAT)-0.8413;
			PHand->pdLspCoeff[9]=(FLOAT)-0.9595;
		}

		PHand->pdLspCoeffQuant=NULL;
		PHand->pdLspCoeffQuant=calloc((lLPCOrder+1),sizeof(FLOAT));
		if(PHand->pdLspCoeffQuant==NULL)ErrorCode=LPCANALYSIS_MEMORY_ALLOCATION_FAILED;

		PHand->fWorkBuffer1=NULL;
		PHand->fWorkBuffer1=calloc((lLPCOrder+1),sizeof(FLOAT));
		if(PHand->fWorkBuffer1==NULL)ErrorCode=LPCANALYSIS_MEMORY_ALLOCATION_FAILED;

		PHand->fWorkBuffer2=NULL;
		PHand->fWorkBuffer2=calloc((lLPCOrder+1),sizeof(FLOAT));
		if(PHand->fWorkBuffer2==NULL)ErrorCode=LPCANALYSIS_MEMORY_ALLOCATION_FAILED;

		PHand->pdLspCoeffOld=NULL;
		PHand->pdLspCoeffOld=calloc((lLPCOrder+1),sizeof(FLOAT));

	}

	if(ErrorCode==LPCANALYSIS_OK)
	{	
		InitQuantizer();  
	}


	if(ErrorCode!=LPCANALYSIS_OK)
	{
		LPCAnalysisDelete(&PHand);
	}
	

	*phLpc=PHand;

	return ErrorCode;
}

/*********************
*
* FUNCTION: LPCAnalysis
*
* DESCRIPTION: main function of LPC analysis
*	
*
***********************/

void LPCAnalysis
(
	LPCANALYSIS_HANDLE hLpc,	
	FLOAT *InputTimeSignals,
	INT32 *lpLpcOrder,
	FLOAT *pfParcor,
	FLOAT *pfReflect,
	FLOAT *pfQuantParcor	
		
)
{
 
	
	FLOAT PredError;
	INT32 i;

	Autocorr(InputTimeSignals,hLpc->lWindowSize,hLpc->pdAutoCorrVector,hLpc->lLPCOrder);
	LagWindow(hLpc->pdAutoCorrVector,hLpc->pdLagWindow,hLpc->lLPCOrder);
	PredError=Levinson(hLpc->pdAutoCorrVector,hLpc->lLPCOrder,hLpc->pdLpcCoeff,hLpc->pdReflectCoeff);


	for(i=0;i<hLpc->lLPCOrder;i++)hLpc->pdLspCoeffOld[i]=hLpc->pdLspCoeff[i];

	ConvertLpcToLsp(hLpc->pdLpcCoeff,hLpc->pdLspCoeff,hLpc->pdLspCoeffOld,hLpc->lLPCOrder,hLpc->fWorkBuffer1,hLpc->fWorkBuffer2);
	
	QuantizeLSP(hLpc->pdLspCoeff, hLpc->pdLspCoeffQuant);

	ConvertLspToLpc(hLpc->pdLspCoeffQuant,hLpc->pdLpcCoeffQuant,hLpc->lLPCOrder,hLpc->fWorkBuffer1,hLpc->fWorkBuffer2);

	if(lpLpcOrder != 0) *lpLpcOrder=hLpc->lLPCOrder;

	if(pfParcor != 0)for(i=0;i<hLpc->lLPCOrder+1;i++)pfParcor[i]=hLpc->pdLpcCoeff[i];
	if(pfReflect != 0)for(i=0;i<hLpc->lLPCOrder+1;i++)pfReflect[i]=hLpc->pdReflectCoeff[i];
	if(pfQuantParcor != 0)for(i=0;i<hLpc->lLPCOrder+1;i++)pfQuantParcor[i]=hLpc->pdLpcCoeffQuant[i];

}

