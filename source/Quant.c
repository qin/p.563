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


#include <math.h>
#include <float.h>
#include "defines.h"

#include "Quant.h"
#include "QuantTab.h"


#define PREDSIZE 4
static FLOAT cfPredictor1[PREDSIZE][LPCORDER] = 
{{(FLOAT)0.2570,(FLOAT)0.2780,(FLOAT)0.2800,(FLOAT)0.2736,(FLOAT)0.2757,(FLOAT)0.2764,(FLOAT)0.2675,(FLOAT)0.2678,(FLOAT)0.2779,(FLOAT)0.2647},
{(FLOAT)0.2142,(FLOAT)0.2194,(FLOAT)0.2331,(FLOAT)0.2230,(FLOAT)0.2272,(FLOAT)0.2252,(FLOAT)0.2148,(FLOAT)0.2123,(FLOAT)0.2115,(FLOAT)0.2096},
{(FLOAT)0.1670,(FLOAT)0.1523,(FLOAT)0.1567,(FLOAT)0.1580,(FLOAT)0.1601,(FLOAT)0.1569,(FLOAT)0.1589,(FLOAT)0.1555,(FLOAT)0.1474,(FLOAT)0.1571},
{(FLOAT)0.1238,(FLOAT)0.0925,(FLOAT)0.0798,(FLOAT)0.0923,(FLOAT)0.0890,(FLOAT)0.0828,(FLOAT)0.1010,(FLOAT)0.0988,(FLOAT)0.0872,(FLOAT)0.1060}};

static FLOAT cfPredictor2[PREDSIZE][LPCORDER] = 
{{(FLOAT)0.2360,(FLOAT)0.2405,(FLOAT)0.2499,(FLOAT)0.2495,(FLOAT)0.2517,(FLOAT)0.2591,(FLOAT)0.2636,(FLOAT)0.2625,(FLOAT)0.2551,(FLOAT)0.2310},
{(FLOAT)0.1285,(FLOAT)0.0925,(FLOAT)0.0779,(FLOAT)0.1060,(FLOAT)0.1183,(FLOAT)0.1176,(FLOAT)0.1277,(FLOAT)0.1268,(FLOAT)0.1193,(FLOAT)0.1211},
{(FLOAT)0.0981,(FLOAT)0.0589,(FLOAT)0.0401,(FLOAT)0.0654,(FLOAT)0.0761,(FLOAT)0.0728,(FLOAT)0.0841,(FLOAT)0.0826,(FLOAT)0.0776,(FLOAT)0.0891},
{(FLOAT)0.0923,(FLOAT)0.0486,(FLOAT)0.0287,(FLOAT)0.0498,(FLOAT)0.0526,(FLOAT)0.0482,(FLOAT)0.0621,(FLOAT)0.0636,(FLOAT)0.0584,(FLOAT)0.0794}};

static FLOAT fpPredHistory[PREDSIZE][LPCORDER];

#define PI04            PI*(FLOAT)0.04   
#define PI92            PI*(FLOAT)0.92
#define CONST12         (FLOAT)1.2


/*********************
*
* FUNCTION: CalcWeight
*
* DESCRIPTION:  Calculate distortion wheigths.					 
*					
*
***********************/

static void CalcWeight
(
 FLOAT  fLsf[],        
 FLOAT  fWeight[] 
)
{
   INT32  i;
   FLOAT  tmp;

   tmp = fLsf[1] - PI04 - (FLOAT)1.0;
   if (tmp > (FLOAT)0.0)fWeight[0] = (FLOAT)1.0;
   else         fWeight[0] = tmp * tmp * (FLOAT)10. + (FLOAT)1.0;

   for ( i=1; i<LPCORDER-1; i++ ) {
      tmp = fLsf[i+1] - fLsf[i-1] - (FLOAT)1.0;
      if (tmp > (FLOAT)0.0)    fWeight[i] = (FLOAT)1.0;
      else              fWeight[i] = tmp * tmp * (FLOAT)10. + (FLOAT)1.0;
   }

   tmp = PI92 - fLsf[LPCORDER-2] - (FLOAT)1.0;
   if (tmp > (FLOAT)0.0)       fWeight[LPCORDER-1] = (FLOAT)1.0;
   else         fWeight[LPCORDER-1] = tmp * tmp * (FLOAT)10. + (FLOAT)1.0;

   fWeight[4] *= CONST12;
   fWeight[5] *= CONST12;
   return;
}


/*********************
*
* FUNCTION: InitQuantizer
*
* DESCRIPTION:  Predictor initialisation 
*	
*
***********************/
    
void InitQuantizer()
{
	
  INT32 i;
  INT32 j;


  for(j=0;j<LPCORDER;j++)
  {
		for(i=0; i<PREDSIZE; i++)
		{
			fpPredHistory[i][j]=(PI/(LPCORDER+1))*(j+1);
		}
  }

}

/*********************
*
* FUNCTION: UpdateHistory
*
* DESCRIPTION:  Predictor History update 
*	
***********************/

static void UpdateHistory
(
 FLOAT fLsf[]
)
{
  INT32 i;
  INT32 j;

  for(i=PREDSIZE-1;i>0;i--)
  {
		for(j=0;j<LPCORDER;j++)
		{
			fpPredHistory[i][j]=fpPredHistory[i-1][j];
		}
  }
  for(j=0;j<LPCORDER;j++)
  {
			fpPredHistory[0][j]=fLsf[j];
  }

  return;
}


/*********************
*
* FUNCTION: MAPredictor
*
* DESCRIPTION: MA Predictor
*	
***********************/

static void MAPredictor(
  FLOAT fLsp[],                
  FLOAT fMaCoeff[][LPCORDER],        
  FLOAT fHistory[][LPCORDER],  
  FLOAT fLspRes[]
)
{
  INT32 j;
  INT32 k;
  FLOAT fNorm=0;

  for( j = 0 ; j < LPCORDER ; j++ ) 
  {
	   fNorm=0;
      fLspRes[j]=fLsp[j];
      for ( k = 0 ; k < PREDSIZE ; k++ )
		{
			fLspRes[j] -= fHistory[k][j] * fMaCoeff[k][j];
			fNorm += fMaCoeff[k][j];
		}
		fLspRes[j] /= (1.0F-fNorm);
   }

   return;
}


/*********************
*
* FUNCTION: MAPredictor
*
* DESCRIPTION: Inverse MA Predictor
*	
***********************/

static void InvMAPredictor
(
  FLOAT fLspRes[],                
  FLOAT fMaCoeff[][LPCORDER],        
  FLOAT fHistory[][LPCORDER],  
  FLOAT fLsp[]
)
{
  INT32 j;
  INT32 k;
  FLOAT fNorm=0;

  for( j = 0 ; j < LPCORDER ; j++ ) 
  {
	   fNorm=0;
      fLsp[j]=0;
      for ( k = 0 ; k < PREDSIZE ; k++ )
		{
       
			fLsp[j] += fHistory[k][j] * fMaCoeff[k][j];
			fNorm += fMaCoeff[k][j];
		}
		fLsp[j]+=fLspRes[j]*(1.0F-fNorm);
   }

   return;
}

/*********************
*
* FUNCTION: QuantizeLSF
*
* DESCRIPTION:  Line Spectrum Frequencies Quanitze routine 
*	
***********************/

void QuantizeLSF(
 FLOAT  fLsf[],                
 FLOAT  fQLsf[]
)
{

	INT32 j;

   FLOAT fDist1;
	FLOAT	fDist2;

   FLOAT   fRes[LPCORDER];
	
	FLOAT	fLsfQTemp1[LPCORDER];
	FLOAT	fLsfQTemp2[LPCORDER];

	FLOAT fQuantBuff1[LPCORDER];
	FLOAT fQuantBuff2[LPCORDER];

	FLOAT fWeight[LPCORDER];


	 CalcWeight( fLsf, fWeight );

	MAPredictor(fLsf,  cfPredictor1, fpPredHistory, fRes);
	VQQuant(fRes,fWeight,fQuantBuff1);
	InvMAPredictor(fQuantBuff1, cfPredictor1, fpPredHistory, fLsfQTemp1);
	fDist1=GetDistortion(fLsf,fLsfQTemp1,fWeight,LPCORDER);


   MAPredictor(fLsf,  cfPredictor2, fpPredHistory, fRes);
	VQQuant(fRes,fWeight,fQuantBuff2);
	InvMAPredictor(fQuantBuff2, cfPredictor2, fpPredHistory, fLsfQTemp2);
	fDist2=GetDistortion(fLsf,fLsfQTemp2,fWeight,LPCORDER);


	if(fDist1>fDist2)
	{
			for( j = 0 ; j < LPCORDER ; j++ ) 
			{
				fQLsf[j]=fLsfQTemp2[j];
			}
			UpdateHistory(fQuantBuff2);
	}
	else
	{
			for( j = 0 ; j < LPCORDER ; j++ ) 
			{
				fQLsf[j]=fLsfQTemp1[j];
			}
			UpdateHistory(fQuantBuff1);
	}

    return;
}


/*********************
*
* FUNCTION: QuantizeLSP
*
* DESCRIPTION:  Line Spectrum Paire Quanitze routine 
*	
***********************/
void QuantizeLSP
(
 FLOAT  fLsp[],    
 FLOAT  fQLsp[]
)
{
   INT32 i;
	FLOAT fLsf[LPCORDER];
	FLOAT fQLsf[LPCORDER];

	
  for (i=0; i<LPCORDER; i++ )
  {
     fLsf[i] = (FLOAT)acos(fLsp[i]);
  }

  QuantizeLSF(fLsf, fQLsf);


  for (i=0; i<LPCORDER; i++ )
  {
     fQLsp[i] = (FLOAT)cos(fQLsf[i]);
  }

}

/*********************
*
* FUNCTION: GetDistortion
*
* DESCRIPTION: Evaluate wheigted distortion
*	
***********************/

FLOAT GetDistortion(FLOAT fOrig[],FLOAT fQuant[],FLOAT fWeight[],INT32 iSize)
{
	INT32 i;
	FLOAT fDist=0;
	FLOAT fTemp=0;

	for(i=0;i<iSize;i++)
	{
		fTemp=fOrig[i]-fQuant[i];
		fDist +=fTemp*fTemp*fWeight[i]; 
	}

	return fDist;

}

/*********************
*
* FUNCTION: VQQuant
*
* DESCRIPTION: Vector Quantizer
*	
***********************/

void VQQuant( FLOAT fLsp[],
				  FLOAT fWeight[],
				  FLOAT fLspQuant[]
				)
{

   INT32  i;
	INT32  j;
	INT32  k;
	INT32 iIndex1;
	INT32 iIndex2;

   FLOAT fMinDist;
	FLOAT fActualDist1;
	FLOAT fActualDist;
	FLOAT fTemp;
	INT32 lSizeBook1;
	INT32 lSizeBook2;

   iIndex1 = 0;
   fMinDist= FLT_MAX;



	lSizeBook1=sizeof(pfCodebook1)/sizeof(pfCodebook1[0]);
	lSizeBook2=sizeof(pfCodebook2)/sizeof(pfCodebook2[0]);


   for(i=0; i<lSizeBook1; i++) 
	{
      fActualDist =(FLOAT)0.;
      for(j=0; j<LPCORDER; j++)
		{
        fTemp = fLsp[j]-pfCodebook1[i][j];
        fActualDist +=fTemp * fTemp;
      }
      if(fActualDist<fMinDist)
      {
        fMinDist=fActualDist;
        iIndex1=i;
      }
    }
  
	for(i=0;i<LPCORDER;i++)
	{
		fLspQuant[i]=pfCodebook1[iIndex1][i];
	}

   iIndex1 = 0;
	iIndex2 = 0;
   fMinDist= FLT_MAX;
  
	for(i=0; i<lSizeBook2; i++) 
	{
      fActualDist1 =(FLOAT)0.;
      for(j=0; j<LPCORDER/2; j++)
		{
        fTemp = fLsp[j]-fLspQuant[j]-pfCodebook2[i][j];
        fActualDist1 += fWeight[j]*fTemp * fTemp;
      }

		for(k=0; k<lSizeBook2; k++) 
		{

		   fActualDist=fActualDist1;
			for(j=LPCORDER/2; j<LPCORDER; j++)
			{
				fTemp = fLsp[j]-fLspQuant[j]-pfCodebook2[k][j];
				fActualDist += fWeight[j]*fTemp * fTemp;
			}

			if(fActualDist<fMinDist)
			{
				fMinDist=fActualDist;
				iIndex1=i;
				iIndex2=k;
			}

		}
   
    }
  

	for(i=0;i<LPCORDER/2;i++)
	{
		fLspQuant[i]=fLspQuant[i]+pfCodebook2[iIndex1][i];
	}
	for(i=LPCORDER/2;i<LPCORDER;i++)
	{
		fLspQuant[i]=fLspQuant[i]+pfCodebook2[iIndex2][i];
	}
}

