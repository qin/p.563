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


#include <stdio.h>
#include <math.h>
// #include <malloc.h>
#include <stdlib.h>
#include "defines.h"
#include "Statistics.h"


/*********************
*
* FUNCTION: SiStatisticsCreate
*
* DESCRIPTION: Cration and Memory allocation of a buffer
*              for statistical analysis of an data vector
*
***********************/


SISTATISTIC_ERRORCODE SiStatisticsCreate(SISTATISTIC_HANDLE *phStatistic,char *szName)
{

	SISTATISTIC_ERRORCODE ErrorCode=STATISTIC_OK;
	SISTATISTIC_HANDLE PHand=*phStatistic;

	if(PHand != NULL)
	{
		SiStatisticsDelete(&PHand);
	}

	if(PHand == NULL)
	{
		PHand = (SISTATISTIC_HANDLE) calloc( 1,sizeof(SISTATISTIC_DATA) );

		if( PHand == NULL )
		{
			 ErrorCode=STATISTIC_MEMORY_ALLOCATION_FAILED;
		}
	}

	if(PHand)
	{
		PHand->szName=szName;
		PHand->lNrEntries=0;			/* */
		PHand->dSummValues=0;
		PHand->dSummSquareValues=0;
		PHand->dMin=0;
		PHand->lMinIndex=0;
		PHand->dMax=0;
		PHand->lMaxIndex=0;
		PHand->lNextIndex=0;
		PHand->lMaxEntries = MAX_ELEMENTS_IN_STAT;
	}

	if(ErrorCode!=STATISTIC_OK)
	{
		SiStatisticsDelete(&PHand);
	}

	*phStatistic=PHand;

	return(ErrorCode);

}

/*********************
*
* FUNCTION: SiStatisticsDelete
*
* DESCRIPTION: Done delete the memory used for statistics
*
*
***********************/

void SiStatisticsDelete(SISTATISTIC_HANDLE *phStatistic)
{
SISTATISTIC_HANDLE PHand=*phStatistic;

	if(PHand != NULL)
	{
		free(PHand);
	}
	*phStatistic=NULL;
}



/*********************
*
* FUNCTION: SetStartStopRange
*
* DESCRIPTION: this function is used to check the start and stop point in the analysis vector.
*					the intention is to avoid multipel use of values if the ranges overlap
*
*
***********************/

INT32 SetStartStopRange(SISTATISTIC_HANDLE hStatistic, UINT32 *ulStart, UINT32 *ulStop)
{
	INT32 i;
	for (i=0; i<hStatistic->lNextIndex; i++)
	{

		if (*ulStart<hStatistic->pulStopIndices[i] && *ulStart>hStatistic->pulStartIndices[i])
		{

			*ulStart = hStatistic->pulStopIndices[i];
			if (*ulStart>=*ulStop)
			{
				*ulStop = *ulStart;
				break;
			}
		}


		if (*ulStop<hStatistic->pulStopIndices[i] && *ulStop>hStatistic->pulStartIndices[i])
		{
			*ulStop = hStatistic->pulStartIndices[i];
			if (*ulStart>=*ulStop)
			{
				*ulStop = *ulStart;
				break;
			}
		}
	}

	if (*ulStop>*ulStart) return 1;
	else						 return 0;
}


/*********************
*
* FUNCTION: SiStatisticsSetDefaultVector
*
* DESCRIPTION: Update statistics with multiple data
*
*
***********************/

void SiStatisticsSetDefaultVector(SISTATISTIC_HANDLE hStatistic,const FLOAT *Vector, UINT32  ulStart, UINT32 ulStop,INT32 iConvertToEnergy)
{

	UINT32 i;
	FLOAT dValue=0;

	if(Vector == NULL) return;

	if (hStatistic->lNextIndex < hStatistic->lMaxEntries)
	{
		if (SetStartStopRange(hStatistic, &ulStart, &ulStop))
		{

			hStatistic->pulStartIndices[hStatistic->lNextIndex] = ulStart;
			hStatistic->pulStopIndices [hStatistic->lNextIndex] = ulStop;
			hStatistic->lNextIndex++;


			for(i=ulStart;i<ulStop;i++)
			{
				dValue=Vector[i];
				if(iConvertToEnergy)dValue = dValue*dValue;
				SiStatisticsSet(hStatistic, dValue);
			}
		}
	}
}

/*********************
*
* FUNCTION: SiStatisticsSet
*
* DESCRIPTION: Update statistics with one value
*
*
***********************/

void SiStatisticsSet(SISTATISTIC_HANDLE hStatistic, FLOAT dValue)
{
	hStatistic->dSummValues += dValue;
	hStatistic->dSummSquareValues += dValue*dValue;

	if(hStatistic->lNrEntries == 0)
	{
		hStatistic->dMax = dValue;
		hStatistic->lMaxIndex=0;
		hStatistic->dMin = dValue;
		hStatistic->lMinIndex=0;
	}

	if(hStatistic->dMax < dValue)
	{
		hStatistic->dMax = dValue;
		hStatistic->lMaxIndex=hStatistic->lNrEntries;
	}
	if(hStatistic->dMin > dValue)
	{
		hStatistic->dMin = dValue;
		hStatistic->lMinIndex=hStatistic->lNrEntries;
	}

	hStatistic->lNrEntries++;
}

/*********************
*
* FUNCTION: SiStatisticsGetMoments
*
* DESCRIPTION: now since data aquisition is done,
*					get first and second moment of the distribution
*
*
***********************/

void SiStatisticsGetMoments(SISTATISTIC_HANDLE hStatistic,FLOAT *dMean,FLOAT *dStd,INT32 *lNrEntries)
{
	FLOAT	cv;
	FLOAT dMeanSquare;
	INT32 lCount;


	*dMean=0;
	*dStd=0;
	*lNrEntries=hStatistic->lNrEntries;

	lCount=hStatistic->lNrEntries;

	if(lCount>0)
	{
		*dMean= hStatistic->dSummValues/(FLOAT)lCount;

       dMeanSquare=hStatistic->dSummSquareValues/(FLOAT)lCount;

		if(lCount > 1)
		{
			*dStd= ((FLOAT)lCount/(FLOAT)(lCount-1) )*( dMeanSquare-(*dMean)*(*dMean) );
         if(*dStd > 0.0) *dStd= (FLOAT)sqrt((FLOAT)(*dStd));

	      cv= (*dStd) / (*dMean);

		}
	}

}
