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


#ifndef __STATISTICS_CLASS__
#define __STATISTICS_CLASS__


#define MAX_ELEMENTS_IN_STAT		500

/* Main handle to STATISTIC Modul */
struct _SISTATISTIC_DATA;
typedef struct _SISTATISTIC_DATA *SISTATISTIC_HANDLE;

/* Errorcode  of STATISTIC modules */
typedef enum
{
	STATISTIC_OK,
	STATISTIC_MEMORY_ALLOCATION_FAILED,
	STATISTIC_SIZEERROR,
	STATISTIC_NOTINITIALIZED

}SISTATISTIC_ERRORCODE;


typedef struct _SISTATISTIC_DATA
{ 	
	char *szName;
	INT32 lNrEntries;			/* */

	FLOAT dSummValues;
	FLOAT dSummSquareValues;
	FLOAT dMin;
	INT32 lMinIndex;
	FLOAT dMax;
	INT32 lMaxIndex;
	UINT32 pulStartIndices[MAX_ELEMENTS_IN_STAT];
	UINT32 pulStopIndices [MAX_ELEMENTS_IN_STAT];
	INT32 lNextIndex;
	INT32 lMaxEntries;
}SISTATISTIC_DATA;
 
SISTATISTIC_ERRORCODE SiStatisticsCreate(SISTATISTIC_HANDLE *phStatistic,char *szName);
 
void SiStatisticsDelete(SISTATISTIC_HANDLE *phStatistic);
void SiStatisticsSet(SISTATISTIC_HANDLE hStatistic, FLOAT dValue);

void SiStatisticsGetMoments(SISTATISTIC_HANDLE hStatistic,FLOAT *Mean,FLOAT *Std,INT32 *NrEntries);

void SiStatisticsSetDefaultVector(SISTATISTIC_HANDLE hStatistic,const FLOAT *Vector, UINT32  ulStart, UINT32 ulStop,INT32 iConvertToEnergy);

#endif /* __STATISTICS_CLASS__*/

