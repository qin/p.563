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


#ifndef __SPEECHLIB_H_
#define __SPEECHLIB_H_


typedef enum
{
	FILTER_OK,
	FILTER_MEMORY_ALLOCATION_FAILED,
	FILTER_SIZEERROR,
	FILTER_NOTINITIALIZED

}FILTER_ERRORCODE;


struct _ALLPOLEFILTER_DATA;
typedef struct _ALLPOLEFILTER_DATA *ALLPOLEFILTER_HANDLE;

typedef struct _ALLPOLEFILTER_DATA
{

	INT32 lOrder;
	INT32 lMaxStepsize;

	FLOAT *pdFilterBuffer;
	FLOAT *pdFilterMemory;
	FLOAT *pdFilterCoeff;

}ALLPOLEFILTER_DATA;

FILTER_ERRORCODE AllPoleFilterCreate(ALLPOLEFILTER_HANDLE *phAllPoleFilter, INT32 lFilterOrder,INT32 MaxStepsize);
void AllPoleFilterDelete(ALLPOLEFILTER_HANDLE *phAllPoleFilter);

struct _ALLZEROFILTER_DATA;
typedef struct _ALLZEROFILTER_DATA *ALLZEROFILTER_HANDLE;

typedef struct _ALLZEROFILTER_DATA
{
	INT32 lOrder;
	INT32 lMaxStepsize;

	FLOAT *pdFilterBuffer;
	FLOAT *pdFilterMemory;
	FLOAT *pdFilterCoeff;

}ALLZEROFILTER_DATA;


struct _IIRFILTER_DATA;
typedef struct _IIRFILTER_DATA *IIRFILTER_HANDLE;

FILTER_ERRORCODE AllZeroFilterCreate(ALLZEROFILTER_HANDLE *phAllZeroFilter, INT32 lFilterOrder,INT32 MaxStepsize);
void AllZeroFilterDelete(ALLZEROFILTER_HANDLE *phAllZeroFilter);


typedef struct _IIRFILTER_DATA
{
	INT32 lOrderEnum;
	INT32 lOrderDenum;
	INT32 lMaxStepsize;

	ALLZEROFILTER_HANDLE hAlZeroFilter;
	ALLPOLEFILTER_HANDLE hAlPoleFilter;

}IIRFILTER_DATA;

FILTER_ERRORCODE IIRFilterCreate(IIRFILTER_HANDLE *phIIRFilter, INT32 lOrderEnum,INT32 lOrderDenum,INT32 MaxStepsize);
void IIRFilterDelete(IIRFILTER_HANDLE *phIIRFilter);


void AllZeroFilterSet(ALLZEROFILTER_HANDLE hAllZeroFilter,  
						
 FLOAT *pCoeff,			
 INT32 iOrder,		
 INT32  Clear			
);

void AllZeroFilter(ALLZEROFILTER_HANDLE hAllZeroFilter,    
				
 FLOAT *pInputSignal,			
 FLOAT *pOutputSignal,			                        
 INT32  iNrSmples			
);

void AllPoleFilterSet(ALLPOLEFILTER_HANDLE hAllPoleFilter, 
						
 FLOAT *pCoeff,			
 INT32 iOrder,		
 INT32  Clear			
);

void AllPoleFilter(ALLPOLEFILTER_HANDLE hAllPoleFilter,
 FLOAT x[],     
 FLOAT y[],   
 INT32  iNrSmples   
);

void IIRFilterSet(IIRFILTER_HANDLE hIIRFilter,  
 FLOAT *pCoeffEnum,			
 INT32 iOrderEnum,		
  FLOAT *pCoeffDenom,			
 INT32 iOrderDenom,		
 INT32  Clear				
);
void IIRFilter(IIRFILTER_HANDLE hIIRFilter,
 FLOAT x[],     
 FLOAT y[],    
 INT32  iNrSmples   
);



void GetEnergy
(
	FLOAT *pInputSignal,
	INT32 SizeSignal,
	FLOAT *fpEnergy
);

#endif

