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


#ifndef __VECTOR_LIB_H
#define __VECTOR_LIB_H

#include "hosm.h"

#ifndef  PI
 #define PI  3.141592653589793238462643383
#endif

#ifdef __cplusplus
 #define CPP "C"
#else
 #define CPP
#endif

static FLOAT fInputArg;
#define SQR(a) ((fInputArg=(a)) == 0.0 ? 0.0 : fInputArg*fInputArg)


typedef struct COMPLEX 
{
  FLOAT r,i;
} complex;


extern CPP void maxv	(FLOAT *a, INT16 as, FLOAT *max, UINT32 *maxo, UINT32 n);
extern CPP void minv	(FLOAT *a, INT16 as, FLOAT *min, UINT32 *mino, UINT32 n);
extern CPP void mini	(INT32	*a, INT16 as, INT32 *min,  UINT32 *mino, UINT32 n);
extern CPP void maxi	(INT32	*a, INT16 as, INT32 *max,  UINT32 *maxo, UINT32 n);
extern CPP void sve		(FLOAT *inp, INT16 astr, FLOAT *res, UINT32 len);
extern CPP void svesq	(FLOAT *inp, INT16 astr, FLOAT *res, UINT32 len);
extern CPP void vsmul	(FLOAT *inp, INT16 astr, FLOAT alpha,FLOAT *out, INT16 bstr, UINT32 len);
extern CPP void vsadd   (FLOAT *inp, INT16 astr, FLOAT alpha,FLOAT *out, INT16 bstr, UINT32 len);
extern CPP void vsub	(FLOAT *inpa,INT16 astr, FLOAT *inpb,INT16 bstr, FLOAT *out, INT16 ostr, UINT32 len);
extern CPP void mcorrel (FLOAT *r,   FLOAT *s, FLOAT *res, UINT32 n);
extern CPP void vclr	(FLOAT *inp, INT16 astr, UINT32 len);
extern CPP void vabs	(FLOAT *inp, INT16 astr, FLOAT *out, INT16 bstr, UINT32 len);
extern CPP void mvemg	(FLOAT *inp, INT16 astr, FLOAT *mean,UINT32 len);
extern CPP void vmov	(FLOAT *inp, INT16 astr, FLOAT *out, INT16 bstr, UINT32 len);
extern CPP void mve		(FLOAT *inp, INT16 astr, FLOAT *mean,UINT32 len);
extern CPP void mvesq	(FLOAT *inp, INT16 astr, FLOAT *mean,UINT32 len);
extern CPP void rmvesq	(FLOAT *inp, INT16 astr, FLOAT *mean,UINT32 len);
extern CPP void vsdiv   (FLOAT *inp, INT16 astr, FLOAT alpha, FLOAT *out, INT16 bstr, UINT32 len);
extern CPP void vfill	(FLOAT fConst,FLOAT *inp,INT16 astr, UINT32 len);
extern CPP void vfloat	(INT16 *iIn, INT16 astr, FLOAT *fOut,INT16 bstr, UINT32 len);
extern CPP void hamm	(FLOAT *in,  INT16 astr, FLOAT *out, INT16 bstr, INT32 length);
extern CPP void cvabs   (complex *inp,INT16 astr,FLOAT *out,INT16 ostr,UINT32 len);
extern CPP void mstdev  (FLOAT *a, INT16 astr, FLOAT *res, UINT32 n);
extern CPP void cepstrum(FLOAT *in, FLOAT *out, INT32 len);
extern CPP void Free    (void ** ptr);
extern CPP void vlogz (FLOAT *inp, INT16 astr, FLOAT alpha, FLOAT *out, INT16 bstr, UINT32 len);

INT32 EvalHigherMoments (FLOAT data[], FLOAT *average, FLOAT *MeanAbsDeviation, FLOAT *StandardDeviation,
						   FLOAT	*var, FLOAT *skewnes, FLOAT *kurtosis, INT32 n);

INT32 LPC_Burg(FLOAT *fInputData, INT32 LPC_order, FLOAT *LPC_coefficients, 
			   FLOAT *MeanSquaredError, INT32 DataLength);


#endif
