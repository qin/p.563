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



/* Standard prototypes */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #include <malloc.h>
#include <limits.h>
#include <assert.h>

/***** local functions *******/
#include "defines.h"
#include "dsp.h"
#include "vector_lib.h"
#include "hosm.h"


#ifndef SWAP
 #define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#endif




/*********************
*
* FUNCTION: EvalHigherMoments
*
* DESCRIPTION:
*	Returns mean average, average deviation MeanAbsDeviation,
*	standard deviation StandardDeviation, variance variance,
*	skewness skewnes, and kurtosis kurtosis of the array data.
*
***********************/

INT32 EvalHigherMoments( FLOAT data[],
						 FLOAT *average,
						 FLOAT *MeanAbsDeviation,
						 FLOAT *StandardDeviation,
						 FLOAT *variance,
						 FLOAT *skewnes,
						 FLOAT *kurtosis,
						 INT32   DataLength )
{
	INT32 j=0;
	FLOAT InterSum=0;
	FLOAT sum=0;
	FLOAT prod=0;

	if (DataLength <= 1)
		return -2;

	sum=0.0;
	for (j=0;j<DataLength;j++)
		sum += data[j];

	*average = sum / DataLength;
	*MeanAbsDeviation = 0;
	*variance	= 0;
	*skewnes	= 0;
	*kurtosis	= 0;

	for (j=0; j < DataLength; j++)
	{
		*MeanAbsDeviation += (FLOAT)(fabs(sum=data[j]-(*average)));
		InterSum += sum;
		prod=sum*sum;
		*variance += prod;

		prod *= sum;
		*skewnes += prod;

		prod *= sum;
		*kurtosis += prod;
	}

	*MeanAbsDeviation /= DataLength;
	*variance = (*variance - InterSum * InterSum / DataLength) / (DataLength-1);
	*StandardDeviation=(FLOAT)(sqrt(*variance));
	if (*variance && *StandardDeviation) {
		*skewnes /= (DataLength*(*variance)*(*StandardDeviation));
		*kurtosis=(FLOAT)((*kurtosis)/(DataLength*(*variance)*(*variance))-3.0);
	} else
		return -1;
	return 0;
}


/*********************
*
* FUNCTION: LPC_Burg
*
* DESCRIPTION:
*	Burg's method for evaluation of the linear prediction coefficients, which involves a
*	recursive procedure for increasing the LPC order per time unit.
*   It also returns the mean square discrepancy as MeanSquaredError.
*
***********************/

INT32 LPC_Burg(FLOAT *fInputData,
			   INT32 LPC_order,
			   FLOAT *LPC_coefficients,
			   FLOAT *MeanSquaredError,
			   INT32 DataLength)
{
	INT32 iCounterA=0;
	INT32 iCounterB=0;
	INT32 iCounterC=0;
	FLOAT fNominator=0.0;
	FLOAT fDenominator=0.0;
	FLOAT coeffs[30];

	FLOAT SqrData=0.0;
	FLOAT *fTempDataA=NULL;
	FLOAT *fTempDataB=NULL;
	FLOAT *fTempDataC=NULL;

	fTempDataA = (FLOAT *)calloc(DataLength, sizeof(FLOAT));
	fTempDataB = (FLOAT *)calloc(DataLength, sizeof(FLOAT));
	fTempDataC = (FLOAT *)calloc(LPC_order,  sizeof(FLOAT));

    if(fTempDataA==NULL || fTempDataB==NULL || fTempDataC==NULL)
		return -1;

	for (iCounterB=1; iCounterB<=DataLength; iCounterB++)
		SqrData += (FLOAT)SQR((FLOAT)(fInputData[iCounterB]));

	*MeanSquaredError = SqrData/DataLength;
	fTempDataA[1]=fInputData[1];
	fTempDataB[DataLength-1]=fInputData[DataLength];

	for (iCounterB=2;iCounterB<=DataLength-1;iCounterB++)
	{
		fTempDataA[iCounterB]=fInputData[iCounterB];
		fTempDataB[iCounterB-1]=fInputData[iCounterB];
	}

	for (iCounterA=1; iCounterA<=LPC_order; iCounterA++)
	{
		fNominator=0.0;
		fDenominator=0.0;

		for (iCounterB=1;iCounterB<=(DataLength-iCounterA);iCounterB++)  {
			fNominator += fTempDataA[iCounterB]*fTempDataB[iCounterB];
			fDenominator += (FLOAT)(SQR((FLOAT)(fTempDataA[iCounterB]))+SQR((FLOAT)(fTempDataB[iCounterB])));
		}


		LPC_coefficients[iCounterA]=(FLOAT)(2.0*fNominator/fDenominator);
		*MeanSquaredError *= (FLOAT)(1-SQR((FLOAT)(LPC_coefficients[iCounterA])));

		for (iCounterC=1;iCounterC<=(iCounterA-1);iCounterC++)
			LPC_coefficients[iCounterC]=fTempDataC[iCounterC]-LPC_coefficients[iCounterA]*fTempDataC[iCounterA-iCounterC];

		if (iCounterA == LPC_order)
		{
			Free((void**)&fTempDataA);
			Free((void**)&fTempDataB);
			Free((void**)&fTempDataC);
			vmov(LPC_coefficients, 1, coeffs, 1, LPC_order);
			return 0;
		}

		for (iCounterC=1;iCounterC<=iCounterA;iCounterC++)
			fTempDataC[iCounterC]=LPC_coefficients[iCounterC];

		for (iCounterB=1;iCounterB<=(DataLength-iCounterA-1);iCounterB++)
		{
			fTempDataA[iCounterB] -= fTempDataC[iCounterA]*fTempDataB[iCounterB];
			fTempDataB[iCounterB]=fTempDataB[iCounterB+1]-fTempDataC[iCounterA]*fTempDataA[iCounterB+1];
		}
	}
	return -2;
}



/*********************
*
* FUNCTION: maxi
*
* DESCRIPTION:
*	Calculates maximum of integer vector elements
*
***********************/

void maxi ( INT32 *a, INT16 as, INT32 *max, UINT32 *maxo, UINT32 n)
  {
  INT32  *data;
  UINT32 i;

  data = a;
  *max = data[0];
  *maxo= 0;


  data += as;
  for (i=1; i<n; i++)
  {
      if(*data > *max)
      {
         *max = *data;
         *maxo = i;
      }
      data += as;
  }
  return;
  }



/*********************
*
* FUNCTION: mini
*
* DESCRIPTION:
*	Calculates minimum of integer vector elements
*
***********************/

void mini ( INT32 *a, INT16 as, INT32 *min, UINT32 *mino, UINT32 n)
  {
  INT32 *data;
  UINT32 cnt;

  data = a;
  *min = data[0];
  *mino= 0;

  data += as;
  for (cnt=1; cnt<n; cnt++)
  {
      if(*data < *min)
      {
         *min = *data;
         *mino = cnt;
      }
      data += as;
  }
  return;
  }




/*********************
*
* FUNCTION: minv
*
* DESCRIPTION:
*	calculate min of FLOAT vector elements
*   [min mino] = minv (a)      as   ..... stride of a
*                              n ........ length of a
*                              min ...... min of a
*                              mino ..... position of min
*
***********************/

void minv ( FLOAT *a, INT16 as, FLOAT *min, UINT32 *mino, UINT32 n)
  {
  UINT32 i;
  FLOAT *data;

  data = a;
  *min = data[0];
  *mino= 0;

  data += as;
  for (i=1; i<n; i++)
  {
      if(*data < *min)
      {
         *min = *data;
         *mino = i;
      }
      data += as;
  }
  return;
  }


/*********************
*
* FUNCTION: maxv
*
* DESCRIPTION:
*	calculate max of FLOAT vector elements
*   [max maxo] = maxv (a)      as   ..... stride of a
*                              n ........ length of a
*                              min ...... min of a
*                              mino ..... position of min
*
***********************/

void maxv ( FLOAT *a, INT16 as, FLOAT *max, UINT32 *maxo, UINT32 n)
  {
  UINT32 i;
  FLOAT *data;

  data = a;
  *max = data[0];
  *maxo= 0;

  data += as;
  for (i=1; i<n; i++)
  {
      if(*data > *max)
      {
         *max = *data;
         *maxo = i;
      }
      data += as;
  }
  return;
  }



/*********************
*
* FUNCTION: mvemg
*
* DESCRIPTION:
*   mean value of vector element magnitudes
*
***********************/

void mvemg(FLOAT *inp, INT16 astr, FLOAT *mean, UINT32 len)
{
  FLOAT sum;
  UINT32 cnt;

  sum = 0.0;
  *mean = *inp;
  if(len==0)
  {
	  return;
  }

  for(cnt=0; cnt<len; cnt++)
  {
  	sum += (FLOAT)fabs(*(inp+cnt*astr));
  }
  sum /= len;
  *mean = sum;
}


/*********************
*
* FUNCTION: mvesq
*
* DESCRIPTION:
*   mean value of vector element squares
*
***********************/

void mvesq(FLOAT *inp, INT16 astr, FLOAT *mean, UINT32 len)
{
  FLOAT sum;
  UINT32 cnt;

  sum = 0.0;
  *mean = *inp;
  if(len==0)
  {
	  return;
  }

  for(cnt=0; cnt<len; cnt++)
  {
  	sum += (*(inp+cnt*astr) * *(inp+cnt*astr));
  }
  sum /= len;
  *mean = sum;
}




/*********************
*
* FUNCTION: mve
*
* DESCRIPTION:
*   mean value of vector elements
*
***********************/

void mve(FLOAT *inp, INT16 astr, FLOAT *mean, UINT32 len)
{
  FLOAT sum;
  UINT32 cnt;

  sum = 0.0;
  *mean = *inp;
  if(len==0)
  {
	  return;
  }

  for(cnt=0; cnt < len; cnt++)
  {
  	sum += *(inp+cnt*astr);
  }
  sum /= len;
  *mean = sum;
}



/*********************
*
* FUNCTION: rmvesq
*
* DESCRIPTION:
*	Root mean square of vector elements
*
***********************/

void rmvesq(FLOAT *inp, INT16 astr, FLOAT *mean, UINT32 len)
{
  FLOAT sum;
  UINT32 cnt;

  sum = 0.0;
  *mean = *inp;
  if(len==0)
  {
	  return;
  }

  for(cnt=0; cnt<len; cnt++)
  {
  	sum += (*(inp+cnt*astr) * *(inp+cnt*astr));
  }
  sum /= len;

  sum = (FLOAT)sqrt(sum);
  *mean = sum;
}




/*********************
*
* FUNCTION: svesq
*
* DESCRIPTION:
*   sum of vector element squares
*
***********************/

void svesq(FLOAT *inp, INT16 astr, FLOAT *res, UINT32 len)
{
  FLOAT sum;
  UINT32 cnt;

  sum = 0.0;
  *res = *inp;
  if(len==0)
  {
	  return;
  }

  for(cnt=0; cnt<len; cnt++)
  {
  	sum += (*(inp+cnt*astr) * *(inp+cnt*astr));
  }
  *res = sum;
}




/*********************
*
* FUNCTION: vsadd
*
* DESCRIPTION:
*   Vector scalar add operation
*	A scalar value "alpha" is added to vector "inp" and a result is saved into "out"
*
***********************/

void vsadd (FLOAT *inp, INT16 astr, FLOAT alpha, FLOAT *out, INT16 bstr, UINT32 len)
{
	UINT32 cnt;
    FLOAT a;
    a = alpha;

	if(len==0)
	{
	  return;
	}

	for(cnt=0; cnt<len; cnt++)
	{
		*(out+cnt*bstr) = (a + *(inp+cnt*astr));
	}
}


/*********************
*
* FUNCTION: vsmul
*
* DESCRIPTION:
*   Vector scalar multiply
*	A scalar value "alpha" is multiplied with a vector "inp" and a result is saved into "out"
*
***********************/

void vsmul (FLOAT *inp, INT16 astr, FLOAT alpha, FLOAT *out, INT16 bstr, UINT32 len)
{
	UINT32 cnt;
    FLOAT a;
    FLOAT mulval;

    a = alpha;

	if(len==0)
	{
	  return;
	}

	for(cnt=0; cnt<len; cnt++)
	{
		mulval = a * *(inp+cnt*astr);
		*(out+cnt*bstr) = mulval;
	}
}



/*********************
*
* FUNCTION: vsub
*
* DESCRIPTION:
*   Vector substraction
*	out=inpa - inpb
*
***********************/

void vsub (FLOAT *inpa, INT16 astr, FLOAT *inpb, INT16 bstr, FLOAT *out, INT16 ostr, UINT32 len)
{
	UINT32 cnt;

	if(len==0)
	{
	  return;
	}

	for(cnt=0; cnt<len; cnt++)
	{
		*(out+cnt*ostr) = *(inpa+cnt*astr) - *(inpb+cnt*bstr);
	}
}



/*********************
*
* FUNCTION: vfill
*
* DESCRIPTION:
*	fill "inp" vector with the scalar values "fConst"
*	inp=fConst
*
***********************/

void vfill ( FLOAT fConst, FLOAT *inp, INT16 astr, UINT32 len)
{
	UINT32 cnt;
	if(len==0)
	{
	  return;
	}

	for(cnt=0; cnt<len; cnt++)
	{
		*(inp+cnt*astr) = fConst;
	}
	return;
}



/*********************
*
* FUNCTION: vlogz
*
* DESCRIPTION:
*	Vector logarithm with zero test
*   if a >  0.0 r=log(a) else r=alpha
*
***********************/

void vlogz (FLOAT *inp, INT16 astr, FLOAT alpha, FLOAT *out, INT16 bstr, UINT32 len)
{
	UINT32 cnt;
	for(cnt=0; cnt<len; cnt++)
	{
      if(*(inp+cnt*bstr) > 0.0)
		*(out+cnt*bstr) = (FLOAT)log(*(inp+cnt*astr));
      else
		*(out+cnt*bstr) = alpha;
	}
}




/*********************
*
* FUNCTION: vlog10z
*
* DESCRIPTION:
*   Vector logarithm, base 10 with zero test
*   if a >  0.0 r=log10(a) else r=alpha
*
***********************/

void vlog10z (FLOAT *inp, INT16 astr, FLOAT alpha, FLOAT *out, INT16 bstr, UINT32 len)
{
	UINT32 cnt;
	for(cnt=0; cnt<len; cnt++)
	{
      if(*(inp+cnt*bstr) > 0.0)
		*(out+cnt*bstr) = (FLOAT)log10(*(inp+cnt*astr));
      else
		*(out+cnt*bstr) = alpha;
	}
}




/*********************
*
* FUNCTION: vmov
*
* DESCRIPTION:
*	vector copy
*	out=inp
*
***********************/

void vmov (FLOAT *inp, INT16 astr, FLOAT *out, INT16 bstr, UINT32 len)
{
	UINT32 cnt;
	for(cnt=0; cnt<len; cnt++)
	{
		*(out+cnt*bstr) = *(inp+cnt*astr);
	}
}



/*********************
*
* FUNCTION: vclr
*
* DESCRIPTION:
*   Vector clear
*	inp=0
*
***********************/

void vclr (FLOAT *inp, INT16 astr, UINT32 len)
{
	UINT32 cnt;
	for(cnt=0; cnt<len; cnt++)
	{
		*(inp+cnt*astr)=0.0;
	}
}




/*********************
*
* FUNCTION: vabs
*
* DESCRIPTION:
*   Vector absolute value
*	out=fabs(inp)
*
***********************/

void vabs (FLOAT *inp, INT16 astr, FLOAT *out, INT16 bstr, UINT32 len)
{
	UINT32 cnt;
	for(cnt=0; cnt<len; cnt++)
	{
		*(out+cnt*bstr) = (FLOAT)fabs(*(inp+cnt*astr));
	}
}




/*********************
*
* FUNCTION: sve
*
* DESCRIPTION:
*   sum of vector elements
*	res=sum(inp)
*
***********************/

void sve(FLOAT *inp, INT16 astr, FLOAT *res, UINT32 len)
{
  FLOAT sum;
  UINT32 cnt;

  sum = 0.0;
  for(cnt=0; cnt<len; cnt++)
  {
  	sum += *(inp+cnt*astr);
  }
  *res = sum;
}



/*********************
*
* FUNCTION: vsdiv
*
* DESCRIPTION:
*	vector scalar divide
*	out=inp/alpha
*
***********************/

void vsdiv (FLOAT *inp, INT16 astr, FLOAT alpha, FLOAT *out, INT16 bstr, UINT32 len)
{
	UINT32 cnt;

	if(alpha==0)
	{
	  return;
	}

	for(cnt=0; cnt<len; cnt++)
	{
		*(out+cnt*bstr) = *(inp+cnt*astr)/alpha;
	}
}



/*********************
*
* FUNCTION: vfloat
*
* DESCRIPTION:
*	integer 16bit to float conversion
*
***********************/

void vfloat(INT16 *iIn, INT16 astr, FLOAT *fOut, INT16 bstr, UINT32 len)
{
	UINT32 cnt;
	for(cnt=0; cnt<len; cnt++)
	{
		*(fOut+cnt*bstr) = (FLOAT)(*(iIn+cnt*astr));
	}
}



/*********************
*
* FUNCTION: hamm
*
* DESCRIPTION:
*	multiplication with a hamming function
*	H=0.54 + 0.46*cos((2*PI*i)/length)
*
***********************/

void hamm (FLOAT *in, INT16 astr, FLOAT *out, INT16 bstr, INT32 length)
  {
    INT32 cntr;
    FLOAT window;
	if(length==0)
	{
	  return;
	}

	for ( cntr = 0; cntr < length; cntr++)
	{
        window = (FLOAT)(0.54 - 0.46*(cos((2.0 * PI *cntr)/length)));
        out[cntr*bstr] = in[cntr*astr] * window;
	}
    return;
  }





/*********************
*
* FUNCTION: mstdev
*
* DESCRIPTION:
*   standard deviation of vector elements
*	s=sqrt((sum(a-amean))/n)
*
***********************/

/*-------------------------------------------------------------------
   standard deviation of vector elements
  -------------------------------------------------------------------*/
void mstdev (FLOAT *a, INT16 astr, FLOAT *res, UINT32 n)
{
   FLOAT sqsum, sum;

   *res = -1.0;
   if(n==0)
   {
	  return;
   }

   sve   (a, astr, &sum,   n);
   svesq (a, astr, &sqsum, n);
   *res = n*sqsum - sum*sum;

   if(*res <= 0.0)
	   *res = 0.0;
   else
	   *res = (FLOAT)sqrt(*res)/n;
}



/*********************
*
* FUNCTION: mcovar
*
* DESCRIPTION:
*   empiric covarianc of vector elements
*
***********************/

void mcovar (FLOAT *x, INT16 xstr, FLOAT *y, INT16 ystr, FLOAT *r, UINT32 n)
{
	FLOAT xmve, ymve;
	UINT32 i;

	*r=0.0;
	if(n==0)
    {
	  return;
    }
	mve  (x, xstr, &xmve, n);
	mve  (y, ystr, &ymve, n);
	for(i=0; i<n; i++)
	{
	   *r += ((*(x+i)-xmve) * (*(y+i)-ymve));
	}
	*r /= n;
}


/*********************
*
* FUNCTION: mcorrel
*
* DESCRIPTION:
*	Empiric correlation coefficient of two dimension set.
*	range:  -1 <= corr <= 1
*
***********************/

void mcorrel (FLOAT *r, FLOAT *s, FLOAT *res, UINT32 n)
{
	FLOAT cov, sdr, sds;

	if(n==0)
    {
	  return;
    }

	mcovar (r, 1, s, 1, &cov, n);
	mstdev (r, 1, &sdr, n);
	mstdev (s, 1, &sds, n);

	if((sdr <= 0) || (sds <= 0))
	{
		*res = -1.0;
	}
	else
	{
		*res = cov/(sdr*sds);
	}
}




/*********************
*
* FUNCTION: cvabs
*
* DESCRIPTION:
*   Complex vector absolute value
*
***********************/

void cvabs(complex *inp,INT16 astr,FLOAT *out,INT16 ostr,UINT32 len)
{
	UINT32 cnt=0;
    complex *in=NULL;
    FLOAT *ou=NULL;
    FLOAT x=0;
	FLOAT y=0;
	FLOAT ans=0;
	FLOAT temp=0;

    in = inp;
    ou = out;

    for(cnt=0; cnt<len; cnt++)
	{
	  x=(FLOAT)fabs(in->r);
	  y=(FLOAT)fabs(in->i);

      if (x == 0.0)
		ans=y;
	  else if (y == 0.0)
		ans=x;
	  else if (x > y)
      {
		temp=y/x;
		ans=x*(FLOAT)sqrt(1.0+temp*temp);
	  } else
      {
		temp=x/y;
		ans=y*(FLOAT)sqrt(1.0+temp*temp);
	  }
      *ou = ans;

      in  += astr;
      ou  += ostr;
	}
}



/*********************
*
* FUNCTION:  cepstrum
*
* DESCRIPTION:
*	calculates the complex cepstrum of the real sequence "in".
*	The Cepstrum separates the glottal frequency from the vocal tract resonances.
*
*	First a logarithmic power spectrum is calculated and then on that
*   an inverse FFT is performed. The result is a signal with a time axis:
*	cepstrum = real(ifft(log2(abs(fft(in))))
*
***********************/

void cepstrum(FLOAT *in, FLOAT *out,INT32 len)
{
	FLOAT *Values=NULL;

	Values = (FLOAT *)calloc(len*2, sizeof(FLOAT));
	if(Values == 0)
	{
		printf("cepstrum: can\'t alocate memory for \'Values\' \n");
		return;
	}

	vclr(Values, 1, len*2);			/* set all to zeros */

	vmov(in, 1, Values, 1, len);	/* save input in a temporary vector Values */

	/* convert real vector into a complex 'out' vector */
	/* imaginar values are zeros. */
	vclr(out, 1, 2*len);		/* set all to zeros */
	vmov(Values, 1, out, 2, len);	/* copy real values */

	FFT(out, len);

	/* calculate absolute values */
   	cvabs ((complex*)(out), 1, Values, 1, len);

	vlogz (Values, 1, 0.0, out, 1, len);

	/* convert real 'out' vector into complex one */
	/* imaginar values are zeros. */
	vclr(Values, 1, 2*len);			/* set all to zeros */

	vmov(out, 1, Values, 2, len);	/* copy real values */

	/* calculate IFFT (inverse Fourier transform) */
	IFFT(Values, len);

	vmov(Values, 2, out, 1, len);

	/* set the second half of 'out' to zero */
	vclr (out+len, 1, len);

	free(Values);

	return;
}




/*********************
*
* FUNCTION: Free
*
* DESCRIPTION:
*	memory free of the blocks allocated by malloc() function
*
***********************/

void Free(void ** ptr)
{
  if (*ptr==0)
    return;
  else
  {
    free(*ptr);
    *ptr = 0;
  }
}
