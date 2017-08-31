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
  
Psytechnics, SwissQual or Opticom can provide licences and further information. 
  
Authors: 
      Ludovic Malfait ludovic.malfait@psytechnics.com 
      Roland Bitto rb@opticom.de 
      Pero Juric pero.juric@swissqual.com
 
********************************************************************/


#include <stdlib.h>
#include <math.h>

#include "defines.h"
#include "lpc.h"

/*********************
*
*	FUNCTION: calc_gen_coef
*
*	DESCRIPTION:
*		Calculates the vocal tract parameter frame by frame.
*
***********************/
void calc_gen_coef(FLOAT in[], INT32 indices[], INT32 Nframes,
		FLOAT coef[], INT32 order )
{
    INT32 frame;
    FLOAT * coefp = coef;

    for( frame = 0; frame < Nframes; frame++ ) 
	{
        vocal_tract_param( in + indices[frame*2],
            indices[frame*2+1] - indices[frame*2],
            coefp, order );
        coefp += order;
    }
}

/*********************
*
*	FUNCTION: poly_to_rc
*
*	DESCRIPTION:
*		Converts LPC coefficient into reflection coefficients
*
***********************/
void poly_to_rc( FLOAT A[], FLOAT K[], INT32 order )
{
    INT32 L = order+1;
	FLOAT * B ; 
    INT32 i, j;
    FLOAT r;

	B = (FLOAT *) calloc( order+1, sizeof(FLOAT) );

    for( i = 0; i < L; i++ ) A[i] /= A[0];
    for( i = L-2; i >= 0; i-- )
    {
        K[i] = -A[i+1];
        if( i > 0 )
        {
            for( j = 0; j < L; j++ ) B[j] = A[j];
            r = (FLOAT) (1.0 / (1.0 - K[i]*K[i]));

            for( j = 1; j <= i; j++ )
                A[j] = (A[j] + K[i] * B[i-j+1]) * r;
        }
    }
	free(B) ;
}

/*********************
*
*	FUNCTION: vocal_tract_param
*
*	DESCRIPTION:
*		Calculates the vocal tract parameter from the given frame.
*		
***********************/
void vocal_tract_param( FLOAT in[], INT32 length, FLOAT vtp[], INT32 order )
{
    FLOAT * lpc_c ;
	FLOAT * coef ;
	lpc_c = (FLOAT *) calloc( order+1, sizeof(FLOAT) );
    coef = (FLOAT *) calloc( order, sizeof(FLOAT) );

    LPC_coef( in, length, lpc_c, order );
    poly_to_rc( lpc_c, coef, order );
    rc_to_vtp( coef, vtp, order );

    free( lpc_c );
    free( coef );
}

/*********************
*
*	FUNCTION: rc_to_vtp
*
*	DESCRIPTION:
*		Converts reflection coefficient to tube sections.
*	
***********************/
void rc_to_vtp( FLOAT mu[], FLOAT S[], INT32 order )
{
    FLOAT SM;
    INT32 M;

    SM = 1.0;
    for( M = order-1; M >= 0; M-- )
    {
        SM = (FLOAT) (SM * (1.0+mu[M])/(1.0-mu[M]));
        S[M] = SM;
    }
}

/*********************
*
*	FUNCTION: LPC_coef
*
*	DESCRIPTION:
*		Calculates the LPC using Shur's method
*
***********************/
void LPC_coef( FLOAT in[], INT32 length, FLOAT coef[], INT32 order )
{
	FLOAT * tmp1 ;
	FLOAT * tmp2 ;

    tmp1 = (FLOAT *) calloc( length, sizeof(FLOAT) );
    tmp2 = (FLOAT *) calloc( length, sizeof(FLOAT) );

    hammingWindow( in, tmp1, length );
    acfn( tmp1, length, tmp2, order );
    schur( tmp2, tmp1, order );
    stepUp( tmp1, coef, order );

    free( tmp1 );
    free( tmp2 );
}

/*********************
*
*	FUNCTION: hammingWindow
*
*	DESCRIPTION:
*		Applies and hamming window to the input frame.
*
***********************/
void hammingWindow( FLOAT in[], FLOAT out[], INT32 length )
{
   INT32 i;
   FLOAT n = TWOPI/(length-1);
   for (i=0; i<length; i++) out[i] = (FLOAT) (0.54 - 0.46*cos((FLOAT)i * n)) * in[i];
}

/*********************
*
*	FUNCTION: acfn
*
*	DESCRIPTION:
*		Autocorrelation function
*
***********************/
void acfn( FLOAT in[], INT32 length, FLOAT acf[], INT32 order )
{
   INT32 i, j;
   for( j=0; j<order+1; j++ )
   {
      acf[j]=0.0;
      for( i=0; i<length-j; i++ ) acf[j] += in[i]*in[i+j];
   }
   if( acf[0] == 0.0 ) acf[0]=1.0;
}

/*********************
*
*	FUNCTION: schur
*
*	DESCRIPTION:
*		Schur recursion
*
***********************/
void schur( FLOAT acf[], FLOAT parcor[], INT32 order )
{
	FLOAT *pp, *kk; 
	INT32 i, m, n;
	pp = (FLOAT *) calloc( order+1, sizeof(FLOAT) );
	kk = (FLOAT *) calloc( order+1, sizeof(FLOAT) );

	for( i = 1; i <= order-1; i++ ) kk[order+1-i] = acf[i];

	for( i = 0; i <= order; i++ ) pp[i] = acf[i];
	for( n = 0; n < order; n++ )
	{
		if( pp[0] < fabs( pp[1] ))
		{
			for( i = n; i < order; i++ ) parcor[i] = 0.0;
			free(pp);
			free(kk);
			return;
		}
		parcor[n] = (FLOAT)fabs( pp[1] ) / pp[0];
		if( pp[1] > 0.0 ) parcor[n] = -parcor[n];
		if( n == order-1 )
		{
			free(pp);
			free(kk);
			return;
		}
		pp[0] += pp[1] * parcor[n];
		for( m = 1; m <=( order-1-n ); m++ )
		{
			pp[m] = pp[1+m] + kk[order+1-m] * parcor[n];
			kk[order+1-m] += pp[1+m] * parcor[n];
		}
	}

	free(pp);
	free(kk);
}

/*********************
*
*	FUNCTION: stepUp
*
*	DESCRIPTION:
*		Step up function
*
***********************/
void stepUp( FLOAT parcor[], FLOAT coef[], INT32 order )
{
   FLOAT *w;
   INT32 i, m;
   w = (FLOAT *) calloc( order+1, sizeof(FLOAT) );

   coef[0] = 1.0;
   coef[1] = parcor[0];
   for( m = 2; m <= order; m++ )
   {
      for( i = 1; i < m; i++ ) w[i] = coef[i] + parcor[m-1] * coef[m-i];
      for( i = 1; i < m; i++ ) coef[i] = w[i];
      coef[m] = parcor[m-1];
   }
   free(w);
}
