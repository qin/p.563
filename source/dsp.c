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
#include <stdlib.h>
#include <limits.h>

#undef VTUNE

#ifdef VTUNE
#define nsp_UsesAll                       /* nsp library functions  */
#include <nsp.h>
#endif 


#include "defines.h"
#include "dsp.h"

/* Variables local to this module - not in header. */
UINT32   FFTSwapInitialised = 0;
UINT32   FFTLog2N;
UINT32 * FFTButter;
UINT32 * FFTBitSwap;
FLOAT  * FFTPhi;


/************************************************************
*    Function definitions                                   *
************************************************************/

UINT32 nextpow2(UINT32 X)
/************************************************************
*   Returns the next power of 2 greater than X,
*   or equal to X if X is a power of 2.
*
*   Checked 22/10/1997 - Antony Rix.
************************************************************/
{
    UINT32 C = 1;
    while( (C < ULONG_MAX) && (C < X) )
        C <<= 1;

    return C;
}

INT ispow2(UINT32 X)
/************************************************************
*   Returns non-zero if X is a power of 2, 0 otherwise.
*
*   Checked 22/10/1997 - Antony Rix.
************************************************************/
{
    UINT32 C = 1;
    while( (C < ULONG_MAX) && (C < X) )
        C <<= 1;
        
    return (C == X);
}

INT intlog2(UINT32 X)
/************************************************************
*   Returns log of X to the power of 2, rounded to the
*   nearest integer.  Calls log().
*
*   Checked 22/10/1997 - Antony Rix.
*   Modified to use log() 13/01/1998 - Antony Rix
************************************************************/
{
    return (INT) floor( log( 1.0 * X ) / log( 2.0 ) + 0.5 );
}

void FFTInit(UINT32 N)
/************************************************************
*   Initialises arrays FFTButter and FFTPhi for use by FFT.
*
*   FFTButter contains the order of each butterfly element in
*   the FFT, given by the bit-reversal of 0 to N/2-1.
*
*   FFTPhi holds the complex angle exp(j * theta) for N/2
*   points spaced by 2*pi/N in the range [0, pi).
*
*   Gives correct results through FFT; checked for basic
*   logic. 22/10/1997 - Antony Rix.  Note that contrary to
*   Watcom documentation, the default matherr() does nothing.
*   When exceptions are raised the program simply terminates
*   with return code 1 if compiled in Watcom C 10.6.  For
*   portability, the exception code has therefore been
*   disabled (13/01/1998 - Antony Rix).
************************************************************/
{
    UINT32   C, L, K;
    FLOAT           Theta;
    FLOAT         * PFFTPhi;
    
    /* Check and clear if space already allocated */
    if( (FFTSwapInitialised != N) && (FFTSwapInitialised != 0) )
        FFTFree();

    if( FFTSwapInitialised == N )
    /* Values already calculated, so return. */
    {
        return;
    }
    else
    {
        /* Compute log to base 2 */
        C = N;
        for( FFTLog2N = 0; C > 1; C >>= 1 )
            FFTLog2N++;

        /* Check that N is a power of 2 */
        C = 1;
        C <<= FFTLog2N;
        if( N == C )
            FFTSwapInitialised = N;
        /* exception raise disabled.
        **else
        **{
        **    dsp_exception.type = DOMAIN;
        **    dsp_exception.name = "FFTInit";
        **    dsp_exception.arg1 = N;
        **    dsp_exception.arg2 = 0;
        **    dsp_exception.retval = N;
        **    matherr( &dsp_exception );
        **    exit( 1 );
        **}
        */
        /* Allocate space */
        FFTButter = (UINT32 *) malloc( sizeof(UINT32) * (N >> 1) );
        FFTBitSwap = (UINT32 *) malloc( sizeof(UINT32) * N );
        FFTPhi = (FLOAT *) malloc( 2 * sizeof(FLOAT) * (N >> 1) );
    
        /* Initialise complex angles */
        PFFTPhi = FFTPhi;
        for( C = 0; C < (N >> 1); C++ )
        {
            Theta = (TWOPI * C) / N;
            (*(PFFTPhi++)) =(FLOAT) cos( Theta );
            (*(PFFTPhi++)) =(FLOAT) sin( Theta );
        }
    
        /* Initialise array containing bit-reversed indices */
        FFTButter[0] = 0;
        L = 1;
        K = N >> 2;
        while( K >= 1 )
        {
            for( C = 0; C < L; C++ )
                FFTButter[C+L] = FFTButter[C] + K;
            L <<= 1;
            K >>= 1;
        }
    }
}

void FFTFree(void)
/************************************************************
*   Clears arrays FFTButter and FFTPhi if already allocated.
*
*   Checked 22/10/1997 - Antony Rix.
************************************************************/
{
    if( FFTSwapInitialised != 0 )
    {
        free( FFTButter );
        free( FFTBitSwap );
        free( FFTPhi );
        FFTSwapInitialised = 0;
    }
}

void FFT(FLOAT * x, UINT32 N)
/************************************************************
*   FFT - Fast Fourier Transform.
*   Perform the FFT of data in x, which is returned in situ
*   in order of increasing frequency.  x should contain N
*   complex pairs, real part first.
*
*   N must be a power of 2.
*
*   FFTFree must be called after the last use of FFT, IFFT,
*   or any FFT-based routine in this library.
*
*   Checked to give correct answers within rounding error
*   using Matlab; fixed N = 1 bug.  22/10/1997 - Antony Rix.
************************************************************/
{
    UINT32   Cycle, C, S, NC;
    UINT32   Step    = N >> 1;
    UINT32   K1, K2;
    register FLOAT  R1, I1, R2, I2;
    FLOAT           ReFFTPhi, ImFFTPhi;

    /* If N = 1 there is nothing to do and FFTInit must not be called */
    if( N > 1 )
    {
        /* Initialise the arrays */
        FFTInit( N );
    
        /* Perform FFT returning frequency points in bit-swapped
           order - they are then sorted at the end. */
        for( Cycle = 1; Cycle < N; Cycle <<= 1, Step >>= 1 )
        {
            /* K1 and K2 point to the two complex pairs that
               are processed at each step. */
            K1 = 0;
            K2 = Step << 1;
    
            for( C = 0; C < Cycle; C++ )
            {
                if( C == 0 )
                {
                    /* This special case needs no multiplying as
                       the value of Theta is exactly zero. */
                    for( S = 0; S < Step; S++ )
                    {
                        R1 = x[K1];
                        I1 = x[K1+1];
                        R2 = x[K2];
                        I2 = x[K2+1];
                        
                        x[K1++] = R1 + R2;
                        x[K1++] = I1 + I2;
                        x[K2++] = R1 - R2;
                        x[K2++] = I1 - I2;
                    }
                }
                else if( C == 1 )
                    /* This case also needs no multiplying as
                       Theta is pi/2, a rotation of 90 degrees. */
                    for( S = 0; S < Step; S++ )
                    {
                        R1 = x[K1];
                        I1 = x[K1+1];
                        R2 = x[K2];
                        I2 = x[K2+1];
                        
                        x[K1++] = R1 + I2;
                        x[K1++] = I1 - R2;
                        x[K2++] = R1 - I2;
                        x[K2++] = I1 + R2;
                    }
                else
                {
                    NC = FFTButter[C] << 1;
                    ReFFTPhi = FFTPhi[NC];
                    ImFFTPhi = FFTPhi[NC+1];
                    for( S = 0; S < Step; S++ )
                    {
                        R1 = x[K1];
                        I1 = x[K1+1];
                        R2 = x[K2];
                        I2 = x[K2+1];
                        
                        x[K1++] = R1 + ReFFTPhi * R2 + ImFFTPhi * I2;
                        x[K1++] = I1 - ImFFTPhi * R2 + ReFFTPhi * I2;
                        x[K2++] = R1 - ReFFTPhi * R2 - ImFFTPhi * I2;
                        x[K2++] = I1 + ImFFTPhi * R2 - ReFFTPhi * I2;
                    }
                }
    
                /* Update K1 and K2 for the next cycle of C. */
                K1 = K2;
                K2 = K1 + (Step << 1);
            }
        }
    
        /* Untangle bit-reversal */
        NC = N >> 1;
        for( C = 0; C < NC; C++ )
        {
            FFTBitSwap[C] = FFTButter[C] << 1;
            FFTBitSwap[C+NC] = 1 + FFTBitSwap[C];
        }
        for( C = 0; C < N; C++ )
            if( (S = FFTBitSwap[C]) != C )
            {
                FFTBitSwap[S] = S;     /* Ensures that only one swap is made */
    
                K1 = C << 1;
                K2 = S << 1;
                R1 = x[K1];
                x[K1++] = x[K2];
                x[K2++] = R1;
                R1 = x[K1];
                x[K1] = x[K2];
                x[K2] = R1;
            }
        /* FFT Completed. */
    }
}

void IFFT(FLOAT * x, UINT32 N)
/************************************************************
*   IFFT - Inverse Fast Fourier Transform.
*   Perform the IFFT of data in x, which is returned in situ
*   in order of increasing time.  x should contain N complex
*   pairs, real part first, in order of increasing frequency.
*
*   N must be a power of 2.
*
*   FFTFree must be called after the last use of FFT, IFFT,
*   or any FFT-based routine in this library.
*
*   Checked to give correct answers within rounding error
*   using Matlab.  22/10/1997 - Antony Rix.
************************************************************/
{
    UINT32   Cycle, C, S, NC;
    UINT32   Step    = N >> 1;
    UINT32   K1, K2;
    register FLOAT  R1, I1, R2, I2;
    FLOAT           ReFFTPhi, ImFFTPhi;

    /* If N = 1 there is nothing to do and FFTInit must not be called */
    if( N > 1 )
    {
        /* Initialise the arrays */
        FFTInit( N );
    
        /* Perform IFFT returning frequency points in bit-swapped
           order - they are then sorted at the end. */
        for( Cycle = 1; Cycle < N; Cycle <<= 1, Step >>= 1 )
        {
            /* K1 and K2 point to the two complex pairs that
               are processed at each step. */
            K1 = 0;
            K2 = Step << 1;
    
            for( C = 0; C < Cycle; C++ )
            {
                if( C == 0 )
                {
                    /* This special case needs no multiplying as
                       the value of Theta is exactly zero. */
                    for( S = 0; S < Step; S++ )
                    {
                        R1 = x[K1];
                        I1 = x[K1+1];
                        R2 = x[K2];
                        I2 = x[K2+1];
                        
                        x[K1++] = R1 + R2;
                        x[K1++] = I1 + I2;
                        x[K2++] = R1 - R2;
                        x[K2++] = I1 - I2;
                    }
                }
                else if( C == 1 )
                    /* This case also needs no multiplying as
                       Theta is pi/2, a rotation of 90 degrees. */
                    for( S = 0; S < Step; S++ )
                    {
                        R1 = x[K1];
                        I1 = x[K1+1];
                        R2 = x[K2];
                        I2 = x[K2+1];
                        
                        x[K1++] = R1 - I2;
                        x[K1++] = I1 + R2;
                        x[K2++] = R1 + I2;
                        x[K2++] = I1 - R2;
                    }
                else
                {
                    NC = FFTButter[C] << 1;
                    ReFFTPhi = FFTPhi[NC];
                    ImFFTPhi = FFTPhi[NC+1];
                    for( S = 0; S < Step; S++ )
                    {
                        R1 = x[K1];
                        I1 = x[K1+1];
                        R2 = x[K2];
                        I2 = x[K2+1];
                        
                        x[K1++] = R1 + ReFFTPhi * R2 - ImFFTPhi * I2;
                        x[K1++] = I1 + ImFFTPhi * R2 + ReFFTPhi * I2;
                        x[K2++] = R1 - ReFFTPhi * R2 + ImFFTPhi * I2;
                        x[K2++] = I1 - ImFFTPhi * R2 - ReFFTPhi * I2;
                    }
                }
    
                /* Update K1 and K2 for the next cycle of C. */
                K1 = K2;
                K2 = K1 + (Step << 1);
            }
        }
    
        /* Untangle bit-reversal */
        NC = N >> 1;
        for( C = 0; C < NC; C++ )
        {
            FFTBitSwap[C] = FFTButter[C] << 1;
            FFTBitSwap[C+NC] = 1 + FFTBitSwap[C];
        }
        for( C = 0; C < N; C++ )
            if( (S = FFTBitSwap[C]) != C )
            {
                FFTBitSwap[S] = S;     /* Ensures that only one swap is made */
    
                K1 = C << 1;
                K2 = S << 1;
                R1 = x[K1];
                x[K1++] = x[K2];
                x[K2++] = R1;
                R1 = x[K1];
                x[K1] = x[K2];
                x[K2] = R1;
            }
    
        /* Divide through by N. */
        NC = N << 1;
        for( C = 0; C < NC; )
            x[C++] /= N;
    
        /* IFFT Completed. */
    }
}

UINT32 FFTNXCorr(
  FLOAT * x1, UINT32 n1,
  FLOAT * x2, UINT32 n2,
  FLOAT * y )
/************************************************************
*   FFTNXCorr - perform cross-correlation using FFT.
*   Returns Conv(reverse(x1), x2) for general n1, n2.
*
*   x1 and x2 contain the two vectors to be correlated, of
*   length n1 and n2 samples.  The correlation product is
*   returned in y, which must hold space for n1+n2-1 samples.
*   The actual length of the correlation product, n1+n2-1,
*   is returned by the function.
*
*   For delay identification it is recommended that any DC
*   level be removed before calling this function.  Identical
*   x1 and x2 will give a maximum in y[Nx-1]; delay of x2
*   with respect to x1 will move this maximum back (to a
*   higher index).
*
*   FFTFree must be called after the last use of FFT, IFFT,
*   or any FFT-based routine in this library.
*
*   Checked to give correct answers within rounding error
*   using Matlab.  22/10/1997 - Antony Rix.
************************************************************/
{
    register FLOAT  r1, i1;
    FLOAT         * tmp1;
    FLOAT         * tmp2;
    INT32            C, D, Nx, Ny;

    /* Allocate temporary storage */
    Nx = nextpow2( max(n1, n2) );
    tmp1 = (FLOAT *) malloc( 2 * sizeof(FLOAT) * 2 * Nx );
    tmp2 = (FLOAT *) malloc( 2 * sizeof(FLOAT) * 2 * Nx );

    /* Copy data, pad with zeros, and perform FFTs. */
    for( C = n1 - 1; C >= 0; C-- )
    {                                   /* Time-reversal is */
        tmp1[C << 1] = *(x1++);         /* implemented here */
        tmp1[1 + (C << 1)] = 0.0;
    }
    for( C = n1 << 1; C < (Nx << 2); C++ )
        tmp1[C] = 0.0;
    FFT( tmp1, 2*Nx );
    
    for( C = 0; C < (INT32) n2; C++ )
    {
        tmp2[C << 1] = x2[C];
        tmp2[1 + (C << 1)] = 0.0;
    }
    for( C = n2 << 1; C < (Nx << 2); C++ )
        tmp2[C] = 0.0;
    FFT( tmp2, 2*Nx );

    /* Convolve in frequency domain. */
    for( C = 0; C < (Nx << 1); C++ )
    {
        D = C << 1; r1 = tmp1[D]; i1 = tmp1[1 + D];
        tmp1[D] = r1 * tmp2[D] - i1 * tmp2[1 + D];
        tmp1[1 + D] = r1 * tmp2[1 + D] + i1 * tmp2[D];
    }

    /* Inverse FFT and copy to position. */
    IFFT( tmp1, 2*Nx );
    Ny = n1 + n2 - 1;
    for( C = 0; C < Ny; C++ )
        y[C] = tmp1[C << 1];
    
    /* Free up temporary memory. */
    free( tmp1 );
    free( tmp2 );

    return Ny;
}

void IIRsos(
    FLOAT * x, UINT32 Nx,
    FLOAT b0, FLOAT b1, FLOAT b2, FLOAT a1, FLOAT a2,
    FLOAT * tz1, FLOAT * tz2 )
/************************************************************
*   IIRsos - fast second-order section IIR filter.
*   Filters Nx samples in x with IIR filter defined by
*   (b0 z^2 + b1 z + b2)/(z^2 + a1 z + a2).  Data is returned
*   in situ in x (sample for sample).  If tz1 or tz2 are NULL
*   the filter registers are initially set to zero.
*   Otherwise z1 is initially set to *tz1; after filtering
*   its value is stored back to *tz1. (Likewise z2 and *tz2.)
*   This allows the filter to be applied in stages.
*
*   Checked to give correct answers within rounding error
*   using Matlab.  22/10/1997 - Antony Rix.
************************************************************/
{
    /* Registers */
    register FLOAT z0;
    register FLOAT z1;
    register FLOAT z2;

    /* Incorporate initial conditions */
    if( tz1 == NULL ) z1 = 0.0f; else z1 = *tz1;
    if( tz2 == NULL ) z2 = 0.0f; else z2 = *tz2;
    
    /* Implement the IIR filter in situ, optimising for redundant taps */
    if( (a1 != 0.0f) || (a2 != 0.0f) )
    {
        /* Use the a taps */
        if( (b1 != 0.0f) || (b2 != 0.0f) )
        {
            /* Implement the full filter */
            while( (Nx) > 0 )
            {
                Nx--;
                                z0 = (*x) - a1 * z1 - a2 * z2;
                *(x++) = b0 * z0 + b1 * z1 + b2 * z2;
                z2 = z1;
                z1 = z0;
            }
        }
        else
        {
            if( b0 != 1.0f )
            {
                /* Omit taps b1 and b2 only */
                while( (Nx) > 0 )
                                {
                                        Nx--;
                    z0 = (*x) - a1 * z1 - a2 * z2;
                    *(x++) = b0 * z0;
                    z2 = z1;
                    z1 = z0;
                }
            }
            else
            {
                /* Omit all b taps */
                while( (Nx) > 0 )
                                {
                                        Nx--;
                    z0 = (*x) - a1 * z1 - a2 * z2;
                    *(x++) = z0;
                    z2 = z1;
                    z1 = z0;
                }
            }
        }
    }
    else
    {
        /* Omit the a taps */
        if( (b1 != 0.0f) || (b2 != 0.0f) )
        {
            /* Implement the filter without the a taps */
            while( (Nx) > 0 )
                        {
                                Nx--;
                z0 = (*x);
                *(x++) = b0 * z0 + b1 * z1 + b2 * z2;
                z2 = z1;
                z1 = z0;
            }
        }
        else
        {
            /* Omit taps a1, a2, b1 and b2 */
            if( b0 != 1.0f )
            {
                /* Implement the filter with b0 only */
                while( (Nx) > 0 )
                                {
                                        Nx--;
                    *x = b0 * (*x);
                    x++;
                }
            }
            /* Otherwise a1=a2=b1=b2=0 and b0=1 so there is nothing to do! */
        }
    }

    /*  */
    if( tz1 != NULL ) (*tz1) = z1;
    if( tz2 != NULL ) (*tz2) = z2;
}

void IIRFilt(
    FLOAT * h, UINT32 Nsos, FLOAT *z,
    FLOAT * x, UINT32 Nx, FLOAT *y )
/************************************************************
*   IIRFilt - cascaded second order sections IIR filter
*   h contains Nsos second order sections in order
*   b0 b1 b2 a1 a2.  2*Nsos initial registers are read from,
*   and stored to, z, if not set to NULL.
*   Nx output values are written to array y if it is not NULL,
*   otherwise Nx data points are processed in situ in x.
*
*   Checked to give correct answers within rounding error
*   using Matlab.  22/10/1997 - Antony Rix.
************************************************************/
{
    UINT32 C;

    if( y == NULL )
        y = x;
    else
    {
        /* Copy input data to output for processing in situ */
        for( C = 0; C < Nx; C++ )
            y[C] = x[C];
    }

    /* Apply second order sections */
    for( C = 0; C < Nsos; C++ )
    {
        if( z != NULL )
        {
            IIRsos( y, Nx, h[0], h[1], h[2], h[3], h[4], z, z+1 );
            z += 2;
        }
        else
            IIRsos( y, Nx, h[0], h[1], h[2], h[3], h[4], NULL, NULL );
        h += 5;
    }
}


#ifdef VTUNE
/***************************************************
*  RealFFT - compute the FFT of a real signal.
*  Inputs:
*   x contains real input data of length N samples
*   N must be a power of 2
*  Modifies:
*   On return, x contains (N/2)+1 complex pairs giving
*   the FFT of the signal in frequency bands 0 .. N/2.
*  Returns:
*   None
*
***************************************************/
void RealFFT(FLOAT * x, UINT32 N) 
{
	 int             log2N = 0;
    unsigned long   p = 1;

    while (p < N) {
        p <<= 1;
        log2N++;
    }

	if(sizeof(x[0]) == sizeof(double))
	{
		nspdRealFft ((double*)x, log2N, NSP_Forw | NSP_NoBitRev);
	}
	if(sizeof(x[0]) == sizeof(float))
	{
		nspsRealFft ((float *)x, log2N, NSP_Forw | NSP_NoBitRev);
	}
}

/***************************************************
*  RealIFFT - compute the inverse FFT of a real (time-domain) signal.
*  using the nsp library
*  Inputs:
*   x contains (N/2)+1 complex pairs giving
*   the FFT of the signal in frequency bands 0 .. N/2.
*   N must be a power of 2
*  Modifies:
*   On return, x contains the N samples of the IFFT of the input.
*  Returns:
*   None
*
***************************************************/
void RealIFFT(FLOAT * x, UINT32 N) 
{

	int              log2N = 0;
	unsigned long    p = 1;

    while (p < N) {
        p <<= 1;
        log2N++;
    }

	if(sizeof(x[0]) == sizeof(double))
	{
		nspdCcsFft ((double*)x, log2N, NSP_Inv | NSP_NoBitRev);
	}
	if(sizeof(x[0]) == sizeof(float))
	{
		 nspsCcsFft ((float*)x, log2N, NSP_Inv | NSP_NoBitRev);
	}



}


#else
 

/***************************************************
*  RealFFT - compute the FFT of a real signal.
*  using the nsp library
*  Inputs:
*   x contains real input data of length N samples
*   N must be a power of 2
*  Modifies:
*   On return, x contains (N/2)+1 complex pairs giving
*   the FFT of the signal in frequency bands 0 .. N/2.
*  Returns:
*   None
*
***************************************************/
void RealFFT(FLOAT * x, UINT32 N) {
    FLOAT * y = malloc(N*2*sizeof(FLOAT));
    UINT32 count;

    for( count = 0; count < N; count++ ) {
        y[count<<1] = x[count];
        y[(count<<1)+1] = 0.0f;
    }

    FFT(y, N);

    for( count = 0; count < (N+2); count++ )
        x[count] = y[count];

    free(y);
}

/***************************************************
*  RealIFFT - compute the inverse FFT of a real (time-domain) signal.
*  Inputs:
*   x contains (N/2)+1 complex pairs giving
*   the FFT of the signal in frequency bands 0 .. N/2.
*   N must be a power of 2
*  Modifies:
*   On return, x contains the N samples of the IFFT of the input.
*  Returns:
*   None
*
***************************************************/
void RealIFFT(FLOAT * x, UINT32 N) {
    FLOAT * y = malloc(N*2*sizeof(FLOAT));
    UINT32 count;

    for( count = 0; count <= (N>>1); count++ ) {
        y[count<<1] = x[count<<1];
        y[(count<<1)+1] = x[(count<<1)+1];
    }
    for( count = 1; count < (N>>1); count++ ) {
        y[(N-count)<<1] = x[count<<1];
        y[((N-count)<<1)+1] = -x[(count<<1)+1];
    }

    IFFT(y, N);

    for( count = 0; count < N; count++ )
        x[count] = y[count<<1];

    free(y);
}

#endif 
