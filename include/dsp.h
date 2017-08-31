/***********************************************************
*
*  DSP.H   General Digital Signal Processing functions.
*
*  Copyright (C) British Telecommunications plc, 1998.
*                Psytechnics Limited, 2003.
*  All rights reserved.
*
*  Author: Antony Rix    antony.rix@bt.com
*
************************************************************
*
*  Declares the functions in DSP.C.
*
***********************************************************/


#ifndef DSP_INCLUDED
  #define DSP_INCLUDED

#ifndef min
  #define min(a,b)  (((a) < (b)) ? (a) : (b))
#endif

#ifndef max
  #define max(a,b)  (((a) > (b)) ? (a) : (b))
#endif


  UINT32 nextpow2(UINT32 X);
  INT ispow2(UINT32 X);
  INT intlog2(UINT32 X);
  void FFTInit(UINT32 N);
  void FFTFree(void);
  void FFT(FLOAT * x, UINT32 N);
  void IFFT(FLOAT * x, UINT32 N);
  UINT32 FFTNXCorr(
    FLOAT * x1, UINT32 n1, FLOAT * x2, UINT32 n2, FLOAT * y );
  void IIRsos(
    FLOAT * x, UINT32 Nx,
    FLOAT b0, FLOAT b1, FLOAT b2, FLOAT a1, FLOAT a2,
    FLOAT * tz1, FLOAT * tz2 );
  void IIRFilt(
    FLOAT * h, UINT32 Nsos, FLOAT * z,
    FLOAT * x, UINT32 Nx, FLOAT * y );
	void RealFFT(FLOAT * x, UINT32 N);
	void RealIFFT(FLOAT * x, UINT32 N);



#endif

