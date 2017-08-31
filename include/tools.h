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

#ifndef TOOLS_H
	#define TOOLS_H

	#define SHORT_TYPE 1
    #define LONG_TYPE  2
    #define FLOAT_TYPE 4

	UINT32 find_max(void *search_array, UINT32 start, UINT32 length, char in_type);
	void multi_remove(INT32 *index, INT32 index_Nelements, INT32 *list_array, INT32 *Nelements);
	INT32 bin_search( INT32 value, INT32 * list_array, INT32 Nelements );
	INT32 bin_insert( INT32 value, INT32 * list_array, INT32 * Nelements, INT32 * Nallocated ) ;
	INT32 bin_remove( INT32 value, INT32 * list_array, INT32 * Nelements );
	INT32 search(INT32 search_element, INT32  *search_array, INT32 search_array_length);
	FLOAT cross_power_ratio(INT16 *x, INT32 x_Nelements, INT16 *y, INT32 y_Nelements);
	void reverse_array(void *array, INT32 array_length, char in_type);
	INT32 search(INT32 search_element, INT32 *search_array, INT32 search_array_length);
	INT32 bin_query( INT32 value, INT32 * list_array, INT32 Nelements );
	FLOAT *zeropad(INT16 *in, INT32 *in_Nelements);
	FLOAT find_min(const FLOAT * array, UINT32 search_length);
	INT32 find_value(const INT32 search_element, const INT32 *search_array, const INT32 search_array_length);
	void zero_removerf(FLOAT *in, INT32 *Nelements);
	FLOAT smooth_averagesf(FLOAT *in, INT32 in_Nelements, FLOAT threshold);
	FLOAT *smooth_sectionsf(FLOAT *in, INT32 in_Nelements, FLOAT threshold, INT32 *out_Nelements);

	FLOAT rms_nosqrt(const void *in, const UINT32 Nelements, const int in_type);
	FLOAT rms(const void *in, const UINT32 Nelements, const int in_type);
	void derivative(const void *in, const UINT32 in_Nelements, FLOAT *deriv);
	FLOAT standard_deviation(const void *in, const UINT32 Nelements, const int in_type);
	FLOAT variance(const void *in, const UINT32 Nelements, const int in_type, FLOAT *Mean);
	FLOAT mean(const void *in, const UINT32 Nelements, const int in_type);
	INT32 round2(FLOAT num);
	int blackmanHarris_window(FLOAT *window, UINT32 window_Nelements);
	void cross_power(INT16 *x, INT32 x_Nelements, INT16 *y, INT32 y_Nelements, FLOAT **mag, INT32 *mag_Nelements);
	FLOAT *cross_correl(INT16 *x, INT32 x_Nelements, INT16 *y, INT32 y_Nelements, INT32 *correl_Nelements);
#endif
