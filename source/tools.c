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
#include <string.h>
#include <math.h>

#include "defines.h"
#include "dsp.h"
#include "tools.h"

/*********************
*
*	FUNCTION: find_max
*
*	DESCRIPTION:
*		Finds the maximum value in the input array and returns its position.
*
***********************/
UINT32 find_max(void *search_array, UINT32 start, UINT32 length, char in_type)
{

	UINT32 max_index ;
	INT16  max_value_INT16;
	INT32  max_value_INT32 ;
	FLOAT  max_value_float;

	UINT32 max_loop ;

	max_index = start ;

	switch(in_type)
	{
		case SHORT_TYPE:
			max_value_INT16  = ((INT16*)search_array) [start] ;
			for(max_loop = start; max_loop < length + start; max_loop++)
			{
				if (((INT16*)search_array)[max_loop] > max_value_INT16)
				{
					max_index = max_loop ;
					max_value_INT16 = ((INT16*)search_array)[max_loop] ;
				}
			}
			break ;

		case LONG_TYPE:
			max_value_INT32  = ((INT32*)search_array) [start] ;
			for(max_loop = start; max_loop < length + start; max_loop++)
			{
				if (((INT32*)search_array)[max_loop] > max_value_INT32)
				{
					max_index = max_loop ;
					max_value_INT32 = ((INT32*)search_array)[max_loop] ;
				}
			}
			break ;

		case FLOAT_TYPE:
			max_value_float  = ((FLOAT*)search_array) [start] ;
			for(max_loop = start; max_loop < length + start; max_loop++)
			{
				if (((FLOAT*)search_array)[max_loop] > max_value_float)
				{
					max_index = max_loop ;
					max_value_float = ((FLOAT*)search_array)[max_loop] ;
				}
			}
			break ;
	}
	return max_index;
}

/*********************
*
*	FUNCTION: multi_remove
*
*	DESCRIPTION:
*		Remove from list_array[] all indexed element from index_array[].
*		list_array[] values are overwritten by the new values.
*
***********************/
void multi_remove(INT32 *index, INT32 index_Nelements, INT32 *list_array, INT32 *Nelements)
{
	INT32 cnt,pos,cnt1;

	for (cnt=0,pos=0; cnt<index_Nelements; cnt++,pos++)
		for (cnt1=pos-cnt; pos<index[cnt]; cnt1++,pos++) list_array[cnt1] = list_array[pos];

	for (cnt1=pos-cnt; pos<(*Nelements); cnt1++,pos++) list_array[cnt1] = list_array[pos];

	*Nelements = cnt1;
}

/*********************
*
*	FUNCTION: bin_search
*
*	DESCRIPTION:
*		Embedded of bin_search() function.
*		Search for a specified value in a sorted array (use bin_query).
*
***********************/
INT32 bin_search( INT32 value, INT32 * list_array, INT32 Nelements )
{
    INT32 ptr = bin_query( value, list_array, Nelements );
    if( ptr == -1L ) return -1L;
    if( list_array[ptr] == value ) return ptr;
    return -1L;
}

/*********************
*
*	FUNCTION: bin_insert
*
*	DESCRIPTION:
*		Insert a value within a specific array (ordered).
*		If the value already exist, then it does not insert it.
*
***********************/
INT32 bin_insert( INT32 value, INT32 * list_array, INT32 * Nelements, INT32 * Nallocated )
{
    INT32 ptr;

    if( (list_array == NULL) || (Nelements == NULL) || (Nallocated == NULL) ) return -1L;
    if( *Nelements >= *Nallocated ) return -1L;

    if( *Nelements == 0 ) {
        list_array[0] = value;
        *Nelements = 1L;
        return 0L;
    }

    ptr = bin_query( value, list_array, *Nelements );
    if( ptr == -1L ) return -1L;
    if( list_array[ptr] == value ) return ptr;

    if( list_array[ptr] < value ) ptr++;

    if( ptr < *Nelements )
        memmove( (list_array + ptr + 1),
            (list_array + ptr),
            (*Nelements - ptr) * sizeof(INT32) );
    *Nelements += 1;

    list_array[ptr] = value;

    return ptr;
}

/*********************
*
*	FUNCTION: bin_remove
*
*	DESCRIPTION:
*		Remove a specific value in the input array.
*		Returns -1 if the value has not been found.
*
***********************/
INT32 bin_remove( INT32 value, INT32 * list_array, INT32 * Nelements )
{

    INT32 ptr ;

	ptr = bin_search( value, list_array, *Nelements );
    if( ptr == -1L ) return -1L;

    if( ptr < (*Nelements - 1) )
        memmove( (list_array + ptr),
            (list_array + ptr + 1),
            (*Nelements - ptr - 1) * sizeof(INT32) );
    *Nelements -= 1;

    return ptr;
}

/*********************
*
*	FUNCTION: cross_power_ratio
*
*	DESCRIPTION:
*		Compute the cross power correlation	and a ratio at the end.
*
***********************/
FLOAT cross_power_ratio(INT16 *x, INT32 x_Nelements, INT16 *y, INT32 y_Nelements)
{
	FLOAT *tempx;
	FLOAT *tempy;
	FLOAT a=0,b=0,tmp1,tmp2;
	INT32 count;

	tempx = zeropad(x,&x_Nelements);
	tempy = zeropad(y,&y_Nelements);

	if (x_Nelements>y_Nelements) x_Nelements = y_Nelements;

	RealFFT(tempx,x_Nelements);
	RealFFT(tempy,y_Nelements);

	for (count=0; count<(x_Nelements*0.15); count+=2)
	{
		tmp1 = tempx[count]*tempy[count]+tempx[count+1]*tempy[count+1];
		tmp2 = tempx[count]*tempy[count+1]-tempx[count+1]*tempy[count];
		a+= (FLOAT) sqrt(tmp1*tmp1 + tmp2*tmp2);
	}

	for (count=(INT32)(x_Nelements*0.15); count<x_Nelements; count+=2)
	{
		tmp1 = tempx[count]*tempy[count]+tempx[count+1]*tempy[count+1];
		tmp2 = tempx[count]*tempy[count+1]-tempx[count+1]*tempy[count];
		b+= (FLOAT) sqrt(tmp1*tmp1 + tmp2*tmp2);
	}

	free(tempx);
	free(tempy);

	if (b!=0) return a/b;
	else return 0;
}

/*********************
*
*	FUNCTION: reverse_array
*
*	DESCRIPTION:
*		Function to reverse the order of an array
*
***********************/
void reverse_array(void *array, INT32 array_length, char in_type)
{
	INT32 reverse_loop ;
	INT32 INT32_tmp ;
	INT16 INT16_tmp ;
	FLOAT double_tmp ;
	div_t div_result ;

	if(array_length == 0) return;
	div_result = div(array_length, 2) ;
	if(!div_result.rem) div_result.quot = div_result.quot - 1;

	switch(in_type)
	{
		case SHORT_TYPE:
			for(reverse_loop=0;reverse_loop<=div_result.quot; reverse_loop++)
			{
				INT16_tmp = ((INT16 *)array)[reverse_loop] ;
				((INT16 *)array)[reverse_loop] = ((INT16 *)array)[(array_length - 1) - reverse_loop] ;
				((INT16 *)array)[(array_length - 1) - reverse_loop] = INT16_tmp ;
			}
			break ;
		case LONG_TYPE:
			for(reverse_loop=0;reverse_loop<=div_result.quot; reverse_loop++)
			{
				INT32_tmp = ((INT32 *)array)[reverse_loop] ;
				((INT32 *)array)[reverse_loop] = ((INT32 *)array)[(array_length - 1) - reverse_loop] ;
				((INT32 *)array)[(array_length - 1) - reverse_loop] = INT32_tmp ;
			}
			break ;
		case FLOAT_TYPE:
			for(reverse_loop=0;reverse_loop<=div_result.quot; reverse_loop++)
			{
				double_tmp = ((FLOAT *)array)[reverse_loop] ;
				((FLOAT *)array)[reverse_loop] = ((FLOAT *)array)[(array_length - 1) - reverse_loop] ;
				((FLOAT *)array)[(array_length - 1) - reverse_loop] = double_tmp ;
			}
	}
}


/*********************
*
*	FUNCTION: search
*
*	DESCRIPTION:
*		Search within an array for a specific value.
*
***********************/
INT32 search(INT32 search_element, INT32 *search_array,
				INT32 search_array_length)
{
	INT32 search_loop ;

	for(search_loop = 0; search_loop < search_array_length; search_loop++)
		if(search_array[search_loop] == search_element)
			return (search_loop);

	return -1 ;
}

/*********************
*
*	FUNCTION: bin_query
*
*	DESCRIPTION:
*		Search a specific value in a sorted array
*
***********************/
INT32 bin_query( INT32 value, INT32 * list_array, INT32 Nelements )
{
    INT32 step_size;
    INT32 ptr;

    if( Nelements <= 0L ) return -1L;
    if( list_array == NULL ) return -1L;
    if( Nelements == 1L ) return 0L;

    step_size = 1L;
    while( step_size < Nelements ) step_size *= 2L;
    step_size /= 2L;

    ptr = step_size;
    while( step_size > 0L )
	{
        if( list_array[ptr] > value )
		{
            if( (ptr - step_size) >= 0L ) ptr -= step_size;
        }
		else if( list_array[ptr] < value )
		{
            if( (ptr + step_size) < Nelements ) ptr += step_size;
        }
		else return ptr;
		step_size /= 2L;
    }
    return ptr;
}

/*********************
*
*	FUNCTION: *zeropad
*
*	DESCRIPTION:
*		Pad the input array until its size is a power of 2.
*
***********************/
FLOAT *zeropad(INT16 *in, INT32 *in_Nelements)
{
	FLOAT *temp;
	INT32 cpt,cpt1;

	cpt = nextpow2(*in_Nelements);

	temp = (FLOAT *) malloc((cpt+2)*sizeof(FLOAT));
	for (cpt1=0;cpt1<*in_Nelements;cpt1++) temp[cpt1]=in[cpt1];
	for (cpt1=*in_Nelements;cpt1<cpt;cpt1++) temp[cpt1]=0;
	*in_Nelements = cpt;

	return temp;
}

/*********************
*
*	FUNCTION: find_min
*
*	DESCRIPTION:
*		Find the minimum value in the given array
*
***********************/
FLOAT find_min(const FLOAT * array, UINT32 search_length)
{
	UINT32 min_loop ;
	FLOAT tmp_min;

	tmp_min = array[0] ;

	for(min_loop = 0; min_loop < search_length; min_loop++)
		if(array[min_loop] < tmp_min) tmp_min = array[min_loop] ;

	return(tmp_min) ;

}

/*********************
*
*	FUNCTION: find_value
*
*	DESCRIPTION:
*		Search within an array for a specific value.
*
***********************/
INT32 find_value(const INT32 search_element, const INT32 *search_array, const INT32 search_array_length)
{
	INT32 search_loop ;

	for(search_loop = 0; search_loop < search_array_length; search_loop++)
		if(search_array[search_loop] == search_element) return search_loop;

	return (-1) ;
}

/*********************
*
*	FUNCTION: zero_removerf
*
*	DESCRIPTION:
*		Remove 0 values contained in the input array.
*		This function overwrite the inpute array and does not reallocate memory.
*
***********************/
void zero_removerf(FLOAT *in, INT32 *Nelements)
{
	INT32 i, count=0;
	for (i=0; i<*Nelements; i++)
	{
		if (in[i]) in[i-count] = in[i];
		else count++;
	}
	*Nelements -= count ;
}

/*********************
*
*	FUNCTION: smooth_averagesf
*
*	DESCRIPTION:
*		This function is computing a power average of the smooth section function
*		and the mean of this smooth section.
*
***********************/
FLOAT smooth_averagesf(FLOAT *in, INT32 in_Nelements, FLOAT threshold)

{
	INT32 i, smooth_Nelements;
	FLOAT *smooth;
	FLOAT output=0;

	smooth = smooth_sectionsf(in, in_Nelements, threshold, &smooth_Nelements);

	for (i = 0; i<smooth_Nelements; i++) output+= smooth[i]*smooth[i];
	output/=smooth_Nelements;
	output = (FLOAT) sqrt(output);

	free(smooth);

	return output;
}

/*********************
*
*	FUNCTION: *smooth_sectionsf
*
*	DESCRIPTION:
*		This function is deriving the input.
*		Smooth part are defined by sections without changes in the input.
*
***********************/
FLOAT *smooth_sectionsf(FLOAT *in, INT32 in_Nelements, FLOAT threshold, INT32 *out_Nelements)

{
	INT32 i, section_size=0, output_Nelements=0;
	FLOAT *deriv, *temp,*out;
	char in_smooth=0;

	temp  = (FLOAT*) malloc( in_Nelements * sizeof(FLOAT) );
	deriv = (FLOAT*) malloc( in_Nelements * sizeof(FLOAT) );

	for (i=0; i<in_Nelements; i++)
		if (in[i]>0) temp[i] = (FLOAT)((INT32) (100*in[i] + 0.5))/100;
		else temp[i] = (FLOAT)((INT32) (100*in[i] - 0.5))/100;

	derivative(temp,in_Nelements,deriv);

	for (i=0;i<in_Nelements;i++)
	{
		if ((-threshold<deriv[i]) && (deriv[i]<=threshold))
		{
			if (in_smooth) section_size++;
			else
			{
				section_size=1;
				in_smooth=1;
			}
		}
		else if (in_smooth)
		{
			deriv[output_Nelements] = (FLOAT) section_size;
			in_smooth = 0;
			output_Nelements++;
		}
	}
	out = (FLOAT *) malloc(output_Nelements * sizeof(FLOAT));

	for (i=0; i<output_Nelements; i++) out[i] = deriv[i];

	*out_Nelements = output_Nelements;

	free( deriv );
	free( temp );

	return out;
}

/*********************
*
*	FUNCTION: rms
*
*	DESCRIPTION:
*		Calculate the Root Mean Square value of the input array
*
***********************/
FLOAT rms(const void *in, const UINT32 Nelements, const int in_type)
{
	return (FLOAT) sqrt(rms_nosqrt(in, Nelements,in_type)/Nelements);
}

/*********************
*
*	FUNCTION: rms_nosqrt
*
*	DESCRIPTION:
*		Calculate the Mean Square value of the input array
*
***********************/
FLOAT rms_nosqrt(const void *in, const UINT32 Nelements, const int in_type)
{
	UINT32 lc_cnt=0;
	INT16 s_calc = 0;
	FLOAT f_calc = 0;
    FLOAT result=0;

	switch(in_type)
	{
	case SHORT_TYPE:
		for (lc_cnt=0; lc_cnt<Nelements; lc_cnt++)
		{
			s_calc = ((INT16*)in)[lc_cnt];
			result+=(FLOAT) s_calc * s_calc;
		}
		break;
	case FLOAT_TYPE:
		for (lc_cnt=0; lc_cnt<Nelements; lc_cnt++)
		{
			f_calc = ((FLOAT *)in)[lc_cnt];
			result+= f_calc * f_calc;
		}
    }
	return result;
}

/*********************
*
*	FUNCTION: derivative
*
*	DESCRIPTION:
*		Calculate the derivative on the input array
*
***********************/
void derivative(const void *in, const UINT32 in_Nelements, FLOAT *deriv)
{
	UINT32 cpt;

	if (in_Nelements>2)
	{
		deriv[0] = ((FLOAT *)in)[1]/2;
		deriv[in_Nelements-1] = ((FLOAT *)in)[in_Nelements-2]/2;
		for (cpt=1; cpt<(in_Nelements-1); cpt++)
			deriv[cpt]= (((FLOAT *)in)[cpt-1]-((FLOAT *)in)[cpt+1])/2;
	}
}

/*********************
*
*	FUNCTION: standard_deviation
*
*	DESCRIPTION:
*		Calculate the standard deviation of the input array
*
***********************/
FLOAT standard_deviation(const void *in, const UINT32 Nelements, const int in_type)
{
	FLOAT mean;
	return (FLOAT) sqrt(variance(in,Nelements,in_type,&mean));
}

/*********************
*
*	FUNCTION: variance
*
*	DESCRIPTION:
*		Calculate the variance of the input array
*
***********************/
FLOAT variance(const void *in, const UINT32 Nelements, const int in_type, FLOAT *Mean)
{
	UINT32 cpt=0;
	FLOAT  count = 0;
	FLOAT  m;
	FLOAT  result=0;

	m= mean(in,Nelements,in_type);
	switch (in_type)
	{
		case SHORT_TYPE:
			for (cpt=0; cpt<Nelements; cpt++)
			{
				count=( ((INT16*)in)[cpt] - m);
				result+=count*count;
			}
			break;
		case FLOAT_TYPE:
			for (cpt=0; cpt<Nelements; cpt++)
			{
				count=( ((FLOAT*)in)[cpt] - m);
				result+=count*count;
			}
			break;
	}
	*Mean = m;

	if (Nelements>1) return result/(Nelements-1);
	else return 0;
}

/*********************
*
*	FUNCTION: mean
*
*	DESCRIPTION:
*		Calculate the average of the input array
*
***********************/
FLOAT mean(const void *in, const UINT32 Nelements, const int in_type)
{
	UINT32 cpt=0;
	FLOAT  count = 0;

	switch (in_type)
	{
		case SHORT_TYPE:
			for (cpt=0; cpt<Nelements; cpt++) count+= ((INT16*)in)[cpt]; break;
		case FLOAT_TYPE:
			for (cpt=0; cpt<Nelements; cpt++) count+= ((FLOAT*)in)[cpt]; break;
		case LONG_TYPE:
			for (cpt=0; cpt<Nelements; cpt++) count+= ((INT32*)in)[cpt];
	}

	if (Nelements)	return (count/Nelements);
	else return 0;

}

/*********************
*
*	FUNCTION: round2
*
*	DESCRIPTION:
*		Round the input value to the closest integer value
*
*
***********************/
INT32 round2(FLOAT num)
{
	if (num<0)
	{
		if ((num - (INT32) num)<-0.5) return (INT32) floor(num);
		else return (INT32) ceil(num);
	}
	else
	{
		if ((num - (INT32) num)<=0.5) return (INT32) floor(num);
		else return (INT32) ceil(num);
	}
}

/*********************
*
*	FUNCTION: blackmanHarris_window
*
*	DESCRIPTION:
*		Calculate a Blackman Harris window
*
***********************/
int blackmanHarris_window(FLOAT *window, UINT32 window_Nelements)
{
	UINT32 loop = 0;
	const FLOAT cst1 = 0.42323f;
	const FLOAT cst2 = 0.49755f;
	const FLOAT cst3 = 0.07922f;

	FLOAT w = 0, wi = 0;

	if (window)
	{
		w = TWOPI / window_Nelements;

		for (loop=0; loop<window_Nelements; loop++)
		{
			wi = loop * w;
			window[loop] = (FLOAT) (cst1 - cst2 * cos (wi) + cst3 * cos (2*wi));
		}
		return 1;
	}
	else return 0;
}

/*********************
*
*	FUNCTION: cross_power
*
*	DESCRIPTION:
*		Calculate the cross power between the 2 input arrays.
*		If the size of both input arrays is not a power of two they are padded to the next power of 2.
*		If input arrays have different size, the minimum is kept as being the size of both arrays.
*
***********************/
void cross_power(INT16 *x, INT32 x_Nelements, INT16 *y, INT32 y_Nelements, FLOAT **mag, INT32 *mag_Nelements)
{
	FLOAT *tempx = NULL;
	FLOAT *tempy = NULL;
	FLOAT tmp1,tmp2;
	INT32 count,count1;

	tempx = zeropad(x,&x_Nelements);
	tempy = zeropad(y,&y_Nelements);

	if (x_Nelements>y_Nelements) x_Nelements = y_Nelements;
	RealFFT(tempy,y_Nelements);
	RealFFT(tempx,x_Nelements);

	if ((*mag)==NULL) (*mag) = (FLOAT*) malloc( (x_Nelements/2) * sizeof(FLOAT));

	for (count=0, count1=0; count<x_Nelements-1; count+=2,count1++)
	{
		tmp1 = tempx[count]*tempy[count]+tempx[count+1]*tempy[count+1];
		tmp2 = tempx[count]*tempy[count+1]-tempx[count+1]*tempy[count];
		(*mag)[count1] = (FLOAT) sqrt(tmp1*tmp1 + tmp2*tmp2);

		(*mag)[count1]/= (x_Nelements*x_Nelements);
	}

	*mag_Nelements = count1;

	free(tempx);
	free(tempy);

}

/*********************
*
*	FUNCTION: *cross_correl
*
*	DESCRIPTION:
*		Calculate the cross correlation between both arrays.
*
***********************/
FLOAT *cross_correl(INT16 *x, INT32 x_Nelements, INT16 *y, INT32 y_Nelements, INT32 *correl_Nelements)
  {
	INT32 i;
	FLOAT *x1 = NULL, *y1 = NULL;
	FLOAT *correl;

	x1 = (FLOAT *) malloc ( x_Nelements * sizeof(FLOAT) );
	y1 = (FLOAT *) malloc ( y_Nelements * sizeof(FLOAT) );

	for (i=0; i<x_Nelements;i++) x1[i] = x[i];
	for (i=0; i<y_Nelements;i++) y1[i] = y[i];

	*correl_Nelements = x_Nelements + y_Nelements - 1;

	correl = (FLOAT *) malloc( (*correl_Nelements+1) * sizeof(FLOAT) );

	FFTNXCorr(x1, x_Nelements, y1, y_Nelements, correl);

	free(x1 );
	free(y1 );

	return correl;
}
