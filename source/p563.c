
/********************************************************/
/*             P.563 Implementation code                */
/********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "defines.h"
#include "module1.h"
#include "module2.h"
#include "module3.h"
#include "mapping.h"

static void usagepseam (void);
static INT16 PrintResults(char *pcProcessFilename,p563Results_struct *ptResults,char *pcResultFilename,char *pcAdditionalInfo);
static int PostProcessMovs(int *PartitionNumber, FLOAT *PredictedMos, 	p563Results_struct *ptResults);

/********************************************************************************************************/
/* Syntax: arg[1]: speech file (no header 16 bit, 8Khz, mono) */
/*	       arg[2]: output file */
/*	       arg[3]: optional parameters passed through to the result file (enclosed in qotes) */
/*	       arg[4]: optional names of parameters passed through to the result file (enclosed in qotes) */

int main( int argc, char * argv[] )
{
	p563Results_struct tResults = {-1}; 

	INT32	i;
	short int * speech_samples;     /* Array of speech samples                               */
	INT32 file_length;           /* Number of samples in speech file                      */
	INT32 num_read;				/* number of elements read from file                     */

	char* PassThroughParameters= "";
	char pSpeechFilename[100];
	char *pResultFileName="stdout";



	FILE *speech_file;

	bool RefFileSet=FALSE;
	bool PassTitlesSet = FALSE;



	int PartitionNumber;
	FLOAT PredictedMos;


	/* Open speech file and results file               */
	
/*	_controlfp(_EM_INEXACT|_EM_DENORMAL|_EM_UNDERFLOW,_MCW_EM); */
 	if (argc < 2)	{
		usagepseam();
		exit(-1);
	}
	
	strcpy(pSpeechFilename,argv[1]);
	RefFileSet=TRUE;

	for (i=2; i<argc; i++)
	{
		if (!strcmp(argv[i], "-out"))
		{
			i++;
			if (i<argc) pResultFileName = argv[i];
		}
		else 
		{
			if      (!RefFileSet)    { strcpy(pSpeechFilename,argv[i]); RefFileSet=TRUE;}
			else if (!PassTitlesSet) { PassThroughParameters = argv[i]; PassTitlesSet = TRUE;}
		}
	}

	/* open speech file */
	if((speech_file = fopen(pSpeechFilename, "rb")) == NULL )
	{
		printf( "Cannot open file %s\n", pSpeechFilename ) ;
		exit (-1) ;
	}

	/* Read speech samples from file */
	fseek(speech_file, 0 , SEEK_END) ;
	file_length = (ftell(speech_file))/(sizeof(short int)) ;
	fseek(speech_file, 0 , SEEK_SET) ;

	speech_samples  = (short int *)calloc(file_length, sizeof(short int)) ;
	num_read = fread(speech_samples, sizeof(short int), file_length, speech_file);

	if(num_read != file_length)
	{
	    printf( "Unable to read all the data from file\n" );
	    exit(-1);
	}
	
/*************************
*     PROCESS MODEL      *
*************************/

	module1(speech_samples, file_length, &tResults);

	module2(speech_samples, file_length, &tResults);

	module3(speech_samples, file_length, &tResults);
		
    PostProcessMovs( &PartitionNumber, &PredictedMos, &tResults);

    PrintResults(pSpeechFilename,&tResults,pResultFileName,PassThroughParameters);

	/* Clean up*/
	fclose(speech_file) ;
	free(speech_samples);

	return(0);
}


static INT16 PrintResults(char *pcProcessFilename,p563Results_struct *ptResults,char *pcResultFilename,char *pcAdditionalInfo)
{
	FILE *hFile=NULL;
	INT32 lWriteHeaderFlag=0;


	if (pcResultFilename!=NULL)
	{
		if(!strcmp(pcResultFilename,"stdout"))
		{
			hFile=stdout;
			lWriteHeaderFlag=1;
		}
		else
		{
			if((hFile = fopen(pcResultFilename, "r")) == NULL )
			{
				lWriteHeaderFlag=1;
			}
			else fclose(hFile);
			if ( (hFile = fopen(pcResultFilename, "a")) == NULL) /* If the file cannot be found */
			{
				printf("\n'%s' cannot be opened for writing!\n\n",pcResultFilename); 
				return -1;
			}
		}

		if (lWriteHeaderFlag)
		{
			fprintf(hFile,"Filename\t");
			fprintf(hFile,"MOS\n");
		}
		fprintf(hFile, "%s\t",pcProcessFilename);

		fprintf(hFile,"%f\t",ptResults->fPredictedMos);
	
		if (pcAdditionalInfo!=NULL) fprintf(hFile, pcAdditionalInfo);
		fprintf(hFile,"\n");
		if(strcmp(pcResultFilename,"stdout"))
		{
			fclose(hFile);
		}
	}
	return 0;
}


static void usagepseam (void)
{
	printf ("Usage:  p563 <SpeechFile> [-out <ResultFile>] <Parameters passed through to resultfile>\n");
}



int PostProcessMovs(int *PartitionNumber, FLOAT *PredictedMos, 	p563Results_struct *ptResults)
{

	GetPartitionNumber( PartitionNumber, ptResults);
	CalculateOverallMapping( PartitionNumber,ptResults, PredictedMos);
	
	ptResults->lPartition=*PartitionNumber;
	ptResults->fPredictedMos=*PredictedMos;

	return 1;
}
