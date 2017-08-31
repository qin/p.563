
#ifndef PSEAM_DEFINITIONS
	#define PSEAM_DEFINITIONS


/*******************************
*     Constant definitions     *
*******************************/

	/* General use definitions */
	#define MIN_SIGNAL_LEN		40000		/* 5 seconds at 8kHz */
	#define MAX_SIGNAL_LEN		160000		/* 20 seconds at 8kHz */
	#define INVALID_CODE		-999
	#define SAMPLE_FREQUENCY	8000
	#define TARGET_SPEECH_LEVEL	-26

	#define PI		3.141592653589793f
	#define TWOPI	6.283185307179586f

	#define VAD_FRAME_LEN				32		/* Frame length for VAD calculation */
	#define VAD_MIN_SPEECH_SECTION_LEN	4		/* Mimimum length for speech sections */
	#define VAD_JOIN_SPEECH_SECTION_LEN	50		/* Miminum length between 2 consecutives speech sections */
	#define REF_DBOV					90.3f	/* Maximum level (dB) an INT16 signal could reach */

	#ifndef FALSE
		#define     FALSE                           0
	#endif

	#ifndef TRUE
		#define     TRUE                            1
	#endif

/*******************************
*        Type definitions      *
*******************************/

	typedef char				INT8;         
	typedef unsigned char	UINT8; 
	
	typedef short				INT16;
	typedef unsigned short	UINT16;
	typedef long				INT32;
	typedef unsigned long	UINT32;

	typedef int					INT;
	typedef unsigned int		UINT;

	typedef char bool ;


#define FLOAT_64  1
#if FLOAT_64
	typedef  double		FLOAT;
#else
	typedef  float			FLOAT;
#endif

	typedef short				BOOL;


/*******************************
*    Structure definitions     *
*******************************/

	/* Structure holding preprocessed signal and VAD output */
	typedef struct
	{
		FLOAT	*pfSignal;
		FLOAT	*pfSignal26dB;
		FLOAT	*pfSignal26dB_IRS;
		INT32	lSignalLen;
		FLOAT	fLevelNorm;
		FLOAT	*pfVadProfile;
		INT32	lVadProfileLen;
	} SignalInfo_struct;

	/* Structure storing basic speech descriptors */
	typedef struct
	{
		FLOAT	fSpeechLevel;
		FLOAT	fPitchAverage;
		FLOAT	fSpeechSectionLevelVar;
	} basic_desc_struct;

	/* Structure storing noise analysis parameters */
	typedef struct
	{

		FLOAT fEstBGNoise;
      FLOAT fEstSegSNR;
      FLOAT fSpecLevelDev;
      FLOAT fSpecLevelRange;
      FLOAT fRelNoiseFloor;



		FLOAT	fNoiseLevel;
		FLOAT	fSnr;
		FLOAT	fHiFreqVar;
		FLOAT	fSpectralClarity;

		FLOAT fLocalBGNoise;
		FLOAT fLocalBGNoiseMean;
		FLOAT fLocalBGNoiseLog;
		FLOAT fLocalMeanDistSamp; 

		FLOAT fGlobalBGNoise;	

	} noise_analysis_struct;

	/* Structure storing mute parameters */
	typedef struct
	{
		FLOAT	fMuteLength;
		FLOAT fSpeechInterruptions;
	   FLOAT fSharpDeclines;
		FLOAT fUnnaturalSilenceMean;


	} mutes_struct;

	/* Structure storing unnatural speech parameters */
	typedef struct
	{
		FLOAT fConsistentArtTracker;
		FLOAT fVtpMaxTubeSection;
		FLOAT fFinalVtpAverage;
		FLOAT fVtpPeakTracker;
		FLOAT fArtAverage;
		FLOAT fVtpVadOverlap;
		FLOAT fPitchCrossCorrelOffset;
		FLOAT fPitchCrossPower;

		FLOAT fFrameRepeats;
		FLOAT fFrameRepeatsMean;     
		FLOAT fUBeeps;	
		FLOAT fUBeepsMean;		
		FLOAT fUnBeepsMeanDistSamp;  
		FLOAT fRobotisation;

		FLOAT fCepADev;			   
		FLOAT fCepSkew;			           
		FLOAT fCepCurt;			           
		FLOAT fLPCCurt;			           
		FLOAT fLPCSkew;			           
		FLOAT fLPCSkewAbs;	

	} unnatural_struct;

	/* Structure storing speech extraction parameters */
	typedef struct
	{

		FLOAT fBasicVoiceQualityAsym;
		FLOAT fBasicVoiceQuality;

	} speech_extract_struct;

	/* Structure grouping all parameters */
	typedef struct
	{
		mutes_struct			tMutes;
		noise_analysis_struct	tNoise;
		unnatural_struct		tUnnatural;
		basic_desc_struct		tBasicDesc;
		speech_extract_struct	tSpeechExtract;

		INT32 lPartition;
		FLOAT fPredictedMos;

	} p563Results_struct;

#endif
