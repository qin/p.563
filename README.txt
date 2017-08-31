--------------------------------------------------------------
*** README file for ITU-T P.563 Version 1.0 software
*** Reference implementation of P.563 (Single Ended 
*** Assessment Model) and conformance data
--------------------------------------------------------------

TITLE
-----
 ITU-T P.563 (05/2004) - Single-ended method for objective
                speech quality assessment in narrow-band
                telephony applications
                
This electronic material is provided as part of ITU-T Rec. P.563.
The package contains the following files:

SOURCE CODE:
------------

README.TXT   				This file

C files
=======
	back_noise.c	  		Local and Global Background Noise      
	beeprob.c 				Evaluation of UnnaturalBeeps Robotization UnnaturalSilences
                            SharpDeclines and FrameRepeats
	dsp.c					Signal processing functions
	Enhance.c				Speech enhancement modules
	EvalQual.c				Basic Quality evaluation 
	hosm.c					Speech statistics functions
	inter_detect.c			Signal interruption functions
	lpc.c					LPC calculation functions
	LpcAnalysis.c			Linear Prediction modules for Speech enhancement
	mapping.c				Perceptual mapping
	module1.c				Pitch, vocal tract, noise and mutes analysis functions
	module2.c				Unnatural speech, noise and mutes analysis functions
	module3.c				Speech statistics interruptions and segmental SNR functions 
	tools1.c				Basic Arithmetic functions
	p563.c					Main program
	pitch.c					Pitch extraction functions
	Quant.c					Speech enhancement modules
	SignalsPercept.c		Basic Voice quality functions
	SpeechLib.c				Basic filter functions
	Statistics.c			Basic Statistic functions
	tools.c					Arithmetic functions
	vector_lib.c			Basic vector arithmetic functions 

Header files
============
	back_noise.h			Header file for back_noise.c
	beeprob.h				Header file for beeprob.c
	defines.h				P.563 model definitions
	dsp.h					Header file for dsp.c
	Enhance.h				Header file for Enhance.c
	EvalQual.h				Header file for EvalQual.c
	generic_typedefs.h
	hosm.h					Header file for hosm.c
	interr_detect.h			Header file for inter_detect.c
	lpc.h					Header file for lpc.c
	LpcAnalysis.h			Header file for LpcAnalysis.c
	mapping.h				Header file for mapping.c
	module1.h				Header file for module1.c
	module2.h				Header file for module2.c
	module3.h				Header file for module3.c
	tools1.h				Header file for tools1.c
	pitch.h					Header file for pitch.c
	Quant.h					Header file for Quant.c
	QuantTab.h		 		Constants for Quant.c
	resource.h
	SignalsPercept.h		Header file for SignalsPercept.c
	SpeechLib.h				Header file for SpeechLib.c
	Statistics.h			Header file for Statistics.c
	tools.h					Header file for tools.c
	vector_lib.h			Header file for vector_lib.c


COMPILING P.563:
-----------------
1. Windows platform
	- cd to the distribution directory
	- nmake /f make_win prepare
	- nmake /f make_win clean
	- nmake /f make_win

2. Linux platform
	- cd to the distribution directory
	- gmake prepare
	- gmake clean
	- gmake

Executables are generated in bin directory.
	

CONFORMANCE DATA:
-----------------
conform\conformance_results_supp23.txt 	Conf. data, ITU-T P-series Suppl.23
conform\process.bat						Batch script


COPYRIGHT AND INTELLECTUAL PROPERTY
-----------------------------------

AUTHOR'S NOTICE
---------------
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
  
  
ITU-T NOTICE
------------

Notwithstanding the above, it is prohibited to place this software
on a network or to copy the software for sale, license or transfer
to third parties.

The ITU does not take any responsibility for difficulties or the 
impossibility of reading out the contents of these files. The ITU 
shall not be held liable for any direct, indirect, consequential or 
incidental damages arising out of the use of or inability to use the 
software.

The present notice must be included in any copy of the software or
the information contained therein.


ITU CONTACTS
------------

For distribution of update software, please contact:
Sales Department
ITU
Place des Nations
CH-1211 Geneve 20
SUISSE
email: sales@itu.int

For reporting problems, please contact TSB helpdesk service at:
TSB Helpdesk service
ITU
Place des Nations
CH-1211 Geneve 20
SUISSE
fax: +41 22 730 5853
email: tsbedh@itu.int
