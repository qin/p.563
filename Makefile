#********************************************************************
# ITU-T Draft Recommendation P.563
# Version 1.0 - 23 March 2004
#
# NOTICE
#
# The Single Ended Assessment Model P.563 algorithm and the copyright therein
# is the joint property of Psytechnics Limited, OPTICOM GmbH and SwissQual AG
# and is protected by UK, US and other patents, either applied for or
# registered.
# Permission is granted to use this source code solely for the purpose of
# evaluation of ITU-T recommendation P.563.
# Any other use of this software requires a licence, which may be obtained
# from:
#
# OPTICOM GmbH
# Am Weichselgarten 7, D- 91058 Erlangen, Germany
# Phone: +49 9131 691 160   Fax: +49 9131 691 325
# E-mail: info@opticom.de         www.3sqm.com
#
# Psytechnics Limited
# Fraser House, 23 Museum Street, Ipswich, IP1 1HN, UK
# Phone: +44 1 473 261 800  Fax: +44 1 473 261 880
# E-mail: info@psytechnics.com    www.psytechnics.com
#
# SwissQual AG
# Gewerbestrasse 2 CH-4528 Zuchwil, Switzerland
# Phone: +41 32 685 08 30   Fax: +41 32 685 08 31
# E-mail: sales@swissqual.com     www.swissqual.com
#
# Psytechnics, SwissQual or Opticom can provide licences and further
# information.
#
# Authors:
#      Ludovic Malfait ludovic.malfait@psytechnics.com
#      Roland Bitto rb@opticom.de
#      Pero Juric pero.juric@swissqual.com
#
#********************************************************************/




########################################################################
#for the executable creation
#DEBUG_FLAG= -g

CC=gcc
LINK=gcc

DEBUG=  -Wall
RELEASE= -O3 -s

DEFINES= -D_LINUX_ -ffloat-store -funroll-loops
INCLUDE=-I . -Iinclude
VPATH	= .:source
OBJDIR=./obj


##################################################################
# Target specifications
##################################################################

ifeq ($(DEBUG_FLAG), -g)
   CFLAGS=  $(DEBUG) $(DEFINES) $(INCLUDE) -c
   OUTDIR=Debug
else
   CFLAGS=  $(RELEASE) $(DEFINES) $(INCLUDE) -c
   OUTDIR=./bin
endif

OBJS_P563=\
	$(OBJDIR)/back_noise.o \
	$(OBJDIR)/beeprob.o \
	$(OBJDIR)/dsp.o \
	$(OBJDIR)/Enhance.o \
	$(OBJDIR)/EvalQual.o \
	$(OBJDIR)/hosm.o \
	$(OBJDIR)/inter_detect.o \
	$(OBJDIR)/LpcAnalysis.o \
	$(OBJDIR)/lpc.o \
	$(OBJDIR)/mapping.o \
	$(OBJDIR)/module1.o \
	$(OBJDIR)/module2.o \
	$(OBJDIR)/module3.o \
	$(OBJDIR)/tools1.o \
	$(OBJDIR)/p563.o \
	$(OBJDIR)/pitch.o \
	$(OBJDIR)/Quant.o \
	$(OBJDIR)/SignalsPercept.o \
	$(OBJDIR)/SpeechLib.o \
	$(OBJDIR)/Statistics.o \
	$(OBJDIR)/tools.o \
	$(OBJDIR)/vector_lib.o \


#################################################################
# compile files
#################################################################

$(OBJDIR)/%.o : %.c
	$(CC) $(CFLAGS) $< -o $@


##################################################################
# targets
##################################################################
binaries=\
	$(OUTDIR)/p563

.PHONY: prepare clean
all: $(binaries)


all: $(binaries)

clean:
	rm -f $(OBJDIR)/*.o
	rm -f $(OUTDIR)/P563

prepare:
	mkdir -p $(OUTDIR)
	mkdir -p $(OBJDIR)

install:
	cp $(OUTDIR)/p563	/usr/local/bin

uninstall:
	rm -f /usr/local/bin/p563


$(OUTDIR)/p563:$(OUTDIR) $(OBJS_P563)
	@echo Making p563
	$(LINK) $(OBJS_P563) -lm -o $(OUTDIR)/p563
