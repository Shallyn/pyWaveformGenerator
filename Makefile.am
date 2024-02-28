# Makefile.am
AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I m4
#pkgpyexecdir = ${HOME}/usr/local

lib_LTLIBRARIES = libEOB.la
#bin_PROGRAMS = MAIN

libEOB_la_SOURCES = \
myUtils.c mySpherical.c myFileIO.c myOptparser.c myLog.c myAlloc.c \
pHamiltonian.c pRK4pdeIntegrator.c pEvolve.c pUtils.c \
pFactorizedWaveform.c pCore.c pInitialCondition.c \
pEnergyFlux.c pNQCcorrection.c pBHRingdown.c newFactorizedWaveform.c \
newFactorizedWaveformPrec.c pInitialConditionExact.c \
pPrec.c pPrecHam.c pPrecUtils.c pPrecRRForce.c pPrecWaveform.c pDebug.c interface.c


libEOB_la_LIBADD = -lgsl -lgslcblas #-lgthread-2.0
libEOB_la_LDFLAGS = -module -avoid-version
#MAIN_SOURCES = main.c
#MAIN_LDADD = libEOB.la

#all:$(bin_PROGRAMS)
#$(bin_PROGRAMS): $(lib_LTLIBRARIES) 	
all:$(lib_LTLIBRARIES) 