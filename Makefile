LIBRARY := SUHH2ZPrimeTotTPrime
DICT := include/ZPrimeTotTPrimeReconstructionHypothesis.h include/SUHH2ZPrime_LinkDef.h
USERLDFLAGS := -lSUHH2core -lSUHH2common -lGenVector
# enable par creation; this is necessary for all packages containing AnalysisModules
# to be loaded from by AnalysisModuleRunner.
PAR := 1
include ../Makefile.common
