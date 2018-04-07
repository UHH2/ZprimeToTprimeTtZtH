LIBRARY := SUHH2ZprimeToTprimeTtZtH
DICT := include/ZPrimeTotTPrimeReconstructionHypothesis.h include/SUHH2ZPrime_LinkDef.h
USERLDFLAGS := -lSUHH2core -lSUHH2common -lGenVector -L/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/lhapdf/6.1.6-ikhhed/lib -lLHAPDF 
USERCXXFLAGS := -I/cvmfs/cms.cern.ch/slc6_amd64_gcc530/external/lhapdf/6.1.6-ikhhed/include/ 
# enable par creation; this is necessary for all packages containing AnalysisModules
# to be loaded from by AnalysisModuleRunner.
PAR := 1
include ../Makefile.common
