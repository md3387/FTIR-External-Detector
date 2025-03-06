# FTIR-External-Detector Data conversion
Main Program File: FTIR_External_Detector_Data_Conversion.m

%Uses subfunctions:
%.....uipickfiles.m
%.....FindMirrorTrigger.m
%.....TDMS_readTDMSFile.m
%..........TDMS_handleGetDataOption
%..........TDMS_preprocessFile
%..........TDMS_getGroupChanNames
%..........TDMS_handleGetDataOption
%.....TDMS_readChannelOrGroup

%MDH 02/10/2019
%Converts a single voltage interferogram, mirror trigger, and HeNe Fringe
%into wavenumber vs transmission

%References
%1. Griffiths, P.R., deHaseth, J.A., "Fourier Transform Infrared Spectrometry 2nd Ed." Wiley 2007
%2. Thermo Fisher Scientific "FT-IR Spectrometers Nicolet iS50 Spectrometer User Guide" 2013
