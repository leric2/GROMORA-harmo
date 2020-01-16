function retrievalTool = import_default_retrievalTool(instrumentName)
%==========================================================================
% NAME          | import_default_retrievalTool(instrumentName)
% TYPE          | Function
% AUTHOR(S)     | Eric Sauvageat
% CREATION      | 01.2020
%               |
% ABSTRACT      | 
%               | 
%               |
%               |
% ARGUMENTS     | INPUTS: * instrumentName: name of the instrument as a
%               | string (ex: 'GROMOS')
%               |
%               | OUTPUTS: * retrievalTool: the default toolbox for
%               | launching a retrieval for this instrument.
%               | 
%               |
% CALLS         |
%               |
%               |
%               |

%==========================================================================
retrievalTool = struct();

% Name of the instrument
retrievalTool.instrumentName=instrumentName;

% Valid properties for all instruments
retrievalTool.bytesPerValue=4;
retrievalTool.binaryType='ieee-be';
retrievalTool.indiceCold=0;
retrievalTool.indiceAntenna=1;
retrievalTool.indiceHot=2;

retrievalTool.hotTemperatureStdThreshold=10;

% minimum number of indices (h-a-c) we want in a calibration cycle for it
% to be valid
retrievalTool.minNumberOfIndicePerCycle=15;

retrievalTool.numberOfAquisitionSpectraHot=60000;
retrievalTool.numberOfAquisitionSpectraAntenna=120000;
retrievalTool.numberOfAquisitionSpectraCold=120000;

retrievalTool.systemTempMaxStd=10;

switch instrumentName
    case 'GROMOS'
        retrievalTool.dataLocation='BERN';
        retrievalTool.PI_NAME='Murk;Axel';
        retrievalTool.PI_AFFILIATION='Universtiy of Bern;UBERN';
        retrievalTool.PI_ADDRESS='Sidlerstrasse 5, University of Bern;3012 Bern;Switzerland';
        retrievalTool.PI_EMAIL='axel.murk@iap.unibe.ch'; 
        retrievalTool.dataSource='MWR.O3_UBERN';
        
        retrievalTool.numberOfChannels=32768;
        
        retrievalTool.lon=7.44;
        retrievalTool.lat=46.95;
        retrievalTool.altitude=560;
        retrievalTool.numberOfTippingCurveExpected=48;
        retrievalTool.toleranceTippingCurves=2;
        retrievalTool.elevationAngleAntenna=40;
        retrievalTool.elevationAngleCold=-84;
        retrievalTool.elevationAngleHot=160;
        retrievalTool.elevationAngleTolerance=5;
        % Considering the expected number of tipping curve:
        retrievalTool.numberOfCyclesExpected=1500;
        retrievalTool.toleranceNumberCycles=15;
        retrievalTool.tippingSize=27;
        retrievalTool.flipped_spectra=true;
        retrievalTool.flip_spectra=@(rawSpectra) flip_spectra_gromos(rawSpectra);
        retrievalTool.THotUnit='degreeC';
        % Harmonizing GROMOS log (for units)
        retrievalTool.harmonize_log=@(log) harmonize_log_gromos(log);  
        
        % Calibration outlier management
        retrievalTool.threshNumRawSpectraHot=0.05*retrievalTool.numberOfChannels;
        retrievalTool.threshNumRawSpectraCold=0.05*retrievalTool.numberOfChannels;
    case 'SOMORA'
        retrievalTool.dataLocation='PAYERNE';
        retrievalTool.PI_NAME='';
        retrievalTool.PI_AFFILIATION='Swiss Meteorological Institute;MCH';
        retrievalTool.PI_ADDRESS='';
        retrievalTool.PI_EMAIL=''; 
        retrievalTool.dataSource='MWR.O3_MCH';
        
        retrievalTool.lon=6.95;
        retrievalTool.lat=46.82;
        retrievalTool.altitude=491;
        retrievalTool.numberOfChannels=16384;
        retrievalTool.numberOfTippingCurveExpected=48;
        retrievalTool.toleranceTippingCurves=2;
        retrievalTool.tippingSize=5;
        retrievalTool.elevationAngleAntenna=38;
        retrievalTool.elevationAngleCold=-90;
        retrievalTool.elevationAngleHot=180;
        retrievalTool.elevationAngleTolerance=5;
        retrievalTool.numberOfCyclesExpected=NaN;
        retrievalTool.toleranceNumberCycles=NaN;
        retrievalTool.flipped_spectra=false;
        % The log needs to be harmonized, for now taking gromos as basis
        retrievalTool.harmonize_log=@(log) harmonize_log_somora(log);
        retrievalTool.THotUnit='K';
        
        % Calibration outlier management
        retrievalTool.threshNumRawSpectraHot=0.05*retrievalTool.numberOfChannels;
        retrievalTool.threshNumRawSpectraCold=0.05*retrievalTool.numberOfChannels;
end

