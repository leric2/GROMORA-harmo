function calibrationTool = import_spectrometer_mopi5(calibrationTool, modelFFTS)
%==========================================================================
% NAME          | import_default_calibrationTool(instrumentName)
% TYPE          | Function
% AUTHOR(S)     | Eric Sauvageat
% CREATION      | 01.2020
%               |
% ABSTRACT      | 
%               | 
%               |
%               |
% ARGUMENTS     | INPUTS: name of the instrument as a
%               | string (ex: 'GROMOS')
%               |         - calibrationTool: the default toolbox
%               |
%               | OUTPUTS: - calibrationTool: the default toolbox for
%               | launching a retrieval for this instrument.
%               | 
%               |
% CALLS         |
%               |
%               |
%               |

%==========================================================================
calibrationTool.ffts_model = modelFFTS;
spectrometerTypes  = {'USRP-A', 'USRP-B','U5303', 'AC240'};
samplingRateFFTS = [200 20  3200 2000]; % sampling rates in MHz 
BW = [200e6 20e6 1.6e9 1e9];

calibrationTool.numberOfChannels=16384;

if modelFFTS<3
    calibrationTool.IQProcessing = true;
    calibrationTool.LOFreqTot = calibrationTool.LOFreq1 + calibrationTool.LOFreq3;
    calibrationTool.badChannels = horzcat([1:1024],[8193],[calibrationTool.numberOfChannels-1024:calibrationTool.numberOfChannels]);
elseif modelFFTS == 3
    calibrationTool.IQProcessing = false;
    calibrationTool.LOFreqTot = calibrationTool.LOFreq1 + calibrationTool.LOFreq2;
    calibrationTool.badChannels = horzcat([1:64],[11000:calibrationTool.numberOfChannels]);
elseif modelFFTS == 4
    calibrationTool.IQProcessing = false;
    calibrationTool.LOFreqTot = calibrationTool.LOFreq1 + calibrationTool.LOFreq2;
    calibrationTool.badChannels = horzcat([1:64],[calibrationTool.numberOfChannels-64:calibrationTool.numberOfChannels]);
end


calibrationTool.instrumentBandwidth = BW(modelFFTS);
calibrationTool.spectrometer=spectrometerTypes{modelFFTS};
calibrationTool.samplingRateFFTS=samplingRateFFTS(modelFFTS);
calibrationTool.filenameLevel1a=['/scratch/MOPI5/Level1/mopi5_level1a_' calibrationTool.spectrometer '_' calibrationTool.dateStr '.nc'];

% Tsys threshold:
TSysCenter = [550 0 550 500];
TSysThresh = [25 0 25 25];
stdTSysThresh = [20 0 15 15];

frequencyBandAroundCenterTSys = 1e6*[20 10 200 200];
calibrationTool.frequencyBandAroundCenterTSys = frequencyBandAroundCenterTSys(modelFFTS);

calibrationTool.TSysCenterTh = TSysCenter(modelFFTS);
calibrationTool.TSysThresh = TSysThresh(modelFFTS);
calibrationTool.stdTSysThresh=stdTSysThresh(modelFFTS);



end