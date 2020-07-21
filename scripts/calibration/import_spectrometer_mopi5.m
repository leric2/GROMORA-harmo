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



calibrationTool.badChannels = [];

calibrationTool.instrumentBandwidth = BW(modelFFTS);
calibrationTool.spectrometer=spectrometerTypes{modelFFTS};
calibrationTool.samplingRateFFTS=samplingRateFFTS(modelFFTS);
calibrationTool.filenameLevel1a=['/scratch/MOPI5/Level1/mopi5_level1a_' calibrationTool.spectrometer '_' calibrationTool.dateStr '.nc'];

% Tsys threshold:
TSysCenter = [780 0 545 505];
TSysThresh = [30 0 20 20];
stdTSysThresh = [15 0 10 10];

frequencyBandAroundCenterTSys = 1e6*[20 10 200 200];
calibrationTool.frequencyBandAroundCenterTSys = frequencyBandAroundCenterTSys(modelFFTS);

calibrationTool.TSysCenterTh = TSysCenter(modelFFTS);
calibrationTool.TSysThresh = TSysThresh(modelFFTS);
calibrationTool.stdTSysThresh=stdTSysThresh(modelFFTS);


end