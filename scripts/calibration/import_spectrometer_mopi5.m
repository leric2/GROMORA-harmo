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

calibrationTool.instrumentBandwidth = BW(modelFFTS);
calibrationTool.spectrometer=spectrometerTypes{modelFFTS};
calibrationTool.samplingRateFFTS=samplingRateFFTS(modelFFTS);
calibrationTool.filenameLevel1a=['/scratch/MOPI5/Level1/mopi5_level1a_' calibrationTool.spectrometer '_' calibrationTool.dateStr '.nc'];

end