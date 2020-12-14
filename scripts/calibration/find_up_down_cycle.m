function [firstIndHalfUp,firstIndHalfDown] = find_up_down_cycle(logFile,calibrationTool)
%==========================================================================
% NAME          | find_up_down_cycle.m
% TYPE          | Nested function in calibrate_generic
% AUTHOR(S)     | Eric Sauvageat
% CREATION      | 01.2020
%               |
% ABSTRACT      | Find up and down cycle in a standard log file for debug
%               | mode
%               | 
% ARGUMENTS     | INPUTS: 1. logFile: standardized log file 
%               |         2. calibrationTool:
%               |               - indiceCold
%               |               - indiceHot
%               |               - indiceAntenna
%               |
%               | OUTPUTS: 1. firstIndHalfUp
%               |          2. firstIndHalfDown
%               |
%==========================================================================
indCold=calibrationTool.indiceCold;
indHot=calibrationTool.indiceHot;
indAntenna=calibrationTool.indiceAntenna;

extendedPos=[logFile.Position' -9999 -9999];

% Finding half cycle 'Up '(cah)
coldIndUp=find(logFile.Tipping_Curve_active==0 & (logFile.Position==indCold));
firstIndHalfUp=coldIndUp((extendedPos(coldIndUp+1)==indAntenna & extendedPos(coldIndUp+2)==indHot));

% Finding half cycle 'Down '(hac)
hotInd=find(logFile.Tipping_Curve_active==0 & (logFile.Position==indHot));
firstIndHalfDown=hotInd((extendedPos(hotInd+1)==indAntenna & extendedPos(hotInd+2)==indCold));
end
