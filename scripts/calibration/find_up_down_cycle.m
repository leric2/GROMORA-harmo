% Nested function in calibrate_generic (copied here for convenience)
function [firstIndHalfUp,firstIndHalfDown] = find_up_down_cycle(logFile,calibrationTool)
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
