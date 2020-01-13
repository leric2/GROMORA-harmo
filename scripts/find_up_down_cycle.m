% Nested function in calibrate_generic (copied here for convenience)
function [firstIndHalfUp,firstIndHalfDown] = find_up_down_cycle(log,retrievalTool)
indCold=retrievalTool.indiceCold;
indHot=retrievalTool.indiceHot;
indAntenna=retrievalTool.indiceAntenna;

extendedPos=[log.Position -9999 -9999];

% Finding half cycle 'Up '(cah)
coldIndUp=find(log.Tipping_Curve_active==0 & (log.Position==indCold));
firstIndHalfUp=coldIndUp((extendedPos(coldIndUp+1)==indAntenna & extendedPos(coldIndUp+2)==indHot));

% Finding half cycle 'Down '(hac)
hotInd=find(log.Tipping_Curve_active==0 & (log.Position==indHot));
firstIndHalfDown=hotInd((extendedPos(hotInd+1)==indAntenna & extendedPos(hotInd+2)==indCold));
end
