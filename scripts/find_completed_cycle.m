% Nested function in calibrate_generic (copied here for convenience)
function firstIndCompleteCycle = find_completed_cycle(log,retrievalTool)
indCold=retrievalTool.indiceCold;
indHot=retrievalTool.indiceHot;
indAntenna=retrievalTool.indiceAntenna;
hotInd=find(log.Tipping_Curve_active==0 & (log.Position==indHot));
extendedPos=[log.Position -9999 -9999 -9999 -9999 -9999];
firstIndCompleteCycle=hotInd((extendedPos(hotInd+1)==indAntenna & extendedPos(hotInd+2)==indCold & extendedPos(hotInd+3)==indCold & extendedPos(hotInd+4)==indAntenna & extendedPos(hotInd+5)==indHot));
end

