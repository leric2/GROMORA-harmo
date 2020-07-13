function integratedSpectra = check_integrated_generic(calibrationTool,integratedSpectra)
%==========================================================================
% NAME      | check_integrated_generic.m
% TYPE      | function
% AUTHOR(S) | Eric Sauvageat
% CREATION  | 2020
%           |
% ABSTRACT  | Function performing some checks on the calibration and
%           | adding meta information to the integratedSpectra structure.
%           | It builts also the flags vector for the level1a data for each
%           | calibration cycle and add every new information to the
%           | calibrated spectra structure (IN/OUT).
%           | 
%           |
% ARGUMENTS | INPUTS:   - standardLog: harmonized GROSOM log file 
%           |           - calibrationTool
%           |           - integratedSpectra
%           |
%           |
%           | OUTPUTS: - integratedSpectra
%           |
%           |
% CALLS     | 
%           | 
%           | 
%           | 
%           |
%==========================================================================
% Checking all calibration cycle
for i = 1:size(integratedSpectra,2)
    
    integratedSpectra(i).meanDatetime=datenum(integratedSpectra(i).dateTime)-datenum(1970,1,1);
    
    calibrationTool.minNumberOfAvgSpectra = 2;
    
    
    %%%%%%%%%%% Flag 1 %%%%%%%%%%%
    % The number of indices for the 3 positions:
    if (integratedSpectra(i).numberOfAveragedSpectra > calibrationTool.minNumberOfAvgSpectra)
        sufficientNumberOfAvgSpectra=1;
    else
        sufficientNumberOfAvgSpectra=0;
        warning('Low number of avg spectra for this integration');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Antenna angle
    %integratedSpectra(i).stdAngleAntenna=std(standardLog.Elevation_Angle(ia));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % System Temperature
    % Computing TN around the line center (approximately +- 200 MHz)
    
    
    %%%%%%%%%%% Flag  %%%%%%%%%%
    % Rain
    calibrationTool.rainAccumulationThreshold = 0.1;
    if (integratedSpectra(i).rainAccumulation <= calibrationTool.rainAccumulationThreshold)
        rain_Accumulation_OK=1;
    else
        rain_Accumulation_OK=0;
    end
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    integratedSpectra(i).numberOfIndices=[
        integratedSpectra.numHotSpectra,...
        integratedSpectra.numColdSpectra,...
        integratedSpectra.numAntSpectra];
    
    % Error vector for this calibration cycle
    integratedSpectra(i).errorVector=[
        sufficientNumberOfAvgSpectra,...
        %systemTemperatureOK,...
        %LN2SensorsOK,...
        %LN2LevelOK,...
        %hotLoadOK,...
        rain_Accumulation_OK];
    
    % Error vector description:
    integratedSpectra(i).errorVectorDescription=[
        "sufficientNumberOfAvgSpectra",...
        %"systemTemperatureOK",...
        %"LN2SensorsOK",...
        %"LN2LevelOK",...
        %"hotLoadOK",...
        "rain_Accumulation_OK"];
    
%     if (sum(integratedSpectra(i).errorVector)<6)
%         errorV=num2str(integratedSpectra(i).errorVector);
%         disp(['Calibration Cycle number ' num2str(i) ', TOD: ' num2str(integratedSpectra(i).timeOfDay)])
%         warning(['Problem with this calibration, error code : ' errorV]);
%         disp(integratedSpectra(i).errorVectorDescription(~integratedSpectra(i).errorVector))
%     end
   
end

end

