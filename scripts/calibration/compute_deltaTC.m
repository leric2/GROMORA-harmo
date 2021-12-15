function dTc = compute_deltaTC(calibrationTool,calibratedSpectra, ind_before, ind_after)
%==========================================================================
% NAME      | planck_Tb.m
% TYPE      | function
% AUTHOR(S) | Eric Sauvageat
% CREATION  | 2021
%           |
% ABSTRACT  | Function computing the Planck brightness temperature. It
%           |
%           |
% ARGUMENTS | INPUTS:   1. calibrationTool:
%           |               - h, Planck's constant
%           |               - kb, Boltzmann’s constant
%           |               - lightSpeed
%           |           2. intensity:
%           |               The intensity of radiation
%           |           3. frequency
%           |
%           | OUTPUTS:  1. TbPlanck: Planck brightness temperature.
%           | Depending on the inputs, it can be a vector.
%           |
%==========================================================================

errorMatrix = vertcat(calibratedSpectra.errorVector);
first_refill = find(errorMatrix(:,3)==0 | errorMatrix(:,4)==0);
if nargin < 3
    if length(first_refill) > 1
        first_refill = first_refill(2);
        ind_before = first_refill - 3;
        ind_after = first_refill + 3;
    elseif isempty(first_refill)
        error('no refill found)') 
    end
end
dTc = (calibrationTool.TCold-calibratedSpectra(ind_before).THot)*(calibratedSpectra(ind_before).mean_Yfactor-calibratedSpectra(ind_after).mean_Yfactor)/ ...
    (calibratedSpectra(ind_before).mean_Yfactor*(calibratedSpectra(ind_after).mean_Yfactor-1));
disp(['Tcold change by ' num2str(dTc) 'K during refill.']);
end