function calibratedSpectra = window_correction_generic(calibratedSpectra)
%==========================================================================
% NAME          | 
% TYPE          |
% AUTHOR(S)     | Susana Fernandez (adapted to GROSOM by ES)
% CREATION      | 09.2014
%               |
% ABSTRACT      |
%               | 
%               |
%               |
% ARGUMENTS     | INPUTS:
%               |
%               | OUTPUTS:
%               |
% CALLS         |
%               |
%               |
%               |

%==========================================================================

% Window correction
for i = 1:length(calibratedSpectra)
    calibratedSpectra(i).Tbw  = (calibratedSpectra(i).Tb - calibratedSpectra(i).TWindow * 0.0012) ./ 0.9988;
end
