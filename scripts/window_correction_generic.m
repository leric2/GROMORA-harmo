function calibratedSpectra = window_correction_generic(calibratedSpectra,retrievalTool)
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

Twin = mean(log.TCAB7(is(1):is(end)));      % Temperature read by sensor close to the window  (Compare TCAB7 instead of TCAB5)
Tbw  = (Tb - Twin * 0.0012) ./ 0.9988;

