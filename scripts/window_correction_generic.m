function level1b = window_correction_generic(retrievalTool,level1b)
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
for t = 1:length(level1b.integration)
    freq=level1b.integration(t).freq;
    Tb=level1b.integration(t).Tb;
    TWindow=level1b.integration(t).TWindow;
    
    % Planck:
    TbWindowP= (2*retrievalTool.h*freq.^2)/(retrievalTool.lightSpeed^2)*(1)./(exp((t* freq)./(retrievalTool.kb*TWindow))-1);
    
    % Railey-Jeans ??
    % TbWindowRJ = (retrievalTool.lightSpeed^2 ./ (2*retrievalTool.kb*freq.^2) ) .* TbWindowP; 

    level1b.integration(t).Tbcorr  = (Tb - TbWindowP*(1-retrievalTool.tWindow))./retrievalTool.tWindow;
end
end
