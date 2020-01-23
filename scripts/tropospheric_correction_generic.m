function log = tropospheric_correction_generic(calibratedSpectra,deltaT)
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

% Mean tropospheric temperature

Tmean = calibratedSpectra.Tair - deltaT;          

Tbg = 2.7;                   % Microwave background    


% Linear fit of spectrum's wings 
   
f_trop_corr  = [calibratedSpectra.freq(100:5000)   calibratedSpectra.f(end-5000:end-100)  ];
Tb_trop_corr = [calibratedSpectra.Tbw(100:5000) calibratedSpectra.Tbw(end-5000:end-100)];
    
[p,s,mu] = polyfit(f_trop_corr,Tb_trop_corr,1 );  % linear fit
Twing    = polyval(p, spectrum.f, [], mu);        % polynomial evaluation


% Transmitance calculated (Ingold) 

transmittance = (Tmean - Twing)./(Tmean - Tbg); 

% Troposph. corr. 

spectrum.Tbw_corr_trop       = ( spectrum.Tbw - Tmean*(1-transmittance) ) ./ transmittance; 
spectrum.trop_transmit       = transmittance;
spectrum.trop_transmit_mean  = mean(transmittance);