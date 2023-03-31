function  Tb_tipping = two_angle_tipping(TC_data, calibrationTool,s_tipping,s_hot,s_cold)   

% NAME      | two_angle_tipping.m
% TYPE      | Script
% AUTHOR(S) | Alistair Bell
% CREATION  | 11.2022
%           |
% ABSTRACT  | Tipping curve calibration from two angles, one low zenith
%           | angle (cold) and one high zenith angle (line). The purpose 
%           | of this file is to validate the new USRP spectrometer. 
%           |
% ARGUMENTS | INPUTS: - ADD ALL INPUTS
%           |          
%           |
%           | OUTPUTS: ADD ALL OUTPUTS 
%           | 
%           |
% CALLS     | ADD ALL CALLS
%           | 
%           |
%           |
%==========================================================================

T_0 = 2.7254;

c1 = 0.69;
c0 = 266.3;
Teff = c1 * (TC_data.Tamb-273) + c0;% mean tropospheric temperature

T0 = 2.7; % cosmic background
Z  = 10e3; %thickness of a relatively thin layer of earths atmosphere (troposphere)
RE = 6371000;%radius of earth used for calculation of airmass factor


if nargin<3
    s_tipping = TC_data.s_tipping([1,end],:);
    za_tipping = TC_data.za_tipping( [1,end]);
    s_hot     = TC_data.s_hot;
    s_cold    = TC_data.s_cold;
    za_cold = TC_data.za_cold;
end

numer = 1 + Z/RE;
denom =sqrt(sin((90-za_tipping)./180*pi).^2 + 2*(Z/RE)+(Z^2/RE^2)) ;
AM = numer./denom;

denom_c =sqrt(sin((90-za_cold)./180*pi).^2 + 2*(Z/RE)+(Z^2/RE^2)) ;
am_cold = numer./denom_c;

%%initialise
Tb = nan(size(s_tipping,2),2);
offset = 1;
niter  = 0;
tau    = .3*ones(size(s_tipping,2),1);
niter_max = 20;

disp('starting sub-iteration')
while mean(offset, 'omitnan') > 20e-3  &&  niter < niter_max

    Tcold = T0 * exp(-tau.*am_cold) + Teff * (1 - exp(-tau.*am_cold));
    za_cold    = TC_data.za_cold;
    za_tipping = TC_data.za_tipping;
    Thot       = TC_data.Thot;
    
    disp('calculating average tbs')
    for i=1:size(s_tipping,1)
        Tb(:,i) = ((s_tipping(i,:)-s_hot)./(s_hot-s_cold).*(Thot-Tcold')+Thot);%AB
    end
    
    tau_slant = log((Teff-T0)./(Teff-Tb));
    for i = 1:size(tau)
        [p,s] = polyfit(AM, tau_slant(i,:), 1);
        tau(i)      = p(1);
        offset(i)   = p(2);
    end 
    niter = niter+1;
    
    if niter >= niter_max
      disp(['niter: ' num2str(niter) ' !!'])
      disp('calculating opacity: iteration failed!');
      tau         = nan;
      Trec        = nan;
      Trec_median = nan;
      return
    end
end
disp('end of iteration')
Tb_tipping = nan(size(s_tipping));
for i=1:size(s_tipping,1)
    Tb_tipping(i,:) = ((s_tipping(i,:)-s_hot)./(s_hot-s_cold).*(Thot-Tcold')+Thot);
end

end



