

% for two polarisations (OMT)
% take antenna pattern into account and calculate airmass with geometrical formula instead of secans, corinne 2010-07-28



function  [tau, Teff, Trec_median, quality, offset, niter, time] = find_tau_iteratively(TC_data, calibrationTool,s_tipping,s_hot,s_cold)   
         %[tau, Teff, error, Trec_median, quality, offset, niter]=
              %S,Shot,za1,Tamb,Thot,spectrometer_offset, error)


if nargin<3
    s_tipping = TC_data.s_tipping;
    s_hot     = TC_data.s_hot;
    s_cold    = TC_data.s_cold;
end 

za_cold    = TC_data.za_cold;
za_tipping = TC_data.za_tipping;
Thot       = TC_data.Thot;
   
time = TC_data.time;

% mean tropospheric temperature

c1 = 0.69;
c0 = 266.3;
Teff = c1 * (TC_data.Tamb-273) + c0;


% cosmic background

T0 = 2.7;
H  = 90e3;
h  = 10e3;
RE = 6371000;


% airmass factor

if isfield(calibrationTool,'antenna_file')
    
    antenna     = load(calibrationTool.antenna_file);
    norm        = sum(antenna(:,2));
    antenna_pat = antenna(:,2)/norm;

    za_cold  = za_cold + antenna(:,1);
    s_cold2  = sqrt((RE+h)^2-RE^2*sin(za_cold/180*pi).^2)-RE*cos(za_cold/180*pi);
    AM_cold1 = s_cold2/h.*antenna_pat;
    am_cold  = sum(AM_cold1);

    za  = ones(length(antenna(:,1)),1) * za_tipping + antenna(:,1) * ones(1,length(za_tipping));
    s   = sqrt((RE+h)^2-RE^2*sin(za/180*pi).^2)-RE*cos(za/180*pi);
    AM1 = s/h.* (antenna_pat *ones(1,length(za_tipping)) );
    am  = sum(AM1)';
else
    am_cold = [];
    am      = [];  
end


%===== some setup values

offset = 1;
niter  = 0;
tau    = .3;

niter_max = 10;


%===== preallocate Tb matrix

Tb = nan(size(s_tipping,1),1);


%===== iteration loop



while abs(offset) > 20e-3  &&  niter < niter_max

  % calculate the cold load brightnes temperature
  
  Tcold = T0 * exp(-tau*am_cold) + Teff * (1 - exp(-tau*am_cold));

  
  % calibrate the signals S
  
  for i=1:size(s_tipping,1)
     Tb(i) = nanmean( ((s_tipping(i,:)-s_hot)./(s_hot-s_cold).*(Thot-Tcold)+Thot)')';
  end

  
  % calculate the slant opacities
  
  tau_slant = log((Teff-T0)./(Teff-Tb));

  
  % fit the airmass-slant opacity data pairs
  
  [p,s] = polyfit (am, tau_slant, 1);
  quality = s.normr/size(s_tipping,1);

  
  % the slope is the new tau
  
  tau      = p(1);
  offset   = p(2);
  % tau_plot = polyval(p,am);

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


%===== calculate Trec to decide on the quality of tau

%  
   Tcold = T0 * exp(-tau*am_cold) + Teff * (1-exp(-tau*am_cold));
   y     = (s_hot)./(s_cold);
   Trec  = (Thot-y.*Tcold)./(y-1);
   Trec_median = median(Trec);
%  
   if Trec_median<110 || Trec_median>200
      disp('calculating opacity: Trec out of range!')
   end
   
% TC.tau         = tau;
% TC.Teff        = Teff;
% TC.Trec_median = Trec_median;
% TC.quality     = quality;
% TC.offset      = offset;
% TC.niter       = niter;
   
end
   