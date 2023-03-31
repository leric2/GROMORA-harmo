function [Tb,T_COLD,T_rec,A,a, deltaTb]=LN2_counts2kelvin(data)
error=0;
%===== initialize output
T_rec=nan(size(data.S_line));

%===== calculate T_COLD
T_COLD=data.T0*exp(-data.tau*AM_cold_trop)+data.Teff*(1-exp(-data.tau*AM_cold_trop));

%===== calculate the brightness temperature of the reference signal by means of calibration
Tb_ref(nonzero_index)=(data.S_ref-data.S_hot)./(data.S_hot-data.S_cold).*(data.Thot-T_COLD)+data.Thot;
Tb_ref(zero_index)=nan;

%===== calculate the brightness temperature of the reference signal by means of the RTE
Tbref=data.T0*exp(-data.tau*AM_ref_trop)+data.Teff*(1-exp(-data.tau*AM_ref_trop));

%===== calculate the equivalent transmission of the ref absorber
A=mean((Tb_ref-data.Tref)./(Tbref-data.Tref), 'omitnan');%AB

%===== calibrate the difference spectrum
T_COLD(nonzero_index)=(data.S_cold-data.S_LN2)./(data.S_hot-data.S_LN2).*(data.Thot-data.T_LN2)+data.T_LN2;
deltaTb(zero_index)=nan;

%===== finally, calculate the factor for tropospheric correction
a=1/(AM_line_ma*exp(-data.tau*AM_line_trop)-A*AM_ref_ma*exp(-data.tau*AM_ref_trop));
Tb=deltaTb*a;

%===== check, wether Tb is well balanced (median(Tb)=[-5 5]) and set process flag accordingly
nonnan_index=find(isnan(Tb)==0);
if median(Tb(nonnan_index))<-5 || median(Tb(nonnan_index))>5
  error=error+8;
  disp('WARNING: spectrum balance is bad!')
end

if isfield(data, 'errorMax')
    if error>=data.errorMax
        ErrMaxTrue = true;
    else
        ErrMaxTrue = false;
    end
else
    ErrMaxTrue = false;
end

if ErrMaxTrue
    disp('WARNING: data excluded from calibration cycle')
    %if errors more than specified, do not include data in calibration
    Tb = ones(size(Tb)) * nan;
    T_COLD = nan;
    T_rec = ones(size(Tb)) * nan;
    a =  nan;
    A = nan;
    deltaTb = ones(size(Tb)) * nan;
else 
    if (~isfield(data, 'errorMax') ||  error<data.errorMax)
        %===== calculate the receiver noise temperature using the y-method
        nonzero_index=find(data.S_cold~=0);
        y=(data.S_hot(nonzero_index))./(data.S_cold(nonzero_index));
        T_rec(nonzero_index)=(data.Thot-y*T_COLD)./(y-1);
    end
end
end
