% for two polarisations (OMT)
% counts2kelvin calibrates the balanced signal (data.S_line-data.S_ref) and applies the 
% correction for the reference absorber and the troposphere.
% 
% input:
% 
% data.S_line: the line signal
% data.S_ref:  the reference signal
% data.S_hot:  the hotload signal
% data.S_cold: the coldload signal
% data.Thot:  the mean hotload temperature
% T_COLD: the average sky temperature corresponding to the elevation angle
%         of data.S_cold
% data.Tref:  the mean reference absorber temperature
% data.Teff: the effective temperature of the troposphere
% data.T0:     the temperature of the cosmic background radiation
% data.tau:    the zenith opacity of the atmosphere
% x_x_elevation: the mirror elevation ot the line/ref signal
% 
% output:
% 
% Tb:     the calibrated balanced signal corrected for the troposphere and ref absorber
% T_rec:  the receiver temperature
% error:  an error variable -> +1 if 0<A<1 is wrong
%                           -> +1 if -10<median(Tb)<10 is wrong
% A:      the equivalent transmission of the reference absorber
% a:      the correction coefficient for the troposphere and ref absorber
% 
% 
% history:
% 
% - v3-1: accounts for the temperature of the reference absorber, 
% January 2008, A.Haefele
% - v4: bases on v3-1 and accounts for tropospheric and middle atmospheric
% airmass factors. July 2008, A.Haefele
% - v5: bases on v4 and accounts for antenna pattern
% - v6: bases on v5, for two polarisations (OMT)

function [Tb,T_COLD,T_rec,A,a, deltaTb]=counts2kelvin(data)
error=0;
%===== initialize output
T_rec=nan(size(data.S_line));

%===== load physical constants
%  physical_constants
RE = 6371000;

%===== 'divide by zero' has to be avoided. thus find indices where (data.S_hot-data.S_cold)~=0
nonzero_index=find(data.S_hot-data.S_cold~=0);
zero_index=find(data.S_hot-data.S_cold==0);
data.S_line=data.S_line(nonzero_index);
data.S_ref=data.S_ref(nonzero_index);
data.S_hot=data.S_hot(nonzero_index);
data.S_cold=data.S_cold(nonzero_index);
%  OFFSET=OFFSET(nonzero_index);


%===== calculate the airmass factor
antenna= load('antenna.txt');
norm=sum(antenna(:,2));
antenna_pat=antenna(:,2)/norm;


% zenith angles
za_line=(90-data.mirror_el_line)+antenna(:,1);
za_ref=90-data.mirror_el_ref+antenna(:,1);
if data.mirror_el_ref>170
za_ref= (180-data.mirror_el_ref)*2+antenna(:,1);
end
if data.mirror_el_ref<10
za_ref= (-data.mirror_el_ref)*2+antenna(:,1);
end

za_cold=90-data.mirror_el_cold+antenna(:,1);

% thickness of the troposphere and the middle atmosphere [m]
H=90e3;
h=10e3;

% calculate airmass
% line
s2_line=sqrt((RE+H+h)^2-RE^2*sin(za_line/180*pi).^2)-RE*cos(za_line/180*pi);
s1_line=sqrt(  (RE+h)^2-RE^2*sin(za_line/180*pi).^2)-RE*cos(za_line/180*pi);
AM_line_trop1=s1_line/h.*antenna_pat;
AM_line_ma1=(s2_line-s1_line)/H.*antenna_pat;
%  figure(1),hold on,plot(za_line,AM_line_trop1)
AM_line_trop=sum(AM_line_trop1);
AM_line_ma=sum(AM_line_ma1);

% ref
s2_ref=sqrt((RE+H+h)^2-RE^2*sin(za_ref/180*pi).^2)-RE*cos(za_ref/180*pi);
s1_ref=sqrt(  (RE+h)^2-RE^2*sin(za_ref/180*pi).^2)-RE*cos(za_ref/180*pi);
AM_ref_trop1=s1_ref/h.*antenna_pat;
AM_ref_ma1=(s2_ref-s1_ref)/H.*antenna_pat;
%  figure(2),hold on,plot(za_ref,AM_ref_trop1)
AM_ref_trop=sum(AM_ref_trop1);
AM_ref_ma=sum(AM_ref_ma1);

% cold
s1_cold=sqrt(  (RE+h)^2-RE^2*sin(za_cold/180*pi).^2)-RE*cos(za_cold/180*pi);
AM_cold_trop1=s1_cold/h.*antenna_pat;
AM_cold_trop=sum(AM_cold_trop1);

%===== calculate T_COLD
T_COLD=data.T0*exp(-data.tau*AM_cold_trop)+data.Teff*(1-exp(-data.tau*AM_cold_trop));


%===== calculate the brightness temperature of the reference signal by means of calibration
Tb_ref(nonzero_index)=(data.S_ref-data.S_hot)./(data.S_hot-data.S_cold).*(data.Thot-T_COLD)+data.Thot;
Tb_ref(zero_index)=nan;

%===== calculate the brightness temperature of the reference signal by means of the RTE
Tbref=data.T0*exp(-data.tau*AM_ref_trop)+data.Teff*(1-exp(-data.tau*AM_ref_trop));

%===== calculate the equivalent transmission of the ref absorber
A=nanmean((Tb_ref-data.Tref)./(Tbref-data.Tref));

%===== calibrate the difference spectrum
deltaTb(nonzero_index)=(data.S_line-data.S_ref)./(data.S_hot-data.S_cold).*(data.Thot-T_COLD);
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

%===== increase the error if the ref abs transmission is out of range
if A>1 || A<0
   error=error+16;
   disp('WARNING: A is out of range!')
end

%===== calculate the receiver temperature using the y-method
nonzero_index=find(data.S_cold~=0);
y=(data.S_hot(nonzero_index))./(data.S_cold(nonzero_index));
T_rec(nonzero_index)=(data.Thot-y*T_COLD)./(y-1);
