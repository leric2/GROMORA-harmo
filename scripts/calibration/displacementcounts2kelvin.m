function [Tb_ref_diff,Tb_hot_diff, Tb_cold_diff]=displacementcounts2kelvin(data)
error=0;
%===== initialize output
T_return=nan(size(data.S_ref_disp1));

%===== load physical constants
%  physical_constants
RE = 6371000;

%===== calculate the airmass factor
antenna= load('/home/alistair/MIAWARA_ret/MIAWARA_pyarts/GROMORA-harmo/files/miawarac_antenna.txt');
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

za_cold=90-data.mirror_el_cold_disp1 + antenna(:,1);

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


% cold
s1_cold=sqrt( (RE+h)^2-RE^2*sin(za_cold/180*pi).^2)-RE*cos(za_cold/180*pi);
AM_cold_trop1=s1_cold/h.*antenna_pat;
AM_cold_trop=sum(AM_cold_trop1);

%===== calculate T_COLD
T_COLD=data.T0*exp(-data.tau*AM_cold_trop)+data.Teff*(1-exp(-data.tau*AM_cold_trop));

%===== calculate the brightness temperature of the reference signal by means of calibration
%Tb_ref(nonzero_index)=(data.S_ref-data.S_hot)./(data.S_hot-data.S_cold).*(data.Thot-T_COLD)+data.Thot;
%Tb_ref(zero_index)=nan;

%===== Calculate differences in brightness temperature of differences in
%reference signal

SmeanHot = (data.S_hot_disp1+data.S_hot_disp2)/2;
SmeanCold = (data.S_cold_disp1+data.S_cold_disp2)/2;

Thotmean = (data.Thot_disp1 + data.Thot_disp2)/2;

Tb_ref_diff=(data.S_ref_disp1-data.S_ref_disp2)./(SmeanHot-SmeanCold).*(Thotmean-T_COLD);
Tb_hot_diff=(data.S_hot_disp1-data.S_hot_disp2)./(SmeanHot-SmeanCold).*(Thotmean-T_COLD);
Tb_cold_diff=(data.S_cold_disp1-data.S_cold_disp2)./(SmeanHot-SmeanCold).*(Thotmean-T_COLD);

