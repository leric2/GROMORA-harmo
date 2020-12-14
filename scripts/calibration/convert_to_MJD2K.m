function mjd2k_date = convert_to_MJD2K(YYYY, MM, DD, hh, mm, ss)
%==========================================================================
% NAME      | convert_to_MJD2K(YYYY, MM, DD, hh, mm, ss)
% TYPE      | Function
% AUTHOR(S) | Eric Sauvageat
% CREATION  | 12.2020
%           |
% ABSTRACT  | Function to compute the MJD2K time. Taken from GEOMS report 1.0:
%           | https://avdc.gsfc.nasa.gov/PDF/GEOMS/geoms-1.0.pdf
%           |
%           |
%           |
% ARGUMENTS | INPUTS: YYYY,MM, DD, hh, mm, ss in UTC time.
%           |
%           | OUTPUTS: - mjd2k_date: a double number corresponding to the
%           | Modified Julian Datetime 2000. It corresponds to the number
%           | of day from 01.01.2000. 
%           | 
%==========================================================================

if MM>2
    y = double(YYYY);
    m = double(MM-3);
    d = double(DD);
else
    y = double(YYYY-1);
    m = double(MM + 9);
    d = double(DD);
end

j = round(365.25*(y+4712)) + round(30.6*m+0.5) + 58.5 + d;

if j < 2299159.5
    jd = j;
else
    gn = 38 - round(3*round(49 + y/100)/4);
    jd = j+gn;
end

df = (hh*3600.0 + mm*60.0 + ss)/86400;

mjd2k_date = jd + df - 2451544.5;
