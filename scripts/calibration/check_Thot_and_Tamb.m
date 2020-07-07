function [Thot, Tamb]= check_Thot_and_Tamb(Thot_time, Thot_in, Tamb_time, Tamb_in, take_tamb) 

if nargin<5
    take_tamb=0;
end

if take_tamb
    
    %===== take Tamb as hotload temp
    
    Thot = interp1(Tamb_time,Tamb_in,Thot_time);
    Tamb = Thot;
	disp('no Thot found!')
    
else
    
    %===== hot load temp
    
    % check Thot 
    ind1 = find (Thot_in>200 & Thot_in~=273.15 & Thot_in~=0);
    if isempty(ind1)
    	Thot=[];
    else
    	Thot_mean1 = mean(Thot_in(ind1));
    	ind2 = find (abs(Thot_in-Thot_mean1)<2);
    	if isempty(ind2)
    		Thot = [];
    	else
    		Thot = mean (Thot_in(ind2));
    	end
    end

    % check Tamb
    ind2 = find (Tamb_in>200 & Tamb_in~=273.15);
    
    % find closest Tamb for given Thot
    [value, meteo_index] = min (abs (Thot_time-Tamb_time(ind2)));

    if value < 0.007      % ~10min
       	Tamb = Tamb_in(ind2(meteo_index));
    	if isempty(Thot)
           Thot=Tamb;
           disp('no Thot found!')
        end 
    elseif isempty(ind2)
        Tamb = Thot;
        disp('no Tamb found!')
    elseif  value > 0.007 && ~isempty(Thot)
        Tamb = Thot;
        disp('no Tamb found!')
    else
        Thot=0;
        Tamb=0;
        disp('no temperature data found!')
    end
    
end