function log = harmonize_log_gromos(log)
% 
% Function to harmonize the log of Gromos
log.T_Hot=log.T_Hot+273.15;
end
