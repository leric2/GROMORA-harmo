function log = harmonize_log_somora(log)
% 
% Function to harmonize the somora log with the gromos basis

log.Tipping_Curve_active=(log.Position==6);
end
