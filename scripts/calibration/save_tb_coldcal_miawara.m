function save_tb_coldcal_miawara(calibrationTool,  coldcal_Tb, xtra_fileID)
%script for saving the basic tipping calibrations to a txt file
%Author: Alistair Bell
%Date: 11/22

disp('starting write file')
filename = ['cold_calibration_' calibrationTool.dateStr '_' xtra_fileID '.dat' ];
fullfile = append(calibrationTool.level1Folder ,filename);

title_string = 'Tb \n';

fileID = fopen(fullfile,'w');
fprintf(fileID, title_string );

fileformat = '%6.6f\n';

disp('writing to file')
fprintf(fileID,fileformat, coldcal_Tb);
disp('closing file')
fclose(fileID);
disp('file closed')

%end