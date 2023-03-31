function save_tipping_data_miawara(calibrationTool, freq, angle, Tb_tipped)
%script for saving the basic tipping calibrations to a txt file
%Author: Alistair Bell
%Date: 11/22

disp('starting write file')
fullfile = append(calibrationTool.level1Folder ,calibrationTool.TippingFilename);

title_string = 'Frequency';
for i=1:length(angle)
    temp_string = append('TB_', string(angle(i)));
    title_string = append(title_string, ' ');
    title_string = append(title_string, temp_string);
end
title_string = append(title_string, '\n');
fileID = fopen(fullfile,'w');
fprintf(fileID, title_string );

write_data(1,:) = freq;
write_data(2:length(angle)+1,:) = Tb_tipped;

disp('data formatted')

fileformat = '%6.2f';
for i=1:length(angle)
    fileformat = append(fileformat, '%6.2f ');
end
fileformat = append(fileformat, '\n');

disp('writing to file')
fprintf(fileID,fileformat, write_data);
disp('closing file')
fclose(fileID);
disp('file closed')

%end