function save_difference_data_miawara(calibrationTool,  hot_differences, ref_differences, cold_differences, xtra_fileID)
%script for saving the basic tipping calibrations to a txt file
%Author: Alistair Bell
%Date: 11/22

disp('starting write file')
filename = ['position_diff_' calibrationTool.dateStr '_' xtra_fileID '.dat' ];
fullfile = append(calibrationTool.level1Folder ,filename);

title_string = 'hot_diff ret_diff cold_diff\n';

fileID = fopen(fullfile,'w');
fprintf(fileID, title_string );

write_data(1,:) = hot_differences;
write_data(2,:) = ref_differences;
write_data(3,:) = cold_differences;

disp('data formatted')

fileformat = '%6.6f %6.6f %6.6f\n';

disp('writing to file')
fprintf(fileID, fileformat, write_data);
disp('closing file')
fclose(fileID);
disp('file closed')

%end