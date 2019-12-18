function [log,rawSpectra,readingLevel0Error] = read_level0_generic(file)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


fid = fopen( [file '.txt'], 'r');
s = fgetl(fid);
header = textscan(s, '%s'); 
header = header{1}; % cell array with all header parameters
N = length(header); % number of header parameters
% x = fscanf(fid, '%f;', [N, inf]);  % data array
x = fscanf(fid, '%f ', [N, inf]);  % data array
M = size(x,2);     % number of data entries
fclose(fid);

D = dir([file '.bin']);
channels=D.bytes/4 /M;

log = inputArg1;

rawSpectra = inputArg2;

readingLevel0Error=2;
end

