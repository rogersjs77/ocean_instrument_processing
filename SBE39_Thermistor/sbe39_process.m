% function []=sbe56_process(filter_time,start_time,end_time)
% Process SBE 39 Thermistor Data
% Justin Rogers, Stanford EFML, 2014
%
% this file takes SBE39 .cnv files, with either
% column 1 = time in julian days
% column 2 = temperature
% column 3 = flag
%
% file name format: 'SBE39_XXMAB_...asc', where XX is vertical placement
%
% INPUT:
% filter_time in minutes
% start_time in matlab format (optional)
% end_time in matlab format (optional)
% start/end times can be single values applied to all files
% or vectors for each individual file
% OUTPUT:
% .mat file containing variables
%%

function []=sbe39_process(filter_time,start_time,end_time)

%% load data files
files=dir('SBE39*.asc');
if isempty(files)
   disp('no files found, check format SBE39_XXMAB_...asc')
   return
end

for i=1:size(files,1); % get order for files in increasing mab
filename=files(i).name;
mab(i) = sscanf(filename,'SBE39_%f MAB_%*f');
end
[mab,ii]=sort(mab,2);
sbe56.mab=mab;

%% load data
for i=1:length(files);
filename=files(ii(i)).name;
delimiter = ',';
startRow = 30;

% Read columns of data as strings:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
textscan(fileID, '%[^\n\r]', startRow-1, 'ReturnOnError', false);
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN,'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

% Allocate imported array to column variable names
sbe56.temp{i} = dataArray{:, 1};
for jj=1:length(dataArray{:, 2})
sbe56.mtime_raw{i}(jj)=datenum([cell2mat(dataArray{2}(jj)) ' ' cell2mat(dataArray{3}(jj))],'dd mmm yyyy HH:MM:SS');
end
% Clear temporary variables
clearvars filename delimiter startRow formatSpec fileID dataArray ans jj;

%%


% clean up mtime if it has bad time readings
sbe56.mtime{i} = sbe56.mtime_raw{i};
sbe56.mtime{i}(sbe56.mtime{i}<sbe56.mtime{i}(1))=nan;
if isnan(sbe56.mtime{i}(end)) % find last good time and interpolate out
    indx = max(find(~isnan(sbe56.mtime{i})));
    sbe56.mtime{i}(end) = sbe56.mtime{i}(indx)+...
        (sbe56.mtime{i}(2)-sbe56.mtime{i}(1))*(length(sbe56.mtime{i})-indx);
end
% interpolate out nans
sbe56.mtime{i} = naninterp(sbe56.mtime{i});

if nargin==3 % user specified start/end
    if length(start_time)==1 % only one start end
    time_start{i}=start_time;
    time_end{i} = end_time; 
    else % vetor of start end times
    time_start{i}=start_time(i);
    time_end{i} = end_time(i);
    end
else % find trim time from temp within 1std of mean
indx = find(sbe56.temp{i}>median(sbe56.temp{i})-std(sbe56.temp{i}) & ...
    sbe56.temp{i}<median(sbe56.temp{i})+std(sbe56.temp{i}));
% trim mtime to time underwater
time_start{i}=sbe56.mtime{i}(indx(1));
time_end{i}=sbe56.mtime{i}(indx(end));
end

time_id{i}= sbe56.mtime{i}>=time_start{i} & ...
    sbe56.mtime{i}<=time_end{i};

% time time and temp to limits
mtime{i}=sbe56.mtime{i}(time_id{i});
dt{i} = (mtime{i}(end)-mtime{i}(1))./length(mtime{i}); % data time step (days)
temp{i} = sbe56.temp{i}(time_id{i});

% low pass filter temperature data
[temp_fft{i},FFT_temp{i},FFT_s{i}] = low_pass_filter(temp{i},dt{i}*24*3600,filter_time*60);
clear DATA indx
end

%% downsample filtered temp data to a common time
dt_filt=filter_time*60*(1/86400); % interpolate to 10 min window
mtime_filt=min(cell2mat(time_start)):dt_filt:max(cell2mat(time_end));

temp_filt=nan(length(mab),length(mtime_filt));
for i=1:length(mab)
    temp_filt(i,:)=interp1(mtime{i}, temp_fft{i},mtime_filt,'linear',NaN);
end

%% clean up and save
clear temp_fft DATA i ii sn_v start_year indx filename start_year temp_fft
save('Processed_SBE56_data.mat');

end

