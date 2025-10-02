% function []=sbe56_process(filter_time,start_time,end_time)
% Process SBE 56 Thermistor Data
% Justin Rogers, Stanford EFML, 2014
%
% this file takes SBE56 .cnv files, with either
% column 1 = time in julian days
% column 2 = temperature
% column 3 = flag
%  OR
% colum 1 = record
% colum 2 = time in julian days
% column 3 = temperature
% column 4 = flag
%
% file name format: 'SBE56_XXMAB_YYY_...cnv', where XX is vertical placement
% and YYY is site name
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

function []=sbe56_process_multiple_sites(filter_time,start_time,end_time)

%% load data files
files=dir('SBE56*.cnv');
if isempty(files)
   disp('no files found, check format SBE56_XXMAB_...cnv')
   return
end

for i=1:size(files,1); % get order for files in increasing mab
filename=files(i).name;
C = strsplit(filename,'_');
mab(i) = sscanf(C{2},'%f MAB');
site{i} = C{3};
sn(i) = sscanf(C{4},'SN %f .cnv');
end

sbe56.mab=mab;
sbe56.site = site;
sbe56.sn = sn;
clear C 
%% load data
for i=1:length(files);
filename=files(i).name;
DATA = importdata(filename,' ',90);
sbe56.hdr{i}=DATA.textdata(:,1);

% find start year for date
j=1;start_year=[];
while isempty(start_year)
%     # start_time = Oct 01 2012 07:19:19   
    start_year=sscanf(sbe56.hdr{i}{j},'# start_time = %*s %*f %f %*f');
    j=j+1;
end
if size(DATA.data,2)==3 % no record number, 3 columns
    sbe56.temp{i}=DATA.data(:,2)';
    sbe56.mtime_raw{i}=jday2matlab(DATA.data(:,1),start_year)'; % convert julian dats to mtime
else % record number, 4 columns
    sbe56.temp{i}=DATA.data(:,3)';
    sbe56.mtime_raw{i}=jday2matlab(DATA.data(:,2),start_year)'; % convert julian dats to mtime
end

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

if ~isempty(start_time) % user specified start/end
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
if ~isempty(mtime{i})
    blank(i)=0;
    dt{i} = (mtime{i}(end)-mtime{i}(1))./length(mtime{i}); % data time step (days)
    temp{i} = sbe56.temp{i}(time_id{i});
else
    blank(i)=1;
    dt{i}=dt{i-1};
    temp{i} = nan+temp{i-1};
    mtime{i} = mtime{i-1};
end

% low pass filter temperature data
if blank(i)==0
    [temp_fft{i},FFT_temp{i},FFT_s{i}] = low_pass_filter(temp{i},dt{i}*24*3600,filter_time*60);
else
   temp_fft{i}=nan;
   FFT_temp{i}=nan;
   FFT_s{i}=nan;
end
    clear DATA indx
end

%% downsample filtered temp data to a common time
dt_filt=filter_time*60*(1/86400); % interpolate to 10 min window
mtime_filt=min(cell2mat(time_start)):dt_filt:max(cell2mat(time_end));

temp_filt=nan(length(mab),length(mtime_filt));
for i=1:length(mab)
    if blank(i)==0
        temp_filt(i,:)=interp1(mtime{i}, temp_fft{i},mtime_filt,'linear',NaN);
    else
        temp_filt(i,:) = nan+mtime_filt;
    end
end

%% clean up and save
clear sbe56 temp_fft DATA i ii sn_v start_year indx filename start_year temp_fft
save('Processed_SBE56_data.mat','-v7.3');

end

