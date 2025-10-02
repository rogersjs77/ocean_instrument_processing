% function []=sbe37_process(filter_time,start_time,end_time)
% Process SBE 37 CTD Data
% Justin Rogers, Stanford EFML, 8/2013
%
% this file takes SBE37 .cnv files, with 
% column 1 = Julian days
% column 2 = Tempuature (C)
% column 3 = Conductivity (S/m)
% column 4 = Salinity (psu)
% column 5 = Pressure (db)
% column 6 = O2 (mg/l)
% column 7 = O2 (% sat)
% file name format: 'SBE37_XXMAB_...cnv', where XX is vertical placement
%                   or 'SBE37...cnv' and MAB is assumed to be zero
% INPUT:
% filter_time in minutes
% start_time in matlab format (optional)
% end_time in matlab format (optional)

function []=sbe37_process(filter_time,start_time,end_time,headerlines)
g=9.81;

%% load data files
if nargin<4
    headerlines=250;
end

files=dir('SBE37*.cnv');
if isempty(files)
   disp('no files found, check format SBE37_XXMAB_...cnv')
   return
end

for i=1:size(files,1); % get order for files in increasing mab
filename=files(i).name;
m = sscanf(filename,'SBE37_%f MAB_%*f');

if isempty(m)
    mab(i)=0;
else
    mab(i)=m;   
end
end
[mab,ii]=sort(mab,2);
sbe37.mab=mab;

%% load data

for i=1:length(files);
filename=files(ii(i)).name;
DATA = importdata(filename,' ',headerlines);

sbe37.hdr{i}=DATA.textdata(:,1);

j=1;
start_year=[];
while isempty(start_year)
    start_year=sscanf(sbe37.hdr{i}{j},'# start_time = %*s %*f %f %*f');
    j=j+1;
end
sbe37.temp{i}=DATA.data(:,2)';
sbe37.cond{i}=DATA.data(:,3)';
sbe37.sal{i}=DATA.data(:,4)';
sbe37.press{i}=DATA.data(:,5)';
sbe37.press_Pa{i} = sbe37.press{i}*1E4;
sbe37.O2_pcnt{i} = DATA.data(:,7)';
sbe37.O2_mgl{i} = DATA.data(:,6)';

sbe37.mtime{i}=jday2matlab(DATA.data(:,1),start_year)'; % convert julian dats to mtime

if ~isempty(start_time) % user specified start/end
    time_start{i}=start_time;
    time_end{i} = end_time;        
else % find trim time from salinity within 1std of mean of middle 1/2 data
indx = find(sbe37.sal{i}>median(sbe37.sal{i}(round(0.25*end:0.75*end)))...
    -std(sbe37.sal{i}(round(0.25*end:0.75*end))));
% trim mtime to time underwater
time_start{i}=sbe37.mtime{i}(indx(1));
time_end{i}=sbe37.mtime{i}(indx(end));
end

time_id{i}= sbe37.mtime{i}>=time_start{i} & ...
    sbe37.mtime{i}<=time_end{i};

% time time and temp to limits
mtime{i}=sbe37.mtime{i}(time_id{i});

dt{i} = (mtime{i}(end)-mtime{i}(1))./length(mtime{i}); % data time step (days)
temp{i} = sbe37.temp{i}(time_id{i});
% temp{i}(indx)=nan;

cond{i} = sbe37.cond{i}(time_id{i});
sal{i} = sbe37.sal{i}(time_id{i});
press_Pa{i} = sbe37.press_Pa{i}(time_id{i});
O2_pcnt{i} = sbe37.O2_pcnt{i}(time_id{i});
O2_mgl{i} = sbe37.O2_mgl{i}(time_id{i});

[~,sigma] = swstate(sal{i},nanmean(temp{i}),nanmean(press_Pa{i}/1E4));
rho{i} = 1000+sigma;
depth{i} = press_Pa{i}./(rho{i}*g); % hydrostatic depth above bottom
Depth(i) = nanmean(depth{i});

% low pass filter temperature data
[temp_fft{i},FFT_temp{i},FFT_s{i}] = low_pass_filter(temp{i},dt{i}*24*3600,filter_time*60);
[cond_fft{i},FFT_cond{i},~] = low_pass_filter(cond{i},dt{i}*24*3600,filter_time*60);
[sal_fft{i},FFT_sal{i},~] = low_pass_filter(sal{i},dt{i}*24*3600,filter_time*60);
[depth_fft{i},FFT_depth{i},~] = low_pass_filter(depth{i},dt{i}*24*3600,filter_time*60);
[O2_mgl_fft{i},FFT_O2{i},~] = low_pass_filter(O2_mgl{i},dt{i}*24*3600,filter_time*60);
[O2_pcnt_fft{i},FFT_O2{i},~] = low_pass_filter(O2_pcnt{i},dt{i}*24*3600,filter_time*60);

clear DATA indx
end

%% interpolate filtered temp data to a common time
dt_filt=filter_time*60*(1/86400); % interpolate to 10 min window
mtime_filt=min(cell2mat(time_start)):dt_filt:max(cell2mat(time_end));

% temp_filt=nan(length(mab),length(mtime_filt));
for i=1:length(mab)
    temp_filt(i,:)=interp1(mtime{i}, temp_fft{i},mtime_filt,'linear',NaN);
    cond_filt(i,:)=interp1(mtime{i}, cond_fft{i},mtime_filt,'linear',NaN);
    sal_filt(i,:)=interp1(mtime{i}, sal_fft{i},mtime_filt,'linear',NaN);
    depth_filt(i,:)=interp1(mtime{i}, depth_fft{i},mtime_filt,'linear',NaN);
    eta_filt(i,:) = depth_filt(i,:)-nanmean(depth_filt(i,:));
    O2_pcnt_filt(i,:) = interp1(mtime{i}, O2_pcnt_fft{i},mtime_filt,'linear',NaN);
    O2_mgl_filt(i,:) = interp1(mtime{i}, O2_mgl_fft{i},mtime_filt,'linear',NaN);
end


    %% tidal analysis
    
[TIDES.nameu,TIDES.fu,TIDES.tidecon,eta_tides]=...
    t_tide(eta_filt(1,:),'interval',24*dt_filt); 

    eta_notides = eta_filt(i,:)-eta_tides;  

%% clean up and save

clear temp_fft DATA i ii sn_v start_year indx filename start_year temp_fft...
    cond_fft sal_fft press_fft depth_fft

save('Processed_SBE37_data.mat');

end

