% function []=sbe19_process()
% Process SBE 19 CTD Data
% Justin Rogers, Stanford EFML, 8/2013
%
% this file takes SBE19 .cnv files, with 
% column 1 = depth
% column 2 = Tempuature
% column 3 = Salinity
% column 4 = oxygen
% file name format: 'SBE19_NAME_...cnv', 

function []=sbe19_process()

%% load data files

files=dir('*.cnv');

for i=1:size(files,1); % get order for files in increasing mab
filename=files(i).name;
m = textscan(filename,'SBE19plus_%s_');

if isempty(m)
    name(i)='';
else
    name(i)=m;   
end
end

sbe19.name=name;

%% load data

for i=1:length(files);
filename=files(i).name;
DATA = importdata(filename,' ',307);

sbe19.hdr{i}=DATA.textdata(:,1);

j=1;
timescan=[];
while isempty(timescan)
%     # datcnv_date = Sep 07 2013 20:09:46, 7.22.5 [datcnv_vars = 4]
    timescan=sscanf(sbe19.hdr{i}{j},'# datcnv_date = %*s %f %f %f:%f:%f,');
    mo=sscanf(sbe19.hdr{i}{j},'# datcnv_date = %s %*f %*f %*f:%*f:%*f,');
    j=j+1;
end

% cast_time(i) = datenum(timescan(2),month(['01-',char(mo'),'-2013']),...
%     timescan(1),timescan(3),timescan(4),timescan(5));
% # name 0 = depSM: Depth [salt water, m], lat = 0.00
% # name 1 = sal00: Salinity, Practical [PSU]
% # name 2 = tv290C: Temperature [ITS-90, deg C]
% # name 3 = sbeox0Mg/L: Oxygen, SBE 43 [mg/l]
% sbe19.depth{i}=DATA.data(:,1)';
% sbe19.sal{i}=DATA.data(:,2)';
% sbe19.temp{i}=DATA.data(:,3)';
% sbe19.oxygen{i}=DATA.data(:,4)';

% # name 0 = temp
% # name 1 = cond
% # name 2 = sal
% # name 3 = depth
% # name 4 = flag
sbe19.depth{i}=DATA.data(:,4)';
sbe19.sal{i}=DATA.data(:,3)';
sbe19.temp{i}=DATA.data(:,1)';
sbe19.oxygen{i}=DATA.data(:,5)';

indx = sbe19.depth{i}>1; % trim to time in water
depth{i} = sbe19.depth{i}(indx);
sal{i} = sbe19.sal{i}(indx);
temp{i} = sbe19.temp{i}(indx);
oxygen{i} = sbe19.oxygen{i}(indx);

indx_dn=0*depth{1};
jj=1;
dz=1;
while dz>0
   indx_dn(jj)=1;
   dz=mean(diff(depth{i}(jj:jj+5)));
   jj=jj+1;   
end
indx_dn=logical(indx_dn);
indx_up = logical(-indx_dn+1);

depth_dn{i} = depth{i}(indx_dn);
temp_dn{i} = temp{i}(indx_dn);
sal_dn{i}=sal{i}(indx_dn);
oxygen_dn{i}=oxygen{i}(indx_dn);

end

%% clean up and save
clear j jj i indx dz m mo timescan filename DATA
save('Processed_SBE19_data.mat');


end

