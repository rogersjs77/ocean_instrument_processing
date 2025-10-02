% this file processes .000 ADCP data using rdradcp process,
% Justin Rogers, stanford efml, 2/2013

function [] = process_rdrADCP()
% run presentation_plot
clear all, close all;

%% process ADCP files

files=dir('*.000');    
for ii=1:size(files,1);

if files(ii).bytes>0 % only process files with data
filename=files(ii).name;  
%[..]=RDRADCP(NAME,NUMAV,NENS,'despike',[ 'no' | 'yes' | 3-element vector ])
% [adcp,cfg,ens,hdr]=rdradcp(filename,1); 
[adcp,cfg,ens,hdr]=rdradcp(filename,1,-1,'despike','yes'); 
if isempty(adcp)==0 % check for no data output
    
if size(adcp.mtime,2)>2 % only process files with sufficient time data

clear i j ii
matname = ['rdrADCP_',filename,'.mat'];
save(matname,'-v7.3');

end
end
end
end

