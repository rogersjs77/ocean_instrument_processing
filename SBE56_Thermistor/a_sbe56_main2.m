% master file for sbe56
clear all, close all;

iprocess=1;
iplot=1;

DATA = importdata('./a_termistor_times.xlsx');
start_time = DATA.data(:,2)+datenum(1900,1,1)-2;
end_time = DATA.data(:,3)+datenum(1900,1,1)-2;

%     start_time = datenum(2017,3,10,18,0,0)*ones(32,1);
%     end_time = datenum(2017,3,28,10,30,0)*ones(32,1);

if iprocess
    filter_time=15;

    sbe56_process_multiple_sites(filter_time,start_time,end_time)  

end

if iplot
    sbe56_plot    
end