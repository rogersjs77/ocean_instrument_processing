% master file for sbe37
clear, close all;

process=1;
plot=1;

if process
    filter_time=30;
    start_time = datenum(2013,9,9,12,0,0);
    end_time = datenum(2014,5,23,1,0,0);
    sbe37_process(filter_time,start_time,end_time)  

end

if plot
    sbe37_plot    
end