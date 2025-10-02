% master file for rbr1050
clear, close all;

iprocess=1;
iplot=1;

if iprocess
    gauge_ht = 0*0.3048;
    filter_time=15;
    P_atm = 1.035E5;
    Waves=1;
    nfft=256*2;
    fcutoff = 1/4;
    time_start = [];%datenum(2013,9,3,16,30,0);
    time_end = [];%datenum(2013,9,8,10,00,0);
    rbr_process(gauge_ht,filter_time,P_atm,Waves,nfft,fcutoff,...
    time_start,time_end)
    
end

if iplot
    rbr_plot    
end