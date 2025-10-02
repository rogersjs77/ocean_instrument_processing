% adp process and plot master script
close all, clear global

iprocess=1;
iplot=1;

if iprocess
    bottom_trim=1.5;
    depth_max=6;
    filter_time=15;
    P_atm= 220;
    Waves=1;
    WavesBurst = 1;
    nfft=200;
    fcutoff = 1/4;
    dirmethod = 1;
    theta=[]; % this adp seems to have the compass off by about +20 deg
    head_height=0.3048*3/12;
    extrap=0;
    z0 = 0.01;
    d = 0.01;
    time_start = [];%datenum(2016,7,25,10,00,0);
    time_end = [];%datenum(2016,8,5,15,00,0);
    % beams were blocked with debris after 10/14/13
    time_bad_vel = [];%[datenum(2013,10,14,12,0,0), datenum(2014,7,12,0,0,0)];
    bins_from_surf = 3;
aquadopp_process_data(bottom_trim,depth_max,filter_time,...
    P_atm,Waves,WavesBurst,nfft,fcutoff,dirmethod,theta,head_height,extrap,z0,d,...
    time_start,time_end,time_bad_vel,bins_from_surf)


end

if iplot
    aquadopp_plot
end