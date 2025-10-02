% adp process and plot master script

clear, close all;

iprocess=1;
iplot=1;

if iprocess
    filter_time=30;
    P_atm=0;
    Waves_Spectra=0;
    Waves_Burst=1;
    nfft=400;
    fcutoff = 1/4;
    dirmethod = 1;
    theta=[]; % 
    head_height=0.6; 
    case_height = 0; % about 1ft from bed
    vertical_orientation='up';
    x_heading = 90; % x= East, deg CW from N
    time_start = [];%datenum(2013,9,11,18,0,0);
    time_end = [];%datenum(2014,7,14,8,00,0);
    extrap = 0;
    sample_size = 1;
    corr_min=45;
    % head rotated on 1/23/14, and then slowly in 3/14, vel data seems ok
    % before 3/2/14, after that corr values are very low (beam interference?).
    time_bad_vel = [];%[datenum(2014,3,2,00,00,00), datenum(2014,7,15,00,00,00)];
vector_process_data(filter_time,P_atm,Waves_Spectra,Waves_Burst,nfft,...
    fcutoff,theta,dirmethod,head_height,case_height,vertical_orientation,...
    x_heading,time_start,time_end,extrap,sample_size,corr_min,time_bad_vel)
end

if iplot
    vector_plot
end