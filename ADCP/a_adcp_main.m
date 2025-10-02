% ADCP master 
clear, close all

iprocess=1;
ipost_process=1;
iplot=1;

if iprocess
   adcp_process_rdrADCP 
end
close 
if ipost_process
   bottom_trim = 1.6; 
   depth_max = 14.5;
   filter_time = 30;
   bad_beam=4;
   P_atm = [];%1.251E4;
   Waves = 0;
   nfft = 100;
   fcutoff = 1/5.5;
   dirmethod = 2;
   corr_min = 0;
   theta = 6.0941+180; % force positive onshore flow
   head_height = 0.3048*16.75/12;
   extrap=0;
   z0 = 0.01;
   d = 0.01;
   time_start = [];%datenum(2012,9,16,11,1,0);
   time_end = [];%datenum(2012,9,26,10,50,0);
   sample_size = 1;
   time_bad_vel = [];
adcp_post_process(bottom_trim,depth_max,filter_time,...
    bad_beam,P_atm,Waves,nfft,fcutoff,dirmethod,corr_min,theta,head_height,...
    extrap,z0,d,time_start,time_end,sample_size,time_bad_vel)
end

if iplot
   adcp_plot;
end