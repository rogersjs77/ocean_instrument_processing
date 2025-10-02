% function [] = rbr1050_plot()
%
% this plots standard RBR data from rbr1050_process
% 
% Justin Rogers, stanford efml, 10/2013

function []=rbr1050_plot()
set(0,'defaulttextinterpreter','latex');
set(0,'defaulttextfontname','times');
set(0,'defaultaxesfontname','times');
set(0,'defaultaxesfontsize',12);
set(0,'defaulttextfontsize',12);

close all

load Processed_RBR1050_data

%% basic plot
p=figure(1);
subplot(2,1,1);
var = RBR.pressure*1E4;
plot(RBR.mtime,var,RBR.mtime,0*RBR.mtime+P_atm);
datetickzoom('x','keeplimits');
legend('P raw','P atm');
ylabel('P (Pa)');
title('presure adjustment');

subplot(2,1,2);
var = RBR.press_Pa;
plot(RBR.mtime,var);
datetickzoom('x','keeplimits');
ylabel('P (Pa)');
legend('P adjusted')

hgsave('figure1_raw_ph_data')

%%
figure(2)
subplot(3,1,1),plot(mtime_filt,pressure_filt)
datetickzoom('x')
ylabel('P (Pa)');
title('Filtered/Trimmed RBR1050 Data')

subplot(3,1,2),plot(mtime_filt,depth_filt)
datetickzoom('x')
ylabel('depth (m)');

subplot(3,1,3), plot(mtime_filt,eta_filt,mtime_filt,eta_notides)
datetickzoom('x')
ylabel('$\eta$ (m)')
legend('lp','detide')

hgsave('figure2_filtered_data')

%% spectra
figure(3)
subplot(2,1,1),plot(1./Fourier.s_vel/3600,2*abs(Fourier.depth))
xlim([0 30])
ylabel('$2|F(d)|$')
xlabel('Period (hr)');

subplot(2,1,2)
loglog(Fourier.s_vel,2*(abs(Fourier.depth)).^2)
axis tight
ylabel('$\Phi(d)$')
xlabel('Hz');

hgsave('figure3_spectra');

%% Wave plots
if Waves==1    
     %% plot swell band
    figure(4), 
    subplot(3,1,1), plot(WaveData.mtime,WaveData.Hsig_ss),
    ylabel('Hsig [m]'); datetick('x','keeplimits');
    title('swell wave properties')
    
    subplot(3,1,2), plot(WaveData.mtime,WaveData.Tm_ss), hold all;
    datetick('x','keeplimits');
    ylabel('Tm (s)');    
    
    subplot(3,1,3), plot(WaveData.mtime,WaveData.Us), hold all;
    datetick('x','keeplimits');
    ylabel('Us (m/s)');
    xlabel('time')
    
    hgsave('figure4_swell_band.fig');
    
    %% ig band
    figure(5), 
    subplot(2,1,1), plot(WaveData.mtime,WaveData.Hsig_ig), hold all;
    datetickzoom('x');
    ylabel('H (m)');
    title('infragravity wave properties')
    
    subplot(2,1,2), plot(WaveData.mtime,WaveData.Tm_ig), hold all;
    datetickzoom('x');
    ylabel('Tm (s)');    

    hgsave('figure5_ig_band.fig');
    %%
    
    p=figure(6);
    imagesc(WaveData.mtime,WaveData.fmt,WaveData.SSE')
    set(gca,'YDir','normal'); 
    q=colorbar; xlabel(q,'Snn'); %caxis([0 5])
    ylabel('f [Hz]'), datetick('x','keeplimits')
    title('wave power spectra');

    hgsave(p,'figure6_power_spectra.fig')
%%    
    p=figure(7);
    plot(WaveData.mtime,WaveData.ztest_ss)
    ylabel('ztest'), datetick('x','keeplimits')
    hgsave(p,'figure7_ztest_Cpu.fig');
    
        %% ig band
    figure(8), 
    subplot(2,1,1), plot(WaveData.mtime,WaveData.Hrms_ss,...
        WaveData.mtime,WaveData.Hrms_all,...
        WaveData.mtime,WaveData.Hrms_ig), hold all;
    datetickzoom('x');
    ylabel('Hrms (m)');
    title('all wave properties')
    legend('swell','all','ig')
    
    subplot(2,1,2), plot(WaveData.mtime,WaveData.Tm_ss,...
        WaveData.mtime,WaveData.Tm_all,...
        WaveData.mtime,WaveData.Tm_ig), hold all;
    datetickzoom('x');
    ylabel('Tm (s)');    

    hgsave('figure8_allband.fig');

%%
end
  
    
end