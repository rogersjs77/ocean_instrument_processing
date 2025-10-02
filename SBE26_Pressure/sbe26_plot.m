% function [] = SBE26_plot()
%
% this plots standard RBR data from SBE26_process
% 
% Justin Rogers, stanford efml, 10/2013

function []=sbe26_plot()
set(0,'defaulttextinterpreter','latex');
set(0,'defaulttextfontname','times');
set(0,'defaultaxesfontname','times');
set(0,'defaultaxesfontsize',12);
set(0,'defaulttextfontsize',12);

close all

load Processed_SBE26_data

%% basic plot
p=figure(1);
subplot(2,1,1);
var = SBE.pressure;
plot(SBE.mtime,var,SBE.mtime,0*SBE.mtime+P_atm);
datetickzoom('x','keeplimits');
legend('P raw','P atm');
ylabel('P (Pa)');
title('presure adjustment');

subplot(2,1,2);
var = SBE.press_Pa;
plot(SBE.mtime,var);
datetickzoom('x','keeplimits');
ylabel('P (Pa)');
legend('P adjusted')

% hgsave('figure1_raw_ph_data')
print -djpeg -r300 figure1_raw_ph_data
%%
figure(2)
subplot(3,1,1),plot(mtime_filt,pressure_filt)
datetickzoom('x')
ylabel('P (Pa)');
title('Filtered/Trimmed SBE26 Data')

subplot(3,1,2),plot(mtime_filt,depth_filt)
datetickzoom('x')
ylabel('depth (m)');

subplot(3,1,3), plot(mtime_filt,eta_filt,mtime_filt,eta_notides)
datetickzoom('x')
ylabel('$\eta$ (m)')
legend('lp','detide')

% hgsave('figure2_filtered_data')
print -djpeg -r300 figure2_filtered_data
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

% hgsave('figure3_spectra');
print -djpeg -r300 figure3_spectra
%% Wave plots
if Waves==1    
     %% plot swell band
    figure(4), 
    subplot(2,1,1), plot(WaveData.mtime,WaveData.Hsig_ss,WaveData.mtime,WaveData.Hsig_ig),
    ylabel('Hsig [m]'); datetick('x','keeplimits');
    title('swell wave properties')
    legend('sw','ig')
    
    subplot(2,1,2), plot(WaveData.mtime,WaveData.Tm_ss), hold all;
    datetick('x','keeplimits');
    ylabel('Tm (s)');    
    
      
%     hgsave('figure4_swell_band.fig');
    print -djpeg -r300 figure4_swell_band
   p=figure(6);
    z = WaveData.SSE';
    zz = reshape(z,[],1);
    c1 = [0 nanmean(zz)+2*nanstd(zz)];
    pcolorjw(WaveData.mtime,WaveData.fmt,WaveData.SSE') 
    q=colorbar; xlabel(q,'S_{\eta\eta}'); 
    caxis(c1)
    ylabel('f [Hz]'), datetick('x','keeplimits')
    title('wave power spectra');

%     hgsave(p,'figure6_power_spectra.fig')
    print -djpeg -r300 figure6_power_spectra
    %% avg power spectra
    
    p=figure(7);
    x = nanmean(WaveData.SSE,1);
     xb = nanstd(WaveData.SSE,1);
     plot(WaveData.fmt,x,'-xk');
     hold on
     plot(WaveData.fmt,x+2*xb,'-b')
     legend('avg','avg+2\sigma')
     
   
     xlabel('f (Hz)')
     ylabel('Avg SSE')
%     hgsave(p,'figure6_power_spectra.fig')
    print -djpeg -r300 figure6_avg_power_spectra
%%    

    
 
%%
end
  
    
end