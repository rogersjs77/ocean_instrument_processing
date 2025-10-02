% function [] = vector_plot()
%
% this plots standard ADV data from vector_process_data
% 
% Justin Rogers, stanford efml,
% June 2017

function [] = vector_plot()
set(0,'defaulttextinterpreter','latex');
set(0,'defaulttextfontname','times');
set(0,'defaultaxesfontname','times');
set(0,'defaultaxesfontsize',12);
set(0,'defaulttextfontsize',12);

% run presentation_plot
close all; clear all;

load('ADV_PostProcess.mat');

%% scatter plots
p=figure(1);
subplot(1,2,1);
[a,b]= adcp_rotation([-1 1],[0 0], theta);
[c,d]= adcp_rotation([0 0],[-1 1], theta);

plot([-1 1],[0 0],'--k',[0 0],[-1 1],'--k');hold on;
plot(b,a,'--m',d,c,'--m');hold on;
scatter(Vel_East,Vel_North)
title('Depth-Avg Velocity');
xlabel('Depth-Avg East (U) (m/s)');
ylabel('Depth-Avg North (V) (m/s)');
axis square

subplot(1,2,2);
plot([-1 1],[0 0],'--k',[0 0],[-1 1],'--k');
hold on;
scatter(Vel_Maj,Vel_Min)
title('Depth Avg Rotated Velocity');
xlabel('Depth-Avg Major (U) (m/s)');
ylabel('Depth-Aveg Minor (V) (m/s)');
axis square
name=['adv_figure1_UV.fig'];
hgsave(p,name);

%% depth averaged earth coordinates
p=figure(2);
subplot(3,1,1); plot(mtime_filt,depth_filt)
datetick('x')
ylabel('depth filt (m)')
subplot(3,1,2); plot(mtime_filt,Vel_East)
datetick('x')
ylabel('U-east filt (m/s)')
subplot(3,1,3); plot(mtime_filt,Vel_North)
datetick('x')
ylabel('V-north filt (m/s)')
xlabel('date')
name=['adv_figure2_UV.fig'];
hgsave(p,name);

%% depth averaged rotated coordinates
p=figure(3);
subplot(3,1,1); plot(mtime_filt,Vel_Maj,mtime_filt,Vel_Maj_notides)
datetick('x')
ylabel('U-maj (m/s)')
legend('lp','detide')
subplot(3,1,2); plot(mtime_filt,Vel_Min,mtime_filt,Vel_Min_notides)
datetick('x')
ylabel('V-min (m/s)')
subplot(3,1,3); plot(mtime_filt,Vel_Vert)
datetick('x')
ylabel('W (m/s)')
xlabel('date')
name=['adv_figure3_UmVm.fig'];
hgsave(p,name);

%% plot correlation
p=figure(4);

subplot(4,1,1); plot(mtime_filt,corr_avg_filt)
datetick('x')
ylabel('avg corr filt')

subplot(4,1,2); plot(adv.mtime,adv.corr1)
datetick('x')
ylabel('corr 1')

subplot(4,1,3); plot(adv.mtime,adv.corr3)
datetick('x')
ylabel('corr 2')

subplot(4,1,4); plot(adv.mtime,adv.corr3)
datetick('x')
ylabel('corr 3')
xlabel('Date')

name=['adv_figure4_corr.fig'];
hgsave(p,name);

%% plot velocities

p=figure(5);

subplot(3,1,1); plot(adv.mtime,adv.vel_x)
datetick('x')
ylabel('raw u(x) (m/s)')

subplot(3,1,2); plot(adv.mtime,adv.vel_y)
datetick('x')
ylabel('raw v(y) (m/s)')

subplot(3,1,3); plot(adv.mtime,adv.vel_z)
datetick('x')
ylabel('raw w(z) (m/s)')
xlabel('Date')

name=['adv_figure5_xyz_raw.fig'];
hgsave(p,name);

%% plot velocities

p=figure(6);

subplot(3,1,1); plot(mtime,vel_east,mtime_filt,vel_east_filt)
datetick('x')
ylabel('u(E) (m/s)')

subplot(3,1,2); plot(mtime,vel_north,mtime_filt,vel_north_filt)
datetick('x')
ylabel('v(N) (m/s)')

subplot(3,1,3); plot(mtime,vel_vert,mtime_filt,vel_vert_filt)
datetick('x')
ylabel('w (m/s)')
xlabel('Date')

name=['adv_figure6_uvw.fig'];
hgsave(p,name);

%% temperature pressure depth

p=figure(7);
subplot(3,1,1);
var = temperature;
plot(mtime,var,mtime_filt,temperature_filt);
datetickzoom('x','keeplimits');
legend('raw','filtered');
ylabel('T (C)');
title('temperature pressure and depth');

subplot(3,1,2);
var = depth;
plot(mtime,var,mtime_filt,depth_filt);
datetickzoom('x','keeplimits');
ylabel('depth (m)');

subplot(3,1,3);
plot(mtime_filt,eta_2hrfilt,mtime_filt,eta_notides);
datetickzoom('x','keeplimits');
xlabel('date'); ylabel('$\eta$ (m)');
legend('2hr lp','detide')

name=['adv_figure7_temp_press_depth.fig'];
hgsave(p,name);

%% pressures

p=figure(8);
subplot(2,1,1);
var = adv.pressure*1E4;
plot(adv.mtime,var,adv.mtime,0*adv.mtime+P_atm);
datetickzoom('x','keeplimits');
legend('P raw','P atm');
ylabel('P (Pa)');
title('presure adjustment');

subplot(2,1,2);
var = adv.press_Pa;
plot(adv.mtime,var);
datetickzoom('x','keeplimits');
ylabel('P (Pa)');
legend('P adjusted')

name=['adv_figure8_pressure_adjustment.fig'];
hgsave(p,name);

%% heading pitch roll voltage
p=figure(9);
subplot(4,1,1);
var = adv.heading+0*adv.mtime;
plot(adv.mtime,var);
var = adv.heading_body+0*adv.mtime;
hold on, plot(adv.mtime,var);
% c = [-3*std(var) 3*std(var)]+mean(var);
% ylim(c);
datetick('x','keeplimits');
ylabel('heading (deg)');
title('heading pitch and roll');
legend('head','body')

subplot(4,1,2);
var = adv.pitch+0*adv.mtime;
plot(adv.mtime,var);
var = adv.pitch_body+0*adv.mtime;
hold on, plot(adv.mtime,var);
% c = [-3*std(var) 3*std(var)]+mean(var);
% ylim(c);
datetick('x','keeplimits');
ylabel('pitch (deg)');

subplot(4,1,3);
var = adv.roll+0*adv.mtime;
plot(adv.mtime,var);
var = adv.roll_body+0*adv.mtime;
hold on, plot(adv.mtime,var);
% c = [-3*std(var) 3*std(var)]+mean(var);
% ylim(c);
datetick('x','keeplimits');
ylabel('roll (deg)');

subplot(4,1,4);
var = adv.voltage;
plot(adv.mtime,var);
% c = [-3*std(var) 3*std(var)]+mean(var);
% ylim(c);
datetick('x','keeplimits');
xlabel('date'); ylabel('voltage');

name=['adv_figure9_head_ptich_roll_volt.fig'];
hgsave(p,name);


%% plot frequency space
p=figure(10);
subplot(3,1,1), 
plot(1./Fourier.s_vel/3600,2*abs(Fourier.Vel_Maj))
xlim([0 30])
ylabel('$2|F(Umaj)|$')
title('tidal amplitudes of depth avg velocity and depth');

subplot(3,1,2), 
plot(1./Fourier.s_vel/3600,2*abs(Fourier.Vel_Min))
xlim([0 30])
ylabel('$2|F(Umin)|$')

subplot(3,1,3), 
plot(1./Fourier.s_vel/3600,2*abs(Fourier.depth))
xlim([0 30])
ylabel('$2|F(d)|$')
xlabel('Period (hr)');
name=['adv_figure10_fftUVd.fig'];
hgsave(p,name);

%% plot power spectrum
p=figure(11);
M=round(length(Fourier.s_vel)/2)+2:length(Fourier.s_vel);
subplot(2,2,1), 
loglog(Fourier.s_vel,2*(abs(Fourier.Vel_Maj)).^2)
ylabel('$\Phi(Umaj)$')
xlabel('Hz');
title('power spectrum of depth avg velocity and depth');

subplot(2,2,2), 
loglog(Fourier.s_vel,2*(abs(Fourier.Vel_Min)).^2)
ylabel('$\Phi(Umin)$')
xlabel('Hz');

subplot(2,2,3), 
loglog(Fourier.s_vel,2*(abs(Fourier.depth)).^2)
axis tight
ylabel('$\Phi(d)$')
xlabel('Hz');

name=['adv_figure141_fftUVd.fig'];
hgsave(p,name);

if input.Waves_Spectra || input.Waves_Burst
% plot wave data
p=figure(20);
subplot(4,1,1),plot(WaveData.mtime,WaveData.Hsig_ss)
if input.Waves_Spectra
   hold on, plot(WaveData.mtime,WaveData.Hsig_ig),
end
ylabel('Hsig (m)'); datetick('x','keeplimits');
legend('swell','infragravity');
title('wave properties')

subplot(4,1,2), plot(WaveData.mtime,WaveData.Tm_ss);
ylabel('Tm (s)'); datetick('x','keeplimits');

subplot(4,1,3),plot(WaveData.mtime,WaveData.dir_heading_ss),
ylabel('dir going'); datetick('x','keeplimits');

subplot(4,1,4),plot(WaveData.mtime,WaveData.spread_ss),
ylabel('spread (deg)'); datetick('x','keeplimits');

hgsave(p,'adv_figure20_Wave_H_dir_sprd.fig');
end
%
if input.Waves_Spectra || input.Waves_Burst
p=figure(21);
subplot(3,1,1),imagesc(WaveData.mtime,WaveData.fmt,WaveData.SSE')
set(gca,'YDir','normal'); 
q=colorbar; xlabel(q,'Snn');
ylabel('f [Hz]'), datetick('x','keeplimits')
title('wave power spectra, direction and spread');

subplot(3,1,2),imagesc(WaveData.mtime,WaveData.fmt,WaveData.DIR_heading')
set(gca,'YDir','normal'); %caxis([-30 30]);
q=colorbar; xlabel(q,'Dir going');
ylabel('f [Hz]'), datetick('x','keeplimits')

subplot(3,1,3),imagesc(WaveData.mtime,WaveData.fmt,WaveData.SPREAD')
set(gca,'YDir','normal'); 
q=colorbar; xlabel(q,'Sprd');
ylabel('f [Hz]'), datetick('x','keeplimits')
hgsave(p,'adv_figure21_SSE_DIR_SPREAD.fig')

%%
p=figure(22);
subplot(3,1,1),imagesc(WaveData.mtime,WaveData.fmt,WaveData.SSE')
set(gca,'YDir','normal'); 
q=colorbar; xlabel(q,'Snn');
ylabel('f [Hz]'), datetick('x','keeplimits')
title('Wave power spectra, U Stoke, V Stoke');

subplot(3,1,2),imagesc(WaveData.mtime,WaveData.fmt,WaveData.USt')
set(gca,'YDir','normal'); 
q=colorbar; xlabel(q,'Ustoke');
ylabel('f [Hz]'), datetick('x','keeplimits')

subplot(3,1,3),imagesc(WaveData.mtime,WaveData.fmt,WaveData.VSt')
set(gca,'YDir','normal'); 
q=colorbar; xlabel(q,'Vstoke');
ylabel('f [Hz]'), datetick('x','keeplimits')

hgsave(p,'adv_figure22_SSE_US_VU.fig');
end
%%
if input.Waves_Spectra || input.Waves_Burst
p=figure(23);
plot(WaveData.mtime,WaveData.Us,...
    WaveData.mtime,WaveData.Vs)
ylabel('Stokes Drift (m/s)'); datetick('x','keeplimits');
legend('Us','Vs')    
hgsave(p,'adv_figure23_Stokes_Drift.fig');
end
%%
if input.Waves_Spectra || input.Waves_Burst
p=figure(24);
subplot(2,1,1), plot(WaveData.mtime,WaveData.ztest_ss)
ylabel('ztest'), datetick('x','keeplimits')
subplot(2,1,2), plot(WaveData.mtime,WaveData.Cpu_ss)
ylabel('Cpu'), datetick('x','keeplimits')
hgsave(p,'adv_figure24_ztest_Cpu.fig');
end

if input.Waves_Spectra || input.Waves_Burst
%% lagrangian flows
p=figure(25);
plot([-1 1],[0 0],'--k',[0 0],[-1 1],'--k'); 
hold on;
scatter(Vel_East_Lagr,Vel_North_Lagr)
title('Depth-Averaged Lagrangian Velocity');
xlabel('Depth-Avg East (Ul) (m/s)');
ylabel('Depth-Avg North (Vl) (m/s)');
axis([-0.5 0.5 -0.5 0.5]);
name=['adv_figure25_UV.fig'];
hgsave(p,name);

%% lagrangian vel depth averaged earth coordinates
p=figure(26);
subplot(2,1,1); plot(mtime_filt,Vel_East_Lagr,mtime_filt,Vel_East_Stokes,...
    mtime_filt,Vel_East)
datetick('x')
legend('Ul','Us','Ue')
ylabel('Ul-east filt (m/s)')
subplot(2,1,2); plot(mtime_filt,Vel_North_Lagr,mtime_filt,Vel_North_Stokes,...
    mtime_filt,Vel_North)
datetick('x')
ylabel('Vl-north filt (m/s)')
xlabel('date')
name=['adv_figure26_UV.fig'];
hgsave(p,name);


end
if input.Waves_Spectra || input.Waves_Burst
%% RMS velocitites
p = figure(27);
plot(WaveData.mtime,WaveData.up_rms,...
    WaveData.mtime,WaveData.vp_rms,...
    WaveData.mtime,WaveData.wp_rms,...
    WaveData.mtime,WaveData.uvwp_rms)
legend('maj','min','vert','mag')
datetickzoom('x')
ylabel('$(u'')_{rms}$')
hgsave(p,'adv_figure27_up_RMS')

%% bottom stress

p=figure(28);
plot(WaveData.mtime,WaveData.Tau_bx,WaveData.mtime,WaveData.Tau_by,...
    WaveData.mtime,WaveData.Tau)
ylabel('$\tau$ (Pa)'), datetick('x','keeplimits')
legend('\tau_{bx}','\tau_{by}','|\tau_{b}|')

hgsave(p,'adv_figure28_bottomT.fig');
end

end