% function [] = aquadopp_plot()
%
% this plots standard ADP data from aquadopp_process_data
% 
% Justin Rogers, stanford efml,
% June 2017

function [] = aquadopp_plot()
set(0,'defaulttextinterpreter','latex');
set(0,'defaulttextfontname','times');
set(0,'defaultaxesfontname','times');
set(0,'defaultaxesfontsize',12);
set(0,'defaulttextfontsize',12);

% run presentation_plot
close all; clear all;



load('ADP_PostProcess.mat');

%% scatter plots
p=figure(1);
p=figure(1);
subplot(1,2,1);
[a,b]= adcp_rotation([-1 1],[0 0], theta);
[c,d]= adcp_rotation([0 0],[-1 1], theta);
c1 = max(reshape(sqrt(Vel_East.^2+Vel_North.^2),[],1));

plot([-1 1],[0 0],'--k',[0 0],[-1 1],'--k');hold on;
plot(b,a,'--m',d,c,'--g');hold on;
plot(Vel_East,Vel_North,'.')
hold on
plot(nanmean(Vel_East),nanmean(Vel_North),'sk')
xlim(c1*[-1 1])
ylim(c1*[-1 1])
title('Depth-Avg Velocity');
xlabel('Depth-Avg East (U) (m/s)');
ylabel('Depth-Avg North (V) (m/s)');
axis square

subplot(1,2,2);
plot([-1 1],[0 0],'--g',[0 0],[-1 1],'--m');
hold on;
plot(Vel_Maj,Vel_Min,'.')
hold on
plot(nanmean(Vel_Maj),nanmean(Vel_Min),'sk')
xlim(c1*[-1 1])
ylim(c1*[-1 1])
title('Depth Avg Rotated Velocity');
xlabel('Depth-Avg Major (U) (m/s)');
ylabel('Depth-Aveg Minor (V) (m/s)');
axis square
% name=['adcp_figure1_UV.fig'];
name=['adp_figure1_UV.fig'];
hgsave(p,name);

%% plot vertical profiles
p=figure(2);
subplot(4,1,1),
var = vel_maj_filt;
pcolorjw(mtime_filt,bin_MAB,var); colormap('french');shading flat;
set(gca,'ydir','normal');
hold on; plot(mtime_filt,depth_filt,'-k');
ylim([min(bin_MAB) max(depth_filt)]);
q=colorbar; %
c = -2*nanmean(nanstd(var,0,2))-abs(nanmean(nanmean(var)));
d = 2*nanmean(nanstd(var,0,2))+abs(nanmean(nanmean(var)));
caxis([c d]) % set limits to 2 stdev of mean
datetickzoom('x');
ylabel(q,'u maj (m/s)'); 
ylabel('MAB (m)');
title('Velocity Profiles (filtered)');

subplot(4,1,2),
var = vel_min_filt;
pcolorjw(mtime_filt,bin_MAB,var); colormap('french');shading flat;
set(gca,'ydir','normal');
hold on; plot(mtime_filt,depth_filt,'-k');
q=colorbar; 
c = -2*nanmean(nanstd(var,0,2))-abs(nanmean(nanmean(var)));
d = 2*nanmean(nanstd(var,0,2))+abs(nanmean(nanmean(var)));
caxis([c d]) % set limits to 2 stdev of mean
datetickzoom('x');
ylabel(q,'u min (m/s)'); 
ylabel('MAB (m)');

subplot(4,1,3),
var = vel_vert_filt;
pcolorjw(mtime_filt,bin_MAB,var); colormap('french');shading flat;
set(gca,'ydir','normal');
hold on; plot(mtime_filt,depth_filt,'-k');
q=colorbar; 
c = -2*nanmean(nanstd(var,0,2))-abs(nanmean(nanmean(var)));
d = 2*nanmean(nanstd(var,0,2))+abs(nanmean(nanmean(var)));
caxis([c d]) % set limits to 2 stdev of mean
datetickzoom('x');
ylabel(q,'w (m/s)'); 

subplot(4,1,4),
var = vel_err_filt;
pcolorjw(mtime_filt,bin_MAB,var); colormap('french');shading flat;
set(gca,'ydir','normal');
q=colorbar; 
c = -2*nanmean(nanstd(var,0,2))-abs(nanmean(nanmean(var)));
d = 2*nanmean(nanstd(var,0,2))+abs(nanmean(nanmean(var)));
caxis([c d]) % set limits to 2 stdev of mean
datetickzoom('x');
ylabel(q,'err (m/s)'); 

xlabel('Date'); ylabel('MAB (m)');


name=['adp_figure2_umvmwm.fig'];
hgsave(p,name);

%% depth averaged earth coordinates
p=figure(3);
subplot(2,1,1); plot(mtime_filt,Vel_East)
datetick('x')
ylabel('U-east filt (m/s)')
subplot(2,1,2); plot(mtime_filt,Vel_North)
datetick('x')
ylabel('V-north filt (m/s)')
xlabel('date')
name=['adp_figure3_UV.fig'];
hgsave(p,name);

%% depth averaged rotated coordinates
p=figure(4);
subplot(2,1,1); plot(mtime_filt,Vel_Maj)
if mtime(end)-mtime(1)>2
    hold all,plot(mtime_filt,Vel_Maj_notides)
end
datetick('x')
ylabel('U-major x filt (m/s)')
subplot(2,1,2); plot(mtime_filt,Vel_Min)
if mtime(end)-mtime(1)>2
    hold all,plot(mtime_filt,Vel_Min_notides)
end
datetick('x')
ylabel('V-minor y filt (m/s)')
xlabel('date')
legend('lp','detide')
name=['adp_figure4_UmVm.fig'];
hgsave(p,name);

%% plot intensity
ttl={'bm 1 amp','bm 2 amp','bm 3 amp'};

p=figure(5);
for j=1:3
q=subplot(3,1,j); 
if j==1
    var=adp.amp1;
    pcolorjw(adp.mtime,bin_MAB_orig,var);
    c = 2*nanmean(nanstd(var,0,2));
    caxis([-c c]+nanmean(nanmean(var))), qq=colorbar('peer',q); % set limits to 2 stdev    
    hold all, plot([mtime(1) mtime(end)],bin_MAB(1)*[1 1],'k') %#ok<*COLND>
    hold all, plot([mtime(1) mtime(end)],bin_MAB(end)*[1 1],'k')
    xlabel(qq,[ttl{j}]);

elseif j==2
    var=adp.amp2;
    pcolorjw(adp.mtime,bin_MAB_orig,var);set(gca,'YDir','normal'),
    c = 2*nanmean(nanstd(var,0,2));
    caxis([-c c]+nanmean(nanmean(var))), qq=colorbar('peer',q); % set limits to 2 stdev    
    hold all, plot([mtime(1) mtime(end)],bin_MAB(1)*[1 1],'k')
    hold all, plot([mtime(1) mtime(end)],bin_MAB(end)*[1 1],'k')
    xlabel(qq,[ttl{j}]);

elseif j==3
    var=adp.amp3;
    pcolorjw(adp.mtime,bin_MAB_orig,var);set(gca,'YDir','normal'),
    c = 2*nanmean(nanstd(var,0,2));
    caxis([-c c]+nanmean(nanmean(var))), qq=colorbar('peer',q); % set limits to 2 stdev    
    hold all, plot([mtime(1) mtime(end)],bin_MAB(1)*[1 1],'k')
    hold all, plot([mtime(1) mtime(end)],bin_MAB(end)*[1 1],'k')
    xlabel(qq,[ttl{j}]);


end

hold all, plot(mtime_filt,depth_filt,'k')
datetick('x','keeplimits')
ylabel('z (m)'); 

end
xlabel('time');
name=['adp_figure5_intens.fig'];
hgsave(p,name);


%% plot velocities
ttl={'east vel (m/s)','north vel (m/s)','vert vel (m/s)','error vel (m/s)'};

p=figure(6);
for j=1:3
q=subplot(3,1,j); 
if j==1
    var=adp.vel_east;
    pcolorjw(adp.mtime,bin_MAB_orig,var);set(gca,'YDir','normal'),
    c = 2*nanmean(nanstd(var,0,2));
    caxis([-c c]), qq=colorbar('peer',q); % set limits to 2 stdev    
    hold all, plot([mtime(1) mtime(end)],bin_MAB(1)*[1 1],'k')
    hold all, plot([mtime(1) mtime(end)],bin_MAB(end)*[1 1],'k')
    xlabel(qq,[ttl{j}]);

elseif j==2
    var=adp.vel_north;
    pcolorjw(adp.mtime,bin_MAB_orig,var);set(gca,'YDir','normal'),
    c = 2*nanmean(nanstd(var,0,2));
    caxis([-c c]),qq= colorbar('peer',q); % set limits to 2 stdev    
    hold all, plot([mtime(1) mtime(end)],bin_MAB(1)*[1 1],'k')
    hold all, plot([mtime(1) mtime(end)],bin_MAB(end)*[1 1],'k')
    xlabel(qq,[ttl{j}]);

elseif j==3
    var=adp.vel_vert;
    pcolorjw(adp.mtime,bin_MAB_orig,var);set(gca,'YDir','normal'),
    c = 2*nanmean(nanstd(var,0,2));
    caxis([-c c]), qq=colorbar('peer',q); % set limits to 2 stdev    
    hold all, plot([mtime(1) mtime(end)],bin_MAB(1)*[1 1],'k')
    hold all, plot([mtime(1) mtime(end)],bin_MAB(end)*[1 1],'k')
    xlabel(qq,[ttl{j}]);


end

hold all, plot(mtime_filt,depth_filt,'k')
datetick('x','keeplimits')
ylabel('z (m)'); 

end
xlabel('time');
name=['adp_figure6_uvwerr.fig'];
hgsave(p,name);

%% temperature pressure depth

p=figure(7);
subplot(3,1,1);
var = temperature;
plot(mtime,var,mtime_filt,temperature_filt);
datetickzoom('x','keeplimits');
legend('raw','filtered');
ylabel('T (C)');
title('figure 10, temperature pressure and depth');

subplot(3,1,2);
var = depth;
plot(mtime,var,mtime_filt,depth_filt);
datetickzoom('x','keeplimits');
ylabel('depth (m)');

subplot(3,1,3);
plot(mtime_filt,eta_filt)
if mtime(end)-mtime(1)>2
    hold all,plot(mtime_filt,eta_notides)
end
datetickzoom('x','keeplimits');
xlabel('date'); ylabel('$\eta$ (m)');
legend('2hr lp','detide')

name=['adp_figure7_temp_press_depth.fig'];
hgsave(p,name);

%% pressures

p=figure(8);
subplot(2,1,1);
var = adp.pressure*1E4;
plot(adp.mtime,var,adp.mtime,0*adp.mtime+P_atm);
datetickzoom('x','keeplimits');
legend('P raw','P atm');
ylabel('P (Pa)');
title('figure 10, presure adjustment');

subplot(2,1,2);
var = adp.press_Pa;
plot(adp.mtime,var);
datetickzoom('x','keeplimits');
ylabel('P (Pa)');
legend('P adjusted')

name=['adp_figure8_pressure_adjustment.fig'];
hgsave(p,name);

%% heading pitch roll

p=figure(9);
subplot(4,1,1);
var = heading;
plot(mtime,var);
c = [-3*std(var) 3*std(var)]+mean(var);
ylim(c);
datetickzoom('x','keeplimits');
ylabel('heading (deg)');
title('figure 11, heading pitch and roll');

subplot(4,1,2);
var = pitch;
plot(mtime,var);
c = [-3*std(var) 3*std(var)]+mean(var);
ylim(c);
datetickzoom('x','keeplimits');
ylabel('pitch (deg)');

subplot(4,1,3);
var = roll;
plot(mtime,var);
c = [-3*std(var) 3*std(var)]+mean(var);
ylim(c);
datetickzoom('x','keeplimits');
ylabel('roll (deg)');

subplot(4,1,4);
var = voltage;
plot(mtime,var);
c = [-3*std(var) 3*std(var)]+mean(var);
ylim(c);
datetickzoom('x','keeplimits');
xlabel('date'); ylabel('voltage');

name=['adp_figure9_head_ptich_roll_volt.fig'];
hgsave(p,name);

%% average profiles

p=figure(10);
subplot(1,2,1), plot(Vel_Maj_Prof,bin_MAB);
xlabel('U Maj (m/s)'), ylabel('MAB');
legend('pos','neg');

subplot(1,2,2), plot(Vel_Min_Prof,bin_MAB);
xlabel('U Min (m/s)'), ylabel('MAB');

hgsave(p,'adp_figure10_profile');


%% plot frequency space
p=figure(11);
subplot(3,1,1), 
% plot(1./Fourier.s_vel/3600,2*abs(Fourier.Vel_Maj))
xlim([0 30])
ylabel('$2|F(Umaj)|$')
title('tidal amplitudes of depth avg velocity and depth');

subplot(3,1,2), 
% plot(1./Fourier.s_vel/3600,2*abs(Fourier.Vel_Min))
xlim([0 30])
ylabel('$2|F(Umin)|$')

subplot(3,1,3), 
% plot(1./Fourier.s_vel/3600,2*abs(Fourier.depth))
xlim([0 30])
ylabel('$2|F(d)|$')
xlabel('Period (hr)');
name=['adp_figure11_fftUVd.fig'];
hgsave(p,name);

%% plot power spectrum
% p=figure(12);
% M=round(length(Fourier.s_vel)/2)+2:length(Fourier.s_vel);
% subplot(2,2,1), 
% loglog(Fourier.s_vel,2*(abs(Fourier.Vel_Maj)).^2)
% ylabel('$\Phi(Umaj)$')
% xlabel('Hz');
% title('power spectrum of depth avg velocity and depth');
% 
% subplot(2,2,2), 
% loglog(Fourier.s_vel,2*(abs(Fourier.Vel_Min)).^2)
% ylabel('$\Phi(Umin)$')
% xlabel('Hz');
% 
% subplot(2,2,3), 
% loglog(Fourier.s_vel,2*(abs(Fourier.depth)).^2)
% axis tight
% ylabel('$\Phi(d)$')
% xlabel('Hz');
% 
% name=['adp_figure12_fftUVd.fig'];
% hgsave(p,name);

if input.Waves
% plot wave data
p=figure(20);
subplot(4,1,1),plot(WaveData.mtime,WaveData.Hsig_ss,...
    WaveData.mtime,WaveData.Hsig_ig),
ylabel('Hsig (m)'); datetick('x','keeplimits');
legend('swell','infragravity');
title('wave properties')

subplot(4,1,2), plot(WaveData.mtime,WaveData.Tm_ss);
ylabel('Tm (s)'); datetick('x','keeplimits');

subplot(4,1,3),plot(WaveData.mtime,WaveData.dir_heading_ss),
ylabel('Heading dir going '); datetick('x','keeplimits');

subplot(4,1,4),plot(WaveData.mtime,WaveData.spread_ss),
ylabel('spread (deg)'); datetick('x','keeplimits');

hgsave(p,'adp_figure20_Wave_H_dir_sprd.fig');

%%
p=figure(21);
subplot(3,1,1),imagesc(WaveData.mtime,WaveData.fmt,WaveData.SSE')
set(gca,'YDir','normal'); 
q=colorbar; xlabel(q,'Snn');
ylabel('f [Hz]'), datetick('x','keeplimits')
title('wave power spectra, direction and spread');

subplot(3,1,2),imagesc(WaveData.mtime,WaveData.fmt,WaveData.DIR_heading')
set(gca,'YDir','normal'); %caxis([-30 30]);
q=colorbar; xlabel(q,'Dir Going');
ylabel('f [Hz]'), datetick('x','keeplimits')

subplot(3,1,3),imagesc(WaveData.mtime,WaveData.fmt,WaveData.SPREAD')
set(gca,'YDir','normal'); 
q=colorbar; xlabel(q,'Sprd');
ylabel('f [Hz]'), datetick('x','keeplimits')
hgsave(p,'adp_figure21_SSE_DIR_SPREAD.fig')

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

hgsave(p,'adp_figure22_SSE_US_VU.fig');

%%
p=figure(23);
plot(WaveData.mtime,WaveData.Us,...
    WaveData.mtime,WaveData.Vs)
ylabel('Stokes Drift (m/s)'); datetick('x','keeplimits');
legend('Us','Vs')    
hgsave(p,'adp_figure23_Stokes_Drift.fig');

%%
p=figure(24);
subplot(2,1,1), plot(WaveData.mtime,WaveData.ztest_ss)
ylabel('ztest'), datetick('x','keeplimits')
subplot(2,1,2), plot(WaveData.mtime,WaveData.Cpu_ss)
ylabel('Cpu'), datetick('x','keeplimits')
hgsave(p,'adp_figure24_ztest_Cpu.fig');



%% lagrangian flows
p=figure(25);
plot([-1 1],[0 0],'--k',[0 0],[-1 1],'--k'); 
hold on;
scatter(Vel_East_Lagr,Vel_North_Lagr)
title('Figure 1, Depth-Averaged Lagrangian Velocity');
xlabel('Depth-Avg East (Ul) (m/s)');
ylabel('Depth-Avg North (Vl) (m/s)');
axis([-0.5 0.5 -0.5 0.5]);
name=['adp_figure25_UV.fig'];
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
name=['adp_figure26_UV.fig'];
hgsave(p,name);

%% plot avg vertical profiles
if extrap
    
p=figure(30);
subplot(1,2,1), 
plot(Vel_Maj_Prof,bin_MAB,'x');
hold all, plot(Vel_Maj_Extrap_Prof,bin_MAB_extrap);
legend('measured','extrapolated')
ylabel('MAB (m)'), xlabel('U Maj (m/s)');

subplot(1,2,2), 
plot(Vel_Min_Prof,bin_MAB,'x');
hold all, plot(Vel_Min_Extrap_Prof,bin_MAB_extrap);
% legend('measured','extrapolated')
ylabel('MAB (m)'), xlabel('U Min (m/s)');

hgsave(p,'adp_figure30_vertical_profiles.fig');

%% plot extrapolated velocities

p=figure(31);
subplot(3,1,1),
var = vel_maj_filt_extrap;
indx2 = isnan(var);
var(indx2) = 0;
imagesc(mtime,bin_MAB_extrap,var); colormap('french');shading flat;
set(gca,'ydir','normal');
hold on; plot(mtime_filt,depth_filt,'-b');
% ylim([min(bin_MAB) max(depth_filt)]);
q=colorbar; %
c = -2*nanmean(nanstd(var,0,2))-abs(nanmean(nanmean(var)));
d = 2*nanmean(nanstd(var,0,2))+abs(nanmean(nanmean(var)));
caxis([c d]) % set limits to 2 stdev of mean
datetickzoom('x');
ylabel(q,'u Maj (m/s)'); 
ylabel('MAB (m)');
title('Velocity Profiles (filtered) extrapolated');

subplot(3,1,2),
var = vel_min_filt_extrap;
indx2 = isnan(var);
var(indx2) = 0;
imagesc(mtime_filt,bin_MAB_extrap,var); colormap('french');shading flat;
set(gca,'ydir','normal');
hold on; plot(mtime_filt,depth_filt,'-b');
q=colorbar; 
c = -2*nanmean(nanstd(var,0,2))-abs(nanmean(nanmean(var)));
d = 2*nanmean(nanstd(var,0,2))+abs(nanmean(nanmean(var)));
caxis([c d]) % set limits to 2 stdev of mean
datetickzoom('x');
ylabel(q,'u Min (m/s)'); 
ylabel('MAB (m)');

subplot(3,1,3),
var = vel_vert_filt_extrap;
indx2 = isnan(var);
var(indx2) = 0;
imagesc(mtime_filt,bin_MAB_extrap,var); colormap('french');shading flat;
set(gca,'ydir','normal');
hold on; plot(mtime_filt,depth_filt,'-b');
q=colorbar; 
c = -2*nanmean(nanstd(var,0,2))-abs(nanmean(nanmean(var)));
d = 2*nanmean(nanstd(var,0,2))+abs(nanmean(nanmean(var)));
caxis([c d]) % set limits to 2 stdev of mean
datetickzoom('x');
ylabel(q,'w (m/s)'); 
ylabel('MAB (m)');

name=['adp_figure31_vel_filt_extrap.fig'];
hgsave(p,name);

end

if input.Waves
%% RMS velocitites
p = figure(40);
subplot(4,1,1)
imagesc(WaveData.mtime,bin_MAB,WaveData.up_rms),shading flat
set(gca,'ydir','normal');
cb = colorbar; xlabel(cb,'u'' rms')
axis tight
datetickzoom('x','keeplimits')
ylabel('z (m)')

subplot(4,1,2)
imagesc(WaveData.mtime,bin_MAB,WaveData.vp_rms),shading flat
set(gca,'ydir','normal');
cb = colorbar; xlabel(cb,'v'' rms')
axis tight
datetickzoom('x','keeplimits')
ylabel('z (m)')

subplot(4,1,3)
imagesc(WaveData.mtime,bin_MAB,WaveData.wp_rms),shading flat
set(gca,'ydir','normal');
cb = colorbar; xlabel(cb,'w'' rms')
axis tight
datetickzoom('x','keeplimits')
ylabel('z (m)')

subplot(4,1,4)
imagesc(WaveData.mtime,bin_MAB,WaveData.uvwp_rms),shading flat
set(gca,'ydir','normal');
cb = colorbar; xlabel(cb,'abs(u'') rms')
axis tight
datetickzoom('x','keeplimits')
ylabel('z (m)')

hgsave(p,'adp_figure40_up_RMS')
end
end