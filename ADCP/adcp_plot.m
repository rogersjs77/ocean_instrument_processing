% function [] = plot_ADCP()
%
% this plots standard ADCP data from a_process_ADCP
% 
% Justin Rogers, stanford efml
% June 2017

function [] = plot_ADCP()
set(0,'defaulttextinterpreter','latex');
set(0,'defaulttextfontname','times');
set(0,'defaultaxesfontname','times');
set(0,'defaultaxesfontsize',12);
set(0,'defaulttextfontsize',12);

% run presentation_plot
close all; clear all;

files=dir('*PostProcess.mat');  

load(files(1).name);

%% scatter plots
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

print -djpeg -r300 adcp_figure1_UV
% close
% hgsave(p,name);

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
ylabel(q,'<u maj> (m/s)'); 
ylabel('MAB (m)');
title('Velocity Profiles (filtered)');

subplot(4,1,2),
var = vel_min_filt;
pcolorjw(mtime_filt,bin_MAB,var); colormap('french');shading flat;
set(gca,'ydir','normal');
hold on; plot(mtime_filt,depth_filt,'-k');
ylim([min(bin_MAB) max(depth_filt)]);
q=colorbar; 
c = -2*nanmean(nanstd(var,0,2))-abs(nanmean(nanmean(var)));
d = 2*nanmean(nanstd(var,0,2))+abs(nanmean(nanmean(var)));
caxis([c d]) % set limits to 2 stdev of mean
datetickzoom('x');
ylabel(q,'<u min> (m/s)'); 
ylabel('MAB (m)');

subplot(4,1,3),
var = vel_vert_filt;
pcolorjw(mtime_filt,bin_MAB,var); colormap('french');shading flat;
set(gca,'ydir','normal');
hold on; plot(mtime_filt,depth_filt,'-k');
ylim([min(bin_MAB) max(depth_filt)]);
q=colorbar; 
c = -2*nanmean(nanstd(var,0,2))-abs(nanmean(nanmean(var)));
d = 2*nanmean(nanstd(var,0,2))+abs(nanmean(nanmean(var)));
caxis([c d]) % set limits to 2 stdev of mean
datetickzoom('x');
ylabel(q,'<w> (m/s)'); 

subplot(4,1,4),
if exist('vel_err_filt')
var = vel_err_filt;
pcolorjw(mtime_filt,bin_MAB,var); colormap('french');shading flat;
set(gca,'ydir','normal');
q=colorbar; 
c = -2*nanmean(nanstd(var,0,2))-abs(nanmean(nanmean(var)));
d = 2*nanmean(nanstd(var,0,2))+abs(nanmean(nanmean(var)));
caxis([c d]) % set limits to 2 stdev of mean
datetickzoom('x');
ylabel(q,'<err> (m/s)'); 
end
xlabel('Date'); ylabel('MAB (m)');

print -djpeg -r300 adcp_figure3_umvmwm
% hgsave(p,name);
% close
%% depth averaged earth coordinates
p=figure(4);
subplot(2,1,1); plot(mtime_filt,Vel_Maj)
if mtime(end)-mtime(1) > 2
hold all, plot(mtime_filt,Vel_Maj_notides)
end
datetickzoom('x')
ylabel('U-major x filt (m/s)')
subplot(2,1,2); plot(mtime_filt,Vel_Min)
if mtime(end)-mtime(1) > 2
hold all, plot(mtime_filt,Vel_Min_notides)
end
datetickzoom('x')
ylabel('V-minor y filt (m/s)')
xlabel('date')
legend('lp','detide')

print -djpeg -r300 adcp_figure4_UV
close

%% depth averaged rotated coordinates
p=figure(5);
subplot(2,1,1); plot(mtime_filt,Vel_Maj)
datetickzoom('x')
ylabel('U-major filt (m/s)')
subplot(2,1,2); plot(mtime_filt,Vel_Min)
datetickzoom('x')
ylabel('V-minor filt (m/s)')
xlabel('date')
print -djpeg -r300 adcp_figure5_UmVm
% hgsave(p,name);
% close

%% plot intensity
var = adcp.intens;
ttl={'adcp.intens beam1','beam2','beam3','beam4'};

p=figure(6);
ii = 1:size(var,3);
for j=1:size(var,2)
    temp=0*adcp.east_vel;    
    for i=1:length(bin_MAB_orig)
        temp(i,ii)=nanmax(var(i,j,ii),2);
    end
    subplot(2,2,j); 
    imagesc(adcp.mtime,bin_MAB_orig,temp),set(gca,'YDir','normal'),
    hold all, plot([mtime(1) mtime(end)],bin_MAB(1)*[1 1],'k')
    hold all, plot([mtime(1) mtime(end)],bin_MAB(end)*[1 1],'k')
    qq=colorbar; 
    caxis([0 150])
    hold all, plot(mtime_filt,depth_filt,'k')
    datetick('x','keeplimits')
    ylabel('z (m)'); %xlabel('time');
    % datetick('x','keeplimits');
    xlabel(qq,[ttl{j}]);
end
 print -djpeg -r300 adcp_figure6_intens

close

%% plot correlation
var = adcp.corr;
ttl={'adcp.corr beam1','beam2','beam3','beam4'};

p=figure(7);
ii = 1:size(var,3);
for j=1:size(var,2)
    temp=0*adcp.east_vel;    
    for i=1:length(bin_MAB_orig)
        temp(i,ii)=nanmax(var(i,j,ii),2);
    end
    subplot(2,2,j); 
    imagesc(adcp.mtime,bin_MAB_orig,temp),set(gca,'YDir','normal'),
    hold all, plot([mtime(1) mtime(end)],bin_MAB(1)*[1 1],'k')
    hold all, plot([mtime(1) mtime(end)],bin_MAB(end)*[1 1],'k')
    qq=colorbar; 
    caxis([0 150])
    hold all, plot(mtime_filt,depth_filt,'k')
    datetick('x','keeplimits')
    ylabel('z (m)'); %xlabel('time');
    % datetick('x','keeplimits');
    xlabel(qq,[ttl{j}]);
end
 print -djpeg -r300 adcp_figure7_corr
close


%% plot perc_good
% If data is collected in beam coordinates, then the four percent 
% good values represent the percentage of the pings collected by each 
% beam for that depth cell whose correlation exceeded a low correlation 
% threshold. In the other coordinate frames (ADCP, Ship and Earth Coordinates),
% the four Percent Good values represent (in order): 
% 1) The percentage of good three beam solutions (one beam rejected); 
% 2) The percentage of good transformations (error velocity threshold not exceeded); 
% 3) The percentage of measurements where more than one beam was bad; and 
% 4) The percentage of measurements with four beam solutions.

var = adcp.perc_good;
ttl={'perc good 3bm sol','perc good error not exceeded','perc meas w 1 bm bad','perc meas w 4bm sol'};

p=figure(8);
ii = 1:size(var,3);
for j=1:size(var,2)
    temp=0*adcp.east_vel;    
    for i=1:length(bin_MAB_orig)
        temp(i,ii)=nanmax(var(i,j,ii),2);
    end
    subplot(2,2,j); 
    imagesc(adcp.mtime,bin_MAB_orig,temp),set(gca,'YDir','normal'),
    hold all, plot([mtime(1) mtime(end)],bin_MAB(1)*[1 1],'k')
    hold all, plot([mtime(1) mtime(end)],bin_MAB(end)*[1 1],'k')
    qq=colorbar; 
    caxis([0 150])
    hold all, plot(mtime_filt,depth_filt,'k')
    datetick('x','keeplimits')
    ylabel('z (m)'); %xlabel('time');
    % datetick('x','keeplimits');
    xlabel(qq,[ttl{j}]);
end
 print -djpeg -r300 adcp_figure8_perc_good
close
%% plot velocities
ttl={'east vel (m/s)','north vel (m/s)','vert vel (m/s)','error vel (m/s)'};

p=figure(9);
for j=1:4
subplot(2,2,j); 
if j==1
    var=adcp.east_vel;
    imagesc(adcp.mtime,bin_MAB_orig,var),set(gca,'YDir','normal'),
    hold all, plot([mtime(1) mtime(end)],bin_MAB(1)*[1 1],'k')
    hold all, plot([mtime(1) mtime(end)],bin_MAB(end)*[1 1],'k')
    qq=colorbar;   c = 2*nanmean(nanstd(var,0,2));
    caxis([-c c]) % set limits to 2 stdev
    xlabel(qq,[ttl{j}]);
elseif j==2
    var=adcp.north_vel;
    imagesc(adcp.mtime,bin_MAB_orig,var),set(gca,'YDir','normal'),
    hold all, plot([mtime(1) mtime(end)],bin_MAB(1)*[1 1],'k')
    hold all, plot([mtime(1) mtime(end)],bin_MAB(end)*[1 1],'k')
    qq=colorbar;   c = 2*nanmean(nanstd(var,0,2));
    caxis([-c c]) % set limits to 2 stdev
    xlabel(qq,[ttl{j}]);
elseif j==3
    var=adcp.vert_vel;
    imagesc(adcp.mtime,bin_MAB_orig,var),set(gca,'YDir','normal'),
    hold all, plot([mtime(1) mtime(end)],bin_MAB(1)*[1 1],'k')
    hold all, plot([mtime(1) mtime(end)],bin_MAB(end)*[1 1],'k')
    qq=colorbar;   c = 2*nanmean(nanstd(var,0,2));
    caxis([-c c]) % set limits to 2 stdev
    xlabel(qq,[ttl{j}]);
elseif j==4 && isfield(adcp,'adcp.error_vel')
    if ~isempty(input.bad_beam)
    var=adcp.error_vel_orig;
    else
    var=adcp.error_vel;
    end
    imagesc(adcp.mtime,bin_MAB_orig,var),set(gca,'YDir','normal'),
    hold all, plot([mtime(1) mtime(end)],bin_MAB(1)*[1 1],'k')
    hold all, plot([mtime(1) mtime(end)],bin_MAB(end)*[1 1],'k')
    qq = colorbar;   c = 2*nanmean(nanstd(var,0,2));
    caxis([-c c]) % set limits to 2 stdev
    xlabel(qq,[ttl{j}]);
end
hold all, plot(mtime_filt,depth_filt,'k')
datetick('x','keeplimits')
ylabel('z (m)'); %xlabel('time');
% title([ttl{j}]);
end

print -djpeg -r300 adcp_figure9_uvwerr

close
%% temperature pressure depth

p=figure(10);
subplot(3,1,1);
if exist('temperature')
var = temperature;
plot(mtime,var,mtime_filt,temperature_filt);
datetickzoom('x','keeplimits');
legend('raw','filtered');
ylabel('T (C)');
title('temperature pressure and depth');
end
subplot(3,1,2);
var = depth;
plot(mtime,var,mtime_filt,depth_filt);
datetickzoom('x','keeplimits');
ylabel('depth (m)');

subplot(3,1,3);
plot(mtime_filt,eta_filt)
if mtime(end)-mtime(1) > 2
hold all, plot(mtime_filt,eta_notides)
end
datetickzoom('x','keeplimits');
xlabel('date'); ylabel('$\eta$ (m)');
legend('lp','detide')

print -djpeg -r300 adcp_figure10_temp_press_depth
close


%% pressures

p=figure(11);
subplot(2,1,1);
if ~isempty(P_atm)
var = adcp.pressure*10;
plot(adcp.mtime,var,adcp.mtime,0*adcp.mtime+P_atm);
datetickzoom('x','keeplimits');
legend('P raw','P atm');
ylabel('P (Pa)');
title('presure adjustment');
end

subplot(2,1,2);
var = adcp.press_Pa;
plot(adcp.mtime,var);
datetickzoom('x','keeplimits');
ylabel('P (Pa)');
legend('P adjusted')

print -djpeg -r300 adcp_figure11_pressure_adjustment
close

%% heading pitch roll
if exist('heading')
    
p=figure(12);
subplot(3,1,1);
var = heading;
plot(mtime,var);
c = [-3*nanstd(var) 3*nanstd(var)]+nanmean(var);
ylim(c);
datetickzoom('x','keeplimits');
ylabel('heading (deg)');
title('heading pitch and roll');

subplot(3,1,2);
var = pitch;
plot(mtime,var);
c = [-3*nanstd(var) 3*nanstd(var)]+nanmean(var);
ylim(c);
datetickzoom('x','keeplimits');
ylabel('pitch (deg)');

subplot(3,1,3);
var = roll;
plot(mtime,var);
c = [-3*nanstd(var) 3*nanstd(var)]+nanmean(var);
ylim(c);
datetickzoom('x','keeplimits');
xlabel('date'); ylabel('roll (deg)');

print -djpeg -r300 adcp_figure12_head_ptich_roll
close
end

%% mean errors
if exist('Vel_Error')
    
p=figure(13);
plot(Vel_Error,bin_MAB_orig,...
    [min(Vel_Error) max(Vel_Error)],Rmax*[1 1],...
    [min(Vel_Error) max(Vel_Error)],bin_MAB(1)*[1 1],...
    [min(Vel_Error) max(Vel_Error)], mean(depth)*[1 1],'-k',...
    [min(Vel_Error) max(Vel_Error)], max(depth)*[1 1],'--k',...
    [min(Vel_Error) max(Vel_Error)], min(depth)*[1 1],'--k');
xlabel('RMS Vel Error (m/s');
ylabel('z (m)');
legend('RMS Uerr','Rmax','Bot trim','Depth','location','nw');
title('RMS error velocity with depth');
print -djpeg -r300 adcp_figure13_rms_error
close
end
%% average profiles

p=figure(14);
subplot(1,2,1), plot(Vel_Maj_Prof,bin_MAB);
xlabel('U Maj (m/s)'), ylabel('MAB');
legend('pos','neg');

subplot(1,2,2), plot(Vel_Min_Prof,bin_MAB);
xlabel('U Min (m/s)'), ylabel('MAB');

print -djpeg -r300 adcp_figure14_profile
close
%% plot frequency space
p=figure(15);
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

print -djpeg -r300 adcp_figure15_fftUVd
close

%% plot power spectrum
p=figure(16);
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

print -djpeg -r300 adcp_figure16_fftUVd
close

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
ylabel('dir going'); datetick('x','keeplimits');

subplot(4,1,4),plot(WaveData.mtime,WaveData.spread_ss),
ylabel('spread (deg)'); datetick('x','keeplimits');

print -djpeg -r300 adcp_figure20_Wave_H_dir_sprd
close

%%
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

print -djpeg -r300 adcp_figure21_SSE_DIR_SPREAD
close

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

print -djpeg -r300 adcp_figure22_SSE_US_VU
close

%%
p=figure(23);
plot(WaveData.mtime,WaveData.Us,...
    WaveData.mtime,WaveData.Vs)
ylabel('Stokes Drift (m/s)'); datetick('x','keeplimits');
legend('Us','Vs')    

print -djpeg -r300 adcp_figure23_Stokes_Drift
close

%%
p=figure(24);
subplot(2,1,1), plot(WaveData.mtime,WaveData.ztest_ss)
ylabel('ztest'), datetick('x','keeplimits')
subplot(2,1,2), plot(WaveData.mtime,WaveData.Cpu_ss)
ylabel('Cpu'), datetick('x','keeplimits')
hgsave(p,'adcp_figure24_ztest_Cpu.fig');

%% lagrangian flows
p=figure(25);
plot([-1 1],[0 0],'--k',[0 0],[-1 1],'--k'); 
hold on;
scatter(Vel_East_Lagr,Vel_North_Lagr)
title('Depth-Averaged Lagrangian Velocity');
xlabel('Depth-Avg East (Ul) (m/s)');
ylabel('Depth-Avg North (Vl) (m/s)');
axis([-0.5 0.5 -0.5 0.5]);

print -djpeg -r300 adcp_figure25_UV
close


%% lagrangian vel depth averaged earth coordinates
p=figure(26);
subplot(2,1,1); plot(mtime_filt,Vel_East_Lagr)
datetickzoom('x')
ylabel('Ul-east filt (m/s)')
subplot(2,1,2); plot(mtime_filt,Vel_North_Lagr)
datetickzoom('x')
ylabel('Vl-north filt (m/s)')
xlabel('date')

print -djpeg -r300 adcp_figure26_UV
close


end


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

hgsave(p,'adcp_figure30_vertical_profiles.fig');
%%

%plot extrapolated velocities

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

print -djpeg -r300 adcp_figure31_vel_filt_extrap
close

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

print -djpeg -r300 adcp_figure40_up_RMS
close
end
end