% plot sbe37 data
% Justin Rogers
% Stanford EFML 8/2013

function []=sbe37_plot()
set(0,'defaulttextinterpreter','latex');
set(0,'defaulttextfontname','times');
set(0,'defaultaxesfontname','times');
set(0,'defaultaxesfontsize',12);
set(0,'defaulttextfontsize',12);

close all; clear all;
load Processed_SBE37_data

%% raw data
p=figure(1);
ha = tight_subplot(4,1,0.01,[0.1 0.01],[0.1 0.01]);
for i=1:length(mab)
axes(ha(1))   
plot(sbe37.mtime{i},sbe37.temp{i}); hold all;

axes(ha(2))    
plot(sbe37.mtime{i},sbe37.sal{i}); hold all;

axes(ha(3))    
plot(sbe37.mtime{i},sbe37.press{i}); hold all;

if isfield(sbe37,'O2_mgl')
axes(ha(4)) 
plot(sbe37.mtime{i},sbe37.O2_mgl{i}); hold all;
end
leg{i}=sprintf('%gMAB',mab(i));
end

axes(ha(1));datetickzoom('x');ylabel('Temp ($^{\circ}C$)');
legend((leg),'location','nw')
title('SBE37 Raw CTD Data');
set(gca,'xticklabel',[]);

axes(ha(2)); datetickzoom('x'); ylabel('Sal');
set(gca,'xticklabel',[]);
axes(ha(3)); datetickzoom('x');  ylabel('Press');
set(gca,'xticklabel',[]);
axes(ha(4)); datetickzoom('x');  ylabel('O_2 (mg/l)');
xlabel('Date'),

hgsave(p,'figure_sbe37_raw_data.fig');

%% trimmed time data
p=figure(2);
ha = tight_subplot(4,1,0.01,[0.1 0.01],[0.1 0.01]);
for i=1:length(mab)
axes(ha(1))  
plot(mtime{i},temp{i}); hold all;

axes(ha(2))     
plot(mtime{i},sal{i}); hold all;

axes(ha(3))     
plot(mtime{i},depth{i}); hold all;
if exist('O2_mgl')
axes(ha(4))  
plot(mtime{i},O2_mgl{i}); hold all;
end

leg{i}=sprintf('%gMAB',mab(i));
end

axes(ha(1));datetickzoom('x');ylabel('Temp ($^{\circ}C$)');
legend((leg),'location','nw')
title('SBE37 Trimmed CTD Data');
set(gca,'xticklabel',[]);
axes(ha(2)); datetickzoom('x'); ylabel('Sal');
set(gca,'xticklabel',[]);
axes(ha(3)); datetickzoom('x');  ylabel('Depth');
set(gca,'xticklabel',[]);
axes(ha(4)); datetickzoom('x');  ylabel('O_2 (mg/l)');
xlabel('Date'),

hgsave(p,'figure_sbe37_trimmed_data.fig');

%% filtered data
p=figure(3);
ha = tight_subplot(4,1,0.01,[0.1 0.01],[0.1 0.01]);
axes(ha(1))
plot(mtime_filt,temp_filt);
for i=1:length(mab)
leg{i}=sprintf('%gMAB',mab(i));
end
datetickzoom('x');
legend((leg),'location','nw')
ylabel('Temp ($^{\circ}C$)');
title('Filtered Data');
set(gca,'xticklabel',[]);

axes(ha(2))
plot(mtime_filt,sal_filt);
datetickzoom('x'); ylabel('Sal');
set(gca,'xticklabel',[]);

axes(ha(3))
plot(mtime_filt,depth_filt);
datetickzoom('x'); ylabel('Depth');
set(gca,'xticklabel',[]);

axes(ha(4))
plot(mtime_filt,eta_filt,mtime_filt,eta_notides)
datetickzoom('x')
ylabel('$\eta$ (m)')
legend('lp','detide')
xlabel('Date')

hgsave(p,'figure_sbe37_filtered_data.fig');


%% power spectrum
p=figure(5);

subplot(2,2,1);
for i=1:length(mab)
loglog(FFT_s{i},2*(abs(FFT_temp{i})).^2), hold all
leg{i}=sprintf('%gMAB',mab(i));
end
axis square
ylabel('$\Phi(T)$')
xlabel('Hz');
title('Power Spectrum of T, Cond, Sal, Press');
legend((leg),'location','nw')


subplot(2,2,2);
if exist('FFT_O2')
for i=1:length(mab)
loglog(FFT_s{i},2*(abs(FFT_O2{i})).^2), hold all
end
end
axis square
ylabel('$\Phi(O_2)$')
xlabel('Hz');

subplot(2,2,3);
for i=1:length(mab)
loglog(FFT_s{i},2*(abs(FFT_sal{i})).^2), hold all
end
axis square
ylabel('$\Phi(Sal)$')
xlabel('Hz');

subplot(2,2,4);
for i=1:length(mab)
loglog(FFT_s{i},2*(abs(FFT_depth{i})).^2), hold all
end
axis square
ylabel('$\Phi(Press)$')
xlabel('Hz');

hgsave(p,'figure_sbe37_fft_temp.fig');

%% TS diagram

p=figure(6);
for j=1:length(mab)
C = mtime{j}-min(mtime{j});    
scatter(temp{j},sal{j},5,C), hold all
end
xlabel('Temp ($^{\circ}C)$'),ylabel('Sal');
title('SBE37 TS Diagram, colors are time (days blue to red)')
colorbar
hgsave(p,'figure_sbe37_TS_plot.fig');
%%
if exist('O2_mgl')
p=figure(7);
for j=1:length(mab)
C = O2_mgl{j};    
scatter(temp{j},sal{j},5,C), hold all
end
xlabel('Temp ($^{\circ}C)$'),ylabel('Sal');
title('SBE37 TS Diagram, colors are O_2 (mg/l)')
colorbar
hgsave(p,'figure_sbe37_TS_O2_plot.fig');
end
end