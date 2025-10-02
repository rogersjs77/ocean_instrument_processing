% plot sbe56 data
% Justin Rogers
% Stanford EFML 8/2013

function []=sbe56_plot()
set(0,'defaulttextinterpreter','latex');
set(0,'defaulttextfontname','times');
set(0,'defaultaxesfontname','times');
set(0,'defaultaxesfontsize',12);
set(0,'defaulttextfontsize',12);

close all; clear all;
load Processed_SBE56_data

%% raw data
if exist('sbe65')
p=figure(1);
for i=1:length(mab)
plot(sbe56.mtime{i},sbe56.temp{i}); hold all;
leg{i}=sprintf([site{i} ' %gMAB'],mab(i));
end
datetickzoom('x');
legend((leg),'location','nw')
xlabel('Date'), ylabel('Temperature ($^{\circ}C$)');
title('Raw Data');
print -djpeg -r300 figure_sbe56_raw_data
% hgsave(p,'figure_sbe56_raw_data.fig');
end
%% trimmed time data
p=figure(2);
for i=1:length(mab)
plot(mtime{i},temp{i}); hold all;
leg{i}=sprintf([site{i} ' %gMAB'],mab(i));
end
datetickzoom('x');
legend((leg),'location','nw')
xlabel('Date'), ylabel('Temperature ($^{\circ}C$)');
title('Raw data trimmed time');
print -djpeg -r300 figure_sbe56_trimmed_data
% hgsave(p,'figure_sbe56_trimmed_data.fig');

%% filtered data
p=figure(3);
plot(mtime_filt,temp_filt);
for i=1:length(mab)
leg{i}=sprintf([site{i} ' %gMAB'],mab(i));
end
datetickzoom('x');
legend((leg),'location','nw')
xlabel('Date'), ylabel('Temperature ($^{\circ}C$)');
title('Filtered Data');
hgsave(p,'figure_sbe56_filtered_data.fig');

%% filtered profile
if length(mab)>1
    if max(mab)>0
       y = mab;
       ylab = 'MAB';
    else
        y = 1:length(mab);
        ylab = 'Site No';
    end

p=figure(4);
imagesc(mtime_filt,mab,temp_filt),
set(gca,'Ydir','normal')
colorbar;
datetick('x')
xlabel('Date'), 
ylabel(ylab)
title('Filtered Data Temperature ($^{\circ}C$)')
hgsave(p,'figure_sbe56_filtered_profile.fig')
end

%% power spectrum
p=figure(5);

for i=1:length(mab)
loglog(FFT_s{i},2*(abs(FFT_temp{i})).^2), hold all
leg{i}=sprintf([site{i} ' %gMAB'],mab(i));
end
axis square
ylabel('$\Phi(T)$')
xlabel('Hz');
title('Power Spectrum of Temperature ($^{\circ}C$)');
legend((leg),'location','nw')
print -djpeg -r300 figure_sbe56_fft_temp
% hgsave(p,'figure_sbe56_fft_temp.fig');

%%

end