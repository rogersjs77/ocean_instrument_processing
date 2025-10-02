% plot sbe37 data
% Justin Rogers
% Stanford EFML 8/2013

function []=sbe19_plot()
set(0,'defaulttextinterpreter','latex');
set(0,'defaulttextfontname','times');
set(0,'defaultaxesfontname','times');
set(0,'defaultaxesfontsize',12);
set(0,'defaulttextfontsize',12);

close all; clear all;
load Processed_SBE19_data

%% raw data
p=figure(1);
for i=1:length(files)
subplot(1,3,1)    
plot(sbe19.temp{i},-sbe19.depth{i}); hold all;
% set('Ydir','reverse');

subplot(1,3,2)    
plot(sbe19.sal{i},-sbe19.depth{i}); hold all;
% set('Ydir','reverse');

subplot(1,3,3)    
plot(sbe19.oxygen{i},-sbe19.depth{i}); hold all;
% set('Ydir','reverse');
% leg(i)=cellstr(name{i});
end

% legend(cellstr(name),'location','nw')
subplot(1,3,1); xlabel('temp C'), ylabel('depth (m)');
title('SBE19 Raw CTD Data');

subplot(1,3,2); xlabel('sal'), ylabel('depth (m)');
subplot(1,3,3); xlabel('oxygen'), ylabel('depth (m)');

hgsave(p,'figure_sbe19_raw_data.fig');

%% downcast data
p=figure(2);
for i=1:length(files)
subplot(1,3,1)    
plot(temp_dn{i},-depth_dn{i}); hold all;
% set('Ydir','reverse');

subplot(1,3,2)    
plot(sal_dn{i},-depth_dn{i}); hold all;
% set('Ydir','reverse');

subplot(1,3,3)    
plot(oxygen_dn{i},-depth_dn{i}); hold all;
% set('Ydir','reverse');
% leg(i)=cellstr(name{i});
end

% legend(cellstr(name),'location','nw')
subplot(1,3,1); xlabel('temp C'), ylabel('depth (m)');
title('SBE19 Downcast CTD Data');

subplot(1,3,2); xlabel('sal'), ylabel('depth (m)');
subplot(1,3,3); xlabel('oxygen'), ylabel('depth (m)');

hgsave(p,'figure_sbe19_downcast_data.fig');


end