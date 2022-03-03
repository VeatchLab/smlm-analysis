%% Resolution estimation example for a nuclear pore complex dataset

%load data
load('sampledata_npc.mat'); % The data is stored in a struct called 'data'
data

%% Optionally draw a new spatial window / ROI. Skip this step to use the one from the paper.
% spacewin_gui is a helper gui for drawing spatial windows that may have
% holes or multiple disjoint segments. Press 'save to base workspace' when
% you are done. That will save the spatial window to the matlab base workspace 
% as a struct with the appropriate format.
% spacewin_gui(data, 'PixelSize', 10) % use 10nm pixels

% then save into the 'data' struct manually
% data.spacewin = spacewin;

%% Run the resolution estimation routine
% Here we supply the data as a struct with fields x,y,t,spacewin,timewin.
% The function also accepts these arguments separately (in that order).
[corrdata, params] = spacetime_resolution(data, 'NTauBin', 10, 'Bootstrap', true);

%% Plot the correlation functions for each tau bin
figure;
plot(corrdata.r, corrdata.cWA)
tau = corrdata.taubincenters;
lh = legend(arrayfun(@num2str, tau, 'UniformOutput', false));
title(lh,'\tau (s)');
xlabel 'r (nm)'
ylabel 'g(r, \tau)'
set(gca, 'YScale', 'log');

%% Plot the normalized correlation function differences
figure;
plot(corrdata.r, corrdata.nDg);
lh = legend(arrayfun(@num2str, tau(1:end-1), 'UniformOutput', false));
title(lh,'\tau (s)');
xlabel 'r (nm)'
ylabel 'g(r, \tau)'

%% Plot the estimated resolution as a function of tau, and report average resolution
figure;
errorbar(tau(1:end-1),corrdata.s,corrdata.confint);
title(sprintf('Average resolution is %.1f nm ', corrdata.S))
set(gca, 'XScale', 'log');
xlabel '\tau (s)'
ylabel 'resolution estimate (nm)'