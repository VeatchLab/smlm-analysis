%% Run stacors and stxcors on some data from Jenny
clear stacors_from_imagestructs
clear stxcors_from_imagestructs
clear all

datadir = '/lipid/group/data/Jenny/2021_09_03_CH27-CD4tmdGFP/PM-mEos3.2_Atto655-streptavidin/cell2'

%%
record = load(fullfile(datadir,'record.mat'));
cdata = load(fullfile(datadir,'cordata.mat'));
is = cdata.is;

r = 12.5:25:1000; % or 25:25:1000 (I forget whether this is bin center or edge...)
frame_time = record.metadata.frame_time;
tau = frame_time:frame_time*50:frame_time*1000;
[c, ~, ~, cN, cNorm] = stxcors_from_imagestructs(is, r, tau);
[g, ~, ~, gN, gNorm] = stacors_from_imagestructs(is, r, tau);

%% prep a normalized struct of answers
newres = struct('c',squeeze(c),'cN',cN, 'cNorm', cNorm, 'g1', squeeze(g{1}), 'g2', squeeze(g{2}), 'gN', gN, 'gNorm',gNorm);
%% if for the first time, do the following
% save testresults.mat -struct newres


oldres = load('testresults.mat');

fprintf('Checking new results for equality with old ones by field\n');
fn = fieldnames(oldres);
for i = 1:numel(fn)
    ndif = sum(1e-15 < abs(1 - (newres.(fn{i})(:) ./ oldres.(fn{i})(:))));
    
    fprintf('%8s: %6d differences\n', fn{i}, ndif);
end