clear;
load traceData.mat;
addpath utilities
% Y = squeeze(traceData.traces(129,7,:));   % pick a particular trace (low SNR)
Y = mean(squeeze(traceData.traces(:,7,:))); % average over ROI (high SNR)
Yr = (Y(:) - min(Y(:)))/(max(Y(:)) - min(Y(:)));   % normalize data
T = length(Yr);
%%
SAMP = cont_ca_sampler(Yr);
plot_continuous_samples(SAMP,Yr);