function spikeRaster = samples_cell2mat(sampleCell,T,Dt)
    if nargin == 2
        Dt = 1;
    end
    bins = 0:Dt:(T-Dt);
    nsamples = length(sampleCell);
    spikeRaster = zeros(nsamples,length(bins));
    for i = 1:nsamples
        tmp = histc([sampleCell{i}(:); inf],[bins, (T+1)]);
        %size(spikeRaster(i,:))
        %size(tmp(1:(end-1))')
        spikeRaster(i,:) = tmp(1:(end-1))';
    end