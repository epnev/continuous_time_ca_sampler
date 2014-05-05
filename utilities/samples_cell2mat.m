function spikeRaster = samples_cell2mat(sampleCell,T)
    nsamples = length(sampleCell);
    spikeRaster = zeros(nsamples,T);
    for i = 1:nsamples
        tmp = histc([sampleCell{i} inf],[0:(T-1) (T+1)]);
        %size(spikeRaster(i,:))
        %size(tmp(1:(end-1))')
        spikeRaster(i,:) = tmp(1:(end-1))';
    end