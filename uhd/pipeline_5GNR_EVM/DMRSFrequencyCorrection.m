function [correctedWaveform,appliedFrequencyCorrection] = DMRSFrequencyCorrection(waveform,sampleRate,frequencyCorrectionRange,tmwavegen,resourcesInfo)
    % waveform - Waveform to be corrected. Needs to be a Nx1 column vector.
    % sampleRate - Sample rate of waveform
    % frequencyCorrectioRange - Range and granularity at which frequency
    % correction is inspected
    % tmwavegen and resourcesInfo - Outputs of generateWaveform method 
    [pdschArray,~,carrier] = hListTargetPDSCHs(tmwavegen.Config,resourcesInfo.WaveformResources);
    bwpCfg = tmwavegen.Config.BandwidthParts{1,1};
    nSlots = carrier.SlotsPerFrame;
    
    % Generate a reference grid spanning 10 ms (one frame). This grid
    % contains only the DM-RS and is used for synchronization.
    refGrid = referenceGrid(carrier,bwpCfg,pdschArray,nSlots);
    
    % Apply frequency offsets to the waveform as specified by
    % freuqnecyCorrectionRange.
    nSamples = (0:length(waveform)-1)';
    frequencyShift = (2*pi*frequencyCorrectionRange.*nSamples)./sampleRate;
    
    % Each column represents an offset waveform.
    offsetWaveforms = waveform.*exp(1j*frequencyShift);
    
    [~,mag] = nrTimingEstimate(offsetWaveforms,carrier.NSizeGrid,...
        carrier.SubcarrierSpacing,nSlots,refGrid, ...
        "SampleRate",sampleRate);
    
    % Find the frequency at which the DM-RS correlation is at a maximum.
    [~,index] = max(max(mag));    
    appliedFrequencyCorrection = frequencyCorrectionRange(index);
    correctedWaveform = offsetWaveforms(:,index);
end

function refGrid = referenceGrid(carrier,bwpCfg,pdschArray,nSlots)
    % Create a reference grid for the required number of slots. The grid
    % contains the DM-RS symbols specified in pdschArray. The function
    % returns REFGRID of dimensions K-by-S-by-L, where K is the number of
    % subcarriers of size carrier.NSizeGrid, S is the number of symbols
    % spanning nSlots, and L is the number of layers.

    nSubcarriers = carrier.NSizeGrid * 12;
    L = carrier.SymbolsPerSlot*nSlots;                           % Number of OFDM symbols in the reference grid
    nLayers = size(pdschArray(1).Resources(1).ChannelIndices,2);
    bwpStart = bwpCfg.NStartBWP;
    bwpLen = bwpCfg.NSizeBWP;
    refGrid = zeros(nSubcarriers,L,nLayers);                     % Empty grid
    bwpGrid = zeros(bwpLen*12,L,nLayers);
    rbsPerSlot = bwpLen*12*carrier.SymbolsPerSlot;

    % Populate the DM-RS symbols in the reference grid for all slots. Place
    % bwpGrid in a carrier grid (at an appropriate location) in case the
    % BWP size is not the same as the carrier grid
    for slotIdx = carrier.NSlot + (0:nSlots-1)
        [~,~,dmrsIndices,dmrsSymbols] = hSlotResources(pdschArray,slotIdx);
        if ~isempty(dmrsIndices)
            for layerIdx = 1:nLayers
                if layerIdx <= size(dmrsIndices,2)
                    dmrsIndices(:,layerIdx) = dmrsIndices(:,layerIdx) - rbsPerSlot*(layerIdx -1) + (L*bwpLen*12*(layerIdx-1));
                    bwpGrid(dmrsIndices(:,layerIdx)+(slotIdx-carrier.NSlot)*rbsPerSlot) = dmrsSymbols(:,layerIdx);
                end
            end
            refGrid(12*bwpStart+1:12*(bwpStart+bwpLen),:,:) = bwpGrid;
        end
    end
end