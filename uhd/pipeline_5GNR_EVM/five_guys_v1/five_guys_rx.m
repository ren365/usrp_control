clear all 
close all

settings = NR5G_setting();
[trBlk,txWaveform,txsetting] = NR5G_send(settings);

rx_fileName = "rx.bin";

[rx_waveforms,SNRs] = NR5G_split_package(rx_fileName,length(txsetting.txWaveform))

function [EVM,BER] = get_EVM_BER()



end

function [pdschEq,BER_ratio] = NR5G_receive(rxWaveform,txsetting)

    % OFDM Demodulation
    rxGrid = nrOFDMDemodulate(txsetting.carrier,rxWaveform);
    % Channel Estimation
    % Perform practical channel estimation between layers and receive antennas.
    [estChGridLayers,noiseEst] = nrChannelEstimate(txsetting.carrier,txsetting.rxGrid,txsetting.dmrsIndices,txsetting.dmrsSymbols,'CDMLengths',txsetting.pdsch.DMRS.CDMLengths);

    % Remove precoding from estChannelGrid before precoding
    % matrix calculation
    estChGridAnts = precodeChannelEstimate(estChGridLayers,conj(txsetting.precodingWeights));

    % Get precoding matrix for next slot
    newPrecodingWeight = getPrecodingMatrix(txsetting.pdsch.PRBSet,txsetting.pdsch.NumLayers,estChGridAnts);

    % Equalization
    [pdschRx,pdschHest] = nrExtractResources(txsetting.pdschIndices,txsetting.rxGrid,estChGridLayers);
    [pdschEq,csi] = nrEqualizeMMSE(pdschRx,pdschHest,noiseEst);
    constPlot.ChannelNames = "Layer "+(txsetting.pdsch.NumLayers:-1:1);
    constPlot.ShowLegend = true;
    % Constellation for the first layer has a higher SNR than that for the
    % last layer. Flip the layers so that the constellations do not mask
    % each other.
    constPlot(fliplr(pdschEq));
    % PDSCH Decoding
    [dlschLLRs,rxSymbols] = nrPDSCHDecode(txsetting.carrier,txsetting.pdsch,pdschEq,noiseEst);
    % Scale LLRs by CSI
    csi = nrLayerDemap(csi);                                    % CSI layer demapping
    for cwIdx = 1:txsetting.pdsch.NumCodewords
        Qm = length(dlschLLRs{cwIdx})/length(rxSymbols{cwIdx}); % Bits per symbol
        csi{cwIdx} = repmat(csi{cwIdx}.',Qm,1);                 % Expand by each bit per symbol
        dlschLLRs{cwIdx} = dlschLLRs{cwIdx} .* csi{cwIdx}(:);   % Scale
    end
    % DL-SCH Decoding
    decodeDLSCH.TransportBlockLength = txsetting.trBlkSizes;
    [decbits,blkerr] = decodeDLSCH(dlschLLRs,txsetting.pdsch.Modulation,txsetting.pdsch.NumLayers, ...
    txsetting.harqEntity.RedundancyVersion,txsetting.harqEntity.HARQProcessID);
    % % HARQ Process Update
    % statusReport = updateAndAdvance(txsetting.harqEntity,blkerr,txsetting.trBlkSizes,txsetting.pdschInfo.G);    
    % disp("Slot "+(nSlot)+". "+statusReport);

    [number,BER_ratio] = biterr(decbits,txsetting.trBlk)

end

% take average across 100 ? 
function [rx_waveforms,SNRs] = NR5G_split_package(rx_fileName,txlen)

    avgNum = 100;
    rx_waveforms = zeros(txlen,avgNum);
    SNRs = [];
    rxWaveform_orignal = File2Wave(rx_fileName);
    rxWaveform_orignal = rxWaveform_orignal(end-txlen*avgNum:end,:);

    rxWaveform_remain = rxWaveform_orignal;
    for indx = 1: avgNum
        % Timing Synchronization
        [t,mag] = nrTimingEstimate(carrier,rxWaveform_remain,dmrsIndices,dmrsSymbols);
        offset = hSkipWeakTimingOffset(offset,t,mag);
        
        rxWaveform = rxWaveform(1+offset:1+offset+txlen-1,:);

        rx_waveforms(:,indx) = rxWaveform;
        rxWaveform_remain = rxWaveform(1+offset+txlen:end,:);

        if offset > txlen
            rxWaveform = rxWaveform_orignal(1+offset:1+offset+txlen-1,:);
            noise_level = rxWaveform_orignal(offset-txlen+1:offset,:);
            snr_calculated = snr(rxWaveform,noise_level);
            SNRs = [SNRs,snr_calculated];
        end
    end

end


function wtx = getPrecodingMatrix(PRBSet,NLayers,hestGrid)
% Calculate precoding matrix given an allocation and a channel estimate
    
    % Allocated subcarrier indices
    allocSc = (1:12)' + 12*PRBSet(:).';
    allocSc = allocSc(:);
    
    % Average channel estimate
    [~,~,R,P] = size(hestGrid);
    estAllocGrid = hestGrid(allocSc,:,:,:);
    Hest = permute(mean(reshape(estAllocGrid,[],R,P)),[2 3 1]);
    
    % SVD decomposition
    [~,~,V] = svd(Hest);
    
    wtx = V(:,1:NLayers).';
    wtx = wtx/sqrt(NLayers); % Normalize by NLayers
end

function estChannelGrid = getInitialChannelEstimate(channel,carrier)
% Obtain an initial channel estimate for calculating the precoding matrix.
% This function assumes a perfect channel estimate

    % Clone of the channel
    chClone = channel.clone();
    chClone.release();

    % No filtering needed to get channel path gains
    chClone.ChannelFiltering = false;    
    
    % Get channel path gains
    [pathGains,sampleTimes] = chClone();
    
    % Perfect timing synchronization
    pathFilters = getPathFilters(chClone);
    offset = nrPerfectTimingEstimate(pathGains,pathFilters);
    
    % Perfect channel estimate
    estChannelGrid = nrPerfectChannelEstimate(carrier,pathGains,pathFilters,offset,sampleTimes);
end

function refPoints = getConstellationRefPoints(mod)
% Calculate the reference constellation points for a given modulation
% scheme.
    switch mod
        case "QPSK"
            nPts = 4;
        case "16QAM"
            nPts = 16;
        case "64QAM"
            nPts = 64;
        case "256QAM"
            nPts = 256;            
    end
    binaryValues = int2bit(0:nPts-1,log2(nPts));
    refPoints = nrSymbolModulate(binaryValues(:),mod);
end

function estChannelGrid = precodeChannelEstimate(estChannelGrid,W)
% Apply precoding matrix W to the last dimension of the channel estimate.

    % Linearize 4-D matrix and reshape after multiplication
    K = size(estChannelGrid,1);
    L = size(estChannelGrid,2);
    R = size(estChannelGrid,3);
    estChannelGrid = reshape(estChannelGrid,K*L*R,[]);
    estChannelGrid = estChannelGrid*W;
    estChannelGrid = reshape(estChannelGrid,K,L,R,[]);

end
