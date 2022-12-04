clear all 
close all

settings = NR5G_setting();
txsetting = NR5G_send(settings);

rx_fileName = "rx.bin";

[rx_waveforms,noise_levels] = NR5G_split_package(rx_fileName,txsetting);

rx_waveform_snr=[];BER_ratios = [];
noise_level_snr=[];pdschEqs=[];pdschSymbols=[];
for i=1:length(rx_waveforms(1,:))
    % snr
    rx_waveform_snr = [rx_waveform_snr;rx_waveforms(:,i)];
    noise_level_snr = [noise_level_snr;noise_levels(:,i)];
    % decode one packet
    [pdschEq,BER_ratio] = NR5G_receive(rx_waveforms(:,i),txsetting);
    % evm
    pdschSymbols = [pdschSymbols;txsetting.pdschSymbols];
    pdschEqs = [pdschEqs;pdschEq];
    % BER
    BER_ratios = [BER_ratios,BER_ratio];
end
SNR = snr(rx_waveform_snr,noise_level_snr);
[EVM,BER] = get_EVM_BER(pdschSymbols,BER_ratios,pdschEqs);
disp("SNR:"+(SNR)+", "+"EVM:"+(EVM)+", "+"BER:"+(BER))


function [EVM,BER] = get_EVM_BER(pdschSymbols,BER_ratio,pdschEq)
    evm = comm.EVM(ReferenceSignalSource="Estimated from reference constellation", ...
    ReferenceConstellation=pdschSymbols);
    EVM = evm(pdschEq);
        
    BER = mean(BER_ratio);
end

function [pdschEq,BER_ratio] = NR5G_receive(rxWaveform,txsetting)

    % OFDM Demodulation
    rxGrid = nrOFDMDemodulate(txsetting.carrier,rxWaveform);
    % Channel Estimation
    % Perform practical channel estimation between layers and receive antennas.
    [estChGridLayers,noiseEst] = nrChannelEstimate(txsetting.carrier,rxGrid,txsetting.dmrsIndices,txsetting.dmrsSymbols,'CDMLengths',txsetting.pdsch.DMRS.CDMLengths);

    % Remove precoding from estChannelGrid before precoding
    % matrix calculation
    estChGridAnts = precodeChannelEstimate(estChGridLayers,conj(txsetting.precodingWeights));

    % Get precoding matrix for next slot
    newPrecodingWeight = getPrecodingMatrix(txsetting.pdsch.PRBSet,txsetting.pdsch.NumLayers,estChGridAnts);

    % Equalization
    [pdschRx,pdschHest] = nrExtractResources(txsetting.pdschIndices,rxGrid,estChGridLayers);
    [pdschEq,csi] = nrEqualizeMMSE(pdschRx,pdschHest,noiseEst);
    
    evm = comm.EVM(ReferenceSignalSource="Estimated from reference constellation", ...
    ReferenceConstellation=txsetting.pdschSymbols);
    EVM = evm(pdschEq)

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
    txsetting.decodeDLSCH.TransportBlockLength = txsetting.trBlkSizes;
    [decbits,blkerr] = txsetting.decodeDLSCH(dlschLLRs,txsetting.pdsch.Modulation,txsetting.pdsch.NumLayers, ...
    txsetting.harqEntity.RedundancyVersion,txsetting.harqEntity.HARQProcessID);
    % % HARQ Process Update
    % statusReport = updateAndAdvance(txsetting.harqEntity,blkerr,txsetting.trBlkSizes,txsetting.pdschInfo.G);    
    % disp("Slot "+(nSlot)+". "+statusReport);

%     [number,BER_ratio] = biterr(decbits,txsetting.trBlk);
    BER_ratio = biterr(decbits,txsetting.trBlk);
    BER_ratio

end

% take average across 100 ? 
function [rx_waveforms,noise_levels] = NR5G_split_package(rx_fileName,txsetting)

    avgNum = 7;
    offset = 0;
    txlen = length(txsetting.txWaveform);
    offset_begin = txlen*20;
    rx_waveforms = zeros(txlen,avgNum);
    noise_levels = zeros(txlen,avgNum);
    SNRs = [];
    rxWaveform_orignal = File2Wave(rx_fileName);
    rxWaveform_orignal = rxWaveform_orignal(end-3*txlen*100-offset_begin:end,:);
%     plot(abs(rxWaveform_orignal))
    rxWaveform_remain = rxWaveform_orignal;
    for indx = 1: avgNum
        % Timing Synchronization
        [t,mag] = nrTimingEstimate(txsetting.carrier,rxWaveform_remain,txsetting.dmrsIndices,txsetting.dmrsSymbols);
        offset = hSkipWeakTimingOffset(offset,t,mag);
indx
        rxWaveform = rxWaveform_remain(1+offset:1+offset+txlen-1,:);

        rx_waveforms(:,indx) = rxWaveform;
        rxWaveform_remain = rxWaveform_remain(1+offset+txlen+1:end,:);

        if offset > txlen
            noise_level = rxWaveform_orignal(offset-txlen+1:offset,:);
            noise_levels(:,indx) = noise_level;
            % snr_calculated = snr(rxWaveform,noise_level);
            % SNRs = [SNRs,snr_calculated];
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
