% Simulation Parameters
SNRdB = 10;                % SNR in dB
totalNoSlots = 1;%20;         % Number of slots to simulate
perfectEstimation = false; % Perfect synchronization and channel estimation
rng("default");            % Set default random number generator for repeatability
carrier = nrCarrierConfig;
% PDSCH and DM-RS Configuration
pdsch = nrPDSCHConfig;
pdsch.Modulation = "16QAM";
pdsch.NumLayers = 1;
pdsch.PRBSet = 0:carrier.NSizeGrid-1;     % Full band allocation
pdsch.DMRS.DMRSAdditionalPosition = 1;
pdsch.DMRS.DMRSConfigurationType = 1;
pdsch.DMRS.DMRSLength = 2;
% DL-SCH Configuration
NHARQProcesses = 16;     % Number of parallel HARQ processes
rvSeq = [0 2 3 1];
% Coding rate
if pdsch.NumCodewords == 1
    codeRate = 490/1024;
else
    codeRate = [490 490]./1024;
end
% Create DL-SCH encoder object
encodeDLSCH = nrDLSCH;
encodeDLSCH.MultipleHARQProcesses = true;
encodeDLSCH.TargetCodeRate = codeRate;

% Create DLSCH decoder object
decodeDLSCH = nrDLSCHDecoder;
decodeDLSCH.MultipleHARQProcesses = true;
decodeDLSCH.TargetCodeRate = codeRate;
decodeDLSCH.LDPCDecodingAlgorithm = "Normalized min-sum";
decodeDLSCH.MaximumLDPCIterationCount = 6;
% HARQ Management
harqEntity = HARQEntity(0:NHARQProcesses-1,rvSeq,pdsch.NumCodewords);

% Channel Configuration
nTxAnts = 1;
nRxAnts = 1;
% Check that the number of layers is valid for the number of antennas
if pdsch.NumLayers > min(nTxAnts,nRxAnts)
    error("The number of layers ("+string(pdsch.NumLayers)+") must be smaller than min(nTxAnts,nRxAnts) ("+string(min(nTxAnts,nRxAnts))+")")
end

channel = nrTDLChannel;
channel.DelayProfile = "TDL-C";
channel.NumTransmitAntennas = nTxAnts;
channel.NumReceiveAntennas = nRxAnts;
ofdmInfo = nrOFDMInfo(carrier);
channel.SampleRate = ofdmInfo.SampleRate;

% Transmission and Reception
constPlot = comm.ConstellationDiagram;                                          % Constellation diagram object
constPlot.ReferenceConstellation = getConstellationRefPoints(pdsch.Modulation); % Reference constellation values
constPlot.EnableMeasurements = 1;                                               % Enable EVM measurements

% Initial timing offset
offset = 0;

estChannelGrid = getInitialChannelEstimate(channel,carrier);
newPrecodingWeight = getPrecodingMatrix(pdsch.PRBSet,pdsch.NumLayers,estChannelGrid);

% for nSlot = 0:totalNoSlots-1
    % New slot
    carrier.NSlot = nSlot;
    % Generate PDSCH indices info, which is needed to calculate the transport
    % block size
    [pdschIndices,pdschInfo] = nrPDSCHIndices(carrier,pdsch);

    % Calculate transport block sizes
    Xoh_PDSCH = 0;
    trBlkSizes = nrTBS(pdsch.Modulation,pdsch.NumLayers,numel(pdsch.PRBSet),pdschInfo.NREPerPRB,codeRate,Xoh_PDSCH);

    % Get new transport blocks and flush decoder soft buffer, as required
    for cwIdx = 1:pdsch.NumCodewords
        if harqEntity.NewData(cwIdx)
            % Create and store a new transport block for transmission
            trBlk = randi([0 1],trBlkSizes(cwIdx),1);
            setTransportBlock(encodeDLSCH,trBlk,cwIdx-1,harqEntity.HARQProcessID);

            % If the previous RV sequence ends without successful
            % decoding, flush the soft buffer
            if harqEntity.SequenceTimeout(cwIdx)
                resetSoftBuffer(decodeDLSCH,cwIdx-1,harqEntity.HARQProcessID);
            end
        end
    end

    codedTrBlock = encodeDLSCH(pdsch.Modulation,pdsch.NumLayers,pdschInfo.G,harqEntity.RedundancyVersion,harqEntity.HARQProcessID);
    pdschSymbols = nrPDSCH(carrier,pdsch,codedTrBlock);
    precodingWeights = newPrecodingWeight;
    pdschSymbolsPrecoded = pdschSymbols*precodingWeights;
    dmrsSymbols = nrPDSCHDMRS(carrier,pdsch);
    dmrsIndices = nrPDSCHDMRSIndices(carrier,pdsch);
    pdschGrid = nrResourceGrid(carrier,nTxAnts);
    [~,pdschAntIndices] = nrExtractResources(pdschIndices,pdschGrid);
    pdschGrid(pdschAntIndices) = pdschSymbolsPrecoded;
    % PDSCH DM-RS precoding and mapping
    for p = 1:size(dmrsSymbols,2)
        [~,dmrsAntIndices] = nrExtractResources(dmrsIndices(:,p),pdschGrid);
        pdschGrid(dmrsAntIndices) = pdschGrid(dmrsAntIndices) + dmrsSymbols(:,p)*precodingWeights(p,:);
    end
    [txWaveform,waveformInfo] = nrOFDMModulate(carrier,pdschGrid);
    padding = zeros(length(txWaveform),1);
    txWaveform_send = [padding;txWaveform;padding];
    plot(abs(txWaveform_send*10))
    Wave2File("tx.bin",txWaveform_send*10)

%     % Propagation Channel
%     chInfo = info(channel);
%     maxChDelay = ceil(max(chInfo.PathDelays*channel.SampleRate)) + chInfo.ChannelFilterDelay;
%     txWaveform = [txWaveform; zeros(maxChDelay,size(txWaveform,2))];
%     
%     [rxWaveform,pathGains,sampleTimes] = channel(txWaveform);
%     noise = generateAWGN(SNRdB,nRxAnts,waveformInfo.Nfft,size(rxWaveform));
%     rxWaveform = rxWaveform + noise;
% 
%     % Timing Synchronization
%     if perfectEstimation
%         % Get path filters for perfect timing estimation
%         pathFilters = getPathFilters(channel); 
%         [offset,mag] = nrPerfectTimingEstimate(pathGains,pathFilters);
%     else
%         [t,mag] = nrTimingEstimate(carrier,rxWaveform,dmrsIndices,dmrsSymbols);
%         offset = hSkipWeakTimingOffset(offset,t,mag);
%     end
%     rxWaveform = rxWaveform(1+offset:end,:);
%     % OFDM Demodulation
%     rxGrid = nrOFDMDemodulate(carrier,rxWaveform);
%     % Channel Estimation
%     if perfectEstimation
%         % Perform perfect channel estimation between transmit and receive
%         % antennas.
%         estChGridAnts = nrPerfectChannelEstimate(carrier,pathGains,pathFilters,offset,sampleTimes);
% 
%         % Get perfect noise estimate (from noise realization)
%         noiseGrid = nrOFDMDemodulate(carrier,noise(1+offset:end ,:));
%         noiseEst = var(noiseGrid(:));
% 
%         % Get precoding matrix for next slot
%         newPrecodingWeight = getPrecodingMatrix(pdsch.PRBSet,pdsch.NumLayers,estChGridAnts);
% 
%         % Apply precoding to estChGridAnts. The resulting estimate is for
%         % the channel estimate between layers and receive antennas.
%         estChGridLayers = precodeChannelEstimate(estChGridAnts,precodingWeights.');
%     else
%         % Perform practical channel estimation between layers and receive
%         % antennas.
%         [estChGridLayers,noiseEst] = nrChannelEstimate(carrier,rxGrid,dmrsIndices,dmrsSymbols,'CDMLengths',pdsch.DMRS.CDMLengths);
% 
%         % Remove precoding from estChannelGrid before precoding
%         % matrix calculation
%         estChGridAnts = precodeChannelEstimate(estChGridLayers,conj(precodingWeights));
% 
%         % Get precoding matrix for next slot
%         newPrecodingWeight = getPrecodingMatrix(pdsch.PRBSet,pdsch.NumLayers,estChGridAnts);
%     end
% 
%     % Equalization
%     [pdschRx,pdschHest] = nrExtractResources(pdschIndices,rxGrid,estChGridLayers);
%     [pdschEq,csi] = nrEqualizeMMSE(pdschRx,pdschHest,noiseEst);
%     constPlot.ChannelNames = "Layer "+(pdsch.NumLayers:-1:1);
%     constPlot.ShowLegend = true;
%     % Constellation for the first layer has a higher SNR than that for the
%     % last layer. Flip the layers so that the constellations do not mask
%     % each other.
%     constPlot(fliplr(pdschEq));
%     % PDSCH Decoding
%     [dlschLLRs,rxSymbols] = nrPDSCHDecode(carrier,pdsch,pdschEq,noiseEst);
%     % Scale LLRs by CSI
%     csi = nrLayerDemap(csi);                                    % CSI layer demapping
%     for cwIdx = 1:pdsch.NumCodewords
%         Qm = length(dlschLLRs{cwIdx})/length(rxSymbols{cwIdx}); % Bits per symbol
%         csi{cwIdx} = repmat(csi{cwIdx}.',Qm,1);                 % Expand by each bit per symbol
%         dlschLLRs{cwIdx} = dlschLLRs{cwIdx} .* csi{cwIdx}(:);   % Scale
%     end
%     % DL-SCH Decoding
%     decodeDLSCH.TransportBlockLength = trBlkSizes;
%     [decbits,blkerr] = decodeDLSCH(dlschLLRs,pdsch.Modulation,pdsch.NumLayers, ...
%         harqEntity.RedundancyVersion,harqEntity.HARQProcessID);
%     % HARQ Process Update
%     statusReport = updateAndAdvance(harqEntity,blkerr,trBlkSizes,pdschInfo.G);    
%     disp("Slot "+(nSlot)+". "+statusReport);
% end % for nSlot = 0:totalNoSlots
% 
% 
% 
% function noise = generateAWGN(SNRdB,nRxAnts,Nfft,sizeRxWaveform)
% % Generate AWGN for a given value of SNR in dB (SNRDB), which is the
% % receiver SNR per RE and antenna, assuming the channel does
% % not affect the power of the signal. NRXANTS is the number of receive
% % antennas. NFFT is the FFT size used in OFDM demodulation. SIZERXWAVEFORM
% % is the size of the receive waveform used to calculate the size of the
% % noise matrix.
% 
%     % Normalize noise power by the IFFT size used in OFDM modulation, as
%     % the OFDM modulator applies this normalization to the transmitted
%     % waveform. Also normalize by the number of receive antennas, as the
%     % channel model applies this normalization to the received waveform by
%     % default. The SNR is defined per RE for each receive antenna (TS
%     % 38.101-4).
%     SNR = 10^(SNRdB/10); % Calculate linear noise gain
%     N0 = 1/sqrt(2.0*nRxAnts*double(Nfft)*SNR);
%     noise = N0*complex(randn(sizeRxWaveform),randn(sizeRxWaveform));
% end
    
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