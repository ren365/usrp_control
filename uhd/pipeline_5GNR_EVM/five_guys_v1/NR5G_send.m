% settings = {};
% settings.Modulation = "16QAM"; % ["QPSK","16QAM","64QAM"]
% settings.CodeRate = 490/1024;
% settings.SubcarrierSpacing = 15; % 15, 30, 60, 120, or 240 kHz
% settings.NSizeGrid= 52; % 24 to 275
% % check https://www.sharetechnote.com/html/5G/5G_ResourceGrid.html 
% % SISO for now, no MIMO support
% settings.NumLayers = 1;
% settings.nTxAnts = 1;
% settings.nRxAnts = 1;
% settings.generateTX = true;
% settings.plotTX = false;

function txsetting = NR5G_send(settings)
    txsetting = {};
    nSlot = 0;
    totalNoSlots = 1;%20;         % Number of slots to simulate
    perfectEstimation = false; % Perfect synchronization and channel estimation
    rng("default");            % Set default random number generator for repeatability
    carrier = nrCarrierConfig;
    carrier.NSizeGrid = settings.NSizeGrid;
    carrier.SubcarrierSpacing = settings.SubcarrierSpacing;
    % PDSCH and DM-RS Configuration
    pdsch = nrPDSCHConfig;
    pdsch.Modulation = settings.Modulation;
    pdsch.NumLayers = settings.NumLayers;
    pdsch.PRBSet = 0:carrier.NSizeGrid-1;     % Full band allocation
    pdsch.DMRS.DMRSAdditionalPosition = 1;
    pdsch.DMRS.DMRSConfigurationType = 1;
    pdsch.DMRS.DMRSLength = 2;
    % DL-SCH Configuration
    NHARQProcesses = 16;     % Number of parallel HARQ processes
    rvSeq = [0 2 3 1];
    % Coding rate
    if pdsch.NumCodewords == 1
        codeRate = settings.CodeRate;
    else
        codeRate = [settings.CodeRate, settings.CodeRate];
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
    nTxAnts = settings.nTxAnts;
    nRxAnts = settings.nRxAnts;
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

    if settings.plotTX
        plot(abs(txWaveform_send*10))
    end

    if settings.generateTX
        Wave2File("tx.bin",txWaveform_send*10)
    end

    % record everything and return 
    txsetting.trBlk = trBlk;
    txsetting.txWaveform = txWaveform;
    txsetting.carrier = carrier;
    txsetting.dmrsIndices = dmrsIndices;
    txsetting.dmrsSymbols = dmrsSymbols;
    txsetting.pdsch = pdsch;
    txsetting.pdschIndices = pdschIndices;
    txsetting.precodingWeights = precodingWeights;
    txsetting.trBlkSizes = trBlkSizes;
    txsetting.harqEntity = harqEntity;
    txsetting.pdschInfo = pdschInfo;
    txsetting.pdschSymbols = pdschSymbols;
    txsetting.decodeDLSCH = decodeDLSCH;

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