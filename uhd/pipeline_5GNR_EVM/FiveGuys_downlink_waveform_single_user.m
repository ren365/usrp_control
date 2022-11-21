% function waveform = FiveGuys_downlink_waveform_single_user()

    carrier = nrSCSCarrierConfig('NSizeGrid',100);
    bwp = nrWavegenBWPConfig('NStartBWP',carrier.NStartGrid+10);
    ssb = nrWavegenSSBurstConfig('BlockPattern','Case A');
    pdcch = nrWavegenPDCCHConfig('AggregationLevel',2,'AllocatedCandidate',4);
    coreset = nrCORESETConfig;
    coreset.FrequencyResources = [1 1 1 1];
    coreset.Duration = 3;
    ss = nrSearchSpaceConfig;
    ss.NumCandidates = [8 4 0 0 0];
    pdsch = nrWavegenPDSCHConfig( ...
        'Modulation','16QAM','TargetCodeRate',658/1024,'EnablePTRS',true);
    dmrs = nrPDSCHDMRSConfig('DMRSTypeAPosition',3);
    pdsch.DMRS = dmrs;
    ptrs = nrPDSCHPTRSConfig('TimeDensity',2);
    pdsch.PTRS = ptrs;
    csirs = nrWavegenCSIRSConfig('RowNumber',4,'RBOffset',10,'NumRB',10,'SymbolLocations',5);
    cfgDL = nrDLCarrierConfig( ...
        'FrequencyRange','FR1', ...
        'ChannelBandwidth',40, ...
        'NumSubframes',20, ...
        'SCSCarriers',{carrier}, ...
        'BandwidthParts',{bwp}, ...
        'SSBurst',ssb, ...
        'CORESET',{coreset}, ...
        'SearchSpaces',{ss}, ...
        'PDCCH',{pdcch}, ...
        'PDSCH',{pdsch}, ...
        'CSIRS',{csirs});
    [waveform, info] = nrWaveformGenerator_wzh(cfgDL);

%     data = info.

%%%%% demode
% Set the NR-TM parameters for the receiver
nrtm = "NR-FR1-TM3.1"; % Reference channel
bw   = "100MHz";        % Channel bandwidth
scs  = "30kHz";        % Subcarrier spacing
dm   = "TDD";          % Duplexing mode

captureSampleRate = 122.88e6;

tmwavegen = hNRReferenceWaveformGenerator(nrtm,bw,scs,dm);
[txWaveform,tmwaveinfo,resourcesInfo] = generateWaveform(tmwavegen);
rxWaveform = waveform;
frequencyCorrectionRange = -100e3:1e3:100e3; 

% rxWaveform = waveStruct.waveform + 0.0001*randn([length(waveStruct.waveform),1]);
[rxWaveform, coarseOffset] = DMRSFrequencyCorrection(rxWaveform,captureSampleRate,...
    frequencyCorrectionRange,tmwavegen,info); 

fprintf("Coarse frequency offset = %.0f Hz", coarseOffset) 

frequencyCorrectionRange = -100:5:100; 

[rxWaveform, fineOffset] = DMRSFrequencyCorrection(rxWaveform,captureSampleRate,...
    frequencyCorrectionRange,tmwavegen,resourcesInfo); 

fprintf("Fine frequency offset = %.1f Hz", fineOffset) 

cfg = struct();
cfg.PlotEVM = true;                 % Plot EVM statistics
cfg.DisplayEVM = true;              % Print EVM statistics
cfg.Label = nrtm;                   % Set to TM name of captured waveform
cfg.SampleRate = captureSampleRate; % Use sample rate during capture

[evmInfo,eqSym,refSym] = hNRPDSCHEVM(tmwavegen.Config,rxWaveform,cfg);

modulation='64QAM';
ncellid=1;
rnti=1;
rxbits = nrPDSCHDecode(eqSym,modulation,ncellid,rnti);



% end