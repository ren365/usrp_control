clear all; close all; clc;

addpath('/home/zehaow/Documents/MATLAB/Examples/R2022b/5g/EVMMeasurementOfNRDownlinkWaveformsExample');

load('waveform_test_100mhz_1frame.mat');

% Set the NR-TM parameters for the receiver
nrtm = "NR-FR1-TM3.1"; % Reference channel
bw   = "100MHz";        % Channel bandwidth
scs  = "30kHz";        % Subcarrier spacing
dm   = "TDD";          % Duplexing mode

captureSampleRate = 122.88e6;

tmwavegen = hNRReferenceWaveformGenerator(nrtm,bw,scs,dm);
[~,tmwaveinfo,resourcesInfo] = generateWaveform(tmwavegen);

frequencyCorrectionRange = -100e3:1e3:100e3; 

rxWaveform = waveStruct.waveform + 0.0001*randn([length(waveStruct.waveform),1]);


[rxWaveform, coarseOffset] = DMRSFrequencyCorrection(rxWaveform,captureSampleRate,...
    frequencyCorrectionRange,tmwavegen,resourcesInfo); 

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