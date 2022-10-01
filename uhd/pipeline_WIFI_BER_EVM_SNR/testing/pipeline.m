% matlab + wlan + UHD executable file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Transmit
cbw = 'CBW20';                     % Channel bandwidth
fs = 50e6;                         % Sample rate (Hz)
ntx = 1;                           % Number of transmit antennas
nsts = 1;                          % Number of space-time streams
nrx = 1;                           % Number of receive antennas

vht = wlanVHTConfig('ChannelBandwidth',cbw,'APEPLength',20, ...
    'NumTransmitAntennas',ntx,'NumSpaceTimeStreams',nsts, ...
    'SpatialMapping','Direct','STBC',false);
% vht.PSDULength = 20;
txPSDU = randi([0 1],vht.PSDULength*8,1);
txPPDU = wlanWaveformGenerator(txPSDU,vht);
Wave2File("Tx/send.bin",txPPDU)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Replace by OTA
% Create a 2x2 TGac channel and an AWGN channel.
tgacChan = wlanTGacChannel('SampleRate',fs,'ChannelBandwidth',cbw, ...
    'NumTransmitAntennas',ntx,'NumReceiveAntennas',nrx, ...
    'LargeScaleFadingEffect','Pathloss and shadowing', ...
    'DelayProfile','Model-C');
awgnChan = comm.AWGNChannel('NoiseMethod','Variance', ...
    'VarianceSource','Input port');
% Calculate the noise variance for a receiver with a 9 dB noise figure. Pass the transmitted waveform through the noisy TGac channel.
pfOffset = comm.PhaseFrequencyOffset('SampleRate',fs,'FrequencyOffsetSource','Input port');
nVar = 10^((-228.6 + 10*log10(290) + 10*log10(fs) + 9)/10);
rxPPDU = awgnChan(tgacChan(txPPDU), nVar);
rxPPDUcfo = pfOffset(rxPPDU,500);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Received Evaluation

% Find the start and stop indices for all component fields of the PPDU.
ind = wlanFieldIndices(vht);
rxLSTF = rxPPDUcfo(ind.LSTF(1):ind.LSTF(2),:);
% Extract the L-STF. Estimate and correct for the carrier frequency offset.
foffset1 = wlanCoarseCFOEstimate(rxLSTF,cbw);
rxPPDUcorr = pfOffset(rxPPDUcfo,-foffset1);

% Extract the L-LTF from the corrected signal. Estimate and correct for the residual frequency offset.
rxLLTF = rxPPDUcorr(ind.LLTF(1):ind.LLTF(2),:);

foffset2 = wlanFineCFOEstimate(rxLLTF,cbw);
rxPPDU2 = pfOffset(rxPPDUcorr,-foffset2);
% Extract and demodulate the VHT-LTF. Estimate the channel coefficients.
rxVHTLTF = rxPPDU2(ind.VHTLTF(1):ind.VHTLTF(2),:);
dLTF = wlanVHTLTFDemodulate(rxVHTLTF,vht);
chEst = wlanVHTLTFChannelEstimate(dLTF,vht);
% Extract the VHT data field from the received and frequency-corrected PPDU. Recover the data field.
rxVHTData = rxPPDU2(ind.VHTData(1):ind.VHTData(2),:);
rxPSDU = wlanVHTDataRecover(rxVHTData,chEst,nVar,vht);
% Calculate the number of bit errors in the received packet.
numErr = biterr(txPSDU,rxPSDU)
