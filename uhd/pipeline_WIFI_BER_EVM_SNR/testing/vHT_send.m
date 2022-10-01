% postprocess to get BER (and EVM)
rng(0);

    sendDatas = [true,false];
    sendData = sendDatas(1);

    scale_alpha = 0.01;
    cbw = 'CBW20';                     % Channel bandwidth
    ntx = 1;                           % Number of transmit antennas
    nsts = 1;                          % Number of space-time streams
%     nrx = 1;                           % Number of receive antennas
    mcs = 0;
    filePathBin = "Tx/";
    filePathMat = "Label/";
    protStr = "VHT";
    psdu = 2000;
    fileName = protStr + '_' + cbw +'_code_'+mcs+ '_length_' + psdu;
    psduData = randi(2, [psdu*8 1]) - 1;

    cfgVHT = wlanVHTConfig( ...
            "ChannelBandwidth", cbw, ...
            "NumUsers", 1, ...
            "NumTransmitAntennas", ntx, ...
            "NumSpaceTimeStreams", nsts, ...
            "MCS", mcs, ...
            "APEPLength", psdu);

    wave  = wlanWaveformGenerator(psduData, cfgVHT);

if sendData
    wave = scale_alpha* wave;
    waveLen=length(wave);
    Wave2File(filePathBin+fileName+".bin", wave);
    save(filePathMat+fileName+".mat", ...
                     "protStr", "cbw", "mcs", "psdu", "psduData", "waveLen");
    return;
end

rxPPDU = File2Wave("Cable/"+fileName+".bin");
% rxPPDU = File2Wave("Cable/"+fileName+".bin");
% pack offset = 0
fieldInd = wlanFieldIndices(cfgVHT);
fs = wlanSampleRate(cfgVHT);

coarsePktOffset = wlanPacketDetect(rxPPDU,cfgVHT.ChannelBandwidth);
lstf = rxPPDU(coarsePktOffset+(fieldInd.LSTF(1):fieldInd.LSTF(2)),:);
coarseFreqOff = wlanCoarseCFOEstimate(lstf,cfgVHT.ChannelBandwidth);
rxPPDU = frequencyOffset(rxPPDU,fs,-coarseFreqOff);

nonhtfields = rxPPDU(coarsePktOffset+(fieldInd.LSTF(1):fieldInd.LSIG(2)),:);
finePktOffset = wlanSymbolTimingEstimate(nonhtfields,...
    cfgVHT.ChannelBandwidth);

return
pktOffset = coarsePktOffset+finePktOffset;
%         % Extract L-STF and perform coarse frequency offset correction
lltf = rxPPDU(pktOffset+(fieldInd.LLTF(1):fieldInd.LLTF(2)),:);
fineFreqOff = wlanFineCFOEstimate(lltf,cbw);
rxPPDU = frequencyOffset(rxPPDU,fs,-fineFreqOff);

% plot the preamble
numSamples = fieldInd.VHTSIGB(2);
time = ([0:double(numSamples)-1]/fs)*1e6;
peak = 1.2*max(abs(rxPPDU(1:numSamples)));
fieldMarkers = zeros(numSamples,1);
fieldMarkers(fieldInd.LSTF(2)-1,1) = peak;
fieldMarkers(fieldInd.LLTF(2)-1,1) = peak;
fieldMarkers(fieldInd.LSIG(2)-1,1) = peak;
fieldMarkers(fieldInd.VHTSIGA(2)-1,1) = peak;
fieldMarkers(fieldInd.VHTSTF(2)-1,1) = peak;
fieldMarkers(fieldInd.VHTLTF(2)-1,1) = peak;
fieldMarkers(fieldInd.VHTSIGB(2)-1,1) = peak;
plot(time,abs(rxPPDU(1:numSamples)),time,fieldMarkers)
xlabel ('Time (microseconds)')
ylabel('Magnitude')
title('VHT Preamble')

%
rxLLTF = rxPPDU(pktOffset+(fieldInd.LLTF(1):fieldInd.LLTF(2)),:);
demodLLTF = wlanLLTFDemodulate(rxLLTF,cfgVHT);
chEstLLTF = wlanLLTFChannelEstimate(demodLLTF,cfgVHT);

%  L-SIG field from the received PPDU
noiseVar = 8e-12;
rxLSIG = rxPPDU(pktOffset+(fieldInd.LSIG(1):fieldInd.LSIG(2)),:);
[recLSIG,failCRC] = wlanLSIGRecover(rxLSIG,chEstLLTF,noiseVar,cfgVHT.ChannelBandwidth);
recLSIG(end-5:end)',recLSIG(1:4)'

% VHT-SIG-A - MCS
rxVHTSIGA = rxPPDU(fieldInd.VHTSIGA(1):fieldInd.VHTSIGA(2),:);
[recVHTSIGA,failCRC] = wlanVHTSIGARecover(rxVHTSIGA, ...
    chEstLLTF,noiseVar,cfgVHT.ChannelBandwidth);
recMCSbits = (recVHTSIGA(29:32))';
recMCS = bi2de(double(recMCSbits));
isequal(recMCS,cfgVHT.MCS) % check MCS

% Extract and demodulate the VHT-LTF.
rxVHTLTF = rxPPDU(fieldInd.VHTLTF(1):fieldInd.VHTLTF(2),:);
demodVHTLTF = wlanVHTLTFDemodulate(rxVHTLTF,cfgVHT);
chEstVHTLTF = wlanVHTLTFChannelEstimate(demodVHTLTF,cfgVHT);
% Extract and recover the VHT-SIG-B.
rxVHTSIGB = rxPPDU(fieldInd.VHTSIGB(1):fieldInd.VHTSIGB(2),:);
recVHTSIGB = wlanVHTSIGBRecover(rxVHTSIGB,chEstVHTLTF,noiseVar,cfgVHT.ChannelBandwidth);
%  the value in the VHT-SIG-B Length field multiplied by 4 is the recovered 
% APEP length for packets carrying data.
% Verify that the APEP length, contained in the first 19 bits of the VHT-SIG-B, 
% corresponds to the specified APEP length.
sigbAPEPbits = recVHTSIGB(1:19)';
sigbAPEPlength = bi2de(double(sigbAPEPbits))*4;
isequal(sigbAPEPlength,cfgVHT.APEPLength)

recPSDU = wlanVHTDataRecover(rxPPDU(fieldInd.VHTData(1):fieldInd.VHTData(2),:),...
    chEstVHTLTF,noiseVar,cfgVHT);

% Determine if any bits are in error, i.e. a packet error
[numErr,ratio] = biterr(psduData,recPSDU);
ratio


