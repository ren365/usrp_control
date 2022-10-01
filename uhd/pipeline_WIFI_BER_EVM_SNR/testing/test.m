% postprocess to get BER (and EVM)
addpath("../../Library")
rng(0);

    cbw = 'CBW20';                     % Channel bandwidth
    ntx = 1;                           % Number of transmit antennas
    nsts = 1;                          % Number of space-time streams
    nrx = 1;                           % Number of receive antennas
    mcs = 4;
    filePathBin = "Tx/";
    filePathMat = "Label/";
    protStr = "Non-HT";
    psdu = 2000;
    fileName = protStr + '_' + cbw + '_' + psdu;
    psduData = randi(2, [psdu*8 1]) - 1;
    [wave,cfgHT] = GetWaveform(protStr, cbw, mcs, psdu, ntx, psduData);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% psduDataStrucrt = load(filePathMat+fileName+".mat","psduData");
% psduData = getfield(psduDataStrucrt,"psduData");
% rxPPDU = File2Wave("Cable/"+fileName+".bin");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get RX data after transmit by USRP/UHD
rxPPDU = File2Wave("Cable/Non-HT_CBW20_2000.bin");
% wave0      = File2Wave(filePathBin+fileName+".bin");

fieldInd = wlanFieldIndices(cfgHT);
fs = wlanSampleRate(cfgHT);

% Extract L-LTF and perform fine frequency offset correction
% coarsePktOffset = wlanPacketDetect(rxPPDU,cbw);
% pktOffset = coarsePktOffset; % already did in Zhihui's python file - coarsePktOffset-500
pktOffset = 0;lltf = rxPPDU(pktOffset+(fieldInd.LLTF(1):fieldInd.LLTF(2)),:);
fineFreqOff = wlanFineCFOEstimate(lltf,cbw);
rxPPDU = frequencyOffset(rxPPDU,fs,-fineFreqOff);

rxLLTF = rxPPDU(pktOffset+(fieldInd.LLTF(1):fieldInd.LLTF(2)),:);
demodLLTF = wlanLLTFDemodulate(rxLLTF,cfgHT);
chEstLLTF = wlanLLTFChannelEstimate(demodLLTF,cfgHT);

noiseVar = 8e-12;
rxLSIG = rxPPDU(pktOffset+(fieldInd.LSIG(1):fieldInd.LSIG(2)),:);
[recLSIG,failCRC] = wlanLSIGRecover(rxLSIG,chEstLLTF,noiseVar,cfgHT.ChannelBandwidth);
recLSIG(end-5:end)'
recLSIG(1:4)'

lltf = rxPPDU(pktOffset+(fieldInd.LLTF(1):fieldInd.LLTF(2)),:);
fineFreqOff = wlanFineCFOEstimate(lltf,cbw);
rxPPDU = frequencyOffset(rxPPDU,fs,-fineFreqOff);

rxLLTF = rxPPDU(pktOffset+(fieldInd.LLTF(1):fieldInd.LLTF(2)),:);
demodLLTF = wlanLLTFDemodulate(rxLLTF,cfgHT);
chEstLLTF = wlanLLTFChannelEstimate(demodLLTF,cfgHT);

noiseVar = 8e-12;
rxLSIG = rxPPDU(pktOffset+(fieldInd.LSIG(1):fieldInd.LSIG(2)),:);
[recLSIG,failCRC] = wlanLSIGRecover(rxLSIG,chEstLLTF,noiseVar,cfgHT.ChannelBandwidth);
recLSIG(end-5:end)'
recLSIG(1:4)'


rxNonHTData = rxPPDU(fieldInd.NonHTData(1):fieldInd.NonHTData(2),:);
 % Recover the transmitted PSDU in HT Data
rxPSDU = wlanNonHTDataRecover(rxNonHTData,chEstLLTF,noiseVar,ht);

% Determine if any bits are in error, i.e. a packet error
[numErr,ratio] = biterr(psduData,rxPSDU)


return
% rxPPDU = wave;pktOffset = 0;

rxPPDU = File2Wave("Tx/Non-HT_CBW20_1000.bin");

rxPPDU = File2Wave("Cable/Non-HT_CBW20_2000.bin");
% rxPPDU = wave;
rxPlot =abs(rxPPDU(pktOffset+(fieldInd.LSTF(1):fieldInd.LSIG(2)),:)) ;
plot(rxPlot)
minD = min(rxPlot)
maxD = max(rxPlot)
hold on
line([fieldInd.LSTF(2),fieldInd.LSTF(2)],[minD,maxD],'Color','r','LineWidth',1.2)
line([fieldInd.LLTF(2),fieldInd.LLTF(2)],[minD,maxD],'Color','r','LineWidth',1.2)
line([fieldInd.LSIG(2),fieldInd.LSIG(2)],[minD,maxD],'Color','r','LineWidth',1.2)
xline(fieldInd.LSTF(2),'Color','r')
xline(fieldInd.LLTF(2),'Color','r')
xline(fieldInd.LSIG(2),'Color','r')




% Packet detect and determine coarse packet offset
        coarsePktOffset = wlanPacketDetect(rxPPDU,cbw);
        if isempty(coarsePktOffset) % If empty no L-STF detected; packet error
            disp("packet error")
        end

%         % Extract L-STF and perform coarse frequency offset correction
        lstf = rxPPDU(coarsePktOffset+(fieldInd.LSTF(1):fieldInd.LSTF(2)),:);
        coarseFreqOff = wlanCoarseCFOEstimate(lstf,cbw);
        rxPPDU = frequencyOffset(rxPPDU,fs,-coarseFreqOff);

        nonhtfields = rxPPDU(coarsePktOffset+(fieldInd.LSTF(1):fieldInd.LSIG(2)),:);
        finePktOffset = wlanSymbolTimingEstimate(nonhtfields,cbw);
        pktOffset = coarsePktOffset+finePktOffset;

        % Extract L-LTF and perform fine frequency offset correction
        lltf = rxPPDU(pktOffset+(fieldInd.LLTF(1):fieldInd.LLTF(2)),:);
        fineFreqOff = wlanFineCFOEstimate(lltf,cbw);
        rxPPDU = frequencyOffset(rxPPDU,fs,-fineFreqOff);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Extract HT-LTF samples from the waveform, demodulate and perform
        % channel estimation
        htltf = rxPPDU(pktOffset+(fieldInd.HTLTF(1):fieldInd.HTLTF(2)),:);
        htltfDemod = wlanLLTFDemodulate(htltf,ht); % wlanNonHTLTFDemodulate(htltf,ht);
        chanEst = wlanLLTFChannelEstimate(htltfDemod,ht);

        % Extract HT Data samples from the waveform
        htdata = rxPPDU(pktOffset+(fieldInd.HTData(1):fieldInd.HTData(2)),:);

        % Estimate the noise power in HT data field
        nVarHT = htNoiseEstimate(htdata,chanEst,ht);

        % Recover the transmitted PSDU in HT Data
        rxPSDU = wlanHTDataRecover(htdata,chanEst,nVarHT,ht);

        % Determine if any bits are in error, i.e. a packet error
       [numErr,ratio] = biterr(psduData,rxPSDU)