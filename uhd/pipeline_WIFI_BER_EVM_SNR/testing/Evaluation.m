% postprocess to get BER (and EVM)
rng(0);

    cbw = 'CBW5';                     % Channel bandwidth
    ntx = 1;                           % Number of transmit antennas
    nsts = 1;                          % Number of space-time streams
    nrx = 1;                           % Number of receive antennas
    mcs = 7;
    filePathBin = "Tx/";
    filePathMat = "Label/";
    protStr = "Non-HT";
    psdu = 3337;
    fileName = protStr + '_' + cbw + '_' + psdu;
    psduData = randi(2, [psdu*8 1]) - 1;
    cfgHT = wlanNonHTConfig('ChannelBandwidth',cbw, ... 
        'MCS',mcs,...
        'PSDULength',psdu);
    wave = wlanWaveformGenerator(psduData,cfgHT);
%     [wave,cfgHT] = GetWaveform(protStr, cbw, mcs, psdu, ntx, psduData);

if true
    wave = 0.01* wave;
    waveLen=length(wave);
    Wave2File(filePathBin+fileName+".bin", wave);
    save(filePathMat+fileName+".mat", ...
                     "protStr", "cbw", "mcs", "psdu", "psduData", "waveLen");
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
psduDataStrucrt = load(filePathMat+fileName+".mat","psduData");
psduData = getfield(psduDataStrucrt,"psduData");
rxPPDU = File2Wave("Cable/"+fileName+".bin");
% data1 = File2Wave("Tx/HT_CBW20_16_4_3337_1.bin");
% biterr(data0,data1)
% Find the start and stop indices for all component fields of the PPDU.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fieldInd = wlanFieldIndices(cfgHT);
fs = wlanSampleRate(cfgHT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
rxLLTF = rxPPDU(fieldInd.LLTF(1):fieldInd.LLTF(2),:);
demodLLTF = wlanLLTFDemodulate(rxLLTF,cfgHT);
chEstLLTF = wlanLLTFChannelEstimate(demodLLTF,cfgHT);

noiseVar = 0.1;
rxLSIG = rxPPDU(fieldInd.LSIG(1):fieldInd.LSIG(2),:);
[recLSIG,failCRC] = wlanLSIGRecover(rxLSIG,chEstLLTF,noiseVar,cfgHT.ChannelBandwidth);
failCRC
return


% Packet detect and determine coarse packet offset
        coarsePktOffset = wlanPacketDetect(rx,cbw);
        if isempty(coarsePktOffset) % If empty no L-STF detected; packet error
            disp("packet error")
        end

%         % Extract L-STF and perform coarse frequency offset correction
        lstf = rx(coarsePktOffset+(ind.LSTF(1):ind.LSTF(2)),:);
        coarseFreqOff = wlanCoarseCFOEstimate(lstf,cbw);
        rx = frequencyOffset(rx,fs,-coarseFreqOff);

% 
%         % Extract the non-HT fields and determine fine packet offset
        nonhtfields = rx(coarsePktOffset+(ind.LSTF(1):ind.LSIG(2)),:);
        finePktOffset = wlanSymbolTimingEstimate(nonhtfields,cbw);
% 
%         % Determine final packet offset
        pktOffset = coarsePktOffset+finePktOffset;
% 
%         % If packet detected outwith the range of expected delays from the
%         % channel modeling; packet error
%         if pktOffset>15
%             disp("packet error 2");
%         end
       

        % Extract L-LTF and perform fine frequency offset correction
        lltf = rx(pktOffset+(ind.LLTF(1):ind.LLTF(2)),:);
        fineFreqOff = wlanFineCFOEstimate(lltf,cbw);
        rx = frequencyOffset(rx,fs,-fineFreqOff);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Extract HT-LTF samples from the waveform, demodulate and perform
        % channel estimation
        htltf = rx(pktOffset+(ind.HTLTF(1):ind.HTLTF(2)),:);
        htltfDemod = wlanLLTFDemodulate(htltf,ht); % wlanNonHTLTFDemodulate(htltf,ht);
        chanEst = wlanLLTFChannelEstimate(htltfDemod,ht);

        % Extract HT Data samples from the waveform
        htdata = rx(pktOffset+(ind.HTData(1):ind.HTData(2)),:);

        % Estimate the noise power in HT data field
        nVarHT = htNoiseEstimate(htdata,chanEst,ht);

        % Recover the transmitted PSDU in HT Data
        rxPSDU = wlanHTDataRecover(htdata,chanEst,nVarHT,ht);

        % Determine if any bits are in error, i.e. a packet error
       [numErr,ratio] = biterr(psduData,rxPSDU)
        
%        waveTemp = rx(pktOffset+(ind.LSIG(1): ind.LSIG(2)));
%        recBits = wlanLSIGRecover(waveTemp, chanEst, nVarHT, "CBW20");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Calculate the number of bit errors in the received packet.
% [numErr,ratio] = biterr(psduData,rxPSDU)

