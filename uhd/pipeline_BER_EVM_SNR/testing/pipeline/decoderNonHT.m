function [rxPSDU,beforeDemap]=decoderNonHT(fileName,cfgNonHT)
rxPPDU = File2Wave(fileName);

fieldInd = wlanFieldIndices(cfgNonHT);
fs = wlanSampleRate(cfgNonHT);
cbw = cfgNonHT.ChannelBandwidth;

pktOffset = 0;lltf = rxPPDU(pktOffset+(fieldInd.LLTF(1):fieldInd.LLTF(2)),:);
fineFreqOff = wlanFineCFOEstimate(lltf,cbw);
rxPPDU = frequencyOffset(rxPPDU,fs,-fineFreqOff);

rxLLTF = rxPPDU(pktOffset+(fieldInd.LLTF(1):fieldInd.LLTF(2)),:);
demodLLTF = wlanLLTFDemodulate(rxLLTF,cfgNonHT);
chEstLLTF = wlanLLTFChannelEstimate(demodLLTF,cfgNonHT);

noiseVar = 8e-12;
rxLSIG = rxPPDU(pktOffset+(fieldInd.LSIG(1):fieldInd.LSIG(2)),:);
[recLSIG,failCRC] = wlanLSIGRecover(rxLSIG,chEstLLTF,noiseVar,cfgNonHT.ChannelBandwidth);
recLSIG(end-5:end)'
recLSIG(1:4)'

rxNonHTData = rxPPDU(fieldInd.NonHTData(1):fieldInd.NonHTData(2),:);
 % Recover the transmitted PSDU in HT Data
[rxPSDU, eqDataSym, varargout,beforeDemap] = wlanNonHTDataRecover_Debug(rxNonHTData,chEstLLTF,noiseVar,cfgNonHT);

end