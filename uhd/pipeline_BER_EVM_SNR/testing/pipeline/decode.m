% postprocess to get BER (and EVM)
addpath("../../Library")
rng(0);

cbw = 'CBW20';                     % Channel bandwidth
mcs = 2;
protStr = "Non-HT";
psdu = 3337;
refSym = qammod((0:(mcs - 1))',mcs,UnitAveragePower=1);

cfgNonHT = wlanNonHTConfig( ...
            "Modulation", 'OFDM', ...
            "ChannelBandwidth", cbw, ...
            "MCS", mcs, ...
            "PSDULength", psdu);
psduData = randi(2, [psdu*8 1]) - 1;
% [txWaveform,afterMap] = wlanWaveformGenerator_Debug(psduData,cfgNonHT);

rxFileName = "Rx/Non-HT_CBW20_2_2_3337_1.bin";
txFileName = "Tx/Non-HT_CBW20_2_2_3337_1.bin";
[rxPSDU,rxSym]=decoderNonHT(rxFileName,cfgNonHT);
[txPSDU,txSym]=decoderNonHT(txFileName,cfgNonHT);

% EVM 
constellation = comm.ConstellationDiagram(ReferenceConstellation=refSym(:));
constellation(rxSym(:))

% Determine if any bits are in error, i.e. a packet error
[numErr,ratio] = biterr(txPSDU,rxPSDU)
