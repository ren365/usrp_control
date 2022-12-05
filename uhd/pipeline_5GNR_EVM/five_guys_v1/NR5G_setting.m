function settings = NR5G_setting()

    settings = {};
    settings.Modulation = "16QAM"; % ["QPSK","16QAM","64QAM"]
    settings.CodeRate = 490/1024;
    settings.SubcarrierSpacing = 15; % 15, 30, 60, 120, or 240 kHz
    settings.NSizeGrid= 52; % 24 to 275
    % check https://www.sharetechnote.com/html/5G/5G_ResourceGrid.html 
    % SISO for now, no MIMO support
    settings.NumLayers = 1;
    settings.nTxAnts = 1;
    settings.nRxAnts = 1;
    settings.generateTX = true;
    settings.plotTX = false;
    
    % setting for decoding
    settings.packetNum = 50;
    settings.fastEVM = true;
end