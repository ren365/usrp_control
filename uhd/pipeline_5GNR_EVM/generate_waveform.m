zeropadding = zeros(1e5,1);
waveform = waveStruct.waveform;
waveform_new = [zeropadding;waveform;zeropadding];
Wave2File("5GNR_sent_withPadding.bin",waveform_new*15);