close all

waveform_rx = File2Wave("rx.bin");
% waveform_rx = File2Wave("5GNR_rx.bin");
% waveform_rx = File2Wave("rx_file.bin");

% waveform_rx_select = fftshift(fft(waveform_rx));

% waveform_rx = File2Wave("rx_file.bin");
waveform_rx_select = waveform_rx;%(end-4.05e4:end,1);
% waveform_rx_select = waveform_rx(end-2.5e5:end-1e5);

figure(1);
plot(abs(waveform_rx_select));
% plot(20*log10(abs(waveform_rx_select)));

figure(2);
% plot(real(fft(waveform_rx_select)));

waveform_tx = File2Wave("tx.bin");
plot(abs(waveform_tx));