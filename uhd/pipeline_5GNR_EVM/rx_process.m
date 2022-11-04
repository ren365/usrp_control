close all

waveform_rx = File2Wave("5GNR_sent_withPadding.bin");
% waveform_rx = File2Wave("5GNR_rx.bin");
% waveform_rx = File2Wave("rx_file.bin");

% waveform_rx_select = fftshift(fft(waveform_rx));

% waveform_rx = File2Wave("rx_file.bin");
waveform_rx_select = waveform_rx(end-2.7e5:end,1);
% waveform_rx_select = waveform_rx(end-2.5e5:end-1e5);

figure(1);
plot(log(abs(waveform_rx_select)));

figure(2);
plot(real(waveform_rx_select));