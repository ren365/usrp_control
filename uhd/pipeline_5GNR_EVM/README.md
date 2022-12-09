# Introduction to Matlab-based 5GNR pipeline [internal use]

## Highlight
1. Compatible with Universal Software Radio Peripheral (USRP)
2. Providing SNR, EVM and BER averaged on multiple packets

## How to run
1. make proper configurations in "NR5G_setting.m"
2. run "five_guys_tx.m" with matlab, which generates "tx.bin" waveform
3. [with usrp connected] run "tx2rx.sh" with proper bandwidth and number of samples configured. After transmission there will be "rx.bin" in the folder.
4. [option] using "rx_test.m" to quickly check with the received waveform
5. run "five_guys_rx.m" to get SNR/EVM/BER

## Folder structure
- "File2Wave.m" helper file to convert usrp received waveform into matlab readable waveform
- "Wave2File.m" helper file to convert matlab readable waveform into usrp received waveform
- "NR5G_send.m" helper file to generate tx waveform
- "NR5G_setting.m" helper file to configure tx waveform