./rx_samples_to_file --help
serialRx="321FCC3"
fileRx="Tx/output.bin"
rateRx=20000000
freqRx=1000000000
gainRx=40
antRx="RX2"
channelRX=1
bandRx=20000000
./rx_samples_to_file --args serial=$serialRx --file $fileRx --type float --nsamps 0 --rate $rateRx --freq $freqRx --gain $gainRx --ant $antRx --channel $channelRX --bw $bandRx

# ./rx --args serial=3231EFC --file output.bin --type float --nsamps 0 --rate 20000000 --freq 1000000000 --gain 40 --ant RX2 --subdev A:A --bw 20000000
