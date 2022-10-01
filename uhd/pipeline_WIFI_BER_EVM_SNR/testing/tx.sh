./tx_samples_from_file --help
serialTx="322A004"
fileTx="send.bin"
rateTx=50000000
freqTx=1000000000
gainTx=20
antTx="TX/RX"
channelTX=0
bandTx=2000000
./tx_samples_from_file --args serial=$serialTx --file $fileTx --rate $rateTx --freq $freqTx --gain $gainTx --ant $antTx --channel $channelTX --bw $bandTx

# ./tx_samples_from_file --args serial=321FCC3 --file send.bin --type float --rate 20000000 --freq 1000000000 --gain 40 --ant TX/RX --subdev A:A --bw 20000000
