#!/bin/bash

# Parameter Setup
downRate=1

serialRx="322A004"
serialTx="322A005"
channelRx=0
channelTx=0
antRx="TX/RX"
antTx="TX/RX"
bandTx=50000000
bandRx=50000000
rateTx=$bandTx
rateRx=$bandRx
gain=10
freq=1000000000
fileTx="tx.bin"
fileRx="rx.bin"


# Transmission
./tx_samples_from_file --args serial=$serialTx --channel $channelTx --ant $antTx --gain $gain --file $fileTx --type float --freq $freq --rate $rateTx --bw $bandTx --repeat &
txpid=$!
sleep 7s
./rx_samples_to_file --args serial=$serialRx --channel $channelRx --ant $antRx --gain $gain --file $fileRx --type float --freq $freq --rate $rateRx --bw $bandRx --nsamps 3000000
kill -9 $txpid
