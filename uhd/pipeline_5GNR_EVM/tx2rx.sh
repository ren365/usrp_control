#!/bin/bash

# Parameter Setup
downRate=1

serialRx="325A876"
serialTx="321FCC3"
channelRx=0
channelTx=0
antRx="TX/RX"
antTx="TX/RX"
bandTx=100000000
bandRx=100000000
rateTx=$bandTx
rateRx=$bandRx
gain=30
freq=3500000000
fileTx="5GNR_sent_withPadding.bin"
fileRx="5GNR_rx.bin"


# Transmission
./tx_samples_from_file --args serial=$serialTx --channel $channelTx --ant $antTx --gain $gain --file $fileTx --type float --freq $freq --rate $rateTx --bw $bandTx --repeat &
txpid=$!
sleep 7s
./rx_samples_to_file --args serial=$serialRx --channel $channelRx --ant $antRx --gain $gain --file $fileRx --type float --freq $freq -- rate $rateRx --bw $bandRx &
rxpid=$!
sleep 9s
kill -9 $rxpid
kill -9 $txpid
