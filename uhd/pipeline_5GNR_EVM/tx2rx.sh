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
gain=15
freq=3000000000
fileTx="5GNR_sent_withPadding.bin"
fileRx="5GNR_rx.bin"


# Transmission
./tx --args serial=$serialTx --channel $channelTx --ant $antTx --gain $gain --file $fileTx --type float --freq $freq --rate $rateTx --bw $bandTx --repeat
pid=$txpid
sleep 1s
./rx --args serial=$serialRx --channel $channelRx --ant $antRx --gain $gain --file $fileRx --type float --freq $freq -- rate $rateRx --bw $bandRx 
pid=$rxpid
sleep 0.25s
kill -9 $rxpid
kill -9 $txpid
