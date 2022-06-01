#!/bin/bash

# Parameter Setup
downRate=1
# -1 for High [MHz]
setup=-1

filePathTx=./Tx/
if [[ $setup -gt 0 ]]
then
	filePathRx=./Cable_$setup/
else
	filePathRx=./Cable/
fi
mkdir $filePathRx

serialRx="321FCC3"
serialTx="322A004"
channelRx=0
channelTx=0
antRx="RX2"
antTx="TX/RX"
gain=10
freq=1000000000

# Transmission
for fileTx in $filePathTx* 
do
	if [ -d $fileTx ];
	then
		continue
	fi
	
	file=$(basename $fileTx)
	fileRx=$filePathRx$file
	if [[ $file == *"CBW5"* ]]
	then
		rateTx=$((5000000/downRate))
	elif [[ $file == *"CBW10"* ]]
	then
		rateTx=$((10000000/downRate))
	elif [[ $file == *"CBW20"* ]]
	then
		rateTx=$((20000000/downRate))
	elif [[ $file == *"CBW40"* ]]
	then
		rateTx=$((40000000/downRate))
	elif [[ $file == *"CBW80"* ]]
	then
		rateTx=$((80000000/downRate))
	elif [[ $file == *"CBW160"* ]]
	then
		rateTx=$((160000000/downRate))
	fi
	
	if [[ $setup -gt 0 ]]
	then
		rateRx=$(($setup*1000000/downRate))
	else
		rateRx=$rateTx
	fi
	bandTx=$rateTx
	bandRx=$rateTx
	
	echo $file $rateTx

	./rx --args serial=$serialRx --channel $channelRx --ant $antRx --gain $gain --file $fileRx --type float --freq $freq -- rate $rateRx --bw $bandRx &
	pid=$!
	sleep 3s
	./tx --args serial=$serialTx --channel $channelTx --ant $antTx --gain $gain --file $fileTx --type float --freq $freq --rate $rateTx --bw $bandTx
	sleep 0.25s
	kill -9 $!
	sleep 0.25s
	python PacketDetection_3.py --fileTx $fileTx --rateTx $rateTx --fileRx $fileRx --rateRx $rateRx --down $downRate --out $fileRx # $filePathRx"_"$file # $fileRx
	sleep 0.5s
done
