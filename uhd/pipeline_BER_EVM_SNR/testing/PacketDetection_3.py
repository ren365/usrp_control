import argparse
import numpy as np
import math
import os
import time
import matplotlib.pyplot as plt
from scipy.io import savemat
from tqdm import tqdm



parser = argparse.ArgumentParser(description='Packet Detection')
parser.add_argument('--fileTx', default='Tx/input.bin', type=str,
    help='Input Tx file')
parser.add_argument('--rateTx', default=20000000, type=int,
    help='Tx Sampling Rate')
parser.add_argument('--fileRx', default='Tx/output.bin', type=str,
    help='Input Rx file')
parser.add_argument('--rateRx', default=5000000, type=int,
    help='Rx Sampling Rate')
parser.add_argument('--down', default=1, type=int,
    help='DownSampling Rate')
parser.add_argument('--out', default='Tx/_output.bin', type=str,
    help='Output Rx file')
parser.add_argument('--threshold', default=0.5, type=float,
    help='Packet Detection Threshold')
parser.add_argument('--field', default="both", type=str, choices=["LLTF", "LSTF", "both"],
    help='Packet Detection Field')



def Correlation(wave, winLen, intLen, deltaLen, mask=None, showProgress=False):
    wave_1 = wave[:-deltaLen]
    wave_2 = wave[deltaLen:]
    waveLen = np.shape(wave_1)[0]

    offsetList = []
    corrList = []
    offset = 0

    offsetBar = tqdm(range(0, waveLen-winLen, intLen), desc="Coarse-Grained Detection...") if showProgress \
        else range(0, waveLen-winLen, intLen)
    for offset in offsetBar:
        offsetList.append(offset)
        waveNow_1 = wave_1[offset: offset+winLen]
        waveNow_2 = wave_2[offset: offset+winLen]
        if mask is not None:
            waveNow_1 = waveNow_1[mask]
            waveNow_2 = waveNow_2[mask]
        waveNow_1 = waveNow_1 - np.mean(waveNow_1)
        waveNow_2 = waveNow_2 - np.mean(waveNow_2)

        C = np.abs(np.dot(waveNow_1, np.conj(waveNow_2)))
        # P = np.sum((np.abs(waveNow_1) ** 2 + np.abs(waveNow_2) ** 2) / 2)
        P = np.sqrt(np.sum(np.abs(waveNow_1) ** 2)) * np.sqrt(np.sum(np.abs(waveNow_2) ** 2))
        corr = C / (P +1e-30)
        corrList.append(corr)

    return offsetList, corrList



def PacketDetect(wave, bandTx, bandRx, threshold=0.9, fieldName="LLTF"):
    if bandTx == 5000000:
        symTime = 16e-6
    elif bandTx == 10000000:
        symTime = 8e-6
    elif bandTx in [20000000, 40000000, 80000000, 160000000]:
        symTime = 4e-6
    else:
        print("Warning: Transmitting Band Not Found!")
    symLen = int(symTime * bandRx)
    fftLen = int(0.8 * symLen)

    if fieldName == "LSTF":
        fieldOffset = 0
        fieldLen = symLen * 2
        deltaLen = int(0.2 * symLen)
        mask = np.arange(fieldLen-deltaLen)

        winLen = fieldLen
        intLen = fieldLen // 2
    elif fieldName == "LLTF":
        fieldOffset = symLen * 2
        fieldLen = symLen * 2
        deltaLen = fftLen
        mask = np.arange(fieldLen-deltaLen)

        winLen = fftLen
        intLen = fftLen // 2
    elif fieldName == "both":
        fieldOffset = 0
        fieldLen = symLen * 4
        deltaLen = fftLen
        mask = np.concatenate((np.arange(0, 2*symLen-fftLen), np.arange(2*symLen, 4*symLen-fftLen)))
        
        winLen = fieldLen - fftLen
        intLen = winLen // 2
    else:
        print("Warning: Field Not Found!")

    offsetCoarseList, corrCoarseList = Correlation(wave, winLen=winLen, intLen=intLen, deltaLen=deltaLen, showProgress=True)
    packetList = [i for i in range(len(corrCoarseList)) if corrCoarseList[i]>threshold]

    offsetFineList = []
    corrFineList = []
    for packet in tqdm(packetList, desc="Fine-Grained Detection..."):
        offset = offsetCoarseList[packet]
        waveNow = wave[offset-winLen: offset+winLen+fieldLen]
        offsetNowList, corrNowList = Correlation(waveNow, winLen=fieldLen-deltaLen, intLen=1, deltaLen=deltaLen, mask=mask, showProgress=False)

        corrMax = max(corrNowList)
        indexMax = corrNowList.index(corrMax)
        offsetMax = offsetNowList[indexMax] + offset - winLen -fieldOffset

        if len(offsetFineList)==0 or offsetMax!=offsetFineList[-1]:
            offsetFineList.append(offsetMax)
            corrFineList.append(corrMax)
    
    return offsetFineList, corrCoarseList



if __name__=='__main__':
    opt = parser.parse_args()
    FILE_TX = opt.fileTx
    RATE_TX = opt.rateTx
    FILE_RX = opt.fileRx
    RATE_RX = opt.rateRx
    DOWN = opt.down
    OUT = opt.out
    THRESHOLD = opt.threshold
    FIELD = opt.field

    waveTx = np.fromfile(open(FILE_TX, 'r'), dtype=np.complex64)
    waveLenTx = np.shape(waveTx)[0]
    sampleRateTx = RATE_TX
    waveRx = np.fromfile(open(FILE_RX, 'r'), dtype=np.complex64)
    sampleRateRx = RATE_RX
    waveLenRx = np.shape(waveRx)[0]

    startTime = time.time()
    offsetList, corrList = PacketDetect(waveRx, bandTx=sampleRateTx, bandRx=sampleRateRx, threshold=THRESHOLD, fieldName=FIELD)
    endTime = time.time()
    # print(endTime - startTime)
    packet = offsetList[-1]

    plt.clf()
    plt.plot(np.arange(len(corrList)), corrList, linewidth=0.1)
    plt.savefig(FILE_RX.replace(".bin", "_1.png"), dpi=200)

    plt.clf()
    plt.plot(
        np.arange(int(0.2*sampleRateRx)), 
        np.log(np.abs(waveRx[packet-int(0.1*sampleRateRx): packet+int(0.1*sampleRateRx)])+1e-5),
        linewidth=0.1)
    plt.plot([int(0.1*sampleRateRx), int(0.1*sampleRateRx)], [-10, 0], linewidth=0.2, color='r', linestyle='dashed')
    plt.savefig(FILE_RX.replace(".bin", "_2.png"), dpi=200)

    waveLen = waveLenTx * sampleRateRx // sampleRateTx
    if len(offsetList) == 0:
        print(FILE_TX + ": Packet NOT Found!")
        os.remove(FILE_TX)
        os.remove(FILE_RX)
        exit()
    # print(FILE_TX + ": Packet Found at " + str(offsetList))
    offset = packet
    wave = waveRx[offset: offset+waveLen]
    wave.tofile(open(OUT, 'w'))

    dict = {"packet": offsetList, "corr": corrList}
    savemat(FILE_RX.replace(".bin", ".mat"), dict)
