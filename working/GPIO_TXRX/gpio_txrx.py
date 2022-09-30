#!/usr/bin/python3
import numpy as np
import uhd
import time

class args():
    # RX 
    freq = 2.4e9
    duration = 1e-6
    rate = 1e7
    channels = [0]
    gain = 0
    # TX
    wave_ampl = 1
    waveform  = "const"
    wave_freq = 3e9
    # GPIO
    all_one = 0xFF
    all_zero = 0x00
    gpio_line = 0xF


WAVEFORMS = {
    "sine": lambda n, tone_offset, rate: np.exp(n * 2j * np.pi * tone_offset / rate),
    "square": lambda n, tone_offset, rate: np.sign(WAVEFORMS["sine"](n, tone_offset, rate)),
    "const": lambda n, tone_offset, rate: 1 + 1j,
    "ramp": lambda n, tone_offset, rate:
            2*(n*(tone_offset/rate) - np.floor(float(0.5 + n*(tone_offset/rate))))
}

data_send = np.array(
        list(map(lambda n: args.wave_ampl * WAVEFORMS[args.waveform](n, args.wave_freq, args.rate),
                 np.arange(
                     int(10 * np.floor(args.rate / args.wave_freq)),
                     dtype=np.complex64))),
        dtype=np.complex64)  # One period


def main():
    """RX TX switch testing for GPIO .... """
    usrp = uhd.usrp.MultiUSRP("addr=192.168.30.2")
    num_samps = int(np.ceil(args.duration*args.rate))
    usrp.set_gpio_attr("FP0", "DDR", args.all_one, args.gpio_line)
    usrp.set_gpio_attr("FP0", "CTRL", args.all_zero, args.gpio_line)
    usrp.set_gpio_attr("FP0", "ATR_RX", args.all_one, args.gpio_line)
    while True:
        usrp.recv_num_samps(num_samps, args.freq, args.rate, args.channels, args.gain)
        usrp.send_waveform(data_send, args.duration, args.freq, args.rate,
                        args.channels, args.gain)
        time.sleep(1/1000000.0)
    # if not isinstance(args.channels, list):
    #     args.channels = [args.channels]
    # samps = usrp.recv_num_samps(num_samps, args.freq, args.rate, args.channels, args.gain)
    # with open(args.output_file, 'wb') as out_file:
    #     if args.numpy:
    #         np.save(out_file, samps, allow_pickle=False, fix_imports=False)
    #     else:
    #         samps.tofile(out_file)

if __name__ == "__main__":
    main()
