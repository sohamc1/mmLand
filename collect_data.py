from AppNotes_WinLinux.Class import TinyRad
import time as time
import numpy as np
from scipy.io import savemat
from scipy.constants import speed_of_light
import pickle

DISTANCE_FROM_TAG = 5.0

try:
    # Change these two parameters for data collection
    description = "test_15cm_across_2m_away"
    duration = 10 # in seconds, for more data collected

    filePrefix = "./"
    filename = time.strftime(filePrefix + description + "-24GHz.pickle", time.localtime())
    

    period = 450
    chirpDur = 150        # Chirp duration, in microseconds (We change this between 200 and 150)
    inter_chirp_delay = period - chirpDur
    c0 = speed_of_light         # 2.998e8
    switchPer = 0.001       # Tag switching period

    modF = 1 / (2 * switchPer)  # Tag modulation frequency

    # Setup Connection
    Brd = TinyRad.TinyRad('Usb', '127.0.0.1')
    Brd.BrdRst()

    # Software Version  
    Brd.BrdDispSwVers()

    # Configure Receiver
    Brd.RfRxEna()
    TxPwr = 100

    # Configure Transmitter (Antenna 0 - 2, Pwr 0 - 100)
    Brd.RfTxEna(0, TxPwr)

    CalDat = Brd.BrdGetCalDat()['Dat']

    # Configure Measurements
    Cfg = dict()
    Cfg['fStrt'] = 24.00e9                  # Start frequency
    Cfg['fStop'] = 24.25e9                  # Stop frequency
    Cfg['TRampUp'] = chirpDur * 1e-6        # Chirp duration (in seconds)
    Cfg['Perd'] = Cfg['TRampUp'] + inter_chirp_delay*1e-6   # Period, chirp duration + between-chirp wait
    Cfg['N'] = 128                          # Samples per chirp
    Cfg['Seq'] = [1]                        # Tx sequence, unused
    Cfg['CycSiz'] = 2                       # ?
    Cfg['FrmSiz'] = 256                     # Chirps per frame, don't use
    Cfg['FrmMeasSiz'] = 256                 # Sampled chirps per frame, use this instead
    frmMeasSiz = Cfg['FrmMeasSiz']

    Brd.RfMeas(Cfg)
    time.sleep(1)

    # Read actual configuration
    N = int(Brd.Get('N'))
    NrChn = int(Brd.Get('NrChn'))
    fs = Brd.Get('fs')
    Perd = Cfg['Perd']

    # I have no idea why there are duplicate lines
    N = int(Brd.Get('N'))

    # Container for raw radar data and timestamps
    DataAll = None
    Timestamps = np.array([])

    # Measure and calculate Range Doppler Map
    # frameSize = int(duration / (frmMeasSiz * Perd))
    # Instead of using a fixed frameSize, use time.time() compared to duration
    print("[Rad24GHz] Started measuring")
    startTime = time.time()
    while True:
        if time.time() - startTime > duration:
            break

        DataFrame = Brd.BrdGetData()

        curr = time.time()
        strcurr = time.strftime("%Y-%m-%d %H:%M:%S:", time.localtime(curr)) + str(int(curr % 1 * 1000))
        Timestamps = np.append(Timestamps, [strcurr])
        reshaped = DataFrame.reshape(DataFrame.shape[0], 1, DataFrame.shape[1])
        DataAll = np.array(reshaped, dtype='float64') if DataAll is None else np.append(DataAll, reshaped, axis=1)

    print('[Rad24GHz] Finished measuring after {0} seconds'.format(time.time() - startTime))
    print('[Rad24GHz] Collected data in shape {0}'.format(DataAll.shape))

    save = {
        "Data": DataAll,
        "dtime": Timestamps,
        # "Brd": Brd, # Brd not used in processing code, but in collectData
        "Cfg": Cfg,
        "N": N,
        "NrChn": NrChn,
        "fs": fs,
        "measRounds": DataAll.shape[1],
        "CalDat": CalDat,
        "switchPer": switchPer,
        "chirpDur": chirpDur,
    }
    # print(save['Data'].dtype)
    with open(filename, 'wb') as f:
        pickle.dump(save, f)

finally:
    del Brd
    time.sleep(1)
