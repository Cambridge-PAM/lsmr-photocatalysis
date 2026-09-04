from lsmr_photocatalysis import ContinuousFlow, FLOW_PARAMS
import numpy as np
import serial
import time
from datetime import datetime, timedelta

"""
Automated pump control script to initiate flow rate changes according to the continuous flow scheme 

tbuffer:    time in s after the solution fully fills the reactor chamber before recording
vrxn:       reaction volume under LED/film in ul
vfull:      full reactor chamber volume in ul
"""

tbuffer, vrxn, vfull = FLOW_PARAMS['250ul'] # change this if setup is different

tRs = np.arange(60, 7*60+1, 60) # residence times in s for 1-7 min in 1 min intervals
flowrates = vrxn/(tRs/60) # flow rates in ul/min for each tR
points = ContinuousFlow.flowpoints(tRs, tbuffer, vrxn, vfull) # time points in s to record absorbances

print()
print(f'tRs (min): {tRs/60}')
print(f'Flowrates (ul/min): {flowrates}')
print(f'Flowrate switch time intervals (s): {points}')
print()

com = 'COM4' # change this if necessary, check device manager

# pump controls using rs-232 serial commands
def rs232():
    ser = serial.Serial(com, 19200, parity=serial.PARITY_NONE, bytesize=serial.EIGHTBITS, stopbits=serial.STOPBITS_ONE, 
                        timeout=0.05, xonxoff=0, rtscts=0)
    print('Pump connected')

    try:
        def send(cmd):
            ser.write((cmd + '\r').encode()) # send commands as encoded strings
            time.sleep(0.05)
            response = ser.readline().decode().strip()
            if response:
                print(f'{cmd} -> {response}')
            else:
                print(f'{cmd} -> (no response)')

        send('RESET')
        send('DIA 9.000')
        send('DIR INF')
        send('PHN 1')
        send('FUN RAT')

        for i in range(len(tRs)):
            print()
            if i == 0:
                send(f'RAT {flowrates[i]} UM')
                print()
                print("Remember to have already started OceanView data collection beforehand!")
                input("Press ENTER and turn LED on simultaneously to start...")
                send('RUN')
                start = time.perf_counter()
            else:
                send(f'RAT {flowrates[i]}')
            print()
            print(f'New flow rate = {flowrates[i]:.2f} ul/min for tR = {tRs[i]/60:.1f} min (Run {i+1}/{len(tRs)})')
            print(f'Elapsed time: {time.perf_counter()-start:.3f} s, Timestamp: {datetime.now().time()}')
            print(f'Next flow rate switch at: {(datetime.now()+timedelta(seconds=points[i])).time()}')
            time.sleep(points[i]-0.1)

        print("Experiment finished!")
        print(f'Elapsed time: {time.perf_counter()-start:.3f} s')
        
    except Exception as e:
        print()
        print(f'Pump error: {e}')

    finally:
        print()
        send('STP')
        print(f'Pump stopped')
        print(f'Elapsed time: {time.perf_counter()-start:.3f} s, Timestamp: {datetime.now().time()}')
        print()
        ser.close()

rs232()