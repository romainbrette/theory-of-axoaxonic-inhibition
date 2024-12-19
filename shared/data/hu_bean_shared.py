'''
Shared analysis functions for Hu and Bean's data (2018)
https://zenodo.org/record/3539297
'''
from __future__ import print_function
import os
from shared.data.abf import ABF
from numpy import zeros,array,median
import re

__all__ = ['load_all_files','load_file']

def load_file(filename):
    '''
    Returns distance, time, somatic voltage, axonal voltage, somatic current, axonal current
    Voltage in mV, current in nA
    '''
    # Distance
    m = re.search('\s(\d+\.?\d*)\s*um', filename)
    try:
        distance = float(m.group(1))
    except:
        print("problem with", filename)

    abf = ABF(filename)
    abf.setSweep(0,channel=0)
    Vs = abf.dataY
    abf.setSweep(0,channel=1)
    Is = abf.dataY * 0.020 # 1 = 20 pA
    abf.setSweep(0,channel=2)
    Va = abf.dataY
    abf.setSweep(0,channel=3)
    Ia = abf.dataY * 0.020
    t = abf.dataX
    if filename == '140107 1-1 48um.abf': # axon and soma were reversed
        Vs,Va = Va,Vs
        Is,Ia = Ia,Is

    # Find the current zero and fix it
    Is0 = median(Is)
    Ia0 = median(Ia)
    Is=Is-Is0
    Ia=Ia-Ia0

    return distance,t,Vs,Va,Is,Ia

def load_all_files(path = "shared/data/"):
    '''
    Returns the data of each file one by one,
    using generator syntax. To be used in a for loop.

    The output is a list of distance, t, Vs, Va, Is, Ia
    '''
    files = [f for f in os.listdir(path)]

    for i,filename in enumerate(files):
        if filename[-4:] == '.abf':
            #print("Loading",filename)
            try:
                yield load_file(path + '/' + filename)
            except:
                print("Could not load",filename)

if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import os

    path = os.path.expanduser('~/Downloads/')

    distance, t, Vs, Va, Is, Ia = load_file(path+'100420 1-2 22um.abf')
    print("dt=",t[1]-t[0])

    fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, ncols=1, sharex=True)

    ax1.plot(t,Vs,'r')
    ax1.plot(t,Va,'b')
    ax2.plot(t,Va-Vs,'k')
    ax3.plot(t,Is,'r')
    ax3.plot(t,Ia,'b')
    plt.show()
