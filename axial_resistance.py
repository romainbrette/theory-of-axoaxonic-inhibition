'''
Calculation of axial resistance based on simultaneous soma-bleb patch recordings.
Data from Hu and Bean (2018): https://zenodo.org/record/3539297
'''
from shared.data.hu_bean_shared import *
import numpy as np

distances, Rap, Ran, Rinput, risetime = [], [], [], [], []

dt = 0.02 # in ms

n_before = int(5/dt) # number of time steps before the pulse (5 ms)

for distance,t,Vs,Va,Is,Ia in load_all_files():
    # Find places where there's an axonal but not somatic current pulse
    i = 1*((abs(Ia)>0.008) & (abs(Is)<0.008))
    start = np.where(np.diff(i) == 1)[0] +1
    end = np.where(np.diff(i) == -1)[0]

    # Make a list of traces
    vs = []
    va = []
    current = []
    for i,j in zip(start,end):
        if j-i>int(480/dt): # long pulse (>480 ms)
            current.append(np.median(Ia[i:j])) # average current
            vs.append(Vs[i-n_before:j])
            va.append(Va[i-n_before:j])

    current = np.array(current)

    try:
        # Select the smallest positive and negative currents
        ind = np.where(current>0)[0]

        i=ind[np.argmin(current[ind])]
        current_pos = current[i]
        vs_pos = vs[i]
        va_pos = va[i]

        ind = np.where(current < 0)[0]
        i=ind[np.argmin(-current[ind])]
        current_neg = current[i]
        vs_neg = vs[i]
        va_neg = va[i]

    except ValueError:
        print("Distance ",distance,"um", "No current pulse found")
        continue

    n = len(vs_pos)

    # Resting potential calculated just before the pulse
    vs_rest = np.median(vs_pos[:n_before])
    va_rest = np.median(va_pos[:n_before])
    # Resistances (at axon, at soma, and axial)
    R_input = (np.median(va_pos[n//2:])-va_rest)/current_pos # input resistance at AIS
    Rm_long = (np.median(vs_pos[n//2:])-vs_rest)/current_pos # long term somatic resistance for AIS injection
    R_axial = R_input-Rm_long
    #print(R_axial, R_input)
    distances.append(distance)
    Rap.append(R_axial)
    Rinput.append(R_input)

    # Rise time (0-90%)
    t75 = np.where(va_pos[n_before:]-va_rest-vs_pos[n_before:]+vs_rest>0.75*R_axial*current_pos)[0][0]
    #print("Rise time = ",t75*dt)
    risetime.append(t75*dt)

    ## Checking with the negative pulse (very similar)
    # Rest
    vs_rest = np.median(vs_neg[:n_before])
    va_rest = np.median(va_neg[:n_before])
    # Resistances
    R_input = np.median(va_neg[n//2:]-va_rest)/current_neg # input resistance at AIS
    Rm_long = np.median(vs_neg[n//2:]-vs_rest)/current_neg # long term somatic resistance for AIS injection
    R_axial = R_input-Rm_long
    Ran.append(R_axial)

distances = np.array(distances)

# Print sorted list
for d,Rp, Rn in sorted(zip(distances,Rap,Ran)):
    print('{} um: positive pulse, {} MOhm; negative pulse, {} MOhm'.format(d,Rp, Rn))

def resistance_stats(R):
    R = np.array(R)
    slope = np.mean(R/distances)
    sd = np.std(R/distances)
    n = len(R)
    print('Axial resistance per unit length = {} +- {} MOhm/um (n={})'.format(slope, sd, n))

print('Negative pulses:')
resistance_stats(Ran)
