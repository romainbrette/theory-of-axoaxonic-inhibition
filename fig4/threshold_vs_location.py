'''
Threshold change vs. synapse position for two AIS positions
'''
from __future__ import print_function
from brian2 import *
from shared.models import params_model_description, model_Na_Kv1
from shared.analysis.trace_analysis import *
from shared.plotting import *
import matplotlib.patches as patches

defaultclock.dt = dt = 0.001*ms
ra = 4*params_model_description.Ri/(pi*params_model_description.axon_diam**2)

### AIS geometry and parameters
params = params_model_description
start1, start2 = 5*um, 25*um # start position of AIS

all_x = linspace(0, 70, 15) * um
I0 = -100*pA

def calculate_threshold(start):
    ### Model
    # The conductance densities at the AIS are fixed:
    morpho = params.morpho
    end = start+30.*um # AIS end position
    neuron = model_Na_Kv1(params, -75.*mV, start, end, morpho=morpho)

    M = StateMonitor(neuron, ('v'), record = 0) # recording at the soma
    run(20*ms)
    store()

    threshold = []
    # Measure threshold
    for x in all_x:
        print(x)
        restore()
        neuron.I[morpho.axon[x]] = I0
        run(5 * ms)
        neuron.I_CC[0] = 1 * nA
        run(5 * ms)
        neuron.I_CC[0] = 0 * amp
        run(10 * ms)
        v = M.v[0]
        #thresholds = v[spike_onsets(v, criterion=40*volt/second*dt, v_peak=-20*mV)] # a common way to measure threshold
        thresholds = v[spike_onsets_dv2(v, v_peak=-20*mV)] # this one is more robust
        threshold.append(thresholds[0])
    threshold = array(threshold)*volt
    return threshold

threshold1 = calculate_threshold(start1)
threshold2 = calculate_threshold(start2)

# Plotting
f1 = figure(1, figsize=(5, 2))
ax1, ax2 = f1.subplots(1, 2, sharey=True)
prepare_panel(ax1, ax2)

ax1.plot(all_x/um, threshold1/mV, 'k')
# ohmic prediction
ax1.plot(all_x/um, (threshold1[0]-all_x*ra*I0)/mV, '--k')
ax2.plot(all_x/um, threshold2/mV, 'k')
# ohmic prediction
ax2.plot(all_x/um, (threshold2[0]-all_x*ra*I0)/mV, '--k')

ax1.set_ylim(top=-55)
ax1.set_ylabel(r'Threshold (mV)')
ax1.set_xlabel(r'Synapse position (µm)')
ax2.set_xlabel(r'Synapse position (µm)')

# Draw AIS
bottom, _ = plt.gca().get_ylim()
rectangle = patches.Rectangle((start1/um, bottom), 30, .3, color='k', alpha=0.5)
rectangle2 = patches.Rectangle((start2/um, bottom), 30, .3, color='k', alpha=0.5)
ax1.add_patch(rectangle)
ax2.add_patch(rectangle2)

tight_layout()

savefig('threshold_locations.png', dpi=300)

show()
