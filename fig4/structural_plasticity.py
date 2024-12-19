'''
Effect of AIS displacement on threshold.
Fixed axoaxonic synapses vs. moving along with the AIS.
'''
from __future__ import print_function
from brian2 import *
from shared.models import params_model_description, model_Na_Kv1_extended_GABA
from shared.analysis.trace_analysis import *
from shared.plotting import *
import matplotlib.patches as patches

defaultclock.dt = dt = 0.001*ms
ra = 4*params_model_description.Ri/(pi*params_model_description.axon_diam**2)

### AIS geometry and parameters
params = params_model_description
morpho = params.morpho

all_x = linspace(0, 45, 20) * um # AIS start position
g_GABA_tot = 5*nS
EGABA = -90*mV

GABA_start = 15*um
GABA_length = 30*um

g_GABA = g_GABA_tot / (pi * GABA_length * morpho.axon.diameter[0])

def calculate_threshold(x1=None, move_along=False): # x1 and x2 are the start and end of the GABA
    ### Model

    if not move_along:
        GABA_segment = morpho.axon[x1:x1+GABA_length]
    threshold = []
    # Measure threshold
    for x in all_x: # AIS position
        print(x)
        if move_along:
            x1 = x #+ 15*um
            GABA_segment = morpho.axon[x1:x1+GABA_length]
        end = x + 30. * um  # AIS end position
        neuron = model_Na_Kv1_extended_GABA(params, -75. * mV, x, end, morpho=morpho)

        M = StateMonitor(neuron, ('v'), record=0)  # recording at the soma

        run(50 * ms)
        neuron.I_CC[0] = 1. * nA
        run(5 * ms)
        neuron.I_CC[0] = 0 * amp
        run(50 * ms)
        neuron.gGABA[GABA_segment] = g_GABA
        run(10 * ms)
        neuron.I_CC[0] = 1. * nA
        run(5 * ms)
        neuron.I_CC[0] = 0 * amp
        run(10 * ms)
        v = M.v[0]
        #thresholds = v[spike_onsets(v, criterion=40*volt/second*dt, v_peak=-20*mV)]
        thresholds = v[spike_onsets_dv2(v, v_peak=-20*mV)]
        threshold.append(thresholds[-1]-thresholds[0])
    threshold = array(threshold)*volt
    return threshold

threshold1 = calculate_threshold(move_along=True)
threshold2 = calculate_threshold(GABA_start, move_along=False)

# Plotting
f1 = figure(1, figsize=(5, 2))
ax1, ax2 = f1.subplots(1, 2, sharey=True)
prepare_panel(ax1, ax2)

ax1.plot(all_x/um + 15, threshold1/mV, 'k') # horizontal: middle position of AIS
ax2.plot(all_x/um + 15, threshold2/mV, 'k') # horizontal: middle position of AIS

bottom, top = plt.gca().get_ylim()

ax1.set_ylabel(r'$\Delta V^*$ (mV)')
ax1.set_xlabel(r'AIS center (µm)')

rectangle2 = patches.Rectangle((GABA_start/um, bottom), GABA_length/um, top-bottom, color='k', alpha=0.2) # 30 um around the synapse
ax2.add_patch(rectangle2)
ax2.set_xlabel(r'AIS center (µm)')

tight_layout()

savefig('threshold_AIS.png', dpi=300)

show()
