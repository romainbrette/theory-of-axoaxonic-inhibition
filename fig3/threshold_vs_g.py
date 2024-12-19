'''
Threshold change vs. synaptic conductance g.
'''
from __future__ import print_function
from brian2 import *
from shared.models import params_model_description, model_Na_Kv1_GABA
from shared.analysis.trace_analysis import *
from shared.plotting import *

defaultclock.dt = dt = 0.001*ms

### AIS geometry and parameters
start =5.*um # AIS start position
end = 35.*um # AIS end position
middle = .5*(start+end)
params = params_model_description

### Model
# The conductance densities at the AIS are fixed:
morpho = params.morpho
neuron = model_Na_Kv1_GABA(params, -75.*mV, start, end, morpho=morpho)
EGABA = -70*mV

M = StateMonitor(neuron, ('v'), record = 0) # recording at the soma
store()

all_g = linspace(0,20,20)*nsiemens
threshold = []
# Measure threshold
for g in all_g:
    print(g)

    ### Simulation
    restore()
    neuron.I_CC[0] = 0*amp
    neuron.gGABA[morpho.axon[middle]] = g
    run(5*ms)
    neuron.I_CC[0] = 1*nA
    run(5*ms)
    neuron.I_CC[0] = 0*amp
    run(10*ms)

    v = M.v[0]
    #thresholds = v[spike_onsets(v, criterion=40 * volt / second * dt, v_peak=-20 * mV)]
    thresholds = v[spike_onsets_dv2(v, v_peak=-20 * mV)]
    threshold.append(thresholds[0])
threshold = array(threshold) * volt


## Plotting
f1 = figure(1, figsize=(5, 2))
ax1 = subplot(122)
prepare_panel(ax1)

ax1.plot(all_g / nS, threshold / mV, 'k')
# Prediction
Ra = middle*4*params_model_description.Ri/(pi*params_model_description.axon_diam**2)
print("Ra =", Ra/Mohm, "MOhm")
prediction = threshold[0]+all_g*Ra*(threshold[0]-EGABA)
ax1.plot(all_g / nS, prediction / mV, '--k')

#ax1.set_ylabel(r'Threshold (mV)')
ax1.set_xlabel(r'$g$ (nS)')

tight_layout()

savefig('threshold_vs_g.png', dpi=300)

show()
