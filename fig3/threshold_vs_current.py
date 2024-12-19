'''
Threshold vs. axosomatic voltage gradient and injected current at the AIS.
'''
from __future__ import print_function
from brian2 import *
from shared.models import params_model_description, model_Na_Kv1
from shared.analysis.trace_analysis import *
from shared.plotting import *

figsize = (5, 2) # width, height

defaultclock.dt = dt = 0.001*ms

### AIS geometry and parameters
start = 5.*um # AIS start position
end = 35.*um # AIS end position
injection = .5*(start+end)
params = params_model_description

### Model
# The conductance densities at the AIS are fixed:
morpho = params.morpho
neuron = model_Na_Kv1(params, -75.*mV, start, end, morpho=morpho)

M = StateMonitor(neuron, ('v'), record = 0) # recording at the soma
run(50*ms)
store()

threshold = []
gradient = []
all_I = linspace(0,200,10)*pA
# Measure threshold
for I in all_I:
    print(I)
    restore()
    neuron.I[morpho.axon[injection]] = -I
    run(5 * ms)
    gradient.append(neuron.v[0] - mean(neuron.v[morpho.axon[start:end]]))
    neuron.I_CC[0] = 1 * nA
    run(5 * ms)
    neuron.I_CC[0] = 0 * amp
    run(10 * ms)
    v = M.v[0]
    #thresholds = v[spike_onsets(v, criterion=40*volt/second*dt, v_peak=-20*mV)]
    thresholds = v[spike_onsets_dv2(v, v_peak=-20*mV)]
    print(v[spike_onsets(v, criterion=40*volt/second*dt, v_peak=-20*mV)][0], v[spike_onsets_dv2(v, v_peak=-20*mV)][0])
    threshold.append(thresholds[0])
threshold = array(threshold)*volt
gradient = array(gradient)*volt

# Plotting
f1 = figure(1, figsize=figsize)
ax1, ax2 = f1.subplots(1, 2)
prepare_panel(ax1)
prepare_panel(ax2)

ax1.plot(all_I/pA, threshold/mV, 'k')
ax1.set_ylabel(r'Threshold (mV)')
ax1.set_xlabel('I (pA)')
ax1.set_yticks([-58, -56, -54])

tight_layout()

print('Threshold :', threshold[0]/mV, 'mV')
ax2.plot(gradient/mV, threshold/mV, 'k')
ax2.plot(gradient/mV, (threshold[0]+gradient-gradient[0])/mV, 'k--')
#ax2.set_ylabel(r'Threshold (mV)')
ax2.set_xlabel('Gradient (mV)')
ax2.set_yticks([-58, -56, -54])
#ax2.set_yticks([-56, -54, -52])
#ax2.set_yticklabels([])

tight_layout()

savefig('threshold_current.png', dpi=300)

show()
