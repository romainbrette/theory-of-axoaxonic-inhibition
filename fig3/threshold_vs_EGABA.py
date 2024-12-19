'''
Threshold change vs. EGABA.
'''
from __future__ import print_function
from brian2 import *
from shared.models import params_model_description, model_Na_Kv1_GABA
from shared.analysis.trace_analysis import *
from shared.plotting import *
from matplotlib.cm import get_cmap
cmap = get_cmap('viridis')

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

M = StateMonitor(neuron, ('v'), record = 0) # recording at the soma
EGABA = -70*mV
run(50*ms)
store()
neuron.I_CC[0] = 1 * nA
run(5 * ms)
neuron.I_CC[0] = 0 * amp
run(10 * ms)
v = M.v[0]
#thresholds = v[spike_onsets(v, criterion=40 * volt / second * dt, v_peak=-20 * mV)]
thresholds = v[spike_onsets_dv2(v, v_peak=-20 * mV)]
threshold0 = thresholds[0]

#f1 = figure(1, figsize=(3, 3))

threshold = []
threshold_soma = []
all_EGABA = linspace(-90,-20, 20)*mV
# Measure threshold
for i, EGABA in enumerate(all_EGABA):
    print(EGABA)
    ### Simulation
    restore()
    neuron.gGABA[morpho.axon[middle]] = 5*nsiemens
    run(5*ms)
    neuron.I_CC[0] = 1*nA
    run(5*ms)
    neuron.I_CC[0] = 0*amp
    run(10*ms)

    v = M.v[0]
    #thresholds = v[spike_onsets(v, criterion=40 * volt / second * dt, v_peak=-20 * mV)]
    thresholds = v[spike_onsets_dv2(v, v_peak=-20 * mV)]
    threshold.append(thresholds[0])
    #plot(M.t / ms, M.v[0]/mV, color = cmap(i/len(all_EGABA)))
threshold = array(threshold) * volt

# Plotting
f2 = figure(2, figsize=(5, 2))
ax1 = subplot(121)
prepare_panel(ax1)
ax1.plot(all_EGABA / mV, threshold / mV, 'k')
ax1.plot(all_EGABA / mV, (threshold0+all_EGABA*0)/mV, 'k--')
ax1.set_ylabel(r'Threshold (mV)')
ax1.set_xlabel(r'$E_{GABA}$ (mV)')
#ax1.set_yticks([-58, -56, -54])

tight_layout()

savefig('threshold_vs_EGABA.png', dpi=300)

show()
