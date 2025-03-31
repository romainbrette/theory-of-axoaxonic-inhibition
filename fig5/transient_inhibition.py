'''
Time-varying threshold modulation by transient axo-axonic inhibition.
'''
from __future__ import print_function
from brian2 import *
from shared.models import params_model_description, model_Na_Kv1_decaying_current
from shared.analysis.trace_analysis import *
from shared.plotting import *

figsize = (5, 6) # width, height

defaultclock.dt = dt = 0.001*ms

### AIS geometry and parameters
start =5.*um # AIS start position
end = 35.*um # AIS end position
middle = .5*(start+end)
params = params_model_description
I0 = -50*pA

### Model
# The conductance densities at the AIS are fixed:
morpho = params.morpho
neuron = model_Na_Kv1_decaying_current(params, -75.*mV, start, end, morpho=morpho)

M = StateMonitor(neuron, ('v'), record = 0) # recording at the soma
M_AIS = StateMonitor(neuron, ('v', 'I'), record = neuron.morphology.axon[middle]) # recording at AIS end

# Plotting
f1 = figure(1, figsize=figsize)
ax1a = subplot(411)
ax1b = subplot(412)
ax1c = subplot(413)
ax1d = subplot(414)
prepare_panel(ax1a, ax1b, ax1c, ax1d)

### Simulation
run(50*ms) # for convergence
store()

### Plotting the effect on soma and AIS
run(5*ms)
neuron.I[morpho.axon[middle]] = I0
run(50*ms)
ax1a.plot(M.t/ms - 50, M_AIS.I[0]/pA, 'k')
ax1b.plot(M.t/ms - 50, M.v[0]/mV, 'k')
ax1c.plot(M.t/ms - 50, (M_AIS.v[0]-M.v[0])/mV, 'k')
ax1a.set_xlim(0,55)
ax1b.set_xlim(0,55)
ax1c.set_xlim(0,55)
ax1a.set_xticklabels([])
ax1b.set_xticklabels([])
ax1c.set_xticklabels([])
ax1a.set_ylabel('I (pA)')
ax1b.set_ylabel('V (mV)')
ax1c.set_ylabel(r'$\Delta$V (mV)')

### Measure threshold change

all_t = linspace(-5,50,50)*ms
threshold = []
spiketime = []
for t_stim in all_t:
    print(t_stim)

    ### Simulation
    restore()
    if (t_stim<0*second):
        run(5*ms+t_stim)
        neuron.I_CC[0] = .7*nA
        run(-t_stim)
        neuron.I[morpho.axon[middle]] = I0
        run(50*ms)
    else:
        run(5*ms)
        neuron.I[morpho.axon[middle]] = I0
        run(t_stim)
        neuron.I_CC[0] = .7*nA
        run(50*ms-t_stim)

    v = M.v[0]
    try:
        first_spike = spike_onsets_dv2(v, v_peak=-20 * mV)[0]
        threshold.append(v[first_spike])
        spiketime.append(first_spike)
        print(first_spike*dt, v[first_spike])
    except:
        print("no spike")
threshold = array(threshold) * volt
spiketime = array(spiketime) * dt

ax1d.plot(spiketime/ms - 50, threshold/mV, 'k')
ax1d.set_xlim(0,55)
ax1d.set_xlabel('Time (ms)')
ax1d.set_ylabel('Threshold (mV)')
f1.tight_layout()
f1.savefig('transient_inhibition.png', dpi=300)
show()
