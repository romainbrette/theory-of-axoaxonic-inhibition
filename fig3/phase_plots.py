'''
Phase plots of APs with varying hyperpolarizing currents at the AIS.
'''
from __future__ import print_function
from brian2 import *
from shared.models import params_model_description, model_Na_Kv1
from shared.analysis.trace_analysis import *
from shared.plotting import *
from matplotlib.cm import get_cmap
cmap = get_cmap('viridis')

figsize = (5, 2) # width, height

defaultclock.dt = dt = 0.001*ms

### AIS geometry and parameters
start =5.*um # AIS start position
end = 35.*um # AIS end position
injection = .5*(start+end)
params = params_model_description

### Model
# The conductance densities at the AIS are fixed:
morpho = params.morpho
neuron = model_Na_Kv1(params, -75.*mV, start, end, morpho=morpho)

M = StateMonitor(neuron, ('v'), record = 0) # recording at the soma
M_AIS = StateMonitor(neuron, ('v'), record = neuron.morphology.axon[end-1*um]) # recording at AIS end

# Plotting
f1 = figure(1, figsize=figsize)
ax1 = subplot(111)
prepare_panel(ax1)

f2 = figure(2, figsize=figsize)
ax3 = subplot(111)
prepare_panel(ax3)

### Simulation
run(50*ms) # for convergence
store()

all_I = linspace(0,200,5)*pA
colors = cmap(np.linspace(0, 1, len(all_I)))

for i,I in enumerate(all_I):
    print(I)
    restore()
    neuron.I_CC[0] = 0 * amp
    neuron.I[0] = 0 * pamp
    neuron.I[morpho.axon[injection]] = -I
    run(5 * ms)
    neuron.I_CC[0] = 1 * nA
    run(5 * ms)
    neuron.I_CC[0] = 0 * amp
    run(10 * ms)
    v = M.v[0]
    v_AIS = M_AIS.v[0]
    onsets = spike_onsets(v, criterion=40*volt/second*dt, v_peak=-20*mV)
    print(v[onsets])

    # Plotting
    dv = diff(v)/dt
    dv_AIS = diff(v_AIS) / dt
    ax1.plot(v[:-1]/mV, dv, color = cmap(i/len(all_I)))

    ax3.plot(M.t/ms, M.v[0]/mV, color_soma_Iaxon, label='soma', color = cmap(i/len(all_I)))

ax1.set_ylabel('dV (V/s)')
ax1.set_xlabel('V (mV)')

f1.tight_layout()
f1.savefig('phase_plots.png', dpi=300)
f2.tight_layout()

show()
