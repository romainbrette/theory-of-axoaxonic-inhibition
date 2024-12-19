'''
Somatic response to a hyperpolarizing current injected at AIS vs. soma.
'''
from brian2 import *
from shared.models import params_model_description, model_Na_Kv1
from shared.plotting import *

figsize = (5, 2.5)

defaultclock.dt = 0.005*ms

params = params_model_description

### Model
# Nav and Kv are turned off
morpho = params.morpho
params.gk_dend = 0* (siemens / meter ** 2)
params.gna_dend = 0* (siemens / meter ** 2)
params.gk_soma = 0* (siemens / meter ** 2)
params.gna_soma = 0* (siemens / meter ** 2)
neuron = model_Na_Kv1(params, -75.*mV, morpho=morpho,
                      gna_density = 0.* (siemens / meter ** 2), gk_density = 0.* (siemens / meter ** 2))

M = StateMonitor(neuron, ('v', 'I'), record = 0) # recording at the soma

run(50*ms)
neuron.I[0] = -100*pA
run(100*ms)
neuron.I[0] = 0*pA
run(150*ms)
neuron.I[morpho.axon[40*um]] = -100*pA
run(100*ms)
neuron.I[morpho.axon[40*um]] = 0*pA
run(100*ms)

# Plotting
fig = figure(1, figsize=figsize)
gs = GridSpec(3, 1, figure=fig)

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[1:, 0])
prepare_panel(ax1, ax2)

ax1.plot(M.t/ms, M.I[0]/pA, 'k')
ax1.set_ylabel('I (pA)')
ax1.set_xlim(0, 250)

ax2.plot(M.t/ms, M.v[0]/mV, 'k')
ax2.plot(M.t/ms -250, M.v[0]/mV, 'grey')
ax2.set_xlabel(r'Time (ms)')
ax2.set_ylabel('Vm (mV)')
ax2.set_xlim(0, 250)

tight_layout()

savefig('current_sink.png', dpi=300)

show()
