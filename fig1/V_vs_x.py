'''
Membrane potential vs. distance from the soma, for a hyperpolarizing current at the AIS.
'''
from brian2 import *
from shared.models import params_model_description, model_Na_Kv1
from shared.plotting import *

figsize = (5,2)

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

compartments = morpho.axon[0*um:80*um]

run(100*ms)
neuron.I[morpho.axon[20*um]] = -100*pA
run(5*ms)
profile1 = neuron.v[compartments]

neuron.I[morpho.axon[20*um]] = 0*pA
run(100*ms)
neuron.I[morpho.axon[40*um]] = -100*pA
run(5*ms)
profile2 = neuron.v[compartments]

# Plotting
fig = figure(1, figsize=figsize)
ax = subplot(111)
prepare_panel(ax)

ax.plot(neuron.distance[compartments]/um, profile1/mV, 'orange')
ax.plot(neuron.distance[compartments]/um, profile2/mV, 'r')
ax.set_xticks([0,20,40,60,80])
ax.set_xlabel(r'Distance from soma ($\mu$m)')
ax.set_ylabel('Vm (mV)')

tight_layout()

savefig('V_vs_x.png', dpi=300)

show()
