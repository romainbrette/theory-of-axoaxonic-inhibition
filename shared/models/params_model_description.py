"""
Parameters of the biophysical model.
"""
from brian2 import *

### Morphology
morpho = Soma(30. * um)
dend_diam = 6.*um 
dend_length = 1000.*um 
axon_diam = 1. * um
axon_length = 500. * um
dendrite = Cylinder(diameter=dend_diam, length=dend_length, n=500)
axon = Cylinder(diameter=axon_diam, length=axon_length, n=500)
morpho.dendrite = dendrite
morpho.axon = axon

# Passive parameters
EL = -75. * mV
Cm = 0.9* uF / cm ** 2
gL = 1./(15000 * ohm * cm**2) 
Ri = 100. * ohm * cm

# Na channels parameters
ENa = 70. * mV

# K channels parameters
EK = -90. * mV

# Channels kinetics
T = 33.
factor = (1. / 2.8) ** ((T - 23.) / 10.)

# Na+
Va = -35. * mV  # Kole: -31 mV
Ka = 5. * mV  
Taum_max = factor * 0.15 * ms
Vh = -60.*mV # value from Kole et al. 2008
Kh = 5. * mV  
Tauh_max = factor * 5.* ms  

# K+:
Vn = -70.*mV 
Kn = 20.*mV 
Taun_max = 1. * ms 

### Soma

# Na channels parameters
gna_soma = 250. * (siemens / meter ** 2) 
gna_dend = 50. * (siemens / meter ** 2) 

# K channels parameters
gk_soma = 250. * (siemens / meter ** 2) 
gk_dend = 50. * (siemens / meter ** 2) 

## Channels kinetics
# Na+:
Va_soma = -30.*mV
Vh_soma = -55.*mV # # value from Kole et al. 2008
