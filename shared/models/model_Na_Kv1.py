'''
Biophysical model of an AP.

This is a model with a spherical isopotential soma, a large dendrite and an unmyelinated axon. 
The AIS is located in the proximal axon and has a spatial extent.

The AIS contains a higher density of Na and K channels than the soma, the dendrite and the rest of the axon.
'''
from brian2 import *

__all__ = ['model_Na_Kv1']

def model_Na_Kv1(params,  resting_vm = -75.*mV, Na_start = 5.*um, Na_end = 35.*um,
                 gL_soma = 1./(15000. * ohm * cm**2) ,
                 gna_density = 4000.* (siemens / meter ** 2), gk_density = 1500.* (siemens / meter ** 2), morpho = None):
    '''
    params: passive and channels parameters,
    neuron.I_CC is a current-clamp stimulation.
    '''

    ### Passive parameters
    EL = resting_vm 
    Cm = params.Cm 
    gL = params.gL 
    Ri = params.Ri 

    ### AIS

    # Na channels parameters
    ENa = params.ENa 

    # K channels parameters
    EK = params.EK 

    Va = params.Va 
    Ka = params.Ka 
    Taum_max = params.Taum_max 
    Vh = params.Vh 
    Kh = params.Kh 
    Tauh_max = params.Tauh_max

    # K+:
    Vn = params.Vn 
    Kn = params.Kn 
    Taun_max = params.Taun_max 

    ### Soma

    # Na channels parameters
    gna_soma = params.gna_soma 
    gna_dend = params.gna_dend 

    # K channels parameters
    gk_soma = params.gk_soma 
    gk_dend = params.gk_dend 

    ## Channels kinetics
    # Na+:
    Va_soma = params.Va_soma 
    Vh_soma = params.Vh_soma 

    # Equations
    eqs = '''
    Im = (gL*(EL-v) + gNa*m*h*(ENa-v) + gK*n**8*(EK-v) ) : amp/meter**2
    INa = gNa*m*h*(ENa-v) : amp/meter**2
    IK = gK*n**8*(EK-v) : amp/meter**2

    dm/dt = alpham*(1-m) - betam*m : 1
    dh/dt = alphah*(1-h) - betah*h : 1
    dn/dt = alphan*(1-n) - betan*n : 1

    alpham = (1/ka)*(v-va) / (1-exp(-(v-va)/ka)) /(2*taum_max) : Hz
    betam = -(1/ka)*(v-va) / (1-exp((v-va)/ka)) /(2*taum_max) : Hz

    alphah = -(1/kh)*(v-vh) / (1-exp((v-vh)/kh)) /(2*tauh_max) : Hz
    betah = (1/kh)*(v-vh) / (1-exp(-(v-vh)/kh)) /(2*tauh_max) : Hz

    alphan = (1/kn)*(v-vn) / (1-exp(-(v-vn)/kn)) /(2*taun_max) : Hz
    betan = -(1/kn)*(v-vn) / (1-exp((v-vn)/kn)) /(2*taun_max): Hz
    
    gL: siemens/meter**2
    gNa : siemens/meter**2
    gK : siemens/meter**2
    va : volt
    vh : volt
    vn : volt
    ka : volt
    kh : volt
    kn : volt
    taum_max : second
    tauh_max : second
    taun_max : second

    I : amp (point current)
    I_CC : amp (point current) 
    '''
    
    neuron = SpatialNeuron(morphology=morpho, model=eqs, Cm=Cm, Ri=Ri,
                                    namespace=dict(EL=EL, ENa=ENa, EK=EK,
                                                   Ka=Ka, Va=Va, Taum_max=Taum_max,
                                                   Kh=Kh, Vh=Vh, Tauh_max=Tauh_max,
                                                   Kn=Kn, Vn=Vn, Taun_max=Taun_max,
                                                   Va_soma=Va_soma, 
                                                   Vh_soma=Vh_soma),
                                                   method="exponential_euler")

    # Parameters of the soma are assigned in the entire neuron for integration purpose, 
    # but gNa and gK at the soma and AIS are then modified
    neuron.va = Va_soma
    neuron.ka = Ka
    neuron.taum_max = Taum_max

    neuron.vh = Vh_soma
    neuron.kh = Kh
    neuron.tauh_max = Tauh_max

    neuron.vn = Vn
    neuron.kn = Kn
    neuron.taun_max = Taun_max

    #Dendrites and axon
    neuron.gNa = gna_dend
    neuron.gK = gk_dend
    neuron.gL = gL

    # Soma
    neuron.gNa[0] = gna_soma
    neuron.gK[0] = gk_soma
    neuron.gL[0] = gL_soma
    
    # Initial segment
    initial_segment = morpho.axon[Na_start:Na_end]
    neuron.gNa[initial_segment] = gna_density
    neuron.gK[initial_segment] = gk_density
    
    neuron.va[initial_segment] = Va
    neuron.vh[initial_segment] = Vh
    
    # Initialisation
    neuron.v = EL
    neuron.h = 1

    return neuron
