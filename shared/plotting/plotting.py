'''
Plotting tools
'''
__all__ = ['prepare_panel', 'color_axon', 'color_soma', 'color_gradient', 'color_soma_Iaxon']

# colors
color_soma = 'k'
color_soma_Iaxon = 'b'
color_axon = 'r'
color_gradient = 'g'

def prepare_panel(*ax):
    for axis in ax:
        axis.spines['right'].set_visible(False)
        axis.spines['top'].set_visible(False)
