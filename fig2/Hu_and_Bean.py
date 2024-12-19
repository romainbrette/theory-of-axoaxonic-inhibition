'''
Injection of a -50 pA hyperpolarizing current at the soma and 75 um away in the axon.
Data from Hu and Bean (2018): https://zenodo.org/record/3539297
'''
from shared.data.hu_bean_shared import *
from shared.plotting import *
from pylab import *
from matplotlib import gridspec
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

figsize = (5,4) # width, height

path = '../shared/data/'

distance, t, Vs, Va, Is, Ia = load_file(path + '100420 1-1 75um.abf')
t *= 1000 # in ms

before, after = 200, 500+200
t0 = [19369, 110853]
t_at_soma = (t>t0[0]-before) & (t<t0[0]+after)
t_at_axon = (t>t0[1]-before) & (t<t0[1]+after)

gs = gridspec.GridSpec(5, 2)
fig = figure('Data', figsize)

for k, interval in enumerate([t_at_soma, t_at_axon]):
    ax_I = fig.add_subplot(gs[0:1, k])
    ax_V = fig.add_subplot(gs[1:3, k])
    ax_dV = fig.add_subplot(gs[3:5, k])
    prepare_panel(ax_I, ax_V, ax_dV)

    ax_I.plot(t[interval]-t0[k], Is[interval]*1000, color_soma, label='soma')
    ax_I.plot(t[interval]-t0[k], Ia[interval]*1000, color_axon, label='axon, 75 µm')
    ax_V.plot(t[interval]-t0[k], Vs[interval], color_soma)
    ax_V.plot(t[interval]-t0[k], Va[interval], color_axon)
    ax_dV.plot(t[interval]-t0[k], Va[interval] - Vs[interval], color_gradient, label='difference')

    ax_I.set_xticklabels([])
    ax_V.set_ylim(-70, -58)
    ax_V.set_xticklabels([])
    ax_dV.set_ylim(-4, 2)
    ax_dV.set_xlabel('t (ms)')
    if k==0:
        ax_I.set_ylabel('I (pA)')
        ax_V.set_ylabel('V (mV)')
        ax_dV.set_ylabel(r'$\Delta$V (mV)')

        ## Legend
        lines, labels = ax_I.get_legend_handles_labels()
        lines2, labels2 = ax_dV.get_legend_handles_labels()
        ax_dV.legend(lines + lines2, labels + labels2, loc='best')

## Inset
''''''
ax_inset = ax_dV.inset_axes([0.33, 0.4, 0.3, 0.5])
ax_inset.plot(t[interval]-t0[k], Va[interval] - Vs[interval], 'g')
ax_inset.set_xlim(-5, 10)
ax_inset.set_ylim(-4, 0)
ax_inset.patch.set_edgecolor('black')
ax_inset.patch.set_linewidth(1)
ax_inset.set_yticklabels([])
mark_inset(ax_dV, ax_inset, loc1=2, loc2=4, fc="none", ec="0.5")

tight_layout()

savefig('Hu_and_Bean.png', dpi=300)

show()
