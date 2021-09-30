#%%
import matplotlib.pylab as plt
from nice_plotting import savefig, set_rcParams
import simulation_tools as tools
import os


#%%

set_rcParams()

dirs = ['sigma_24', 'sigma_12', 'sigma_3']
fig, axs = plt.subplots(len(dirs), 4, figsize=(15, 8), sharey='row', sharex='col')

shift=40

for j,tempdir in enumerate(dirs):
    plot_data = tools.pickle_load(f"plot_data_{tempdir}.pkl", f"Figures/{tempdir}")
    for i,signal in enumerate(plot_data):
        data = plot_data[signal]
        axs[j,i].plot(data['SlopeOP']+shift, label='SlopeOP', lw=4)
        axs[j,i].plot(data['signal'], 'k--', label='true signal', lw=2)
        segments = tools.get_linear_segments(data['OP2D']-shift)
        tools.plot_picewise_linear(segments, axs[j,i], 'C2', label='OP2D', lw=4)
        # axs[j,i].set_ylim(-50,160)
        if i==0:
            axs[j,i].set_ylabel(tempdir.replace('_', ' = ').replace("sigma","$\\sigma$"))
        if j==0:
            axs[j,i].set_title(signal.replace('_', ' '))

axs[0,0].legend(loc ='upper left')
fig.suptitle('Segmentations with best MSE', fontsize=16)

fig_dir = f"Figures"
if not os.path.exists(fig_dir):
    os.mkdir(fig_dir)

savefig(f"{fig_dir}/best_MSE_segmentation.pdf")
# %%
