# %% IMPORTS
from multiprocessing import Pool
import itertools as itt
import slopeOP as sop
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pickle
from sklearn.metrics import adjusted_rand_score
import os

# local imports
import simulation_tools as tools
from nice_plotting import savefig, set_standard_yrange, set_rcParams

read_csv = pd.read_csv

set_rcParams()

tempdir = 'simulation'


with open(f"{tempdir}/simulation_metadata.pickle", 'rb') as f:
    metadata = pickle.load(f)

# sigma = metadata['sigma']
sigma = 12
betas = metadata['betas']
seeds = metadata['seeds']
REPETITIONS = metadata['repetitions']

signals = ["signal_1", "signal_2", "signal_3", "signal_4"]

# %% Calc MSE and simulation run time

mse_vs_beta = dict()
exec_times = dict()

print("calc MSE and simulation run times...", end='')


for signal in signals:

    # ====== average log(MSE) vs Beta ====== #

    mse_vs_beta[signal] = {  # data_col1
        'SlopeOP': tools.get_mse_vs_beta(signal, "SlopeOP", betas, sigma, seeds, tempdir),
        'CPOP': tools.get_mse_vs_beta(signal, "CPOP", betas, sigma, seeds, tempdir),
        'OP2D': tools.get_mse_vs_beta(signal, "OP2D", betas, sigma, seeds, tempdir)
    }

    algos = mse_vs_beta[signal].keys()

    data_col2 = {}
    minMSEs = {}

    # ====== Simulation run times ====== #
    exec_times[signal] = []
    for beta in betas:
        for algo in algos:
            time_all_seeds = tools.get_time(signal, algo, beta, tempdir)
            exec_times[signal].append(
                {
                    'algo': algo, 'beta': beta,
                    'mean': time_all_seeds.mean(),
                    'std': time_all_seeds.std()
                })

print("Done")

# %% calc ARI


algos = ['SlopeOP', 'CPOP', 'OP2D']

print("Calc ARI...", end='')


def work(signal, algo, beta):
    dfs = (pd.read_csv(f'{tempdir}/simu_{signal}_{algo}_{beta:.3f}.csv')
           .rename(columns={'Unnamed: 0': 'seed'})
           .set_index('seed')
           )
    dfs_labels = dfs.apply(tools.get_labels, axis=0).astype(int)

    out = []

    for x in dfs:
        out.append({
            'algo': algo,
            'seed': x,
            'signal': signal,
            'beta': beta,
            'ari': adjusted_rand_score(dfs_labels[x], dfs_labels.X0_0)
        })

    return out


arglist = itt.product(signals, algos, betas)

with Pool(40) as pool:
    result = pool.starmap(work, arglist)

ari = list(itt.chain(*result))

ari_mean = pd.DataFrame(ari).groupby(['signal', 'algo', 'beta']).mean()
ari_std = pd.DataFrame(ari).groupby(['signal', 'algo', 'beta']).std()

print("Done")

# %%

fig, axs = plt.subplots(4, 3, figsize=(15, 10), sharex='col')

minMSEs = dict()
for i, signal in enumerate(signals):
    for algo in mse_vs_beta[signal]:

        # Column 1
        MSE_argmin = np.argmin(mse_vs_beta[signal][algo])
        MSE_min = mse_vs_beta[signal][algo][MSE_argmin]
        MSE_min_beta = betas[MSE_argmin]
        minMSEs[signal,algo] = (round(MSE_min_beta, 2), MSE_min)

        axs[i, 0].scatter(betas, mse_vs_beta[signal][algo],
                          marker=tools.marker_algo[algo], label=algo)
        if i == 0 and algo == "OP2D":
            axs[i, 0].scatter(MSE_min_beta, MSE_min, s=50,
                              c='k', marker='x', label="min MSE")
        else:
            axs[i, 0].scatter(MSE_min_beta, MSE_min, s=50, c='k', marker='x')

    # Column 2
        axs[i, 1].errorbar(
            betas,
            ari_mean.loc[signals[i], algo].values,
            # yerr=ari_std.loc[signals[i], algo].values.ravel(),
            marker=tools.marker_algo[algo], label=algo)

    # Column 3
        df = pd.DataFrame(exec_times[signal]).set_index(['algo'])

        axs[i, 2].errorbar(
            df.loc[algo]['beta'],
            df.loc[algo]['mean'].values,
            yerr=df.loc[algo]['std'].values,
            ls='-',
            marker=tools.marker_algo[algo],
            label=algo)

for i in range(4):
    axs[i, 1].set_ylim(0, 1.1)


# axs[0,0].get_shared_x_axes().join((axs[i,0] for i in range(0,4)))
# for i in range(4):
#     axs[i,2].set_yticklabels([])

# set optimal yrange
yrange_cols = [[0, 0], [0, 0], [0, 0]]
for i in range(3):
    for ax in axs[:, i]:
        temp = ax.get_ybound()
        if temp[0] < yrange_cols[i][0]:
            yrange_cols[i][0] = temp[0]
        if temp[1] > yrange_cols[i][1]:
            yrange_cols[i][1] = temp[1]


set_standard_yrange(axs[:, 2], yrange=yrange_cols[2], xrange=None)
set_standard_yrange(axs[:, 0], yrange=yrange_cols[0], xrange=None)
# set_standard_yrange(axs[:,1], yrange=yrange_cols[1], xrange=None)

#  TITLES & LABELS
axs[0, 0].set_title("log(MSE) vs b")
axs[0, 1].set_title(f"ARI over {REPETITIONS} simulations")
axs[0, 2].set_title("exec time")
axs[-1, 2].set_ylabel("time [s]")

axs[-1, 0].set_xlabel("b")
axs[-1, 1].set_xlabel("b")
axs[-1, 2].set_xlabel("b")


axs[0, 0].legend(loc="upper right", prop={'size': 12}, ncol=2)
# axs[0,2].legend(loc='lower left', prop={'size': 12}, ncol=2)

fig_dir = f"Figures/{tempdir}"
if not os.path.exists(fig_dir):
    os.mkdir(fig_dir)

savefig(f"{fig_dir}/beta-mse-boxplots_sigma{sigma}.pdf")
# %%

### Segmentation

plot_data = {}

for i,signal in enumerate(signals):
    data_slopeOP = read_csv(f"{tempdir}/simu_{signal}_SlopeOP_{minMSEs[(signal,'SlopeOP')][0]:.3f}.csv")[f"X{sigma}_1"]
    data_signal = read_csv(f"{tempdir}/{signal}.csv")["0_0"]
    data_OP2D = read_csv(f"{tempdir}/simu_{signal}_OP2D_{minMSEs[(signal,'OP2D')][0]:.3f}.csv")[f"X{sigma}_1"]
    # plt.title(signal.replace('_',' '))
    plot_data[signal] = {
        'SlopeOP':data_slopeOP,
        'signal':data_signal,
        'OP2D':data_OP2D
    }

tools.pickle_dump(plot_data, f"plot_data_{tempdir}.pkl", f"Figures/{tempdir}")
# %%

# %%
