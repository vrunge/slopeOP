import numpy as np
import os
import multiprocessing as mp
from matplotlib.pyplot import plot
import time

import pickle
import time
from subprocess import Popen, PIPE

from sklearn.metrics import mean_squared_error as mse
import pandas as pd

script_name = 'Rinterface.R'

################ SIGNAL GENERATION ###########################

def signal_1(n=500,noise=0, seed=42):
    np.random.seed(seed)
    return np.interp(np.arange(n), (0,n//2,n-1), (0,60,0)) + np.random.normal(0,noise,n)

def signal_1_at_50pt(seed=42, noise=0):
    return signal_1(n=50, seed=seed, noise=noise)

def signal_1_at_300pt(seed=42, noise=0):
    return signal_1(n=300, seed=seed, noise=noise)

def signal_2(noise=0, seed=42):
    
    data_a = [1] * 125
    data_b = [15] * 125
    data_c = [30] * 125
    data_d = [60] * 125
    data = data_a + data_b + data_c + data_d
    
    np.random.seed(seed)
    return np.array(data + np.random.normal(0,noise,len(data)))

def signal_3(noise=0, seed=42):
    
    a = list(np.linspace(1,30,91))
    b = list(np.linspace(30,1,61))
    c = [x for x in [30]*111]
    data = a[1:] + b[1:] + a[1:] + c[1:] + list(np.array(a)+29) + list(np.array(b)+30)[1:-1]
    
    np.random.seed(seed)
    return np.array(data + np.random.normal(0,noise,len(data)))

def signal_4(noise=0, seed=42):
    
    data_a = list(np.linspace(1,10,65))
    data_b = list(np.linspace(10,30,65))
    data_c = list(np.linspace(30,70,65))
    data_d = list(np.linspace(70,60,65))
    data_e = list(np.linspace(60,30,65))
    data_f = list(np.linspace(30,1,26))
    data_g = list(np.linspace(1,10,96))
    data_h = list(np.linspace(10,40,60))
    data = data_a + data_b[1:] + data_c[1:] + data_d[1:] + data_e[1:] + data_f[1:] + data_g[1:] + data_h[1:]
    
    np.random.seed(seed)
    return np.array(data + np.random.normal(0,noise,len(data)))


def get_run_list(algos, signals, betas, working_directory=os.getcwd(),
                 first_column=None, last_column=None):
    """Return a list of simulation parameters in all possible combinations"""
    run_list = []
    for algo in algos:
        for signal in signals:
            for beta in betas:
                if (first_column is not None) and (last_column is not None):
                    run_list.append([
                        'Rscript', script_name, algo, signal,
                        str(beta), working_directory, str(first_column), str(last_column)])
                else:
                    run_list.append(
                        ['Rscript', script_name, algo, signal, str(beta), working_directory])

    return run_list


def run_command(cmd):
    """run cmd in a subshell and time it"""
    start = time.time()
    # print(f'start {cmd}:', start)
    # err = os.system(' '.join(cmd))
    proc = Popen(cmd, stdout=PIPE, stderr=PIPE)
    proc.wait()
    elapsed = (time.time() - start)
    # print(f'finished {cmd}:', elapsed)
    out = dict(
        algo=cmd[2],
        signal=cmd[3],
        beta=cmd[4],
        # start_time = start,
        elapsed_time=elapsed)

    return out


# PLOTTING FUNCTIONS HELPERS
def get_cpts_from_signal(y, min_d2=0.1):
    """return the changepoints given a signal.
    The signal is assumed to be peachewise linear.

    min_d2: minimum '1st derivative' difference to
    be conisdered a discontinuity"""

    d2y = np.diff(y, 2)
    cpts = [0]
    cp_detected = False
    for i, y_this in enumerate(d2y[1:]):
        if cp_detected:
            cp_detected = False
            continue
        if y_this > min_d2:
            cpts.append(i+3)
            cp_detected = True
    cpts.append(len(y)-1)
    return cpts


def get_linear_segments(y):
    """Transform a piece-wise linear signal into a list of
    its linear segments."""
    cpts = get_cpts_from_signal(y)
    assert(len(cpts) > 0)
    start = cpts[0]
    segments = []
    for cpt in cpts[1:]:
        end = cpt-1
        segments.append(y[start:end])
        start = cpt
    return segments


def plot_picewise_linear(segments, ax, *args, **kwargs):
    """plot all segments in same plot"""
    label = kwargs.pop('label', None)

    for s in segments[:-1]:
        ax.plot(s, *args, **kwargs)

    ax.plot(segments[-1], *args, **kwargs, label=label)


def pickle_dump(obj, fname, dir_path='.'):
    """pickle obj"""
    with open(os.path.join(dir_path, fname), 'wb') as f:
        pickle.dump(obj, f)

def pickle_load(fname, dir_path='.'):
    """pickle obj"""
    with open(os.path.join(dir_path, fname), 'br') as f:
        obj = pickle.load(f)
    return obj

def get_changepoints_and_labels(v, epsilon=1e-10):
    """Return a vector containing the changepoints and segment labels
    correspoinding to each segment of a picewise linear signal

    Args:
        v (np.ndarray): a picewise linear signal
    """

    d = np.diff(v, 2)
    # threshold
    d[abs(d) < epsilon] = 0

    # changepoints
    p = (np.argwhere(d).ravel()+2).astype(int)

    labels = np.zeros(len(v))

    for i in range(0, len(p)-1):
        labels[p[i]:p[i+1]] = i+2

    # last segment
    if len(p) > 0:
        labels[p[-1]:] = len(p)

    return p, labels


def get_labels(v, epsilon=1e-10):
    _, labels = get_changepoints_and_labels(v, epsilon)
    return labels


if __name__ == '__main__':
    coords = np.arange(1, 10)
    for i in range(1000):
        xs = sorted(np.random.choice(coords, size=5, replace=False)*10)
        ys = np.random.choice(coords, size=5, replace=False)+np.random.randn(5)
        x = np.interp(np.arange(0, 100), xs, ys)
        cp, labels = get_changepoints_and_labels(x)

        # check that the discovered changepoints are correct
        if len(xs) != len(cp) or any(xs != cp-1):
            plt.plot(x)
            plt.plot(labels, label='labels')
            plt.legend()
            print(xs, cp, ys)
            raise Exception("wrong labels?")

    cp, labels = get_changepoints_and_labels(np.ones(100))
    assert len(cp) == 0
    assert all(labels == 0)

    print("all test passed")


# ================== plotting =========================== #
def get_mse_vs_beta(signal_name, algo, betas, sigma, seeds, path):
    """Return the log(MSE) vs Beta of a given segmentation algorithm,
    averaged over the different noise realizations"""
    mses = np.zeros(len(betas))
    
    for i, beta in enumerate(betas):
        signal_data = pd.read_csv(
            f"{path}/simu_{signal_name}_{algo}_{beta:.3f}.csv")
        mse_temp = []
        for seed in seeds:
            mse_temp.append(
                mse(signal_data["X0_0"], signal_data[f"X{sigma}_{seed}"]))
        mses[i] = np.mean(mse_temp)
    return np.log(mses)


def get_mses(signal_name, algo, beta, sigma, seeds, path):
    """Return the RMSEs of a given segmentation algorithm,
    for all different noise realizations"""
    mses = []
    signal_data = pd.read_csv(
        f"{path}/simu_{signal_name}_{algo}_{beta:.3f}.csv")
    for seed in seeds:
        mses.append(mse(signal_data["X0_0"], signal_data[f"X{sigma}_{seed}"]))
    return mses


def get_time(signal_name, algo, beta, path):
    """Return the execution times for the specified simulation"""
    time_data = pd.read_csv(
        f"{path}/simu_{signal_name}_{algo}_{beta:.3f}_time.csv").iloc[1:].drop("Unnamed: 0", axis=1).values
    return time_data.flatten()


marker_algo = {
    'SlopeOP': '.',
    'CPOP': '+',
    'OP2D': '*',
    'FPOP': '+',
    'RFPOP': '*'
}
