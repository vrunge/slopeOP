import numpy as np
import os
import asyncio
from matplotlib.pyplot import plot


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
    data_c = list(np.linspace(30,60,65))
    data_d = list(np.linspace(60,50,65))
    data_e = list(np.linspace(50,30,65))
    data_f = list(np.linspace(30,1,26))
    data_g = list(np.linspace(1,10,96))
    data_h = list(np.linspace(10,40,60))
    data = data_a + data_b[1:] + data_c[1:] + data_d[1:] + data_e[1:] + data_f[1:] + data_g[1:] + data_h[1:]
    
    np.random.seed(seed)
    return np.array(data + np.random.normal(0,noise,len(data)))


########################## simulation
def get_run_list(algos, signals, betas, working_directory=os.getcwd(),
                 first_column=None,last_column=None):
    """Return a list of simulation parameters in all possible combinations"""
    script_name = 'algorithmsDataFrame05_10.R'
    run_list = []
    for algo in algos:
        for signal in signals:
            for beta in betas:
                if (first_column is not None) and (last_column is not None):
                    run_list.append(['Rscript', script_name, algo, signal, str(beta), working_directory, str(first_column), str(last_column)])
                else:
                    run_list.append(['Rscript', script_name, algo, signal, str(beta), working_directory])

    return run_list


async def process_run(cmd,verbose=False):
    proc = await asyncio.create_subprocess_shell(
        cmd,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE)

    stdout, stderr = await proc.communicate()

    if verbose:
        print(f'[{cmd!r} exited with {proc.returncode}]')
        #if stdout:
        #    print(f'[stdout]\n{stdout.decode()}')
    if stderr:
        print(f'[stderr]\n{stderr.decode()}')

async def simulation_run(run_list):
    await asyncio.gather(
        *[process_run (" ".join(x)) for x in run_list]
    )
    
## PLOTTING FUNCTIONS HELPERS
def get_cpts_from_signal(y, min_d2=0.1):
    """return the changepoints given a signal.
    The signal is assumed to be peachewise linear.
    
    min_d2: minimum '1st derivative' difference to
    be conisdered a discontinuity"""
    
    d2y = np.diff(y,2)
    cpts=[0]
    cp_detected=False
    for i, y_this in enumerate(d2y[1:]):
        if cp_detected:
            cp_detected=False
            continue
        if y_this>min_d2:
            cpts.append(i+3)
            cp_detected=True
    cpts.append(len(y)-1)
    return cpts

def get_linear_segments(y):
    """Transform a piece-wise linear signal into a list of
    its linear segments."""
    cpts=get_cpts_from_signal(y)
    assert(len(cpts)>0)
    start=cpts[0]
    segments=[]
    for cpt in cpts[1:]:
        end=cpt-1
        segments.append(y[start:end])
        start=cpt
    return segments

def plot_picewise_linear(y,*args, **kwargs):
    """plot a piecewise linear signal as disconnected segments"""
    segs=get_linear_segments(y)
    label=kwargs.pop('label',None)
    
    for s in segs[:-1]:
        plot(s,*args,**kwargs)
    
    plot(segs[-1],*args,**kwargs,label=label)
