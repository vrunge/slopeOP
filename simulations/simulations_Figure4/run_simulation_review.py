# %%
import numpy as np
import pandas as pd

# global imports
import os
import shutil
from multiprocessing import Pool, cpu_count
import pandas as pd
from tqdm import tqdm
from random import shuffle
from glob import glob

# local imports
import simulation_tools as tools

FIRST_EXECUTION = True
read_csv = pd.read_csv
tempdir = 'sigma_3'

# %% FIRST EXECUTION SETUP
if FIRST_EXECUTION:
    # ask confirmation
    r = input(
        "WARNING: this will delete old simulation files. Press ENTER to continue...")
    print(r)
    if r != '':
        raise(KeyboardInterrupt("Execution aborted"))

    print("Simulation directory:", tempdir)
    if not os.path.exists(tempdir):
        os.mkdir(tempdir)
    else:
        print("remove old simulation files")
        os.system(f"rm -f {tempdir}/*")

    print("copy CPOP code")
    filenames = glob(os.path.join("CPOP_code", "*"))
    for filename in filenames:
        shutil.copy(filename, tempdir)

    FIRST_EXECUTION = False
else:
    print("Already executed, doing nothing.")


# %% SIMULATION PARAMETERS
# number of repetition of each simulation for averaging
REPETITIONS = 150
beta_min = 0
beta_step = 0.1
beta_max = 4.5

betas = np.arange(beta_min+beta_step, beta_max, beta_step)
seeds = np.arange(1, REPETITIONS+1)
sigma = 3

print("sigma", sigma)
print("seeds", REPETITIONS)

metadata = {'tempdir': tempdir, 'betas': betas, 'sigma': sigma, 'seeds': seeds, 'repetitions':REPETITIONS}
tools.pickle_dump(metadata, 'simulation_metadata.pickle', dir_path=tempdir)
print("metadata saved")

#%% generate signals
for f in [tools.signal_1, tools.signal_2, tools.signal_3, tools.signal_4,]:
    
    signal={}

    column_format = "{sigma}_{seed}"
    column_name = column_format.format(sigma=0, seed=0)
    signal[column_name] = f(noise=0)
    

    for seed in seeds:
        column_name = column_format.format(sigma=sigma, seed=seed)
        signal[column_name] = f(noise=sigma,seed=seed)

    # Write data
    file_name = os.path.join(tempdir, f"{f.__name__}.csv")
    pd.DataFrame(signal).to_csv(file_name, index=False)


# %% Preparationo of the run list
run_list = tools.get_run_list(
    algos=['CPOP', 'SlopeOP', 'OP2D'],
    signals=['signal_1.csv', 'signal_2.csv', 'signal_3.csv', 'signal_4.csv'],
    betas=betas,
    working_directory=tempdir
)

print(len(run_list), "simulation commands generated")
tools.pickle_dump(run_list, 'run_list.pickle', dir_path=tempdir)


# %% RUN SIMULATION

shuffle(run_list)

cores = 10
chunksize = len(run_list)//cores

# # WITH TQDM
# with Pool(cores) as pool:
#     work_results = list(tqdm(
#         pool.imap(
#             tools.run_command, run_list, chunksize=1),
#         total=len(run_list)))

# BARE
with Pool(cores) as pool:
    work_results = pool.map(tools.run_command, run_list, chunksize=chunksize)


print("Done")
