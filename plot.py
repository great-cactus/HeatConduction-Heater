import pandas as pd
import glob
import matplotlib.pyplot as plt
#from mpi4py import MPI

def process_csv_file(csv_file):
    df = pd.read_csv(csv_file)
    timestep = csv_file.split('_')[1][:4]

    fig, ax = plt.subplots(1,1)
    secax = ax.twinx()
    ax.set_title(f'{float( timestep ):.2e} ms')
    ax.plot(   df[' x[m]'], df['T[degC]'], color='black')
    secax.plot(df[' x[m]'], df['conductivity[-]'], label='cond')
    secax.plot(df[' x[m]'], df['radiation[-]']   , label='rad')
    secax.plot(df[' x[m]'], df['HeatTransfer[-]'], label='Heat Transfer')
    secax.plot(df[' x[m]'], df['HeatGain[-]']    , label='Heat Gain')
    secax.set_ylabel('Transport')
    secax.legend(loc='upper left')
    ax.set_xlabel('x [m]')
    ax.set_ylabel('T [degC]')
    ax.set_ylim(200, 1500)
    ax.set_xlim(0,2.45e-2)

    plt.tight_layout()
    plt.savefig(f'PNG/output_{timestep}.png', dpi=300)
    plt.cla()
    plt.close()
    del fig, ax

if __name__ == '__main__':
    # Initialize MPI
    #comm = MPI.COMM_WORLD
    #rank = comm.Get_rank()
    #size = comm.Get_size()

    # List all csv files
    total_files = sum(1 for _ in glob.iglob('DATA/*.csv'))

    # Divide the tasks among processes
    #tasks_per_process = total_files // size
    #start_index = rank * tasks_per_process
    #end_index = (rank + 1) * tasks_per_process if rank != size - 1 else total_files

    # Process the allocated csv files
    for file_index, csv_file in enumerate(glob.iglob('DATA/*.csv')):
        #if start_index <= file_index < end_index:
        process_csv_file(csv_file)
