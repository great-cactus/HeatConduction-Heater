import pandas as pd
import glob
import matplotlib.pyplot as plt

csv_files = glob.glob('DATA/*.csv')
for idx, csv_file in enumerate(csv_files):
    df = pd.read_csv(csv_file)
    timestep = csv_file.split('_')[1][:4]

    fig, ax = plt.subplots(1,1)
    secax = ax.twinx()
    ax.set_title(f'{float( timestep ):.2e} ms')
    ax.plot(df[' x[cm]'], df['T[K]'], color='black')
    secax.plot(df[' x[cm]'], df['conductivity[-]'], label='cond')
    secax.plot(df[' x[cm]'], df['HeatTransfer[-]'], label='Heat Transfer')
    secax.plot(df[' x[cm]'], df['HeatGain[-]']    , label='Heat Gain')
    secax.set_ylabel('Transport')
    secax.legend(loc='upper left')
    ax.set_xlabel('x [m]')
    ax.set_ylabel('T [K]')
    ax.set_ylim(0, 1000)
    ax.set_xlim(0,0.1)

    plt.tight_layout()
    plt.savefig(f'PNG/output_{timestep}.png', dpi=300)
    plt.cla()
    plt.close()
    del fig, ax
