import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt

whole_file = 'whole.csv'
with open(whole_file, 'w') as ww:
    ww.writelines('t[s],maxT[degC],Tatquat[degC],HeatGain[-]\n')

csv_files = glob.glob('DATA/*.csv')
for idx, csv_file in enumerate(sorted( csv_files )):
    df = pd.read_csv(csv_file)
    timestep = csv_file.split('_')[1][:4]

    maxT = df['T[degC]'].max()
    HG = df['HeatGain[-]'].max()
    quadIdx = np.argmin( abs( df[' x[m]'].values - 1.5e-2 ) )
    quadT = df['T[degC]'].iloc[quadIdx]

    with open(whole_file, 'a') as wa:
        wa.writelines(f'{float(timestep):.6e},{maxT:.6e},{quadT:.6e},{HG:.6e}\n')
