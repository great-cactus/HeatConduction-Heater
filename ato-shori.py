import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt

whole_file = 'whole.csv'
with open(whole_file, 'w') as ww:
    ww.writelines('t[s],maxT[K],Tatquat[K]\n')

csv_files = glob.glob('DATA/*.csv')
for idx, csv_file in enumerate(csv_files):
    df = pd.read_csv(csv_file)
    timestep = csv_file.split('_')[1][:4]

    maxT = df['T[K]'].max()
    quadIdx = np.argmin( abs( df[' x[cm]'].values - 0.04 ) )
    quadT = df['T[K]'].iloc[quadIdx]

    with open(whole_file, 'a') as wa:
        wa.writelines(f'{float(timestep):.6e},{maxT:.6e},{quadT:.6e}\n')
