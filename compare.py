import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ref_df = pd.read_csv('ref.csv')
whole_df = pd.read_csv('whole.csv')

ref_df = ref_df[ref_df['time'] < 80]
whole_df = whole_df[whole_df['t[s]'] < 80]

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(15, 10))

#.... original temperature time history
ax1.set_title('Original Temperature Time History')
ax1.set_xlabel('Time (s)')
ax1.set_ylabel('Temperature (degC)')
ax1.plot(whole_df['t[s]'], whole_df['maxT[degC]']   , label='maxT'        , color='tab:red', linestyle='dashed')
ax1.plot(whole_df['t[s]'], whole_df['Tatquat[degC]'], label='TC1 location', color='tab:blue', linestyle='dashed')
ax1.plot(ref_df['time']  , ref_df['Tctip[degC]']    , label='TCtip'       , color='tab:red')
ax1.plot(ref_df['time']  , ref_df['TC1[degC]']      , label='TC1'         , color='tab:blue')
ax1.legend()

#.... normalized temperature time history
normed_maxT  = whole_df['maxT[degC]'].values/whole_df['maxT[degC]'].max()
normed_quadT = whole_df['Tatquat[degC]'].values/whole_df['maxT[degC]'].max()
normed_TCtip = ref_df['Tctip[degC]'].values/ref_df['Tctip[degC]'].max()
normed_TC1   = ref_df['TC1[degC]'].values/ref_df['Tctip[degC]'].max()

ax2.set_title('Normalized Temperature Time History')
ax2.set_xlabel('Time (s)')
ax2.set_ylabel('Normalized Temperature (-)')
ax2.plot(whole_df['t[s]'], normed_maxT , label='maxT'        , color='tab:red' , linestyle='dashed')
ax2.plot(whole_df['t[s]'], normed_quadT, label='TC1 location', color='tab:blue', linestyle='dashed')
ax2.plot(ref_df['time']  , normed_TCtip, label='TCtip'       , color='tab:red')
ax2.plot(ref_df['time']  , normed_TC1  , label='TC1'         , color='tab:blue')
ax2.legend()

#.... Temperature hysteresis
ax3.set_title('Temperature hysteresis')
ax3.set_xlabel('High temperature (degC)')
ax3.set_ylabel('Low temperature (degC)')
ax3.scatter(ref_df['Tctip[degC]'] , ref_df['TC1[degC]']      , label='original', color='tab:orange')
ax3.scatter(whole_df['maxT[degC]'], whole_df['Tatquat[degC]'], label='model'   , color='tab:blue')
ax3.legend()

#.... normalized temperature hysteresis
ax4.set_title('Normalized temperature hysteresis')
ax4.set_xlabel('Normalized high temperature (degC)')
ax4.set_ylabel('Normalized low temperature (degC)')
ax4.scatter(normed_TCtip, normed_TC1  , label='original', color='tab:orange')
ax4.scatter(normed_maxT , normed_quadT, label='model'   , color='tab:blue')
ax4.legend()

ax1.set_xlim(0., 80)
ax2.set_xlim(0., 80)
fig.tight_layout()
plt.show()

fig.savefig('FigCompare.png', dpi=300)

