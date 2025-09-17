import os
from datetime import datetime
import datetime as dt
from pyproj import Transformer
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import fnmatch
from pyproj import Transformer
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from pyproj import Transformer
from scipy.stats import pearsonr
import sys
import warnings

parent_dir = os.path.abspath('C:/Users/Will.Rust/OneDrive - Cranfield University/postdoc/Environment/Projects/RESTRECO/remote sensing WP/sweep_code/source')
sys.path.append(parent_dir)

from dlm import run_dlm
from helper_funcs import (
    calculate_bounds,
    to_df,
    get_api_response,
    wave_variance,
    plot_dlm_results,
    plot_gcc_comparison,
    plot_sweep_location,
    plot_sweep_extract,
    ews
)
from sweep import (
    sweep_locate,
    sweep_infil,
    sweep_extract,
    plot_landscape,
    plot_landscape_gif,
    synth_max_wp
)

#----------------------
#IMPORT AND CLEAN NDVI
#----------------------

ndvi_in = pd.read_csv('HarzForest_raw_ndvi.csv', parse_dates=[1], date_parser=lambda x: pd.to_datetime(x, format="%Y_%m_%d", errors="coerce"))
#ndvi_in = pd.read_csv('DinkeyCreekForest_raw_ndvi.csv', parse_dates=[1], date_parser=lambda x: pd.to_datetime(x, format="%Y_%m_%d", errors="coerce"))
ndvi_in = ndvi_in.iloc[:, [1, 4]]

ndvi_in.head()

plt.plot(ndvi_in['ndvi'])
plt.show()

#-----------------------------
#IMPORT AND CLEAN CLIMATE DATA
#-----------------------------

climate_in = pd.read_csv('Rainfall_Anom_Harz.csv', parse_dates=[2], date_parser=lambda x: pd.to_datetime(x, format="%Y%m%d", errors="coerce"))
#climate_in = pd.read_csv('Rain_Anom_Dinkey_Fresno.csv', parse_dates=[2], date_parser=lambda x: pd.to_datetime(x, format="%Y%m%d", errors="coerce"))

#DINKY climate data has multiple stations (i assume?) Take only Freiburg
#climate_in = climate_in[climate_in["Name"] == "Freiburg"]

climate_in = climate_in.iloc[:, [2, 3]]

climate_in['month'] = climate_in['Date'].dt.month
monthly_mean_precip = climate_in.groupby(['month'])['Rainfall_mm'].transform('mean')
climate_in['precip_an'] = climate_in['Rainfall_mm'] - monthly_mean_precip
climate_in = climate_in.fillna({'precip_an':0})

#------------------------------
#MERGE NDVI AND CLIAMTE
#------------------------------

data_merged = pd.merge(ndvi_in, climate_in, left_on="id", right_on="Date", how="inner")

x = data_merged['ndvi'].values
anCLM = data_merged['precip_an']

#------------------------------
#CALCUALTE SWEEP
#------------------------------

#sweep params
nseas = [1, 2]
rseas = nseas
wav = ('morlet', {'mu': 6})
vid = 2

l_period =  128 
u_period = 365.25 * 1.5
min_thresh = 0.3
fs = 1
buff = 0.1
sig_lvl = 0.05
cost = 0.5
cost_i = 0
quantile = 0.9     

exp = 1
mean_step = 0.05
range_step = mean_step/2

win_mean = np.arange(mean_step, 1 + mean_step, mean_step)
win_range = np.arange(range_step, 1 + range_step, range_step)

#calculating significant wavelet power takes a while
#n = 1000
#sorted_synth_wp = synth_max_wp(win_mean, win_range, wav, fs, n, l_period, u_period)
#sig_ind = int(np.ceil((1 - sig_lvl) * n)) - 1
#sig_wp = sorted_synth_wp[sig_ind]
sig_wp = 0.13964520394802094


# NaNs are gaps for the DLM to handle
with warnings.catch_warnings():
    warnings.filterwarnings('ignore', r'All-NaN (slice|axis) encountered')
    warnings.filterwarnings('ignore', r'Mean of empty slice')
    
    signal_loc = sweep_locate(
        x,
        win_mean,
        win_range,
        wav,
        fs,
        l_period,
        u_period,
        min_thresh,
        sig_wp,
        cost_i,
        buff
    )

    signal_ind = sweep_infil(
        signal_loc['signal_index'],
        signal_loc['signal_location'],
        signal_loc['best_power'],
        signal_loc['sig_vector'],
        signal_loc['av_power_mat'],
        win_mean, win_range, buff
    )

    x_sweep, wx_og, wx_filt, factor_mask = sweep_extract(
        x,
        signal_ind['signal_location_infil'],
        wav,
        fs,
        l_period,
        exp
    )

sweep_loc_plot = plot_sweep_location(
    data_merged["id"],
    x,
    x_sweep,
    signal_loc['sig_vector'],
    signal_ind['signal_location_infil'],
    signal_ind['signal_index_infil'],
    signal_loc['center_mat_w']
)

sweep_loc_plot.show()
sweep_loc_plot.savefig("Harz_sweep_loc.png", dpi=300, bbox_inches="tight")

sweep_extract_plot = plot_sweep_extract(
    data_merged["id"],
    np.nan_to_num(x),
    x_sweep,
    signal_ind['signal_location_infil'],
    wx_og,
    wx_filt,
    l_period,
    wav
)

sweep_extract_plot.show()
sweep_extract_plot.savefig("Harz_sweep_extract.png", dpi=300, bbox_inches="tight")

#--------------------------
#RUN DLM
#--------------------------

fs1_deltas = np.ones(4) * 0.9967

sm_sweep, sC_sweep, snu_sweep, FF_sweep, *_ = run_dlm(
    x_sweep.values,
    anCLM,
    vid,
    1,
    rseas,
    fs1_deltas
)

warmup = int(365.25 * 1.5)

point_size = 0.7

sm_sweep = FF_sweep['sm'][2,:]
sC_sweep = FF_sweep['sC'][2,2,:]
snu_sweep = FF_sweep['snu']

sweep_lbounds, sweep_ubounds = calculate_bounds(
    sm_sweep, sC_sweep, snu_sweep, quantile)

sweep_ews_ind = ews(quantile, sm_sweep, sC_sweep, snu_sweep,
                    warmup, 1, 365)

sweep_ews = np.full(len(x_sweep), np.nan)
if len(sweep_ews_ind) != 0:
    sweep_ews[sweep_ews_ind] = sm_sweep[sweep_ews_ind]

fig, ax = plt.subplots(2, 1, figsize=(4, 4),
                        gridspec_kw={'width_ratios': [1],
                                    'height_ratios': [1, 2]})

daily_dates = data_merged['id']
sig_vector = signal_loc['sig_vector']
cosmos_vwc = data_merged['precip_an']

ax[0].scatter(daily_dates[sig_vector], x_sweep[sig_vector],
                color='black', label='SWEEP (sig.)', s=point_size)
ax[0].scatter(daily_dates[~sig_vector], x_sweep[~sig_vector],
                color='black', label='SWEEP (non-sig.)', marker='+',
                s=point_size + 5, linewidths=0.3)
ax[0].legend(loc='lower right', fontsize=6, ncol=3)
ax[0].set_ylabel('NDVI', fontsize=6)
ax[0].tick_params(axis='x', labelsize=6)
ax[0].tick_params(axis='y', labelsize=6)

ax[1].fill_between(daily_dates, -sweep_lbounds, -sweep_ubounds,
                    facecolor='lightgrey', alpha=0.8, zorder=4)
ax[1].plot(daily_dates, -sm_sweep, color='black', linestyle='dashed',
            label='SWEEP Speed', linewidth=0.6, zorder=5)
ax[1].plot(daily_dates, -sweep_ews, color='black', label='SWEEP CSD',
            linewidth=2, zorder=6)

ax[1].set_ylim(-1, 1)
ax[1].legend(loc='upper right', fontsize=6, ncol=2)
ax[1].set_ylabel('System Speed (-AC1)', fontsize=6)
ax[1].set_xlabel('Date', fontsize=6)
ax[1].tick_params(axis='x', labelsize=6)
ax[1].tick_params(axis='y', labelsize=6)

for ax_i in [ax[0], ax[1]]:
    plt.setp(ax_i.get_xticklabels(), visible=False)
    ax_i.tick_params(axis='x', which='both', bottom=False, top=False)

fig.show()
fig.savefig("Harz_dlm.png", dpi=300, bbox_inches="tight")
