__author__ = 'marf'

import pandas as pd
from parameter import *
from dspike_formulas import *
import numpy as np
import flatdict
import collections

spike1 = [["116"],["117", "120", "122"]]

# Define Model Parameter#
#*** Natural Fractionation ***#
fnat_sim = -1
#*** Instrumental Fractionation ***#
fins_sim =  2.2
#*** Sample/Spike Ratio ***#
mix = 0.5

    # Spike Calculation#

sim1 = calc_dspike_samples(Sn_meas_obj,df_new,spike_obj,Sn_mass_obj,spike1,"120")
sim1.spike_sim(fnat_sim,fins_sim,mix,3,6,-0.1,-2,'z')

log_dict = {}
counter = 0

for calc_dspike_object in calc_dspike:
    log_dict[counter] = calc_dspike_object.log_dict
    counter += 1

log_df = pd.DataFrame.from_dict(log_dict, orient='index')
print log_df
log_df.to_csv(path + "Sn117-122_116_1367.csv")