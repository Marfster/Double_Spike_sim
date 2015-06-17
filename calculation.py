__author__ = 'marf'

import pandas as pd
from parameter import *
from dspike_formulas import *
import numpy as np

spike1 = [["116"],["117", "120", "122"]]

# Define Model Parameter#
#*** Natural Fractionation ***#
fnat = -1
#*** Instrumental Fractionation ***#
fins =  2.2
#*** Sample/Spike Ratio ***#
mix = 0.5


# Mix Sample-Spike
def fract(x, pos, ma1, ma2, frac):
        X = x[pos] * (ma1/ma2) ** frac
        return X

Sn_ratios_meas_sim_dict = {}
for ratio in Sn_ratios_meas.get_all_ratios("120"):
    Sn_ratios_meas_sim_dict[ratio] = fract(Sn_ratios_meas.get_all_ratios("120"), ratio,
                                             Sn_masses.get_Isotope_mass(ratio),
                                             Sn_masses.get_Isotope_mass("120"), fnat)

Sn_ratios_meas_sim = Isotope_Ratios()
Sn_ratios_meas_sim.add_ratios_dict("120", Sn_ratios_meas_sim_dict)
Sn_abund_meas_sim = Isotope_Abundances()
Sn_abund_meas_sim.add_abundances_dict(Sn_ratios_meas_sim.get_all_abundances("120"))


sample_spike_mix_abund_dict = {}

for ratio in Sn_abund_spike.get_all_abundances():
    sample_spike_mix_abund_dict[ratio] = mix * Sn_abund_meas_sim.get_all_abundances()[ratio] + (1 - mix) * Sn_abund_spike.get_all_abundances()[ratio]

sample_spike_mix_abund = Isotope_Abundances()
sample_spike_mix_abund.add_abundances_dict(sample_spike_mix_abund_dict)
sample_spike_mix_ratio = sample_spike_mix_abund.get_all_ratios("120")

sample_spike_mix_ratio_sim = {}
for ratio in sample_spike_mix_ratio:
    sample_spike_mix_ratio_sim[ratio] = fract(sample_spike_mix_ratio, ratio,
                                                 Sn_masses.get_Isotope_mass(ratio),
                                                 Sn_masses.get_Isotope_mass("120"),fins)
#simulate mixture

nat_frac = []
for value in range(len(df_new)):
    mix = {}
    for isotope in df_new:
         mix[isotope] = (sample_spike_mix_ratio_sim[isotope]/df_new[isotope][value])*df_new[isotope].mean()

    mix_ratios = Isotope_Ratios()
    mix_ratios.add_ratios_dict("120", mix)
    mix_abund = Isotope_Abundances()
    mix_abund.add_abundances_dict(mix_ratios.get_all_abundances("120"))

    # Spike Calculation#

    std = dspike_formulas("std", Sn_meas_obj, spike_obj, Sn_mass_obj, spike1)
    sample = dspike_formulas("mix1", mix_abund, spike_obj, Sn_mass_obj, spike1)
    dspike_single = calc_dspike()

    nat_frac.append(dspike_single.dspike_calc(std,sample,3,6,-1,-2.2,-0.1,-2))

print nat_frac
print np.mean(nat_frac)
print np.std(nat_frac)
print (np.std(nat_frac)/np.mean(nat_frac))*1000000
print ((np.mean(nat_frac)/fnat)-1)*1000000

for calc_dspike_object in calc_dspike:
    print calc_dspike_object.N
    print calc_dspike_object.frac_ins
    print calc_dspike_object.frac_nat