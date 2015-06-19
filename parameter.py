__author__ = 'marf'

from iso_properties import *
import pandas as pd

# Masses of Isotopes
#*** Sn after IUPAC - De Laeter (2003)***#
Sn_masses = {"112" : 111.904823, "114" : 113.902783, "115" : 114.903347, "116" : 115.901745, "117" : 116.902955,
             "118" : 117.901608, "119" : 118.903311, "120" : 119.9021985, "122" : 121.9034411, "124" : 123.9052745}

# True Ratios
#*** Sn after Lee et al. (1995)***#
Sn_ratios_Lee = {"112" : 0.029812, "114" : 0.020195, "115" : 0.010366, "116" : 0.4460, "117" : 0.235313,
                 "118" : 0.742935, "119" : 0.26343, "122" : 0.142086, "124" : 0.177588}
Lee_iso_denom = "120"

#*** Import composition of Sn-Std (intern. norm) or unspiked measurement

meas_dict = {"112" : 0.029825, "114" : 0.020191, "115" : 0.010363, "116" : 0.4460, "117" : 0.235303,
             "118" : 0.742926, "119" : 0.263445, "122" : 0.142079, "124" : 0.177547}
meas_denom = "120"

# Spike - Isotope Abundances
#*** Import composition of Sn Spike ***#

spike_dict_117_122 = {"112" : 0.0003468, "114" : 0.0003468, "115": 0.0003958,"116" : 0.01300079, "117": 0.43863839,
              "118": 0.02698559, "119" : 0.01077660, "120" : 0.030053, "122" : 0.4700396, "124" : 0.00941660}

spike_117 = {"112" : 0.0005, "114" : 0.0005, "115" : 0.0006,"116" : 0.0231, "117": 0.8911,
              "118": 0.0452, "119" : 0.0112, "120" : 0.0218, "122" : 0.0029, "124" : 0.0021}

spike_122 = {"112" : 0.0002, "114" : 0.0002, "115" : 0.0002,"116" : 0.00033, "117": 0.0036,
              "118": 0.0095, "119" : 0.0104, "120" : 0.0381, "122" : 0.9219, "124" : 0.0165}

# Data - measured mixture
path = "/Users/marf/Desktop/PhD Temp/Double Spike/Planning Sn 117-122/real data -intern norm/"
df = pd.read_csv(str(path + "2015_03_25_1367_raw_ratio_corr.csv"))
df_new = df.ix[:,"112":"124"]
df_new = df_new.drop(["117_2", "118_2", "119_2", "122_2"], axis=1)

# Masses object used for calc
Sn_mass_obj = load_mass_dict(Sn_masses)

# Abundance object used for calc
Sn_Lee_obj = load_ratio_dict(Sn_ratios_Lee, Lee_iso_denom)
Sn_meas_obj = load_ratio_dict(meas_dict, meas_denom)
spike_obj = load_abundance_dict(spike_dict_117_122)

