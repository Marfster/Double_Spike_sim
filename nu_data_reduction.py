__author__ = 'Matthias Friebel'

import pandas as pd
import re
from collections import Counter
from math import log
import numpy as np
import math



class NU_data_read(object):

    def __init__(self, path, datafile_list, cup_configuration):

        '''******************************************************************************
           Values for Reading in the data and storing it accordingly to the cup settings
           ****************************************************************************** '''
        self.path = path #Path of Datafiles
        self.files = datafile_list  # List of datafile numbers
        #self.cycles = no_cycles # isotope measurement lines
        self.cups = cup_configuration # "dictionary with cup configuration
        self.zero_cycles = [k for k in cup_configuration.keys() if 'zero' in k] # cycles for one measurements
    '''***********************************
       Functions used for Reading the data
       ***********************************'''

    # Extracts sample name, date, ... of NU csv file
    def extract_metadata(self, filenumber, metadata_string):
        datafile = "Data_" + str(filenumber) + ".csv"

        df = pd.read_table(self.path + datafile, dtype=str, header=None)
        for row in range(len(df.index)):
            line = str(df[0][row])
            if metadata_string in line:
                metadata = re.search(r""+metadata_string+",(.*)", line).group(1)
        return metadata

    # Reads in the data from NU csv files
    def data_read(self, filenumber):
        datafile = "Data_" + str(filenumber) + ".csv"


        data = pd.read_csv(self.path + datafile, skiprows=58, index_col='Cycle') # before line 58 only headers
        return data

    # Return the measured signals in the cycle of each sample in form of dictionary
    def data_signals(self, filenumber):
        data = self.data_read(filenumber)
        columns = data.columns.values.tolist()
        dict_signals = {}
        for cycle in self.cups: # reads through cup configuration and writes signals per cup into a dictionary
            dict_signals[cycle] = {}
            for cup in self.cups[cycle]:
                if (cup in columns):
                    dict_signals[cycle][cup] = data[cup].to_dict()
        return dict_signals

    # Applies the Zero correction of each sample (csv file) - Signal of measurement cycle - zero cycle
    def data_zero_corr(self, filenumber):
        data = self.data_signals(filenumber)
        # Outlier correction 2SD
        #cyclex_signal = cyclex_signal[np.abs(cyclex_signal-cyclex_signal.mean())<=(2*cyclex_signal.std())]
        # Zero correction
        zero_signals = {}
        counter = 0
        for cycle in self.cups:
            if (cycle not in self.zero_cycles):
                zero_signals[cycle] = {}
                if (len(self.zero_cycles) == 1):
                    for cup in self.cups[cycle]:
                        match = [k for k in self.cups[self.zero_cycles[counter]].keys() if cup[:-4] in k] #screen throug zero1 dict and match cycle cups
                        if (len(match) > 0):
                            cup_cycle = Counter(data[cycle][cup])
                            cup_zero = Counter(data[self.zero_cycles[counter]][match[0]])
                            cup_cycle.subtract(cup_zero)
                            zero_signals[cycle][cup] = cup_cycle
                else:
                    counter += 1
                    for cup in self.cups[cycle]:
                        match = [k for k in self.cups[self.zero_cycles[counter]].keys() if cup[:-4] in k] #screen throug zero1 dict and match cycle cups
                        if (len(match) > 0):
                            cup_cycle = Counter(data[cycle][cup])
                            cup_zero = Counter(data[self.zero_cycles[counter]][match[0]])
                            cup_cycle.subtract(cup_zero)
                            zero_signals[cycle][cup] = cup_cycle
        return zero_signals


    # Background correction for one sample with one background - not used
    def data_bgd_corr(self, filenumber, filenumber_bgd):
        zero_data = self.data_zero_corr(filenumber)
        zero_data_bgd = self.data_zero_corr(filenumber_bgd)
        bgd_signals = {}
        for cycle in zero_data:
            bgd_signals[cycle] = {}
            for cup in zero_data[cycle]:
                if (cup in zero_data_bgd[cycle]):
                    cup_cycle = Counter(zero_data[cycle][cup])
                    cup_cycle_bgd = Counter(zero_data_bgd[cycle][cup])
                    cup_cycle.subtract(cup_cycle_bgd)
                    bgd_signals[cycle][cup] = cup_cycle
        return bgd_signals

class Neptune_data_read(object):

    def __init__(self,path, datafile_list, cup_configuration):

        '''******************************************************************************
        Values for Reading in the data and storing it accordingly to the cup settings
        ****************************************************************************** '''
        self.path = path  # Path of Datafiles
        self.files = datafile_list  # List of datafile numbers
        #self.cycles = no_cycles # isotope measurement lines
        self.cups = cup_configuration # "dictionary with cup configuration
        '''***********************************
        Functions used for Reading the data
        ***********************************'''

    # Extracts sample name
    def extract_metadata(self, filex, metadata_string):

        df = pd.read_table(filex, dtype=str, header=None)
        for row in range(len(df.index)):
            line = str(df[0][row])
            if metadata_string in line:
                metadata = re.search(r"" + metadata_string + ":(.*)", line).group(1)
        return metadata

    # Reads in the data
    def data_read(self, filex):
        data = pd.read_table(filex, skiprows=16, index_col = "Cycle")
        data = data.ix[:,:"124Sn"]
        data = data.drop(["***"],axis=0)
        data.ix[:,"117Sn":] = data.ix[:,"117Sn":].astype(np.float64)
        data = data.dropna(axis=0)
        data.index = data.index.astype(np.int64)
        return data

    # Return the measured signals in the cycle of each sample
    def data_signals(self, filex):
        data = self.data_read(filex)
        columns = data.columns.values.tolist()
        dict_signals = {}
        for cycle in self.cups:
            dict_signals[cycle] = {}
            for cup in self.cups[cycle]:
                if (cup in columns):
                    dict_signals[cycle][cup] = data[cup].to_dict()
        return dict_signals


    # Background correction for one sample with one background
    def data_bgd_corr(self, filenumber, filenumber_bgd):
        zero_data = self.data_zero_corr(filenumber)
        zero_data_bgd = self.data_zero_corr(filenumber_bgd)
        bgd_signals = {}
        for cycle in zero_data:
            bgd_signals[cycle] = {}
            for cup in zero_data[cycle]:
                if (cup in zero_data_bgd[cycle]):
                    cup_cycle = Counter(zero_data[cycle][cup])
                    cup_cycle_bgd = Counter(zero_data_bgd[cycle][cup])
                    cup_cycle.subtract(cup_cycle_bgd)
                    bgd_signals[cycle][cup] = cup_cycle
        return bgd_signals

class normalisation(object):

    def __init__(self, data_dict, cycle_no, cup_configuration, database, mass_range, isotopes_for_corr, denom_corr_ratio=None, law="exp", n_GPL=None):

        self.data_dict = data_dict                                                      #dictionary with data (usually zero and/or bgd corrected)
        self.cycle_no = cycle_no                                                        # number of cycles per measurement
        self.cups = cup_configuration                                                   # cup configuration used for measurement (see "sn_config.py")
        self.zero_cycles = [k for k in self.cups.keys() if 'zero' in k]                 # number of cycles per zero measurement
        self.database = database                                                        # refers to database in "sn_config.py", which contains the used Isotope masses and ratios
        self.mass_range = mass_range                                                    # refers to defined mass range of isotopes for calculation (see "sn_config.py")
        self.isotopes_for_corr = isotopes_for_corr                                      # defines Isotopes used for Interferences correction (see "sn_config.py or ipython notebook")
        self.graph_of_corr = self.mass_range.get_graph_of_corr(self.isotopes_for_corr)  # creates a graph of dependencies for topologic sorting
        self.order_of_corr = self.mass_range.get_order_of_corr(self.isotopes_for_corr)  # creates a order for topologic sorting
        self.denom_corr_ratio = denom_corr_ratio                                        # alternative (interference-free) Isotope used for Interference correction on denominator isotope
        self.law_mass_frac = law
        self.n_GPL = n_GPL
    '''*********************************************************************
       Functions used for internal normalisation and interference correction
       *********************************************************************'''

    # Lookup Isotope Signal
    def lookup_signal(self, isotope, isotope_from_line1 = True):
        invers_dict = {}
        cup_list = []
        for line in self.data_dict:
            invers_dict[line] = dict(zip(self.cups[line].values(), self.cups[line]))
            if (isotope in invers_dict[line]):
                cup_list.append((line, invers_dict[line][isotope]))
        if isotope_from_line1 == False and len(cup_list) > 1:
            return self.data_dict[cup_list[1][0]][cup_list[1][1]][self.cycle_no]
        else:
            return self.data_dict[cup_list[0][0]][cup_list[0][1]][self.cycle_no]

    # Raw Isotope Ratios
    def raw_ratio(self, isotope_nom, isotope_denom, isotope_from_line1 = True):
        return self.lookup_signal(isotope_nom, isotope_from_line1)/self.lookup_signal(isotope_denom, isotope_from_line1)

    # Mass fractionation correction - using exponential law
    def mass_frac_exp_law(self, element_nom, element_denom, isotope_nom, isotope_denom, beta, isotope_from_line1 = True, corr_isotope_denom = None):
        signal_isotope_nom = self.lookup_signal(isotope_nom, isotope_from_line1)
        if corr_isotope_denom == None:
            signal_isotope_denom = self.lookup_signal(isotope_denom, isotope_from_line1 = True)
        else:
            signal_isotope_denom = corr_isotope_denom
        mass_isotope_nom = self.database[element_nom]["Masses"].get_Isotope_mass(isotope_nom)
        mass_isotope_denom = self.database[element_denom]["Masses"].get_Isotope_mass(isotope_denom)

        return (signal_isotope_nom/signal_isotope_denom) * (mass_isotope_nom/mass_isotope_denom)**beta

    # beta - using exponential law
    def beta_exp_law(self, true_ratio, ratio_raw, mass_isotope_nom, mass_isotope_denom):
        return log(true_ratio / (ratio_raw)) / log(mass_isotope_nom/mass_isotope_denom) #using natural logarithm

    # Mass fractionation correction - using exponential law
    def mass_frac_GPL(self, element_nom, element_denom, isotope_nom, isotope_denom, beta, isotope_from_line1 = True, corr_isotope_denom = None):
        signal_isotope_nom = self.lookup_signal(isotope_nom, isotope_from_line1)
        if corr_isotope_denom == None:
            signal_isotope_denom = self.lookup_signal(isotope_denom, isotope_from_line1 = True)
        else:
            signal_isotope_denom = corr_isotope_denom
        mass_isotope_nom = self.database[element_nom]["Masses"].get_Isotope_mass(isotope_nom)
        mass_isotope_denom = self.database[element_denom]["Masses"].get_Isotope_mass(isotope_denom)

        return (signal_isotope_nom/signal_isotope_denom) * beta**(mass_isotope_nom ** self.n_GPL - mass_isotope_denom ** self.n_GPL)

    # beta - using GPL
    def beta_GPL(self, true_ratio, ratio_raw, mass_isotope_nom, mass_isotope_denom):
        return (true_ratio / ratio_raw) ** (1/ (mass_isotope_nom ** self.n_GPL - mass_isotope_denom ** self.n_GPL))

    # calculate - n for GPL
    def n_inf_GPL(self, element_nom, element_denom, isotope_nom, isotope_denom, norm_ratio, isotope_from_line1 = True, corr_isotope_denom = None):
        signal_isotope_nom = self.lookup_signal(isotope_nom, isotope_from_line1)
        signal_isotope_ref = self.lookup_signal(norm_ratio[0], isotope_from_line1)
        if corr_isotope_denom == None:
            signal_isotope_denom = self.lookup_signal(isotope_denom, isotope_from_line1=True)
        else:
            signal_isotope_denom = corr_isotope_denom

        mass_isotope_nom = self.database[element_nom]["Masses"].get_Isotope_mass(isotope_nom)
        mass_isotope_denom = self.database[element_denom]["Masses"].get_Isotope_mass(isotope_denom)

        mass_isotope_nom_ref = self.database[element_nom]["Masses"].get_Isotope_mass(norm_ratio[0])

        true_ratio = self.database[element_nom]["Ratios"].get_ratio(isotope_nom, isotope_denom)
        true_ratio_ref = self.database[element_nom]["Ratios"].get_ratio(norm_ratio[0], isotope_denom)

        beta_kin = log(mass_isotope_denom/mass_isotope_nom)/log(mass_isotope_denom/mass_isotope_nom_ref)
        beta_eq = (1/mass_isotope_denom - 1/mass_isotope_nom)/(1/mass_isotope_denom - 1/mass_isotope_nom_ref)
        beta_inf = log((signal_isotope_nom/signal_isotope_denom)/true_ratio) / log((signal_isotope_ref/signal_isotope_denom)/true_ratio_ref)
        return (beta_inf - beta_kin) / (beta_kin - beta_eq)

    # chooses formula based on defined mass-fractionation law
    def mass_frac_law(self, element_nom, element_denom, isotope, isotope_denom, beta, isotope_from_line1, corr_isotope_denom):
        if self.law_mass_frac == "exp":
            return self.mass_frac_exp_law(element_nom, element_denom, isotope, isotope_denom, beta, isotope_from_line1, corr_isotope_denom)
        elif self.law_mass_frac == "GPL":
            return self.mass_frac_GPL(element_nom, element_denom, isotope, isotope_denom, beta, isotope_from_line1, corr_isotope_denom)
        else:
            print 'wrong input for mass-fractionation law, use "exp" or "GPL"'

    # chooses formula (invers caluculation to get raw ratios) based on defined mass-fractionation law
    def mass_frac_law_inv(self, isotope_ratio, mass_isotope_nom, mass_isotope_denom, beta):
        if self.law_mass_frac == "exp":
            return isotope_ratio * (mass_isotope_nom/mass_isotope_denom)**-beta
        elif self.law_mass_frac == "GPL":
            return isotope_ratio * beta**-(mass_isotope_nom ** self.n_GPL - mass_isotope_denom ** self.n_GPL)
        else:
            print 'wrong input for mass-fractionation law, use "exp" or "GPL"'

    # chooses formula based on defined mass-fractionation law
    def beta_law(self, true_ratio, ratio_raw, mass_isotope_nom, mass_isotope_denom):
        if self.law_mass_frac == "exp":
            return self.beta_exp_law(true_ratio, ratio_raw, mass_isotope_nom, mass_isotope_denom)

        elif self.law_mass_frac == "GPL":
            return self.beta_GPL(true_ratio, ratio_raw, mass_isotope_nom, mass_isotope_denom)

    # Interference correction for interference on the nominator isotopes
    def interference_corr_all(self, element_denom, isotope_denom, beta, corr_isotope_denom = None, isotope_from_line1 = True):
        corr_dict = {}
        counter = 0
        for level in self.order_of_corr:
            if counter == 0: # isotope masses with no isobaric interference - only internal normalisation
                for isotope in level:
                    element_nom = self.mass_range.get_isotopes(isotope)[0]
                    corr_dict[isotope] = self.mass_frac_law(element_nom, element_denom, isotope, isotope_denom, beta, isotope_from_line1, corr_isotope_denom)

                counter += 1
            elif counter > 0:
                for isotope in level: # isotope masses with one or more interferences - search for interference element on isotope mass and the corresponding isotope mass used interference correction
                    corr_isotopes = list(self.graph_of_corr[isotope])


                    for corr_isotope in corr_isotopes: # search for element of isobaric interference (e.g. "Te" on 124Sn) &
                        element_nom = self.mass_range.get_isotopes(isotope)
                        element_nom = list(element_nom)[0]
                        element_corr = set(self.mass_range.get_isotopes(corr_isotope))
                        element_corr = element_corr.intersection(self.mass_range.get_isotopes(isotope))
                        element_corr = list(element_corr)

                        if isotope in corr_dict:
                            None
                        else:
                            corr_dict[isotope] = self.mass_frac_law(element_nom, element_denom, isotope, isotope_denom, beta, isotope_from_line1, corr_isotope_denom)

                        true_ratio_corr = self.database[element_corr[0]]["Ratios"].get_ratio(isotope, corr_isotope) # get true ratio (corr_isotope/isotope_nominator - e.g 125Te/124Sn) from database used for correction
                        corr_dict[isotope] -= corr_dict[corr_isotope] * true_ratio_corr # substract interfering ratio of ratio of interest

                counter += 1

        return corr_dict

    # Interference correction allows to correct also interference on denominator isotope
    def interference_corr_ratio(self, element, isotope_nom, isotope_denom, beta, isotope_denom_corr = True, isotope_from_line1 = True):

        if (isotope_nom in self.order_of_corr[0] == True) and (isotope_denom in self.order_of_corr[0] == True): # interference correction for case that nominator and denominator isotope is free of interferences
            isotope_ratio = self.interference_corr_all(element, isotope_denom, beta, isotope_from_line1 = isotope_from_line1)[isotope_nom]

        elif (isotope_denom in self.order_of_corr[0]) == False and (isotope_denom_corr == True): # interference correction for case that nominator and denominator isotope contain interferences and interference on denominator isotope is corrected by using an interference free isotope as new denominator isotope
            new_isotope_denom = self.denom_corr_ratio['isotope_denom']
            isotope_ratio = self.interference_corr_all(element, new_isotope_denom, beta, isotope_from_line1 = isotope_from_line1)[isotope_denom] # eg(Sn120/Sn119)
            mass_isotope_denom = self.database[element]["Masses"].get_Isotope_mass(isotope_denom)
            mass_new_isotope_denom = self.database[element]["Masses"].get_Isotope_mass(new_isotope_denom)
            isotope_ratio_raw = isotope_ratio * (mass_isotope_denom/mass_new_isotope_denom)**(-beta)
            isotope_signal = isotope_ratio_raw * self.lookup_signal(new_isotope_denom, isotope_from_line1 = isotope_from_line1)
            isotope_ratio = self.interference_corr_all(element, isotope_denom, beta, isotope_signal, isotope_from_line1 = isotope_from_line1)[isotope_nom]

        elif (isotope_denom in self.order_of_corr[0]) == False:
            isotope_ratio = self.interference_corr_all(element, isotope_denom, beta, isotope_from_line1 = isotope_from_line1)[isotope_nom]

        else:
            isotope_ratio = self.interference_corr_all(element, isotope_denom, beta, isotope_from_line1 = isotope_from_line1)[isotope_nom]

        return isotope_ratio

    # beta calculation by applying interference correction if necessary
    def beta(self, iterations, element, isotope_nom, isotope_denom, isotope_denom_corr = True, isotope_from_line1 = True):
        # beta used in exponential law for normalisation (uses true_ratio from database)

        if not list(self.graph_of_corr[isotope_nom]):
            # no interference correction necessary to determine beta for elements and isotopes masses used for normalisation (e.g 116Sn/119Sn or 123Sb/121Sb)

            signal_isotope_nom = self.lookup_signal(isotope_nom, isotope_from_line1) # Looks up signal of nominator isotope
            signal_isotope_denom = self.lookup_signal(isotope_denom, isotope_from_line1) # Looks up signal for denominator isotope
            true_ratio = self.database[element]["Ratios"].get_ratio(isotope_nom, isotope_denom) # Looks up true ratio of nominator/denominator from database
            mass_isotope_nom = self.database[element]["Masses"].get_Isotope_mass(isotope_nom) # Looks up mass of nominator isotope
            mass_isotope_denom = self.database[element]["Masses"].get_Isotope_mass(isotope_denom) # Looks up mass of denominator isotope


            raw_ratio = signal_isotope_nom/signal_isotope_denom
            beta_temp = self.beta_law(true_ratio, raw_ratio, mass_isotope_nom, mass_isotope_denom)
            #beta_temp = log(true_ratio / (signal_isotope_nom/ signal_isotope_denom)) / log(mass_isotope_nom/mass_isotope_denom) # calculates beta using exponential law

        else:

            if isotope_denom_corr == False:
                # interference correction necessary to determine beta for elements and isotopes masses used for normalisation, but only on nominator isotopes (e.g. 116Sn/120Sn - 120Sn not interference corrected before detemine beta)

                signal_isotope_nom = self.lookup_signal(isotope_nom, isotope_from_line1) # Looks up signal of nominator isotope
                signal_isotope_denom = self.lookup_signal(isotope_denom, isotope_from_line1) # Looks up signal of denominator isotope
                true_ratio = self.database[element]["Ratios"].get_ratio(isotope_nom, isotope_denom) # Looks up true ratio of nominator/denominator from database
                mass_isotope_nom = self.database[element]["Masses"].get_Isotope_mass(isotope_nom) # Looks up mass of nominator isotope
                mass_isotope_denom = self.database[element]["Masses"].get_Isotope_mass(isotope_denom) # Looks up mass of nominator isotope

                corr_isotopes = list(self.graph_of_corr[isotope_nom]) # Looks up correction isotopes to correct isobaric interferences - e.g for 116Sn -- 111Cd

                for i in corr_isotopes: # Iterates over corr_isotopes list e.g. ["111", ..] - "111" then "..."
                    element_corr = set(self.mass_range.get_isotopes(i)) # Check which elements are on mass "111" - "Cd"
                    element_corr = element_corr.intersection(self.mass_range.get_isotopes(isotope_nom))  # Select corr element based on intersection of corr_element and elements for isotope_nom (e.g. "111" - "Cd", "116" - "Sn,Cd" --> "Cd")
                    element_corr = list(element_corr)[0] # element_corr = "Te" or "Xe" ...
                    true_ratio_corr = self.database[element_corr]["Ratios"].get_ratio(isotope_nom, i) # Looks up true ratio of nominator_isotope/corr_isotope for element_corr (e.g. "116"/"111" for Cd)
                    signal_corr_isotope = self.lookup_signal(i, isotope_from_line1) # Looks up signal of corr_isotope - e.g. signal on  "111Cd"
                    signal_isotope_nom -= (signal_corr_isotope * true_ratio_corr) # subtract e.g "111Cd" * true_ratio("116Cd"/"111Cd")

                raw_ratio = signal_isotope_nom / signal_isotope_denom #assumed raw ratio after interference correction
                beta_temp = self.beta_law(true_ratio, raw_ratio, mass_isotope_nom, mass_isotope_denom) #calculate first beta

                # iteratively solve beta
                for i in range(iterations):
                    isotope_ratio = self.interference_corr_all(element, isotope_denom, beta_temp)[isotope_nom] # interference correction

                    if isotope_ratio > 0:
                        raw_ratio = self.mass_frac_law_inv(isotope_ratio, mass_isotope_nom, mass_isotope_denom, beta_temp) # calculate raw_ratio from interference corr ratio

                        beta_update = self.beta_law(true_ratio, raw_ratio, mass_isotope_nom, mass_isotope_denom) # use raw_ratio to calculate a refined beta value

                        convergence = (beta_temp - beta_update) / beta_temp
                        beta_temp = beta_update # update new beta value for start of new loop
                    else:
                        # jump out of loop if a beta value creates negative isotope values and use beta of previous iteration step
                        print "Iterative Beta Correction Failed! --> beta_temp before iteration used"
                        #beta_temp = log(true_ratio / (signal_isotope_nom/ signal_isotope_denom)) / log(mass_isotope_nom/mass_isotope_denom)
                        beta_temp = self.beta_law(true_ratio, raw_ratio, mass_isotope_nom, mass_isotope_denom)
                        break

            elif isotope_denom_corr == True:
                # nominator and denominator isotopes are interference corrected (e.g 116Sn/120Sn - beta determined by using 117Sn/119Sn for interference corr of 120Sn)

                mass_isotope_nom = self.database[element]["Masses"].get_Isotope_mass(isotope_nom) #e.g. M(116Sn)
                mass_isotope_denom = self.database[element]["Masses"].get_Isotope_mass(isotope_denom)# e.g. M(120Sn)
                true_ratio = self.database[element]["Ratios"].get_ratio(isotope_nom, isotope_denom) # e.g 116Sn/120Sn

                new_isotope_nom = self.denom_corr_ratio['isotope_nom'] #e.g. 117Sn
                signal_new_isotope_nom = self.lookup_signal(new_isotope_nom, isotope_from_line1) #e.g. Signal on 117Sn
                mass_new_isotope_nom = self.database[element]["Masses"].get_Isotope_mass(new_isotope_nom) #e.g. M(117Sn)

                new_isotope_denom = self.denom_corr_ratio['isotope_denom'] #e.g. 119Sn
                signal_new_isotope_denom = self.lookup_signal(new_isotope_denom, isotope_from_line1) #e.g. Signal on 119Sn
                mass_new_isotope_denom = self.database[element]["Masses"].get_Isotope_mass(new_isotope_denom) #e.g. M(119Sn)

                true_ratio_new = self.database[element]["Ratios"].get_ratio(new_isotope_nom, new_isotope_denom) # e.g 117Sn/119Sn

                raw_ratio_new = signal_new_isotope_nom / signal_new_isotope_denom #e.g. raw ratio 117Sn/119Sn
                beta_temp = self.beta_law(true_ratio_new, raw_ratio_new, mass_new_isotope_nom, mass_new_isotope_denom) # determine initial beta for interference corr

                # iteratively solve beta
                for i in range(iterations):
                    new_isotope_ratio_nom = self.interference_corr_all(element, new_isotope_denom, beta_temp)[isotope_nom]  #e.g. 116Sn/119Sn - interference correction
                    new_isotope_ratio_denom = self.interference_corr_all(element, new_isotope_denom, beta_temp)[isotope_denom] # e.g. 120Sn/119Sn - interferences correction

                    new_isotope_ratio_nom_raw = self.mass_frac_law_inv(new_isotope_ratio_nom, mass_isotope_nom, mass_new_isotope_denom, beta_temp) # determine raw ratio of e.g. 116Sn/119Sn - interference corrected
                    new_isotope_ratio_denom_raw = self.mass_frac_law_inv(new_isotope_ratio_denom, mass_isotope_denom, mass_new_isotope_denom, beta_temp) # determine raw ratio of e.g. 120Sn/119Sn - interference corrected

                    raw_ratio = new_isotope_ratio_nom_raw / new_isotope_ratio_denom_raw # raw ratio of e.g. 116Sn/120Sn by 116Sn/119Sn / 120Sn/119Sn

                    beta_update = self.beta_law(true_ratio, raw_ratio, mass_isotope_nom, mass_isotope_denom) #calculate new beta value using interference corrected raw ratio
                    convergence = (beta_temp - beta_update) / beta_temp
                    beta_temp = beta_update # update beta_temp for new iteration step

        return beta_temp

    # beta calculation by applying interference correction if necessary, but with a given true_ratio to use
    def beta_true_change(self, iterations, element, isotope_nom, isotope_denom, true_ratio, isotope_denom_corr = False, isotope_from_line1 = True):
        # beta used in exponential law but with self-defined true_ratio for normalisation

        if not list(self.graph_of_corr[isotope_nom]):
            # no interference correction necessary to determine beta for elements and isotopes masses used for normalisation (e.g 116Sn/119Sn or 123Sb/121Sb)

            signal_isotope_nom = self.lookup_signal(isotope_nom, isotope_from_line1) # Looks up signal of nominator isotope
            signal_isotope_denom = self.lookup_signal(isotope_denom, isotope_from_line1) # Looks up signal for denominator isotope
            true_ratio = true_ratio # Choosen true ratio of nominator/denominator from database
            mass_isotope_nom = self.database[element]["Masses"].get_Isotope_mass(isotope_nom) # Looks up mass of nominator isotope
            mass_isotope_denom = self.database[element]["Masses"].get_Isotope_mass(isotope_denom) # Looks up mass of denominator isotope


            raw_ratio = signal_isotope_nom/signal_isotope_denom
            beta_temp = self.beta_law(true_ratio, raw_ratio, mass_isotope_nom, mass_isotope_denom)
            #beta_temp = log(true_ratio / (signal_isotope_nom/ signal_isotope_denom)) / log(mass_isotope_nom/mass_isotope_denom) # calculates beta using exponential law

        else:

            if isotope_denom_corr == False:
                # interference correction necessary to determine beta for elements and isotopes masses used for normalisation, but only on nominator isotopes (e.g. 116Sn/120Sn - 120Sn not interference corrected before detemine beta)

                signal_isotope_nom = self.lookup_signal(isotope_nom, isotope_from_line1) # Looks up signal of nominator isotope
                signal_isotope_denom = self.lookup_signal(isotope_denom, isotope_from_line1) # Looks up signal of denominator isotope
                true_ratio = true_ratio # Choosen true ratio of nominator/denominator from database
                mass_isotope_nom = self.database[element]["Masses"].get_Isotope_mass(isotope_nom) # Looks up mass of nominator isotope
                mass_isotope_denom = self.database[element]["Masses"].get_Isotope_mass(isotope_denom) # Looks up mass of nominator isotope

                corr_isotopes = list(self.graph_of_corr[isotope_nom]) # Looks up correction isotopes to correct isobaric interferences - e.g for 116Sn -- 111Cd

                for i in corr_isotopes: # Iterates over corr_isotopes list e.g. ["111", ..] - "111" then "..."
                    element_corr = set(self.mass_range.get_isotopes(i)) # Check which elements are on mass "111" - "Cd"
                    element_corr = element_corr.intersection(self.mass_range.get_isotopes(isotope_nom))  # Select corr element based on intersection of corr_element and elements for isotope_nom (e.g. "111" - "Cd", "116" - "Sn,Cd" --> "Cd")
                    element_corr = list(element_corr)[0] # element_corr = "Te" or "Xe" ...
                    true_ratio_corr = self.database[element_corr]["Ratios"].get_ratio(isotope_nom, i) # Looks up true ratio of nominator_isotope/corr_isotope for element_corr (e.g. "116"/"111" for Cd)
                    signal_corr_isotope = self.lookup_signal(i, isotope_from_line1) # Looks up signal of corr_isotope - e.g. signal on  "111Cd"
                    signal_isotope_nom -= (signal_corr_isotope * true_ratio_corr) # subtract e.g "111Cd" * true_ratio("116Cd"/"111Cd")

                raw_ratio = signal_isotope_nom / signal_isotope_denom #assumed raw ratio after interference correction
                beta_temp = self.beta_law(true_ratio, raw_ratio, mass_isotope_nom, mass_isotope_denom) #calculate first beta

                # iteratively solve beta
                for i in range(iterations):
                    isotope_ratio = self.interference_corr_all(element, isotope_denom, beta_temp)[isotope_nom] # interference correction

                    if isotope_ratio > 0:
                        raw_ratio = self.mass_frac_law_inv(isotope_ratio, mass_isotope_nom, mass_isotope_denom, beta_temp) # calculate raw_ratio from interference corr ratio

                        beta_update = self.beta_law(true_ratio, raw_ratio, mass_isotope_nom, mass_isotope_denom) # use raw_ratio to calculate a refined beta value

                        convergence = (beta_temp - beta_update) / beta_temp
                        beta_temp = beta_update # update new beta value for start of new loop
                    else:
                        # jump out of loop if a beta value creates negative isotope values and use beta of previous iteration step
                        print "Iterative Beta Correction Failed! --> beta_temp before iteration used"
                        #beta_temp = log(true_ratio / (signal_isotope_nom/ signal_isotope_denom)) / log(mass_isotope_nom/mass_isotope_denom)
                        beta_temp = self.beta_law(true_ratio, raw_ratio, mass_isotope_nom, mass_isotope_denom)
                        break

            elif isotope_denom_corr == True:
                # nominator and denominator isotopes are interference corrected (e.g 116Sn/120Sn - beta determined by using 117Sn/119Sn for interference corr of 120Sn)

                mass_isotope_nom = self.database[element]["Masses"].get_Isotope_mass(isotope_nom) #e.g. M(116Sn)
                mass_isotope_denom = self.database[element]["Masses"].get_Isotope_mass(isotope_denom)# e.g. M(120Sn)
                true_ratio = true_ratio # Choosen true ratio of nominator/denominator from database

                new_isotope_nom = self.denom_corr_ratio['isotope_nom'] #e.g. 117Sn
                signal_new_isotope_nom = self.lookup_signal(new_isotope_nom, isotope_from_line1) #e.g. Signal on 117Sn
                mass_new_isotope_nom = self.database[element]["Masses"].get_Isotope_mass(new_isotope_nom) #e.g. M(117Sn)

                new_isotope_denom = self.denom_corr_ratio['isotope_denom'] #e.g. 119Sn
                signal_new_isotope_denom = self.lookup_signal(new_isotope_denom, isotope_from_line1) #e.g. Signal on 119Sn
                mass_new_isotope_denom = self.database[element]["Masses"].get_Isotope_mass(new_isotope_denom) #e.g. M(119Sn)

                true_ratio_new = self.database[element]["Ratios"].get_ratio(new_isotope_nom, new_isotope_denom) # e.g 117Sn/119Sn

                raw_ratio_new = signal_new_isotope_nom / signal_new_isotope_denom #e.g. raw ratio 117Sn/119Sn
                beta_temp = self.beta_law(true_ratio_new, raw_ratio_new, mass_new_isotope_nom, mass_new_isotope_denom) # determine initial beta for interference corr

                # iteratively solve beta
                for i in range(iterations):
                    new_isotope_ratio_nom = self.interference_corr_all(element, new_isotope_denom, beta_temp)[isotope_nom]  #e.g. 116Sn/119Sn - interference correction
                    new_isotope_ratio_denom = self.interference_corr_all(element, new_isotope_denom, beta_temp)[isotope_denom] # e.g. 120Sn/119Sn - interferences correction

                    new_isotope_ratio_nom_raw = self.mass_frac_law_inv(new_isotope_ratio_nom, mass_isotope_nom, mass_new_isotope_denom, beta_temp) # determine raw ratio of e.g. 116Sn/119Sn - interference corrected
                    new_isotope_ratio_denom_raw = self.mass_frac_law_inv(new_isotope_ratio_denom, mass_isotope_denom, mass_new_isotope_denom, beta_temp) # determine raw ratio of e.g. 120Sn/119Sn - interference corrected

                    raw_ratio = new_isotope_ratio_nom_raw / new_isotope_ratio_denom_raw # raw ratio of e.g. 116Sn/120Sn by 116Sn/119Sn / 120Sn/119Sn

                    beta_update = self.beta_law(true_ratio, raw_ratio, mass_isotope_nom, mass_isotope_denom) #calculate new beta value using interference corrected raw ratio
                    convergence = (beta_temp - beta_update) / beta_temp
                    beta_temp = beta_update # update beta_temp for new iteration step

        return beta_temp

    # Interference corr on single raw ratio
    def interference_corr_all_raw(self, element_denom, isotope_denom, corr_isotope_denom = None, isotope_from_line1 = True):
        # Does interference correction with no mass-bias corrected ratios (beta calculation)
        corr_dict = {}
        counter = 0
        for level in self.order_of_corr:
            if counter == 0:
                for isotope in level:
                    corr_dict[isotope] = self.raw_ratio(isotope, isotope_denom, isotope_from_line1 = isotope_from_line1)
                counter += 1
            elif counter > 0:
                for isotope in level:
                    corr_isotopes = list(self.graph_of_corr[isotope])

                    for corr_isotope in corr_isotopes:  # search for element of isobaric interference (e.g. "Te" on 124Sn) &
                        element_nom = self.mass_range.get_isotopes(isotope)
                        element_nom = list(element_nom)[0]
                        element_corr = set(self.mass_range.get_isotopes(corr_isotope))
                        element_corr = element_corr.intersection(self.mass_range.get_isotopes(isotope))
                        element_corr = list(element_corr)

                        if isotope in corr_dict:
                            None
                        else:
                            corr_dict[isotope] = self.raw_ratio(isotope, isotope_denom, isotope_from_line1 = isotope_from_line1)

                        true_ratio_corr = self.database[element_corr[0]]["Ratios"].get_ratio(isotope, corr_isotope)  # get true ratio (corr_isotope/isotope_nominator - e.g 125Te/124Sn) from database used for correction
                        corr_dict[isotope] -= corr_dict[corr_isotope] * true_ratio_corr  # substract interfering ratio of ratio of interest

                counter += 1

        return corr_dict

    # Interference correction on raw ratios
    def interference_corr_ratio_raw(self, element, isotope_nom, isotope_denom, isotope_denom_corr = True, isotope_from_line1 = True):

        if (isotope_nom in self.order_of_corr[0] == True) and (isotope_denom in self.order_of_corr[0] == True):
            isotope_ratio = self.interference_corr_all_raw(element, isotope_denom, isotope_from_line1 = isotope_from_line1)[isotope_nom]

        elif (isotope_denom in self.order_of_corr[0]) == False and (isotope_denom_corr == True):
            new_isotope_denom = self.denom_corr_ratio['isotope_denom']
            isotope_ratio = self.interference_corr_all_raw(element, new_isotope_denom, isotope_from_line1 = isotope_from_line1)[isotope_denom] # eg(Sn120/Sn119)
            isotope_signal = isotope_ratio * self.lookup_signal(new_isotope_denom, isotope_from_line1 = isotope_from_line1)
            isotope_ratio = self.interference_corr_all_raw(element, isotope_denom, isotope_signal, isotope_from_line1 = isotope_from_line1)[isotope_nom]

        elif (isotope_denom in self.order_of_corr[0]) == False:
            isotope_ratio = self.interference_corr_all_raw(element, isotope_denom, isotope_from_line1 = isotope_from_line1)[isotope_nom]

        else:
            isotope_ratio = self.interference_corr_all_raw(element, isotope_denom, isotope_from_line1 = isotope_from_line1)[isotope_nom]

        return isotope_ratio

# class contains all methods to do different evaluation types on one sample measurement with n cycles
class evaluation(object):

    def __init__(self, data_dict, cycles, isotopes ,cup_configuration, database, mass_range, isotopes_for_corr, denom_corr_ratio=None, law="exp", n_GPL=None):
        self.data_dict = data_dict #load data dictionary
        self.cycles = cycles #cycles in
        self.isotopes = isotopes # isotopes used for calculation
        self.cups = cup_configuration #cup configuration of instrument
        self.database = database # database with isotope masses, true ratios, ...
        self.mass_range = mass_range # isotope masses link to elements on that mass
        self.isotopes_for_corr = isotopes_for_corr # isotope masses used for interference correction
        self.denom_corr_ratio = denom_corr_ratio # corr isotope ratio (e.g. 119Sn/117Sn) for interference corr on denominator isotope
        self.law_mass_frac = law # mass-fractionation law used for mass-bias correction
        self.n_GPL = n_GPL # n values used for GPL law

    # return data dictionary stored
    def get_df(self):
        return self.data_dict

    # adjust signal from line 1 to line 2
    def line2_corr(self, df, isotope):
        invers_dict = {}
        if "cycle2" in df:
            for cycle in df:
                invers_dict[cycle] = dict(zip(self.cups[cycle].values(), self.cups[cycle]))
                if cycle == "cycle2":
                    keys = df[cycle].keys()
                    keys.append(keys.pop(keys.index(invers_dict["cycle2"][isotope])))
                    for cup in keys:
                        for meas_point in df[cycle][cup]:
                                if (df["cycle1"][invers_dict["cycle1"][isotope]][meas_point]) and  df["cycle2"][invers_dict["cycle2"][isotope]][meas_point]:
                                    line2_corr = df["cycle1"][invers_dict["cycle1"][isotope]][meas_point]/df["cycle2"][invers_dict["cycle2"][isotope]][meas_point]
                                    df[cycle][cup][meas_point] = df[cycle][cup][meas_point] * line2_corr
            self.data_dict = df
        return self.data_dict

    #background correction
    def data_bgd_corr(self, df_bgd_1, df_bgd_2):
        if df_bgd_2:
            bgd_signals = {}
            for cycle in df_bgd_1:
                bgd_signals[cycle] = {}
                for cup in df_bgd_1[cycle]:
                    if (cup in df_bgd_2[cycle]):
                        bgd_signals[cycle][cup] = {}
                        for meas_point in df_bgd_1[cycle][cup]:
                            if (meas_point in df_bgd_2[cycle][cup]):
                                avg_cup_cycle_bgd = np.nanmean([df_bgd_1[cycle][cup][meas_point], df_bgd_2[cycle][cup][meas_point]])
                                #avg_cup_cycle_bgd = np.divide((np.add(df_bgd_1[cycle][cup][meas_point], df_bgd_2[cycle][cup][meas_point])),2)
                                bgd_signals[cycle][cup][meas_point] = avg_cup_cycle_bgd
        else:
            bgd_signals = df_bgd_1

        df_bgd_corr = {}
        for cycle in self.data_dict:
            df_bgd_corr[cycle] = {}
            for cup in self.data_dict[cycle]:
                if (cup in bgd_signals[cycle]):
                        names = ['id','data']
                        formats = ['float','float']
                        dtype = dict(names = names, formats=formats)
                        cup_cycle = np.array(self.data_dict[cycle][cup].items(), dtype=dtype)
                        cup_cycle_bgd = np.array(bgd_signals[cycle][cup].items(), dtype=dtype)
                        mean = np.nanmean(cup_cycle_bgd["data"])
                        x2 = np.full(len(self.data_dict[cycle][cup]),(mean))
                        cup_cycle = np.array(cup_cycle)
                        cup_cycle["data"] = cup_cycle["data"] - x2
                        df_bgd_corr[cycle][cup] = Counter(dict(enumerate(cup_cycle["data"], 1)))
        self.data_dict = df_bgd_corr
        return self.data_dict

    def data_bgd_corr_2(self, df_bgd_1, df_bgd_2):
        if df_bgd_2:
            bgd_signals = {}
            for cycle in df_bgd_1:
                bgd_signals[cycle] = {}
                for cup in df_bgd_1[cycle]:
                    if (cup in df_bgd_2[cycle]):
                        avg_cup_cycle_bgd_1 = []
                        avg_cup_cycle_bgd_2 = []
                        for meas_point in df_bgd_1[cycle][cup]:
                            if (meas_point in df_bgd_2[cycle][cup]):
                                avg_cup_cycle_bgd_1.append(df_bgd_1[cycle][cup][meas_point])
                                avg_cup_cycle_bgd_2.append(df_bgd_2[cycle][cup][meas_point])

                        bgd_signals[cycle][cup] = np.nanmean([np.nanmean(avg_cup_cycle_bgd_1), np.nanmean(avg_cup_cycle_bgd_2)])

        else:
            bgd_signals = {}
            for cycle in df_bgd_1:
                bgd_signals[cycle] = {}
                for cup in df_bgd_1[cycle]:
                     avg_cup_cycle_bgd_1 = []
                     for meas_point in df_bgd_1[cycle][cup]:
                        avg_cup_cycle_bgd_1.append(df_bgd_1[cycle][cup][meas_point])

                     bgd_signals[cycle][cup] = np.nanmean(avg_cup_cycle_bgd_1)

        df_bgd_corr = {}
        for cycle in self.data_dict:
            df_bgd_corr[cycle] = {}
            for cup in self.data_dict[cycle]:
                if (cup in bgd_signals[cycle]):
                    names = ['id', 'data']
                    formats = ['float', 'float']
                    dtype = dict(names=names, formats=formats)
                    cup_cycle = np.array(self.data_dict[cycle][cup].items(), dtype=dtype)
                    avg_bgd = bgd_signals[cycle][cup]
                    x2 = np.full(len(self.data_dict[cycle][cup]), (avg_bgd))
                    #cup_cycle = np.array(cup_cycle)
                    cup_cycle["data"] = cup_cycle["data"] - x2
                    df_bgd_corr[cycle][cup] = Counter(dict(enumerate(cup_cycle["data"], 1)))
        self.data_dict = df_bgd_corr
        return self.data_dict

    #returns raw signals
    def raw_signals_all(self):
        data_sample = {}

        for n in self.cycles:
            corr = normalisation(self.data_dict, n, self.cups, self.database, self.mass_range, self.isotopes_for_corr, self.denom_corr_ratio)
            data_sample[n] = {}

            for key in self.mass_range.get_mass_range():
                data_sample[n][key] = corr.lookup_signal(key, isotope_from_line1=True)

            if len(self.isotopes) > 1:
                for isotope in self.isotopes[1]:
                    data_sample[n][key + "_2"] = corr.lookup_signal(key, isotope_from_line1=False)

        return data_sample

    # return raw signals including isotope denominator
    def raw_signals(self, isotope_denom):
        data_sample = {}
        isotopes_0 = self.isotopes[0]
        isotopes_0.append(isotope_denom)

        for n in self.cycles:
            corr = normalisation(self.data_dict, n, self.cups, self.database, self.mass_range, self.isotopes_for_corr, self.denom_corr_ratio)
            data_sample[n] = {}

            for isotope in isotopes_0:
                data_sample[n][isotope] = corr.lookup_signal(isotope, isotope_from_line1=True)

            if len(self.isotopes) > 1:
                for isotope in self.isotopes[1]:
                    data_sample[n][isotope + "_2"] = corr.lookup_signal(isotope, isotope_from_line1=False)

        return data_sample

    #raw ratios
    def raw_ratios(self, isotope_denom):
        data_sample = {}
        for n in self.cycles:
            corr = normalisation(self.data_dict, n, self.cups, self.database, self.mass_range, self.isotopes_for_corr, self.denom_corr_ratio)
            data_sample[n] = {}

            for isotope in self.isotopes[0]:
                data_sample[n][isotope] = corr.raw_ratio(isotope, isotope_denom)

            if len(self.isotopes) > 1:
                for isotope in self.isotopes[1]:
                    data_sample[n][isotope + "_2"] = corr.raw_ratio(isotope, isotope_denom)
        return data_sample

    # raw ratios interference corrected
    def raw_ratios_corr(self, isotope_denom):
        data_sample = {}
        for n in self.cycles:
            corr = normalisation(self.data_dict, n, self.cups, self.database, self.mass_range, self.isotopes_for_corr, self.denom_corr_ratio)

            data_sample[n] = {}
            for isotope in self.isotopes[0]:
                data_sample[n][isotope] = corr.interference_corr_ratio_raw("Sn", isotope , isotope_denom, isotope_denom_corr = True, isotope_from_line1 = True)

            if len(self.isotopes) > 1:
                for isotope in self.isotopes[1]:
                    data_sample[n][isotope + "_2"] = corr.interference_corr_ratio_raw("Sn", isotope , isotope_denom, isotope_denom_corr = True, isotope_from_line1 = False)

        return data_sample

    # internal normalisation

    # 1 only mass fractionation
    def mass_fractionation(self, norm_ratio, isotope_denom, iter):
        data_sample = {}
        for n in self.cycles:
            corr = normalisation(self.data_dict, n, self.cups, self.database, self.mass_range, self.isotopes_for_corr, self.denom_corr_ratio, self.law_mass_frac, self.n_GPL)
            beta = corr.beta(iter, "Sn", norm_ratio[0], norm_ratio[1], isotope_denom_corr = False)
            data_sample[n] = {}

            for isotope in self.isotopes[0]:
                data_sample[n][isotope] = corr.mass_frac_law("Sn", "Sn", isotope, isotope_denom, beta)

            if len(self.isotopes) > 1:
                for isotope in self.isotopes[1]:
                    data_sample[n][isotope + "_2"] = corr.mass_frac_law("Sn", "Sn", isotope, isotope_denom, beta, isotope_from_line1 = False)

        return data_sample
    # 2 no corr on isotope_denom,
    def internal_norm_1(self, norm_ratio, isotope_denom, iter):
        data_sample = {}
        for n in self.cycles:
            corr = normalisation(self.data_dict, n, self.cups, self.database, self.mass_range, self.isotopes_for_corr, self.denom_corr_ratio, self.law_mass_frac, self.n_GPL)
            beta = corr.beta(iter, "Sn", norm_ratio[0], norm_ratio[1], isotope_denom_corr = False)

            data_sample[n] = {}
            for isotope in self.isotopes[0]:
                data_sample[n][isotope] = corr.interference_corr_ratio("Sn", isotope , isotope_denom, beta, isotope_denom_corr = False, isotope_from_line1 = True)

            if len(self.isotopes) > 1:
                for isotope in self.isotopes[1]:
                    data_sample[n][isotope + "_2"] = corr.interference_corr_ratio("Sn", isotope , isotope_denom, beta, isotope_denom_corr = False, isotope_from_line1 = False)
        return data_sample


    # 3 corr on isotope_denom
    def internal_norm_2(self, norm_ratio, isotope_denom, iter):
        data_sample = {}
        for n in self.cycles:
            corr = normalisation(self.data_dict, n, self.cups, self.database, self.mass_range, self.isotopes_for_corr, self.denom_corr_ratio, self.law_mass_frac, self.n_GPL)
            beta = corr.beta(iter, "Sn", norm_ratio[0], norm_ratio[1], isotope_denom_corr = True)

            data_sample[n] = {}
            for isotope in self.isotopes[0]:
                data_sample[n][isotope] = corr.interference_corr_ratio("Sn", isotope , isotope_denom, beta, isotope_denom_corr = True, isotope_from_line1 = True)

            if len(self.isotopes) > 1:
                for isotope in self.isotopes[1]:
                    data_sample[n][isotope + "_2"] = corr.interference_corr_ratio("Sn", isotope , isotope_denom, beta, isotope_denom_corr = True, isotope_from_line1 = False)

        return data_sample

    # External normalisation with Sb
    def external_norm_Sb(self, norm_ratio, isotope_denom, iter):
        data_sample = {}
        for n in self.cycles:
            corr = normalisation(self.data_dict, n, self.cups, self.database, self.mass_range, self.isotopes_for_corr, self.denom_corr_ratio, self.law_mass_frac, self.n_GPL)
            beta = corr.beta(iter, "Sb", norm_ratio[0], norm_ratio[1] , isotope_denom_corr = False)

            data_sample[n] = {}
            for isotope in self.isotopes[0]:
                data_sample[n][isotope] = corr.interference_corr_ratio("Sn", isotope , isotope_denom, beta, isotope_denom_corr = True, isotope_from_line1 = True)

            if len(self.isotopes) > 1:
                for isotope in self.isotopes[1]:
                    data_sample[n][isotope + "_2"] = corr.interference_corr_ratio("Sn", isotope , isotope_denom, beta, isotope_denom_corr = True, isotope_from_line1 = False)

        return data_sample

    # normalisation with given beta
    def norm_beta(self, element, isotope_denom, beta):
        data_sample = {}
        for n in self.cycles:
            corr = normalisation(self.data_dict, n, self.cups, self.database, self.mass_range, self.isotopes_for_corr,
                                     self.denom_corr_ratio, self.law_mass_frac, self.n_GPL)

            data_sample[n] = {}
            for isotope in self.isotopes[0]:
                data_sample[n][isotope] = corr.interference_corr_ratio(element, isotope, isotope_denom, beta[n-1],
                                                                       isotope_denom_corr=True,
                                                                       isotope_from_line1=True)

            if len(self.isotopes) > 1:
                for isotope in self.isotopes[1]:
                    data_sample[n][isotope + "_2"] = corr.interference_corr_ratio(element, isotope, isotope_denom,
                                                                                  beta[n-1], isotope_denom_corr=True,
                                                                                  isotope_from_line1=False)
        return data_sample

    # converts normalised isotope ratios into raw_ratios
    def norm_beta_to_raw(self, element, isotope_denom, beta):
        data_sample = {}
        for n in self.cycles:
            corr = normalisation(self.data_dict, n, self.cups, self.database, self.mass_range, self.isotopes_for_corr,
                                     self.denom_corr_ratio, self.law_mass_frac, self.n_GPL)

            data_sample[n] = {}
            for isotope in self.isotopes[0]:
                data_sample[n][isotope] = corr.interference_corr_ratio(element, isotope, isotope_denom, beta[n-1],
                                                                       isotope_denom_corr=True,
                                                                       isotope_from_line1=True)

                mass_isotope_nom = self.database[element]["Masses"].get_Isotope_mass(isotope)
                mass_isotope_denom = self.database[element]["Masses"].get_Isotope_mass(isotope_denom)

                data_sample[n][isotope] = corr.mass_frac_law_inv(data_sample[n][isotope], mass_isotope_nom, mass_isotope_denom, beta[n-1])


            if len(self.isotopes) > 1:
                for isotope in self.isotopes[1]:
                    data_sample[n][isotope + "_2"] = corr.interference_corr_ratio(element, isotope, isotope_denom,
                                                                                  beta[n-1], isotope_denom_corr=True,
                                                                                  isotope_from_line1=False)
                    mass_isotope_nom = self.database[element]["Masses"].get_Isotope_mass(isotope)
                    mass_isotope_denom = self.database[element]["Masses"].get_Isotope_mass(isotope_denom)

                    data_sample[n][isotope + "_2"] = corr.mass_frac_law_inv(data_sample[n][isotope], mass_isotope_nom, mass_isotope_denom, beta[n-1])

        return data_sample

    # Sb normalisation with change true_ratio used for beta_calculation
    def external_norm_Sb_change_true(self, norm_ratio, isotope_denom, iter, changed_true_ratio):
        data_sample = {}
        for n in self.cycles:
            corr = normalisation(self.data_dict, n, self.cups, self.database, self.mass_range, self.isotopes_for_corr, self.denom_corr_ratio, self.law_mass_frac, self.n_GPL)
            beta = corr.beta_true_change(iter, "Sb", norm_ratio[0], norm_ratio[1], changed_true_ratio , isotope_denom_corr = False)

            data_sample[n] = {}
            for isotope in self.isotopes[0]:
                data_sample[n][isotope] = corr.interference_corr_ratio("Sn", isotope , isotope_denom, beta, isotope_denom_corr = True, isotope_from_line1 = True)

            if len(self.isotopes) > 1:
                for isotope in self.isotopes[1]:
                    data_sample[n][isotope + "_2"] = corr.interference_corr_ratio("Sn", isotope , isotope_denom, beta, isotope_denom_corr = True, isotope_from_line1 = False)

        return data_sample

    # normalisation with change true_ratio used for beta_calculation
    def internal_norm_Sb_change_true(self, norm_ratio, isotope_denom, iter, changed_true_ratio):
        data_sample = {}
        for n in self.cycles:
            corr = normalisation(self.data_dict, n, self.cups, self.database, self.mass_range, self.isotopes_for_corr, self.denom_corr_ratio, self.law_mass_frac, self.n_GPL)
            beta = corr.beta_true_change(iter, "Sn", norm_ratio[0], norm_ratio[1], changed_true_ratio , isotope_denom_corr = False)

            data_sample[n] = {}
            for isotope in self.isotopes[0]:
                data_sample[n][isotope] = corr.interference_corr_ratio("Sn", isotope , isotope_denom, beta, isotope_denom_corr = True, isotope_from_line1 = True)

            if len(self.isotopes) > 1:
                for isotope in self.isotopes[1]:
                    data_sample[n][isotope + "_2"] = corr.interference_corr_ratio("Sn", isotope , isotope_denom, beta, isotope_denom_corr = True, isotope_from_line1 = False)

        return data_sample

    # calculate beta and return it
    def beta(self, element, norm_ratio, iter):
        beta = {}
        for n in self.cycles:
            corr = normalisation(self.data_dict, n, self.cups, self.database, self.mass_range, self.isotopes_for_corr,
                                     self.denom_corr_ratio, self.law_mass_frac, self.n_GPL)
            beta[n] = corr.beta(iter, element, norm_ratio[0], norm_ratio[1], isotope_denom_corr=False)

        return beta

    # calculate beta with change true ratio
    def beta_true_change(self, element, norm_ratio, iter, true_ratio):
        beta = {}
        for n in self.cycles:
            corr = normalisation(self.data_dict, n, self.cups, self.database, self.mass_range, self.isotopes_for_corr,
                                     self.denom_corr_ratio, self.law_mass_frac, self.n_GPL)
            beta[n] = corr.beta_true_change(iter, element, norm_ratio[0], norm_ratio[1], true_ratio, isotope_denom_corr=False)

        return beta

    # determines n_inferred to see which law the mass_bias correction should follow
    def n_inf(self, element, norm_ratio, isotope_denom):
        n_inf = {}
        for n in self.cycles:
            corr = normalisation(self.data_dict, n, self.cups, self.database, self.mass_range, self.isotopes_for_corr, self.denom_corr_ratio, self.law_mass_frac, self.n_GPL)
            n_inf[n] = {}

            for isotope in self.isotopes[0]:
                if isotope == norm_ratio[0]:
                    None
                else:
                    n_inf[n][isotope] = corr.n_isotope_flux(element, element, isotope, isotope_denom, norm_ratio)

        return n_inf

    #outlier detection from https://github.com/joferkington/oost_paper_code/blob/master/utilities.py - based on modified z-score
    def mad_based_outlier(self, points, thresh=3.5):
        """
        Returns a boolean array with True if points are outliers and False
        otherwise.

        Parameters:
        -----------
            points : An numobservations by numdimensions array of observations
            thresh : The modified z-score to use as a threshold. Observations with
                a modified z-score (based on the median absolute deviation) greater
                than this value will be classified as outliers.

        Returns:
        --------
            mask : A numobservations-length boolean array.

        References:
        ----------
            Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
            Handle Outliers", The ASQC Basic References in Quality Control:
            Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.
        """
        if len(points.shape) == 1:
            points = points[:,None]
        median = np.median(points, axis=0)
        diff = np.sum((points - median)**2, axis=-1)
        diff = np.sqrt(diff)
        med_abs_deviation = np.median(diff)
     
        modified_z_score = 0.6745 * diff / med_abs_deviation

        return modified_z_score > thresh

    #outlier rejection in dataframe
    def mad_outlier_rejection(self, df):
        data_sample_outlier = pd.DataFrame()
        columns = df.columns.tolist()
        for column in columns:
            data_sample_outlier[column] = df[column].where(~self.mad_based_outlier(df[column]), other=np.NaN)
        return data_sample_outlier

    # outlier rejection in dictionary
    def mad_outlier_rejection_dict(self, dictionary):
        data_sample_outlier = {}
        for cycle in dictionary:
            temp = pd.DataFrame.from_dict(dictionary[cycle], orient = 'columns')
            temp = self.mad_outlier_rejection(temp)
            data_sample_outlier[cycle] = temp.to_dict()
        for cup in data_sample_outlier[cycle]:
            data_sample_outlier[cycle][cup] = Counter(data_sample_outlier[cycle][cup])
        return data_sample_outlier

    # normal z-score outlier rejection
    def z_score_outlier_rejection(self, df):
        df = df[np.abs(df - df.mean())<=(2*df.std())]

        return df


    #data processing wrap up - not working yet
    def data_process(self, path, files, cup_config, isotopes, mass_range, corr_isotopes, denom_corr_ratio, line2_corr, isotope_line2_corr, bgd_corr, option, iter_beta, isotope_denom):
        # Empty dataframes
        data_sample = {}
        sample_names = {}
        avg_ratio_sample_all = {}
        sd2_ratio_sample_all = {}

        if bgd_corr:
            for i in files:
                files[i] = [i, [i-1, i+1]]


        for sample in files:
            df = NU_data_read(path, sample, cup_config)
            cycles = range(1, len(df.data_read(sample).index)+1)

            if bgd_corr == True:
                df_zero = df.data_zero_corr(sample)
                blk_1 = NU_data_read(path, files[sample][1][0], cup_config)
                blk_2 = NU_data_read(path, files[sample][1][1], cup_config)
                df_bgd_1 = blk_1.data_zero_corr(files[sample][1][0])
                df_bgd_2 = blk_2.data_zero_corr(files[sample][1][1])
                new_corr = evaluation(df_zero, cycles, isotopes, cup_config, database, mass_range, corr_isotopes, denom_corr_ratio)
                new_corr.data_bgd_corr(df_bgd_1, df_bgd_2)

            df_zero = df.data_zero_corr(sample)
            new_corr = evaluation(df_zero, cycles, isotopes, cup_config, database, mass_range, corr_isotopes, denom_corr_ratio)

            if line2_corr == True:
                new_corr.line2_corr(df_zero, isotope_line2_corr)

            methods = {1 : new_corr.raw_signals, 2: new_corr.raw_ratios, 3: new_corr.raw_ratios_corr,
                       4: new_corr.mass_fractionation, 5: new_corr.internal_norm_1, 6: new_corr.internal_norm_2,
                       7: new_corr.external_norm_Sb}

            if option < 4:
                data_sample[sample] = methods[option](isotope_denom)

            else:
                data_sample[sample] = methods[option](isotope_denom, iter_beta)

            avg_ratio_sample_all[sample] = new_corr.avg_to_df(data_sample, sample)
            sd2_ratio_sample_all[sample] = new_corr.SD_to_df(data_sample, sample)
            sample_names[sample] = df.extract_metadata(sample, "Sample Name")

            data_all = new_corr.to_df_all(sample_names, avg_ratio_sample_all, sd2_ratio_sample_all, ratios = False)

        return data_all
