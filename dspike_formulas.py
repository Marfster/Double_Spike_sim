__author__ = 'Matthias Friebel'

from math import log10, log1p
from iso_properties import *
import collections
import flatdict
import numpy as np
import pandas as pd

class dspike_formulas():
    """ Contains all formulas used for Double Spike calculation after Siebert et al. (2001)

    abundances_sample - isotope composition of N (Isotopes_Abundances object)
    abundances_spike - isotope composition of SP (Isotopes_Abundances object)
    isotope_masses - Isotope masses of the Element (Isotopes_Masses object)
    list_spike_isotopes - four isotopes used for Double spike calculation
    [[denominator], nominator1, nominator2, nominator3] """

    def __init__(self, abundances_sample, abundances_spike, isotope_masses, list_spike_isotopes, law="exp", n_GPL=None):
        # input values defined in sn_config.py/parameter.py or are read in from ipython notebook

        self.sample_ratios = abundances_sample.get_all_ratios(list_spike_isotopes[0][0]) #calculate isotope ratios of sample/std relative to denominator isotope
        self.spike_ratios = abundances_spike.get_all_ratios(list_spike_isotopes[0][0]) #calculate isotope ratios of spike relative to denominator isotope
        self.columns_name = {'x' : str(list_spike_isotopes[1][0] + "/" + list_spike_isotopes[0][0]), # names of isotope columns
                             'y' : str(list_spike_isotopes[1][1] + "/" + list_spike_isotopes[0][0]),
                             'z' : str(list_spike_isotopes[1][2] + "/" + list_spike_isotopes[0][0])}
        self.x = {'x' : self.sample_ratios[list_spike_isotopes[1][0]], #isotope ratios of Std / natural sample
                  'y' : self.sample_ratios[list_spike_isotopes[1][1]],
                  'z' : self.sample_ratios[list_spike_isotopes[1][2]]}
        self.SP = {'x' : self.spike_ratios[list_spike_isotopes[1][0]], #isotope ratios of Sn Spike
                  'y' : self.spike_ratios[list_spike_isotopes[1][1]],
                  'z' : self.spike_ratios[list_spike_isotopes[1][2]]}
        self.ma1 = {'x' : isotope_masses.get_Isotope_mass(list_spike_isotopes[1][0]), # nominator isotope "x" for DS inv
                    'y' : isotope_masses.get_Isotope_mass(list_spike_isotopes[1][1]), # nominator isotope "y" for DS inv
                    'z' : isotope_masses.get_Isotope_mass(list_spike_isotopes[1][2])} # nominator isotope "z" for DS inv
        self.ma2 = isotope_masses.get_Isotope_mass(list_spike_isotopes[0][0]) # denominator isotope for DS inv
        self.law_mass_frac = law # law used for mass_bias correction, by default uses "exp" - exponential law
        self.n_GPL = n_GPL # n value for GPL, if not defined, n_GPL is by default a NONE-Typ


    def X(self, x, pos, frac):
        #*** Natural Fractionation & Instrumental Fractionation***#
        # x = n or m, X = N or M
        if self.law_mass_frac == "exp": # uses exponential law
            X = x[pos] * (self.ma1[pos]/self.ma2) ** frac
        elif self.law_mass_frac == "GPL": # uses GPL
            X = x[pos] * frac**(self.ma1[pos]**self.n_GPL - self.ma2**self.n_GPL)
        else:
            print 'wrong input for mass-fractionation law, use "exp" or "GPL"'
            X = None
        return X

    #describe plane x - X - SP with a, b, c:
    def a(self, x, X): # x = n or m, X = N or M, dictionaries : x = {'x' : 0.5, 'y' : 0.3, 'z' = 0.4}
        a = (x['y'] * (X['z']- self.SP['z']) + X['y'] * (self.SP['z'] - x['z']) + self.SP['y'] * (x['z'] - X['z'])) / (x['y'] * (X['x']- self.SP['x']) + X['y'] * (self.SP['x'] - x['x']) + self.SP['y'] * (x['x'] - X['x']))
        return a

    def b(self, x, X): # x = n or m, X = N or M, dictionaries : x = {'x' : 0.5, 'y' : 0.3, 'z' = 0.4}
        b = (x['x'] * (X['z']- self.SP['z']) + X['x'] * (self.SP['z'] - x['z']) + self.SP['x'] * (x['z'] - X['z'])) / (x['x'] * (X['y']- self.SP['y']) + X['x'] * (self.SP['y'] - x['y']) + self.SP['x'] * (x['y'] - X['y']))
        return b

    def c(self, x, a, b): # x = n or m, X = N or M, dictionaries : x = {'x' : 0.5, 'y' : 0.3, 'z' = 0.4}
        c = x['z'] - a * x['x'] - b * x['y']
        return c

    # describe line x-X with d, e, f, g:
    def d(self, x,  X):
        d = (x['z'] - X['z']) / (x['x'] - X['x'])
        return d

    def e(self, x, d):
        e = x['z'] - d * x['x']
        return e

    def f(self, x, X):
        f = (x['z'] - X['z']) / (x['y'] - X['y'])
        return f

    def g(self, x, f):
        g = x['z'] - f * x['y']
        return g

    # coordinates of intersection of line m-M with n-N-SP or n-N with m-M-SP

    def xint(self, a, b, c, d, e, f, g):
        xint = (b * g - b * e + e * f - c * f) / (a * f + b * d - d * f)
        return xint

    def yint(self, a, b, c, d, e, f, g):
        yint = (a * e - a * g + d * g - c * d) / (a * f + b * d - d * f)
        return yint

    def zint(self, a, b, c, xint, yint):
        zint = a * xint + b * yint + c
        return zint

    def frac(self, x, X, pos): # pos = 'x', 'y' or 'z'
        if self.law_mass_frac == "exp": # uses exponential law
            frac = np.log(X[pos]/x[pos])/np.log(self.ma1[pos]/self.ma2)
        elif self.law_mass_frac == "GPL": # uses GPL
            frac = (X[pos] / x[pos]) ** (1/(self.ma1[pos]**self.n_GPL - self.ma2**self.n_GPL))
        else:
            print 'wrong input for mass-fractionation law, use "exp" or "GPL"'
            frac = None
        return frac

class IterRegistry(type): # used to store single values from iteration of DS inversion
    def __iter__(cls):
        return iter(cls._registry)

class calc_dspike(object):
    """ Contains method for double spike calculation for one measurement of a sample """
    __metaclass__ = IterRegistry
    _registry = []

    def __init__(self):
        # empty dictionaries to store single values for each variable and step of iteration
        self._registry.append(self)
        self.n = {}
        self.N = {}
        self.a_nat = {}
        self.b_nat = {}
        self.c_nat = {}
        self.d_nat = {}
        self.e_nat = {}
        self.f_nat = {}
        self.g_nat = {}
        self.frac_nat = {}

        self.m = {}
        self.M = {}
        self.a_ins = {}
        self.b_ins = {}
        self.c_ins = {}
        self.d_ins = {}
        self.e_ins = {}
        self.f_ins = {}
        self.g_ins = {}
        self.frac_ins = {}
        self.log_dict = collections.OrderedDict()

    def dspike_calc(self, cls_nat, cls_ins, iter_nat, iter_ins, frac_nat, frac_ins, frac_ratio):
        """ Performs double spike calculation for one measurement of a sample
        # cls_nat - sample object of class calc_dspike
        # cls_ins - mixture object of class calc_dspike
        # iter_nat - number of iterations used for calculation of natural fractionation
        # iter_ins - number of iterations used for calculation of instrumental fractionation
        # frac_nat - assumed initial natural fractionation
        # frac_ins - assumed initial instrumental fractionation
        # frac_ratio - which ratio: 'x', 'y' or 'z' for fractionation calculation should be used
        """
        dict_log_inner = collections.OrderedDict() # Dictionary of values of Inner Inversion Loop
        dict_log_outer_1 = collections.OrderedDict() # Dictionary of values of Outer Inversion Loop - #Plane n-N-SP
        dict_log_outer_2 = collections.OrderedDict()

        n = cls_nat.x
        m = cls_ins.x
        for i in range(iter_nat):
            # STEP 1# - First Outer Loop - start by calculation natural isotope ratios assuming a certain natural fraction
            N = {}
            for ratio in cls_nat.x:
                N[ratio] = cls_nat.X(n, ratio, frac_nat)
            #Plane n-N-SP
            a_nat = cls_nat.a(n, N)
            b_nat = cls_nat.b(n, N)
            c_nat = cls_nat.c(n, a_nat, b_nat)

            dict_log_outer_1.update(sorted(flatdict.FlatDict({'n'+str(i):n}).items())
                                    + sorted(flatdict.FlatDict({'N'+str(i):N}).items())
                                    + {'a_nat'+str(i):a_nat}.items() + {'b_nat'+str(i):b_nat}.items()
                                    + {'c_nat'+str(i):c_nat}.items())

            #STEP 2# - Inner Loop - calculate mass-bias corrected Mixture - Line m-M and coordinations of intersection with plane n-N-SP
            for s in range(iter_ins):
                M = {}
                for ratio in cls_ins.x:
                    M[ratio] = cls_ins.X(m, ratio, frac_ins)
                #Line m-M
                d_ins = cls_ins.d(m, M)
                e_ins = cls_ins.e(m, d_ins)
                f_ins = cls_ins.f(m, M)
                g_ins = cls_ins.g(m, f_ins)

                xint_ins = cls_ins.xint(a_nat, b_nat, c_nat, d_ins, e_ins, f_ins, g_ins)
                yint_ins = cls_ins.yint(a_nat, b_nat, c_nat, d_ins, e_ins, f_ins, g_ins)
                zint_ins = cls_ins.zint(a_nat, b_nat, c_nat, xint_ins, yint_ins)

                # Log variables
                list_log_inner = [self.m.update({"m"+str(i)+'.'+str(s):m}), self.M.update({"M"+str(i)+'.'+str(s):M}),
                                  self.d_ins.update({'d_ins'+str(i)+'.'+str(s):d_ins}),
                                  self.e_ins.update({'e_ins'+str(i)+'.'+str(s):e_ins}),
                                  self.f_ins.update({'f_ins'+str(i)+'.'+str(s):f_ins}),
                                  self.g_ins.update({'g_ins'+str(i)+'.'+str(s):g_ins})]
                dict_log_inner.update(sorted(flatdict.FlatDict(self.m).items()) + sorted(flatdict.FlatDict(self.M).items())
                                      + self.d_ins.items() + self.e_ins.items() + self.f_ins.items() + self.g_ins.items())

                # Update M
                Mr = {}
                Mr['x'] = xint_ins
                Mr['y'] = yint_ins
                Mr['z'] = zint_ins

                frac_ins = cls_ins.frac(m, Mr, frac_ratio)
                self.frac_ins['frac_ins_'+frac_ratio+str(i)+'.'+str(s)] = frac_ins

                dict_log_inner.update(sorted(flatdict.FlatDict({"Mr"+str(i)+'.'+str(s):Mr}).items()) + self.frac_ins.items())

            #STEP 3# - second outer loop using mass-bias corrected M to calculate line n-N & plane m-M-SP and find intersection point to get correct N
            #Line n-N
            d_nat = cls_nat.d(n, N)
            e_nat = cls_nat.e(n, d_nat)
            f_nat = cls_nat.f(n, N)
            g_nat = cls_nat.g(n, f_nat)

            #Plane m-M-SP
            a_ins = cls_ins.a(m, M)
            b_ins = cls_ins.b(m, M)
            c_ins = cls_ins.c(m, a_ins, b_ins)

            xint_nat = cls_nat.xint(a_ins, b_ins, c_ins, d_nat, e_nat, f_nat, g_nat)
            yint_nat = cls_ins.yint(a_ins, b_ins, c_ins, d_nat, e_nat, f_nat, g_nat)
            zint_nat = cls_ins.zint(a_ins, b_ins, c_ins, xint_nat, yint_nat)

            # Log variables
            list_log_outer = [self.n.update({'n'+str(i):n}), self.N.update({'N'+str(i):N}),
                              self.a_nat.update({'a_nat'+str(i):a_nat}), self.b_nat.update({'b_nat'+str(i):b_nat}),
                              self.c_nat.update({'c_nat'+str(i):c_nat}), self.d_nat.update({'d_nat'+str(i):d_nat}),
                              self.e_nat.update({'e_nat'+str(i):e_nat}), self.f_nat.update({'f_nat'+str(i):f_nat}),
                              self.g_nat.update({'g_nat'+str(i):g_nat}), self.a_ins.update({'a_ins'+str(i):a_ins}),
                              self.b_ins.update({'b_ins'+str(i):b_ins}), self.c_ins.update({'c_ins'+str(i):c_ins})]

            dict_log_outer_2.update(self.d_nat.items() + self.e_nat.items() + self.f_nat.items()
                                    + self.g_nat.items() + self.a_ins.items() + self.b_ins.items()
                                    + self.c_ins.items())
            # Update N

            N['x'] = xint_nat
            N['y'] = yint_nat
            N['z'] = zint_nat

            frac_nat = cls_nat.frac(n, N, frac_ratio)
            self.frac_nat['frac_nat_'+frac_ratio+str(i)] = frac_nat

            dict_log_outer_2.update(sorted(flatdict.FlatDict({"Nr"+str(i):N}).items()) + self.frac_nat.items())
            self.log_dict.update(dict_log_outer_1.items() + dict_log_inner.items() + dict_log_outer_2.items())
        return frac_nat

class calc_dspike_sample(object):
    """ Contains methods for Spike Simulation and Double Spike correction

        # Sn_meas_obj - composition of the sample/standard (internal norm) (Isotopes_Abundances object)
        # data - dataframe containing the measurement (m for dspike correction, n for dspike simulation)
        # spike_obj = composition of the double spike (internal norm) (Isotopes_Abundances object)
        # Sn_mass_obj - Isotope masses of the Element (Isotopes_Masses object)
        # spike_lists - Isotope used for spike calculation [[denominator],[nominator1, nominator2, nominator3]]
        # data_isotope_denom - "Denominator isotope in measured data"
        # law_mass_frac - law used for mass-bias correction - "exp"- exponential law / "GPL" - general power law
        # n_GPL_ins - fixed assumed n for GPL of instrumental fractionation
        # n_GPL_nat - fixed assumed n for GPL of natural fractionation
        """
    def __init__(self, Sn_meas_obj, data, spike_obj, Sn_mass_obj, spike_list, data_isotope_denom, law, n_GPL_ins, n_GPL_nat):
        self.Sn_std = Sn_meas_obj
        self.Sn_data = data
        self.Sn_spike = spike_obj
        self.mix = {}
        self.Sn_masses = Sn_mass_obj
        self.spike_list = spike_list
        self.data_denom = data_isotope_denom
        self.std = dspike_formulas(Sn_meas_obj, spike_obj, Sn_mass_obj, spike_list, law, n_GPL_nat) # class with formulas and values of Sample-Spike-Mix used for DS Inversion - outer loop
        self.log_file_rang = {}
        self.law_mass_frac = law
        self.n_GPL_ins = n_GPL_ins
        self.n_GPL_nat = n_GPL_nat

    def mix_sim(self, fnat_sim, fins_sim, mix_ratio):
        """ Simulates a Sample-Spike-Mix with natural and instrumental fractionation
            # fnat_sim - simulated natural fractionation
            # fins_sim - simulated instrumental fractionation
            # mix_ratio - p (Ratio of Sample/Spike-Mix)
        """
        def fract(x, pos, ma1, ma2, frac):
                X = x[pos] * (ma1/ma2) ** frac
                return X

        Sn_std_sim_dict = {}
        for ratio in self.Sn_std.get_all_ratios(self.data_denom):
            Sn_std_sim_dict[ratio] = fract(self.Sn_std.get_all_ratios(self.data_denom), ratio,
                                                     self.Sn_masses.get_Isotope_mass(ratio),
                                                     self.Sn_masses.get_Isotope_mass(self.data_denom), fnat_sim)

        Sn_std_sim = load_ratio_dict(Sn_std_sim_dict,self.data_denom)

        sample_spike_mix_abund_dict = {}
        for isotope in self.Sn_spike.get_all_abundances():
            sample_spike_mix_abund_dict[isotope] = mix_ratio * Sn_std_sim.get_all_abundances()[isotope] + (1 - mix_ratio) * self.Sn_spike.get_all_abundances()[isotope]

        sample_spike_mix_abund = load_abundance_dict(sample_spike_mix_abund_dict)
        sample_spike_mix_ratio = sample_spike_mix_abund.get_all_ratios(self.data_denom)

        sample_spike_mix_ratio_sim = {}
        for ratio in sample_spike_mix_ratio:
            sample_spike_mix_ratio_sim[ratio] = fract(sample_spike_mix_ratio, ratio,
                                                         self.Sn_masses.get_Isotope_mass(ratio),
                                                         self.Sn_masses.get_Isotope_mass(self.data_denom),fins_sim)
        self.mix = sample_spike_mix_ratio_sim
        return self.mix
#       
    def spike_sim(self, fnat_sim, fins_sim, mix_ratio, dampening,iter_nat, iter_ins, frac_nat, frac_ins, frac_ratio):
        """ Performs Double-Spike Simulation

            # fnat_sim - simulated natural fractionation
            # fins_sim - simulated instrumental fractionation
            # mix_ratio - p (Ratio of Sample/Spike-Mix)
            # dampening - adds noise on "m" with variation of the measured n
            # iter_nat - number of iterations used for calculation of natural fractionation
            # iter_ins - number of iterations used for calculation of instrumental fractionation
            # frac_nat - assumed initial natural fractionation
            # frac_ins - assumed initial instrumental fractionation
            # frac_ratio - which ratio: 'x', 'y' or 'z' for fractionation calculation should be used"""

        # Creates an "m" for Sample-Spike-Mix with the data of the measurement for n
        mix_ini = self.mix_sim(fnat_sim, fins_sim, mix_ratio)
        mix_sim = {}
        mix_sim_mean = {}
        mix_sim_up = {}
        for value in range(len(self.Sn_data)):
            mix_sim[value] = {}
            for isotope in self.Sn_data:
                 mix_sim[value][isotope] = (mix_ini[isotope]/self.Sn_data[isotope].mean()) * self.Sn_data[isotope][value]


        for isotope in mix_sim[0]:
            mix_sim_mean[isotope] = {}
            sum = 0
            for value in range(len(self.Sn_data)):
                mix_sim_mean[isotope][value] = (mix_ini[isotope]/self.Sn_data[isotope].mean()) * self.Sn_data[isotope][value]
                sum += mix_sim_mean[isotope][value]
            mix_sim_mean[isotope] = sum/len(mix_sim_mean[isotope])

        for value in range(len(self.Sn_data)):
            mix_sim_up[value] = {}
            for isotope in mix_sim[value]:
                mix_sim_up[value][isotope] = (dampening * (mix_sim[value][isotope] - mix_sim_mean[isotope])) + mix_sim_mean[isotope]

         # Spike Calculation for all measurements of an sample#
        for value in range(len(mix_sim_up)):
            mix_sim_abund = load_ratio_dict(mix_sim_up[value],self.data_denom)
            mix = dspike_formulas(mix_sim_abund, self.Sn_spike, self.Sn_masses, self.spike_list)
            dspike_single = calc_dspike()

            dspike_single.dspike_calc(self.std, mix, iter_nat, iter_ins, frac_nat, frac_ins, frac_ratio)

        return self.log_file() # Return a log_file (dataframe) containing all parameters used Double-Spike calculation

    def spike_sim_p_range(self, mix_range, fnat_sim, fins_sim, dampening, iter_nat, iter_ins, frac_nat, frac_ins, frac_ratio):
        """ Double-Spike Simulation with varying Sample-Spike-Ratio (p)

        # mix_range - list of p-values (Ratio of Sample/Spike-Mix) [0.1, 0.2, 0.3, ...]
        # fnat_sim - simulated natural fractionation
        # fins_sim - simulated instrumental fractionation
        # dampening - adds noise on "m" with variation of the measured n
        # iter_nat - number of iterations used for calculation of natural fractionation
        # iter_ins - number of iterations used for calculation of instrumental fractionation
        # frac_nat - assumed initial natural fractionation
        # frac_ins - assumed initial instrumental fractionation
        # frac_ratio - which ratio: 'x', 'y' or 'z' for fractionation calculation should be used"""

        log_file_range = collections.OrderedDict()

        for mixed in mix_range:
            log_file_range.update({mixed : self.spike_sim(fnat_sim,fins_sim,mixed,dampening,iter_nat,iter_ins,frac_nat,frac_ins,frac_ratio)})

        return log_file_range # Returns a log-file for each p containing all parameters


    def dspike_corr(self, iter_nat, iter_ins, frac_nat, frac_ins, frac_ratio):
        """ Double Spike correction for all measurement lines of a measured Sample-Spike-Mix (self.Sn_data)

            # iter_nat - number of iterations used for calculation of natural fractionation
            # iter_ins - number of iterations used for calculation of instrumental fractionation
            # frac_nat - assumed initial natural fractionation
            # frac_ins - assumed initial instrumental fractionation
            # frac_ratio - which ratio: 'x', 'y' or 'z' for fractionation calculation should be used"""

        for index, row in self.Sn_data.iterrows(): # iterates over all measurement lines
            m_abund = load_ratio_dict(self.Sn_data.ix[index,:].to_dict(),self.data_denom) # calculates Isotope abundances for measured Sample-Spike-Mix
            m_cls = dspike_formulas(m_abund, self.Sn_spike, self.Sn_masses, self.spike_list, self.law_mass_frac, self.n_GPL_ins) # class with formulas and values of Sample-Spike-Mix used for DS Inversion - inner loop
            dspike_single = calc_dspike() #Loads in class with complete DS inversion for one sample and one measurement line

            dspike_single.dspike_calc(self.std, m_cls, iter_nat, iter_ins, frac_nat, frac_ins, frac_ratio) # Performs DS inversion for one measurement line

        return self.log_file() # Return a log_file (dataframe) containing all parameters used for Double-Spike Inversion

    def log_file(self):
        # creates a dataframe from the log-file dictionary
        counter = 0
        log_dict = {}

        for calc_dspike_object in calc_dspike:
            log_dict[counter] = calc_dspike_object.log_dict
            counter += 1

        log_df = pd.DataFrame.from_dict(log_dict, orient='index')
        calc_dspike._registry = []

        return log_df

    def error_vs_p(self, log_file_range, frac_ratio):
        """ Calculates the Error of fnat (natural fractionation)
            log-file_range returned from "spike_sim_p_range" method
            frac_ratio - which ratio: 'x', 'y' or 'z' for fractionation calculation should be used to calculate the error"""
        frac_nat_ppm = []
        S_SP_ratio = []
        for p in log_file_range:
            frac_nat_ppm.append(np.abs(log_file_range[p]['frac_nat_'+frac_ratio+"2"].std()))
            #frac_nat_ppm.append(np.abs(((log_file_range[p]['frac_nat_'+frac_ratio+"2"].std()/
                                           #np.sqrt(len(log_file_range[p]['frac_nat_'+frac_ratio+"2"])))/log_file_range[p]['frac_nat_'+frac_ratio+"2"].mean())*10**6))
            S_SP_ratio.append((1-p)/p)

        return S_SP_ratio, frac_nat_ppm

def spike_sim_q_range(q_range, spike1, spike2, Sn_meas_obj, df_new, Sn_mass_obj, spike_ls, mix, fnat_sim, fins_sim, dampening, iter_nat, iter_ins, frac_nat, frac_ins, frac_ratio):
    """ Double-Spike Simulation with varying Spike1-Spike2-Ratio (q)

        # q_range - list of q-values (Ratio of Spike1/Spike2-Mix) [0.1, 0.2, 0.3, ...]
        # spike1 - composition of spike1 (Isotopes_Abundances object)
        # spike2 - composition of spike2 (Isotopes_Abundances object)
        # Sn_meas_obj - composition of the sample/standard (internal norm) (Isotopes_Abundances object)
        # df_new - dataframe containing the measurement (m for dspike correction, n for dspike simulation)
        # Sn_mass_obj - Isotope masses of the Element (Isotopes_Masses object)
        # spike_ls - Isotope used for spike calculation [[denominator],[nominator1, nominator2, nominator3]]
        # mix - p (Ratio of Sample/Spike-Mix)
        # fnat_sim - simulated natural fractionation
        # fins_sim - simulated instrumental fractionation
        # dampening - adds noise on "m" with variation of the measured n
        # iter_nat - number of iterations used for calculation of natural fractionation
        # iter_ins - number of iterations used for calculation of instrumental fractionation
        # frac_nat - assumed initial natural fractionation
        # frac_ins - assumed initial instrumental fractionation
        # frac_ratio - which ratio: 'x', 'y' or 'z' for fractionation calculation should be used"""

    log_file_range = collections.OrderedDict()

    for q in q_range:
        spike_mix = {}
        for isotope in spike1:
             spike_mix[isotope] = spike1[isotope] * q + spike2[isotope] * (1-q)
        spike_obj = load_abundance_dict(spike_mix)
        q_sim = calc_dspike_sample(Sn_meas_obj,df_new,spike_obj,Sn_mass_obj,spike_ls,"120")
        log_file_range.update({q : q_sim.spike_sim(fnat_sim,fins_sim,mix,dampening,iter_nat,iter_ins,frac_nat,frac_ins,frac_ratio)})
#
    return log_file_range

def error_vs_q(log_file_range, frac_ratio):
    """ Calculates the Error of fnat (natural fractionation)
        log-file_range returned from "spike_sim_q_range" method
        frac_ratio - which ratio: 'x', 'y' or 'z' for fractionation calculation should be used to calculate the error"""
    frac_nat_ppm = []
    q_ls = []
    for q in log_file_range:
        frac_nat_ppm.append(np.abs(((log_file_range[q]['frac_nat_'+frac_ratio+"2"].std()/
                                           np.sqrt(len(log_file_range[q]['frac_nat_'+frac_ratio+"2"])))/log_file_range[q]['frac_nat_'+frac_ratio+"2"].mean())*10**6))
        #frac_nat_ppm.append(np.abs(log_file_range[q]['frac_nat_'+frac_ratio+"2"].std()/
        #                                   np.sqrt(len(log_file_range[q]['frac_nat_'+frac_ratio+"2"]))))
        q_ls.append(q)

    return q_ls, frac_nat_ppm













