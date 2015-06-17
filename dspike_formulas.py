__author__ = 'marf'

from math import log10
from iso_properties import *
# Formulas for Double spike calculations#

#*** Natural Fractionation & Instrumental Fractionation***#
    # x = n or m, X = N or M
class dspike_formulas():

    def __init__(self, abundances_sample, abundances_spike, isotope_masses, list_spike_isotopes):

        # list_spike_isotopes = ['116']['117', '120', '122']

        self.sample_ratios = abundances_sample.get_all_ratios(list_spike_isotopes[0][0])
        self.spike_ratios = abundances_spike.get_all_ratios(list_spike_isotopes[0][0])
        self.columns_name = {'x' : str(list_spike_isotopes[1][0] + "/" + list_spike_isotopes[0][0]),
                             'y' : str(list_spike_isotopes[1][1] + "/" + list_spike_isotopes[0][0]),
                             'z' : str(list_spike_isotopes[1][2] + "/" + list_spike_isotopes[0][0])}
        self.x = {'x' : self.sample_ratios[list_spike_isotopes[1][0]],
                  'y' : self.sample_ratios[list_spike_isotopes[1][1]],
                  'z' : self.sample_ratios[list_spike_isotopes[1][2]]}
        self.SP = {'x' : self.spike_ratios[list_spike_isotopes[1][0]],
                  'y' : self.spike_ratios[list_spike_isotopes[1][1]],
                  'z' : self.spike_ratios[list_spike_isotopes[1][2]]}
        self.ma1 = {'x' : isotope_masses.get_Isotope_mass(list_spike_isotopes[1][0]),
                    'y' : isotope_masses.get_Isotope_mass(list_spike_isotopes[1][1]),
                    'z' : isotope_masses.get_Isotope_mass(list_spike_isotopes[1][2])}
        self.ma2 = isotope_masses.get_Isotope_mass(list_spike_isotopes[0][0])


    def X(self, x, pos, frac):
        X = x[pos] * (self.ma1[pos]/self.ma2) ** frac
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
        frac = log10(X[pos]/x[pos])/log10(self.ma1[pos]/self.ma2)
        return frac

class IterRegistry(type):
    def __iter__(cls):
        return iter(cls._registry)

class calc_dspike(object):
    __metaclass__ = IterRegistry
    _registry = []

    def __init__(self):
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

    def dspike_calc(self, cls_nat, cls_ins, iter_nat, iter_ins, frac_nat, frac_ins, frac_ratio):
        # cls_nat - sample object of class calc_dspike
        # cls_ins - mixture object of class calc_dspike
        # iter_nat - number of iterations used for calculation of natural fractionation
        # iter_ins - number of iterations used for calculation of instrumental fractionation
        # frac_nat - assumed initial natural fractionation
        # frac_ins - assumed initial instrumental fractionation
        # frac_ratio - which ratio: 'x', 'y' or 'z' for fractionation calculation should be used
        n = cls_nat.x
        m = cls_ins.x
        for i in range(iter_nat):
            # STEP 1#
            N = {}
            for ratio in cls_nat.x:
                N[ratio] = cls_nat.X(n, ratio, frac_nat)
            #Plane n-N-SP
            a_nat = cls_nat.a(n, N)
            b_nat = cls_nat.b(n, N)
            c_nat = cls_nat.c(n, a_nat, b_nat)

            #STEP 2#
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
                list_log_inner = [self.m.update({s:m}), self.M.update({s:M}), self.d_ins.update({s:d_ins}), self.e_ins.update({s:e_ins}),
                                  self.f_ins.update({s:f_ins}), self.g_ins.update({s:g_ins})]
                # Update M
                M['x'] = xint_ins
                M['y'] = yint_ins
                M['z'] = zint_ins

                frac_ins = cls_ins.frac(m, M, frac_ratio)
                self.frac_ins[s] = frac_ins

            #STEP 3#
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
            list_log_outer = [self.n.update({i:n}), self.N.update({i:N}), self.a_nat.update({i:a_nat}),
                              self.b_nat.update({i:b_nat}), self.c_nat.update({i:c_nat}), self.d_nat.update({i:d_nat}),
                              self.e_nat.update({i:e_nat}), self.f_nat.update({i:f_nat}), self.g_nat.update({i:g_nat}),
                              self.a_ins.update({i:a_ins}), self.b_ins.update({i:b_ins}), self.c_ins.update({i:c_ins})]
            # Update N

            N['x'] = xint_nat
            N['y'] = yint_nat
            N['z'] = zint_nat

            frac_nat = cls_nat.frac(n, N, frac_ratio)
            self.frac_nat[i] = frac_nat

        return frac_nat



class calc_dspike_samples(object):

    def __init__(self, Sn_meas_obj, data, spike_obj, Sn_mass_obj, spike_list, data_isotope_denom):
        self.Sn_std = Sn_meas_obj
        self.Sn_data = data
        self.Sn_spike = spike_obj
        self.mix = {}
        self.Sn_masses = Sn_meas_obj
        self.data_denom = data_isotope_denom
        self.std = dspike_formulas(Sn_meas_obj, spike_obj, Sn_mass_obj, spike_list)

    def mix_sim(self, fnat, fins, mix):
        # Mix Sample-Spike
        def fract(x, pos, ma1, ma2, frac):
                X = x[pos] * (ma1/ma2) ** frac
                return X

        def load_ratio_dict(dict_r, isotope_denom):
            ratio_dict = Isotope_Ratios()
            ratio_dict.add_ratios_dict(isotope_denom, dict_r)
            abund = Isotope_Abundances()
            abund = abund.add_abundances_dict(ratio_dict.get_all_abundances(isotope_denom))
            return abund

        def load_abundance_dict(dict_a):
            abund = Isotope_Abundances()
            abund.add_abundances_dict(dict_a)
            return abund

        Sn_std_sim_dict = {}
        for ratio in self.Sn_std.get_all_ratios(self.data_denom):
            Sn_std_sim_dict[ratio] = fract(self.Sn_std.get_all_ratios(self.data_denom), ratio,
                                                     self.Sn_masses.get_Isotope_mass(ratio),
                                                     self.Sn_masses.get_Isotope_mass(self.data_denom), fnat)
        Sn_std_sim = load_ratio_dict(Sn_std_sim_dict,self.data_denom)

        sample_spike_mix_abund_dict = {}
        for isotope in self.Sn_spike.get_all_abundances():
            sample_spike_mix_abund_dict[isotope] = mix * Sn_std_sim.get_all_abundances()[isotope] + (1 - mix) * self.Sn_spike.get_all_abundances()[isotope]

        sample_spike_mix_abund = load_abundance_dict(sample_spike_mix_abund_dict)
        sample_spike_mix_ratio = sample_spike_mix_abund.get_all_ratios(self.data_denom)

        sample_spike_mix_ratio_sim = {}
        for ratio in sample_spike_mix_ratio:
            sample_spike_mix_ratio_sim[ratio] = fract(sample_spike_mix_ratio, ratio,
                                                         self.Sn_masses.get_Isotope_mass(ratio),
                                                         self.Sn_masses.get_Isotope_mass(self.data_denom),fins)
        self.mix = sample_spike_mix_ratio_sim
        return self.mix
#       
    def spike_sim(self, fnat, fins, mix):
        nat_frac = []

        for value in range(len(self.data)):
            mix = {}
            for isotope in self.data:
                 mix[isotope] = (sample_spike_mix_ratio_sim[isotope]/df_new[isotope][value])*df_new[isotope].mean()

            mix_ratios = Isotope_Ratios()
            mix_ratios.add_ratios_dict(self.data_denom, mix)
            mix_abund = Isotope_Abundances()
            mix_abund.add_abundances_dict(mix_ratios.get_all_abundances(self.data_denom))

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

#   def dspike_calc:


#    def log_file:


#    def vis_sim:









