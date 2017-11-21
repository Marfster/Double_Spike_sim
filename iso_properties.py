__author__ = 'Matthias Friebel'



from toposort import toposort, toposort_flatten

# collection of dictionary and methods to import isotope masses into a mass database and read them out again
class Isotope_Masses():

    def __init__(self):
        self.isotope_masses = {} # dictionary e.g. {"124" : 123.9052745} - {isotope_mass_num : isotope mass}

    def add_masses_dict(self, masses_dict): #store dictionary of isotope masses in mass database
        self.isotope_masses = masses_dict

    def add_single_isotope(self, isotope_mass_no, isotope_mass): #add single isotope mass to database
        self.isotope_masses[isotope_mass_no] = isotope_mass

    def get_Isotope_mass(self, isotope_mass_no): # read out single isotope mass by mass number
        return self.isotope_masses[isotope_mass_no]

    def get_all_Isotope_masses(self): # read out all isotope masses at once
        return self.isotope_masses

# collection of dictionary and methods to import isotope ratios into a isotope ratio database and read them out again or convert them
class Isotope_Ratios():

    def __init__(self):
        self.isotope_ratios = {} # dictionary e.g. {"124" : 0.05} - {isotope_mass_num : isotope abundance}
        self.isotope_abundances = {} #dictionary e.g. {"120" : {"112" : 0.029182}} - {denominator_isotope : {nominator isotope : isotope_ratio}}

    def add_isotope_ratio(self, isotope_nom, isotope_denom, value): #add single isotope ratio to database
        if isotope_denom in self.isotope_ratios:
            self.isotope_ratios[isotope_denom][isotope_nom] = value
        else:
            self.isotope_ratios[isotope_denom] = {isotope_nom : value}

        self.isotope_ratios[isotope_denom][isotope_nom] = value

    def add_ratios_dict(self, isotope_denom, ratios_dict): #store dictionary of isotope ratio in ratio database
        self.isotope_ratios[isotope_denom] = ratios_dict

    def remove_ratio(self, isotope_nom, isotope_denom): #remove single isotope ratio from ratio database, necessary for conversion
        self.isotope_ratios[isotope_denom].pop(isotope_nom)

    def transf_ratio(self, isotope_denom_old, isotope_denom_new): #convert isotope ratios from one denominator to another
        for isotope in self.isotope_ratios[isotope_denom_old]:
            trans_ratio = self.isotope_ratios[isotope_denom_old][isotope] / self.isotope_ratios[isotope_denom_old][isotope_denom_new]
            self.add_isotope_ratio(isotope, isotope_denom_new, trans_ratio)

        invers_ratio = 1 / self.isotope_ratios[isotope_denom_old][isotope_denom_new]
        self.add_isotope_ratio(isotope_denom_old, isotope_denom_new, invers_ratio)

        self.remove_ratio(isotope_denom_new, isotope_denom_new)

    def calc_abundances(self, isotope_denom): #calculate isotope abundances from isotope ratios, needs denominator isotope of isotope ratios
        sum_ratios = 0
        for isotope in self.get_all_ratios(isotope_denom):
            sum_ratios += self.get_ratio(isotope, isotope_denom)
        for isotope in self.get_all_ratios(isotope_denom):
            self.isotope_abundances[isotope] = 100/(sum_ratios + 1) * self.get_ratio(isotope, isotope_denom)

        abundances_isotope_denom = 100
        for isotope in self.isotope_abundances:
            abundances_isotope_denom -= self.isotope_abundances[isotope]
        self.isotope_abundances[isotope_denom] = abundances_isotope_denom

    def get_ratio(self, isotope_nom, isotope_denom): # read out single isotope ratio from database based on nominator and denominator
        return self.isotope_ratios[isotope_denom][isotope_nom]

    def get_all_ratios(self, isotope_denom): # read out all isotope ratios from database at once
        return self.isotope_ratios[isotope_denom]

    def get_abundance(self, isotope, isotope_denom): #reads back isotope abundance from isotope ratio, needs nominator and denominator isotope
        self.calc_abundances(isotope_denom)
        return self.isotope_abundances[isotope]

    def get_all_abundances(self, isotope_denom): #reads back isotope abundances from database, needs denominator isotope of isotope ratios
        self.calc_abundances(isotope_denom)
        return self.isotope_abundances

# collection of dictionary and methods to import isotope abundances directly into a isotope ratio database and read them out again or convert them
class Isotope_Abundances(object):

    def __init__(self):
        self.isotope_abundances = {} # dictionary e.g. {"124" : 0.05} - {isotope_mass_num : isotope abundance}
        self.isotope_ratios = {} #dictionary e.g. {"120" : {"112" : 0.029182}} - {denominator_isotope : {nominator isotope : isotope_ratio}}

    def add_abundance(self, isotope_nom, value): # add single isotope abundance
        self.isotope_abundances[isotope_nom] = value

    def add_abundances_dict(self, abundances_dict): # read in dictionary of isotope abundances and scale them to 1
        sum_abundances = 0
        self.isotope_abundances = abundances_dict
        for isotope in self.isotope_abundances:
            sum_abundances += self.isotope_abundances[isotope]
        scale_factor = 1/sum_abundances
        for isotope in self.isotope_abundances:
            self.isotope_abundances[isotope] = self.isotope_abundances[isotope] * scale_factor # Scale to 1

    def calc_ratios(self, isotope_denom): # calculate ratios from isotope abundances by isotope denominator
        sum_abundances = 0
        self.isotope_ratios[isotope_denom] = {}
        for isotope in self.isotope_abundances:
            sum_abundances += self.isotope_abundances[isotope]
        scale_factor = 1/sum_abundances
        for isotope in self.isotope_abundances:
            self.isotope_abundances[isotope] = self.isotope_abundances[isotope] * scale_factor # Scale to 1

        for isotope in self.isotope_abundances:
            self.isotope_ratios[isotope_denom][isotope] = self.isotope_abundances[isotope]/self.isotope_abundances[isotope_denom]
        del self.isotope_ratios[isotope_denom][isotope_denom]

    def get_abundance(self, isotope): # read back a single isotope abundance based on isotope mass number
        return self.isotope_abundances[isotope]

    def get_all_abundances(self): # read back all isotope abundances at once
        return self.isotope_abundances

    def get_ratio(self, isotope_nom, isotope_denom): # read a single isotope ratio based on mass numbers for nominator and denominator
        return self.isotope_ratios[isotope_denom][isotope_nom]

    def get_all_ratios(self, isotope_denom): # read back all isotope ratios at once based on mass number for denominator
        self.calc_ratios(isotope_denom)
        return self.isotope_ratios[isotope_denom]

# collections of dictionary and methods to import the Elements associated with certain isotope mass number
class Isotopes_Mass_Range():

    def __init__(self):
        self.isotopes_mass_range = {} # dictionary e.g. {"124" : ["Sn", "Te", "Xe"]} - {isotope_mass_num : [Element_of_calculation, Interference element 1, ...]}

    def add_mass_range_dict(self, mass_range_dict):
        self.isotopes_mass_range = mass_range_dict

    def add_single_mass_isotope(self, isotope_mass_no, isotopes):
        self.isotopes_mass_range[isotope_mass_no] = isotopes

    def get_isotopes(self, isotope_mass_no):
        return self.isotopes_mass_range[isotope_mass_no]

    def get_mass_range(self):
        return self.isotopes_mass_range

    def get_interferences(self, isotope_mass_no, Element):
        isotopes = self.isotopes_mass_range[isotope_mass_no]
        isotopes.remove(Element)
        return isotopes

    # Graph of dependencies - dictionary with mass number of isotope used for calculation (e.g. 120 - Sn) and the associated mass number of isotope used for interference correction (e.g. 125 - Te)
    def get_graph_of_corr(self, corr_isotopes):
        graph = {} #{'117': set(), '118': set(), '119': set(), '120': {'125'}, '121': set(), '122': {'125'}, '123': {'125'}, '124': {'125', '129'}, ...}
        for mass_no in self.isotopes_mass_range:
            isotopes = self.isotopes_mass_range[mass_no]
            graph[mass_no] = set()
            for isotope in isotopes:
                if isotope in corr_isotopes and set(isotopes).issubset(corr_isotopes) == False:
                    if isotope == self.isotopes_mass_range[mass_no][0]:
                        None
                    else:
                        graph[mass_no].add(corr_isotopes[isotope])
                elif isotope not in corr_isotopes:
                    None
                elif len(self.isotopes_mass_range[corr_isotopes[isotope]]) == 1 and len(self.isotopes_mass_range[mass_no]) > 1:
                    if isotope == self.isotopes_mass_range[mass_no][0]:
                        None
                    else:
                        graph[mass_no].add(corr_isotopes[isotope])
        return graph
    # Order the dependencies - directed topology / creates scheme to follow for interference correction start with isotope with no dependencies,
    # than with isotopes with one dependency, than with two dependencies ...
    def get_order_of_corr(self, corr_isotopes):
        # [{'117', '118', '119', '121', '125', '127', '129', '131'}, {'120', '122', '123', '124', '126'}]
        # - first calculate with all isotopes in first dictionary (free of interferences), then calculate with second dictionary (one interference/dependency)
        return list(toposort(self.get_graph_of_corr(corr_isotopes)))


# load in isotope ratios into database and create a database of isotope abundances from it
def load_ratio_dict(dict_r, isotope_denom):
    ratio_dict = Isotope_Ratios()
    ratio_dict.add_ratios_dict(isotope_denom, dict_r)
    abund = Isotope_Abundances()
    abund.add_abundances_dict(ratio_dict.get_all_abundances(isotope_denom))
    return abund

# load in dictionary of isotope abundances & create database
def load_abundance_dict(dict_a):
    abund = Isotope_Abundances()
    abund.add_abundances_dict(dict_a)
    return abund

# load in dictionary for isotope masses & create database
def load_mass_dict(dict_m):
    mass_dict = Isotope_Masses()
    mass_dict.add_masses_dict(dict_m)
    return mass_dict
