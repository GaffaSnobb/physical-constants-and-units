from numpy import pi, e, array
import codecs
import urllib.request
import os.path
import os
import multiprocessing
from urllib.error import HTTPError

def download_html(data_directory, which="all"):
    """
    Download the html files with mass data.

    Parameters
    ----------
    data_directory : string
        Name of the folder of downloaded data.

    which : list, string
        A list of element names by their abbreviations, or a single
        element as a string.
    """
    element_names_list = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F",
        "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
        "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga",
        "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo",
        "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I",
        "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu",
        "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W",
        "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po",
        "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am",
        "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db",
        "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc",
        "Lv", "Ts", "Og", "all"]

    element_names_dict = {"H":1, "He":2, "Li":3, "Be":4, "B":5, "C":6, "N":7,
        "O":8, "F":9, "Ne":10, "Na":11, "Mg":12, "Al":13, "Si":14, "P":15,
        "S":16, "Cl":17, "Ar":18, "K":19, "Ca":20, "Sc":21, "Ti":22, "V":23,
        "Cr":24, "Mn":25, "Fe":26, "Co":27, "Ni":28, "Cu":29, "Zn":30, "Ga":31,
        "Ge":32, "As":33, "Se":34, "Br":35, "Kr":36, "Rb":37, "Sr":38, "Y":39,
        "Zr":40, "Nb":41, "Mo":42, "Tc":43, "Ru":44, "Rh":45, "Pd":46, "Ag":47,
        "Cd":48, "In":49, "Sn":50, "Sb":51, "Te":52, "I":53, "Xe":54, "Cs":55,
        "Ba":56, "La":57, "Ce":58, "Pr":59, "Nd":60, "Pm":61, "Sm":62, "Eu":63,
        "Gd":64, "Tb":65, "Dy":66, "Ho":67, "Er":68, "Tm":69, "Yb":70, "Lu":71,
        "Hf":72, "Ta":73, "W":74, "Re":75, "Os":76, "Ir":77, "Pt":78, "Au":79,
        "Hg":80, "Tl":81, "Pb":82, "Bi":83, "Po":84, "At":85, "Rn":86, "Fr":87,
        "Ra":88, "Ac":89, "Th":90, "Pa":91, "U":92, "Np":93, "Pu":94, "Am":95,
        "Cm":96, "Bk":97, "Cf":98, "Es":99, "Fm":100, "Md":101, "No":102,
        "Lr":103, "Rf":104, "Db":105, "Sg":106, "Bh":107, "Hs":108, "Mt":109,
        "Ds":110, "Rg":111, "Cn":112, "Nh":113, "Fl":114, "Mc":115, "Lv":116,
        "Ts":117, "Og":118}

    if which == "all": which = element_names_list[:-1]
    
    if isinstance(which, str): which = [which]
    elif isinstance(which, list):
        if not set(which) < set(element_names_list):
            print("You have specified invalid elements. " +
            " Please use the two letter abbreviations, either in a list" + 
            " or a single element as a string.")
            print(element_names_list)
            return
    else:
        print(f"Valid input types are {list} and {str}, got {type(which)}.")
        return

    if not os.path.isdir(data_directory): os.mkdir(data_directory)
    
    pool = multiprocessing.Pool()
    pool.map(fetch_data_in_parallel,
        ((i, which, element_names_dict, data_directory) for i in range(len(which))))


def fetch_data_in_parallel(index_list_dict_folder):
    """
    The delay between request and recieve makes fetching data for many
    elements slow. This function is parallelised with the built-in
    multiprocessing module for fetching mass data in parallel.

    Parameters
    ----------
    index_list_dict_folder : tuple
        Since pool.map only supports a single function argument, all
        objects needed inside this function is passed as a tuple.

        i : int
            Current index of the input list of elements to fetch.

        which : list
            A list of elements to fetch.

        element_names_dict : dictionary
            A dictionary where key = element name, value = element
            number.

        data_directory : string
            Name of the sub-directory where the mass data are stored.
    """
    i, which, element_names_dict, data_directory = index_list_dict_folder
    
    element_number = f"{element_names_dict[which[i]]:03d}"
    filepath = data_directory + "/" + element_number + "-" + which[i] + "-masses.html"
    element_number = f"{element_names_dict[which[i]]:02d}"
    
    if not os.path.isfile(filepath):
        urlname = "https://wwwndc.jaea.go.jp/cgi-bin/nucltab14?" + element_number
        try:
            urllib.request.urlretrieve(urlname, filepath)
            print(filepath + " downloaded")
        except HTTPError:
            print(f"Unable to download {urlname}. Invalid URL.")


def extract_mass_from_html(filename):
    """
    Read (pre-downloaded) html files from 'wwwndc.jaea.go.jp' with mass
    data.
    """
    data = {}
    names = []


    with codecs.open(filename, "r") as infile:

        for i in range(3): infile.readline()    # Skip non-interesting lines.
        infile.readline()[4:14]  # Extract the name of the element.
        infile.readline()
        categories = infile.readline()  # Extract the titles of the table.
        infile.readline()

        max_iterations = 1000   # Avoid infinite loop.
        iterations = 0

        while True: # Iterate through the file, line by line.
            content = infile.readline()
            content = content.split()

            if (content[0][0:4] == "----"): break   # End of file.
            if (iterations > max_iterations):
                print("Reached maximum number of iterations. The program did"
                    + " most likely not find the end of the file.")
                break

            name = ""
            reference_point = 0 # Index of list element reference point where to find value and name data.
            name_location_1 = 0 # Substring index of (partial or entire) name of the isotope.
            
            for elem in content:
                """
                Use '</a>' as a point of reference.  Mass value is
                always the next element in the list 'content' from
                '</a>'.
                """
                if (elem[-4:] == "</a>"):
                    if (len(elem) > 4):
                        """
                        The name of the element is sometimes (partially
                        or entirely) in the same string as '</a>', eg:
                        '98</a>', so we extract all characters before
                        '<', which has index -4, if the length of elem
                        is greater than 4.
                        """
                        name = elem[0:-4]
                        
                        for char in elem[:-1]:
                            """
                            Sometimes the entire name is in the same
                            string as '</a>', eg: '...>Tc-102</a>', and
                            the name is always preceeded by '>'.  We
                            look for '>' to extract only the name of the
                            isotope. If the name is only partially
                            contained with '</a>' like '98</a>', then
                            there is no nothing to strip away in front,
                            but we must make sure that the final '>' is
                            not mistaken as the character preceeding the
                            name, hence 'elem[:-1].
                            """
                            name_location_1 += 1
                            
                            if (char == ">"):
                                name = elem[name_location_1:-4]
                                break
                            
                    break   # Break the loop if we find '</a>'.
                
                reference_point += 1

            try:
                """
                If the next element of 'content' can be converted to
                float, then it is the mass of the isotope.
                """
                value = float(content[reference_point+1])
            except ValueError:
                """
                If the next element of 'content' cannot be converted to
                float, then the line is not interesting and can be
                skipped entirely.
                """
                continue

            name_location_2 = 0 # Substring index of (partial or entire) name of the isotope.

            try:
                """
                If the previous element of 'content' can be converted to
                int, then that element is a part of the name of the
                isotope, ex. [..., ...H-, 1, </a>], and the rest of the
                name is located 2 steps back in the list.
                """
                name = int(content[reference_point-1])
                name = str(name)
                reference_point -= 2
            except ValueError:
                """
                If the previous element of 'content' cannot be converted
                to int, then the rest of the name is located 1 step back
                in the list.
                """
                reference_point -= 1
                
            for char in content[reference_point]:
                """
                Using '</a>' to find the start of the name of the
                isotope.
                """
                name_location_2 += 1
                if (char == ">"): break

            name = content[reference_point][name_location_2:] + name
            name = name.replace("-", "")
            names.append(name)
            data[name] = value

            iterations += 1

        # for key in data:
        #     print(key, data[key])

        return data, names



class MassesClass(dict):
    def __init__(self, which="all"):
        """
        Parameters
        ----------
        which : string, list
            A single string, or a list of strings, containing the
            abbreviated names of the elements you wish to download and
            load the mass data for. Note that the download only happens
            the first time you run the program, so it is probably the
            best to just leave this to the default value of 'all'. Input
            is case sensitive. Arbitrary ordering of input elements is
            allowed.

            Examples:
                which = 'all' - download and load all elements.
                which = 'U' - download and load mass data for uranium.
                which = ['H', 'He', 'Au'] - download and load mass data
                    for hydrogen, helium, and gold.
        """
        self.amu_to_kg = 1.66053904e-27     # 1 amu in kg.
        self.u_to_kg = self.amu_to_kg       # Alias.
        self.kg_to_amu = 1/self.amu_to_kg   # 1 kg in amu.
        self.kg_to_u = self.kg_to_amu       # Alias.
        
        self.amu_to_eV = 931.49410242e6     # 1 amu in eV/c**2.
        self.amu_to_ev = self.amu_to_eV     # Alias.
        self.u_to_eV = self.amu_to_eV       # Alias.
        self.u_to_ev = self.amu_to_eV       # Alias.
        self.eV_to_amu = 1/self.amu_to_eV   # 1 eV/c**2 in amu.
        self.ev_to_amu = self.eV_to_amu     # Alias.
        self.eV_to_u = self.eV_to_amu       # Alias.
        self.ev_to_u = self.eV_to_amu       # Alias.

        self.eV_to_kg = 1.782662e-36        # 1 eV/c**2 in kg.
        self.ev_to_kg = self.eV_to_kg       # Alias.
        self.kg_to_eV = 1/self.eV_to_kg  # 1 kg in eV/c**2.
        self.kg_to_ev = self.kg_to_eV       # Alias.

        self.mass_data_in_amu = {}
        self.mass_data_in_eV = {}
        self.mass_data_in_kg = {}
        self.names = []

        self.current_unit = ""

        # NUCLEAR DATA PART IS CURRENTLY NOT WORKING. MUST BE DEBUGGED.
        # mass_data_tmp = {}
        # names_tmp = []
        
        # data_directory = "nuclear_data/"
        
        # if __name__ == "__main__":
        #     """
        #     Does not work in an interactive session. Run the script
        #     directly before interactive usage.
        #     """
        #     download_html(data_directory, which)
        # else:
        #     if not os.path.isdir(data_directory):
        #         print("Mass data files from https://wwwndc.jaea.go.jp/NuC/"
        #          + " must be downloaded once first by directly running this script.")
        #         print(f"Script located at {os.path.dirname(os.path.realpath(__file__))}")
        #         return
        #     elif (len(os.listdir(data_directory)) == 0):
        #         print("Mass data files from https://wwwndc.jaea.go.jp/NuC/"
        #          + " must be downloaded once first by directly running this script.")
        #         print(f"Script located at {os.path.dirname(os.path.realpath(__file__))}")
        #         return
        
        # for filename in os.listdir(data_directory):
        #     if filename.endswith(".html"):
        #         mass_data_tmp, names_tmp = extract_mass_from_html(data_directory + filename)
        #         self.mass_data_in_amu.update(mass_data_tmp)
        #         self.names += names_tmp

        # for key in self.mass_data_in_amu:
        #     self.mass_data_in_eV[key] = self.mass_data_in_amu[key]*self.amu_to_eV
        #     self.mass_data_in_kg[key] = self.mass_data_in_amu[key]*self.amu_to_kg

        # self.names = sorted(self.names)
        self.set_unit()

    def set_unit(self, unit="ev"):
        """
        Set the mass unit.

        Parameters
        ----------
        unit : string
            Valid inputs are: 'ev', 'eV', 'ev/c^2', 'ev/c**2', 'eV/c^2',
            'eV/c**2', 'amu', 'u', 'kg'. Defaults to 'ev'.
        """

        if (unit == "ev") or (unit == "eV") or (unit == "ev/c^2") or (unit == "ev/c**2") or (unit == "eV/c^2") or (unit == "eV/c**2"):
            self.current_unit = "eV"
            self.set_attributes(self.mass_data_in_eV)
            fac = self.kg_to_eV # Base unit of all other masses is kg.
            print(f"Unit set to 'eV/c**2'.")
        
        elif (unit == "amu") or (unit == "u"):
            self.current_unit = "amu"
            self.set_attributes(self.mass_data_in_amu)
            fac = self.kg_to_amu    # Base unit of all other masses is kg.
            print(f"Unit set to 'amu'.")
       
        elif (unit == "kg"):
            self.current_unit = "kg"
            self.set_attributes(self.mass_data_in_kg)
            fac = 1 # Base unit of all other masses is kg.
            print(f"Unit set to 'kg'.")
        
        else:
            print(f"Invalid unit. Use 'eV', 'amu', or 'kg'. Got {unit}.")
            return
        
        
        # Nuclear elements (mass of nuclei).
        self.alpha = 6.644657230e-27*fac        # Mass of alpha particle, base [kg].

        # Hadrons.
        self.p = 1.6726219e-27*fac              # mass of proton, base [kg].
        self.proton = self.p
        self.n = 1.67492749804e-27*fac          # mass of neutron, base [kg].
        self.neutron = self.n

        # Leptons.
        self.e   = 9.10938356e-31*fac           # mass of electron, base [kg].
        self.electron = self.e
        self.mu  = 1.8835315972221714e-28*fac   # mass of muon, base [kg].
        self.muon = self.mu
        self.tau = 3.16754080132e-27*fac        # mass of tau, base [kg].
        self.tauon = self.tau

        # Mesons.
        self.pion = array([2.4880613534e-28, 2.406176557092e-28, 2.4880613534e-28])*fac # pi minus, pi zero, pi plus masses, base [kg].
        self.pi = self.pion
        self.pi_m = self.pion[0]    # Alias.
        self.pim = self.pion[0]     # Alias.
        self.pi_z = self.pion[1]    # Alias.
        self.piz = self.pion[1]     # Alias.
        self.pi_0 = self.pion[1]    # Alias.
        self.pi0 = self.pion[1]     # Alias.
        self.pi_n = self.pion[1]    # Alias.
        self.pin = self.pion[1]     # Alias.
        self.pi_p = self.pion[2]    # Alias.
        self.pip = self.pion[2]     # Alias.

        self.kaon = array([8.80059228174e-28, 8.87138178976e-28, 8.80059228174e-28])*fac # K minus, K zero, K plus masses, base [kg].
        self.K_m = self.kaon[0]    # Alias.
        self.Km = self.kaon[0]     # Alias.
        self.K_z = self.kaon[1]    # Alias.
        self.Kz = self.kaon[1]     # Alias.
        self.K_0 = self.kaon[1]    # Alias.
        self.K0 = self.kaon[1]     # Alias.
        self.K_n = self.kaon[1]    # Alias.
        self.Kn = self.kaon[1]     # Alias.
        self.K_p = self.kaon[2]    # Alias.
        self.Kp = self.kaon[2]     # Alias.
        
        ## Astronomy.
        self.sol = 1.9891e30*fac               # mass of sun, base [kg]
        self.sun = self.sol
        self.earth = 5.972e24*fac              # mass of earth, base [kg]
    
    def print_isotopes(self, which="all"):
        """
        Print a given set of available isotopes.

        Parameters
        ----------
        which : string
            Choose element for which to list all available isotopes and
            masses. Defaults to all.
        """
        if (self.current_unit == "eV"):
            current_mass_data = self.mass_data_in_eV
        elif (self.current_unit == "amu"):
            current_mass_data = self.mass_data_in_amu
        elif (self.current_unit == "kg"):
            current_mass_data = self.mass_data_in_kg
            
        for name in self.names:
            if (name.startswith(which)) or (which == "all"):
                print(f"{name}: {current_mass_data[name]}")


    def set_attributes(self, *args, **kwargs):
        """
        Update class attributes with entries in input dictionary. Each
        key will be the name of the attribute and each value will be
        the value of the attribute.
        """
        super(MassesClass, self).__init__(*args, **kwargs)
        self.__dict__.update(self)
        
        
class LifetimeClass:

    # Leptons.
    mu = 2.1969811e-6   # Muon lifetime, [s].
    muon = mu

    tau = 2.903e-13     # Tau lifetime, [s].
    tauon = tau

    # Mesons.
    pion = array([2.6e-8, 8.4e-17, 2.6e-8])    # pi minus, pi zero, pi plus lifetimes, [s].
    pi = pion
    pi_m = pion[0]
    pi_z = pion[1]
    pi_n = pi_z
    pi_p = pion[2]
    
    def __init__(self):
        pass

class HalfLifeClass:

    Mo99 = 237513.6         # Half-life of molybdenum, [s].
    mo99 = Mo99             # Alias.
    Tc99 = 6661667078600    # Half-life of technetium, [s].
    tc99 = Tc99             # Alias.
    Tc99m = 216241200       # Half-life of first excited state of technetium, [s].
    tc99m = Tc99m           # Alias.
    
    def __init__(self):
        pass

class UnitsClass:
    """
    Class of non-SI units converted to the base SI equivalent.
    """
    curie = 3.7e10  # Unit of radioactivity, [1/s].
    ci = curie
    Ci = curie

    barn = 1e-28    # Unit for expressing the cross-sectional area of nuclear reactions, [m**2].
    b = barn

    def __init__(self):
        pass

m = MassesClass()
t = LifetimeClass()
hl = HalfLifeClass()
u = UnitsClass()

man = "constants.m for masses. constants.t for lifetimes. constants.u for units. constants.hl for half-lives."
manual = man

# constants
c        = 299792458                     # Speed of light in vacuum, [m/s]
h        = 6.62607004e-34                # Planck's constant, [m^2*kg/s]
hbar     = h/(2*pi)                      # reduced Planck's constant, [m^2*kg/s]
h_ev     = 4.135667696e-15               # Planck's constant, [eV*s]
hbar_ev  = h_ev/(2*pi)                   # Reduced Planck's constant, [eV*s]
kb       = 1.38064852e-23                # Boltzmann's constant, [m^2*kg/(k*s^2)]
kb_ev    = 8.6173303e-5                  # Boltzmann's constant, [ev/K]
G        = 6.67408e-11                   # gravitational constant, [m^3/(kg*s^2)]
g        = 9.807                         # gravitational acceleration, [m/s**2]
R        = 8.3144598                     # gas constant, [J/(K mol)]
pi       = pi                            # Ratio between circumference and diameter of a circle, [unitless].
e        = e                             # Eulers number, [unitless].
sigma    = 2*pi**5*kb**4/(15*c**2*h**3)  # Stefan-Boltzmanns constant, [W/(m^2*K^4)].
# constants, units
c_unit     = "m/s"
h_unit     = "m^2*kg/s"
hbar_unit  = "m^2*kg/s"
kb_unit    = "m^2*kg/(K*s^2) = J/K"
kb_ev_unit = "ev/K"
G_unit     = "m^3/(kg*s^2)"
R_unit     = "J/(K mol)"
pi_unit    = "unitless"
e_unit     = "unitless"
sigma_unit = "W/(m^2*K^4)"


# lengths
AU = 149597871                  # astronomical unit, [km]
ly = 9.4605284e15               # light year, [m]
pc = 3.08567758e16              # parsec, [m]
r_sol = 695508e3                # solar radius, [m]
r_earth = 6371e3                # earth radius, [m]

# lengths, units
AU_unit  = "km"
ly_unit = "m"
pc_unit = "m"
r_sol_unit = "m"
r_earth_unit = "m"

# energies
eV_J = 1.60217662e-19           # electron volt to joule, [J]
ev_j = eV_J # Alias.
ev_J = eV_J # Alias.
ev_joule = eV_J # Alias.

# energies, units
eV_J_unit = "J"

# unsorted
atm         = 101325                  # atmospheric pressure, [Pa]
mol         = 6.0221415e23            # mole
T_sol       = 1.57e7                  # core temp of the sun, [K]
L_sol       = 3.828e26                # solar luminosity, [W]
yr          = 31556926                # seconds in a year, [s]
ec          = 1.60217662e-19          # elementary charge, [coulomb]
alpha       = 0.0072973525693         # fine structure constant
# unsorted, units
atm_unit   = "Pa"
mol_unit   = "unitless"
T_sol_unit = "K"
L_sol_unit = "W"
yr_unit    = "s"
ec_unit    = "C"
