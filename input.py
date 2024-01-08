import numpy as np
import os


## Constants ##
kboltz = 1.38064852e-16    # Boltzmann's constant
amu = 1.660539040e-24      # atomic mass unit
gamma = 0.57721
rjup = 7.1492e9            # equatorial radius of Jupiter
rsun = 6.9566e10           # solar radius
rearth = 6.378e8           # earth radius
pressure_probed = 1e-2     # probed pressure in bars
# pressure_cia = 1e-2      # pressure for cia in bars
# m = 2.4*amu              # assummed hydrogen-dominated atmosphere
m_water = 18.0*amu         # mean molecular mass of any molecules you want to consider
m_cyanide = 27.0*amu
m_ammonia = 17.0*amu
m_methane = 16.0*amu
m_carbon_dioxide = 44.0*amu
m_carbon_monoxide = 28.0*amu
solar_h2 = 0.5
solar_he = 0.085114

## Planet Data ##

planet_name = 'WASP-31b'

g = 456
g_uncertainty = 38
rstar = 1.252
rstar_uncertainty = 0.033
r0 = 1.379   # Rp = [1.499,1.599]
r0_uncertainty = 0.24   # R0 = [1.139,1.619]
teq = 1576
teq_uncertainty = 58

wavelength_bins = np.array([1.1107999999999998,1.1416,1.1709,1.1987999999999999,1.2257,1.2522,1.2791,1.3058,1.3321,1.3586,1.3860000000000001,1.414,1.4425,1.4718999999999998,1.5027,1.5345,1.5682,1.6042,1.6431999999999998])
transit_depth = np.array([1.5743159859215052,1.5859818163806723,1.588476101905992,1.5914625257359747,1.567865944837875,1.5929431931144669,1.5775325731036487,1.5726596045274135,1.5891307759796662,1.608101131959885,1.608512921004833,1.600552088151131,1.598115742763793,1.6190470153106553,1.5748452916752858,1.5717277695155019,1.5707941096962712,1.5865932382856405])
transit_depth_error = np.array([0.01734615604227587,0.016342276079246452,0.017926175418191116,0.013964552619703587,0.015159275024948123,0.014735339597991035,0.018407329258117697,0.015308681749044109,0.015073431681620895,0.013140862222708867,0.014375899511235704,0.01778640008048767,0.012593402150575894,0.015337523202063033,0.015881840919207208,0.014886186400454618,0.015236612861458759,0.014181554994181176])


## Retrieval info ##
model_name = 'greycloud'
approach_name = 'non_isobaric'

molecules = ['1H2-16O__POKAZATEL_e2b']
parameters = ['T', 'log_xh2o', 'log_P_cloudtop', 'R0', 'Rstar', 'G']
res = 2         # resolution used for opacities
live = 1000     # live points used in nested sampling
wavenumber = True     # True if opacity given in terms of wavenumber, False if wavelength
num_levels = 200

priors = {'T': [6*teq_uncertainty, teq-(3*teq_uncertainty)], 'log_xh2o': [13,-13], 'log_xnh3': [13,-13], 'log_xch4': [13,-13], 'log_xco2': [13,-13], 'log_xco': [13,-13],
          'log_P0': [4,-1], 'R0': [2*r0_uncertainty, r0-r0_uncertainty], 'log_tau_ref': [7,-5], 'Q0': [99,1], 'a': [10,3],
          'log_r_c': [6,-9], 'log_p_cia': [3,-3], 'log_P_cloudtop': [5,-4], 'log_cloud_depth': [2,0],
          'Rstar': [2*rstar_uncertainty, rstar-rstar_uncertainty],
          'G': [2*g_uncertainty, g-g_uncertainty], 'line': [5,0]} # priors for all possible parameters

pmin = 1e-6



## Info for all possible parameters ##
molecular_name_dict = {'1H2-16O__POKAZATEL_e2b': 'water', '14N-1H3__CoYuTe_e2b': 'ammonia', '12C-1H4__YT34to10_e2b': 'methane', '12C-16O2__CDSD_4000_e2b': 'carbon_dioxide', '12C-16O__Li2015_e2b': 'carbon_monoxide'}  # dictionary list of all possible molecules and corresponding names
molecular_abundance_dict = {'1H2-16O__POKAZATEL_e2b': 'log_xh2o', '14N-1H3__CoYuTe_e2b': 'log_xnh3', '12C-1H4__YT34to10_e2b': 'log_xch4', '12C-16O2__CDSD_4000_e2b': 'log_xco2', '12C-16O__Li2015_e2b': 'log_xco'}  # dictionary list of all possible molecules and corresponding abundance names

parameter_dict = {'T': 1000, 'log_xh2o': 'Off', 'log_xnh3': 'Off', 'log_xch4': 'Off', 'log_xco2': 'Off', 'log_xco': 'Off', 
                  'log_P_cloudtop': 'Off', 'log_P0': 1, 'R0': r0, 'Q0': 'Off', 'a': 'Off', 'log_r_c': 'Off', 'log_p_cia': -2,
                  'log_tau_ref': 'Off', 'log_cloud_depth': 'Off', 'Rstar': rstar, 'G': g, 'line': 'Off'}    # default parameter values used if not retrieved

molecular_mass_dict = {'1H2-16O__POKAZATEL_e2b': m_water, '14N-1H3__CoYuTe_e2b': m_ammonia, '12C-1H4__YT34to10_e2b': m_methane, '12C-16O2__CDSD_4000_e2b': m_carbon_dioxide, '12C-16O__Li2015_e2b': m_carbon_monoxide}   # dictionary of molecules and their mean molecular masses

wavenumber_dict = {'1H2-16O__POKAZATEL_e2b': '42000', '14N-1H3__CoYuTe_e2b': '20000', '12C-1H4__YT34to10_e2b': '12000', '12C-16O2__CDSD_4000_e2b': '09000', '12C-16O__Li2015_e2b': '22000'}

pressure_array_opacities = 10**np.array([-8, -7.66, -7.33, -7, -6.66, -6.33, -6, -5.66, -5.33, -5, -4.66, -4.33, -4, -3.66,
                                         -3.33, -3, -2.66, -2.33, -2, -1.66, -1.33, -1, -0.66, -0.33, 0, 0.33, 0.66, 1.0])  # full log pressure array for opacities in bars
temperature_array = np.r_[50:700:50, 700:1500:100, 1500:3100:200]

temp_dict = {'1H2-16O__POKAZATEL_e2b': temperature_array, '14N-1H3__CoYuTe_e2b': temperature_array[:22], '12C-1H4__YT34to10_e2b': temperature_array,
             '12C-16O2__CDSD_4000_e2b': temperature_array, '12C-16O__Li2015_e2b': temperature_array}   # temperature values for corresponding opacity tables
temperature_array_cia = np.r_[200:3025:25]          # temperature array for CIA table

drive_path = '/Users/astroaline/Library/CloudStorage/GoogleDrive-alinovais2@gmail.com/My Drive/ALINES DRIVE/'
opacity_path = drive_path + 'PhD/OPACITY_DATA/SPECIES/'   # path to opacity binary files
cia_path = drive_path + 'PhD/OPACITY_DATA/CIA/HITRAN/'    # path to CIA files

planet_file = planet_name + "_" + approach_name + "_" + model_name

this_path = '/Users/astroaline/Library/CloudStorage/GoogleDrive-alinovais2@gmail.com/My Drive/ALINES DRIVE/PhD/NEW/HELIOS-T_non_isobaric/'
planet_path = this_path + "planets/" + planet_name + "/"