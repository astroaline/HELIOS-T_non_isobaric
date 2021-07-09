import numpy as np
import os

## Constants ##

kboltz = 1.38064852e-16    # Boltzmann's constant
amu = 1.660539040e-24      # atomic mass unit
gamma = 0.57721
rjup = 7.1492e9            # equatorial radius of Jupiter
rsun = 6.9566e10           # solar radius
rearth = 6.378e8            # earth radius
pressure_probed = 1e-2      # probed pressure in bars
# pressure_cia = 1e-2         # pressure for cia in bars
# m = 2.4*amu                 # assummed hydrogen-dominated atmosphere
m_water = 18.0*amu          # mean molecular mass of any molecules you want to consider
m_cyanide = 27.0*amu
m_ammonia = 17.0*amu
m_methane = 16.0*amu
m_carbon_monoxide = 28.0*amu


## Planet Data ##

planet_name = 'WASP-52b'

g = 646
g_uperr = 45
g_loerr = 45
g_uncertainty = (g_uperr + g_loerr)/2
rstar = 0.79
rstar_uperr = 0.02
rstar_loerr = 0.02
rstar_uncertainty = (rstar_uperr + rstar_loerr)/2
r0 = 1.199
r0_uncertainty = 0.09

wavelength_bins = np.array([1.1107999999999998,1.1416,1.1709,1.1987999999999999,1.2257,1.2522,1.2791,1.3058,1.3321,1.3586,1.3860000000000001,1.414,1.4425,1.4718999999999998,1.5027,1.5345,1.5682,1.6042,1.6431999999999998])
transit_depth = np.array([2.718477117350561,2.720461986046827,2.692395426210081,2.71633103482348,2.699371522240168,2.6854591342633687,2.685953812973567,2.7113144731241237,2.734573447702161,2.7754339978059117,2.7465526304776047,2.7206363146982318,2.7154493398982633,2.7041020387754977,2.708223989100348,2.6867324036147693,2.681259162249782,2.687464478204251])
transit_depth_error = np.array([0.012273944764401792,0.009560440189967807,0.012308393828975021,0.010007946965824297,0.012573731038705637,0.010990328530552821,0.010862833736521506,0.012686953933709765,0.014366851120553534,0.012302853462452685,0.011207573135051597,0.012050238696995636,0.013560627518689413,0.013664400775793901,0.0113032052384349,0.014227050111014258,0.013246974684411633,0.016463991969084563])





## Retrieval info ##

model_name = 'greycloud'

molecules = ['1H2-16O__POKAZATEL_e2b']  # list of molecules (determines which opacity tables are loaded)
parameters = ["T", "log_xh2o", "log_kappa_cloud", "log_P0", "Rstar", "G"]   # parameters you wish to retrieve (MUST MATCH MOLECULES)
res = 2         # resolution used for opacities
live = 1000     # live points used in nested sampling
wavenumber=True     # True if opacity given in terms of wavenumber, False if wavelength

priors = {"T": [2800, 100], "log_xh2o": [13,-13], "log_xch4": [13,-13], "log_xco": [13,-13], "log_kappa_cloud": [14,-12],
          "log_P0": [4,-1], "log_kappa_0": [9,-10], "Q0": [99,1], "a": [3,3],
          "log_r_c": [6,-7], "log_p_cia": [3,-3], "Rstar": [rstar_uperr+rstar_loerr,rstar-rstar_loerr], "G": [g_uperr+g_loerr,g-g_loerr], "line": [5,0]} # priors for all possible parameters



## info for all possible parameters ##
molecular_abundance_dict = {'1H2-16O__POKAZATEL_e2b': 'log_xh2o', '12C-1H4__YT10to10_e2b': 'log_xch4', '12C-16O__HITEMP2010_e2b': 'log_xco'}  # dictionary list of all possible molecules and corresponding abundance names

parameter_dict = {"T": 1000, "log_xh2o": "Off", "log_xch4": "Off", "log_xco": "Off", "log_kappa_cloud": "Off", "R0": r0, "Rstar": rstar,
                  "log_P0": 1, "log_kappa_0": "Off", "Q0": "Off", "a": "Off", "log_r_c": "Off", "log_p_cia": -2, "G": g, "line": "Off"}    # default parameter values used if not retrieved

molecular_mass_dict = {'1H2-16O__POKAZATEL_e2b': m_water, '12C-1H4__YT10to10_e2b': m_methane, '12C-16O__HITEMP2010_e2b': m_carbon_monoxide}   # dictionary of molecules and their mean molecular masses
temperature_array = np.r_[50:700:50, 700:1500:100, 1500:2900:200]
temp_dict = {'01': temperature_array, '12C-1H4__YT10to10_e2b': temperature_array, '12C-16O__HITEMP2010_e2b': temperature_array}   # temperature values for corresponding opacity tables
temperature_array_cia = np.r_[200:3025:25]          # temperature array for CIA table
opacity_path = os.environ['HOME'] + "/Desktop/PhD/OPACITIES/"  # path to opacity binary files
cia_path = os.environ['HOME'] + "/Desktop/PhD/HITRAN/"      # path to CIA files
