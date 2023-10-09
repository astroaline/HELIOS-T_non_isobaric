from input import *

## info for plotting ##

parameter_ranges = {}
for parameter, prior in priors.items():
    parameter_ranges[parameter] = [prior[1], prior[0]+prior[1]]

parameter_labels = {'T': r'$T$ (K)', 'log_xh2o': r'log $X_{\rm H_2O}$', 'log_xnh3': r'log $X_{\rm NH_3}$',
                    'log_xco2': r'log $X_{\rm CO_2}$', 'log_xco': r'log $X_{\rm CO}$', 'log_xch4': r'log $X_{\rm CH_4}$',
                    'log_P_cloudtop': r'log $P_{\rm cloudtop}$ (bar)', 'log_P0': r'log $P_0$ (bar)', 'R0': r'$R_0$ ($R_{\rm Jup}$)', 
                    'Q0': r'$Q_0$', 'a': r'$a$', 'log_r_c': r'log $r_{\rm c}$ (cm)', 'log_p_cia': r'log $P_{\rm CIA}$ (bar)',
                    'log_tau_ref': r'log $\tau_{\rm ref}$', 'log_cloud_depth': r'$b_{\rm c}$',
                    'Rstar': r'$R_{\rm star}$ ($R_\odot$)', 'G': r'$g$ (cm s$^{-2}$)', 'line': r'flat-line'}      # labels for all possible parameters

parameter_colors = {'T': ['Reds', 0.4], 'log_xh2o': ['Blues', 0.4], 'log_xnh3': ['YlGn', 0.3],
                    'log_xco2': ['Oranges', 0.3], 'log_xco': ['GnBu', 0.3], 'log_xch4': ['RdPu', 0.3],
                    'log_P_cloudtop': ['YlOrRd', 0.3], 'log_P0': ['GnBu', 0.4], 'R0': ['PuRd', 0.4], 
                    'Q0': ['BuPu', 0.5], 'a': ['YlGnBu', 0.4], 'log_r_c': ['PuBu', 0.3], 'log_p_cia': ['YlOrBr', 0.4],
                    'log_tau_ref': ['RdPu', 0.4], 'log_cloud_depth': ['PuBu', 0.6],
                    'Rstar': ['YlOrRd', 0.4], 'G': ['YlGn', 0.4], 'line': ['Blues', 0.7]}      # colours for plots
