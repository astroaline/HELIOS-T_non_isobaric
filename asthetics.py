from input import *

## info for plotting ##

parameter_ranges = {}
for parameter, prior in priors.items():
    parameter_ranges[parameter] = [prior[1], prior[0]+prior[1]]

parameter_labels = {'T': r'$T$ (K)', 'log_xh2o': r'log $X_{\rm H_2O}$', 'log_xhcn': r'log $X_{\rm HCN}$', 'log_xnh3': r'log $X_{\rm NH_3}$',
                    'log_xco2': r'log $X_{\rm CO_2}$', 'log_xco': r'log $X_{\rm CO}$', 'log_xch4': r'log $X_{\rm CH_4}$',
                    'log_P0': r'log $P_0$ (bar)', 'R0': r'$R_0$ ($R_{\rm Jup}$)',
                    'log_kappa_0': r'log $\kappa_0$ (cm$^2$ g$^{-1}$)', 'log_tau_ref': r'log $\tau_{\rm ref}$', 'log_cloud_depth': r'$b_{\rm c}$',                    
                    'log_kappa_cloud': r'log $\kappa_{\rm cloud}$ (cm$^2$ g$^{-1}$)', 'log_P_cloudtop': r'log $P_{\rm cloudtop}$ (bar)',
                    'Q0': r'$Q_0$', 'a': r'$a$', 'log_r_c': r'log $r_{\rm c}$ (cm)', 'log_p_cia': r'log $P_{\rm CIA}$ (bar)',
                    'Rstar': r'$R_{\rm star}$ ($R_\odot$)', 'G': r'$g$ (cm s$^{-2}$)', 'line': r'flat-line'}

parameter_units = {'T': r'$\rm{K}$', 'log_xh2o': r'', 'log_xhcn': r'', 'log_xnh3': r'',
                    'log_xco2': r'', 'log_xco': r'', 'log_xch4': r'',
                    'log_P0': r'$\rm{bar}$', 'R0': r'$R_{\rm Jup}$',
                    'log_kappa_0': r'$\rm{cm}^2 \rm{ g}^{-1}$', 'log_tau_ref': r'', 'log_cloud_depth': r'',                    
                    'log_kappa_cloud': r'$\rm{cm}^2\rm{\,g}^{-1}$', 'log_P_cloudtop': r'$\rm{bar}$',
                    'Q0': r'', 'a': r'', 'log_r_c': r'$\rm{cm}$', 'log_p_cia': r'$\rm{bar}$',
                    'Rstar': r'$R_\odot$', 'G': r'$\rm{cm\,s}^{-2}$', 'line': r''}

parameter_format = {'T': '.0f', 'log_xh2o': '.2f', 'log_xhcn': '.2f', 'log_xnh3': '.2f',
                    'log_xco2': '.2f', 'log_xco': '.2f', 'log_xch4': '.2f',
                    'log_P0': '.2f', 'R0': '.2f',
                    'log_kappa_0': '.2f', 'log_tau_ref': '.2f', 'log_cloud_depth': '.2f',                    
                    'log_kappa_cloud': '.2f', 'log_P_cloudtop': '.2f',
                    'Q0': '.2f', 'a': '.2f', 'log_r_c': '.2f', 'log_p_cia': '.2f',
                    'Rstar': '.2f', 'G': '.0f', 'line': '.2f'}
                                        
parameter_colors = {'T': ['Reds', 0.4], 'log_xh2o': ['Blues', 0.4], 'log_xhcn': ['YlOrRd', 0.15], 'log_xnh3': ['YlGn', 0.3],
                    'log_xco2': ['Oranges', 0.3], 'log_xco': ['GnBu', 0.3], 'log_xch4': ['RdPu', 0.3],
                    'log_P0': ['GnBu', 0.4], 'R0': ['PuRd', 0.4],
                    'log_kappa_0': ['copper_r', 0.0], 'log_tau_ref': ['RdPu', 0.4], 'log_cloud_depth': ['PuBu', 0.6],                    
                    'log_kappa_cloud': ['Purples', 0.4], 'log_P_cloudtop': ['YlOrRd', 0.3],
                    'Q0': ['BuPu', 0.5], 'a': ['YlGnBu', 0.4], 'log_r_c': ['PuBu', 0.3], 'log_p_cia': ['YlOrBr', 0.4],
                    'Rstar': ['YlOrRd', 0.4], 'G': ['YlGn', 0.4], 'line': ['Blues', 0.7]}
