import pymultinest
import numpy as np
from input import *
from asthetics import *
import data_setup
import ns_setup
import model
from matplotlib import pyplot as plt
import cornerplot
import spectrum
import time
import json
import os
import csv
import pdb
import pandas as pd

start = time.time()

print('\n\n\n' + 'n_live_points = ' + str(live) + '        num_levels = ' + str(num_levels))
print('\n\n\n' + 'now running .................... ' + planet_file + '\n\n\n')

#planet_file = planet_name + '_' + approach_name + '_' + model_name

n_params = len(parameters)

output_directory = 'out/'
os.makedirs(os.path.dirname(output_directory), exist_ok=True)

x, x_full, integral_dict, bin_indices, ydata, yerr, wavelength_centre, wavelength_err = data_setup.data()
len_x = len(x)


#pdb.set_trace()

# Run PyMultinest ##
b = ns_setup.Priors(1, n_params)
pymultinest.run(b.loglike, b.prior, n_params,
		        loglike_args=[len_x, x_full, bin_indices, integral_dict, ydata, yerr],
                outputfiles_basename=output_directory + planet_name + '_',
                resume=False, verbose=False,
                n_live_points=live)   # verbose is False if you do not want to print information while running code

json.dump(parameters, open(output_directory+planet_name+'_params.json', 'w'))  # save parameter names

a = pymultinest.Analyzer(outputfiles_basename=output_directory+planet_name+'_', n_params=n_params)
samples = a.get_equal_weighted_posterior()[:, :-1]
bestfit_params = a.get_best_fit()
#stats = a.get_stats()


## set up results ##
retrieved_results = list(map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]), zip(*np.percentile(samples, [16, 50, 84], axis=0))))

plot_percentiles = []
new_param_dict = parameter_dict
retrieved_parameters_full = {}
retrieved_parameters_list = []
for i in range(n_params):
    retrieved_parameters_full[parameters[i]] = retrieved_results[i]
    new_param_dict[parameters[i]] = bestfit_params['parameters'][i]
    # new_param_dict[parameters[i]] = retrieved_results[i][0]
    retrieved_parameters_list.append(retrieved_results[i][0])
    plot_percentiles.append(retrieved_results[i])


print("""PyMultinest result:""")
for param in parameters:
    print(param,' = ',retrieved_parameters_full[param][0],' +',retrieved_parameters_full[param][1],' -',retrieved_parameters_full[param][2])

a_lnZ = a.get_stats()['global evidence']
log_evidence = a_lnZ/np.log(10)

print
print ('************************')
print ('MAIN RESULT: Evidence Z ')
print ('************************')
print ('  log Z for model = %.1f' % (log_evidence))
print


## save retrieved results ##
with open(planet_path + 'RETRIEVED RESULTS/retrieved_results_' + planet_file + '.txt', 'w') as save_retrieved_results:
    for param in parameters:
        save_retrieved_results.write(str(param)+' = '+str(retrieved_parameters_full[param][0])+' + '+str(retrieved_parameters_full[param][1])+' - '+str(retrieved_parameters_full[param][2])+'\n')
    save_retrieved_results.write('\n')
    save_retrieved_results.write('log Z for model = %.1f' % (log_evidence))


if model_name == 'flat_line':
    pass

else:
    ## save posterior samples ##
    with open(planet_path + 'POSTERIOR SAMPLES/posterior_samples_' + planet_file + '.csv', 'w') as save_posterior_samples:
        for param in parameters:
            save_posterior_samples.write(str(param) + ' ')
        save_posterior_samples.write('\n')
        for i in range(len(samples)):
            for value in samples[i]:
                save_posterior_samples.write(format(value, '.6e') + ' ')
            save_posterior_samples.write('\n')


    ## compute model for retrieved results ##
    y = model.Model(len_x, x_full, bin_indices, new_param_dict, integral_dict)
    yfit = y.transit_depth()
    yfit_binned = y.binned_model()

    wavenumber_min = int(1e4 / wavelength_bins[-1])
    wavenumber_max = int(1e4 / wavelength_bins[0])

    index_min = int((wavenumber_min) / res)
    if res == 2:
        index_max = int((wavenumber_max) / res) - 1
    else:
        index_max = int((wavenumber_max) / res)
    step_size = int(res/0.01)


    ## change to wavelength ##
    if wavenumber is True:
        x_full_plot = 1.0 / (x_full * 1e-4)  # convert to wavelength
    else:
        x_full_plot = x_full


    ## save modeled spectrum ##
    df1 = pd.DataFrame({'wavelength': x_full_plot, 'transit_depth': yfit})
    df1.to_csv(planet_path + 'MODELED SPECTRA RESULTS/modeled_spectrum_results_' + planet_file + '.csv', index=False, sep=' ', header=True)

    df2 = pd.DataFrame({'wavelength_centre': wavelength_centre, 'yfit_binned': yfit_binned, 'ydata': ydata, 'wavelength_err': wavelength_err, 'yerr': yerr})
    df2.to_csv(planet_path + 'MODELED SPECTRA RESULTS/binned_spectrum_results_' + planet_file + '.csv', index=False, sep=' ', header=True)


    ## plot posteriors ##
    posterior_ranges, retrieved_parameter_labels, labels_units, labels_format, color_list = [], [], [], [], []
    
    for param in parameters:
        retrieved_parameter_labels.append(parameter_labels[param])
        color_list.append(parameter_colors[param])    
        labels_units.append(parameter_units[param])
        labels_format.append(parameter_format[param])
        #posterior_ranges.append(parameter_ranges[param])

    for i, param in enumerate(parameters):
        new_param_ranges = [float(np.percentile(samples[i], [1])), float(np.percentile(samples[i], [100]))] 
        posterior_ranges.append(new_param_ranges)

    line_color = 'b'
    
    fig = cornerplot.corner(samples, x_full_plot, wavelength_centre, yfit, yfit_binned, ydata, wavelength_err, yerr, posterior_ranges, line_color, color_list,
                            labels=retrieved_parameter_labels, labels_units=labels_units, truths=plot_percentiles, title_fmt=labels_format, max_n_ticks=3, show_titles=True)       
    fig2 = spectrum.spec(samples, x_full_plot, wavelength_centre, yfit, yfit_binned, ydata, wavelength_err, yerr, posterior_ranges, line_color, color_list,
                         labels=retrieved_parameter_labels, truths=plot_percentiles, max_n_ticks=3)  

    fig_name = 'cornerplot_' + planet_file
    fig.savefig(planet_path + 'CORNERPLOTS/' + fig_name + '.png', format='png')
    fig.savefig(planet_path + 'CORNERPLOTS/' + fig_name + '.pdf', format='pdf')

    fig2_name = 'spectrum_' + planet_file
    fig2.savefig(planet_path + 'SPECTRA/' + fig2_name + '.png', bbox_inches='tight', pad_inches=0.2, format='png')
    fig2.savefig(planet_path + 'SPECTRA/' + fig2_name + '.pdf', bbox_inches='tight', pad_inches=0.2, format='pdf')


end = time.time()

print('time in secs =', end-start)
print('time in mins =', (end-start)/60)
print('time in hours =', (end-start)/3600)