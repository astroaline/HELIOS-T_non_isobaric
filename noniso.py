from input import *
from load_files import *
import numpy as np


def load_opacity(temperature, pressure):
    res = 2     # resolution for the opacities
    step_size = int(res/0.01)

    wavenumber_min = int(1e4/wavelength_bins[-1])
    wavenumber_max = int(1e4/wavelength_bins[0])

    index_min = int((wavenumber_min)/res)
    if res == 2:
        index_max = int((wavenumber_max)/res) - 1
    else:
        index_max = int((wavenumber_max)/res)

    temp_str = str(temperature).zfill(5)     # temperature as in opacity filename

    pressure_load = int(np.log10(pressure) * 100)

    if pressure_load < 0:
        pressure_str = 'n' + str(abs(pressure_load)).rjust(3, '0')  # pressure as in opacity filename
    else:
        pressure_str = 'p' + str(abs(pressure_load)).rjust(3, '0')
    
    filename = '1H2-16O__POKAZATEL_e2b/Out_00000_42000_' +temp_str +'_' +pressure_str +'.bin'

    data = []
    with open(opacity_path + filename, "rb") as f:
        byte = f.read(4)
        while byte:
            data.extend(struct.unpack("f", byte))
            byte = f.read(4)

    x_full = np.r_[0:42000:0.01]
    x_full = x_full[index_min * step_size:index_max * step_size:step_size]

    data = np.array(data[index_min * step_size:index_max * step_size:step_size])

    return data, x_full


def tau(temperature, pmin, p0):
    # Compute tau for all pressures

    pressure_array = 10 ** np.array([-8, -7.66, -7.33, -7, -6.66, -6.33, -6, -5.66, -5.33, -5, -4.66, -4.33, -4, -3.66,
                                     -3.33, -3, -2.66, -2.33, -2, -1.66, -1.33, -1, -0.66, -0.33, 0, 0.33, 0.66, 1.0])

    pressure_array_pmin = pressure_array[np.where(pressure_array == pmin)[0][0]:]  # remove everything below pmin

    wavenumber_min = int(1e4/wavelength_bins[-1])
    wavenumber_max = int(1e4/wavelength_bins[0])
    opacity_line_length = int((wavenumber_max - wavenumber_min) / res)
 
    if (opacity_line_length % 2) == 0:
        opacity_line_length = int((wavenumber_max - wavenumber_min) / res)
    else:
        opacity_line_length = int((wavenumber_max - wavenumber_min) / res) - 1

    integrand_grid = np.zeros((len(pressure_array_pmin), opacity_line_length))  # This will be the integrand
    # we will integrate over pressure,
    # for one temperature, for all wavelengths

    # Load integrands for all pressures
    for i, p in enumerate(pressure_array_pmin):
        opacity, x_full = load_opacity(temperature, p)  # load opacity for this temperature
        print(len(opacity), len(x_full))
        integrand_grid[i] = opacity / np.sqrt(np.log(p0 / p))  # compute kappa/sqrt(ln(P0/P))

    integral_grid = np.zeros((len(pressure_array_pmin), opacity_line_length))


    # Integrate for each pressure p, from pmin to p
    for i, p in enumerate(pressure_array_pmin):
        pressure_sliced = pressure_array_pmin[:i + 1]  # pass in pressure values and integrand values for all pressures
        integrand_grid_sliced = integrand_grid[:i + 1]  # below p, above pmin

        integral_value = np.trapz(integrand_grid_sliced, pressure_sliced,
                                  axis=0)  # calculate integral using trapezoid approximation

        integral_grid[i] = integral_value

    return pressure_array_pmin, integral_grid, x_full
