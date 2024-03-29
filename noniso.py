from input import *
import load_files
import numpy as np
import struct
import pdb
from scipy.interpolate import RegularGridInterpolator, interp1d


def load_opacity(temperature, pressure, molecule):
    # res = 0.01   # resolution for the opacities
    step_size = int(res/0.01)

    wavenumber_min = int(1e4/wavelength_bins[-1])
    wavenumber_max = int(1e4/wavelength_bins[0])

    index_min = int((wavenumber_min)/res)
    index_max = int((wavenumber_max)/res)

    if res == 2:
        index_max = index_max - 1

    temp_str = str(temperature).zfill(5)   # temperature as in opacity filename
    pressure_load = int(np.log10(pressure) * 100)

    if pressure_load < 0:
        pressure_str = 'n' + str(abs(pressure_load)).rjust(3, '0')   # pressure as in opacity filename
    else:
        pressure_str = 'p' + str(abs(pressure_load)).rjust(3, '0')

    filename = molecule + '/Out_00000_' + wavenumber_dict[molecule] + '_' + temp_str + '_' + pressure_str + '.bin'

    data = []
    with open(opacity_path + filename, "rb") as f:
        byte = f.read(4)
        while byte:
            data.extend(struct.unpack("f", byte))
            byte = f.read(4)

    x_full = np.r_[0:42000:0.01]

    x_full = x_full[index_min * step_size:index_max * step_size:step_size]
    data = np.array(data[index_min * step_size:index_max * step_size:step_size])

    if len(data) < len(x_full):
        print('padding opacities...')
        data = np.pad(data, (0, len(x_full) - len(data)), mode='constant', constant_values=1e-6)

    return data, x_full



def load_cia(x_full):
    sigma_h2h2 = load_files.load_sigma('H2', 'H2', x_full)
    sigma_h2he = load_files.load_sigma('H2', 'He', x_full)

    sigma_cia = sigma_h2h2 + (solar_he/solar_h2)*sigma_h2he

    return sigma_cia



def interpolate_opacity(my_pressure, pressure_arr, x_full, opacity):
    fn = RegularGridInterpolator((pressure_arr, x_full), opacity, bounds_error=False, fill_value=None)
    pt = (my_pressure, x_full)
    y = fn(pt)
    return y



def tau(p0_bar):
    # Compute tau for all pressures
    pressure_array_pmin_opacities = pressure_array_opacities[np.where(pressure_array_opacities == pmin)[0][0]:]   # remove everything below pmin (this is in bars)

    pressure_levels_pmin_log = np.linspace(np.log10(pmin), np.log10(p0_bar), num_levels)   # log pressure array with num_levels (bars)
    pressure_levels_pmin = 10**pressure_levels_pmin_log
    p0_cgs = p0_bar * 1e6   # convert to cgs

    integral_dict = {}


    for molecule in molecules:
        _, x_full = load_opacity(temp_dict[molecule][0], pressure_levels_pmin[0], molecule)  # load one to get x_full
        integral_grid_molecule = np.zeros((len(temp_dict[molecule]), len(pressure_levels_pmin), len(x_full)))

        # we will integrate over pressure, for each temperature, for all wavelengths

        # Load integrands for all pressures
        for i, t in enumerate(temperature_array):

            # load opacities for all available pressures
            opacity_grid = np.zeros((len(pressure_array_pmin_opacities), len(x_full)))
            for j, p in enumerate(pressure_array_pmin_opacities):
                opacity_grid[j],_ = load_opacity(t, p, molecule)

            opacity_interpolator = interp1d(np.log10(pressure_array_pmin_opacities), opacity_grid, axis=0,
                                            bounds_error=False, fill_value=(opacity_grid[0], opacity_grid[-1]),
                                            assume_sorted=True)
            opacity_grid_all_levels = opacity_interpolator(pressure_levels_pmin_log)


            for j, p in enumerate(pressure_levels_pmin):
                p_cgs = p * 1e6

                pressure_sliced = pressure_levels_pmin[:j+1] * 1e6   # slice up to p and convert to cgs
                opacity_grid_sliced = opacity_grid_all_levels[:j+1]

                sigma = opacity_grid_sliced * molecular_mass_dict[molecule]   # array of (len(pressure_sliced),1458)
                integrand = (sigma.T * pressure_sliced).T   # array of (len(pressure_sliced),1458)
                y_integral = np.sqrt(np.log(p_cgs / pressure_sliced))   # array of len(pressure_sliced)

                integral_value = -np.trapz(integrand, y_integral, axis=0)   # negative because we're doing the integral upside down
                integral_grid_molecule[i, j] = integral_value

        integral_dict[molecule] = interp1d(temperature_array, integral_grid_molecule, axis=0,
                                           bounds_error=False,
                                           fill_value=(integral_grid_molecule[0], integral_grid_molecule[-1]),
                                           assume_sorted=True)  # save interpolator instead

    # now we do the same for CIA
    integral_cia_grid = np.zeros((len(temperature_array_cia), len(pressure_levels_pmin), len(x_full)))
    sigma_cia_full = load_cia(x_full)

    for i, t in enumerate(temperature_array_cia):

        sigma_cia_0 = np.array([sigma_cia_full[i]])   # reshape to 2d
        for j, p in enumerate(pressure_levels_pmin):

            p_cgs = p*1e6

            ntot_factor = 1 / kboltz / t   # extra number density since there are two species in each CIA
            sigma_cia = ntot_factor * sigma_cia_0

            pressure_sliced = pressure_levels_pmin[:j + 1] * 1e6   # pass in pressure values and integrand values for all pressures
            pressure_sliced_sq = pressure_sliced*pressure_sliced   # squared to account for P in ntot missing from ntot_factor
            y_integral = np.sqrt(np.log(p_cgs / pressure_sliced))
            integrand_cia = (sigma_cia.T * pressure_sliced_sq).T   # array of (len(pressure_sliced),1458)
            integral_value = -np.trapz(integrand_cia, y_integral, axis=0)   # calculate integral using trapezoid approximation

            integral_cia_grid[i, j] = integral_value


    # And Rayleigh scattering, but this one is independent of temperature
    sigma_rayleigh = np.array([8.4909e-45 * (x_full ** 4)])
    integral_rayleigh_grid = np.zeros((len(pressure_levels_pmin), len(x_full)))

    for j, p in enumerate(pressure_levels_pmin):
        p_cgs = p*1e6

        pressure_sliced = pressure_levels_pmin[:j + 1]*1e6   # pass in pressure values and integrand values for all pressures
        y_integral = np.sqrt(np.log(p_cgs / pressure_sliced))
        integrand_rayleigh = (sigma_rayleigh.T * pressure_sliced).T   # array of (len(pressure_sliced),1458)
        integral_value = -np.trapz(integrand_rayleigh, y_integral, axis=0)   # calculate integral using trapezoid approximation

        integral_rayleigh_grid[j] = integral_value


    # Here is a dictionary of integral grids, each with different associated temperature arrays !!!
    integral_dict['cia'] = interp1d(temperature_array_cia, integral_cia_grid, axis=0,
                                    bounds_error=False, fill_value=(integral_cia_grid[0], integral_cia_grid[-1]),
                                    assume_sorted=True) # save interpolator again

    integral_dict['rayleigh'] = integral_rayleigh_grid

    return integral_dict, x_full
