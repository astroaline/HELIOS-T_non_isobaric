from input import *
from load_files import *
import numpy as np
from scipy.interpolate import RegularGridInterpolator



class Model:

    def __init__(self, len_x, x_full, bin_indices, parameter_dict, opacity_grid):

        self.len_x = len_x
        self.x_full = x_full
        self.bin_indices = bin_indices
        self.parameter_dict = parameter_dict
        self.opacity_grid = opacity_grid


    def molecular_opacity(self, my_temp, opacity_table):

        fn = RegularGridInterpolator((self.x_full, temperature_array), opacity_table)
        pt = (self.x_full, my_temp)
        y = fn(pt)
        return y


    def cia_cross_section(self, my_temp, cia_table):

        fn = RegularGridInterpolator((temperature_array_cia, self.x_full), cia_table)
        pt = (my_temp, self.x_full)
        y = fn(pt)
        return y


    def transit_depth(self):
        ## calculates transit depth ##

        my_temp = self.parameter_dict["T"]
        R0 = self.parameter_dict["R0"]*rjup
        P0 = (10**self.parameter_dict["log_P0"])*1e6   # convert to cgs
        Rstar = self.parameter_dict["Rstar"]*rsun
        G = self.parameter_dict["G"]
        pressure_cia = 10**self.parameter_dict["log_p_cia"]

        res = 2     # resolution for the opacities
        step_size = int(res/0.01)

        wavenumber_min = int(1e4/wavelength_bins[-1])
        wavenumber_max = int(1e4/wavelength_bins[0])

        index_min = int((wavenumber_min)/res)
        if res == 2:
            index_max = int((wavenumber_max)/res) - 1
        else:
            index_max = int((wavenumber_max)/res)

        epsilon = 0.000001
        p0 = 10 + epsilon   # add some tiny value to p0 to avoid infinities in integration
        


        def load_opacity(temperature, pressure):

            temp_str = str(temperature).zfill(5)     # temperature as in opacity filename

            pressure_load = int(np.log10(pressure) * 100)

            if pressure_load < 0:
                pressure_str = 'n' + str(abs(pressure_load)).rjust(3, '0')  # pressure as in opacity filename
            else:
                pressure_str = 'p' + str(abs(pressure_load)).rjust(3, '0')

            filename = '1H2-16O__POKAZATEL_e2b/Out_00000_42000_'+temp_str+'_'+pressure_str+'.bin'

            data = []
            opacity_table = []
            with open(opacity_path + filename, "rb") as f:
                byte = f.read(4)
                while byte:
                    data.extend(struct.unpack("f", byte))
                    byte = f.read(4)

            data = np.array(data[index_min*step_size:index_max*step_size:step_size])

            return data




        mass_fraction = []
        molecular_mass = []

        for molecule in molecules:
            abundance_name = molecular_abundance_dict[molecule]
            mass_fraction.append([10**self.parameter_dict[abundance_name]])
            molecular_mass.append([molecular_mass_dict[molecule]])

        mass_fraction = np.array(mass_fraction)
        molecular_mass = np.array(molecular_mass)

        xh2 = (1 - np.sum(mass_fraction))/1.1   # calculate abundance of H2
        if 'm' not in globals():                # set mean molecular weight if not given in input
            m = 2.4*xh2*amu + np.sum(mass_fraction*molecular_mass)
            
        scale_height = kboltz*my_temp/m/G



        def tau(temperature, pmin, p0):
        # Compute tau for all pressures

            pressure_array_pmin = pressure_array[np.where(pressure_array==pmin)[0][0]:] # remove everything below pmin

            opacity_line_length = int((wavenumber_max-wavenumber_min)/res) 
            integrand_grid = np.zeros((len(pressure_array_pmin), opacity_line_length))  # This will be the integrand
                                                                                        # we will integrate over pressure,
                                                                                        # for one temperature, for all wavelengths

            # Load integrands for all pressures
            for i, p in enumerate(pressure_array_pmin):
                opacity = load_opacity(temperature, p)   # load opacity for this temperature
                
                integrand_grid[i] = opacity/np.sqrt(np.log(p0/p))   # compute kappa/sqrt(ln(P0/P))

            kappa_grid = np.zeros((len(pressure_array_pmin), opacity_line_length))

            factor = np.sqrt(2*scale_height*R0)/(kboltz*temperature)

            # Integrate for each pressure p, from pmin to p
            for i, p in enumerate(pressure_array_pmin):

                pressure_sliced = pressure_array_pmin[:i+1]     # pass in pressure values and integrand values for all pressures
                integrand_grid_sliced = integrand_grid[:i+1]    # below p, above pmin

                integral_value = np.trapz(integrand_grid_sliced, pressure_sliced, axis=0)   # calculate integral using trapezoid approximation

                kappa_grid[i] = factor*integral_value

            tau_grid = kappa_grid*m_water     # convert opacity to cross-section

            return pressure_array_pmin, tau_grid

        pressure_values, tau_values = tau(my_temp, 1e-6, p0)
        
        


        def h(temperature, pmin, p0):

            h_values = np.zeros(len(tau_values[0]))

            for i in range(len(tau_values[0])):        # This should be 1458

                new_integrand = np.zeros(len(pressure_values))
 
                for j in range(len(pressure_values)):       # This should be 22

                    new_integrand[j] = (1 - np.exp(-tau_values[j,i]))/pressure_values[j]*(R0 + scale_height*np.log(p0/pressure_values[j]))
                    print(new_integrand)
                integral_value = np.trapz(new_integrand, pressure_values)     # This should be a scalar
                print(integral_value)
           
                h_values[i] = (scale_height/R0)*integral_value          

            return pressure_values, h_values

        pressure_values, h_values = h(my_temp, 1e-6, p0)
        #print(h_values)



        r = R0 + h_values
        
        result = 100.0 * (r / Rstar) ** 2  # return percentage transit depth
        #print(result)
        return result
        



    def binned_model(self):
    ## calculates average transit depth in given bins ##


        if self.parameter_dict['line'] == 'Off':
            y_full = self.transit_depth()
            #print(y_full)
            y_mean = np.zeros(self.len_x)
            for i in range(self.len_x):
                j = int(self.bin_indices[i])
                k = int(self.bin_indices[i + 1])
                y_mean[i] = np.mean(y_full[k:j])  # bin transit depth      
        else:
            y_mean = np.full(self.len_x, self.parameter_dict['line'])
            
        return y_mean
