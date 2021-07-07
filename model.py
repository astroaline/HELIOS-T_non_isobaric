from input import *
from load_files import *
import numpy as np
from scipy.interpolate import RegularGridInterpolator
import pdb


class Model:

    def __init__(self, len_x, x_full, bin_indices, parameter_dict, integral_grid):

        self.len_x = len_x
        self.x_full = x_full
        self.bin_indices = bin_indices
        self.parameter_dict = parameter_dict
        # self.opacity_grid = opacity_grid
        self.integral_grid = integral_grid


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

        epsilon = 0.000001
        p0 = 10 + epsilon   # add some tiny value to p0 to avoid infinities in integration


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

        pmin = 1e-6
        pressure_array = 10 ** np.array([-8, -7.66, -7.33, -7, -6.66, -6.33, -6, -5.66, -5.33, -5, -4.66, -4.33, -4, -3.66, -3.33, -3, -2.66, -2.33, -2, -1.66, -1.33, -1, -0.66, -0.33, 0, 0.33, 0.66, 1.0])

        pressure_values = pressure_array[np.where(pressure_array == pmin)[0][0]:]  # remove everything below pmin

        factor = np.sqrt(2 * scale_height * R0) / (kboltz * my_temp)

        kappa_grid = factor * self.integral_grid
        tau_values = kappa_grid * m_water  # convert opacity to cross-section

        def h(p0):

            h_values = np.zeros(len(tau_values[0]))

            for i in range(len(tau_values[0])):        # This should be 1458

                new_integrand = np.zeros(len(pressure_values))
 
                for j in range(len(pressure_values)):       # This should be 22

                    new_integrand[j] = (1 - np.exp(-tau_values[j,i]))/pressure_values[j]*(R0 + scale_height*np.log(p0/pressure_values[j]))

                    # print(new_integrand)


                integral_value = np.trapz(new_integrand, pressure_values)     # This should be a scalar
                # print(integral_value)
           
                h_values[i] = (scale_height/R0)*integral_value

            return pressure_values, h_values

        pressure_values, h_values = h(p0)
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
