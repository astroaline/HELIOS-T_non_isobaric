import numpy as np
from input import *
import struct


def load_opacity_data(molecule):
    ## Loads opacity tables from binary files ##

    wavenumber_min = 1e4/wavelength_bins[-1]
    wavenumber_max = 1e4/wavelength_bins[0]

    index_min = int((wavenumber_min)/res)
    index_max = int((wavenumber_max)/res) - 1

    for param in parameters:
        if param == 'log_xh2o':
            max_wavenumber_database = '42000'
            x_full = np.r_[0:42000:res]
            # wavenumber = np.r_[0:43000:1000]
        if param == 'log_xch4':
            max_wavenumber_database = '13000'
            x_full = np.r_[0:13000:res]
            #cwavenumber = np.r_[0:14000:1000]
        if param == 'log_xco':
            max_wavenumber_database = '09000'
            x_full = np.r_[0:9000:res]
            # wavenumber = np.r_[0:10000:1000]   
    
    pressure = int(np.log10(pressure_probed) * 100)
    if pressure < 0:
        pressure_str = 'n' + str(abs(pressure)).rjust(3, '0')
    else:
        pressure_str = 'p' + str(abs(pressure)).rjust(3, '0')

    temp = temp_dict[molecule]
    n_temp = len(temp)
    opacity_table = []
    for i in range(n_temp):
        data = []
        temp_str = str(temp[i]).zfill(5)
        filename = 'Out_00000_' + max_wavenumber_database + '_' + temp_str + '_' + pressure_str + '.bin'    
        with open(opacity_path + molecule + '/' + filename, "rb") as f:
            byte = f.read(4)
            while byte:
                data.extend(struct.unpack("f", byte))
                byte = f.read(4)
        step = int(res / 0.01)
        data = data[::step]
        opacity_table.append(data)
    opacity_table = np.array(opacity_table).T
    x_full = x_full[index_min:index_max]
    opacity_table = opacity_table[index_min:index_max, :]
    return x_full, opacity_table
    print(x_full, opacity_table)


def load_sigma(molecule1, molecule2, x_full): 
    if os.path.exists(cia_path+molecule1+'-'+molecule2+'_2011.cia'):
        file_name = molecule1+'-'+molecule2+'_2011.cia'
    elif os.path.exists(cia_path+molecule1+'-'+molecule2+'_2018.cia'):
        file_name = molecule1+'-'+molecule2+'_2018.cia'
    elif os.path.exists(cia_path+molecule1+'-'+molecule2+'_2018b.cia'):
        file_name = molecule1+'-'+molecule2+'_2018b.cia'
    elif os.path.exists(cia_path+molecule1+'-'+molecule2+'_norm_2011.cia'):
        file_name = molecule1+'-'+molecule2+'_norm_2011.cia'
    elif os.path.exists(cia_path+molecule1+'-'+molecule2+'_eq_2011.cia'):
        file_name = molecule1+'-'+molecule2+'_eq_2011.cia'
    else:
        print('{}-{} not found'.format(molecule1,molecule2))
        import sys
        sys.exit(1)

    with open(cia_path+file_name) as file:
        data = file.readlines()

    header = str(data[0])
    header_list = header.split(' ')
    header_list = [x for x in header_list if x]

    lines_per_temp = float(header_list[3])

    num_temps = len(temperature_array_cia)

    cia_data = []
    if molecule2 == 'H2':
        wavenumber_array = np.r_[x_full[0]:10001:res]
    else:
        wavenumber_array = np.r_[x_full[0]:17001:res]

    for i in range(int(num_temps)):
        j = int(i*(lines_per_temp+1))
        header = data[j]
        header_list = header.split(' ')
        header_list = [x for x in header_list if x]


        start_v = float(header_list[1])

        min_v = int(wavenumber_array[0] - start_v)
        max_v = int(wavenumber_array[-1] - start_v)

        cia_line = []

        for k in range(min_v, max_v+1):           # For CIA, temperature goes down rows, wavenumber goes across columns
            data_line = data[j+k].split(' ')
            data_line = [x for x in data_line if x]
            cia_line.append(float(data_line[1][:-1]))

        cia_line = cia_line[::res]

        cia_data.append(cia_line)

    cia_data = np.array(cia_data)
    wavenumber_array = np.array(wavenumber_array)

    pad_end = int((x_full[-1] - wavenumber_array[-1]) / res)

    if pad_end < 0:
        cia_data = cia_data[:, :pad_end]
    else:
        cia_data = np.pad(cia_data, ((0, 0), (0, pad_end)), 'constant')

    return cia_data



