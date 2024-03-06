import sys, os
import astropy.units as u

def check(id_check, config):
    """
    Returns True if id_check says 'y' or 'yes', False otherwise.
    
    Parameters:
    --
    id_check: :class:`~str` key to check;

    config: :class:`~dictionary` dictionary containing the key; 
    """
    if id_check in config:
        config[id_check] = str(config[id_check].lower())
        if config[id_check] in ['y', 'yes']:
            return True
        elif config[id_check] in ['n', 'no']:
            return False
        else:
            print("Found "+ str(id_check) + " but can't understand it. I am assuming it is a 'no'...\n")
            return False
    else:
        return False

def get_parameters(list_of_parameters, label):
    """
    Get min and max information for a certain parameter.

    Parameters:
    - list_of_parameters (list): list of dictionaries containing parameters information
    - label (str): Name of the label to look up

    Returns:
    - tuple: (min, max) if the parameter is found, else (None, None)

    """
    for parameter in list_of_parameters:
        if parameter['name'].strip().lower() == label.strip("'").lower():
            return parameter['min'], parameter['max']
    return None, None

def recap(config):
    """
    Produces a name depending on the graphs conditions (histograms, contours) and prints a summary as doublecheck. 

    Parameters:
    ---
    config: :class:`~dictionary` dictionary containing the values needed;
    """
    name = config['file_path']
    name = os.path.basename(name)
    name, file_extension = os.path.splitext(name)

    print('Initialization complete. Here a short summary:')
    print("------- Summary -------")
    print("File name:\n- ", name+file_extension)
    intervals = 'Intervals:'
    if config['kind'] == 'histogram': 
        saving_strings_message = "The resulting plot(s) will be called:\n"
        name_dict = {}
        for i, label in enumerate(config['label_list']):
            label_min, label_max = get_parameters(config['list_of_parameters'], label)
            intervals += '\n-- '+ label + ' [' + str(label_min.value) + ', ' + str(label_max.value) +']' + str(label_min.unit)
            savename = name + '_histogram_' + label
            name_dict[label] = savename
            saving_strings_message = saving_strings_message + "                             -"+savename+'.png\n'
        config['name'] = name_dict
        saving_strings_message += 'They will be saved in: ' + config['file_saving_path']
    elif config['kind'] == 'contour':
        x_min, x_max = get_parameters(config['list_of_parameters'], config['x_label'])
        y_min, y_max = get_parameters(config['list_of_parameters'], config['y_label'])
        intervals = '-- ' + config['x_label'] + ':['+ str(x_min) + ', ' + str(x_max) + ']\n'
        intervals = intervals + '-- ' + config['y_label'] + ':[' + str(y_min) + ', ' + str(y_max) + ']'       
        name = name + '_contour_'+'_'+config['x_label']+'_'+config['y_label']
        config['name'] = name
        saving_strings_message = "The resulting plot will be called: " +name+".png\nIt will be saved in: "+config['file_saving_path']
    
    print("Kind of graph selectioned:", config['kind'])
    print(intervals)
    print(saving_strings_message)

    if not os.path.exists(config['file_saving_path']):
        os.makedirs(config['file_saving_path'])
    return config

def check_quantity(config, string, unit):
    """
    Check if the string is an astropy.unit or not.

    Parameters:
    ---
    config: :class:`~dictionary` dictionary containing the key;

    string: :class:`~str` string to check;

    unit:  :class:`~str` unit to associate.
    """

    try:
        config[string] = u.Unit(config[string]).to(u.Unit(unit))*u.Unit(unit)
    except ValueError:
        sys.stderr.write("Error: " + string + " must be an astropy.unit.\n")
        sys.exit(1)
    except KeyError:
        sys.stderr.write("Error: " + string + " missing.\n")
        sys.exit(1)

def check_str(config, string) -> dict:
    """
    Check if the string is present. Otherwise it returns an error.
    """
    try:
        config[string] = str(config[string])
        return config
    except KeyError:
        sys.stderr.write("Error: "+ string + " missing.\n")
        sys.exit()

def check_int(config, value_name) -> dict:
    """
    Check if the value is an int. Otherwise it returns an error.
    """
    try:
        config[value_name] = int(config[value_name])
        return config
    except KeyError:
        sys.stderr.write("Error: " + value_name + " missing.\n")
        sys.exit(1)
    except ValueError:
        sys.stderr.write("Error: " + value_name + " must be an int.\n")
        sys.exit(1)

def reader():
    """
    Reads the 'config.txt' files and check the main conditions. Returns a dictionary.
    """


    current_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    config_name = 'config_plot.txt'
    config_path = os.path.join(current_dir, config_name)
    save_dir = 'plots/graphs'
    config = {}
    try:
        with open(config_path, 'r') as config_file:
            for line in config_file:
                if line.startswith('#'):
                    # Skip comments
                    continue
                # Elaborate data
                key, value = line.strip().split(' = ')
                config[key] = value
    except (FileNotFoundError):
        config_file.close()
        sys.stderr.write('Error: file "config_plot.txt" not found.\n')
        sys.exit(1)
    except (AttributeError):
        config_file.close()
        sys.stderr.write('Error: "config_plot.txt" format file is not valid.\n')
        sys.stderr.write('Last valid input: ' + key + '\n')
        sys.exit(1)
    except ValueError:
        config_file.close()
        sys.stderr.write('Error: expected two values, found one.\n')
        sys.stderr.write('       Either one or more entries are empty or the format is "variable= number" instead of "variable = number".\n')
        sys.stderr.write('Last valid input: '+ key +'.\n')
        sys.exit(1)

    config_file.close()

    try:
            config['file_path'] = str(config['file_path'])
    except KeyError:
            sys.stderr.write("Error, 'file_path' missing.\n")
            sys.exit(1)


    try:
        if config['file_saving_path'] == "''" or '""':
            config['file_saving_path'] = os.path.join(current_dir, save_dir)
        else:
            config['file_saving_path'] = str(config['file_saving_path'])
    except KeyError:
            sys.stderr.write("Error: 'file_saving_path' value not specified.\n")
            sys.exit(1)


    quantities = {
        'pres_Min': 'cm-1 g s-2',
        'pres_Max': 'cm-1 g s-2',
        'dens_Min': 'g cm-3',
        'dens_Max': 'g cm-3',
        'energy_Min': 'erg',
        'energy_Max': 'erg',
        'T_Min': 'K',
        'T_Max': 'K',
        'DistanceUnit': 'cm',
        'TimeUnit': 's',
        'MassUnit': 'g',
    }

    for key in quantities.keys():
        check_quantity(config, key, quantities[key])


    list_of_parameters = [
    {'name': 'pressure', 'min': config['pres_Min'], 'max': config['pres_Max']},
    {'name': 'density', 'min': config['dens_Min'], 'max': config['dens_Max']},
    {'name': 'temperature', 'min': config['T_Min'], 'max': config['T_Max']},
    {'name': 'energy', 'min': config['energy_Min'], 'max': config['energy_Max']},
    ]
    config['list_of_parameters'] = list_of_parameters

    try:
        config['kind'] == str(config['kind'])
        if config['kind'] == 'histogram':
            label_list = []
            i = 1
            while i != None:
                try:
                    label_list.append(config['label'+str(i)])
                    i += 1
                except KeyError:
                    i = None
            config['label_list'] = label_list
            check_str(config, 'bin_rule')
            if not config['bin_rule'] in ['sqrt', 'rice', 'sturges']:
                sys.stderr.write('Error: "bin_rule" does not have the option', config['bin_rule'])
                sys.exit() 
        elif config['kind'] == 'contour':

            check_str(config,'x_label')
            check_str(config,'y_label')
            check_int(config, 'clim')
            check_int(config, 'n')

    except KeyError:
            sys.stderr.write("Error, 'kind' missing.\n")
            sys.exit(1)
    recap(config)

    return config
