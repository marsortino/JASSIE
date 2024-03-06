import sys
import os
import astropy.units as u
from astropy.constants import G, M_sun, c
from agnpy.targets import SSDisk
import numpy as np
import lib.conversions as cs

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
            print("Found "+ str(id_check) + " can't understand it. I am assuming it is a 'no'...\n")
            return False
    else:
        return False

def check_emissions(id_check, config):
    """
    Checks if the emissions region needs to be computed or not. Returns config[id_check] = True if 'y' or 'yes', False otherwise.
    Main difference between check() is that this function modifies the entry of the dictionary. While check() main purpose is to serve as a simple check.
   
    Parameters:
    --
    id_check: :class:`~str` key to check;

    config: :class:`~dictionary` dictionary containing the key; 
    """


    try:
        config[id_check] = (config[id_check].lower())
        if config[id_check] in ['y', 'yes']:
            config[id_check] = True
            return config[id_check]
        elif config[id_check] in ['n', 'no']:
            config[id_check] = False
        else:
            print("Found " + id_check + " can't understand it. I am assuming it is a 'no'...\n")
            config[id_check] = False
        return config[id_check]
    except KeyError:
        sys.stderr.write("Error: "+ id_check + " not specified.\n")
        sys.exit(1)

def check_float(config, string):
    """
    Check if the string is a float or not.

    Parameters:
    --
    config: :class:`~dictionary` dictionary containing the key;

    string: :class:`~str` string to check; 
    """

    try:
        config[string] = float(config[string])
    except ValueError:
        sys.stderr.write("Error: " + string + " must be a float.\n")
        sys.exit(1)
    except KeyError:
        sys.stderr.write("Error: " + string + " is missing.\n")
        sys.exit(1)

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

def disc_reader(config):
    """
    Defines an agnpy.target.SSDisk and appends it to config['target_list'] if it exists, creates the key otherwise.

    'a' and 'P_jet' are needed if no 'M_BH' is stated. The black hole mass is then computed as:

            M_BH = (P_jet*1e8*M_sun/(a**2*2*1e45*u.Unit('erg s-1'))).to(u.g) 

        (A.R. King, J.E. Pringle 2021 [3], DOI 10.3847/2041-8213/ac19a1)
            
        where:

        'a' is the dimensionless spin factor

        'P_jet' is the jet kinetic power (erg s-1)

    Entry needed:
    -
    'eta': accretion efficiency

    'L_disc': Disc luminosity (erg s-1)

    'R_g_unit': Gravitational radius - R_g = (2*G*Mass_BlackHole/c**2), can be yes or no. 

                If yes R_in and R_out are stated as multiple of R_g. E.g. R_in = 6 -> R_in = 6 R_G

    'R_in': Internal disc radius. If 'R_g_unit' is 'y'/'yes' then it must be a float.

    'R_out': External disc radius. If 'R_g_unit' is 'y'/'yes' then it must be a float. 

    Parameters:
    --
    config: :class:`~dictionary` dictionary containing the values needed;
    """
    removelist = ['a', 'P_jet', 'eta', 'L_disc', 'R_g_unit', 'R_in', 'R_out']
    if 'M_BH' not in config:
        try:
            a = int(config['a']) 
        except ValueError:
            sys.stderr.write("Error: 'a' must be an int.\n")
            sys.exit(1)
        except KeyError:
            sys.stderr.write("Error: 'a' is missing.\n")
            sys.exit(1)   

        try:
            P_jet = u.Unit(config['P_jet']).to(u.Unit('erg s-1'))*u.Unit('erg s-1') 
        except ValueError:
            sys.stderr.write("Error: 'P_jet' must be an astropy.unit.\n")
            sys.exit(1)
        except KeyError:
            sys.stderr.write("Error: 'P_jet' is missing.\n")
            sys.exit(1)
   
        Mass_BlackHole = (P_jet*1e8*M_sun/(a**2*2*1e45*u.Unit('erg s-1'))).to(u.g)
    else:
        try:
            Mass_BlackHole = u.Unit(config['M_BH']).to(u.Unit('g'))*u.Unit('g')
        except ValueError:
            sys.stderr.write("Error: 'M_BH' must be an astropy.unit.\n")
            sys.exit(1)
        except KeyError:
            sys.stderr.write("Error: 'M_BH' is missing.\n")
            sys.exit(1)                   

    try:
        eta = float(config['eta']) 
    except ValueError:
        sys.stderr.write("Error: 'eta' must be a float.\n")
        sys.exit(1)
    except KeyError:
        sys.stderr.write("Error: 'eta' is missing.\n")
        sys.exit(1)   
    try:
        L_disc = u.Unit(config['L_disc']).to(u.Unit('erg s-1'))*u.Unit('erg s-1')
    except ValueError:
        sys.stderr.write("Error: 'L_disc' must be an astropy.unit.\n")
        sys.exit(1)
    except KeyError:
            sys.stderr.write("Error: 'L_disc' is missing.\n")
            sys.exit(1)   
    
    if check('R_g_unit', config) == True:
        R_g = (2*G*Mass_BlackHole/c**2).to(u.cm)
        try:
            R_in = float(config['R_in'])*R_g
        except ValueError:
            sys.stderr.write("Error: since 'R_g_unit' is active, 'R_in' must be a float.\n")
            sys.exit(1)
        except KeyError:
            sys.stderr.write("Error: 'R_in' is missing.\n")
            sys.exit(1)

        try:
            R_out = float(config['R_out'])*R_g
        except ValueError:
            sys.stderr.write("Error: since 'R_g_unit' is active, 'R_out' must be a float.\n")
            sys.exit(1)
        except KeyError:
            sys.stderr.write("Error: 'R_out' is missing.\n")
            sys.exit(1)

    else:
        try:
            R_in = u.Unit(config['R_in']).to(u.Unit('cm'))*u.Unit('cm')
        except ValueError:
            sys.stderr.write("Error: 'R_in' must be an astropy.unit.\n")
            sys.exit(1)
        except KeyError:
            sys.stderr.write("Error: 'R_in' is missing.\n")
            sys.exit(1)
        try:
            R_in = u.Unit(config['R_out']).to(u.Unit('cm'))*u.Unit('cm')
        except ValueError:
            sys.stderr.write("Error: 'R_out' must be an astropy.unit.\n")
            sys.exit(1)
        except KeyError:
            sys.stderr.write("Error: 'R_out' is missing.\n")
            sys.exit(1)
    
    for key in removelist:
        del config[key]
    if 'target_list' in config:
        config['target_list'].append(SSDisk(Mass_BlackHole, L_disc, eta, R_in, R_out))
    else:
        config['target_list'] = []
        config['target_list'].append(SSDisk(Mass_BlackHole, L_disc, eta, R_in, R_out))

    return config

def recap(config):
    """
    Produces a name depending on the plotting conditions (e.g. bremmstrahlung, synchrotron, disc etc) and prints a summary as doublecheck. 

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
    name = name +''
    
    print("Emission processes:")
    if config['id_brem']:
        name = name +'_B'
        print("- Bremmstrahlung")
    if config['id_syn']:
        name = name +'_S'
        print("- Synchrotron")
        if config['id_ssa']:
            print("- - self-absorption")
            name = name +'_a'
        if config['id_ssc']:
            name = name +'_c'
            print("- - self-Compton")
    if config['id_ec']:
        name= name +'_EC'
        print("- External Compton")
        if config['id_cmb']:
            name = name+'_cmb'
            print("- - CMB")
        if check('id_disc', config):
            name = name+'_disc'
            print("- - Disc")
        if check('id_dusttorus', config):
            name = name+'_dt'
            print("- - Dust Torus")
        if check('id_blr', config):
            name = name+'_blr'
            print("- - BLR")
    elif config['id_ec'] == False:
        warning = '\nAttention, External Compton is set as off but the following processes are still set to be computed:'
        i = 0
        if config['id_cmb']:
            warning = warning + '\n- CMB'
            i += 1
        if check('id_disc', config):
            warning = warning +'\n- Disc'
            i += 1
        if check('id_dusttorus', config):
            warning = warning +'\n- Dust Torus'
            i += 1
        if check('id_blr', config):
            warning = warning+'\n- BLR'
            i += 1
        if i != 0:
            print(warning)
            answer = input('\nThis may result in long computational times without any contribute to the final spectrum. Do you still want to proceed? (y/yes)\n')
            if answer != 'y' and answer != 'yes':
                print("Interrupting.")
                sys.exit(1)
    
    print("The selected power law is:", config['n_e'])

    config['name'] = name
    print("The resulting plot will be called: "+ name+ ".png")
    print("It will be saved in:", config['file_saving_path'])
    if not os.path.exists(config['file_saving_path']):
        os.makedirs(config['file_saving_path'])
    return config

def polishing(config):
    """
    Defines the key 'nu' as np.logspace(nu_min, nu_max) and 'ObsCoords' in cartesian coordinates. 
    Removes unused key
    
    Parameters:
    --
    config: :class:`~dictionary` dictionary containing the values needed;  
    """
    removelist = ['Obs_theta', 'Obs_phi', 'nu_Min', 'nu_Max']

    config['nu'] = np.logspace(config['nu_Min'], config['nu_Max'])*u.Hz
    config['ObsCoords'] = cs.spherical_to_cartesian(config['ObserverDistance'], config['Obs_theta'], config['Obs_phi'])

    for key in removelist:
        del config[key]

    config = recap(config)
    return config

def reader():
    """
    Reads the 'config.txt' files and check the main conditions. Returns a dictionary.
    """


    current_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    config_name = 'config.txt'
    config_path = os.path.join(current_dir, config_name)
    save_dir = 'plots/seds'
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
        sys.stderr.write('Error: file "config.txt" not found.\n')
        sys.exit(1)
    except (AttributeError):
        config_file.close()
        sys.stderr.write('Error: "config.txt" format file is not valid.\n')
        sys.stderr.write('Last valid input: ' + key + '\n')
        sys.exit(1)
    except ValueError:
        config_file.close()
        sys.stderr.write('Error: expected two values, found one.\n')
        sys.stderr.write('       Either one or more entries are empty or the format is "variable= number" instead of "variable = number".\n')
        sys.stderr.write('Last valid input: '+ key +'.\n')
        sys.exit(1)

    config_file.close()

    # In case the user want to launch multiple using the same config file.
    if check('id_launch', config) == True:
        try:
            id_run = sys.argv[1]
            fname = open(config['file_list_path'], 'r')
            lines = fname.read().splitlines()
            fname.close()

            id_run = int(id_run)
            id_run = lines[id_run-1]
            config['file_path'] = os.path.join(config['file_path'],id_run)
            
        except IndexError:
            sys.stderr.write('Error: user input not provided. The correct syntax is: python3 main.py i\n')
            sys.exit(1)
        except KeyError:
            sys.stderr.write("Error: 'file_list_path' not provided.\n")
            sys.exit(1)
        except:
            sys.stderr.write("Error: 'id_launch' not correctly stated.\n")                    
            sys.exit(1)

    else:
        try:
            config['file_path'] = str(config['file_path'])
        except KeyError:
            sys.stderr.write("Error, 'file_path' missing.\n")
            sys.exit(1)


    try:
        if config['file_saving_path'] == "''":
            config['file_saving_path'] = os.path.join(current_dir, save_dir)
        else:
            config['file_saving_path'] = str(config['file_saving_path'])
    except KeyError:
            sys.stderr.write("Error: 'file_saving_path' value not specified.\n")
            sys.exit(1)

    floats = ['Z', 'mu_0', 'E_released', 'f_delta', 'nu_Min', 'nu_Max', 'gamma_Min']

    quantities = {
        'pres_Min': 'cm-1 g s-2',
        'pres_Max': 'cm-1 g s-2',
        'dens_Min': 'g cm-3',
        'dens_Max': 'g cm-3',
        'energy_Min': 'erg',
        'energy_Max': 'erg',
        'T_Min': 'K',
        'T_Max': 'K',
        'B': 'G',
        'ObserverDistance': 'cm',
        'Obs_theta': 'deg',
        'Obs_phi': 'deg',
        'DistanceUnit': 'cm',
        'TimeUnit': 's',
        'MassUnit': 'g',
        'GridUnit': 'cm'
    }

    for float_value in floats:
        check_float(config, float_value)

    for key in quantities.keys():
        check_quantity(config, key, quantities[key])

    try:
        if config['n_e'] == 'PowerLaw':
            check_float(config, 'p')
        elif config['n_e'] == 'BrokenPowerLaw':
            check_float(config, 'p1')
            check_float(config, 'p2')
            check_float(config, 'gamma_B')
        else:
            sys.stderr.write("Error: No power-law distribution function has been stated.\n")
            sys.stderr.write("       Available options are 'n_e = PowerLaw' -> 'p = ..' and 'n_e = BrokenPowerLaw' -> 'p1 = ..', 'p2 = ..'\n")
            sys.exit(1)
    except KeyError:
        sys.stderr.write("Error: 'n_e' is missing.\n")
        sys.exit(1)

    if check('id_disc', config) == True:
        disc_reader(config)

    id_list = ['id_cmb', 'id_brem', 'id_syn', 'id_ssc', 'id_ssa', 'id_ec']

    for id_value in id_list:
        check_emissions(id_value, config)

    if config['id_ec']:
        if 'target_list' not in config and config['id_cmb']==False:
            sys.stderr.write("Error: id_ec is set to 'yes' but no emitting region had been declared.\n")
            sys.exit(1)

    config = polishing(config)

    return config

