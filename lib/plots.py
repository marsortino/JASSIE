import sys
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from lib.init_plot import get_parameters

def histogram(config, data, label):
    """
    Computes the histograms of the labeled quantity.
    
    Parameters:
    config (dictionary): dictionary containing the values of config_plot.txt
    data (misc.datamap): data from the input file
    label (str): quantity to compute histogram on.
    """

    x_min, x_max = get_parameters(config['list_of_parameters'], label)
    label = label.strip("'")
    if label =='density':
        quantity = data.density
    elif label =='pressure':
        quantity = data.pressure
    elif label =='temperature':
        quantity = data.temp
    elif label =='energy':
        quantity = data.energy
    else:
        sys.stderr.write('Error in labeling.. I cannot find anything')
        sys.exit()

    if config['bin_rule'] == 'sqrt':
        n = round(np.power(len(quantity),1/2))
    elif config['bin_rule'] == 'rice':
        n = 2*round(np.power(len(quantity), 1/3))
    elif config['bin_rule'] == 'sturges':
        n = round(np.log2(len(quantity)))+1
    
    try:
        display_arrays = str(config['display_arrays'])
        if display_arrays in ['yes', 'y']:
            print(quantity)
    except KeyError:
        print('display_arrays not found.')
    log_bin_edges = np.logspace(np.log10(x_min.value), np.log10(x_max.value), n)
    counts, bins = np.histogram(quantity, bins = log_bin_edges)
    norm_counts = np.zeros(len(counts))
    for i in range(len(bins)-1):
        norm_counts[i] = counts[i]/(np.sum(counts)*(bins[i+1]-bins[i]))
    norm_counts = norm_counts/np.sum(norm_counts)

    plt.bar((bins[1:]+bins[:-1])/2, norm_counts, width=np.diff(bins))

    try:
        xscale = config['xscale']
    except KeyError:
        print('xscale not stated. Assuming xscale = log')
        xscale = 'log'

    try:
        yscale = config['yscale']
    except KeyError:
        print('yscale not stated. Assuming yscale = linear')
        yscale = 'linear'

    name = config['name'][label]

    plt.xscale(xscale)
    plt.yscale(yscale)
    plt.xticks()
    plt.xlabel(label)
    plt.ylabel('Frequency')
    plt.title(label + ' histogram')
    plt.savefig(config['file_saving_path']+'/'+ name +'.png')
    plt.close()

def contour(config, data):#x, y, n, xlabel = '', ylabel = ''):
    """
    Produce the heatmaps of x and y with a number of bins n 
    Parameters:
    x: astropy quantity
    y: astropy quantity
    n: int, number of bins to produce
    xlabel: string, label of x axis
    ylabel: string, label of y axis
    """
    labels = [config['x_label'], config['y_label']]

    quantity = [] 
    print('plotting:')
    for label in labels:
        if label =='density':
            print('--density')
            quantity.append(data.density)                  
        elif label =='pressure':
            print('--pressure')
            quantity.append(data.pressure)
        elif label =='temperature':
            print('--temperature')
            quantity.append(data.temp)
        elif label =='energy':
            print('--energy')
            quantity.append(data.energy)
        else:
            print('can t find anything')
            exit()
    
    try:
        display_arrays = str(config['display_arrays'])
        if display_arrays in ['yes', 'y']:
            print(quantity)
    except KeyError:
        print('display_arrays not found.')
        
    n = config['n']

    height, xedges, yedges = np.histogram2d(quantity[0], quantity[1], n)

    if not config['clim']:
        config['clim'] = len(quantity[0])

    X = np.logspace(np.log10(xedges[0]), np.log10(xedges[-1]), n )
    Y = np.logspace(np.log10(yedges[0]), np.log10(yedges[-1]), n )

    X, Y = np.meshgrid(X, Y)
    cp = plt.contourf(X, Y, height.T)
    cp.set_clim(0, config['clim'])
    cbar = plt.colorbar(cp)

    x_label = config['x_label']
    y_label = config['y_label']

    if x_label == 'density':
        plt.xlabel('density [g cm-3]')
    elif x_label =='pressure':
        plt.xlabel('pressure [g cm-2 s-2]')
    elif x_label == 'energy':
        plt.xlabel('energy [erg]')
    elif x_label == 'temperature':
        plt.xlabel('temperature [K]')
    else:
        plt.xlabel(x_label)

    if y_label == 'pressure':
        plt.ylabel('pressure [g cm-2 s-2]')
    elif y_label == 'density':
        plt.ylabel('density [g cm-3]')
    elif y_label == 'energy':
        plt.ylabel('energy [erg]')
    elif y_label == 'temperature':
        plt.ylabel('temperature [K]')
    else:
        plt.ylabel(y_label)

    plt.yscale('log')
    plt.xscale('log')
    plt.savefig(config['file_saving_path']+'/'+config['name']+'.png')
    plt.close()