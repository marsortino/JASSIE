"""
Plots the histograms of the energies temperatures pressure and densities of the blobs.

"""

import lib.misc as m
import lib.plots as plots

import lib.init_plot as init_plot
import lib.indexer as ind

config = init_plot.reader()

data = m.tinymap(config)
dataset = ind.PlotIndex(data, config)


if config['kind'] == 'histogram':
    for label in config['label_list']:
        plots.histogram(config, dataset, label)

elif config['kind'] == 'contour':
    plots.contour(config, dataset)

else:
    print('No correct "kind" stated.\nIn config_plot.txt I have found: kind =', config['kind'],'\nAvailable "kind" settings are "histogram" and "contour".')
