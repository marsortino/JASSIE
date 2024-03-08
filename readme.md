## Introduction



This tool computes the Spectral Energy Distribution that a FLASH simulated AGN Jet would generate. It allows a certain degree of customability, providing the user the possibility to decide:
* the distance;
* the inclination with respect the line of view;
* which blobs take part in the emission;
* which radiative processes to consider;

![example sed](/examples/run1_S_c.png)

It comes with two codes:
- main.py
- plotting.py

Where the first one computes the SED and produces a plot, the latter can be used to show blobs distribution via histograms and contour plots.

For further informations check the [User's Manual](usermanual.pdf)

## Requirement and installation

It currently works on **python3.7** or above. The following libraries need to be installed:
* numpy
* matplotlib
* astropy
* agnpy
* sherpa
* gammapy

It has been correctly ran and tested on **Ubuntu 23.10**. It does **not** work on **Windows** since _sherpa_ is not supported. I was not able to test it on **Mac** so I can't provide any info on that.

Before using it, it is necessary to compile the _C code_ found in lib/c_routines/codes using:
> gcc -g - Wall raytrace . c -o raytrace - lm

and then move the output 'raytrace' in lib/c_routines/exec folder.

Alternatevily one can launch _compile.sh_ via terminal making sure to issue the command

> sh compile.sh

from the same folder in which _main.py plotter.py_ and _lib_ are located.

### Input file structure

The tool expects HDF5 files produced via FLASH. It looks for the following parameters with the respective structure:
* temp
    * 8x8x8
* dens
    * 8x8x8
* pres
    * 8x8x8
* ener
    * 8x8x8
* velx
    * 8x8x8
* vely
    * 8x8x8
* velz
    * 8x8x8
* block size
    * 1x3
* coordinates
    * 1x3

It also looks for the time between the current input simulation time and the next one in _real runtime parameters_ section and it assumes to be labeled as _dt_.

## How to Use

Both _main.py_ and _plotting.py_ come with a txt file which serves the purpose of providing the launching conditions, respectively *config.txt* and *config_plot.txt*.

### _main_
In order to run _main.py_, one must indicate in _config.txt_:
* data file path
* output file directory path
    If no directory is indicated, the tool will automatically create a directory _plots/seds/_ and store the output there. The output name is automatically produced by the tool.
* Blobs physical delimiters
    The intervals of temp, density, energy and pressure must be specified. Only blobs inside the interval will be selected    
* Jet physical quantities
* Emitting Regions and Process settings
    User can choose whether to consider synchrotron radiation, inverse compton, bremmstrahlung as well as the presence of a disk etc
* Units of measure and grid units dimensions

An example of _config.txt_ is provided. The [User's Manual](usermanual.pdf) provide an in detail explanation of the features and how to set it to provide maximum results.

### _plotting_

Similarly, in order to launch _plotting.py_ it is necessary to set _config\_plot.txt_.
One must indicate:

* data file path
* output file directory path
    If no directory is indicated, the tool will automatically create a directory _plots/graph/_ and store the output there. The output name is automatically produced by the tool.
* Blobs physical delimiters
    The intervals of temp, density, energy and pressure must be specified. Only blobs inside the interval will be selected    
* Graph Kind
* Units of measure and grid units dimensions

There are two graphs structure available: *histogram* and *counter* plots. In order to launch one or the other, the _config\_plot_ structure must be changed.

#### _histogram_

![histogram example](/examples/run1_histogram_energy.png)

In this case it is necessary to specify:
* label1 - the physical quantity to consider
    * The user can indicate to compute more than 1 histogram per once. This can be done by entering more than one parameter. i.e.:
    > label1 = temp
    >
    > label2 = energy 
    >
    > label3 = pres
    
    etc. 
* bin_rule - the rule to use when determining the number of bins.
    * Three different options available: *sqrt*, *rice*, *sturges*. More info on [User's Manual](usermanual.pdf).

#### _counter_

![counter example](/examples/run1_contour__temperature_energy.png)

In this case it is necessary to specify:
* x_label - which physical quantity to be put on the x-axis
* y_label - which physical quantity to be put on the y-axis
* n - number of bins
* cmap - max value of the color map to display on the colour bar
