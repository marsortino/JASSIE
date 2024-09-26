## Introduction

<img src="/examples/logo_sfondo_bianco.png" alt="drawing" width="200"/>

**Jetted AGN Spectrum Simulator InterfacE** (**JASSIE**) is a module built on top of the python library [agnpy](https://agnpy.readthedocs.io/en/latest/index.html). This tool is currently in *beta* version.

JASSIE computes the Spectral Energy Distribution (SED) that a FLASH simulated AGN Jet would generate. It allows a certain degree of customability, enabling users to adjust:
* the distance;
* the inclination with respect the line of view;
* which blobs take part in the emission;
* which radiative processes to consider;

![example sed](/examples/run1_S_c.png)

It comes with two scripts:
- main.py
- plotting.py

_main.py_ computes the SED and produces a plot, while _plotting.py_ visualizes blobs distribution using histograms and contour plots.
For more details, please refer to the [User's Manual](usermanual.pdf)

## Requirement and installation

It currently works on **python3.7** or above. The following libraries need to be installed:
* numpy
* matplotlib
* astropy
* h5py
* agnpy
* sherpa
* gammapy

It has been successfully tested on **Ubuntu 23.10**. Note that it does **not** work on **Windows** as _sherpa_ is unsupported. Compability with **Mac OS** is untested.

Before using it, it is necessary to compile the _C code_ found in lib/c_routines/codes using:
> gcc -g - Wall raytrace . c -o raytrace - lm

and then move the output 'raytrace' in lib/c_routines/exec folder.

Alternatevily one can launch _compile.sh_ via terminal making sure to issue the command

> sh compile.sh

Make sure you run this command from the directory where _main.py plotting.py_ and _lib/_ folder are located.
from the same folder in which _main.py plotter.py_ and _lib_ are located.

### Input file structure

JASSIE expects HDF5 files produced via FLASH and looks for the following parameters with these structures:
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

It also expects the time interval between the current and next simulation steps to be labeled _dt_ in the _real runtime parameters_ section.

## How to Use

Both _main.py_ and _plotting.py_ come with a corresponding conffiguration file  (*config.txt* for _main.py_ and *config_plot.txt* for _plotting.py_)

### _main_
To run _main.py_, specify the following in _config.txt_:
* launch settings (file path, number of cores to use)
* output file directory path
    If not specified, the tool will automatically create a directory _plots/seds/_ and store the output there. The output name is generated automatically
* Blobs delimiters
    They can be based on spatial coordinates or on physical parameters. In such case the intervals of temp, density, energy and pressure must be specified.
* Jet physical quantities
* Emitting Regions and Process settings
    Choose whether to consider synchrotron radiation, inverse compton, bremmstrahlung as well as the presence of a disk etc
* Plots settings
     range of frequencies, angle of view, distance of the source, y axis plot interval
* Units of measure and grid units dimensions

An example of _config.txt_ is provided. The [User's Manual](usermanual.pdf) provide an in detail explanation of the features and how to set it to provide maximum results.

### _plotting_

Similarly, in order to launch _plotting.py_ it is necessary to set _config\_plot.txt_.
One must indicate:

* data file path
* output file directory path
    If not specified, JASSIE will automatically create a directory _plots/graph/_ and store the output there. The output name is generated automatically.
* Blobs physical delimiters
    Specify intervals for temperature, density, energy and pressure.
* Graph type
* Units of measure and grid units dimensions

There are two available graphs type: *histogram* and *counter* plots. Specify the desired type in the _config\_plot_ file.

#### _histogram_

![histogram example](/examples/run1_histogram_energy.png)

In this case it is necessary to specify:
* label1 - the physical quantity to consider
    * You can plot multiple histograms by specifying additional labels, e.g.:
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
* x_label - physical quantity for the x-axis
* y_label - physical quantity for the x-axis
* n - number of bins
* cmap - maximum color map value for the color bar

### Output name and output log
JASSIE names output files based on the input file name, with suffixes added according to the emission processes considered. For example:
* _run1.hdf5 → run1_B_Sac_Ecd_t0000_p0000_0000.hdf5_


Suffixes includes:
* Bremsstrahlung -> *_B*
* Synchrotron -> *_S*
* Synchrotron Self Absorption -> *a*
* Synchrotron Self Compton -> *c*
* External Compton -> *_E*
* disk -> *d*
* CMB -> *c*

  
* Angle of view:
  
    theta -> *_tnnnn* where *nnnn* stands for the angle in degree
  
    phi -> *_pnnnn* as above.
  
* a numbered label if the input file with same configuration already exists.

The tool also generates a log file with configuration details, routine execution times, and the final flux value, accessible in the _/log_ directory.
