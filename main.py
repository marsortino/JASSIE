import numpy as np
from agnpy.utils.plot import plot_sed, load_mpl_rc
import matplotlib.pyplot as plt

import time

import lib.init as init
import lib.indexer as ind
import lib.misc as m

import lib.listbuilder as listbuilder
import lib.energy as e
import lib.rayline as rl
import lib.objects as o
import lib.emission as em
import lib.parallelize as pl

# //*************************************
# DATA READER
sets = init.reader()
time.sleep(2)

start_time = time.time()
data = m.datamap(sets)
end_time = time.time()
print('Data read in', "{:.2f}".format(end_time-start_time), 'seconds.')

# ************************************//

# # SPEED UP AND SCALABILITY TEST PURPOSES
speed_up_test = False
if speed_up_test:
    import os
    timing_file_pathfile = os.path.dirname(os.path.abspath(__file__))
    timing_file_pathfile = os.path.join(timing_file_pathfile, sets['name'] + '_n' + str(sets['num_cores']) + '_log_timer.txt')

    if os.path.exists(timing_file_pathfile):
        ftime = open(timing_file_pathfile, 'a')
    else:
        ftime = open(timing_file_pathfile, 'x')

    ftime.write("num_cores: " + str(sets['num_cores']) + "\n")
# ************************************//

time_table = {}
initial_time = time.time()

# //*************************************
# INDEX COMPILER
start_time = time.time()
index = ind.MainIndex(data, sets)
end_time = time.time()
print('Index complete in', "{:.2f}".format(end_time-start_time), 'seconds.')
time_table['index'] = end_time-start_time
# *************************************//

# //*************************************
# LIST BUILDER
start_time = time.time()
builder = listbuilder.builder(sets)
pointlist = builder.pointsbuilder(index, data.coordSet)

index_positive_blocks = np.where(pointlist[:,2]>0)[0]
index_negative_blocks = np.where(pointlist[:,2]<0)[0]

positivelist = pointlist[index_positive_blocks]
negativelist = pointlist[index_negative_blocks]

blocklist = builder.blockbuilder(index, data)
end_time = time.time()
print('Blocks built in', "{:.2f}".format(end_time-start_time), 'seconds.')
time_table['blocks'] = end_time-start_time
if speed_up_test:
    ftime.write("num_cores: " + str(len(blocklist)) + "\n")

# *************************************//

# //*************************************
# QUICKHULL
start_time = time.time()
qh_checker = 0

if sets['id_cmb']:
    m.sphere(positivelist, blocklist, sets['qhull_depth'])
    m.sphere(negativelist, blocklist, sets['qhull_depth'])
    qh_checker += 1

if sets['id_disc'] == 'y' or sets['id_disc'] == 'yes':
    m.lower_sphere(positivelist, blocklist, sets['qhull_depth']) 
    m.lower_sphere(negativelist, blocklist, sets['qhull_depth'])
    qh_checker += 1

end_time = time.time()

if qh_checker != 0:
    print("QuickHull in:", "{:.2f}".format(end_time-start_time), 'seconds.')
    time_table['quickhull'] = end_time-start_time
else:
    time_table['quickhull'] = 0
# *************************************//

# //*************************************
# RAYTRACING
start_time = time.time()
m.BlocklistTXT(blocklist, sets['name'])
src = o.target('source', np.array([0.0,0.0,0.0]))
obs = o.target('observer', sets['ObsCoords'])
rl.C_raytracing((src , obs), blocklist, sets['name'])
end_time = time.time()
print('Raytracing complete in', "{:.2f}".format(end_time-start_time), 'seconds.')
time_table['raytracing'] = end_time-start_time

# *************************************

# //*************************************
# ENERGY BLOCKS 
start_time = time.time()
energies = e.Block_Energy(sets)
if 'target_list' in sets:
    energies.EnergyComputer(blocklist, data.dtmin, sets['target_list'], sets['id_cmb'])
else:
    energies.EnergyComputer(blocklist, data.dtmin, cmb = sets['id_cmb'])
end_time = time.time()
print('Energy computed in', "{:.2f}".format(end_time-start_time), 'seconds.')
time_table['energy'] = end_time-start_time

# *************************************//

# //*************************************
# BLOBLIST BUILDER
start_time = time.time()
builder.blobbuilder(blocklist)
end_time = time.time()
print('Bloblist built in', "{:.2f}".format(end_time-start_time), 'seconds.')
time_table['blobs_built'] = end_time-start_time

# *************************************//


##  SED
sed_per_emission_process = {}

if sets['id_brem']:
    # //*************************************
    # BREMMSTRAHLUNG
    start_time = time.time()
    sed_value = em.Bremsstrahlung(blocklist, sets).total_sed
    sed_per_emission_process['id_brem'] = sed_value
    end_time = time.time()
    print('Bremmstrahlung computed in', "{:.2f}".format(end_time-start_time), 'seconds.')
    time_table['Bremmstrahlung'] = end_time-start_time

    # *************************************//


if sets['id_syn']:
#     # //*************************************
#     # SYNCHROTRON    
    start_time = time.time()
    sed_value = pl.synchro(blocklist, sets)
    sed_per_emission_process['id_syn'] = sed_value
    end_time = time.time()
    print('Synchrotron computed in',"{:.2f}".format(end_time-start_time), 'seconds.')
    time_table['Synchrotron'] = end_time-start_time

    if speed_up_test:
        ftime.write('SYN: ' + str(end_time-start_time) + '\n')

#     # *************************************//

if sets['id_ec']:
    # //*************************************
    # EXTERNAL COMPTON
    start_time = time.time()
    sed_value = pl.external_compton(blocklist, sets)
    sed_per_emission_process['id_ec'] = sed_value
    end_time = time.time()
    print('External Compton computed in', "{:.2f}".format(end_time-start_time), 'seconds.')
    time_table['External_Compton'] = end_time-start_time
    if speed_up_test:
        ftime.write('ExC: ' + str(end_time-start_time))

# #     # *************************************//
if speed_up_test:
    ftime.close()

## PLOTTING
load_mpl_rc()

final_sed = 0

for key in sed_per_emission_process.keys():
    plot_sed(sets['nu'], sed_per_emission_process[key], label=key[3:])
    final_sed += sed_per_emission_process[key]
    print(final_sed)
# SED taking into account emission process
plot_sed(sets['nu'], final_sed, label='total flux')
if sets['yplot_max'] == 'auto':
    y_max = (max(final_sed)*10).value
y_min = sets['yplot_min']
if y_max < y_min:
    print('Error, the minimum threshold provided for the output is too high.\n Max value found is higher than threshold value.')
    print('Max value of the flux: ', y_max)
    print('Lower threshold provided: ', y_min)
    print('Manually putting as lower threshold Max*1e-15: ', y_max*1e-15)
    y_min = y_max*1e-15
plt.ylim([y_min, y_max])
plt.title(sets['name'])
plt.savefig(sets['file_saving_path']+'/'+sets['name']+'.png')
plt.close()
# SED with only the final flux
plot_sed(sets['nu'], final_sed, label='total_flux')
plt.ylim([y_min, y_max])
plt.title(sets['name'] + ' final')
plt.savefig(sets['file_saving_path']+'/'+sets['name']+'_final.png')
plt.close()


# Final Output
total_time = time.time()
elapsed_time = total_time - initial_time
minutes = int(elapsed_time // 60)
seconds = int(elapsed_time % 60)

time_table['final_time'] = elapsed_time
m.write_log(sets, time_table, len(blocklist), final_sed)

m.clean(sets['name'])

print("Done!")
print('In: ', minutes, 'minutes and', seconds, 'seconds.')

