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

initial_time = time.time()

# //*************************************
# DATA READER
sets = init.reader()
time.sleep(2)

start_time = time.time()
data = m.datamap(sets)
end_time = time.time()
print('Data read in', "{:.2f}".format(end_time-start_time), 'seconds.')
# ************************************//

# //*************************************
# INDEX COMPILER
start_time = time.time()
index = ind.MainIndex(data, sets)
end_time = time.time()
print('Index complete in', "{:.2f}".format(end_time-start_time), 'seconds.')
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
# *************************************//

# //*************************************
# QUICKHULL
start_time = time.time()
qh_checker = 0
if sets['id_cmb']:
    m.sphere(positivelist, blocklist)
    m.sphere(negativelist, blocklist)
    qh_checker += 1

if sets['id_disc'] == 'y' or sets['id_disc'] == 'yes':
    m.lower_sphere(positivelist, blocklist) 
    m.lower_sphere(negativelist, blocklist)
    qh_checker += 1

end_time = time.time()

if qh_checker != 0:
    print("QH in:", "{:.2f}".format(end_time-start_time), 'seconds.')
# *************************************//

# //*************************************
# RAYTRACING
start_time = time.time()
m.BlocklistTXT(blocklist)
src = o.target('source', np.array([0.0,0.0,0.0]))
obs = o.target('observer', sets['ObsCoords'])
rl.C_raytracing(src, blocklist)
rl.C_raytracing(obs, blocklist)
end_time = time.time()   
print('Raytracing complete in', "{:.2f}".format(end_time-start_time), 'seconds.')
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
# *************************************//

# //*************************************
# BLOBLIST BUILDER
start_time = time.time()
builder.blobbuilder(blocklist)
end_time = time.time()
print('Bloblist built in', "{:.2f}".format(end_time-start_time), 'seconds.')
# *************************************//


##  SED
seds = []

if sets['id_brem']:
    # //*************************************
    # BREMMSTRAHLUNG
    start_time = time.time()
    seds.append(em.Bremsstrahlung(blocklist, sets))
    end_time = time.time()
    print(end_time-start_time)
    print('Bremmstrahlung computed in', "{:.2f}".format(end_time-start_time), 'seconds.')
    # *************************************//

if sets['id_syn']:
    # //*************************************
    # SYNCHROTRON
    start_time = time.time()
    seds.append(em.synchrotron(blocklist, sets))
    end_time = time.time()
    print(end_time-start_time)
    print('Synchrotron computed in',"{:.2f}".format(end_time-start_time), 'seconds.')
    # *************************************//

if sets['id_ec']:
    # //*************************************
    # EXTERNAL COMPTON
    start_time = time.time()
    if 'target_list' in sets:
        seds.append(em.ExtCompton(blocklist, sets['nu'], sets['target_list'], id_cmb = sets['id_cmb']))
    else:
        seds.append(em.ExtCompton(blocklist, sets['nu'], id_cmb = sets['id_cmb']))
    end_time = time.time()
    print('External Compton computed in', "{:.2f}".format(end_time-start_time), 'seconds.')
    # *************************************//
total_sed = 0
for sed in seds:
    total_sed += sed.total_sed

# //*************************************
## PLOTTING
load_mpl_rc()

plot_sed(sets['nu'], total_sed)
plt.title(sets['name'])
#plt.xlim([1e9, 1e17])
plt.ylim([1e-5, 1e5])
plt.savefig(sets['file_saving_path']+'/'+sets['name']+'.png')
total_time = time.time()
print('Done in: ', "{:.2f}".format(total_time-initial_time), 'seconds.')

elapsed_time = total_time - initial_time
minutes = int(elapsed_time // 60)
seconds = int(elapsed_time % 60)

print('i.e.:', minutes, 'minutes and', seconds, 'seconds.')
