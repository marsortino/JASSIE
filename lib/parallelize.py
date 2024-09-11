import lib.emission as em
import lib.emission_debug as em_debug
import multiprocessing
import os

def check_condition(num_cores):
    """
    Check the number of cores requested by the user and if they are available.
    If such value is more than those available, the code returns the maximum value possible.

    Works on SLURM and Linux atm.
    """

    try:
        num_cores_available = int(os.environ['SLURM_CPUS_PER_TASK'])
        if num_cores > num_cores_available:
            print("Number of cores issued is more than those available by the machine.")
            print("Number of cores requested: {}. Number of cores available: {}".format(num_cores, num_cores_available))
            print("Continuing with only {} cores.".format(num_cores_available))
            return num_cores_available
        else:
            return num_cores
    except (KeyError):
        num_cores_available = multiprocessing.cpu_count()
        if num_cores > num_cores_available:
            print("Number of cores issued is more than those available by the machine.")
            print("Number of cores requested: {}. Number of cores available: {}".format(num_cores, num_cores_available))
            print("Continuing with only {} cores.".format(num_cores_available))
            return num_cores_available
        else:
            return num_cores
        
def compute_synchro(blocklist, config, q, whole_blocklist):
    """
    Calls a core to compute the synchrotron emission process.
    """
    q.put(em.synchrotron(blocklist, config, whole_blocklist).total_sed)


def synchro(blocklist, config):
    """
    This function parallelize the synchrotron emission computation. 
    If num_cores is == 1, then the code will go on without using multiprocessing.
    If num_cores is != 1, it issues several cores to run calling the compute_synchro() function.

    For each core stated in config.txt, it put in queue a multiprocess.Process instance with a slice of the blocklist to compute.
    Once they all have finished, it is gathered by process 0 which sums the flux and returns the final result.
    """

    num_cores = check_condition(config['num_cores'])
    len_blocklist = len(blocklist)

    block_per_core = len_blocklist // num_cores
    whole_blocklist = blocklist
    if block_per_core == 0 or num_cores == 1:
        total_synchro_sed = em_debug.synchrotron(blocklist, config, blocklist).total_sed
    else:
        q = multiprocessing.Queue()
        for i in range(0, num_cores-1):
            p = multiprocessing.Process(target=compute_synchro, args=(blocklist[i*block_per_core : (i+1)*block_per_core], config, q, whole_blocklist))
            p.start()
        p = multiprocessing.Process(target=compute_synchro, args=(blocklist[(num_cores-1)*block_per_core:], config, q, whole_blocklist))    
        p.start()
        total_synchro_sed = 0
        p.join()

        for i in range(num_cores):
            total_synchro_sed += q.get()
        return total_synchro_sed
        

def blocklist_reducer(blocklist):
    """
    Reduces the list of blocks, returning a list of only those which take part into EC emission.
    This is done in order to properly split the computation load on each core.
    """
    min_blocklist = []

    for block in blocklist:
        if block.EC or block.CMB:
            min_blocklist.append(block)
    return min_blocklist

def compute_EC_targetlist(blocklist, config, q):
    """
    Calls a core to compute the EC emission process.
    """
    q.put(em.ExtCompton(blocklist, config['nu'], config['target_list'], id_cmb = config['id_cmb']).total_sed, config['mu_s'])

def compute_EC(blocklist, config, q):
    """
    Calls a core to compute the EC emission in case there is only the CMB.
    """
    q.put(em.ExtCompton(blocklist, config['nu'], id_cmb = config['id_cmb'], mu_s = config['mu_s']).total_sed)

def external_compton(blocklist, config):
    """
    This function parallelize the External Compton emission computation. 
    If num_cores is == 1, then the code will go on without using multiprocessing.
    If num_cores is != 1, it issues several cores to run calling the compute_synchro() function.

    For each core stated in config.txt, it put in queue a multiprocess.Process instance with a slice of the blocklist to compute.
    Once they all have finished, it is gathered by process 0 which sums the flux and returns the final result.
    """
    min_blocklist = blocklist_reducer(blocklist)
    num_cores = check_condition(config['num_cores'])
    len_blocklist = len(min_blocklist)

    block_per_core = len_blocklist // num_cores
    q = multiprocessing.Queue()
    total_EC_sed = 0
    print("Starting External Compton computation.")
    if 'target_list' in config:
        if block_per_core == 0 or num_cores == 1:
            total_EC_sed += em.ExtCompton(blocklist, config['nu'], config['target_list'], id_cmb = config['id_cmb'], mu_s = config['mu_s']).total_sed
        else:
            for i in range(0, num_cores-1):
                p = multiprocessing.Process(target=compute_EC_targetlist, args=(min_blocklist[i*block_per_core : (i+1)*block_per_core], config, q,))
                p.start()
            p = multiprocessing.Process(target=compute_EC_targetlist, args=(min_blocklist[(num_cores-1)*block_per_core:], config, q,))    
            p.start()
            p.join()

            for i in range(num_cores):
                total_EC_sed += q.get()


    else:
        if block_per_core == 0  or num_cores == 1:
            total_EC_sed += em.ExtCompton(blocklist, config['nu'], id_cmb = config['id_cmb'], mu_s = config['mu_s']).total_sed

        else:
            for i in range(0, num_cores-1):
                p = multiprocessing.Process(target=compute_EC, args=(min_blocklist[i*block_per_core : (i+1)*block_per_core], config, q,))
                p.start()
            p = multiprocessing.Process(target=compute_EC, args=(min_blocklist[(num_cores-1)*block_per_core:], config, q,))    
            p.start()
            p.join()
            
            for i in range(num_cores):
                total_EC_sed += q.get()


    return total_EC_sed
        
