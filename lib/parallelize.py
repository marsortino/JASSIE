import lib.emission as em
import multiprocessing


def compute_synchro(blocklist, config, q):
    q.put(em.synchrotron(blocklist, config).total_sed)


def synchro(blocklist, config):
    """
    
    """

    num_cores = config['num_cores']
    len_blocklist = len(blocklist)

    block_per_core = len_blocklist // num_cores

    q = multiprocessing.Queue()
    if block_per_core == 0:
        p = multiprocessing.Process(target=compute_synchro, args=(blocklist, config, q,))
        p.start()
        num_cores = 1
    else:
        for i in range(0, num_cores-1):
            p = multiprocessing.Process(target=compute_synchro, args=(blocklist[i*block_per_core : (i+1)*block_per_core], config, q,))
            p.start()
        p = multiprocessing.Process(target=compute_synchro, args=(blocklist[(num_cores-1)*block_per_core:], config, q,))    
        p.start()
        total_synchro_sed = 0
        p.join()

    for i in range(num_cores):
        total_synchro_sed += q.get()
    return total_synchro_sed
        

def blocklist_reducer(blocklist):
    min_blocklist = []

    for block in blocklist:
        if block.EC or block.CMB:
            min_blocklist.append(block)
    return min_blocklist

def compute_EC_targetlist(blocklist, config, q):
    q.put(em.ExtCompton(blocklist, config['nu'], config['target_list'], id_cmb = config['id_cmb']).total_sed)

def compute_EC(blocklist, config, q):
    q.put(em.ExtCompton(blocklist, config['nu'], id_cmb = config['id_cmb']).total_sed)

def external_compton(blocklist, config):
    """
    
    """
    min_blocklist = blocklist_reducer(blocklist)
    num_cores = config['num_cores']
    len_blocklist = len(min_blocklist)

    block_per_core = len_blocklist // num_cores
    q = multiprocessing.Queue()
    total_EC_sed = 0
    print("Starting External Compton computation.")
    if 'target_list' in config:
        if block_per_core == 0:
            p = multiprocessing.Process(target=compute_EC_targetlist, args=(min_blocklist, config, q,))
            p.start()
            num_cores = 1
        else:
            for i in range(0, num_cores-1):
                p = multiprocessing.Process(target=compute_EC_targetlist, args=(min_blocklist[i*block_per_core : (i+1)*block_per_core], config, q,))
                p.start()
            p = multiprocessing.Process(target=compute_EC_targetlist, args=(min_blocklist[(num_cores-1)*block_per_core:], config, q,))    
            p.start()
            p.join()
    else:
        if block_per_core == 0:
            p = multiprocessing.Process(target=compute_EC, args=(min_blocklist, config, q,))
            p.start()
            num_cores = 1
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
        
