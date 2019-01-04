from helper_funcs import *

def compute_af(vcf_file, names_index):

    variants = []
    afs = []

    file_test_result = file_test(vcf_file)
    open_func = file_test_result[0]
    mode = file_test_result[1]
    with open_func(vcf_file, mode) as f:
        for l in f:
            line = l.rstrip('\n')
            if not line.startswith('#'):
                items = line.split('\t')
                count = 0
                if items[8] == 'GT':
                    for i in names_index:
                        genotype = items[i]
                        if genotype == '0|1' or genotype == '1|0' or genotype == '0/1' or genotype == '1/0':
                            count += 1
                        if genotype == '1|1' or genotype == '1/1':
                            count += 2
                    variants.append(int(items[1])-1)
                    afs.append(float(count) / (len(names_index)*2))

                else:
                    genotypes = [items[i].split(':')[0] for i in names_index]
                    if not any(i == './.' for i in genotypes): #TODO: think about how to deal with missing data
                        for i in names_index:
                            genotype = items[i].split(':')[0]
                            if genotype == '0|1' or genotype == '1|0' or genotype == '0/1' or genotype == '1/0':
                                count += 1
                            if genotype == '1|1' or genotype == '1/1':
                                count += 2
                    variants.append(int(items[1]) - 1)
                    afs.append(float(count) / (len(names_index)*2))
    return variants, afs

def compute_pi_all(afs, num_seq):
    total_pi = 0
    for i in afs:
        pi = 2 * i * (1 - i)
        total_pi += pi
    total_pi_adjusted = (num_seq / (num_seq - 1)) * total_pi
    return total_pi_adjusted

def calculate(num):
    return 2*num*(1-num)

def compute_pi_target(targets, variants, afs, num_seq):
    vals_in = zip(*(variants, afs))
    out_dict = {}
    beds_pi = {}
    # with open(bed_file, 'r') as f:
    #     for line in f:
    #         chrName, s, e = line.split('\t')
    #         se_list = []
    for target in targets:
        se_list = []
        for i in range(target[1], target[2]):
            out_dict[i] = se_list
        beds_pi[(target[0], target[1], target[2])] = se_list

    list(map(lambda v: out_dict.get(int(v[0]), []).append(calculate(float(v[1]))), vals_in))

    return {k: (sum(v))*(num_seq/(num_seq-1)) for (k, v) in beds_pi.items()}, {k: len(v) for (k, v) in beds_pi.items()}

def place_target_into_window(windows, targets):
    new_targets = []
    for window in windows:
        for target in targets:
            if target[1] >= window[0] and target[1] < window[1]:
                if target[2] <= window[1]:
                    new_targets.append(target)
                else:
                    new_targets.append((target[0], target[1], window[1]))
            elif target[1] < window[0] and target[2] > window[0]:
                if target[2] < window[1]:
                    new_targets.append((target[0], window[0], target[2]))
                elif (target[2]-1) > (window[1]-1):
                    new_targets.append((target[0], window[0], window[1]))
    return new_targets

def calc_total_sites_pi_per_window(windows, new_targets, pi):
    windows_total_sites = {}
    windows_pi = {}
    for window in windows:
        total_sites = 0
        total_pi = 0
        for i in range(len(new_targets)):
            if new_targets[i][1] >= window[0] and new_targets[i][2] <= window[1]:
                total_sites += new_targets[i][2] - new_targets[i][1]
                total_pi += pi[i]
        windows_total_sites[window] = total_sites
        windows_pi[window] = total_pi
    return windows_total_sites, windows_pi