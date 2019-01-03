from helper_funcs import *
def compute_af(vcf_file, names_index):

    variants = []
    afs = []

    open_func = file_test(vcf_file)
    with open_func(vcf_file, 'r') as f:
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
# def compute_pi(bed_file, variants, afs, num_seq):
#
#     beds_pi = {}
#     with open(bed_file, 'r') as f:
#         for line in f:
#             # s, e = line.split('\t'[1:])
#             line = line.rstrip('\n')
#             line = line.split('\t')
#             start = int(line[1])
#             end = int(line[2])
#
#             total_pi = 0
#             for i in range(len(variants)):
#                 if start <= variants[i] and variants[i] < end:
#                     pi = 2* afs[i] * (1-afs[i])
#                     total_pi += pi
#                 elif variants[i] >= end:
#                     break
#
#             total_pi_adjusted = (num_seq/(num_seq-1))*total_pi
#             beds_pi[(start, end)] = total_pi_adjusted
#
#     return beds_pi

def compute_pi(bed_file, variants, afs, num_seq):
    vals_in = zip(*(variants, afs))
    out_dict = {}
    beds_pi = {}
    with open(bed_file, 'r') as f:
        for line in f:
            _, s, e = line.split('\t')
            se_list = []
            for i in range(int(s), int(e)):
                out_dict[i] = se_list
            beds_pi[(int(s), int(e))] = se_list

    list(map(lambda v: out_dict.get(int(v[0]), []).append(calculate(float(v[1]))), vals_in))

    return {k: (sum(v))*(num_seq/(num_seq-1)) for (k, v) in beds_pi.items()}, {k: len(v) for (k, v) in beds_pi.items()}


def tabulate_callable_sites_each_window(windows, targets):
    window_callable_sites = {}
    for each_window in windows:
        callable_sites = 0
        for each_target in targets:
            if each_target[0] >= each_window[0] and each_target[0] < each_window[1]:
                if (each_target[1]-1) <= (each_window[1]-1):
                    callable_sites += each_target[1] - each_target[0]
                elif (each_target[1]-1) > (each_window[1]-1):
                    callable_sites += each_window[1] - each_target[0]
            elif each_target[0] < each_window[0] and each_target[1] > each_window[0]:
                if (each_target[1]-1) < (each_window[1]-1):
                    callable_sites += each_target[1] - each_window[0]
                elif (each_target[1]-1) > (each_window[1]-1):
                    callable_sites += each_window[1] - each_window[0]
        window_callable_sites[each_window] = callable_sites
    return window_callable_sites