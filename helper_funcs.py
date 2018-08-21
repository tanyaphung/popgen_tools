import collections
import gzip


def file_test(vcf_file):
    if vcf_file.endswith('.gz'):
        return gzip.open
    else:
        return open

def find_index(vcf_file, names_file):
    index = {}

    open_func = file_test(vcf_file)
    with open_func(vcf_file, 'r') as f:
        for l in f:
            line = l.rstrip('\n')
            if line.startswith('#CHROM'):
                items = line.split('\t')
                for i in range(9, len(items)):
                    index[items[i]] = i
                break

    name_index = []
    with open(names_file, 'r') as f:
        for l in f:
            line = l.rstrip('\n')
            name = line.split('\t')[0]
            name_index.append(index[name])
    return name_index

def find_index_all(vcf_file):
    name_index = []
    open_func = file_test(vcf_file)
    with open_func(vcf_file, 'r') as f:
        for l in f:
            line = l.rstrip('\n')
            if line.startswith('#CHROM'):
                items = line.split('\t')
                for i in range(9, len(items)):
                    name_index.append(i)
                break

    return name_index

def count_alt_allele(vcf_file, names_index, bed_file):
    ranges = []
    with open(bed_file, 'r') as f:
        for l in f:
            line = l.rstrip('\n')
            items = line.split('\t')
            ranges.append((int(items[1]), int(items[2])))

    alt_allele_count = []

    open_func = file_test(vcf_file)
    with open_func(vcf_file, 'r') as f:
        for l in f:
            line = l.rstrip('\n')
            if not line.startswith('#'):
                items = line.split('\t')
                if any(lower <= int(items[1]) - 1 < upper for (lower, upper) in ranges):
                    count = 0
                    if items[8] == 'GT':
                        for i in names_index:
                            genotype = items[i]
                            if genotype == '0|1' or genotype == '1|0' or genotype == '0/1' \
                                    or genotype == '1/0':
                                count += 1
                            if genotype == '1|1' or genotype == '1/1':
                                count += 2
                        alt_allele_count.append(count)

                    else:
                        for i in names_index:
                            genotype = items[i].split(':')[0]
                            if genotype == '0|1' or genotype == '1|0' or genotype == '0/1' \
                                    or genotype == '1/0':
                                count += 1
                            if genotype == '1|1' or genotype == '1/1':
                                count += 2
                        alt_allele_count.append(count)
    return alt_allele_count

def count_alt_allele_no_target_bed(vcf_file, names_index):

    alt_allele_count = []

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
                        if genotype == '0|1' or genotype == '1|0' or genotype == '0/1' \
                                or genotype == '1/0':
                            count += 1
                        if genotype == '1|1' or genotype == '1/1':
                            count += 2
                    alt_allele_count.append(count)

                else:
                    for i in names_index:
                        genotype = items[i].split(':')[0]
                        if genotype == '0|1' or genotype == '1|0' or genotype == '0/1' \
                                or genotype == '1/0':
                            count += 1
                        if genotype == '1|1' or genotype == '1/1':
                            count += 2
                    alt_allele_count.append(count)
    return alt_allele_count


def make_sfs(num_seq, alt_allele_count):
    """
    This function makes a site frequency spectrum following equation 1.2 in John Wakely's
    Coalescent book.

    :param num_seq: number of sequences in the VCF file.
    :param alt_allele_count: a list where each item in the list is the count of the
    alternate alleles for each variant
    :return: a dictionary where keys are frequency and values are the count of variants
    at that frequency
    """

    # frequency_count is a dictionary where keys indicate the frequency and values
    # indicate how many variants there are at that particular frequency.
    frequency_count = collections.Counter(alt_allele_count)

    # Generate a folded site frequency spectrum
    sfs = {}
    for i in range(1, (num_seq / 2) + 1):
        if i != (num_seq - i):
            count = float(frequency_count[i] + frequency_count[num_seq - i])
            sfs[i] = int(count)
        else:
            count = float(frequency_count[i] + frequency_count[num_seq - i]) / 2
            sfs[i] = int(count)
    return sfs


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
                    if not any(i == './.' for i in genotypes):
                        for genotype in genotypes:
                    # for i in names_index:
                    #     genotype = items[i].split(':')[0]
                            if genotype == '0|1' or genotype == '1|0' or genotype == '0/1' or genotype == '1/0':
                                count += 1
                            if genotype == '1|1' or genotype == '1/1':
                                count += 2
                    variants.append(int(items[1]) - 1)
                    afs.append(float(count) / (len(names_index)*2))
    return variants, afs

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
