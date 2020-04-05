import argparse
import pandas as pd
import gzip
import collections
import time

#TODO: Add error checking for inputs
#TODO: When generating an SFS, allow for the requirement of target bed to be optional
#TODO: Have a check so that the chrName is concordant between target and window

def parse_args():
    """
    Parsing command-line arguments
    """
    parser = argparse.ArgumentParser(description='This program: (1) generates a 1-D site frequency '
                                                 'spectrum, and (2) calculates genetic diversity, pi. Pi is defined as '
                                                 'the average number of DNA differences between all pairs of sequence).')

    parser.add_argument('--vcf_file', required=True,
                        help='REQUIRED. Input the path to a VCF file. Either gzipped or '
                             'not gzipped file works.')

    parser.add_argument('--names_list', required=False,
                        help='If you only want to generate the SFS and/or calculate genetic diversity for a subset of '
                             'individuals from the VCF file, input the path to that file here. This file lists each '
                             'individual on each line. If you want to perform the calculations for all of the individuals,'
                             'do not use this flag.')

    # Options for calculating pi
    parser.add_argument('--pi', action='store_true', default=False,
                        help='Turn on this flag if you want to calculate pi.')

    parser.add_argument('--ploidy', help='Input either haploid or diploid')

    parser.add_argument('--pi_all', action='store_true', default=False,
                        help='Turn on this flag if you want to calculate pi using all of the sites in the VCF file.')

    parser.add_argument('--pi_target', action='store_true', default=False,
                        help='Turn on this flag if you want to calculate pi in regions of the genome specified by a BED file.')

    parser.add_argument('--pi_window', action='store_true', default=False,
                        help='Turn on this flag if you want to calculate pi in windows in the presence of target_bed file.')

    # Options for calculating SFS
    parser.add_argument('--sfs_all', action='store_true', default=False,
                        help='Turn on this flag if you want to calculate sfs using all of the sites in the VCF file.')

    # Extra required files for different options of calculating pi
    parser.add_argument('--target_bed', required=False,
                        help='If you want to calculate pi in regions of the genome specified by a BED file, input the '
                             'path here.')

    parser.add_argument('--window_bed', required=False,
                        help='If you want to calculate genetic diversity in window in the presence of target_bed file, '
                             'input the path here.')

    # Output files
    parser.add_argument('--sfs_all_out', required=False,
                        help='Input the file to the output of sfs_all.')

    parser.add_argument('--pi_target_out', required=False,
                        help='Path to the output pi_target_out.')

    parser.add_argument('--pi_target_per_site_out', required=False,
                        help='Path to the output pi_target_per_site_out.')

    args = parser.parse_args()
    return args

def file_test(vcf_file):
    """
    This function checks if the input VCF file is gzip or not.
    """
    if vcf_file.endswith('.gz'):
        return gzip.open, 'rt'
    else:
        return open, 'r'

def find_index(vcf_file, names_list=None):
    """
    This function finds the index of the individuals where you want to calculate pi or SFS on.
    :param vcf_file:
    :param names_list:
    :return: a list where each item in the list is the index of the individual.
    """

    name_index = []

    open_func, mode = file_test(vcf_file)
    with open_func(vcf_file, mode) as f:
        for line in f:
            if line.startswith('#CHROM'):
                items = line.rstrip('\n').split('\t')
                if names_list is None:
                    name_index = [i for i in range(9, len(items)) if names_list is None]
                    break
                else:
                    index = {items[i] : i for i in range(9, len(items))}
                    with open(names_list, 'r') as f:
                        individuals = [l.rstrip('\n') for l in f]
                        name_index = [index[i] for i in individuals]
                    break
    return name_index

def compute_af(vcf_file, names_index, ploidy):

    variants = []
    afs = []

    open_func, mode = file_test(vcf_file)
    with open_func(vcf_file, mode) as f:
        for l in f:
            line = l.rstrip('\n')
            if not line.startswith('#'):
                items = line.split('\t')
                count = 0
                if items[8] == 'GT':
                    for i in names_index:
                        genotype = items[i]
                        if genotype == '0|1' or genotype == '1|0' or genotype == '0/1' or genotype == '1/0' or genotype == '1':
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
                            if genotype == '0|1' or genotype == '1|0' or genotype == '0/1' or genotype == '1/0' or genotype == '1':
                                count += 1
                            if genotype == '1|1' or genotype == '1/1':
                                count += 2
                    variants.append(int(items[1]) - 1)
                    if ploidy == "diploid":
                        afs.append(float(count) / (len(names_index)*2))
                    elif ploidy == "haploid":
                        afs.append(float(count) / len(names_index))
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


def count_alt_allele_from_row(line, indices, ploidy):
    if line.startswith('#'):
        return
    pop_subset = [line.split('\t')[i] for i in indices]
    if ploidy == 'diploid':
        return sum([int(s[0]) + int(s[2]) for s in pop_subset if s[0] != '.' and s[2] != '.' and s[0] != '/' and s[2] != '/' and s[0] != '|' and s[2] != '|'])
    elif ploidy == 'haploid':
        return sum([int(s[0]) for s in pop_subset if s[0] != '.'])


def count_alt_allele_all(vcf, names_index, ploidy):
    open_func, mode = file_test(vcf)
    with open_func(vcf, mode) as f:
        alt_allele_count = list(map(lambda l: count_alt_allele_from_row(l, names_index, ploidy), f))
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
    frequency_count = collections.Counter(alt_allele_count)


    sfs = {i: frequency_count[i] + frequency_count[num_seq-i] for i in range(1, int(num_seq/2) + 1)}
    sfs[int(num_seq/2)] = int(sfs[int(num_seq/2)]/2)

    return sfs


def main():

    args = parse_args()

    # Check if the user wants to perform the calculations for all of the individuals in the VCF or just a subset of
    # individuals, which is specified by the presence of the argument names_list
    if args.names_list:
        # Find the index of the individuals whose genetic diversity you want to calculate
        names_index = find_index(args.vcf_file, args.names_list)

    else:
        names_index = find_index(args.vcf_file)

    if args.pi:
        # For pi, get the variants and afs
        variants_afs = compute_af(args.vcf_file, names_index, args.ploidy)

    # Check if the user wants to calculate pi using all of the sites in the VCF file.
    if args.pi_all:
        # Calculate allele frequency
        afs = variants_afs[1]
        # Compute pi_all
        if args.ploidy == "diploid":
            pi_all = compute_pi_all(afs, len(names_index)*2)
            with open('pi_all.out', 'w') as out:
                out.write('pi_all' + '\n' + str(pi_all))

            print ('Pi for all of the variants in this population is ', str(pi_all))

        elif args.ploidy == "haploid":
            pi_all = compute_pi_all(afs, len(names_index))
            with open('pi_all.out', 'w') as out:
                out.write('pi_all' + '\n' + str(pi_all))

            print ('Pi for all of the variants in this population is ', str(pi_all))




    # Check if the user wants to calculate pi in regions of the genome specified by the BED file
    if args.pi_target:
        targets = []
        with open(args.target_bed, 'r') as f:
            for line in f:
                chrName, s, e = line.split('\t')
                targets.append((chrName, int(s), int(e)))
        # Calculate adjusted pi
        if args.ploidy == "diploid":
            adjusted_pi = compute_pi_target(targets, variants_afs[0], variants_afs[1], len(names_index) * 2)
        elif args.ploidy == "haploid":
            adjusted_pi = compute_pi_target(targets, variants_afs[0], variants_afs[1], len(names_index))
        pi_outfile = open(args.pi_target_out, 'w')
        header = ['chr', 'start', 'end', 'pi']
        pi_outfile.write('\t'.join(header) + '\n')
        for k, v in adjusted_pi[0].items():
            out = [str(k[0]), str(k[1]), str(k[2]), str(v)]
            pi_outfile.write('\t'.join(out) + '\n')
        pi_outfile.close()

        # Calculate pi per site
        pi_target = pd.read_table(args.pi_target_out)
        pi_total = pi_target['pi'].sum()
        callable_sites = (pi_target['end'] - pi_target['start']).sum()
        pi_per_site = pi_total/callable_sites

        with open(args.pi_target_per_site_out, 'w') as out:
            out.write('pi_per_site' + '\n' + str(pi_per_site))

    # Check if the user wants to calculate pi in regions of the genome specified by the BED file but also in windows
    if args.pi_window:
        targets = []
        windows = []
        with open(args.target_bed, 'r') as f:
            for line in f:
                chrName, s, e = line.split('\t')
                targets.append((chrName, int(s), int(e)))
        with open(args.window_bed, 'r') as f:
            for line in f:
                _, s, e = line.split('\t')
                windows.append((int(s), int(e)))

        # Make new targets
        new_targets = place_target_into_window(windows, targets)

        # Calculate adjusted pi
        adjusted_pi = compute_pi_target(new_targets, variants_afs[0], variants_afs[1], len(names_index) * 2)

        pi = []
        for k, v in adjusted_pi[0].items():
            pi.append(v)

        results = calc_total_sites_pi_per_window(windows, new_targets, pi)

        pi_outfile = open("pi_window.out", 'w')
        for window in results[0]:
            out = [str(window[0]), str(window[1]), str(results[0][window]), str(results[1][window])]
            pi_outfile.write('\t'.join(out) + '\n')
        pi_outfile.close()

    # Check if the user wants to calculate the sfs using all of the sites in the vcf file
    if args.sfs_all:

        s = time.time()

        # Calculate the number of alternate alleles for each variant.
        alt_allele_count_all = count_alt_allele_all(args.vcf_file, names_index, args.ploidy)

        # Generate a folded site frequency spectrum
        if args.ploidy == 'diploid':
            sfs = make_sfs(len(names_index) * 2, alt_allele_count_all)
        elif args.ploidy == 'haploid':
            sfs = make_sfs(len(names_index), alt_allele_count_all)

        sfs_outfile = open(args.sfs_all_out, "w")
        header = ['af_bin', 'count']
        sfs_outfile.write('\t'.join(header) + '\n')
        out = [(str(k) , str(sfs[k])) for k in sfs]
        sfs_outfile.write('\n'.join(['\t'.join(s) for s in out]))

        e = time.time()
        print ("Time %d ms" %int(e-s))

if __name__ == '__main__':
    main()
