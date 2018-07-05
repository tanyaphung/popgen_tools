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

    name_index = []
    with open(names_file, 'r') as f:
        for l in f:
            line = l.rstrip('\n')
            name = line.split('\t')[0]
            name_index.append(index[name])
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

    variants_af = {}

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
                    variants_af[int(items[1])-1] = float(count)/num_seq

                else:
                    for i in names_index:
                        genotype = items[i].split(':')[0]
                        if genotype == '0|1' or genotype == '1|0' or genotype == '0/1' or genotype == '1/0':
                            count += 1
                        if genotype == '1|1' or genotype == '1/1':
                            count += 2
                    variants_af[int(items[1])-1] = float(count) / (len(names_index)*2)
    return variants_af


def compute_pi(bed_file, variants_af, num_seq):
    """
    This function computes pi in each nonoverlapping bed

    :param bed_file: a BED file specifying the coordinates for each nonoverlapping bed
    :param variants_af: a dictionary where key is the position of the variant and value
    is the allele frequency of that variant. This is the output from the function compute_af.
    :param num_seq: number of sequences in the VCF file
    :return: a dictionary where key is the (start, end) coordinates and value is pi for that bed.
    """
    beds_pi = {}
    with open(bed_file, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            line = line.split('\t')
            start = int(line[1])
            end = int(line[2])
            total_pi = 0
            for variant in variants_af.keys():
                if start <= variant and variant < end:
                    af = variants_af[variant]
                    pi = 2* af * (1-af)
                    total_pi += pi
            total_pi_adjusted = (num_seq/(num_seq-1))*total_pi
            beds_pi[(start, end)] = total_pi_adjusted


    return beds_pi

