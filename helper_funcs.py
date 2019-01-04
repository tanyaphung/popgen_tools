import collections
import gzip


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
    index = {}

    file_test_result = file_test(vcf_file)
    open_func = file_test_result[0]
    mode = file_test_result[1]
    with open_func(vcf_file, mode) as f:
        for l in f:
            line = l.rstrip('\n')
            if line.startswith('#CHROM'):
                items = line.split('\t')
                if names_list is None:
                    for i in range(9, len(items)): #TODO: have to check if the genotype actually starts at position 9.
                        name_index.append(i)
                    break
                else:
                    for i in range(9, len(items)):
                        index[items[i]] = i
                    with open(names_list, 'r') as f:
                        for l in f:
                            line = l.rstrip('\n')
                            name = line.split('\t')[0]
                            name_index.append(index[name])
                    break
    return name_index

def count_alt_allele(vcf_file, names_index, bed_file):
    """
    This function counts the number of alternate allele for each variant
    :param vcf_file:
    :param names_index:
    :param bed_file:
    :return: a list where each item in the list is the count of the alternate allele.
    """
    # Generate a set of positions of interest (defined by the bed file)
    sites_set = set()
    with open(bed_file, 'r') as f:
        for l in f:
            line = l.rstrip('\n')
            items = line.split('\t')
            for i in range(int(items[1]), int(items[2])):
                sites_set.add(i)

    alt_allele_count = []
    open_func = file_test(vcf_file)
    with open_func(vcf_file, 'r') as f:
        for l in f:
            line = l.rstrip('\n')
            if not line.startswith('#'):
                items = line.split('\t')
                if int(items[1])-1 in sites_set: #subtract 1 because VCF file is 1-based whereas BED file is 0-based.
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
    # print (snps)
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
    # print (len(alt_allele_count))
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



