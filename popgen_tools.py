from helper_funcs import *
import argparse

#TODO: Add error checking for inputs
#TODO: When generating an SFS, allow for the requirement of target bed to be optional
#TODO: Instead of inputting a list of individuals, the default should be to calculate pi and to generate SFS for all of the individuals in the vcf file

def parse_args():
	"""
	Parsing command-line arguments
	"""
	parser = argparse.ArgumentParser(description='This script generates a site frequency '
						'spectrum and calculates genetic diversity, pi. Pi is defined as '
						'the average number of DNA differences between all pairs of sequence).')

	parser.add_argument('--vcf_file', required=True,
						help='REQUIRED. Input the path to a VCF file. Either gzipped or '
							 'not gzipped file works.')

	parser.add_argument('--target_bed', required=True,
						help='REQUIRED. Input the path to the BED file that specifies the coordinates '
							 'to generate SFS or to calculate pi. For example, this file '
							 'specificies the coordinates for the putatively neutral regions')

	parser.add_argument('--names_list', required=True,
						help='REQUIRED. Input the path to the file that lists the individuals from '
							 'the VCF file that you want to calculate genetic diversity '
							 'or to generate the SFS.')

	parser.add_argument('--sfs_out', required=False,
						help='Input the path for the output file for the site frequency spectrum. '
							 'If this parameter is not specified, an output file called '
							 'sfs.out will be outputted in the current directory.')

	parser.add_argument('--pi_out', required=False,
						help='Input the path for the output file for pi. If this parameter '
							 'is not specified, an output file called pi.out will be outputted '
							 'in the current directory.')

	# Options for turning on/off parts of the pipeline
	parser.add_argument('--no_sfs', action='store_true', default=False,
						help='Turn on this flag if you do not want '
							 'to generate a site frequency spectrum.')

	parser.add_argument('--no_pi', action='store_true', default=False,
						help='Turn on this flag if you do not want '
							 'to calculate genetic diversity.')

	args = parser.parse_args()
	return args


def main():

	args = parse_args()

	if args.no_sfs is not True:

		# Find the index of the individuals whose genetic diversity you want to calculate
		names_index = find_index(args.vcf_file, args.names_list)

		# Calculate the number of alternate alleles for each variant.
		alt_allele_count = count_alt_allele(args.vcf_file, names_index, args.target_bed)

		# Generate a folded site frequency spectrum
		sfs = make_sfs(len(names_index)*2, alt_allele_count)

		if args.sfs_out:
			outfile = open(args.sfs_out, 'w')
			header = ['frequency', 'num_variants']
			outfile.write('\t'.join(header) + '\n')

			for k, v in sfs.items():
				out = [str(k), str(v)]
				outfile.write('\t'.join(out) + '\n')

			outfile.close()

		else:
			outfile = open('sfs.out', 'w')
			header = ['frequency', 'num_variants']
			outfile.write('\t'.join(header) + '\n')

			for k, v in sfs.items():
				out = [str(k), str(v)]
				outfile.write('\t'.join(out) + '\n')

			outfile.close()

	if args.no_pi is not True:

		# Find the index of the individuals whose genetic diversity you want to calculate
		names_index = find_index(args.vcf_file, args.names_list)

		# Calculate genetic diversity, pi
		variants_af = compute_af(args.vcf_file, names_index)

		# Calculate adjusted pi
		total_pi_adjusted = compute_pi(args.target_bed, variants_af[0], variants_af[1], len(names_index)*2)

		# for k, v in total_pi_adjusted.items():
		# 	print (k[0], k[1], v)

		if args.pi_out:
			outfile = open(args.pi_out, 'w')
			header = ['start', 'end', 'pi']
			outfile.write('\t'.join(header) + '\n')

			for k, v in total_pi_adjusted.items():
				out = [str(k[0]), str(k[1]), str(v)]
				outfile.write('\t'.join(out) + '\n')

			outfile.close()

		else:
			outfile = open('pi.out', 'w')
			header = ['start', 'end', 'pi']
			outfile.write('\t'.join(header) + '\n')

			for k, v in total_pi_adjusted.items():
				out = [str(k[0]), str(k[1]), str(v)]
				outfile.write('\t'.join(out) + '\n')

			outfile.close()


if __name__ == '__main__':
	main()
