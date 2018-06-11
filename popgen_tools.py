from helper_funcs import *
import argparse


def parse_args():
	"""
	Parsing command-line arguments
	"""
	parser = argparse.ArgumentParser(description='This script generates a site frequency spectrum and calculates '
						'genetic diversity, pi in nonoverlapping windows (average number of '
						'DNA differences between all pairs of sequence).')

	parser.add_argument('--vcf_file', required=True,
						help='REQUIRED. Input the path to a VCF file')

	parser.add_argument('--num_chr', required=True,
						help='REQUIRED. Input the number of chromosomes in your VCF file. The number of alleles '
							 'is equal to the number of individuals multiplying by 2.')

	parser.add_argument('--window_bed', required=False,
						help='Input the path to the BED file that specifies the coordinates for each '
							 'nonoverlapping window. This argument is required for the analysis of '
							 'genetic diversity.')

	parser.add_argument('--sfs_out', required=False,
						help='Input the path for the output file for the site frequency spectrum. '
							 'If this parameter is not specified, an output file called sfs.out will be output '
							 'in the current directory.')

	parser.add_argument('--pi_out', required=False,
						help='Input the path for the output file for pi. If this parameter is not specified, '
							 'an output file called pi.out will be output in the current directory.')

	# Options for turning on/off parts of the pipeline
	parser.add_argument('--no_sfs', action='store_true', default=False, help='Turn on this flag if you do not want '
																			 'to generate a site frequency spectrum.')

	parser.add_argument('--no_pi', action='store_true', default=False, help='Turn on this flag if you do not want '
																			'to calculate genetic diversity.')

	args = parser.parse_args()
	return args


def main():

	args = parse_args()

	if args.no_sfs is not True:
		# Calculate the number of alternate alleles for each variant. This is a list where each item in the list is
		# the count of the alternate allele for each variant.
		alt_allele_count = count_alt_allele(args.vcf_file)

		# Generate a folded site frequency spectrum
		sfs = make_sfs(int(args.num_chr), alt_allele_count)

		if args.sfs_out:
			outfile = open(args.sfs_out, 'w')
			header = ['bin', 'count']
			outfile.write('\t'.join(header) + '\n')

			for k, v in sfs.items():
				out = [str(k), str(v)]
				outfile.write('\t'.join(out) + '\n')

			outfile.close()

		else:
			outfile = open('sfs.out', 'w')
			header = ['bin', 'count']
			outfile.write('\t'.join(header) + '\n')

			for k, v in sfs.items():
				out = [str(k), str(v)]
				outfile.write('\t'.join(out) + '\n')

			outfile.close()

	if args.no_pi is not True:
		# Calculate genetic diversity, pi which is the average number of DNA differences across all pairs of sequences
		variants_af = compute_af(args.vcf_file, float(args.num_chr))
		total_pi_adjusted = compute_pi(args.window_bed, variants_af, float(args.num_chr))

		if args.pi_out:
			outfile = open(args.pi_out, 'w')
			header = ['start', 'end', 'pi', 'pi_per_site']
			outfile.write('\t'.join(header) + '\n')

			for k, v in total_pi_adjusted.items():
				out = [str(k[0]), str(k[1]), str(v), str(float(v)/(k[1]-k[0]))]
				outfile.write('\t'.join(out) + '\n')

			outfile.close()

		else:
			outfile = open('pi.out', 'w')
			header = ['start', 'end', 'pi', 'pi_per_site']
			outfile.write('\t'.join(header) + '\n')

			for k, v in total_pi_adjusted.items():
				out = [str(k[0]), str(k[1]), str(v), str(float(v) / (k[1] - k[0]))]
				outfile.write('\t'.join(out) + '\n')

			outfile.close()


if __name__ == '__main__':
	main()
