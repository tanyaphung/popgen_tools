from helper_funcs import *
import argparse

#TODO: Add error checking for inputs
#TODO: When generating an SFS, allow for the requirement of target bed to be optional

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

	parser.add_argument('--names_list', required=False,
						help='Input the path to the file that lists the individuals from '
							 'the VCF file that you want to calculate genetic diversity '
							 'or to generate the SFS.')

	parser.add_argument('--target_bed', required=False,
						help='Input the path to the BED file that specifies the coordinates '
							 'to generate SFS or to calculate pi. For example, this file '
							 'specificies the coordinates for the putatively neutral regions.')

	parser.add_argument('--window_bed', required=False,
						help='Input the path the the BED file that specifies the coordinates '
							 'for nonoverlapping windows.')

	# Providing the output files
	parser.add_argument('--sfs_out', required=False,
						help='Input the path for the output file for the site frequency spectrum. '
							 'If this parameter is not specified, an output file called '
							 'sfs.out will be outputted in the current directory.')

	parser.add_argument('--pi_out', required=False,
						help='Input the path for the output file for pi. If this parameter '
							 'is not specified, an output file called pi.out will be outputted '
							 'in the current directory.')

	parser.add_argument('--total_SNPs', required=False,
						help='Input the path to the output file for the number of SNPs.') #TODO: make outputing the number of SNPs an option.

	# Options for turning on/off parts of the pipeline
	parser.add_argument('--no_sfs', action='store_true', default=False,
						help='Turn on this flag if you do not want '
							 'to generate a site frequency spectrum.')

	parser.add_argument('--no_pi', action='store_true', default=False,
						help='Turn on this flag if you do not want '
							 'to calculate genetic diversity.')

	parser.add_argument('--window', action='store_true', default=False,
						help='Turn on this flag if you want to calculate genetic diversity'
							 ' in windows. The default when calculating pi is to calculate'
							 ' pi in the target regions.')

	args = parser.parse_args()
	return args


def main():

	args = parse_args()

	if args.names_list:
		# Find the index of the individuals whose genetic diversity you want to calculate
		names_index = find_index(args.vcf_file, args.names_list)

	else:
		names_index = find_index_all(args.vcf_file)

	if args.no_sfs is not True:

		# Calculate the number of alternate alleles for each variant.
		alt_allele_count = count_alt_allele(args.vcf_file, names_index, args.target_bed)

		# Generate a folded site frequency spectrum
		sfs = make_sfs(len(names_index)*2, alt_allele_count)

		outfile = open(args.sfs_out, 'w')
		header = ['frequency', 'num_variants']
		outfile.write('\t'.join(header) + '\n')

		for k, v in sfs.items():
			out = [str(k), str(v)]
			outfile.write('\t'.join(out) + '\n')

		outfile.close()

	if args.no_pi is not True:

		if args.window is not True: #calculate pi not in non-overlapping windows

			# Calculate allele frequency
			variants_af = compute_af(args.vcf_file, names_index)

			# Calculate adjusted pi
			adjusted_pi = compute_pi(args.target_bed, variants_af[0], variants_af[1], len(names_index)*2)

			# Save output for genetic diversity in a file
			pi_outfile = open(args.pi_out, 'w')
			for k, v in adjusted_pi[0].items():
				out = [str(k[0]), str(k[1]), str(v)]
				pi_outfile.write('\t'.join(out) + '\n')
			pi_outfile.close()

			# Save output for the number of SNPs in a file
			total_SNPs_outfile = open(args.total_SNPs, 'w')
			for k, v in adjusted_pi[1].items():
				out = [str(k[0]), str(k[1]), str(v)]
				total_SNPs_outfile.write('\t'.join(out) + '\n')
			total_SNPs_outfile.close()

		else:

			# Calculate allele frequency
			variants_af = compute_af(args.vcf_file, names_index)

			# Calculate adjusted pi
			adjusted_pi = compute_pi(args.window_bed, variants_af[0], variants_af[1], len(names_index)*2)

			# Calculate the total callable site in each window
			windows = []
			with open(args.window_bed, 'r') as f:
				for line in f:
					_, s, e = line.split('\t')
					windows.append((int(s), int(e)))

			targets = []
			with open(args.target_bed, 'r') as f:
				for line in f:
					_, s, e = line.split('\t')
					targets.append((int(s), int(e)))

			total_callable = tabulate_callable_sites_each_window(windows, targets)

			pi_outfile = open(args.pi_out, 'w')

			for k, v in adjusted_pi[0].items():
				if total_callable[k] != 0:
					out = [str(k[0]), str(k[1]), str(total_callable[k]), str(v)]
					pi_outfile.write('\t'.join(out) + '\n')

			pi_outfile.close()


if __name__ == '__main__':
	main()
