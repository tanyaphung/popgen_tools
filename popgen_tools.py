from helper_funcs import *
from calc_pi import *
import argparse
import pandas as pd

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
	parser.add_argument('--pi_all', action='store_true', default=False,
						help='Turn on this flag if you want to calculate pi using all of the sites in the VCF file.')

	parser.add_argument('--pi_target', action='store_true', default=False,
						help='Turn on this flag if you want to calculate pi in regions of the genome specified by a BED file.')

	parser.add_argument('--pi_window', action='store_true', default=False,
						help='Turn on this flag if you want to calculate pi in windows in the presence of target_bed file.')

	# Extra required files for different options of calculating pi
	parser.add_argument('--target_bed', required=False,
						help='If you want to calculate pi in regions of the genome specified by a BED file, input the '
							 'path here.')

	parser.add_argument('--window_bed', required=False,
						help='If you want to calculate genetic diversity in window in the presence of target_bed file, '
							 'input the path here.')

	# # Providing the output files
	# parser.add_argument('--sfs_out', required=False,
	# 					help='Input the path for the output file for the site frequency spectrum. '
	# 						 'If this parameter is not specified, an output file called '
	# 						 'sfs.out will be outputted in the current directory.')
	#
	# parser.add_argument('--pi_out', required=False,
	# 					help='Input the path for the output file for pi. If this parameter '
	# 						 'is not specified, an output file called pi.out will be outputted '
	# 						 'in the current directory.')
	#
	# parser.add_argument('--total_SNPs', required=False,
	# 					help='Input the path to the output file for the number of SNPs.') #TODO: make outputing the number of SNPs an option.
	#
	# # Options for turning on/off parts of the pipeline
	# parser.add_argument('--no_sfs', action='store_true', default=False,
	# 					help='Turn on this flag if you do not want '
	# 						 'to generate a site frequency spectrum.')
	#
	# parser.add_argument('--no_pi', action='store_true', default=False,
	# 					help='Turn on this flag if you do not want '
	# 						 'to calculate genetic diversity.')
	#
	# parser.add_argument('--window', action='store_true', default=False,
	# 					help='Turn on this flag if you want to calculate genetic diversity'
	# 						 ' in windows. The default when calculating pi is to calculate'
	# 						 ' pi in the target regions.')
	#
	# parser.add_argument('--sfs_no_target_bed', action='store_true', default=False,
	# 					help='Turn on this flag if you are making the SFS using all of the '
	# 						 'variants from the VCF file.')

	args = parser.parse_args()
	return args


def main():

	args = parse_args()

	# Check if the user wants to perform the calculations for all of the individuals in the VCF or just a subset of
	# individuals, which is specified by the presence of the argument names_list
	if args.names_list:
		# Find the index of the individuals whose genetic diversity you want to calculate
		names_index = find_index(args.vcf_file, args.names_list)

	else:
		names_index = find_index(args.vcf_file)

	# For pi, get the variants and afs
	variants_afs = compute_af(args.vcf_file, names_index)

	# Check if the user wants to calculate pi using all of the sites in the VCF file.
	if args.pi_all:
		# Calculate allele frequency
		afs = variants_afs[1]
		# Compute pi_all
		pi_all = compute_pi_all(afs, len(names_index)*2)

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
		adjusted_pi = compute_pi_target(targets, variants_afs[0], variants_afs[1], len(names_index) * 2)
		pi_outfile = open("pi_target.out", 'w')
		header = ['chr', 'start', 'end', 'pi']
		pi_outfile.write('\t'.join(header) + '\n')
		for k, v in adjusted_pi[0].items():
			out = [str(k[0]), str(k[1]), str(k[2]), str(v)]
			pi_outfile.write('\t'.join(out) + '\n')
		pi_outfile.close()

		# Calculate pi per site
		pi_target = pd.read_table("pi_target.out")
		pi_total = pi_target['pi'].sum()
		callable_sites = (pi_target['end'] - pi_target['start']).sum()
		pi_per_site = pi_total/callable_sites

		with open('pi_target_per_site.out', 'w') as out:
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

		# # Calculate the total callable site in each window
		# windows = []
		# with open(args.window_bed, 'r') as f:
		# 	for line in f:
		# 		_, s, e = line.split('\t')
		# 		windows.append((int(s), int(e)))
		#
		# targets = []
		# with open(args.target_bed, 'r') as f:
		# 	for line in f:
		# 		_, s, e = line.split('\t')
		# 		targets.append((int(s), int(e)))
		#
		# total_callable = tabulate_callable_sites_each_window(windows, targets)
		#
		# pi_outfile = open(args.pi_out, 'w')
		#
		# for k, v in adjusted_pi[0].items():
		# 	if total_callable[k] != 0:
		# 		out = [str(k[0]), str(k[1]), str(total_callable[k]), str(v)]
		# 		pi_outfile.write('\t'.join(out) + '\n')
		#
		# pi_outfile.close()


# if args.no_sfs is not True:
	#
	# 	if args.sfs_no_target_bed: #when this flag is turned on, use all of the variants from the
	# 	#  VCF file
	# 	# Calculate the number of alternate alleles for each variant.
	# 		alt_allele_count_no_target_bed = count_alt_allele_no_target_bed(args.vcf_file, names_index)
	#
	# 		# Generate a folded site frequency spectrum
	# 		sfs = make_sfs(len(names_index) * 2, alt_allele_count_no_target_bed)
	#
	# 		outfile = open(args.sfs_out, 'w')
	# 		header = ['frequency', 'num_variants']
	# 		outfile.write('\t'.join(header) + '\n')
	#
	# 		for k, v in sfs.items():
	# 			out = [str(k), str(v)]
	# 			outfile.write('\t'.join(out) + '\n')
	#
	# 		outfile.close()
	#
	# 	else:
	#
	# 		# Calculate the number of alternate alleles for each variant.
	# 		alt_allele_count = count_alt_allele(args.vcf_file, names_index, args.target_bed)
	#
	# 		# Generate a folded site frequency spectrum
	# 		sfs = make_sfs(len(names_index)*2, alt_allele_count)
	#
	# 		outfile = open(args.sfs_out, 'w')
	# 		header = ['frequency', 'num_variants']
	# 		outfile.write('\t'.join(header) + '\n')
	#
	# 		for k, v in sfs.items():
	# 			out = [str(k), str(v)]
	# 			outfile.write('\t'.join(out) + '\n')
	#
	# 		outfile.close()

	# if args.no_pi is not True:
	#
	# 	# If the --no_pi flag is not turned on, then the user wants to calculate genetic diversity
	#
	# 	if args.window is not True: #calculate pi not in non-overlapping windows
	#
	# 		# Calculate allele frequency
	# 		variants_af = compute_af(args.vcf_file, names_index)
	#
	# 		# Calculate adjusted pi
	# 		adjusted_pi = compute_pi(args.target_bed, variants_af[0], variants_af[1], len(names_index)*2)
	#
	# 		# Save output for genetic diversity in a file
	# 		pi_outfile = open(args.pi_out, 'w')
	# 		for k, v in adjusted_pi[0].items():
	# 			out = [str(k[0]), str(k[1]), str(v)]
	# 			pi_outfile.write('\t'.join(out) + '\n')
	# 		pi_outfile.close()
	#
	# 		# Save output for the number of SNPs in a file
	# 		total_SNPs_outfile = open(args.total_SNPs, 'w')
	# 		for k, v in adjusted_pi[1].items():
	# 			out = [str(k[0]), str(k[1]), str(v)]
	# 			total_SNPs_outfile.write('\t'.join(out) + '\n')
	# 		total_SNPs_outfile.close()
	#
	# 	else:
	#
	# 		# Calculate allele frequency
	# 		variants_af = compute_af(args.vcf_file, names_index)
	#
	# 		# Calculate adjusted pi
	# 		adjusted_pi = compute_pi(args.window_bed, variants_af[0], variants_af[1], len(names_index)*2)
	#
	# 		# Calculate the total callable site in each window
	# 		windows = []
	# 		with open(args.window_bed, 'r') as f:
	# 			for line in f:
	# 				_, s, e = line.split('\t')
	# 				windows.append((int(s), int(e)))
	#
	# 		targets = []
	# 		with open(args.target_bed, 'r') as f:
	# 			for line in f:
	# 				_, s, e = line.split('\t')
	# 				targets.append((int(s), int(e)))
	#
	# 		total_callable = tabulate_callable_sites_each_window(windows, targets)
	#
	# 		pi_outfile = open(args.pi_out, 'w')
	#
	# 		for k, v in adjusted_pi[0].items():
	# 			if total_callable[k] != 0:
	# 				out = [str(k[0]), str(k[1]), str(total_callable[k]), str(v)]
	# 				pi_outfile.write('\t'.join(out) + '\n')
	#
	# 		pi_outfile.close()


if __name__ == '__main__':
	main()
