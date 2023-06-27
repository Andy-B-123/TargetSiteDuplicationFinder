### Written by Andy Bachler, uploaded 16/06/2023
### Please ensure that the .bam file to be processed is <500MB, significant speed decreases occur for larger .bam files
### Consider subsetting your bam file to your region of interesting using samtools:
###	samtools view -b -o subset.bam original.bam contig:start-end
###	samtools index subset.bam
### Candidate TSD sites are initially identified using cluster_identifier from Scramble (https://github.com/GeneDx/scramble). Please run this on the .bam file you are looking to process (eg your subset file above)
### 	cluster_identifier subset.bam > subset.potential_sites.out
### Run this script:
### python TSD_parser.py --input_bam subset.bam --input_cluster_identifier_file subset.potential_sites.out --output_base subset.TSD_sites.bed --input_type SR

import pysam
from icecream import ic
import statistics
import csv
import pandas as pd
from tqdm import tqdm
import argparse
import numpy as np
import shutil
import subprocess

parser = argparse.ArgumentParser(description='Process BAM files and identify potential Target Site Duplication events.')
parser.add_argument('--input_bam', help='Full path to input .BAM file to process. Assumed sorted and indexed.')
parser.add_argument('--output_base', help='Output base name. Default current working directory and "TSD_parse_out"', default = "TSD_parse_out" )
parser.add_argument('--input_type', help = 'One of SR, LR, AS for Short-Read, Long-Read or Assembly for the type of bam to process. Default is Short-read.', default = "SR")
parser.add_argument('--cluster_identifer_path', help='The path to the cluster_identifier executable if not availbale in current PATH. Assumed to be in current PATH.', default = "PATH")

args = parser.parse_args()

input_path_bam = args.input_bam
output_base = args.output_base
input_bam_type = args.input_type
cluster_identifer_path = args.cluster_identifer_path
ic(input_path_bam)
ic(output_base)
ic(input_bam_type)
ic(cluster_identifer_path)

def process_cluster_identifier_data(cluster_identifier_path):

	candidate_TSD_sites = []

	tsv_file = open(cluster_identifier_path)
	ic(tsv_file)
	tsv_reader = csv.reader(tsv_file, delimiter="\t")
	ic(tsv_reader)

	previous_contig = False
	previous_coord = False
	previous_orientation = False
	previous_anchor_seq = False


	for line in tsv_reader:
		current_contig = line[0].split(":")[0]
		current_coord = int(line[0].split(":")[1])
		current_orientation = line[1]
		current_anchor_seq = line[4]

		if previous_contig == False:
			previous_contig = current_contig
			previous_coord = current_coord
			previous_orientation = current_orientation
			previous_anchor_seq = current_anchor_seq
			continue
		elif ((current_coord - previous_coord) < (window_size*2)) and (previous_orientation != current_orientation):
			if check_for_ref_seq_overlap_from_cluster_consensus(previous_anchor_seq,current_anchor_seq):
				candidate_TSD_sites.append([current_contig,previous_coord,current_coord])
			previous_contig = current_contig
			previous_coord = current_coord
			previous_orientation = current_orientation
			previous_anchor_seq = current_anchor_seq
		else:
			previous_contig = current_contig
			previous_coord = current_coord
			previous_orientation = current_orientation
			previous_anchor_seq = current_anchor_seq
			continue
	return candidate_TSD_sites

def check_for_ref_seq_overlap_from_cluster_consensus(anchor_seq_first, anchor_seq_second):
	for i in range(5,40):
		if anchor_seq_first[:i] == anchor_seq_second[-i:]:
			return True
	return False

def get_alignment_characteristics(bam_file):
	# Initilise coordinates as empty
	alignment_contigs_and_lengths = []

	# Open bam and read in header information to get names and end lengths of contigs
	input_bam = pysam.AlignmentFile(bam_file, "rb")
	header = input_bam.header

	contigs = header.get("SQ", [])
	for contig in contigs:
			contig_name = contig.get("SN")
			contig_length = contig.get("LN")
			alignment_contigs_and_lengths.append([contig_name,contig_length])
	input_bam.close()
	# Return the header info in a nice list
	return alignment_contigs_and_lengths

def generate_coordinates(bam_characteristics,cluster_identifier_hits,step_size,window_size):
	# Take the hits from cluster_identifier and generate windows around them for finer processing
	# ic(bam_characteristics)
	# ic(cluster_identifier_hits[0])

	coordinates = []
	for potential_TSD in cluster_identifier_hits:
		for j in range( (potential_TSD[1]-(window_size*2)) , (potential_TSD[1]+(window_size*2)) , step_size):
			coordinates.append([potential_TSD[0],j,j + window_size, potential_TSD[1]])
	# ic(coordinates)
	return coordinates

def read_BAM_file(bam_file, bam_file_coordinates):
	# Open the BAM file
	bam = pysam.AlignmentFile(bam_file, "rb")

	window_evaluation = []	
	ic(len(bam_file_coordinates))
	# ic(bam_file_coordinates)
	# iterate over contigs with window approach, for contigs that have reads
	for potential_site in tqdm(bam_file_coordinates):
		if check_for_TSD(bam, potential_site[0], potential_site[1],potential_site[2]) != None:
			TSD_scaffold, TSD_start, TSD_stop = check_for_TSD(bam, potential_site[0], potential_site[1],potential_site[2])
			if check_sufficient_coverage(bam, TSD_scaffold, TSD_start, TSD_stop):
				if check_for_clipped_bases(bam,TSD_scaffold, TSD_start, TSD_stop):
					if check_for_mismapped_reads(bam,TSD_scaffold, TSD_start, TSD_stop):
						window_evaluation.append([TSD_scaffold, TSD_start, TSD_stop, potential_site[0] +":"+ str(potential_site[3])])
		else:
			continue
	bam.close()
	# ic(window_evaluation)
	window_evaluation_deduped = []
	[window_evaluation_deduped.append(item) for item in window_evaluation if item not in window_evaluation_deduped]
	# ic(window_evaluation_deduped)
	return window_evaluation_deduped


def check_sufficient_coverage(bam_iterator,scaffold,start,stop):
	coverage_lists = bam_iterator.count_coverage(contig=scaffold, start=start, stop=stop, read_callback="nofilter")
	total_coverage = np.sum(coverage_lists)
	if total_coverage > coverage_threshold:
		return True
	else:
		return False

def check_for_TSD(bam_iterator,scaffold,start,stop):
	pileup_window = []
	bumps_identified = 0
	bump_start_coordinate = 0
	bump_end_coordinate = 0

	for pileupcolumn in bam_iterator.pileup(scaffold,start,stop,stepper = "nofilter",truncate =True):
		pileup_window.append(len(list(filter(None,pileupcolumn.get_query_sequences()))))
	
	#Filter out windows with deletions present which mess up assessment of TSD:
	if len(pileup_window) == 0:
		return None
	count_under_threshold = sum(1 for value in pileup_window if value < coverage_threshold)
	if count_under_threshold / len(pileup_window) > 0.15:
		return None
	
	average_coverage = np.mean(pileup_window)
	
	stdev_threshold = 5
	stdev_coverage = np.std(pileup_window)

	start_bump = False

	bump_difference_threshold = 0.25

	for i in range(len(pileup_window)):
		if i > 1:
			current_position_normalisedCov = pileup_window[i]/average_coverage
			previous_position_normalisedCov = pileup_window[i-1]/ average_coverage
			if (start_bump == True ) and abs(current_position_normalisedCov - previous_position_normalisedCov) > bump_difference_threshold:
				start_bump = False
				bump_end_coordinate = i
				bumps_identified += 1
			elif (current_position_normalisedCov - previous_position_normalisedCov) > bump_difference_threshold :
				start_bump = True
				bump_start_coordinate = i
				bumps_identified += 1
			else:
				continue
		else:
			continue
	if (bumps_identified == 2):
		return scaffold, (start+bump_start_coordinate), (start+bump_end_coordinate)
	else:
		return None

def check_for_clipped_bases(bam_iterator,scaffold,start,stop):
	total_reads = 0
	number_clipped_bases = 0
	
	for read in bam_iterator.fetch(scaffold,start,stop):
		total_reads = total_reads + 1
		if read.cigartuples is not None:
			for operation, length in read.cigartuples:
				if operation == 4:
					number_clipped_bases = number_clipped_bases + length
				elif operation == 5:
					number_clipped_bases = number_clipped_bases + length
	if number_clipped_bases == 0:
		return False
	elif total_reads == 0:
		return False
	elif (number_clipped_bases / (total_reads * input_read_length) ) > 0.25:
		return True
	else:
		return False

def check_for_mismapped_reads(bam_iterator,scaffold,start,stop):
	total_reads = 0
	mismapped_count = 0
	primary_count = 0
	non_primary_count = 0
	for read in bam_iterator.fetch(scaffold,start,stop):
		total_reads += 1
		if read.reference_name != read.next_reference_name:
			mismapped_count +=  1
		elif (read.reference_name == read.next_reference_name) and (read.template_length > 500 or read.template_length < 0):
			mismapped_count +=  1
		elif (read.is_secondary or read.is_supplementary):
			non_primary_count += 1
		elif read.mapping_quality < 40:
			non_primary_count += 1
		else:
			primary_count +=  1
	if ((mismapped_count / total_reads) > 0.75):
		return True
	else:
		return False

def main():
	# Example usage
	bam_file = input_path_bam
	global window_size 
	window_size= 30
	step_size = 6
	global coverage_threshold 
	coverage_threshold= 15
	global input_read_length 
	input_read_length = 150		 # Assumed to be 150

	# Check for cluster_identifier presence and establish path
	custom_path = False
	if cluster_identifer_path == "PATH":
		print(" 'cluster_identifer' from SCRAMBLE assumed to be in PATH, checking...")
		path = shutil.which("cluster_identifier") 
		if path is None:
			print("No 'cluster_identifier' executable present in path. Exiting...")
			exit()
		else:
			print("'cluster_identifier' present. Proceeding...")
	else:
		print("Custom path for 'cluster_identifer' from SCRAMBLE provided, checking...")
		path = shutil.which(cluster_identifer_path)
		if path is None:
			print("No 'cluster_identifier' executable present in path. Exiting...")
			exit()
		else:
			print("'cluster_identifier' present. Proceeding...")
			custom_path = True

	if custom_path == False:
		output_cluster_file_path = bam_file.replace(".bam",'.cluster_identifier.out')
		cluster_identifer_command = "cluster_identifier " + bam_file + " > " + output_cluster_file_path
		process = subprocess.Popen(cluster_identifer_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
	else:
		output_cluster_file_path = bam_file.replace(".bam",'.cluster_identifier.out')
		cluster_identifer_command = cluster_identifer_path + " " + bam_file + " > " + output_cluster_file_path
		process = subprocess.Popen(cluster_identifer_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)

	bam_file_parameters = get_alignment_characteristics(bam_file)
	ic(len(bam_file_parameters))
	

	if (input_bam_type == "SR"):
		ic("Short read sample indicated")
		cluster_identifier_file_TSD_coordinates = process_cluster_identifier_data(output_cluster_file_path)
		ic(len(cluster_identifier_file_TSD_coordinates))
		bam_file_coordinates = generate_coordinates(bam_file_parameters,cluster_identifier_file_TSD_coordinates,step_size,window_size)
		ic(len(bam_file_coordinates))
		Potential_TSD_list = read_BAM_file(bam_file, bam_file_coordinates)
		pd.DataFrame(Potential_TSD_list).to_csv(output_base + ".WindowSize_" + str(window_size) + ".bed", index=False, sep = "\t", header = False)
	else:
		print("Working on long-read and assembly alignments at the moment! Please stay tuned :-)")
		exit()

if __name__ == "__main__":
	main()
